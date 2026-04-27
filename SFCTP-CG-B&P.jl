# Stochastic fixed-charge transportation problem (SFCTP) with Branch-and-Price

using JuMP
using HiGHS
using Printf

const NSD_COST = 100.0
const COLUMN_TOLERANCE = 1e-6
const INTEGRALITY_TOLERANCE = 1e-6

struct ScenarioColumn
    y::Matrix{Float64}
    x::Matrix{Float64}
    z::Vector{Float64}
    cost::Float64
end

struct BranchNode
    id::Int
    depth::Int
    fixed_zero::Set{Tuple{Int,Int}}
    fixed_one::Set{Tuple{Int,Int}}
    columns_by_scenario::Vector{Vector{ScenarioColumn}}
end

function create_extensive_form_model(I, J, S, A, B, C, F, P)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, 0 <= X[s = 1:S, i = 1:I, j = 1:J] <= min(A[i], B[s, j]))
    @variable(model, Y[i = 1:I, j = 1:J], Bin)
    @variable(model, Z[s = 1:S, j = 1:J] >= 0)

    @objective(
        model,
        Min,
        sum(F[i, j] * Y[i, j] for i = 1:I, j = 1:J) +
        sum(C[i, j] * P[s] * X[s, i, j] for s = 1:S, i = 1:I, j = 1:J) +
        sum(NSD_COST * P[s] * Z[s, j] for s = 1:S, j = 1:J)
    )

    @constraint(model, [s = 1:S, i = 1:I], sum(X[s, i, j] for j = 1:J) <= A[i])
    @constraint(model, [s = 1:S, j = 1:J], sum(X[s, i, j] for i = 1:I) >= B[s, j] - Z[s, j])
    @constraint(
        model,
        [s = 1:S, i = 1:I, j = 1:J],
        X[s, i, j] <= min(A[i], B[s, j]) * Y[i, j]
    )

    return model
end

function build_artificial_column(I, J, demand)
    return ScenarioColumn(
        zeros(I, J),
        zeros(I, J),
        Float64.(vec(demand)),
        NSD_COST * sum(demand),
    )
end

function clone_columns(columns_by_scenario)
    return [copy(columns) for columns in columns_by_scenario]
end

function column_already_exists(columns, candidate; atol = 1e-6)
    for column in columns
        if isapprox(column.cost, candidate.cost; atol = atol) &&
           all(isapprox.(column.y, candidate.y; atol = atol))
            return true
        end
    end
    return false
end

function build_restricted_master(I, J, S, F, P, columns_by_scenario, fixed_zero, fixed_one)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, 0 <= Y[1:I, 1:J] <= 1)
    for (i, j) in fixed_zero
        fix(Y[i, j], 0.0; force = true)
    end
    for (i, j) in fixed_one
        fix(Y[i, j], 1.0; force = true)
    end

    lambda = [VariableRef[] for _ = 1:S]
    for s = 1:S
        for _ in eachindex(columns_by_scenario[s])
            push!(lambda[s], @variable(model, lower_bound = 0))
        end
    end

    @objective(
        model,
        Min,
        sum(F[i, j] * Y[i, j] for i = 1:I, j = 1:J) + sum(
            P[s] * columns_by_scenario[s][k].cost * lambda[s][k] for s = 1:S,
            k in eachindex(columns_by_scenario[s])
        )
    )

    convexity = @constraint(model, [s = 1:S], sum(lambda[s]) == 1)
    link = @constraint(
        model,
        [s = 1:S, i = 1:I, j = 1:J],
        Y[i, j] >= sum(
            columns_by_scenario[s][k].y[i, j] * lambda[s][k] for
            k in eachindex(columns_by_scenario[s])
        )
    )

    return (; model, Y, lambda, convexity, link)
end

function solve_pricing_problem(
    s,
    I,
    J,
    A,
    B,
    C,
    P,
    convexity_dual,
    linking_dual,
    fixed_zero,
)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, 0 <= X[i = 1:I, j = 1:J] <= min(A[i], B[s, j]))
    @variable(model, Y[i = 1:I, j = 1:J], Bin)
    @variable(model, Z[j = 1:J] >= 0)

    for (i, j) in fixed_zero
        fix(Y[i, j], 0.0; force = true)
    end

    @expression(
        model,
        operating_cost,
        sum(C[i, j] * X[i, j] for i = 1:I, j = 1:J) + sum(NSD_COST * Z[j] for j = 1:J)
    )

    @objective(
        model,
        Min,
        P[s] * operating_cost + sum(linking_dual[i, j] * Y[i, j] for i = 1:I, j = 1:J) -
        convexity_dual
    )

    @constraint(model, [i = 1:I], sum(X[i, j] for j = 1:J) <= A[i])
    @constraint(model, [j = 1:J], sum(X[i, j] for i = 1:I) >= B[s, j] - Z[j])
    @constraint(model, [i = 1:I, j = 1:J], X[i, j] <= min(A[i], B[s, j]) * Y[i, j])

    optimize!(model)
    @assert is_solved_and_feasible(model)

    reduced_cost = objective_value(model)
    if reduced_cost < -COLUMN_TOLERANCE
        return ScenarioColumn(value.(Y), value.(X), value.(Z), value(operating_cost)),
        reduced_cost
    end

    return nothing, reduced_cost
end

function is_integral_solution(Y)
    return all(abs(value - round(value)) <= INTEGRALITY_TOLERANCE for value in Y)
end

function select_branch_arc(Y)
    best_arc = nothing
    best_score = Inf
    for i in axes(Y, 1), j in axes(Y, 2)
        value = Y[i, j]
        if INTEGRALITY_TOLERANCE < value < 1.0 - INTEGRALITY_TOLERANCE
            score = abs(value - 0.5)
            if score < best_score
                best_score = score
                best_arc = (i, j, value)
            end
        end
    end
    return best_arc
end

function print_node_header(node)
    println(
        "Node ",
        lpad(node.id, 4),
        " | depth = ",
        node.depth,
        " | fixed_zero = ",
        length(node.fixed_zero),
        " | fixed_one = ",
        length(node.fixed_one),
    )
    return
end

function print_node_iteration(
    node_id,
    iteration,
    node_bound,
    best_reduced_cost,
    total_columns,
    new_columns,
)
    println(
        "  node ",
        lpad(node_id, 4),
        " iter ",
        lpad(iteration, 3),
        " | LP = ",
        @sprintf("%10.4f", node_bound),
        " | best rc = ",
        @sprintf("%10.4e", best_reduced_cost),
        " | columns = ",
        lpad(total_columns, 4),
        " | new = ",
        lpad(new_columns, 2),
    )
    return
end

function solve_node_relaxation!(node, I, J, S, A, B, C, F, P; max_iterations = 1_000)
    columns_by_scenario = clone_columns(node.columns_by_scenario)
    iteration_history = NamedTuple[]
    master = nothing

    for iteration = 1:max_iterations
        master = build_restricted_master(
            I,
            J,
            S,
            F,
            P,
            columns_by_scenario,
            node.fixed_zero,
            node.fixed_one,
        )
        optimize!(master.model)
        @assert is_solved_and_feasible(master.model; dual = true)

        node_bound = objective_value(master.model)
        best_reduced_cost = Inf
        new_columns = 0

        for s = 1:S
            column, reduced_cost = solve_pricing_problem(
                s,
                I,
                J,
                A,
                B,
                C,
                P,
                dual(master.convexity[s]),
                dual.(master.link[s, :, :]),
                node.fixed_zero,
            )
            best_reduced_cost = min(best_reduced_cost, reduced_cost)
            if column !== nothing && !column_already_exists(columns_by_scenario[s], column)
                push!(columns_by_scenario[s], column)
                new_columns += 1
            end
        end

        total_columns = sum(length(columns) for columns in columns_by_scenario)
        push!(
            iteration_history,
            (
                iteration = iteration,
                objective = node_bound,
                best_reduced_cost = best_reduced_cost,
                columns = total_columns,
                new_columns = new_columns,
            ),
        )
        print_node_iteration(
            node.id,
            iteration,
            node_bound,
            best_reduced_cost,
            total_columns,
            new_columns,
        )

        if new_columns == 0
            Y_value = value.(master.Y)
            return (
                bound = node_bound,
                Y = Y_value,
                master = master,
                columns_by_scenario = columns_by_scenario,
                iteration_history = iteration_history,
            )
        end
    end

    error(
        "Column generation did not converge at node $(node.id) within the iteration limit.",
    )
end

function branch_and_price!(
    I,
    J,
    S,
    A,
    B,
    C,
    F,
    P;
    max_nodes = 10_000,
    max_iterations_per_node = 1_000,
)
    for s = 1:S
        if sum(B[s, :]) >= sum(A)
            @warn "The total demand for scenario $s exceeds the total supply. Non-supplied demand will absorb the deficit."
        end
    end

    root_columns = [[build_artificial_column(I, J, B[s, :])] for s = 1:S]
    open_nodes =
        [BranchNode(1, 0, Set{Tuple{Int,Int}}(), Set{Tuple{Int,Int}}(), root_columns)]
    node_history = NamedTuple[]
    next_node_id = 2
    best_objective = Inf
    best_Y = nothing
    best_master = nothing
    processed_nodes = 0

    while !isempty(open_nodes)
        if processed_nodes >= max_nodes
            error(
                "Branch-and-price terminated because the node limit was reached before proving optimality.",
            )
        end

        node = pop!(open_nodes)
        processed_nodes += 1
        print_node_header(node)
        result = solve_node_relaxation!(
            node,
            I,
            J,
            S,
            A,
            B,
            C,
            F,
            P;
            max_iterations = max_iterations_per_node,
        )

        node_bound = result.bound
        push!(
            node_history,
            (
                node_id = node.id,
                depth = node.depth,
                bound = node_bound,
                fixed_zero = length(node.fixed_zero),
                fixed_one = length(node.fixed_one),
            ),
        )

        if node_bound >= best_objective - INTEGRALITY_TOLERANCE
            println("  pruned by bound")
            continue
        end

        if is_integral_solution(result.Y)
            best_objective = node_bound
            best_Y = round.(Int, result.Y)
            best_master = result.master
            println("  incumbent updated to $(round(best_objective; digits = 4))")
            continue
        end

        branch_arc = select_branch_arc(result.Y)
        @assert branch_arc !== nothing
        (i_branch, j_branch, value_branch) = branch_arc
        println(
            "  branching on Y[",
            i_branch,
            ", ",
            j_branch,
            "] = ",
            @sprintf("%.4f", value_branch),
        )

        child_columns = result.columns_by_scenario
        child_zero = BranchNode(
            next_node_id,
            node.depth + 1,
            union(node.fixed_zero, Set([(i_branch, j_branch)])),
            copy(node.fixed_one),
            clone_columns(child_columns),
        )
        next_node_id += 1

        child_one = BranchNode(
            next_node_id,
            node.depth + 1,
            copy(node.fixed_zero),
            union(node.fixed_one, Set([(i_branch, j_branch)])),
            clone_columns(child_columns),
        )
        next_node_id += 1

        push!(open_nodes, child_one)
        push!(open_nodes, child_zero)
    end

    if best_Y === nothing
        error("Branch-and-price finished without finding an integer incumbent.")
    end

    return (
        objective = best_objective,
        Y = best_Y,
        master = best_master,
        node_history = node_history,
        processed_nodes = processed_nodes,
    )
end

## Define data

### system dimensions
I = 4
J = 3
S = 3

### parameters
A = [20, 30, 40, 20]
P = [0.5, 0.3, 0.2]

B = [
    21 51 31
    32 22 52
    53 33 23
]

C = [
    1 2 3
    3 2 1
    2 3 4
    4 3 2
]

F = [
    10 20 30
    20 30 40
    30 40 50
    40 50 60
]

## Solve the problem with branch-and-price
bp_result = branch_and_price!(I, J, S, A, B, C, F, P)
@info "Branch-and-price objective: $(bp_result.objective)"
@show bp_result.Y
@info "Processed nodes: $(bp_result.processed_nodes)"

## Solve the extensive form for validation
extensive_form = create_extensive_form_model(I, J, S, A, B, C, F, P)
optimize!(extensive_form)
@assert is_solved_and_feasible(extensive_form)

extensive_objective = objective_value(extensive_form)
Y_extensive = round.(Int, value.(extensive_form[:Y]))
@info "Extensive-form objective: $extensive_objective"
@show Y_extensive

@assert isapprox(bp_result.objective, extensive_objective; atol = 1e-6)
@assert bp_result.Y == Y_extensive
