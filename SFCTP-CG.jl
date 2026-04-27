# Stochastic fixed-charge transportation problem (SFCTP) with Column Generation

using JuMP
using HiGHS
using Plots
using Printf

const NSD_COST = 100.0
const COLUMN_TOLERANCE = 1e-6

struct ScenarioColumn
    y::Matrix{Float64}
    x::Matrix{Float64}
    z::Vector{Float64}
    cost::Float64
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

function column_already_exists(columns, candidate; atol = 1e-6)
    for column in columns
        if isapprox(column.cost, candidate.cost; atol = atol) &&
           all(isapprox.(column.y, candidate.y; atol = atol))
            return true
        end
    end
    return false
end

function build_restricted_master(
    I,
    J,
    S,
    F,
    P,
    columns_by_scenario;
    relax_integrality = true,
)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    if relax_integrality
        @variable(model, 0 <= Y[1:I, 1:J] <= 1)
    else
        @variable(model, Y[1:I, 1:J], Bin)
    end

    lambda = [VariableRef[] for _ = 1:S]
    for s = 1:S
        for _ in eachindex(columns_by_scenario[s])
            variable =
                relax_integrality ? @variable(model, lower_bound = 0) :
                @variable(model, binary = true)
            push!(lambda[s], variable)
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

function solve_pricing_problem(s, I, J, A, B, C, P, convexity_dual, linking_dual)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, 0 <= X[i = 1:I, j = 1:J] <= min(A[i], B[s, j]))
    @variable(model, Y[i = 1:I, j = 1:J], Bin)
    @variable(model, Z[j = 1:J] >= 0)

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

function print_iteration(
    iteration,
    objective_value_master,
    best_reduced_cost,
    total_columns,
    new_columns,
)
    println(
        lpad(iteration, 9),
        " ",
        @sprintf("%12.4f", objective_value_master),
        " ",
        @sprintf("%12.4e", best_reduced_cost),
        " ",
        lpad(total_columns, 8),
        " ",
        lpad(new_columns, 11),
    )
    return
end

function column_generation!(I, J, S, A, B, C, F, P; max_iterations = 100)
    for s = 1:S
        if sum(B[s, :]) >= sum(A)
            @warn "The total demand for scenario $s exceeds the total supply. Non-supplied demand will absorb the deficit."
        end
    end

    columns_by_scenario = [[build_artificial_column(I, J, B[s, :])] for s = 1:S]
    iteration_history = NamedTuple[]

    println("Iteration  LP Objective  Best Red. Cost  Columns  New Columns")
    for iteration = 1:max_iterations
        master = build_restricted_master(
            I,
            J,
            S,
            F,
            P,
            columns_by_scenario;
            relax_integrality = true,
        )
        optimize!(master.model)
        @assert is_solved_and_feasible(master.model; dual = true)

        best_reduced_cost = Inf
        new_columns = 0
        objective_value_master = objective_value(master.model)

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
                objective = objective_value_master,
                best_reduced_cost = best_reduced_cost,
                columns = total_columns,
                new_columns = new_columns,
            ),
        )
        print_iteration(
            iteration,
            objective_value_master,
            best_reduced_cost,
            total_columns,
            new_columns,
        )

        if new_columns == 0
            @info "No improving columns found. Terminating column generation."
            break
        end
    end

    final_master = build_restricted_master(
        I,
        J,
        S,
        F,
        P,
        columns_by_scenario;
        relax_integrality = false,
    )
    optimize!(final_master.model)
    @assert is_solved_and_feasible(final_master.model)

    return iteration_history, columns_by_scenario, final_master
end

function describe_selected_columns(columns_by_scenario, master)
    for s = 1:length(columns_by_scenario)
        println("Scenario $s selected columns:")
        any_selected = false
        for k in eachindex(columns_by_scenario[s])
            lambda_value = value(master.lambda[s][k])
            if lambda_value > 1e-6
                any_selected = true
                println(
                    "  column ",
                    lpad(k, 3),
                    " with lambda = ",
                    @sprintf("%.2f", lambda_value),
                    ", operating cost = ",
                    @sprintf("%.2f", columns_by_scenario[s][k].cost),
                )
            end
        end
        if !any_selected
            println("  no selected columns")
        end
    end
    return
end

function plot_column_generation_history(iteration_history)
    iterations = [entry.iteration for entry in iteration_history]
    lp_objective = [entry.objective for entry in iteration_history]
    best_reduced_cost = [entry.best_reduced_cost for entry in iteration_history]
    total_columns = [entry.columns for entry in iteration_history]
    new_columns = [entry.new_columns for entry in iteration_history]

    objective_plot = plot(
        iterations,
        lp_objective;
        label = "LP objective",
        xlabel = "Iteration",
        ylabel = "Objective value",
        title = "Restricted Master Convergence",
        lw = 2,
        marker = :circle,
    )

    reduced_cost_plot = plot(
        iterations,
        best_reduced_cost;
        label = "Best reduced cost",
        xlabel = "Iteration",
        ylabel = "Reduced cost",
        title = "Pricing Progress",
        lw = 2,
        marker = :diamond,
        color = :firebrick,
    )
    hline!(
        reduced_cost_plot,
        [0.0];
        linestyle = :dash,
        color = :black,
        label = "Zero threshold",
    )

    columns_plot = plot(
        iterations,
        [total_columns new_columns];
        label = ["Total columns" "New columns"],
        xlabel = "Iteration",
        ylabel = "Columns",
        title = "Column Pool Growth",
        lw = 2,
        marker = [:square :utriangle],
    )

    return plot(
        objective_plot,
        reduced_cost_plot,
        columns_plot;
        layout = (3, 1),
        size = (900, 900),
    )
end

function plot_investment_heatmaps(Y_cg, Y_extensive)
    cg_plot = heatmap(
        1:size(Y_cg, 2),
        1:size(Y_cg, 1),
        Y_cg;
        xlabel = "Destination j",
        ylabel = "Origin i",
        title = "Column Generation Investment",
        color = cgrad(:blues),
        clims = (0, 1),
    )

    extensive_plot = heatmap(
        1:size(Y_extensive, 2),
        1:size(Y_extensive, 1),
        Y_extensive;
        xlabel = "Destination j",
        ylabel = "Origin i",
        title = "Extensive-Form Investment",
        color = cgrad(:greens),
        clims = (0, 1),
    )

    return plot(cg_plot, extensive_plot; layout = (1, 2), size = (900, 400))
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

## Solve the problem with column generation
iteration_history, columns_by_scenario, final_master =
    column_generation!(I, J, S, A, B, C, F, P)

cg_objective = objective_value(final_master.model)
Y_cg = value.(final_master.Y)
@info "Restricted master objective after column generation: $cg_objective"
@show Y_cg
describe_selected_columns(columns_by_scenario, final_master)

## Solve the extensive form for validation
extensive_form = create_extensive_form_model(I, J, S, A, B, C, F, P)
optimize!(extensive_form)
@assert is_solved_and_feasible(extensive_form)

extensive_objective = objective_value(extensive_form)
Y_extensive = value.(extensive_form[:Y])
@info "Extensive-form objective: $extensive_objective"
@info "Restricted-master gap against extensive form: $(cg_objective - extensive_objective)"
@show Y_extensive

if !isapprox(cg_objective, extensive_objective; atol = 1e-6)
    @warn "The restricted master does not match the extensive-form optimum. This is expected without branch-and-price if some needed columns were not generated."
end

history_plot = plot_column_generation_history(iteration_history)
display(history_plot)

investment_plot = plot_investment_heatmaps(Y_cg, Y_extensive)
display(investment_plot)
