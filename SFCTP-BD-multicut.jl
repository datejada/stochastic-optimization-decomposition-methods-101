# Stochastic fixed-charge transportation problem (SFCTP) with Bender's decomposition and Multi cuts

## Define the packages
using JuMP
using HiGHS
using DataFrames
using Printf
using Plots

## Define functions
function create_first_stage_model(I, J, S, F, P)
    # Scalar values
    M = -10000 # to avoid the problem is unbounded 

    first_stage = Model(HiGHS.Optimizer)
    set_silent(first_stage)

    @variable(first_stage, Y[i = 1:I, j = 1:J], Bin) # arc investment decision
    @variable(first_stage, theta[s = 1:S] ≥ M)       # Benders' cut

    @objective(
        first_stage,
        Min,
        sum(F[i, j] * Y[i, j] for i = 1:I, j = 1:J) + sum(P[s] * theta[s] for s = 1:S)
    )

    return first_stage
end

function create_and_solve_subproblem(I, J, A, B, C, Y_proposal)

    NSD_COST = 100 # non-supplied demand cost

    subproblem = Model(HiGHS.Optimizer)
    set_silent(subproblem)

    @variable(subproblem, 0 <= X[i = 1:I, j = 1:J] <= min(A[i], B[j])) # arc flow
    @variable(subproblem, 0 <= Z[j = 1:J]) # slack variable for non-supplied demand

    @objective(
        subproblem,
        Min,
        sum(C[i, j] * X[i, j] for i = 1:I, j = 1:J) + sum(NSD_COST * Z[j] for j = 1:J)
    )

    @constraint(subproblem, Offer[i = 1:I], sum(X[i, j] for j = 1:J) ≤ A[i])
    @constraint(subproblem, Demand[j = 1:J], sum(X[i, j] for i = 1:I) ≥ B[j] - Z[j])
    @constraint(
        subproblem,
        FlowLimit[i = 1:I, j = 1:J],
        X[i, j] <= min(A[i], B[j]) * Y_proposal[i, j]
    )

    optimize!(subproblem)

    # Check if the model is optimal
    @assert is_solved_and_feasible(subproblem; dual = true) # Check if the model is optimal and there is a dual solution

    return subproblem
end

function add_cut!(first_stage_model, subproblem, I, J, S, A, B, Y_proposal, iteration)
    # Add Benders' cut
    for s = 1:S
        p_dual = dual.(subproblem[s][:FlowLimit])
        @constraint(
            first_stage_model,
            base_name = "cut_iter_$(iteration)_scenario_$(s)",
            first_stage_model[:theta][s] >=
            objective_value(subproblem[s]) + sum(
                -p_dual[i, j] *
                min(A[i], B[s, j]) *
                (Y_proposal[i, j] - first_stage_model[:Y][i, j]) for i = 1:I, j = 1:J
            )
        )
    end

    return nothing
end

"""
    print_iteration(iteration)

Prints the current iteration number.

# Arguments
- `iteration`: An integer representing the current iteration number.
"""
function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end

function benders_decomposition!(I, J, S, A, B, C, F, P)
    for s = 1:S
        if sum(B[s, :]) >= sum(A)
            @warn "The total demand for scenario $s exceeds the total supply. Consider increasing the supply or reducing the demand to avoid non-supplied demand."
        end
    end

    # Define constants
    MAXIMUM_ITERATIONS = 100       # maximum number of iterations
    ABSOLUTE_OPTIMALITY_GAP = 1e-6 # gap between lower and upper bounds to stop the algorithm

    # Define outputs
    df_Benders_interations = DataFrame(;
        iteration = Int[],
        lower_bound = Float64[],
        upper_bound = Float64[],
        gap = Float64[],
    )

    df_investment_per_iteration =
        DataFrame(; iteration = Int[], connection = String[], investment = Float64[])

    # Create first stage model
    first_stage_model = create_first_stage_model(I, J, S, F, P)

    # Start Benders' Decomposition
    upper_bound = Inf
    println("Iteration  Lower Bound  Upper Bound          Gap")
    for iteration = 1:MAXIMUM_ITERATIONS
        optimize!(first_stage_model)
        @assert is_solved_and_feasible(first_stage_model)
        lower_bound = objective_value(first_stage_model)
        Y_proposal = value.(first_stage_model[:Y])
        subproblem = Dict{Int,JuMP.Model}()
        for s = 1:S
            B_sc = reshape(B[s, :], (1, 3))
            subproblem[s] = create_and_solve_subproblem(I, J, A, B_sc, C, Y_proposal)
        end
        upper_bound = minimum([
            upper_bound,
            objective_value(first_stage_model) -
            sum(P[s] * value.(first_stage_model[:theta][s]) for s = 1:S) +
            sum(P[s] * objective_value(subproblem[s]) for s = 1:S),
        ])
        gap = abs(upper_bound - lower_bound) / abs(upper_bound)
        push!(df_Benders_interations, (iteration, lower_bound, upper_bound, gap))
        print_iteration(iteration, lower_bound, upper_bound, gap)

        if gap < ABSOLUTE_OPTIMALITY_GAP
            println("Terminating with the optimal solution")
            break
        end
        add_cut!(first_stage_model, subproblem, I, J, S, A, B, Y_proposal, iteration)

        for i = 1:I, j = 1:J
            ij = "($i, $j)"
            push!(df_investment_per_iteration, (iteration, ij, Y_proposal[i, j]))
        end
    end

    return df_Benders_interations, df_investment_per_iteration, first_stage_model
end

## Define data

### system dimensions
I = 4  # origins
J = 3  # destinations
S = 3  # scenarios

###  parameters
A = [20, 30, 40, 20] # A(i) product offer 
P = [0.5, 0.3, 0.2]  # P(s) probability of each scenario 

### - B(s,j) product demand per scenario 
B = [
    21 51 31
    32 22 52
    53 33 23
]

### - C(i,j) per unit variable transportation cost
C = [
    1 2 3
    3 2 1
    2 3 4
    4 3 2
]

### - F(i,j) fixed transportation cost
F = [
    10 20 30
    20 30 40
    30 40 50
    40 50 60
]

## Solve the stochastic model using Benders' decomposition
df_Benders_interations, df_investment_per_iteration, first_stage_model =
    benders_decomposition!(I, J, S, A, B, C, F, P)
@info "Objective value at the final iteration: $(objective_value(first_stage_model))"
stochastic_solution = value.(first_stage_model[:Y])
## Plot upper and lower bounds per iteration
plot(
    df_Benders_interations.iteration,
    [df_Benders_interations.lower_bound, df_Benders_interations.upper_bound],
    label = ["Lower Bound" "Upper Bound"],
    xlabel = "Iteration",
    ylabel = "Bound",
    title = "Benders' Decomposition Convergence",
)
## Plot gap per iteration
plot(
    df_Benders_interations.iteration,
    df_Benders_interations.gap,
    label = "Gap",
    xlabel = "Iteration",
    ylabel = "Gap",
    title = "Benders' Decomposition Gap",
)

## Plot investment decisions per iteration for a given connection (i,j)
connection_to_plot = "(2, 3)" # specify the connection to plot
df_filtered =
    filter(row -> row.connection == connection_to_plot, df_investment_per_iteration)
plot(
    df_filtered.iteration,
    df_filtered.investment,
    label = connection_to_plot,
    xlabel = "Iteration",
    ylabel = "Investment Decision",
    title = "Investment Decisions per Iteration",
)
## Print the first stage problem .lp
write_to_file(first_stage_model, "first_stage_model.lp")
print(first_stage_model)
