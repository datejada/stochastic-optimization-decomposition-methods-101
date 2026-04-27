# Stochastic fixed-charge transportation problem (SFCTP)

## Define the packages
using JuMP    # used for mathematical programming
using HiGHS   # used as the solver

## Define model function
function create_and_optimize_sfctp_model(I, J, S, A, B, C, F, P)
    for s = 1:S
        if sum(B[s, :]) >= sum(A)
            @warn "The total demand for scenario $s exceeds the total supply. Consider increasing the supply or reducing the demand to avoid non-supplied demand."
        end
    end

    NSD_COST = 100 # non-supplied demand cost

    FCTP = Model(HiGHS.Optimizer)

    @variable(FCTP, 0 <= X[s = 1:S, i = 1:I, j = 1:J] <= min(A[i], B[s, j])) # arc flow
    @variable(FCTP, Y[i = 1:I, j = 1:J], Bin) # arc investment decision
    @variable(FCTP, 0 <= Z[s = 1:S, j = 1:J]) # slack variable for non-supplied demand

    @objective(
        FCTP,
        Min,
        sum(F[i, j] * Y[i, j] for i = 1:I, j = 1:J) +
        sum(C[i, j] * P[s] * X[s, i, j] for s = 1:S, i = 1:I, j = 1:J) +
        sum(NSD_COST * P[s] * Z[s, j] for s = 1:S, j = 1:J)
    )

    @constraint(FCTP, Offer[s = 1:S, i = 1:I], sum(X[s, i, j] for j = 1:J) ≤ A[i])
    @constraint(
        FCTP,
        Demand[s = 1:S, j = 1:J],
        sum(X[s, i, j] for i = 1:I) ≥ B[s, j] - Z[s, j]
    )
    @constraint(
        FCTP,
        FlowLimit[s = 1:S, i = 1:I, j = 1:J],
        X[s, i, j] <= min(A[i], B[s, j]) * Y[i, j]
    )

    optimize!(FCTP)
    obj_fun = objective_value(FCTP)
    @info "Status of the problem is $(termination_status(FCTP)) with o.f. $obj_fun"

    return FCTP
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

## Solve the stochastic model
stochastic = create_and_optimize_sfctp_model(I, J, S, A, B, C, F, P)
Y_stochastic = value.(stochastic[:Y])
@show Y_stochastic

## Solve deterministic scenarios
determ_sc = Dict()

### solve deterministic scenario 1
B_sc1 = reshape(B[1, :], (1, 3)) # reshape to make it a 1x3 matrix, as the function expects a 2D array for B (other option is to use: similar(B[1,:], Int64, 1, 3))
determ_sc[1] = create_and_optimize_sfctp_model(I, J, 1, A, B_sc1, C, F, [1])
Y_sc1 = value.(determ_sc[1][:Y])
@info "Optimal solution for scenario 1: $Y_sc1"

### solve deterministic scenario 2
B_sc2 = reshape(B[2, :], (1, 3))
determ_sc[2] = create_and_optimize_sfctp_model(I, J, 1, A, B_sc2, C, F, [1.0])
Y_sc2 = value.(determ_sc[2][:Y])
@info "Optimal solution for scenario 2: $Y_sc2"

### solve deterministic scenario 3
B_sc3 = reshape(B[3, :], (1, 3))
determ_sc[3] = create_and_optimize_sfctp_model(I, J, 1, A, B_sc3, C, F, [1.0])
Y_sc3 = value.(determ_sc[3][:Y])
@info "Optimal solution for scenario 3: $Y_sc3"

## Stochastic measures

### 1. save the solution from the stochastic model
recourse_problem_solution = objective_value(stochastic)
@info "Recourse problem solution: $recourse_problem_solution"

### 2. Determine the wait and see solution
wait_and_see_solution = sum(P[s] * objective_value(determ_sc[s]) for s = 1:S)
@info "Wait and see solution: $wait_and_see_solution"

### 3. Determine the expected value of perfect information (EVPI) or mean regret
expected_value_perfect_information = recourse_problem_solution - wait_and_see_solution
@info "Expected value of perfect information: $expected_value_perfect_information"

### 4. solve create the expected demand scenario
B_exp = sum(P .* B; dims = 1)
determ_expected_demand = create_and_optimize_sfctp_model(I, J, 1, A, B_exp, C, F, [1.0])

### 4.1 Save the solution of the expected value
investment_expected_scenario = value.(determ_expected_demand[:Y])
@info "Optimal solution for expected demand scenario: $investment_expected_scenario"

### 4.2 Fix the investment decision for each scenario 
for s = 1:S, i = 1:I, j = 1:J
    fix(determ_sc[s][:Y][i, j], investment_expected_scenario[i, j]; force = true)
end

### 4.3 solve each scenario again with the fixed investment decision
for s = 1:S
    optimize!(determ_sc[s])
end

### 4.4 Determine the expectation of EV Problem (EEV)
expectation_expected_value_solution = sum(P[s] * objective_value(determ_sc[s]) for s = 1:S)
@info "Expectation of expected value solution: $expectation_expected_value_solution"

### 5. Determine the Value of the stochastic solution (VSS) or Expected Value of Including Uncertainty (EVIU)
# Always positive for minimization problems.
value_stochastic_solution = expectation_expected_value_solution - recourse_problem_solution
@info "Value of the stochastic solution: $value_stochastic_solution"

## Simulate and out of sample evaluation
B_ofs = [20 20 70] # out-of-sample demand scenario

### 1. optimize the model with the out-of-sample demand
determ_ofs = create_and_optimize_sfctp_model(I, J, 1, A, B_ofs, C, F, [1.0])
investment_ofs = value.(determ_ofs[:Y])
@info "Optimal solution for out-of-sample demand scenario: $investment_ofs"

### 2. fix the investment decision from the stochastic solution for the out-of-sample demand scenario
for i = 1:I, j = 1:J
    fix(determ_ofs[:Y][i, j], Y_stochastic[i, j]; force = true)
end
optimize!(determ_ofs)
objective_value_ofs = objective_value(determ_ofs)
@info "Objective value for out-of-sample demand scenario with stochastic solution: $objective_value_ofs"
@info "Total non-supplied demand for out-of-sample scenario with stochastic solution: $(sum(value.(determ_ofs[:Z])))"

### 3. fix the investment decision from the scenario 1 solution for the out-of-sample demand scenario
for i = 1:I, j = 1:J
    fix(determ_ofs[:Y][i, j], Y_sc1[i, j]; force = true)
end
optimize!(determ_ofs)
objective_value_ofs = objective_value(determ_ofs)
@info "Objective value for out-of-sample demand scenario with scenario 1 solution: $objective_value_ofs"
@info "Total non-supplied demand for out-of-sample scenario with scenario 1 solution: $(sum(value.(determ_ofs[:Z])))"
