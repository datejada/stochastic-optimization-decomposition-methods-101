# Fixed-charge transportation problem (FCTP)

# Define the packages
using JuMP    # used for mathematical programming
using HiGHS   # used as the solver

# Define model function
function create_and_optimize_fctp_model(I, J, A, B, C, F)

    @assert sum(B) < sum(A) "Infeasible problem"

    FCTP = Model(HiGHS.Optimizer)

    @variable(FCTP, 0 <= X[i = 1:I, j = 1:J] <= min(A[i], B[j])) # arc flow
    @variable(FCTP, Y[i = 1:I, j = 1:J], Bin) # arc investment decision

    @objective(
        FCTP,
        Min,
        sum(F[i, j] * Y[i, j] for i = 1:I, j = 1:J) +
        sum(C[i, j] * X[i, j] for i = 1:I, j = 1:J)
    )

    @constraint(FCTP, Offer[i = 1:I], sum(X[i, j] for j = 1:J) ≤ A[i])
    @constraint(FCTP, Demand[j = 1:J], sum(X[i, j] for i = 1:I) ≥ B[j])
    @constraint(FCTP, FlowLimit[i = 1:I, j = 1:J], X[i, j] <= min(A[i], B[j]) * Y[i, j])

    optimize!(FCTP)
    Z1 = objective_value(FCTP)
    println("Status of the problem is ", termination_status(FCTP), " with o.f. ", Z1)

    return FCTP
end

# Define data
# system dimensions
I = 4  # origins
J = 3  # destinations

#  parameters
A = [20, 30, 40, 20] # A(i) product offer 
B = [21, 51, 31]     # B(j) product demand   

# - C(i,j) per unit variable transportation cost
C = [
    1 2 3
    3 2 1
    2 3 4
    4 3 2
]

# - F(i,j) fixed transportation cost
F = [
    10 20 30
    20 30 40
    30 40 50
    40 50 60
]

# deterministic scenario 1
determ_sc1 = create_and_optimize_fctp_model(I, J, A, B, C, F)

# demand scenario 2
B = [32, 22, 52]
determ_sc2 = create_and_optimize_fctp_model(I, J, A, B, C, F)

# demand scenario 3
B = [53, 33, 23]
determ_sc3 = create_and_optimize_fctp_model(I, J, A, B, C, F)

# mean demand scenario
B = [sum([21, 32, 53]) / 3, sum([51, 22, 33]) / 3, sum([31, 52, 23]) / 3]
determ_mean = create_and_optimize_fctp_model(I, J, A, B, C, F)
