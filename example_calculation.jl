# Demonstration of RQA calculation using sampling approach
# N. Marwan: Energy-efficient recurrence quantification analysis, 
#    European Physical Journal – Special Topics (in press).
#    DOI:10.1140/epjs/s11734-025-02121-w

# Import required packages
using DifferentialEquations
using DelayEmbeddings
using Statistics
include("RPLineLengths.jl")

# Create exemplary data (Roessler)
function roessler!(du, u, p, t)
    du[1] = -u[2] - u[3]
    du[2] = u[1] + a * u[2]
    du[3] = 0.2 + u[3] * (u[1] - 5.7)
end

N = 10000                 # length of final time series
N_trans = 100             # transients to be removed from time series
a = 0.2                   # control parameter in the Roessler system
deltaT = 0.2              # sampling time
u0 = [-6.2668, -1.3413, 0.0166]   # initial conditions

prob = ODEProblem(roessler!, u0, (0.0, (N+N_trans-1)*deltaT)) # problem statement for the solver
sol = solve(prob, Tsit5(), saveat=deltaT) # solving ODEs
x = transpose(hcat(sol.u...))     # get solution as a matrix
x = x[N_trans+1:end, :]   # remove transients

# RQA estimation
ε = 0.1 * (maximum(x) - minimum(x))    # recurrence threshold value
M = Int(round(0.2 * N))   # number of random subsamples

# Get histograms of line lengths
@time histL1 = get_hist_diagonal_woRP(x, ε); # standard calculation (using efficient woRP approach)
@time histL2, N_rp = get_hist_diagonal_sampled(x, ε, M); # calculation using sampling schema

# Get RQA measures
r1 = rqa(histL1, N*(N-1)/2) # RQA using standard approach
r2 = rqa(histL2, N_rp)       # RQA using sampling approach
