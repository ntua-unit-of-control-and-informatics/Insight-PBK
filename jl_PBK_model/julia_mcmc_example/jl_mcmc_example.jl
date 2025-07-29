using Turing
using DifferentialEquations
# Load StatsPlots for visualizations and diagnostics.
using StatsPlots
using LinearAlgebra
using Distributions
# Set a seed for reproducibility.
using Random
Random.seed!(14);

# Define Lotka–Volterra model.
function lotka_volterra(du, u, p, t)
    # Model parameters.
    α, β, γ, δ = p
    # Current state.
    x, y = u

    # Evaluate differential equations.
    du[1] = (α - β * y) * x # prey
    du[2] = (δ * x - γ) * y # predator

    return nothing
end

# Define initial-value problem.
u0 = [1.0, 1.0]
p = [1.5, 1.0, 3.0, 1.0]
tspan = (0.0, 10.0)
prob = ODEProblem(lotka_volterra, u0, tspan, p)

# Plot simulation.
plot(solve(prob, Tsit5()))

sol = solve(prob, Tsit5(); saveat=0.1)
q = 1.7
odedata = rand.(Poisson.(q * Array(sol)))

# Plot simulation and noisy observations.
plot(sol, label=["Prey" "Predator"])
scatter!(sol.t, odedata'; color=[1 2], label="")

@model function fitlv(data, prob)
    # Prior distributions.
    α ~ truncated(Normal(1.5, 0.2); lower=0.5, upper=2.5)
    β ~ truncated(Normal(1.1, 0.2); lower=0, upper=2)
    γ ~ truncated(Normal(3.0, 0.2); lower=1, upper=4)
    δ ~ truncated(Normal(1.0, 0.2); lower=0, upper=2)
    q ~ truncated(Normal(1.7, 0.2); lower=0, upper=3)

    # Simulate Lotka–Volterra model. 
    p = [α, β, γ, δ]
    predicted = solve(prob, Tsit5(); p=p, saveat=0.1, abstol=1e-6, reltol=1e-6)
    ϵ = 1e-5
    
    # Observations.
    for i in eachindex(predicted)
        data[:, i] ~ arraydist(Poisson.(q .* predicted[i] .+ ϵ))
    end

    return nothing
end

model = fitlv(odedata, prob)

# Sample 3 independent chains with forward-mode automatic differentiation (the default).
chain = sample(model, NUTS(), MCMCSerial(), 1000, 3; progress=false)
# Print summary statistics of the chain
println("Summary statistics of the MCMC chains:")
describe(chain)
# Print quantiles of the posterior distributions

# Create the plot
plot(chain)
# Save the plot to a file
savefig(p, "jl_mcmc_example/mcmc_chains_plot.png")
println("Plot saved as mcmc_chains_plot.png")

