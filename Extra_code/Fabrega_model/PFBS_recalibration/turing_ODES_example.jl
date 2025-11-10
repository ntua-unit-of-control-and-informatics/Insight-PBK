using Turing
using DifferentialEquations
using Sundials
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

# Define prior parameters - CHANGE THESE VALUES TO UPDATE ALL PRIORS
const PRIORS = Dict(
    :α => (mean=1.5, std=0.2, lower=0.5, upper=2.5),
    :β => (mean=1.1, std=0.2, lower=0, upper=2),
    :γ => (mean=3.0, std=0.2, lower=1, upper=4),
    :δ => (mean=1.0, std=0.2, lower=0, upper=2),
    :q => (mean=1.7, std=0.2, lower=0, upper=3)
)

# Define initial-value problem.
u0 = [1.0, 1.0]
p = [1.5, 1.0, 3.0, 1.0]
tspan = (0.0, 10.0)
prob = ODEProblem(lotka_volterra, u0, tspan, p)

sol = solve(prob, Tsit5(); saveat=0.1)
q = 1.7
odedata = rand.(Poisson.(q * Array(sol)))

@model function fitlv(data, prob)
    # Prior distributions - automatically generated from PRIORS dict
    α ~ truncated(Normal($(PRIORS[:α].mean), $(PRIORS[:α].std)); lower=$(PRIORS[:α].lower), upper=$(PRIORS[:α].upper))
    β ~ truncated(Normal($(PRIORS[:β].mean), $(PRIORS[:β].std)); lower=$(PRIORS[:β].lower), upper=$(PRIORS[:β].upper))
    γ ~ truncated(Normal($(PRIORS[:γ].mean), $(PRIORS[:γ].std)); lower=$(PRIORS[:γ].lower), upper=$(PRIORS[:γ].upper))
    δ ~ truncated(Normal($(PRIORS[:δ].mean), $(PRIORS[:δ].std)); lower=$(PRIORS[:δ].lower), upper=$(PRIORS[:δ].upper))
    q ~ truncated(Normal($(PRIORS[:q].mean), $(PRIORS[:q].std)); lower=$(PRIORS[:q].lower), upper=$(PRIORS[:q].upper))

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
describe(chain)

# Print prior and posterior distributions
using Statistics
using MCMCChains

println("\n" * "="^60)
println("PRIOR vs POSTERIOR COMPARISON")
println("="^60)

# Generate prior_info automatically from PRIORS dictionary
prior_info = Dict(
    param => (
        distribution = "truncated(Normal($(params.mean), $(params.std)); lower=$(params.lower), upper=$(params.upper))",
        mean = params.mean,
        std = params.std
    ) for (param, params) in PRIORS
)

# Extract posterior statistics
param_names = names(chain)
for param in param_names
    if param in keys(prior_info)
        # Posterior statistics
        chain_values = vec(Array(chain[param]))
        post_mean = mean(chain_values)
        post_std = std(chain_values)
        post_q025 = quantile(chain_values, 0.025)
        post_q975 = quantile(chain_values, 0.975)
        
        # Prior information
        prior = prior_info[param]
        
        println("\n$param:")
        println("  Prior: $(prior.distribution)")
        println("  Prior mean ± SD: $(prior.mean) ± $(prior.std)")
        println("  Posterior mean ± SD: $(round(post_mean, digits=4)) ± $(round(post_std, digits=4))")
        println("  Posterior 95% CI: [$(round(post_q025, digits=4)), $(round(post_q975, digits=4))]")
    end
end

println("\n" * "="^60)
