using Serialization
using Turing
using MCMCChains
using StatsPlots
using Plots
using Statistics

println("="^80)
println("PFBS POPULATION MODEL - CHAIN DIAGNOSTICS")
println("="^80)

# Load MCMC chains
println("\nLoading MCMC chains...")
chains = deserialize("PFAS_models/PFBS/PFBS_population_chains.jls")
println("✓ Chains loaded successfully")

# Display summary statistics
println("\n" * "="^80)
println("CHAIN SUMMARY STATISTICS")
println("="^80)
describe(chains)

# Calculate derived quantities
println("\n" * "="^80)
println("POPULATION PARAMETERS (ORIGINAL SCALE)")
println("="^80)

# Extract posterior samples
mu_pop_samples = Array(chains[:mu_pop])
sigma_pop_samples = Array(chains[:sigma_pop])

# Transform to original scale
mean_RAFapi_samples = exp.(mu_pop_samples .+ sigma_pop_samples.^2 ./ 2)
var_RAFapi_samples = (exp.(sigma_pop_samples.^2) .- 1) .* exp.(2 .* mu_pop_samples .+ sigma_pop_samples.^2)
sd_RAFapi_samples = sqrt.(var_RAFapi_samples)
cv_RAFapi_samples = sqrt.(exp.(sigma_pop_samples.^2) .- 1)

# Calculate posterior summaries
mean_RAFapi_mean = mean(mean_RAFapi_samples)
mean_RAFapi_sd = std(mean_RAFapi_samples)
mean_RAFapi_ci = quantile(mean_RAFapi_samples, [0.025, 0.975])

sd_RAFapi_mean = mean(sd_RAFapi_samples)
sd_RAFapi_sd = std(sd_RAFapi_samples)
sd_RAFapi_ci = quantile(sd_RAFapi_samples, [0.025, 0.975])

cv_RAFapi_mean = mean(cv_RAFapi_samples)
cv_RAFapi_sd = std(cv_RAFapi_samples)
cv_RAFapi_ci = quantile(cv_RAFapi_samples, [0.025, 0.975])

println("\nPopulation Mean RAFapi:")
println("  Mean ± SD: $(round(mean_RAFapi_mean, sigdigits=4)) ± $(round(mean_RAFapi_sd, sigdigits=4))")
println("  95% CI: [$(round(mean_RAFapi_ci[1], sigdigits=4)), $(round(mean_RAFapi_ci[2], sigdigits=4))]")

println("\nPopulation SD RAFapi:")
println("  Mean ± SD: $(round(sd_RAFapi_mean, sigdigits=4)) ± $(round(sd_RAFapi_sd, sigdigits=4))")
println("  95% CI: [$(round(sd_RAFapi_ci[1], sigdigits=4)), $(round(sd_RAFapi_ci[2], sigdigits=4))]")

println("\nPopulation CV:")
println("  Mean ± SD: $(round(cv_RAFapi_mean, digits=3)) ± $(round(cv_RAFapi_sd, digits=3))")
println("  95% CI: [$(round(cv_RAFapi_ci[1], digits=3)), $(round(cv_RAFapi_ci[2], digits=3))]")

Km_apical_samples = Array(chains[:Km_apical])
Km_apical_mean = mean(Km_apical_samples)
Km_apical_sd = std(Km_apical_samples)
Km_apical_ci = quantile(Km_apical_samples, [0.025, 0.975])

println("\nKm_apical (original scale):")
println("  Mean ± SD: $(round(Km_apical_mean, sigdigits=4)) ± $(round(Km_apical_sd, sigdigits=4)) μg/L")
println("  95% CI: [$(round(Km_apical_ci[1], sigdigits=4)), $(round(Km_apical_ci[2], sigdigits=4))] μg/L")

println("\nLog-scale parameters:")
println("  mu_pop (mean of log(RAFapi)): $(round(mean(chains[:mu_pop]), digits=3))")
println("  sigma_pop (SD of log(RAFapi)): $(round(mean(chains[:sigma_pop]), digits=3))")

# Generate diagnostic plots
println("\n" * "="^80)
println("GENERATING DIAGNOSTIC PLOTS")
println("="^80)

chains_plot = plot(chains, size=(1800, 2800), dpi=300,
                    left_margin=25Plots.mm, bottom_margin=8Plots.mm,
                    top_margin=6Plots.mm, right_margin=6Plots.mm)
savefig(chains_plot, "PFAS_models/PFBS/PPC_hierarchical/PFBS_population_diagnostic_plots.png")
println("✓ Saved to: PFAS_models/PFBS/PPC_hierarchical/PFBS_population_diagnostic_plots.png")

println("\n" * "="^80)
println("DIAGNOSTICS COMPLETE")
println("="^80)
