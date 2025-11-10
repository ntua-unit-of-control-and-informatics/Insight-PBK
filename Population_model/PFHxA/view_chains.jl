using Serialization
using Turing
using MCMCChains
using StatsPlots
using Plots
using Statistics

println("="^80)
println("PFHxA SINGLE-DOSE MODEL - CHAIN DIAGNOSTICS")
println("="^80)

# Load MCMC chains (choose which file to load)
println("\nLoading MCMC chains...")
chains = deserialize("Population_model/PFHxA/PFHxA_with_urine_chains.jls")  # or PFHxA_chains.jls
println("✓ Chains loaded successfully")

# Display summary statistics
println("\n" * "="^80)
println("CHAIN SUMMARY STATISTICS")
println("="^80)
describe(chains)

# Calculate parameter summaries
println("\n" * "="^80)
println("PARAMETER ESTIMATES (POSTERIOR SUMMARIES)")
println("="^80)

# RAFapi
RAFapi_samples = Array(chains[:RAFapi])
RAFapi_mean = mean(RAFapi_samples)
RAFapi_sd = std(RAFapi_samples)
RAFapi_ci = quantile(RAFapi_samples, [0.025, 0.975])

println("\nRAFapi:")
println("  Mean ± SD: $(round(RAFapi_mean, sigdigits=4)) ± $(round(RAFapi_sd, sigdigits=4))")
println("  95% CI: [$(round(RAFapi_ci[1], sigdigits=4)), $(round(RAFapi_ci[2], sigdigits=4))]")

# RAFbaso
RAFbaso_samples = Array(chains[:RAFbaso])
RAFbaso_mean = mean(RAFbaso_samples)
RAFbaso_sd = std(RAFbaso_samples)
RAFbaso_ci = quantile(RAFbaso_samples, [0.025, 0.975])

println("\nRAFbaso:")
println("  Mean ± SD: $(round(RAFbaso_mean, sigdigits=4)) ± $(round(RAFbaso_sd, sigdigits=4))")
println("  95% CI: [$(round(RAFbaso_ci[1], sigdigits=4)), $(round(RAFbaso_ci[2], sigdigits=4))]")

# Km_apical
Km_apical_samples = Array(chains[:Km_apical])
Km_apical_mean = mean(Km_apical_samples)
Km_apical_sd = std(Km_apical_samples)
Km_apical_ci = quantile(Km_apical_samples, [0.025, 0.975])

println("\nKm_apical:")
println("  Mean ± SD: $(round(Km_apical_mean, sigdigits=4)) ± $(round(Km_apical_sd, sigdigits=4)) μg/L")
println("  95% CI: [$(round(Km_apical_ci[1], sigdigits=4)), $(round(Km_apical_ci[2], sigdigits=4))] μg/L")

# PC_scale
PC_scale_samples = Array(chains[:PC_scale])
PC_scale_mean = mean(PC_scale_samples)
PC_scale_sd = std(PC_scale_samples)
PC_scale_ci = quantile(PC_scale_samples, [0.025, 0.975])

println("\nPC_scale:")
println("  Mean ± SD: $(round(PC_scale_mean, sigdigits=4)) ± $(round(PC_scale_sd, sigdigits=4))")
println("  95% CI: [$(round(PC_scale_ci[1], sigdigits=4)), $(round(PC_scale_ci[2], sigdigits=4))]")

# Observation errors (check if sigma_plasma and sigma_urine exist, otherwise just sigma)
if :sigma_plasma in keys(chains)
    sigma_plasma_samples = Array(chains[:sigma_plasma])
    sigma_plasma_mean = mean(sigma_plasma_samples)
    sigma_plasma_sd = std(sigma_plasma_samples)
    sigma_plasma_ci = quantile(sigma_plasma_samples, [0.025, 0.975])

    println("\nObservation error - Plasma (sigma_plasma):")
    println("  Mean ± SD: $(round(sigma_plasma_mean, digits=4)) ± $(round(sigma_plasma_sd, digits=4))")
    println("  95% CI: [$(round(sigma_plasma_ci[1], digits=4)), $(round(sigma_plasma_ci[2], digits=4))]")

    sigma_urine_samples = Array(chains[:sigma_urine])
    sigma_urine_mean = mean(sigma_urine_samples)
    sigma_urine_sd = std(sigma_urine_samples)
    sigma_urine_ci = quantile(sigma_urine_samples, [0.025, 0.975])

    println("\nObservation error - Urine (sigma_urine):")
    println("  Mean ± SD: $(round(sigma_urine_mean, digits=4)) ± $(round(sigma_urine_sd, digits=4))")
    println("  95% CI: [$(round(sigma_urine_ci[1], digits=4)), $(round(sigma_urine_ci[2], digits=4))]")
else
    sigma_samples = Array(chains[:sigma])
    sigma_mean = mean(sigma_samples)
    sigma_sd = std(sigma_samples)
    sigma_ci = quantile(sigma_samples, [0.025, 0.975])

    println("\nObservation error (sigma):")
    println("  Mean ± SD: $(round(sigma_mean, digits=4)) ± $(round(sigma_sd, digits=4))")
    println("  95% CI: [$(round(sigma_ci[1], digits=4)), $(round(sigma_ci[2], digits=4))]")
end

# Generate diagnostic plots
println("\n" * "="^80)
println("GENERATING DIAGNOSTIC PLOTS")
println("="^80)

chains_plot = plot(chains, size=(1800, 2800), dpi=300,
                    left_margin=25Plots.mm, bottom_margin=8Plots.mm,
                    top_margin=6Plots.mm, right_margin=6Plots.mm)
savefig(chains_plot, "Population_model/PFHxA/PFHxA_population_diagnostic_plots.png")
println("✓ Saved to: Population_model/PFHxA/PFHxA_population_diagnostic_plots.png")

println("\n" * "="^80)
println("DIAGNOSTICS COMPLETE")
println("="^80)
