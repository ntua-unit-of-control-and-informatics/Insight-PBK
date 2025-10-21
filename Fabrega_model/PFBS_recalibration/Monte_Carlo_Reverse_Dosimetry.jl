using Distributions
using DataFrames
using CSV
using Random
using Statistics

# Set seed for reproducibility
Random.seed!(42)

# PFAS kinetic parameters with 10% CV uncertainty
PFAS_kinetic_params = Dict(
    "PFHxA" => (
        Vd_mean = 119, 
        Vd_cv = 0.1,
        k_el_mean = 0.478, 
        k_el_cv = 0.1
    ),
    "PFBS" => (
        Vd_mean = 137, 
        Vd_cv = 0.05,
        k_el_mean = 0.0137, 
        k_el_cv = 0.1
    )
)

# Steady-state concentration data with placeholder standard deviations
study_data = Dict(
    "Zhou_PFBS" => (css_mean = 19.3, css_sd = 19.7),
    "Zhou_PFHxA" => (css_mean = 0.41, css_sd = 0.25),
    "Olsen_PFBS" => (css_mean = 397, css_sd = 115.0*(6^(1/2))), # sd estimated from standard error and N=6
    "He_PFBS" => (css_mean = 64.7, css_sd = 98.75)
)

function create_parameter_distributions(pfas_type)
    params = PFAS_kinetic_params[pfas_type]
    
    # Convert to LogNormal parameters
    # For LogNormal: μ_log = log(μ) - 0.5*log(1 + (σ/μ)²)
    # σ_log = sqrt(log(1 + (σ/μ)²))
    
    # Vd distribution
    cv_vd = params.Vd_cv
    μ_log_vd = log(params.Vd_mean) - 0.5 * log(1 + cv_vd^2)
    σ_log_vd = sqrt(log(1 + cv_vd^2))
    Vd_dist = LogNormal(μ_log_vd, σ_log_vd)
    
    # k_el distribution
    cv_kel = params.k_el_cv
    μ_log_kel = log(params.k_el_mean) - 0.5 * log(1 + cv_kel^2)
    σ_log_kel = sqrt(log(1 + cv_kel^2))
    k_el_dist = LogNormal(μ_log_kel, σ_log_kel)
    
    return (Vd = Vd_dist, k_el = k_el_dist)
end

function monte_carlo_reverse_dosimetry(pfas_type, study_key, n_simulations = 10000)
    # Get parameter distributions
    param_dists = create_parameter_distributions(pfas_type)
    
    # Get steady-state concentration distribution  
    study = study_data[study_key]
    
    # Convert Css to LogNormal distribution
    cv_css = study.css_sd / study.css_mean
    μ_log_css = log(study.css_mean) - 0.5 * log(1 + cv_css^2)
    σ_log_css = sqrt(log(1 + cv_css^2))
    css_dist = LogNormal(μ_log_css, σ_log_css)
    
    # Monte Carlo sampling
    estimated_intakes = Float64[]
    
    for i in 1:n_simulations
        # Sample parameters from LogNormal distributions (naturally positive)
        Vd_sample = rand(param_dists.Vd)
        k_el_sample = rand(param_dists.k_el)
        css_sample = rand(css_dist)
        
        # Calculate estimated intake: intake = k_el * Css * Vd
        estimated_intake = k_el_sample * css_sample * Vd_sample
        push!(estimated_intakes, estimated_intake)
    end
    
    return estimated_intakes
end

function analyze_results(intakes, study_name)
    results = Dict(
        "study" => study_name,
        "mean" => mean(intakes),
        "median" => median(intakes),
        "std" => std(intakes),
        "q25" => quantile(intakes, 0.25),
        "q75" => quantile(intakes, 0.75),
        "q95" => quantile(intakes, 0.95),
        "min" => minimum(intakes),
        "max" => maximum(intakes)
    )
    return results
end

# Run Monte Carlo simulations for each study
println("Running Monte Carlo Reverse Dosimetry Simulations...")
println("=" ^ 50)

n_sims = 10000
all_results = []

# PFBS studies
for (study_key, study_name) in [
    ("Zhou_PFBS", "Zhou et al. 2014 - PFBS"),
    ("Olsen_PFBS", "Olsen et al. 2009 - PFBS"),
    ("He_PFBS", "He et al. 2023 - PFBS")
]
    intakes = monte_carlo_reverse_dosimetry("PFBS", study_key, n_sims)
    results = analyze_results(intakes, study_name)
    push!(all_results, results)
    
    println("$(study_name):")
    println("  Mean intake: $(round(results["mean"], digits=2)) ng/kg-bw/day")
    println("  Median intake: $(round(results["median"], digits=2)) ng/kg-bw/day")
    println("  95th percentile: $(round(results["q95"], digits=2)) ng/kg-bw/day")
    println()
end

# PFHxA study
intakes_pfhxa = monte_carlo_reverse_dosimetry("PFHxA", "Zhou_PFHxA", n_sims)
results_pfhxa = analyze_results(intakes_pfhxa, "Zhou et al. 2014 - PFHxA")
push!(all_results, results_pfhxa)

println("Zhou et al. 2014 - PFHxA:")
println("  Mean intake: $(round(results_pfhxa["mean"], digits=2)) ng/kg-bw/day")
println("  Median intake: $(round(results_pfhxa["median"], digits=2)) ng/kg-bw/day")
println("  95th percentile: $(round(results_pfhxa["q95"], digits=2)) ng/kg-bw/day")

# Save results to DataFrame and CSV
df_results = DataFrame(
    study = [r["study"] for r in all_results],
    mean = [r["mean"] for r in all_results],
    median = [r["median"] for r in all_results],
    std = [r["std"] for r in all_results],
    q25 = [r["q25"] for r in all_results],
    q75 = [r["q75"] for r in all_results],
    q95 = [r["q95"] for r in all_results],
    min = [r["min"] for r in all_results],
    max = [r["max"] for r in all_results]
)

CSV.write("PFBS_recalibration/MC_reverse_dosimetry_results.csv", df_results)
println("\nResults saved to PFBS_recalibration/MC_reverse_dosimetry_results.csv")