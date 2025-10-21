using Turing
using DifferentialEquations
using DataFrames
using CSV
using MCMCChains
using StatsPlots
using Plots
using Distributions
using Distributions
using LinearAlgebra
using Random
using ReverseDiff
using Sundials

# Include original PBK functions
include("../jl_PBK_model/julia_model/pbk_functions.jl")

"""
PFBS Parameter Calibration using Hierarchical MCMC

This script calibrates Tm (maximum reabsorption rate) and Kt (Michaelis constant) 
parameters of the PBK model using serum and urine data from three studies:
- Olsen et al. (2009): Pharmacokinetic elimination phase
- He et al. (2023): Occupational exposure (20 years)
- Zhou et al. (2014): General population (10 years)

Approach:
- Population-level LogNormal priors for Tm, Kt (CV=200%)
- Proportional error model for serum and urine predictions
- 4 parallel chains, 10,000 iterations each
- New partition coefficients from EDM model
"""

# Set random seed for reproducibility
Random.seed!(42)

println("="^80)
println("PFBS PBK MODEL PARAMETER CALIBRATION")
println("="^80)

# Load all required data
println("Loading data sources...")

# 1. Load preprocessed study data
println("  Loading preprocessed study data...")
olsen_serum = CSV.read("data_preprocessing/olsen_2009_serum_summary.csv", DataFrame)
olsen_urine = CSV.read("data_preprocessing/olsen_2009_urine_summary.csv", DataFrame)
he_serum = CSV.read("data_preprocessing/he_2023_serum_summary.csv", DataFrame)
he_urine = CSV.read("data_preprocessing/he_2023_urine_summary.csv", DataFrame)
zhou_serum = CSV.read("data_preprocessing/zhou_2014_serum_summary.csv", DataFrame)
zhou_urine = CSV.read("data_preprocessing/zhou_2014_urine_summary.csv", DataFrame)

# 2. Load exposure estimates from Abraham model
println("  Loading exposure estimates...")
abraham_results = CSV.read("../Abraham_model/Abraham_model_results.csv", DataFrame)

# Extract PFBS exposures for each study (ng/kg-bw/day)
exposure_zhou = abraham_results[abraham_results.Study .== "Zhou et al., 2014" .&& abraham_results.PFAS .== "PFBS", :Estimated_Intake_ng_kg_bw_day][1]
exposure_olsen = abraham_results[abraham_results.Study .== "Olsen et al., 2009" .&& abraham_results.PFAS .== "PFBS", :Estimated_Intake_ng_kg_bw_day][1]
exposure_he = abraham_results[abraham_results.Study .== "He et al., 2023" .&& abraham_results.PFAS .== "PFBS", :Estimated_Intake_ng_kg_bw_day][1]

println("  Exposure estimates:")
println("    Zhou et al.: $(round(exposure_zhou, digits=1)) ng/kg-bw/day")
println("    Olsen et al.: $(round(exposure_olsen, digits=1)) ng/kg-bw/day") 
println("    He et al.: $(round(exposure_he, digits=1)) ng/kg-bw/day")

# 3. Load new partition coefficients
println("  Loading new partition coefficients...")
partition_coeffs = CSV.read("../EDM_Partition_Coefficients/partition_coefficients_results.csv", DataFrame)
pfbs_coeffs = filter(row -> row.Compound == "PFBS", partition_coeffs)

# Create partition coefficient dictionary
new_partition_coeffs = Dict(
    "liver" => pfbs_coeffs[pfbs_coeffs.Tissue .== "liver", :Partition_Coefficient][1],
    "adipose" => pfbs_coeffs[pfbs_coeffs.Tissue .== "adipose", :Partition_Coefficient][1],
    "kidney" => pfbs_coeffs[pfbs_coeffs.Tissue .== "kidney", :Partition_Coefficient][1],
    "gut" => pfbs_coeffs[pfbs_coeffs.Tissue .== "gut", :Partition_Coefficient][1],
    "lung" => pfbs_coeffs[pfbs_coeffs.Tissue .== "lung", :Partition_Coefficient][1],
    "brain" => pfbs_coeffs[pfbs_coeffs.Tissue .== "brain", :Partition_Coefficient][1],
    "rest_of_body" => pfbs_coeffs[pfbs_coeffs.Tissue .== "rest of body", :Partition_Coefficient][1]
)

println("  New partition coefficients:")
for (tissue, coeff) in new_partition_coeffs
    println("    $tissue: $(round(coeff, digits=3))")
end

# 4. Define exposure times (days) - TO BE FILLED BY USER
exposure_time_zhou = 10.0*365  # User to fill
exposure_time_olsen = 10.0*365  # User to fill  
exposure_time_he = 20*365  # User to fill

println("\n  Exposure times:")
println("    Zhou et al.: $exposure_time_zhou days")
println("    Olsen et al.: $exposure_time_olsen days")
println("    He et al.: $exposure_time_he days")

# Modified PBK parameter creation function with new partition coefficients
function create_params_updated(user_input, new_PCs)
    params = create_params(user_input)
    
    # Update partition coefficients with new values using splatting
    # params = (params..., 
    #           PL = new_PCs["liver"],
    #           PF = new_PCs["adipose"], 
    #           PK = new_PCs["kidney"],
    #           PG = new_PCs["gut"],
    #           PLu = new_PCs["lung"],
    #           PB = new_PCs["brain"],
    #           PR = new_PCs["rest_of_body"])
    
    return params
end

# Modified extract_concentrations function that includes urine concentration
function extract_concentrations_with_urine(sol, parameters)
    @unpack VPlas, VG, VL, VF, VLu, VB, VK, VFil, VR, Free = parameters
    
    # Standard tissue concentrations
    concentrations = Dict(
        "CPlas" => [u[1] for u in sol.u] ./ VPlas ./ Free,
        "CG" => [u[2] for u in sol.u] ./ VG,
        "CL" => [u[3] for u in sol.u] ./ VL,
        "CF" => [u[4] for u in sol.u] ./ VF,
        "CLu" => [u[5] for u in sol.u] ./ VLu,
        "CB" => [u[6] for u in sol.u] ./ VB,
        "CK" => [u[7] for u in sol.u] ./ VK,
        "CFil" => [u[8] for u in sol.u] ./ VFil,
        "CR" => [u[11] for u in sol.u] ./ VR
    )
    
    # Urine concentrations - use filtrate compartment concentration
    # Filtrate compartment (index 8) represents the concentration in kidney filtrate
    # This is more physiologically meaningful than calculating from cumulative excretion
    filtrate_conc = [u[8] for u in sol.u] ./ VFil  # μg/L
    
    # Convert from μg/L to ng/ml (1 μg/L = 1 ng/ml)
    concentrations["CUrine"] = filtrate_conc
    
    return concentrations
end

# Function to run PBK model for a specific study
function run_pbk_study(Tm, Kt, exposure_dose, exposure_time, measurement_times, study_name="", BW=70.0)

    if study_name == "Olsen"
        # Special case for Olsen: Two-phase exposure (exposure + elimination)
        # Phase 1: 10 years continuous exposure (0 to 3650 days)
        # Phase 2: Elimination phase (3650 to 3650+max_measurement_time days)
        
        elimination_start = exposure_time  # 10 years = 3650 days
        total_time = elimination_start + maximum(measurement_times)  # Total simulation time
        
        # Measurement times are already relative to elimination start
        # Convert to absolute times (relative to start of exposure)
        absolute_measurement_times = elimination_start .+ measurement_times
        
        user_input = (
            BW = BW,
            substance = "PFBS",
            admin_dose = [0.0],  # No IV dose
            admin_time = [0.0],
            f_unabs = 0.0,
            ingestion = [exposure_dose * 70.0 , 0.0],  # Exposure on, then off
            ingestion_time = [0.01, elimination_start],  # Start exposure, then stop
            admin_type = "oral",
            exp_type = "biomonitoring"
        )
        
        # Time span covers both exposure and elimination phases
        tspan = (0.0, total_time)
        base_grid = collect(range(0.0, total_time, length=100))
        # Include observation times in simulation grid (converted to absolute times)
        absolute_measurement_times = elimination_start .+ measurement_times
        saveat_times = sort(unique(vcat(base_grid, absolute_measurement_times)))
        
    else
        # Standard continuous exposure for Zhou and He studies
        user_input = (
            BW = BW,
            substance = "PFBS",
            admin_dose = [0.0],  # No IV dose
            admin_time = [0.0],
            f_unabs = 0.0,
            ingestion = [exposure_dose * 70.0],  # Continuous exposure
            ingestion_time = [0.01],
            admin_type = "oral",
            exp_type = "biomonitoring"
        )
        
        tspan = (0.0, maximum(measurement_times))
        base_grid = collect(range(0.0, maximum(measurement_times), length=100))
        # Include observation times in simulation grid
        saveat_times = sort(unique(vcat(base_grid, measurement_times)))
    end
    
    # Create parameters with new partition coefficients and updated Tm, Kt
    constant_params = create_params_updated(user_input, new_partition_coeffs)
    # Update both Tm and Kt, keep PR at EDM value from new_partition_coeffs
    params = (constant_params..., Tm = ReverseDiff.value(Tm), Kt = ReverseDiff.value(Kt))
    # params = (constant_params..., Tm = Tm, Kt = 5.0)
    # Create initial conditions
    inits = create_inits(params)
    
    # Create events for dosing
    events = create_events(params)
    
    # Define ODE problem
    prob = ODEProblem(ode_func, inits, tspan, saveat=saveat_times, params, callback=events)
    
    # Solve ODE
    sol = solve(prob, CVODE_BDF(), saveat=saveat_times, reltol=1e-8, abstol=1e-10)
    
    # Extract concentrations with urine concentration calculation
    concentrations = extract_concentrations_with_urine(sol, params)
    
    # Return all concentrations and simulation times
    return concentrations, sol.t
end

# Prior parameters for Tm and Kt (LogNormal with CV=200%)
# Original values: Tm=6.1, Kt=5.0
function lognormal_params(μ, cv)
    σ = μ * cv  # Standard deviation
    # Convert to log-normal parameters
    μ_log = log(μ) - 0.5 * log(1 + (σ/μ)^2)
    σ_log = sqrt(log(1 + (σ/μ)^2))
    return μ_log, σ_log
end

μ_log_Tm, σ_log_Tm = lognormal_params(146, 0.1)  # CV = 10%
μ_log_Kt, σ_log_Kt = lognormal_params(5.0, 0.1)  # CV = 10%

println("\nPrior parameters:")
println("  Tm: LogNormal($(round(μ_log_Tm, digits=3)), $(round(σ_log_Tm, digits=3)))")
println("  Kt: LogNormal($(round(μ_log_Kt, digits=3)), $(round(σ_log_Kt, digits=3)))")
println("  PR: Fixed at EDM value from partition coefficient calculation")

# Combine all data for MCMC
study_data = [
    (name="Zhou", serum=zhou_serum, urine=zhou_urine, exposure=exposure_zhou, exp_time=exposure_time_zhou),
    (name="Olsen", serum=olsen_serum, urine=olsen_urine, exposure=exposure_olsen, exp_time=exposure_time_olsen), 
    (name="He", serum=he_serum, urine=he_urine, exposure=exposure_he, exp_time=exposure_time_he)
]

println("\nData summary:")
for study in study_data
    println("  $(study.name):")
    println("    Serum: $(nrow(study.serum)) time points")
    println("    Urine: $(nrow(study.urine)) time points")
    println("    Exposure: $(round(study.exposure, digits=1)) ng/kg-bw/day")
end

# MCMC Model Definition
@model function pbk_calibration(study_data)
    
    # Priors for parameters
    Tm ~ LogNormal(μ_log_Tm, σ_log_Tm)  # Maximum reabsorption rate (μg/h)
    Kt ~ LogNormal(μ_log_Kt, σ_log_Kt)  # Michaelis constant (μg/L)
    # PR is fixed at EDM value from partition coefficient calculation
    
    # Proportional error parameter
    # sigma ~ truncated(Normal(0, 0.1), 0, Inf)  # Proportional error for both serum and urine
    sigma ~ truncated(Normal(0.0, 1.0), 0, Inf)
    # Loop through each study
    # for study in study_data
    study = study_data[1]  # Only one study in this example
    # Extract measurement times and convert to arrays
    serum_times = study.serum.time_days
    urine_times = study.urine.time_days
    all_times = sort(unique(vcat(serum_times, urine_times)))
    
    # Run PBK model for this study
    try
        concentrations, sim_times = run_pbk_study(Tm, Kt, study.exposure, study.exp_time, all_times, study.name)
        
        # Find matching predictions for all time points
        serum_pred_vals = Float64[]
        urine_pred_vals = Float64[]
        
        # Match serum observations to predictions
        for i in 1:nrow(study.serum)
            # For Olsen, convert relative elimination times to absolute times
            if study.name == "Olsen"
                abs_time = study.exp_time + study.serum.time_days[i]  # elimination_start + relative_time
                time_idx = findfirst(x -> x ≈ abs_time, sim_times)
            else
                time_idx = findfirst(x -> x ≈ study.serum.time_days[i], sim_times)
            end
            
            # Now guaranteed to find a match since obs times are in simulation grid
            push!(serum_pred_vals, concentrations["CPlas"][time_idx])
        end
        
        # Match urine observations to predictions
        for i in 1:nrow(study.urine)
            # For Olsen, convert relative elimination times to absolute times
            if study.name == "Olsen"
                abs_time = study.exp_time + study.urine.time_days[i]  # elimination_start + relative_time
                time_idx = findfirst(x -> x ≈ abs_time, sim_times)
            else
                time_idx = findfirst(x -> x ≈ study.urine.time_days[i], sim_times)
            end
            
            # Now guaranteed to find a match since obs times are in simulation grid
            push!(urine_pred_vals, concentrations["CUrine"][time_idx])
        end
        
        # Study-specific likelihood with unique variable names
        if study.name == "Zhou"
            obs_serum_Zhou = study.serum.mean_ng_ml
            obs_urine_Zhou = study.urine.mean_ng_ml
            # Additive error model
            obs_serum_Zhou ~ arraydist(Normal.(serum_pred_vals, sigma))
            obs_urine_Zhou ~ arraydist(Normal.(urine_pred_vals, sigma))
            # Proportional error model (commented out)
            # obs_serum_Zhou ~ arraydist(LogNormal.(log.(serum_pred_vals), sigma))
            # obs_urine_Zhou ~ arraydist(LogNormal.(log.(urine_pred_vals), sigma))
        elseif study.name == "Olsen"  
            obs_serum_Olsen = study.serum.mean_ng_ml
            obs_urine_Olsen = study.urine.mean_ng_ml
            # Additive error model
            obs_serum_Olsen ~ arraydist(Normal.(serum_pred_vals, sigma))
            obs_urine_Olsen ~ arraydist(Normal.(urine_pred_vals, sigma))
            # Proportional error model (commented out)
            # obs_serum_Olsen ~ arraydist(LogNormal.(log.(serum_pred_vals), sigma))
            # obs_urine_Olsen ~ arraydist(LogNormal.(log.(urine_pred_vals), sigma))
        elseif study.name == "He"
            obs_serum_He = study.serum.mean_ng_ml
            obs_urine_He = study.urine.mean_ng_ml
            # Additive error model
            obs_serum_He ~ arraydist(Normal.(serum_pred_vals, sigma))
            obs_urine_He ~ arraydist(Normal.(urine_pred_vals, sigma))
            # Proportional error model (commented out)
            # obs_serum_He ~ arraydist(LogNormal.(log.(serum_pred_vals), sigma))
            # obs_urine_He ~ arraydist(LogNormal.(log.(urine_pred_vals), sigma))
        end
        
    catch e
        println("Error in ODE solving for $(study.name): $e")
        # Add penalty for failed ODE solutions
        Turing.@addlogprob! -1000
    end
    # end
end

println("\n" * "="^80)
println("STARTING MCMC SAMPLING")
println("="^80)

# MCMC Settings
n_chains = 4
n_iterations = 500
n_burnin = 250

println("MCMC settings:")
println("  Chains: $n_chains")
println("  Iterations per chain: $n_iterations") 
println("  Burn-in: $n_burnin")

# Sample from the model
println("\nStarting parallel MCMC sampling...")
# chains = sample(pbk_calibration(study_data), NUTS(), MCMCThreads(), n_iterations, n_chains)
chains = sample(pbk_calibration(study_data), NUTS(adtype=AutoReverseDiff(false)), MCMCThreads(), n_iterations, n_chains)
# chains = sample(pbk_calibration(study_data), NUTS(), MCMCThreads(), n_iterations, n_chains)

println("\nMCMC sampling completed!")

# Save results
println("Saving MCMC results...")
using Serialization
serialize("PFBS_calibration_chains.jls", chains)

# Also save parameter samples as CSV for easier access
println("Saving parameter samples as CSV...")
try
    tm_samples = vec(Array(chains[:Tm]))
    kt_samples = vec(Array(chains[:Kt]))
    sigma_samples = vec(Array(chains[:sigma]))

    param_samples_df = DataFrame(
        Tm = tm_samples,
        Kt = kt_samples,
        sigma = sigma_samples
    )

    CSV.write("PFBS_parameter_samples.csv", param_samples_df)
    println("Parameter samples saved successfully.")
catch e
    println("Error saving parameter samples: $e")
end

# Generate diagnostics
println("\n" * "="^80)
println("DIAGNOSTICS")
println("="^80)

# Print detailed summary statistics
println("Detailed MCMC Summary:")
display(chains)

# Calculate and display R-hat values with details
rhat_values = rhat(chains)
println("\nR-hat Convergence Diagnostics:")
rhat_df = DataFrame(rhat_values)
for row in eachrow(rhat_df)
    status = row.rhat < 1.1 ? "✓ Converged" : "⚠ Poor convergence"
    println("  $(row.parameters): $(round(row.rhat, digits=4)) ($status)")
end

# Calculate and display Effective Sample Size (ESS)
ess_values = ess(chains)
println("\nEffective Sample Size (ESS) Diagnostics:")
ess_df = DataFrame(ess_values)
total_samples = n_iterations * n_chains
for row in eachrow(ess_df)
    ess_ratio = row.ess / total_samples
    status = if row.ess > 400
        "✓ Good (>400)"
    elseif row.ess > 100
        "⚠ Moderate (>100)"
    else
        "❌ Poor (<100)"
    end
    println("  $(row.parameters): $(round(row.ess, digits=0)) samples ($(round(ess_ratio*100, digits=1))% efficiency) $status")
end

# Extract and display key parameter statistics
println("\nPrior vs Posterior Comparison:")
key_params = [:Tm, :Kt, :sigma]

# Calculate prior statistics for LogNormal distributions
prior_mean_Tm = exp(μ_log_Tm + 0.5 * σ_log_Tm^2)
prior_var_Tm = (exp(σ_log_Tm^2) - 1) * exp(2*μ_log_Tm + σ_log_Tm^2)
prior_sd_Tm = sqrt(prior_var_Tm)

prior_mean_Kt = exp(μ_log_Kt + 0.5 * σ_log_Kt^2)
prior_var_Kt = (exp(σ_log_Kt^2) - 1) * exp(2*μ_log_Kt + σ_log_Kt^2)
prior_sd_Kt = sqrt(prior_var_Kt)

# Prior information
prior_info = Dict(
    :Tm => (mean = prior_mean_Tm, sd = prior_sd_Tm, description = "LogNormal($(round(μ_log_Tm, digits=3)), $(round(σ_log_Tm, digits=3)))"),
    :Kt => (mean = prior_mean_Kt, sd = prior_sd_Kt, description = "LogNormal($(round(μ_log_Kt, digits=3)), $(round(σ_log_Kt, digits=3)))"),
    :sigma => (mean = "N/A", sd = "N/A", description = "truncated(Normal(0, 1.0), 0, Inf)")
)

for param in key_params
    if param in names(chains)
        chain_values = vec(Array(chains[param]))
        post_mean = mean(chain_values)
        post_sd = std(chain_values)
        q025 = quantile(chain_values, 0.025)
        q975 = quantile(chain_values, 0.975)
        
        println("  $param:")
        println("    Prior: $(prior_info[param].description)")
        if prior_info[param].mean != "N/A"
            println("    Prior mean ± SD: $(round(prior_info[param].mean, digits=3)) ± $(round(prior_info[param].sd, digits=3))")
        end
        println("    Posterior mean ± SD: $(round(post_mean, digits=3)) ± $(round(post_sd, digits=3))")
        println("    Posterior 95% CI: [$(round(q025, digits=3)), $(round(q975, digits=3))]")
        println()
    end
end

# Generate trace plots
println("\nGenerating diagnostic plots...")

# Create individual trace plots with better formatting
p1 = plot(chains[:,:Tm,:], 
          title="Tm Trace Plot", 
          titlefontsize=14, 
          guidefontsize=12, 
          tickfontsize=10,
          ylabel="Tm (μg/h)",
          xlabel="Iteration")

p2 = plot(chains[:,:Kt,:], 
          title="Kt Trace Plot", 
          titlefontsize=14, 
          guidefontsize=12, 
          tickfontsize=10,
          ylabel="Kt (μg/L)",
          xlabel="Iteration")

p3 = plot(chains[:,:sigma,:], 
          title="σ Trace Plot", 
          titlefontsize=14, 
          guidefontsize=12, 
          tickfontsize=10,
          ylabel="σ",
          xlabel="Iteration")

# Combine plots with better sizing and margins
trace_plot = plot(p1, p2, p3, 
                 layout=(1,3), 
                 size=(1400,500), 
                 left_margin=5Plots.mm, 
                 right_margin=5Plots.mm, 
                 top_margin=10Plots.mm, 
                 bottom_margin=10Plots.mm,
                 dpi=300)
savefig(trace_plot, "PFBS_calibration_trace_plots.png")

# Create density plots with better formatting - show individual chains
p1_dens = plot(title="Tm Posterior Distribution", 
               titlefontsize=14, 
               guidefontsize=12, 
               tickfontsize=10,
               xlabel="Tm (μg/h)",
               ylabel="Density")

# Plot density for each chain separately
for chain_id in 1:n_chains
    tm_chain = chains[:, :Tm, chain_id]
    density!(p1_dens, tm_chain, 
            label="Chain $chain_id",
            linewidth=2,
            alpha=0.7)
end

p2_dens = plot(title="Kt Posterior Distribution", 
               titlefontsize=14, 
               guidefontsize=12, 
               tickfontsize=10,
               xlabel="Kt (μg/L)",
               ylabel="Density")

# Plot density for each chain separately  
for chain_id in 1:n_chains
    kt_chain = chains[:, :Kt, chain_id]
    density!(p2_dens, kt_chain,
            label="Chain $chain_id", 
            linewidth=2,
            alpha=0.7)
end

p3_dens = plot(title="σ Posterior Distribution", 
               titlefontsize=14, 
               guidefontsize=12, 
               tickfontsize=10,
               xlabel="σ",
               ylabel="Density")

# Plot density for each chain separately  
for chain_id in 1:n_chains
    sigma_chain = chains[:, :sigma, chain_id]
    density!(p3_dens, sigma_chain,
            label="Chain $chain_id", 
            linewidth=2,
            alpha=0.7)
end

# Combine density plots with better sizing and margins  
density_plot = plot(p1_dens, p2_dens, p3_dens, 
                   layout=(1,3), 
                   size=(1400,500),
                   left_margin=5Plots.mm, 
                   right_margin=5Plots.mm,
                   top_margin=10Plots.mm, 
                   bottom_margin=10Plots.mm,
                   dpi=300)
savefig(density_plot, "PFBS_calibration_density_plots.png")

println("\n" * "="^80)
println("CALIBRATION COMPLETE")
println("="^80)
println("Results saved to:")
println("  - PFBS_calibration_chains.jls (MCMC chains)")
println("  - PFBS_calibration_trace_plots.png (trace plots)")
println("  - PFBS_calibration_density_plots.png (posterior densities)")
println("="^80)