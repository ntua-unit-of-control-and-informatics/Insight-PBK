using Turing
using DifferentialEquations
using DataFrames
using CSV
using MCMCChains
using StatsPlots
using Distributions
using LinearAlgebra
using Random

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
    
    # Update partition coefficients with new values
    params = merge(params, (
        PL = new_PCs["liver"],
        PF = new_PCs["adipose"], 
        PK = new_PCs["kidney"],
        PG = new_PCs["gut"],
        PLu = new_PCs["lung"],
        PB = new_PCs["brain"],
        PR = new_PCs["rest_of_body"]
    ))
    
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
    
    # Calculate urine concentration from cumulative excretion
    # Urine production rate: 60 ml/hour = 1.44 L/day
    urine_production_rate = 1.44  # L/day
    
    # Extract cumulative urine amounts (μg)
    AUrine_cumulative = [u[10] for u in sol.u]  # μg
    
    # Calculate urine concentrations
    urine_concentrations = Float64[]
    
    for i in 1:length(sol.t)
        if i == 1
            # First time point: concentration = cumulative amount / total volume produced
            time_elapsed = sol.t[i]  # days from start
            total_volume = urine_production_rate * time_elapsed  # L
            
            if total_volume > 0
                urine_conc = AUrine_cumulative[i] / total_volume  # μg/L
            else
                urine_conc = 0.0
            end
        else
            # Subsequent time points: incremental concentration
            time_interval = sol.t[i] - sol.t[i-1]  # days
            volume_interval = urine_production_rate * time_interval  # L
            amount_interval = AUrine_cumulative[i] - AUrine_cumulative[i-1]  # μg
            
            if volume_interval > 0
                urine_conc = amount_interval / volume_interval  # μg/L
            else
                urine_conc = 0.0
            end
        end
        
        # Convert from μg/L to ng/ml (1 μg/L = 1 ng/ml)
        push!(urine_concentrations, urine_conc)
    end
    
    # Add urine concentration to results
    concentrations["CUrine"] = urine_concentrations
    
    return concentrations
end

# Function to run PBK model for a specific study
function run_pbk_study(Tm, Kt, exposure_dose, exposure_time, measurement_times, study_name="", BW=70.0)
    
    if study_name == "Olsen"
        # Special case for Olsen: Two-phase exposure (exposure + elimination)
        # Phase 1: 10 years continuous exposure
        # Phase 2: 180 days elimination (exposure = 0)
        
        elimination_start = exposure_time  # 10 years = 3650 days
        total_time = elimination_start + maximum(measurement_times)  # Total simulation time
        
        # Measurement times are relative to start of elimination phase
        # Convert to absolute times (relative to start of exposure)
        absolute_measurement_times = elimination_start .+ measurement_times
        
        user_input = (
            BW = BW,
            substance = "PFBS",
            admin_dose = [0.0],  # No IV dose
            admin_time = [0.0],
            f_unabs = 0.0,
            ingestion = [exposure_dose, 0.0],  # Exposure on, then off
            ingestion_time = [0.0, elimination_start],  # Start exposure, then stop
            admin_type = "oral",
            exp_type = "biomonitoring"
        )
        
        # Time span covers both exposure and elimination phases
        tspan = (0.0, total_time)
        saveat_times = vcat([0.0], absolute_measurement_times)
        
    else
        # Standard continuous exposure for Zhou and He studies
        user_input = (
            BW = BW,
            substance = "PFBS",
            admin_dose = [0.0],  # No IV dose
            admin_time = [0.0],
            f_unabs = 0.0,
            ingestion = [exposure_dose],  # Continuous exposure
            ingestion_time = [0.0],
            admin_type = "oral",
            exp_type = "biomonitoring"
        )
        
        tspan = (0.0, maximum(measurement_times))
        saveat_times = measurement_times
    end
    
    # Create parameters with new partition coefficients and updated Tm, Kt
    params = create_params_updated(user_input, new_partition_coeffs)
    params = merge(params, (Tm = Tm, Kt = Kt))
    
    # Create initial conditions
    inits = create_inits(params)
    
    # Create events for dosing
    events = create_events(params)
    
    # Define ODE problem
    prob = ODEProblem(ode_func, inits, tspan, params, callback=events)
    
    # Solve ODE
    sol = solve(prob, Tsit5(), saveat=saveat_times, reltol=1e-8, abstol=1e-10)
    
    # Extract concentrations with urine concentration calculation
    concentrations = extract_concentrations_with_urine(sol, params)
    
    # For Olsen study, extract only the elimination phase concentrations
    if study_name == "Olsen"
        # Skip the first time point (t=0 before elimination) and return elimination phase
        elimination_concentrations = Dict()
        for (key, values) in concentrations
            elimination_concentrations[key] = values[2:end]  # Skip first point
        end
        return elimination_concentrations
    else
        return concentrations
    end
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

μ_log_Tm, σ_log_Tm = lognormal_params(6.1, 2.0)  # CV = 200%
μ_log_Kt, σ_log_Kt = lognormal_params(5.0, 2.0)  # CV = 200%

println("\nPrior parameters:")
println("  Tm: LogNormal($(round(μ_log_Tm, digits=3)), $(round(σ_log_Tm, digits=3)))")
println("  Kt: LogNormal($(round(μ_log_Kt, digits=3)), $(round(σ_log_Kt, digits=3)))")

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
    
    # Proportional error parameters
    sigma_serum ~ Normal(0.0, sigma_serum^2)  # Proportional error for serum
    sigma_urine ~ Normal(0.0, sigma_urine^2)  # Proportional error for urine
    
    # Loop through each study
    for study in study_data
        
        # Extract measurement times and convert to arrays
        serum_times = study.serum.time_days
        urine_times = study.urine.time_days
        all_times = sort(unique(vcat(serum_times, urine_times)))
        
        # Run PBK model for this study
        try
            concentrations = run_pbk_study(Tm, Kt, study.exposure, study.exp_time, all_times, study.name)
            
            # Compare serum predictions to observations
            for i in 1:nrow(study.serum)
                time_idx = findfirst(x -> x ≈ study.serum.time_days[i], all_times)
                if !isnothing(time_idx)
                    pred_serum = concentrations["CPlas"][time_idx]
                    obs_serum = study.serum.mean_ng_ml[i]
                    # Likelihood with proportional error
                    # Cobs = Cpred*exp(epsilon)
                    # epsilon ~ Normal(0, sigma_serum^2)
                    # log(C_obs) ~ Normal(log(Cpred), sigma_serum^2)
                    log(obs_serum) ~ Normal(log(pred_serum), sigma_serum)
                end
            end
            
            # Compare urine predictions to observations  
            for i in 1:nrow(study.urine)
                time_idx = findfirst(x -> x ≈ study.urine.time_days[i], all_times)
                if !isnothing(time_idx)
                    pred_urine = concentrations["CUrine"][time_idx]  # ng/ml
                    obs_urine = study.urine.mean_ng_ml[i]
                    # Likelihood with proportional error
                    # Curine_obs = Curine_pred*exp(epsilon)
                    # epsilon ~ Normal(0, sigma_urine^2)
                    # log(Curine_obs) ~ Normal(log(Curine_pred), sigma_urine^2)
                    if pred_urine > 0 && obs_urine > 0
                        log(obs_urine) ~ Normal(log(pred_urine), sigma_urine)
                    end
                end
            end
            
        catch e
            println("Error in ODE solving for $(study.name): $e")
            # Add penalty for failed ODE solutions
            Turing.@addlogprob! -1000
        end
    end
end

println("\n" * "="^80)
println("STARTING MCMC SAMPLING")
println("="^80)

# MCMC Settings
n_chains = 4
n_iterations = 10000
n_burnin = 5000

println("MCMC settings:")
println("  Chains: $n_chains")
println("  Iterations per chain: $n_iterations") 
println("  Burn-in: $n_burnin")

# Sample from the model
println("\nStarting parallel MCMC sampling...")
chains = sample(pbk_calibration(study_data), NUTS(), MCMCThreads(), n_iterations, n_chains)

println("\nMCMC sampling completed!")

# Save results
println("Saving MCMC results...")
write("PFBS_calibration_chains.jls", chains)

# Generate diagnostics
println("\n" * "="^80)
println("DIAGNOSTICS")
println("="^80)

# Print summary statistics
println(chains)

# Calculate R-hat values
rhat_values = rhat(chains)
println("\nR-hat values:")
for (param, value) in pairs(rhat_values)
    println("  $param: $(round(value, digits=4))")
end

# Generate trace plots
println("\nGenerating diagnostic plots...")
p1 = plot(chains[:Tm], title="Tm Trace Plot")
p2 = plot(chains[:Kt], title="Kt Trace Plot") 
p3 = plot(chains[:σ_serum], title="σ_serum Trace Plot")
p4 = plot(chains[:σ_urine], title="σ_urine Trace Plot")

trace_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(800,600))
savefig(trace_plot, "PFBS_calibration_trace_plots.png")

# Density plots
p1_dens = density(chains[:Tm], title="Tm Posterior")
p2_dens = density(chains[:Kt], title="Kt Posterior")
p3_dens = density(chains[:σ_serum], title="σ_serum Posterior") 
p4_dens = density(chains[:σ_urine], title="σ_urine Posterior")

density_plot = plot(p1_dens, p2_dens, p3_dens, p4_dens, layout=(2,2), size=(800,600))
savefig(density_plot, "PFBS_calibration_density_plots.png")

# Extract posterior means
posterior_means = mean(chains)
println("\nPosterior means:")
for (param, value) in pairs(posterior_means)
    println("  $param: $(round(value, digits=4))")
end

println("\n" * "="^80)
println("CALIBRATION COMPLETE")
println("="^80)
println("Results saved to:")
println("  - PFBS_calibration_chains.jls (MCMC chains)")
println("  - PFBS_calibration_trace_plots.png (trace plots)")
println("  - PFBS_calibration_density_plots.png (posterior densities)")
println("="^80)