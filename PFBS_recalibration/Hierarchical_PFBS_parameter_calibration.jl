using Turing
using DifferentialEquations
using DataFrames
using CSV
using MCMCChains
using StatsPlots
using Plots
using Distributions
using LinearAlgebra
using Random
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
# Random.seed!(42)

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
exposure_time_he = 10.0*365  # User to fill

println("\n  Exposure times:")
println("    Zhou et al.: $exposure_time_zhou days")
println("    Olsen et al.: $exposure_time_olsen days")
println("    He et al.: $exposure_time_he days")

# Modified PBK parameter creation function with new partition coefficients
function create_params_updated(user_input, new_PCs)
    params = create_params(user_input)
    
    # Use EDM partition coefficients without scaling
    params = (params..., 
              PL = new_PCs["liver"],
              PF = new_PCs["adipose"], 
              PK = new_PCs["kidney"],
              PG = new_PCs["gut"],
              PLu = new_PCs["lung"],
              PB = new_PCs["brain"],
              PR = new_PCs["rest_of_body"])
    
    return params
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
            ingestion = [exposure_dose * 70.0 * 1e-03, 0.0],  # Exposure on, then off
            ingestion_time = [0.01, elimination_start],  # Start exposure, then stop
            admin_type = "oral",
            exp_type = "biomonitoring"
        )
        
        # Time span covers both exposure and elimination phases
        tspan = (0.0, total_time)
        # For Olsen: Only save at measurement times (converted to absolute times)
        # Add elimination start time to ensure proper transition modeling
        absolute_measurement_times = elimination_start .+ measurement_times
        saveat_times = sort(unique(vcat([0.01, elimination_start], absolute_measurement_times)))
        
    else
        # Standard continuous exposure for Zhou and He studies
        user_input = (
            BW = BW,
            substance = "PFBS",
            admin_dose = [0.0],  # No IV dose
            admin_time = [0.0],
            f_unabs = 0.0,
            ingestion = [exposure_dose * 70.0 * 1e-03],  # Continuous exposure
            ingestion_time = [0.01],
            admin_type = "oral",
            exp_type = "biomonitoring"
        )
        
        tspan = (0.0, maximum(measurement_times))
        # For Zhou/He: Only save at measurement times
        # Add initial time point to ensure proper initial conditions
        saveat_times = sort(unique(vcat([0.01], measurement_times)))
    end
    
    # Create parameters with new partition coefficients and updated Tm
    constant_params = create_params_updated(user_input, new_partition_coeffs)
    # Update Tm and Kt parameters
    params = (constant_params..., Tm = Tm, Kt = Kt)
    # Create initial conditions
    inits = create_inits(params)
    
    # Create events for dosing
    events = create_events(params)
    
    # Define ODE problem
    prob = ODEProblem(ode_func, inits, tspan, params, callback=events)
    
    # Solve ODE
    sol = solve(prob, Rodas5P(), saveat=saveat_times, reltol=1e-6, abstol=1e-8)

    # Extract concentrations using main function (includes physiological urine concentration)
    concentrations = extract_concentrations(sol, params)

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

# μ_log_Tm, σ_log_Tm = lognormal_params(7000.0, 1.0)  # CV = 100%
Kt_fixed = 15000.0  # Fixed Kt value

println("\nPrior parameters:")
# println("  Tm: LogNormal($(round(μ_log_Tm, digits=3)), $(round(σ_log_Tm, digits=3)))")
println("  Kt: Fixed at $(Kt_fixed) μg/L")
println("  Partition coefficients: Fixed at EDM calculated values")

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

# Create 3D log-transformed data array
# Dimensions: [observation_index, study_index, compartment_type]
# study_index: 1=Zhou, 2=Olsen, 3=He
# compartment_type: 1=serum, 2=urine
println("\nCreating 3D log-transformed data array...")
max_obs_serum = maximum([nrow(zhou_serum), nrow(olsen_serum), nrow(he_serum)])
max_obs_urine = maximum([nrow(zhou_urine), nrow(olsen_urine), nrow(he_urine)])
max_obs = max(max_obs_serum, max_obs_urine)

log_data = fill(NaN, max_obs, 3, 2)  # (max_obs, 3 studies, 2 compartments)

# Fill with log-transformed observations
log_data[1:nrow(zhou_serum), 1, 1] = log.(zhou_serum.mean_ng_ml)
log_data[1:nrow(zhou_urine), 1, 2] = log.(zhou_urine.mean_ng_ml)
log_data[1:nrow(olsen_serum), 2, 1] = log.(olsen_serum.mean_ng_ml)
log_data[1:nrow(olsen_urine), 2, 2] = log.(olsen_urine.mean_ng_ml)
log_data[1:nrow(he_serum), 3, 1] = log.(he_serum.mean_ng_ml)
log_data[1:nrow(he_urine), 3, 2] = log.(he_urine.mean_ng_ml)

println("  Log-data array shape: $(size(log_data))")
println("  Zhou serum obs: $(nrow(zhou_serum)), urine obs: $(nrow(zhou_urine))")
println("  Olsen serum obs: $(nrow(olsen_serum)), urine obs: $(nrow(olsen_urine))")
println("  He serum obs: $(nrow(he_serum)), urine obs: $(nrow(he_urine))")

# MCMC Model Definition
@model function pbk_calibration(study_data, log_data, Kt_fixed)

    # Priors for parameters
    μ_Tm ~ truncated(Normal(15000.0, 7500.0), 0, Inf)  # Hyperprior for log-mean of Tm
    σ_Tm ~ truncated(Normal(7500.0, 2500.0), 0, Inf)  # Hyperprior for log-sd of Tm
    μ_log_Tm, σ_log_Tm = lognormal_params(μ_Tm, σ_Tm/μ_Tm)
    Tm_zhou ~ LogNormal(μ_log_Tm, σ_log_Tm)  # Maximum reabsorption rate (μg/day)
    Tm_olsen ~ LogNormal(μ_log_Tm, σ_log_Tm)  # Maximum reabsorption rate (μg/day)
    Tm_he ~ LogNormal(μ_log_Tm, σ_log_Tm)  # Maximum reabsorption rate (μg/day)
    Kt = Kt_fixed  # Fixed Michaelis constant (μg/L)

    # Proportional error parameter
    sigma_serum ~ truncated(Normal(0.0, 0.5), 0, Inf)  #  error for serum
    sigma_urine ~ truncated(Normal(0.0, 0.5), 0, Inf)  #  error for urine

    # Loop through each study
    for (study_idx, study) in enumerate(study_data)
        # Extract measurement times and convert to arrays
        serum_times = study.serum.time_days
        urine_times = study.urine.time_days
        all_times = sort(unique(vcat(serum_times, urine_times)))
        if study.name == "Zhou"
            Tm = Tm_zhou
        elseif study.name == "Olsen"
            Tm = Tm_olsen
        elseif study.name == "He"
            Tm = Tm_he
        end
        # Run PBK model for this study
        try
            concentrations, sim_times = run_pbk_study(Tm, Kt, study.exposure, study.exp_time, all_times, study.name)
            # Find matching predictions for all time points (vectorized)
            # Convert time arrays to target times
            if study.name == "Olsen"
                # For Olsen: convert relative elimination times to absolute times
                target_serum_times = study.exp_time .+ study.serum.time_days
                target_urine_times = study.exp_time .+ study.urine.time_days
            else
                # For Zhou/He: use times directly
                target_serum_times = study.serum.time_days
                target_urine_times = study.urine.time_days
            end
            
            # Find all matching indices at once
            serum_indices = [findfirst(x -> x ≈ t, sim_times) for t in target_serum_times]
            urine_indices = [findfirst(x -> x ≈ t, sim_times) for t in target_urine_times]
            
            # Extract all predictions at once
            serum_pred_vals = concentrations["CPlas"][serum_indices]
            urine_pred_vals = concentrations["CUrine"][urine_indices]

            # Get number of observations for this study
            n_serum = length(serum_pred_vals)
            n_urine = length(urine_pred_vals)

            # Use log_data from function argument - Turing recognizes it as observed
            log_data[1:n_serum, study_idx, 1] ~ arraydist(Normal.(log.(serum_pred_vals), sigma_serum))
            log_data[1:n_urine, study_idx, 2] ~ arraydist(Normal.(log.(urine_pred_vals), sigma_urine))
            
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
n_iterations = 2000
n_burnin = 1000

println("MCMC settings:")
println("  Chains: $n_chains")
println("  Iterations per chain: $n_iterations") 
println("  Burn-in: $n_burnin")

# Sample from the model
println("\nStarting parallel MCMC sampling...")
chains = sample(pbk_calibration(study_data, log_data, Kt_fixed), NUTS(), MCMCThreads(), n_iterations, n_chains)


println("\nMCMC sampling completed!")

# Save results
println("Saving MCMC results...")
using Serialization
serialize("PFBS_calibration_chains.jls", chains)


# Generate diagnostics
println("\n" * "="^80)
println("DIAGNOSTICS")
println("="^80)

# Print detailed summary statistics using describe
describe(chains)
# Generate diagnostic plots
println("\nGenerating diagnostic plots...")
chains_plot = plot(chains)
savefig(chains_plot, "PFBS_calibration_diagnostic_plots.png")

println("\n" * "="^80)
println("CALIBRATION COMPLETE")
println("="^80)
println("Results saved to:")
println("  - PFBS_calibration_chains.jls (MCMC chains)")
println("  - PFBS_calibration_diagnostic_plots.png (trace plots and density plots)")
println("="^80)