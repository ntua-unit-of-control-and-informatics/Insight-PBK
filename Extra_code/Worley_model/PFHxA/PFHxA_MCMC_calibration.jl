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
using Printf

# Include Worley PBK model functions
include("../PBK_model/julia_code.jl")

"""
Worley Model Parameter Calibration using MCMC

This script calibrates Vmax_baso (basolateral transporter maximum velocity)
of the Worley PBK model using synthetic steady-state plasma concentration data.

Approach:
- Individual-level calibration for each subject
- Synthetic observation: CA = 100 μg/L at 10 years exposure
- Exposure data from intake estimation (ng/kg-bw/day)
- LogNormal prior for Vmax_baso (CV to be specified)
- Proportional error model for plasma predictions
- 4 parallel chains, iterations to be specified
"""

println("="^80)
println("WORLEY PBK MODEL PARAMETER CALIBRATION - PFHxA")
println("="^80)

# Load intake estimation data
println("\nLoading intake estimation data...")
intake_data = CSV.read("Worley_model/PFHxA/Data/intake_estimation_results.csv", DataFrame)

# Load partition coefficients
println("\nLoading partition coefficients...")
partition_data = CSV.read("Worley_model/Worley_partition_coefficients_results.csv", DataFrame)

# Extract PFHxA partition coefficients
pfhxA_partitions = filter(row -> row.Compound == "PFHxA", partition_data)
PL_worley = pfhxA_partitions[pfhxA_partitions.Tissue .== "liver", :Partition_Coefficient][1]
PK_worley = pfhxA_partitions[pfhxA_partitions.Tissue .== "kidney", :Partition_Coefficient][1]
PR_worley = pfhxA_partitions[pfhxA_partitions.Tissue .== "rest of body", :Partition_Coefficient][1]

println("  PFHxA Partition Coefficients:")
println("    Liver:blood (PL) = $(round(PL_worley, digits=4))")
println("    Kidney:blood (PK) = $(round(PK_worley, digits=4))")
println("    Rest of body:blood (PR) = $(round(PR_worley, digits=4))")

println("  Total subjects: $(nrow(intake_data))")
println("  Studies: $(unique(intake_data.Study))")

# Create synthetic observation data
# Each subject has one observation: CA = 100 μg/L at 10 years
println("\nCreating synthetic observation data...")
observation_time = 10.0  # years
observation_concentration = 10.0  # μg/L

# Create observation dataframe
obs_data = DataFrame(
    Subject_ID = 1:nrow(intake_data),
    Study = intake_data.Study,
    Half_life_days = intake_data.Half_life_days,
    Intake_ng_kg_bw_day = intake_data.Estimated_Intake_ng_kg_bw_day,
    Observation_time_years = fill(observation_time, nrow(intake_data)),
    Observed_CA_ug_L = fill(observation_concentration, nrow(intake_data))
)

println("  Observation time: $observation_time years")
println("  Observed plasma concentration: $observation_concentration μg/L")
println("  Number of subjects: $(nrow(obs_data))")

# Summary by study
println("\nData summary by study:")
for study in unique(obs_data.Study)
    study_data = filter(row -> row.Study == study, obs_data)
    n_subjects = nrow(study_data)
    mean_intake = mean(study_data.Intake_ng_kg_bw_day)
    println("  $study:")
    println("    Subjects: $n_subjects")
    println("    Mean intake: $(round(mean_intake, digits=1)) ng/kg-bw/day")
end

# Function to run Worley PBK model for a specific subject
function run_pbk_worley(RAFbaso, RAFapi, exposure_dose_ng_kg_day, exposure_time_years, measurement_times_years, BW=70.0)
    """
    Run Worley PBK model with continuous exposure

    Parameters:
    - RAFbaso: Relative activity factor for basolateral transporters
    - RAFapi: Relative activity factor for apical transporters
    - exposure_dose_ng_kg_day: Daily intake (ng/kg-bw/day)
    - exposure_time_years: Duration of exposure (years)
    - measurement_times_years: Times to extract predictions (years)
    - BW: Body weight (kg)

    Returns:
    - CA values at measurement times (μg/L)
    - Simulation times (years)
    """

    # Convert intake from ng/kg-bw/day to μg/year for the model
    # ng/kg-bw/day × kg × (1 μg / 1000 ng) × (365 day / year) = μg/year
    intake_ug_per_year = exposure_dose_ng_kg_day * BW * (1.0/1000.0) * 365.0

    # Set up user input for Worley model
    user_input = (
        BW = BW,
        substance = "PFHxA",
        ingestion = [intake_ug_per_year],  # Continuous exposure (μg/year)
        ingestion_time = [0.001],  # Start at time 0
        admin_dose = 0.0,
        admin_time = 0.0,
        admin_type = "oral",
        exp_type = "continuous",
        time_scale = "years"
    )

    # Create base parameters with calibrated RAF values
    base_params = create_params_worley(user_input, RAF_baso=RAFbaso, RAF_api=RAFapi)

    # Update with loaded partition coefficients
    params = (base_params..., PL = PL_worley, PK = PK_worley, PR = PR_worley)

    # Create initial conditions
    inits = create_inits_worley(params)

    # Create events
    events = create_events_worley(params)

    # Time span
    tspan = (0.0, maximum(measurement_times_years))
    saveat_times = sort(unique(vcat([0.0], measurement_times_years)))

    # Define and solve ODE problem
    prob = ODEProblem(ode_func_worley, inits, tspan, params, callback=events)
    sol = solve(prob, Rodas5P(), saveat=saveat_times,
                reltol=1e-6, abstol=1e-8)

    # Extract CA (total plasma concentration) from solution
    # CA is calculated in the ODE system but we need to extract it
    CA_values = []

    for u in sol.u
        Aplas_free = u[18]  # Free plasma amount (μg)
        VPlas = params.VPlas  # Plasma volume (L)
        Free = params.Free  # Free fraction

        CA_free = Aplas_free / VPlas  # Free concentration (μg/L)
        CA = CA_free / Free  # Total concentration (μg/L)

        push!(CA_values, CA)
    end

    return CA_values, sol.t
end

# Helper function to create LogNormal prior parameters
function lognormal_params(μ, cv)
    """
    Convert mean and CV to LogNormal distribution parameters

    Parameters:
    - μ: Mean value
    - cv: Coefficient of variation (e.g., 0.5 for 50% CV)

    Returns:
    - μ_log: LogNormal location parameter
    - σ_log: LogNormal scale parameter
    """
    σ = μ * cv  # Standard deviation
    # Convert to log-normal parameters
    μ_log = log(μ) - 0.5 * log(1 + (σ/μ)^2)
    σ_log = sqrt(log(1 + (σ/μ)^2))
    return μ_log, σ_log
end

# Prior parameters for RAF values
# Default values: RAFbaso = 1.0, RAFapi = 0.0007
RAFbaso_prior_mean = 1.0  # Prior mean for RAFbaso (dimensionless)
RAFbaso_prior_cv = 1.0  # CV for prior (e.g., 0.5 for 50% CV)

RAFapi_prior_mean = 0.0007  # Prior mean for RAFapi (dimensionless)
RAFapi_prior_cv = 1.0  # CV for prior (e.g., 0.5 for 50% CV)

μ_log_RAFbaso, σ_log_RAFbaso = lognormal_params(RAFbaso_prior_mean, RAFbaso_prior_cv)
μ_log_RAFapi, σ_log_RAFapi = lognormal_params(RAFapi_prior_mean, RAFapi_prior_cv)

println("\n" * "="^80)
println("PRIOR PARAMETERS")
println("="^80)
println("RAFbaso prior mean: $RAFbaso_prior_mean (dimensionless)")
println("RAFbaso prior CV: $(RAFbaso_prior_cv*100)%")
println("  LogNormal parameters: μ=$(round(μ_log_RAFbaso, digits=3)), σ=$(round(σ_log_RAFbaso, digits=3))")
println()
println("RAFapi prior mean: $RAFapi_prior_mean (dimensionless)")
println("RAFapi prior CV: $(RAFapi_prior_cv*100)%")
println("  LogNormal parameters: μ=$(round(μ_log_RAFapi, digits=3)), σ=$(round(σ_log_RAFapi, digits=3))")

# Prepare data for MCMC
# Create vectors for vectorized sampling
n_subjects = nrow(obs_data)
log_obs = log.(obs_data.Observed_CA_ug_L)  # Log-transform observations

println("\n" * "="^80)
println("MCMC DATA PREPARATION")
println("="^80)
println("Number of subjects: $n_subjects")
println("Log-transformed observations: $(length(log_obs))")

# MCMC Model Definition
@model function worley_calibration(obs_data, log_obs)
    """
    Turing model for Worley PBK calibration

    Estimates:
    - RAFbaso: Relative activity factor for basolateral transporters
    - RAFapi: Relative activity factor for apical transporters
    - sigma: Proportional error for plasma predictions
    """

    # Priors
    # RAFbaso ~ LogNormal(μ_log_RAFbaso, σ_log_RAFbaso)
    RAFbaso = 1.0 # for identifiability
    RAFapi ~ LogNormal(μ_log_RAFapi, σ_log_RAFapi)
    sigma ~ truncated(Normal(0.0, 0.5), 0, Inf)  # Proportional error for plasma

    # Loop through each subject
    for i in 1:nrow(obs_data)
        subject = obs_data[i, :]

        # Extract subject-specific data
        exposure_dose = subject.Intake_ng_kg_bw_day
        exposure_time = subject.Observation_time_years
        measurement_times = [exposure_time]

        # Run PBK model for this subject
        try
            CA_predictions, sim_times = run_pbk_worley(
                RAFbaso,
                RAFapi,
                exposure_dose,
                exposure_time,
                measurement_times
            )

            # Extract prediction at observation time
            # (should be the last time point since we only have one measurement)
            CA_pred = CA_predictions[end]

            # Likelihood: log-normal error model
            log_obs[i] ~ Normal(log(CA_pred), sigma)

        catch e
            println("Error in ODE solving for subject $i: $e")
            # Add penalty for failed ODE solutions
            Turing.@addlogprob! -1000
        end
    end
end

println("\n" * "="^80)
println("MCMC MODEL DEFINED")
println("="^80)

# MCMC Settings (TO BE CONFIGURED)
n_chains = 4
n_iterations = 1000  # Adjust as needed
n_burnin = 500

println("\n" * "="^80)
println("MCMC SETTINGS")
println("="^80)
println("Chains: $n_chains")
println("Iterations per chain: $n_iterations")
println("Burn-in: $n_burnin")

println("\n" * "="^80)
println("READY TO START MCMC SAMPLING")

# Uncomment below to run MCMC sampling
println("\n" * "="^80)
println("STARTING MCMC SAMPLING")
println("="^80)

chains = sample(worley_calibration(obs_data, log_obs), NUTS(), MCMCThreads(), n_iterations, n_chains, discard_initial=n_burnin)

println("\nMCMC sampling completed!")

# Save results
println("\nSaving MCMC results...")
using Serialization
serialize("Worley_model/PFHxA/Worley_PFHxA_calibration_chains.jls", chains)

# Generate diagnostics
println("\n" * "="^80)
println("DIAGNOSTICS")
println("="^80)

# Extract and display parameter statistics
chain_array = Array(chains)
param_names = String.(names(chains))

# RAFapi statistics
RAFapi_idx = findfirst(==("RAFapi"), param_names)
if !isnothing(RAFapi_idx)
    RAFapi_samples = chain_array[:, RAFapi_idx]
    RAFapi_mean = mean(RAFapi_samples)
    RAFapi_std = std(RAFapi_samples)
    println("\nRAFapi posterior: mean=$(@sprintf("%.2e", RAFapi_mean)), std=$(@sprintf("%.2e", RAFapi_std))")
end

# sigma statistics
sigma_idx = findfirst(==("sigma"), param_names)
if !isnothing(sigma_idx)
    sigma_samples = chain_array[:, sigma_idx]
    sigma_mean = mean(sigma_samples)
    sigma_std = std(sigma_samples)
    println("sigma posterior: mean=$(@sprintf("%.2e", sigma_mean)), std=$(@sprintf("%.2e", sigma_std))")
end

println()
describe(chains)

# Generate diagnostic plots
println("\nGenerating diagnostic plots...")
chains_plot = plot(chains)
savefig(chains_plot, "Worley_model/PFHxA/Worley_PFHxA_calibration_diagnostic_plots.png")

println("\n" * "="^80)
println("CALIBRATION COMPLETE")
println("="^80)
println("Results saved to:")
println("  - Worley_model/PFHxA/Worley_PFHxA_calibration_chains.jls (MCMC chains)")
println("  - Worley_model/PFHxA/Worley_PFHxA_calibration_diagnostic_plots.png (trace plots and density plots)")
println("="^80)
