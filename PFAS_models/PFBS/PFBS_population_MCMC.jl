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
Population-level Worley PBK Model Parameter Calibration using MCMC

This script implements hierarchical Bayesian inference to estimate:
- mu_pop: Mean of log(RAFapi) across population
- sigma_pop: SD of log(RAFapi) across population
- RAFapi[i]: Subject-specific RAFapi values
- Km_apical: Michaelis constant for apical transporters
- sigma: Observation error (log scale)

The model transforms mu_pop and sigma_pop to obtain:
- Population mean RAFapi (original scale)
- Population SD RAFapi (original scale)
- Population CV

Data: Elimination phase concentration-time profiles for 6 subjects
"""

println("="^80)
println("POPULATION WORLEY PBK MODEL - PFBS")
println("="^80)

# Load partition coefficients
println("\nLoading partition coefficients...")
partition_data = CSV.read("PFAS_models/partition_coefficients_results.csv", DataFrame)
pfbs_partitions = filter(row -> row.Compound == "PFBS", partition_data)
PL = pfbs_partitions[pfbs_partitions.Tissue .== "liver", :Partition_Coefficient][1]
PK = pfbs_partitions[pfbs_partitions.Tissue .== "kidney", :Partition_Coefficient][1]
PR = pfbs_partitions[pfbs_partitions.Tissue .== "rest of body", :Partition_Coefficient][1]

println("  PFBS Partition Coefficients:")
println("    Liver:blood (PL) = $(round(PL, digits=4))")
println("    Kidney:blood (PK) = $(round(PK, digits=4))")
println("    Rest of body:blood (PR) = $(round(PR, digits=4))")

# Load data for all subjects
println("\nLoading subject data...")
n_subjects = 6
subject_data = []

for i in 1:n_subjects
    # Read CSV file
    data = CSV.read("PFAS_models/PFBS/digitized_data/subject_$i.csv", DataFrame)

    # Round time to remove decimals
    data.Time_days = round.(data.Time_days)

    # Rename concentration column for consistency
    rename!(data, Symbol("C_serum_ug/L") => :C_serum)

    # Pre-compute log-transformed observations
    data.log_C_serum = log.(data.C_serum)

    push!(subject_data, data)
    println("  Subject $i: $(nrow(data)) observations, time range: $(data.Time_days[1])-$(data.Time_days[end]) days")
end

# Background exposure parameter (to be set later)
background_dose_ng_kg_day = 0.93/70.0  # Low background exposure Haug et al. (2009)

println("\nBackground exposure: $background_dose_ng_kg_day ng/kg-bw/day")

# Function to run Worley PBK model for elimination phase
function run_pbk_elimination(RAFapi, Km_apical, initial_C_serum, observation_times_days,
                             background_dose_ng_kg_day, BW=70.0)
    """
    Run Worley PBK model for elimination phase with initial conditions

    Parameters:
    - RAFapi: Relative activity factor for apical transporters
    - Km_apical: Michaelis constant for apical transporters
    - initial_C_serum: Initial serum concentration (μg/L)
    - observation_times_days: Times to extract predictions (days)
    - background_dose_ng_kg_day: Background exposure (ng/kg-bw/day)
    - BW: Body weight (kg)

    Returns:
    - CA values at observation times (μg/L)
    - Time values (days)
    """

    # Convert background dose from ng/kg-bw/day to μg/day
    background_ug_per_day = background_dose_ng_kg_day * BW * (1.0/1000.0)

    # Set up user input for Worley model
    user_input = (
        BW = BW,
        substance = "PFBS",
        ingestion = [background_ug_per_day],  # Low background exposure
        ingestion_time = [1.0e-10],  # Small delay needed for callbacks
        admin_dose = 0.0,
        admin_time = 0.0,
        admin_type = "oral",
        exp_type = "continuous",
        time_scale = "days"
    )

    # Create base parameters with RAFapi
    base_params = create_params_worley(user_input, RAF_api=RAFapi)

    # Update with partition coefficients
    params = (base_params..., PL = PL, PK = PK, PR = PR, Km_apical = Km_apical)

    # Calculate initial conditions for tissue compartments
    # Using partition coefficients: C_tissue = P_tissue * C_serum
    # Initial amounts = C_tissue * V_tissue

    # Get tissue volumes from params
    VPlas = params.VPlas
    VL = params.VL
    VK = params.VK
    VR = params.VR
    Free = params.Free

    # Calculate initial free serum concentration
    C_serum_free_initial = initial_C_serum * Free

    # Calculate initial tissue concentrations (total)
    C_liver_initial = PL * initial_C_serum
    C_kidney_initial = PK * initial_C_serum
    C_rest_initial = PR * initial_C_serum

    # Calculate initial amounts (free)
    Aplas_free_initial = C_serum_free_initial * VPlas
    AL_initial = C_liver_initial * VL
    AK_initial = C_kidney_initial * VK
    AR_initial = C_rest_initial * VR

    # Create initial conditions vector (all compartments)
    # Order based on Worley model: [gut, liver, kidney, rest, filtrate, storage, plasma, ...]
    inits = create_inits_worley(params)

    # Update relevant compartments with calculated initial conditions
    # Order: [Agut, AL_free, ..., AR_free, ..., APTC (index 8), ..., Aplas_free (index 18)]
    inits[16] = AL_initial   # Liver free
    inits[8] = AK_initial   # PTC (Proximal Tubule Cells) - kidney initial condition
    inits[1] = AR_initial   # Rest of body free
    inits[18] = Aplas_free_initial  # Plasma free

    # Create events
    events = create_events_worley(params)

    # Time span
    tspan = (0.0, maximum(observation_times_days))
    saveat_times = sort(unique(observation_times_days))

    # Solve ODE
    prob = ODEProblem(ode_func_worley, inits, tspan, params, callback=events)
    sol = solve(prob, Rodas5P(), saveat=saveat_times, reltol=1e-6, abstol=1e-8)

    # Extract CA concentrations
    CA_values = []
    for u in sol.u
        Aplas_free = u[18]
        CA_free = Aplas_free / VPlas
        CA = CA_free / Free
        push!(CA_values, CA)
    end

    return CA_values, sol.t
end

# Prepare data for MCMC
println("\n" * "="^80)
println("PREPARING DATA FOR POPULATION MCMC")
println("="^80)

# Create 2D array for log-transformed observations
# Handle variable number of observations per subject by using NaN padding
# Exclude initial time point (t=0) from likelihood data
# Dimensions: [max_observations, n_subjects]
max_obs = maximum([sum(subject_data[i].Time_days .!= 0.0) for i in 1:n_subjects])
log_data_2d = fill(NaN, max_obs, n_subjects)

# Fill with log-transformed observations (excluding t=0)
for i in 1:n_subjects
    # Find indices where time != 0
    non_zero_indices = findall(subject_data[i].Time_days .!= 0.0)
    n_obs_likelihood = length(non_zero_indices)
    log_data_2d[1:n_obs_likelihood, i] = subject_data[i].log_C_serum[non_zero_indices]
end

println("  Log-data 2D array shape: $(size(log_data_2d))")
println("  Maximum observations per subject (excluding t=0): $max_obs")
for i in 1:n_subjects
    n_obs_total = nrow(subject_data[i])
    n_obs_likelihood = sum(subject_data[i].Time_days .!= 0.0)
    println("    Subject $i: $n_obs_total total, $n_obs_likelihood for likelihood (excluding t=0)")
end

# Prior parameters for population-level RAFapi
RAFapi_prior_mean = 0.0007 
RAFapi_prior_cv = 0.5
Km_apical_prior_mean = 52338.44 #1500000.0
Km_apical_prior_cv = 0.5

# Calculate LogNormal parameters for prior
function lognormal_params(mu, cv)
    sigma = mu * cv
    mu_log = log(mu) - 0.5 * log(1 + (sigma/mu)^2)
    sigma_log = sqrt(log(1 + (sigma/mu)^2))
    return mu_log, sigma_log
end

mu_log_prior, sigma_log_prior = lognormal_params(RAFapi_prior_mean, RAFapi_prior_cv)
mu_log_km_apical, sigma_log_km_apical = lognormal_params(Km_apical_prior_mean, Km_apical_prior_cv)

println("\nPopulation prior for RAFapi:")
println("  Prior mean: $RAFapi_prior_mean")
println("  Prior CV: $(RAFapi_prior_cv*100)%")
println("  LogNormal parameters: μ=$(round(mu_log_prior, digits=3)), σ=$(round(sigma_log_prior, digits=3))")

println("\nPopulation prior for Km_apical:")
println("  Prior mean: $Km_apical_prior_mean")
println("  Prior CV: $(Km_apical_prior_cv*100)%")
println("  LogNormal parameters: μ=$(round(mu_log_km_apical, digits=3)), σ=$(round(sigma_log_km_apical, digits=3))")

# Population MCMC Model
@model function population_worley_model(subject_data, log_data_2d)
    """
    Hierarchical Bayesian model for population PBK calibration

    Parameters:
    - subject_data: Vector of DataFrames with subject information
    - log_data_2d: 2D array [max_obs, n_subjects] with log-transformed observations (NaN-padded)

    Estimates (log-scale parameters):
    - mu_pop: Mean of log(RAFapi) - log-scale location parameter
    - sigma_pop: SD of log(RAFapi) - log-scale dispersion parameter
    - RAFapi[i]: Subject-specific RAFapi values ~ LogNormal(mu_pop, sigma_pop)
    - Km_apical: Michaelis constant for apical transporters (μg/L)
    - sigma: Observation error (log scale)
    """

    # Population-level hyperpriors (log-scale parameters)
    mu_pop ~ Normal(mu_log_prior, 1.0)  # Mean of log(RAFapi)
    sigma_pop ~ truncated(Normal(0.0, 0.5), 0.0, Inf)  # SD of log(RAFapi)

    # Population-level Km
    Km_apical ~ LogNormal(mu_log_km_apical, sigma_log_km_apical)

    # Observation error
    sigma ~ truncated(Normal(0.0, 0.5), 0.0, Inf)

    # Subject-level parameters
    RAFapi = Vector{Float64}(undef, n_subjects)

    # Loop through each subject
    for i in 1:n_subjects
        RAFapi[i] ~ LogNormal(mu_pop, sigma_pop)
        # Extract data for subject i
        data = subject_data[i]
        observation_times = data.Time_days
        initial_C_serum = data.C_serum[1]

        try
            # Run PBK model with all time points (including t=0)
            CA_pred_all, time_solution = run_pbk_elimination(
                RAFapi[i],
                Km_apical,
                initial_C_serum,
                observation_times,
                background_dose_ng_kg_day
            )

            # Keep only predictions at observation times and exclude t=0
            # Find which time_solution values are in observation_times
            matching_indices = findall(t -> t in observation_times, time_solution)
            CA_pred_matched = CA_pred_all[matching_indices]
            time_matched = time_solution[matching_indices]

            # Remove initial time point (t=0) from likelihood
            non_zero_indices = findall(time_matched .!= 0.0)
            CA_pred = CA_pred_matched[non_zero_indices] 
            n_obs_likelihood = length(non_zero_indices)

            # Likelihood: log-normal error model
            # Slice log_data_2d to only non-initial observations (excluding t=0)
            log_data_2d[1:n_obs_likelihood, i] ~ arraydist(Normal.(log.(CA_pred), sigma))

        catch e
            println("Error in ODE solving for subject $i: $e")
            # Add penalty for failed ODE solutions
            Turing.@addlogprob! -1000
        end
    end
end

println("\n" * "="^80)
println("POPULATION MCMC MODEL DEFINED")
println("="^80)

# MCMC Settings
n_chains = 4
n_iterations = 1000
n_burnin = 500

println("\nMCMC Settings:")
println("  Chains: $n_chains")
println("  Iterations per chain: $n_iterations")
println("  Burn-in: $n_burnin")

println("\n" * "="^80)
println("STARTING POPULATION MCMC SAMPLING")
println("="^80)

# Run MCMC sampling
chains = sample(population_worley_model(subject_data, log_data_2d),
                NUTS(), MCMCThreads(), n_iterations, n_chains,
                discard_initial=n_burnin)

println("\nMCMC sampling completed!")

# Save results
println("\nSaving MCMC results...")
using Serialization
serialize("PFAS_models/PFBS/PFBS_population_chains.jls", chains)

# Generate diagnostics
println("\n" * "="^80)
println("DIAGNOSTICS")
println("="^80)
println()
describe(chains)

# Calculate derived quantities: Population mean and SD on original scale
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
mean_RAFapi_ci = quantile(mean_RAFapi_samples, [0.025, 0.975])

sd_RAFapi_mean = mean(sd_RAFapi_samples)
sd_RAFapi_ci = quantile(sd_RAFapi_samples, [0.025, 0.975])

cv_RAFapi_mean = mean(cv_RAFapi_samples)
cv_RAFapi_ci = quantile(cv_RAFapi_samples, [0.025, 0.975])

println("\nPopulation Mean RAFapi:")
println("  Estimate: $(round(mean_RAFapi_mean, sigdigits=4))")
println("  95% CI: [$(round(mean_RAFapi_ci[1], sigdigits=4)), $(round(mean_RAFapi_ci[2], sigdigits=4))]")

println("\nPopulation SD RAFapi:")
println("  Estimate: $(round(sd_RAFapi_mean, sigdigits=4))")
println("  95% CI: [$(round(sd_RAFapi_ci[1], sigdigits=4)), $(round(sd_RAFapi_ci[2], sigdigits=4))]")

println("\nPopulation CV:")
println("  Estimate: $(round(cv_RAFapi_mean, digits=3))")
println("  95% CI: [$(round(cv_RAFapi_ci[1], digits=3)), $(round(cv_RAFapi_ci[2], digits=3))]")

println("\nLog-scale parameters:")
println("  mu_pop (mean of log(RAFapi)): $(round(mean(chains[:mu_pop]), digits=3))")
println("  sigma_pop (SD of log(RAFapi)): $(round(mean(chains[:sigma_pop]), digits=3))")
println("  Km_apical: $(round(mean(chains[:Km_apical]), sigdigits=4)) μg/L")

# Generate diagnostic plots
println("\n" * "="^80)
println("GENERATING DIAGNOSTIC PLOTS")
println("="^80)

chains_plot = plot(chains, size=(2000, 2500), dpi=300)
savefig(chains_plot, "PFAS_models/PFBS/PFBS_population_diagnostic_plots.png")
println("  Saved to: PFAS_models/PFBS/PFBS_population_diagnostic_plots.png")

println("\n" * "="^80)
println("POPULATION CALIBRATION COMPLETE")
println("="^80)
println("Results saved to:")
println("  - PFAS_models/PFBS/PFBS_population_chains.jls (MCMC chains)")
println("  - PFAS_models/PFBS/PFBS_population_diagnostic_plots.png (diagnostic plots)")
println("="^80)
