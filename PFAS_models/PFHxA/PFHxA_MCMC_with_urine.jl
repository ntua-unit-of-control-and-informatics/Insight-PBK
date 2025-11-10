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
Single-Dose Worley PBK Model Parameter Calibration using MCMC with Plasma and Urine Data

This script implements Bayesian inference to estimate:
- RAFapi: Relative activity factor for apical transporters
- RAFbaso: Relative activity factor for basolateral transporters
- Km_apical: Michaelis constant for apical transporters (μg/L)
- PC_scale: Partition coefficient scaling factor
- sigma_plasma: Observation error for plasma (log scale)
- sigma_urine: Observation error for urine (log scale)

Data: Single oral dose (3.99 μg) with plasma concentration and cumulative urine excretion time-series
"""

println("="^80)
println("SINGLE-DOSE WORLEY PBK MODEL - PFHxA (PLASMA + URINE)")
println("="^80)

# Load partition coefficients
println("\nLoading partition coefficients...")
partition_data = CSV.read("PFAS_models/partition_coefficients_results.csv", DataFrame)
pfhxa_partitions = filter(row -> row.Compound == "PFHxA", partition_data)
PL = pfhxa_partitions[pfhxa_partitions.Tissue .== "liver", :Partition_Coefficient][1]
PK = pfhxa_partitions[pfhxa_partitions.Tissue .== "kidney", :Partition_Coefficient][1]
PR = pfhxa_partitions[pfhxa_partitions.Tissue .== "rest of body", :Partition_Coefficient][1]

println("  PFHxA Partition Coefficients:")
println("    Liver:blood (PL) = $(round(PL, digits=4))")
println("    Kidney:blood (PK) = $(round(PK, digits=4))")
println("    Rest of body:blood (PR) = $(round(PR, digits=4))")

# Load plasma data
println("\nLoading plasma data...")
plasma_data = CSV.read("PFAS_models/PFHxA/digitized_data/Abraham_2024_plasma_data.csv", DataFrame)
rename!(plasma_data, Symbol("C_plasma_ug/L") => :C_plasma)
plasma_data = dropmissing(plasma_data)
plasma_data.log_C_plasma = log.(plasma_data.C_plasma)

println("  Plasma data points: $(nrow(plasma_data))")
println("  Time range: $(minimum(plasma_data.Time_days)) - $(maximum(plasma_data.Time_days)) days")

# Load urine data
println("\nLoading urine data...")
urine_data = CSV.read("PFAS_models/PFHxA/digitized_data/Abraham_2024_urine_data_processed.csv", DataFrame)
urine_data = dropmissing(urine_data)
urine_data.log_Cumulative_Excretion = log.(urine_data.Cumulative_Excretion_ug)

println("  Urine data points: $(nrow(urine_data))")
println("  Time range: $(minimum(urine_data.Time_days)) - $(maximum(urine_data.Time_days)) days")
println("  Cumulative excretion range: $(round(minimum(urine_data.Cumulative_Excretion_ug), digits=4)) - $(round(maximum(urine_data.Cumulative_Excretion_ug), digits=4)) μg")

# Exposure parameters
oral_dose_ug = 3.99  # Single oral dose (μg)
background_dose_ng_kg_day = 0.93 / 82.0  # Low background exposure (ng/kg-bw/day)
BW = 82.0  # Body weight (kg) - Haug et al. 2009

println("\nExposure parameters:")
println("  Single oral dose: $oral_dose_ug μg")
println("  Background exposure: $(round(background_dose_ng_kg_day, digits=4)) ng/kg-bw/day")
println("  Body weight: $BW kg")

# Function to run Worley PBK model for single oral dose
function run_pbk_single_dose(RAFapi, RAFbaso, Km_apical, PC_scale, observation_times_days,
                             oral_dose_ug, background_dose_ng_kg_day, BW=82.0)
    """
    Run Worley PBK model for single oral dose with background exposure

    Parameters:
    - RAFapi: Relative activity factor for apical transporters
    - RAFbaso: Relative activity factor for basolateral transporters
    - Km_apical: Michaelis constant for apical transporters (μg/L)
    - PC_scale: Partition coefficient scaling factor
    - observation_times_days: Times to extract predictions (days)
    - oral_dose_ug: Single oral dose (μg)
    - background_dose_ng_kg_day: Background exposure (ng/kg-bw/day)
    - BW: Body weight (kg)

    Returns:
    - CA values at observation times (μg/L)
    - Aurine values (cumulative urine excretion, u[10]) at observation times (μg)
    - Time values (days)
    """

    background_ug_per_day = background_dose_ng_kg_day * BW / 1000.0

    user_input = (
        BW = BW,
        substance = "PFHxA",
        ingestion = [background_ug_per_day],
        ingestion_time = [1.0e-10],
        admin_dose = oral_dose_ug,
        admin_time = 0.0,
        admin_type = "oral",
        exp_type = "single",
        time_scale = "days"
    )

    base_params = create_params_worley(user_input, RAF_api=RAFapi, RAF_baso=RAFbaso)
    params = (base_params..., PL = PL, PK = PK, PR = PR * PC_scale, Km_apical = Km_apical)

    inits = create_inits_worley(params)
    inits[11] = oral_dose_ug  # AST (stomach compartment)

    events = create_events_worley(params)
    tspan = (0.0, maximum(observation_times_days))
    saveat_times = sort(unique(observation_times_days))

    prob = ODEProblem(ode_func_worley, inits, tspan, params)#, callback=events)
    sol = solve(prob, Rodas5P(), saveat=saveat_times, reltol=1e-8, abstol=1e-10)

    # Extract CA concentrations (plasma) and Aurine (u[10])
    CA_values = []
    Aurine_values = []
    for u in sol.u
        Aplas_free = u[18]
        CA_free = Aplas_free / params.VPlas
        CA = CA_free / params.Free
        push!(CA_values, CA)
        push!(Aurine_values, u[10])  # Cumulative urine excretion
    end

    return CA_values, Aurine_values, sol.t
end

# Prepare data for MCMC
println("\n" * "="^80)
println("PREPARING DATA FOR MCMC")
println("="^80)

# Combine observation times from both plasma and urine
all_observation_times = sort(unique(vcat(plasma_data.Time_days, urine_data.Time_days)))
println("  Total unique observation times: $(length(all_observation_times))")
println("  Plasma observations: $(nrow(plasma_data))")
println("  Urine observations: $(nrow(urine_data))")

# Prior parameters
RAFapi_prior_mean = 0.0007  
RAFapi_prior_cv = 1.0  # 100% CV (wide prior)

RAFbaso_prior_mean = 0.5  
RAFbaso_prior_cv = 0.5  # 100% CV (wide prior)

Km_apical_prior_mean = 52338.44  # Default from model (μg/L)
Km_apical_prior_cv = 0.2  # 20% CV

PC_scale_prior_mean = 0.5  # Partition coefficient scaling factor
PC_scale_prior_std = 0.5  # Standard deviation

# Calculate LogNormal parameters for priors
function lognormal_params(mu, cv)
    """Convert mean and CV to LogNormal distribution parameters"""
    sigma = mu * cv
    mu_log = log(mu) - 0.5 * log(1 + (sigma/mu)^2)
    sigma_log = sqrt(log(1 + (sigma/mu)^2))
    return mu_log, sigma_log
end

mu_log_RAFapi, sigma_log_RAFapi = lognormal_params(RAFapi_prior_mean, RAFapi_prior_cv)
mu_log_RAFbaso, sigma_log_RAFbaso = lognormal_params(RAFbaso_prior_mean, RAFbaso_prior_cv)
mu_log_Km, sigma_log_Km = lognormal_params(Km_apical_prior_mean, Km_apical_prior_cv)

# Bayesian MCMC Model
@model function single_dose_with_urine_model(plasma_data, urine_data, all_observation_times,
                                              oral_dose_ug, background_dose_ng_kg_day, BW)
    """
    Bayesian model for single-dose PBK calibration with plasma and urine data

    Parameters to estimate:
    - RAFapi, RAFbaso, Km_apical, PC_scale
    - sigma_plasma: Plasma observation error (log scale)
    - sigma_urine: Urine observation error (log scale)
    """

    # Priors
    RAFapi ~ LogNormal(mu_log_RAFapi, sigma_log_RAFapi)
    RAFbaso ~ LogNormal(mu_log_RAFbaso, sigma_log_RAFbaso)
    Km_apical ~ LogNormal(mu_log_Km, sigma_log_Km)
    PC_scale ~ truncated(Normal(PC_scale_prior_mean, PC_scale_prior_std), 0.0, Inf)
    sigma_plasma ~ truncated(Normal(0.0, 0.2), 0.0, Inf)
    sigma_urine ~ truncated(Normal(0.0, 0.2), 0.0, Inf)

    try
        # Run PBK model
        CA_pred_all, Aurine_pred_all, time_solution = run_pbk_single_dose(
            RAFapi, RAFbaso, Km_apical, PC_scale, all_observation_times,
            oral_dose_ug, background_dose_ng_kg_day, BW
        )

        # Match predictions to plasma observation times
        for i in 1:nrow(plasma_data)
            t_obs = plasma_data.Time_days[i]
            idx = findfirst(==(t_obs), time_solution)
            if !isnothing(idx)
                log_CA_pred = log(CA_pred_all[idx])
                plasma_data.log_C_plasma[i] ~ Normal(log_CA_pred, sigma_plasma)
            end
        end

        # Match predictions to urine observation times
        for i in 1:nrow(urine_data)
            t_obs = urine_data.Time_days[i]
            idx = findfirst(==(t_obs), time_solution)
            if !isnothing(idx)
                log_Aurine_pred = log(Aurine_pred_all[idx])
                urine_data.log_Cumulative_Excretion[i] ~ Normal(log_Aurine_pred, sigma_urine)
            end
        end

    catch e
        println("Error in ODE solving: $e")
        Turing.@addlogprob! -1000.0
    end
end

println("\n" * "="^80)
println("BAYESIAN MODEL DEFINED")
println("="^80)

# MCMC Settings
n_chains = 4
n_iterations = 1000
n_burnin = 500

println("\nMCMC Settings:")
println("  Sampler: NUTS")
println("  Chains: $n_chains")
println("  Iterations per chain: $n_iterations")
println("  Burn-in: $n_burnin")

println("\n" * "="^80)
println("STARTING MCMC SAMPLING")
println("="^80)

# Run MCMC sampling
Random.seed!(123)
chains = sample(
    single_dose_with_urine_model(plasma_data, urine_data, all_observation_times,
                                 oral_dose_ug, background_dose_ng_kg_day, BW),
    NUTS(n_burnin, 0.65),
    MCMCThreads(),
    n_iterations,
    n_chains,
    discard_initial=n_burnin
)

println("\n" * "="^80)
println("MCMC SAMPLING COMPLETED!")
println("="^80)

# Save results
println("\nSaving MCMC results...")
using Serialization
output_dir = "PFAS_models/PFHxA"
mkpath(output_dir)  # Create directory if it doesn't exist

serialize("$output_dir/PFHxA_with_urine_chains.jls", chains)
println("  Saved chains to: $output_dir/PFHxA_with_urine_chains.jls")

describe(chains)

# Calculate additional diagnostics
println("\nParameter estimates (posterior means ± SD):")
println("  RAFapi = $(round(mean(chains[:RAFapi]), sigdigits=4)) ± $(round(std(chains[:RAFapi]), sigdigits=4))")
println("  RAFbaso = $(round(mean(chains[:RAFbaso]), sigdigits=4)) ± $(round(std(chains[:RAFbaso]), sigdigits=4))")
println("  Km_apical = $(round(mean(chains[:Km_apical]), sigdigits=4)) ± $(round(std(chains[:Km_apical]), sigdigits=4)) μg/L")
println("  PC_scale = $(round(mean(chains[:PC_scale]), sigdigits=4)) ± $(round(std(chains[:PC_scale]), sigdigits=4))")
println("  sigma_plasma = $(round(mean(chains[:sigma_plasma]), sigdigits=4)) ± $(round(std(chains[:sigma_plasma]), sigdigits=4))")
println("  sigma_urine = $(round(mean(chains[:sigma_urine]), sigdigits=4)) ± $(round(std(chains[:sigma_urine]), sigdigits=4))")

# Generate diagnostic plots
println("\nGenerating diagnostic plots...")
chains_plot = plot(chains, size=(2000, 1200), dpi=300)
savefig(chains_plot, "$output_dir/PFHxA_with_urine_diagnostic_plots.png")
println("  Saved diagnostic plots to: $output_dir/PFHxA_with_urine_diagnostic_plots.png")

println("\n" * "="^80)
println("CALIBRATION COMPLETE")
println("="^80)
