using Serialization
using DifferentialEquations
using DataFrames
using CSV
using MCMCChains
using Turing
using StatsPlots
using Plots
using Statistics
using Random
using Sundials

# Include Worley PBK model functions
include("../PBK_model/julia_code.jl")

println("="^80)
println("POSTERIOR PREDICTIVE CHECK - SINGLE-DOSE PFHxA MODEL (PLASMA + URINE)")
println("="^80)

# Create output directory
output_dir = "Population_model/PFHxA/PPC"
if !isdir(output_dir)
    mkpath(output_dir)
    println("\nCreated output directory: $output_dir")
end

# Load MCMC chains
println("\nLoading MCMC chains...")
chains = deserialize("Population_model/PFHxA/PFHxA_with_urine_chains.jls")

# Load partition coefficients
println("Loading partition coefficients...")
partition_data = CSV.read("Population_model/partition_coefficients_results.csv", DataFrame)
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
plasma_data = CSV.read("Population_model/PFHxA/digitized_data/Abraham_2024_plasma_data.csv", DataFrame)
rename!(plasma_data, Symbol("C_plasma_ug/L") => :C_plasma)
plasma_data = dropmissing(plasma_data)

println("  Plasma data points: $(nrow(plasma_data))")
println("  Time range: $(minimum(plasma_data.Time_days)) - $(maximum(plasma_data.Time_days)) days")

# Load urine data
println("\nLoading urine data...")
urine_data = CSV.read("Population_model/PFHxA/digitized_data/Abraham_2024_urine_data_processed.csv", DataFrame)
urine_data = dropmissing(urine_data)

println("  Urine data points: $(nrow(urine_data))")
println("  Time range: $(minimum(urine_data.Time_days)) - $(maximum(urine_data.Time_days)) days")

# Exposure parameters
oral_dose_ug = 3.99
background_dose_ng_kg_day = 0.93 / 82.0
BW = 82.0

println("\nExposure parameters:")
println("  Single oral dose: $oral_dose_ug μg")
println("  Background exposure: $(round(background_dose_ng_kg_day, digits=4)) ng/kg-bw/day")
println("  Body weight: $BW kg")

# Function to run Worley PBK model for single oral dose
function run_pbk_single_dose(RAFapi, RAFbaso, Km_apical, PC_scale, observation_times_days,
                             oral_dose_ug, background_dose_ng_kg_day, BW=82.0)
    """
    Run Worley PBK model for single oral dose with background exposure

    Returns:
    - CA values at observation times (μg/L)
    - Aurine values (cumulative urine excretion, u[10]) at observation times (μg)
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

    prob = ODEProblem(ode_func_worley, inits, tspan, params)
    sol = solve(prob, Rodas5P(), saveat=saveat_times, reltol=1e-8, abstol=1e-10)

    # Extract CA concentrations (plasma) and Aurine (u[10])
    CA_values = Float64[]
    Aurine_values = Float64[]
    for u in sol.u
        Aplas_free = u[18]
        CA_free = Aplas_free / params.VPlas
        CA = CA_free / params.Free
        push!(CA_values, CA)
        push!(Aurine_values, u[10])  # Cumulative urine excretion
    end

    return CA_values, Aurine_values
end

# Extract posterior samples
println("\n" * "="^80)
println("EXTRACTING POSTERIOR SAMPLES")
println("="^80)

chain_array = Array(chains)
n_samples = size(chain_array, 1)

param_names = String.(names(chains))
println("  Available parameters: $(join(param_names, ", "))")

# Extract parameter samples
RAFapi_idx = findfirst(==(("RAFapi")), param_names)
RAFbaso_idx = findfirst(==(("RAFbaso")), param_names)
Km_apical_idx = findfirst(==(("Km_apical")), param_names)
PC_scale_idx = findfirst(==(("PC_scale")), param_names)
sigma_plasma_idx = findfirst(==(("sigma_plasma")), param_names)
sigma_urine_idx = findfirst(==(("sigma_urine")), param_names)

RAFapi_samples = chain_array[:, RAFapi_idx]
RAFbaso_samples = chain_array[:, RAFbaso_idx]
Km_apical_samples = chain_array[:, Km_apical_idx]
PC_scale_samples = chain_array[:, PC_scale_idx]
sigma_plasma_samples = chain_array[:, sigma_plasma_idx]
sigma_urine_samples = chain_array[:, sigma_urine_idx]

println("  Using all $n_samples posterior samples")

# Generate predictions
println("\n" * "="^80)
println("GENERATING POSTERIOR PREDICTIVE SAMPLES")
println("="^80)

# Create dense time grid for smooth predictions
max_time = max(maximum(plasma_data.Time_days), maximum(urine_data.Time_days))
time_points = collect(range(0.0, max_time, length=300))

println("  Prediction time points: $(length(time_points))")
println("  Time range: 0-$(max_time) days")

# Initialize prediction matrices
predictions_plasma = zeros(n_samples, length(time_points))
predictions_urine = zeros(n_samples, length(time_points))

for i in 1:n_samples
    try
        # Run model
        CA_pred, Aurine_pred = run_pbk_single_dose(
            RAFapi_samples[i],
            RAFbaso_samples[i],
            Km_apical_samples[i],
            PC_scale_samples[i],
            time_points,
            oral_dose_ug,
            background_dose_ng_kg_day,
            BW
        )

        # Add observation error (log-normal)
        error_plasma = rand(Normal(0, sigma_plasma_samples[i]))
        error_urine = rand(Normal(0, sigma_urine_samples[i]))

        for t in 1:length(time_points)
            predictions_plasma[i, t] = CA_pred[t] * exp(error_plasma)
            predictions_urine[i, t] = Aurine_pred[t] * exp(error_urine)
        end

    catch e
        println("  Warning: ODE solve failed for sample $i: $e")
        predictions_plasma[i, :] .= NaN
        predictions_urine[i, :] .= NaN
    end
end

# Remove failed samples
valid_samples = .!isnan.(predictions_plasma[:, 1])
predictions_plasma = predictions_plasma[valid_samples, :]
predictions_urine = predictions_urine[valid_samples, :]
n_valid = sum(valid_samples)

println("  Valid samples: $n_valid / $n_samples")

if n_valid < 100
    error("Too few valid samples ($n_valid). Cannot generate reliable predictions.")
end

# Calculate prediction intervals
println("\n" * "="^80)
println("CALCULATING PREDICTION INTERVALS")
println("="^80)

lower_plasma = Float64[quantile(predictions_plasma[:, t], 0.025) for t in 1:length(time_points)]
upper_plasma = Float64[quantile(predictions_plasma[:, t], 0.975) for t in 1:length(time_points)]
median_plasma = Float64[median(predictions_plasma[:, t]) for t in 1:length(time_points)]

lower_urine = Float64[quantile(predictions_urine[:, t], 0.025) for t in 1:length(time_points)]
upper_urine = Float64[quantile(predictions_urine[:, t], 0.975) for t in 1:length(time_points)]
median_urine = Float64[median(predictions_urine[:, t]) for t in 1:length(time_points)]

# Plot 1: Plasma only
println("\n" * "="^80)
println("CREATING PLASMA PLOT")
println("="^80)

p1 = plot(time_points, median_plasma,
          ribbon=(median_plasma .- lower_plasma, upper_plasma .- median_plasma),
          fillalpha=0.3,
          label="95% Prediction Interval",
          xlabel="Time (days)",
          ylabel="Plasma Concentration (μg/L)",
          title="PFHxA Single-Dose Model - Plasma PPC",
          linewidth=2,
          color=:blue,
          size=(900, 600),
          legend=:topright,
          left_margin=6Plots.mm,
          bottom_margin=6Plots.mm)

scatter!(p1, plasma_data.Time_days, plasma_data.C_plasma,
         label="Observed Data",
         markersize=6,
         markercolor=:red,
         markerstrokewidth=2)

filename_plasma = joinpath(output_dir, "PFHxA_Plasma_PPC_with_urine.png")
savefig(p1, filename_plasma)
println("  Saved plasma plot: $filename_plasma")

# Plot 2: Urine only
println("\n" * "="^80)
println("CREATING URINE PLOT")
println("="^80)

p2 = plot(time_points, median_urine,
          ribbon=(median_urine .- lower_urine, upper_urine .- median_urine),
          fillalpha=0.3,
          label="95% Prediction Interval",
          xlabel="Time (days)",
          ylabel="Cumulative Urine Excretion (μg)",
          title="PFHxA Single-Dose Model - Urine PPC",
          linewidth=2,
          color=:green,
          size=(900, 600),
          legend=:bottomright,
          left_margin=6Plots.mm,
          bottom_margin=6Plots.mm)

scatter!(p2, urine_data.Time_days, urine_data.Cumulative_Excretion_ug,
         label="Observed Data",
         markersize=6,
         markercolor=:red,
         markerstrokewidth=2)

filename_urine = joinpath(output_dir, "PFHxA_Urine_PPC_with_urine.png")
savefig(p2, filename_urine)
println("  Saved urine plot: $filename_urine")

# Plot 3: Combined 2-panel plot
println("\n" * "="^80)
println("CREATING COMBINED PLOT")
println("="^80)

p_combined = plot(layout=(2, 1), size=(1000, 1000),
                  left_margin=8Plots.mm,
                  bottom_margin=6Plots.mm,
                  top_margin=4Plots.mm,
                  right_margin=4Plots.mm)

# Plasma panel
plot!(p_combined[1], time_points, median_plasma,
      ribbon=(median_plasma .- lower_plasma, upper_plasma .- median_plasma),
      fillalpha=0.3,
      label="95% Prediction Interval",
      xlabel="Time (days)",
      ylabel="Plasma Concentration (μg/L)",
      title="PFHxA Single-Dose Model - Plasma",
      linewidth=2,
      color=:blue,
      legend=:topright)

scatter!(p_combined[1], plasma_data.Time_days, plasma_data.C_plasma,
         label="Observed Data",
         markersize=6,
         markercolor=:red,
         markerstrokewidth=2)

# Urine panel
plot!(p_combined[2], time_points, median_urine,
      ribbon=(median_urine .- lower_urine, upper_urine .- median_urine),
      fillalpha=0.3,
      label="95% Prediction Interval",
      xlabel="Time (days)",
      ylabel="Cumulative Urine Excretion (μg)",
      title="PFHxA Single-Dose Model - Urine",
      linewidth=2,
      color=:green,
      legend=:bottomright)

scatter!(p_combined[2], urine_data.Time_days, urine_data.Cumulative_Excretion_ug,
         label="Observed Data",
         markersize=6,
         markercolor=:red,
         markerstrokewidth=2)

filename_combined = joinpath(output_dir, "PFHxA_Combined_PPC_with_urine.png")
savefig(p_combined, filename_combined)
println("  Saved combined plot: $filename_combined")

println("\n" * "="^80)
println("POSTERIOR PREDICTIVE CHECK COMPLETE")
println("="^80)
println("Plots saved to: $output_dir")
println("="^80)
