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
println("POSTERIOR PREDICTIVE CHECK - HIERARCHICAL POPULATION MODEL PFBS")
println("="^80)

# Create output directory
output_dir = "Population_model/PFBS/PPC_hierarchical"
if !isdir(output_dir)
    mkpath(output_dir)
    println("\nCreated output directory: $output_dir")
end

# Load MCMC chains from hierarchical model
println("\nLoading MCMC chains...")
chains = deserialize("Population_model/PFBS/PFBS_population_chains.jls")

# Load partition coefficients
println("Loading partition coefficients...")
partition_data = CSV.read("Population_model/partition_coefficients_results.csv", DataFrame)
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
    data = CSV.read("Population_model/PFBS/digitized_data/subject_$i.csv", DataFrame)

    # Round time to remove decimals
    data.Time_days = round.(data.Time_days)

    # Rename concentration column for consistency
    rename!(data, Symbol("C_serum_ug/L") => :C_serum)

    push!(subject_data, data)
    println("  Subject $i: $(nrow(data)) observations, time range: $(data.Time_days[1])-$(data.Time_days[end]) days")
end

# Background exposure parameter
background_dose_ng_kg_day = 0.93/70.0  # Low background exposure

println("\nBackground exposure: $background_dose_ng_kg_day ng/kg-bw/day")

# Extract posterior samples
println("\nExtracting posterior samples...")
chain_array = Array(chains)
n_samples = size(chain_array, 1)

# Extract parameter samples from chains
param_names = String.(names(chains))
println("  Available parameters: $(join(param_names, ", "))")

# Extract subject-specific RAFapi samples (one for each subject)
RAFapi_samples = Vector{Vector{Float64}}(undef, n_subjects)
for subject_idx in 1:n_subjects
    param_name = "RAFapi[$subject_idx]"
    RAFapi_idx = findfirst(==(param_name), param_names)
    if isnothing(RAFapi_idx)
        error("Parameter '$param_name' not found in chains. Available parameters: $(join(param_names, ", "))")
    end
    RAFapi_samples[subject_idx] = chain_array[:, RAFapi_idx]
    println("  Extracted $param_name: $(length(RAFapi_samples[subject_idx])) samples")
end

# Extract Km_apical samples (population-level)
Km_apical_idx = findfirst(==("Km_apical"), param_names)
if isnothing(Km_apical_idx)
    error("Parameter 'Km_apical' not found in chains. Available parameters: $(join(param_names, ", "))")
end
Km_apical_samples = chain_array[:, Km_apical_idx]

# Extract sigma samples (population-level)
sigma_idx = findfirst(==("sigma"), param_names)
if isnothing(sigma_idx)
    error("Parameter 'sigma' not found in chains. Available parameters: $(join(param_names, ", "))")
end
sigma_samples = chain_array[:, sigma_idx]

println("  Using all $n_samples posterior samples")

# Function to run Worley PBK model for elimination phase
function run_pbk_elimination(RAFapi, Km_apical, initial_C_serum, observation_times_days,
                             background_dose_ng_kg_day, BW=70.0)
    """
    Run Worley PBK model for elimination phase with initial conditions

    Parameters:
    - RAFapi: Relative activity factor for apical transporters
    - Km_apical: Michaelis constant for apical transporters
    - initial_C_serum: Initial serum concentration (ug/L)
    - observation_times_days: Times to extract predictions (days)
    - background_dose_ng_kg_day: Background exposure (ng/kg-bw/day)
    - BW: Body weight (kg)

    Returns:
    - CA values at observation times (ug/L)
    """

    # Convert background dose from ng/kg-bw/day to ug/day
    background_ug_per_day = background_dose_ng_kg_day * BW * (1.0/1000.0)

    # Set up user input for Worley model
    user_input = (
        BW = BW,
        substance = "PFBS",
        ingestion = [background_ug_per_day],
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

    # Create initial conditions vector
    inits = create_inits_worley(params)

    # Update relevant compartments with calculated initial conditions
    inits[16] = AL_initial   # Liver free
    inits[8] = AK_initial   # PTC (Proximal Tubule Cells)
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
    CA_values = Float64[]
    for u in sol.u
        Aplas_free = u[18]
        CA_free = Aplas_free / VPlas
        CA = CA_free / Free
        push!(CA_values, CA)
    end

    return CA_values
end

println("\n" * "="^80)
println("GENERATING POSTERIOR PREDICTIVE SAMPLES")
println("="^80)

# Store predictions for all subjects
all_predictions = []

for subject_idx in 1:n_subjects
    data = subject_data[subject_idx]
    initial_C_serum = data.C_serum[1]
    max_time = data.Time_days[end]

    println("\nProcessing Subject $subject_idx...")
    println("  Initial concentration: $(round(initial_C_serum, digits=2)) ug/L")
    println("  Time range: 0-$(max_time) days")

    # Create dense grid of time points (every 0.1 day)
    time_points = collect(0.0:0.1:max_time)
    println("  Prediction time points: $(length(time_points))")

    # Generate predictions for each posterior sample
    predictions_matrix = zeros(n_samples, length(time_points))

    for i in 1:n_samples
        try
            # Run model with subject-specific RAFapi and population-level Km_apical
            CA_pred = run_pbk_elimination(
                RAFapi_samples[subject_idx][i],  # Use subject-specific RAFapi
                Km_apical_samples[i],
                initial_C_serum,
                time_points,
                background_dose_ng_kg_day
            )

            # Add observation error (log-normal)
            # Sample error once per posterior sample (paired with parameters)
            error = rand(Normal(0, sigma_samples[i]))
            for t in 1:length(time_points)
                predictions_matrix[i, t] = CA_pred[t] * exp(error)
            end

        catch e
            println("  Warning: ODE solve failed for sample $i: $e")
            predictions_matrix[i, :] .= NaN
        end
    end

    # Remove failed samples
    valid_samples = .!isnan.(predictions_matrix[:, 1])
    predictions_matrix = predictions_matrix[valid_samples, :]
    n_valid = sum(valid_samples)

    println("  Valid samples: $n_valid / $n_samples")

    if n_valid < 100
        println("  Warning: Only $n_valid valid samples. Skipping subject $subject_idx")
        continue
    end

    # Calculate prediction intervals
    lower_bound = Float64[quantile(predictions_matrix[:, t], 0.025) for t in 1:length(time_points)]
    upper_bound = Float64[quantile(predictions_matrix[:, t], 0.975) for t in 1:length(time_points)]
    median_pred = Float64[median(predictions_matrix[:, t]) for t in 1:length(time_points)]

    # Store for multi-panel plot
    push!(all_predictions, (
        subject_idx = subject_idx,
        time_points = time_points,
        median_pred = median_pred,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        obs_time = data.Time_days,
        obs_conc = data.C_serum
    ))
end

# Create multi-panel plot (3x2 grid)
println("\n" * "="^80)
println("CREATING MULTI-PANEL PLOT")
println("="^80)

p = plot(layout=(3, 2), size=(1400, 1200), legend=:topright)

for (idx, pred_data) in enumerate(all_predictions)
    subject_idx = pred_data.subject_idx

    # Plot prediction interval
    plot!(p[idx], pred_data.time_points, pred_data.median_pred,
          ribbon=(pred_data.median_pred .- pred_data.lower_bound,
                  pred_data.upper_bound .- pred_data.median_pred),
          fillalpha=0.3,
          label="95% Prediction Interval",
          xlabel="Time (days)",
          ylabel="Serum Concentration (ug/L)",
          title="Subject $subject_idx",
          linewidth=2,
          color=:blue)

    # Add observed data points
    scatter!(p[idx], pred_data.obs_time, pred_data.obs_conc,
             label="Observed Data",
             markersize=6,
             markercolor=:red,
             markerstrokewidth=2)
end

# Save multi-panel plot
filename = joinpath(output_dir, "All_Subjects_Hierarchical_PPC.png")
savefig(p, filename)
println("\nSaved multi-panel plot: $filename")

println("\n" * "="^80)
println("POSTERIOR PREDICTIVE CHECK COMPLETE")
println("="^80)
println("Plot saved to: $output_dir")
println("="^80)
