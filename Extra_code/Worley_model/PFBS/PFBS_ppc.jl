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
using Distributions

# Include Worley PBK model functions
include("../PBK_model/julia_code.jl")

println("="^80)
println("POSTERIOR PREDICTIVE CHECK - WORLEY MODEL")
println("="^80)

# Create output directory
output_dir = "Worley_model/PFBS/PPC"
if !isdir(output_dir)
    mkpath(output_dir)
    println("\nCreated output directory: $output_dir")
end

# Load MCMC chains
println("\nLoading MCMC chains...")
chains = deserialize("Worley_model/PFBS/Worley_PFBS_calibration_chains.jls")

# Load intake estimation data
println("Loading experimental data...")
intake_data = CSV.read("Worley_model/PFBS/Data/intake_estimation_results.csv", DataFrame)

# Load partition coefficients
partition_data = CSV.read("Worley_model/Worley_partition_coefficients_results.csv", DataFrame)
pfbs_partitions = filter(row -> row.Compound == "PFBS", partition_data)
PL_worley = pfbs_partitions[pfbs_partitions.Tissue .== "liver", :Partition_Coefficient][1]
PK_worley = pfbs_partitions[pfbs_partitions.Tissue .== "kidney", :Partition_Coefficient][1]
PR_worley = pfbs_partitions[pfbs_partitions.Tissue .== "rest of body", :Partition_Coefficient][1]

println("\nData loaded:")
println("  Total subjects: $(nrow(intake_data))")
println("  Studies: $(join(unique(intake_data.Study), ", "))")

# Extract posterior samples
println("\nExtracting posterior samples...")
chain_array = Array(chains)
n_samples = size(chain_array, 1)

# Extract parameter samples from chains
param_names = String.(names(chains))
println("  Available parameters: $(join(param_names, ", "))")

# Extract RAFapi samples
RAFapi_idx = findfirst(==("RAFapi"), param_names)
if isnothing(RAFapi_idx)
    error("Parameter 'RAFapi' not found in chains. Available parameters: $(join(param_names, ", "))")
end
RAFapi_samples = chain_array[:, RAFapi_idx]

# Extract sigma samples
sigma_idx = findfirst(==("sigma"), param_names)
if isnothing(sigma_idx)
    error("Parameter 'sigma' not found in chains. Available parameters: $(join(param_names, ", "))")
end
sigma_samples = chain_array[:, sigma_idx]

println("  Using all $n_samples posterior samples for RAFapi and sigma")

# Function to run PBK model and get predictions
function run_pbk_prediction(RAFapi, exposure_dose_ng_kg_day, time_points_years, BW=70.0)
    """
    Run Worley PBK model and get predictions at specified time points

    Parameters:
    - RAFapi: Relative activity factor for apical transporters
    - exposure_dose_ng_kg_day: Daily intake (ng/kg-bw/day)
    - time_points_years: Time points to extract predictions (years)
    - BW: Body weight (kg)

    Returns:
    - CA values at time points (μg/L)
    """

    # Convert intake from ng/kg-bw/day to μg/year for the model
    intake_ug_per_year = exposure_dose_ng_kg_day * BW * (1.0/1000.0) * 365.0

    # Set up user input for Worley model
    user_input = (
        BW = BW,
        substance = "PFBS",
        ingestion = [intake_ug_per_year],
        ingestion_time = [0.000000001],  # Small delay needed for callbacks
        admin_dose = 0.0,
        admin_time = 0.0,
        admin_type = "oral",
        exp_type = "continuous",
        time_scale = "years"
    )

    # Create parameters with RAFapi (RAFbaso uses default value from create_params_worley)
    base_params = create_params_worley(user_input, RAF_api=RAFapi)

    # Update with partition coefficients
    params = (base_params..., PL = PL_worley, PK = PK_worley, PR = PR_worley)

    # Create initial conditions
    inits = create_inits_worley(params)

    # Create events
    events = create_events_worley(params)

    # Time span and save times
    tspan = (0.0, maximum(time_points_years))
    saveat_times = sort(unique(vcat([0.0], time_points_years)))

    # Solve ODE
    prob = ODEProblem(ode_func_worley, inits, tspan, params, callback=events)
    sol = solve(prob, Rodas5P(), saveat=saveat_times, reltol=1e-6, abstol=1e-8)

    # Extract CA concentrations
    CA_values = Float64[]
    for u in sol.u
        Aplas_free = u[18]
        CA_free = Aplas_free / params.VPlas
        CA = CA_free / params.Free
        push!(CA_values, CA)
    end

    return CA_values
end

# Time points for predictions
# Add early time points to capture initial dynamics after ingestion starts
# Then regular 0.1 year intervals
early_times = [0.0, 0.0001, 0.001, 0.01, 0.05]
regular_times = collect(0.1:0.1:10.0)
time_points = sort(unique(vcat(early_times, regular_times)))

println("\n" * "="^80)
println("GENERATING POSTERIOR PREDICTIVE SAMPLES")
println("="^80)
println("Time points: $(length(time_points)) points from $(time_points[1]) to $(time_points[end]) years")

# Process each subject
for subject_idx in 1:nrow(intake_data)
    subject = intake_data[subject_idx, :]
    study = subject.Study
    exposure_dose = subject.Estimated_Intake_ng_kg_bw_day

    println("\nProcessing Subject $subject_idx (Study: $study)...")

    # Generate predictions for each posterior sample
    predictions_matrix = zeros(n_samples, length(time_points))

    for i in 1:n_samples
        try
            # Run model with posterior sample
            CA_pred = run_pbk_prediction(
                RAFapi_samples[i],
                exposure_dose,
                time_points
            )

            # Add observation error (log-normal)
            # log(CA) ~ Normal(log(CA_pred), sigma)
            # Equivalent to: CA = CA_pred * exp(error) where error ~ Normal(0, sigma)
            # Sample error once per posterior sample (paired with RAFapi)
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

    if n_valid < 100
        println("  Warning: Only $n_valid valid samples. Skipping subject $subject_idx")
        continue
    end

    # Calculate prediction intervals (2.5% and 97.5% percentiles)
    lower_bound = Float64[quantile(predictions_matrix[:, t], 0.025) for t in 1:length(time_points)]
    upper_bound = Float64[quantile(predictions_matrix[:, t], 0.975) for t in 1:length(time_points)]
    median_pred = Float64[median(predictions_matrix[:, t]) for t in 1:length(time_points)]

    # Create plot
    p = plot(time_points, median_pred,
             ribbon=(median_pred .- lower_bound, upper_bound .- median_pred),
             fillalpha=0.3,
             label="95% Prediction Interval",
             xlabel="Time (years)",
             ylabel="Serum Concentration (μg/L)",
             title="Subject $subject_idx - Study: $study\nIntake: $(round(exposure_dose, digits=1)) ng/kg-bw/day",
             linewidth=2,
             legend=:bottomright,
             size=(800, 600))

    # Add observed concentration at 10 years
    scatter!(p, [10.0], [100.0],
             label="Observed (10 years)",
             markersize=8,
             markercolor=:red,
             markerstrokewidth=2)

    # Save plot
    filename = joinpath(output_dir, "Subject_$(subject_idx)_$(study).png")
    savefig(p, filename)
    println("  Saved: $filename")
end

println("\n" * "="^80)
println("POSTERIOR PREDICTIVE CHECK COMPLETE")
println("="^80)
println("Plots saved to: $output_dir")
println("="^80)
