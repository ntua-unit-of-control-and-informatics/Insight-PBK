using DifferentialEquations
using DataFrames
using CSV
using Plots
using Sundials

# Include Worley PBK model functions
include("../PBK_model/julia_code.jl")

"""
Check ODE solutions for PFBS elimination phase
Solve ODEs once for each subject with a single RAFapi value
"""

println("="^80)
println("ODE CHECK - PFBS ELIMINATION PHASE")
println("="^80)

# Load partition coefficients
println("\nLoading partition coefficients...")
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

# Define single RAFapi value for all subjects
RAFapi = 0.007  # Prior mean value

println("\nRAFapi value: $RAFapi")

# Function to run Worley PBK model for elimination phase
function run_pbk_elimination(RAFapi, initial_C_serum, observation_times_days,
                             background_dose_ng_kg_day, BW=70.0)
    """
    Run Worley PBK model for elimination phase with initial conditions

    Parameters:
    - RAFapi: Relative activity factor for apical transporters
    - initial_C_serum: Initial serum concentration (ug/L)
    - observation_times_days: Times to extract predictions (days)
    - background_dose_ng_kg_day: Background exposure (ng/kg-bw/day)
    - BW: Body weight (kg)

    Returns:
    - CA values at observation times (ug/L)
    - Time values (days)
    """

    # Convert background dose from ng/kg-bw/day to ug/day
    background_ug_per_day = background_dose_ng_kg_day * BW * (1.0/1000.0)

    # Set up user input for Worley model
    user_input = (
        BW = BW,
        substance = "PFBS",
        ingestion = [background_ug_per_day],  # Low background exposure
        ingestion_time = [0.000000001],  # Small delay needed for callbacks
        admin_dose = 0.0,
        admin_time = 0.0,
        admin_type = "oral",
        exp_type = "continuous",
        time_scale = "days"
    )

    # Create base parameters with RAFapi
    base_params = create_params_worley(user_input, RAF_api=RAFapi)

    # Update with partition coefficients
    params = (base_params..., PL = PL, PK = PK, PR = PR, Km_apical = 10000000000.7905)#1500000) #77500.0)
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

    # Extract CA concentrations and urine amounts
    CA_values = []
    Aurine_values = []
    A8_values = []
    A9_values = []
    for u in sol.u
        Aplas_free = u[18]
        CA_free = Aplas_free / VPlas
        CA = CA_free / Free
        push!(CA_values, CA)

        # Extract urine compartment (index 10 in Worley model)
        Aurine = u[10]
        push!(Aurine_values, Aurine)

        # Extract compartment 8
        A8 = u[8]
        push!(A8_values, A8)

        # Extract compartment 9
        A9 = u[9]
        push!(A9_values, A9)
    end

    return CA_values, Aurine_values, A8_values, A9_values, sol.t
end

# Solve ODEs for each subject and prepare data for plotting
println("\n" * "="^80)
println("SOLVING ODES FOR ALL SUBJECTS")
println("="^80)

plot_data = []

for i in 1:n_subjects
    data = subject_data[i]
    observation_times = data.Time_days
    initial_C_serum = data.C_serum[1]

    println("\nSubject $i:")
    println("  Initial C_serum: $(round(initial_C_serum, digits=2)) ug/L")

    # Run PBK model
    CA_pred, Aurine_pred, A8_pred, A9_pred, time_pred = run_pbk_elimination(
        RAFapi,
        initial_C_serum,
        observation_times,
        background_dose_ng_kg_day
    )

    println("  ODE solved successfully")
    println("  Time points: $(length(time_pred))")
    println("  Final urine amount: $(round(Aurine_pred[end], digits=2)) ug")
    println("  Final A8 value: $(round(A8_pred[end], digits=2)) ug")
    println("  Final A9 value: $(round(A9_pred[end], digits=2)) ug")

    push!(plot_data, (subject=i, time=time_pred, CA=CA_pred, Aurine=Aurine_pred, A8=A8_pred, A9=A9_pred,
                      obs_time=observation_times, obs_CA=data.C_serum))
end

# Create plot
println("\n" * "="^80)
println("CREATING PLOTS")
println("="^80)

p = plot(layout=(3,2), size=(1200, 900),
         xlabel="Time (days)", ylabel="Serum Concentration (ug/L)",
         legend=:topright)

for i in 1:n_subjects
    subject_plot = plot_data[i]

    # Plot model predictions
    plot!(p[i], subject_plot.time, subject_plot.CA,
          label="Model prediction", linewidth=2, color=:blue)

    # Plot observations
    scatter!(p[i], subject_plot.obs_time, subject_plot.obs_CA,
             label="Observed data", markersize=5, color=:red)

    # Add title
    plot!(p[i], title="Subject $i")
end

# Save serum concentration plot
output_file = "Population_model/PFBS/check_odes_serum_plot.png"
savefig(p, output_file)
println("\nSerum concentration plot saved to: $output_file")

# Create urine excretion plot
p_urine = plot(size=(800, 600),
               xlabel="Time (days)",
               ylabel="Cumulative Urine Excretion (ug)",
               title="Cumulative PFBS Excretion in Urine",
               legend=:bottomright)

for i in 1:n_subjects
    subject_plot = plot_data[i]

    # Plot cumulative urine excretion for each subject
    plot!(p_urine, subject_plot.time, subject_plot.Aurine,
          label="Subject $i", linewidth=2)
end

# Save urine plot
output_file_urine = "Population_model/PFBS/check_odes_urine_plot.png"
savefig(p_urine, output_file_urine)
println("Urine excretion plot saved to: $output_file_urine")

# Create A8 compartment plot
p_a8 = plot(size=(800, 600),
            xlabel="Time (days)",
            ylabel="Compartment 8 Amount (ug)",
            title="PFBS Amount in Compartment 8",
            legend=:bottomright)

for i in 1:n_subjects
    subject_plot = plot_data[i]

    # Plot A8 for each subject
    plot!(p_a8, subject_plot.time, subject_plot.A8,
          label="Subject $i", linewidth=2)
end

# Save A8 plot
output_file_a8 = "Population_model/PFBS/check_odes_A8_plot.png"
savefig(p_a8, output_file_a8)
println("Compartment 8 plot saved to: $output_file_a8")

# Create A9 compartment plot
p_a9 = plot(size=(800, 600),
            xlabel="Time (days)",
            ylabel="Compartment 9 Amount (ug)",
            title="PFBS Amount in Compartment 9",
            legend=:bottomright)

for i in 1:n_subjects
    subject_plot = plot_data[i]

    # Plot A9 for each subject
    plot!(p_a9, subject_plot.time, subject_plot.A9,
          label="Subject $i", linewidth=2)
end

# Save A9 plot
output_file_a9 = "Population_model/PFBS/check_odes_A9_plot.png"
savefig(p_a9, output_file_a9)
println("Compartment 9 plot saved to: $output_file_a9")

println("\n" * "="^80)
println("ODE CHECK COMPLETE")
println("="^80)
