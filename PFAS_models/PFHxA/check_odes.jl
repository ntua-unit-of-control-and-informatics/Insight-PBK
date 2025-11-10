using DifferentialEquations
using DataFrames
using CSV
using Plots
using Sundials
using Statistics

# Include Worley PBK model functions
include("../PBK_model/julia_code.jl")

"""
Check ODE solutions for PFHxA single-dose
Solve ODEs with specified RAFapi and Km_apical values
"""

println("="^80)
println("ODE CHECK - PFHxA SINGLE-DOSE")
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
data = CSV.read("PFAS_models/PFHxA/digitized_data/Abraham_2024_plasma_data.csv", DataFrame)

# Rename concentration column for consistency
rename!(data, Symbol("C_plasma_ug/L") => :C_plasma)

# Remove any rows with missing values
data = dropmissing(data)

println("  Data points: $(nrow(data))")
println("  Time range: $(minimum(data.Time_days)) - $(maximum(data.Time_days)) days")
println("  Concentration range: $(round(minimum(data.C_plasma), digits=4)) - $(round(maximum(data.C_plasma), digits=4)) μg/L")

# Exposure parameters
oral_dose_ug = 3.99  # Single oral dose (μg)
background_dose_ng_kg_day = 0.93 / 82.0  # Low background exposure (ng/kg-bw/day)
BW = 82.0  # Body weight (kg)

println("\nExposure parameters:")
println("  Single oral dose: $oral_dose_ug μg")
println("  Background exposure: $(round(background_dose_ng_kg_day, digits=4)) ng/kg-bw/day")
println("  Body weight: $BW kg")

# ============================================================================
# MODIFY THESE PARAMETERS TO TEST DIFFERENT VALUES
# ============================================================================
PC_scale = 0.2472 # Partition coefficient scaling factor
RAFapi = 2.532e-5  # Relative activity factor for apical transporters
RAFbaso = 0.02332  # Relative activity factor for basolateral transporters
Km_apical = 401305.1965  # Michaelis constant for apical transporters (μg/L)
Km_baso = 15859.5  # Michaelis constant for basolateral transporters (μg/L)
# ============================================================================

println("\n" * "="^80)
println("MODEL PARAMETERS (MODIFY THESE TO TEST DIFFERENT VALUES)")
println("="^80)
println("  RAFapi = $RAFapi")
println("  Km_apical = $Km_apical μg/L")
println("="^80)

# Function to run Worley PBK model for single oral dose
function run_pbk_single_dose(RAFapi, Km_apical, observation_times_days,
                             oral_dose_ug, background_dose_ng_kg_day, BW=82.0)
    """
    Run Worley PBK model for single oral dose with background exposure

    Parameters:
    - RAFapi: Relative activity factor for apical transporters
    - Km_apical: Michaelis constant for apical transporters (μg/L)
    - observation_times_days: Times to extract predictions (days)
    - oral_dose_ug: Single oral dose (μg)
    - background_dose_ng_kg_day: Background exposure (ng/kg-bw/day)
    - BW: Body weight (kg)

    Returns:
    - CA values at observation times (μg/L)
    - Time values (days)
    - Full solution for additional compartment analysis
    """

    # Convert background dose from ng/kg-bw/day to μg/day
    background_ug_per_day = background_dose_ng_kg_day * BW / 1000.0

    # Set up user input for Worley model
    user_input = (
        BW = BW,
        substance = "PFHxA",
        ingestion = [background_ug_per_day],  # Background exposure
        ingestion_time = [1.0e-10],  # Small delay for callbacks
        admin_dose = oral_dose_ug,  # Single oral dose
        admin_time = 0.0,  # Dose at t=0
        admin_type = "oral",
        exp_type = "single",
        time_scale = "days"
    )

    # Create base parameters with RAFapi
    base_params = create_params_worley(user_input, RAF_api=RAFapi, RAF_baso=RAFbaso)

    # Update with partition coefficients and Km_apical
    params = (base_params..., PL = PL, PK = PK, PR = PR*PC_scale, Km_apical = Km_apical,
    Km_baso = Km_baso)

    # Create initial conditions (all zeros for single dose)
    inits = create_inits_worley(params)

    # For single oral dose, add the dose to stomach compartment
    inits[11] = oral_dose_ug  # AST (stomach compartment)

    # Create events
    events = create_events_worley(params)

    # Time span - create a fine grid for smooth plotting
    tspan = (0.0, maximum(observation_times_days))

    # Create a fine time grid for plotting
    fine_times = collect(range(0.0, maximum(observation_times_days), length=200))
    all_times = sort(unique(vcat(fine_times, observation_times_days)))

    # Solve ODE
    prob = ODEProblem(ode_func_worley, inits, tspan, params, callback=events)
    sol = solve(prob, Rodas5P(), saveat=all_times, reltol=1e-8, abstol=1e-10)

    # Extract CA concentrations (total plasma concentration)
    CA_values = []
    Aplas_free_values = []  # Free plasma amount (u[18])
    Aplas_free_over_Free_values = []  # u[18]/Free
    AR_values = []         # Rest of body (u[1])
    AKb_values = []        # Kidney blood (u[4])
    APTC_values = []       # Proximal tubule cells (u[8])
    Afil_values = []       # Filtrate (u[9])
    Aurine_values = []     # Urine compartment (u[10])
    AST_values = []        # Stomach compartment (u[11])
    ASI_values = []        # Small intestine compartment (u[13])
    Afeces_values = []     # Feces compartment (u[15])
    AL_values = []         # Liver compartment (u[16])

    for u in sol.u
        Aplas_free = u[18]
        CA_free = Aplas_free / params.VPlas
        CA = CA_free / params.Free  # Convert free to total
        push!(CA_values, CA)
        push!(Aplas_free_values, Aplas_free)
        push!(Aplas_free_over_Free_values, Aplas_free / params.Free)

        # Extract rest of body (index 1)
        AR = u[1]
        push!(AR_values, AR)

        # Extract kidney blood (index 4)
        AKb = u[4]
        push!(AKb_values, AKb)

        # Extract proximal tubule cells (index 8)
        APTC = u[8]
        push!(APTC_values, APTC)

        # Extract filtrate (index 9)
        Afil = u[9]
        push!(Afil_values, Afil)

        # Extract urine compartment (index 10)
        Aurine = u[10]
        push!(Aurine_values, Aurine)

        # Extract stomach (index 11)
        AST = u[11]
        push!(AST_values, AST)

        # Extract small intestine (index 13)
        ASI = u[13]
        push!(ASI_values, ASI)

        # Extract feces (index 15)
        Afeces = u[15]
        push!(Afeces_values, Afeces)

        # Extract liver (index 16)
        AL = u[16]
        push!(AL_values, AL)
    end

    return CA_values, Aplas_free_values, Aplas_free_over_Free_values, AR_values, AKb_values, APTC_values, Afil_values, Aurine_values, AST_values, ASI_values, Afeces_values, AL_values, sol.t
end

# Solve ODE
println("\n" * "="^80)
println("SOLVING ODE")
println("="^80)

observation_times = data.Time_days

CA_pred, Aplas_free_pred, Aplas_free_over_Free_pred, AR_pred, AKb_pred, APTC_pred, Afil_pred, Aurine_pred, AST_pred, ASI_pred, Afeces_pred, AL_pred, time_pred = run_pbk_single_dose(
    RAFapi,
    Km_apical,
    observation_times,
    oral_dose_ug,
    background_dose_ng_kg_day,
    BW
)

println("\n  ODE solved successfully")
println("  Time points in solution: $(length(time_pred))")
println("  Final rest of body amount: $(round(AR_pred[end], digits=4)) μg")
println("  Final kidney blood amount: $(round(AKb_pred[end], digits=6)) μg")
println("  Final PTC amount: $(round(APTC_pred[end], digits=6)) μg")
println("  Final filtrate amount: $(round(Afil_pred[end], digits=6)) μg")
println("  Final urine excretion: $(round(Aurine_pred[end], digits=4)) μg")
println("  Final feces excretion: $(round(Afeces_pred[end], digits=4)) μg")
println("  Final stomach amount: $(round(AST_pred[end], digits=6)) μg")
println("  Final small intestine amount: $(round(ASI_pred[end], digits=6)) μg")
println("  Final liver amount: $(round(AL_pred[end], digits=6)) μg")

# Create plot
println("\n" * "="^80)
println("CREATING PLOT")
println("="^80)

# Plasma concentration plot
p = plot(time_pred, CA_pred,
         label="Model prediction",
         xlabel="Time (days)",
         ylabel="Plasma Concentration (μg/L)",
         title="PFHxA Single-Dose: Model vs Observed",
         linewidth=2,
         color=:blue,
         size=(900, 600),
         legend=:topright)

scatter!(p, observation_times, data.C_plasma,
        label="Observed data",
        markersize=5,
        markercolor=:red,
        markerstrokewidth=0)

output_file = "PFAS_models/PFHxA/check_odes_plasma_plot.png"
savefig(p, output_file)
println("\n  Plasma concentration plot saved to: $output_file")

# Stomach compartment plot (u[11])
p_stomach = plot(time_pred, AST_pred,
                label="Stomach amount (u[11])",
                xlabel="Time (days)",
                ylabel="Amount in Stomach (μg)",
                title="PFHxA Amount in Stomach Compartment (u[11])",
                linewidth=2,
                color=:purple,
                size=(900, 600),
                legend=:topright)

output_file_stomach = "PFAS_models/PFHxA/check_odes_stomach_plot.png"
savefig(p_stomach, output_file_stomach)
println("  Stomach compartment plot saved to: $output_file_stomach")

# Small intestine compartment plot (u[13])
p_intestine = plot(time_pred, ASI_pred,
                  label="Small intestine amount (u[13])",
                  xlabel="Time (days)",
                  ylabel="Amount in Small Intestine (μg)",
                  title="PFHxA Amount in Small Intestine Compartment (u[13])",
                  linewidth=2,
                  color=:orange,
                  size=(900, 600),
                  legend=:topright)

output_file_intestine = "PFAS_models/PFHxA/check_odes_intestine_plot.png"
savefig(p_intestine, output_file_intestine)
println("  Small intestine compartment plot saved to: $output_file_intestine")

# Liver compartment plot (u[16])
p_liver = plot(time_pred, AL_pred,
              label="Liver amount (u[16])",
              xlabel="Time (days)",
              ylabel="Amount in Liver (μg)",
              title="PFHxA Amount in Liver Compartment (u[16])",
              linewidth=2,
              color=:green,
              size=(900, 600),
              legend=:topright)

output_file_liver = "PFAS_models/PFHxA/check_odes_liver_plot.png"
savefig(p_liver, output_file_liver)
println("  Liver compartment plot saved to: $output_file_liver")

# Urine compartment plot (u[10])
p_urine = plot(time_pred, Aurine_pred,
              label="Cumulative urine excretion (u[10])",
              xlabel="Time (days)",
              ylabel="Cumulative Amount in Urine (μg)",
              title="PFHxA Cumulative Urine Excretion (u[10])",
              linewidth=2,
              color=:cyan,
              size=(900, 600),
              legend=:bottomright)

output_file_urine = "PFAS_models/PFHxA/check_odes_urine_plot.png"
savefig(p_urine, output_file_urine)
println("  Urine compartment plot saved to: $output_file_urine")

# Feces compartment plot (u[15])
p_feces = plot(time_pred, Afeces_pred,
              label="Cumulative feces excretion (u[15])",
              xlabel="Time (days)",
              ylabel="Cumulative Amount in Feces (μg)",
              title="PFHxA Cumulative Fecal Excretion (u[15])",
              linewidth=2,
              color=:brown,
              size=(900, 600),
              legend=:bottomright)

output_file_feces = "PFAS_models/PFHxA/check_odes_feces_plot.png"
savefig(p_feces, output_file_feces)
println("  Feces compartment plot saved to: $output_file_feces")

# u[18]/Free plot
p_u18_over_Free = plot(time_pred, Aplas_free_over_Free_pred,
                      label="u[18]/Free",
                      xlabel="Time (days)",
                      ylabel="u[18]/Free (μg)",
                      title="PFHxA: u[18]/Free Over Time",
                      linewidth=2,
                      color=:red,
                      size=(900, 600),
                      legend=:topright)

output_file_u18 = "PFAS_models/PFHxA/check_odes_u18_over_Free_plot.png"
savefig(p_u18_over_Free, output_file_u18)
println("  u[18]/Free plot saved to: $output_file_u18")

# Rest of body compartment plot (u[1])
p_rest = plot(time_pred, AR_pred,
             label="Rest of body amount (u[1])",
             xlabel="Time (days)",
             ylabel="Amount in Rest of Body (μg)",
             title="PFHxA Amount in Rest of Body Compartment (u[1])",
             linewidth=2,
             color=:magenta,
             size=(900, 600),
             legend=:topright)

output_file_rest = "PFAS_models/PFHxA/check_odes_rest_body_plot.png"
savefig(p_rest, output_file_rest)
println("  Rest of body compartment plot saved to: $output_file_rest")

# Kidney blood compartment plot (u[4])
p_kidney_blood = plot(time_pred, AKb_pred,
                     label="Kidney blood amount (u[4])",
                     xlabel="Time (days)",
                     ylabel="Amount in Kidney Blood (μg)",
                     title="PFHxA Amount in Kidney Blood Compartment (u[4])",
                     linewidth=2,
                     color=:darkblue,
                     size=(900, 600),
                     legend=:topright)

output_file_kb = "PFAS_models/PFHxA/check_odes_kidney_blood_plot.png"
savefig(p_kidney_blood, output_file_kb)
println("  Kidney blood compartment plot saved to: $output_file_kb")

# Proximal tubule cells compartment plot (u[8])
p_ptc = plot(time_pred, APTC_pred,
            label="PTC amount (u[8])",
            xlabel="Time (days)",
            ylabel="Amount in Proximal Tubule Cells (μg)",
            title="PFHxA Amount in Proximal Tubule Cells (u[8])",
            linewidth=2,
            color=:darkgreen,
            size=(900, 600),
            legend=:topright)

output_file_ptc = "PFAS_models/PFHxA/check_odes_PTC_plot.png"
savefig(p_ptc, output_file_ptc)
println("  PTC compartment plot saved to: $output_file_ptc")

# Filtrate compartment plot (u[9])
p_filtrate = plot(time_pred, Afil_pred,
                 label="Filtrate amount (u[9])",
                 xlabel="Time (days)",
                 ylabel="Amount in Filtrate (μg)",
                 title="PFHxA Amount in Filtrate Compartment (u[9])",
                 linewidth=2,
                 color=:lightblue,
                 size=(900, 600),
                 legend=:topright)

output_file_fil = "PFAS_models/PFHxA/check_odes_filtrate_plot.png"
savefig(p_filtrate, output_file_fil)
println("  Filtrate compartment plot saved to: $output_file_fil")

println("\n" * "="^80)
println("ODE CHECK COMPLETE")
println("="^80)
println("\nTo test different parameter values:")
println("  1. Edit lines 60-61 in this file")
println("  2. Change RAFapi and/or Km_apical values")
println("  3. Run the script again: julia PFAS_models/PFHxA/check_odes.jl")
println("="^80)