using DataFrames
using CSV
using Plots
using DifferentialEquations
using UnPack
using Sundials

include("../jl_PBK_model/julia_model/pbk_functions.jl")

println("PFBS Scenario Verification")
println("="^40)

# Load study data
olsen_serum = CSV.read("PFBS_recalibration/data_preprocessing/olsen_2009_serum_summary.csv", DataFrame)
olsen_urine = CSV.read("PFBS_recalibration/data_preprocessing/olsen_2009_urine_summary.csv", DataFrame)
he_serum = CSV.read("PFBS_recalibration/data_preprocessing/he_2023_serum_summary.csv", DataFrame)
he_urine = CSV.read("PFBS_recalibration/data_preprocessing/he_2023_urine_summary.csv", DataFrame)
zhou_serum = CSV.read("PFBS_recalibration/data_preprocessing/zhou_2014_serum_summary.csv", DataFrame)
zhou_urine = CSV.read("PFBS_recalibration/data_preprocessing/zhou_2014_urine_summary.csv", DataFrame)

# Load exposures
abraham_results = CSV.read("Abraham_model/Abraham_model_results.csv", DataFrame)
exposure_zhou = abraham_results[abraham_results.Study .== "Zhou et al., 2014" .&& abraham_results.PFAS .== "PFBS", :Estimated_Intake_ng_kg_bw_day][1]
exposure_olsen = abraham_results[abraham_results.Study .== "Olsen et al., 2009" .&& abraham_results.PFAS .== "PFBS", :Estimated_Intake_ng_kg_bw_day][1]
exposure_he = abraham_results[abraham_results.Study .== "He et al., 2023" .&& abraham_results.PFAS .== "PFBS", :Estimated_Intake_ng_kg_bw_day][1]

# Load partition coefficients
partition_coeffs = CSV.read("EDM_Partition_Coefficients/partition_coefficients_results.csv", DataFrame)
pfbs_coeffs = filter(row -> row.Compound == "PFBS", partition_coeffs)

new_partition_coeffs = Dict(
    "liver" => pfbs_coeffs[pfbs_coeffs.Tissue .== "liver", :Partition_Coefficient][1],
    "adipose" => pfbs_coeffs[pfbs_coeffs.Tissue .== "adipose", :Partition_Coefficient][1],
    "kidney" => pfbs_coeffs[pfbs_coeffs.Tissue .== "kidney", :Partition_Coefficient][1],
    "gut" => pfbs_coeffs[pfbs_coeffs.Tissue .== "gut", :Partition_Coefficient][1],
    "lung" => pfbs_coeffs[pfbs_coeffs.Tissue .== "lung", :Partition_Coefficient][1],
    "brain" => pfbs_coeffs[pfbs_coeffs.Tissue .== "brain", :Partition_Coefficient][1],
    "rest_of_body" => pfbs_coeffs[pfbs_coeffs.Tissue .== "rest of body", :Partition_Coefficient][1]
)
println("Using partition coefficients:")
for (tissue, coeff) in new_partition_coeffs
    println("  $tissue: $(round(coeff, digits=3))")
end

# Study data structure
studies = [
    (name="Zhou", serum=zhou_serum, urine=zhou_urine, exposure=exposure_zhou, exp_time=10.0*365),
    (name="Olsen", serum=olsen_serum, urine=olsen_urine, exposure=exposure_olsen, exp_time=10.0*365),
    (name="He", serum=he_serum, urine=he_urine, exposure=exposure_he, exp_time=10.0*365)
]

# Mean parameter values (same as calibration script)
Tm_mean = 13791.1552  # Literature value used as mean
Kt_fixed = 5.0  # Fixed value
PC_scale_large = 1.0
PC_scale_low = 1.0

# Create updated parameters function
function create_params_updated(user_input, new_PCs)
    params = create_params(user_input)
    return (params..., 
            PL = new_PCs["liver"],
            PF = new_PCs["adipose"], 
            PK = new_PCs["kidney"],
            PG = new_PCs["gut"],
            PLu = new_PCs["lung"],
            PB = new_PCs["brain"],
            PR = new_PCs["rest_of_body"])
end


# Run simulation for each study
for study in studies
    println("Simulating $(study.name)...")
    
    if study.name == "Olsen"
        # Olsen: Two-phase exposure (exposure + elimination)
        elimination_start = study.exp_time
        max_obs_time = maximum(vcat(study.serum.time_days, study.urine.time_days))
        total_time = elimination_start + max_obs_time
        study.serum.time_days .+= elimination_start
        study.urine.time_days .+= elimination_start

        # Simulation time points
        time_points = collect(range(0.0, elimination_start + max_obs_time, length=100))
        absolute_times = time_points
        
        user_input = (BW = 70.0, substance = "PFBS", admin_dose = [0.0], admin_time = [0.0],
                     f_unabs = 0.0, ingestion = [study.exposure * 70.0*1e-3, 0.0], 
                     ingestion_time = [0.01, elimination_start], admin_type = "oral", exp_type = "biomonitoring")
        
        tspan = (0.0, total_time)
        saveat_times = absolute_times
        plot_times = time_points  # Relative times for plotting
    else
        # Zhou/He: Continuous exposure
        max_time = maximum(vcat(study.serum.time_days, study.urine.time_days))
        time_points = collect(range(0.0, max_time, length=100))
        
        user_input = (BW = 70.0, substance = "PFBS", admin_dose = [0.0], admin_time = [0.0],
                     f_unabs = 0.0, ingestion = [study.exposure * 70.0*1e-3], ingestion_time = [0.01], 
                     admin_type = "oral", exp_type = "biomonitoring")
        
        tspan = (0.0, max_time)
        saveat_times = time_points
        plot_times = time_points
    end
    
    # Create parameters
    constant_params = create_params_updated(user_input, new_partition_coeffs)
    # constant_params = create_params(user_input)
    params = (constant_params..., Tm = Tm_mean, Kt = Kt_fixed)
    # Apply scaling to all partition coefficients
    # PL = 128.8, #constant_params.PL * PC_scale_large,
    # PF = 0.467,#constant_params.PF * PC_scale_low,
    # PK = 6.27,#constant_params.PK * PC_scale_low,
    # PG = 0.05,#constant_params.PG * PC_scale_low,
    # PLu = 56.11,#constant_params.PLu * PC_scale_large,
    # PB = 201.6,#constant_params.PB * PC_scale_large,
    # PR = 0.12)#constant_params.PR * PC_scale_low)
    # PL = constant_params.PL * PC_scale_large,
    # PF = constant_params.PF * PC_scale_low,
    # PK = constant_params.PK * PC_scale_low,
    # PG = constant_params.PG * PC_scale_low,
    # PLu = constant_params.PLu * PC_scale_large,
    # PB = constant_params.PB * PC_scale_large,
    # PR = constant_params.PR * PC_scale_low)

    # Solve ODE
    inits = create_inits(params)
    events = create_events(params)
    println("Event schedule from callbacks:")
    prob = ODEProblem(ode_func, inits, tspan, params, callback=events)
    sol = solve(prob, Rodas5P(), saveat=saveat_times, reltol=1e-6, abstol=1e-6)

    # Debug: Check solution times
    println("  tspan: $tspan")
    println("  saveat_times length: $(length(saveat_times))")
    println("  saveat_times range: $(saveat_times[1]) to $(saveat_times[end])")
    println("  sol.t length: $(length(sol.t))")
    println("  sol.t range: $(sol.t[1]) to $(sol.t[end])")
    
    # Extract concentrations using main function (now includes physiological urine concentration)
    concentrations = extract_concentrations(sol, params)

    # Create serum plot - show entire simulation for all studies
    p1 = plot(sol.t, concentrations["CPlas"], 
              label="Model prediction", lw=2, color=:blue,
              xlabel="Time (days)", ylabel="Serum concentration (ng/ml)",
              title="$(study.name) - Serum PFBS")
    scatter!(p1, study.serum.time_days, study.serum.mean_ng_ml,
            label="Observed", color=:red, ms=6)
    
    # Create urine plot - show entire simulation for all studies
    p2 = plot(sol.t, concentrations["CUrine"],
              label="Model prediction", lw=2, color=:blue,
              xlabel="Time (days)", ylabel="Urine concentration (ng/ml)",
              title="$(study.name) - Urine PFBS")
    scatter!(p2, study.urine.time_days, study.urine.mean_ng_ml,
             label="Observed", color=:red, ms=6)

    # Extract Astorage (compartment 9) across all time points
    filtrate_values = [u[8] for u in sol.u]

    p_Filtrate = plot(sol.t, filtrate_values,
              label="Model prediction", lw=2, color=:blue,
              xlabel="Time (days)", ylabel="Filtrate ug",
              title="$(study.name) - Filtrate PFBS")

    # Save plots
    savefig(p1, "PFBS_recalibration/Scenario_check/$(study.name)_serum_check.png")
    savefig(p2, "PFBS_recalibration/Scenario_check/$(study.name)_urine_check.png")
    savefig(p_Filtrate, "PFBS_recalibration/Scenario_check/$(study.name)_Filtrate_check.png")

    # Create simple plot with all compartments (different line styles) - add small offset for log scale
    ϵ = 1e-10  # Small offset to avoid log(0)
    p6 = plot(sol.t, [u[1] + ϵ for u in sol.u], lw=4, ls=:solid, label="Serum", size=(800, 600), color=:black, yaxis=:log, ylims=(1e-2, 1e6))
    plot!(p6, sol.t, [u[3] + ϵ for u in sol.u], lw=4, ls=:dash, label="Liver", color=:black)
    plot!(p6, sol.t, [u[7] + ϵ for u in sol.u], lw=4, ls=:dot, label="Kidney", color=:black)
    plot!(p6, sol.t, [u[4] + ϵ for u in sol.u], lw=4, ls=:dashdot, label="Fat", color=:black)
    plot!(p6, sol.t, [u[6] + ϵ for u in sol.u], lw=4, ls=:dashdotdot, label="Brain", color=:black)
    plot!(p6, sol.t, [u[5] + ϵ for u in sol.u], lw=3, ls=:solid, label="Lung", color=:gray)
    plot!(p6, sol.t, [u[11] + ϵ for u in sol.u], lw=3, ls=:dash, label="Rest", color=:gray)
    plot!(p6, sol.t, [u[10] + ϵ for u in sol.u], lw=2, ls=:auto, label="Urine", color=:red)
    plot!(p6, sol.t, [u[9] + ϵ for u in sol.u], lw=2, ls=:solid, label="Storage", color=:blue)
    plot!(p6, sol.t, [u[8] + ϵ for u in sol.u], lw=2, ls=:dash, label="Filtrate", color=:green)

    savefig(p6, "PFBS_recalibration/Scenario_check/$(study.name)_mass_compartments.png")

    # Create simple plot with all compartments (different line styles)
    p3 = plot(sol.t, concentrations["CPlas"], lw=4, ls=:solid, label="Serum", size=(800, 600), color=:black)
    plot!(p3, sol.t, concentrations["CL"], lw=4, ls=:dash, label="Liver", color=:black)
    plot!(p3, sol.t, concentrations["CK"], lw=4, ls=:dot, label="Kidney", color=:black)
    plot!(p3, sol.t, concentrations["CF"], lw=4, ls=:dashdot, label="Fat", color=:black)
    plot!(p3, sol.t, concentrations["CB"], lw=4, ls=:dashdotdot, label="Brain", color=:black)
    plot!(p3, sol.t, concentrations["CLu"], lw=2, ls=:solid, label="Lung", color=:black)
    plot!(p3, sol.t, concentrations["CR"], lw=2, ls=:dash, label="Rest", color=:black)
    plot!(p3, sol.t, concentrations["CUrine"], lw=2, ls=:dot, label="Urine", color=:black)

    savefig(p3, "PFBS_recalibration/Scenario_check/$(study.name)_compartments.png")
    
    # Create cumulative urine excretion plot (mass)
    # Extract cumulative urine amount (compartment 10, AUrine)
    cumulative_urine_mass = [u[10] for u in sol.u]  # μg
    
    p4 = plot(sol.t, cumulative_urine_mass,
              label="Cumulative excretion", lw=2, color=:green,
              xlabel="Time (days)", ylabel="Cumulative urine excretion (μg)",
              title="$(study.name) - Cumulative Urine Excretion")
    
    savefig(p4, "PFBS_recalibration/Scenario_check/$(study.name)_urine_mass.png")
    
    # Create plot for AStore + AFil (kidney compartments)
    # Extract AStore (compartment 9) and AFil (compartment 8)
    astore_mass = [u[9] for u in sol.u]  # μg - storage compartment
    afil_mass = [u[8] for u in sol.u]    # μg - filtrate compartment
    
    p5 = plot(sol.t, astore_mass, label="AStore", lw=2, color=:blue)
    plot!(p5, sol.t, afil_mass, label="AFil", lw=2, color=:red)
    plot!(p5, xlabel="Time (days)", ylabel="Amount (μg)",
          title="$(study.name) - Kidney Storage and Filtrate Amounts")
    
    savefig(p5, "PFBS_recalibration/Scenario_check/$(study.name)_kidney_amounts.png")

    # Print summary
    println("  $(study.name) simulation completed:")
    println("    Exposure: $(round(study.exposure, digits=1)) ng/kg-bw/day")
    println("    Exposure duration: $(round(study.exp_time, digits=0)) days")
    println("    Serum range: $(round(minimum(concentrations["CPlas"]), digits=3)) - $(round(maximum(concentrations["CPlas"]), digits=3)) ng/ml")
    println("    Urine range: $(round(minimum(concentrations["CUrine"]), digits=3)) - $(round(maximum(concentrations["CUrine"]), digits=3)) ng/ml")
    println("    Total urine excretion: $(round(maximum(cumulative_urine_mass), digits=3)) μg")
    println("    Plots saved: $(study.name)_serum_check.png, $(study.name)_urine_check.png, $(study.name)_compartments.png, $(study.name)_urine_mass.png, $(study.name)_kidney_amounts.png")
    println()
end

println("Scenario verification complete!")
println("Check the generated plots to verify simulation scenarios.")