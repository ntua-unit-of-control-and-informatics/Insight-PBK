using DataFrames
using CSV
using Plots
using DifferentialEquations
using UnPack
using Sundials


include("../jl_PBK_model/julia_model/pbk_functions.jl")

println("PFBS Scenario Verification")
println("="^40)

zhou_serum = CSV.read("PFBS_recalibration/data_preprocessing/zhou_2014_serum_summary.csv", DataFrame)
zhou_urine = CSV.read("PFBS_recalibration/data_preprocessing/zhou_2014_urine_summary.csv", DataFrame)

# Load exposures
abraham_results = CSV.read("Abraham_model/Abraham_model_results.csv", DataFrame)
exposure_zhou = abraham_results[abraham_results.Study .== "Zhou et al., 2014" .&& abraham_results.PFAS .== "PFBS", :Estimated_Intake_ng_kg_bw_day][1]

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

# Study data structure
study = [
    (name="Zhou", serum=zhou_serum, urine=zhou_urine, exposure=exposure_zhou, exp_time=10.0*365)
]


# Mean parameter values (same as calibration script)
Tm_mean = 6.1  # Literature value used as mean
Kt_fixed = 5.0  # Fixed value

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

# Extract concentrations with urine
function extract_concentrations_with_urine(sol, parameters)
    @unpack VPlas, Free = parameters
    
    # Serum concentrations
    serum_conc = [u[1] for u in sol.u] ./ VPlas ./ Free
    
    # Urine concentrations
    urine_production_rate = 1.44  # L/day
    AUrine_cumulative = [u[10] for u in sol.u]
    
    urine_concentrations = Float64[]
    for i in 1:length(sol.t)
        if i == 1
            time_elapsed = sol.t[i]
            total_volume = urine_production_rate * time_elapsed
            urine_conc = total_volume > 0 ? AUrine_cumulative[i] / total_volume : 0.0
        else
            time_interval = sol.t[i] - sol.t[i-1]
            volume_interval = urine_production_rate * time_interval
            amount_interval = AUrine_cumulative[i] - AUrine_cumulative[i-1]
            urine_conc = volume_interval > 0 ? amount_interval / volume_interval : 0.0
        end
        push!(urine_concentrations, urine_conc)
    end
    
    return Dict("CPlas" => serum_conc, "CUrine" => urine_concentrations)
end

study = study[1]  # Only one study in this example
max_time = maximum(vcat(study.serum.time_days, study.urine.time_days))
time_points = collect(range(0.0, max_time, length=100))

user_input = (BW = 70.0, substance = "PFBS", admin_dose = [0.0], admin_time = [0.0],
                f_unabs = 0.0, ingestion = [study.exposure * 70.0 * 1e-03], ingestion_time = [0.01], 
                admin_type = "oral", exp_type = "biomonitoring")

# tspan = (0.0, max_time)
tspan = (0.0, max_time)
saveat_times = time_points
plot_times = time_points

# Create parameters
constant_params = create_params_updated(user_input, new_partition_coeffs)
params = (constant_params..., Tm = Tm_mean, Kt = Kt_fixed)

# Solve ODE
inits = create_inits(params)
inits[12] = 1.148287
events = create_events(params)
println("Event schedule from callbacks:")
prob = ODEProblem(ode_func, inits, tspan, params)
sol = solve(prob, CVODE_BDF(), saveat=saveat_times, reltol=1e-6, abstol=1e-6)

# Debug: Check solution times
println("  tspan: $tspan")
println("  saveat_times length: $(length(saveat_times))")
println("  saveat_times range: $(saveat_times[1]) to $(saveat_times[end])")
println("  sol.t length: $(length(sol.t))")
println("  sol.t range: $(sol.t[1]) to $(sol.t[end])")