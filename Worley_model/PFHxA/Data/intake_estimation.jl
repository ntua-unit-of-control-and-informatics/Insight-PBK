using UnPack
using DataFrames
using CSV
using Plots

# Abraham kinetic model for PFAS
# Vd = Volume of distribution (ml/kg BW)
# k_el = Elimination rate constant (1/day)
PFAS_kinetic_params = Dict(
    "PFHxA" => (Vd = 119, k_el = 0.478),
    "PFBS" => (Vd = 137, k_el = 0.0137)
)

function inverse_abraham_kinetic_model(pfas_type, params, Css, half_life)
    @unpack Vd, k_el= params
    estimated_intake = log(2.0)/half_life * Css *Vd
    return estimated_intake
end

function run_inverse_abraham_kinetic_model(pfas_type, Css, half_life)
    params = PFAS_kinetic_params[pfas_type]
    if isnothing(params)
        error("PFAS type not recognized. Available types: $(keys(PFAS_kinetic_params))")
    end
    estimated_intake = inverse_abraham_kinetic_model(pfas_type, params, Css, half_life)
    return estimated_intake
end

half_life_data = CSV.read("Worley_model/PFHxA/Data/PFHxA_half_life_data.csv", DataFrame)
# Pre-allocate the column with the correct type
half_life_data.Estimated_Intake_ng_kg_bw_day = zeros(nrow(half_life_data))
Css = 10.0 # ug/L

for (i, row) in enumerate(eachrow(half_life_data))
    pfas_type = "PFHxA"
    half_life = row.Half_life_days
    estimated_intake = run_inverse_abraham_kinetic_model(pfas_type, Css, half_life)
    half_life_data.Estimated_Intake_ng_kg_bw_day[i] = estimated_intake
end

# Write results to CSV
CSV.write("Worley_model/PFHxA/Data/intake_estimation_results.csv", half_life_data)


