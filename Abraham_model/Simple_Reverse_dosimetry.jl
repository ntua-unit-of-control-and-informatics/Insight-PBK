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

function inverse_abraham_kinetic_model(pfas_type, params, Css)
    @unpack Vd, k_el= params
    estimated_intake = k_el * Css *Vd
    return estimated_intake
end

function run_inverse_abraham_kinetic_model(pfas_type, Css)
    params = PFAS_kinetic_params[pfas_type]
    if isnothing(params)
        error("PFAS type not recognized. Available types: $(keys(PFAS_kinetic_params))")
    end
    estimated_intake = inverse_abraham_kinetic_model(pfas_type, params, Css)
    return estimated_intake
end

# Simulation: PFBS (Zhou et al., 2014)
pfas_type = "PFBS"
Css = 8.74 # ng/ml
estimated_intake_zhou_PFBS = run_inverse_abraham_kinetic_model(pfas_type, Css)
println("Estimated daily intake of $pfas_type based on Css of $Css ng/ml: $estimated_intake_zhou_PFBS ng/kg_bw/day")

# Simulation: PFBS (Olsen et al., 2009)
pfas_type = "PFBS"
Css = 397 # ng/ml 3 average concentration at study onset
estimated_intake_olsen_PFBS = run_inverse_abraham_kinetic_model(pfas_type, Css)
println("Estimated daily intake of $pfas_type based on Css of $Css ng/ml: $estimated_intake_olsen_PFBS ng/kg_bw/day")

# Simulation: PFBS (He et al., 2023)
pfas_type = "PFBS"
Css = 64.7 # ng/ml
estimated_intake_he_pfbs = run_inverse_abraham_kinetic_model(pfas_type, Css)
println("Estimated daily intake of $pfas_type based on Css of $Css ng/ml: $estimated_intake_he_pfbs ng/kg_bw/day")

# Collect all results in a DataFrame
results = DataFrame(
    Study = ["Zhou et al., 2014", "Olsen et al., 2009", "He et al., 2023"],
    PFAS = ["PFBS", "PFBS", "PFBS"],
    Css_ng_ml = [8.74, 397, 64.7],
    Estimated_Intake_ng_kg_bw_day = [
        estimated_intake_zhou_PFBS,
        estimated_intake_olsen_PFBS,
        estimated_intake_he_pfbs
    ]
)

# Write results to CSV
CSV.write("Abraham_model/Abraham_model_results.csv", results)
println("Results saved to Abraham_model/Abraham_model_results.csv")

