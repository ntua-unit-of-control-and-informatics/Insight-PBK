using DifferentialEquations
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

function abraham_kinetic_model!(du, u, p, t)
    @unpack Vd, k_el, I = p
    C = u[1]  # Concentration of PFAS in the system
    du[1] = I/Vd - k_el*C  # Rate of change of concentration
end

function run_abraham_kinetic_model(pfas_type, intake, tspan, u0)
    params = PFAS_kinetic_params[pfas_type]
    if isnothing(params)
        error("PFAS type not recognized. Available types: $(keys(PFAS_kinetic_params))")
    end
    params = merge(params, (I = intake,))  # Initial concentration as input
    prob = ODEProblem(abraham_kinetic_model!, u0, tspan, params)
    sol = solve(prob, Tsit5())
    display(plot(sol.t, sol[1,:], xlabel="Time (days)", ylabel="Concentration (ng/ml)", 
        title="$pfas_type Concentration vs Time", linewidth=2))
    Css = sol[end]
    return Css
end

# example = run_abraham_kinetic_model("PFHxA", 100.0, (0.0, 100), [0.0])  # Example usage
# print(example)

# Simulation: PFBS (Zhou et al., 2014)
pfas_type = "PFBS"
intake = 8.33 # ng/kg_bw/day
tspan = (0.0, 6*365)  # 6 years in days
u0 = [0.0]  # Initial concentration
Css = run_abraham_kinetic_model(pfas_type, intake, tspan, u0)
println("Steady-state concentration of $pfas_type after 6 years: $Css ng/ml")

# Simulation: PFHxA (Zhou et al., 2014)
pfas_type = "PFHxA"
intake = 0.01 # ng/kg_bw/day
tspan = (0.0, 6*365)  # 6 years in days
u0 = [0.0]  # Initial concentration
Css = run_abraham_kinetic_model(pfas_type, intake, tspan, u0)
println("Steady-state concentration of $pfas_type after 6 years: $Css ng/ml")

