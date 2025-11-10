using DifferentialEquations
using BenchmarkTools
using CSV
using DataFrames
using Statistics
using Sundials
using UnPack

# Define the create.params function (adapted from original)
function create_params(user_input)
    @unpack BW, substance, admin_dose, admin_time, f_unabs, ingestion, ingestion_time, admin_type, exp_type = user_input

    QCC = 12.5		  # Cardiac blood output (L/h/kg^0.75)
    QFC = 0.052		  # Fraction cardiac output going to fat
    QLC = 0.19 #0.25 from Loccisano 2011
    # Fraction cardiac output going to liver
    QKC = 0.175		  # Fraction cardiac output going to kiney
    #QFilC = 0.035		# Fraction cardiac output to the filtrate
    # compartment (20% of kiney blood flow)
    QGC = 0.181		  # Fraction cardiac output going to gut
    QLuC = 0.034    # Fraction cardiac output going to lungs
    QBC = 0.117    # Fraction cardiac output going to brain
    
    #fractional tissue volumes
    VLC = 0.023		  # Fraction liver volume
    VFC = 0.02 #0.214		  # Fraction fat volume
    VKC = 0.004		  # Fraction kidney volume
    VFilC = 0.0004	# Fraction filtrate compartment volume (10% of kidney volume)
    VGC = 0.0171		# Fraction gut volume
    VPlasC = 0.0428	# Fraction plasma volume
    VLuC = 0.014	# Fraction lungs volume
    VBC = 0.021	# Fraction lungs volume

    # Scaling parameters
    QC = QCC*BW^0.75*24	#Cardiac output (L/day)
    Htc = 0.44      #hematocrit
    QCP = QC*(1-Htc)	# Plasma flow
    QL = QLC*QCP			# Plasma flow to liver (L/day)
    QF = QFC*QCP			# Plasma flow to fat (L/day)
    QK = QKC*QCP	    # Plasma flow to kiney (L/day)
    QFil = 0.2*QK   	# Plasma flow to filtrate compartment (L/day)# 20% of QK
    QG = QGC*QCP	    # Plasma flow to gut (L/day)
    QLu = QLuC*QCP	  # Plasma flow to lungs (L/day)
    QB = QBC*QCP	    # Plasma flow to brain (L/day)
    QR = QCP - QL - QF - QK - QFil - QG - QLu - QB
    # Plasma flow to rest of the boy (L/day)
    
    Qbal = QCP - (QL+QF+QK+QFil+QG+QR+QLu+QB)        # balance check 
    
    VL = VLC*BW			    # Liver volume (L)
    VF = VFC*BW			    # Fat volume (L)
    VK = VKC*BW			    # Kiney volume (L)
    VFil = VFilC*BW	    # Fitrate compartment volume (L)
    VG = VGC*BW			    # Gut volume (L)
    VPlas = VPlasC*BW		# Plasma volume (L)
    VLu = VLuC*BW			  # Lungs volume (L)
    VB = VBC*BW			    # Brain volume (L)
    VR = 1*BW - VL - VF - VK - VFil - VG -VLu -VB - VPlas
    # Rest of the boy volume (L) # Loccisano 2011 uses 0.84*BW
    
    Vbal = (1*BW)-(VL+VF+VK+VFil+VG+VPlas+VR+VLu+VB) 

    kinetic_parameters = Dict(
      "PFBS" => (Liver=128.8, Brain=201.6, Lung=56.11, Kidney=6.27, Tm=6.1, Kt=5, Free=0.001),
      "PFHxA" => (Liver=0.001, Brain=43.68, Lung=31.92, Kidney=11.57, Tm=245.6, Kt=0.6, Free=0.01)
    )
    keep_params = kinetic_parameters[substance]

    Tm = keep_params.Tm #ug/h Fabrega (2015) Table 2
    Kt = keep_params.Kt #ug/L Fabrega (2015) Table 2
    Free = keep_params.Free #unitless Fabrega (2015) Table 2
    PL = keep_params.Liver
    # Partition Coefficient for Liver: Fabrega (2015) Table 2
    PF = 0.467 # Partition Coefficient for Fat: Fabrega (2014) Table 1
    PB = keep_params.Brain
    # Partition Coefficient for Brain: Fabrega (2015) Table 2
    PLu = keep_params.Lung
    # Partition Coefficient for Lung: Fabrega (2015) Table 2
    PK = keep_params.Kidney
    # Partition Coefficient for Kidney: Fabrega (2015) Table 2
    PG = 0.05 # Partition Coefficient for Gut: Loccisano (2011) Table 1
    PR = 0.12 # Partition Coefficient for Rest of body: Loccisano (2011) Table 1
    kurinec = 3e-04	#urinary elimination rate constant  (/h/kg^-0.25); 
    kurine = kurinec*BW^(-0.25)*24 # Elimination rate (1/day)

    return(QC=QC, QCP=QCP, QL=QL, QF=QF, QK=QK, 
            QFil=QFil, QG=QG, QLu=QLu, QB=QB, QR=QR,
            VPlas=VPlas, VL=VL, VF=VF, VK=VK, 
            VFil=VFil, VG=VG, VLu=VLu, VB=VB, VR=VR,
            PL=PL, PF=PF, PB=PB, PLu=PLu, PK=PK, 
            PG=PG, PR=PR,
            Tm=Tm, Kt=Kt, Free=Free, kurine=kurine,
            ingestion=ingestion, ingestion_time=ingestion_time,
            admin_dose=admin_dose, admin_time=admin_time,
            admin_type=admin_type, exp_type=exp_type,
            f_unabs=f_unabs)
end

function create_inits(parameters)
    # APlas, AG, AL, AF, ALu, AB, AK, AFil, AStore, AUrine, AR, ingestion
    return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
end

function create_events(parameters)
    @unpack admin_dose, admin_time, admin_type, exp_type, ingestion, ingestion_time, Free = parameters
    
    events = []
    
    if admin_type == "iv"
        # Calculate number of administrated doses and corresponding administration time for IV
        ldose = length(admin_dose)
        ltimes = length(admin_time)
        # If not equal, then stop 
        if ltimes != ldose
            error("The times of administration should be equal in number to the doses")
        else
            # IV: Add dose*Free to plasma compartment (index 1)
            for (dose, time) in zip(admin_dose, admin_time)
                affect!(integrator) = integrator.u[1] += dose * Free
                push!(events, PresetTimeCallback(time, affect!))
            end
        end
        
    elseif admin_type == "oral"
        lingest = length(ingestion)
        lingesttimes = length(ingestion_time)
        # If not equal, then stop 
        if lingest != lingesttimes
            error("The times of ingestion rate change should be equal to the ingestion time vector")
        else
            if exp_type == "pharmacokinetics"
                # For pharmacokinetic studies, add the dose directly to the gut compartment (index 2)
                for (dose, time) in zip(ingestion, ingestion_time)
                    affect!(integrator) = integrator.u[2] += dose
                    push!(events, PresetTimeCallback(time, affect!))
                end
            elseif exp_type == "biomonitoring"
                # For continuous exposure studies, set the ingestion rate (index 12)
                for (rate, time) in zip(ingestion, ingestion_time)
                    affect!(integrator) = integrator.u[12] = rate
                    push!(events, PresetTimeCallback(time, affect!))
                end
            end
        end
    end
    
    # Return empty callback set if no events, otherwise combine all events
    return length(events) > 0 ? CallbackSet(events...) : nothing
end

function ode_func(du, u, p, t)
    @unpack QC, QCP, QL, QF, QK, QFil, QG, QLu, QB, QR,
            VPlas, VL, VF, VK, VFil, VG, VLu, VB, VR,
            PL, PF, PB, PLu, PK, PG, PR,
            Tm, Kt, Free, kurine = p
    
    # Amount of PFAS in each compartment (ug)
    APlas = u[1] # Plasma mass (ug)
    AG = u[2] # Gut mass (ug)
    AL = u[3] # Liver mass (ug)
    AF = u[4] # Fat mass (ug)
    ALu = u[5] # Lung mass (ug)
    AB = u[6] # Brain mass (ug)
    AK = u[7] # Kidney mass (ug)
    AFil = u[8] # Filtrate compartment mass (ug)
    AStore = u[9] # Storage compartment mass (ug)
    AUrine = u[10] # Urine mass (ug)
    AR = u[11] # Rest of body mass (ug)
    ingestion = u[12] # Ingestion mass (ug)

    # Concentrations 
    CPlas_free = APlas / VPlas # Free plasma concentration (ug/L)
    CPlas = CPlas_free / Free # Total plasma concentration (ug/L)
    CG = AG / VG # Gut concentration (ug/L)
    CL = AL / VL # Liver concentration (ug/L)
    CF = AF / VF # Fat concentration (ug/L)
    CLu = ALu / VLu # Lung concentration (ug/L)
    CB = AB / VB # Brain concentration (ug/L)
    CK = AK / VK # Kidney concentration (ug/L)
    CFil = AFil / VFil # Filtrate compartment concentration (ug/L)
    CR = AR / VR # Rest of body concentration (ug/L)

    # Differential equations
    # Plasma compartment
    du[1] = QF*Free*CF/PF + (QL+QG)*Free*CL/PL + QR*Free*CR/PR + QK*Free*CK/PK + 
    QB*Free*CB/PB + QLu*Free*CLu/PLu - QCP*Free*CPlas

    # Gut compartment
    du[2] = QG*Free*(CPlas - CG/PG) + ingestion

    # Liver compartment
    du[3] = QL*Free*CPlas + QG*Free*CG/PG - (QL+QG)*Free*CL/PL

    # Fat compartment
    du[4] = QF*Free*(CPlas - CF/PF)

    # Lung compartment
    du[5] = QLu*Free*(CPlas - CLu/PLu)

    # Brain compartment
    du[6] = QB*Free*(CPlas - CB/PB)

    # Kidney compartment
    du[7] = QK*Free*(CPlas - CK/PK) + (Tm*CFil)/(Kt+CFil)

    # Filtrate compartment
    du[8] = QFil*(Free*CPlas - CFil) - (Tm*CFil)/(Kt+CFil)

    # Storage compartment for Urine
    du[9] = QFil*CFil - kurine*AStore

    # Urine compartment
    du[10] = kurine*AStore

    # Rest of body compartment
    du[11] = QR*Free*(CPlas - CR/PR)

    # Ingestion compartment
    du[12] = 0 # Ingestion is not modeled as a differential equation,
    # but rather as a source term in the gut compartment
end

# Post-process to get Concentrations
function extract_concentrations(sol, parameters)
    @unpack VPlas, VG, VL, VF, VLu, VB, VK, VFil, VR, Free = parameters
    concentrations = Dict(
        "CPlas" => [u[1] for u in sol.u] ./ VPlas ./ Free,
        "CG" => [u[2] for u in sol.u] ./ VG,
        "CL" => [u[3] for u in sol.u] ./ VL,
        "CF" => [u[4] for u in sol.u] ./ VF,
        "CLu" => [u[5] for u in sol.u] ./ VLu,
        "CB" => [u[6] for u in sol.u] ./ VB,
        "CK" => [u[7] for u in sol.u] ./ VK,
        "CFil" => [u[8] for u in sol.u] ./ VFil,
        "CR" => [u[11] for u in sol.u] ./ VR
    )
    return concentrations
end

# Benchmark configuration
struct BenchmarkScenario
    name::String
    substance::String
    tspan::Tuple{Float64, Float64}
    saveat::Float64
    ingestion::Vector{Float64}
    ingestion_time::Vector{Float64}
    exp_type::String
    description::String
end

# Define the four scenarios
scenarios = [
    BenchmarkScenario("Scenario1", "PFHxA", (0.0, 10.0), 1.0, [100.0], [0.01], "biomonitoring", "10 days, timestep 1 day, PFHxA"),
    BenchmarkScenario("Scenario2", "PFHxA", (0.0, 10.0), 0.1, [100.0], [0.01], "biomonitoring", "10 days, timestep 0.1 day, PFHxA"), 
    BenchmarkScenario("Scenario3", "PFHxA", (0.0, 5*365.0), 1.0, [100.0], [0.01], "biomonitoring", "5 years, timestep 1 day, PFHxA"),
    BenchmarkScenario("Scenario4", "PFBS", (0.0, 365.0), 1.0, collect(700.0 * 70 for _ in 0:30:365), collect(0:30.0:365), "pharmacokinetics", "1 year, monthly dosing, PFBS")
]

# Solvers to test
solvers = [
    (Tsit5(), "Tsit5"),
    (CVODE_Adams(), "CVODE_Adams"),
    (CVODE_BDF(), "CVODE_BDF"),
    (Rodas4(), "Rodas4"), 
    (Vern7(), "Vern7"),
    (TRBDF2(), "TRBDF2")
]

# Model parameters
BW = 70.0  # Body weight in kg
f_unabs = 0.5

function create_user_input(scenario)
    return Dict(
        "BW" => BW,
        "substance" => scenario.substance,
        "f_unabs" => f_unabs,
        "admin_dose" => 0,
        "admin_time" => 0,
        "ingestion" => scenario.ingestion,
        "ingestion_time" => scenario.ingestion_time,
        "admin_type" => "oral",
        "exp_type" => scenario.exp_type
    )
end

function run_benchmark(scenario, solver_alg, solver_name)
    println("Running $(scenario.name) with $(solver_name)...")
    
    # Set up user input and parameters
    user_input = create_user_input(scenario)
    params = create_params(user_input)
    inits = create_inits(params)
    events = create_events(params)
    
    # Set up ODE problem
    prob = ODEProblem(ode_func, inits, scenario.tspan, params, callback=events)
    
    # Benchmark the solver
    benchmark_result = @benchmark solve($prob, $solver_alg, saveat=$(scenario.saveat), reltol=1e-8, abstol=1e-10, maxiters=1000000) samples=3 seconds=30
    
    # Extract timing statistics
    times = benchmark_result.times / 1e6  # Convert to milliseconds
    
    return Dict(
        "scenario" => scenario.name,
        "substance" => scenario.substance,
        "solver" => solver_name,
        "description" => scenario.description,
        "mean_time_ms" => mean(times),
        "median_time_ms" => median(times),
        "std_time_ms" => std(times),
        "min_time_ms" => minimum(times),
        "max_time_ms" => maximum(times),
        "memory_mb" => benchmark_result.memory / 1024^2,
        "allocations" => benchmark_result.allocs,
        "n_samples" => length(times)
    )
end

function verify_solution_accuracy(scenario, solver_alg)
    user_input = create_user_input(scenario)
    params = create_params(user_input)
    inits = create_inits(params)
    events = create_events(params)
    
    prob = ODEProblem(ode_func, inits, scenario.tspan, params, callback=events)
    sol = solve(prob, solver_alg, saveat=scenario.saveat, reltol=1e-8, abstol=1e-10, maxiters=1000000)
    
    # Extract final plasma concentration
    concentrations = extract_concentrations(sol, params)
    final_cplas = concentrations["CPlas"][end]
    
    return Dict(
        "scenario" => scenario.name,
        "substance" => scenario.substance,
        "final_plasma_conc" => final_cplas,
        "n_timepoints" => length(sol.t),
        "final_time" => sol.t[end]
    )
end

# Main benchmarking function
function run_all_benchmarks()
    println("Starting Julia PBK benchmarking...")
    println("Model: 12-compartment PFAS PBK model")
    println("=" ^ 60)
    
    results = []
    accuracy_results = []
    
    # Run benchmarks for all scenario-solver combinations
    for scenario in scenarios
        println("\n--- $(scenario.name): $(scenario.description) ---")
        
        for (solver_alg, solver_name) in solvers
            # Skip explicit solvers for stiff Scenario 4
            if scenario.name == "Scenario4" && solver_name in ["Tsit5", "Vern7"]
                println("  Skipping $(solver_name) for $(scenario.name) (explicit solver not suitable for stiff problem)")
                continue
            end
            
            try
                # Benchmark performance
                result = run_benchmark(scenario, solver_alg, solver_name)
                push!(results, result)
                
                println("  $(solver_name): $(round(result["mean_time_ms"], digits=3)) ms (mean)")
                
                # Test accuracy (only once per scenario, using appropriate solver)
                accuracy_solver = (scenario.name == "Scenario4") ? "Rodas4" : "Tsit5"
                if solver_name == accuracy_solver
                    accuracy = verify_solution_accuracy(scenario, solver_alg)
                    accuracy["solver"] = solver_name
                    push!(accuracy_results, accuracy)
                    println("  Final plasma conc: $(round(accuracy["final_plasma_conc"], digits=6)) ng/mL")
                end
                
            catch e
                println("  ERROR with $(solver_name): $e")
                push!(results, Dict(
                    "scenario" => scenario.name,
                    "substance" => scenario.substance,
                    "solver" => solver_name,
                    "description" => scenario.description,
                    "error" => string(e)
                ))
            end
        end
    end
    
    # Save results
    results_df = DataFrame(results)
    accuracy_df = DataFrame(accuracy_results)
    
    CSV.write("julia_pbk_benchmark_results.csv", results_df)
    CSV.write("julia_pbk_accuracy_results.csv", accuracy_df)
    
    println("\n" * "=" ^ 60)
    println("PBK Benchmarking complete!")
    println("Results saved to:")
    println("  - julia_pbk_benchmark_results.csv")
    println("  - julia_pbk_accuracy_results.csv")
    
    return results_df, accuracy_df
end

# Run the benchmarks
if abspath(PROGRAM_FILE) == @__FILE__
    results_df, accuracy_df = run_all_benchmarks()
end