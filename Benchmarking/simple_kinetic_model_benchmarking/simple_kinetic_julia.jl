using DifferentialEquations
using BenchmarkTools
using CSV
using DataFrames
using Statistics
using Sundials

# Simple kinetic model: dC/dt = -k_el * C
function simple_kinetic!(du, u, p, t)
    k_el = p
    C = u[1]
    du[1] = -k_el * C
end

# Benchmark configuration
struct BenchmarkScenario
    name::String
    tspan::Tuple{Float64, Float64}
    saveat::Float64
    description::String
end

# Define the three scenarios
scenarios = [
    BenchmarkScenario("Scenario1", (0.0, 10.0), 1.0, "10 days, timestep 1 day"),
    BenchmarkScenario("Scenario2", (0.0, 10.0), 0.1, "10 days, timestep 0.1 day"), 
    BenchmarkScenario("Scenario3", (0.0, 5*365.0), 1.0, "5 years, timestep 1 day")
]

# Solvers to test
solvers = [
    (Tsit5(), "Tsit5"),
    (CVODE_Adams(), "CVODE_Adams"),
    (Rodas4(), "Rodas4"), 
    (RadauIIA5(), "RadauIIA5")
]

# Model parameters
k_el = 0.1  # elimination rate constant (1/day)
C0 = [100.0]  # initial concentration

function run_benchmark(scenario, solver_alg, solver_name)
    println("Running $(scenario.name) with $(solver_name)...")
    
    # Set up ODE problem
    prob = ODEProblem(simple_kinetic!, C0, scenario.tspan, k_el)
    
    # Benchmark the solver
    benchmark_result = @benchmark solve($prob, $solver_alg, saveat=$(scenario.saveat), reltol=1e-8, abstol=1e-10) samples=10 seconds=10
    
    # Extract timing statistics
    times = benchmark_result.times / 1e6  # Convert to milliseconds
    
    return Dict(
        "scenario" => scenario.name,
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
    prob = ODEProblem(simple_kinetic!, C0, scenario.tspan, k_el)
    sol = solve(prob, solver_alg, saveat=scenario.saveat, reltol=1e-8, abstol=1e-10)
    
    # Analytical solution: C(t) = C0 * exp(-k_el * t)
    analytical_final = C0[1] * exp(-k_el * scenario.tspan[2])
    numerical_final = sol.u[end][1]
    relative_error = abs(numerical_final - analytical_final) / analytical_final
    
    return Dict(
        "scenario" => scenario.name,
        "analytical_final" => analytical_final,
        "numerical_final" => numerical_final,
        "relative_error" => relative_error,
        "n_timepoints" => length(sol.t)
    )
end

# Main benchmarking function
function run_all_benchmarks()
    println("Starting Julia benchmarking for simple kinetic model...")
    println("Model: dC/dt = -k_el * C, where k_el = $k_el /day, C0 = $(C0[1])")
    println("=" ^ 60)
    
    results = []
    accuracy_results = []
    
    # Run benchmarks for all scenario-solver combinations
    for scenario in scenarios
        println("\n--- $(scenario.name): $(scenario.description) ---")
        
        for (solver_alg, solver_name) in solvers
            try
                # Benchmark performance
                result = run_benchmark(scenario, solver_alg, solver_name)
                push!(results, result)
                
                println("  $(solver_name): $(round(result["mean_time_ms"], digits=3)) ms (mean)")
                
                # Test accuracy (only once per scenario, using Tsit5)
                if solver_name == "Tsit5"
                    accuracy = verify_solution_accuracy(scenario, solver_alg)
                    accuracy["solver"] = solver_name
                    push!(accuracy_results, accuracy)
                    println("  Accuracy: $(round(accuracy["relative_error"]*100, digits=6))% relative error")
                end
                
            catch e
                println("  ERROR with $(solver_name): $e")
                push!(results, Dict(
                    "scenario" => scenario.name,
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
    
    CSV.write("julia_benchmark_results.csv", results_df)
    CSV.write("julia_accuracy_results.csv", accuracy_df)
    
    println("\n" * "=" ^ 60)
    println("Benchmarking complete!")
    println("Results saved to:")
    println("  - julia_benchmark_results.csv")
    println("  - julia_accuracy_results.csv")
    
    return results_df, accuracy_df
end

# Run the benchmarks
if abspath(PROGRAM_FILE) == @__FILE__
    results_df, accuracy_df = run_all_benchmarks()
end