library(deSolve)
library(microbenchmark)

# Simple kinetic model: dC/dt = -k_el * C
simple_kinetic_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dC <- -k_el * C
    return(list(c(dC)))
  })
}

# Model parameters
k_el <- 0.1  # elimination rate constant (1/day)
C0 <- c(C = 100.0)  # initial concentration

# R solvers to test
solvers <- list(
  list(method = "lsoda", name = "lsoda")
)

# Benchmark scenarios
scenarios <- list(
  list(name = "Scenario1", tspan = c(0, 10), timestep = 1.0, 
       description = "10 days, timestep 1 day"),
  list(name = "Scenario2", tspan = c(0, 10), timestep = 0.1, 
       description = "10 days, timestep 0.1 day"),
  list(name = "Scenario3", tspan = c(0, 5*365), timestep = 1.0, 
       description = "5 years, timestep 1 day")
)

# Function to run benchmark for a single scenario and solver
run_benchmark <- function(scenario, solver) {
  cat("Running", scenario$name, "with", solver$name, "...\n")
  
  # Create time sequence
  times <- seq(scenario$tspan[1], scenario$tspan[2], by = scenario$timestep)
  
  # Benchmark the solver
  benchmark_result <- microbenchmark(
    result = ode(y = C0, times = times, func = simple_kinetic_model, 
                 parms = list(k_el = k_el), method = solver$method,
                 rtol = 1e-8, atol = 1e-10),
    times = 10,
    unit = "ms"
  )
  
  # Extract timing statistics
  times_ms <- benchmark_result$time / 1e6  # Convert to milliseconds
  
  return(list(
    scenario = scenario$name,
    solver = solver$name,
    description = scenario$description,
    mean_time_ms = mean(times_ms),
    median_time_ms = median(times_ms),
    std_time_ms = sd(times_ms),
    min_time_ms = min(times_ms),
    max_time_ms = max(times_ms),
    n_samples = length(times_ms),
    n_timepoints = length(times)
  ))
}

# Function to verify solution accuracy
verify_solution_accuracy <- function(scenario, solver) {
  times <- seq(scenario$tspan[1], scenario$tspan[2], by = scenario$timestep)
  
  # Solve ODE
  sol <- ode(y = C0, times = times, func = simple_kinetic_model, 
             parms = list(k_el = k_el), method = solver$method,
             rtol = 1e-8, atol = 1e-10)
  
  # Analytical solution: C(t) = C0 * exp(-k_el * t)
  final_time <- scenario$tspan[2]
  analytical_final <- C0[1] * exp(-k_el * final_time)
  numerical_final <- tail(sol[, "C"], 1)
  relative_error <- abs(numerical_final - analytical_final) / analytical_final
  
  return(list(
    scenario = scenario$name,
    solver = solver$name,
    analytical_final = analytical_final,
    numerical_final = numerical_final,
    relative_error = relative_error,
    n_timepoints = nrow(sol)
  ))
}

# Main benchmarking function
run_all_benchmarks <- function() {
  cat("Starting R benchmarking for simple kinetic model...\n")
  cat("Model: dC/dt = -k_el * C, where k_el =", k_el, "/day, C0 =", C0[1], "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  results <- list()
  accuracy_results <- list()
  
  # Run benchmarks for all scenario-solver combinations
  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    cat("\n--- ", scenario$name, ": ", scenario$description, " ---\n", sep = "")
    
    for (j in seq_along(solvers)) {
      solver <- solvers[[j]]
      
      tryCatch({
        # Benchmark performance
        result <- run_benchmark(scenario, solver)
        results[[length(results) + 1]] <- result
        
        cat("  ", solver$name, ": ", round(result$mean_time_ms, 3), " ms (mean)\n", sep = "")
        
        # Test accuracy (for all solvers)
        accuracy <- verify_solution_accuracy(scenario, solver)
        accuracy_results[[length(accuracy_results) + 1]] <- accuracy
        cat("  ", solver$name, " Accuracy: ", round(accuracy$relative_error * 100, 6), "% relative error\n", sep = "")
        
      }, error = function(e) {
        cat("  ERROR with ", solver$name, ": ", e$message, "\n", sep = "")
        results[[length(results) + 1]] <- list(
          scenario = scenario$name,
          solver = solver$name,
          description = scenario$description,
          error = e$message
        )
      })
    }
  }
  
  # Convert lists to data frames
  results_df <- do.call(rbind, lapply(results, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  accuracy_df <- do.call(rbind, lapply(accuracy_results, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  
  # Save results
  write.csv(results_df, "r_benchmark_results.csv", row.names = FALSE)
  write.csv(accuracy_df, "r_accuracy_results.csv", row.names = FALSE)
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("Benchmarking complete!\n")
  cat("Results saved to:\n")
  cat("  - r_benchmark_results.csv\n")
  cat("  - r_accuracy_results.csv\n")
  
  return(list(results = results_df, accuracy = accuracy_df))
}

# Run the benchmarks
benchmark_results <- run_all_benchmarks()