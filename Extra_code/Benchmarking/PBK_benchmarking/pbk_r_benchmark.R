library(deSolve)
library(microbenchmark)

create.params.Fabrega <- function(user_input){
  with(as.list(user_input),{
    # Physiological parameters (from Brown, et al)
    # fractional blood flows to tissues
    QCC = 12.5		  # Cardiac blood output (L/h/kg^0.75)
    QFC = 0.052		  # Fraction cardiac output going to fat
    QLC = 0.19 #0.25 from Loccisano 2011 		  # Fraction cardiac output going to liver
    QKC = 0.175		  # Fraction cardiac output going to kidney
    #QFilC = 0.035		# Fraction cardiac output to the filtrate compartment (20% of kidney blood flow)
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
    VBC = 0.021	# Fraction brain volume
    
    # Scaling parameters
    QC = QCC*(BW^0.75)*24	#Cardiac output (L/day)
    Htc = 0.44      #hematocrit
    QCP = QC*(1-Htc)	# Plasma flow
    QL = QLC*QCP			# Plasma flow to liver (L/day)
    QF = QFC*QCP			# Plasma flow to fat (L/day)
    QK = QKC*QCP	    # Plasma flow to kidney (L/day)
    QFil = 0.2*QK   	# Plasma flow to filtrate compartment (L/day)# 20% of QK
    QG = QGC*QCP	    # Plasma flow to gut (L/day)
    QLu = QLuC*QCP	  # Plasma flow to lungs (L/day)
    QB = QBC*QCP	    # Plasma flow to brain (L/day)
    QR = QCP - QL - QF - QK - QFil - QG - QLu - QB	 # Plasma flow to rest of the body (L/day)
    
    Qbal = QCP - (QL+QF+QK+QFil+QG+QR+QLu+QB)        # balance check 
    
    VL = VLC*BW			    # Liver volume (L)
    VF = VFC*BW			    # Fat volume (L)
    VK = VKC*BW			    # Kidney volume (L)
    VFil = VFilC*BW	    # Filtrate compartment volume (L)
    VG = VGC*BW			    # Gut volume (L)
    VPlas = VPlasC*BW		# Plasma volume (L)
    VLu = VLuC*BW			  # Lungs volume (L)
    VB = VBC*BW			    # Brain volume (L)
    VR = 1*BW - VL - VF - VK - VFil - VG -VLu -VB - VPlas		# Rest of the body volume (L) # Loccisano 2011 uses 0.84*BW
    
    Vbal = (1*BW)-(VL+VF+VK+VFil+VG+VPlas+VR+VLu+VB)         # Balance check
    
    
    if(substance == 'PFOA'){
      Tm = 147.4*24 #ug/day
      Kt = 0.116 #ug/L
      Free = 0.03 # Fabrega (2014) Table 1
      PL = 1.03 # Partition Coefficient for Liver: Fabrega (2014) Table 1
      PF = 0.467 # Partition Coefficient for Fat: Fabrega (2014) Table 1
      PB = 0.17 # Partition Coefficient for Brain: Fabrega (2014) Table 1
      PLu = 1.27 # Partition Coefficient for Lungs: Fabrega (2014) Table 1
      PK = 1.17 # Partition Coefficient for Kidney: Fabrega (2014) Table 1
      PG = 0.05 # Partition Coefficient for Gut: Loccisano (2011) Table 1
      PR = 0.12 # Partition Coefficient for Rest of body: Loccisano (2011) Table 1
      kurinec = 3e-04	#urinary elimination rate constant  (/h/kg^-0.25); 
      kurine = kurinec*(BW^-0.25)*24 # Elimination rate (1/day)
    }else if(substance == 'PFOS'){
      Tm = 24*86 #ug/day
      Kt = 0.0176 #ug/L
      Free = 0.03 # Fabrega (2014) Table 1
      PL = 2.67 # Partition Coefficient for Liver: Fabrega (2014) Table 1
      PF = 0.33 # Partition Coefficient for Fat: Fabrega (2014) Table 1
      PB = 0.255 # Partition Coefficient for Brain: Fabrega (2014) Table 1
      PLu = 0.155 # Partition Coefficient for Lungs: Fabrega (2014) Table 1
      PK = 1.26 # Partition Coefficient for Kidney: Fabrega (2014) Table 1
      PG = 0.57 # Partition Coefficient for Gut: Loccisano (2011) Table 1
      PR = 0.2 # Partition Coefficient for Rest of body: Loccisano (2011) Table 1
      kurinec = 1e-3	#urinary elimination rate constant  (/h/kg^-0.25); estimated from Harada, 
      kurine = 24*kurinec*BW^(-0.25) # Elimination rate (1/day)
    }else if(substance == 'PFHxA'){
      Tm = 245.6 # PFHxA Tmax
      Kt = 0.6 # Estimated value for PFHxA
      Free = 0.01 # Fabrega (2014) Table 1
      PL = 0.001 # Partition Coefficient for Liver: Estimated
      PF = 0.467 # Partition Coefficient for Fat: Estimated
      PB = 43.68 # Partition Coefficient for Brain: Estimated
      PLu = 31.92 # Partition Coefficient for Lungs: Estimated
      PK = 11.57 # Partition Coefficient for Kidney: Estimated
      PG = 0.05 # Partition Coefficient for Gut: Estimated
      PR = 0.12 # Partition Coefficient for Rest of body: Estimated
      kurinec = 3e-04	#urinary elimination rate constant (/h/kg^-0.25); estimated from Harada, 
      kurine = kurinec*(BW^-0.25)*24 # Elimination rate (1/day)
    }else if(substance == 'PFBS'){
      Tm = 6.1*24 #ug/day - converted from ug/h
      Kt = 5 #ug/L
      Free = 0.001 # Fabrega (2015) Table 2
      PL = 128.8 # Partition Coefficient for Liver: Fabrega (2015) Table 2
      PF = 0.467 # Partition Coefficient for Fat: Fabrega (2014) Table 1
      PB = 201.6 # Partition Coefficient for Brain: Fabrega (2015) Table 2
      PLu = 56.11 # Partition Coefficient for Lungs: Fabrega (2015) Table 2
      PK = 6.27 # Partition Coefficient for Kidney: Fabrega (2015) Table 2
      PG = 0.05 # Partition Coefficient for Gut: Loccisano (2011) Table 1
      PR = 0.12 # Partition Coefficient for Rest of body: Loccisano (2011) Table 1
      kurinec = 3e-04	#urinary elimination rate constant (/h/kg^-0.25); 
      kurine = kurinec*(BW^-0.25)*24 # Elimination rate (1/day)
    }
    return(list('QC'=QC, 'QCP'=QCP, 'QL'=QL, 'QF'=QF, 'QK'=QK, 
                'QFil'=QFil, 'QG'=QG, 'QLu'=QLu, 'QB'=QB, 'QR'=QR,
                'VPlas'=VPlas, 'VL'=VL, 'VF'=VF, 'VK'=VK, 
                'VFil'=VFil, 'VG'=VG, 'VLu'=VLu, 'VB'=VB, 'VR'=VR,
                'PL'=PL, 'PF'=PF, 'PB'=PB, 'PLu'=PLu, 'PK'=PK, 
                'PG'=PG, 'PR'=PR,
                'Tm'=Tm, 'Kt'=Kt, 'Free'=Free, 'kurine'=kurine, 
                "ingestion" = ingestion, "ingestion_time" = ingestion_time,
                "admin_dose" = admin_dose, "admin_time" = admin_time,
                "admin_type" = admin_type, "exp_type" = exp_type))
  })
}

create.inits.Fabrega <- function(parameters){
  with(as.list(parameters),{
    
    APlas<-0; AG<-0; AL<-0; AF<-0; ALu<-0; AB<-0
    AK<-0; AFil<-0; AStore<-0; AUrine<-0
    AR<-0; ingestion<-0; 
    
    return(c("APlas"=APlas, "AG"=AG, "AL"=AL, "AF"=AF, "ALu"=ALu, "AB"=AB,
             "AK"=AK, "AFil"=AFil, "AStore"=AStore, "AUrine"=AUrine,
             "AR"=AR, 'ingestion'=ingestion
             ))
  })
}


create.events.Fabrega <- function(parameters){
  with(as.list(parameters), {
    if (admin_type == "iv"){
      # Calculate number of administrated doses and corresponding administration time for IV
      ldose <- length(admin_dose)
      ltimes <- length(admin_time)
      # If not equal, then stop 
      if (ltimes != ldose){
        stop("The times of administration should be equal in number to the doses")
      }else{
        events <- list(data = rbind(data.frame(var = c("APlas"),  time = admin_time, 
                                             value = admin_dose*Free, method = c("add")) ))
      }
    }else if(admin_type == "oral"){
      lingest <- length(ingestion)
      lingesttimes <- length(ingestion_time)
      # If not equal, then stop 
      if (lingest != lingesttimes){
        stop("The times of ingestion rate change should be equal to the ingestion time vector")
      }else{
        if(exp_type == "pharmacokinetics"){
          # For pharmacokinetic studies, add the dose directly to the gut compartment
          events <- list(data = data.frame(var = c("AG"),  time = ingestion_time, 
                                           value = ingestion, method = c("add")))
        }else if(exp_type == "biomonitoring"){
          # For continuous exposure studies, set the ingestion rate
          events <- list(data = data.frame(var = c("ingestion"),  time = ingestion_time, 
                                         value = ingestion, method = c("rep")))
        }
      }
    }
    
    return(events)
  })
}

custom.func <- function(){
  return()
}

ode.func.Fabrega <- function(time, inits, params, custom.func){
  with(as.list(c(inits,params)),{
    
    # Units
    # C_i: Concentration in tissue i (ug/L) 
    # ingestion: hourly ingestion of PFASs (ug/h)
    
    # Concentrations (ug/L)
    CPlas_free <- APlas/VPlas # free concentration in plasma 
    CPlas <- CPlas_free/Free # Concentration in plasma
    CG <- AG/VG # Concentration in Gut
    CL <- AL/VL # Concentration in Liver
    CF <- AF/VF # Concentration in Fat
    CLu <- ALu/VLu # Concentration in Lungs
    CB <- AB/VB # Concentration in Brain
    CK <- AK/VK # Concentration in Kidney
    CFil <- AFil/VFil # Concentration in Filtration compartment
    CR <- AR/VR  # Concentration in Rest of the body compartment
    
    # Plasma compartment (ug)
    dAPlas <- QF*Free*CF/PF + (QL+QG)*Free*CL/PL + QR*Free*CR/PR + QK*Free*CK/PK + 
      QB*Free*CB/PB + QLu*Free*CLu/PLu - QCP*Free*CPlas
    
    # Gut compartment
    dAG <- QG*Free*(CPlas - CG/PG) + ingestion
    
    # Liver compartment
    dAL <- QL*Free*CPlas + QG*Free*CG/PG - (QL+QG)*Free*CL/PL
    
    # Fat compartment
    dAF <- QF*Free*(CPlas - CF/PF)
    
    # Lungs Compartment
    dALu <- QLu*Free*(CPlas - CLu/PLu)
    
    # Brain Compartment
    dAB <- QB*Free*(CPlas - CB/PB)
    
    # Kidney compartment
    dAK <- QK*Free*(CPlas - CK/PK) + (Tm*CFil)/(Kt+CFil)
    
    # Filtrate compartment
    dAFil = QFil*(Free*CPlas - CFil) - (Tm*CFil)/(Kt+CFil)
    
    # Storage compartment for urine
    dAStore <- QFil*CFil - kurine*AStore
    
    # Amount excrete via urine
    dAUrine <- kurine*AStore
    
    # Rest of the body
    dAR <- QR*Free*(CPlas - CR/PR)
    
    # PFAS concentration
    dingestion <- 0
    
    
    return(list(c("dAPlas"=dAPlas, "dAG"=dAG, "dAL"=dAL, "dAF"=dAF, "dALu"=dALu, 
                  "dAB"=dAB, "dAK"=dAK, "dAFil"=dAFil, "dAStore"=dAStore, "dAUrine"=dAUrine,
                  "dAR"=dAR, 'dingestion'=dingestion),
                "CPlas"=CPlas, "CG"=CG, "CL"=CL, "CF"=CF, "CLu"=CLu, "CB"=CB,
                "CK"=CK, "CFil"=CFil, "CR"=CR,   "Cserum" = CPlas, "Cliver" = CL) ) 
  })
}

# Benchmark configuration
create_scenario <- function(name, substance, tspan, timestep, ingestion, ingestion_time, exp_type, description) {
  list(
    name = name,
    substance = substance,
    tspan = tspan,
    timestep = timestep,
    ingestion = ingestion,
    ingestion_time = ingestion_time,
    exp_type = exp_type,
    description = description
  )
}

# Define the four scenarios
scenarios <- list(
  create_scenario("Scenario1", "PFHxA", c(0, 10), 1.0, c(100.0), c(0.01), "biomonitoring", "10 days, timestep 1 day, PFHxA"),
  create_scenario("Scenario2", "PFHxA", c(0, 10), 0.1, c(100.0), c(0.01), "biomonitoring", "10 days, timestep 0.1 day, PFHxA"),
  create_scenario("Scenario3", "PFHxA", c(0, 5*365), 1.0, c(100.0), c(0.01), "biomonitoring", "5 years, timestep 1 day, PFHxA"),
  create_scenario("Scenario4", "PFBS", c(0, 365), 1.0, rep(700.0 * 70, length(seq(0, 365, by=30))), seq(0, 365, by=30), "pharmacokinetics", "1 year, monthly dosing, PFBS")
)

# R solvers to test
solvers <- list(
  list(method = "lsoda", name = "lsoda"),
  list(method = "lsodes", name = "lsodes"),
  list(method = "bdf", name = "bdf")
)

# Model parameters
BW <- 70.0  # Body weight in kg
f_unabs <- 0.5

create_user_input <- function(scenario) {
  list(
    BW = BW,
    substance = scenario$substance,
    f_unabs = f_unabs,
    admin_dose = 0,
    admin_time = 0,
    ingestion = scenario$ingestion,
    ingestion_time = scenario$ingestion_time,
    admin_type = "oral",
    exp_type = scenario$exp_type
  )
}

# Function to run benchmark for a single scenario and solver
run_benchmark <- function(scenario, solver) {
  cat("Running", scenario$name, "with", solver$name, "...\n")
  
  # Create time sequence
  times <- seq(scenario$tspan[1], scenario$tspan[2], by = scenario$timestep)
  
  # Set up user input and parameters
  user_input <- create_user_input(scenario)
  params <- create.params.Fabrega(user_input)
  inits <- create.inits.Fabrega(params)
  events <- create.events.Fabrega(params)
  
  # Benchmark the solver
  benchmark_result <- microbenchmark(
    result = ode(y = inits, times = times, func = ode.func.Fabrega, 
                 parms = params, events = events, method = solver$method,
                 rtol = 1e-8, atol = 1e-10),
    times = 3,
    unit = "ms"
  )
  
  # Extract timing statistics
  times_ms <- benchmark_result$time / 1e6  # Convert to milliseconds
  
  return(list(
    scenario = scenario$name,
    substance = scenario$substance,
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
  
  user_input <- create_user_input(scenario)
  params <- create.params.Fabrega(user_input)
  inits <- create.inits.Fabrega(params)
  events <- create.events.Fabrega(params)
  
  # Solve ODE
  sol <- ode(y = inits, times = times, func = ode.func.Fabrega, 
             parms = params, events = events, method = solver$method,
             rtol = 1e-8, atol = 1e-10)
  
  # Extract final plasma concentration
  final_cplas <- tail(sol[, "CPlas"], 1)
  
  return(list(
    scenario = scenario$name,
    substance = scenario$substance,
    solver = solver$name,
    final_plasma_conc = final_cplas,
    n_timepoints = nrow(sol),
    final_time = max(times)
  ))
}

# Main benchmarking function
run_all_benchmarks <- function() {
  cat("Starting R PBK benchmarking...\n")
  cat("Model: 12-compartment PFAS PBK model\n")
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
        
        # Test accuracy (only once per scenario, using lsoda)
        if (solver$name == "lsoda") {
          accuracy <- verify_solution_accuracy(scenario, solver)
          accuracy_results[[length(accuracy_results) + 1]] <- accuracy
          cat("  Final plasma conc: ", round(accuracy$final_plasma_conc, 6), " ng/mL\n", sep = "")
        }
        
      }, error = function(e) {
        cat("  ERROR with ", solver$name, ": ", e$message, "\n", sep = "")
        results[[length(results) + 1]] <- list(
          scenario = scenario$name,
          substance = scenario$substance,
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
  write.csv(results_df, "r_pbk_benchmark_results.csv", row.names = FALSE)
  write.csv(accuracy_df, "r_pbk_accuracy_results.csv", row.names = FALSE)
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("PBK Benchmarking complete!\n")
  cat("Results saved to:\n")
  cat("  - r_pbk_benchmark_results.csv\n")
  cat("  - r_pbk_accuracy_results.csv\n")
  
  return(list(results = results_df, accuracy = accuracy_df))
}

# Run the benchmarks
benchmark_results <- run_all_benchmarks()