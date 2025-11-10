# Acetone PBK model after Mork & Johanson (2006)
# ---------------------------------------------------------------
# This script implements the human physiologically-based kinetic

setwd("Mork and Johanson,2006")
  library(deSolve)
  library(ggplot2)
  library(jaqpotr)

#=========================
#1. Parameters of the model
#=========================

create.params <- function(user.input) {
  with(as.list(user.input), {

#Table 1 independent parameters
 #Partition coefficients
    
    P_blood_air   = 196    # blood:air
    P_rich_air    = 137      # richly perfused tissue:air
    P_muscle_air  = 151   # muscle:air
    P_fat_air     = 86       # fat:air
    P_water_air   = 263   # water (mucosa):air
    P_liver_air   = 113     # liver:air

 #Volumes for respiratory tract

    V_ds  = 0.15   # anatomical dead space (L)
    V_bro = 0.10   # bronchioles air space (L)
    V_muc = 0.005  # mucosa (L)
    V_alv = 4.10   # alveolar air space (L)
   
    
#Metabolic parameters

    Vmax = 18.6/60*BW^0.75      # mg / (kg^0.75 · min)
    Km   = 48.4      # mg / L

    #Ventilation and flows (L/min)
    
    if (workload == 0) {
      Qp <- 9.1; VPR <- 1.4; RF <- 15.4;
      fQfat <- 0.03; fQrich <- 0.44; fQrest <- 0.105;
      fQwork <- 0.105; fQliver <- 0.32
    } else if (workload == 50) {
      Qp <- 23.8; VPR  <- 2.3; RF <- 17.1;
      fQfat <- 0.05; fQrich <- 0.29;
      fQrest <- 0.050; fQwork <- 0.450; fQliver <- 0.16
    } else if (workload == 100) {
      Qp  <- 35.3; VPR  <- 2.5; RF  <- 19.9;
      fQfat <- 0.05; fQrich <- 0.21;
      fQrest <- 0.040; fQwork <- 0.610; fQliver <- 0.09
    } else if (workload == 150) {
      Qp  <- 50.1; VPR  <- 2.5; RF  <- 24.0;
      fQfat <- 0.03; fQrich <- 0.19;
      fQrest  <- 0.03; fQwork <- 0.70; fQliver <- 0.05
    }
    
    #(L/min)
    
    Q_alv  <- Qp  - (V_ds * RF ) 
    Qc  <- Q_alv  / VPR
    Q_fat  <- Qc  * fQfat 
    Q_rich  <- Qc  * fQrich
    Q_rest  <- Qc  * fQrest
    Q_work  <- Qc  * fQwork
    Q_liver  <- Qc  * fQliver
    
#Table 3   
#Body composition  
  
  TBW <- -12.86 + 0.1757 * height + 0.331 * BW #total body water (males/cm,kg)
  FFM <- TBW / 0.72 #fat free mass (kg)
  LBV <- FFM / 1.1 #lean body volume (L)
  
#Compartment volumes (L)
  
  V_fat = (BW - FFM)/0.92                       #fat
  V_art_WW = (0.01933 + 0.00907) * LBV          #Washin–washout model
  V_art_IT = (0.0256 + 0.01933 + 0.00907) * LBV #Inert-tube model
  V_rich = (0.00532 + 0.0103) * LBV            #Richly perfused tissue
  V_rest = 0.344 * LBV                   #Resting muscle and skin
  V_work = 0.344 * LBV                #working muscle
  V_liver = 0.0285 * LBV                        #liver

  #===============================
  # Steady-state helpers (Eqs. 2,5,6 & Eq. 3)
  #===============================
  
  compute_ss_washin <- function(params, Cendo) {
    Pb_a <- params$P_blood_air
    Pw_a <- params$P_water_air
    Pliv_a <- params$P_liver_air
    Qalv <- params$Q_alv   # L/min
    Qc   <- params$Qc      # L/min
    
    A <- Cendo / Pb_a
    r <- (Qc/Qalv) * Pw_a
    
    Cbro_ss <- (r / (2 + 3*r)) * A
    Calv_ss <- ( (2 + 4*r) / (2 + 3*r) ) * A
    Cmuc_ss <- (3*r / (2 + 3*r)) * Pw_a * A
    Cliv_ss <- (Pliv_a/Pb_a) * Cendo  # well-stirred equilibrium in liver
    
    list(C_bro_SS = Cbro_ss, C_alv_SS = Calv_ss,
         C_muc_SS = Cmuc_ss, C_liv_SS = Cliv_ss)
  }
  
  compute_PR <- function(params, Cendo, model_type = c("washin–washout","inert-tube")) {
    model_type <- match.arg(model_type)
    
    Pb_a <- params$P_blood_air
    Pw_a <- params$P_water_air
    Pliv_a <- params$P_liver_air
    Qalv <- params$Q_alv   # L/min
    Qc   <- params$Qc      # L/min
    Vmax <- params$Vmax    # mg/(kg^0.75·min)
    Km   <- params$Km      # mg/L
    
    if (model_type == "washin–washout") {
      ss <- compute_ss_washin(params, Cendo)
      term_exh_in   <- Qalv * (Cendo / Pb_a)
      term_metab    <- (Vmax * ss$C_liv_SS * Pb_a/Pliv_a) / (Km + ss$C_liv_SS * Pb_a/Pliv_a)
      term_muc_xchg <- Qc   * (Cendo*Pb_a/Pw_a - ss$C_muc_SS)
      term_exh_out  <- Qalv * ss$C_alv_SS
      PR <- term_exh_in + term_metab + term_muc_xchg - term_exh_out  # mg/min
      
      c(PR = PR, ss)
      
    } else { # inert-tube: omit mucosa/alveoli exchange terms

      C_liv_SS <- (Pliv_a/Pb_a) * Cendo
      term_exh_in <- Qalv * (Cendo / Pb_a)
      term_metab  <- (Vmax * C_liv_SS * Pb_a/Pliv_a) / (Km + C_liv_SS * Pb_a/Pliv_a)
      PR <- term_exh_in + term_metab  # mg/min
      
      c(PR = PR,
        C_bro_SS = NA_real_, C_alv_SS = NA_real_, C_muc_SS = NA_real_, C_liv_SS = C_liv_SS)
    }
  }

  #-------------------------------
  # Fixed endogenous arterial acetone
  #-------------------------------
  
  
  core <- list(P_blood_air = P_blood_air, P_water_air = P_water_air, P_liver_air = P_liver_air,
               Q_alv = Q_alv, Qc = Qc, Vmax = Vmax, Km = Km)
  
  # model_type should be provided in user.input (e.g., "inert-tube" or "washin–washout")
  ss_pr <- compute_PR(core, Cendo, model_type)
  
  
  
  return(list('BW' = BW, 'height' = height,
    'P_blood_air' = P_blood_air, 'P_rich_air' = P_rich_air, 'P_muscle_air' = P_muscle_air,
    'P_fat_air' = P_fat_air, 'P_water_air' = P_water_air, 'P_liver_air' = P_liver_air,

    'V_ds' = V_ds, 'V_bro' = V_bro, 'V_muc' = V_muc, 'V_alv' = V_alv,
    'Vmax' = Vmax, 'Km' = Km,

    'workload' = workload, 'Qp' = Qp, 'VPR' = VPR, 'RF' = RF,
    'fQfat' = fQfat, 'fQrich' = fQrich, 'fQrest' = fQrest,
    'fQwork' = fQwork, 'fQliver' = fQliver,
    'Q_alv' = Q_alv, 'Qc' = Qc, 'Q_fat' = Q_fat, 'Q_rich' = Q_rich,
    'Q_rest' = Q_rest, 'Q_work' = Q_work, 'Q_liver' = Q_liver,
    
    'TBW' = TBW, 'FFM' = FFM, 'LBV' = LBV,
    'V_fat' = V_fat, 'V_art_WW' = V_art_WW, 'V_art_IT' = V_art_IT,
    'V_rich' = V_rich, 'V_rest' = V_rest,
    'V_work' = V_work, 'V_liver' = V_liver,
    
    'Cendo' = Cendo, 'PR' = unname(ss_pr['PR']),
    'C_bro_SS' = unname(ss_pr['C_bro_SS']),
    'C_muc_SS' = unname(ss_pr['C_muc_SS']),
    'C_alv_SS' = unname(ss_pr['C_alv_SS']),
    'C_liv_SS' = unname(ss_pr['C_liv_SS']),
    
    'inhalation'=inhalation, 'inhalation_time'=inhalation_time, 'model_type'=model_type

    
  ))
  })
}

#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    A_fat <-0 ; A_rich <-0 ; A_rest <-0 ; A_work <-0;
    A_liv <-0 ; A_bro <-0 ; A_muc <-0 ; A_alv <-0;
    A_art <-0 ; A_met <-0 ; A_exh <-0; PR <-0; inhalation <- 0
    
    
    return(c(  'A_fat' =  A_fat, 'A_rich' =  A_rich, 'A_rest' =  A_rest, 'A_work' =  A_work,
               'A_liv' =  A_liv, 'A_bro' =  A_bro, 'A_muc' =  A_muc, 'A_alv' =  A_alv,
               'A_art' =  A_art,'A_met' =  A_met, 'A_exh' =  A_exh, 'PR'=PR, 
               'inhalation'=inhalation
               
               
    ))
  })
}


#================================================================================================
#===================
#3. Events function
#===================

create.events <- function(parameters){
  with(as.list(parameters), {
    ldose <- length(inhalation)
    ltimes <- length(inhalation_time)
    # If not equal, then stop 
    if (ltimes != ldose){
      stop("The times of administration should be equal in number to the doses")
    }else{
      if (model_type == "washin–washout"){
      events <- list(data = data.frame(var = c("inhalation"),  time = inhalation_time, 
                                       value = inhalation, method = c("rep")))
    }else if (model_type == "inert-tube"){
      events <- list(data = data.frame(var = c("inhalation"),  time = inhalation_time, 
                                       value = inhalation, method = c("rep")))
      }
    
    }
    return(events)
  })
}


#================================================================================================
#==================
#4. Custom function 
#==================

custom.func <- function(){
  return()
}

# ============================================================================
#5. ODEs System
# ============================================================================

ode.func <- function(time, inits, params, custom.func){
  with(as.list(c(inits, params)),{

# Calculate tissue concentrations (mg/L)
#========================================================================================
   
    C_bro  <- A_bro  / V_bro
    C_muc  <- A_muc  / V_muc
    C_alv  <- A_alv  / V_alv
    C_liv  <- A_liv  / V_liver
    C_rich <- A_rich / V_rich
    C_rest <- A_rest / V_rest #rest muscle and skin
    C_work <- A_work/ V_work
    C_fat  <- A_fat  / V_fat
    
    if (model_type == "washin–washout") { 
      C_art  <- A_art  / V_art_WW
      
    } else if (model_type == "inert-tube")
      C_art  <- A_art  / V_art_IT
      

#ODEs
#======================================================================================== 

    #non-metabolizing compartments
    #fat, richly perfused tissue, resting muscle and skin and working muscle
    
    dA_fat  <- Q_fat*C_art - Q_fat*C_fat*P_blood_air/P_fat_air
    dA_rich <- Q_rich*C_art - Q_rich*C_rich*P_blood_air/P_rich_air
    dA_rest <- Q_rest*C_art - Q_rest*C_rest*P_blood_air/P_muscle_air
    dA_work <- Q_work*C_art - Q_work*C_work*P_blood_air/P_muscle_air
    
    #metabolizing compartment (liver)
    
    metab_rate <- (Vmax * C_liv * P_blood_air / P_liver_air) /
                  (Km + C_liv * P_blood_air / P_liver_air)
   
    
    dA_liv  <- Q_liver*C_art - Q_liver*C_liv*P_blood_air/P_liver_air + PR - metab_rate
    
   #respiratory tract (washin–washout model)
    
    if (model_type == "washin–washout") {  
    dA_bro <- inhalation - 2*Q_alv*C_bro - Q_alv*C_bro + Q_alv*C_muc/P_water_air
    dA_muc <- Q_alv*C_bro - Q_alv*C_muc/P_water_air + Qc*C_art*P_water_air/P_blood_air - Qc*C_muc
    dA_alv <- Q_alv*C_bro + Q_alv*C_art/P_blood_air - Q_alv*C_alv
   
    tissue_out <- Q_liver*C_liv*P_blood_air/P_liver_air + Q_rich*C_rich*P_blood_air/P_rich_air +
                  Q_rest*C_rest*P_blood_air/P_muscle_air + Q_work*C_work*P_blood_air/P_muscle_air +
                  Q_fat*C_fat*P_blood_air/P_fat_air
    
    dA_art <- Q_alv*C_alv - Q_alv*C_art/P_blood_air + Qc*C_muc - Qc*C_art*P_water_air/P_blood_air +
              tissue_out - Q_fat*C_art - Q_rich*C_art - Q_rest*C_art -
              Q_work*C_art - Q_liver*C_art
    
    dA_exh <-  Q_alv*C_bro
    
#I deleted Qc*C_art leaving from arterial blood CMT
    
    #respiratory tract (inert-tube model)
    
    } else if (model_type == "inert-tube") {
    dA_bro <-0
    dA_muc <-0
    dA_alv <-0
    
    tissue_out <- Q_liver*C_liv*P_blood_air/P_liver_air + Q_rich*C_rich*P_blood_air/P_rich_air +
                  Q_rest*C_rest*P_blood_air/P_muscle_air + Q_work*C_work*P_blood_air/P_muscle_air +
                  Q_fat*C_fat*P_blood_air/P_fat_air
    
    dA_art <- inhalation + tissue_out - Q_alv*C_art/P_blood_air - Q_fat*C_art - Q_rich*C_art - 
                Q_rest*C_art - Q_work*C_art - Q_liver*C_art
  
    dA_exh <- Q_alv*C_art/P_blood_air
    
    }
    
    
    dA_met <- metab_rate
    dPR <- - PR
    
    dinhalation = 0

    list(c( 'dA_fat' = dA_fat, 'dA_rich' = dA_rich, 'dA_rest' = dA_rest, 'dA_work' = dA_work,
            'dA_liv' = dA_liv, 'dA_bro' = dA_bro, 'dA_muc' = dA_muc, 'dA_alv' = dA_alv,
            'dA_art' = dA_art,'dA_met' = dA_met, 'dA_exh' = dA_exh,
            'dPR'=dPR, 'dinhalation' = dinhalation),
        
           'C_fat' = C_fat, 'C_rich' = C_rich, 'C_rest' = C_rest, 'C_work' = C_work,
           'C_liv' = C_liv, 'C_bro' = C_bro, 'C_muc' = C_muc, 'C_alv' = C_alv,
           'C_art' = C_art)
  })
}


#=============
#6. User input 
#=============

################################################################################

setwd("Mork and Johanson,2006/Data")


# Read data
Acetone_arterial_blood_IT_model_1 <- openxlsx::read.xlsx("Acetone_arterial blood_fit_Inert_tube_model_1.xlsx")
Acetone_arterial_blood_WW_model_1 <- openxlsx::read.xlsx("Acetone_arterial blood_fit_Washin-washout_model_1.xlsx")


setwd("Mork and Johanson,2006")

dataset <- list("df1" = Acetone_arterial_blood_IT_model_1, "df2" = Acetone_arterial_blood_WW_model_1)
                

#################################
#--------------------------------
# Wigaeus et al., 1981 - inert-tube model - series 1
#--------------------------------
#################################

# Set up parameters for event-based exposure
BW <- 71 #kg
height <- 181 #cm

workload <- 0
model_type <- "inert-tube"
Cendo <- 1.4  # mg/L (fixed value)
exposure <- 1.309 #mg/L
duration <- 120 #min
inhalation <- 1
inhalation_time <-1

sex <- "M"

initial_input <- list('BW'=BW,
                   'height' = height,
                   'Cendo'=Cendo,
                   'inhalation'=inhalation,
                   'inhalation_time'=inhalation_time,
                   "workload"=workload,
                   "model_type"=model_type,
                   "sex" = sex)

params_i <- create.params(initial_input)

total_admin.dose <- exposure* params_i$Q_alv* duration
n_doses <- floor(duration/1) + 1 
interval_min <- 1
eps_min <- 0.01
inj_times_min <- seq(0, by = interval_min, length.out = n_doses)

dose_per_event <- total_admin.dose / n_doses
dose_rate <- dose_per_event/interval_min

times_raw  <- as.vector(rbind(inj_times_min, inj_times_min + interval_min - eps_min))
values_raw <- as.vector(rbind(rep(dose_per_event, n_doses), rep(0, n_doses)))

# Round to match your printed format
inhalation_time <- round(times_raw, 2)
inhalation      <- round(values_raw, 2)

user_input <- list('BW'=BW,
                   'height' = height,
                   'Cendo'=Cendo,
                   "workload"=workload,
                   "model_type"=model_type,
                   "inhalation"=inhalation,
                   "inhalation_time" = inhalation_time,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 240, 0.1)  #min

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))

#======================================df1=========================================================

exp_data <- dataset$df1 
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
column_names <- c("C_art")
preds_Acetone_arterial_blood_IT_model_1 <- list()

for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_Acetone_arterial_blood_IT_model_1 [[i]] <- solution[solution$time %in% exp_time, column_names[i]]} #mg/L

obs_Acetone_arterial_blood_IT_model_1 <- list(exp_data[exp_data$Tissue == "Arterial_Blood", "concentration"]) #mg/L


#################################
#--------------------------------
# Wigaeus et al., 1981 - Washin-washout model - series 1
#--------------------------------
#################################

# Set up parameters for event-based exposure
BW <- 71 #kg
height <- 181 #cm

workload <- 0
model_type <- "inert-tube"
Cendo <- 1.4  # mg/L (fixed value)
exposure <- 1.309 #mg/L
duration <- 120 #min
inhalation <- 1
inhalation_time <-1

sex <- "M"

initial_input <- list('BW'=BW,
                      'height' = height,
                      'Cendo'=Cendo,
                      'inhalation'=inhalation,
                      'inhalation_time'=inhalation_time,
                      "workload"=workload,
                      "model_type"=model_type,
                      "sex" = sex)

params_i <- create.params(initial_input)

total_admin.dose <- exposure* params_i$Q_alv* duration
n_doses <- floor(duration/1) + 1 
interval_min <- 1
eps_min <- 0.01
inj_times_min <- seq(0, by = interval_min, length.out = n_doses)

dose_per_event <- total_admin.dose / n_doses
dose_rate <- dose_per_event/interval_min

times_raw  <- as.vector(rbind(inj_times_min, inj_times_min + interval_min - eps_min))
values_raw <- as.vector(rbind(rep(dose_per_event, n_doses), rep(0, n_doses)))

# Round to match your printed format
inhalation_time <- round(times_raw, 2)
inhalation      <- round(values_raw, 2)

user_input <- list('BW'=BW,
                   'height' = height,
                   'Cendo'=Cendo,
                   "workload"=workload,
                   "model_type"=model_type,
                   "inhalation"=inhalation,
                   "inhalation_time" = inhalation_time,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 240, 0.1)  #min

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))

#======================================df2=========================================================

exp_data <- dataset$df2 
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
column_names <- c("C_art")
preds_Acetone_arterial_blood_WW_model_1 <- list()

for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_Acetone_arterial_blood_WW_model_1 [[i]] <- solution[solution$time %in% exp_time, column_names[i]]} #mg/L

obs_Acetone_arterial_blood_WW_model_1 <- list(exp_data[exp_data$Tissue == "Arterial_Blood", "concentration"]) #mg/L


#=====================================================================================================================================
#################################
#For simulations and plots

#################################
#--------------------------------
# Wigaeus et al., 1981 - inert-tube model - series 1
#--------------------------------
#################################

# Set up parameters for event-based exposure
BW <- 71 #kg
height <- 181 #cm

workload <- 0
model_type <- "inert-tube"
Cendo <- 1.4  # mg/L (fixed value)
exposure <- 1.309 #mg/L
duration <- 120 #min
inhalation <- 1
inhalation_time <-1

sex <- "M"

initial_input <- list('BW'=BW,
                      'height' = height,
                      'Cendo'=Cendo,
                      'inhalation'=inhalation,
                      'inhalation_time'=inhalation_time,
                      "workload"=workload,
                      "model_type"=model_type,
                      "sex" = sex)

params_i <- create.params(initial_input)

total_admin.dose <- exposure* params_i$Q_alv* duration
n_doses <- floor(duration/1) + 1 
interval_min <- 1
eps_min <- 0.01
inj_times_min <- seq(0, by = interval_min, length.out = n_doses)

dose_per_event <- total_admin.dose / n_doses
dose_rate <- dose_per_event/interval_min

times_raw  <- as.vector(rbind(inj_times_min, inj_times_min + interval_min - eps_min))
values_raw <- as.vector(rbind(rep(dose_per_event, n_doses), rep(0, n_doses)))

# Round to match your printed format
inhalation_time <- round(times_raw, 2)
inhalation      <- round(values_raw, 2)

user_input <- list('BW'=BW,
                   'height' = height,
                   'Cendo'=Cendo,
                   "workload"=workload,
                   "model_type"=model_type,
                   "inhalation"=inhalation,
                   "inhalation_time" = inhalation_time,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 240, 0.1)  #min

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))

preds_Acetone_arterial_blood_IT_model_1 <-  solution[, c("time","C_art")]

#################################
#--------------------------------
# Wigaeus et al., 1981 - Washin-washout model - series 1
#--------------------------------
#################################

# Set up parameters for event-based exposure
BW <- 71 #kg
height <- 181 #cm

workload <- 0
model_type <- "inert-tube"
Cendo <- 1.4  # mg/L (fixed value)
exposure <- 1.309 #mg/L
duration <- 120 #min
inhalation <- 1
inhalation_time <-1

sex <- "M"

initial_input <- list('BW'=BW,
                      'height' = height,
                      'Cendo'=Cendo,
                      'inhalation'=inhalation,
                      'inhalation_time'=inhalation_time,
                      "workload"=workload,
                      "model_type"=model_type,
                      "sex" = sex)

params_i <- create.params(initial_input)

total_admin.dose <- exposure* params_i$Q_alv* duration
n_doses <- floor(duration/1) + 1 
interval_min <- 1
eps_min <- 0.01
inj_times_min <- seq(0, by = interval_min, length.out = n_doses)

dose_per_event <- total_admin.dose / n_doses
dose_rate <- dose_per_event/interval_min

times_raw  <- as.vector(rbind(inj_times_min, inj_times_min + interval_min - eps_min))
values_raw <- as.vector(rbind(rep(dose_per_event, n_doses), rep(0, n_doses)))

# Round to match your printed format
inhalation_time <- round(times_raw, 2)
inhalation      <- round(values_raw, 2)

user_input <- list('BW'=BW,
                   'height' = height,
                   'Cendo'=Cendo,
                   "workload"=workload,
                   "model_type"=model_type,
                   "inhalation"=inhalation,
                   "inhalation_time" = inhalation_time,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 240, 0.1)  #min

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))

preds_Acetone_arterial_blood_WW_model_1 <-  solution[, c("time","C_art")]


# ######################################################################################
#Plot the predictions against the observations
library(ggplot2) 

# Function that creates a plot given a compartment name and the respective predictions and observations
create.plots <- function(predictions, observations, compartment){  
  #Colours of observations and predictions
  cls <-  c("predictions" = "#56B4E9",  "Observations" = "#D55E00")
  
  ggplot(data = predictions)+
    geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                          color = '"predictions"'),  size=1.5,alpha = 0.7) +
    geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                             color='"Observations"'), size=4)+
    labs(title = rlang::expr(!!compartment), 
         y = expression("Acetone concentration (mg/L)" ),
         x = "Time (min)")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual("", values=cls,
                       guide = guide_legend(override.aes =
                                              list(shape = c(16,NA),
                                                   linetype = c(0,1))))+
    theme_light() + 
    theme(legend.position=c(1,1), 
          legend.justification=c(0, 1), 
          legend.key.size = unit(1.5, 'cm'),  
          legend.title = element_text(size=14),
          axis.title=element_text(size=14),
          legend.text = element_text(size=14)
    )
  
}


# Convert TDI_urine_data from long to wide format using reshape
experiment1 <- reshape(Acetone_arterial_blood_IT_model_1[c("Tissue" ,"Time_(min)", 
                                        "Concentration_(mg/L)")], 
                       idvar = "Time_(min)", timevar = "Tissue", direction = "wide")
colnames(experiment1) <- c("Time",unique(Acetone_arterial_blood_IT_model_1$Tissue))

experiment2 <- reshape(Acetone_arterial_blood_WW_model_1[c("Tissue" ,"Time_(min)", 
                                                         "Concentration_(mg/L)")], 
                       idvar = "Time_(min)", timevar = "Tissue", direction = "wide")
colnames(experiment2) <- c("Time",unique(Acetone_arterial_blood_WW_model_1$Tissue))

# Put the experiments in a list
experiments <- list(experiment1 = experiment1, experiment2 = experiment2)

# Rename predictions so that they share the same name as the names of the experimental data 
colnames(preds_Acetone_arterial_blood_IT_model_1) <- c( "Time", "Arterial_Blood")
colnames(preds_Acetone_arterial_blood_WW_model_1) <- c( "Time", "Arterial_Blood")

# Create a list containing the corresponding predictions
simulations <- list(predictions1 = preds_Acetone_arterial_blood_IT_model_1,
                    predictions2 = preds_Acetone_arterial_blood_WW_model_1)


# Iterate over all existing experiments and create the accompanying plots
for(i in 1:length(experiments)){
  # Retrieve the corresponding observations and simulations
  observations <- experiments[[i]]
  predictions <- simulations[[i]]
  # Extract the compartment names
  compartments <- names(predictions)[2:length(predictions)]
  
  # Use lapply to iterate over the column names and create plots
  plots <- lapply(compartments, function(compartment) {
    create.plots(predictions, observations, compartment )
  })
  if(length(compartments) == 1){
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 1, nrow = 1,
                                               common.legend = TRUE, legend = "right"))
    
  }else{
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 3, nrow = ceiling(length(plots) / 3),
                                               common.legend = TRUE, legend = "right"))
  }
  
  
  plot.margin=unit(c(0,0,0,0), "pt")
  
  
  # Save the plot with dynamically adjusted dimensions
  ggsave(paste0("experiment_jaq_", i,".png"), plot = final_plot,
         device = 'png', dpi = 300,
         width = 13,
         height = 10,
         units = "in")
}

#save.image("Acetone_PBK_model.RData")

# ====================
# Upload on Jaqpot
# ===================


# Subset of features to be displayed on the user interface
predicted.feats <- c("A_fat" , "A_rich", "A_rest" , "A_work" ,
                     "A_liv" , "A_bro" , "A_muc", "A_alv",
                     "A_art" ,"A_met" , "PR", "A_exh" ,

                     "C_fat" , "C_rich" , "C_rest", "C_work",
                     "C_liv" , "C_bro" , "C_muc", "C_alv",
                     "C_art")



setwd("Mork and Johanson,2006")
jaqpotr::deploy.pbpk(user.input = user_input,out.vars = predicted.feats,
                     create.params = create.params,  create.inits = create.inits,
                     create.events = create.events, custom.func = custom.func,
                     envFile = "/Users/eviepapakyriakopoulou/Documents/Documents/NTUA Post-Doc/INSIGHT/jaqpotkeys/keys.env")
