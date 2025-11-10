# PBPK Model for Ethylene Glycol and Glycolic Acid, Corley et al. (2005)
#
# This model simulates the disposition of ethylene glycol (EG) and its metabolite 
# glycolic acid (GA) in humans after various routes of exposure.
#
# Key Features:
# - Multi-compartment PBPK model
# - Saturable metabolism of EG to GA in liver
# - Saturable metabolism of GA to glyoxylic acid in liver  
# - Kidney model with glomerular filtration and tubular reabsorption for GA
# - Multiple exposure routes (oral, IV, inhalation, dermal, subcutaneous)

library(deSolve)
library(ggplot2)

create.params <- function(user.input) {
  with(as.list(user.input), {
    
    # Physiological parameters for humans (Tables 2 and 4)
    BW_ref <- 60  # kg 
    
    QCC <- 15  # cardiac output, L/h/kg, Scalable by (BW)^0.74.
    QC <- QCC  * BW^0.74 # L/h
    QPC <- 15  # alveolar ventilation, L/h/kg
    Q_alv <- QPC * (BW^0.74)  # Alveolar ventilation L/h
    
    # Tissue volumes as fraction of BW (Table 2 - Human values)
    PVplasma <- 3.5/70 #L for a normal 70 kg adult, https://doi.org/10.1016/B978-0-323-03707-5.50065-6
    Vplasma <- PVplasma * BW #plasma volume kg=L
    PVVen <- 0.022
    VVen <- PVVen*BW 	#volume of venous plasma (L); Corley 2011
    PVArt <- 0.045
    VArt <- PVArt*BW	#volume of arterial plasma (L); Corley 2011
    
    PVBC <- 0.059   # Blood fraction
    VBC <- PVBC * BW   # Blood volume (L)
    PVLC <- 0.0314  # Liver fraction  
    VLC <- PVLC * BW   # Liver volume (L)
  
    if (sex == "M"){
      PVKC <- 0.00443  
      VKC <- PVKC * BW  # Kidney volume (L)
    }else if(sex == "F"){
      PVKC <- 0.00474 
      VKC <- PVKC * BW  # Kidney volume (L)
    }
    
    PVluC <- 0.0115   # Lungs fraction
    VluC <- PVluC * BW # Lung volume (L)
    PVGIC <- 0.034    # GI tract fraction
    VGIC <- PVGIC * BW # GI tract volume (L)
    PVFC <- 0.231     # Fat fraction
    VFC <- PVFC * BW   # Fat volume (L)
    PVSKC <- 0.051    # Skin fraction
    VSKC <- PVSKC * BW # Skin volume (L)
    PVRC <- 0.0371    # Richly perfused fraction
    VRC <- PVRC * BW  # Richly perfused volume (L)
    
    # Slowly perfused fraction (remaining tissues)
    PVSC <- 1.0 - (PVBC + PVLC + PVKC + PVluC + PVGIC + PVFC + PVSKC + PVRC)
    VSC <- PVSC * BW  # Slowly perfused volume (L)
    
    # Kidney tubule volume (10% of kidney volume)
    V_tubule <- 0.1 * VKC
    
    # Blood flows (fraction of cardiac output) - Table 2 Human values
    PQLC <- 0.25   # liver
    QLC <- PQLC * QC  # L/h
    PQGIC <- 0.21  # GI tract  
    QGIC <- PQGIC * QC  # L/h
    PQKC <- 0.25   # Kidney
    QKC <- PQKC * QC  # L/h
    PQFC <- 0.05   # Fat
    QFC <- PQFC * QC  # L/h
    PQSKC <- 0.03  # Skin
    QSKC <- PQSKC * QC  # L/h
    PQSC <- 0.17   # Slowly perfused tissues
    QSC <- PQSC * QC  # L/h
    PQRC <- 0.21   # Richly perfused tissues (from table)
    QRC <- PQRC * QC  # L/h
    
    # Renal physiology parameters (Table 4 - Human values, averaged male/female)
    
    if (sex == "M"){
      PGFR <- 24.19  #L/h/kg  
      GFR <- PGFR * VKC #L/h
    }else if(sex == "F"){
      PGFR <- 27.28  #L/h/kg  
      GFR <- PGFR * VKC #L/h
    }
    
    
    if (sex == "M"){
      QURC <- 0.212  #1/h/kg  
      Q_urine <- QURC * VKC  # L/h
    }else if(sex == "F"){
      QURC <- 0.152  #1/h/kg  
      Q_urine <- QURC * VKC  # L/h
    }
    
    
    # Partition coefficients for EG (Table 3 - Human values)
    PB_EG  <- 17542    # blood:air
    PBS_EG  <- 1.14    # blood:saline
    PSKS_EG <- 1.36    # skin:saline
    PSKA_EG <- 17542   # skin:air
    PL_EG  <- 0.96     # liver:blood
    PK_EG  <- 1.22     # kidney:blood
    PLu_EG  <- 0.96    # lung:blood
    PF_EG  <- 0.64     # fat:blood
    PSK_EG  <- 1.19    # skin:blood
    PGI_EG  <- 1.48    # GI:blood
    PR_EG  <- 0.96     # richly perfused:blood
    PS_EG  <- 0.67     # slowly perfused:blood
    
    # Partition coefficients for GA (Table 3 - Human values)
    PBS_GA  <- 3.36    # blood:saline
    PSKS_GA <- 2.51    # skin:saline
    PL_GA  <- 0.97     # liver:blood
    PK_GA  <- 1.40     # kidney:blood
    PF_GA  <- 1.09     # fat:blood
    PSK_GA  <- 0.75    # skin:blood
    PGI_GA  <- 0.95    # GI:blood
    PR_GA  <- 0.97     # richly perfused:blood
    PS_GA  <- 0.70     # slowly perfused:blood
    
    
    # Absorption rate constants (Table 3)
    KAS <- 1.0        # oral absorption rate (h-1)
    #KSC <- 1.0        # subcutaneous absorption rate (h-1) for rats
    
    # Urinary clearance for EG (Table 3)
    CLC <- 0.06    # l/h/kg
    CL_renal <- CLC * (BW^0.7)  # L/h, scaled
    
    # Metabolism parameters (Table 4 - Human values)
    # EG metabolism: EG -> GA
    KM1 <- 23.8       # mM (Michaelis constant, same for humans and rats)
    KM1_mg <- KM1 * 62.07  # mg/L (convert to mg/L, MW EG = 62.07)
    VMAX1C <- 1300    # mg/h/kg (higher than rats)
    VMAX1 <- VMAX1C * (BW^0.7)  # mg/h, scaled
    
    # GA metabolism: GA -> glyoxylic acid  
    KM2 <- 0.19       # mM (lower than rats)
    KM2_mg <- KM2 * 76.05  # mg/L (convert to mg/L, MW GA = 76.05)
    VMAX2C <- 107.6   # mg/h/kg
    VMAX2 <- VMAX2C * (BW^0.7)  # mg/h, scaled
    
    # Renal reabsorption parameters for GA (Table 4 - Human values)
    KT <- 840         # mg/L (Michaelis constant for tubular reabsorption)
    TMAXEC <- 14      # mg/h/kg (maximum reabsorption capacity)
    TMAX <- TMAXEC * (BW^0.7)  # mg/h, scaled
    
    return(list('BW' = BW, 'QC' = QC, 'Q_alv' = Q_alv, 'Vplasma'=Vplasma,
                'VVen'=VVen, 'VArt'=VArt, 'VBC' = VBC,
                'VLC' = VLC, 'VKC' = VKC, 'VluC' = VluC, 'VGIC' = VGIC,
                'VFC' = VFC, 'VSKC' = VSKC, 'VRC' = VRC, 'VSC' = VSC, 'V_tubule' = V_tubule,
                'QLC' = QLC, 'QGIC' = QGIC, 'QKC' = QKC, 'QFC' = QFC, 'QSKC' = QSKC, 
                'QRC' = QRC, 'QSC' = QSC, 'GFR' = GFR, 'Q_urine' = Q_urine,
                
                # EG partition coefficients
                'PB_EG' = PB_EG, 'PBS_EG' = PBS_EG, 'PL_EG' = PL_EG, 'PK_EG' = PK_EG,
                'PLu_EG' = PLu_EG, 'PF_EG' = PF_EG, 'PSK_EG' = PSK_EG, 'PGI_EG' = PGI_EG,
                'PR_EG' = PR_EG, 'PS_EG' = PS_EG,
                
                # GA partition coefficients
                'PBS_GA' = PBS_GA, 'PL_GA' = PL_GA, 'PK_GA' = PK_GA, 'PF_GA' = PF_GA,
                'PSK_GA' = PSK_GA, 'PGI_GA' = PGI_GA, 'PR_GA' = PR_GA, 'PS_GA' = PS_GA,
                
                # Other parameters
                'KAS' = KAS, 'CL_renal' = CL_renal, 'KM1_mg' = KM1_mg, 'VMAX1' = VMAX1,
                'KM2_mg' = KM2_mg, 'VMAX2' = VMAX2, 'KT' = KT, 'TMAX' = TMAX
    ))
  })
}  


ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
  # Calculate tissue concentrations (mg/L)
#======================================================================================== 
    C_lung_EG <- A_lung_EG / VluC
    C_liver_EG <- A_liver_EG / VLC
    C_kidney_EG <- A_kidney_EG / VKC
    C_fat_EG <- A_fat_EG / VFC
    C_skin_EG <- A_skin_EG / VSKC
    C_GI_EG <- A_GI_EG / VGIC
    C_rich_EG <- A_rich_EG / VRC
    C_slow_EG <- A_slow_EG / VSC
    C_ven_EG <- A_ven_EG/VVen
    C_art_EG <- A_art_EG/VArt
    C_plasma_EG <- (A_ven_EG + A_art_EG)/VBC
    
    C_liver_GA <- A_liver_GA / VLC
    C_kidney_GA <- A_kidney_GA / VKC
    C_fat_GA <- A_fat_GA / VFC
    C_skin_GA <- A_skin_GA / VSKC
    C_GI_GA <- A_GI_GA / VGIC
    C_rich_GA <- A_rich_GA / VRC
    C_slow_GA <- A_slow_GA / VSC
    C_tubule_GA <- A_tubule_GA / V_tubule
    C_ven_GA <- A_ven_GA/VVen
    C_art_GA <- A_art_GA/VArt
    C_plasma_GA <- (A_ven_GA + A_art_GA)/VBC
    
    
 #ODEs
#======================================================================================== 
    #Venous blood
    dA_ven_EG <- (QLC+QGIC) * C_liver_EG / PL_EG + QKC * C_kidney_EG / PK_EG + QFC * C_fat_EG / PF_EG +
                 QSKC * C_skin_EG / PSK_EG + QRC * C_rich_EG / PR_EG + QSC * C_slow_EG / PS_EG -
                 QC * C_ven_EG
    
    dA_ven_GA <- (QLC+QGIC) * C_liver_GA / PL_GA + QKC * C_kidney_GA / PK_GA + QFC * C_fat_GA / PF_GA +
                 QSKC * C_skin_GA / PSK_GA + QRC * C_rich_GA / PR_GA + QSC * C_slow_GA / PS_GA -
                 QC * C_ven_GA
    
    #Arterial blood
    dA_art_EG <- QC * C_lung_EG / PLu_EG - QLC * C_art_EG - QKC * C_art_EG - QFC * C_art_EG -
                 QSKC * C_art_EG - QGIC * C_art_EG - QRC * C_art_EG - QSC * C_art_EG - 
                 CL_renal * C_art_EG - GFR * C_art_GA
    
    dA_art_GA <- QC * C_ven_GA - QLC * C_art_GA - QKC * C_art_GA - QFC * C_art_GA -
                 QSKC * C_art_GA - QGIC * C_art_GA - QRC * C_art_GA - QSC * C_art_GA 
    
    
    
    # Differential equations for EG compartments
    dA_lung_EG <- QC * C_ven_EG - QC * C_lung_EG / PLu_EG
    dA_liver_EG <- QLC * C_art_EG - (QLC+QGIC) * C_liver_EG / PL_EG - 
                   (VMAX1 * C_liver_EG) / (KM1_mg + C_liver_EG) + QGIC * C_GI_EG / PGI_EG
    dA_kidney_EG <- QKC * C_art_EG - QKC * C_kidney_EG / PK_EG
    dA_fat_EG <- QFC * C_art_EG - QFC * C_fat_EG / PF_EG
    dA_skin_EG <- QSKC * C_art_EG - QSKC * C_skin_EG / PSK_EG
    dA_GI_EG <- QGIC * C_art_EG - QGIC * C_GI_EG / PGI_EG
    dA_rich_EG <- QRC * C_art_EG - QRC * C_rich_EG / PR_EG
    dA_slow_EG <- QSC * C_art_EG - QSC * C_slow_EG / PS_EG
    
    
    # Differential equations for GA compartments
    dA_liver_GA <- QLC * C_art_GA - (QLC+QGIC) * C_liver_GA / PL_GA + QGIC * C_GI_GA / PGI_GA +
                   (VMAX1 * C_liver_EG) / (KM1_mg + C_liver_EG) - (VMAX2 * C_liver_GA) / (KM2_mg + C_liver_GA)
    dA_kidney_GA <- QKC * C_art_GA - QKC * C_kidney_GA / PK_GA + (TMAX * C_tubule_GA) / (KT + C_tubule_GA)
    dA_fat_GA <- QFC * C_art_GA - QFC * C_fat_GA / PF_GA
    dA_skin_GA <- QSKC * C_art_GA - QSKC * C_skin_GA / PSK_GA
    dA_GI_GA <- QGIC * C_art_GA - QGIC * C_GI_GA / PGI_GA
    dA_rich_GA <- QRC * C_art_GA - QRC * C_rich_GA / PR_GA
    dA_slow_GA <- QSC * C_art_GA - QSC * C_slow_GA / PS_GA
    
    
    # Tubule urine GA
    dA_tubule_GA <- GFR * C_art_GA - (TMAX * C_tubule_GA) / (KT + C_tubule_GA) - Q_urine * C_tubule_GA
    
    # Cumulative amounts
    dCum_urine_EG <- CL_renal * C_art_EG
    dCum_urine_GA <- Q_urine * C_tubule_GA
    dCum_met_GA <- (VMAX2 * C_liver_GA) / (KM2_mg + C_liver_GA)
    
    # Return derivatives and concentrations
    return(list(c('dA_ven_EG'=dA_ven_EG, 'dA_ven_GA'=dA_ven_GA,
                  'dA_art_EG'=dA_art_EG, 'dA_art_GA'=dA_art_GA,
                  'dA_lung_EG'=dA_lung_EG, 'dA_liver_EG'=dA_liver_EG,
                  'dA_kidney_EG'=dA_kidney_EG, 'dA_fat_EG'=dA_fat_EG, 'dA_skin_EG'=dA_skin_EG,
                  'dA_GI_EG'=dA_GI_EG, 'dA_rich_EG'=dA_rich_EG, 'dA_slow_EG'=dA_slow_EG,
                  'dA_liver_GA'=dA_liver_GA, 'dA_kidney_GA'=dA_kidney_GA,
                  'dA_fat_GA'=dA_fat_GA, 'dA_skin_GA'=dA_skin_GA, 
                  'dA_GI_GA'=dA_GI_GA, 'dA_rich_GA'=dA_rich_GA, 'dA_slow_GA'=dA_slow_GA,
                  'dA_tubule_GA'=dA_tubule_GA,'dCum_urine_EG'=dCum_urine_EG,
                  'dCum_urine_GA'=dCum_urine_GA, 'dCum_met_GA'=dCum_met_GA),
                
                'C_lung_EG' = C_lung_EG, 'C_liver_EG' = C_liver_EG, 'C_kidney_EG' = C_kidney_EG,
                'C_fat_EG' = C_fat_EG, 'C_skin_EG' = C_skin_EG, 'C_GI_EG' = C_GI_EG,
                'C_rich_EG' = C_rich_EG, 'C_slow_EG' = C_slow_EG,
                'C_liver_GA' = C_liver_GA, 'C_kidney_GA' = C_kidney_GA, 'C_fat_GA' = C_fat_GA,
                'C_skin_GA' = C_skin_GA, 'C_GI_GA' = C_GI_GA, 'C_rich_GA' = C_rich_GA,
                'C_slow_GA' = C_slow_GA, 'C_tubule_GA' = C_tubule_GA,
                'C_art_EG' = C_art_EG, 'C_art_GA' = C_art_GA, 'C_ven_EG' = C_ven_EG,
                'C_ven_GA' = C_ven_GA, 'C_plasma_EG'=C_plasma_EG, 'C_plasma_GA'=C_plasma_GA
                
                
    ))
    
  })
}


create.inits <- function(parameters){
  with(as.list(parameters),{
    

    A_ven_EG <- 0; A_ven_GA <- 0; A_art_EG <-0; A_art_GA <-0;
    A_lung_EG <- 0 ; A_liver_EG <- 0; A_kidney_EG <- 0; A_fat_EG <- 0; A_skin_EG <- 0;
    A_GI_EG <- 0; A_rich_EG <- 0; A_slow_EG <- 0; A_liver_GA <- 0;
    A_kidney_GA <- 0; A_fat_GA <- 0; A_skin_GA <- 0; A_GI_GA <- 0; A_rich_GA <- 0;
    A_slow_GA <- 0; A_tubule_GA <- 0; Cum_urine_EG <- 0; 
    Cum_urine_GA <- 0; Cum_met_GA <-0
    
    return(c(  'A_ven_EG'= A_ven_EG, 'A_ven_GA'= A_ven_GA,
               'A_art_EG'= A_art_EG, 'A_art_GA'= A_art_GA,
               'A_lung_EG'= A_lung_EG, 'A_liver_EG'= A_liver_EG,
               'A_kidney_EG'= A_kidney_EG, 'A_fat_EG'= A_fat_EG, 'A_skin_EG'= A_skin_EG,
               'A_GI_EG'= A_GI_EG, 'A_rich_EG'= A_rich_EG, 'A_slow_EG'= A_slow_EG,
               'A_liver_GA'= A_liver_GA, 'A_kidney_GA'= A_kidney_GA,
               'A_fat_GA'= A_fat_GA, 'A_skin_GA'= A_skin_GA, 
               'A_GI_GA'= A_GI_GA, 'A_rich_GA'= A_rich_GA, 'A_slow_GA'= A_slow_GA,
               'A_tubule_GA'= A_tubule_GA,'Cum_urine_EG'= Cum_urine_EG,
               'Cum_urine_GA'= Cum_urine_GA, 'Cum_met_GA'= Cum_met_GA
      
      
    ))
  })
}

create.events <- function(parameters){
  with(as.list(parameters), {
    
    # Calculate number of administrated doses and corresponding administration time
    ldose <- length(admin.dose)
    ltimes <- length(admin.time)
    # If not equal, then stop 
    if (ltimes != ldose){
      stop("The times of administration should be equal in number to the doses")
    }else{
      if (admin.type == "ip"){
        events <- list(data = rbind(data.frame(var = c("A_liver_EG"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "oral"){
        events <- list(data = rbind(data.frame(var = c("A_GI_EG"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "dermal"){
        events <- list(data = rbind(data.frame(var = c("A_skin_EG"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "inh"){
        events <- list(data = rbind(data.frame(var = c("A_lung_EG"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "SC"){
        events <- list(data = rbind(data.frame(var = c("A_ven_EG"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      } 
    }
    return(events)
  })
}


################################################################################

setwd("Corley et al.2005/Data")

# Read data
EG_vol1_plasma <- openxlsx::read.xlsx("EGvol1_plasma_fit.xlsx")
EG_vol2_plasma <- openxlsx::read.xlsx("EGvol2_plasma_fit.xlsx")
GA_vol1_plasma <- openxlsx::read.xlsx("GAvol1_plasma_fit.xlsx")
GA_vol2_plasma <- openxlsx::read.xlsx("GAvol2_plasma_fit.xlsx")
EG_vol1_urine <- openxlsx::read.xlsx("EGvol1_urine_fit.xlsx")
EG_vol2_urine <- openxlsx::read.xlsx("EGvol2_urine_fit.xlsx")
GA_vol1_urine <- openxlsx::read.xlsx("GAvol1_urine_fit.xlsx")
GA_vol2_urine <- openxlsx::read.xlsx("GAvol2_urine_fit.xlsx")


setwd("Corley et al.2005")

dataset <- list("df1" = EG_vol1_plasma, "df2" = EG_vol2_plasma, "df3" = GA_vol1_plasma, "df4" = GA_vol2_plasma,
                "df5" = EG_vol1_urine, "df6" = EG_vol2_urine, "df7" = GA_vol1_urine, "df8" = GA_vol2_urine)


#################################
#--------------------------------
# EG inhalation scenario - volunteer 1
#--------------------------------
#################################

# Set up parameters for event-based exposure
BW <- 96

admin.dose_peg_kg <- 0.96 #mg/kg BW
admin.dose_total <- admin.dose_peg_kg*BW

n_inj <- 16
interval_min <- 15
inj_times_h <- seq(0, by = interval_min/60, length.out = n_inj) # 0..3.75 h
dose_per_inj_mg <- admin.dose_total/n_inj

admin.time <- inj_times_h
admin.dose <- rep(dose_per_inj_mg, n_inj) # mg per pulse
admin.type <- "inh"
sex <- "M"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 16, 0.1)  #h

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df1=========================================================

exp_data <- dataset$df1 
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
column_names <- c("C_plasma_EG")
preds_EG_vol1_plasma <- list()

for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_EG_vol1_plasma [[i]] <- solution[solution$time %in% exp_time, column_names[i]]} #mg/L

obs_EG_vol1_plasma <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"]) #mg/L

#======================================df3=========================================================

exp_data <- dataset$df3 
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
column_names <- c("C_plasma_GA")
preds_GA_vol1_plasma <- list()

for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_GA_vol1_plasma [[i]] <- solution[solution$time %in% exp_time, column_names[i]]} #mg/L

obs_GA_vol1_plasma <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"]) #mg/L


#======================================df5=========================================================
# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 20, 0.1)  #h

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))

exp_data <- dataset$df5 
colnames(exp_data)[c(2,3)] <- c("time", "amount")
column_names <- c("Cum_urine_EG")
preds_EG_vol1_urine <- list()

for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_EG_vol1_urine [[i]] <- solution[solution$time %in% exp_time, column_names[i]]} #mg/L

obs_EG_vol1_urine <- list(exp_data[exp_data$Tissue == "Urine", "amount"]) #mg/L


#======================================df7=========================================================

exp_data <- dataset$df7 
colnames(exp_data)[c(2,3)] <- c("time", "amount")
column_names <- c("Cum_urine_GA")
preds_GA_vol1_urine <- list()

for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_GA_vol1_urine [[i]] <- solution[solution$time %in% exp_time, column_names[i]]} #mg/L

obs_GA_vol1_urine <- list(exp_data[exp_data$Tissue == "Urine", "amount"]) #mg/L



#################################
#--------------------------------
# EG inhalation scenario - volunteer 2
#--------------------------------
#################################

# Set up parameters for event-based exposure
BW <- 57

admin.dose_peg_kg <- 1.51 #mg/kg BW
admin.dose_total <- admin.dose_peg_kg*BW

n_inj <- 16
interval_min <- 15
inj_times_h <- seq(0, by = interval_min/60, length.out = n_inj) # 0..3.75 h
dose_per_inj_mg <- admin.dose_total/n_inj

admin.time <- inj_times_h
admin.dose <- rep(dose_per_inj_mg, n_inj) # mg per pulse
admin.type <- "inh"
sex <- "M"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 16, 0.1)  #h

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df2=========================================================

exp_data <- dataset$df2 
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
column_names <- c("C_plasma_EG")
preds_EG_vol2_plasma <- list()

for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_EG_vol2_plasma [[i]] <- solution[solution$time %in% exp_time, column_names[i]]} #mg/L

obs_EG_vol2_plasma <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"]) #mg/L

#======================================df4=========================================================

exp_data <- dataset$df4 
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
column_names <- c("C_plasma_GA")
preds_GA_vol2_plasma <- list()

for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_GA_vol2_plasma [[i]] <- solution[solution$time %in% exp_time, column_names[i]]} #mg/L

obs_GA_vol2_plasma <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"]) #mg/L

#======================================df6=========================================================
# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 20, 0.1)  #h

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))


exp_data <- dataset$df6 
colnames(exp_data)[c(2,3)] <- c("time", "amount")
column_names <- c("Cum_urine_EG")
preds_EG_vol2_urine <- list()

for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_EG_vol2_urine [[i]] <- solution[solution$time %in% exp_time, column_names[i]]} #mg/L

obs_EG_vol2_urine <- list(exp_data[exp_data$Tissue == "Urine", "amount"]) #mg/L


#======================================df8=========================================================

exp_data <- dataset$df8 
colnames(exp_data)[c(2,3)] <- c("time", "amount")
column_names <- c("Cum_urine_GA")
preds_GA_vol2_urine <- list()

for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_GA_vol2_urine [[i]] <- solution[solution$time %in% exp_time, column_names[i]]} #mg/L

obs_GA_vol2_urine <- list(exp_data[exp_data$Tissue == "Urine", "amount"]) #mg/L


#Simulations and plots
################################################################################
#EG inhalation scenario - volunteer 1
################################################################################

# Set up parameters for event-based exposure
BW <- 96

admin.dose_peg_kg <- 0.96 #mg/kg BW
admin.dose_total <- admin.dose_peg_kg*BW

n_inj <- 16
interval_min <- 15
inj_times_h <- seq(0, by = interval_min/60, length.out = n_inj) # 0..3.75 h
dose_per_inj_mg <- admin.dose_total/n_inj

admin.time <- inj_times_h
admin.dose <- rep(dose_per_inj_mg, n_inj) # mg per pulse
admin.type <- "inh"
sex <- "M"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 16, 0.1)  #h

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))


preds_EG_vol1_plasma <- solution[, c("time","C_plasma_EG")]
preds_GA_vol1_plasma <- solution[, c("time","C_plasma_GA")]

# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 20, 0.1)  #h

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))

preds_EG_vol1_urine <- solution[, c("time","Cum_urine_EG")]
preds_GA_vol1_urine <- solution[, c("time","Cum_urine_GA")]

################################################################################
#EG inhalation scenario - volunteer 2
################################################################################

# Set up parameters for event-based exposure
BW <- 57

admin.dose_peg_kg <- 1.51 #mg/kg BW
admin.dose_total <- admin.dose_peg_kg*BW

n_inj <- 16
interval_min <- 15
inj_times_h <- seq(0, by = interval_min/60, length.out = n_inj) # 0..3.75 h
dose_per_inj_mg <- admin.dose_total/n_inj

admin.time <- inj_times_h
admin.dose <- rep(dose_per_inj_mg, n_inj) # mg per pulse
admin.type <- "inh"
sex <- "M"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 16, 0.1)  #h

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))


preds_EG_vol2_plasma <- solution[, c("time","C_plasma_EG")]
preds_GA_vol2_plasma <- solution[, c("time","C_plasma_GA")]

# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 20, 0.1)  #h

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                           y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-07, atol = 1e-07))

preds_EG_vol2_urine <- solution[, c("time","Cum_urine_EG")]
preds_GA_vol2_urine <- solution[, c("time","Cum_urine_GA")]

# ######################################################################################
#Plot the predictions against the observations
library(ggplot2) 

# Function that creates a plot given a compartment name and the respective predictions and observations
create.plots_concentration <- function(predictions, observations, compartment){  
  #Colours of observations and predictions
  cls <-  c("predictions" = "#56B4E9",  "Observations" = "#D55E00")
  
  ggplot(data = predictions)+
    geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                          color = '"predictions"'),  size=1.5,alpha = 0.7) +
    geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                             color='"Observations"'), size=4)+
    labs(title = rlang::expr(!!compartment), 
         y = expression("Plasma concentration (mg/L)" ),
         x = "Time (hours)")+
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
experiment1 <- reshape(EG_vol1_plasma[c("Tissue" ,"Time_(h)", 
                                        "Concentration_(mg/L)")], 
                       idvar = "Time_(h)", timevar = "Tissue", direction = "wide")
colnames(experiment1) <- c("Time",unique(EG_vol1_plasma$Tissue))

experiment2 <- reshape(EG_vol2_plasma[c("Tissue" ,"Time_(h)", 
                                        "Concentration_(mg/L)")], 
                       idvar = "Time_(h)", timevar = "Tissue", direction = "wide")
colnames(experiment2) <- c("Time",unique(EG_vol2_plasma$Tissue))
                                         
experiment3 <- reshape(GA_vol1_plasma[c("Tissue" ,"Time_(h)", 
                                        "Concentration_(mg/L)")], 
                       idvar = "Time_(h)", timevar = "Tissue", direction = "wide")
colnames(experiment3) <- c("Time",unique(GA_vol1_plasma$Tissue))

experiment4 <- reshape(GA_vol2_plasma[c("Tissue" ,"Time_(h)", 
                                        "Concentration_(mg/L)")], 
                       idvar = "Time_(h)", timevar = "Tissue", direction = "wide")
colnames(experiment4) <- c("Time",unique(GA_vol2_plasma$Tissue))


# Put the experiments in a list
experiments <- list(experiment1 = experiment1, experiment2 = experiment2, experiment3 = experiment3,
                    experiment4 = experiment4)

# Rename predictions so that they share the same name as the names of the experimental data 
colnames(preds_EG_vol1_plasma) <- c( "Time", "Plasma")
colnames(preds_EG_vol2_plasma) <- c( "Time", "Plasma")
colnames(preds_GA_vol1_plasma) <- c( "Time", "Plasma")
colnames(preds_GA_vol2_plasma) <- c( "Time", "Plasma")

# Create a list containing the corresponding predictions
simulations <- list(predictions1 = preds_EG_vol1_plasma, predictions2 = preds_EG_vol2_plasma,
                    predictions3 = preds_GA_vol1_plasma, predictions4 = preds_GA_vol2_plasma)


# Iterate over all existing experiments and create the accompanying plots
for(i in 1:length(experiments)){
  # Retrieve the corresponding observations and simulations
  observations <- experiments[[i]]
  predictions <- simulations[[i]]
  # Extract the compartment names
  compartments <- names(predictions)[2:length(predictions)]
  
  # Use lapply to iterate over the column names and create plots
  plots <- lapply(compartments, function(compartment) {
    create.plots_concentration(predictions, observations, compartment )
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
  ggsave(paste0("experiment_plasma", i,".png"), plot = final_plot,
         device = 'png', dpi = 300,
         width = 13,
         height = 10,
         units = "in")
}

create.plots_amount <- function(predictions, observations, compartment){  
  #Colours of observations and predictions
  cls <-  c("predictions" = "#56B4E9",  "Observations" = "#D55E00")
  
  ggplot(data = predictions)+
    geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                          color = '"predictions"'),  size=1.5,alpha = 0.7) +
    geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                             color='"Observations"'), size=4)+
    labs(title = rlang::expr(!!compartment), 
         y = expression("Urine cumulative amount (mg)" ),
         x = "Time (hours)")+
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
experiment5 <- reshape(EG_vol1_urine[c("Tissue" ,"Time_(h)", 
                                        "Amount_mg")], 
                       idvar = "Time_(h)", timevar = "Tissue", direction = "wide")
colnames(experiment5) <- c("Time",unique(EG_vol1_urine$Tissue))

experiment6 <- reshape(EG_vol2_urine[c("Tissue" ,"Time_(h)", 
                                        "Amount_mg")], 
                       idvar = "Time_(h)", timevar = "Tissue", direction = "wide")
colnames(experiment6) <- c("Time",unique(EG_vol2_urine$Tissue))

experiment7 <- reshape(GA_vol1_urine[c("Tissue" ,"Time_(h)", 
                                        "Amount_mg")], 
                       idvar = "Time_(h)", timevar = "Tissue", direction = "wide")
colnames(experiment7) <- c("Time",unique(GA_vol1_urine$Tissue))

experiment8 <- reshape(GA_vol2_urine[c("Tissue" ,"Time_(h)", 
                                        "Amount_mg")], 
                       idvar = "Time_(h)", timevar = "Tissue", direction = "wide")
colnames(experiment8) <- c("Time",unique(GA_vol2_urine$Tissue))


# Put the experiments in a list
experiments <- list(experiment5 = experiment5, experiment6 = experiment6, experiment7 = experiment7,
                    experiment8 = experiment8)

# Rename predictions so that they share the same name as the names of the experimental data 
colnames(preds_EG_vol1_urine) <- c( "Time", "Urine")
colnames(preds_EG_vol2_urine) <- c( "Time", "Urine")
colnames(preds_GA_vol1_urine) <- c( "Time", "Urine")
colnames(preds_GA_vol2_urine) <- c( "Time", "Urine")


# Create a list containing the corresponding predictions
simulations <- list(predictions5 = preds_EG_vol1_urine, predictions6 = preds_EG_vol2_urine,
                    predictions7 = preds_GA_vol1_urine, predictions8 = preds_GA_vol2_urine)


# Iterate over all existing experiments and create the accompanying plots
for(i in 1:length(experiments)){
  # Retrieve the corresponding observations and simulations
  observations <- experiments[[i]]
  predictions <- simulations[[i]]
  # Extract the compartment names
  compartments <- names(predictions)[2:length(predictions)]
  
  # Use lapply to iterate over the column names and create plots
  plots <- lapply(compartments, function(compartment) {
    create.plots_amount(predictions, observations, compartment )
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
  ggsave(paste0("experiment_urine", i,".png"), plot = final_plot,
         device = 'png', dpi = 300,
         width = 13,
         height = 10,
         units = "in")
}

save.image("EG-GA_model.RData")