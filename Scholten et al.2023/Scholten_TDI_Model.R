# ============================================================================
# Scholten et al. 2023 - PBK Model for TDI Exposure
# ============================================================================
#amount in ug
#concentration in ug/L
#Kg=L

library(deSolve)

create.params <- function(user.input){
  with(as.list(user.input),{
    
# ============================================================================
#Physiological Parameters
# ============================================================================

#Chemical specific parameters
# ============================================================================
#Molecular weights (g/mol)
MW <- 174.16
MWmet <- 122.17
LogP <- 3.74
VP <- 2.7 #vapor pressure, Pa
Talbmcr_TDI <- 0.167 #TDI Half-life of binding to albumin and other macromolecules, h
Talbmcr_MDI <- 1.5 #MDI Half-life of binding to albumin and other macromolecules, h
Falb <- 0.2 #Proportion bound to albumin (%)
Tvas_TDI <- 0.083 #TDI Half-life transfer from lung interstitial to lung vascular, h
Tvas_MDI <- 14 #MDI Half-life transfer from lung interstitial to lung vascular, h

# ============================================================================

# Body weight (kg)
BW <- 85
#Hematocrit (%)
Fhem <- 0.45

# ============================================================================

# Tissue fractional volumes
VFbld <- 0.073 #Blood
VFlun <- 0.015 # Liver
VFskn <- 0.045 #Skin
VFgut <- 0.016 #Gut
VFlumen <- 0.9 #Gut lumen
VFkidney <- 0.004  #Kidney volume

#Excretion parameters
Tfeces <- 2.3 #1/h, Half-life transfer lumen to feces
UPR <- 0.00125 #urine production, L/kg/h
Ccreat <- 1 #creatinine in urine, g/L
Telim <- 456 #Half-life elimination from albumin, h

#Skin parameters
SAderm <- 100 #Exposed surface area of skin, cm2
Hscve <- 0.012 #Height of skin barrier, cm
Hsurf <- 0.01 #Height of skin surface, cm
Fatepi <- 0.02 #Fraction of fat in epidermis
Fatbld <- 0.07 #Fraction of fat in blood
Talb_scve <- 2.8*10-6 #Half-life of binding to albumin in scve, H

#============================================================================================
#Flow rates (L/h)
#============================================================================================

Qc <- 390 #cardiac output
Qgfr <- 5.5 #glomerular filtration rate

# Fractional blood flows to tissue blood flow
QFskn <- 0.05   #Skin
QFgut <- 0.15  #Gut

#Blood flows
Qupr <- UPR*BW #Urine production rate, L/h
Qgut <- QFgut*Qc #Blood flow to the gut, L/h 
Qskn <- QFskn*Qc #Blood flow to the skin, L/h

#Fat-content adjusted stratum corneum permeability coefficient transferred into a clearance value, L/h
Qscve <- 100.74*LogP - 0.006*MW - 2.8 *SAderm/1000 

# ============================================================================
# Calculated Parameters
# ============================================================================

#Binding constants
#=============================================================================

Kvas_TDI <- log(2)/Tvas_TDI #TDI Transfer rate constant for transfer from the interstitial to the vascular lung compartment, 1/h
Kvas_MDI <- log(2)/Tvas_MDI #MDI Transfer rate constant for transfer from the interstitial to the vascular lung compartment, 1/h

Kalb_TDI <- Falb*log(2)/Talbmcr_TDI #TDI Binding rate constants for binding of unbound diisocyanates to albumin, 1/h
Kalb_MDI <- Falb*log(2)/Talbmcr_MDI #MDI Binding rate constants for binding of unbound diisocyanates to albumin, 1/h

Kmcr_TDI <- (1-Falb)*(log(2)/Talbmcr_TDI) #TDI Binding rate constants for binding of unbound diisocyanates to macromolecules, 1/h
Kmcr_MDI <- (1-Falb)*(log(2)/Talbmcr_MDI) #MDI Binding rate constants for binding of unbound diisocyanates to macromolecules, 1/h

Kelim <- log(2)/Telim #Elimination rate constant for the degradation of albumin into (smaller) macromolecules), 1/h
Thydro <- 1.0  # Half-life for hydrolysis (assumed value), h
khydro <- log(2)/Thydro #Rate constant for hydrolysis of unbound diisocyanates, 1/h
Kfeces <- log(2)/Tfeces #Transfer rate constant from the gut lumen to the feces, 1/h

Kalb_scve <- log(2)/Talb_scve #Binding rate constant of binding to albumin in the skin
Pbld_scve <- (1-Fatepi + Fatepi*10^LogP)/(1-Fatbld + Fatbld* 10^LogP) #Skin barrier (stratum corneum + viable epidermis) to blood partition coefficient


# Tissue volumes (L)
#=============================================================================
Vbld <- VFbld * BW
Vlun <- VFlun * BW  
Vskn <- VFskn * BW
Vgut <- VFgut * BW
Vkidney <- VFkidney * BW
Vscve <- SAderm * Hscve / 1000  # Skin barrier volume (L)
Vsurf <- SAderm * Hsurf / 1000  # Skin surface volume (L)
Vrem <- BW - Vbld - Vlun - Vskn - Vgut - Vkidney  # Remaining tissue volume


#==========================================================================================================================
#Other parameters
#==========================================================================================================================

#Evaporation parameters
#kevap = cevap * (MW*8.7*VP*SAderm)/8.314*303*10

# Additional parameters needed for simulation
Dderm_ub <- 0  # Dermal deposition rate (mg/h/cm2)
Qoral <- 0  # Oral dose rate (mg/h)
Kabs <- 1   # Absorption rate constant from gut lumen (1/h)
Qkid <- 0.2 * Qc  # Kidney blood flow (L/h)
Vscve <- SAderm * Hscve / 1000  # Skin barrier volume (L)
Vsurf <- SAderm * Hsurf / 1000  # Skin surface volume (L)


return(list(
  'MW'=MW, 'MWmet'=MWmet, 'LogP'=LogP, 'VP'=VP,
  'Talbmcr_TDI'=Talbmcr_TDI, 'Talbmcr_MDI'=Talbmcr_MDI, 'Falb'=Falb,
  'Tvas_TDI'=Tvas_TDI, 'Tvas_MDI'=Tvas_MDI, 'BW'=BW,
  'VFbld'=VFbld, 'VFlun'=VFlun, 'VFskn'=VFskn, 'VFgut'=VFgut, 'VFlumen'=VFlumen,
  'VFkidney'=VFkidney, 'Qc'=Qc, 'Qgfr'=Qgfr, 'QFskn'=QFskn, 'QFgut'=QFgut,
  'Tfeces'=Tfeces, 'UPR'=UPR, 'Ccreat'=Ccreat, 'Telim'=Telim,
  'SAderm'=SAderm, 'Hscve'=Hscve, 'Hsurf'=Hsurf, 'Fatepi'=Fatepi, 'Fatbld'=Fatbld,
  'Talb_scve'=Talb_scve, 'Kvas_TDI'=Kvas_TDI, 'Kvas_MDI'=Kvas_MDI, 'Kalb_TDI'=Kalb_TDI,
  'Kalb_MDI'=Kalb_MDI, 'Kmcr_TDI'=Kmcr_TDI, 'Kmcr_MDI'=Kmcr_MDI, 'Kelim'=Kelim,
  'Thydro'=Thydro, 'khydro'=khydro, 'Kfeces'=Kfeces, 'Kalb_scve'=Kalb_scve,
  'Pbld_scve'=Pbld_scve, 'Qscve'=Qscve, 'Qupr'=Qupr, 'Qgut'=Qgut,
  'Qskn'=Qskn, 'Dderm_ub'=Dderm_ub, 'Qoral'=Qoral, 'Kabs'=Kabs,
  #'kevap'=kevap,
  'Vbld'=Vbld, 'Vlun'=Vlun, 'Vskn'=Vskn, 'Vgut'=Vgut, 'Vkidney'=Vkidney,
  'Vscve'=Vscve, 'Vsurf'=Vsurf, 'Vrem'=Vrem, 
  'Qkid'=Qkid, 'Vscve'=Vscve, 'Vsurf'=Vsurf, 'Cair'=Cair, 'Qbr'=Qbr,
  'Fabs'=Fabs, 'admin.dose'=admin.dose, 'admin.dose_gut'=admin.dose_gut,
  'admin.time' = admin.time, 'admin.type' = admin.type

))

  })
}  


# ============================================================================
# Differential Equations for PBK Model
# ============================================================================

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    
    # Concentrations (mg/L) for all fractions
    Clung_int_ub <- Alung_int_ub / Vlun
    Clung_int_alb <- Alung_int_alb / Vlun
    Clung_int_mcr <- Alung_int_mcr / Vlun
    Clung_vas_ub <- Alung_vas_ub / Vlun
    Clung_vas_alb <- Alung_vas_alb / Vlun
    Clung_vas_mcr <- Alung_vas_mcr / Vlun
    Cart_ub <- Aart_ub / Vbld
    Cart_alb <- Aart_alb / Vbld
    Cart_mcr <- Aart_mcr / Vbld
    Cven_ub <- Aven_ub / Vbld
    Cven_alb <- Aven_alb / Vbld
    Cven_mcr <- Aven_mcr / Vbld
    Cskn_ub <- Askn_ub / Vskn
    Cskn_mcr <- Askn_mcr / Vskn
    Clum <- Alum / Vgut
    Cscve_ub <- Ascve_ub / Vscve
    Cscve_mcr <- Ascve_mcr / Vscve
    Csurf_ub <- Asurf_ub / Vsurf
    Ckid_ub <- Akid_ub / Vkidney
    Ckid_mcr <- Akid_mcr / Vkidney
    Crem_ub <- Arem_ub / Vrem
    Curine <- Aurine / Vuri
    Curine_amine <- Auri_amine / Vuri
  
    
    # ========== LUNG COMPARTMENT ==========
    
    
    # Lung Interstitial - Unbound
    dAlung_int_ub <-  - Kalb_TDI * Alung_int_ub - Kmcr_TDI * Alung_int_ub -
                     Kvas_TDI * Alung_int_ub

    # Lung Interstitial - Albumin-bound
    dAlung_int_alb <- Kalb_TDI * Alung_int_ub - Kelim * Alung_int_alb - Kvas_TDI * Alung_int_alb
    
    # Lung Interstitial - Macromolecule-bound
    dAlung_int_mcr <- Kmcr_TDI * Alung_int_ub + Kelim * Alung_int_alb - Kvas_TDI * Alung_int_mcr
     
    # Lung Vascular - Unbound
    dAlung_vas_ub <- Kvas_TDI * Alung_int_ub - Kalb_TDI * Alung_vas_ub - Kmcr_TDI * Alung_vas_ub +
                     Qc * Cven_ub - Qc * Clung_vas_ub
    
    # Lung Vascular - Albumin-bound
    dAlung_vas_alb <- Kvas_TDI * Alung_int_alb + Kalb_TDI * Alung_vas_ub - Kelim * Alung_vas_alb +
                      Qc * Cven_alb - Qc * Clung_vas_alb
    
    # Lung Vascular - Macromolecule-bound
    dAlung_vas_mcr <- Kvas_TDI * Alung_int_mcr + Kmcr_TDI * Alung_vas_ub + Kelim * Alung_vas_alb + 
                      Qc * Cven_mcr - Qc * Clung_vas_mcr
    
    
    # ========== BLOOD COMPARTMENTS ==========
    
    # Arterial Blood - Unbound
    dAart_ub <- Qc * Clung_vas_ub - Qc * Cart_ub - Qskn * Cart_ub - Qkid * Cart_ub -
                Kalb_TDI*Aart_ub - Kmcr_TDI*Aart_ub
    
    # Arterial Blood - Albumin-bound  
    dAart_alb <- Qc * Clung_vas_alb - Qc * Cart_alb + Kalb_TDI*Aart_ub - Kelim*Aart_alb
    
    # Arterial Blood - Macromolecule-bound
    dAart_mcr <- Qc * Clung_vas_mcr - Qc * Cart_mcr - Qskn * Cart_mcr - Qkid * Cart_mcr +
                 Kmcr_TDI*Aart_ub + Kelim*Aart_alb 
    
    # Venous Blood - Unbound
    dAven_ub <- Qskn * Cskn_ub + Qc * Crem_ub + Qkid * Ckid_ub - Qc * Cven_ub - Kalb_TDI*Aven_ub - Kmcr_TDI*Aven_ub
    
    # Venous Blood - Albumin-bound
    dAven_alb <- Qc * Cart_alb - Qc * Cven_alb + Kalb_TDI*Aven_ub - Kelim*Aven_alb
    
    # Venous Blood - Macromolecule-bound
    dAven_mcr <- Qskn * Cskn_mcr + Qkid * Ckid_mcr  + Qc * Cart_mcr  - Qc * Cven_mcr   + 
                 Kmcr_TDI*Aven_ub + Kelim*Aven_alb
    
      
    # ========== SKIN COMPARTMENT ==========
    
    #Dermal exposure
    dDerm <- - SAderm*Dderm_ub
    
    #Evaporation
    #devap <- kevap*Asurf_ub
    
    # Skin Surface - Unbound
    dAsurf_ub <- SAderm*Dderm_ub - Qscve * Csurf_ub + Qscve * Cscve_ub #- kevap*Asurf_ub
    
    # Skin Barrier - Unbound
    dAscve_ub <- Qscve * Csurf_ub - Qscve * Cscve_ub - Kalb_scve * Ascve_ub - Pbld_scve*Qscve*Cscve_ub
    
    # Skin Barrier - Albumin-bound
    dAscve_alb <- Kalb_scve * Ascve_ub - Kelim * Ascve_alb
    
    # Skin Barrier - Macromolecule-bound
    dAscve_mcr <- Kelim * Ascve_alb - Pbld_scve*Qscve*Cscve_mcr
    
    # Skin Tissue - Unbound
    dAskn_ub <- Qskn * Cart_ub - Qskn * Cskn_ub + Pbld_scve*Qscve*Cscve_ub
    
    # Skin Tissue - Macromolecule-bound
    dAskn_mcr <- - Qskn * Cskn_mcr + Pbld_scve*Qscve*Cscve_mcr + Qskn * Cart_mcr
    
    # ========== GUT COMPARTMENT ==========
    
    # Gut Lumen (for oral exposure)
    dAlum <-  - Kfeces*Alum
    
    
    # ========== ELIMINATION PATHWAYS ==========
    
    #Kidney Tissue - Unbound
    dAkid_ub <- Qkid * Cart_ub - Qkid * Ckid_ub
    
    #Kidney Tissue - Macromolecule-bound
    dAkid_mcr <- Qkid * Cart_mcr - Qkid * Ckid_mcr - Qgfr*Cart_mcr
    
    #Amine metabolite
    dAuri_amine <- MWmet/MW * Qgfr*Ckid_mcr
    
    # Urine elimination (albumin-bound + free in blood)
    dAurine <- Qgfr*Cart_mcr - MWmet/MW * Qgfr*Ckid_mcr
    
    # Feces elimination
    dAfeces <- Kfeces*Alum
    
    # ========== REMAINING COMPARTMENTS ==========
    
    # Remaining tissues - Unbound
    dArem_ub <- Qc * Cart_ub - Qc * Crem_ub
    
    dVuri <- Qupr
    
    # Return derivatives grouped by organ
    list(c('dAlung_int_ub' = dAlung_int_ub, 'dAlung_int_alb' = dAlung_int_alb, 'dAlung_int_mcr'=dAlung_int_mcr,
           'dAlung_vas_ub'=dAlung_vas_ub, 'dAlung_vas_alb'=dAlung_vas_alb, 'dAlung_vas_mcr'=dAlung_vas_mcr,
           'dAart_ub'=dAart_ub, 'dAart_alb'=dAart_alb, 'dAart_mcr'=dAart_mcr, 'dAven_ub'=dAven_ub,
           'dAven_alb'=dAven_alb, 'dAven_mcr'=dAven_mcr, 'dDerm'=dDerm, 'dAsurf_ub'=dAsurf_ub, 'dAscve_ub'=dAscve_ub,
           'dAscve_alb'=dAscve_alb, 'dAscve_mcr'=dAscve_mcr, 'dAskn_ub'=dAskn_ub, 'dAskn_mcr'=dAskn_mcr,
           'dAlum'=dAlum, 'dAkid_ub'=dAkid_ub, 'dAkid_mcr'=dAkid_mcr, 'dAuri_amine'=dAuri_amine, 'dAurine'=dAurine,
           'dAfeces'=dAfeces, 'dArem_ub'=dArem_ub, 'dVuri'=dVuri 
     
    ),
    
    # Concentrations for output
    'Clung_int_ub'=Clung_int_ub, 'Clung_int_alb'=Clung_int_alb, 'Clung_int_mcr'=Clung_int_mcr,
    'Clung_vas_ub'=Clung_vas_ub, 'Clung_vas_alb'=Clung_vas_alb, 'Clung_vas_mcr'=dAlung_vas_mcr,
    'Cart_ub'=Cart_ub, 'Cart_alb'=Cart_alb, 'Cart_mcr'=Cart_mcr,
    'Cven_ub'=Cven_ub, 'Cven_alb'=Cven_alb, 'Cven_mcr'=Cven_mcr,
    'Curine'=Curine,  'Curine_amine'=Curine_amine
    )
    
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    
    Alung_int_ub <-0 ; Alung_int_alb <-0  ; Alung_int_mcr <-0 ;
    Alung_vas_ub <-0 ; Alung_vas_alb <-0 ; Alung_vas_mcr <-0 ;
    Aart_ub <-0 ; Aart_alb <-0 ; Aart_mcr <-0 ; Aven_ub <-0;
    Aven_alb <-0 ; Aven_mcr <-0 ; Derm <-0; Asurf_ub <-0 ; Ascve_ub <-0 ;
    Ascve_alb <-0 ; Ascve_mcr <-0 ; Askn_ub <-0 ; Askn_mcr <-0 ;
    Alum <-0 ; Akid_ub <-0 ; Akid_mcr <-0 ; Auri_amine <-0; Aurine <-0 ;
    Afeces <-0 ; Arem_ub <-0 ; Vuri <-0
  
  return(c(
    'Alung_int_ub' =  Alung_int_ub, 'Alung_int_alb' =  Alung_int_alb, 'Alung_int_mcr'= Alung_int_mcr,
    'Alung_vas_ub'= Alung_vas_ub, 'Alung_vas_alb'= Alung_vas_alb, 'Alung_vas_mcr'= Alung_vas_mcr,
    'Aart_ub'= Aart_ub, 'Aart_alb'= Aart_alb, 'Aart_mcr'= Aart_mcr, 'Aven_ub'= Aven_ub,
    'Aven_alb'= Aven_alb, 'Aven_mcr'= Aven_mcr, 'Derm'=Derm, 'Asurf_ub'= Asurf_ub, 'Ascve_ub'= Ascve_ub,
    'Ascve_alb'= Ascve_alb, 'Ascve_mcr'= Ascve_mcr, 'Askn_ub'= Askn_ub, 'Askn_mcr'= Askn_mcr,
    'Alum'= Alum, 'Akid_ub'= Akid_ub, 'Akid_mcr'= Akid_mcr, 'Auri_amine'=Auri_amine, 'Aurine'= Aurine,
    'Afeces'= Afeces, 'Arem_ub'= Arem_ub, 'Vuri'= Vuri
    
  ))
  
  
  })
}
  


#================================================================================================

create.events <- function(parameters){
  with(as.list(parameters), {
    
    # Calculate number of administered doses and corresponding administration time
    ldose <- length(admin.dose)
    ltimes <- length(admin.time)
    # If not equal, then stop 
    if (ltimes != ldose){
      stop("The times of administration should be equal in number to the doses")
    }else{
      if (admin.type == "single_inh"){ 
        exposure_events <- rbind(data.frame(var = c("Alung_int_ub"),  time = admin.time, 
                                               value = admin.dose, method = c("add")),
                                    data.frame(var = c("Alum"),  time = admin.time, 
                                               value = admin.dose_gut, method = c("add"))
        
        )
        
      }else if (admin.type == "chronic_inh"){
        #chronic exposure schedule: 8h/day, 5 days/week, 8 weeks
        exposure_events <- data.frame()
        
        # Loop through 26 weeks
        for (week in 0:25) {
          # Loop through 5 working days (Monday to Friday)
          for (day in 0:4) {
            current_day <- week * 25 + day  # Day number from start
            
            # Start exposure at 9:00 (hour 9 + 24*current_day)
            start_time <- 9 + 24 * current_day
            # End exposure at 17:00 (hour 17 + 24*current_day) 
            end_time <- 17 + 24 * current_day
            
            # Add exposure start event
            exposure_events <- rbind(exposure_events,
                                     data.frame(var = c("Alung_int_ub"), time = start_time, 
                                                value = admin.dose, method = c("add")),
                                     data.frame(var = c("Alum"), time = start_time, 
                                                value = admin.dose_gut, method = c("add"))
            )
            
            # during the 8-hour workday, adding doses every hour
            for (hour in 1:7) {  # Hours 1-7 after start (total 8 hours)
              dose_time <- start_time + hour
              exposure_events <- rbind(exposure_events,
                                       data.frame(var = c("Alung_int_ub"), time = dose_time, 
                                                  value = admin.dose, method = c("add")),
                                       data.frame(var = c("Alum"), time = dose_time, 
                                                  value = admin.dose_gut, method = c("add"))
              )
            }
          }
        }
        
        # Create urination schedule: 8:00, 13:00, 18:00, 23:00 daily
        urination_events <- data.frame()
        
        # Loop through all days (8 weeks = 56 days)
        for (day in 0:55) {
          urination_times_daily <- c(8, 13, 18, 23) + 24 * day
          
          for (urin_time in urination_times_daily) {
            urination_events <- rbind(urination_events,
                                      data.frame(var = c("Vuri"), time = urin_time, 
                                                 value = 0, method = c("replace")),
                                      data.frame(var = c("Auri_amine"), time = urin_time, 
                                                 value = 0, method = c("replace"))
            )
          }
        }
        
        # Combine exposure and urination events
        all_events <- rbind(exposure_events, urination_events)
        # Sort by time
        all_events <- all_events[order(all_events$time), ]
        
        events <- list(data = all_events)
      
      }else if (admin.type == "dermal"){
        exposure_events <- rbind(data.frame(var = c("Asurf_ub"),  time = admin.time, 
                                               value = admin.dose, method = c("add"))
        )
      }
      
      urination_times <- seq(4, max(4368), by = 4)  # Every 4 hours
      
      urination_events <- data.frame(
        var = rep(c("Vuri", "Auri_amine"), each = length(urination_times)),
        time = rep(urination_times, 2),
        value = rep(c(0, 0), each = length(urination_times)),  # Reset to 0
        method = rep(c("replace", "replace"), each = length(urination_times))
      )
      
      # Combine exposure and urination events
      all_events <- rbind(exposure_events, urination_events)
      
      events <- list(data = all_events)
    }
    return(events)
  })
}

################################################################################

setwd("Scholten et al.2023/Data")

# Read data
TDA_urine_data <- openxlsx::read.xlsx("Figure S1_urine TDA.xlsx")
TDA_urine_data_chronic <- openxlsx::read.xlsx("Figure 3A_urine.xlsx")

setwd("Scholten et al.2023")
source("Goodness-of-fit-metrics.R")

dataset <- list("df1" = TDA_urine_data, "df2" = TDA_urine_data_chronic)


  #################################
  #--------------------------------
  # TDI single inhalation scenario
  #--------------------------------
  #################################
  
  # Set up parameters for event-based exposure
  BW <- 85
  Cair <- 0.040 #ug/L, 40 ug/m^3 TDI = 40 ug/1000L TDI
  duration <- 4 #h
  Qbr <- 1000 #breathing rate L/h
  Fabs <- 0.2 #fraction absorbed
  admin.dose <- Qbr*Fabs*Cair*duration
  admin.dose_gut <- (1-Fabs)*Qbr*Cair*duration
  admin.time <- 0  
  admin.type <- "single_inh"
  Ccreat <- 1 # g/L
  sex <- "M"
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.dose_gut"=admin.dose_gut,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs for full week
  sample_time=seq(0, 72, 0.5)  #h
  
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
  column_names <- c("Curine_amine")
  preds_TDA_in <- list()
  
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_TDA_in [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/Ccreat} #mg/L
  
    obs_TDA_in <- list(exp_data[exp_data$Tissue == "Urine", "concentration"]) #ug/g creatinine
  

#   ########################################################################################

    #################################
    #--------------------------------
    # TDI chronic exposure scenario
    #--------------------------------
    #################################
    
    # Set up parameters for event-based exposure
    BW <- 85
    Cair <- 0.020 #ug/L, 20 ug/m^3 TDI = 20 ug/1000L TDI
    duration <- 1 #h, 1 h dose rate
    Qbr <- 1000 #breathing rate L/h
    Fabs <- 0.2 #fraction absorbed
    admin.dose <- Qbr*Fabs*Cair*duration
    admin.dose_gut <- (1-Fabs)*Qbr*Cair*duration
    admin.time <- 0  
    admin.type <- "chronic_inh"
    Ccreat <- 1 # g/L
    sex <- "M"
    
    user_input <- list('BW'=BW,
                       "admin.dose"= admin.dose,
                       "admin.dose_gut"=admin.dose_gut,
                       "admin.time" = admin.time, 
                       "admin.type" = admin.type,
                       "sex" = sex)
    
    params <- create.params(user_input)
    inits <- create.inits(params)
    events <- create.events(params)
    
    # sample_time: a vector of time points to solve the ODEs for full week
    sample_time=seq(0, 4368, 1)  #h
    
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
    column_names <- c("Curine_amine")
    preds_TDA_in_chronic <- list()
    
    for (i in 1:length(unique(exp_data$Tissue))) {
      compartment <- unique(exp_data$Tissue)[i]
      #Retrieve time points at which measurements are available for compartment i
      exp_time <- exp_data[exp_data$Tissue == compartment, 2]
      
      preds_TDA_in_chronic [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/Ccreat} #mg/L
    
    obs_TDA_in_chronic <- list(exp_data[exp_data$Tissue == "Urine", "concentration"]) #ug/g creatinine
    
    



################################################################################
#Single exposure scenario
################################################################################

BW <- 85
Cair <- 0.040 #ug/L, 40 ug/m^3 TDI = 40 ug/1000L TDI
duration <- 4 #h
Qbr <- 1000 #breathing rate L/h
Fabs <- 0.2 #fraction absorbed
admin.dose <- Qbr*Fabs*Cair*duration
admin.dose_gut <- (1-Fabs)*admin.dose
admin.time <- 0  
admin.type <- "single_inh"
Ccreat <- 1 # g/L
sex <- "M"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.dose_gut"=admin.dose_gut,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 72, 0.5)  #h

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-07, atol = 1e-07))


preds_TDA_in <-  solution[, c("time","Curine_amine")]

################################################################################
#Chronic exposure scenario
################################################################################

# Set up parameters for event-based exposure
BW <- 85
Cair <- 0.020 #ug/L, 20 ug/m^3 TDI = 20 ug/1000L TDI
duration <- 1 #h, 1 hour dose rate
Qbr <- 1000 #breathing rate L/h
Fabs <- 0.2 #fraction absorbed
admin.dose <- Qbr*Fabs*Cair*duration
admin.dose_gut <- (1-Fabs)*admin.dose
admin.time <- 0  
admin.type <- "chronic_inh"
Ccreat <- 1 # g/L
sex <- "M"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.dose_gut"=admin.dose_gut,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs for full week
sample_time=seq(0, 4368, 1)  #h

# ode(): The solver of the ODEs
solution <- data.frame(ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-07, atol = 1e-07))

preds_TDA_in_chronic <-  solution[, c("time","Curine_amine")]


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
         y = expression("TDA urinary concentration (" * mu* "ug/g creatinin)" ),
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
experiment1 <- reshape(TDA_urine_data[c("Tissue" ,"Time_(h)", 
                                        "Concentration_(ug/gcreatinine)")], 
                       idvar = "Time_(h)", timevar = "Tissue", direction = "wide")
colnames(experiment1) <- c("Time",unique(TDA_urine_data$Tissue))

experiment2 <- reshape(TDA_urine_data_chronic[c("Tissue" ,"Time_(h)", 
                                        "Concentration_(ug/gcreatinine)")], 
                       idvar = "Time_(h)", timevar = "Tissue", direction = "wide")
colnames(experiment2) <- c("Time",unique(TDA_urine_data_chronic$Tissue))

experiment2_filtered <- experiment2[experiment2$Time >= 4200 & experiment2$Time <= 4368, ]

# Put the experiments in a list
experiments <- list(experiment1 = experiment1, experiment2 = experiment2, experiment3 = experiment2_filtered)

# Rename predictions so that they share the same name as the names of the experimental data 
colnames(preds_TDA_in) <- c( "Time", "Urine")
colnames(preds_TDA_in_chronic) <- c( "Time", "Urine")

preds_TDA_in_chronic_filtered <- preds_TDA_in_chronic[preds_TDA_in_chronic$Time >= 4200 & preds_TDA_in_chronic$Time <= 4368, ]

# Create a list containing the corresponding predictions
simulations <- list(predictions1 = preds_TDA_in, predictions2 = preds_TDA_in_chronic, predictions3 = preds_TDA_in_chronic_filtered)


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
  ggsave(paste0("experiment", i,".png"), plot = final_plot,
         device = 'png', dpi = 300,
         width = 13,
         height = 10,
         units = "in")
}
save.image("TDI_model.RData")