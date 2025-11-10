library(deSolve)

#===============================================
# PFBS Population PBK Model (Deterministic)
#===============================================
# This R implementation uses mean posterior estimates from the
# hierarchical Bayesian PFBS population model.
# Based on the Worley et al. PBK model structure.
#===============================================

create.params <- function(user_input){
  with( as.list(user_input),{
    # User input: BW(kg)

    #Initial time: hours
    if  (missing(time_scale) || is.na(time_scale)){
      time_scale <- 1/24
    }else if (time_scale == "minutes"){
      time_scale <- 60
    }else if (time_scale == "hours"){
      time_scale <- 1
    }else if (time_scale == "days"){
      time_scale <- 1/24
    }else if (time_scale == "weeks"){
      time_scale <- (1/24)/7
    }else if (time_scale == "months"){
      time_scale <- ((1/24)/30)
    }else if (time_scale == "years"){
      time_scale <- ((1/24)/365)
    }
    inv_time_scale <- 1/time_scale

    #Cardiac Output and Bloodflow (as fraction of cardiac output)
    QCC = 12.5*inv_time_scale #cardiac output in L/time_scale/kg^0.75; Brown 1997
    QLC = 0.25 #fraction blood flow to liver; Brown 1997
    QKC = 0.175 #fraction blood flow to kidney; Brown 1997.
    Htc = 0.44 #hematocrit for the rat; Davies 1993

    #Tissue Volumes
    VplasC = 0.0428 #fraction vol. of plasma (L/kg BW); Davies 1993
    VLC = 0.026 #fraction vol. of liver (L/kg BW); Brown 1997
    VKC = 0.004 #fraction vol. of kidney (L/kg BW); Brown 1997
    VfilC = 4e-4	#fraction vol. of filtrate (L/kg BW)
    VPTCC = 1.35e-4 #vol. of proximal tubule cells (L/g kidney) (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)

    #Chemical Specific Parameters
    MW = 300.10	#PFBS molecular mass (g/mol)

    Free = 0.02 #free fraction in plasma

    #Kidney Transport Parameters
    # Basolateral transporter parameters (from Worley model - OAT1/OAT3)
    Vmax_baso_invitro = 439.2 #Vmax of basolateral transporter (pmol/mg protein/min)
    Km_baso = 20100 #Km of basolateral transporter (ug/L)

    # Apical transporter parameters
    Vmax_apical_invitro = 37400 #Vmax of apical transporter (pmol/mg protein/min); base value from Worley model

    # ====================================================================
    # POSTERIOR MEAN ESTIMATES - TO BE FILLED IN BY USER
    # ====================================================================
    Km_apical = 109100.0  # Mean Km_apical from MCMC (ug/L) - REPLACE WITH POSTERIOR MEAN
    RAFbaso = 1.0	# Fixed (for identifiability)
    RAFapi = 7.141e-5	# Mean population RAFapi from MCMC - REPLACE WITH POSTERIOR MEAN

    protein = 2.0e-6	#amount of protein in proximal tubule cells (mg protein/proximal tubule cell)
    GFRC = 24.19*inv_time_scale	#glomerular filtration rate (L/time_scale/kg kidney)

    # ====================================================================
    # PFBS PARTITION COEFFICIENTS - TO BE FILLED IN BY USER
    # ====================================================================
    # Values from Population_model/partition_coefficients_results.csv
    PL = 0.4587 #liver:blood partition coefficient - REPLACE WITH ACTUAL VALUE
    PK = 0.4539 #kidney:blood partition coefficient - REPLACE WITH ACTUAL VALUE
    PR = 0.4104 #rest of body:blood partition coefficient - REPLACE WITH ACTUAL VALUE

    #Rate constants (from Worley model)
    kdif = 0.001*inv_time_scale	#diffusion rate from proximal tubule cells (L/time_scale)
    kabsc = 2.12*inv_time_scale	#rate of absorption of chemical from small intestine to liver (1/(time_scale*BW^-0.25))
    kunabsc = 7.06e-5*inv_time_scale	#rate of unabsorbed dose to appear in feces (1/(time_scale*BW^-0.25))
    GEC = 3.5*inv_time_scale #gastric emptying time (1/(time_scale*BW^-0.25))
    k0C = 1.0*inv_time_scale	#rate of uptake from the stomach into the liver (1/(time_scale*BW^-0.25))

    keffluxc = 0.1*inv_time_scale #rate of clearance from proximal tubule cells into blood (1/(time_scale*BW^-0.25))
    kbilec = 0.0001*inv_time_scale #biliary elimination rate; liver to feces storage (1/(time_scale*BW^-0.25))
    kurinec = 0.063*inv_time_scale #rate of urine elimination from urine storage (1/(time_scale*BW^-0.25))
    kvoid = 0.06974*inv_time_scale  #daily urine volume rate (L/time_scale)

    #Scaled Parameters
    #Cardiac output and blood flows
    QC = QCC*(BW^0.75)*(1-Htc)	#cardiac output in L/time_scale; adjusted for plasma
    QK = (QKC*QC)	#plasma flow to kidney (L/time_scale)
    QL = (QLC*QC)	#plasma flow to liver (L/time_scale)
    QR = QC - QK - QL 	#plasma flow to rest of body (L/time_scale)
    QBal = QC - (QK + QL + QR) #Balance check of blood flows; should equal zero

    #Tissue Volumes
    VPlas = VplasC*BW 	#volume of plasma (L)
    VK = VKC*BW 	#volume of kidney (L)
    MK = VK*1.0*1000	#mass of the kidney (g)
    VKb = VK*0.16	#volume of blood in the kidney (L); fraction blood volume of kidney (0.16) from Brown, 1997
    Vfil = VfilC*BW	#volume of filtrate (L)
    VL = VLC*BW	#volume of liver (L)
    ML = VL*1.05*1000	#mass of the liver (g)

    #Kidney Parameters
    PTC = VKC*1000*6e7	#number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney, assuming density of 1 kg/L)
    VPTC = VK*1000*VPTCC	#volume of proximal tubule cells (L)
    MPTC = VPTC*1000 #mass of the proximal tubule cells (g) (assuming density 1 kg/L)
    VR = (0.93*BW) - VPlas - VPTC - Vfil - VL	#volume of remaining tissue (L);
    VBal = (0.93*BW) - (VR + VL + VPTC + Vfil + VPlas)	#Balance check of tissue volumes; should equal zero

    # Calculate Vmax with PFBS molecular weight and estimated RAFapi
    Vmax_basoC = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1e6)*inv_time_scale #Vmax of basolateral transporters (ug/time_scale/kg BW)
    Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1e6)*inv_time_scale #Vmax of apical transporters (ug/time_scale/kg BW)
    Vmax_baso = Vmax_basoC*BW^0.75	#(ug/time_scale)
    Vmax_apical = Vmax_apicalC*BW^0.75	#(ug/time_scale)

    kbile = kbilec*BW^(-0.25)	#biliary elimination; liver to feces storage (/time_scale)
    kurine = kurinec*BW^(-0.25)	#urinary elimination, from filtrate (/time_scale)
    kefflux = keffluxc*BW^(-0.25)	#efflux clearance rate, from PTC to blood (/time_scale)
    GFR = GFRC*VK	#glomerular filtration rate, scaled to mass of kidney(in kg)(L/time_scale)

    #GI Tract Parameters
    kabs = kabsc*BW^(-0.25)	#rate of absorption of chemical from small intestine to liver (/time_scale)
    kunabs = kunabsc*BW^(-0.25)	#rate of unabsorbed dose to appear in feces (/time_scale)
    GE = GEC*BW^(-0.25)	#gastric emptying time (/time_scale)
    k0 = k0C*BW^(-0.25) 	#rate of uptake from the stomach into the liver (/time_scale)

    water_consumption <- 1.36 # L/time_scale

    return(list( "Free" = Free, "QC" = QC, "QK" = QK, "QL" = QL, "QR" = QR,
                 "VPlas" = VPlas,
                 "VKb" = VKb, "Vfil" = Vfil, "VL" = VL, "VR" = VR, "ML" = ML,
                 "VPTC" = VPTC, "Vmax_baso" = Vmax_baso, "Vmax_apical" = Vmax_apical,
                 'kdif' = kdif, "Km_baso" = Km_baso, "Km_apical" = Km_apical,
                 "kbile" = kbile, "kurine" = kurine, "kefflux" = kefflux,
                 "GFR" = GFR, "kabs" = kabs, "kunabs" = kunabs, "GE" = GE, "k0" = k0,
                 "PL" = PL, "PK" = PK, "PR" = PR, "kvoid" = kvoid,
                 "ingestion" = ingestion, "ingestion_time" = ingestion_time,
                 "admin_dose" = admin_dose, "admin_time" = admin_time,
                 "exp_type" = exp_type
    ))

  })
}

#===============================================
#2. Function to create initial values for ODEs
#===============================================

create.inits <- function(parameters){
  with( as.list(parameters),{
    "AR" = 0; "Adif" = 0; "A_baso" = 0; "AKb" = 0;
    "ACl" = 0; "Aefflux" = 0;
    "A_apical" = 0; "APTC" = 0; "Afil" = 0;
    "Aurine" = 0; "AST" = 0;
    "AabsST" = 0; "ASI" = 0; "AabsSI" = 0; "Afeces" = 0;
    "AL" = 0; "Abile" = 0; "Aplas_free" = 0;
    "ingestion" = 0;

    return(c("AR" = AR, "Adif" = Adif, "A_baso" = A_baso, "AKb" = AKb,
             "ACl" = ACl, "Aefflux" = Aefflux,
             "A_apical" = A_apical, "APTC" = APTC, "Afil" = Afil,
             "Aurine" = Aurine, "AST" = AST,
             "AabsST" = AabsST, "ASI" = ASI, "AabsSI" = AabsSI, "Afeces" = Afeces,
             "AL" = AL, "Abile" = Abile, "Aplas_free" = Aplas_free,
             "ingestion" = ingestion))
  })
}

#===================
#3. Events function
#===================
create.events <- function(parameters){
  with(as.list(parameters), {

    if (exp_type == "continuous") {
      # Continuous exposure: set ingestion rates
      lingest <- length(ingestion)
      lingesttimes <- length(ingestion_time)

      # Check that ingestion and ingestion_time have equal length
      if (lingest != lingesttimes){
        stop("The times of ingestion rate change should be equal to the ingestion time vector")
      }

      # Set the ingestion rate at specified times
      events <- list(data = data.frame(
        var = "ingestion",
        time = ingestion_time,
        value = ingestion,
        method = "rep"
      ))

    } else if (exp_type == "bolus") {
      # Bolus dosing: add admin_dose to stomach (AST) at admin_time
      # This adds the specified dose directly to the stomach compartment
      events <- list(data = data.frame(
        var = "AST",
        time = admin_time,
        value = admin_dose,
        method = "add"
      ))

    } else {
      stop(paste("Unknown exp_type:", exp_type, ". Must be 'continuous' or 'bolus'."))
    }

    return(events)
  })
}

#==================
#4. Custom function
#==================
custom.func <- function(){
  return()
}

#==============
#5. ODEs System
#==============

ode.func <- function(time, inits, params, custom.func){
  with(as.list(c(inits,params)),{

    CR = AR/VR #concentration in rest of body (ug/L)
    CVR = CR/PR	#concentration in venous blood leaving the rest of the body (ug/L)
    CKb = AKb/VKb	#concentration in kidney blood (ug/L)
    CVK = CKb #/PK	#concentration in venous blood leaving kidney (ug/L)
    CPTC = APTC/VPTC	#concentration in PTC (ug/L)
    Cfil = Afil/Vfil	#concentration in filtrate (ug/L)
    CL = AL/VL	#concentration in the liver (ug/L)
    CLiver = AL/ML #	concentration in the liver (ug/g)
    CVL = CL/PL	#concentration in the venous blood leaving the liver (ug/L)
    CA_free = Aplas_free/VPlas		#concentration in plasma (ug/L)
    CA = CA_free/Free	#concentration of total PFBS in plasma (ug/L)
    Curine = Aurine/kvoid

    # Rest of Body (Tis)
    dAR = QR*(CA-CVR)*Free	#rate of change in rest of body (ug/time_scale)

    #Kidney
    #Kidney Blood (Kb)
    dAdif = kdif*(CKb - CPTC)	#rate of diffusion from into the PTC (ug/time_scale)
    dA_baso = (Vmax_baso*CKb)/(Km_baso + CKb)
    dAKb = QK*(CA-CVK)*Free - CA*GFR*Free - dAdif - dA_baso #rate of change in kidney blood (ug/time_scale).
    dACl = CA*GFR*Free	#rate of clearance via glomerular filtration (ug/time_scale)

    #Proximal Tubule Cells (PTC)
    dAefflux = kefflux*APTC
    dA_apical = (Vmax_apical*Cfil)/(Km_apical + Cfil)
    dAPTC =  dAdif + dA_apical + dA_baso - dAefflux #rate of change in PTC(ug/time_scale)

    #Filtrate (Fil)
    dAfil = CA*GFR*Free - dA_apical - Afil*kurine	#rate of change in filtrate (ug/time_scale)

    #Urinary elimination
    dAurine = kurine*Afil	#rate of change in urine (ug/time_scale)

    #GI Tract (Absorption site of oral dose)
    #Stomach
    dAST=  ingestion - k0*AST - GE*AST	#rate of change in the stomach (ug/time_scale)
    dAabsST = k0*AST	#rate of absorption in the stomach (ug/time_scale)

    #Small Intestine
    dASI = GE*AST - kabs*ASI - kunabs*ASI	#rate of change in the small intestine (ug/time_scale)
    dAabsSI = kabs*ASI	#rate of absorption in the small intestine (ug/time_scale)

    total_oral_uptake = AabsSI + AabsST	#total oral uptake in the GI tract (ug)

    #Feces compartment
    dAfeces = kbile*AL + kunabs*ASI #rate of change in the feces compartment (ug/time_scale)

    #Liver
    dAL = QL*(CA-CVL)*Free - kbile*AL + kabs*ASI + k0*AST #rate of change in the liver (ug/time_scale)
    dAbile = kbile*AL
    amount_per_gram_liver = CLiver	#amount of PFBS in liver per gram liver (ug/g)

    #Plasma compartment
    dAplas_free = (QR*CVR*Free) + (QK*CVK*Free) + (QL*CVL*Free) -
      (QC*CA*Free) + dAefflux  #rate of change in the plasma (ug/time_scale)

    dingestion = 0

    #Mass Balance Check
    Atissue = Aplas_free + AR + AKb + Afil + APTC + AL + AST + ASI 	#sum of mass in all compartments (ug)
    Aloss = Aurine + Afeces #sum of mass lost through urinary and fecal excretion (ug)
    Atotal = Atissue + Aloss 	#total mass; should equal total dose

    list(c("dAR" = dAR, "dAdif" = dAdif, "dA_baso" = dA_baso, "dAKb" = dAKb,
           "dACl" = dACl, "dAefflux" = dAefflux,
           "dA_apical" = dA_apical, "dAPTC" = dAPTC, "dAfil" = dAfil,
           "dAurine" = dAurine, "dAST" =dAST,
           "dAabsST" = dAabsST, "dASI" = dASI, "dAabsSI" = dAabsSI, "dAfeces" = dAfeces,
           "dAL" = dAL, "dAbile" = dAbile, "dAplas_free" = dAplas_free,
           "dingestion" = dingestion),
         "total_oral_uptake" = total_oral_uptake, "amount_per_gram_liver" = amount_per_gram_liver,
         "Atissue" = Atissue, "Aloss" = Aloss, "Atotal" = Atotal, "CR" =CR, "CVR" = CVR, "CKb" = CKb,
         "CVK" = CVK, "CPTC" = CPTC,
         "Cfil" = Cfil, "CL" = CL, "CVL" = CVL, "CA_free" = CA_free,
         "CA" = CA, "Cliver" = CL)

  })
}

#======================================
# Example Usage
#======================================
# Uncomment to test:

user_input <- list(
  'BW' = 70,
  "ingestion" = c(10, 0),  # 10 ug/day for period, then 0 (elimination)
  "ingestion_time" = c(0, 5),  # Exposure for 5 years, then elimination
  "admin_dose" = 0,
  "admin_time" = 0,
  "exp_type" = "continuous",  # "continuous" or "bolus"
  "time_scale" = "years"
)

# For bolus dosing example:
# user_input <- list(
#   'BW' = 70,
#   "ingestion" = 0,
#   "ingestion_time" = 0,
#   "admin_dose" = c(100,50),  # 100 ug bolus dose
#   "admin_time" = c(0, 5),    # at time 0
#   "exp_type" = "bolus",
#   "time_scale" = "days"
# )

# params <- create.params(user_input)
# inits <- create.inits(params)
# events <- create.events(params)

# times <- seq(0, 10, 0.1)
# out <- ode(y = inits, times = times, func = ode.func, parms = params,
#            events = events, method = "lsodes", atol = 1e-8, rtol = 1e-8)

# plot(out[,"time"], out[,"Cserum"], type = 'l', log = 'y',
#      xlab = 'Time (years)', ylab = 'Serum Concentration (ug/L)',
#      main = 'PFBS Population Model - 70 kg Human')


# ====================
# Upload on Jaqpot
# ====================
# Subset of features to be displayed on the user interface
# Concentrations (ug/L)
predicted.feats <- c("CR", "CVR", "CVK", "CPTC", 
                      "Cfil", "CL", "CVL", "CA_free", "CA",
                     # Amounts/Masses (ug)
                     "AR", "AKb", "APTC", "Afil", "AL", "Aplas_free", "AST", "ASI",
                     "Aurine", "Afeces", "Abile",
                     # Mass balance
                     "Atissue", "Aloss", "Atotal")

jaqpotr::deploy.pbpk(user.input = user_input,out.vars = predicted.feats,
create.params = create.params,  create.inits = create.inits,
create.events = create.events, custom.func = custom.func,
envFile = "/Users/vassilis/Documents/GitHub/jaqpotpy/.env")