library(deSolve)

create.params.worley  <- function(user_input){
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
    MW = 414.07	#PFOA molecular mass (g/mol)
    
    Free = 0.02#0.006 #free fraction in plasma (male)
    
    #Kidney Transport Parameters
    Vmax_baso_invitro = 439.2 #Vmax of basolateral transporter (pmol/mg protein/min); averaged in vitro value of OAT1 and OAT3 from Nakagawa, 2007
    Km_baso = 20100 #Km of basolateral transporter (ug/L) Average of OAT1 and OAT3 from Nakagawa et. al, 2007
    Vmax_apical_invitro = 37400 #Vmax of apical transporter (pmol/mg protein/min); invitro value for OAT4 from Yang et al, 2010
    Km_apical = 77500#Km of apical transporter (ug/L), in vitro value for OAT4 and URAT1 from Yang et al, 2010.
    RAFbaso = 1	#relative activity factor, basolateral transporters (male) (fit to data)	
    RAFapi = 0.0007	#relative activity factor, apical transporters (male) (fit to data)	
    protein = 2.0e-6	#amount of protein in proximal tubule cells (mg protein/proximal tubule cell)
    GFRC = 24.19*inv_time_scale	#glomerular filtration rate (L/time_scale/kg kidney) (male); Corley, 2005
    
    #Partition Coefficients (from rat tissue data, Kudo et al, 2007)
    PL = 1.03 #liver:blood
    PK = 1.17 #kidney:blood
    PR = 0.11 #rest of body:blood
    
    #rate constants
    kdif = 0.001*inv_time_scale	#diffusion rate from proximal tubule cells (L/time_scale)
    kabsc = 2.12*inv_time_scale	#rate of absorption of chemical from small intestine to liver (1/(time_scale*BW^-0.25))(fit to data)
    kunabsc = 7.06e-5*inv_time_scale	#rate of unabsorbed dose to appear in feces (1/(time_scale*BW^-0.25))(fit to data)
    GEC = 3.5*inv_time_scale#gastric emptying time (1/(time_scale*BW^-0.25)); from Yang, 2013
    k0C = 1.0*inv_time_scale	#rate of uptake from the stomach into the liver (1/(time_scale*BW^-0.25)) (fit to data)
    
    keffluxc =0.1*inv_time_scale #rate of clearance of PFOA from proximal tubule cells into blood (1/(time_scale*BW^-0.25))
    kbilec = 0.0001*inv_time_scale #biliary elimination rate ((male); liver to feces storage (1/(time_scale*BW^-0.25)) (fit to data)
    kurinec = 0.063*inv_time_scale #rate of urine elimination from urine storage (male) (1/(time_scale*BW^-0.25))(fit to data)
    kvoid = 0.06974*inv_time_scale  #daily urine volume rate (L/time_scale); Van Haarst, 2004                                                   
    
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
    
    Vmax_basoC = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1e6)*inv_time_scale#Vmax of basolateral transporters (ug/time_scale/kg BW)
    Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1e6)*inv_time_scale #Vmax of basolateral transporters (ug/time_scale/kg BW)
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
    
    water_consumption <- 1.36# L/time_scale
    
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
                 "admin_type" = admin_type, "exp_type" = exp_type
    ))
    
  })
}

#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits.worley <- function(parameters){
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
create.events.worley <- function(parameters){
  with(as.list(parameters), {
   
      lingest <- length(ingestion)
      lingesttimes <- length(ingestion_time)
      # If not equal, then stop 
      if (lingest != lingesttimes){
        stop("The times of ingestion rate change should be equal to the ingestion time vector")
      }else{
       
          # For continuous exposure studies, set the ingestion rate
          events <- list(data = data.frame(var = c("ingestion"),  time = ingestion_time, 
                                           value = ingestion, method = c("rep")))
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

ode.func.worley <- function(time, inits, params, custom.func){
  with(as.list(c(inits,params)),{
    
    CR = AR/VR #concentration in rest of body (ug/L)
    CVR = CR/PR	#concentration in venous blood leaving the rest of the body (ug/L)
    CKb = AKb/VKb	#concentration in kidney blodd (ug/L) 
    CVK = CKb #/PK	#concentration in venous blood leaving kidney (ug/L)
    CPTC = APTC/VPTC	#concentration in PTC (ug/L)
    Cfil = Afil/Vfil	#concentration in filtrate (ug/L)
    CL = AL/VL	#concentration in the liver (ug/L)
    CLiver = AL/ML #	concentration in the liver (ug/g)
    CVL = CL/PL	#concentration in the venous blood leaving the liver (ug/L)
    CA_free = Aplas_free/VPlas		#concentration in plasma (ug/L)
    CA = CA_free/Free	#concentration of total PFOA in plasma (ug/L)
    Curine = Aurine/kvoid
    
    # Rest of Body (Tis)
    dAR = QR*(CA-CVR)*Free	#rate of change in rest of body (ug/time_scale)
    
    #Kidney 
    #Kidney Blood (Kb)
    dAdif = kdif*(CKb - CPTC)	#rate of diffusion from into the PTC (ug/time_scale)
    dA_baso = (Vmax_baso*CKb)/(Km_baso + CKb)	
    dAKb = QK*(CA-CVK)*Free - CA*GFR*Free - dAdif - dA_baso #rate of change in kidney blood (ug/time_scale).
    dACl = CA*GFR*Free	#rate of clearance via glormerular filtration (ug/time_scale)
    
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
    dAabsSI = kabs*ASI	#rate of absorption inthe small intestine (ug/time_scale)
    
    total_oral_uptake = AabsSI + AabsST	#total oral uptake in the GI tract (ug) 
    
    #Feces compartment
    dAfeces = kbile*AL + kunabs*ASI #rate of change in the feces compartment (ug/time_scale)
    
    #Liver
    dAL = QL*(CA-CVL)*Free - kbile*AL + kabs*ASI + k0*AST #rate of change in the liver (ug/time_scale)
    dAbile = kbile*AL
    amount_per_gram_liver = CLiver	#amount of PFOA in liver per gram liver (ug/g)
    
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
         "Cfil" = Cfil, "CL" = CL, "CVL" = CVL, "CA_free" = CA_free, "CA" = CA,
         "Cserum" = CA, "Cliver" = CL)
    
  })
}

user_input <- list(
  'BW' = 70,
  "exposure_time" = seq(0,30,1),
  "ingestion" = c(10, 0, 12),
  "ingestion_time" = c(0, 15, 20),
  "admin_dose" = 0,
  "admin_time" = 0,
  "admin_type" = "oral",
  "exp_type" = "continuous",
  "time_scale" = "years"
)
params <- create.params.worley(user_input)
inits <- create.inits.worley(params)
events <- create.events.worley(params)

times <- seq(0,30,1)
out <- ode(y = inits, times = times, func = ode.func.worley, parms = params, 
           events = events, method = "lsodes", atol = 1e-8, rtol = 1e-8)

plot(out[,"time"], out[,"Cserum"], type = 'l', log = 'y', xlab = 'Time (years)', ylab = 'Serum Concentration (ug/L)', main = '70 kg Human, 10 ug/day for 15 years')

tail(out)
