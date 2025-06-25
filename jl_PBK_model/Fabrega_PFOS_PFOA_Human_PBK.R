library(deSolve)
create.params.Fabrega <- function(user_input){
  with(as.list(user_input),{
    # Physiological parameters (from Brown, et al)
    # fractional bloo flows to tissues
    QCC = 12.5		  # Cardiac blood output (L/h/kg^0.75)
    QFC = 0.052		  # Fraction cardiac output going to fat
    QLC = 0.19 #0.25 from Loccisano 2011 		  # Fraction cariac output going to liver
    QKC = 0.175		  # Fraction cardiac output going to kiney
    #QFilC = 0.035		# Fraction cardiac output to the filtrate compartment (20% of kiney blood flow)
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
    QC = QCC*(BW^0.75)*24	#Cardiac output (L/day)
    Htc = 0.44      #hematocrit
    QCP = QC*(1-Htc)	# Plasma flow
    QL = QLC*QCP			# Plasma flow to liver (L/day)
    QF = QFC*QCP			# Plasma flow to fat (L/day)
    QK = QKC*QCP	    # Plasma flow to kiney (L/day)
    QFil = 0.2*QK   	# Plasma flow to filtrate compartment (L/day)# 20% of QK
    QG = QGC*QCP	    # Plasma flow to gut (L/day)
    QLu = QLuC*QCP	  # Plasma flow to lungs (L/day)
    QB = QBC*QCP	    # Plasma flow to brain (L/day)
    QR = QCP - QL - QF - QK - QFil - QG - QLu - QB	 # Plasma flow to rest of the boy (L/day)
    
    Qbal = QCP - (QL+QF+QK+QFil+QG+QR+QLu+QB)        # balance check 
    
    VL = VLC*BW			    # Liver volume (L)
    VF = VFC*BW			    # Fat volume (L)
    VK = VKC*BW			    # Kiney volume (L)
    VFil = VFilC*BW	    # Fitrate compartment volume (L)
    VG = VGC*BW			    # Gut volume (L)
    VPlas = VPlasC*BW		# Plasma volume (L)
    VLu = VLuC*BW			  # Lungs volume (L)
    VB = VBC*BW			    # Brain volume (L)
    VR = 1*BW - VL - VF - VK - VFil - VG -VLu -VB - VPlas		# Rest of the boy volume (L) # Loccisano 2011 uses 0.84*BW
    
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
      Tm = 245.6 # PFHxA does not have a Tmax, so set to 0
      Kt = 0.6 # Estimated value for PFHxA
      Free = 0.01 # Fabrega (2014) Table 1
      PL = 0.001 # Partition Coefficient for Liver: Estimated
      PF = 0.467 # Partition Coefficient for Fat: Estimated
      PB = 43.68 # Partition Coefficient for Brain: Estimated
      PLu = 31.92 # Partition Coefficient for Lungs: Estimated
      PK = 11.57 # Partition Coefficient for Kidney: Estimated
      PG = 00.05 # Partition Coefficient for Gut: Estimated
      PR = 0.12 # Partition Coefficient for Rest of body: Estimated
      kurinec = 3e-04	#urinary elimination rate constant (/h/kg^-0.25); estimated from Harada, 
      kurine = kurinec*(BW^-0.25)*24 # Elimination rate (1/day
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
    
    # Concentrations (ng/g)
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
    
    # Plasma compartment (ng)
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

# Test with same parameters as Julia model
BW <- 70
substance <- "PFHxA"  # Note: Need to update function to handle PFHxA
f_unabs <- 0.5
admin_dose <- 0
admin_time <- 0
ingestion <- 100
ingestion_time <- 0.01  # Convert to days (R model uses days)
admin_type <- "oral"
exp_type <- "biomonitoring"

user_input <- list(
  BW = BW,
  substance = substance,
  f_unabs = f_unabs,
  admin_dose = admin_dose,
  admin_time = admin_time,
  ingestion = ingestion,
  ingestion_time = ingestion_time,
  admin_type = admin_type,
  exp_type = exp_type
)

# Create parameters, initial conditions, and events
params <- create.params.Fabrega(user_input)
inits <- create.inits.Fabrega(params)
events <- create.events.Fabrega(params)

# Solve ODE
library(deSolve)
times <- seq(0, 10, by = 0.2)  # 10 days with 0.2 day intervals
solution <- ode(y = inits, times = times, func = ode.func.Fabrega, 
                parms = params, events = events, method = "lsoda")

# Convert to dataframe and show tail
solution_df <- as.data.frame(solution)
print("Tail of R solution:")
print(tail(solution_df))
