using DifferentialEquations
# using UnPack
using DataFrames
using CSV

# Define the create.params function.
function create_params(user_input)
    BW = user_input.BW
    substance = user_input.substance
    admin_dose = user_input.admin_dose
    admin_time = user_input.admin_time
    f_unabs = user_input.f_unabs
    ingestion = user_input.ingestion
    ingestion_time = user_input.ingestion_time
    admin_type = user_input.admin_type
    exp_type = user_input.exp_type

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
    VBC = 0.021	# Fraction brain volume

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
      "PFBS" => (Liver=128.8, Brain=201.6, Lung=56.11, Kidney=6.27, Tm=6.1, Kt=5, Free=0.0365), # Original value Free=0.001
      "PFHxA" => (Liver=0.001, Brain=43.68, Lung=31.92, Kidney=11.57, Tm=245.6, Kt=0.6, Free=0.038) # Original value Free=0.01
    )
    keep_params = kinetic_parameters[substance]

    Tm = keep_params.Tm*24 #ug/h Fabrega (2015) Table 2
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
            # Tm=Tm, 
            Kt=Kt, Free=Free, kurine=kurine,
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
    admin_dose = parameters.admin_dose
    admin_time = parameters.admin_time
    admin_type = parameters.admin_type
    exp_type = parameters.exp_type
    ingestion = parameters.ingestion
    ingestion_time = parameters.ingestion_time
    Free = parameters.Free
    
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
    # @unpack QC, QCP, QL, QF, QK, QFil, QG, QLu, QB, QR,
    #         VPlas, VL, VF, VK, VFil, VG, VLu, VB, VR,
    #         PL, PF, PB, PLu, PK, PG, PR,
    #         Tm, Kt, Free, kurine = p
    
    # Extract parameters explicitly from p tuple
    QC = p.QC
    QCP = p.QCP
    QL = p.QL
    QF = p.QF
    QK = p.QK
    QFil = p.QFil
    QG = p.QG
    QLu = p.QLu
    QB = p.QB
    QR = p.QR
    VPlas = p.VPlas
    VL = p.VL
    VF = p.VF
    VK = p.VK
    VFil = p.VFil
    VG = p.VG
    VLu = p.VLu
    VB = p.VB
    VR = p.VR
    PL = p.PL
    PF = p.PF
    PB = p.PB
    PLu = p.PLu
    PK = p.PK
    PG = p.PG
    PR = p.PR
    Tm = p.Tm
    Kt = p.Kt
    Free = p.Free
    kurine = p.kurine
    
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
    du[9] = 0 #QFil*CFil - kurine*AStore

    # Urine compartment
    du[10] = QFil*CFil # kurine*AStore

    # Rest of body compartment
    du[11] = QR*Free*(CPlas - CR/PR)

    # Ingestion compartment
    du[12] = 0 # Ingestion is not modeled as a differential equation,
    # but rather as a source term in the gut compartment
    
    # Calculate instantaneous urine concentration
    # Excretion rate (μg/day) / Urine flow rate (L/day) = μg/L = ng/mL
    urine_flow_rate = 1.26  # L/day
    instantaneous_excretion_rate = kurine * AStore  # μg/day (same as du[10])
    urine_concentration_instant = instantaneous_excretion_rate / urine_flow_rate  # ng/mL
    
    # Store this as additional output (can be accessed via callback)
    # This will be available in the solution for post-processing
end

# Post-process to get Concentrations
function extract_concentrations(sol, parameters)
    VPlas = parameters.VPlas
    VG = parameters.VG
    VL = parameters.VL
    VF = parameters.VF
    VLu = parameters.VLu
    VB = parameters.VB
    VK = parameters.VK
    VFil = parameters.VFil
    VR = parameters.VR
    Free = parameters.Free
    kurine = parameters.kurine
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
    
    # Calculate instantaneous urine concentration using physiological approach
    urine_flow_rate = 1.26  # L/day
    urine_concentrations = []
    
    for u in sol.u
        AFiltrate = u[8]  # Storage compartment mass (μg)
        # instantaneous_excretion_rate = kurine * AStore  # μg/day
        instantaneous_excretion_rate = AFiltrate / 0.028 * 250  # μg/day
        urine_conc = instantaneous_excretion_rate / urine_flow_rate  # ng/mL
        push!(urine_concentrations, urine_conc)
    end
    
    concentrations["CUrine"] = urine_concentrations
    return concentrations
end