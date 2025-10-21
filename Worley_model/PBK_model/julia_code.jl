using DifferentialEquations
using Plots

#===============================================#
#1. Function to create parameters for Worley model
#===============================================#

function create_params_worley(user_input; RAF_baso=nothing, RAF_api=nothing)
    BW = user_input.BW
    substance = user_input.substance
    time_scale_str = get(user_input, :time_scale, "days")
    ingestion = user_input.ingestion
    ingestion_time = user_input.ingestion_time
    admin_dose = user_input.admin_dose
    admin_time = user_input.admin_time
    admin_type = user_input.admin_type
    exp_type = user_input.exp_type

    # Time scale conversion (hours as base unit)
    if time_scale_str == "minutes"
        time_scale = 60.0
    elseif time_scale_str == "hours"
        time_scale = 1.0
    elseif time_scale_str == "days"
        time_scale = 1.0/24.0
    elseif time_scale_str == "weeks"
        time_scale = (1.0/24.0)/7.0
    elseif time_scale_str == "months"
        time_scale = (1.0/24.0)/30.0
    elseif time_scale_str == "years"
        time_scale = (1.0/24.0)/365.0
    else
        time_scale = 1.0/24.0  # Default to days
    end

    inv_time_scale = 1.0/time_scale

    # Cardiac Output and Bloodflow (as fraction of cardiac output)
    QCC = 12.5 * inv_time_scale  # cardiac output in L/time_scale/kg^0.75; Brown 1997
    QLC = 0.25  # fraction blood flow to liver; Brown 1997
    QKC = 0.175  # fraction blood flow to kidney; Brown 1997
    Htc = 0.44  # hematocrit for the rat; Davies 1993

    # Tissue Volumes
    VplasC = 0.0428  # fraction vol. of plasma (L/kg BW); Davies 1993
    VLC = 0.026  # fraction vol. of liver (L/kg BW); Brown 1997
    VKC = 0.004  # fraction vol. of kidney (L/kg BW); Brown 1997
    VfilC = 4e-4  # fraction vol. of filtrate (L/kg BW)
    VPTCC = 1.35e-4  # vol. of proximal tubule cells (L/g kidney)

    # Chemical Specific Parameters - substance dependent
    substance_parameters = Dict(
        "PFBS" => (MW = 300.1, Free = 0.0365),
        "PFHxA" => (MW = 314.05, Free = 0.038),
        "PFOA" => (MW = 414.07, Free = 0.02)  # Original default
    )

    if !haskey(substance_parameters, substance)
        error("Unknown substance: $substance. Available: PFBS, PFHxA, PFOA")
    end

    substance_params = substance_parameters[substance]
    MW = substance_params.MW  # Molecular mass (g/mol)
    Free = substance_params.Free  # Free fraction in plasma

    # Kidney Transport Parameters
    Vmax_baso_invitro = 439.2  # Vmax of basolateral transporter (pmol/mg protein/min)
    if substance == "PFOA"
        Km_baso = 20100.0  # Km of basolateral transporter (ug/L)
    elseif substance == "PFBS"
        Km_baso = 10515.0  # Km of basolateral transporter (ug/L)
        # Estimated based on Weaver et al. 2010 data and 
        # Ryu et al. 2024 data 
        # as Km_PFBS = Km_PFOA * Uptake_rate_PFBS / Uptake_rate_PFOA * MW_PFBS
        # where Km_PFOA = 65.7 uM, Uptake_rate_PFBS = 4.40 and Uptake_rate_PFOA = 8.25
    elseif substance == "PFHxA"
        Km_baso = 15859.5  # Km of basolateral transporter (ug/L) - Weaver et al. 2010 (Km = 50.5 uM*MW, where MW=314.05 ug/umol)
    end
    Vmax_apical_invitro = 37400.0  # Vmax of apical transporter (pmol/mg protein/min)
    Km_apical = 77500.0  # Km of apical transporter (ug/L)

    # Use provided RAF values or defaults
    RAFbaso = isnothing(RAF_baso) ? 1.0 : RAF_baso  # relative activity factor, basolateral transporters (male)
    RAFapi = isnothing(RAF_api) ? 0.0007 : RAF_api  # relative activity factor, apical transporters (male)
    protein = 2.0e-6  # amount of protein in proximal tubule cells (mg protein/cell)
    GFRC = 24.19 * inv_time_scale  # glomerular filtration rate (L/time_scale/kg kidney)

    # Partition Coefficients (from rat tissue data, Kudo et al, 2007)
    PL = 1.03  # liver:blood
    PK = 1.17  # kidney:blood
    PR = 0.11  # rest of body:blood

    # Rate constants
    kdif = 0.001 * inv_time_scale  # diffusion rate from proximal tubule cells
    kabsc = 2.12 * inv_time_scale  # rate of absorption from small intestine to liver
    kunabsc = 7.06e-5 * inv_time_scale  # rate of unabsorbed dose to appear in feces
    GEC = 3.5 * inv_time_scale  # gastric emptying time
    k0C = 1.0 * inv_time_scale  # rate of uptake from stomach into liver
    keffluxc = 0.1 * inv_time_scale  # rate of clearance from PTC into blood
    kbilec = 0.0001 * inv_time_scale  # biliary elimination rate
    kurinec = 0.063 * inv_time_scale  # rate of urine elimination from urine storage
    kvoid = 0.06974 * inv_time_scale  # daily urine volume rate (L/time_scale)

    # Scaled Parameters
    # Cardiac output and blood flows
    QC = QCC * (BW^0.75) * (1 - Htc)  # cardiac output in L/time_scale; adjusted for plasma
    QK = QKC * QC  # plasma flow to kidney (L/time_scale)
    QL = QLC * QC  # plasma flow to liver (L/time_scale)
    QR = QC - QK - QL  # plasma flow to rest of body (L/time_scale)
    QBal = QC - (QK + QL + QR)  # Balance check of blood flows; should equal zero

    # Tissue Volumes
    VPlas = VplasC * BW  # volume of plasma (L)
    VK = VKC * BW  # volume of kidney (L)
    MK = VK * 1.0 * 1000  # mass of the kidney (g)
    VKb = VK * 0.16  # volume of blood in the kidney (L)
    Vfil = VfilC * BW  # volume of filtrate (L)
    VL = VLC * BW  # volume of liver (L)
    ML = VL * 1.05 * 1000  # mass of the liver (g)

    # Kidney Parameters
    PTC = VKC * 1000 * 6e7  # number of PTC (cells/kg BW)
    VPTC = VK * 1000 * VPTCC  # volume of proximal tubule cells (L)
    MPTC = VPTC * 1000  # mass of the proximal tubule cells (g)
    VR = (0.93 * BW) - VPlas - VPTC - Vfil - VL  # volume of remaining tissue (L)
    VBal = (0.93 * BW) - (VR + VL + VPTC + Vfil + VPlas)  # Balance check

    Vmax_basoC = (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW/1e12) * 1e6) * inv_time_scale
    Vmax_apicalC = (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW/1e12) * 1e6) * inv_time_scale
    Vmax_baso = Vmax_basoC * BW^0.75  # (ug/time_scale)
    Vmax_apical = Vmax_apicalC * BW^0.75  # (ug/time_scale)
    kbile = kbilec * BW^(-0.25)  # biliary elimination
    kurine = kurinec * BW^(-0.25)  # urinary elimination
    kefflux = keffluxc * BW^(-0.25)  # efflux clearance rate
    GFR = GFRC * VK  # glomerular filtration rate

    # GI Tract Parameters
    kabs = kabsc * BW^(-0.25)  # rate of absorption from small intestine to liver
    kunabs = kunabsc * BW^(-0.25)  # rate of unabsorbed dose to appear in feces
    GE = GEC * BW^(-0.25)  # gastric emptying time
    k0 = k0C * BW^(-0.25)  # rate of uptake from stomach into liver

    water_consumption = 1.36  # L/time_scale

    return (Free=Free, QC=QC, QK=QK, QL=QL, QR=QR,
            VPlas=VPlas, VKb=VKb, Vfil=Vfil, VL=VL, VR=VR, ML=ML,
            VPTC=VPTC, Vmax_baso=Vmax_baso, Vmax_apical=Vmax_apical,
            kdif=kdif, Km_baso=Km_baso, Km_apical=Km_apical,
            kbile=kbile, kurine=kurine, kefflux=kefflux,
            GFR=GFR, kabs=kabs, kunabs=kunabs, GE=GE, k0=k0,
            PL=PL, PK=PK, PR=PR, kvoid=kvoid,
            ingestion=ingestion, ingestion_time=ingestion_time,
            admin_dose=admin_dose, admin_time=admin_time,
            admin_type=admin_type, exp_type=exp_type)
end

#===============================================#
#2. Function to create initial values for ODEs
#===============================================#

function create_inits_worley(parameters)
    # AR, Adif, A_baso, AKb, ACl, Aefflux, A_apical, APTC, Afil, Aurine,
    # AST, AabsST, ASI, AabsSI, Afeces, AL, Abile, Aplas_free, ingestion
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
end

#===============================================#
#3. Events function
#===============================================#

function create_events_worley(parameters)
    ingestion = parameters.ingestion
    ingestion_time = parameters.ingestion_time

    lingest = length(ingestion)
    lingesttimes = length(ingestion_time)

    # If not equal, then stop
    if lingest != lingesttimes
        error("The times of ingestion rate change should be equal to the ingestion time vector")
    end

    events = []

    # For continuous exposure studies, set the ingestion rate (index 19)
    for (rate, time) in zip(ingestion, ingestion_time)
        affect!(integrator) = integrator.u[19] = rate
        push!(events, PresetTimeCallback(time, affect!))
    end

    # Return empty callback set if no events, otherwise combine all events
    return length(events) > 0 ? CallbackSet(events...) : nothing
end

#===============================================#
#4. ODEs System
#===============================================#

function ode_func_worley(du, u, p, t)
    # Extract state variables
    AR = u[1]          # Rest of body
    Adif = u[2]        # Diffusion (cumulative)
    A_baso = u[3]      # Basolateral transport (cumulative)
    AKb = u[4]         # Kidney blood
    ACl = u[5]         # Clearance (cumulative)
    Aefflux = u[6]     # Efflux (cumulative)
    A_apical = u[7]    # Apical transport (cumulative)
    APTC = u[8]        # Proximal tubule cells
    Afil = u[9]        # Filtrate
    Aurine = u[10]     # Urine
    AST = u[11]        # Stomach
    AabsST = u[12]     # Absorption from stomach (cumulative)
    ASI = u[13]        # Small intestine
    AabsSI = u[14]     # Absorption from small intestine (cumulative)
    Afeces = u[15]     # Feces
    AL = u[16]         # Liver
    Abile = u[17]      # Bile (cumulative)
    Aplas_free = u[18] # Plasma (free)
    ingestion = u[19]  # Ingestion rate

    # Extract parameters
    Free = p.Free
    QC = p.QC
    QK = p.QK
    QL = p.QL
    QR = p.QR
    VPlas = p.VPlas
    VKb = p.VKb
    Vfil = p.Vfil
    VL = p.VL
    VR = p.VR
    ML = p.ML
    VPTC = p.VPTC
    Vmax_baso = p.Vmax_baso
    Vmax_apical = p.Vmax_apical
    kdif = p.kdif
    Km_baso = p.Km_baso
    Km_apical = p.Km_apical
    kbile = p.kbile
    kurine = p.kurine
    kefflux = p.kefflux
    GFR = p.GFR
    kabs = p.kabs
    kunabs = p.kunabs
    GE = p.GE
    k0 = p.k0
    PL = p.PL
    PK = p.PK
    PR = p.PR
    kvoid = p.kvoid

    # Calculate concentrations
    CR = AR / VR                    # concentration in rest of body (ug/L)
    CVR = CR / PR                   # concentration in venous blood leaving rest of body (ug/L)
    CKb = AKb / VKb                 # concentration in kidney blood (ug/L)
    CVK = CKb                       # concentration in venous blood leaving kidney (ug/L)
    CPTC = APTC / VPTC              # concentration in PTC (ug/L)
    Cfil = Afil / Vfil              # concentration in filtrate (ug/L)
    CL = AL / VL                    # concentration in the liver (ug/L)
    CLiver = AL / ML                # concentration in the liver (ug/g)
    CVL = CL / PL                   # concentration in venous blood leaving liver (ug/L)
    CA_free = Aplas_free / VPlas    # concentration in plasma (ug/L)
    CA = CA_free / Free             # concentration of total PFOA in plasma (ug/L)
    Curine = Aurine / kvoid

    # Rest of Body
    du[1] = QR * (CA - CVR) * Free  # dAR

    # Kidney
    # Kidney Blood (Kb)
    dAdif_rate = kdif * (CKb - CPTC)                    # rate of diffusion into PTC
    dA_baso_rate = (Vmax_baso * CKb) / (Km_baso + CKb)  # basolateral transport rate
    du[2] = dAdif_rate                                  # dAdif (cumulative)
    du[3] = dA_baso_rate                                # dA_baso (cumulative)
    du[4] = QK * (CA - CVK) * Free - CA * GFR * Free - dAdif_rate - dA_baso_rate  # dAKb
    du[5] = CA * GFR * Free                             # dACl (cumulative clearance)

    # Proximal Tubule Cells (PTC)
    dAefflux_rate = kefflux * APTC
    dA_apical_rate = (Vmax_apical * Cfil) / (Km_apical + Cfil)
    du[6] = dAefflux_rate                               # dAefflux (cumulative)
    du[7] = dA_apical_rate                              # dA_apical (cumulative)
    du[8] = dAdif_rate + dA_apical_rate + dA_baso_rate - dAefflux_rate  # dAPTC

    # Filtrate
    du[9] = CA * GFR * Free - dA_apical_rate - Afil * kurine  # dAfil

    # Urinary elimination
    du[10] = kurine * Afil  # dAurine

    # GI Tract
    # Stomach
    du[11] = ingestion - k0 * AST - GE * AST  # dAST
    du[12] = k0 * AST                          # dAabsST (cumulative)

    # Small Intestine
    du[13] = GE * AST - kabs * ASI - kunabs * ASI  # dASI
    du[14] = kabs * ASI                             # dAabsSI (cumulative)

    # Feces compartment
    du[15] = kbile * AL + kunabs * ASI  # dAfeces

    # Liver
    du[16] = QL * (CA - CVL) * Free - kbile * AL + kabs * ASI + k0 * AST  # dAL
    du[17] = kbile * AL                                                     # dAbile (cumulative)

    # Plasma compartment
    du[18] = (QR * CVR * Free) + (QK * CVK * Free) + (QL * CVL * Free) -
             (QC * CA * Free) + dAefflux_rate  # dAplas_free

    # Ingestion rate (controlled by events)
    du[19] = 0.0  # dingestion
end

#===============================================#
#5. Post-processing function to extract outputs
#===============================================#

function extract_outputs_worley(sol, parameters)
    VPlas = parameters.VPlas
    VR = parameters.VR
    VKb = parameters.VKb
    VPTC = parameters.VPTC
    Vfil = parameters.Vfil
    VL = parameters.VL
    ML = parameters.ML
    PR = parameters.PR
    PL = parameters.PL
    Free = parameters.Free
    kvoid = parameters.kvoid

    outputs = Dict()

    # Calculate concentrations and metrics for each timepoint
    for (i, u) in enumerate(sol.u)
        AR = u[1]
        AKb = u[4]
        APTC = u[8]
        Afil = u[9]
        Aurine = u[10]
        AST = u[11]
        ASI = u[13]
        AabsSI = u[14]
        AabsST = u[12]
        Afeces = u[15]
        AL = u[16]
        Aplas_free = u[18]

        CR = AR / VR
        CVR = CR / PR
        CKb = AKb / VKb
        CVK = CKb
        CPTC = APTC / VPTC
        Cfil = Afil / Vfil
        CL = AL / VL
        CLiver = AL / ML
        CVL = CL / PL
        CA_free = Aplas_free / VPlas
        CA = CA_free / Free
        Curine = Aurine / kvoid

        # Mass balance
        Atissue = Aplas_free + AR + AKb + Afil + APTC + AL + AST + ASI
        Aloss = Aurine + Afeces
        Atotal = Atissue + Aloss
        total_oral_uptake = AabsSI + AabsST

        # Store outputs
        if i == 1
            outputs["CR"] = Float64[]
            outputs["CVR"] = Float64[]
            outputs["CKb"] = Float64[]
            outputs["CVK"] = Float64[]
            outputs["CPTC"] = Float64[]
            outputs["Cfil"] = Float64[]
            outputs["CL"] = Float64[]
            outputs["CLiver"] = Float64[]
            outputs["CVL"] = Float64[]
            outputs["CA_free"] = Float64[]
            outputs["CA"] = Float64[]
            outputs["Cserum"] = Float64[]
            outputs["Cliver"] = Float64[]
            outputs["Curine"] = Float64[]
            outputs["Atissue"] = Float64[]
            outputs["Aloss"] = Float64[]
            outputs["Atotal"] = Float64[]
            outputs["total_oral_uptake"] = Float64[]
        end

        push!(outputs["CR"], CR)
        push!(outputs["CVR"], CVR)
        push!(outputs["CKb"], CKb)
        push!(outputs["CVK"], CVK)
        push!(outputs["CPTC"], CPTC)
        push!(outputs["Cfil"], Cfil)
        push!(outputs["CL"], CL)
        push!(outputs["CLiver"], CLiver)
        push!(outputs["CVL"], CVL)
        push!(outputs["CA_free"], CA_free)
        push!(outputs["CA"], CA)
        push!(outputs["Cserum"], CA)
        push!(outputs["Cliver"], CL)
        push!(outputs["Curine"], Curine)
        push!(outputs["Atissue"], Atissue)
        push!(outputs["Aloss"], Aloss)
        push!(outputs["Atotal"], Atotal)
        push!(outputs["total_oral_uptake"], total_oral_uptake)
    end

    return outputs
end

#===============================================#
#6. Example simulation (matching R code)
#===============================================#

# User input
# user_input = (
#     BW = 70.0,
#     substance = "PFBS",  # Options: "PFBS", "PFHxA", "PFOA"
#     ingestion = [10.0, 0.0, 12.0],
#     ingestion_time = [0.1, 15.0, 20.0],
#     admin_dose = 0.0,
#     admin_time = 0.0,
#     admin_type = "oral",
#     exp_type = "continuous",
#     time_scale = "years"
# )

# # Create parameters
# params = create_params_worley(user_input)

# # Create initial conditions
# inits = create_inits_worley(params)

# # Create events
# events = create_events_worley(params)

# # Time span
# times = 0.0:1.0:30.0

# # Solve ODE system
# prob = ODEProblem(ode_func_worley, inits, (times[1], times[end]), params)
# sol = solve(prob, Rodas5P(), callback=events, saveat=times,
#             abstol=1e-8, reltol=1e-8)

# # Extract outputs
# outputs = extract_outputs_worley(sol, params)
# println("Last value of Cserum: ", outputs["Cserum"][end])
# # Plot results (matching R plot)
# plot(sol.t, outputs["Cserum"],
#      xlabel="Time (years)",
#      ylabel="Serum Concentration (ug/L)",
#      title="70 kg Human, 10 ug/day for 15 years",
#      linewidth=2,
#      legend=false,
#      ylims=())

# # Display plot
# display(plot!())

# # Print summary statistics
# println("\n=== Simulation Summary ===")
# println("Final serum concentration: ", round(outputs["Cserum"][end], digits=4), " ug/L")
# println("Final liver concentration: ", round(outputs["CLiver"][end], digits=4), " ug/L")
# println("Final urine concentration: ", round(outputs["Curine"][end], digits=4), " ug/L")
# println("Total mass in tissue: ", round(outputs["Atissue"][end], digits=2), " ug")
# println("Total mass lost: ", round(outputs["Aloss"][end], digits=2), " ug")
# println("Total mass (balance): ", round(outputs["Atotal"][end], digits=2), " ug")
