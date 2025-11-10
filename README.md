# Insight-PBK

Repository of physiologically-based pharmacokinetic (PBK) models for various chemical compounds developed for the Insight EU Project. This collection includes PBK models for acetone, toluene diisocyanate (TDI), ethylene glycol (EG), and PFAS compounds (PFBS and PFHxA).

## Available Models

### 1. Acetone PBK Model (Mork and Johanson, 2006)

Physiologically-based pharmacokinetic model for acetone exposure based on the Mork and Johanson (2006) framework. The model describes acetone distribution and elimination following inhalation or oral exposure.

**Location**: `Mork and Johanson,2006/`

**Implementation**: R (deSolve)

**Deployment**: Includes Jaqpot deployment version (`Acetone_PBK_model_Jaqpot.R`)

### 2. TDI PBK Model (Scholten et al. 2023)

PBK model for toluene diisocyanate (TDI) exposure based on Scholten et al. (2023). The model captures TDI kinetics and metabolism following inhalation exposure.

**Location**: `Scholten et al.2023/`

**Implementation**: R (deSolve)

**Deployment**: Includes Jaqpot deployment version (`Scholten_TDI_Model_Jaqpot upload.R`)

### 3. Ethylene Glycol PBK Model (Corley et al. 2005)

PBK model for ethylene glycol (EG) based on Corley et al. (2005). The model describes EG distribution, metabolism, and elimination following oral exposure.

**Location**: `Corley et al.2005/`

**Implementation**: R (deSolve)

**Deployment**: Includes Jaqpot deployment version (`Corley_et_al_2005_model_jaqpot upload.R`)

### 4. PFAS PBK Models (PFBS and PFHxA)

Physiologically-based pharmacokinetic models for PFAS compounds implementing the Worley model with active renal reabsorption. These models use Bayesian inference to estimate transporter activity parameters from biomonitoring data.

**Location**: `PFAS_models/`

**Implementation**: Julia (DifferentialEquations.jl, Turing.jl) with R deployment versions

**Compounds**:
- **PFBS**: Hierarchical Bayesian model for population-level analysis (6 subjects)
- **PFHxA**: Non-hierarchical Bayesian model for single-subject analysis

**Key Features**:
- Active renal reabsorption mechanisms (basolateral and apical transporters)
- MCMC parameter estimation using NUTS sampler
- Posterior predictive checks for model validation
- Separate observation errors for plasma and urine data (PFHxA)

**Deployment**: Both models include Jaqpot deployment versions in their respective `jaqpot_deployment/` directories

### 5. Partition Coefficient Calculator (EDM)

Implementation of the Allendorf Equilibrium Distribution Model (EDM) for calculating tissue:blood partition coefficients from tissue composition data.

**Location**: `EDM_Partition_Coefficients/`

**Implementation**: Julia

**Purpose**: Generates partition coefficients for liver, kidney, and rest of body compartments used in PBK models

## Project Structure

```
Insight-PBK/
├── Mork and Johanson,2006/      # Acetone PBK model
├── Scholten et al.2023/          # TDI PBK model
├── Corley et al.2005/            # Ethylene glycol PBK model
├── PFAS_models/                  # PFAS PBK models
│   ├── PFBS/                     # PFBS hierarchical model
│   │   └── jaqpot_deployment/    # R code for Jaqpot
│   └── PFHxA/                    # PFHxA non-hierarchical model
│       └── jaqpot_deployment/    # R code for Jaqpot
└── EDM_Partition_Coefficients/   # Partition coefficient calculations
```
