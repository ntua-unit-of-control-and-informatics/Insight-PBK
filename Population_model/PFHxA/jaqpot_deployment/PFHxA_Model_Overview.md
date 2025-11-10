# PFHxA Single-Dose PBK Model Overview

## What is this model?

This document describes a **physiologically-based pharmacokinetic (PBK) model** for PFHxA (perfluorohexanoic acid), a PFAS compound. The model predicts how PFHxA moves through the human body over time - how it's absorbed, distributed to different organs, and eliminated through urine and feces.

Think of it like a detailed map of PFHxA's journey through the body, with mathematical equations describing the speed and direction of travel at each step.

---

## Model Structure

### Compartments (Body Parts Tracked)

The model divides the body into several compartments:

1. **Stomach (AST)** - Where oral doses first arrive
2. **Small Intestine (ASI)** - Where absorption into blood occurs
3. **Liver (AL)** - Major metabolic organ
4. **Kidney Blood (AKb)** - Blood circulating through kidneys
5. **Proximal Tubule Cells (APTC)** - Kidney cells that actively reabsorb PFHxA
6. **Filtrate (Afil)** - Fluid filtered by kidneys
7. **Plasma (Aplas_free)** - Blood plasma (free fraction)
8. **Rest of Body (AR)** - All other tissues (heart, muscle, skin, etc.)
9. **Urine (Aurine)** - Cumulative urinary excretion
10. **Feces (Afeces)** - Cumulative fecal excretion

### Key Processes

- **Absorption**: PFHxA moves from the gut into the liver and bloodstream
- **Distribution**: Blood carries PFHxA to all organs based on blood flow rates
- **Kidney Reabsorption**: Special kidney transporters actively pull PFHxA back from urine into blood
  - **Basolateral transporters** (OAT1/OAT3): Move PFHxA from blood into kidney cells
  - **Apical transporters** (OAT4): Move PFHxA from filtrate back into kidney cells
- **Elimination**: PFHxA leaves the body through urine (main route) and feces (minor route)

---

## The Data

### Study Design

We used data from a **controlled human study** (Abraham et al. 2024) where a single volunteer received an oral dose of PFHxA and had both plasma and urine monitored over time.

- **Subject**: 1 healthy adult male
- **Dose**: 3.99 μg single oral dose
- **Measurements**:
  - Plasma concentration time-series
  - Cumulative urine excretion time-series
- **Time Range**: Several days post-dosing
- **Study Type**: Single-dose pharmacokinetic study
- **Background Exposure**: Low continuous exposure (0.93 ng/kg-bw/day from Haug et al. 2009)

### Why Plasma and Urine Data?

By using both plasma concentration and cumulative urine excretion data simultaneously, we can:
- Better constrain kidney transporter parameters (both basolateral and apical)
- Distinguish between different elimination pathways
- Improve parameter identifiability compared to using plasma data alone
- Validate that the model correctly predicts both blood levels and urinary elimination

---

## Statistical Model: Bayesian Approach

### What Does "Bayesian" Mean?

The Bayesian approach combines prior knowledge with observed data to estimate model parameters:

**Bayesian framework**:
- We start with prior beliefs about plausible parameter values (based on literature)
- We update these beliefs using the observed plasma and urine data
- We end up with posterior probability distributions for each parameter
- These distributions quantify both the best estimate AND the uncertainty

---

### Mathematical Formulation

#### Estimated Parameters

We estimate 4 key parameters from the single-dose study data:

**1. RAFapi:** Apical transporter activity factor
- Represents how active the kidney OAT4 transporters are at reabsorbing PFHxA from filtrate

**2. RAFbaso:** Basolateral transporter activity factor
- Represents how active the kidney OAT1/OAT3 transporters are at uptaking PFHxA from blood

**3. Km_apical:** Michaelis constant for apical transporters (μg/L)
- Binding affinity of PFHxA to OAT4 transporters
- Lower Km = stronger binding

**4. PC_scale:** Partition coefficient scaling factor
- Adjusts the "rest of body" tissue:blood partition coefficient
- Accounts for uncertainty in tissue distribution

**5. Fixed Parameters**
- Km_baso = 8952.96 μg/L (from in vitro studies, PFHxA-specific)
- Partition coefficients for liver and kidney (from EDM calculations)

---

### Prior Distributions

**Transport Parameters:**

```
RAFapi ~ LogNormal(μ_RAFapi, σ_RAFapi)
  where: μ_RAFapi = log(0.0007) - σ²_RAFapi/2
         σ_RAFapi = sqrt(log(1 + CV²))
         CV = 1.0 (100% coefficient of variation)

RAFbaso ~ LogNormal(μ_RAFbaso, σ_RAFbaso)
  where: μ_RAFbaso = log(0.5) - σ²_RAFbaso/2
         σ_RAFbaso = sqrt(log(1 + 0.5²))
         CV = 0.5 (50% coefficient of variation)

Km_apical ~ LogNormal(μ_Km, σ_Km)
  where: μ_Km = log(52338.44) - σ²_Km/2
         σ_Km = sqrt(log(1 + 0.2²))
         CV = 0.2 (20% coefficient of variation)
```

**Partition Coefficient Scaling:**

```
PC_scale ~ truncated(Normal(0.5, 0.5), 0, ∞)
  (constrained to be positive, centered at 0.5)
```

**Observation Errors:**

```
σ_plasma ~ truncated(Normal(0, 0.2), 0, ∞)
  (log-scale error for plasma concentrations)

σ_urine ~ truncated(Normal(0, 0.2), 0, ∞)
  (log-scale error for cumulative urine excretion)
```

---

### Likelihood (How Data Relates to Model)

The model fits both plasma and urine data simultaneously with separate error terms:

**Plasma Data:**
```
For each plasma observation j:
  log(C_plasma,obs[j]) ~ Normal(log(C_plasma,pred[j]), σ_plasma)
```

**Urine Data:**
```
For each urine observation k:
  log(Aurine,obs[k]) ~ Normal(log(Aurine,pred[k]), σ_urine)
```

Where:
- **C_plasma,obs**: Observed plasma concentration (μg/L)
- **C_plasma,pred**: Model-predicted plasma concentration using estimated parameters
- **Aurine,obs**: Observed cumulative urine excretion (μg)
- **Aurine,pred**: Model-predicted cumulative urine excretion
- **σ_plasma**, **σ_urine**: Separate observation errors for each data type

**Why separate errors?**
- Plasma and urine measurements have different sources of variability
- Analytical methods differ between plasma and urine assays
- This prevents one data type from dominating the fit

**Why log-transform?**
- Both concentrations and cumulative amounts span multiple orders of magnitude
- Errors are proportional to measured values (relative error)
- Log-normal distribution ensures predictions remain positive

---

### Summary of Parameter Structure

| Parameter | Type | Description |
|-----------|------|-------------|
| RAFapi | Estimated | Apical transporter activity factor |
| RAFbaso | Estimated | Basolateral transporter activity factor |
| Km_apical | Estimated | Apical transporter binding affinity (μg/L) |
| PC_scale | Estimated | Partition coefficient scaling factor |
| σ_plasma | Estimated | Plasma observation error (log scale) |
| σ_urine | Estimated | Urine observation error (log scale) |
| Km_baso | Fixed | Basolateral binding affinity = 8952.96 μg/L |

**Total estimated parameters:** 6 parameters

---

### MCMC Sampling

We use **NUTS (No-U-Turn Sampler)**, an advanced Markov Chain Monte Carlo algorithm that efficiently explores the parameter space.

- **4 chains** running in parallel (to check convergence)
- **1000 iterations** per chain
- **500 burn-in** iterations (discarded to ensure convergence)

**Output**:
- Posterior distributions for all parameters
- Uncertainty quantification (95% credible intervals)
- Parameter correlations

---

## Partition Coefficients

### What Are Partition Coefficients?

Partition coefficients (PCs) describe how PFHxA distributes between blood and tissues when at equilibrium.

**Example**: If liver:blood PC = 0.53, then at equilibrium, the concentration in liver is 53% of the blood concentration.

### How Were They Calculated?

We used the **Equilibrium Distribution Model (EDM)** developed by Allendorf et al., which calculates partition coefficients based on:

1. **Tissue composition**:
   - Water content
   - Protein content (albumin, phospholipids, neutral lipids)

2. **Chemical properties of PFHxA**:
   - How strongly it binds to different tissue components
   - Its water solubility and protein binding affinity

### PFHxA Partition Coefficients

| Tissue | Partition Coefficient (Tissue:Blood) |
|--------|-------------------------------------|
| Liver | 0.5279 |
| Kidney | 0.5141 |
| Rest of Body (base) | 0.3991 |
| Rest of Body (scaled) | 0.3991 × PC_scale |

**Interpretation**:
- PFHxA concentrations in tissues are about 40-53% of blood concentrations
- PFHxA stays primarily in blood rather than accumulating in tissues
- This is typical for highly water-soluble, protein-bound compounds
- The PC_scale factor accounts for uncertainty in "rest of body" tissue distribution

### Rest of Body Composition

"Rest of body" is a **mass-weighted average** of multiple tissues:
- Heart
- Gonads (reproductive organs)
- Skin
- Spleen
- Muscle

Each tissue's partition coefficient is weighted by its mass fraction in the body. The scaling factor (PC_scale) is estimated from the data to account for uncertainty in this composite tissue.

---

## Model Calibration Results

### Parameter Estimates

From the Bayesian calibration using plasma and urine data, we obtained:

**Transporter Activity Factors:**
- **RAFapi** = [1.41e-6 ± 3.589e-7] (95% CI: [8.137e-7, 2.195e-6])
  - Apical transporter activity

- **RAFbaso** = [0.02883 ± 0.007836] (95% CI: [0.01624, 0.04663])
  - Basolateral transporter activity

**Transporter Binding Affinity:**
- **Km_apical** = [72750.0 ± 14210.0] μg/L (95% CI: [49000.0, 104100.0])
  - Binding affinity for apical transporters (OAT4)

**Distribution Parameter:**
- **PC_scale** = [0.1226 ± 0.02496] (95% CI: [0.07336, 0.1707])
  - Scaling factor for rest of body partition coefficient

**Observation Errors:**
- **σ_plasma** = [0.3636 ± 0.0398] (95% CI: [0.2953, 0.4481])
  - Plasma measurement variability

- **σ_urine** = [0.5592 ± 0.0725] (95% CI: [0.4308, 0.714])
  - Urine measurement variability

---

### Model Diagnostics

The model was validated using:

1. **Posterior Predictive Checks (PPC)**
   - Simulating new plasma and urine data from the fitted model
   - Comparing predictions to observed data for both endpoints
   - Good fit: observed data falls within 95% prediction interval
   - Separate validation for plasma and urine time-courses

2. **MCMC Diagnostics**
   - **R̂ (R-hat)**: Should be < 1.1 (indicates chain convergence)
   - **Effective Sample Size (ESS)**: Should be > 400 (indicates sufficient sampling)
   - **Trace Plots**: Should show good mixing without trends

3. **Parameter Identifiability**
   - Using both plasma and urine data improves identifiability
   - RAFbaso and RAFapi can be estimated separately (unlike models with plasma-only data)
   - Posterior distributions are well-constrained

---

## Software & Reproducibility

### Implementation

- **Language**: Julia (for Bayesian calibration) and R (for deployment)
- **ODE Solver**: High-precision adaptive solvers (Rodas5P in Julia, lsodes in R)
- **Bayesian Inference**: Turing.jl with NUTS sampler
- **Visualization**: Plots.jl and StatsPlots.jl

### Files

- **MCMC Script**: `PFHxA_MCMC_with_urine.jl`
- **Posterior Predictive Checks**: `PFHxA_ppc_with_urine.jl`
- **Chain Diagnostics**: `view_chains.jl`
- **R Implementation**: `jaqpot_deployment/deploy_on_jaqpot.R`

### Reproducibility

All code and data are version-controlled and documented. The model can be re-run to verify results or adapted for new scenarios.

---

## References

### Model Framework
- **Worley et al.**: Worley RR, Yang X, Fisher J. Physiologically based pharmacokinetic modeling of human exposure to perfluorooctanoic acid suggests historical non drinking-water exposures are important for predicting current serum concentrations. Toxicol Appl Pharmacol. 2017 Sep 1;330:9-21. doi: 10.1016/j.taap.2017.07.001. Epub 2017 Jul 3. PMID: 28684146; PMCID: PMC5664934.

### Data Sources
- **Haug et al. (2009)**: Background PFHxA exposure estimates
- **Abraham et al. (2024)**: Single-dose pharmacokinetic study with plasma and urine data
  - Abraham K, Gerofke A, Razzazi-Fazeli E, Pabel U, Pestemer S, Trapp A, Schafft H. Mass Balance and Plasma and Urine Kinetics of Perfluorohexanoic Acid (PFHxA) in Male Volunteers. Chem Res Toxicol. 2024 Jul 15;37(7):1107-1113. doi: 10.1021/acs.chemrestox.4c00032. Epub 2024 Jun 3. PMID: 38829655; PMCID: PMC11261607.

### Partition Coefficients
- **Allendorf EDM**: Allendorf F, Goss KU, Ulrich N. Estimating the Equilibrium Distribution of Perfluoroalkyl Acids and 4 of Their Alternatives in Mammals. Environ Toxicol Chem. 2021 Mar;40(3):910-920. doi: 10.1002/etc.4954. Epub 2021 Feb 9. PMID: 33289938.
- **Tissue Composition Data**: Literature values for human tissue composition

### Transporter Parameters
- **In Vitro Studies**: OAT1, OAT3, OAT4 transporter kinetics from literature
- **Calibrated Values**: RAFapi, RAFbaso, and Km_apical estimated from human data
