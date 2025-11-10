# PFBS Population PBK Model Overview

## What is this model?

This document describes a **physiologically-based pharmacokinetic (PBK) model** for PFBS (perfluorobutanesulfonic acid), a PFAS compound. The model predicts how PFBS moves through the human body over time - how it's absorbed, distributed to different organs, and eliminated through urine and feces.

**Model Structure**: The PBK model follows the same compartmental structure as **Worley et al. (2017)**, which was originally developed for PFOA and includes active renal reabsorption mechanisms.

Think of it like a detailed map of PFBS's journey through the body, with mathematical equations describing the speed and direction of travel at each step.

---

## Model Implementation and Deployment

**Training**: The model parameters were estimated using a **hierarchical Bayesian framework** with MCMC sampling (see Statistical Model section below). This approach provides full uncertainty quantification for all parameters.

**Deployment**: The Jaqpot deployment provides **deterministic predictions** using the mean posterior estimates from the Bayesian calibration. This means:
- Predictions are made with fixed parameter values (no uncertainty propagation)
- Fast and efficient for real-time predictions
- Parameter uncertainty is not reflected in the output
- Suitable for point predictions and scenario testing

---

## Model Structure

### Compartments (Body Parts Tracked)

The model divides the body into several compartments:

1. **Stomach (AST)** - Where oral doses first arrive
2. **Small Intestine (ASI)** - Where absorption into blood occurs
3. **Liver (AL)** - Major metabolic organ
4. **Kidney Blood (AKb)** - Blood circulating through kidneys
5. **Proximal Tubule Cells (APTC)** - Kidney cells that actively reabsorb PFBS
6. **Filtrate (Afil)** - Fluid filtered by kidneys
7. **Plasma (Aplas_free)** - Blood plasma (free fraction)
8. **Rest of Body (AR)** - All other tissues (heart, muscle, skin, etc.)
9. **Urine (Aurine)** - Cumulative urinary excretion
10. **Feces (Afeces)** - Cumulative fecal excretion

### Key Processes

- **Absorption**: PFBS moves from the gut into the liver and bloodstream
- **Distribution**: Blood carries PFBS to all organs based on blood flow rates
- **Kidney Reabsorption**: Special kidney transporters actively pull PFBS back from urine into blood
  - **Basolateral transporters** (OAT1/OAT3): Move PFBS from blood into kidney cells
  - **Apical transporters** (OAT4): Move PFBS from filtrate back into kidney cells
- **Elimination**: PFBS leaves the body through urine (main route) and feces (minor route)

---

## The Data

### Study Population

We used data from **6 human volunteers** who are occupationally exposed to PFBS and had their blood concentrations monitored during the **elimination phase** (after leaving from work).

- **Subjects**: 6 healthy adults
- **Exposure**: Occupational
- **Measurements**: Serum concentration over time (elimination phase only)
- **Time Range**: Days to weeks after dosing
- **Background Exposure**: Low continuous exposure (0.93 ng/kg-bw/day from Haug et al. 2009)

### Why Elimination Phase Data?

By focusing on the elimination phase (when PFBS is leaving the body), we can better estimate the kidney's ability to reabsorb PFBS - the key process that makes PFBS stay in the body for a long time.

---

## Statistical Model: Hierarchical Bayesian Approach

### What Does "Hierarchical Bayesian" Mean?

Instead of estimating separate parameters for each person independently, we recognize that all 6 subjects are part of the same human population.

**Hierarchical structure**:
- There's a **population-level** average (how the typical person responds)
- Each **individual** has their own value, which varies around that average
- We estimate both the population average AND how much individuals vary

**Bayesian approach**:
- We start with prior beliefs about plausible parameter values
- We update these beliefs using the observed data
- We end up with a probability distribution (not just a single number) for each parameter

---

### Mathematical Formulation

#### Population-Level Parameters (Shared Across All Subjects)

**1. Population hyperparameters for RAFapi:**

- μ_pop: Mean of log(RAFapi) across the population
- σ_pop: Standard deviation of log(RAFapi) across the population

**2. Km_apical:** Michaelis constant for apical transporters (μg/L)
- Shared across all subjects

**3. σ:** Observation error (log scale)
- Shared across all subjects

**4. RAFbaso = 1.0:** Fixed (for identifiability)

---

#### Individual-Level Parameters (Subject-Specific)

For each subject i (i = 1, 2, ..., 6):

**RAFapi[i]:** Individual apical transporter activity factor

---

### Hierarchical Structure

The hierarchical model links individual and population parameters:

```
Population level:
  μ_pop ~ Prior_μ
  σ_pop ~ Prior_σ

Individual level (for each subject i):
  log(RAFapi[i]) ~ Normal(μ_pop, σ_pop)

Shared population parameters:
  Km_apical ~ Prior_Km
  σ ~ Prior_σobs
```

This means:
- Each individual's RAFapi is drawn from a population distribution
- The population distribution is characterized by μ_pop and σ_pop
- μ_pop and σ_pop are estimated from the data

---

### Prior Distributions

**Population hyperparameters:**

```
μ_pop ~ Normal(μ_prior, σ_prior)
  where: μ_prior = log(0.0007) - σ²_prior/2
         σ_prior = sqrt(log(1 + CV²))
         CV = 1.0 (100% coefficient of variation)

σ_pop ~ truncated(Normal(0, 1), 0, ∞)
  (constrained to be positive)
```

**Shared population parameters:**

```
Km_apical ~ LogNormal(μ_Km, σ_Km)
  where: μ_Km = log(52338.44) - σ²_Km/2
         σ_Km = sqrt(log(1 + 0.2²))
         (CV = 20%)

σ ~ truncated(Normal(0, 0.2), 0, ∞)
  (constrained to be positive)
```

**Individual parameters:**

```
For each subject i:
  log(RAFapi[i]) ~ Normal(μ_pop, σ_pop)
```

---

### Likelihood (How Data Relates to Model)

For each observation j in subject i:

```
log(C_obs[i,j]) ~ Normal(log(C_pred[i,j]), σ)
```

Where:
- **C_obs[i,j]**: Observed serum concentration for subject i at time j
- **C_pred[i,j]**: Model-predicted concentration using RAFapi[i], Km_apical, and fixed parameters
- **σ**: Observation error (shared across all subjects)

**Why log-normal?**
- Serum concentrations span multiple orders of magnitude
- Errors are proportional to concentration level (relative error)
- Log-normal distribution ensures concentrations remain positive

---

### Summary of Parameter Structure

| Parameter | Type | Dimensionality | Description |
|-----------|------|----------------|-------------|
| μ_pop | Population | Scalar | Mean of log(RAFapi) across population |
| σ_pop | Population | Scalar | SD of log(RAFapi) across population |
| RAFapi[i] | Individual | 6 values | Subject-specific transporter activity |
| Km_apical | Population | Scalar | Transporter binding affinity (shared) |
| σ | Population | Scalar | Observation error (shared) |
| RAFbaso | Fixed | Scalar | Basolateral activity = 1.0 (fixed) |

**Total estimated parameters:** 2 + 6 + 1 + 1 = 10 parameters

### MCMC Sampling

We use **NUTS (No-U-Turn Sampler)**, an advanced Markov Chain Monte Carlo algorithm that efficiently explores the parameter space.

- **4 chains** running in parallel (to check convergence)
- **1000 iterations** per chain
- **500 burn-in** iterations (discarded to ensure convergence)

**Output**:
- Posterior distributions for all parameters
- Uncertainty quantification (95% credible intervals)

---

## Partition Coefficients

### What Are Partition Coefficients?

Partition coefficients (PCs) describe how PFBS distributes between blood and tissues when at equilibrium.

**Example**: If liver:blood PC = 0.46, then at equilibrium, the concentration in liver is 46% of the blood concentration.

### How Were They Calculated?

We used the **Equilibrium Distribution Model (EDM)** developed by Allendorf et al., which calculates partition coefficients based on:

1. **Tissue composition**:
   - Water content
   - Protein content (albumin, phospholipids, neutral lipids)

2. **Chemical properties of PFBS**:
   - How strongly it binds to different tissue components
   - Its water solubility and protein binding affinity

### PFBS Partition Coefficients

| Tissue | Partition Coefficient (Tissue:Blood) |
|--------|-------------------------------------|
| Liver | 0.4587 |
| Kidney | 0.4539 |
| Rest of Body | 0.4104 |

**Interpretation**:
- PFBS concentrations in tissues are about 40-46% of blood concentrations
- PFBS stays primarily in blood rather than accumulating in tissues
- This is typical for highly water-soluble, protein-bound compounds

### Rest of Body Composition

"Rest of body" is a **mass-weighted average** of multiple tissues:
- Heart
- Gonads (reproductive organs)
- Skin
- Spleen
- Muscle

Each tissue's partition coefficient is weighted by its mass fraction in the body.

---

## Model Calibration Results

### Population-Level Estimates

From the hierarchical Bayesian calibration, we obtained:

**Population Hyperparameters:**
- **μ_pop** = -9.80 ± 0.43 (95% CI: [-10.61, -8.91])
  - Mean of log(RAFapi) across the population

- **σ_pop** = 0.47 ± 0.20 (95% CI: [0.19, 0.96])
  - Standard deviation of log(RAFapi) across the population
  - Corresponds to ~50% population coefficient of variation

**Mean Population RAFapi** = 7.14 × 10⁻⁵
- Indicates very low kidney reabsorption activity for PFBS in humans
- Explains why PFBS is eliminated much faster than longer-chain PFAS

**Transporter Binding Affinity:**
- **Km_apical** = 109,083 ± 45,508 μg/L (95% CI: [48,508, 225,681])
  - High Km indicates weak binding affinity

**Observation Error:**
- **σ** = 0.526 ± 0.059 (95% CI: [0.426, 0.662])

### Model Diagnostics

The model was validated using:

1. **Posterior Predictive Checks (PPC)**
   - Simulating new data from the fitted model
   - Comparing predictions to observed data
   - Good fit: observed data falls within 95% prediction interval

2. **MCMC Diagnostics**
   - **R̂ (R-hat)**: Should be < 1.1 (indicates chain convergence)
   - **Effective Sample Size (ESS)**: Should be > 400 (indicates sufficient sampling)
   - **Trace Plots**: Should show good mixing without trends

---

## Software & Reproducibility

### Implementation

- **Language**: Julia (for Bayesian calibration) and R (for deployment)
- **ODE Solver**: High-precision adaptive solvers (Rodas5P in Julia, lsodes in R)
- **Bayesian Inference**: Turing.jl with NUTS sampler
- **Visualization**: Plots.jl and StatsPlots.jl

### Files

- **MCMC Script**: `PFBS_population_MCMC.jl`
- **Posterior Predictive Checks**: `PFBS_population_ppc.jl`
- **Chain Diagnostics**: `view_chains.jl`
- **R Implementation**: `jaqpot_deployment/deploy_on_jaqpot.R`

### Reproducibility

All code and data are version-controlled and documented. The model can be re-run to verify results or adapted for new scenarios.

---

## References

### Model Framework
- **Worley et al.**: OWorley RR, Yang X, Fisher J. Physiologically based pharmacokinetic modeling of human exposure to perfluorooctanoic acid suggests historical non drinking-water exposures are important for predicting current serum concentrations. Toxicol Appl Pharmacol. 2017 Sep 1;330:9-21. doi: 10.1016/j.taap.2017.07.001. Epub 2017 Jul 3. PMID: 28684146; PMCID: PMC5664934.

### Data Sources
- **Haug et al. (2009)**: Background PFBS exposure estimates
- **Experimental data**: Olsen GW, Chang SC, Noker PE, Gorman GS, Ehresman DJ, Lieder PH, Butenhoff JL. A comparison of the pharmacokinetics of perfluorobutanesulfonate (PFBS) in rats, monkeys, and humans. Toxicology. 2009 Feb 4;256(1-2):65-74. doi: 10.1016/j.tox.2008.11.008. Epub 2008 Nov 19. PMID: 19059455.

### Partition Coefficients
- **Allendorf EDM**: Allendorf F, Goss KU, Ulrich N. Estimating the Equilibrium Distribution of Perfluoroalkyl Acids and 4 of Their Alternatives in Mammals. Environ Toxicol Chem. 2021 Mar;40(3):910-920. doi: 10.1002/etc.4954. Epub 2021 Feb 9. PMID: 33289938.
