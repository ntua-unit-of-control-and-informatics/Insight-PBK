# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository implements physiologically-based pharmacokinetic (PBK) modeling for PFAS (Per- and polyfluoroalkyl substances) compounds, specifically PFBS and PFHxA, developed for the Insight EU Project. The work focuses on modeling renal reabsorption mechanisms and estimating transporter activity parameters from biomonitoring data.

## Project Structure

### Worley PBK Model (`Worley_model/`)

The main modeling framework implementing the Worley et al. PBK model with active renal reabsorption:

- **PBK_model/julia_code.jl**: Core PBK model functions including compartmental ODEs and transporter kinetics
- **PFBS/**: PFBS-specific calibration, posterior predictive checks, and data
- **PFHxA/**: PFHxA-specific calibration, posterior predictive checks, and data

### Supporting Models

- **Abraham Model**: One-compartment kinetic model used in reverse dosimetry mode to estimate daily intake from observed serum concentrations
- **EDM Partition Coefficients** (`EDM_Partition_Coefficients/`): Allendorf equilibrium distribution model for calculating tissue:blood partition coefficients from tissue composition data

## Modeling Approach

### 1. Intake Estimation
Using the Abraham one-compartment model, we perform reverse dosimetry to estimate daily PFAS intake (ng/kg-bw/day) from biomonitoring serum concentration data.

### 2. Partition Coefficient Calculation
Tissue:blood partition coefficients (liver, kidney, rest of body) are calculated using the Allendorf EDM based on tissue composition and compound-specific properties.

### 3. Bayesian Parameter Calibration
MCMC sampling (Turing.jl with NUTS) is used to estimate the apical transporter activity factor (RAFapi):
- **Fixed**: RAFbaso = 1.0 (for identifiability)
- **Estimated**: RAFapi and sigma (observation error)
- **Prior**: LogNormal distribution for RAFapi with CV=100%
- **Likelihood**: Log-normal error model for plasma concentrations

### 4. Posterior Predictive Checks
Model validation through posterior predictive checks showing median predictions and 95% prediction intervals across continuous exposure scenarios.

## Key Dependencies

**Julia Packages:**
- `DifferentialEquations.jl` - ODE solving with adaptive time-stepping
- `Turing.jl` - Bayesian inference with NUTS sampler
- `DataFrames.jl`, `CSV.jl` - Data handling
- `Distributions.jl` - Statistical distributions
- `StatsPlots.jl`, `Plots.jl` - Visualization
- `Serialization.jl` - Saving/loading MCMC chains

## Technical Implementation Details

### Model Features
- Continuous oral exposure modeling with callback events
- Log-normal observation error model: log(CA_obs) ~ Normal(log(CA_pred), sigma)
- Tissue:blood partition coefficients calculated using EDM for liver, kidney, and rest of body
- Rest of body compartment represents weighted average of heart, gonads, skin, spleen, and muscle

### Calibration Strategy
- RAFbaso fixed at 1.0 for identifiability (basolateral transporter activity)
- RAFapi estimated from synthetic steady-state data (apical transporter activity)
- Synthetic observations: CA = 100 μg/L (PFBS) or 10 μg/L (PFHxA) at 10 years continuous exposure
- MCMC: 4 parallel chains, 1000 iterations, 500 burn-in

### Posterior Predictive Checks
Model validation using all posterior samples directly from MCMC chains to generate:
- Median prediction trajectories over 10 years
- 95% prediction intervals showing parameter and observation uncertainty
- Comparison with synthetic observations
- Individual plots for each subject saved as PNG files