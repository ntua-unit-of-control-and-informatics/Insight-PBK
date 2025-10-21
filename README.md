# Insight-PBK

Physiologically-based pharmacokinetic (PBK) modeling for PFAS compounds (PFBS and PFHxA) developed for the Insight EU Project. This work implements the Worley PBK model with active renal reabsorption to estimate transporter activity parameters from biomonitoring data.

## Overview

This repository focuses on modeling renal reabsorption mechanisms in PFAS compounds using Bayesian inference to calibrate transporter activity factors from synthetic steady-state observations.

### Key Components

- **Worley PBK Model**: Compartmental model with active renal reabsorption (basolateral and apical transporters)
- **Abraham Model**: One-compartment model for reverse dosimetry to estimate daily intake from serum concentrations
- **EDM Partition Coefficients**: Allendorf equilibrium distribution model for tissue:blood partition coefficients

### Workflow

1. **Intake Estimation**: Reverse dosimetry with Abraham model to estimate daily PFAS intake
2. **Partition Coefficients**: Calculate tissue:blood ratios using EDM from tissue composition data
3. **Bayesian Calibration**: MCMC estimation of apical transporter activity (RAFapi) with Turing.jl
4. **Model Validation**: Posterior predictive checks with 95% prediction intervals

### Implementation

- **Language**: Julia
- **Inference**: NUTS sampler (Turing.jl)
- **Model Structure**: Continuous oral exposure with log-normal observation error
- **Compounds**: PFBS and PFHxA

## Project Structure

```
Worley_model/
├── PBK_model/          # Core PBK model functions
├── PFBS/               # PFBS calibration and validation
└── PFHxA/              # PFHxA calibration and validation

EDM_Partition_Coefficients/  # Partition coefficient calculations
```
