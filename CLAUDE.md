# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

This repository contains physiologically-based pharmacokinetic (PBK) models for PFAS (Per- and polyfluoroalkyl substances) developed for the Insight EU Project. The main goal is to model PFAS exposure and biokinetics in humans for regulatory risk assessment.

## Key Architecture Components

### Core PBK Models
- **Primary Julia Model**: `jl_PBK_model/julia_model/julia_PBK_model.jl` - 12-compartment PBK model for PFBS, PFHxA, PFOA, PFOS
- **Shared Functions**: `jl_PBK_model/julia_model/pbk_functions.jl` - Common PBK modeling utilities
- **Simple Kinetic Model**: `Abraham_model/Abraham_kinetic_model.jl` - One-compartment model for reverse dosimetry
- **Partition Coefficient Calculator**: `EDM_Partition_Coefficients/calculate_partition_coefficients.jl` - Allendorf model for tissue:blood partition coefficients

### Model Structure
The main PBK model uses a 12-compartment system: Plasma, Gut, Liver, Fat, Lungs, Brain, Kidney, Filtrate, Storage, Urine, Rest of body, Ingestion. It supports both IV and oral exposure routes with Michaelis-Menten elimination kinetics.

### Bayesian Inference Workflow
1. **Exposure Estimation**: Use Abraham model with Monte Carlo reverse dosimetry to estimate exposure distributions from biomonitoring data
2. **Partition Coefficient Calculation**: Use Allendorf equilibrium distribution model to calculate tissue:blood partition coefficients from tissue composition data
3. **Parameter Recalibration**: Apply hierarchical MCMC (Turing.jl) to recalibrate PBK parameters using estimated exposures and paired serum-urine data

## Common Commands

### Julia Model Execution
```bash
julia jl_PBK_model/julia_model/julia_PBK_model.jl
```

### Monte Carlo Reverse Dosimetry
```bash
julia PFBS_recalibration/Monte_Carlo_Reverse_Dosimetry.jl
```

### Bayesian Parameter Estimation
```bash
julia PFBS_recalibration/PFBS_recalibration_hierarchical.jl
```

### Benchmarking Analysis
```bash
julia Benchmarking/PBK_julia_benchmarking.jl
```

### Partition Coefficient Calculation
```bash
julia EDM_Partition_Coefficients/calculate_partition_coefficients.jl
```

## Key Dependencies

**Julia Packages:**
- `DifferentialEquations.jl` - ODE solving
- `Turing.jl` - Bayesian inference
- `UnPack.jl`, `DataFrames.jl`, `CSV.jl` - Data handling
- `Distributions.jl` - Statistical distributions
- `StatsPlots.jl` - Visualization

## Data Sources

Experimental data in `Data/Data_files/` includes biomonitoring studies:
- He et al. (2023) - Occupational PFAS exposure
- Zhou et al. (2014) - Background population exposure  
- Olsen et al. (2009) - PFBS pharmacokinetics

## Development Notes

- Julia implementation provides 60-151x performance advantage over R equivalent
- Models are designed for both pharmacokinetic studies and biomonitoring data interpretation
- Hierarchical modeling accounts for between-study heterogeneity in parameter estimation
- LogNormal distributions used for naturally positive parameters (Vd, k_el, exposure doses)
- Partition coefficients calculated using Allendorf equilibrium distribution model with mass-weighted rest of body composition
- Rest of body compartment represents weighted average of Heart, Gonads, Skin, Spleen, and Muscle tissues