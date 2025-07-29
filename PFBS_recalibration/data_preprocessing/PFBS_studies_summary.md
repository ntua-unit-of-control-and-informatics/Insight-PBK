# PFBS Studies Summary for Parameter Recalibration

This document summarizes the three PFBS studies used for hierarchical MCMC parameter recalibration of the PBK model.

## Study Overview

| Study | Exposure Duration | Study Type | Estimated Exposure* | Sample Size |
|-------|------------------|------------|-------------------|-------------|
| Zhou et al. (2014) | 10 years | Background/General Population | 36.2 ng/kg-bw/day | 28 (serum), 36 (urine) |
| Olsen et al. (2009) | Single dose (elimination phase) | Pharmacokinetic | 745.1 ng/kg-bw/day | 6 subjects |
| He et al. (2023) | 20 years | Occupational | 121.4 ng/kg-bw/day | 60 workers |

*Estimated daily intake based on Abraham model reverse dosimetry (from `Abraham_model/Abraham_model_results.csv`)

## Measurement Time Points

### Zhou et al. (2014)
- **Time Points**: 1 (steady-state)
- **Day**: 3653 (10 years)
- **Data Type**: Single cross-sectional measurement
- **Available Data**: Serum and urine concentrations

### Olsen et al. (2009)
- **Time Points**: 4 (grouped)
- **Days**: 0, 5, 60, 180
- **Data Type**: Time-course pharmacokinetic study (**elimination phase after single dose**)
- **Available Data**: Serum and urine concentrations
- **Study Design**: Measurements taken after cessation of exposure to track elimination kinetics
- **Grouping Applied**:
  - Day 5: Combined measurements from days 2, 5, 8
  - Day 60: Combined measurements from days 49, 57, 60, 62
  - Day 180: Combined measurements from days 175, 180

### He et al. (2023)
- **Time Points**: 1 (steady-state)
- **Day**: 7305 (20 years)
- **Data Type**: Occupational exposure cross-sectional study
- **Available Data**: Serum and urine concentrations

## Data Processing Summary

| Study | Serum Data | Urine Data | Statistical Processing |
|-------|------------|------------|----------------------|
| Zhou et al. (2014) | ✓ | ✓ | Mean ± SD (direct) |
| Olsen et al. (2009) | ✓ | ✓ | Individual → Mean ± SD (grouped by time) |
| He et al. (2023) | ✓ | ✓ | 95% CI → SD conversion |

## Exposure Characteristics

- **Lowest Exposure**: Zhou et al. (36.2 ng/kg-bw/day) - General population
- **Moderate Exposure**: He et al. (121.4 ng/kg-bw/day) - Occupational
- **Highest Exposure**: Olsen et al. (745.1 ng/kg-bw/day) - Single pharmacokinetic dose

## Model Calibration Strategy

The hierarchical MCMC approach will use:
1. **Exposure estimates** from Abraham model reverse dosimetry
2. **Serum and urine concentration data** from all three studies
3. **Study-specific parameters** to account for population differences
4. **Population-level parameters** to capture overall PFBS kinetics

This multi-study approach provides:
- **Dose range coverage**: 3 orders of magnitude difference in exposure
- **Time scale diversity**: Acute (days) to chronic (years) exposure
- **Population diversity**: General, occupational, and clinical populations
- **Rich temporal data**: From Olsen elimination phase study for model validation
- **Kinetic insights**: Olsen data captures elimination kinetics after exposure cessation, critical for parameter calibration