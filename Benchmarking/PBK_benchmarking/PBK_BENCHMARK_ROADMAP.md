# PBK Model Benchmarking Roadmap

## **Objective**
Compare execution times for solving complex physiologically-based pharmacokinetic (PBK) models between Julia and R using the existing 12-compartment PFAS models.

## **Model Specifications**

### **PBK Model Details**
- **Compartments**: 12 (Plasma, Gut, Liver, Fat, Lungs, Brain, Kidney, Filtrate, Storage, Urine, Rest of body, Ingestion)
- **Compounds**: PFBS and PFHxA (existing in both Julia and R implementations)
- **Physiological parameters**: Body weight = 70 kg
- **Same tolerances**: rtol=1e-8, atol=1e-10
- **Model complexity**: Multi-compartment with saturable kinetics (Michaelis-Menten)

### **Test Scenarios**
1. **Scenario 1**: 10 days, timestep = 1 day, continuous oral exposure (100 ng/kg/day) - **PFHxA**
2. **Scenario 2**: 10 days, timestep = 0.1 day, continuous oral exposure (100 ng/kg/day) - **PFHxA**
3. **Scenario 3**: 5 years, timestep = 1 day, continuous oral exposure (100 ng/kg/day) - **PFHxA**
4. **Scenario 4**: 5 years, timestep = 1 day, weekly oral doses (700 ng/kg BW per dose) - **PFBS**

## **Solver Matrix**

### **Julia Solvers** (6 total)
- `Tsit5()` - Default fast explicit RK (primary comparison)
- `CVODE_Adams()` - lsoda equivalent (Adams-Moulton)
- `CVODE_BDF()` - lsodes equivalent (external wrapper)
- `Rodas4()` - Native Julia stiff solver
- `Vern7()` - High-order explicit (potential fastest)
- `TRBDF2()` - Julia-native implicit solver

### **R Solvers** (2 total)
- `lsoda` - Automatic stiff/non-stiff switching (baseline)
- `lsodes` - Stiff solver with sparse Jacobian

## **Implementation Plan**

### **Files to Create in `PBK_benchmarking/`**

1. **`PBK_BENCHMARK_ROADMAP.md`** ✅ - This roadmap document
2. **`pbk_julia_benchmark.jl`** ✅ - Julia implementation
   - Adapted existing `julia_PBK_model.jl` with benchmarking framework
   - Test 6 solvers across 4 scenarios (PFHxA for 1-3, PFBS for 4)
   - 3 samples per test (reduced from 10)
   - Automated result collection and CSV export
   - Weekly dosing events for Scenario 4 (261 doses over 5 years)
3. **`pbk_r_benchmark.R`** ✅ - R implementation  
   - Adapted existing `Fabrega_PFOS_PFOA_Human_PBK.R` with benchmarking framework
   - Added PFBS kinetic parameters to match Julia implementation
   - Test 2 solvers across 4 scenarios (PFHxA for 1-3, PFBS for 4)
   - 3 samples per test
   - Weekly dosing events for Scenario 4 (261 doses over 5 years)
   - Automated result collection and CSV export
4. **Result files**: 
   - `julia_pbk_benchmark_results.csv` - Julia timing results
   - `julia_pbk_accuracy_results.csv` - Julia accuracy results
   - `r_pbk_benchmark_results.csv` - R timing results
   - `r_pbk_accuracy_results.csv` - R accuracy results

### **Benchmark Protocol**
- **Repetitions**: 3 samples per solver-scenario-compound combination
- **Output metrics**: Plasma concentration time series
- **Performance metrics**: Mean/median/std timing, memory usage, allocations
- **Accuracy verification**: Compare final plasma concentrations between solvers
- **Error handling**: Robust try-catch blocks for solver failures

### **Events Implementation**
- **Scenarios 1-3**: Continuous exposure via ingestion rate parameter
  - Constant ingestion rate throughout simulation
  - Biomonitoring exposure type
- **Scenario 4**: Weekly oral doses (700 ng/kg BW) using callback/events systems
  - **Julia**: `PresetTimeCallback` for weekly dosing at days 0, 7, 14, 21, ... for 5 years (261 doses total)
  - **R**: `events` parameter in `ode()` for weekly dosing schedule
  - Pharmacokinetic exposure type (bolus dosing)

### **Expected Testing Matrix**
- **Total Julia tests**: 6 solvers × 4 scenarios × 3 samples = **72 tests**
- **Total R tests**: 2 solvers × 4 scenarios × 3 samples = **24 tests**
- **Grand total**: **96 benchmark tests**
- **Estimated runtime**: 30-60 minutes for complete benchmark suite

## **Performance Hypotheses**

### **Expected Julia Advantages**
- **Compilation benefits**: Better performance on longer simulations (Scenarios 3-4)
- **Memory efficiency**: Lower allocation overhead
- **Solver ecosystem**: Native Julia solvers may outperform wrapped libraries

### **Complexity Impact**
- **PBK vs Simple**: Expect larger performance gaps due to model complexity
- **Event handling**: Scenario 4 will test discrete event performance
- **Stiff dynamics**: Some compartments may exhibit stiffness requiring appropriate solvers

## **Deliverables**

### **Primary Outputs**
1. **Comprehensive performance comparison tables** (Julia vs R)
2. **Solver-specific recommendations** for PBK modeling
3. **Scaling analysis** (simple kinetic vs complex PBK performance differences)
4. **Memory efficiency comparison** between languages

### **Analysis Components**
5. **Optimal solver selection guidelines** for PFAS PBK modeling
6. **Performance impact analysis** of dosing regimen complexity (continuous vs pulsed)
7. **Compound-specific performance** (PFBS vs PFHxA differences)
8. **Scenario complexity scaling** (computational cost vs problem size)

## **Success Criteria**
- All solvers successfully complete at least 80% of test scenarios
- Clear performance ranking established for each scenario type
- Reproducible results with statistical confidence (3 samples minimum)
- Practical recommendations for solver selection based on use case

## **Progress Status**

### ✅ Completed
1. **Roadmap planning** - Comprehensive benchmarking strategy defined
2. **Julia benchmark implementation** - `pbk_julia_benchmark.jl` created
   - 6 solvers: Tsit5, CVODE_Adams, CVODE_BDF, Rodas4, Vern7, TRBDF2
   - 4 scenarios with appropriate substance selection
   - Event handling for weekly dosing (Scenario 4)
   - Comprehensive error handling and result collection
3. **R benchmark implementation** - `pbk_r_benchmark.R` created
   - 2 solvers: lsoda, lsodes
   - Added PFBS kinetic parameters to match Julia model
   - Event handling for complex dosing scenarios
   - Same tolerance settings and scenario structure as Julia

### 🔄 Ready for Execution
4. **Benchmark suite ready** - Both implementations complete
   - Total: 96 tests (72 Julia + 24 R)
   - Identical scenarios and substance assignments
   - Same numerical tolerances for fair comparison

### ⏳ Pending
5. **Execution and data collection**
   - Run Julia benchmarks (72 tests)
   - Run R benchmarks (24 tests)
6. **Results analysis and reporting**
   - Performance comparison tables
   - Solver recommendations
   - Scaling analysis vs simple kinetic model