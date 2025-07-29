using CSV
using DataFrames
using Statistics

"""
    preprocess_he_2023_data()

Preprocess He et al. (2023) PFBS occupational exposure data by converting 
95% confidence intervals to standard deviation for serum and urine measurements.

Converts 95% CI to SD using: SD = (CI_upper - CI_lower) × √N / 3.92

Treats data as steady-state measurements after 20 years of occupational exposure.

Output:
- he_2023_serum_summary.csv: Serum concentration with mean ± SD
- he_2023_urine_summary.csv: Urine concentration with mean ± SD
"""
function preprocess_he_2023_data()
    
    println("Loading He et al. (2023) PFBS occupational exposure data...")
    
    # Load raw data
    data = CSV.read("../../Data/Data_files/He_2023.csv", DataFrame)
    
    println("Raw data dimensions: $(nrow(data)) rows × $(ncol(data)) columns")
    
    # Filter for PFBS occupational exposure only
    pfbs_occupational = filter(row -> row.PFAS == "PFBS" && row.Exposure == "Occupational", data)
    
    println("PFBS occupational measurements found:")
    for row in eachrow(pfbs_occupational)
        println("  $(row.Media): Mean = $(row.Mean_ng_ml) ng/ml, 95% CI: $(row.L_CI) - $(row.U_CI), N = $(row.N)")
    end
    
    # Study parameters
    time_years = 20
    
    println("\nStudy parameters:")
    println("  Exposure duration: $time_years years (steady-state)")
    println("  95% CI to SD conversion: SD = (CI_upper - CI_lower) × √N / 3.92")
    
    # Process serum data
    println("\nProcessing serum data...")
    serum_data = filter(row -> row.Media == "Serum", pfbs_occupational)
    
    serum_summary = DataFrame(
        time_days = Float64[],
        n_subjects = Int64[],
        mean_ng_ml = Float64[],
        sd_ng_ml = Float64[]
    )
    
    if nrow(serum_data) > 0
        serum_row = serum_data[1, :]
        mean_val = serum_row.Mean_ng_ml
        ci_lower = serum_row.L_CI
        ci_upper = serum_row.U_CI
        n_subjects = serum_row.N
        
        # Convert 95% CI to SD: SD = (CI_upper - CI_lower) × √N / 3.92
        sd_val = (ci_upper - ci_lower) * sqrt(n_subjects) / 3.92
        
        # Convert years to days for consistency with other studies
        time_days = time_years * 365.25
        
        push!(serum_summary, (time_days, n_subjects, mean_val, sd_val))
        
        println("  Time: $time_years years ($(Int(time_days)) days)")
        println("  Mean: $(round(mean_val, digits=1)) ng/ml")
        println("  95% CI: $(round(ci_lower, digits=1)) - $(round(ci_upper, digits=1)) → SD: $(round(sd_val, digits=1)) ng/ml")
        println("  n = $n_subjects workers")
    end
    
    # Process urine data
    println("\nProcessing urine data...")
    urine_data = filter(row -> row.Media == "Urine", pfbs_occupational)
    
    urine_summary = DataFrame(
        time_days = Float64[],
        n_subjects = Int64[],
        mean_ng_ml = Float64[],
        sd_ng_ml = Float64[]
    )
    
    if nrow(urine_data) > 0
        urine_row = urine_data[1, :]
        mean_val = urine_row.Mean_ng_ml
        ci_lower = urine_row.L_CI
        ci_upper = urine_row.U_CI
        n_subjects = urine_row.N
        
        # Convert 95% CI to SD: SD = (CI_upper - CI_lower) × √N / 3.92
        sd_val = (ci_upper - ci_lower) * sqrt(n_subjects) / 3.92
        
        # Convert years to days for consistency with other studies
        time_days = time_years * 365.25
        
        push!(urine_summary, (time_days, n_subjects, mean_val, sd_val))
        
        println("  Time: $time_years years ($(Int(time_days)) days)")
        println("  Mean: $(round(mean_val, digits=1)) ng/ml")
        println("  95% CI: $(round(ci_lower, digits=1)) - $(round(ci_upper, digits=1)) → SD: $(round(sd_val, digits=1)) ng/ml")
        println("  n = $n_subjects workers")
    end
    
    # Export results
    println("\nExporting preprocessed data...")
    
    CSV.write("he_2023_serum_summary.csv", serum_summary)
    println("Serum summary saved to: he_2023_serum_summary.csv")
    
    CSV.write("he_2023_urine_summary.csv", urine_summary)
    println("Urine summary saved to: he_2023_urine_summary.csv")
    
    # Display summary statistics
    println("\n" * "="^60)
    println("HE ET AL. (2023) PFBS OCCUPATIONAL EXPOSURE SUMMARY")
    println("="^60)
    
    if nrow(serum_summary) > 0
        println("\nSERUM MEASUREMENTS:")
        println("-"^40)
        for row in eachrow(serum_summary)
            println("Day $(Int(row.time_days)): $(row.n_subjects) workers, $(round(row.mean_ng_ml, digits=1)) ± $(round(row.sd_ng_ml, digits=1)) ng/ml")
        end
    end
    
    if nrow(urine_summary) > 0
        println("\nURINE MEASUREMENTS:")
        println("-"^40)
        for row in eachrow(urine_summary)
            println("Day $(Int(row.time_days)): $(row.n_subjects) workers, $(round(row.mean_ng_ml, digits=1)) ± $(round(row.sd_ng_ml, digits=1)) ng/ml")
        end
    end
    
    return (serum = serum_summary, urine = urine_summary)
end

# Run preprocessing if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    results = preprocess_he_2023_data()
    println("\nPreprocessing completed successfully!")
end