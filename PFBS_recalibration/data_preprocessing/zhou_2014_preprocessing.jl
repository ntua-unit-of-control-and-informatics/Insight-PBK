using CSV
using DataFrames
using Statistics

"""
    preprocess_zhou_2014_data()

Preprocess Zhou et al. (2014) PFBS concentration data by age group.
Converts 25th and 75th percentiles (P25, P75) to standard deviation.

Focuses on Age_group 3 data only.
P25 and P75 represent 25th and 75th percentiles respectively.

Converts percentiles to SD using: SD = (P75 - P25) / 1.349
(where 1.349 is the interquartile range factor for normal distribution)

Treats data as steady-state serum measurements.

Output:
- zhou_2014_serum_summary.csv: Serum concentration with mean ± SD for Age_group 3
- zhou_2014_urine_summary.csv: Urine concentration with mean ± SD from original dataset
"""
function preprocess_zhou_2014_data()
    
    println("Loading Zhou et al. (2014) PFBS data...")
    
    # Load serum data by age group (for Age_group 3)
    serum_age_data = CSV.read("../../Data/Data_files/Zhou_2014_PFBS_concentration_by_group_of_years.csv", DataFrame)
    
    # Load original data for urine measurements
    original_data = CSV.read("../../Data/Data_files/Zhou_2014.csv", DataFrame)
    
    println("Serum age data dimensions: $(nrow(serum_age_data)) rows × $(ncol(serum_age_data)) columns")
    println("Original data dimensions: $(nrow(original_data)) rows × $(ncol(original_data)) columns")
    
    # Filter for Age_group 3 from serum data
    age_group_3 = filter(row -> row.Age_group == 3, serum_age_data)
    
    println("Age_group 3 data found:")
    for row in eachrow(age_group_3)
        println("  Age_group $(row.Age_group): Mean = $(row.Mean) ng/ml, 25th percentile = $(row.P25), 75th percentile = $(row.P75), N = $(row.N)")
    end
    
    # Study parameters
    time_years = 6
    time_days = time_years * 365
    
    println("\nStudy parameters:")
    println("  Exposure duration: $time_years years (steady-state)")
    println("  Time point: $(round(Int, time_days)) days")
    println("  Data format: Mean with 25th and 75th percentiles")
    println("  25th-75th percentile to SD conversion: SD = (P75 - P25) / 1.349")
    
    # Process serum data (Age_group 3 represents serum measurements)
    println("\nProcessing serum data...")
    serum_summary = DataFrame(
        time_days = Float64[],
        n_subjects = Int64[],
        mean_ng_ml = Float64[],
        sd_ng_ml = Float64[]
    )
    
    if nrow(age_group_3) > 0
        data_row = age_group_3[1, :]
        mean_val = data_row.Mean
        p25 = data_row.P25  # 25th percentile
        p75 = data_row.P75  # 75th percentile
        n_subjects = data_row.N
        
        # Convert 25th-75th percentile to SD: SD = (P75 - P25) / 1.349
        sd_val = (p75 - p25) / 1.349
        
        push!(serum_summary, (time_days, n_subjects, mean_val, sd_val))
        
        println("  Time: $time_years years ($(round(Int, time_days)) days)")
        println("  Mean: $(round(mean_val, digits=1)) ng/ml")
        println("  25th percentile: $(round(p25, digits=1)) ng/ml")
        println("  75th percentile: $(round(p75, digits=1)) ng/ml")
        println("  → SD: $(round(sd_val, digits=1)) ng/ml")
        println("  n = $n_subjects subjects")
    end
    
    # Process urine data from original Zhou_2014.csv
    println("\nProcessing urine data from original dataset...")
    pfbs_urine = filter(row -> row.PFAS == "PFBS" && row.Media == "Urine", original_data)
    
    urine_summary = DataFrame(
        time_days = Float64[],
        n_subjects = Int64[],
        mean_ng_ml = Float64[],
        sd_ng_ml = Float64[]
    )
    
    if nrow(pfbs_urine) > 0
        urine_row = pfbs_urine[1, :]
        mean_val = urine_row.Mean
        sd_val = urine_row.SD
        n_subjects = urine_row.N
        
        push!(urine_summary, (time_days, n_subjects, mean_val, sd_val))
        
        println("  Time: $time_years years ($(round(Int, time_days)) days)")
        println("  Mean: $(round(mean_val, digits=1)) ng/ml")
        println("  SD: $(round(sd_val, digits=1)) ng/ml")
        println("  n = $n_subjects subjects")
    end
    
    # Export results
    println("\nExporting preprocessed data...")
    
    CSV.write("zhou_2014_serum_summary.csv", serum_summary)
    println("Serum summary saved to: zhou_2014_serum_summary.csv")
    
    CSV.write("zhou_2014_urine_summary.csv", urine_summary)
    println("Urine summary saved to: zhou_2014_urine_summary.csv")
    
    # Display summary statistics
    println("\n" * "="^60)
    println("ZHOU ET AL. (2014) PFBS AGE_GROUP 3 SUMMARY")
    println("="^60)
    
    if nrow(serum_summary) > 0
        println("\nSERUM MEASUREMENTS:")
        println("-"^40)
        for row in eachrow(serum_summary)
            println("Day $(round(Int, row.time_days)): $(row.n_subjects) subjects, $(round(row.mean_ng_ml, digits=1)) ± $(round(row.sd_ng_ml, digits=1)) ng/ml")
        end
    end
    
    if nrow(urine_summary) > 0
        println("\nURINE MEASUREMENTS:")
        println("-"^40)
        for row in eachrow(urine_summary)
            println("Day $(round(Int, row.time_days)): $(row.n_subjects) subjects, $(round(row.mean_ng_ml, digits=1)) ± $(round(row.sd_ng_ml, digits=1)) ng/ml")
        end
    end
    
    return (serum = serum_summary, urine = urine_summary)
end

# Run preprocessing if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    results = preprocess_zhou_2014_data()
    println("\nPreprocessing completed successfully!")
end