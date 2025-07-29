using CSV
using DataFrames
using Statistics

"""
    preprocess_zhou_2014_data()

Preprocess Zhou et al. (2014) PFBS data by extracting mean and standard deviation
for serum and urine measurements.

Data already contains mean and SD, so no statistical conversion is needed.
Background levels are considered negligible and ignored.

Treats data as steady-state measurements after 10 years of exposure.

Output:
- zhou_2014_serum_summary.csv: Serum concentration with mean ± SD
- zhou_2014_urine_summary.csv: Urine concentration with mean ± SD
"""
function preprocess_zhou_2014_data()
    
    println("Loading Zhou et al. (2014) PFBS data...")
    
    # Load raw data
    data = CSV.read("../../Data/Data_files/Zhou_2014.csv", DataFrame)
    
    println("Raw data dimensions: $(nrow(data)) rows × $(ncol(data)) columns")
    
    # Filter for PFBS data only (excluding background)
    pfbs_data = filter(row -> row.PFAS == "PFBS", data)
    
    println("PFBS measurements found:")
    for row in eachrow(pfbs_data)
        println("  $(row.Media): Mean = $(row.Mean) ng/ml, SD = $(row.SD), N = $(row.N)")
    end
    
    # Study parameters
    time_years = 10
    time_days = time_years * 365.25
    
    println("\nStudy parameters:")
    println("  Exposure duration: $time_years years (steady-state)")
    println("  Time point: $(round(Int, time_days)) days")
    println("  Data format: Mean and SD already provided")
    
    # Process serum data
    println("\nProcessing serum data...")
    serum_data = filter(row -> row.Media == "Serum", pfbs_data)
    
    serum_summary = DataFrame(
        time_days = Float64[],
        n_subjects = Int64[],
        mean_ng_ml = Float64[],
        sd_ng_ml = Float64[]
    )
    
    if nrow(serum_data) > 0
        serum_row = serum_data[1, :]
        mean_val = serum_row.Mean
        sd_val = serum_row.SD
        n_subjects = serum_row.N
        
        push!(serum_summary, (time_days, n_subjects, mean_val, sd_val))
        
        println("  Time: $time_years years ($(round(Int, time_days)) days)")
        println("  Mean: $(round(mean_val, digits=1)) ng/ml")
        println("  SD: $(round(sd_val, digits=1)) ng/ml")
        println("  n = $n_subjects subjects")
    end
    
    # Process urine data
    println("\nProcessing urine data...")
    urine_data = filter(row -> row.Media == "Urine", pfbs_data)
    
    urine_summary = DataFrame(
        time_days = Float64[],
        n_subjects = Int64[],
        mean_ng_ml = Float64[],
        sd_ng_ml = Float64[]
    )
    
    if nrow(urine_data) > 0
        urine_row = urine_data[1, :]
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
    println("ZHOU ET AL. (2014) PFBS DATA SUMMARY")
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