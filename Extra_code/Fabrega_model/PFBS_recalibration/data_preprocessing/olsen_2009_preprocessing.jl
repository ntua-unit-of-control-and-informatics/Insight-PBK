using CSV
using DataFrames
using Statistics

"""
    preprocess_olsen_2009_data()

Preprocess Olsen et al. (2009) PFBS data by calculating mean and standard deviation
for serum and urine measurements at each time point across all subjects.

Handles:
- Missing values (empty cells)
- Below limit of quantification values ("LOQ")
- Individual subject measurements aggregated by time point

Output:
- olsen_2009_serum_summary.csv: Time points with mean ± SD serum concentrations
- olsen_2009_urine_summary.csv: Time points with mean ± SD urine concentrations
"""
function preprocess_olsen_2009_data()
    
    println("Loading Olsen et al. (2009) PFBS data...")
    
    # Load raw data
    data = CSV.read("../../Data/Data_files/olsen_2009.csv", DataFrame)
    
    println("Raw data dimensions: $(nrow(data)) rows × $(ncol(data)) columns")
    println("Subjects: $(length(unique(data.subject)))")
    println("Time points: $(sort(unique(data.day)))")
    
    # Function to clean numeric values
    function clean_numeric(value)
        if ismissing(value) || value == "" || value == "LOQ"
            return missing
        else
            return parse(Float64, string(value))
        end
    end
    
    # Clean serum and urine columns
    data.serum_clean = [clean_numeric(x) for x in data.serum]
    data.urine_clean = [clean_numeric(x) for x in data.urine]
    
    # Define time point grouping
    function group_time_point(day)
        if day == 0
            return 0
        elseif day in [2, 5, 8]
            return 5
        elseif day in [49, 57, 60, 62]
            return 60
        elseif day in [175, 180]
            return 180
        else
            return day  # Keep original if not in predefined groups
        end
    end
    
    # Apply time point grouping
    data.grouped_day = [group_time_point(day) for day in data.day]
    
    # Get unique grouped time points
    time_points = sort(unique(data.grouped_day))
    
    println("Time point grouping applied:")
    original_days = sort(unique(data.day))
    for orig_day in original_days
        grouped_day = group_time_point(orig_day)
        if orig_day != grouped_day
            println("  Day $orig_day → Day $grouped_day")
        else
            println("  Day $orig_day → Day $grouped_day (unchanged)")
        end
    end
    
    println("\nProcessing serum data...")
    
    # Process serum data
    serum_summary = DataFrame(
        time_days = Float64[],
        n_subjects = Int64[],
        mean_ng_ml = Float64[],
        sd_ng_ml = Float64[]
    )
    
    for time_point in time_points
        # Get serum measurements at this grouped time point
        time_data = filter(row -> row.grouped_day == time_point, data)
        serum_values = filter(!ismissing, time_data.serum_clean)
        
        if length(serum_values) > 0
            n_subjects = length(serum_values)
            mean_val = mean(serum_values)
            sd_val = length(serum_values) > 1 ? std(serum_values) : 0.0
            
            push!(serum_summary, (time_point, n_subjects, mean_val, sd_val))
            
            # Show which original days contributed to this grouped time point
            original_days_in_group = sort(unique(time_data.day))
            println("  Day $time_point (from days: $(join(original_days_in_group, ", "))): n=$n_subjects, mean=$(round(mean_val, digits=1)) ± $(round(sd_val, digits=1)) ng/ml")
        end
    end
    
    println("\nProcessing urine data...")
    
    # Process urine data
    urine_summary = DataFrame(
        time_days = Float64[],
        n_subjects = Int64[],
        mean_ng_ml = Float64[],
        sd_ng_ml = Float64[]
    )
    
    for time_point in time_points
        # Get urine measurements at this grouped time point
        time_data = filter(row -> row.grouped_day == time_point, data)
        urine_values = filter(!ismissing, time_data.urine_clean)
        
        if length(urine_values) > 0
            n_subjects = length(urine_values)
            mean_val = mean(urine_values)
            sd_val = length(urine_values) > 1 ? std(urine_values) : 0.0
            
            push!(urine_summary, (time_point, n_subjects, mean_val, sd_val))
            
            # Show which original days contributed to this grouped time point
            original_days_in_group = sort(unique(time_data.day))
            println("  Day $time_point (from days: $(join(original_days_in_group, ", "))): n=$n_subjects, mean=$(round(mean_val, digits=1)) ± $(round(sd_val, digits=1)) ng/ml")
        end
    end
    
    # Export results
    println("\nExporting preprocessed data...")
    
    CSV.write("olsen_2009_serum_summary.csv", serum_summary)
    println("Serum summary saved to: olsen_2009_serum_summary.csv")
    
    CSV.write("olsen_2009_urine_summary.csv", urine_summary)
    println("Urine summary saved to: olsen_2009_urine_summary.csv")
    
    # Display summary statistics
    println("\n" * "="^60)
    println("OLSEN ET AL. (2009) PFBS DATA SUMMARY")
    println("="^60)
    
    println("\nSERUM MEASUREMENTS:")
    println("-"^40)
    for row in eachrow(serum_summary)
        println("Day $(Int(row.time_days)): $(row.n_subjects) subjects, $(round(row.mean_ng_ml, digits=1)) ± $(round(row.sd_ng_ml, digits=1)) ng/ml")
    end
    
    println("\nURINE MEASUREMENTS:")
    println("-"^40)
    for row in eachrow(urine_summary)
        println("Day $(Int(row.time_days)): $(row.n_subjects) subjects, $(round(row.mean_ng_ml, digits=1)) ± $(round(row.sd_ng_ml, digits=1)) ng/ml")
    end
    
    return (serum = serum_summary, urine = urine_summary)
end

# Run preprocessing if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    results = preprocess_olsen_2009_data()
    println("\nPreprocessing completed successfully!")
end