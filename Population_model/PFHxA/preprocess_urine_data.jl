using DataFrames
using CSV

println("="^80)
println("PREPROCESS PFHxA URINE DATA")
println("="^80)

# Load digitized urine data
data = CSV.read("Population_model/PFHxA/digitized_data/Abraham_2024_urine_data.csv", DataFrame)
# Use excretion rates (ug/hour) to estimate the absolute amount excreted in each time interval
# transforme excretion rates from ug/hour to ug/day
data.excretion_ug_h = data.excretion_ug_h .* 24.0
#rename columns for clarity 
rename!(data, Symbol("Time_days") => :Time_days, Symbol("excretion_ug_h") => :Excretion_ug_day)

# Calculate cumulative amount excreted over time
data.Cumulative_Excretion_ug = zeros(nrow(data))
for row in 1:nrow(data)
    if row == 1
        data.Cumulative_Excretion_ug[row] = data.Excretion_ug_day[row] * (data.Time_days[row])
    else
        delta_time = data.Time_days[row] - data.Time_days[row - 1]
        data.Cumulative_Excretion_ug[row] = data.Cumulative_Excretion_ug[row - 1] + data.Excretion_ug_day[row] * delta_time
    end
end

# Save preprocessed data
CSV.write("Population_model/PFHxA/digitized_data/Abraham_2024_urine_data_processed.csv", data)