using CSV
using DataFrames
using Statistics

"""
    calculate_rest_of_body_composition(tissue_data, tissue_mass_data, requested_tissues)

Calculate mass-weighted average composition for rest of body compartment.
Automatically includes all available tissues except those explicitly requested by the user.

# Arguments
- `tissue_data`: DataFrame with tissue composition data
- `tissue_mass_data`: DataFrame with tissue mass fractions
- `requested_tissues`: Vector of tissue names that user wants as separate compartments

# Returns
- Named tuple with rest of body composition (sl, ml, alb, sp, w, fabp)
"""
function calculate_rest_of_body_composition(tissue_data, tissue_mass_data, requested_tissues=[])

    # Get all available tissues from the mass fractions data
    all_available_tissues = unique(tissue_mass_data.Compartment)

    # Normalize requested tissue names to match data (capitalize first letter)
    requested_tissues_normalized = [map_tissue_name(t) for t in requested_tissues if lowercase(t) != "rest of body"]

    # Exclude tissues that should never be in rest of body:
    # - Blood (reference compartment)
    # - Plasma (part of blood compartment)
    # - Rest_of_Body (avoid circular definition)
    # - Total_Body (this is a sum, not a tissue)
    # - User-requested tissues (they get their own compartments)
    exclude_tissues = vcat(requested_tissues_normalized, ["Blood", "Plasma", "Rest_of_Body", "Total_Body"])

    # Determine which tissues go into "rest of body"
    rest_body_tissues = [t for t in all_available_tissues if !(t in exclude_tissues)]

    println("\n" * "="^60)
    println("Determining rest of body composition:")
    println("  Available tissues: ", join(all_available_tissues, ", "))
    println("  Requested as separate compartments: ", join(requested_tissues_normalized, ", "))
    println("  Automatically assigned to rest of body: ", join(rest_body_tissues, ", "))
    println("="^60)
    
    # Get mass fractions for rest of body tissues
    rest_mass_fractions = Dict{String, Float64}()
    total_rest_mass = 0.0
    
    for tissue in rest_body_tissues
        mass_row = filter(row -> row.Compartment == tissue, tissue_mass_data)
        if nrow(mass_row) > 0
            mass_fraction = mass_row[1, :Mass_Fraction_Percent]
            rest_mass_fractions[tissue] = mass_fraction
            total_rest_mass += mass_fraction
        else
            @warn "Mass fraction for $tissue not found"
            rest_mass_fractions[tissue] = 0.0
        end
    end
    
    println("Rest of body tissues and mass fractions:")
    for (tissue, mass) in rest_mass_fractions
        println("  $tissue: $(round(mass, digits=2))%")
    end
    println("  Total rest of body mass: $(round(total_rest_mass, digits=2))%")
    
    # Calculate mass-weighted average composition
    rest_composition = Dict("sl" => 0.0, "ml" => 0.0, "alb" => 0.0, "sp" => 0.0, "w" => 0.0, "fabp" => 0.0)
    
    for tissue in rest_body_tissues
        # Get composition for this tissue
        comp_row = filter(row -> row.Tissue == tissue, tissue_data)
        if nrow(comp_row) == 0
            @warn "Composition for $tissue not found in tissue data"
            continue
        end
        
        tissue_comp = comp_row[1, :]
        weight = rest_mass_fractions[tissue] / total_rest_mass  # Normalize to sum to 1
        
        # Add weighted contribution to each matrix
        for matrix in ["sl", "ml", "alb", "sp", "w", "fabp"]
            matrix_symbol = Symbol(matrix)
            rest_composition[matrix] += tissue_comp[matrix_symbol] * weight
        end
    end
    
    println("\nCalculated rest of body composition:")
    for (matrix, value) in rest_composition
        println("  $matrix: $(round(value, digits=2))%")
    end
    
    # Validate that composition sums to approximately 100%
    total_comp = sum(values(rest_composition))
    if abs(total_comp - 100.0) > 1.0
        @warn "Rest of body composition sums to $(round(total_comp, digits=2))%, not 100%"
    end
    
    return (
        sl = rest_composition["sl"],
        ml = rest_composition["ml"], 
        alb = rest_composition["alb"],
        sp = rest_composition["sp"],
        w = rest_composition["w"],
        fabp = rest_composition["fabp"],
        Tissue = "Rest_of_Body"
    )
end

"""
    calculate_partition_coefficients(tissue_composition_file, distribution_coeffs_file, tissue_mass_file, compounds, tissues)

Calculate tissue:blood partition coefficients for PFAS compounds using the Allendorf equilibrium distribution model.

The "rest of body" compartment is automatically calculated as a mass-weighted average of all tissues
NOT explicitly listed in the tissues vector (excluding Blood which is the reference compartment).

# Arguments
- `tissue_composition_file`: Path to CSV file with tissue composition data
- `distribution_coeffs_file`: Path to CSV file with distribution coefficients (log scale)
- `tissue_mass_file`: Path to CSV file with tissue mass fractions for rest of body calculation
- `compounds`: Vector of compound names to calculate (e.g., ["PFBS", "PFHxS"])
- `tissues`: Vector of tissue names to calculate (e.g., ["liver", "kidney", "rest of body"])
  - Include "rest of body" to get the partition coefficient for all remaining tissues
  - All tissues NOT in this list will automatically be lumped into "rest of body"

# Returns
- Dictionary with partition coefficients for each tissue-compound pair

# Examples
```julia
# Example 1: Only liver and kidney as separate compartments
compounds = ["PFBS", "PFHxS"]
tissues = ["liver", "kidney", "rest of body"]
results = calculate_partition_coefficients(
    "volume_fractions.csv",
    "distribution_coefficients.csv",
    "tissue_fractions_allendorf_s1.8.1.csv",
    compounds,
    tissues
)
# Result: liver PC, kidney PC, and rest of body PC (includes adipose, gut, lung, brain, etc.)

# Example 2: Multiple tissues as separate compartments
tissues = ["liver", "adipose", "kidney", "gut", "lung", "brain", "rest of body"]
results = calculate_partition_coefficients(...)
# Result: PCs for each listed tissue, and rest of body PC (includes heart, gonads, skin, etc.)
```
"""
function calculate_partition_coefficients(tissue_composition_file, distribution_coeffs_file, tissue_mass_file, compounds, tissues)
    
    # Load data files
    println("Loading tissue composition data...")
    tissue_data = CSV.read(tissue_composition_file, DataFrame)
    
    println("Loading distribution coefficients data...")
    dist_coeffs = CSV.read(distribution_coeffs_file, DataFrame)
    
    println("Loading tissue mass fractions data...")
    tissue_mass_data = CSV.read(tissue_mass_file, DataFrame)
    
    # Initialize results dictionary
    results = Dict{String, Dict{String, Float64}}()
    
    # Process each compound
    for compound in compounds
        results[compound] = Dict{String, Float64}()
        
        # Get distribution coefficients for this compound
        compound_row = filter(row -> row.Compound == compound, dist_coeffs)
        if nrow(compound_row) == 0
            @warn "Compound $compound not found in distribution coefficients data"
            continue
        end
        
        compound_data = compound_row[1, :]
        
        # Convert log-scale to linear scale, handling missing values
        D_ml_w = ismissing(compound_data.Dml_w) || compound_data.Dml_w == -999 ? missing : 10^compound_data.Dml_w
        D_alb_w = ismissing(compound_data.Dalb_w) || compound_data.Dalb_w == -999 ? missing : 10^compound_data.Dalb_w
        D_sp_w = ismissing(compound_data.Dsp_w) || compound_data.Dsp_w == -999 ? missing : 10^compound_data.Dsp_w
        D_sl_w = ismissing(compound_data.Dsl_w) || compound_data.Dsl_w == -999 ? missing : 10^compound_data.Dsl_w
        D_fabp_w = ismissing(compound_data.Dfapb_w) || compound_data.Dfapb_w == -999 || compound_data.Dfapb_w == "" ? missing : 10^compound_data.Dfapb_w
        
        # Check for essential coefficients
        if ismissing(D_alb_w)
            @warn "Missing albumin distribution coefficient for $compound, skipping"
            continue
        end
        
        # Get blood composition for final calculation
        blood_row = filter(row -> lowercase(row.Tissue) == "blood", tissue_data)
        if nrow(blood_row) == 0
            @error "Blood composition not found in tissue data"
            continue
        end
        blood_comp = blood_row[1, :]
        
        # Calculate blood organ-water distribution coefficient
        D_blood_w = calculate_organ_water_distribution(
            blood_comp, D_ml_w, D_alb_w, D_sp_w, D_sl_w, D_fabp_w, compound
        )
        
        if ismissing(D_blood_w)
            @warn "Could not calculate blood distribution coefficient for $compound"
            continue
        end
        
        # Process each tissue
        for tissue in tissues
            # Map tissue names to match data
            tissue_name = map_tissue_name(tissue)
            
            # Handle rest of body calculation
            if tissue_name == "Rest_of_Body"
                println("\nCalculating rest of body composition for $compound...")
                tissue_comp = calculate_rest_of_body_composition(tissue_data, tissue_mass_data, tissues)
            else
                # Get tissue composition from standard data
                tissue_row = filter(row -> lowercase(row.Tissue) == lowercase(tissue_name), tissue_data)
                if nrow(tissue_row) == 0
                    @warn "Tissue $tissue_name not found in composition data"
                    continue
                end
                tissue_comp = tissue_row[1, :]
            end
            
            # Calculate tissue organ-water distribution coefficient
            D_tissue_w = calculate_organ_water_distribution(
                tissue_comp, D_ml_w, D_alb_w, D_sp_w, D_sl_w, D_fabp_w, compound
            )
            
            if ismissing(D_tissue_w)
                @warn "Could not calculate tissue distribution coefficient for $compound in $tissue"
                continue
            end
            
            # Calculate partition coefficient: PC = D_tissue_w / D_blood_w
            PC = D_tissue_w / D_blood_w
            results[compound][tissue] = PC
            
            println("$compound - $tissue: PC = $(round(PC, digits=3))")
        end
    end
    
    return results
end

"""
    calculate_organ_water_distribution(tissue_comp, D_ml_w, D_alb_w, D_sp_w, D_sl_w, D_fabp_w, compound)

Calculate organ-water distribution coefficient using tissue composition and distribution coefficients.
"""
function calculate_organ_water_distribution(tissue_comp, D_ml_w, D_alb_w, D_sp_w, D_sl_w, D_fabp_w, compound)
    
    # Extract volume fractions (convert percentages to fractions)
    V_sl = tissue_comp.sl / 100.0
    V_ml = tissue_comp.ml / 100.0
    V_alb = tissue_comp.alb / 100.0
    V_sp = tissue_comp.sp / 100.0
    V_w = tissue_comp.w / 100.0
    V_fabp = tissue_comp.fabp / 100.0
    
    # Validate that fractions sum to approximately 1
    total_fraction = V_sl + V_ml + V_alb + V_sp + V_w + V_fabp
    if abs(total_fraction - 1.0) > 0.01
        @warn "Volume fractions for $(tissue_comp.Tissue) sum to $total_fraction, not 1.0"
    end
    
    # Calculate relative distribution coefficients (relative to albumin)
    # Handle missing values by setting relative coefficient to 0 (no partitioning to that matrix)
    D_ml_alb = ismissing(D_ml_w) || ismissing(D_alb_w) ? 0.0 : D_ml_w / D_alb_w
    D_sp_alb = ismissing(D_sp_w) || ismissing(D_alb_w) ? 0.0 : D_sp_w / D_alb_w
    D_sl_alb = ismissing(D_sl_w) || ismissing(D_alb_w) ? 0.0 : D_sl_w / D_alb_w
    D_fabp_alb = ismissing(D_fabp_w) || ismissing(D_alb_w) ? 0.0 : D_fabp_w / D_alb_w
    D_w_alb = ismissing(D_alb_w) ? missing : 1.0 / D_alb_w  # Water relative to albumin
    
    if ismissing(D_w_alb)
        return missing
    end
    
    # Calculate sorption fractions
    # fi,alb = 1 / (1 + D_w/alb × Vw/Valb + D_sp/alb × Vsp/Valb + D_ml/alb × Vml/Valb + D_sl/alb × Vsl/Valb + D_FABP/alb × VFABP/Valb)
    
    denominator = 1.0
    if V_alb > 0
        denominator += D_w_alb * (V_w / V_alb)
        denominator += D_sp_alb * (V_sp / V_alb)
        denominator += D_ml_alb * (V_ml / V_alb)
        denominator += D_sl_alb * (V_sl / V_alb)
        denominator += D_fabp_alb * (V_fabp / V_alb)
    else
        # Handle case where V_alb = 0 (no albumin in tissue)
        if V_w > 0 || V_sp > 0 || V_ml > 0 || V_sl > 0 || V_fabp > 0
            # Distribute among available matrices proportionally
            total_other = V_w + V_sp + V_ml + V_sl + V_fabp
            f_w = V_w / total_other
            f_sp = V_sp / total_other
            f_ml = V_ml / total_other
            f_sl = V_sl / total_other
            f_fabp = V_fabp / total_other
            f_alb = 0.0
        else
            @warn "No available matrices for partitioning in $(tissue_comp.Tissue)"
            return missing
        end
    end
    
    if V_alb > 0
        f_alb = 1.0 / denominator
        f_w = f_alb * D_w_alb * (V_w / V_alb)
        f_sp = f_alb * D_sp_alb * (V_sp / V_alb)
        f_ml = f_alb * D_ml_alb * (V_ml / V_alb)
        f_sl = f_alb * D_sl_alb * (V_sl / V_alb)
        f_fabp = f_alb * D_fabp_alb * (V_fabp / V_alb)
    end
    
    # Validate that fractions sum to 1
    total_f = f_alb + f_w + f_sp + f_ml + f_sl + f_fabp
    if abs(total_f - 1.0) > 0.01
        @warn "Sorption fractions for $compound in $(tissue_comp.Tissue) sum to $total_f, not 1.0"
    end
    
    # Calculate organ-water distribution coefficient
    # D_organ_w = fi,alb × Di,alb/w + fi,sp × Di,sp/w + fi,ml × Di,ml/w + fi,sl × Di,sl/w + fi,FABP × Di,FABP/w + fi,w × 1.0
    
    D_organ_w = f_w * 1.0  # Water contribution
    
    if !ismissing(D_alb_w)
        D_organ_w += f_alb * D_alb_w
    end
    
    if !ismissing(D_sp_w)
        D_organ_w += f_sp * D_sp_w
    end
    
    if !ismissing(D_ml_w)
        D_organ_w += f_ml * D_ml_w
    end
    
    if !ismissing(D_sl_w)
        D_organ_w += f_sl * D_sl_w
    end
    
    if !ismissing(D_fabp_w)
        D_organ_w += f_fabp * D_fabp_w
    end
    
    return D_organ_w
end

"""
    map_tissue_name(tissue)

Map common tissue names to those used in the data file.
"""
function map_tissue_name(tissue)
    tissue_mapping = Dict(
        "liver" => "Liver",
        "adipose" => "Adipose", 
        "kidney" => "Kidney",
        "gut" => "Gut",
        "lung" => "Lung",
        "lungs" => "Lung",
        "brain" => "Brain",
        "muscle" => "Muscle",
        "rest of body" => "Rest_of_Body"  # Will be calculated from mass-weighted average
    )
    
    return get(tissue_mapping, lowercase(tissue), tissue)
end

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    # Define target compounds and tissues
    compounds = ["PFBS", "PFHxA"]
    # tissues = ["liver", "adipose", "kidney", "gut", "lung", "brain", "rest of body"]
    tissues = ["liver", "kidney", "rest of body"]

    # Calculate partition coefficients
    results = calculate_partition_coefficients(
        "volume_fractions.csv",
        "distribution_coefficients.csv",
        "tissue_fractions_allendorf_s1.8.1.csv",
        compounds,
        tissues
    )
    
    # Display results in a formatted table
    println("\n" * "="^60)
    println("TISSUE:BLOOD PARTITION COEFFICIENTS")
    println("="^60)
    
    for compound in compounds
        if haskey(results, compound)
            println("\n$compound:")
            println("-"^40)
            for tissue in tissues
                if haskey(results[compound], tissue)
                    pc_value = round(results[compound][tissue], digits=3)
                    println("  $(rpad(tissue, 15)): $pc_value")
                end
            end
        end
    end
    
    # Save results to CSV
    println("\nSaving results to partition_coefficients_results.csv...")
    
    # Create DataFrame for export
    export_data = []
    for compound in compounds
        if haskey(results, compound)
            for tissue in tissues
                if haskey(results[compound], tissue)
                    push!(export_data, (
                        Compound = compound,
                        Tissue = tissue,
                        Partition_Coefficient = results[compound][tissue]
                    ))
                end
            end
        end
    end
    
    if !isempty(export_data)
        df_results = DataFrame(export_data)
        CSV.write("../Worley_model/Worley_partition_coefficients_results.csv", df_results)
        println("Results saved successfully!")
    end
end