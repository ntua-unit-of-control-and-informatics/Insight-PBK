using CSV
using DataFrames
using Statistics

# Include the original functions
include("calculate_partition_coefficients.jl")

"""
    calculate_partition_coefficients_plasma(tissue_composition_file, distribution_coeffs_file, tissue_mass_file, compounds, tissues)

Calculate tissue:plasma partition coefficients for PFAS compounds using the Allendorf equilibrium distribution model.

This is identical to calculate_partition_coefficients() but uses PLASMA as the reference compartment
instead of BLOOD.

# Arguments
- `tissue_composition_file`: Path to CSV file with tissue composition data
- `distribution_coeffs_file`: Path to CSV file with distribution coefficients (log scale)
- `tissue_mass_file`: Path to CSV file with tissue mass fractions for rest of body calculation
- `compounds`: Vector of compound names to calculate (e.g., ["PFBS", "PFHxA"])
- `tissues`: Vector of tissue names to calculate (e.g., ["liver", "kidney", "rest of body"])

# Returns
- Dictionary with partition coefficients for each tissue-compound pair
"""
function calculate_partition_coefficients_plasma(tissue_composition_file, distribution_coeffs_file, tissue_mass_file, compounds, tissues)

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

        # Get PLASMA composition for final calculation (KEY CHANGE)
        plasma_row = filter(row -> lowercase(row.Tissue) == "plasma", tissue_data)
        if nrow(plasma_row) == 0
            @error "Plasma composition not found in tissue data"
            continue
        end
        plasma_comp = plasma_row[1, :]

        # Calculate plasma organ-water distribution coefficient (KEY CHANGE)
        D_plasma_w = calculate_organ_water_distribution(
            plasma_comp, D_ml_w, D_alb_w, D_sp_w, D_sl_w, D_fabp_w, compound
        )

        if ismissing(D_plasma_w)
            @warn "Could not calculate plasma distribution coefficient for $compound"
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

            # Calculate partition coefficient: PC = D_tissue_w / D_plasma_w (KEY CHANGE)
            PC = D_tissue_w / D_plasma_w
            results[compound][tissue] = PC

            println("$compound - $tissue: PC (tissue:plasma) = $(round(PC, digits=3))")
        end
    end

    return results
end

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    # Define target compounds and tissues
    compounds = ["PFBS", "PFHxA"]
    tissues = ["liver", "kidney", "rest of body"]

    # Calculate tissue:plasma partition coefficients
    results = calculate_partition_coefficients_plasma(
        "volume_fractions.csv",
        "distribution_coefficients.csv",
        "tissue_fractions_allendorf_s1.8.1.csv",
        compounds,
        tissues
    )

    # Display results in a formatted table
    println("\n" * "="^60)
    println("TISSUE:PLASMA PARTITION COEFFICIENTS")
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
    println("\nSaving results to partition_coefficients_plasma_results.csv...")

    # Create DataFrame for export
    export_data = []
    for compound in compounds
        if haskey(results, compound)
            for tissue in tissues
                if haskey(results[compound], tissue)
                    push!(export_data, (
                        Compound = compound,
                        Tissue = tissue,
                        Partition_Coefficient_Plasma = results[compound][tissue]
                    ))
                end
            end
        end
    end

    if !isempty(export_data)
        df_results = DataFrame(export_data)
        CSV.write("../Worley_model/Worley_partition_coefficients_plasma_results.csv", df_results)
        println("Results saved successfully!")
    end
end
