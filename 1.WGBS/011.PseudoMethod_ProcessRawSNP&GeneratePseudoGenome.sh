#!/usr/bin/env bash
################################################################################
# Script: process_snp_and_pseudogenome.sh
# Description: Process raw SNP data and generate pseudo-genome for foxtail millet
# Author: Your Name
# Date: 2026-07-22
# Version: 1.0
#
# Usage: 
#   ./process_snp_and_pseudogenome.sh [OPTIONS]
#
# Options:
#   -i, --input-dir      Directory containing raw .gz SNP files [default: /public1/home/miyj/Project/00.SNP/chen]
#   -o, --output-dir     Output directory for processed files [default: /public1/home/miyj/projects/sitalica/snp]
#   -r, --ref-genome     Reference genome file path [default: /public1/home/chenjf/Share/Foxtailmillet/00.Reference/Sitalica/v2.2/assembly/Sitalica_312_v2_ptld.fa]
#   -s, --script-dir     Directory containing R and Python scripts [default: /public1/home/miyj/projects/sitalica/snp]
#   -t, --threads        Number of CPU threads to use [default: 4]
#   -h, --help           Show this help message
#
# Dependencies:
#   - R (with required packages)
#   - Python (with snp_filter.py script)
#   - awk, sed, less, gunzip
#
################################################################################

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

#===============================================================================
# Configuration variables
#===============================================================================

# Default paths
INPUT_DIR="${INPUT_DIR:-/public1/home/miyj/Project/00.SNP/chen}"
OUTPUT_DIR="${OUTPUT_DIR:-/public1/home/miyj/projects/sitalica/snp}"
REF_GENOME="${REF_GENOME:-/public1/home/chenjf/Share/Foxtailmillet/00.Reference/Sitalica/v2.2/assembly/Sitalica_312_v2_ptld.fa}"
SCRIPT_DIR="${SCRIPT_DIR:-/public1/home/miyj/projects/sitalica/snp}"
THREADS="${THREADS:-4}"

# Program paths
PYTHON_FILTER="${SCRIPT_DIR}/snp_filter.py"
R_SCRIPT="${SCRIPT_DIR}/pseudoRef_test.R"
HEADER_FILE="${OUTPUT_DIR}/header.txt"

#===============================================================================
# Functions
#===============================================================================

print_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Process raw SNP data and generate pseudo-genome for foxtail millet.

OPTIONS:
    -i, --input-dir DIR     Input directory with raw .gz SNP files
                            [default: ${INPUT_DIR}]
    -o, --output-dir DIR    Output directory for processed files
                            [default: ${OUTPUT_DIR}]
    -r, --ref-genome FILE   Reference genome FASTA file
                            [default: ${REF_GENOME}]
    -s, --script-dir DIR    Directory with R and Python scripts
                            [default: ${SCRIPT_DIR}]
    -t, --threads NUM       Number of CPU threads [default: ${THREADS}]
    -h, --help              Show this help message

DEPENDENCIES:
    Requires the following scripts:
        - ${PYTHON_FILTER}
        - ${R_SCRIPT}
    
    Requires R with packages and Python 3.

EXAMPLES:
    # Run with default settings
    $0
    
    # Specify custom directories
    $0 -i /path/to/input -o /path/to/output -r /path/to/reference.fa
    
    # Run with 8 threads
    $0 -t 8

EOF
    exit 0
}

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_dependencies() {
    local missing_deps=0
    
    # Check required scripts
    for script in "${PYTHON_FILTER}" "${R_SCRIPT}"; do
        if [[ ! -f "${script}" ]]; then
            log_message "ERROR: Required script not found: ${script}"
            missing_deps=1
        fi
    done
    
    # Check reference genome
    if [[ ! -f "${REF_GENOME}" ]]; then
        log_message "ERROR: Reference genome not found: ${REF_GENOME}"
        missing_deps=1
    fi
    
    # Check input directory
    if [[ ! -d "${INPUT_DIR}" ]]; then
        log_message "ERROR: Input directory not found: ${INPUT_DIR}"
        missing_deps=1
    fi
    
    # Create output directory if it doesn't exist
    mkdir -p "${OUTPUT_DIR}"
    
    if [[ ${missing_deps} -eq 1 ]]; then
        log_message "ERROR: Missing dependencies. Please check the paths."
        exit 1
    fi
}

#===============================================================================
# Parse command line arguments
#===============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -r|--ref-genome)
            REF_GENOME="$2"
            shift 2
            ;;
        -s|--script-dir)
            SCRIPT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            print_help
            ;;
        *)
            log_message "ERROR: Unknown option: $1"
            print_help
            ;;
    esac
done

#===============================================================================
# Main pipeline
#===============================================================================

main() {
    local start_time=$(date +%s)
    
    log_message "=========================================="
    log_message "Starting SNP processing pipeline"
    log_message "=========================================="
    log_message "Input directory:  ${INPUT_DIR}"
    log_message "Output directory: ${OUTPUT_DIR}"
    log_message "Reference genome: ${REF_GENOME}"
    log_message "Script directory: ${SCRIPT_DIR}"
    log_message "Threads:          ${THREADS}"
    log_message "=========================================="
    
    # Check dependencies
    check_dependencies
    
    #===========================================================================
    # Step 1: Uncompress .gz files and extract SNP data
    #===========================================================================
    
    log_message "STEP 1: Processing compressed SNP files"
    
    # Change to input directory
    cd "${INPUT_DIR}" || { log_message "ERROR: Cannot change to input directory"; exit 1; }
    
    # Process all .gz files
    local gz_count=0
    for gz_file in *.gz; do
        if [[ -f "${gz_file}" ]]; then
            gz_count=$((gz_count + 1))
            local base_name="${gz_file%.gz}"
            local output_file="${OUTPUT_DIR}/${base_name}"
            
            log_message "  Uncompressing: ${gz_file} -> ${output_file}"
            
            # Uncompress and save to output directory
            gunzip -c "${gz_file}" > "${output_file}" || {
                log_message "  ERROR: Failed to uncompress ${gz_file}"
                continue
            }
        fi
    done
    
    log_message "  Processed ${gz_count} .gz files"
    
    #===========================================================================
    # Step 2: Remove header lines (first 4 lines) from .snps files
    #===========================================================================
    
    log_message "STEP 2: Removing headers from SNP files"
    
    cd "${OUTPUT_DIR}" || { log_message "ERROR: Cannot change to output directory"; exit 1; }
    
    local snps_count=0
    for snps_file in *.snps; do
        if [[ -f "${snps_file}" ]]; then
            snps_count=$((snps_count + 1))
            log_message "  Removing header from: ${snps_file}"
            
            # Remove first 4 lines (in-place)
            sed -i '1,4d' "${snps_file}" || {
                log_message "  ERROR: Failed to process ${snps_file}"
                continue
            }
        fi
    done
    
    log_message "  Processed ${snps_count} .snps files"
    
    #===========================================================================
    # Step 3: Filter SNPs using Python script
    #===========================================================================
    
    log_message "STEP 3: Filtering SNPs with Python"
    
    local txt_count=0
    for snps_file in *.snps; do
        if [[ -f "${snps_file}" ]]; then
            txt_count=$((txt_count + 1))
            local base_name="${snps_file%.snps}"
            local txt_output="${base_name}.txt"
            
            log_message "  Filtering: ${snps_file} -> ${txt_output}"
            
            python "${PYTHON_FILTER}" "${snps_file}" "${txt_output}" || {
                log_message "  ERROR: Python filtering failed for ${snps_file}"
                continue
            }
        fi
    done
    
    log_message "  Processed ${txt_count} .snps files"
    
    #===========================================================================
    # Step 4: Format SNP files (extract specific columns)
    #===========================================================================
    
    log_message "STEP 4: Formatting SNP data"
    
    local format_count=0
    for txt_file in *.txt; do
        if [[ -f "${txt_file}" ]]; then
            format_count=$((format_count + 1))
            local base_name="${txt_file%.txt}"
            local snp_output="${base_name}.snp"
            
            log_message "  Formatting: ${txt_file} -> ${snp_output}"
            
            # Extract columns: $13, $1, $2, $3, $3 (sample, chr, pos, ref, alt)
            awk -v FS="\t" -v OFS="\t" '{print $13, $1, $2, $3, $3}' "${txt_file}" > "${snp_output}" || {
                log_message "  ERROR: Formatting failed for ${txt_file}"
                continue
            }
        fi
    done
    
    log_message "  Processed ${format_count} .txt files"
    
    #===========================================================================
    # Step 5: Rename .snp files back to .txt
    #===========================================================================
    
    log_message "STEP 5: Renaming .snp files to .txt"
    
    local rename_count=0
    for snp_file in *.snp; do
        if [[ -f "${snp_file}" ]]; then
            rename_count=$((rename_count + 1))
            local base_name="${snp_file%.snp}"
            
            log_message "  Renaming: ${snp_file} -> ${base_name}.txt"
            mv "${snp_file}" "${base_name}.txt" || {
                log_message "  ERROR: Rename failed for ${snp_file}"
                continue
            }
        fi
    done
    
    log_message "  Renamed ${rename_count} .snp files"
    
    # Optional: Clean up intermediate files
    # rm -f *.snps  # Uncomment to remove intermediate .snps files
    
    #===========================================================================
    # Step 6: Generate pseudo-genome using R
    #===========================================================================
    
    log_message "STEP 6: Generating pseudo-genome"
    
    # Get the first .txt file for processing (array job emulation)
    local first_txt=$(ls *.txt 2>/dev/null | head -n 1)
    
    if [[ -z "${first_txt}" ]]; then
        log_message "ERROR: No .txt files found for pseudo-genome generation"
        exit 1
    fi
    
    local prefix="${first_txt%.txt}"
    local result_file="${OUTPUT_DIR}/${prefix}.result"
    local header_file="${OUTPUT_DIR}/header.txt"
    
    log_message "  Processing: ${first_txt}"
    log_message "  Prefix: ${prefix}"
    log_message "  Result file: ${result_file}"
    
    # Create header for result file
    echo -e "chr\tpos\tref\talt\t${prefix}" > "${header_file}"
    
    # Activate R environment and run script
    log_message "  Running R script..."
    
    # Check if conda is available
    if command -v conda &> /dev/null; then
        source activate
        conda deactivate 2>/dev/null || true
        conda activate R 2>/dev/null || {
            log_message "WARNING: Cannot activate R conda environment, trying to use system R"
        }
    fi
    
    # Run R script
    Rscript "${R_SCRIPT}" "${result_file}" "${first_txt}" "${REF_GENOME}" || {
        log_message "ERROR: R script failed"
        exit 1
    }
    
    log_message "  Pseudo-genome generation completed"
    
    #===========================================================================
    # Cleanup and summary
    #===========================================================================
    
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    local h=$((runtime / 3600))
    local m=$(((runtime % 3600) / 60))
    local s=$((runtime % 60))
    
    log_message "=========================================="
    log_message "PIPELINE COMPLETED SUCCESSFULLY"
    log_message "=========================================="
    log_message "Total runtime: ${h}h ${m}m ${s}s"
    log_message "Output directory: ${OUTPUT_DIR}"
    log_message "=========================================="
    
    # Generate summary file
    cat > "${OUTPUT_DIR}/pipeline_summary.txt" << EOF
Pipeline Summary
================
Date: $(date)
Input directory: ${INPUT_DIR}
Output directory: ${OUTPUT_DIR}
Reference genome: ${REF_GENOME}
Threads: ${THREADS}
Total runtime: ${h}h ${m}m ${s}s

Processed files:
- GZ files: ${gz_count}
- SNP files: ${snps_count}
- TXT files (filtered): ${txt_count}
- Formatted files: ${format_count}
- Renamed files: ${rename_count}

Dependencies:
- Python script: ${PYTHON_FILTER}
- R script: ${R_SCRIPT}
EOF
    
    log_message "Summary written to: ${OUTPUT_DIR}/pipeline_summary.txt"
}

#===============================================================================
# Execute main function
#===============================================================================

# Trap errors for better error reporting
trap 'log_message "ERROR: Pipeline failed at line $LINENO"' ERR

# Run the main pipeline
main "$@"

exit 0