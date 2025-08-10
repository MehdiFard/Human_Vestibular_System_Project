## Concatenate all phenotype tables produced by Pantry into one, along with a
## groups file that specifies all phenotypes per gene as one group. This allows
## for special downstream analyses such as QTL mapping that reports
## conditionally independent QTLs per gene even across modalities.
##
## Usage: bash scripts/combine_modalities.sh alt_polyA alt_TSS expression isoforms splicing stability
##
## Outputs: output/cross_modality.bed.gz, output/cross_modality.bed.gz.tbi, output/cross_modality.phenotype_groups.txt
##
## The phenotype names are modified by prepending '{modality}:' so they are
## unique and so their modalities can be recovered. After running, add 'cross_modality' as
## a modality in the Pheast config and generate QTLs for it.


#!/bin/bash

# Default value
output_path="output"

# Parse arguments
modalities=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --path_to_output)
            output_path="$2"
            shift 2
            ;;
        *)
            modalities+=("$1")
            shift
            ;;
    esac
done

# Make sure at least one modality was provided
if [ "${#modalities[@]}" -lt 1 ]; then
    echo "Usage: bash scripts/combine_modalities.sh [--path_to_output DIR] alt_polyA alt_TSS ..."
    exit 1
fi

# Confirm headers match
first="${modalities[0]}"
for modality in "${modalities[@]}"; do
    diff <(zcat "$output_path/$first.bed.gz" | head -n 1) <(zcat "$output_path/$modality.bed.gz" | head -n 1)
done

# Combine body rows and prepend modality name
for modality in "${modalities[@]}"; do
    echo "$modality"
    zcat "$output_path/$modality.bed.gz" | tail -n+2 | awk -v OFS='\t' -v m="$modality" '{$4 = m ":" $4; print}' >> "$output_path/cross_modality_tmp.bed"
done

# Create final BED file
zcat "$output_path/$first.bed.gz" | head -n 1 > "$output_path/cross_modality.bed"
sort -k1,1 -k2,2n "$output_path/cross_modality_tmp.bed" >> "$output_path/cross_modality.bed"
rm "$output_path/cross_modality_tmp.bed"

# Compress and index
bgzip -f "$output_path/cross_modality.bed"
tabix -p bed "$output_path/cross_modality.bed.gz"

# Make phenotype_groups file
zcat "$output_path/cross_modality.bed.gz" | tail -n+2 | cut -f4 | awk -F':' '{ g=$2; sub(/__.*$/, "", g); print $0 "\t" g }' > "$output_path/cross_modality.phenotype_groups.txt"

### END