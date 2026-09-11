#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=200G
#SBATCH -t 100:00:00
#SBATCH -p amd_m9_768

start=`date +%s`

# Environment setup
echo "Loading ENV ..."
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate Ixchel

# Configuration
Ixchel=/publicfs10/fs10-m9/home/m9s003656/program/Ixchel/SourceCode/Ixchel.py
prefix=Sit60.d3
graph_input_dir=GraphMethyl_input
graph_output_dir=GraphMethyl_output

echo "Step 1: Prepare graph files ..."
python3 $Ixchel prepareGraphFiles ${prefix}.gfa --reference_name Yugu01

echo "Step 2: Precompute conversion files (array job) ..."
mkdir -p split_annotations

# Process annotation files in parallel
for ann_file in split_annotations/Annotations.Segments.${prefix}__*; do
    if [[ -f "$ann_file" ]]; then
        echo "Processing: $ann_file"
        python3 $Ixchel precompute_conversion \
            $ann_file \
            RefOnly.Segments.${prefix}.pkl \
            QueryOnly.Segments.${prefix}.pkl \
            FilteredLinks.Links.${prefix}.pkl \
            UpstreamArray.RefOnly.Segments.${prefix}.pkl \
            DownstreamArray.RefOnly.Segments.${prefix}.pkl \
            DoubleAnchored.FilteredLinks.Links.${prefix}.pkl
    fi
done

echo "Step 3: Merge converted files ..."
if [ ! -e Annotations.Segments.${prefix}__merged.converted ]; then
    cat split_annotations/Annotations.Segments.${prefix}__?????.converted > Annotations.Segments.${prefix}__merged.converted
fi

echo "Step 4: Build surjection database ..."
python3 $Ixchel build_db \
    Annotations.Segments.${prefix}__merged.converted \
    Annotations.Segments.${prefix}__merged.converted.db

echo "Step 5: Convert GraphMethyl to MethylC ..."
mkdir -p ${graph_output_dir}

for graph_file in ${graph_input_dir}/*.graph.methyl; do
    if [[ -f "$graph_file" ]]; then
        sample=$(basename $graph_file .graph.methyl)
        echo "Processing: $sample"
        python3 $Ixchel convertGraphMethylToMethylC \
            $graph_file \
            Annotations.Segments.${prefix}__merged.converted.db \
            ${graph_output_dir}/${sample}.graph.methyl.methylC
    fi
done

# Runtime statistics
end=`date +%s`
runtime=$((end-start))
h=$(($runtime/3600))
m=$((($runtime%3600)/60))
s=$(($runtime%60))

echo "Start: $start"
echo "End: $end"
echo "Run time: ${h}h ${m}m ${s}s"
echo "Done!"