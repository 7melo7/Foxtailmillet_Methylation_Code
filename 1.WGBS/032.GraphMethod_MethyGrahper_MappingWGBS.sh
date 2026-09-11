#!/bin/bash
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=100G
##SBATCH --mem-per-cpu=100G
#SBATCH -t 100:00:00
#SBATCH -p amd_m9_768
#SBATCH --chdir=./

start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then CPU=2; fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then N=1; fi

#----------------
fq=`ls /publicfs10/fs10-m9/home/m9s003656/Project_SetariaMe/0.RawLibraryReads/*_R1.fq.gz | head -n $N | tail -n 1`
fq_dir=$(dirname $fq)
prefix=$(basename $fq _R1.fq.gz)

input_gfa=Sit60.d3.gfa
genome_prefix=Sit60.d3

lambda=lambda.fa

work_dir=out_$prefix
mkdir $work_dir

#----------------
echo "Loading ENV ..."
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate MethylGrapher
#which vg: ~/prgram/vg
 
#----------------
if [ ! -e $genome_prefix.wl.G2A.dist ]; then
  echo "<<<<< Prepare genome ..."
  # Adds lambda phage genome to genome graph, converts a GFA file into fully G->A and C->T converted GFA file, and indexes it for vg giraffe alignment.
  methylGrapher PrepareGenome -t 1 -gfa $input_gfa -lp $lambda -prefix $genome_prefix
fi

echo "<<<<< (just in case) link genome ..."
# Change the path !
ln -s /publicfs10/fs10-m9/home/m9s003656/Project_SetariaMe/3A.methylGrapher/Sit60.d3.* $work_dir

#----------------
echo "<<<<< (additional) fastp ..."
if [ ! -e $work_dir/${prefix}_R2.clean.fq.gz ]; then
  /publicfs10/fs10-m9/home/m9s003656/program/fastp -i $fq_dir/${prefix}_R1.fq.gz -I $fq_dir/${prefix}_R2.fq.gz \
    -o $work_dir/${prefix}_R1.clean.fq.gz -O $work_dir/${prefix}_R2.clean.fq.gz \
    -h $work_dir/$prefix.fastp.html -j $work_dir/$prefix.fastp.json \
    -z 4 -q 20 -u 30
fi

#----------------
echo "<<<<< Alignment ..."
# VG Giraffe alignment, please provide work directory and index prefix.
methylGrapher Align -t $CPU -fq1 $work_dir/${prefix}_R1.clean.fq.gz -fq2 $work_dir/${prefix}_R2.clean.fq.gz \
  -index_prefix $genome_prefix -work_dir $work_dir -compress Y

echo "<<<<< Methylation Extraction ..."
# Methylation call from vg giraffe alignment result.
methylGrapher MethylCall -t $CPU -index_prefix $genome_prefix -work_dir $work_dir -cg_only N

echo "<<<<< Merge CpG..."
# Merge cytosine methylation call (graph.methyl) into CpG methylation call.
# During graph indexing, all CpG locations are identified and stored in a separate TSV file.
# The graph CpG locations are stored under {index_prefix}cpg.tsv, with CpG id and both cytosine location on graph coordinate.
# MergeCpG function will merge the cytosine methylation call (graph.methyl) into CpG methylation call using graph CpG id.
methylGrapher MergeCpG -index_prefix $genome_prefix -work_dir $work_dir

echo "<<<<< Simulate bisulfite conversion rate ..."
# Estimate conversion rate from methylation extraction. MethylCall must be executed before this step.
methylGrapher ConversionRate -index_prefix $genome_prefix -work_dir $work_dir

echo "<<<<< (OPTIONAL) sorting output ..."
sort -k1n -k2n $work_dir/graph.methyl -o $work_dir/graph.methyl.sort

echo "<<<<< (OPTIONAL) Clean-up ..."
rm $work_dir/${prefix}_R1.clean.fq.gz
rm $work_dir/${prefix}_R2.clean.fq.gz
rm $work_dir/*.gaf
rm $work_dir/graph.methyl


#----------------
end=`date +%s`
runtime=$((end-start))
h=$(($runtime/3600))
hh=$(($runtime%3600))
m=$(($hh/60))
s=$(($hh%60))

echo "Start= $start"
echo "End= $end"
echo "Run time= $h:$m:$s"
echo "Done!"

