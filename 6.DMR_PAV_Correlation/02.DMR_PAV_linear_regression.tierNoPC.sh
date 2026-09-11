#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
##SBATCH --mem=20G
##SBATCH --mem-per-cpu=20G
#SBATCH -t 100:00:00
#SBATCH -p amd_m9_768
#SBATCH --chdir=./

start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then CPU=2; fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then N=1; fi

#----------------
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate base

#----------------

PART="10"
python dmr_pav_analysis_ols_tier_noPC.py --pairwise-file DMR_overlap_PAV.10kbp.pairwise --output-dir results_1000_permu_tierNoPC --n-permutations 1000 --n-threads 20

#----------------
#echo "Run samtools stats ..."
#samtools stats -@ $CPU $bam > $bam.stats

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
