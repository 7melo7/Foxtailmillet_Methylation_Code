#!/bin/bash
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=600G
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
#echo "Loading ENV ..."
#CONDA_BASE=$(conda info --base)
#source ${CONDA_BASE}/etc/profile.d/conda.sh
#conda activate Ixchel
 
#----------------
gfa=Sit60.d3.gfa

prefix=Sit60.d3.vg_srna

mkdir vg_tmp

#----------------
echo "Run vg autoindex ..."
# Create spliced pangenome graph and indexes for vg map
vg autoindex -t $CPU --target-mem 550G \
  --tmp-dir ./vg_tmp \
  --prefix $prefix \
  --workflow map \
  --gfa $gfa


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

