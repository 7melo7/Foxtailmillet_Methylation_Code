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
#echo "Loading ENV ..."
#CONDA_BASE=$(conda info --base)
#source ${CONDA_BASE}/etc/profile.d/conda.sh
#conda activate Ixchel
 
#----------------
fq=`ls /publicfs10/fs10-m9/home/m9s003656/Project_SetariaMe/0.clean_fastq/*.clean.fq.gz | head -n $N | tail -n 1`
fq_dir=$(dirname $fq)
fq_prefix=$(basename $fq .clean.fq.gz)

gfa=Sit60.d3.gfa
index_prefix=Sit60.d3.vg_srna

if [ ! -d output ]; then mkdir output; fi

#----------------
echo "CPU= $CPU"
echo "READS_PREFIX= $fq_prefix"

#----------------
if [ ! -e ./output/$fq_prefix.srna.bam ]; then
  echo "Run vg map ..."
  # Map small RNA reads using vg map
  vg map -t $CPU \
    -f $fq \
    -x $index_prefix.xg \
    -g $index_prefix.gcsa \
    --min-mem 11 --min-ident 1.0 --max-multimaps 2 \
    --surject-to bam \
    > ./output/$fq_prefix.srna.bam
fi


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

