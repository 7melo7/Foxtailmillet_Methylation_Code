#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=2000G
##SBATCH --mem-per-cpu=100G
#SBATCH -t 1000:00:00
#SBATCH -p amd_2T

start=`date +%s`

#----------------
# create work_dir
if [ ! -e ./work_dir ]; then
  mkdir ./work_dir
fi

# create out_dir
if [ ! -e ./out_dir ]; then
  mkdir ./out_dir
fi

# remove job_store
#if [ -e ./job_store ]; then
#  rm -rf ./job_store
#fi

#if [ ! -e ${prefix} ]; then
  echo "Loading ENV ..."
#  source activate
#  conda deactivate   
#  conda activate cactus
   
  CONDA_BASE=$(conda info --base)
  source ${CONDA_BASE}/etc/profile.d/conda.sh
  conda activate cactus
  export PATH=/public3/home/a6s001194/program/cactus-bin-v3.1.2/bin:$PATH
  export PYTHONPATH=/public3/home/a6s001194/program/cactus-bin-v3.1.2/lib:$PYTHONPATH
  export LD_LIBRARY_PATH=/public3/home/a6s001194/program/cactus-bin-v3.1.2/lib:$LD_LIBRARY_PATH
 
  echo "Run cactus-pangenome ..."
#  cactus-pangenome ./job_store ./genome.txt \
#    --outDir ./out_dir --outName Sit --reference Yugu01 \
#    --vcf --giraffe --gfa --gbz \
#    --maxCores 20 --maxMemory 2000G \
#    --workDir ./work_dir \
#    --binariesMode local

OLD_PATH=/public3/data/Chenlab_Data/Share_Bigdata/Liuyang/Project_SetariaMethyl/17.Minigraph-Cactus
NEW_PATH=$(pwd)

proot -b ${NEW_PATH}:${OLD_PATH} \
  cactus-pangenome ./job_store ./genome.txt --restart \
    --outDir ./out_dir --outName Sit --reference Yugu01 \
    --gfa --vcf \
    --maxCores 20 --maxMemory 1900G \
    --mgCores 20 --mgMemory 1900G \
    --mapCores 20 \
    --consCores 20 --consMemory 1900G \
    --indexCores 19 --indexMemory 1900G \
    --workDir ./work_dir \
    --binariesMode local
#fi


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

