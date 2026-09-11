#:wq/usr/bin/bash
#!/bin/bash
#PBS -N EDTA
#PBS -l nodes=1:ppn=30
#PBS -l mem=100gb
#PBS -l walltime=9999:00:00
#PBS -d ./
#PBS -j oe
#PBS -V

cd $PBS_O_WORKDIR

start=`date +%s`

CPU=$PBS_NP
if [ ! $CPU ]; then
   CPU=2
fi

N=$PBS_ARRAYID
if [ ! $N ]; then
    N=1
fi

echo "$CPU"
echo "$N"

source activate EDTA

edta="/public1/home/guowl/miniconda3/envs/EDTA/share/EDTA/EDTA.pl"

perl ${edta} --genome ../Sitalica_Yugu1.genomeV2.fa --cds ../Sitalica_Yugu1.CDS.fa --sensitive 1 \
             --anno 1 --evaluate 0 --overwrite 1 -t 30 

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

