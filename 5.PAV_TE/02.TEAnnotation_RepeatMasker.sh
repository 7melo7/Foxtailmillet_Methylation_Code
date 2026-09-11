##/usr/bin/bash
#!/bin/bash
#PBS -N makeRepDB
#PBS -l nodes=1:ppn=5
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

path="/public3/data/Chenlab_Data/Share_Bigdata/Guowenlei/SitalicaData/RepMask"

${acc}="C12"

BuildDatabase -name ${acc} -engine ncbi ${acc}.genome.fa
RepeatModeler -pa 5 -database ${acc} -engine ncbi

RepeatMasker -lib ${acc}/${acc}-families.fa -gff -html -poly -pa 30 -dir ${acc} ${acc}.genome.fa

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"