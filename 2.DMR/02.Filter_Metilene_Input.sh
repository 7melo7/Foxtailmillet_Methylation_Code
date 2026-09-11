##/usr/bin/bash
#!/bin/bash
#PBS -N FilterC
#PBS -l nodes=1:ppn=1
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


for prefix in `ls *CHH.input`;do
    echo $prefix
    python ../../script/filter_matrix_NA.py ${prefix} ${prefix%.input}
done



end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"