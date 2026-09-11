##/usr/bin/bash
#!/bin/bash
#PBS -N Cluster
#PBS -l nodes=1:ppn=20
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

cd-hit-est -i te.fa -o te.cluster.txt -sc 1 -c 0.8 -aL 0.8 -aS 0.8 \
           -n 5 -M 0 -T 20

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"