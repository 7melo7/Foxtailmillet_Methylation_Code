##/usr/bin/bash
#!/bin/bash
#PBS -N DMR
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

program1="/public1/home/guowl/software/metilene_v0.2-9/metilene"
program2="/public1/home/guowl/software/metilene_v0.2-9/metilene_output.pl"


context="CG"

$program1 -M 300 -m 8 -d 0.4 -t 20 -f 1 -X 12 -Y 12 -a Q_${context} -b L_${context} Metilene.Q_vs_L.${context}.fixed.txt > Metilene.Q_vs_L.${context}.fixed.output
perl $program2 -q Metilene.Q_vs_L.${context}.fixed.output -p 0.01 -c 8 -a Q_${context} -b L_${context} -o DMR.${context}.Q_vs_L.fixed

$program1 -M 300 -m 8 -d 0.4 -t 20 -f 1 -a Q_${context} -b L_${context} Metilene.Q_vs_L.${context}.onegroup.txt > Metilene.Q_vs_L.${context}.onegroup.output
perl $program2 -q Metilene.Q_vs_L.${context}.onegroup.output -p 0.01 -c 8 -a Q_${context} -b L_${context} -o DMR.${context}.Q_vs_L.onegroup

context="CHG"

$program1 -M 300 -m 8 -d 0.4 -t 20 -f 1 -X 12 -Y 12 -a Q_${context} -b L_${context} Metilene.Q_vs_L.${context}.fixed.txt > Metilene.Q_vs_L.${context}.fixed.output
perl $program2 -q Metilene.Q_vs_L.${context}.fixed.output -p 0.01 -c 8 -a Q_${context} -b L_${context} -o DMR.${context}.Q_vs_L.fixed

$program1 -M 300 -m 8 -d 0.4 -t 20 -f 1 -a Q_${context} -b L_${context} Metilene.Q_vs_L.${context}.onegroup.txt > Metilene.Q_vs_L.${context}.onegroup.output
perl $program2 -q Metilene.Q_vs_L.${context}.onegroup.output -p 0.01 -c 8 -a Q_${context} -b L_${context} -o DMR.${context}.Q_vs_L.onegroup


context="CHH"

$program1 -M 300 -m 8 -d 0.2 -t 20 -f 1 -X 12 -Y 12 -a Q_${context} -b L_${context} Metilene.Q_vs_L.${context}.fixed.txt > Metilene.Q_vs_L.${context}.fixed.output
perl $program2 -q Metilene.Q_vs_L.${context}.fixed.output -p 0.01 -c 8 -a Q_${context} -b L_${context} -o DMR.${context}.Q_vs_L.fixed

$program1 -M 300 -m 8 -d 0.2 -t 20 -f 1 -a Q_${context} -b L_${context} Metilene.Q_vs_L.${context}.onegroup.txt > Metilene.Q_vs_L.${context}.onegroup.output
perl $program2 -q Metilene.Q_vs_L.${context}.onegroup.output -p 0.01 -c 8 -a Q_${context} -b L_${context} -o DMR.${context}.Q_vs_L.onegroup



end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"