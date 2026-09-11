##/usr/bin/bash
#!/bin/bash
#PBS -N MergeC
#PBS -l nodes=1:ppn=1
#PBS -l walltime=9999:00:00
#PBS -d ./
#PBS -j oe
#PBS -V

start=`date +%s`

CPU=$PBS_NP
if [ ! $CPU ]; then
   CPU=2
fi

N=$PBS_ARRAYID
if [ ! $N ]; then
    N=1
fi

program="/public1/home/guowl/software/metilene_v0.2-9/metilene_input.pl"
bed_dir="/public1/home/guowl/Foxtailmillet_population/02.DNA_methylation_calls/CytosineSites_MethyGraph"

typ="CG"

group_Q=Q14_${typ}.bed,Q16_${typ}.bed,Q17_${typ}.bed,Q18_${typ}.bed,Q1_${typ}.bed,Q20_${typ}.bed,Q21_${typ}.bed,Q23_${typ}.bed,Q24_${typ}.bed,Q25_${typ}.bed,Q26_${typ}.bed,Q28_${typ}.bed,Q2_${typ}.bed,Q31_${typ}.bed,Q33_${typ}.bed,Q34_${typ}.bed,Q35_${typ}.bed,Q37_${typ}.bed,Q4_${typ}.bed,Q7_${typ}.bed

group_L=L12_${typ}.bed,L14_${typ}.bed,L16_${typ}.bed,L18_${typ}.bed,L1_${typ}.bed,L20_${typ}.bed,L22_${typ}.bed,L24_${typ}.bed,L26_${typ}.bed,L27_${typ}.bed,L29_${typ}.bed,L32_${typ}.bed,L33_${typ}.bed,L34_${typ}.bed,L36_${typ}.bed,L37_${typ}.bed,L38_${typ}.bed,L5_${typ}.bed,L6_${typ}.bed,L9_${typ}.bed

group_C=C12_${typ}.bed,C13_${typ}.bed,C14_${typ}.bed,C16_${typ}.bed,C19_${typ}.bed,C22_${typ}.bed,C25_${typ}.bed,C27_${typ}.bed,C28_${typ}.bed,C29_${typ}.bed,C2_${typ}.bed,C31_${typ}.bed,C32_${typ}.bed,C34_${typ}.bed,C35_${typ}.bed,C3_${typ}.bed,C5_${typ}.bed,C6_${typ}.bed,C8_${typ}.bed,C9_${typ}.bed


perl $program --in1 $group_Q --in2 $group_L --h1 Q_${typ} --h2 L_${typ} --out Metilene.Q_vs_L.${typ}.input
#sed 's/NA/0/g' ${dir}/Q_${typ}_L_${typ}.input >${dir}/metilene_Q_${typ}_L_${typ}.input
#rm ${dir}/Q_${typ}_L_${typ}.input

perl $program --in1 $group_L --in2 $group_C --h1 L_${typ} --h2 C_${typ} --out Metilene.L_vs_C.${typ}.input
#sed 's/NA/0/g'  ${dir}/L_${typ}_C_${typ}.input>${dir}/metilene_L_${typ}_C_${typ}.input
#rm ${dir}/L_${typ}_C_${typ}.input


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