#/usr/bin/bash
#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=20b
#PBS -l walltime=40:00:00
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

vcf=merged.vcf.gz
vcftools=/public1/home/guowl/software/vcftools_0.1.13/bin/vcftools


pops=("w" "l" "c")
window_size=(10000 20000 50000 100000)
window_step=(1000 2000 5000 10000)


# calculate the pi with different levels of pop, window size
for i in ${pops[*]};do
    for j in {0..3};do
        $vcftools --gzvcf $vcf --window-pi ${window_size[${j}]} --window-pi-step ${window_step[${j}]} \
                  --keep acc.lst.${i} --out pi/${i}.win_${window_size[${j}]}
    done
done  

# calculate the Fst between the groups for each chromosome
for pair in w,l l,c ;do
    str=(${pair//,/ })
    pop1=${str[0]}
    pop2=${str[1]}
    for i in {0..3};do
        $vcftools --gzvcf $vcf --weir-fst-pop acc.lst.${pop1} \
        --weir-fst-pop acc.lst.${pop2} \
        --out fst/${pop1}_vs_${pop2}.win_${window_size[${i}]} \
        --fst-window-size ${window_size[${i}]} --fst-window-step ${window_step[${i}]}
    done
done 

# calculate the tajimaD for each group

for i in ${pops[*]};do
    for j in {0..3};do
        $vcftools --gzvcf $vcf --TajimaD ${window_size[${j}]} --keep acc.lst.${i} --out tajimaD/${i}.win_${window_size[${j}]}
    done
done


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

