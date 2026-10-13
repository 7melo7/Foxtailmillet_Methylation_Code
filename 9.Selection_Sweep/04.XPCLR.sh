##/usr/bin/bash
#!/bin/bash
#PBS -N XPCLR
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

VCF=merged.vcf.gz

pop1=w
pop2=l
chr=${N}

window_size=(10000 20000 50000 100000)
window_step=(1000 2000 5000 10000)

for pair in w,l l,c ;do
    str=(${pair//,/ })
    pop1=${str[0]}
    pop2=${str[1]}
    for i in {0..3};do
        xpclr --out xpclr/${pop2}_vs_${pop1}.chr${chr}.size${window_size[${i}]}.xpclr.txt --input $VCF \
              --samplesA acc.lst.${pop1} \
              --samplesB acc.lst.${pop2} \
              --chr scaffold_${chr} --maxsnps 300 --size ${window_size[${i}]} --step  ${window_step[${i}]}
    done
done

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

