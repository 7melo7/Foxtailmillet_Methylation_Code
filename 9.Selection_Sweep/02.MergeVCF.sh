##/usr/bin/bash
#!/bin/bash
#PBS -N mergeVcf
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

for i in `cat acc.lst`;do 
    python SNPs2VCF.py snps/$i.qr.snps $i
    mv $i.vcf vcf_files/$i.vcf
    bgzip vcf_files/$i.vcf
    bcftools index -t vcf_files/$i.vcf.gz
done 

bcftools merge vcf_files/*.vcf.gz -Oz -o merged.vcf.gz


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

