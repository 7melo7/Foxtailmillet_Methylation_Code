##/usr/bin/bash
#!/bin/bash
#PBS -N MUMMER_SNP_Vcf
#PBS -l nodes=1:ppn=30
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

for i in `cat acc.lst.l`;do
    nucmer -t 30 -p $i Sitalica.Yugu1.fa ../supplementary_data_01_pangenome_sequences/$i.fa
    delta-filter -1 -q -r $i.delta > $i.qr.filter
    show-snps -Clr -T -x 1 delta_files/$i.qr.filter > mummer_snps/$i.qr.snps
done

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

