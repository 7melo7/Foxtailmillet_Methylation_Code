#!/bin/bash
#PBS -N bismark_statistic
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

echo "CPU: $CPU"
echo "N: $N"

path="/public1/home/guowl/Foxtailmillet_population/00.WGBS"

cytosine_statistic() {
    prefix=$1
    result_dir="Methylation_Results_Ref_Yugu1_Extract/${prefix}"
    base_name="${prefix}_R1.clean_bismark_bt2_pe.deduplicated"
    
    # Extract and filter chromosome data
    gunzip -c ${result_dir}/${base_name}.CX_report.txt.gz | grep -v contig > ${result_dir}/${prefix}_chr1-9.txt
    
    # Calculate p-values
    python scripts/methylated_p.py ${path}/nc_rates/C13.nc.txt \
           ${result_dir}/${prefix}_chr1-9.txt \
           ${result_dir}/${prefix}_chr1-9.pvalue.txt 4
    
    # FDR correction
    Rscript scripts/methylated_fdr.R ${result_dir}/${prefix}_chr1-9.pvalue.txt \
           ${result_dir}/${prefix}_chr1-9.pvalue.fdr.txt
    
    # Label cytosines and compress
    python scripts/label_cytosine.py ${result_dir}/${prefix}_chr1-9.pvalue.fdr.txt \
           ${result_dir}/${prefix}_cytosine.txt
    gzip ${result_dir}/${prefix}_cytosine.txt
    
    # Clean up intermediate files
    rm -f ${result_dir}/${prefix}_chr1-9.txt \
          ${result_dir}/${prefix}_chr1-9.pvalue.txt \
          ${result_dir}/${prefix}_chr1-9.pvalue.fdr.txt
}

# Process samples
samples="L29_10X L29_12X L29_14X L29_16X L29_18X"
for i in $samples; do
    cytosine_statistic $i
done

# Runtime statistics
end=`date +%s`
runtime=$((end-start))
h=$(($runtime/3600))
m=$((($runtime%3600)/60))
s=$(($runtime%60))

echo "Start: $start"
echo "End: $end"
echo "Run time: ${h}h ${m}m ${s}s"
echo "Done"