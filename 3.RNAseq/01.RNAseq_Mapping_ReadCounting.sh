#!/bin/bash
#PBS -N Sitalica_RNAseq
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
echo "Array ID: $N"

# Directories
raw_lp="/public3/data/Chenlab_Data/Foxtailmillet_RNAseq/Si_RNAseq_add1repfor9acc_fastq_cleandata"
output_lp="/public3/data/Chenlab_Data/Foxtailmillet_Replace"
bam_dir="${output_lp}/bam_data_0412"
clean_dir="${output_lp}/clean_fastq_0412"
summary_dir="/public1/home/guowl/methy/13_RNA-seq_0407/align_summary"


# Create output directories if needed
mkdir -p ${bam_dir} ${clean_dir} ${summary_dir}

# Get sample ID from array task
acc_lt=($(awk '{print $1}' rna_seq.ids))
acc_id=${acc_lt[${N}]}

echo "Start processing ${acc_id} RNA-seq data ..."

# QC and filtering (commented out)
fastp -i ${raw_lp}/${acc_id}_R1.fq.gz \
      -I ${raw_lp}/${acc_id}_R2.fq.gz \
      -o ${clean_dir}/${acc_id}_R1.clean.fq.gz \
      -O ${clean_dir}/${acc_id}_R2.clean.fq.gz \
      -z 4

# Alignment with bowtie2 (commented out)
bowtie2 -x /public1/home/guowl/dataset/Sitalica/Sitalica_312_v2.index \
        -1 ${clean_dir}/${acc_id}_R1.clean.fq.gz \
        -2 ${clean_dir}/${acc_id}_R2.clean.fq.gz \
        -p 5 -S ${bam_dir}/${acc_id}.sam \
        2> ${summary_dir}/${acc_id}.align.summary

# Convert SAM to BAM
echo "Converting SAM to BAM for ${acc_id} ..."
samtools sort -O bam -o ${bam_dir}/${acc_id}.bam ${bam_dir}/${acc_id}.sam

featureCounts -t exon -g gene_id -p -a test.gtf -T 30 -o read_counts.0412.gtf.txt ${output_lp}/bam_data_0412/*.sam




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