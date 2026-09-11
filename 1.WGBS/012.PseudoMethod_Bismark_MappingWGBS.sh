#!/bin/bash
#PBS -N BISMARK
#PBS -l nodes=1:ppn=30
#PBS -l walltime=360:00:00
#PBS -d ./
#PBS -j oe

start=`date +%s`

CPU=$PBS_NP
if [ ! $CPU ]; then
   CPU=2
fi

N=$PBS_ARRAYID
if [ ! $N ]; then
    N=1
fi

source activate python37

path="/public1/home/guowl/Foxtailmillet_population/00.WGBS"
pseudo_path="/public1/home/guowl/Foxtailmillet_population/01.Pseudogenome/Pseudo_Reference_Genomes4Bismark"

while read -r prefix; do
    # Quality control and trimming
    fastp -i ${path}/RawLibraryReads_Extract/${prefix}_R1.fq.gz \
          -I ${path}/RawLibraryReads_Extract/${prefix}_R2.fq.gz \
          -o ${path}/CleanReads_Extract/${prefix}_R1.clean.fq.gz \
          -O ${path}/CleanReads_Extract/${prefix}_R2.clean.fq.gz \
          -z 4 -q 20 -u 30

    # Create output directory
    mkdir -p Methylation_Results_Ref_Yugu1_Extract/${prefix}

    # Bismark alignment
    bismark --genome ${pseudo_path}/${prefix%_*} \
            --path_to_bowtie2 /public1/home/guowl/miniconda3/envs/python37/bin/ \
            -1 ${path}/CleanReads_Extract/${prefix}_R1.clean.fq.gz \
            -2 ${path}/CleanReads_Extract/${prefix}_R2.clean.fq.gz \
            -p 30 -o Methylation_Results_Ref_Yugu1_Extract/${prefix}

    # Deduplicate
    deduplicate_bismark -p --bam Methylation_Results_Ref_Yugu1_Extract/${prefix}/${prefix}_R1.clean_bismark_bt2_pe.bam \
                        --output_dir Methylation_Results_Ref_Yugu1_Extract/${prefix}

    # Methylation extraction
    bismark_methylation_extractor -p --no_overlap --comprehensive --counts --bedGraph \
                                  --CX_context --cytosine_report --buffer_size 10G \
                                  --gzip --parallel 10 \
                                  --genome_folder ${pseudo_path}/${prefix%_*} \
                                  Methylation_Results_Ref_Yugu1_Extract/${prefix}/${prefix}_R1.clean_bismark_bt2_pe.deduplicated.bam \
                                  -o Methylation_Results_Ref_Yugu1_Extract/${prefix}
done

# Runtime statistics
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