#!/bin/bash
#PBS -N pangenome_methy
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

echo "CPU: $CPU"
echo "N: $N"

path="/public1/home/guowl/Public3_GuoWL/Foxtailmillet_PanMethylation"
bowtie2_path="/public1/home/guowl/software/bowtie2-2.5.3-linux-x86_64"

# Process samples
for acc in `acc.lst`; do
    echo "Processing ${acc} ..."
    
    # Build genome index (commented out)
    # echo "Build index for ${acc} genome ..."
    bismark_genome_preparation --path_to_aligner ${bowtie2_path} \
                               --parallel 10 --verbose ${path}/00.Reference/${acc}

    # Filter raw reads (commented out)
    # echo "Filter raw WGBS reads for ${acc} ..."
    fastp -i ${path}/01.CleanReads/${acc}_R1.fq.gz \
          -I ${path}/01.CleanReads/${acc}_R2.fq.gz \
          -o ${path}/01.CleanReads/${acc}_1.clean.fq.gz \
          -O ${path}/01.CleanReads/${acc}_2.clean.fq.gz \
          -z 4 -q 20 -u 30

    # Bismark alignment (commented out)
    # echo "Mapping WGBS reads to ${acc} genome ..."
    bismark --genome ${path}/00.Reference/${acc} \
            --path_to_bowtie2 ${bowtie2_path} \
            -1 ${path}/01.CleanReads/${acc}_1.clean.fq.gz \
            -2 ${path}/01.CleanReads/${acc}_2.clean.fq.gz \
            -p 30 -o ${path}/02.BismarkOutput/${acc}

    # Deduplicate (commented out)
    # echo "Deduplicate Bismark results of ${acc} ..."
    deduplicate_bismark -p --bam ${path}/02.BismarkOutput/${acc}/${acc}_1.clean_bismark_bt2_pe.bam \
                        --output_dir ${path}/02.BismarkOutput/${acc}

    # Methylation extraction
    echo "Extract cytosine sites of ${acc} ..."
    bismark_methylation_extractor -p --no_overlap --comprehensive --counts --bedGraph \
                                  --CX_context --cytosine_report --buffer_size 20G \
                                  --gzip --multicore 30 \
                                  --genome_folder ${path}/00.Reference/${acc} \
                                  ${path}/02.BismarkOutput/${acc}/${acc}_1.clean_bismark_bt2_pe.deduplicated.bam \
                                  -o ${path}/02.BismarkOutput/${acc}
    
    # Clean up (commented out)
    rm -f ${path}/01.CleanReads/${acc}_1.clean.fq.gz ${path}/01.CleanReads/${acc}_2.clean.fq.gz
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