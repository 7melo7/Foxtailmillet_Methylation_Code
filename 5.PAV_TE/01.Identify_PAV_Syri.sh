##/usr/bin/bash
#!/bin/bash
#PBS -N SyriCompar
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

local_path="/public1/home/guowl/Foxtailmillet_population/00.PAV2Graphgenome/VariantsFromSyri"
ref_path="/public1/home/guowl/Foxtailmillet_population/00.Reference/Sitalica_Yugu1_Genome"
alt_path="/public1/home/guowl/Foxtailmillet_population/00.Pangenome/PangenomeSequence"


syri="/public1/home/guowl/miniconda3/envs/syri/bin/syri"
plotsr="/public1/home/guowl/miniconda3/envs/syri/bin/plotsr"

source activate syri

for alt in `cat 110.exlC1`;do
    if [ ! -d ${local_path}/${alt} ]; then
        mkdir ${local_path}/${alt}
    else
        rm ${local_path}/${alt}/*
    fi
    
    minimap2 -t 20 -ax asm5 --eqx ${ref_path}/Sitalica_Yugu1.genomeV2.fa ${alt_path}/${alt}.fa > ${local_path}/${alt}/${alt}.sam
    syri -c ${local_path}/${alt}/${alt}.sam -r ${ref_path}/Sitalica_Yugu1.genomeV2.fa -q ${alt_path}/${alt}.fa -k -F S \
        --prefix ${alt}.
    mv ${alt}.* ${local_path}/${alt}/
done


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"