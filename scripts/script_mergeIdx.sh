date

REF_PATH=/netscratch/dep_psl/grp_rgo/yniu/ref/Kathrin_SynCom

KALLISTO_PATH=/home/yniu/Biotools/kallisto_v0.46.1
HISAT2_PATH=/home/yniu/Biotools/hisat2-2.1.0
SAMTOOLS_PATH=/home/yniu/Biotools/samtools_1.9/bin

CORENUM=40

cd ${REF_PATH}

condi=($(ls))

for i in ${condi[@]}; do

    echo "================================"
    echo "Kallisto index for ${i}"
    ${KALLISTO_PATH}/kallisto index \
                    -i ${REF_PATH}/${i}/${i}.kindex \
                    ${REF_PATH}/${i}/${i}_cDNA.fasta

    echo "HISAT2 index for ${i}"
    ${SAMTOOLS_PATH}/samtools faidx ${REF_PATH}/${i}/${i}_genome.fasta

    mkdir ${REF_PATH}/${i}/${i}ht2index
    cd ${REF_PATH}/${i}/${i}ht2index
    ${HISAT2_PATH}/hisat2-build -f \
                  ${REF_PATH}/${i}/${i}_genome.fasta\
                  genome

done
