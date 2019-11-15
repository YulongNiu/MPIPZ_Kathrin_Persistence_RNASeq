#########################alignment##############################
date

REF_PATH=/netscratch/dep_psl/grp_rgo/yniu/ref/Kathrin_SynCom
CLEAN_PATH=/netscratch/dep_psl/grp_rgo/yniu/KathrinPersistence/clean_data
ALIGN_PATH=/netscratch/dep_psl/grp_rgo/yniu/KathrinPersistence/align_data_plantsyncom

KALLISTO_PATH=/home/yniu/Biotools/kallisto_v0.46.1
HISAT2_PATH=/home/yniu/Biotools/hisat2-2.1.0
SAMTOOL_PATH=/home/yniu/Biotools/samtools_1.9/bin
RM_PATH=/bin/rm

CORENUM=120

cd ${REF_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plantsyncom~~~~~~~~~~~~~~~~~~~~~~~~~~
fq=($(ls ${CLEAN_PATH} | grep .fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    SPECIES=${i%_*}
    REF_INDEX_PATH=${REF_PATH}/${SPECIES}

    echo "===================================="
    echo "Kallisto using ${SPECIES} cDNA for ${i}."
    ${KALLISTO_PATH}/kallisto quant \
                    -t ${CORENUM} \
                    -i ${REF_INDEX_PATH}/${SPECIES}.kindex \
                    -o ${ALIGN_PATH}/${i}_plantsyncom_kallisto \
                    ${CLEAN_PATH}/${i}_R1.fq.gz ${CLEAN_PATH}/${i}_R2.fq.gz

    echo "HISAT2 using ${SPECIES} genome for ${i}."
    ${HISAT2_PATH}/hisat2 -p ${CORENUM} \
                  -x ${REF_INDEX_PATH}/${SPECIES}ht2index/genome \
                  -1 ${CLEAN_PATH}/${i}_R1.fq.gz \
                  -2 ${CLEAN_PATH}/${i}_R2.fq.gz \
                  -S ${ALIGN_PATH}/${i}_plantsyncom_hisat2.sam

    ${SAMTOOL_PATH}/samtools sort \
                   -@ ${CORENUM} \
                   -o ${ALIGN_PATH}/${i}_plantsyncom_hisat2.bam \
                   ${ALIGN_PATH}/${i}_plantsyncom_hisat2.sam

    ${SAMTOOL_PATH}/samtools index ${ALIGN_PATH}/${i}_plantsyncom_hisat2.bam

    ${RM_PATH} ${ALIGN_PATH}/${i}_plantsyncom_hisat2.sam
    echo "====================================="

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################

