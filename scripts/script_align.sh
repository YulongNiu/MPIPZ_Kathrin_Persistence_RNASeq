#########################alignment##############################
date

REF_PATH=/netscratch/dep_psl/grp_rgo/yniu/ref
CLEAN_PATH=/netscratch/dep_psl/grp_rgo/yniu/KathrinPersistence/clean_data
ALIGN_PATH=/netscratch/dep_psl/grp_rgo/yniu/KathrinPersistence/align_data

KALLISTO_PATH=/home/yniu/Biotools/kallisto_v0.46.1
HISAT2_PATH=/home/yniu/Biotools/hisat2-2.1.0
SAMTOOL_PATH=/opt/share/software/bin
RM_PATH=/bin/rm
MOVE_PATH=/bin/mv

CORENUM=120

cd ${REF_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1stadd~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fq=($(ls ${CLEAN_PATH} | grep .fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    if [ ${i%%_*} = "C" ]; then
        SPECIES='ath'
    elif [ ${i%%_*} = "L" ]; then
        SPECIES='lotus'
    fi

    echo "===================================="
    echo "Kallisto using ${SPECIES} cDNA for ${i}."
    ${KALLISTO_PATH}/kallisto quant \
                    -t ${CORENUM} \
                    -i ${REF_PATH}/${SPECIES}/${SPECIES}.kindex \
                    -o ${ALIGN_PATH}/${i}_${SPECIES}_kallisto \
                    ${CLEAN_PATH}/${i}_R1.fq.gz ${CLEAN_PATH}/${i}_R2.fq.gz

    echo "HISAT2 using ${SPECIES} genome for ${i}."
    ${HISAT2_PATH}/hisat2 -p ${CORENUM} \
                  -x ${REF_PATH}/${SPECIES}/${SPECIES}ht2index/genome \
                  -1 ${CLEAN_PATH}/${i}_R1.fq.gz \
                  -2 ${CLEAN_PATH}/${i}_R2.fq.gz \
                  -S ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools sort \
                   -@ ${CORENUM} \
                   -o ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam \
                   ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools index ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam

    ${RM_PATH} ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam
    echo "====================================="

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################

