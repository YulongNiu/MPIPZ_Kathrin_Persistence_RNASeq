
## originally by Yulong Niu
## yulong.niu@hotmail.com

date

RNAMMER_PATH=/opt/share/software/bin
SEQTK_PATH=/home/yniu/Biotools/seqtk

LOTUSGENOME_PATH=/netscratch/dep_psl/grp_rgo/yniu/ref/lotus_gifu_collaborator_v1.2
OUTPUT_PATH=/netscratch/dep_psl/grp_rgo/yniu/KathrinPersistence/results/RNAPred

cd ${OUTPUT_PATH}

## extract sequences
${SEQTK_PATH}/seqtk subseq ${LOTUSGENOME_PATH}/LjGifu1.1_pseudomol.fa nucleolus.lst > nucleolus.fa

${SEQTK_PATH}/seqtk subseq ${LOTUSGENOME_PATH}/LjGifu1.1_pseudomol.fa mitochondria.lst > mitochondria.fa

${SEQTK_PATH}/seqtk subseq ${LOTUSGENOME_PATH}/LjGifu1.1_pseudomol.fa chloroplast.lst > chloroplast.fa

## predict ribosomal RNA
${RNAMMER_PATH}/rnammer -S euk -m lsu,ssu,tsu -xml nucleolus.xml -gff nucleolus.gff -h nucleolus.hmmreport < nucleolus.fa

${RNAMMER_PATH}/rnammer -S bac -m lsu,ssu,tsu -xml mitochondria.xml -gff mitochondria.gff -h mitochondria.hmmreport < mitochondria.fa

${RNAMMER_PATH}/rnammer -S bac -m lsu,ssu,tsu -xml chloroplast.xml -gff chloroplast.gff -h chloroplast.hmmreport < chloroplast.fa

date
