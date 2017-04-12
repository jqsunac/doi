#!/bin/bash

source ../bin/global.sh


JOB_ROOT=${PROJECT_ROOT}/qc


FASTQ_REPORTS_BEFORE=${FASTQDIR}/reports
FASTQ_REPORTS_AFTER=${FASTQQCDIR}/reports


mkdir -p ${FASTQ_REPORTS_BEFORE}
mkdir -p ${FASTQ_REPORTS_AFTER}




cd ${FASTQDIR}

for rnaseq_lib in ${RNASEQ_LIBS[@]}; do
    ${BIN}/fastqc --nogroup --noextract -o ${FASTQ_REPORTS_BEFORE} ${FASTQDIR}/${rnaseq_lib}_R1.fq.gz
    ${BIN}/fastqc --nogroup --noextract -o ${FASTQ_REPORTS_BEFORE} ${FASTQDIR}/${rnaseq_lib}_R2.fq.gz
done



for rnaseq_lib in ${RNASEQ_LIBS[@]}; do
    java -Xmx4g -jar ${BIN}/Trimmomatic-0.36/trimmomatic-0.36.jar \
        PE -threads ${THREADS} -phred33 \
        ${FASTQDIR}/${rnaseq_lib}_R1.fq.gz ${FASTQDIR}/${rnaseq_lib}_R2.fq.gz \
        ${FASTQQCDIR}/${rnaseq_lib}_R1.fq.gz \
        ${FASTQQCDIR}/${rnaseq_lib}_unpaired_R1.fq.gz \
        ${FASTQQCDIR}/${rnaseq_lib}_R2.fq.gz \
        ${FASTQQCDIR}/${rnaseq_lib}_unpaired_R2.fq.gz \
        ILLUMINACLIP:${CWD}/adapters.txt:2:30:10 \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
done



cd ${FASTQQCDIR}

for rnaseq_lib in ${RNASEQ_LIBS[@]}; do
    ${BIN}/fastqc --nogroup --noextract -o ${FASTQ_REPORTS_AFTER} ${FASTQQCDIR}/${rnaseq_lib}_R1.fq.gz
    ${BIN}/fastqc --nogroup --noextract -o ${FASTQ_REPORTS_AFTER} ${FASTQQCDIR}/${rnaseq_lib}_R2.fq.gz
done



