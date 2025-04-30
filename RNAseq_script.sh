#!/bin/sh

# Set DIR_WORKING, PATH_GENOME, PATH_ANNOTATION, DIR_INDEX, and SMAPLE_LIST (Choose one set when you run)

## P. formosa
DIR_WORKING=~/main/240827_RNAseq_Pformosa
PATH_GENOME=~/genome/Poecilia_formosa/GCF_000485575.1_Poecilia_formosa-5.1.2_genomic.fna
PATH_ANNOTATION=~/genome/Poecilia_formosa/genomic.gtf
DIR_INDEX=~/genome/index/240827Pfor_RNAseq
SMAPLE_LIST=(SRR13349690 SRR13349693 SRR13349684 SRR13349687 SRR13349692 SRR13349689 SRR13349683 SRR13349686)

# ## P. mexicana
# DIR_WORKING=/home/nakaharu/HDD5/machii/240827RNAseq_Pmex_Ov
# PATH_GENOME=/home/nakaharu/HDD5/machii/240827RNAseq_Pmex_Ov/genome/GCF_001443325.1_P_mexicana-1.0_genomic.fna
# PATH_ANNOTATION=/home/nakaharu/HDD5/machii/240827RNAseq_Pmex_Ov/genome/GCF_001443325.1_P_mexicana-1.0_genomic.gtf
# DIR_INDEX=/home/nakaharu/HDD5/machii/240827RNAseq_Pmex_Ov/genome/index
# SMAPLE_LIST=(SRR5224072 SRR5224073 SRR5224074 SRR5224075)

# ## B. splendens
# DIR_WORKING=~/main/240827_RNAseq_Bsplendens
# PATH_GENOME=~/genome/Betta_splendens/GCF_900634795.4_fBetSpl5.4_genomic.fna
# PATH_ANNOTATION=~/genome/Betta_splendens/genomic.gtf
# DIR_INDEX=~/genome/index/240827RNAseq_Bsplendens
# SMAPLE_LIST=(SRR18098736 SRR18098737 SRR18098738 SRR18098739 SRR18098740)

# ## O. latipes
# DIR_WORKING=~/main/240827_RNAseq_Olatipes
# PATH_GENOME=~/genome/Oryzias_latipes/GCF_002234675.1_ASM223467v1_genomic.fna
# PATH_ANNOTATION=~/genome/Oryzias_latipes/genomic.gtf
# DIR_INDEX=~/genome/index/240827Olat_RNAseq
# SMAPLE_LIST=(SRR23648299 SRR23648300 SRR23648301 SRR27475147 SRR27475148 SRR27475149)

# ## A. testudineus
# DIR_WORKING=~/main/240827_RNAseq_Atestudineus
# PATH_GENOME=~/genome/Anabas_testudineus/GCF_900324465.2_fAnaTes1.2_genomic.fna
# PATH_ANNOTATION=~/genome/Anabas_testudineus/genomic.gtf
# DIR_INDEX=~/genome/index/240827RNAseq_Atestudineus
# SMAPLE_LIST=(SRR8885647 SRR8885644 SRR8885645 SRR8885653 SRR8885654 SRR8885656 SRR8885657 SRR8885640)

# ## H. abdominalis
# DIR_WORKING=~/main/241110_RNAseq_Habd
# PATH_GENOME=~/main/241110_RNAseq_Habd/GCA_018466805.1_ZJU1.0_genomic.fna
# PATH_ANNOTATION=~/main/241110_RNAseq_Habd/Habdominalis_gene_annotation_substituted.gtf
# DIR_INDEX=~/main/241110_RNAseq_Habd/index
# SMAPLE_LIST=(SRR15058210 SRR15058211)

THREADS="16"
FASTP_THREADS="16"
FASTERQ_THREADS="16"
PIGZ_THREADS="16"

mkdir -p ${DIR_WORKING}
mkdir -p ${DIR_WORKING}/FASTP
mkdir -p ${DIR_WORKING}/STAR
mkdir -p ${DIR_WORKING}/RSEM


### STEP1: Download short reads and quality control
FASTP_start_time=$(date "+%Y-%m-%d %H:%M:%S")
echo "${FASTP_start_time}: Quality Check and Trimming with FASTP"
cd ${DIR_WORKING}/FASTP || exit

for sample in "${SMAPLE_LIST[@]}"; do
    prefetch \
        -O ./ \
        --max-size u \
        "${sample}"
    
    fasterq-dump \
        ./"${sample}" \
        -e "${FASTERQ_THREADS}" \
        --split-files
    
    pigz -p "${PIGZ_THREADS}" *.fastq

    fastp \
        -i ${DIR_WORKING}/FASTP/${sample}_1.fastq.gz \
        -I ${DIR_WORKING}/FASTP/${sample}_2.fastq.gz \
        -o ${sample}_trim_1.fastq.gz \
        -O ${sample}_trim_2.fastq.gz \
        -h ${sample}_fastp.html \
        -w ${FASTP_THREADS} 
done

### STEP2: Making star index and run star
mkdir -p ${DIR_INDEX}
STAR_Index_start_time=$(date "+%Y-%m-%d %H:%M:%S")
echo "${STAR_Index_start_time}: make index of STAR"
cd "${DIR_INDEX}" || exit

STAR \
    --runThreadN ${THREADS} \
    --runMode genomeGenerate \
    --genomeFastaFiles ${PATH_GENOME} \
    --genomeDir ${DIR_INDEX} \
    --limitGenomeGenerateRAM 100000000000 \
    --genomeSAindexNbases 13 \
    --sjdbGTFfile ${PATH_ANNOTATION}

STAR_start_time=$(date "+%Y-%m-%d %H:%M:%S")
echo "${STAR_start_time}: run STAR"

for sample in "${SMAPLE_LIST[@]}"; do
    mkdir -p "${DIR_WORKING}"/STAR/"${sample}"
    cd "${DIR_WORKING}"/STAR/"${sample}"

    STAR \
        --runThreadN ${THREADS} \
        --genomeDir ${DIR_INDEX} \
        --readFilesIn ${DIR_WORKING}/FASTP/${sample}_trim_1.fastq.gz ${DIR_WORKING}/FASTP/${sample}_trim_2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM
done

## STEP3: Making rsem index and run rsem-calculate-expression
rsem-prepare-reference \
    -p ${THREADS} \
    --gtf ${PATH_ANNOTATION} \
    ${PATH_GENOME} \
    ${DIR_INDEX}/RSEM

for sample in "${SMAPLE_LIST[@]}"; do
    mkdir -p "${DIR_WORKING}"/RSEM/"${sample}"
    cd "${DIR_WORKING}"/RSEM/"${sample}"
    
    rsem-calculate-expression \
        -p ${THREADS} \
        --alignments \
        --no-bam-output \
        --append-names \
        --bam ${DIR_WORKING}/STAR/${sample}/${sample}Aligned.toTranscriptome.out.bam \
        --strandedness reverse \
        --paired-end \
        ${DIR_INDEX}/RSEM \
        ${sample}_rsem
done
