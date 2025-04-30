#!/bin/bash

# set SRR_LIST and WORK_DIR (Choose one set of WORK_DIR and SRR_LIST)
## P. formosa Liver
WORK_DIR=/home/nakaharu/HDD5/machii/241107_Trinity_Pformosa_Lv
SRR_LIST=(SRR13349683 SRR13349686 SRR13349689 SRR13349692)

# ## P. formosa Ovary
# WORK_DIR=/home/nakaharu/HDD5/machii/241107_Trinity_Pformosa_Ov
# SRR_LIST=(SRR13349684 SRR13349687 SRR13349690 SRR13349693)

# ## P. mexicana Liver
# WORK_DIR=/home/nakaharu/HDD5/machii/241107_Trinity_Pformosa_Lv
# SRR_LIST=(SRR5224072 SRR5224073 SRR5224074 SRR5224075)

# ## O. latipes Liver
# WORK_DIR=/home/nakaharu/HDD5/machii/241107_Trinity_Olatipes_Lv
# SRR_LIST=(SRR27475147 SRR27475148 SRR27475149)

# ## O. latipes Ovary
# WORK_DIR=/home/nakaharu/HDD5/machii/241107_Trinity_Olatipes_Ov
# SRR_LIST=(SRR23648299 SRR23648300 SRR23648301)

# ## H. abdominalis Ovary
# WORK_DIR=/home/nakaharu/HDD5/machii/241107_Trinity_Habdominalis_Ov
# SRR_LIST=(SRR15058210 SRR15058211)

# ## B. splendens Ovary
# WORK_DIR=/home/nakaharu/HDD5/machii/241107_Trinity_Bsplendens_Ov
# SRR_LIST=(SRR18098736 SRR18098737 SRR18098738 SRR18098739 SRR18098740)

# ## A. testudineus Liver
# WORK_DIR=/home/nakaharu/HDD5/machii/241107_Trinity_Atestudineus_Lv
# SRR_LIST=(SRR8885640 SRR8885644 SRR8885645 SRR8885647 SRR8885653 SRR8885654 SRR8885656 SRR8885657)

OUTPUT_DIR=${WORK_DIR}/out

THREADS="4"

STRAND_SPECIFICITY="RF"

date=`date '+%y%m%d_%H%M'`

mkdir -p ${WORK_DIR}
mkdir -p ${WORK_DIR}/fastq
mkdir -p ${WORK_DIR}/fastp
mkdir -p ${OUTPUT_DIR}


### STEP1:
echo "Downloading and converting SRA data..."
cd ${WORK_DIR}/fastq

for srr in "${SRR_LIST[@]}"; do
    prefetch \
        -O ./ \
        --max-size u \
        "${srr}"
    
    fasterq-dump \
        ./"${srr}" \
        -e "${THREADS}" \
        --split-files
    
    pigz \
        -p "${THREADS}" \
        *.fastq
done

### STEP2
echo "Quality checking and trimming with fastp..."
cd ${WORK_DIR}/fastp

for srr in "${SRR_LIST[@]}"; do
    fastp \
        -i ${WORK_DIR}/fastq/${srr}_1.fastq.gz \
        -I ${WORK_DIR}/fastq/${srr}_2.fastq.gz \
        -o ${srr}_trim_1.fastq.gz \
        -O ${srr}_trim_2.fastq.gz \
        -h ${srr}_fastp.html \
        -w ${THREADS} 
done

### STEP3
echo "Running Trinity for de novo assembly..."

left=""
right=""
for srr in "${SRR_LIST[@]}"; do
    if [ "$left" = "" ]; then
        left="${WORK_DIR}/fastp/${srr}_trim_1.fastq.gz"
        right="${WORK_DIR}/fastp/${srr}_trim_2.fastq.gz"
    else
        left="${left},${WORK_DIR}/fastp/${srr}_trim_1.fastq.gz"
        right="${right},${WORK_DIR}/fastp/${srr}_trim_2.fastq.gz"
    fi
done

# Run Trinity
Trinity \
    --seqType fq \
    --left ${left} \
    --right ${right} \
    --CPU ${THREADS} \
    --max_memory 100G \
    --SS_lib_type ${STRAND_SPECIFICITY} \
    --output ${OUTPUT_DIR}/trinity_${date} \
    --full_cleanup \
    --workdir ${OUTPUT_DIR}/trinity_tmp

echo "Pipeline completed successfully!"

