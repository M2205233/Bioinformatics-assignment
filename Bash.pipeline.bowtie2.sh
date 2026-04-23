#!/usr/bin/env bash
# =============================================================================
#  WES Pipeline for NGS0001
# =============================================================================

set -euo pipefail

SAMPLE="NGS0001"
THREADS=4
BASE_DIR="/home/ubuntu/Bioinf_assignment/dnaseq_pipeline"
ANNOVAR_DIR="/home/ubuntu/annovar"

RAW_R1="${BASE_DIR}/data/untrimmed_fastq/${SAMPLE}.R1.fastq.gz"
RAW_R2="${BASE_DIR}/data/untrimmed_fastq/${SAMPLE}.R2.fastq.gz"
REFERENCE="${BASE_DIR}/data/reference/hg19.fa.gz"

echo "=== Starting Final WES Pipeline for ${SAMPLE} ==="

# Check ANNOVAR exists
if [ ! -f "${ANNOVAR_DIR}/convert2annovar.pl" ]; then
    echo "ERROR: convert2annovar.pl not found in ${ANNOVAR_DIR}"
    exit 1
fi

mkdir -p "${BASE_DIR}/results/fastqc"
mkdir -p "${BASE_DIR}/results/annotation"

# ========================= TRIMMING =========================
echo "[1/8] Trimmomatic..."
mkdir -p "${BASE_DIR}/data/trimmed_fastq"
trimmomatic PE -threads ${THREADS} -phred33 \
    "${RAW_R1}" "${RAW_R2}" \
    "${BASE_DIR}/data/trimmed_fastq/${SAMPLE}_R1_paired.fastq.gz" \
    "${BASE_DIR}/data/trimmed_fastq/${SAMPLE}_R1_unpaired.fastq.gz" \
    "${BASE_DIR}/data/trimmed_fastq/${SAMPLE}_R2_paired.fastq.gz" \
    "${BASE_DIR}/data/trimmed_fastq/${SAMPLE}_R2_unpaired.fastq.gz" \
    ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 TRAILING:25 MINLEN:50

# ========================= FASTQC =========================
echo "[2/8] FastQC on trimmed reads..."
fastqc -t ${THREADS} \
    "${BASE_DIR}/data/trimmed_fastq/${SAMPLE}_R1_paired.fastq.gz" \
    "${BASE_DIR}/data/trimmed_fastq/${SAMPLE}_R2_paired.fastq.gz" \
    -o "${BASE_DIR}/results/fastqc"

# ========================= ALIGNMENT =========================
echo "[3/8] Building Bowtie2 index..."

bowtie2-build "${REFERENCE}" "${BASE_DIR}/data/reference/hg19_bowtie2"

echo "[3.5/8] Bowtie2 alignment + sorting..."

bowtie2 -p ${THREADS} \
    --rg-id NGS0001 \
    --rg "SM:NGS0001\tPL:ILLUMINA\tLB:NGS0001\tPU:unit1" \
    -x "${BASE_DIR}/data/reference/hg19_bowtie2" \
    -1 "${BASE_DIR}/data/trimmed_fastq/${SAMPLE}_R1_paired.fastq.gz" \
    -2 "${BASE_DIR}/data/trimmed_fastq/${SAMPLE}_R2_paired.fastq.gz" \
    | samtools sort -@ ${THREADS} -o "${SAMPLE}_sorted.bam" -

samtools index "${SAMPLE}_sorted.bam"

rm -rf "${BASE_DIR}/data/trimmed_fastq"   # Free space

# ========================= MARK DUPLICATES =========================
echo "[4/8] Marking duplicates..."
picard MarkDuplicates I="${SAMPLE}_sorted.bam" O="${SAMPLE}_sorted_marked.bam" M=marked_dup_metrics.txt
rm -f "${SAMPLE}_sorted.bam"
samtools index "${SAMPLE}_sorted_marked.bam"

# ========================= QUALITY FILTER =========================
echo "[5/8] Quality filtering..."
samtools view -F 1796 -q 20 -@ ${THREADS} -b "${SAMPLE}_sorted_marked.bam" > "${SAMPLE}_sorted_filtered.bam"
rm -f "${SAMPLE}_sorted_marked.bam"
samtools index "${SAMPLE}_sorted_filtered.bam"

# ========================= VARIANT CALLING =========================
echo "[6/8] FreeBayes variant calling..."
samtools faidx "${REFERENCE%.gz}"
freebayes -f "${REFERENCE%.gz}" "${SAMPLE}_sorted_filtered.bam" > "${SAMPLE}.vcf"
bgzip "${SAMPLE}.vcf"
tabix -p vcf "${SAMPLE}.vcf.gz"
rm -f "${SAMPLE}_sorted_filtered.bam"

# ========================= FILTERING =========================
echo "[7/8] Filtering..."
vcffilter -f "QUAL > 20 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
    "${SAMPLE}.vcf.gz" > "${SAMPLE}_filtered.vcf"

bedtools intersect -header -wa \
    -a "${SAMPLE}_filtered.vcf" \
    -b "${BASE_DIR}/data/annotation.bed" > "${SAMPLE}_Chr22.vcf"

bgzip "${SAMPLE}_Chr22.vcf"
tabix -p vcf "${SAMPLE}_Chr22.vcf.gz"

# ========================= ANNOTATION (Fixed Path) =========================
echo "[8/8] Running ANNOVAR + snpEff..."

# ANNOVAR with correct path
${ANNOVAR_DIR}/convert2annovar.pl -format vcf4 "${SAMPLE}_Chr22.vcf.gz" > "${SAMPLE}_Chr22.avinput"

${ANNOVAR_DIR}/table_annovar.pl "${SAMPLE}_Chr22.avinput" ${ANNOVAR_DIR}/humandb/ \
    -buildver hg19 \
    -out "${SAMPLE}_annotated" \
    -remove \
    -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro \
    -operation g,g,f,f,f \
    -otherinfo -nastring . -csvout

# snpEff
snpEff hg19 "${SAMPLE}_filtered.vcf" > "${SAMPLE}_snpeff.vcf"
#snpEff annotation - hg19 v5 database unavailable due to network restrictions
# Command would be: snpEff hg19 ${SAMPLE}_filtered.vcf > ${SAMPLE}_snpeff.vcf
# snpEff hg19 "${SAMPLE}_filtered.vcf" > "${SAMPLE}_snpeff.vcf"

# Prioritisation
grep -i "exonic" "${SAMPLE}_annotated.hg19_multianno.csv" > "${BASE_DIR}/results/${SAMPLE}_exonic.csv" 2>/dev/null || true

echo "=== Pipeline completed successfully for ${SAMPLE} ==="

