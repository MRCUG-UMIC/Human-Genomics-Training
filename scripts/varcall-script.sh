#!/bin/bash -l
#SBATCH -J varcall
#SBATCH -e errors_%j.txt
#SBATCH --account=roche
#SBATCH --cpus-per-task=2
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8

# Load modules

module load bwa
module load samtools
module load bcftools
module load gatk


# Create bam folder 
mkdir -p bam
mkdir -p qc
mkdir -p gvcf
mkdir -p vcf

REF=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.fasta

DBSNP=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.dbsnp138.vcf.gz

KNOWN_INDELS=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.known_indels.vcf.gz

BRCA1=chr17:43041295-43173327 # Added a 3kb padding upstream and downstream

THREADS=8

while read SAMPLE; do
  R1=/mnt/roche/resources/datasets/human-genomics/short-read/${SAMPLE}_R1.fastq.gz
  R2=/mnt/roche/resources/datasets/human-genomics/short-read/${SAMPLE}_R2.fastq.gz
  RG="@RG\tID:${SAMPLE}.L1\tSM:${SAMPLE}\tLB:${SAMPLE}_lib1\tPL:ILLUMINA\tPU:FLOWCELL.L1"
	
  # Read alignment
  bwa mem -t ${THREADS} -R "${RG}" ${REF} ${R1} ${R2} \
    | samtools sort -@ ${THREADS} -o bam/${SAMPLE}.sorted.bam -

  samtools index bam/${SAMPLE}.sorted.bam

  # Mark duplicates
  gatk MarkDuplicates -I bam/${SAMPLE}.sorted.bam \
    -O bam/${SAMPLE}.dedup.bam -M qc/${SAMPLE}.dup_metrics.txt 

  samtools index bam/${SAMPLE}.dedup.bam

  # Base recalibration
  gatk BaseRecalibrator -I bam/${SAMPLE}.dedup.bam -R ${REF} \
    --known-sites ${DBSNP} --known-sites ${KNOWN_INDELS} \
    -L ${BRCA1} -O bam/${SAMPLE}.recal.table

  gatk ApplyBQSR -I bam/${SAMPLE}.dedup.bam -R ${REF} \
    --bqsr-recal-file bam/${SAMPLE}.recal.table -L ${BRCA1} \
    -O bam/${SAMPLE}.recal.bam

  # samtools index bam/${SAMPLE}.recal.bam

  # Determine variant sites
  gatk HaplotypeCaller -R ${REF} -I bam/${SAMPLE}.recal.bam \
    -L ${BRCA1} -ERC GVCF -O gvcf/${SAMPLE}.g.vcf.gz

done < samples.txt

# Single sample genotyping (e.g. HG02562)
gatk GenotypeGVCFs \
  -R ${REF} \
  -V gvcf/HG02562.g.vcf.gz \
  -L ${BRCA1} \
  -O vcf/HG02562.raw.vcf.gz

# Joint calling

## Create sample-gvcf map: sample-name <TAB> path-to-gvcf
awk '{print $1 "\tgvcf/" $1 ".g.vcf.gz"}' samples.txt > cohort.sample_map

## Option 1: Using GenomicsDBImport

gatk GenomicsDBImport \
  --genomicsdb-workspace-path genomicsdb_brca1 \
  --sample-name-map cohort.sample_map \
  -L ${BRCA1} \
  --reader-threads ${THREADS}


gatk GenotypeGVCFs \
 -R ${REF} \
 -V gendb://genomicsdb_brca1 \
 -O vcf/cohort.joint1.vcf.gz \
 --genomicsdb-shared-posixfs-optimizations

## Remove the database to save space
rm -rf genomicsdb_brca1

## Option 2: Using CombineGVCFs

awk '{print "-V", $2}' cohort.sample_map > gvcf_inputs.args


gatk CombineGVCFs \
    -R ${REF} \
    --arguments_file gvcf_inputs.args \
    -L ${BRCA1} \
    -O gvcf/cohort.g.vcf.gz

gatk GenotypeGVCFs \
    -R ${REF} \
    -V gvcf/cohort.g.vcf.gz \
    -O vcf/cohort.joint2.vcf.gz
