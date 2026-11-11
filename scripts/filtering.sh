#!/bin/bash -l
#SBATCH -J var-filter
#SBATCH -e errors_%j.txt
#SBATCH --account=roche
#SBATCH --cpus-per-task=2
#SBATCH --mem=10GB

REF=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.fasta

DBSNP=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.dbsnp138.vcf.gz

KNOWN_INDELS=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.known_indels.vcf.gz

BRCA1=chr17:43041295-43173327 # Added a 3kb padding upstream and downstream

THREADS=8

# SNVs
gatk SelectVariants -R ${REF} -V vcf/cohort.raw.vcf.gz \
  --select-type-to-include SNP -O vcf/cohort.snps.vcf.gz

gatk VariantFiltration -R ${REF} -V vcf/cohort.snps.vcf.gz \
  --filter-expression "QD < 2.0"              --filter-name "QD2" \
  --filter-expression "FS > 60.0"             --filter-name "FS60" \
  --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
  -O vcf/cohort.snps.filtered.vcf.gz

# Indels
gatk SelectVariants -R ${REF} -V vcf/cohort.raw.vcf.gz \
  --select-type-to-include INDEL -O vcf/cohort.indels.vcf.gz

gatk VariantFiltration -R ${REF} -V vcf/cohort.indels.vcf.gz \
  --filter-expression "QD < 2.0"               --filter-name "QD2" \
  --filter-expression "FS > 200.0"             --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  --filter-expression "SOR > 10.0"             --filter-name "SOR10" \
  -O vcf/cohort.indels.filtered.vcf.gz

# Merge + keep PASS
gatk MergeVcfs \
  -I vcf/cohort.snps.filtered.vcf.gz \
  -I vcf/cohort.indels.filtered.vcf.gz \
  -O vcf/cohort.filtered.vcf.gz

bcftools view -f PASS vcf/cohort.filtered.vcf.gz -Oz -o vcf/cohort.pass.vcf.gz

bcftools index -t vcf/cohort.pass.vcf.gz

