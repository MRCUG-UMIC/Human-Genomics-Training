REF=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.fasta

DBSNP=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.dbsnp138.vcf.gz

KNOWN_INDELS=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.known_indels.vcf.gz

BRCA1=chr17:43041295-43173327 # Added a 3kb padding upstream and downstream

THREADS=8


# Stats
bcftools stats -s - vcf/cohort.pass.vcf.gz > qc/cohort.bcftools_stats.txt


