
process ALIGNMENT {

    tag "$sample"

    input:
    tuple val(sample), path(reads)
    path ref_dir

    output:
    tuple val(sample), path("${sample}.bam")

    script:
    def (r1, r2) = reads

    """
    set -euo pipefail

    bwa mem -t ${task.cpus} \
        -R '@RG\\tID:HG02562.L1\\tSM:HG02562\\tLB:HG02562_lib1\\tPL:ILLUMINA\\tPU:FLOWCELL.L1' \
        ${ref_dir}/Homo_sapiens_assembly38.fasta \
        ${r1} ${r2} \
        | samtools view -@ ${task.cpus} -b -o ${sample}.bam -
    """
}


process SORT_BAM {

    tag "$sample"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample),
          path("${sample}.sorted.bam"),
          path("${sample}.sorted.bam.bai")

    script:
    """
    set -euo pipefail

    samtools sort -@ ${task.cpus} \
        -o ${sample}.sorted.bam ${bam}

    samtools index -@ ${task.cpus} \
        ${sample}.sorted.bam
    """
}


process MARK_DUPLICATES {

    tag "$sample"

    publishDir "${params.outdir}/qc",
        pattern: '*.dup_metrics.txt',
        mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample),
          path("${sample}.dedup.bam"),
          path("${sample}.dedup.bam.bai"),
          emit: bam

    path("${sample}.dup_metrics.txt"), emit: metrics

    script:
    """
    set -euo pipefail

    gatk MarkDuplicates \
        -I ${bam} \
        -O ${sample}.dedup.bam \
        -M ${sample}.dup_metrics.txt

    samtools index -@ ${task.cpus} \
        ${sample}.dedup.bam
    """
}


process BASE_RECALIBRATOR {

    tag "$sample"

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_dir

    output:
    tuple val(sample),
          path("${sample}.dedup.bam"),
          path("${sample}.dedup.bam.bai"),
          path("${sample}.recal.table")

    script:
    """
    set -euo pipefail

    gatk BaseRecalibrator \
        -I ${bam} \
        -R ${ref_dir}/Homo_sapiens_assembly38.fasta \
        --known-sites ${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
        --known-sites ${ref_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
        -L '${params.interval}' \
        -O ${sample}.recal.table
    """
}


process APPLY_BQSR {

    tag "$sample"

    publishDir "${params.outdir}/bam",
        pattern: '*.recal.bam*',
        mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai), path(recal_table)
    path ref_dir

    output:
    tuple val(sample),
          path("${sample}.recal.bam"),
          path("${sample}.recal.bam.bai")

    script:
    """
    set -euo pipefail

    gatk ApplyBQSR \
        -I ${bam} \
        -R ${ref_dir}/Homo_sapiens_assembly38.fasta \
        --bqsr-recal-file ${recal_table} \
        -L '${params.interval}' \
        -O ${sample}.recal.bam

    samtools index -@ ${task.cpus} \
        ${sample}.recal.bam
    """
}

