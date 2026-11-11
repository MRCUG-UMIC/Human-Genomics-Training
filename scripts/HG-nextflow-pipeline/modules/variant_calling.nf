// This file contains the complete variant-calling module: per-sample calling and
// the two alternative cohort joint-genotyping processes.
process VARIANT_CALLING {
    tag "$sample"
    publishDir "${params.outdir}/gvcf", pattern: '*.g.vcf.gz*', mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_dir

    output:
    tuple val(sample), path("${sample}.g.vcf.gz"), path("${sample}.g.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    gatk HaplotypeCaller -R ${ref_dir}/Homo_sapiens_assembly38.fasta -I ${bam} \\
      -L '${params.interval}' -ERC GVCF -O ${sample}.g.vcf.gz
    """
}

process JOINT_COMBINE {
    publishDir "${params.outdir}/vcf", pattern: 'cohort.joint.vcf.gz*', mode: 'copy'
    input:
    path gvcfs
    path ref_dir

    output:
    path 'cohort.joint.vcf.gz'
    path 'cohort.joint.vcf.gz.tbi'

    script:
    """
    set -euo pipefail
    : > gvcf_inputs.args
    for f in *.g.vcf.gz; do printf -- '-V %s\\n' "\$f" >> gvcf_inputs.args; done
    gatk CombineGVCFs -R ${ref_dir}/Homo_sapiens_assembly38.fasta --arguments_file gvcf_inputs.args \\
      -L '${params.interval}' -O cohort.combined.g.vcf.gz
    gatk GenotypeGVCFs -R ${ref_dir}/Homo_sapiens_assembly38.fasta -V cohort.combined.g.vcf.gz \\
      -L '${params.interval}' -O cohort.joint.vcf.gz
    """
}

process JOINT_GENOMICSDB {
    // Keep import and genotyping in one task, on node-local scratch.
    scratch true
    publishDir "${params.outdir}/vcf", pattern: 'cohort.joint.vcf.gz*', mode: 'copy'
    input:
    path gvcfs
    path ref_dir

    output:
    path 'cohort.joint.vcf.gz'
    path 'cohort.joint.vcf.gz.tbi'

    script:
    """
    set -euo pipefail
    : > cohort.sample_map
    for f in *.g.vcf.gz; do
      sample="\${f%.g.vcf.gz}"
      printf '%s\\t%s\\n' "\$sample" "\$f" >> cohort.sample_map
    done
    gatk GenomicsDBImport --sample-name-map cohort.sample_map \\
      --genomicsdb-workspace-path genomicsdb_brca1 -L '${params.interval}' \\
      --reader-threads ${task.cpus}
    gatk GenotypeGVCFs -R ${ref_dir}/Homo_sapiens_assembly38.fasta -V gendb://genomicsdb_brca1 \\
      -L '${params.interval}' -O cohort.joint.vcf.gz
    """
}
