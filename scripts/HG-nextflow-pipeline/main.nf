nextflow.enable.dsl = 2

include { 
    ALIGNMENT; 
    SORT_BAM;
    MARK_DUPLICATES;
    BASE_RECALIBRATOR;
    APPLY_BQSR 
    } from './modules/alignment'
include { VARIANT_CALLING; JOINT_COMBINE; JOINT_GENOMICSDB } from './modules/variant_calling'

workflow {
    if (!(params.joint_method in ['combine', 'genomicsdb'])) {
        error "--joint_method must be 'combine' or 'genomicsdb'"
    }
    def samplesFile = file(params.samples)
    if (!samplesFile.exists()) error "Missing samples file: ${params.samples}"

    def refFile = file(params.ref, checkIfExists: true)
    def refDir = refFile.parent

    samples_ch = Channel.fromPath(params.samples).splitText()
        .map { it.trim() }
        .filter { it && !it.startsWith('#') }
        .map { sample ->
            def r1 = file("${params.reads}/${sample}_R1.fastq.gz")
            def r2 = file("${params.reads}/${sample}_R2.fastq.gz")
            if (!r1.exists() || !r2.exists()) error "FASTQs missing for ${sample}"
            tuple(sample, [r1, r2])
        }

    ref_ch = Channel.value(refDir)   

    aligned = ALIGNMENT(samples_ch, ref_ch)
    sorted = SORT_BAM(aligned)
    dedup = MARK_DUPLICATES(sorted)
    recal = BASE_RECALIBRATOR(dedup.bam, ref_ch)
    final_bams = APPLY_BQSR(recal, ref_ch)
    called = VARIANT_CALLING(final_bams, ref_ch)
    gvcfs = called.map { sample, gvcf, index -> [gvcf, index] }.collect().map { it.flatten() }

    if (params.joint_method == 'combine') {
        JOINT_COMBINE(gvcfs, ref_ch)
    } else {
        JOINT_GENOMICSDB(gvcfs, ref_ch)
    }
}
