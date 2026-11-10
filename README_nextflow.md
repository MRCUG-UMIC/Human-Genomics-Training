# Read alignment and variant calling with Nextflow DSL2

In Part 1, we used a Bash script submitted with `sbatch` to align paired-end Illumina reads to GRCh38, preprocess BAMs, call variants in a BRCA1 interval, and jointly genotype samples. This part introduces **Nextflow DSL2**, as an alternative approach. It covers the same analysis but tracks individual steps, manages dependencies, and can schedule samples independently through SLURM.

> **Scope:** A teaching example for `chr17:43041295-43173327` (GRCh38).

## 1. What changes from SBATCH?

| Part 1: Bash + SBATCH | Part 2: Nextflow DSL2 |
|---|---|
| A shell `while read SAMPLE` loop | A channel emits one sample and its FASTQs at a time |
| Commands run sequentially in a single batch script | Separate processes declare inputs and outputs |
| You manually manage intermediate files | Nextflow stages files in task-specific `work/` directories |
| Restarting may repeat completed steps | `-resume` can reuse successful tasks |
| SLURM resources apply to the batch job | CPU, memory, and time can be configured per process |

## 2. Workflow overview

```text
samples.txt + paired FASTQs
           |
           v
       ALIGNMENT (BWA-MEM)
           |
           v
       SORT_BAM (SAMtools)
           |
           v
       MARK_DUPLICATES (GATK)
           |
           v
       BASE_RECALIBRATOR (GATK)
           |
           v
       APPLY_BQSR (GATK)
           |
           v
       VARIANT_CALLING (HaplotypeCaller, per sample)
           |
           v
       Collect cohort GVCFs
           |
           +---- JOINT_COMBINE (CombineGVCFs + GenotypeGVCFs)
           |
           +---- JOINT_GENOMICSDB (GenomicsDBImport + GenotypeGVCFs)
```

The **two module files** group related processes; they are not limited to one process each.

## 3. Project layout

Clone the repository to your project directory:

```bash
cd /mnt/roche/${USER}_roche_project
git clone HG-nextflow-pipeline
cd HG-nextflow-pipeline
```

The supplied archive contains:

```text
HG-nextflow-pipeline/
├── main.nf
├── nextflow.config
└── modules/
    ├── alignment.nf
    └── variant_calling.nf
```

- `main.nf`: defines the workflow and connects processes.
- `nextflow.config`: sets default input paths, resources, and execution profiles.
- `modules/alignment.nf`: contains `ALIGNMENT`, `SORT_BAM`, `MARK_DUPLICATES`, `BASE_RECALIBRATOR`, and `APPLY_BQSR`.
- `modules/variant_calling.nf`: contains `VARIANT_CALLING`, `JOINT_COMBINE`, and `JOINT_GENOMICSDB`.

## 4. Requirements and input data

You need Nextflow, Java compatible with your installed Nextflow release, and environment modules for BWA, SAMtools, and GATK. The `slurm` profile also requires access to SLURM commands (`sbatch`, `squeue`).

Check the installations:

```bash
module load nextflow
nextflow -version
module avail bwa samtools gatk
sinfo
```

The example input sample list is one identifier per line:

```text
HG02562
HG02573
```

For each sample, Nextflow expects paired FASTQs in the directory given by `--reads`:

```text
HG02562_R1.fastq.gz
HG02562_R2.fastq.gz
HG02573_R1.fastq.gz
HG02573_R2.fastq.gz
```

The supplied configuration uses:

```text
Reference: /mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.fasta
Reads (Directory):     /mnt/roche/resources/datasets/human-genomics/short-read
Interval:  chr17:43041295-43173327
Output:    results/
```

The reference directory must contain the FASTA, its BWA index, FASTA index (`.fai`), sequence dictionary (`.dict`), and the two known-sites VCFs with their indexes. 

## 5. How `main.nf` connects the processes

The supplied workflow first reads `samples.txt`, constructs a pair of FASTQ paths for each sample, and emits tuples resembling:

```groovy
['HG02562', [file('HG02562_R1.fastq.gz'), file('HG02562_R2.fastq.gz')]]
```

It then passes the reference directory as a reusable value channel:

```groovy
def refFile = file(params.ref, checkIfExists: true)
def refDir = refFile.parent
ref_ch = Channel.value(refDir)
```

The process chain in the supplied `main.nf` is:

```groovy
aligned = ALIGNMENT(samples_ch, ref_ch)
sorted = SORT_BAM(aligned)
dedup = MARK_DUPLICATES(sorted)
recal = BASE_RECALIBRATOR(dedup.bam, ref_ch)
final_bams = APPLY_BQSR(recal, ref_ch)
called = VARIANT_CALLING(final_bams, ref_ch)
```

A **tuple** keeps a sample identifier together with its BAM and index. For example, `APPLY_BQSR` emits:

```groovy
tuple val(sample), path("${sample}.recal.bam"), path("${sample}.recal.bam.bai")
```

The variant-calling process accepts that same three-element tuple. After individual GVCFs are produced, the workflow collects them and selects either joint-calling method according to `--joint_method`.

## 6. What happens in each module?

### Module A — `modules/alignment.nf`

1. **`ALIGNMENT`:** BWA-MEM maps paired FASTQs to GRCh38 and writes an unsorted BAM through `samtools view`.
2. **`SORT_BAM`:** sorts and indexes the BAM with SAMtools.
3. **`MARK_DUPLICATES`:** identifies duplicate reads with GATK and publishes duplication metrics to `results/qc/`.
4. **`BASE_RECALIBRATOR`:** builds a recalibration table using dbSNP and known-indel sites in the BRCA1 interval.
5. **`APPLY_BQSR`:** applies the recalibration, indexes the resulting BAM, and publishes it to `results/bam/`.

### Module B — `modules/variant_calling.nf`

1. **`VARIANT_CALLING`:** runs GATK HaplotypeCaller with `-ERC GVCF` for each sample and publishes per-sample GVCFs to `results/gvcf/`.
2. **`JOINT_COMBINE`:** combines GVCFs with `CombineGVCFs` and runs `GenotypeGVCFs`.
3. **`JOINT_GENOMICSDB`:** imports GVCFs into a GenomicsDB workspace and runs `GenotypeGVCFs` in the same task. The supplied process has `scratch true`.

Only one joint-genotyping route is chosen per run.

## 7. Configure and run on SLURM

The supplied `nextflow.config` defines two relevant profiles:

```groovy
profiles {
    modules {
        process.beforeScript = 'module load bwa samtools gatk'
    }
    slurm {
        process.executor = 'slurm'
        process.queue = 'short'
    }
}
```

Replace `short` with a valid partition if rerunning this on another HPC. Check available partitions with:

```bash
sinfo -o '%P %a %l %D'
```

Run with both profiles:

```bash
mkdir /tmp/${USER} 
nextflow run main.nf -profile modules,slurm --joint_method combine -w /tmp/${USER}/nextflow-work
```

To select GenomicsDB instead:

```bash
nextflow run main.nf -profile modules,slurm --joint_method genomicsdb -w /tmp/${USER}/nextflow-work
```

The defaults in `nextflow.config` can be overriden on the command line. 


> **Important:** Node-local `/tmp` may be deleted or unavailable from another node, so it is not ideal for a long-lived `-resume` cache. Prefer persistent storage that supports locking when available.

## 8. Results and monitoring

Expected published outputs include:

```text
results/
├── bam/
│   ├── HG02562.recal.bam
│   └── HG02562.recal.bam.bai
├── qc/
│   └── HG02562.dup_metrics.txt
├── gvcf/
│   ├── HG02562.g.vcf.gz
│   └── HG02562.g.vcf.gz.tbi
└── vcf/
    ├── cohort.joint.vcf.gz
    └── cohort.joint.vcf.gz.tbi
```

The `work/` directory contains intermediate task outputs. Check SLURM jobs and inspect task logs:

```bash
squeue -u "$USER"
nextflow log
find work -name .command.err | head
```

To rerun after fixing a failed process, use `-resume` **from the same launch directory**, retaining the Nextflow cache and work directory.


## 9. Important checks and corrections

The following are **recommendations based on inspection of the supplied archive**, not changes already present in that archive.

**A. Task memory and runtime.** The supplied config uses two CPUs and two-hour limits for most processes. These are example settings, not guaranteed to be sufficient for real WGS data. Review SLURM limits and adjust as needed.

**B. Joint calling and sample counts.** Confirm that the number of per-sample GVCFs matches the expected sample list before interpreting the cohort VCF. Joint calling should not proceed with a silently incomplete cohort.

## 10. Key takeaways

- **DSL2 modules** organise reusable processes; a single module file can contain several processes.
- **Channels and tuples** connect outputs from one process to the next.
- **Profiles** allow the same workflow to use environment modules and SLURM (`-profile modules,slurm`).
- **`-resume`** reuses completed tasks when the execution cache and work directory are retained.
- **`--joint_method`** selects between CombineGVCFs and GenomicsDBImport in the supplied workflow.

**Suggested exercise:** Run two samples with `--joint_method combine`, inspect the published BAMs/GVCFs/VCF, then rerun with `-resume` and observe which tasks Nextflow reuses.

---

