
# Read alignment and variant calling tutorial using SBATCH on the UMIC HPC

This tutorial demonstrates alignment of paired-end Illumina FASTQ reads to **GRCh38**, preprocessing of BAM files, and germline variant calling across a **BRCA1 region on chromosome 17**. The paths for this training are tailored to UMIC HPC users under the UMIC-Roche Bursary Initiative.

> **Scope:** This is a teaching workflow for an interval (BRCA1 region), not a fully validated clinical pipeline. The tutorial assumes single-library paired-end reads per sample.

## 1. Workflow overview

```text
Paired-end FASTQ (R1, R2)
         |
         v
BWA-MEM alignment + SAMtools sorting
         |
         v
Sorted BAM + BAM index
         |
         v
GATK MarkDuplicates
         |
         v
Deduplicated BAM + duplicate metrics
         |
         v
GATK BaseRecalibrator -> ApplyBQSR
         |
         v
Recalibrated BAM + BAM index
         |
         v
GATK HaplotypeCaller (-ERC GVCF)
         |
         v
Per-sample GVCFs
         |
         +--> GenotypeGVCFs (one sample)
         |
         +--> Joint genotyping (choose ONE):
                  A. GenomicsDBImport -> GenotypeGVCFs
                  B. CombineGVCFs -> GenotypeGVCFs
```

## 2. Requirements

- Access to a Linux HPC with **SLURM** (`sbatch`, `squeue`, `sacct`).
- Environment modules for **BWA**, **SAMtools**, **GATK**, and optionally **BCFtools**.
- Paired-end FASTQ files named `<SAMPLE>_R1.fastq.gz` and `<SAMPLE>_R2.fastq.gz`.
- GRCh38 reference FASTA with its `.fai`, sequence dictionary (`.dict`), and BWA index files.
- GRCh38 dbSNP and known-indel VCFs with compatible indexes.

Key paths

```bash
REF=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.fasta
DBSNP=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.dbsnp138.vcf.gz
KNOWN_INDELS=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.known_indels.vcf.gz
READS_DIR=/mnt/roche/resources/datasets/human-genomics/short-read
BRCA1=chr17:43041295-43173327
```


## 3. Prepare a sample list


Create a `samples.txt` file with **one sample ID per line**. For example:

```text
HG02562
HG02573
HG02585
HG02678
HG02703
```

Additional samples can be found at `/mnt/roche/resources/datasets/human-genomics/short-read`

For each sample, the script expects files such as:

```text
/mnt/roche/resources/datasets/human-genomics/short-read/HG02562_R1.fastq.gz
/mnt/roche/resources/datasets/human-genomics/short-read/HG02562_R2.fastq.gz
```

(Optional) Check that the input files exist:

```bash
while read -r SAMPLE; do
    for READ in R1 R2; do
        test -r "${READS_DIR}/${SAMPLE}_${READ}.fastq.gz" || echo "Missing: ${SAMPLE}_${READ}.fastq.gz"
    done
done < samples.txt
```

Run this check after setting `READS_DIR` in your shell.

## 4. Create the SBATCH script

Save the following as `varcall.sh`. 

```bash
#!/bin/bash -l
#SBATCH --job-name=varcall
#SBATCH --output=varcall_%j.out
#SBATCH --error=varcall_%j.err
#SBATCH --account=roche
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G

set -euo pipefail

module load bwa samtools bcftools gatk

mkdir -p bam qc gvcf vcf

REF=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.fasta
DBSNP=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.dbsnp138.vcf.gz
KNOWN_INDELS=/mnt/roche/resources/hg38_ref_genome/Homo_sapiens_assembly38.known_indels.vcf.gz
READS_DIR=/mnt/roche/resources/datasets/human-genomics/short-read
BRCA1=chr17:43041295-43173327
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Select one joint-genotyping method: combine or genomicsdb
JOINT_METHOD=combine

while read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue
    R1="${READS_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${READS_DIR}/${SAMPLE}_R2.fastq.gz"
    RG="@RG\\tID:${SAMPLE}.L1\\tSM:${SAMPLE}\\tLB:${SAMPLE}_lib1\\tPL:ILLUMINA\\tPU:FLOWCELL.L1"

    echo "Processing ${SAMPLE}"

    # 1. Align and coordinate-sort
    bwa mem -t "$THREADS" -R "$RG" "$REF" "$R1" "$R2" \
        | samtools sort -@ "$THREADS" -o "bam/${SAMPLE}.sorted.bam" -
    samtools index "bam/${SAMPLE}.sorted.bam"

    # 2. Mark duplicates and index
    gatk MarkDuplicates \
        -I "bam/${SAMPLE}.sorted.bam" \
        -O "bam/${SAMPLE}.dedup.bam" \
        -M "qc/${SAMPLE}.dup_metrics.txt"
    samtools index "bam/${SAMPLE}.dedup.bam"

    # 3. Build BQSR table for the target interval
    gatk BaseRecalibrator \
        -I "bam/${SAMPLE}.dedup.bam" -R "$REF" \
        --known-sites "$DBSNP" --known-sites "$KNOWN_INDELS" \
        -L "$BRCA1" -O "bam/${SAMPLE}.recal.table"

    # 4. Apply BQSR over the target interval
    gatk ApplyBQSR \
        -I "bam/${SAMPLE}.dedup.bam" -R "$REF" \
        --bqsr-recal-file "bam/${SAMPLE}.recal.table" \
        -L "$BRCA1" -O "bam/${SAMPLE}.recal.bam"
    samtools index "bam/${SAMPLE}.recal.bam"

    # 5. Call variants in GVCF mode
    gatk HaplotypeCaller \
        -R "$REF" -I "bam/${SAMPLE}.recal.bam" \
        -L "$BRCA1" -ERC GVCF \
        -O "gvcf/${SAMPLE}.g.vcf.gz"
done < /mnt/roche/resources/datasets/human-genomics/short-read/samples.txt

# 6. Generate a tab-delimited sample-to-GVCF map
awk '{print $1 "\tgvcf/" $1 ".g.vcf.gz"}' samples.txt > cohort.sample_map

# 7. Choose ONE joint-genotyping approach
if [[ "$JOINT_METHOD" == "genomicsdb" ]]; then
    gatk GenomicsDBImport \
        --genomicsdb-workspace-path genomicsdb_brca1 \
        --sample-name-map cohort.sample_map \
        -L "$BRCA1" --reader-threads "$THREADS"

    gatk GenotypeGVCFs \
        -R "$REF" -V gendb://genomicsdb_brca1 \
        -L "$BRCA1" -O vcf/cohort.joint_genomicsdb.vcf.gz

elif [[ "$JOINT_METHOD" == "combine" ]]; then
    awk '{print "-V", $2}' cohort.sample_map > gvcf_inputs.args

    gatk CombineGVCFs \
        -R "$REF" --arguments_file gvcf_inputs.args \
        -L "$BRCA1" -O gvcf/cohort.g.vcf.gz

    gatk GenotypeGVCFs \
        -R "$REF" -V gvcf/cohort.g.vcf.gz \
        -L "$BRCA1" -O vcf/cohort.joint_combine.vcf.gz
else
    echo "Unknown JOINT_METHOD: $JOINT_METHOD" >&2
    exit 1
fi

echo "Variant calling complete"
```

> **Resource note:** To optimise runtime, tweak `--cpus-per-task`, `--tasks-per-node` and `THREADS`. 

> **BQSR note:** Because `ApplyBQSR` uses `-L "$BRCA1"`, the recalibrated BAM contains only the selected interval, rather than the entire genome. Remove that `-L` restriction if you need a genome-wide recalibrated BAM. The known-sites resources should be suitable for the reference build.

## 5. Submit and monitor the job

Submit from the directory containing `samples.txt`:

```bash
sbatch varcall.sh
```

Monitor your jobs:

```bash
squeue -u "$USER"
```

Inspect completed-job accounting:

```bash
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode
```

Inspect logs:

```bash
tail -n 50 varcall_JOB_ID.err
tail -n 50 varcall_JOB_ID.out
```

Replace `JOB_ID` with the actual SLURM job ID.

## 6. What each step does

| Step | Tool | Purpose | Main output |
|---|---|---|---|
| Alignment | `bwa mem` | Map paired-end reads to GRCh38 | Alignment stream |
| Sorting and indexing | `samtools sort`, `samtools index` | Coordinate-sort reads and index BAM | `*.sorted.bam` + `.bai` |
| Duplicate marking | `gatk MarkDuplicates` | Flag duplicate reads and report metrics | `*.dedup.bam`, `*.dup_metrics.txt` |
| Base recalibration | `gatk BaseRecalibrator` | Estimate systematic base-quality errors using known sites | `*.recal.table` |
| Apply recalibration | `gatk ApplyBQSR` | Adjust base quality scores | `*.recal.bam` + `.bai` |
| Per-sample variant calling | `gatk HaplotypeCaller -ERC GVCF` | Generate GVCF records for later genotyping | `*.g.vcf.gz` |
| Joint genotyping | GATK | Produce cohort-level genotypes | `cohort.joint_*.vcf.gz` |

### Read groups matter

The `-R` argument includes an `ID` and sample name (`SM`) for every sample. The `\\t` in the Bash variable assignment represents the **literal backslash-t sequence** that BWA interprets as a tab separator. Do not replace it with actual tab characters in the argument passed to BWA.

## 7. Optional: genotype a single sample

The supplied script also demonstrates genotyping one GVCF (example: `HG02562`):

```bash
gatk GenotypeGVCFs \
    -R "$REF" \
    -V gvcf/HG02562.g.vcf.gz \
    -L "$BRCA1" \
    -O vcf/HG02562.raw.vcf.gz
```

Run this **after** the per-sample GVCF has been created, in an environment where `REF` and `BRCA1` are defined. The complete example script above focuses on cohort joint genotyping; this single-sample command is optional.

## 8. Choosing a joint-genotyping approach

**Option A — GenomicsDBImport:** Builds a GenomicsDB workspace and jointly genotypes the imported GVCFs. This is generally preferred for larger GATK cohorts, but its database operations can be sensitive to filesystem behaviour. On some Lustre configurations, a node-local workspace may be more reliable; ensure both import and genotyping can access the same workspace on the same node.

**Option B — CombineGVCFs:** Combines GVCFs into a multi-sample GVCF before `GenotypeGVCFs`. It avoids GenomicsDB and is simpler for smaller cohorts, but generally scales less well as sample count increases.

The script defaults to `JOINT_METHOD=combine`. Change it to `genomicsdb` to select Option A. **Do not delete the GenomicsDB workspace until genotyping has succeeded and the outputs are verified.**

## 9. Inspect the results

Expected directory structure:

```text
.
├── samples.txt
├── varcall.sbatch
├── cohort.sample_map
├── bam/
│   ├── HG02562.sorted.bam
│   ├── HG02562.dedup.bam
│   ├── HG02562.recal.table
│   └── HG02562.recal.bam
├── qc/
│   └── HG02562.dup_metrics.txt
├── gvcf/
│   ├── HG02562.g.vcf.gz
│   └── cohort.g.vcf.gz        # CombineGVCFs option only
└── vcf/
    └── cohort.joint_combine.vcf.gz
```

For the GenomicsDB option, the final VCF is named `vcf/cohort.joint_genomicsdb.vcf.gz` instead.

Check the BAM and VCF:

```bash
samtools quickcheck -v bam/HG02562.recal.bam
samtools flagstat bam/HG02562.recal.bam
bcftools view -h vcf/cohort.joint_combine.vcf.gz | tail
bcftools query -l vcf/cohort.joint_combine.vcf.gz
```

For the GenomicsDB option, substitute the appropriate final VCF filename.

## 10. Common troubleshooting

| Symptom | What to check |
|---|---|
| `bwa_set_rg` errors | Read-group `-R` formatting; use escaped `\\t`, not literal tabs |
| `samtools sort: failed to read header` | Check upstream `bwa mem` error output |
| Reference/index errors | FASTA `.fai`, `.dict`, and BWA indexes; matching contig names |
| GATK cannot read BAM | Ensure BAM exists and is indexed |
| GenomicsDB cannot open array | Workspace path, filesystem compatibility, and whether import completed |
| Job runs out of memory | Inspect `sacct` and adjust `#SBATCH --mem` |
| Job takes too long | The original script loops over samples **sequentially** in one SLURM job; use job arrays or Nextflow for sample-level parallelism |

## 11. What comes next?

**Part 2** will implement the same analysis with **Nextflow DSL2**, separating alignment/preprocessing and variant calling into reusable modules, with independent tasks, SLURM scheduling, and `-resume` support.

---

