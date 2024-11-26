import glob
from pathlib import Path
import os
import logging

# Define base directories
SCRATCH_DIR = "/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta"
GENOME_DIR = SCRATCH_DIR
RESULTS_DIR = os.path.join(SCRATCH_DIR, "results")

# Create necessary directories
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(os.path.join(RESULTS_DIR, "fastqc"), exist_ok=True)
os.makedirs(os.path.join(RESULTS_DIR, "multiqc"), exist_ok=True)
os.makedirs(os.path.join(RESULTS_DIR, "bowtie2"), exist_ok=True)
os.makedirs(os.path.join(RESULTS_DIR, "peaks"), exist_ok=True)
os.makedirs(os.path.join(RESULTS_DIR, "bigwig"), exist_ok=True)
os.makedirs(os.path.join(SCRATCH_DIR, "logs"), exist_ok=True)
os.makedirs(os.path.join(RESULTS_DIR, "peak_analysis"), exist_ok=True)

# Load config file
configfile: "config.yaml"

# Get all samples and their directories
SAMPLE_DIRS = {}
for fastq in glob.glob("90-1102945428/00_fastq/*_R1_001.fastq.gz"):
    sample = Path(fastq).stem.replace("_R1_001.fastq", "")
    if sample.startswith(('BG', 'BM')):  # Only process BG and BM samples
        sample_dir = "90-1102945428"  # Hardcode the sample directory
        SAMPLE_DIRS[sample] = sample_dir

SAMPLES = list(SAMPLE_DIRS.keys())

# Add after SAMPLES definition to verify sample discovery
print("Found samples:", SAMPLES)

# Debug information
print("\nSample paths:")
for sample, directory in SAMPLE_DIRS.items():
    r1 = f"{directory}/00_fastq/{sample}_R1_001.fastq.gz"
    print(f"{sample}: {r1} (exists: {os.path.exists(r1)})")

# Add this after SAMPLES definition
print("\nTarget outputs:")
print("\nFastQC outputs:")
print(expand(os.path.join(RESULTS_DIR, "fastqc/{sample}_R{read}_001_fastqc.html"), sample=SAMPLES, read=[1,2]))
print("\nBAM outputs:")
print(expand(os.path.join(RESULTS_DIR, "bowtie2/{sample}.sorted.bam"), sample=SAMPLES))
print("\nPeak outputs:")
print(expand(os.path.join(RESULTS_DIR, "peaks/{sample}_peaks.narrowPeak"), sample=SAMPLES))

rule verify_index:
    output:
        touch(temp(os.path.join(config["log_dir"], "index_verified.done")))
    params:
        index = config["genome"]["index"],
        log_file = os.path.join(config["log_dir"], "index_verification.log")
    run:
        with open(params.log_file, "w") as log:
            log.write(f"Verifying index files at: {params.index}\n")
            required_files = [f"{params.index}.{i}.bt2" for i in range(1,5)] + \
                           [f"{params.index}.rev.{i}.bt2" for i in range(1,3)]
            
            for f in required_files:
                exists = os.path.exists(f)
                log.write(f"Checking {f}: {'EXISTS' if exists else 'MISSING'}\n")
                if exists:
                    log.write(f"File size: {os.path.getsize(f)} bytes\n")
            
            missing = [f for f in required_files if not os.path.exists(f)]
            if missing:
                error_msg = f"Missing index files: {missing}"
                log.write(f"ERROR: {error_msg}\n")
                raise ValueError(error_msg)
            else:
                log.write("All index files found successfully!\n")

rule bowtie2_align:
    input:
        r1 = lambda wildcards: f"{SAMPLE_DIRS[wildcards.sample]}/00_fastq/{wildcards.sample}_R1_001.fastq.gz",
        r2 = lambda wildcards: f"{SAMPLE_DIRS[wildcards.sample]}/00_fastq/{wildcards.sample}_R2_001.fastq.gz",
        index_verified = os.path.join(config["log_dir"], "index_verified.done")
    output:
        bam = os.path.join(RESULTS_DIR, "bowtie2/{sample}.sorted.bam"),
        bai = os.path.join(RESULTS_DIR, "bowtie2/{sample}.sorted.bam.bai")
    threads: 8
    params:
        index = config["genome"]["index"]
    log:
        os.path.join(config["log_dir"], "bowtie2/{sample}.log")
    shell:
        """
        echo "Starting alignment for {wildcards.sample}" > {log}
        echo "Using index: {params.index}" >> {log}
        bowtie2 -p {threads} -x {params.index} \
            -1 {input.r1} -2 {input.r2} 2>> {log} | \
        samtools sort -@ {threads} -o {output.bam} - && \
        samtools index {output.bam}
        echo "Finished alignment for {wildcards.sample}" >> {log}
        """

rule all:
    input:
        # QC outputs
        expand(os.path.join(RESULTS_DIR, "fastqc/{sample}_R{read}_001_fastqc.html"), 
               sample=SAMPLES, read=[1,2]),
        os.path.join(RESULTS_DIR, "multiqc/multiqc_report.html"),
        # Analysis outputs
        expand(os.path.join(RESULTS_DIR, "bowtie2/{sample}.sorted.bam"), 
               sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "bowtie2/{sample}.sorted.bam.bai"), 
               sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "peaks/{sample}_peaks.narrowPeak"), 
               sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "bigwig/{sample}.bw"), sample=SAMPLES),
        # os.path.join(RESULTS_DIR, "peak_analysis/peak_size_comparison.txt")

rule fastqc:
    input:
        r1 = lambda wildcards: f"{SAMPLE_DIRS[wildcards.sample]}/00_fastq/{wildcards.sample}_R1_001.fastq.gz",
        r2 = lambda wildcards: f"{SAMPLE_DIRS[wildcards.sample]}/00_fastq/{wildcards.sample}_R2_001.fastq.gz"
    output:
        html_r1 = os.path.join(RESULTS_DIR, "fastqc/{sample}_R1_001_fastqc.html"),
        html_r2 = os.path.join(RESULTS_DIR, "fastqc/{sample}_R2_001_fastqc.html"),
        zip_r1 = os.path.join(RESULTS_DIR, "fastqc/{sample}_R1_001_fastqc.zip"),
        zip_r2 = os.path.join(RESULTS_DIR, "fastqc/{sample}_R2_001_fastqc.zip")
    threads: 2
    shell:
        """
        fastqc -t {threads} -o {RESULTS_DIR}/fastqc {input.r1} {input.r2}
        """

rule multiqc:
    input:
        expand(os.path.join(RESULTS_DIR, "fastqc/{sample}_R{read}_001_fastqc.zip"), 
               sample=SAMPLES, read=[1,2])
    output:
        os.path.join(RESULTS_DIR, "multiqc/multiqc_report.html")
    shell:
        """
        multiqc {RESULTS_DIR}/fastqc -o {RESULTS_DIR}/multiqc
        """

rule call_peaks:
    input:
        treatment = os.path.join(RESULTS_DIR, "bowtie2/{sample}.sorted.bam")
    output:
        os.path.join(RESULTS_DIR, "peaks/{sample}_peaks.narrowPeak")
    shell:
        """
        macs2 callpeak -t {input.treatment} \
            -f BAMPE -g mm -n {wildcards.sample} \
            --outdir {RESULTS_DIR}/peaks \
            --nomodel --extsize 200
        """

rule create_bigwig:
    input:
        bam = os.path.join(RESULTS_DIR, "bowtie2/{sample}.sorted.bam"),
        bai = os.path.join(RESULTS_DIR, "bowtie2/{sample}.sorted.bam.bai")
    output:
        os.path.join(RESULTS_DIR, "bigwig/{sample}.bw")
    threads: 4
    shell:
        """
        bamCoverage -b {input.bam} -o {output} \
            --binSize 10 --normalizeUsing RPKM \
            --numberOfProcessors {threads}
        """

# rule merge_peak_groups:
#     input:
#         bg_peaks = expand(os.path.join(RESULTS_DIR, "peaks/{sample}_peaks.narrowPeak"),
#                          sample=[s for s in SAMPLES if s.startswith('BG')]),
#         bm_peaks = expand(os.path.join(RESULTS_DIR, "peaks/{sample}_peaks.narrowPeak"),
#                          sample=[s for s in SAMPLES if s.startswith('BM')])
#     output:
#         bg_merged = os.path.join(RESULTS_DIR, "peak_analysis/BG_merged.bed"),
#         bm_merged = os.path.join(RESULTS_DIR, "peak_analysis/BM_merged.bed")
#     shell:
#         """
#         cat {input.bg_peaks} | sort -k1,1 -k2,2n | bedtools merge > {output.bg_merged}
#         cat {input.bm_peaks} | sort -k1,1 -k2,2n | bedtools merge > {output.bm_merged}
#         """

# rule find_common_cpg_peaks:
#     input:
#         bg_peaks = os.path.join(RESULTS_DIR, "peak_analysis/BG_merged.bed"),
#         bm_peaks = os.path.join(RESULTS_DIR, "peak_analysis/BM_merged.bed"),
#         cpg_islands = "cpg_islands.bed"
#     output:
#         common_peaks = os.path.join(RESULTS_DIR, "peak_analysis/common_cpg_peaks.bed")
#     shell:
#         """
#         # Find peaks present in both BG and BM
#         bedtools intersect -a {input.bg_peaks} -b {input.bm_peaks} -wa | \
#         # Then intersect with CpG islands
#         bedtools intersect -a - -b {input.cpg_islands} -wa | \
#         sort -k1,1 -k2,2n > {output.common_peaks}
#         """

# rule count_reads_in_peaks:
#     input:
#         peaks = os.path.join(RESULTS_DIR, "peak_analysis/common_cpg_peaks.bed"),
#         bam = os.path.join(RESULTS_DIR, "bowtie2/{sample}.sorted.bam")
#     output:
#         counts = os.path.join(RESULTS_DIR, "peak_analysis/{sample}_peak_counts.txt")
#     shell:
#         """
#         bedtools coverage -a {input.peaks} -b {input.bam} -counts > {output.counts}
#         """

# rule compare_peak_sizes:
#     input:
#         peak_counts = expand(os.path.join(RESULTS_DIR, "peak_analysis/{sample}_peak_counts.txt"),
#                            sample=SAMPLES)
#     output:
#         comparison = os.path.join(RESULTS_DIR, "peak_analysis/peak_size_comparison.txt")
#     params:
#         min_bg_samples = 3,  # We expect 3 BG samples
#         min_bm_samples = 1   # We expect 1 BM sample
#     script:
#         "scripts/compare_peaks.py"

# Debug information
for sample, directory in SAMPLE_DIRS.items():
    r1 = f"{directory}/00_fastq/{sample}_R1_001.fastq.gz"
    r2 = f"{directory}/00_fastq/{sample}_R2_001.fastq.gz"
    print(f"Checking {sample}:")
    print(f"R1 exists: {os.path.exists(r1)}")
    print(f"R2 exists: {os.path.exists(r2)}")

print("\nDEBUG - All found FASTQ files:")
for fastq in glob.glob("*/00_fastq/*_R1_001.fastq.gz"):
    print(f"Found FASTQ: {fastq}")

print("\nDEBUG - Sample directories:")
for sample, directory in SAMPLE_DIRS.items():
    print(f"Sample: {sample}, Directory: {directory}") 