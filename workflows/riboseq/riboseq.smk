from snakemake.utils import min_version
import pandas as pd
min_version("6.6")


configfile: "riboseq.config.yaml"
samples = pd.read_table(config['sample_list'], delim_whitespace=True).set_index('sample', drop=False)

def get_fastq(wildcards):
    return samples.loc[wildcards.sample].loc['file']

rule all:
    input:
        expand("mapped_reads/{sample}_Aligned.sortedByCoord.out.bam", sample=samples['sample']),
        expand("genome_cov/{sample}_fw.bw", sample=samples['sample']),
        expand("genome_cov/{sample}_rc.bw", sample=samples['sample']),
        expand("quant/{sample}_gene.txt", sample=samples['sample']),
        expand("quant/{sample}_tx.txt", sample=samples['sample']),


rule pseudotx:
    input:
        gtf=config['gtf'],
        genome=config['genome']
    output:
        "psite/pseudotx.fa"
    threads: 1
    run:
        shell("python scripts/GTFtools.py convert2bed -g {genome} -t cds -e 30 | "
              "bedtools getfasta -nameOnly -fi {genome} -bed stdin >{output}")

rule cutadapt_trim:
    input:
        get_fastq
    output:
        "clean_reads/{sample}_clean.fq.gz"
    log:
        "logs/trim_{sample}.log"
    threads: 4
    params:
        adaptor=lambda w: samples.loc[w.sample].loc['adaptor']
    shell:
        "cutadapt -j {threads} -m 18 --trim-n -a {params.adaptor} -o {output} {input} 2>&1 >{log}"

rule bowtie2_filter:
    input:
        "clean_reads/{sample}_clean.fq.gz"
    output:
        "clean_reads/{sample}_filter.fq.gz"
    log:
        "logs/filter_{sample}.log"
    threads: 8
    params:
        bt2_idx=config['bt2_idx']
    shell:
        "bowtie2 -p {threads} --local --un-gz {output} -x {params.bt2_idx} -U {input} >/dev/null 2>{log}"

rule star_map:
    input:
        "clean_reads/{sample}_filter.fq.gz"
    output:
        "mapped_reads/{sample}_Aligned.sortedByCoord.out.bam"
    log:
        "logs/map_{sample}.log"
    threads: 32
    params:
        star_idx=config['star_idx'],
        output_prefix="mapped_reads/{sample}_"
    resources: mem_mb=16384
    shell:
        "STAR --runMode alignReads --genomeDir {params.star_idx} --readFilesIn {input} --readFilesCommand zcat "
        "--outFileNamePrefix {params.output_prefix} --runThreadN {threads} --outSAMtype BAM SortedByCoordinate "
        "--outSAMattributes All --limitBAMsortRAM 16000000000 2>&1 >{log}"

rule index_bam:
    input:
        "mapped_reads/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "mapped_reads/{sample}_Aligned.sortedByCoord.out.bam.bai"
    threads: 1
    shell:
        "samtools index {input}"


rule quant_gene:
    input:
        bam="mapped_reads/{sample}_Aligned.sortedByCoord.out.bam",
        bai="mapped_reads/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        "quant/{sample}_gene.txt"
    log:
        "logs/quant_{sample}_gene.log"
    threads: 1
    params:
        gtf=config['gtf'],
        stranded=config['stranded'],
        rlen_lower=config['rlen_lower'],
        rlen_upper=config['rlen_upper']
    shell:
        "python scripts/bam_filter.py -l {params.rlen_lower} -u {params.rlen_upper} -q 10 {input.bam} | "
        "featureCounts -t CDS -Q 10 -s 1 -a {params.gtf} -o {output} 2>{log}"

rule quant_tx:
    input:
        bam="mapped_reads/{sample}_Aligned.sortedByCoord.out.bam",
        bai="mapped_reads/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        "quant/{sample}_tx.txt"
    log:
        "logs/quant_{sample}_tx.log"
    threads: 1
    params:
        gtf=config['gtf'],
        stranded=config['stranded'],
        rlen_lower=config['rlen_lower'],
        rlen_upper=config['rlen_upper']
    shell:
        "python scripts/bam_filter.py -l {params.rlen_lower} -u {params.rlen_upper} -q 10 {input.bam} | "
        "featureCounts -t CDS -g transcript_id -O -Q 10 -s 1 -a {params.gtf} -o {output} 2>{log}"

rule genome_cov:
    input:
        bam="mapped_reads/{sample}_Aligned.sortedByCoord.out.bam",
        bai="mapped_reads/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        "genome_cov/{sample}_fw.bw",
        "genome_cov/{sample}_rc.bw"
    log:
        "logs/genome_cov_{sample}.log"
    threads: 4
    params:
        outprefix="genome_cov/{sample}",
        poffset=config['poffset'],
        rlen_lower=config['rlen_lower'],
        rlen_upper=config['rlen_upper']
    shell:
        "python scripts/bam_coverage.py -p {params.poffset} -l {params.rlen_lower} -u {params.rlen_upper} {input.bam} {params.outprefix}"

