import pandas as pd


configfile: "example_config.yaml"
samples = pd.read_table(config['samples']).set_index('sample', drop=False)

def get_fastq(wildcards):
    return samples.loc[wildcards.sample].loc['file']

rule all:
    input:
        "plots/quals.svg"


rule bwa_map:
    input:
        config["genome"],
        get_fastq
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule bcftools_call:
    input:
        fa=config["genome"],
        bam=expand("sorted_reads/{sample}.bam", sample=samples['sample'].to_list()),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=samples['sample'].to_list())
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"

