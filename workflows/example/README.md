Example snakemake workflow
--------------------------
The example workflow is build based on the tutorial of [snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

To see what it does:
```
$ snakemake -np -s example.smk
Building DAG of jobs...
Job stats:
job               count    min threads    max threads
--------------  -------  -------------  -------------
all                   1              1              1
bcftools_call         1              1              1
bwa_map               2              1              1
plot_quals            1              1              1
samtools_index        2              1              1
samtools_sort         2              1              1
total                 9              1              1


[Mon Aug 23 11:53:41 2021]
rule bwa_map:
    input: genome.fa, B.fastq
    output: mapped_reads/B.bam
    jobid: 6
    wildcards: sample=B
    resources: tmpdir=/tmp

bwa mem -t 1 genome.fa B.fastq | samtools view -Sb - > mapped_reads/B.bam

[Mon Aug 23 11:53:41 2021]
rule bwa_map:
    input: genome.fa, A.fastq
    output: mapped_reads/A.bam
    jobid: 4
    wildcards: sample=A
    resources: tmpdir=/tmp

bwa mem -t 1 genome.fa A.fastq | samtools view -Sb - > mapped_reads/A.bam

[Mon Aug 23 11:53:41 2021]
rule samtools_sort:
    input: mapped_reads/B.bam
    output: sorted_reads/B.bam
    jobid: 5
    wildcards: sample=B
    resources: tmpdir=/tmp

samtools sort -T sorted_reads/B -O bam mapped_reads/B.bam > sorted_reads/B.bam

[Mon Aug 23 11:53:41 2021]
rule samtools_sort:
    input: mapped_reads/A.bam
    output: sorted_reads/A.bam
    jobid: 3
    wildcards: sample=A
    resources: tmpdir=/tmp

samtools sort -T sorted_reads/A -O bam mapped_reads/A.bam > sorted_reads/A.bam

[Mon Aug 23 11:53:41 2021]
rule samtools_index:
    input: sorted_reads/A.bam
    output: sorted_reads/A.bam.bai
    jobid: 7
    wildcards: sample=A
    resources: tmpdir=/tmp

samtools index sorted_reads/A.bam

[Mon Aug 23 11:53:41 2021]
rule samtools_index:
    input: sorted_reads/B.bam
    output: sorted_reads/B.bam.bai
    jobid: 8
    wildcards: sample=B
    resources: tmpdir=/tmp

samtools index sorted_reads/B.bam

[Mon Aug 23 11:53:41 2021]
rule bcftools_call:
    input: genome.fa, sorted_reads/A.bam, sorted_reads/B.bam, sorted_reads/A.bam.bai, sorted_reads/B.bam.bai
    output: calls/all.vcf
    jobid: 2
    resources: tmpdir=/tmp

samtools mpileup -g -f genome.fa sorted_reads/A.bam sorted_reads/B.bam | bcftools call -mv - > calls/all.vcf

[Mon Aug 23 11:53:41 2021]
rule plot_quals:
    input: calls/all.vcf
    output: plots/quals.svg
    jobid: 1
    resources: tmpdir=/tmp

[Mon Aug 23 11:53:41 2021]
localrule all:
    input: plots/quals.svg
    jobid: 0
    resources: tmpdir=/tmp

Job stats:
job               count    min threads    max threads
--------------  -------  -------------  -------------
all                   1              1              1
bcftools_call         1              1              1
bwa_map               2              1              1
plot_quals            1              1              1
samtools_index        2              1              1
samtools_sort         2              1              1
total                 9              1              1

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

