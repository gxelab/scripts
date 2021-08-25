##### snakemake workflow for RNA-Seq (single-end)
To use this workflow, put sample information in the `rnaseq_se.samples.tsv` following the example file.
If you do not need some of the output files, you can comment out the corresponding lines in rule `all`.

##### How to run this worflow
```bash
# dry-run: show what will be done without execution
snakemake -s riboseq_se.smk -np

# run the workflow if the dry-run seems OK
snakemake -s riboseq_se.smk

# use `--cores` for multithreading on a local server
snakemake -s riboseq_se.smk --cores 16

# use `slurm_profile` to submit jobs with `sbatch` on a cluster
snakemake -s riboseq_se.smk --cores 16 --profile slurm_profile
```
