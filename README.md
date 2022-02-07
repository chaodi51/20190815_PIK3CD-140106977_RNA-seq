# Snakemake workflow: rna-seq_pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/ribo-seq_pipeline.svg?branch=master)](https://travis-ci.org/snakemake-workflows/ribo-seq_pipeline)

This is a RNA-seq pipeline structured by snakemake template. It is for the PI3K project stand on the standard RNA-seq pipeline.
Insert your code into the respective folders, i.e. `scripts`, `rules`, and `envs`. Define the entry point of the workflow in the `Snakefile` and the main configuration in the `config.yaml` file.

## Authors

* chaodi (@dic)

## Usage
Running on respublica by:
snakemake --use-conda -c "qsub -l h_vmem={params.mem} -l mem_free={params.mem} -pe smp {threads} -V -cwd -e qsub/{params.jobName}.e -o qsub/{params.jobName}.o" -j -p

