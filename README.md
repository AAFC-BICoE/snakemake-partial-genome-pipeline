# Snakemake Partial Genome Sequencing Pipeline

Pipeline for processing Illumina sequencing data generated by target enrichment via hybrid capture experiments. Heavily follows the Phyluce methodology outlined in [Tutorial I: UCE Phylogenomics](https://phyluce.readthedocs.io/en/latest/tutorial-one.html) 

1) Trims adapters and bases below <20 quality score [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
2) Assembles trimmed reads [SPAdes, rnaSPAdes](http://cab.spbu.ru/software/spades/)
3) Detects and extracts target contigs [Phyluce](https://phyluce.readthedocs.io/en/latest/index.html) 
4) Summary statistics on targets and assemblies [BBTools Stats](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/statistics-guide/)
5) Optional Script to proceed with Phylogeny

### Prerequisites

* [Conda Installation](https://conda.io/docs/user-guide/install/index.html)
* [Snakemake Installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
* Git


```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install -c bioconda -c conda-forge snakemake
```

## Getting Started

Within a working directory:

```
git clone https://github.com/AAFC-BICoE/snakemake-partial-genome-pipeline.git .
```
* Create a folder named "fastq" that contains Illumina based raw reads in fastq 
* Create a folder named "probes" that contains a probe fasta file with fasta headers in [Phyluce UCE format](https://phyluce.readthedocs.io/en/latest/uce-processing.html#match-contigs-to-probes)

Pipeline can be invoked from within working directory with conda enviroment containing snakemake active
```
snakemake --use-conda -k
```

### Example

Incomplete...
* Short reads from a Coleoptera target enrichment via hybrid capture experiment:  [PRJNA379181](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA379181)
    * [UCE raw reads of Platynectes sp. SLE956](https://www.ncbi.nlm.nih.gov/sra/SRX2727692[accn])
* Probe set available from [Ultra Conserved Elements](http://ultraconserved.org/)
  * [Coleoptera-UCE-1.1K-v1.fasta](https://ndownloader.figshare.com/files/6042081)

```
conda install -c bioconda sra-tools 
mkdir fastq
fastq-dump --split-files SRR5437755
```

## Built With

* [Python](https://www.python.org/doc/) - Programming language
* [Conda](https://conda.io/docs/index.html) - Package, dependency and environment management
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) - Workflow management system
* [Phyluce](https://phyluce.readthedocs.io/en/latest/index.html) - Target enrichment data analysis
* [BioPython](https://biopython.org/) - Tools for biological computation


## Copyright
Government of Canada, Agriculture & Agri-Food Canada

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
