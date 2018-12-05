# Snakemake file to process partial genome sequences generated from UCE target enrichment experiments
# Author: Jackson Eyres jackson.eyres@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

import glob
import os
from shutil import copyfile

# Configuration Settings

# Location of fastq folder, default "fastq". Phyluce requires files to not mix "-" and "_", so fastq files renamed
for f in glob.glob('fastq/*.fastq.gz'):
    basename = os.path.basename(f)
    new_basename = ""

    new_basename = basename.replace("-","_")
    if new_basename[0].isdigit(): # Phyluce doesn't work with samples starting with digits
        new_basename = "Sample_" + new_basename
    new_path = f.replace(basename, new_basename)
    os.rename(f, new_path)

SAMPLES = set([os.path.basename(f).replace("_L001_R1_001.fastq.gz","").replace("_L001_R2_001.fastq.gz","") for f in glob.glob('fastq/*.fastq.gz')])

# Location of adaptor.fa for trimming
adaptors = "pipeline_files/adapters.fa"


rule all:
    input:
        fastqc_dir = directory("fastqc"),

        r1_trimmed = expand("trimmed/{sample}_trimmed_L001_R1_001.fastq.gz", sample=SAMPLES),
        r2_trimmed = expand("trimmed/{sample}_trimmed_L001_R2_001.fastq.gz", sample=SAMPLES),

        fastqc_trimmed_dir = directory("fastqc_trimmed"),

        multiqc_report = "multiqc/multiqc_report.html",
        multiqc_report_trimmed = "multiqc/multiqc_report_trimmed.html",

        spades_assemblies = expand("spades_assemblies/{sample}/contigs.fasta", sample=SAMPLES),
        rna_assemblies = expand("rnaspades_assemblies/{sample}/transcripts.fasta", sample=SAMPLES),

        phyluce_spades_assemblies = expand("phyluce-spades/assemblies/{sample}_S.fasta", sample=SAMPLES),
        phyluce_rnaspades_assemblies = expand("phyluce-rnaspades/assemblies/{sample}_R.fasta", sample=SAMPLES),

        spades_taxon = "phyluce-spades/taxon.conf",
        rnaspades_taxon = "phyluce-rnaspades/taxon.conf",

        spades_db = "phyluce-spades/uce-search-results/probe.matches.sqlite",
        spades_log = "phyluce-spades/phyluce_assembly_match_contigs_to_probes.log",
        spades_taxon_sets = "phyluce-spades/taxon-sets/all/all-taxa-incomplete.conf",
        spades_all_taxa = "phyluce-spades/taxon-sets/all/all-taxa-incomplete.fasta",
        spades_exploded_fastas = directory("phyluce-spades/taxon-sets/all/exploded-fastas"),
        spades_exploded_locus = directory("phyluce-spades/taxon-sets/all/exploded-locus"),

        rnaspades_db = "phyluce-rnaspades/uce-search-results/probe.matches.sqlite",
        rnaspades_log = "phyluce-rnaspades/phyluce_assembly_match_contigs_to_probes.log",
        rnaspades_taxon_sets = "phyluce-rnaspades/taxon-sets/all/all-taxa-incomplete.conf",
        rnaspades_all_taxa = "phyluce-rnaspades/taxon-sets/all/all-taxa-incomplete.fasta",
        rnaspades_exploded_fastas = directory("phyluce-rnaspades/taxon-sets/all/exploded-fastas"),
        rnaspades_exploded_locus = directory("phyluce-rnaspades/taxon-sets/all/exploded-locus"),

        merged_uces = directory("merged_uces"),

        spades_assembly_metrics = "spades_assembly_metrics.tsv",
        rnaspades_assembly_metrics = "rnaspades_assembly_metrics.tsv",

        phyluce_spades_uce_metrics = "phyluce_spades_uce_metrics.tsv",
        phyluce_rnaspades_uce_metrics = "phyluce_rnaspades_uce_metrics.tsv",


        spades_uce_summary = "spades_uce_summary.csv",
        rnaspades_uce_summary = "rnaspades_uce_summary.csv"

rule fastqc:
    # Quality Control check on raw data before adaptor trimming
    input:
        directory("fastq")
    output:
         fastqc_dir = directory("fastqc")
    log: "logs/fastqc.log"
    conda: "pipeline_files/pg_assembly.yml"
    shell:
        "mkdir fastqc; fastqc -o fastqc fastq/*.fastq.gz 2>{log} 2>&1"


rule bbduk:
    # Sequencing Adaptor trimming
    input:
        r1 = 'fastq/{sample}_L001_R1_001.fastq.gz',
        r2 = 'fastq/{sample}_L001_R2_001.fastq.gz'
    output:
        out1 = "trimmed/{sample}_trimmed_L001_R1_001.fastq.gz",
        out2 = "trimmed/{sample}_trimmed_L001_R2_001.fastq.gz",
    log: "logs/bbduk.{sample}.log"
    conda: "pipeline_files/pg_assembly.yml"
    shell: "bbduk.sh in1={input.r1} out1={output.out1} in2={input.r2} out2={output.out2} ref={adaptors} qtrim=rl trimq=20 ktrim=r k=23 mink=11 hdist=1 tpe tbo 2>{log} 2>&1; touch {output.out1} {output.out2}"


rule fastqc_trimmed:
    # Quality Control check after adaptor trimming
    input: r1=expand("trimmed/{sample}_trimmed_L001_R1_001.fastq.gz", sample=SAMPLES),
        r2 = expand("trimmed/{sample}_trimmed_L001_R2_001.fastq.gz", sample=SAMPLES)
    output:
        fastqc_trimmed_dir = directory("fastqc_trimmed")
    log: "logs/fastqc_trimmed.log"
    conda: "pipeline_files/pg_assembly.yml"
    shell:
        "mkdir fastqc_trimmed; fastqc -o fastqc_trimmed trimmed/*.fastq.gz 2>{log} 2>&1"


rule multiqc:
    # Consolidates all QC files into single report pre/post trimming
    input:
        dir1 = directory("fastqc"),
        dir2 = directory("fastqc_trimmed")
    output:
        r1_report = "multiqc/multiqc_report.html",
        r2_report = "multiqc/multiqc_report_trimmed.html",
    conda: "pipeline_files/multiqcenv.yml"
    shell:
        "multiqc -n multiqc_report.html -o multiqc fastqc; multiqc -n multiqc_report_trimmed.html -o multiqc fastqc_trimmed;"


rule spades:
    # Assembles fastq files using default settings
    input:
        r1 = "trimmed/{sample}_trimmed_L001_R1_001.fastq.gz",
        r2 = "trimmed/{sample}_trimmed_L001_R2_001.fastq.gz"
    output:
        "spades_assemblies/{sample}/contigs.fasta"
    log: "logs/spades.{sample}.log"
    conda: "pipeline_files/pg_assembly.yml"
    shell:
        "spades.py -1 {input.r1} -2 {input.r2} -o spades_assemblies/{wildcards.sample} 2>{log} 2>&1"


rule rnaspades:
    # Variation of SPAdes that assembles fastq files using default settings
    input:
        r1 = "trimmed/{sample}_trimmed_L001_R1_001.fastq.gz",
        r2 = "trimmed/{sample}_trimmed_L001_R2_001.fastq.gz"
    output:
        "rnaspades_assemblies/{sample}/transcripts.fasta"
    log: "logs/rnaspades.{sample}.log"
    conda: "pipeline_files/pg_assembly.yml"
    shell:
        "rnaspades.py -1 {input.r1} -2 {input.r2} -o rnaspades_assemblies/{wildcards.sample} 2>{log} 2>&1"


rule gather_assemblies:
    # Rename all spades assemblies and copy to a folder for further analysis
    input:
        assembly = "spades_assemblies/{sample}/contigs.fasta"
    output:
        renamed_assembly = "phyluce-spades/assemblies/{sample}_S.fasta"
    run:
        if os.path.exists(input.assembly):
            if os.path.exists("phyluce-spades/assemblies"):
                pass
            else:
                os.path.mkdir("phyluce-spades/assemblies")
            copyfile(input.assembly,output.renamed_assembly)


rule gather_rna_assemblies:
    # Rename all rnaspades assemblies and copy to a folder for further analysis
    input:
        assembly = "rnaspades_assemblies/{sample}/transcripts.fasta"
    output:
        renamed_assembly = "phyluce-rnaspades/assemblies/{sample}_R.fasta"
    run:
        if os.path.exists(input.assembly):
            if os.path.exists("phyluce-rnaspades/assemblies"):
                pass
            else:
                os.path.mkdir("phyluce-rnaspades/assemblies")
            copyfile(input.assembly,output.renamed_assembly)


rule generate_taxons_conf:
    # List of assembly names required for Phyluce processing
    output: w1="phyluce-spades/taxon.conf", w2="phyluce-rnaspades/taxon.conf"
    run:
        with open (output.w1, "w") as f:
            f.write("[all]\n")
            for item in SAMPLES:
                f.write(item + "_S\n")
        with open (output.w2, "w") as f:
            f.write("[all]\n")
            for item in SAMPLES:
                f.write(item + "_R\n")


# The following scripts are derived from https://phyluce.readthedocs.io/en/latest/tutorial-one.html
# This tutorial outlines the key steps of following Phyluce to derive UCEs from probes and assemblies
rule phyluce_spades:
    input: "phyluce-rnaspades/taxon.conf"
    output: db="phyluce-spades/uce-search-results/probe.matches.sqlite", log="phyluce-spades/phyluce_assembly_match_contigs_to_probes.log"
    conda: "pipeline_files/phyenv.yml"
    #Must remove the auto generated output directory before running script
    shell: "rm -r phyluce-spades/uce-search-results; cd phyluce-spades; phyluce_assembly_match_contigs_to_probes --keep-duplicates KEEP_DUPLICATES --contigs assemblies --output uce-search-results --probes ../probes/*.fasta"

rule phyluce_assembly_get_match_counts_spades:
    input: conf="phyluce-spades/taxon.conf", db="phyluce-spades/uce-search-results/probe.matches.sqlite"
    output: "phyluce-spades/taxon-sets/all/all-taxa-incomplete.conf"
    conda: "pipeline_files/phyenv.yml"
    shell: "cd phyluce-spades; phyluce_assembly_get_match_counts --locus-db uce-search-results/probe.matches.sqlite --taxon-list-config taxon.conf --taxon-group 'all' --incomplete-matrix --output taxon-sets/all/all-taxa-incomplete.conf"

rule phyluce_assembly_get_fastas_from_match_counts_spades:
    input: "phyluce-spades/uce-search-results/probe.matches.sqlite"
    output: "phyluce-spades/taxon-sets/all/all-taxa-incomplete.fasta"
    conda: "pipeline_files/phyenv.yml"
    shell: "cd phyluce-spades/taxon-sets/all; mkdir log; phyluce_assembly_get_fastas_from_match_counts --contigs ../../assemblies --locus-db ../../uce-search-results/probe.matches.sqlite --match-count-output all-taxa-incomplete.conf --output all-taxa-incomplete.fasta --incomplete-matrix all-taxa-incomplete.incomplete --log-path log"

rule phyluce_assembly_explode_get_fastas_file_spades:
    input: alignments="phyluce-spades/taxon-sets/all/all-taxa-incomplete.fasta"
    output: exploded_fasta=directory("phyluce-spades/taxon-sets/all/exploded-fastas"), exploded_locus=directory("phyluce-spades/taxon-sets/all/exploded-locus")
    conda: "pipeline_files/phyenv.yml"
    # Command requires --input and not --alignments
    shell: "cd phyluce-spades/taxon-sets/all; phyluce_assembly_explode_get_fastas_file --input all-taxa-incomplete.fasta --output exploded-fastas --by-taxon; phyluce_assembly_explode_get_fastas_file --input all-taxa-incomplete.fasta --output exploded-locus"

rule phyluce_rnaspades:
    input: "phyluce-rnaspades/taxon.conf"
    output: db="phyluce-rnaspades/uce-search-results/probe.matches.sqlite", log="phyluce-rnaspades/phyluce_assembly_match_contigs_to_probes.log"
    conda: "pipeline_files/phyenv.yml"
    #Must remove the auto generated output directory before running script
    shell: "rm -r phyluce-rnaspades/uce-search-results; cd phyluce-rnaspades; phyluce_assembly_match_contigs_to_probes --keep-duplicates KEEP_DUPLICATES --contigs assemblies --output uce-search-results --probes ../probes/*.fasta"

rule phyluce_assembly_get_match_counts_rnaspades:
    input: conf="phyluce-rnaspades/taxon.conf", db="phyluce-rnaspades/uce-search-results/probe.matches.sqlite"
    output: "phyluce-rnaspades/taxon-sets/all/all-taxa-incomplete.conf"
    conda: "pipeline_files/phyenv.yml"
    shell: "cd phyluce-rnaspades; phyluce_assembly_get_match_counts --locus-db uce-search-results/probe.matches.sqlite --taxon-list-config taxon.conf --taxon-group 'all' --incomplete-matrix --output taxon-sets/all/all-taxa-incomplete.conf"

rule phyluce_assembly_get_fastas_from_match_counts_rnaspades:
    input: "phyluce-rnaspades/uce-search-results/probe.matches.sqlite"
    output: "phyluce-rnaspades/taxon-sets/all/all-taxa-incomplete.fasta"
    conda: "pipeline_files/phyenv.yml"
    shell: "cd phyluce-rnaspades/taxon-sets/all; mkdir log; phyluce_assembly_get_fastas_from_match_counts --contigs ../../assemblies --locus-db ../../uce-search-results/probe.matches.sqlite --match-count-output all-taxa-incomplete.conf --output all-taxa-incomplete.fasta --incomplete-matrix all-taxa-incomplete.incomplete --log-path log"

rule phyluce_assembly_explode_get_fastas_file_rnaspades:
    input: alignments="phyluce-rnaspades/taxon-sets/all/all-taxa-incomplete.fasta"
    output: exploded_fasta=directory("phyluce-rnaspades/taxon-sets/all/exploded-fastas"), exploded_locus=directory("phyluce-rnaspades/taxon-sets/all/exploded-locus")
    conda: "pipeline_files/phyenv.yml"
    # Command requires --input and not --alignments
    shell: "cd phyluce-rnaspades/taxon-sets/all; phyluce_assembly_explode_get_fastas_file --input all-taxa-incomplete.fasta --output exploded-fastas --by-taxon; phyluce_assembly_explode_get_fastas_file --input all-taxa-incomplete.fasta --output exploded-locus"


rule combine_uces:
    # Combines SPAdes and rnaSPAdes derived UCEs into singular files for further analysis
    input: exploded_spades_fastas=directory("phyluce-spades/taxon-sets/all/exploded-fastas"), exploded_rnaspades_fastas=directory("phyluce-rnaspades/taxon-sets/all/exploded-fastas")
    output: directory("merged_uces")
    conda: "pipeline_files/pg_assembly.yml"
    shell: "python pipeline_files/merge_uces.py -o {output} -s {input.exploded_spades_fastas} -r {input.exploded_rnaspades_fastas}"

rule spades_quality_metrics:
    # BBMap's Stats.sh assembly metrics for spades assemblies
    input: "phyluce-spades/taxon.conf"
    output: "spades_assembly_metrics.tsv"
    conda: "pipeline_files/pg_assembly.yml"
    shell: "statswrapper.sh phyluce-spades/assemblies/*.fasta > {output}"

rule rnaspades_quality_metrics:
    # BBMap's Stats.sh assembly metrics for rnaspades assemblies
    input: "phyluce-rnaspades/taxon.conf"
    output: "rnaspades_assembly_metrics.tsv"
    conda: "pipeline_files/pg_assembly.yml"
    shell: "statswrapper.sh phyluce-rnaspades/assemblies/*.fasta > {output}"

rule fastq_quality_metrics:
    # BBMap's Stats.sh assembly metrics for fastq files
    input: directory("fastq")
    output: "fastq_metrics.tsv"
    conda: "pipeline_files/pg_assembly.yml"
    shell: "statswrapper.sh {input}/*.fastq.gz > {output}"

rule phyluce_rnaspades_quality_metrics:
# BBMap's Stats.sh assembly metrics for rnaspades assemblies
    input: directory("phyluce-rnaspades/taxon-sets/all/exploded-fastas")
    output: "phyluce_rnaspades_uce_metrics.tsv"
    conda: "pipeline_files/pg_assembly.yml"
    shell: "statswrapper.sh {input}/*.fasta > {output}"

rule phyluce_spades_quality_metrics:
# BBMap's Stats.sh assembly metrics for rnaspades assemblies
    input: directory("phyluce-spades/taxon-sets/all/exploded-fastas")
    output: "phyluce_spades_uce_metrics.tsv"
    conda: "pipeline_files/pg_assembly.yml"
    shell: "statswrapper.sh {input}/*.fasta > {output}"

rule summarize_uces:
    input: r1="phyluce-rnaspades/phyluce_assembly_match_contigs_to_probes.log", s1="phyluce-spades/phyluce_assembly_match_contigs_to_probes.log"
    output: r2="rnaspades_uce_summary.csv", s2="spades_uce_summary.csv"
    conda: "pipeline_files/pg_assembly.yml"
    shell: "python pipeline_files/evaluate.py -r1 {input.r1} -r2 {output.r2} -s1 {input.s1} -s2 {output.s2} "