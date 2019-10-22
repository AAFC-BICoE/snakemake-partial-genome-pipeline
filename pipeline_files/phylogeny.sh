# Code to produce a default Phyluce phylogeny using Mafft and RAxML
# Adapted from https://phyluce.readthedocs.io/en/latest/tutorial-one.html, Copyright 2012-2015, Brant C. Faircloth

# Conda enviroment phyenv.yml and iqtree.yml must be installed before this code will work. Provided primarily
# as an example for generated a simple phylogeny with target enrichment data
# Phyluce runs in a Python 2.7 environment, such as the one provided with phyenv.yml

CORES=80
TAXA=106
conda activate phyenv
mkdir log
# Run Mafft
phyluce_align_seqcap_align \
    --fasta all-taxa-incomplete.fasta \
    --output mafft-nexus-internal-trimmed \
    --taxa $TAXA \
    --aligner mafft \
    --cores $CORES \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path log

# Run gblocks trimming on the alignments
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments mafft-nexus-internal-trimmed \
    --output mafft-nexus-internal-trimmed-gblocks \
    --cores $CORES \
    --log log

phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-internal-trimmed-gblocks \
    --cores $CORES \
    --log-path log

phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments mafft-nexus-internal-trimmed-gblocks \
    --output mafft-nexus-internal-trimmed-gblocks-clean \
    --cores $CORES \
    --log-path log

# Create a 75p Completeness Matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean \
    --taxa $TAXA \
    --percent 0.75 \
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --cores $CORES \
    --log-path log

phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean \
    --taxa $TAXA \
    --percent 0.50 \
    --output mafft-nexus-internal-trimmed-gblocks-clean-50p \
    --cores $CORES \
    --log-path log

# Prepare files for RAxML
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml \
    --charsets \
    --log-path log

phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean-50p \
    --output mafft-nexus-internal-trimmed-gblocks-clean-50p-raxml \
    --charsets \
    --log-path log

conda activate iqtree
cd mafft-nexus-internal-trimmed-gblocks-clean-50p-raxml
iqtree -s mafft-nexus-internal-trimmed-gblocks-clean-50p.phylip -m MFP -bb 1000 -alrt 1000 -nt AUTO

cd ../mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml
iqtree -s mafft-nexus-internal-trimmed-gblocks-clean-75p.phylip -m MFP -bb 1000 -alrt 1000 -nt AUTO
