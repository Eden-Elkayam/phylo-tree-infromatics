#!/bin/bash
# Example: Running the GTDB tree pruning pipeline
#
# This script demonstrates how to run the pruning pipeline
# with your data files.

# Basic usage with default column names
python prune_gtdb_tree.py \
    --csv carbon_by_electron_donor_phylum.csv \
    --tree bac120.tree \
    --taxonomy bac120_taxonomy.tsv \
    --outdir results

# Creating iTOL annotation files for the pruned tree
python generate_iTOL_annotation_files.py \
  --tree results/phylum_representative.tree \
  --metadata results/representative_tip_metadata.tsv \
  --taxonomy bac120_taxonomy.tsv \
  --represented-level phylum \
  --label-rank phylum \
  --strip-ranks domain,phylum \
  --outdir results/itol_annotations


