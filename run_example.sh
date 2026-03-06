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

# After running, check the results:
echo "Pipeline completed! Check the following files:"
echo "  results/pruned_bac120.tree          - Pruned tree in Newick format"
echo "  results/pruned_tip_metadata.tsv     - Per-tip metadata (long format)"
echo "  results/retained_genomes.tsv        - Genome to phylum mapping"
echo "  results/target_phyla.txt            - List of target phyla"
echo "  results/pruning_report.txt          - Summary statistics"
