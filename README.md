# Phylo Tree Curation

Small utilities for curating phylogenetic trees and preparing them for visualization.

This repo currently includes tools for:
- pruning large GTDB trees
- selecting representative genomes (e.g. one genome per phylum)
- attaching taxonomy and metadata
- exporting annotations for visualization tools like **iTOL** and **ggtree**

The goal is to go from large genome phylogenies to clean, interpretable trees that are easy to visualize and explore

Files required for this repo to run:
- bac120.tree
- bac120_taxonomy.tsv (both are found on GTDB website (gtdb.ecogenomic.org)
- csv file for genome count by electron source, according to C1-genomics output (A. Flamholz) 
