#!/bin/bash
# Example: Generate iTOL annotation files for a phylum-level representative tree

# Basic usage - phylum-level tree with domain and phylum annotations
python generate_iTOL_annotation_files.py \
  --tree results/phylum_representative.tree \
  --metadata results/representative_tip_metadata.tsv \
  --taxonomy bac120_taxonomy.tsv \
  --represented-level phylum \
  --label-rank phylum \
  --strip-ranks domain,phylum \
  --outdir itol_annotations

# With debug taxonomy export
# python generate_iTOL_annotation_files.py \
#   --tree results/phylum_representative.tree \
#   --metadata results/representative_tip_metadata.tsv \
#   --taxonomy bac120_taxonomy.tsv \
#   --represented-level phylum \
#   --label-rank phylum \
#   --strip-ranks domain,phylum \
#   --debug-taxonomy \
#   --outdir itol_annotations

# Example for a genome-level tree (all ranks available)
# python generate_iTOL_annotation_files.py \
#   --tree genome_level.tree \
#   --metadata genome_metadata.tsv \
#   --taxonomy bac120_taxonomy.tsv \
#   --represented-level genome \
#   --label-rank species \
#   --strip-ranks domain,phylum,class,order \
#   --outdir itol_annotations_genome

# Example for a class-level tree (domain, phylum, class available)
# python generate_iTOL_annotation_files.py \
#   --tree class_representative.tree \
#   --metadata class_metadata.tsv \
#   --taxonomy bac120_taxonomy.tsv \
#   --represented-level class \
#   --label-rank class \
#   --strip-ranks domain,phylum,class \
#   --outdir itol_annotations_class
