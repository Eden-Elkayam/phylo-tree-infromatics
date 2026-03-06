# Example R script for working with phylum-level representative tree
# This script demonstrates how to load and visualize the representative tree and metadata

library(ape)
library(ggtree)
library(tidyverse)

# =============================================================================
# 1. Load the phylum-level representative tree
# =============================================================================

tree <- read.tree("results/phylum_representative.tree")

# Check tree structure
cat("Tree has", length(tree$tip.label), "tips\n")
cat("Tree info:\n")
print(tree)

# =============================================================================
# 2. Load metadata - link representative genomes to phylum-level data
# =============================================================================

# Load full metadata (may have multiple rows per genome if CSV had multiple rows per phylum)
metadata_full <- read_tsv("results/representative_tip_metadata.tsv")

cat("Full metadata rows:", nrow(metadata_full), "\n")
cat("Unique genomes:", n_distinct(metadata_full$genome_id), "\n")

# For simpler analysis, get one row per genome (first occurrence)
metadata <- metadata_full %>%
  distinct(genome_id, .keep_all = TRUE)

cat("Simple metadata rows:", nrow(metadata), "\n")

# Check what columns we have
print(colnames(metadata))

# =============================================================================
# 3. Basic tree visualization
# =============================================================================

# Simple tree plot
(p_simple <- ggtree(tree) + theme_tree2())

# Tree with small points at tips
(p_dots <- ggtree(tree) +
  geom_tippoint(size = 2) +
  theme_tree2())

# =============================================================================
# 4. Plot with metadata - size by n_genomes, color by electron_donor
# =============================================================================

# Join metadata to tree for visualization
(p_colored <- ggtree(tree) %<+% metadata +
  geom_tippoint(aes(color = electron_donor, size = n_genomes), alpha = 0.6) +
  scale_size_continuous(trans = "log10", name = "n_genomes") +
  theme_tree2() +
  theme(legend.position = "right"))

# =============================================================================
# 5. Add phylum labels at tips
# =============================================================================

# This may be crowded if tree is large - adjust size as needed
(p_labeled <- p_colored +
   geom_tiplab(aes(label = phylum), size = 1, offset = 0.1))

# =============================================================================
# 6. Faceted visualization by electron donor
# =============================================================================

# Show how representatives are distributed across electron donors
(p_faceted <- ggtree(tree) %<+% metadata +
   geom_tippoint(aes(color = phylum, size = n_genomes), alpha = 0.7) +
   scale_size_continuous(trans = "log10") +
   facet_wrap(~electron_donor) +
   theme_tree2() +
   theme(legend.position = "bottom"))

# =============================================================================
# 7. Summary statistics by phylum
# =============================================================================

# Number of genomes per phylum (from metadata)
summary_by_phylum <- metadata %>%
  group_by(phylum) %>%
  summarize(
    representative_genome = first(genome_id),
    total_genomes = first(n_genomes),
    domain = first(domain),
    .groups = "drop"
  ) %>%
  arrange(desc(total_genomes))

cat("\nTop 20 phyla by genome count:\n")
print(summary_by_phylum %>% head(20))

# Summary by electron donor
summary_by_donor <- metadata %>%
  group_by(electron_donor) %>%
  summarize(
    phyla = n_distinct(phylum),
    total_genomes = sum(n_genomes),
    .groups = "drop"
  ) %>%
  arrange(desc(phyla))

cat("\nPhyla per electron donor:\n")
print(summary_by_donor)

# =============================================================================
# 8. Check representative genomes selected
# =============================================================================

# Load the mapping of which genome was selected for each phylum
representatives <- read_tsv("results/representative_genomes.tsv")

cat("\nRepresentative genomes per phylum:\n")
print(representatives %>% head(10))

# =============================================================================
# 9. View all genomes per phylum (before selection)
# =============================================================================

# This shows the pool from which representatives were selected
phylum_genomes <- read_tsv("results/phylum_to_genomes.tsv")

# Example: show all genomes for Pseudomonadota
pseudomonadota_genomes <- phylum_genomes %>%
  filter(phylum == "Pseudomonadota") %>%
  pull(all_genome_ids) %>%
  str_split(";") %>%
  unlist()

cat("\nPseudomonadota has", length(pseudomonadota_genomes), "genomes total\n")
cat("Representative selected:", 
    representatives %>% filter(phylum == "Pseudomonadota") %>% pull(representative_genome_id),
    "\n")
cat("First 5 alternatives:\n")
print(head(setdiff(pseudomonadota_genomes, 
                   representatives %>% filter(phylum == "Pseudomonadota") %>% pull(representative_genome_id)), 5))

# =============================================================================
# 10. Save a high-quality plot
# =============================================================================

# Save the main visualization
pdf("phylum_representative_tree.pdf", width = 12, height = 8)
print(p_colored + ggtitle("GTDB Phylum-Level Representative Tree"))
dev.off()

cat("\nPlot saved to: phylum_representative_tree.pdf\n")

# =============================================================================
# 11. Example: Filter to specific electron donors
# =============================================================================

# If you only want to analyze H2-oxidizing phyla
h2_metadata <- metadata %>%
  filter(electron_donor == "H2 oxidation") %>%
  distinct(genome_id, .keep_all = TRUE)

cat("\nH2-oxidizing phyla:", nrow(h2_metadata), "\n")

# Filter tree to only these genomes
h2_genomes <- h2_metadata %>% pull(genome_id)
h2_tree <- keep.tip(tree, h2_genomes)

cat("H2-oxidizing tree tips:", length(h2_tree$tip.label), "\n")

# Plot H2-oxidizing subgroup
(p_h2 <- ggtree(h2_tree) %<+% h2_metadata +
   geom_tippoint(aes(size = n_genomes), color = "darkblue", alpha = 0.6) +
   geom_tiplab(aes(label = phylum), size = 2, offset = 0.1) +
   scale_size_continuous(trans = "log10") +
   theme_tree2() +
   ggtitle("H2-oxidizing phyla only"))

# =============================================================================

cat("\n✓ Example analysis complete!\n")
cat("Try modifying the visualizations to explore your data:\n")
cat("  - Change 'electron_donor' to other metadata columns\n")
cat("  - Filter by domain (Bacteria vs Archaea)\n")
cat("  - Look at specific phyla or clades\n")
