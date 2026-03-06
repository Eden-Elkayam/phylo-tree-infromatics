#!/usr/bin/env python3
"""
Create a Phylum-Level Representative Tree from GTDB

This script creates a **phylum-level representative tree** by:
1. Extracting target phyla from a phylum-level metadata CSV
2. Selecting ONE representative genome per phylum from the full GTDB tree
3. Pruning the tree to retain only these representative genomes
4. Creating metadata linking representative genomes to phylum-level information

The resulting tree has ~N_phyla tips (genome IDs used as phylogenetic proxies)
rather than millions of genome tips.

Usage:
    python prune_gtdb_tree.py --csv electron_donor_by_phylum.csv \\
                               --tree bac120.tree \\
                               --taxonomy bac120_taxonomy.tsv \\
                               --outdir results

Required inputs:
    - CSV file with phylum-level metadata (must have 'phylum' column)
    - GTDB tree in Newick format
    - GTDB taxonomy TSV file

Outputs (in outdir):
    - phylum_representative.tree: Tree with one representative genome per phylum
    - representative_tip_metadata.tsv: Phylum-level metadata for each representative
    - target_phyla.txt: List of target phyla from CSV
    - phylum_to_genomes.tsv: All genomes per phylum (before selection)
    - pruning_report.txt: Summary statistics
"""

import argparse
import sys
import os
from pathlib import Path
import pandas as pd
from collections import defaultdict

# Configuration - easily changeable column names
DEFAULT_COLUMNS = {
    'phylum': 'phylum',
    'n_genomes': 'n_genomes',
    'electron_donor': 'electron_donor',
    'domain': 'domain'
}


def normalize_phylum_name(phylum_str):
    """
    Normalize phylum name by removing GTDB prefix (p__) and stripping whitespace.
    
    Args:
        phylum_str: Raw phylum name (e.g., "p__Acidobacteriota" or "Acidobacteriota")
    
    Returns:
        Normalized phylum name (e.g., "Acidobacteriota")
    """
    if pd.isna(phylum_str):
        return None
    
    phylum = str(phylum_str).strip()
    
    # Remove GTDB prefix if present
    if phylum.startswith('p__'):
        phylum = phylum[3:]
    
    return phylum


def extract_phylum_from_taxonomy(taxonomy_str):
    """
    Extract phylum from GTDB taxonomy string.
    
    Args:
        taxonomy_str: GTDB taxonomy string (e.g., "d__Bacteria;p__Pseudomonadota;c__...")
    
    Returns:
        Normalized phylum name or None
    """
    if pd.isna(taxonomy_str):
        return None
    
    # Split by semicolon and find phylum field
    for field in str(taxonomy_str).split(';'):
        if field.startswith('p__'):
            return normalize_phylum_name(field)
    
    return None


def read_target_phyla(csv_path, phylum_col='phylum'):
    """
    Read phylum-level CSV and extract unique target phyla.
    
    Args:
        csv_path: Path to CSV file
        phylum_col: Name of phylum column
    
    Returns:
        tuple: (set of normalized phylum names, full DataFrame)
    """
    print(f"\n{'='*60}")
    print("STEP 1: Reading phylum-level CSV")
    print(f"{'='*60}")
    
    # Read CSV
    df = pd.read_csv(csv_path)
    print(f"✓ Read CSV with {len(df)} rows and {len(df.columns)} columns")
    
    # Check if phylum column exists
    if phylum_col not in df.columns:
        print(f"\n✗ ERROR: Column '{phylum_col}' not found in CSV")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    # Extract and normalize phyla
    phyla_raw = df[phylum_col].dropna().unique()
    phyla_normalized = {normalize_phylum_name(p) for p in phyla_raw}
    phyla_normalized.discard(None)
    
    print(f"✓ Found {len(phyla_normalized)} unique phyla")
    print(f"  Examples: {', '.join(list(phyla_normalized)[:5])}")
    
    return phyla_normalized, df


def read_taxonomy_and_map_genomes(taxonomy_path, target_phyla):
    """
    Read GTDB taxonomy and map target phyla to ALL genome IDs.
    
    Args:
        taxonomy_path: Path to taxonomy TSV file
        target_phyla: Set of target phylum names
    
    Returns:
        dict: Maps phylum -> list of genome IDs for that phylum
    """
    print(f"\n{'='*60}")
    print("STEP 2: Reading taxonomy and mapping genomes to phyla")
    print(f"{'='*60}")
    
    phylum_genomes_map = {}  # phylum -> [genome_ids]
    total_genomes = 0
    matched_genomes = 0
    
    print("Processing taxonomy file (may take a moment for large files)...")
    
    with open(taxonomy_path, 'r') as f:
        for line in f:
            total_genomes += 1
            parts = line.strip().split('\t')
            
            if len(parts) >= 2:
                genome_id = parts[0]
                taxonomy_str = parts[1]
                phylum = extract_phylum_from_taxonomy(taxonomy_str)
                
                if phylum in target_phyla:
                    if phylum not in phylum_genomes_map:
                        phylum_genomes_map[phylum] = []
                    phylum_genomes_map[phylum].append(genome_id)
                    matched_genomes += 1
            
            # Progress indicator
            if total_genomes % 100000 == 0:
                print(f"  Processed {total_genomes:,} genomes, found {matched_genomes:,} matches...")
    
    print(f"\n✓ Processed {total_genomes:,} total genomes")
    print(f"✓ Found {matched_genomes:,} genomes from {len(phylum_genomes_map)} target phyla")
    
    # Check for phyla that weren't found
    unmatched_phyla = target_phyla - set(phylum_genomes_map.keys())
    if unmatched_phyla:
        print(f"\n⚠ WARNING: {len(unmatched_phyla)} phyla from CSV not found in taxonomy:")
        for phylum in sorted(unmatched_phyla)[:10]:
            print(f"  - {phylum}")
        if len(unmatched_phyla) > 10:
            print(f"  ... and {len(unmatched_phyla) - 10} more")
    
    if len(phylum_genomes_map) == 0:
        print("\n✗ ERROR: No genomes matched! Check phylum name formatting.")
        sys.exit(1)
    
    return phylum_genomes_map


def get_tree_tip_ids(tree_path):
    """
    Extract the set of tip IDs (genome IDs) from the Newick tree file.
    
    Args:
        tree_path: Path to Newick format tree
    
    Returns:
        Set of genome IDs that are tips in the tree
    """
    # Try ete3 first
    try:
        from ete3 import Tree
        tree = Tree(tree_path, format=1)
        return {leaf.name for leaf in tree.get_leaves()}
    except ImportError:
        pass
    
    # Try BioPython
    try:
        from Bio import Phylo
        tree = Phylo.read(tree_path, 'newick')
        return {term.name for term in tree.get_terminals()}
    except ImportError:
        pass
    
    # Fallback: simple regex parsing of Newick
    print("⚠ Warning: Tree libraries not available, using simple Newick parsing...")
    tip_ids = set()
    with open(tree_path, 'r') as f:
        content = f.read()
        # Simple regex to find genome IDs (assumes they don't contain special Newick chars)
        import re
        # Match patterns like RS_GCF_XXXXXXX or similar before :, ), or ,
        matches = re.findall(r'([A-Za-z0-9_\.]+)(?:[:\),]|$)', content)
        tip_ids = set(matches)
    return tip_ids


def select_representative_genomes(phylum_genomes_map, tree_tip_ids):
    """
    Select one representative genome per phylum (the first one present in tree).
    
    Args:
        phylum_genomes_map: dict mapping phylum -> [genome_ids]
        tree_tip_ids: Set of genome IDs that are actually in the tree
    
    Returns:
        tuple: (dict of phylum -> representative_genome_id, 
                DataFrame with phylum -> all_genomes mapping)
    """
    print(f"\n{'='*60}")
    print("STEP 2b: Selecting representative genomes per phylum")
    print(f"{'='*60}")
    
    representative_genomes = {}  # phylum -> genome_id
    phyla_found_in_tree = []
    phyla_not_found_in_tree = []
    
    # Create comprehensive mapping before filtering
    all_mappings = []
    
    for phylum in sorted(phylum_genomes_map.keys()):
        genomes = phylum_genomes_map[phylum]
        all_mappings.append({
            'phylum': phylum,
            'total_genomes_for_phylum': len(genomes),
            'all_genome_ids': ';'.join(sorted(genomes))
        })
    
    # Now select representatives from those present in tree
    for phylum in sorted(phylum_genomes_map.keys()):
        genomes_for_phylum = sorted(phylum_genomes_map[phylum])
        
        # Find first genome that's actually in the tree
        representative = None
        for genome_id in genomes_for_phylum:
            if genome_id in tree_tip_ids:
                representative = genome_id
                break
        
        if representative:
            representative_genomes[phylum] = representative
            phyla_found_in_tree.append(phylum)
        else:
            phyla_not_found_in_tree.append(phylum)
    
    print(f"✓ Selected {len(representative_genomes)} representative genomes")
    print(f"  Found in tree: {len(phyla_found_in_tree)} phyla")
    
    if phyla_not_found_in_tree:
        print(f"\n⚠ WARNING: {len(phyla_not_found_in_tree)} phyla have no genomes in tree:")
        for phylum in phyla_not_found_in_tree[:5]:
            print(f"  - {phylum}")
        if len(phyla_not_found_in_tree) > 5:
            print(f"  ... and {len(phyla_not_found_in_tree) - 5} more")
    
    all_mappings_df = pd.DataFrame(all_mappings)
    
    return representative_genomes, all_mappings_df


def prune_tree_to_representatives(tree_path, representative_genome_ids, output_path):
    """
    Prune tree to retain only representative genomes.
    
    Args:
        tree_path: Path to input tree
        representative_genome_ids: Set of representative genome IDs to retain
        output_path: Path to output tree
    
    Returns:
        tuple: (original_tip_count, pruned_tip_count)
    """
    print(f"\n{'='*60}")
    print("STEP 3: Pruning tree to representative genomes")
    print(f"{'='*60}")
    
    print(f"Reading tree from: {tree_path}")
    
    # Try to use ete3 first (most robust)
    try:
        from ete3 import Tree
        
        tree = Tree(tree_path, format=1)
        original_tips = len(tree.get_leaves())
        print(f"✓ Tree loaded with {original_tips:,} tips")
        
        # Prune tree to keep only representative genomes
        tree.prune(representative_genome_ids, preserve_branch_length=True)
        pruned_tips = len(tree.get_leaves())
        
        # Write pruned tree
        tree.write(format=1, outfile=output_path)
        print(f"✓ Pruned tree saved with {pruned_tips:,} tips (one per phylum)")
        
        return original_tips, pruned_tips
        
    except ImportError:
        print("⚠ ete3 not available, trying BioPython...")
        
        # Try BioPython
        try:
            from Bio import Phylo
            
            tree = Phylo.read(tree_path, 'newick')
            original_tips = len(tree.get_terminals())
            print(f"✓ Tree loaded with {original_tips:,} tips")
            
            # Get terminals to remove
            terminals_to_remove = [term for term in tree.get_terminals() 
                                 if term.name not in representative_genome_ids]
            
            # Prune
            for term in terminals_to_remove:
                tree.prune(term)
            
            pruned_tips = len(tree.get_terminals())
            
            # Write
            Phylo.write(tree, str(output_path), 'newick')
            print(f"✓ Pruned tree saved with {pruned_tips:,} tips (one per phylum)")
            
            return original_tips, pruned_tips
            
        except ImportError:
            print("✗ ERROR: Neither ete3 nor BioPython is installed!")
            print("\nPlease install one of these packages:")
            print("  pip install ete3")
            print("  pip install biopython")
            sys.exit(1)


def prune_tree_simple(tree_path, retained_genome_ids, output_path):
    """
    Prune tree using manual Newick parsing (no external dependencies).
    
    This is a simple implementation that works for basic pruning.
    For complex trees, consider using BioPython or ete3.
    
    Args:
        tree_path: Path to input tree
        retained_genome_ids: Set of genome IDs to retain
        output_path: Path to output tree
    
    Returns:
        tuple: (original_tip_count, pruned_tip_count, set of actual tip IDs in pruned tree)
    """
    print(f"\n{'='*60}")
    print("STEP 3: Pruning tree")
    print(f"{'='*60}")
    
    print(f"Reading tree from: {tree_path}")
    
    # Try to use ete3 first (most robust)
    try:
        from ete3 import Tree
        
        tree = Tree(tree_path, format=1)
        original_tips = len(tree.get_leaves())
        print(f"✓ Tree loaded with {original_tips:,} tips")
        
        # Prune tree
        tree.prune(retained_genome_ids, preserve_branch_length=True)
        pruned_tips = len(tree.get_leaves())
        
        # Get actual tip IDs
        actual_tip_ids = {leaf.name for leaf in tree.get_leaves()}
        
        # Write pruned tree
        tree.write(format=1, outfile=output_path)
        print(f"✓ Pruned tree saved with {pruned_tips:,} tips")
        
        return original_tips, pruned_tips, actual_tip_ids
        
    except ImportError:
        print("⚠ ete3 not available, trying BioPython...")
        
        # Try BioPython
        try:
            from Bio import Phylo
            
            tree = Phylo.read(tree_path, 'newick')
            original_tips = len(tree.get_terminals())
            print(f"✓ Tree loaded with {original_tips:,} tips")
            
            # Get terminals to remove
            terminals_to_remove = [term for term in tree.get_terminals() 
                                 if term.name not in retained_genome_ids]
            
            # Prune
            for term in terminals_to_remove:
                tree.prune(term)
            
            pruned_tips = len(tree.get_terminals())
            
            # Get actual tip IDs
            actual_tip_ids = {term.name for term in tree.get_terminals()}
            
            # Write
            Phylo.write(tree, str(output_path), 'newick')
            print(f"✓ Pruned tree saved with {pruned_tips:,} tips")
            
            return original_tips, pruned_tips, actual_tip_ids
            
        except ImportError:
            print("✗ ERROR: Neither ete3 nor BioPython is installed!")
            print("\nPlease install one of these packages:")
            print("  pip install ete3")
            print("  pip install biopython")
            sys.exit(1)


def create_tip_metadata(representative_genomes, phylum_data_df, phylum_col='phylum'):
    """
    Create per-tip metadata for representative genomes by joining with phylum-level data.
    
    Args:
        representative_genomes: dict mapping phylum -> representative_genome_id
        phylum_data_df: DataFrame with phylum-level metadata from CSV
        phylum_col: Name of phylum column
    
    Returns:
        DataFrame with per-tip metadata
    """
    print(f"\n{'='*60}")
    print("STEP 4: Creating representative tip metadata")
    print(f"{'='*60}")
    
    # Create a simple dataframe with representative genomes
    rep_genomes_df = pd.DataFrame([
        {'genome_id': genome_id, 'phylum': phylum}
        for phylum, genome_id in representative_genomes.items()
    ])
    
    # Normalize phylum names in CSV for joining
    phylum_data_normalized = phylum_data_df.copy()
    phylum_data_normalized[phylum_col] = phylum_data_normalized[phylum_col].apply(normalize_phylum_name)
    
    # Join on phylum
    tip_metadata = rep_genomes_df.merge(
        phylum_data_normalized,
        left_on='phylum',
        right_on=phylum_col,
        how='left'
    )
    
    # Remove duplicate phylum column if it exists
    if f'{phylum_col}_x' in tip_metadata.columns:
        tip_metadata = tip_metadata.drop(columns=[f'{phylum_col}_y'])
        tip_metadata = tip_metadata.rename(columns={f'{phylum_col}_x': phylum_col})
    
    print(f"✓ Created metadata table with {len(tip_metadata)} rows")
    print(f"  Phyla represented: {len(tip_metadata['phylum'].drop_duplicates())}")
    print(f"  Columns: {', '.join(tip_metadata.columns[:8])}")
    if len(tip_metadata.columns) > 8:
        print(f"  ... and {len(tip_metadata.columns) - 8} more")
    
    return tip_metadata


def build_unmatched_phyla_table(phylum_data_df, target_phyla, phylum_genomes_map,
                                phylum_col='phylum', domain_col='domain', n_genomes_col='n_genomes'):
    """
    Build a table of phyla that are in the CSV but not found in taxonomy mapping.

    Args:
        phylum_data_df: Original phylum-level DataFrame from CSV
        target_phyla: Set of normalized target phyla from CSV
        phylum_genomes_map: Dict mapping matched phylum -> genome IDs
        phylum_col: Name of phylum column in CSV
        domain_col: Name of domain column in CSV

    Returns:
        DataFrame with columns: phylum, domain, total_n_genomes
    """
    unmatched_phyla = set(target_phyla) - set(phylum_genomes_map.keys())

    if not unmatched_phyla:
        return pd.DataFrame(columns=['phylum', 'domain', 'total_n_genomes'])

    csv_normalized = phylum_data_df.copy()
    csv_normalized['phylum_normalized'] = csv_normalized[phylum_col].apply(normalize_phylum_name)

    # Ensure n_genomes is numeric before summing
    csv_normalized[n_genomes_col] = pd.to_numeric(csv_normalized[n_genomes_col], errors='coerce').fillna(0)

    unmatched_df = (
        csv_normalized[
            csv_normalized['phylum_normalized'].isin(unmatched_phyla)
        ]
        .dropna(subset=['phylum_normalized'])
        .groupby('phylum_normalized', as_index=False)
        .agg(
            domain=(domain_col, 'first'),
            total_n_genomes=(n_genomes_col, 'sum')
        )
        .rename(columns={'phylum_normalized': 'phylum'})
        [['phylum', 'domain', 'total_n_genomes']]
        .sort_values(['total_n_genomes', 'phylum'], ascending=[False, True])
        .reset_index(drop=True)
    )

    return unmatched_df


def save_outputs(outdir, target_phyla, representative_genomes, tip_metadata_df,
                all_mappings_df, unmatched_phyla_df, original_tips, pruned_tips):
    """
    Save all output files.
    
    Args:
        outdir: Output directory
        target_phyla: Set of target phylum names
        representative_genomes: dict mapping phylum -> representative_genome_id
        tip_metadata_df: DataFrame with per-tip metadata
        all_mappings_df: DataFrame with phylum -> all genomes mapping
        unmatched_phyla_df: DataFrame with unmatched phyla and CSV domain
        original_tips: Number of tips in original tree
        pruned_tips: Number of tips in pruned tree
    """
    print(f"\n{'='*60}")
    print("STEP 5: Saving outputs")
    print(f"{'='*60}")
    
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)
    
    # 1. Save target phyla list
    phyla_file = outdir / 'target_phyla.txt'
    with open(phyla_file, 'w') as f:
        for phylum in sorted(target_phyla):
            f.write(f"{phylum}\n")
    print(f"✓ Saved {phyla_file}")
    
    # 2. Save representative genomes mapping
    rep_genomes_file = outdir / 'representative_genomes.tsv'
    rep_df = pd.DataFrame([
        {'phylum': phylum, 'representative_genome_id': genome_id}
        for phylum, genome_id in representative_genomes.items()
    ]).sort_values('phylum')
    rep_df.to_csv(rep_genomes_file, sep='\t', index=False)
    print(f"✓ Saved {rep_genomes_file}")
    
    # 3. Save all phylum-to-genomes mapping
    mapping_file = outdir / 'phylum_to_genomes.tsv'
    all_mappings_df.to_csv(mapping_file, sep='\t', index=False)
    print(f"✓ Saved {mapping_file}")
    
    # 4. Save tip metadata
    metadata_file = outdir / 'representative_tip_metadata.tsv'
    tip_metadata_df.to_csv(metadata_file, sep='\t', index=False)
    print(f"✓ Saved {metadata_file}")

    # 5. Save unmatched phyla table
    unmatched_file = outdir / 'unmatched_phyla.tsv'
    unmatched_phyla_df.to_csv(unmatched_file, sep='\t', index=False)
    if len(unmatched_phyla_df) > 0:
        print(f"⚠ {len(unmatched_phyla_df)} phyla from the CSV were not found in bac120 taxonomy")
        print(f"Details saved to {unmatched_file}")
    else:
        print(f"✓ Saved {unmatched_file} (0 unmatched phyla)")
    
    # 6. Create summary report
    report_file = outdir / 'pruning_report.txt'
    with open(report_file, 'w') as f:
        f.write("GTDB Phylum-Level Representative Tree Report\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"Target phyla: {len(target_phyla)}\n")
        f.write(f"Representative genomes selected: {len(representative_genomes)}\n")
        f.write(f"Original tree tips (genomes): {original_tips:,}\n")
        f.write(f"Pruned tree tips (representatives): {pruned_tips}\n")
        f.write(f"Reduction: {original_tips:,} -> {pruned_tips} tips\n\n")
        
        f.write("Representative genomes per phylum:\n")
        for phylum in sorted(representative_genomes.keys()):
            f.write(f"  {phylum}: {representative_genomes[phylum]}\n")
    
    print(f"✓ Saved {report_file}")
    print(f"\n✓ All outputs saved to: {outdir}")


def main():
    """Main pipeline function."""
    parser = argparse.ArgumentParser(
        description='Prune GTDB tree by phyla from a phylum-level CSV file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Required arguments
    parser.add_argument('--csv', required=True,
                       help='Path to phylum-level CSV file')
    parser.add_argument('--tree', required=True,
                       help='Path to GTDB tree in Newick format')
    parser.add_argument('--taxonomy', required=True,
                       help='Path to GTDB taxonomy TSV file')
    parser.add_argument('--outdir', required=True,
                       help='Output directory for results')
    
    # Optional arguments for column customization
    parser.add_argument('--phylum-col', default='phylum',
                       help='Name of phylum column in CSV (default: phylum)')
    parser.add_argument('--n-genomes-col', default='n_genomes',
                       help='Name of n_genomes column in CSV (default: n_genomes)')
    parser.add_argument('--electron-donor-col', default='electron_donor',
                       help='Name of electron_donor column in CSV (default: electron_donor)')
    
    args = parser.parse_args()
    
    # Validate input files exist
    for path, name in [(args.csv, 'CSV'), (args.tree, 'Tree'), (args.taxonomy, 'Taxonomy')]:
        if not os.path.exists(path):
            print(f"✗ ERROR: {name} file not found: {path}")
            sys.exit(1)
    
    print("\n" + "="*60)
    print("GTDB Phylum-Level Representative Tree Pipeline")
    print("="*60)
    print(f"Input CSV: {args.csv}")
    print(f"Input tree: {args.tree}")
    print(f"Input taxonomy: {args.taxonomy}")
    print(f"Output directory: {args.outdir}")
    
    # Ensure output directory exists early
    outdir_path = Path(args.outdir)
    outdir_path.mkdir(exist_ok=True, parents=True)
    
    # Step 1: Read target phyla from CSV
    target_phyla, phylum_data_df = read_target_phyla(args.csv, args.phylum_col)
    
    # Step 2: Map all genomes to phyla using taxonomy
    phylum_genomes_map = read_taxonomy_and_map_genomes(args.taxonomy, target_phyla)

    # Build unmatched phyla table (CSV-only phyla not found in bac120 taxonomy)
    unmatched_phyla_df = build_unmatched_phyla_table(
        phylum_data_df,
        target_phyla,
        phylum_genomes_map,
        phylum_col=args.phylum_col,
        domain_col=DEFAULT_COLUMNS['domain'],
        n_genomes_col=args.n_genomes_col
    )
    
    # Get tree tip IDs to check which phyla are present in the tree
    print(f"\n{'='*60}")
    print("STEP 2.5: Getting tree tip IDs")
    print(f"{'='*60}")
    tree_tip_ids = get_tree_tip_ids(args.tree)
    print(f"✓ Tree has {len(tree_tip_ids):,} tips")
    
    # Step 2b: Select representative genomes (one per phylum)
    representative_genomes, all_mappings_df = select_representative_genomes(
        phylum_genomes_map, 
        tree_tip_ids
    )
    
    # Step 3: Prune tree to representatives
    representative_genome_ids = set(representative_genomes.values())
    pruned_tree_path = outdir_path / 'phylum_representative.tree'
    original_tips, pruned_tips = prune_tree_to_representatives(
        args.tree, 
        representative_genome_ids, 
        pruned_tree_path
    )
    
    # Step 4: Create tip metadata
    tip_metadata_df = create_tip_metadata(
        representative_genomes, 
        phylum_data_df,
        args.phylum_col
    )
    
    # Step 5: Save all outputs
    save_outputs(
        args.outdir,
        target_phyla,
        representative_genomes,
        tip_metadata_df,
        all_mappings_df,
        unmatched_phyla_df,
        original_tips,
        pruned_tips
    )
    
    print(f"\n{'='*60}")
    print("✓ Pipeline completed successfully!")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
