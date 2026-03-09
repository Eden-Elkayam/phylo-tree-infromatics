#!/usr/bin/env python3
"""Build a phylum-representative GTDB tree and companion metadata tables."""

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd

DEFAULT_COLUMNS = {
    'phylum': 'phylum',
    'n_genomes': 'n_genomes',
    'electron_donor': 'electron_donor',
    'domain': 'domain',
}

TAXONOMIC_LEVELS = [
    'species',
    'genus',
    'family',
    'order',
    'class',
    'phylum',
    'domain',
]

RANK_PREFIX = {
    'domain': 'd__',
    'phylum': 'p__',
    'class': 'c__',
    'order': 'o__',
    'family': 'f__',
    'genus': 'g__',
    'species': 's__',
}


def step(title):
    print(f"\n{'=' * 60}")
    print(title)
    print(f"{'=' * 60}")


def fail(msg):
    print(msg)
    sys.exit(1)


def validate_input_file(path, label):
    if not os.path.exists(path):
        fail(f"✗ ERROR: {label} file not found: {path}")


def write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False)
    print(f"✓ Saved {path}")


def write_text_lines(lines, path):
    with open(path, 'w') as f:
        for line in lines:
            f.write(f"{line}\n")
    print(f"✓ Saved {path}")


def normalize_taxon_name(value, rank):
    if pd.isna(value):
        return None
    s = str(value).strip()
    prefix = RANK_PREFIX[rank]
    return s[3:] if s.startswith(prefix) else s


def extract_rank_from_taxonomy(taxonomy, rank):
    if pd.isna(taxonomy):
        return None
    prefix = RANK_PREFIX[rank]
    for field in str(taxonomy).split(';'):
        if field.startswith(prefix):
            return normalize_taxon_name(field, rank)
    return None


def build_representative_genomes_df(representative_genomes, target_rank='phylum'):
    return pd.DataFrame(
        [{target_rank: p, 'representative_genome_id': g} for p, g in representative_genomes.items()]
    )


def _load_tree_with_backend(tree_path, fail_hard=True):
    try:
        from ete3 import Tree
        return 'ete3', Tree(tree_path, format=1)
    except ImportError:
        pass

    try:
        from Bio import Phylo
        return 'biopython', Phylo.read(tree_path, 'newick')
    except ImportError:
        if not fail_hard:
            return None, None
        fail(
            "✗ ERROR: Neither ete3 nor BioPython is installed!\n"
            "\nPlease install one of these packages:\n"
            "  pip install ete3\n"
            "  pip install biopython"
        )


def _tip_ids_from_tree(backend, tree):
    return {leaf.name for leaf in tree.get_leaves()} if backend == 'ete3' else {t.name for t in tree.get_terminals()}


def _prune_loaded_tree(backend, tree, keep_ids):
    if backend == 'ete3':
        tree.prune(keep_ids, preserve_branch_length=True)
        return len(tree.get_leaves())

    terminals_to_remove = [t for t in tree.get_terminals() if t.name not in keep_ids]
    for t in terminals_to_remove:
        tree.prune(t)
    return len(tree.get_terminals())


def _clear_internal_node_labels(backend, tree):
    if backend == 'ete3':
        for node in tree.traverse():
            if not node.is_leaf():
                node.name = ''
        return

    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = None


def _write_loaded_tree(backend, tree, out_path):
    if backend == 'ete3':
        tree.write(format=1, outfile=out_path)
    else:
        from Bio import Phylo
        Phylo.write(tree, str(out_path), 'newick')


def _prune_tree(tree_path, keep_ids, out_path, step_title, suffix='', include_tip_ids=False):
    step(step_title)
    print(f"Reading tree from: {tree_path}")

    backend, tree = _load_tree_with_backend(tree_path)
    if backend == 'biopython':
        print("⚠ ete3 not available, trying BioPython...")

    original_tips = len(_tip_ids_from_tree(backend, tree))
    print(f"✓ Tree loaded with {original_tips:,} tips")

    pruned_tips = _prune_loaded_tree(backend, tree, keep_ids)
    _clear_internal_node_labels(backend, tree)
    tip_ids = _tip_ids_from_tree(backend, tree)
    _write_loaded_tree(backend, tree, out_path)

    print(f"✓ Pruned tree saved with {pruned_tips:,} tips{suffix}")
    return (original_tips, pruned_tips, tip_ids) if include_tip_ids else (original_tips, pruned_tips)


def get_tree_tip_ids(tree_path):
    backend, tree = _load_tree_with_backend(tree_path, fail_hard=False)
    if backend is not None:
        return _tip_ids_from_tree(backend, tree)

    print("⚠ Warning: Tree libraries not available, using simple Newick parsing...")
    with open(tree_path, 'r') as f:
        content = f.read()
    return set(re.findall(r'([A-Za-z0-9_\.]+)(?:[:\),]|$)', content))


def read_target_phyla(csv_path, phylum_col='phylum'):
    step("STEP 1: Reading phylum-level CSV")

    df = pd.read_csv(csv_path)
    print(f"✓ Read CSV with {len(df)} rows and {len(df.columns)} columns")

    if phylum_col in df.columns:
        target_rank = phylum_col
    else:
        target_rank = next((lvl for lvl in TAXONOMIC_LEVELS if lvl in df.columns), None)
        if target_rank is None:
            fail(
                "\n✗ ERROR: No supported taxonomy level column found in CSV. "
                f"Expected one of: {', '.join(TAXONOMIC_LEVELS)}\n"
                f"Available columns: {', '.join(df.columns)}"
            )

    phyla = {normalize_taxon_name(p, target_rank) for p in df[target_rank].dropna().unique()}
    phyla.discard(None)

    print(f"✓ Using taxonomy level: {target_rank}")
    print(f"✓ Found {len(phyla)} unique targets")
    print(f"  Examples: {', '.join(list(phyla)[:5])}")
    return phyla, df, target_rank


def read_taxonomy_and_map_genomes(taxonomy_path, target_phyla, target_rank):
    step("STEP 2: Reading taxonomy and mapping genomes to phyla")

    phylum_to_genomes = {}
    total_genomes = 0
    matched_genomes = 0

    print("Processing taxonomy file (may take a moment for large files)...")
    with open(taxonomy_path, 'r') as f:
        for line in f:
            total_genomes += 1
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            genome_id, taxonomy = parts[0], parts[1]
            phylum = extract_rank_from_taxonomy(taxonomy, target_rank)
            if phylum in target_phyla:
                phylum_to_genomes.setdefault(phylum, []).append(genome_id)
                matched_genomes += 1

            if total_genomes % 100000 == 0:
                print(f"  Processed {total_genomes:,} genomes, found {matched_genomes:,} matches...")

    print(f"\n✓ Processed {total_genomes:,} total genomes")
    print(f"✓ Found {matched_genomes:,} genomes from {len(phylum_to_genomes)} target phyla")

    unmatched = target_phyla - set(phylum_to_genomes.keys())
    if unmatched:
        print(f"\n⚠ WARNING: {len(unmatched)} phyla from CSV not found in taxonomy:")
        for p in sorted(unmatched)[:10]:
            print(f"  - {p}")
        if len(unmatched) > 10:
            print(f"  ... and {len(unmatched) - 10} more")

    if not phylum_to_genomes:
        fail("\n✗ ERROR: No genomes matched! Check phylum name formatting.")

    return phylum_to_genomes


def select_representative_genomes(phylum_genomes_map, tree_tip_ids):
    step("STEP 2b: Selecting representative genomes per phylum")

    representative_genomes = {}
    missing_in_tree = []
    all_mappings = []

    for phylum in sorted(phylum_genomes_map):
        genomes = sorted(phylum_genomes_map[phylum])
        all_mappings.append({
            'phylum': phylum,
            'total_genomes_for_phylum': len(genomes),
            'all_genome_ids': ';'.join(genomes),
        })

        rep = next((g for g in genomes if g in tree_tip_ids), None)
        if rep:
            representative_genomes[phylum] = rep
        else:
            missing_in_tree.append(phylum)

    print(f"✓ Selected {len(representative_genomes)} representative genomes")
    print(f"  Found in tree: {len(representative_genomes)} phyla")

    if missing_in_tree:
        print(f"\n⚠ WARNING: {len(missing_in_tree)} phyla have no genomes in tree:")
        for p in missing_in_tree[:5]:
            print(f"  - {p}")
        if len(missing_in_tree) > 5:
            print(f"  ... and {len(missing_in_tree) - 5} more")

    return representative_genomes, pd.DataFrame(all_mappings)


def prune_tree_to_representatives(tree_path, representative_genome_ids, output_path):
    return _prune_tree(
        tree_path,
        representative_genome_ids,
        output_path,
        step_title="STEP 3: Pruning tree to representative genomes",
        suffix=" (one per phylum)",
    )


def prune_tree_simple(tree_path, retained_genome_ids, output_path):
    return _prune_tree(
        tree_path,
        retained_genome_ids,
        output_path,
        step_title="STEP 3: Pruning tree",
        include_tip_ids=True,
    )


def create_tip_metadata(representative_genomes, phylum_data_df, phylum_col='phylum'):
    step("STEP 4: Creating representative tip metadata")

    rep_df = build_representative_genomes_df(representative_genomes, phylum_col).rename(
        columns={'representative_genome_id': 'genome_id'}
    )[['genome_id', phylum_col]]

    data = phylum_data_df.copy()
    data[phylum_col] = data[phylum_col].apply(lambda x: normalize_taxon_name(x, phylum_col))

    tip_metadata = rep_df.merge(data, on=phylum_col, how='left')
    if f'{phylum_col}_x' in tip_metadata.columns:
        tip_metadata = tip_metadata.drop(columns=[f'{phylum_col}_y']).rename(
            columns={f'{phylum_col}_x': phylum_col}
        )

    print(f"✓ Created metadata table with {len(tip_metadata)} rows")
    print(f"  Phyla represented: {tip_metadata[phylum_col].nunique()}")
    print(f"  Columns: {', '.join(tip_metadata.columns[:8])}")
    if len(tip_metadata.columns) > 8:
        print(f"  ... and {len(tip_metadata.columns) - 8} more")

    return tip_metadata


def build_unmatched_phyla_table(phylum_data_df, target_phyla, phylum_genomes_map,
                                phylum_col='phylum', domain_col='domain', n_genomes_col='n_genomes'):
    unmatched = set(target_phyla) - set(phylum_genomes_map.keys())
    if not unmatched:
        return pd.DataFrame(columns=[phylum_col, 'domain', 'total_n_genomes'])

    df = phylum_data_df.copy()
    df['phylum_normalized'] = df[phylum_col].apply(lambda x: normalize_taxon_name(x, phylum_col))
    df[n_genomes_col] = pd.to_numeric(df[n_genomes_col], errors='coerce').fillna(0)

    return (
        df[df['phylum_normalized'].isin(unmatched)]
        .dropna(subset=['phylum_normalized'])
        .groupby('phylum_normalized', as_index=False)
        .agg(domain=(domain_col, 'first'), total_n_genomes=(n_genomes_col, 'sum'))
        .rename(columns={'phylum_normalized': phylum_col})
        [[phylum_col, 'domain', 'total_n_genomes']]
        .sort_values(['total_n_genomes', phylum_col], ascending=[False, True])
        .reset_index(drop=True)
    )


def build_report_text(target_phyla, representative_genomes, original_tips, pruned_tips):
    lines = [
        "GTDB Phylum-Level Representative Tree Report",
        "=" * 60,
        "",
        f"Target phyla: {len(target_phyla)}",
        f"Representative genomes selected: {len(representative_genomes)}",
        f"Original tree tips (genomes): {original_tips:,}",
        f"Pruned tree tips (representatives): {pruned_tips}",
        f"Reduction: {original_tips:,} -> {pruned_tips} tips",
        "",
        "Representative genomes per phylum:",
    ]
    lines.extend([f"  {p}: {representative_genomes[p]}" for p in sorted(representative_genomes)])
    return "\n".join(lines) + "\n"


def save_outputs(outdir, target_phyla, representative_genomes, tip_metadata_df,
                 all_mappings_df, unmatched_phyla_df, original_tips, pruned_tips, target_rank='phylum'):
    step("STEP 5: Saving outputs")

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    write_text_lines(sorted(target_phyla), outdir / 'target_phyla.txt')
    write_tsv(build_representative_genomes_df(representative_genomes, target_rank).sort_values(target_rank),
              outdir / 'representative_genomes.tsv')
    write_tsv(all_mappings_df, outdir / f'{target_rank}_to_genomes.tsv')
    write_tsv(tip_metadata_df, outdir / 'representative_tip_metadata.tsv')

    unmatched_file = outdir / 'unmatched_phyla.tsv'
    write_tsv(unmatched_phyla_df, unmatched_file)
    if len(unmatched_phyla_df) > 0:
        print(f"⚠ {len(unmatched_phyla_df)} phyla from the CSV were not found in bac120 taxonomy")
        print(f"Details saved to {unmatched_file}")

    report_file = outdir / 'pruning_report.txt'
    with open(report_file, 'w') as f:
        f.write(build_report_text(target_phyla, representative_genomes, original_tips, pruned_tips))
    print(f"✓ Saved {report_file}")

    print(f"\n✓ All outputs saved to: {outdir}")


def main():
    parser = argparse.ArgumentParser(
        description='Prune GTDB tree by phyla from a phylum-level CSV file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument('--csv', required=True, help='Path to phylum-level CSV file')
    parser.add_argument('--tree', required=True, help='Path to GTDB tree in Newick format')
    parser.add_argument('--taxonomy', required=True, help='Path to GTDB taxonomy TSV file')
    parser.add_argument('--outdir', required=True, help='Output directory for results')

    parser.add_argument('--phylum-col', default='phylum', help='Name of phylum column in CSV (default: phylum)')
    parser.add_argument('--n-genomes-col', default='n_genomes', help='Name of n_genomes column in CSV (default: n_genomes)')
    parser.add_argument('--electron-donor-col', default='electron_donor', help='Name of electron_donor column in CSV (default: electron_donor)')

    args = parser.parse_args()

    for path, label in [(args.csv, 'CSV'), (args.tree, 'Tree'), (args.taxonomy, 'Taxonomy')]:
        validate_input_file(path, label)

    print("\n" + "=" * 60)
    print("GTDB Phylum-Level Representative Tree Pipeline")
    print("=" * 60)
    print(f"Input CSV: {args.csv}")
    print(f"Input tree: {args.tree}")
    print(f"Input taxonomy: {args.taxonomy}")
    print(f"Output directory: {args.outdir}")

    outdir_path = Path(args.outdir)
    outdir_path.mkdir(exist_ok=True, parents=True)

    target_phyla, phylum_data_df, target_rank = read_target_phyla(args.csv, args.phylum_col)
    phylum_genomes_map = read_taxonomy_and_map_genomes(args.taxonomy, target_phyla, target_rank)

    unmatched_phyla_df = build_unmatched_phyla_table(
        phylum_data_df,
        target_phyla,
        phylum_genomes_map,
        phylum_col=target_rank,
        domain_col=DEFAULT_COLUMNS['domain'],
        n_genomes_col=args.n_genomes_col,
    )

    step("STEP 2.5: Getting tree tip IDs")
    tree_tip_ids = get_tree_tip_ids(args.tree)
    print(f"✓ Tree has {len(tree_tip_ids):,} tips")

    representative_genomes, all_mappings_df = select_representative_genomes(phylum_genomes_map, tree_tip_ids)

    representative_genome_ids = set(representative_genomes.values())
    pruned_tree_path = outdir_path / f'{target_rank}_representative.tree'
    original_tips, pruned_tips = prune_tree_to_representatives(
        args.tree,
        representative_genome_ids,
        pruned_tree_path,
    )

    tip_metadata_df = create_tip_metadata(representative_genomes, phylum_data_df, target_rank)

    save_outputs(
        args.outdir,
        target_phyla,
        representative_genomes,
        tip_metadata_df,
        all_mappings_df,
        unmatched_phyla_df,
        original_tips,
        pruned_tips,
        target_rank,
    )

    print(f"\n{'=' * 60}")
    print("✓ Pipeline completed successfully!")
    print(f"{'=' * 60}\n")


if __name__ == '__main__':
    main()
