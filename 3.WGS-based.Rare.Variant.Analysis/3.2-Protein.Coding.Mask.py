#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge VEP annotations with precomputed SpliceAI scores and generate
BRaVa-style variant annotations together with a SAIGE group file.

Purpose
-------
1. Read VEP annotations for protein-coding variants.
2. Merge per-variant SpliceAI scores from Illumina precomputed SNV and INDEL VCFs.
3. Optionally supplement missing CADD scores for indels.
4. Classify variants into BRaVa-style functional groups.(reference:https://github.com/BRaVa-genetics/variant-annotation?tab=readme-ov-file)
5. Export:
   - long-format annotation table
   - SAIGE group file

Notes
-----
- Supports site-level variant annotation merging.
- Uses tabix streaming by contig for memory-efficient SpliceAI parsing.
- Allows optional parallel parsing across contigs.

Author: Shelley
Date: 2026-03-18
"""

from __future__ import annotations
import argparse
import re
import shlex
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count

# -----------------------------
# Utils and logging
# -----------------------------

def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def eprint(*args, **kwargs):
    print(f"[{ts()}]", *args, file=sys.stderr, **kwargs)

def extract_tail_snpid(s: str) -> str:
    """From 'DRAGEN:chr20:60001:T:A' (or any prefix:xxx:xxx:xxx:xxx) -> 'chr20:60001:T:A'."""
    if not isinstance(s, str):
        return s
    parts = s.split(':')
    if len(parts) >= 4:
        tail = ':'.join(parts[-4:])
        if not tail.split(':')[0].startswith('chr'):
            t = tail.split(':')
            t[0] = 'chr' + t[0]
            tail = ':'.join(t)
        return tail
    return s

def add_chr_prefix(s: str) -> str:
    return s if s.startswith('chr') else f'chr{s}'

def run_and_capture(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True).strip()

def list_contigs(vcf_path: str) -> list[str]:
    """List contigs from a tabix-indexed VCF; raises if no index."""
    try:
        out = run_and_capture(['tabix', '-l', vcf_path])
    except subprocess.CalledProcessError as ex:
        raise RuntimeError(f"tabix failed on {vcf_path}. Is it bgzipped and indexed (.tbi present)? {ex}")
    contigs = [c for c in out.splitlines() if c]
    if not contigs:
        raise RuntimeError(f"No contigs found via tabix -l for {vcf_path}.")
    return contigs

def normalize_region_arg(contigs: list[str], contig: str) -> str:
    """Return a contig name that exists in the file, accepting '22' or 'chr22'."""
    if contig in contigs:
        return contig
    if contig.startswith('chr') and contig[3:] in contigs:
        return contig[3:]
    if ('chr' + contig) in contigs:
        return 'chr' + contig
    raise RuntimeError(f"Contig '{contig}' not in VCF index. Available: {', '.join(contigs[:10])} ...")

# -----------------------------
# SpliceAI parsing
# -----------------------------

def parse_spliceai_info(info: str) -> List[Tuple[str, float]]:
    """
    Parse SpliceAI INFO -> list of (gene_symbol_or_ENSG, DSmax).
    Per transcript: DSmax = max(DS_AG, DS_AL, DS_DG, DS_DL)
    Collapse to max per gene. Supports Illumina v1.3 and older '---' layouts.
    """
    if not isinstance(info, str) or info == '.' or not info.strip():
        return []
    if info.startswith('SpliceAI='):
        info = info.split('=', 1)[1]

    per_gene: defaultdict[str, float] = defaultdict(float)

    for rec in info.split(','):
        parts = rec.split('|')
        if len(parts) < 6:
            continue

        ds_vals = []
        for x in parts[2:6]:
            try:
                ds_vals.append(float(x))
            except Exception:
                pass
        if not ds_vals:
            continue
        ds_max = float(np.max(ds_vals))

        gene_field = parts[1]
        if '---' in gene_field:
            gene_id = None
            for chunk in gene_field.split('---'):
                if chunk.startswith('ENSG'):
                    gene_id = chunk.split('.')[0]
                    break
            gene = gene_id if gene_id else gene_field.split('---')[0]
        else:
            gene = gene_field

        if ds_max > per_gene[gene]:
            per_gene[gene] = ds_max

    return list(per_gene.items())

_re_splice_payload = re.compile(r'SpliceAI=([^;\t]+)')

def read_spliceai_region(path: str, region: str) -> pd.DataFrame:
    """
    Stream a single contig or region from a tabix-indexed SpliceAI VCF.
    Returns DataFrame ['ID','SP_GENE','SpliceAI_max'] collapsed per (variant, gene).
    """
    cmd = f"tabix {shlex.quote(path)} {shlex.quote(region)}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, errors='replace')
    if proc.stdout is None:
        raise RuntimeError(f"Failed to open {path} region {region} with tabix")

    ids, genes, dsvals = [], [], []

    for line in proc.stdout:
        if line.startswith('#'):
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 8:
            continue

        chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]
        chrom = chrom if chrom.startswith('chr') else ('chr' + chrom)
        var_id = f"{chrom}:{pos}:{ref}:{alt}"

        info_field = parts[7]
        m = _re_splice_payload.search(info_field)
        if not m:
            continue

        for (g, ds) in parse_spliceai_info(m.group(1)):
            if isinstance(g, str) and g.startswith('ENSG'):
                g = g.split('.')[0]
            ids.append(var_id)
            genes.append(g)
            dsvals.append(ds)

    proc.stdout.close()
    proc.wait()

    if not ids:
        return pd.DataFrame(columns=['ID', 'SP_GENE', 'SpliceAI_max'])

    sp = pd.DataFrame({'ID': ids, 'SP_GENE': genes, 'SpliceAI_max': dsvals})
    sp['SpliceAI_max'] = pd.to_numeric(sp['SpliceAI_max'], errors='coerce')
    sp = sp.groupby(['ID', 'SP_GENE'], as_index=False, sort=False)['SpliceAI_max'].max()
    return sp

def read_spliceai_all_by_contig(path: str, n_jobs: int = 1) -> pd.DataFrame:
    """Iterate all contigs via tabix; optionally parallelize by contig."""
    contigs = list_contigs(path)
    eprint(f"[INFO] SpliceAI contigs: {len(contigs)}")

    if n_jobs and n_jobs > 1:
        n_jobs = min(n_jobs, cpu_count())
        eprint(f"[INFO] Parallel parsing with {n_jobs} workers ...")
        with Pool(processes=n_jobs) as pool:
            parts = pool.starmap(read_spliceai_region, [(path, c) for c in contigs])
    else:
        parts = []
        for c in contigs:
            eprint(f"[INFO] Parsing contig {c} ...")
            parts.append(read_spliceai_region(path, c))

    if not parts:
        return pd.DataFrame(columns=['ID', 'SP_GENE', 'SpliceAI_max'])

    sp = pd.concat(parts, ignore_index=True)
    sp = sp.groupby(['ID', 'SP_GENE'], as_index=False, sort=False)['SpliceAI_max'].max()
    return sp

# -----------------------------
# Classification (BRaVa-like)
# -----------------------------

PLOF_CSQS = {
    'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
    'stop_gained', 'frameshift_variant'
}
MISSENSE_CSQS = {
    'stop_lost', 'start_lost', 'transcript_amplification',
    'inframe_insertion', 'inframe_deletion', 'missense_variant',
    'protein_altering_variant'
}
SYNONYMOUS_CSQS = {'stop_retained_variant', 'synonymous_variant'}
OTHER_CSQS = {
    'mature_miRNA_variant', '5_prime_UTR_variant', '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant', 'intron_variant', 'NMD_transcript_variant',
    'non_coding_transcript_variant', 'upstream_gene_variant', 'downstream_gene_variant',
    'TFBS_ablation', 'TFBS_amplification', 'TF_binding_site_variant',
    'regulatory_region_ablation', 'regulatory_region_amplification', 'feature_elongation',
    'regulatory_region_variant', 'feature_truncation', 'intergenic_variant'
}
INFRAME_CSQS = {'inframe_deletion', 'inframe_insertion'}

def classify_variant(row: pd.Series, revel_cut: float, cadd_cut: float, spliceai_cut: float) -> str:
    csq = str(row.get('CSQ', '')) if pd.notna(row.get('CSQ')) else ''
    consequences = set(csq.split('&')) if csq else set()
    is_plof = row.get('LOF', '') == 'HC'
    is_lof_low = row.get('LOF', '') == 'LC'
    missense = len(consequences & MISSENSE_CSQS) > 0
    synonymous = len(consequences & SYNONYMOUS_CSQS) > 0
    other = len(consequences & OTHER_CSQS) > 0
    inframe = len(consequences & INFRAME_CSQS) > 0

    revel = row.get('REVEL_SCORE', np.nan)
    cadd = row.get('CADD_PHRED', np.nan)
    spmax = row.get('SpliceAI_max', np.nan)

    if is_plof:
        return 'pLoF'
    if missense and ((pd.notna(revel) and revel >= revel_cut) or (pd.notna(cadd) and cadd >= cadd_cut)):
        return 'damaging_missense_or_protein_altering'
    if pd.notna(spmax) and spmax >= spliceai_cut:
        return 'damaging_missense_or_protein_altering'
    if is_lof_low:
        return 'damaging_missense_or_protein_altering'
    if missense or inframe:
        return 'other_missense_or_protein_altering'
    if other:
        return 'non_coding'
    if synonymous:
        return 'synonymous'
    return ''

# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description='Merge VEP and SpliceAI to generate annotations and a SAIGE group file')
    ap.add_argument('--vep', '-v', required=True, help='VEP table (whitespace-delimited). SNP_ID may be DRAGEN:chr:pos:ref:alt')
    ap.add_argument('--spliceai-snv', required=True, help='Illumina precomputed SNV SpliceAI VCF (.vcf.gz)')
    ap.add_argument('--spliceai-indel', required=True, help='Illumina precomputed INDEL SpliceAI VCF (.vcf.gz)')
    ap.add_argument('--out', '-o', required=True, help='Output prefix')
    ap.add_argument('--cadd-indels', default=None, help='Combined CADD indel file (.tsv or .tsv.gz): Chrom Pos Ref Alt RawScore PHRED')

    ap.add_argument('--contig', default=None, help='Chromosome or contig to parse only (e.g., 22 or chr22)')
    ap.add_argument('--n-jobs', type=int, default=1, help='Parallel workers when scanning all contigs')

    ap.add_argument('--col-snpid', default='SNP_ID')
    ap.add_argument('--col-gene', default='GENE')
    ap.add_argument('--col-lof', default='LOF')
    ap.add_argument('--col-revel', default='REVEL_SCORE')
    ap.add_argument('--col-cadd', default='CADD_PHRED')
    ap.add_argument('--col-csq', default='CSQ')
    ap.add_argument('--col-canonical', default='CANONICAL')
    ap.add_argument('--col-biotype', default='BIOTYPE')
    ap.add_argument('--col-mane', default='MANE_SELECT')

    ap.add_argument('--revel-cut', type=float, default=0.773)
    ap.add_argument('--cadd-cut', type=float, default=28.1)
    ap.add_argument('--spliceai-cut', type=float, default=0.20)

    args = ap.parse_args()

    CN = dict(
        SNP=args.col_snpid,
        GENE=args.col_gene,
        LOF=args.col_lof,
        REVEL=args.col_revel,
        CADD=args.col_cadd,
        CSQ=args.col_csq,
        CANON=args.col_canonical,
        BIOTYPE=args.col_biotype,
        MANE=args.col_mane
    )

    # 1) Read VEP
    usecols = [CN['SNP'], CN['GENE'], CN['LOF'], CN['REVEL'], CN['CADD'], CN['CSQ'], CN['CANON'], CN['BIOTYPE'], CN['MANE']]
    eprint("[INFO] Reading VEP ...")
    vep_df = pd.read_csv(
        args.vep, sep=r"\s+", usecols=usecols, na_values='.', dtype=str, engine="python", encoding='cp1252'
    )
    eprint(f"[INFO] VEP rows: {len(vep_df):,}")

    vep_df['SNP_ID_final'] = vep_df[CN['SNP']].astype(str)
    vep_df['SNP_ID_merge'] = vep_df['SNP_ID_final'].map(extract_tail_snpid)

    total_rows = len(vep_df)
    keep = (vep_df[CN['BIOTYPE']] == 'protein_coding') & (
        vep_df[CN['MANE']].notna() | ((vep_df[CN['MANE']].isna()) & (vep_df[CN['CANON']] == 'YES'))
    )
    vep_df = vep_df.loc[keep].copy()
    eprint(f"[INFO] protein_coding + MANE/CANONICAL: kept {len(vep_df):,} / {total_rows:,}")

    def parse_revel(cell: Optional[str]) -> float:
        if cell is None or (isinstance(cell, float) and np.isnan(cell)):
            return np.nan
        for p in str(cell).split('&'):
            if p and p != '.':
                try:
                    return float(p)
                except ValueError:
                    pass
        return np.nan

    def to_float(cell):
        try:
            return float(cell)
        except Exception:
            return np.nan

    eprint("[INFO] Parsing REVEL and CADD fields ...")
    vep_df[CN['REVEL']] = vep_df[CN['REVEL']].map(parse_revel)
    vep_df[CN['CADD']] = vep_df[CN['CADD']].map(to_float)
    eprint(f"[INFO] REVEL non-NA: {vep_df[CN['REVEL']].notna().sum():,} "
           f"(>=cut: {(vep_df[CN['REVEL']] >= args.revel_cut).sum():,})")

    # 2) Read SpliceAI
    eprint("[INFO] Reading SpliceAI (region-aware via tabix) ...")

    if args.contig:
        contigs = list_contigs(args.spliceai_snv)
        region = normalize_region_arg(contigs, args.contig)
        eprint(f"[INFO] SNV: using region {region}")
        sp_snv = read_spliceai_region(args.spliceai_snv, region)
    else:
        sp_snv = read_spliceai_all_by_contig(args.spliceai_snv, n_jobs=max(1, args.n_jobs))

    if args.contig:
        contigs_i = list_contigs(args.spliceai_indel)
        region_i = normalize_region_arg(contigs_i, args.contig)
        eprint(f"[INFO] INDEL: using region {region_i}")
        sp_indel = read_spliceai_region(args.spliceai_indel, region_i)
    else:
        sp_indel = read_spliceai_all_by_contig(args.spliceai_indel, n_jobs=max(1, args.n_jobs))

    sp_all = pd.concat([sp_snv, sp_indel], ignore_index=True)
    if sp_all.empty:
        raise RuntimeError("Parsed 0 SpliceAI rows from SNV and INDEL files. Check tabix index and INFO parsing.")
    eprint(f"[INFO] SpliceAI combined rows (per gene, collapsed): {len(sp_all):,}")

    # 3) Merge VEP and SpliceAI
    eprint("[INFO] Merging VEP with SpliceAI on SNP only (per-variant DSmax across genes) ...")
    per_variant = sp_all.groupby('ID', as_index=False)['SpliceAI_max'].max()
    merged = vep_df.merge(per_variant, left_on='SNP_ID_merge', right_on='ID', how='left')

    n_sp = merged['SpliceAI_max'].notna().sum()
    eprint(f"[CHECK] After SpliceAI merge, rows with SpliceAI_max: {n_sp:,} / {len(merged):,}")
    if n_sp == 0:
        raise RuntimeError("SpliceAI coverage is 0 after merge. ID normalization or contig selection may be off.")

    # 4) Optional CADD indel merge
    if args.cadd_indels:
        before = merged[CN['CADD']].notna().sum()
        eprint(f"[INFO] Reading combined CADD indels from {args.cadd_indels} ...")
        try:
            cadd = pd.read_csv(
                args.cadd_indels, sep=r'\t', comment='#', header=None,
                names=['Chrom', 'Pos', 'Ref', 'Alt', 'RawScore', 'PHRED'], dtype=str, engine='python'
            )
        except Exception as ex:
            raise RuntimeError(f"Failed reading {args.cadd_indels}: {ex}")
        if not cadd.empty:
            cadd['Chrom'] = cadd['Chrom'].astype(str).str.replace(r'^chr', '', regex=True)
            cadd['snpid'] = 'chr' + cadd['Chrom'] + ':' + cadd['Pos'].astype(str) + ':' + cadd['Ref'] + ':' + cadd['Alt']
            cadd = cadd[['snpid', 'PHRED']].dropna().drop_duplicates()
            cadd = cadd.rename(columns={'PHRED': 'PHRED_indel'})
            merged = merged.merge(cadd, left_on='SNP_ID_merge', right_on='snpid', how='left')
            merged[CN['CADD']] = merged[CN['CADD']].where(
                merged[CN['CADD']].notna(),
                pd.to_numeric(merged['PHRED_indel'], errors='coerce')
            )
            merged.drop(columns=['snpid', 'PHRED_indel'], inplace=True, errors='ignore')
        after = merged[CN['CADD']].notna().sum()
        eprint(f"[INFO] CADD PHRED non-NA: before {before:,} -> after {after:,}")

    # 5) Classify
    eprint("[INFO] Classifying variants ...")
    merged['SpliceAI_max'] = pd.to_numeric(merged['SpliceAI_max'], errors='coerce')
    merged['annotation'] = merged.apply(
        lambda r: classify_variant(
            r,
            args.revel_cut,
            28.1 if np.isnan(r.get('CADD_PHRED', np.nan)) else args.cadd_cut,
            args.spliceai_cut
        ),
        axis=1
    )

    kept = merged[merged['annotation'] != ''].copy()
    eprint(f"[INFO] Empty-annotation filter: kept {len(kept):,} / {len(merged):,}")
    eprint(kept['annotation'].value_counts())

    # 6) Outputs
    out_csv = args.out + '.long.csv.gz'
    eprint(f"[INFO] Writing long CSV: {out_csv}")
    cols_out = [
        'SNP_ID_final', 'SNP_ID_merge',
        CN['GENE'], CN['LOF'], CN['REVEL'], CN['CADD'], CN['CSQ'],
        'SpliceAI_max', 'annotation'
    ]
    cols_out = [c for c in cols_out if c in kept.columns]
    kept[cols_out].to_csv(out_csv, index=False, compression='gzip')

    out_grp = args.out
    eprint(f"[INFO] Writing SAIGE group file (using SNP_ID_final): {out_grp}")
    with open(out_grp, 'w') as fout:
        for gene, sub in kept.groupby(CN['GENE']):
            snps = ' '.join(sub['SNP_ID_final'].astype(str).tolist())
            annos = ' '.join(sub['annotation'].astype(str).tolist())
            fout.write(f"{gene} var {snps}\n")
            fout.write(f"{gene} anno {annos}\n")

    eprint("[INFO] Done.")

if __name__ == '__main__':
    main()