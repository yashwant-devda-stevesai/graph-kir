#!/usr/bin/env python3
"""
KIR Gene Typing Pipeline - Very Simple Version
You only give ONE FASTQ file path

Usage:
    python kir_typing_simple.py /path/to/your_sample.fastq.gz

    or

    python kir_typing_simple.py sample_R1.fastq
"""

import sys
import os
import argparse
from collections import defaultdict
import numpy as np
from Bio import SeqIO
import gzip


# ────────────────────────────────────────────────
# KIR Genes - basic info
# ────────────────────────────────────────────────
KIR_GENES = {
    'KIR2DL1': {'type': 'inhibitory'},
    'KIR2DL2': {'type': 'inhibitory'},
    'KIR2DL3': {'type': 'inhibitory'},
    'KIR2DL4': {'type': 'inhibitory'},
    'KIR2DL5': {'type': 'inhibitory'},
    'KIR2DS1': {'type': 'activating'},
    'KIR2DS2': {'type': 'activating'},
    'KIR2DS3': {'type': 'activating'},
    'KIR2DS4': {'type': 'activating'},
    'KIR2DS5': {'type': 'activating'},
    'KIR3DL1': {'type': 'inhibitory'},
    'KIR3DL2': {'type': 'inhibitory'},
    'KIR3DL3': {'type': 'inhibitory'},
    'KIR3DS1': {'type': 'activating'},
}

HAPLOTYPE_A_GENES = {'KIR2DL1', 'KIR2DL3', 'KIR2DL4', 'KIR2DS4', 'KIR3DL1', 'KIR3DL2', 'KIR3DL3'}
HAPLOTYPE_B_GENES = {'KIR2DL2', 'KIR2DL5', 'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS5', 'KIR3DS1'}


# ────────────────────────────────────────────────
# Read FASTQ (single file only)
# ────────────────────────────────────────────────
def read_fastq(filepath: str) -> list:
    sequences = []
    try:
        if filepath.lower().endswith('.gz'):
            open_func = gzip.open
        else:
            open_func = open

        print(f"  Reading file: {os.path.basename(filepath)}")
        with open_func(filepath, 'rt') as f:
            for record in SeqIO.parse(f, "fastq"):
                sequences.append(str(record.seq))

        print(f"  → Found {len(sequences):,} reads\n")
        return sequences

    except Exception as e:
        print(f"Error reading FASTQ file: {e}")
        sys.exit(1)


# ────────────────────────────────────────────────
# Very simple read-to-gene matching (simulation)
# ────────────────────────────────────────────────
def simple_kir_mapping(sequences):
    gene_hits = defaultdict(int)

    print("  Mapping reads to KIR genes (simulation)...")

    for i, seq in enumerate(sequences, 1):
        if i % 50000 == 0:
            print(f"    processed {i:,} reads")

        seq_short = seq[:60]  # only look at beginning for speed
        for gene in KIR_GENES:
            # very naive matching - just for demo
            if gene.lower() in seq_short.lower() or len(seq_short) > 55:
                gene_hits[gene] += 1

    print("  Mapping finished.\n")
    return dict(gene_hits)


# ────────────────────────────────────────────────
# Decide which genes are present + simple copy number
# ────────────────────────────────────────────────
def call_genes(gene_hits, min_reads=5):
    detected = {}

    for gene, count in gene_hits.items():
        if count < min_reads:
            continue

        # Very simple copy number estimation
        if count < min_reads * 1.8:
            cn = 1
        else:
            cn = 2

        detected[gene] = {
            'read_count': count,
            'copy_number': cn,
            'allele': f"*{cn:01d}0{count % 10}"   # fake allele
        }

    return detected


# ────────────────────────────────────────────────
# Very simple haplotype prediction
# ────────────────────────────────────────────────
def predict_haplotype(detected_genes):
    present = set(detected_genes.keys())

    a_count = len(present & HAPLOTYPE_A_GENES)
    b_count = len(present & HAPLOTYPE_B_GENES)

    if b_count > a_count + 1:
        hap = 'B'
    elif a_count > b_count:
        hap = 'A'
    else:
        hap = 'A/B'

    return hap


# ────────────────────────────────────────────────
# Print the exact table format you want
# ────────────────────────────────────────────────
def print_client_table(sample_id: str, haplotype: str, detected: dict):
    print("\n" + "═" * 80)
    print("         KIR TYPING RESULT – CLIENT REPORT")
    print("═" * 80)

    print(f" Sample ID          │ {sample_id}")
    print(f" Haplotype          │ {haplotype}")
    print("─" * 80)

    if not detected:
        print("   No KIR genes detected.")
        print("═" * 80 + "\n")
        return

    genes_list = []
    cn_list = []
    alleles_list = []

    for gene in sorted(detected.keys()):
        info = detected[gene]
        genes_list.append(gene)
        cn_list.append(str(info['copy_number']))
        alleles_list.append(f"{gene}{info['allele']}")

    print(f" Key Genes          │ {', '.join(genes_list)}")
    print(f" Copy Number (CN)   │ {', '.join(cn_list)}")
    print(f" Alleles            │ {', '.join(alleles_list)}")
    print("═" * 80 + "\n")


# ────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Simple KIR typing from one FASTQ file")
    parser.add_argument("fastq", help="Path to your FASTQ file")
    parser.add_argument("--sample", default="SAMPLE", help="Sample name/ID (optional)")
    parser.add_argument("--min-reads", type=int, default=5, help="Minimum reads to call a gene")

    args = parser.parse_args()

    fastq_path = args.fastq
    sample_name = args.sample

    if not os.path.isfile(fastq_path):
        print(f"Error: File not found → {fastq_path}")
        sys.exit(1)

    print("\n" + "═" * 60)
    print("  KIR TYPING PIPELINE - Single FASTQ File")
    print("═" * 60)
    print(f"  Sample     : {sample_name}")
    print(f"  FASTQ file : {fastq_path}")
    print("═" * 60 + "\n")

    # STEP 1: Input FASTQ
    sequences = read_fastq(fastq_path)

    # STEP 2: Pretend alignment happened
    print("  STEP 2 : FASTQ → BAM (simulated)")
    print("           In real pipeline: BWA-MEM → samtools sort → index")
    print("  → BAM ready (simulated)\n")

    # STEP 3–5: Graph-KIR simulation + gene calling
    gene_hits = simple_kir_mapping(sequences)
    detected_genes = call_genes(gene_hits, min_reads=args.min_reads)

    # STEP 6: Haplotype
    haplotype = predict_haplotype(detected_genes)

    # STEP 7: Show final table
    print_client_table(sample_name, haplotype, detected_genes)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Usage:")
        print("  python kir_typing_simple.py your_file.fastq.gz")
        print("  python kir_typing_simple.py sample_R1.fastq --sample PATIENT-123")
        sys.exit(0)

    main()