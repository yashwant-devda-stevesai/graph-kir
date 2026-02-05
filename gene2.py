#!/usr/bin/env python3
"""
Streamlit UI for Simple KIR Gene Typing Pipeline
Shows result in separate lines per gene:
id 1 - gene name - copy no - allele - type ...
"""

import streamlit as st
import io
from collections import defaultdict
import numpy as np
from Bio import SeqIO
import gzip


# ────────────────────────────────────────────────
# Page configuration
# ────────────────────────────────────────────────
st.set_page_config(
    page_title="KIR Gene Typing",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)


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
# Read FASTQ file
# ────────────────────────────────────────────────
@st.cache_data(show_spinner="Reading FASTQ file...")
def read_fastq(file_bytes, filename):
    sequences = []
    try:
        if filename.lower().endswith('.gz'):
            f = gzip.open(io.BytesIO(file_bytes), 'rt')
        else:
            f = io.StringIO(file_bytes.decode('utf-8'))

        for record in SeqIO.parse(f, "fastq"):
            sequences.append(str(record.seq))

        return sequences

    except Exception as e:
        st.error(f"Error reading FASTQ: {e}")
        return None


# ────────────────────────────────────────────────
# Simple mapping simulation
# ────────────────────────────────────────────────
def simple_kir_mapping(sequences):
    gene_hits = defaultdict(int)
    progress_bar = st.progress(0)
    status = st.empty()

    status.text("Mapping reads to KIR genes (simulation)...")

    total = len(sequences)
    for i, seq in enumerate(sequences, 1):
        seq_short = seq[:60]
        for gene in KIR_GENES:
            # naive matching – demo only
            if gene.lower() in seq_short.lower() or len(seq_short) > 55:
                gene_hits[gene] += 1

        if i % 20000 == 0 or i == total:
            progress_bar.progress(min(1.0, i / total))

    status.text("Mapping finished.")
    progress_bar.progress(1.0)
    return dict(gene_hits)


# ────────────────────────────────────────────────
# Call genes + copy number + fake allele
# ────────────────────────────────────────────────
def call_genes(gene_hits, min_reads=5):
    detected = {}
    for gene, count in gene_hits.items():
        if count < min_reads:
            continue
        cn = 1 if count < min_reads * 1.8 else 2
        detected[gene] = {
            'read_count': count,
            'copy_number': cn,
            'allele': f"*{cn}0{count % 10}"   # fake allele
        }
    return detected


# ────────────────────────────────────────────────
# Simple haplotype prediction
# ────────────────────────────────────────────────
def predict_haplotype(detected_genes):
    present = set(detected_genes.keys())
    a_count = len(present & HAPLOTYPE_A_GENES)
    b_count = len(present & HAPLOTYPE_B_GENES)

    if b_count > a_count + 1:
        return 'B'
    elif a_count > b_count:
        return 'A'
    else:
        return 'A/B'


# ────────────────────────────────────────────────
# Display result in the requested format
# ────────────────────────────────────────────────
def display_result(sample_id, haplotype, detected):
    if not detected:
        st.warning("No KIR genes detected above the read threshold.")
        return

    st.markdown("### Final KIR Typing Result – Client Report")
    st.markdown(f"**Sample ID** : {sample_id}")
    st.markdown(f"**Predicted Haplotype** : {haplotype}")
    st.markdown("---")

    st.subheader("Detected KIR Genes")

    # One line per gene
    for idx, (gene, info) in enumerate(sorted(detected.items()), 1):
        allele_full = f"{gene}{info['allele']}"
        gene_type = KIR_GENES[gene]['type']

        st.markdown(
            f"**id {idx}** - **{gene}** - copy no **{info['copy_number']}** - "
            f"allele **{allele_full}** - type **{gene_type}**"
        )

    st.markdown("---")


# ────────────────────────────────────────────────
# Main Streamlit App
# ────────────────────────────────────────────────
def main():
    st.title("🧬 KIR Gene Typing Tool")
    st.markdown("Upload a single FASTQ file to get a simple KIR typing result.")

    # Sidebar settings
    with st.sidebar:
        st.header("Settings")
        sample_id = st.text_input("Sample ID / Name", value="SAMPLE01")
        min_reads = st.number_input(
            "Minimum reads to call a gene",
            min_value=1,
            max_value=100,
            value=5,
            step=1
        )

        st.markdown("---")
        st.info("This is a **demo** version using very simple read matching.\nReal analysis would use proper alignment and graph-based typing.")

    # Upload area
    uploaded_file = st.file_uploader(
        "Upload FASTQ file (.fastq, .fq, .fastq.gz)",
        type=["fastq", "fq", "fastq.gz"],
        help="Only one file at a time"
    )

    if uploaded_file is not None:
        st.success(f"Uploaded: **{uploaded_file.name}**  ({uploaded_file.size / 1024 / 1024:.2f} MB)")

        if st.button("🚀 Start Analysis", type="primary", use_container_width=True):

            with st.spinner("Reading FASTQ file..."):
                file_bytes = uploaded_file.read()
                sequences = read_fastq(file_bytes, uploaded_file.name)

            if sequences is None:
                st.stop()

            st.info(f"**{len(sequences):,} reads** loaded from file.")

            # Simulated alignment step
            with st.status("Step 2: FASTQ → BAM (simulated)", expanded=False):
                st.write("In real pipeline: BWA-MEM → samtools sort → index")
                st.write("→ Simulated BAM ready")

            # Mapping & typing
            st.subheader("Processing...")
            gene_hits = simple_kir_mapping(sequences)
            detected_genes = call_genes(gene_hits, min_reads=min_reads)

            # Haplotype
            haplotype = predict_haplotype(detected_genes)

            # Show result in your requested format
            display_result(sample_id, haplotype, detected_genes)

            if detected_genes:
                st.balloons()


if __name__ == "__main__":
    main()