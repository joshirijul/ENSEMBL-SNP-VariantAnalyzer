# Ensembl SNP VariantAnalyzer 🧬

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![Ensembl API](https://img.shields.io/badge/Data-Ensembl_REST-green.svg)
![License](https://img.shields.io/badge/License-MIT-purple.svg)

A pure-Python pipeline for functional genomics and high-throughput variant effect analysis. `Ensembl-SNP-VariantAnalyzer` interfaces directly with the live Ensembl REST API to evaluate the biological and biophysical consequences of Single Nucleotide Polymorphisms (SNPs) on canonical mRNA transcripts and primary protein structures.

The tool features an **Automated Discovery Mode** that dynamically mines the Ensembl database for pathogenic missense and nonsense mutations across entire genes, automatically filtering by Sequence Ontology (SO) terms and extracting structural damage profiles.

## ✨ Core Features
*   **Automated Target Discovery:** Input a gene symbol (e.g., `BRCA1`, `XPA`), and the pipeline will map its genomic coordinates, extract all overlapping regional variants, and isolate pathogenic targets automatically.
*   **Live Ensembl Integration:** Automated fetching of Variant Effect Predictor (VEP) data, canonical transcripts, and wild-type cDNA/protein sequences.
*   **Robust API Client:** Built-in network error handling, string-safe type casting, and HTTP 429 rate-limit backoff.
*   **Biophysical Shift Profiling:** Calculates and quantifies changes in:
    *   Kyte-Doolittle Hydropathy Index
    *   Molecular Weight (Da)
    *   Isoelectric Charge (at pH 7.4)
    *   Steric Volume (Å³)
*   **Academic Export:** Generates terminal summaries and exports JSON/Markdown reports with local sequence alignments.

## 🔬 Downstream Applications
This pipeline serves as a critical pre-processing and feature-extraction engine for:
*   **Deep Learning Feature Generation:** Automating the extraction of physical property deltas (e.g., Hydropathy, Volume) to feed numerical datasets into custom frameworks like 3D U-Nets for oncology classification.
*   **Pre-Docking Target Validation:** Identifying steric clashes or charge flips in therapeutic targets (e.g., BACE1, APP) before running molecular docking simulations.
*   **High-Throughput Triage:** Rapidly triaging actionable genomic variants from sequencing panels.

## 🚀 Installation

Clone the repository and install the required dependencies:

```bash
git clone [https://github.com/yourusername/Ensembl-SNP-VariantAnalyzer.git](https://github.com/yourusername/Ensembl-SNP-VariantAnalyzer.git)
cd Ensembl-SNP-VariantAnalyzer
pip install -r requirements.txt