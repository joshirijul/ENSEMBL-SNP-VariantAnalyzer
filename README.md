# Ensembl SNP VariantAnalyzer 🧬

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![Ensembl API](https://img.shields.io/badge/Data-Ensembl_REST-green.svg)
![License](https://img.shields.io/badge/License-MIT-purple.svg)

A lightweight, pure-Python pipeline for functional genomics and variant effect analysis. `Ensembl-SNP-VariantAnalyzer` interfaces directly with the live Ensembl REST API to evaluate the biological and biophysical consequences of Single Nucleotide Polymorphisms (SNPs) on canonical mRNA transcripts and primary protein structures.

This tool is designed to automate the preprocessing of variant data for downstream structural biology workflows, such as identifying steric shifts or charge alterations prior to molecular docking simulations.

## ✨ Core Features
*   **Live Ensembl Integration:** Automated fetching of Variant Effect Predictor (VEP) data, genomic coordinates, and wild-type cDNA/protein sequences.
*   **Robust API Client:** Built-in network error handling and HTTP 429 rate-limit backoff.
*   **In-Silico Translation:** Accurately reconstructs mutated coding sequences and isolates the exact codon/amino acid shift.
*   **Biophysical Shift Profiling:** Calculates changes in:
    *   Kyte-Doolittle Hydropathy Index
    *   Molecular Weight (Da)
    *   Isoelectric Charge (at pH 7.4)
    *   Steric Volume (Å³)
*   **Academic Export:** Generates terminal summaries and exports JSON/Markdown reports with local sequence alignments.

## 🚀 Installation

Clone the repository and install the required dependencies:

```bash
git clone [https://github.com/yourusername/Ensembl-SNP-VariantAnalyzer.git](https://github.com/yourusername/Ensembl-SNP-VariantAnalyzer.git)
cd Ensembl-SNP-VariantAnalyzer
pip install -r requirements.txt