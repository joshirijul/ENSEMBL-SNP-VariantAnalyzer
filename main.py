#Execution script for the Ensembl-SNP-VariantAnalyzer.

import argparse
from ensembl_client import EnsemblClient
from variant_analyzer import VariantAnalyzer
from report_generator import ReportGenerator

def main():
    parser = argparse.ArgumentParser(description="Ensembl SNP Variant Analyzer")
    parser.add_argument("--rsid", type=str, default="rs63750847", help="dbSNP rsID to analyze")
    args = parser.parse_args()

    print(f"[*] Initializing Ensembl REST pipeline for {args.rsid}...")
    
    # Initialize modules
    client = EnsemblClient()
    analyzer = VariantAnalyzer(client)
    reporter = ReportGenerator()

    # Execute pipeline
    try:
        print("[*] Fetching genetic, transcript, and VEP data...")
        result_data = analyzer.analyze_snp(args.rsid)
        
        print("[*] Generating reports and exporting files...")
        reporter.print_terminal(result_data)
        reporter.export(result_data)
        
        print("[+] Pipeline execution complete. Files saved in /results.")
        
    except Exception as e:
        print(f"[!] Pipeline failed: {str(e)}")

if __name__ == "__main__":
    main()