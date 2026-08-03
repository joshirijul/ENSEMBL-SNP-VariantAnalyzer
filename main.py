#Execution script for the Ensembl-SNP-VariantAnalyzer.

import argparse
from ensembl_client import EnsemblClient
from variant_analyzer import VariantAnalyzer
from report_generator import ReportGenerator

def process_variant(rsid: str, analyzer: VariantAnalyzer, reporter: ReportGenerator) -> bool:
    """Helper function to process a single variant and handle local errors."""
    try:
        result_data = analyzer.analyze_snp(rsid)
        # Skip if no coding sequence changes were found (e.g., intronic variants)
        if 'wt_aa' not in result_data:
            print(f"[-] {rsid} is non-coding or lacks amino acid shifts. Skipping report.")
            return False
            
        reporter.print_terminal(result_data)
        reporter.export(result_data)
        return True
    except Exception as e:
        print(f"[!] Failed to process {rsid}: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Ensembl SNP Variant Analyzer")
    
    # Create a mutually exclusive group: User must provide EITHER rsid OR gene
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--rsid", type=str, help="dbSNP rsID to analyze directly (e.g., rs63750847)")
    group.add_argument("--gene", type=str, help="Gene symbol for pathogenic variant discovery (e.g., XPA)")
    
    # Optional limit for batch processing
    parser.add_argument("--limit", type=int, default=3, help="Max variants to process in gene mode (default: 3)")
    
    args = parser.parse_args()

    # Initialize modules
    client = EnsemblClient()
    analyzer = VariantAnalyzer(client)
    reporter = ReportGenerator()

    # Route 1: Single rsID execution
    if args.rsid:
        print(f"\n[*] Initializing direct lookup for {args.rsid}...")
        process_variant(args.rsid, analyzer, reporter)
        print("\n[+] Pipeline execution complete. Files saved in /results.")
        
    # Route 2: Automated Gene Discovery execution
    elif args.gene:
        print(f"\n[*] Initializing discovery mode for gene {args.gene}...")
        rsids = client.discover_pathogenic_snps(args.gene)
        
        if not rsids:
            print(f"[-] No pathogenic variants found for {args.gene}.")
            return
            
        limit = min(args.limit, len(rsids))
        print(f"[*] Processing the first {limit} variants...\n")
        
        successful = 0
        for rsid in rsids[:limit]:
            if process_variant(rsid, analyzer, reporter):
                successful += 1
                
        print(f"\n[+] Discovery execution complete. Successfully mapped {successful}/{limit} variants. Files saved in /results.")

if __name__ == "__main__":
    main()