#Handles all interactions with the Ensembl REST API.

import time
import requests
from typing import Dict, Any, Optional

class EnsemblClient:
    BASE_URL = "https://rest.ensembl.org"

    def __init__(self, species: str = "human"):
        self.species = species
        self.headers = {"Content-Type": "application/json"}

    def _request(self, endpoint: str) -> Optional[Any]:
        """Internal method to handle requests and HTTP 429 rate limiting."""
        url = f"{self.BASE_URL}{endpoint}"
        max_retries = 3
        
        for attempt in range(max_retries):
            response = requests.get(url, headers=self.headers)
            
            if response.status_code == 200:
                return response.json()
            elif response.status_code == 429:
                # Ensembl specific rate limit header
                retry_after = int(response.headers.get("Retry-After", 1))
                print(f"[Warning] Rate limited by Ensembl. Retrying in {retry_after}s...")
                time.sleep(retry_after)
            else:
                response.raise_for_status()
                
        raise Exception(f"Failed to fetch {url} after {max_retries} retries.")

    def get_vep_by_rsid(self, rsid: str) -> list:
        """Fetch Variant Effect Predictor (VEP) data to map variant to transcripts."""
        # Added ?canonical=1 to ensure the 'canonical' flag is returned in the JSON
        return self._request(f"/vep/{self.species}/id/{rsid}?canonical=1")

    def get_variation_details(self, rsid: str) -> Dict[str, Any]:
        """Fetch ancestral allele, synonyms, and clinical significance."""
        return self._request(f"/variation/{self.species}/{rsid}")

    def get_sequence(self, ensembl_id: str, seq_type: str = "protein") -> str:
        """Fetch canonical transcript (cDNA) or protein sequence by ID."""
        # Added ?type=protein so Ensembl automatically translates transcripts
        data = self._request(f"/sequence/id/{ensembl_id}?type={seq_type}")
        return data.get("seq", "") if data else ""

    def discover_pathogenic_snps(self, gene_symbol: str) -> list:
        """Finds all clinically significant rsIDs for a given gene symbol."""
        print(f"[*] Discovering variants for gene: {gene_symbol}...")
        
        # Step 1: Look up the Ensembl Gene ID & Genomic Coordinates
        gene_data = self._request(f"/lookup/symbol/{self.species}/{gene_symbol}")
        if not gene_data or 'seq_region_name' not in gene_data:
            print(f"[!] Could not find Ensembl coordinates for gene {gene_symbol}")
            return []
            
        chrom = gene_data.get("seq_region_name")
        start = gene_data.get("start")
        end = gene_data.get("end")
        
        print(f"[*] Mapped {gene_symbol} to chr{chrom}:{start}-{end}. Fetching regional variants (this may take a moment)...")
        
        # Step 2: Query all variants in this physical genomic region
        # Note: feature=variation ensures we get the raw variant data with ClinVar tags
        variants = self._request(f"/overlap/region/{self.species}/{chrom}:{start}-{end}?feature=variation")
        
        pathogenic_rsids = set()
        
        # Step 3: Filter for mutations explicitly flagged as pathogenic
        if variants:
            for v in variants:
                clin_sigs = v.get('clinical_significance', [])
                
                # Ensure we are dealing with a list before iterating
                if isinstance(clin_sigs, list):
                    # Check if 'pathogenic' is in any of the significance tags
                    if any('pathogenic' in str(cs).lower() for cs in clin_sigs):
                        rsid = v.get('id', '')
                        if rsid.startswith('rs'):
                            pathogenic_rsids.add(rsid)
                
        # Convert set to a sorted list for consistent processing
        final_list = sorted(list(pathogenic_rsids))
        print(f"[+] Found {len(final_list)} pathogenic variants for {gene_symbol}.")
        return final_list