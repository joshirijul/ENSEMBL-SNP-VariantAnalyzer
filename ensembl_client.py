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