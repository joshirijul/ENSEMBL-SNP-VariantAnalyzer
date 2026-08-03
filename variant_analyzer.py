#Performs in-silico translation and calculates biophysical property shifts.

from typing import Dict, Any

# Biophysical properties mapping (Values: Kyte-Doolittle, MW (Da), pI charge roughly at pH 7.4, Volume (Å³))
AA_PROPS = {
    'A': {'hydro': 1.8, 'mw': 89.1, 'charge': 0, 'vol': 67},
    'R': {'hydro': -4.5, 'mw': 174.2, 'charge': 1, 'vol': 148},
    'N': {'hydro': -3.5, 'mw': 132.1, 'charge': 0, 'vol': 96},
    'D': {'hydro': -3.5, 'mw': 133.1, 'charge': -1, 'vol': 91},
    'C': {'hydro': 2.5, 'mw': 121.2, 'charge': 0, 'vol': 86},
    'Q': {'hydro': -3.5, 'mw': 146.1, 'charge': 0, 'vol': 114},
    'E': {'hydro': -3.5, 'mw': 147.1, 'charge': -1, 'vol': 109},
    'G': {'hydro': -0.4, 'mw': 75.1, 'charge': 0, 'vol': 48},
    'H': {'hydro': -3.2, 'mw': 155.2, 'charge': 0.1, 'vol': 118},
    'I': {'hydro': 4.5, 'mw': 131.2, 'charge': 0, 'vol': 124},
    'L': {'hydro': 3.8, 'mw': 131.2, 'charge': 0, 'vol': 124},
    'K': {'hydro': -3.9, 'mw': 146.2, 'charge': 1, 'vol': 135},
    'M': {'hydro': 1.9, 'mw': 149.2, 'charge': 0, 'vol': 124},
    'F': {'hydro': 2.8, 'mw': 165.2, 'charge': 0, 'vol': 135},
    'P': {'hydro': -1.6, 'mw': 115.1, 'charge': 0, 'vol': 90},
    'S': {'hydro': -0.8, 'mw': 105.1, 'charge': 0, 'vol': 73},
    'T': {'hydro': -0.7, 'mw': 119.1, 'charge': 0, 'vol': 93},
    'W': {'hydro': -0.9, 'mw': 204.2, 'charge': 0, 'vol': 163},
    'Y': {'hydro': -1.3, 'mw': 181.2, 'charge': 0, 'vol': 141},
    'V': {'hydro': 4.2, 'mw': 117.1, 'charge': 0, 'vol': 105},
    '*': {'hydro': 0.0, 'mw': 0.0, 'charge': 0, 'vol': 0} # Stop codon
}

class VariantAnalyzer:
    def __init__(self, client):
        self.client = client

    def analyze_snp(self, rsid: str) -> Dict[str, Any]:
        """Main orchestration method to analyze the SNP and return a structured report."""
        vep_data = self.client.get_vep_by_rsid(rsid)[0]
        var_data = self.client.get_variation_details(rsid)
        
        # Isolate canonical transcript or fallback to the first coding transcript
        transcripts = vep_data.get('transcript_consequences', [])
        target_tx = next((tx for tx in transcripts if tx.get('canonical') == 1 and 'protein_start' in tx), None)
        if not target_tx:
            target_tx = next((tx for tx in transcripts if 'protein_start' in tx), transcripts[0])

        result = {
            'rsid': rsid,
            'gene': target_tx.get('gene_symbol', 'Unknown'),
            'transcript_id': target_tx.get('transcript_id'),
            'protein_id': target_tx.get('protein_id'),
            'consequence': ", ".join(target_tx.get('consequence_terms', [])),
            'clinical_significance': ", ".join(var_data.get('clinical_significance', ['Unknown']))
        }

        if 'missense_variant' in result['consequence'] or 'nonsense_variant' in result['consequence']:
            result.update(self._process_coding_variant(target_tx))

        return result

    def _process_coding_variant(self, tx: Dict) -> Dict[str, Any]:
        """Processes the exact sequence and biophysical changes for a coding variant."""
        aa_change = tx.get('amino_acids', 'X/X').split('/')
        codon_change = tx.get('codons', 'xxx/xxx').split('/')
        pos = tx.get('protein_start', 1)
        
        wt_aa = aa_change[0]
        mut_aa = aa_change[1] if len(aa_change) > 1 else aa_change[0]
        
        # Fetch Sequences using transcript_id and explicit protein translation
        transcript_id = tx.get('transcript_id')
        wt_prot_seq = self.client.get_sequence(transcript_id, seq_type="protein")
        
        # Construct Alignment
        start_idx = max(0, pos - 11)
        end_idx = min(len(wt_prot_seq), pos + 10)
        
        wt_snippet = wt_prot_seq[start_idx:end_idx]
        mut_snippet = wt_snippet[: (pos - 1 - start_idx)] + mut_aa + wt_snippet[(pos - start_idx):]
        
        # Calculate Shifts
        wt_props = AA_PROPS.get(wt_aa, AA_PROPS['A'])
        mut_props = AA_PROPS.get(mut_aa, AA_PROPS['A'])
        
        return {
            'position': pos,
            'wt_codon': codon_change[0],
            'mut_codon': codon_change[1] if len(codon_change) > 1 else codon_change[0],
            'wt_aa': wt_aa,
            'mut_aa': mut_aa,
            'alignment': {
                'start_pos': start_idx + 1,
                'wt': wt_snippet,
                'mut': mut_snippet,
                'indicator': ' ' * (pos - 1 - start_idx) + '^'
            },
            'biophysics': {
                'delta_hydropathy': round(mut_props['hydro'] - wt_props['hydro'], 2),
                'delta_mw_da': round(mut_props['mw'] - wt_props['mw'], 2),
                'delta_charge': round(mut_props['charge'] - wt_props['charge'], 2),
                'delta_volume': round(mut_props['vol'] - wt_props['vol'], 2)
            }
        }