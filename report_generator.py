#Generates terminal outputs and exports results to disk.
import os
import json

class ReportGenerator:
    def __init__(self, output_dir: str = "results"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def print_terminal(self, data: dict):
        """Prints a highly formatted biological summary to the terminal."""
        print(f"\n{'='*55}")
        print(f" VARIANT ANALYSIS REPORT: {data['rsid']} ({data['gene']})")
        print(f"{'='*55}")
        print(f"Transcript ID : {data.get('transcript_id')}")
        print(f"Consequence   : {data.get('consequence')}")
        print(f"Clinical Sig  : {data.get('clinical_significance').capitalize()}")
        
        if 'wt_aa' in data:
            print(f"\n[MOLECULAR CHANGES]")
            print(f"Codon         : {data['wt_codon']} -> {data['mut_codon']}")
            print(f"Amino Acid    : {data['wt_aa']} -> {data['mut_aa']} (Pos: {data['position']})")
            
            bp = data['biophysics']
            print(f"\n[BIOPHYSICAL SHIFTS]")
            print(f"Hydropathy (Kyte-Doolittle) : {bp['delta_hydropathy']:+0.2f}")
            print(f"Molecular Weight (Da)       : {bp['delta_mw_da']:+0.2f}")
            print(f"Charge (pH 7.4)             : {bp['delta_charge']:+0.2f}")
            print(f"Steric Volume (Å³)          : {bp['delta_volume']:+0.2f}")
            
            aln = data['alignment']
            print(f"\n[LOCAL ALIGNMENT]")
            print(f"Pos {aln['start_pos']:<4} WT:  {aln['wt']}")
            print(f"         MUT: {aln['mut']}")
            print(f"              {aln['indicator']}")
        print(f"{'='*55}\n")

    def export(self, data: dict):
        """Exports data to required submission formats."""
        # 1. JSON Summary
        with open(f"{self.output_dir}/variant_summary.json", 'w') as f:
            json.dump(data, f, indent=4)
            
        if 'wt_aa' in data:
            aln = data['alignment']
            # 2. Alignment txt
            with open(f"{self.output_dir}/alignment.txt", 'w') as f:
                f.write(f"Alignment for {data['rsid']} in {data['gene']}\n\n")
                f.write(f"WT  : {aln['wt']}\n")
                f.write(f"MUT : {aln['mut']}\n")
                f.write(f"      {aln['indicator']}\n")
                
            # 3. Markdown Report
            with open(f"{self.output_dir}/report.md", 'w') as f:
                f.write(f"# Variant Analysis: {data['rsid']} ({data['gene']})\n\n")
                f.write(f"**Consequence:** {data['consequence']}\n")
                f.write(f"**Clinical Significance:** {data['clinical_significance']}\n\n")
                f.write("## Sequence Details\n")
                f.write(f"- **Codon Shift:** {data['wt_codon']} -> {data['mut_codon']}\n")
                f.write(f"- **Amino Acid Shift:** {data['wt_aa']} -> {data['mut_aa']} at position {data['position']}\n\n")
                f.write("## Biophysical Property Shifts\n")
                f.write("| Property | Delta |\n|---|---|\n")
                f.write(f"| Hydropathy | {data['biophysics']['delta_hydropathy']:+.2f} |\n")
                f.write(f"| Mol. Weight (Da) | {data['biophysics']['delta_mw_da']:+.2f} |\n")
                f.write(f"| Volume (Å³) | {data['biophysics']['delta_volume']:+.2f} |\n")