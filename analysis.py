from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
from data_retrieval import fetch_wt_sequence, CYP2D6_MUTATIONS

def analyze_mutations():
    # 1. Get the Wild-Type sequence
    wt_record = fetch_wt_sequence()
    if not wt_record:
        return
    wt_seq_str = str(wt_record.seq)
    
    # Initialize a list to store the data for Pandas
    results = []
    
    # 2. Analyze the Wild-Type baseline
    wt_analyzer = ProteinAnalysis(wt_seq_str)
    wt_instability = wt_analyzer.instability_index()
    wt_gravy = wt_analyzer.gravy()
    wt_mw = wt_analyzer.molecular_weight()
    wt_pi = wt_analyzer.isoelectric_point()
    
    # 3. Loop through each mutation and calculate properties
    for allele, pos, wt_aa, mut_aa, rs_id in CYP2D6_MUTATIONS:
        # ProtParam uses 0-based indexing, so we subtract 1 from the position
        idx = pos - 1
        
        # Verify the wild-type amino acid matches
        if wt_seq_str[idx] != wt_aa:
            print(f"[WARNING] {allele}: Expected {wt_aa} at pos {pos}, but found {wt_seq_str[idx]}")

        # Create the mutant sequence
        # Handle the stop codon (X) by truncating the sequence
        if mut_aa == "X":
            mutant_seq_str = wt_seq_str[:idx] 
        else:
            mutant_seq_str = wt_seq_str[:idx] + mut_aa + wt_seq_str[idx+1:]
            
        # Analyze the mutant
        mut_analyzer = ProteinAnalysis(mutant_seq_str)
        mut_instability = mut_analyzer.instability_index()
        mut_gravy = mut_analyzer.gravy()
        mut_mw = mut_analyzer.molecular_weight()
        mut_pi = mut_analyzer.isoelectric_point()
        
        # Calculate the DELTA (Change) between Mutant and Wild-Type
        # Note: A positive change in instability means the protein is MORE unstable
        delta_instability = mut_instability - wt_instability
        delta_gravy = mut_gravy - wt_gravy
        delta_mw = mut_mw - wt_mw
        
        effect = "Destabilizing" if delta_instability > 0 else "Stabilizing"
        
        # Save the results
        results.append({
            "Allele": allele,
            "Mutation": f"{wt_aa}{pos}{mut_aa}",
            "rsID": rs_id,
            "WT_Instability": round(wt_instability, 2),
            "Mut_Instability": round(mut_instability, 2),
            "Delta_Instability": round(delta_instability, 2),
            "Effect_on_Stability": effect,
            "WT_GRAVY": round(wt_gravy, 3),
            "Mut_GRAVY": round(mut_gravy, 3),
            "Delta_GRAVY": round(delta_gravy, 3),
            "WT_MW": round(wt_mw, 2),
            "Mut_MW": round(mut_mw, 2),
            "Delta_MW": round(delta_mw, 2),
            "WT_pI": round(wt_pi, 2),
            "Mut_pI": round(mut_pi, 2)
        })

    # 4. Convert to a Pandas DataFrame and save
    df = pd.DataFrame(results)
    df.to_csv("cyp2d6_analysis_results.csv", index=False)
    print("\n[SUCCESS] Analysis complete! Results saved to cyp2d6_analysis_results.csv")
    
    # Print a preview of the data
    print("\n--- PREVIEW OF RESULTS ---")
    print(df[["Allele", "Mutation", "Delta_Instability", "Effect_on_Stability"]])

if __name__ == "__main__":
    analyze_mutations()