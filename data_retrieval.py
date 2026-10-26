from Bio import SeqIO
import os

# ============================================================
# Step 2: Data Retrieval (Offline Mode)
# ============================================================

# We are now reading directly from the FASTA file you created.
WT_FASTA_FILE = "cyp2d6_wt.fasta"

def fetch_wt_sequence():
    """Load the WT sequence directly from the local FASTA file."""
    if not os.path.exists(WT_FASTA_FILE):
        print(f"[ERROR] Could not find {WT_FASTA_FILE}. Please make sure it exists.")
        return None
        
    print(f"[INFO] Loading WT sequence from {WT_FASTA_FILE}")
    record = SeqIO.read(WT_FASTA_FILE, "fasta")
    
    print(f"[INFO] WT sequence: {record.id}, length = {len(record.seq)} aa")
    return record

# --- 2.2 Hardcoded clinical nsSNP dataset ---
# Format: (allele_name, position_1based, wt_aa, mutant_aa, rs_id)
CYP2D6_MUTATIONS = [
    ("CYP2D6*7",  324, "H", "P", "rs5030867"),   # His324Pro
    ("CYP2D6*8",  296, "R", "C", "rs16947"),      # Cys296Arg
    ("CYP2D6*10",  34, "P", "S", "rs1065852"),    # Pro34Ser
    ("CYP2D6*17", 107, "T", "I", "rs28371706"),   # Thr107Ile
    ("CYP2D6*12",  42, "G", "R", "rs5030862"),    # Gly42Arg
    ("CYP2D6*49", 120, "F", "I", "rs1135822"),    # Phe120Ile
    ("CYP2D6*56", 344, "R", "X", "rs72549347"),   # Arg344Ter (nonsense)
]

if __name__ == "__main__":
    print("\n[MUTATION DATASET]")
    for m in CYP2D6_MUTATIONS:
        print(f"  {m[0]}: {m[2]}{m[1]}{m[3]}  ({m[4]})")

    # Run the fetch function
    wt_record = fetch_wt_sequence()
    
    if wt_record:
        print("\n[SUCCESS] Step 2 complete. Ready for Step 3.")