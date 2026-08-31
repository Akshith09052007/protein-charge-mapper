"""
Protein Sequence Charge Distribution Mapper — Backend Core Logic

Responsibility: Backend core logic development, including FASTA sequence 
parsing and amino acid charge calculation algorithms.
"""

from collections import Counter

# Standard 20 single-letter amino acid codes
VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# Standard pKa values for titratable groups (Lehninger scale)
PKA_VALUES = {
    "N_TERM": 8.6,
    "C_TERM": 3.6,
    "K": 10.5,  # Lysine (+)
    "R": 12.5,  # Arginine (+)
    "H": 6.0,   # Histidine (+)
    "D": 3.9,   # Aspartate (-)
    "E": 4.1,   # Glutamate (-)
    "C": 8.3,   # Cysteine (-)
    "Y": 10.9,  # Tyrosine (-)
}


def clean_and_validate_protein(raw_sequence: str) -> str:
    """Strips whitespace, forces uppercase, and validates 20 standard amino acids."""
    sequence = raw_sequence.strip().upper().replace("\n", "").replace(" ", "")
    
    invalid_chars = set(sequence) - VALID_AMINO_ACIDS
    if invalid_chars:
        raise ValueError(
            f"Sequence contains invalid amino acid character(s): {sorted(invalid_chars)}. "
            f"Only standard 20 amino acids are allowed."
        )
    
    if len(sequence) == 0:
        raise ValueError("Sequence cannot be empty.")
        
    return sequence


def compute_residue_charge(residue: str, ph: float) -> float:
    """Calculates partial charge of a side chain at target pH (Henderson-Hasselbalch)."""
    if residue in ("K", "R", "H"):
        pka = PKA_VALUES[residue]
        return 1.0 / (1.0 + 10 ** (ph - pka))
    elif residue in ("D", "E", "C", "Y"):
        pka = PKA_VALUES[residue]
        return -1.0 / (1.0 + 10 ** (pka - ph))
    return 0.0


def compute_total_net_charge(sequence: str, ph: float = 7.0) -> float:
    """Computes total net charge of a sequence including N/C termini at given pH."""
    n_term_charge = 1.0 / (1.0 + 10 ** (ph - PKA_VALUES["N_TERM"]))
    c_term_charge = -1.0 / (1.0 + 10 ** (PKA_VALUES["C_TERM"] - ph))
    
    sidechain_charge = sum(compute_residue_charge(aa, ph) for aa in sequence)
    return n_term_charge + c_term_charge + sidechain_charge


def calculate_sliding_window_charge(
    sequence: str, ph: float = 7.0, window_size: int = 5
) -> list[float]:
    """Maps local charge density across the protein using a sliding window."""
    if window_size % 2 == 0:
        raise ValueError("Window size must be an odd integer.")
        
    half_win = window_size // 2
    per_residue_charges = [compute_residue_charge(aa, ph) for aa in sequence]
    window_profile = []
    
    for i in range(len(sequence)):
        start = max(0, i - half_win)
        end = min(len(sequence), i + half_win + 1)
        window_charge = sum(per_residue_charges[start:end])
        window_profile.append(round(window_charge, 3))
        
    return window_profile


def analyze_protein_sequence(
    raw_sequence: str, ph: float = 7.0, window_size: int = 5
) -> dict:
    """Primary entry point returning complete charge analytics dictionary."""
    sequence = clean_and_validate_protein(raw_sequence)
    total_charge = compute_total_net_charge(sequence, ph)
    window_profile = calculate_sliding_window_charge(sequence, ph, window_size)
    
    return {
        "sequence_length": len(sequence),
        "ph": ph,
        "window_size": window_size,
        "total_net_charge": round(total_charge, 3),
        "amino_acid_counts": dict(Counter(sequence)),
        "sliding_window_profile": window_profile,
    }


if __name__ == "__main__":
    sample_protein = "GIVEQCCTSICSLYQLENYCN"
    res = analyze_protein_sequence(sample_protein, ph=7.0, window_size=5)
    print("Net Charge:", res["total_net_charge"])
    print("Charge Profile:", res["sliding_window_profile"])
