#!/usr/bin/env python3
"""Educational RNA secondary-structure prediction with the Nussinov algorithm."""
from __future__ import annotations

import argparse
import json
import re
from typing import Any

CANONICAL_PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
WOBBLE_PAIRS = {("G", "U"), ("U", "G")}


def validate_sequence(sequence: str) -> str:
    """Remove whitespace, normalize DNA T to RNA U, and validate bases."""
    if not isinstance(sequence, str):
        raise ValueError("Sequence must be provided as text.")
    seq = re.sub(r"\s+", "", sequence).upper().replace("T", "U")
    if not seq:
        raise ValueError("RNA sequence cannot be empty.")
    invalid = sorted(set(seq) - set("ACGU"))
    if invalid:
        raise ValueError(
            f"Invalid RNA base(s): {', '.join(invalid)}. Use only A, C, G, and U."
        )
    return seq


def can_pair(base_a: str, base_b: str, allow_wobble: bool = True) -> bool:
    """Return whether two bases can form a canonical or enabled wobble pair."""
    return ((base_a, base_b) in CANONICAL_PAIRS or
            (allow_wobble and (base_a, base_b) in WOBBLE_PAIRS))


def predict_structure(
    sequence: str, min_loop: int = 3, allow_wobble: bool = True
) -> dict[str, Any]:
    """Predict a maximum-base-pair structure and return 1-based pair positions."""
    seq = validate_sequence(sequence)
    if isinstance(min_loop, bool) or not isinstance(min_loop, int) or min_loop < 0:
        raise ValueError("Minimum loop length must be a non-negative integer.")
    if not isinstance(allow_wobble, bool):
        raise ValueError("allow_wobble must be true or false.")

    n = len(seq)
    dp = [[0] * n for _ in range(n)]

    # dp[i][j] is the maximum number of pairs in the interval i..j.
    for span in range(1, n):
        for i in range(n - span):
            j = i + span
            best = max(dp[i + 1][j], dp[i][j - 1])

            # At least min_loop unpaired nucleotides must lie inside a hairpin.
            if j - i > min_loop and can_pair(seq[i], seq[j], allow_wobble):
                best = max(best, dp[i + 1][j - 1] + 1)

            # Split the interval into two independently folded parts.
            for k in range(i + 1, j):
                best = max(best, dp[i][k] + dp[k + 1][j])
            dp[i][j] = best

    pairs: list[tuple[int, int]] = []

    def traceback(i: int, j: int) -> None:
        if i >= j:
            return
        if dp[i][j] == dp[i + 1][j]:
            traceback(i + 1, j)
        elif dp[i][j] == dp[i][j - 1]:
            traceback(i, j - 1)
        elif (j - i > min_loop
              and can_pair(seq[i], seq[j], allow_wobble)
              and dp[i][j] == dp[i + 1][j - 1] + 1):
            pairs.append((i + 1, j + 1))
            traceback(i + 1, j - 1)
        else:
            for k in range(i + 1, j):
                if dp[i][j] == dp[i][k] + dp[k + 1][j]:
                    traceback(i, k)
                    traceback(k + 1, j)
                    return

    if n > 1:
        traceback(0, n - 1)

    pairs.sort()
    structure = ["."] * n
    base_pairs = []
    for left, right in pairs:
        structure[left - 1] = "("
        structure[right - 1] = ")"
        pair_type = (
            "canonical"
            if (seq[left - 1], seq[right - 1]) in CANONICAL_PAIRS
            else "wobble"
        )
        base_pairs.append({
            "position_1": left,
            "position_2": right,
            "bases": f"{seq[left - 1]}-{seq[right - 1]}",
            "type": pair_type,
        })

    return {
        "sequence": seq,
        "length": n,
        "structure": "".join(structure),
        "base_pair_count": len(pairs),
        "base_pairs": base_pairs,
        "minimum_loop_length": min_loop,
        "wobble_pairs_enabled": allow_wobble,
        "algorithm": "Nussinov dynamic programming",
        "complexity": {"time": "O(n^3)", "space": "O(n^2)"},
        "warning": (
            "Educational model: maximizes the number of base pairs. "
            "It does not calculate thermodynamic free energy or pseudoknots."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Predict RNA secondary structure in dot-bracket notation."
    )
    parser.add_argument("sequence", nargs="?", help="RNA sequence (A, C, G, U)")
    parser.add_argument("--min-loop", type=int, default=3)
    parser.add_argument("--no-wobble", action="store_true")
    parser.add_argument("--json", action="store_true", help="Print JSON output")
    args = parser.parse_args()

    try:
        sequence = args.sequence if args.sequence is not None else input(
            "Enter RNA sequence (A, C, G, U): "
        )
        result = predict_structure(sequence, args.min_loop, not args.no_wobble)
    except ValueError as error:
        parser.error(str(error))

    if args.json:
        print(json.dumps(result, indent=2))
        return

    print("\n=== RNA Secondary Structure Basic Predictor ===")
    print("Authors: R. Akshith Narasimha (SE25UBIT048); "
          "D. Anjana Reddy (SE25UBIT008)")
    print(f"Sequence : {result['sequence']}")
    print(f"Length   : {result['length']}")
    print(f"Structure: {result['structure']}")
    print(f"Pairs    : {result['base_pair_count']}")
    for pair in result["base_pairs"]:
        print(f"  {pair['position_1']}-{pair['position_2']} "
              f"{pair['bases']} ({pair['type']})")
    print(result["warning"])


if __name__ == "__main__":
    main()
