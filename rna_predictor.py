#!/usr/bin/env python3
"""Basic RNA secondary-structure predictor using Nussinov dynamic programming."""
from __future__ import annotations
import argparse
import json
import re
from typing import List, Tuple

CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
WOBBLE = {("G", "U"), ("U", "G")}


def validate_sequence(sequence: str) -> str:
    seq = re.sub(r"\s+", "", sequence).upper().replace("T", "U")
    if not seq:
        raise ValueError("RNA sequence cannot be empty.")
    invalid = sorted(set(seq) - set("ACGU"))
    if invalid:
        raise ValueError(f"Invalid RNA base(s): {', '.join(invalid)}. Use A, C, G, U.")
    return seq


def can_pair(a: str, b: str, allow_wobble: bool = True) -> bool:
    return (a, b) in CANONICAL or (allow_wobble and (a, b) in WOBBLE)


def predict_structure(sequence: str, min_loop: int = 3, allow_wobble: bool = True) -> dict:
    seq = validate_sequence(sequence)
    if min_loop < 0:
        raise ValueError("Minimum loop length cannot be negative.")
    n = len(seq)
    dp = [[0] * n for _ in range(n)]

    for span in range(1, n):
        for i in range(n - span):
            j = i + span
            best = max(dp[i + 1][j], dp[i][j - 1])
            if j - i > min_loop and can_pair(seq[i], seq[j], allow_wobble):
                best = max(best, dp[i + 1][j - 1] + 1)
            for k in range(i + 1, j):
                best = max(best, dp[i][k] + dp[k + 1][j])
            dp[i][j] = best

    pairs: List[Tuple[int, int]] = []

    def traceback(i: int, j: int) -> None:
        if i >= j:
            return
        if dp[i][j] == dp[i + 1][j]:
            traceback(i + 1, j)
        elif dp[i][j] == dp[i][j - 1]:
            traceback(i, j - 1)
        elif (j - i > min_loop and can_pair(seq[i], seq[j], allow_wobble)
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
    for i, j in pairs:
        structure[i - 1] = "("
        structure[j - 1] = ")"

    return {
        "sequence": seq,
        "length": n,
        "structure": "".join(structure),
        "base_pair_count": len(pairs),
        "base_pairs": [
            {"position_1": i, "position_2": j, "bases": f"{seq[i-1]}-{seq[j-1]}",
             "type": "canonical" if (seq[i-1], seq[j-1]) in CANONICAL else "wobble"}
            for i, j in pairs
        ],
        "minimum_loop_length": min_loop,
        "wobble_pairs_enabled": allow_wobble,
        "algorithm": "Nussinov dynamic programming",
        "complexity": {"time": "O(n^3)", "space": "O(n^2)"},
        "warning": "Educational basic predictor; it maximizes base-pair count and does not calculate thermodynamic free energy."
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict RNA secondary structure in dot-bracket notation.")
    parser.add_argument("sequence", nargs="?", help="RNA sequence using A/C/G/U")
    parser.add_argument("--min-loop", type=int, default=3)
    parser.add_argument("--no-wobble", action="store_true", help="Disable G-U wobble pairs")
    parser.add_argument("--json", action="store_true", help="Output machine-readable JSON")
    args = parser.parse_args()
    sequence = args.sequence or input("Enter RNA sequence (A, C, G, U): ")
    try:
        result = predict_structure(sequence, args.min_loop, not args.no_wobble)
    except ValueError as exc:
        parser.error(str(exc))
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print("\n=== RNA Secondary Structure Predictor ===")
        print(f"Sequence : {result['sequence']}")
        print(f"Length   : {result['length']}")
        print(f"Structure: {result['structure']}")
        print(f"Pairs    : {result['base_pair_count']}")
        print(f"Pairs    : {', '.join(f'{p['position_1']}-{p['position_2']}' for p in result['base_pairs']) or 'None'}")
        print(f"Algorithm: {result['algorithm']}")


if __name__ == "__main__":
    main()
