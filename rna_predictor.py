#!/usr/bin/env python3
"""Basic RNA secondary-structure predictor using Nussinov dynamic programming.

This is an educational predictor. It maximizes the number of canonical base pairs
(A-U, U-A, G-C, C-G) and allows G-U wobble pairs. It does not model thermodynamics.
"""

from typing import List, Tuple

MIN_LOOP = 3
PAIR_SCORES = {
    ("A", "U"): 1,
    ("U", "A"): 1,
    ("G", "C"): 1,
    ("C", "G"): 1,
    ("G", "U"): 1,
    ("U", "G"): 1,
}


def can_pair(a: str, b: str) -> bool:
    return (a, b) in PAIR_SCORES


def nussinov(sequence: str) -> Tuple[str, List[Tuple[int, int]]]:
    """Predict a structure in dot-bracket notation."""
    seq = sequence.upper().replace(" ", "").replace("\n", "")
    if not seq:
        raise ValueError("RNA sequence cannot be empty.")
    invalid = set(seq) - set("ACGU")
    if invalid:
        raise ValueError(f"Invalid RNA base(s): {', '.join(sorted(invalid))}")

    n = len(seq)
    dp = [[0] * n for _ in range(n)]

    for span in range(1, n):
        for i in range(n - span):
            j = i + span
            best = max(dp[i + 1][j], dp[i][j - 1])

            if j - i > MIN_LOOP and can_pair(seq[i], seq[j]):
                best = max(best, dp[i + 1][j - 1] + PAIR_SCORES[(seq[i], seq[j])])

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
        elif j - i > MIN_LOOP and can_pair(seq[i], seq[j]) and dp[i][j] == dp[i + 1][j - 1] + 1:
            pairs.append((i, j))
            traceback(i + 1, j - 1)
        else:
            for k in range(i + 1, j):
                if dp[i][j] == dp[i][k] + dp[k + 1][j]:
                    traceback(i, k)
                    traceback(k + 1, j)
                    return

    traceback(0, n - 1)
    pairs.sort()

    structure = ["."] * n
    for i, j in pairs:
        structure[i] = "("
        structure[j] = ")"
    return "".join(structure), pairs


def main() -> None:
    print("Basic RNA Secondary Structure Predictor")
    print("Authors: R. Akshith Narasimha (SE25UBIT048); D. Anjana Reddy (SE25UBIT048)")
    sequence = input("Enter RNA sequence (A, C, G, U): ").strip()
    try:
        structure, pairs = nussinov(sequence)
    except ValueError as exc:
        print(f"Error: {exc}")
        return

    print("\nSequence :", sequence.upper().replace(" ", ""))
    print("Structure:", structure)
    print("Base pairs:", len(pairs))
    if pairs:
        print("Paired positions (1-based):", ", ".join(f"{i + 1}-{j + 1}" for i, j in pairs))
    else:
        print("No compatible base pairs found.")
    print("\nNote: This is a basic educational predictor based on the Nussinov algorithm.")


if __name__ == "__main__":
    main()
