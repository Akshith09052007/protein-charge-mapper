#!/usr/bin/env python3
"""
RNA Secondary Structure Basic Predictor
---------------------------------------
A beginner-friendly, Python-only RNA secondary-structure predictor.

Algorithm: Nussinov dynamic programming
Pairing rules: A-U, U-A, G-C, C-G; optional G-U wobble pairs
Output: dot-bracket notation and 1-based base-pair positions

Run:
    python RNA_Secondary_Structure_Basic_Predictor.py
    python RNA_Secondary_Structure_Basic_Predictor.py GGGAAACCC
    python RNA_Secondary_Structure_Basic_Predictor.py GGGAAACCC --min-loop 3 --no-wobble
    python RNA_Secondary_Structure_Basic_Predictor.py GGGAAACCC --json
"""

import argparse
import json
import sys


CANONICAL_PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
WOBBLE_PAIRS = {("G", "U"), ("U", "G")}


def clean_sequence(sequence):
    """Remove whitespace, convert DNA T to RNA U, and validate the sequence."""
    sequence = "".join(sequence.split()).upper().replace("T", "U")
    if not sequence:
        raise ValueError("Please enter an RNA sequence.")
    invalid = sorted(set(sequence) - set("ACGU"))
    if invalid:
        raise ValueError("Invalid character(s): " + ", ".join(invalid)
                         + ". Use only A, C, G, and U.")
    return sequence


def can_pair(left, right, allow_wobble=True):
    """Return True when two nucleotides are allowed to form a base pair."""
    if (left, right) in CANONICAL_PAIRS:
        return True
    return allow_wobble and (left, right) in WOBBLE_PAIRS


def predict_structure(sequence, min_loop=3, allow_wobble=True):
    """
    Predict a simple structure using the Nussinov algorithm.

    The dynamic-programming table stores the maximum number of base pairs
    possible in each sequence interval. Traceback reconstructs one optimal
    non-crossing set of pairs.
    """
    sequence = clean_sequence(sequence)
    if min_loop < 0:
        raise ValueError("Minimum loop length cannot be negative.")

    n = len(sequence)
    dp = [[0] * n for _ in range(n)]

    # Fill intervals from shortest to longest.
    for span in range(1, n):
        for i in range(n - span):
            j = i + span
            best = max(dp[i + 1][j], dp[i][j - 1])

            if j - i > min_loop and can_pair(
                sequence[i], sequence[j], allow_wobble
            ):
                best = max(best, dp[i + 1][j - 1] + 1)

            # Split the interval into two independent subintervals.
            for k in range(i, j):
                best = max(best, dp[i][k] + dp[k + 1][j])

            dp[i][j] = best

    pairs = []

    def traceback(i, j):
        if i >= j:
            return

        if dp[i][j] == dp[i + 1][j]:
            traceback(i + 1, j)
            return

        if dp[i][j] == dp[i][j - 1]:
            traceback(i, j - 1)
            return

        if (j - i > min_loop
                and can_pair(sequence[i], sequence[j], allow_wobble)
                and dp[i][j] == dp[i + 1][j - 1] + 1):
            pairs.append((i + 1, j + 1))
            traceback(i + 1, j - 1)
            return

        for k in range(i, j):
            if dp[i][j] == dp[i][k] + dp[k + 1][j]:
                traceback(i, k)
                traceback(k + 1, j)
                return

    if n:
        traceback(0, n - 1)

    pairs.sort()
    structure = ["."] * n
    pair_details = []

    for left, right in pairs:
        structure[left - 1] = "("
        structure[right - 1] = ")"
        pair_type = sequence[left - 1] + "-" + sequence[right - 1]
        pair_details.append({
            "left_position": left,
            "right_position": right,
            "pair": pair_type
        })

    return {
        "sequence": sequence,
        "structure": "".join(structure),
        "base_pairs": pair_details,
        "number_of_pairs": len(pairs),
        "minimum_loop_length": min_loop,
        "wobble_pairing": bool(allow_wobble),
        "algorithm": "Nussinov dynamic programming",
        "time_complexity": "O(n^3)",
        "space_complexity": "O(n^2)"
    }


def print_result(result):
    """Print the prediction in a readable command-line format."""
    print("\nRNA SECONDARY STRUCTURE PREDICTION")
    print("-" * 38)
    print("RNA sequence :", result["sequence"])
    print("Structure    :", result["structure"])
    print("Base pairs   :", result["number_of_pairs"])
    print("\nPaired positions (1-based):")
    if result["base_pairs"]:
        for item in result["base_pairs"]:
            print("  Position {left_position} - Position {right_position}: {pair}".format(**item))
    else:
        print("  No base pairs predicted.")
    print("\nAlgorithm    :", result["algorithm"])
    print("Complexity   : Time O(n^3), Space O(n^2)")
    print("\nNote: This educational model maximizes allowed base pairs. It does not")
    print("calculate thermodynamic free energy or predict pseudoknots.")


def run_gui():
    """Launch a simple desktop interface using Python's Tkinter library."""
    try:
        import tkinter as tk
        from tkinter import messagebox, scrolledtext
    except ImportError:
        print("Tkinter is not available. Run the program in command-line mode.")
        return

    root = tk.Tk()
    root.title("RNA Secondary Structure Basic Predictor")
    root.geometry("720x590")
    root.minsize(580, 470)

    title = tk.Label(root, text="RNA Secondary Structure Basic Predictor",
                     font=("Arial", 17, "bold"))
    title.pack(pady=(16, 8))

    tk.Label(root, text="Enter an RNA sequence (A, C, G, U):",
             font=("Arial", 11)).pack(anchor="w", padx=20)
    sequence_box = scrolledtext.ScrolledText(root, height=5, wrap=tk.WORD,
                                              font=("Courier New", 12))
    sequence_box.pack(fill=tk.X, padx=20, pady=7)

    options = tk.Frame(root)
    options.pack(fill=tk.X, padx=20, pady=6)

    tk.Label(options, text="Minimum loop length:").pack(side=tk.LEFT)
    loop_value = tk.StringVar(value="3")
    tk.Spinbox(options, from_=0, to=20, width=5, textvariable=loop_value).pack(
        side=tk.LEFT, padx=(7, 20)
    )
    wobble_value = tk.BooleanVar(value=True)
    tk.Checkbutton(options, text="Allow G-U wobble pairs",
                   variable=wobble_value).pack(side=tk.LEFT)

    output = scrolledtext.ScrolledText(root, height=16, wrap=tk.WORD,
                                       font=("Courier New", 10), state=tk.DISABLED)
    output.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)

    def predict():
        try:
            loop_length = int(loop_value.get())
            result = predict_structure(sequence_box.get("1.0", tk.END),
                                       loop_length, wobble_value.get())
            lines = [
                "RNA sequence : " + result["sequence"],
                "Structure    : " + result["structure"],
                "Number pairs : " + str(result["number_of_pairs"]),
                "",
                "Paired positions (1-based):"
            ]
            if result["base_pairs"]:
                for item in result["base_pairs"]:
                    lines.append("  {left_position} - {right_position}: {pair}".format(**item))
            else:
                lines.append("  No base pairs predicted.")
            lines.extend([
                "",
                "Algorithm: " + result["algorithm"],
                "Time complexity: O(n^3)",
                "Space complexity: O(n^2)",
                "",
                "Educational limitation: this basic model maximizes the number",
                "of allowed pairs. It does not calculate free energy or model",
                "pseudoknots, so it is not a replacement for research tools."
            ])
            output.configure(state=tk.NORMAL)
            output.delete("1.0", tk.END)
            output.insert(tk.END, "\n".join(lines))
            output.configure(state=tk.DISABLED)
        except (ValueError, TypeError) as error:
            messagebox.showerror("Input error", str(error))

    tk.Button(root, text="Predict RNA Structure", command=predict,
              font=("Arial", 11, "bold"), padx=12, pady=6).pack(pady=(0, 14))
    root.mainloop()


def main():
    parser = argparse.ArgumentParser(
        description="Predict a basic RNA secondary structure using the Nussinov algorithm."
    )
    parser.add_argument("sequence", nargs="?",
                        help="RNA sequence (A, C, G, U); DNA T is converted to U.")
    parser.add_argument("--min-loop", type=int, default=3,
                        help="minimum number of unpaired bases in a hairpin loop (default: 3)")
    parser.add_argument("--no-wobble", action="store_true",
                        help="disable G-U wobble base pairing")
    parser.add_argument("--json", action="store_true",
                        help="print the result as JSON")
    args = parser.parse_args()

    if args.sequence is None:
        if args.json:
            parser.error("--json requires an RNA sequence.")
        run_gui()
        return

    try:
        result = predict_structure(args.sequence, args.min_loop,
                                   allow_wobble=not args.no_wobble)
    except ValueError as error:
        print("Error:", error, file=sys.stderr)
        raise SystemExit(2)

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print_result(result)


if __name__ == "__main__":
    main()
