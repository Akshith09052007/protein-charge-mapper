#!/usr/bin/env python3
"""Python-only desktop interface for the RNA secondary-structure predictor."""
import tkinter as tk
from tkinter import messagebox, ttk

from rna_predictor import predict_structure


class RNAPredictorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("RNA Secondary Structure Basic Predictor")
        self.root.geometry("850x700")
        self.root.minsize(650, 550)

        main = ttk.Frame(root, padding=18)
        main.pack(fill="both", expand=True)

        ttk.Label(
            main, text="RNA Secondary Structure Predictor",
            font=("Arial", 19, "bold")
        ).pack(anchor="w")
        ttk.Label(
            main,
            text="Predict a basic RNA secondary structure using the Nussinov algorithm."
        ).pack(anchor="w", pady=(4, 14))

        ttk.Label(main, text="RNA sequence (A, C, G, U; DNA T is also accepted):").pack(anchor="w")
        self.sequence_box = tk.Text(main, height=5, wrap="word", font=("Courier", 12))
        self.sequence_box.pack(fill="x", pady=(5, 12))
        self.sequence_box.insert("1.0", "GGGAAACCC")

        options = ttk.Frame(main)
        options.pack(fill="x", pady=(0, 12))
        ttk.Label(options, text="Minimum loop length:").pack(side="left")
        self.loop_value = tk.StringVar(value="3")
        ttk.Spinbox(
            options, from_=0, to=50, width=5, textvariable=self.loop_value
        ).pack(side="left", padx=(6, 18))
        self.wobble_value = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            options, text="Allow G-U wobble pairs", variable=self.wobble_value
        ).pack(side="left")

        ttk.Button(
            main, text="Predict RNA Structure", command=self.run_prediction
        ).pack(anchor="w", pady=(0, 12))

        ttk.Label(main, text="Prediction results:", font=("Arial", 12, "bold")).pack(anchor="w")
        output_frame = ttk.Frame(main)
        output_frame.pack(fill="both", expand=True, pady=(5, 12))
        self.output = tk.Text(
            output_frame, wrap="word", font=("Courier", 11), state="disabled"
        )
        scrollbar = ttk.Scrollbar(output_frame, command=self.output.yview)
        self.output.configure(yscrollcommand=scrollbar.set)
        self.output.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        ttk.Label(
            main,
            text="R. Akshith Narasimha — SE25UBIT048  |  D. Anjana Reddy — SE25UBIT008"
        ).pack(anchor="center")
        self.show_output(
            "Enter an RNA sequence and click 'Predict RNA Structure' to begin."
        )

    def show_output(self, text):
        self.output.configure(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)
        self.output.configure(state="disabled")

    def run_prediction(self):
        sequence = self.sequence_box.get("1.0", "end").strip()
        try:
            min_loop = int(self.loop_value.get())
            result = predict_structure(
                sequence, min_loop=min_loop, allow_wobble=self.wobble_value.get()
            )
            if result["length"] > 300:
                raise ValueError(
                    "For this basic desktop version, please use sequences of 300 "
                    "nucleotides or fewer."
                )
        except (ValueError, TypeError) as error:
            messagebox.showerror("Input error", str(error))
            return

        lines = [
            "PREDICTION SUMMARY",
            "=" * 55,
            f"RNA sequence : {result['sequence']}",
            f"Length       : {result['length']} nucleotides",
            f"Structure    : {result['structure']}",
            f"Base pairs   : {result['base_pair_count']}",
            f"Algorithm    : {result['algorithm']}",
            f"Minimum loop : {result['minimum_loop_length']}",
            f"G-U wobble   : {'Enabled' if result['wobble_pairs_enabled'] else 'Disabled'}",
            "",
            "BASE-PAIR POSITIONS (1-based)",
            "-" * 55,
        ]
        if result["base_pairs"]:
            for pair in result["base_pairs"]:
                lines.append(
                    f"{pair['position_1']} - {pair['position_2']}  "
                    f"{pair['bases']}  ({pair['type']})"
                )
        else:
            lines.append("No base pairs were predicted.")
        lines.extend([
            "",
            "COMPUTATIONAL COMPLEXITY",
            "-" * 55,
            f"Time: {result['complexity']['time']}",
            f"Space: {result['complexity']['space']}",
            "",
            "SCIENTIFIC LIMITATION",
            "-" * 55,
            result["warning"],
        ])
        self.show_output("\n".join(lines))


def main():
    root = tk.Tk()
    RNAPredictorApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
