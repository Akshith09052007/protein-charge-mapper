# RNA Secondary Structure Basic Predictor

A Python-only educational bioinformatics project that predicts a basic RNA secondary structure using the **Nussinov dynamic programming algorithm**. The project includes a Python desktop interface (Tkinter), a Python prediction engine, a Python command-line interface, and Python unit tests.

## Authors
- R. Akshith Narasimha — SE25UBIT048
- D. Anjana Reddy — SE25UBIT008

## Technologies
- **Python 3.10+** — all application logic, validation, prediction, and tests.
- **Tkinter** — desktop graphical interface (included with most Python installations).
- **Python standard library only** — no Flask, HTML, CSS, JavaScript, or third-party Python packages are required.

## Features
- Accepts RNA sequences containing A, C, G, and U (DNA T is automatically converted to U).
- Predicts a maximum-base-pair structure in dot-bracket notation.
- Supports canonical A-U and G-C pairs, with optional G-U wobble pairs.
- Lets the user choose the minimum hairpin loop length.
- Displays paired nucleotide positions and pair types.
- Includes input validation, a desktop interface, a command-line interface, and automated unit tests.
- Runs locally without an internet connection or external API.

## Run the desktop application

Python 3.10 or newer is recommended. Tkinter should be available in your Python installation.

```bash
# Clone the repository
git clone https://github.com/Akshith09052007/protein-charge-mapper.git
cd protein-charge-mapper
git checkout rna-secondary-structure-predictor

# Launch the desktop application
python app.py
```

A desktop window will open. Enter an RNA sequence, choose the minimum loop length, optionally enable G-U wobble pairing, and click **Predict RNA Structure**.

## Command-line usage

```bash
python rna_predictor.py GGGAAACCC
python rna_predictor.py GGGAAACCC --min-loop 3 --no-wobble
python rna_predictor.py GGGAAACCC --json
```

Omit the sequence to enter it interactively.

## Run tests

```bash
python -m unittest discover -s tests -v
```

## Algorithm

The Nussinov algorithm fills a dynamic-programming table. For each sequence interval it considers:
1. Leaving the left nucleotide unpaired.
2. Leaving the right nucleotide unpaired.
3. Pairing the two ends when the bases are compatible and the loop constraint is satisfied.
4. Splitting the interval into two smaller independent intervals.

A traceback reconstructs the base pairs. Positions in the output are **1-based**.

- Time complexity: **O(n³)**
- Space complexity: **O(n²)**

## Scientific scope and limitations

This is a basic educational predictor, not a production-grade RNA folding package. It maximizes the number of allowed base pairs; it does not use thermodynamic free-energy parameters, model pseudoknots, or guarantee the experimentally observed structure. For real research, compare results with validated tools such as ViennaRNA and experimental evidence.
