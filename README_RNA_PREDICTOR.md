# Basic RNA Secondary Structure Predictor

**Authors**
- R. Akshith Narasimha — SE25UBIT048
- D. Anjana Reddy — SE25UBIT048

## Overview

This project implements a basic RNA secondary-structure predictor using the **Nussinov dynamic programming algorithm**. Given an RNA sequence, it predicts a secondary structure in **dot-bracket notation** by maximizing the number of compatible base pairs.

## Features

- Accepts RNA sequences containing A, C, G and U.
- Predicts canonical A-U and G-C pairs.
- Includes G-U wobble pairs.
- Enforces a minimum loop length of 3 nucleotides.
- Reports predicted dot-bracket structure and paired positions.
- Provides clear input validation and an educational explanation.

## Requirements

- Python 3.8+
- No external packages required.

## Run

```bash
python3 rna_predictor.py
```

Example input:

```text
GGGAAACCC
```

The program prints the input sequence, predicted dot-bracket structure, number of base pairs, and paired nucleotide positions.

## Algorithm

For a sequence of length `n`, the Nussinov dynamic-programming table stores the maximum number of base pairs possible for every subsequence. Each state considers leaving either end unpaired, pairing compatible bases, or splitting the subsequence into two smaller subsequences. A traceback reconstructs the predicted base pairs.

### Complexity

- Time: **O(n^3)**
- Space: **O(n^2)**

## Scope and limitations

This is a **basic educational predictor**, not a production-grade RNA folding package. It does not use experimentally measured constraints, pseudoknot prediction, or thermodynamic free-energy parameters. Consequently, its structure should not be treated as a biological ground truth.
