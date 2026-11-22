# RNA Secondary Structure Basic Predictor

A beginner-friendly bioinformatics project that predicts a basic RNA secondary structure using the **Nussinov dynamic programming algorithm**. It provides a browser interface, a JSON API, and a command-line interface.

## Authors
- R. Akshith Narasimha — SE25UBIT048
- D. Anjana Reddy — SE25UBIT008

## Features
- Accepts RNA sequences containing A, C, G, and U (DNA T is automatically converted to U).
- Predicts a maximum-base-pair structure in dot-bracket notation.
- Supports canonical A-U and G-C pairs, with optional G-U wobble pairs.
- Lets the user choose the minimum hairpin loop length.
- Displays paired nucleotide positions and pair types.
- Includes input validation, a health endpoint, and automated unit tests.
- Runs locally; no external bioinformatics service or API key is required.

## Run locally

Python 3.10 or newer is recommended.

```bash
# 1. Clone the repository and enter this branch
git clone https://github.com/Akshith09052007/protein-charge-mapper.git
cd protein-charge-mapper
git checkout rna-secondary-structure-predictor

# 2. Create and activate a virtual environment
python -m venv .venv
# macOS / Linux:
source .venv/bin/activate
# Windows PowerShell:
# .venv\Scripts\Activate.ps1

# 3. Install dependencies
python -m pip install -r requirements.txt

# 4. Start the web application
python app.py
```

Open **http://127.0.0.1:5000** in your browser.

## Command-line usage

```bash
python rna_predictor.py GGGAAACCC
python rna_predictor.py GGGAAACCC --min-loop 3 --no-wobble
python rna_predictor.py GGGAAACCC --json
```

Omit the sequence to enter it interactively.

## API

### Health check
`GET /api/health`

### Predict
`POST /api/predict`

Example request:

```json
{
  "sequence": "GGGAAACCC",
  "min_loop": 3,
  "allow_wobble": true
}
```

A successful response contains `success: true` and a `result` object with the normalized sequence, dot-bracket structure, base-pair list, settings, algorithm, and complexity. Invalid input returns HTTP 400 with an explanatory error.

## Algorithm

The Nussinov algorithm fills a dynamic-programming table. For each sequence interval it considers:
1. Leaving the left nucleotide unpaired.
2. Leaving the right nucleotide unpaired.
3. Pairing the two ends when the bases are compatible and the loop constraint is satisfied.
4. Splitting the interval into two smaller independent intervals.

A traceback reconstructs the base pairs. Positions in the output are **1-based**.

- Time complexity: **O(n³)**
- Space complexity: **O(n²)**

## Run tests

```bash
python -m unittest discover -s tests -v
```

## Scientific scope and limitations

This is a basic educational predictor, not a production-grade RNA folding package. It maximizes the number of allowed base pairs; it does not use thermodynamic free-energy parameters, model pseudoknots, or guarantee the experimentally observed structure. For real research, compare results with validated tools such as ViennaRNA and experimental evidence.
