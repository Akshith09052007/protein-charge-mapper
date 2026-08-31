# 🧬 RNA Secondary Structure Basic Predictor

**Authors**
- R. Akshith Narasimha — SE25UBIT048
- D. Anjana Reddy — SE25UBIT048

## Overview
A full-stack educational RNA secondary-structure predictor using the **Nussinov dynamic-programming algorithm**. It predicts dot-bracket structure by maximizing compatible base pairs.

## Features
- Clean browser interface
- Flask REST backend
- JSON prediction API
- A/C/G/U validation; T automatically becomes U
- Canonical A-U and G-C pairs
- Optional G-U wobble pairs
- Configurable minimum loop length
- Base-pair positions and pair types
- Health-check endpoint
- CLI mode

## Run locally
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python app.py
```
Then open **http://localhost:5000**.

## CLI
```bash
python rna_predictor.py GGGAAACCC
python rna_predictor.py GGGAAACCC --json
python rna_predictor.py GGGAAACCC --min-loop 4 --no-wobble
```

## API
### `GET /api/health`
Returns service status.

### `POST /api/predict`
Request:
```json
{"sequence":"GGGAAACCC","min_loop":3,"allow_wobble":true}
```
Response includes normalized sequence, length, dot-bracket structure, pair count, paired positions, pair types, algorithm and complexity.

## Algorithm
Nussinov considers leaving either endpoint unpaired, pairing compatible endpoints, or splitting an interval into two subproblems. Complexity is **O(n³) time** and **O(n²) space**.

## Limitations
This is a basic educational predictor. It does not calculate thermodynamic free energy, pseudoknots, experimental constraints, or the full set of biological folding effects. Predictions should not be treated as biological ground truth.
