# nef2tab

Convert NEF files from CCPNMR to `.tab` format for TALOS-N

`nef2tab.py` is a Python script designed to convert NMR-STAR NEF files exported from [CCPNMR Analysis](httpswww.ccpn.ac.uk) (after backbone chemical shift assignment) into `.tab` format compatible with the [TALOS-N Server](httpsspin.niddk.nih.govbaxnmrservertalosn).

This tool facilitates the prediction of protein backbone and sidechain torsion angles from assigned chemical shifts.

---

## Features

- Converts NEF chemical shift data to TALOS-N `.tab` format.
- Filters out unassigned residues.
- Automatically includes `n-1` residue shifts to help with incomplete assignments.
- Extracts protein name from the NEF file and adds it as a REMARK.
- Requires only the NEF file and the 1-letter protein sequence.

---

## Installation

This script requires only Python 3.6

## Usage

Run the script with a NEF file as input:

```bash

python nef2tab.py input.nef output.tab
```
You will be prompted to enter your protein sequence as a one-letter code (no spaces).

## Example:
```
Please enter the full protein sequence (1-letter code, no spaces):
GSSHHHHHHENLYFQG...
```
---

## Output

The script generates a `.tab` file formatted for the TALOS-N server, containing:

- One line per chemical shift
- Columns: Residue Number, Residue Name, Atom Type, Chemical Shift
- REMARK line with the protein name (from NEF)
- Automatic skipping of unassigned residues
- Inclusion of residue n-1 shifts when available
