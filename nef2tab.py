# Python script to convert a NEF file to .tab format for TALOS+
# Corrected version to properly handle n-1 and n values, and use n-1 only to fill missing residues.
# Added: ignores unassigned residues (preceded by @ in sequence_code or chain_code)
# Modified: adds protein name from NEF file to REMARK, and writes DATA SEQUENCE on every sequence line

import re

# 3-letter to 1-letter amino acid conversion dictionary
AA_3_to_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    # Add variants if needed
    '.': '.', '': ''
}

# Atoms of interest for TALOS+ (include 'H' for H->HN mapping)
TALOS_ATOMS = {'HN', 'N', 'CA', 'CB', 'C', 'HA', 'HA2', 'HA3', 'H'}

def parse_nef_chemical_shifts(nef_file):
    """
    Parse the chemical_shift_list block from the NEF file and return a list of shifts.
    Each shift is a dict with keys: chain_code, sequence_code, residue_name, atom_name, value
    Ignores unassigned residues (preceded by @ in sequence_code or chain_code)
    """
    shifts = []
    in_loop = False
    columns = []
    with open(nef_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('save_nef_chemical_shift_list'):
                in_block = True
            if line.startswith('loop_'):
                in_loop = True
                columns = []
                continue
            if in_loop and line.startswith('_nef_chemical_shift.'):
                columns.append(line)
                continue
            if in_loop and columns and (line == 'stop_' or line.startswith('save_')):
                in_loop = False
                continue
            if in_loop and columns and not line.startswith('_'):
                # Data
                fields = re.split(r'\s+', line)
                if len(fields) < len(columns):
                    continue
                entry = {}
                for i, col in enumerate(columns):
                    key = col.replace('_nef_chemical_shift.', '')
                    entry[key] = fields[i]
                # Ignore unassigned residues (preceded by @ in sequence_code or chain_code)
                seq_code = entry.get('sequence_code', '')
                chain_code = entry.get('chain_code', '')
                if seq_code.startswith('@') or chain_code.startswith('@'):
                    continue
                shifts.append(entry)
    return shifts

def get_sequence_from_user():
    print("Please enter the full protein sequence (1-letter code, no spaces):")
    seq = input().strip().upper().replace(" ", "")
    return seq

def extract_protein_name_from_nef(nef_file):
    """
    Extracts the protein name from the NEF file, looking for a line like:
    ; Substance:Sec5.
    or
         Substance:BamE.
    Returns the name (e.g., 'Sec5' or 'BamE'), or None if not found.
    This function is robust to leading/trailing whitespace and to the presence/absence of the leading semicolon.
    It also scans the entire file, including the end.
    """
    protein_name = None
    # Regex matches lines like:
    # ; Substance:BamE.
    #    Substance:BamE.
    # or with extra spaces
    pattern = re.compile(r'^[; ]*\s*Substance:([A-Za-z0-9_.-]+)\.', re.IGNORECASE)
    with open(nef_file, 'r') as f:
        for line in f:
            l = line.strip()
            m = pattern.match(l)
            if m:
                protein_name = m.group(1)
                # Do not break; keep last match in case of multiple
    return protein_name

def convert_shifts_to_talos(shifts, sequence):
    """
    Converts the list of NEF shifts to TALOS+ format (RESID, RESNAME, ATOMNAME, SHIFT)
    Handles cases where the residue is unknown (n-1) using the provided sequence.
    Correctly maps H (NEF) -> HN (TAB) in all cases, including n-1.
    Uses n-1 values only to fill missing residues, without overwriting n values.
    Ignores unassigned residues (preceded by @ in sequence_code or chain_code)
    """
    talos_rows = []
    # Index shifts by (resid, atom) and (resid_str, atom) to handle n-1
    shift_by_resid_atom = {}      # {resid: {atom: {...}}} for normal residues
    shift_by_residstr_atom = {}   # {resid_str: {atom: {...}}} for n-1

    for entry in shifts:
        atom = entry.get('atom_name', '').strip()
        if atom not in TALOS_ATOMS:
            continue
        # Ignore unassigned residues (preceded by @ in sequence_code or chain_code)
        resid_str_raw = entry.get('sequence_code', '')
        chain_code_raw = entry.get('chain_code', '')
        if resid_str_raw.startswith('@') or chain_code_raw.startswith('@'):
            continue
        resid_str = resid_str_raw.replace('@', '')
        # Keep the raw version for n-1 (e.g., '19-1')
        # Extract the main residue number (e.g., '19' from '19-1')
        m = re.match(r'(\d+)(?:-1)?$', resid_str)
        if not m:
            continue
        resid_int = int(m.group(1))
        resname3 = entry.get('residue_name', '').strip()
        try:
            shift = float(entry.get('value', '.'))
        except Exception:
            continue
        # Index by number (for normal residues) and by raw key (for n-1)
        if not resid_str.endswith('-1'):
            shift_by_resid_atom.setdefault(resid_int, {})[atom] = {
                'resname3': resname3,
                'shift': shift,
                'resid_str': resid_str
            }
        else:
            shift_by_residstr_atom.setdefault(resid_str, {})[atom] = {
                'resname3': resname3,
                'shift': shift
            }

    # Iterate over the sequence to generate TALOS lines
    seq_len = len(sequence)
    for i, aa1 in enumerate(sequence):
        resid = i + 1
        atom_dict = shift_by_resid_atom.get(resid, {})
        aa1_forced = aa1

        # 1. Add normal shifts for this residue
        for atom, d in atom_dict.items():
            # TALOS mapping: H (NEF) -> HN (TAB), N -> N, others unchanged
            if atom == 'H':
                atomname = 'HN'
            elif atom == 'HN':
                atomname = 'HN'
            elif atom == 'N':
                atomname = 'N'
            else:
                atomname = atom
            talos_rows.append((resid, aa1_forced, atomname, d['shift']))

        # 2. Add n-1 shifts to fill the previous residue (resid-1), if needed
        if i > 0:
            resid_n1 = resid - 1
            resid_n1_str = f"{resid}-1"
            atom_dict_n1 = shift_by_residstr_atom.get(resid_n1_str, {})
            for atom, d in atom_dict_n1.items():
                # Only take shifts where the residue is unknown (resname3 == '.' or '')
                if d['resname3'] not in ('.', ''):
                    continue
                # Only add the n-1 value if the previous residue does not already have a shift for this atom
                atom_dict_prev = shift_by_resid_atom.get(resid_n1, {})
                if atom in atom_dict_prev:
                    continue
                # TALOS mapping: H (NEF) -> HN (TAB), N -> N, others unchanged
                if atom == 'H':
                    atomname = 'HN'
                elif atom == 'HN':
                    atomname = 'HN'
                elif atom == 'N':
                    atomname = 'N'
                else:
                    atomname = atom
                # The residue name is taken from the sequence at position resid-1
                aa1_prev = sequence[resid_n1 - 1]
                talos_rows.append((resid_n1, aa1_prev, atomname, d['shift']))

    # Sort by residue number then atom name
    talos_rows.sort(key=lambda x: (x[0], x[2]))
    return talos_rows

def write_talos_tab(talos_rows, output_file, sequence, first_resid=1, protein_name=None):
    """
    Write the shifts in TALOS+ .tab format with SEQUENCE and FIRST_RESID headers.
    - Adds REMARK with protein name if provided.
    - Writes DATA SEQUENCE at the start of every sequence line (not just the first).
    """
    with open(output_file, 'w') as f:
        if protein_name:
            f.write(f"REMARK Chemical Shift Table for {protein_name}\n\n")
        else:
            f.write("REMARK Chemical Shift Table\n\n")
        f.write(f"DATA FIRST_RESID {first_resid}\n\n")
        # Write the sequence in lines of 10 aa, 5 lines per block (50 aa per line)
        seq_chunks = [sequence[i:i+10] for i in range(0, len(sequence), 10)]
        # Group chunks into lines of 5 chunks (50 aa per line)
        lines = []
        for i in range(0, len(seq_chunks), 5):
            line = " ".join(seq_chunks[i:i+5])
            lines.append(line)
        for line in lines:
            f.write(f"DATA SEQUENCE {line}\n")
        f.write("\n")
        f.write("VARS   RESID RESNAME ATOMNAME SHIFT\n")
        f.write("FORMAT %4d   %1s     %4s      %8.3f\n\n")
        for resid, resname, atomname, shift in talos_rows:
            f.write(f"{resid:4d} {resname:1s} {atomname:4s} {shift:8.3f}\n")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python convert_nef_to_talos.py <input.nef> <output.tab>")
        sys.exit(1)
    nef_file = sys.argv[1]
    output_file = sys.argv[2]
    # Ask the user for the sequence
    sequence = get_sequence_from_user()
    # Extract protein name from NEF file
    protein_name = extract_protein_name_from_nef(nef_file)
    shifts = parse_nef_chemical_shifts(nef_file)
    talos_rows = convert_shifts_to_talos(shifts, sequence)
    write_talos_tab(talos_rows, output_file, sequence, protein_name=protein_name)
    print(f"Conversion complete: {output_file}")