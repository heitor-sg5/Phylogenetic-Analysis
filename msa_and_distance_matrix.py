from tkinter import filedialog
import tkinter as tk
import subprocess
import tempfile
import math
import os

def muscle_align(sequences, muscle_path="muscle"):
    if not sequences:
        raise ValueError("The sequence list is empty.")
    
    if not os.path.isfile(muscle_path):
        raise FileNotFoundError(f"MUSCLE executable not found at {muscle_path}")

    with tempfile.TemporaryDirectory() as temp_dir:
        input_fasta = os.path.join(temp_dir, "input.fasta")
        output_fasta = os.path.join(temp_dir, "output.fasta")

        with open(input_fasta, 'w') as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq{i+1}\n{seq}\n")

        result = subprocess.run(
            [muscle_path, "-align", input_fasta, "-output", output_fasta],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            raise RuntimeError(f"MUSCLE failed:\n{result.stderr}")

        aligned_sequences = []
        with open(output_fasta, 'r') as f:
            seq = ""
            for line in f:
                if line.startswith(">"):
                    if seq:
                        aligned_sequences.append(seq)
                        seq = ""
                else:
                    seq += line.strip()
            if seq:
                aligned_sequences.append(seq)

        return aligned_sequences
    
def remove_indel_columns(aligned_seqs):
    if not aligned_seqs:
        return []
    transposed = list(zip(*aligned_seqs))
    columns_to_keep = [i for i, col in enumerate(transposed) if '-' not in col]
    new_seqs = [''.join(seq[i] for i in columns_to_keep) for seq in aligned_seqs]
    return new_seqs

def hamming_distance(v, w):
    return sum(1 for a, b in zip(v, w) if a != b)

def kimura_distance(v, w):
    transitions = 0
    transversions = 0
    valid_positions = 0

    transitions_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

    for a, b in zip(v, w):
        if a == '-' or b == '-':
            continue
        if a == b:
            valid_positions += 1
        else:
            valid_positions += 1
            if (a, b) in transitions_pairs:
                transitions += 1
            else:
                transversions += 1

    if valid_positions == 0:
        return 0.0

    P = transitions / valid_positions
    Q = transversions / valid_positions

    try:
        distance = -0.5 * math.log((1 - 2 * P - Q) * math.sqrt(1 - 2 * Q))
    except ValueError:
        distance = float('inf')

    return round (distance, 2)

def read_fasta_sequence(file_path):
    sequence = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            sequence.append(line.strip())
    return ''.join(sequence)

def get_fasta_files():
    root = tk.Tk()
    root.withdraw()

    folder_path = filedialog.askdirectory(title="Select Folder with FASTA Files")
    if not folder_path:
        print("No folder selected.")
        return []

    fasta_files = [f for f in os.listdir(folder_path) if f.lower().endswith(('.fasta', '.fa', '.fna'))]
    if not fasta_files:
        print("No FASTA files found in the folder.")
        return []

    sequences = []
    for file in fasta_files:
        file_path = os.path.join(folder_path, file)
        seq = read_fasta_sequence(file_path)
        sequences.append(seq)

    return sequences

def build_distance_matrix(seq_list, d):
    aligned_seqs = muscle_align(seq_list)
    if d == 2:
        return remove_indel_columns(aligned_seqs)
        
    n = len(aligned_seqs)
    dist_matrix = [[0]*n for _ in range(n)]

    for i in range(n):
        for j in range(i+1, n):
            if d == 0:
                dist = hamming_distance(aligned_seqs[i], aligned_seqs[j])
            else:
                dist = kimura_distance(aligned_seqs[i], aligned_seqs[j])

            dist_matrix[i][j] = dist
            dist_matrix[j][i] = dist

    return dist_matrix

seqs = muscle_align(seqs)
print(build_distance_matrix(seqs, 0))
