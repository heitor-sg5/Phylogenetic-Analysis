# Phylogenetic Tree Analysis

This repository provides a range of scripts for computing evolutionary distances among sequences and reconstructing phylogenetic trees. It includes tools to build distance matrices from nucleotide sequences using multiple alignment and substitution models, and to infer trees via Unweighted Pair Group Method with Arithmetic Mean (UPGMA), Neighbor‑joining, and Nearest‑neighbour‑interchange (NNI) combined with small‑parsimony algorithms.

---

## 🧬 What is a Phylogenetic Tree?

A phylogenetic tree is a branching diagram that models the evolutionary history of a set of organisms or genes by depicting how they diverged from common ancestors based on shared genetic or phenotypic features. It conveys both the order and relative timing of speciation or duplication events, grouping taxa into clades that reflect common descent. 

Because the number of possible tree topologies grows super‑exponentially with the number of taxa, computational algorithms are essential to efficiently search this vast space and identify the tree that best explains the observed sequence data under explicit evolutionary models.

---

## 📁 Files in This Repository

- `msa_and_distance_matrix.py`: Aligns a set of sequences with MUSCLE, then computes either: the indel-free consensus; Hamming's pairwise distance; or Kimura's two-parameter distance.
- `upgma.py`: Implements the UPGMA algorithm to build an ultrametric tree from a distance matrix and output a rooted tree and edge lengths.
- `neighbour_joining.py`: Implements the Neighbor‑joining algorithm to infer an unrooted tree from a distance matrix, outputting the a tree with branch lengths.
- `small_parsimony_and_nni.py`: Applies NNI moves combined with the small‑parsimony score to optimize a tree topology for a set of character strings.

---

## ⚙️ How to Use

### 1. Prepare Input

The input may vary between scripts:

- To build the distance matrix or remove indels from aligned sequences using `msa_and_distance_matrix.py`: place your FASTA files in one directory and the script will prompt you to select that folder via a GUI dialog.
- For `upgma.py` and `neighbour_joining.py` the input will be a distance matrix of length n (e.g. [[0, 2, 1], [2, 0, 1], [1, 2, 0]]).
- For `small_parsimony_and_nni.py` use an indel-free list of (string) aligned sequences of equal length.

### 2. Run the Algorithms

Each script will print:

- A distance matrix using Kimura's or Hamming's distance model, or an indel-free sequence alignment, if using `msa_and_distance_matrix.py`.
- All other scripts will output the produced phylogenetic tree, clustering similar sequences close to one another.

---

#### MSA and Distance Matrix Builder 

  bash
```msa_and_distance_matrix.py```

#### UPGMA Tree Constructor 

  bash
```upgma.py```

#### Neighbour-joining Tree Constructor 

  bash
```neighbour_joining.py```

#### NNI with Small-parsimony Tree Constructor 

  bash
```small_parsimony_and_nni.py```

The variables D (distance matrix) and seqs (sequence list) is defined at the bottom of the tree-building scripts, and can be replaced by your own values.

---

## 🧠 Algorithm Overviews

### MSA and Distance Matrix Builder

- Applies multiple sequence alignment using MUSCLE.
- Calculates the (Kimura's or Hamming's) distance between every pair of sequences to build a 2D distance matrix.
- Removes columns containing indels from aligned sequences, producing indel-free sequences of equal length.
- Time complexity: O(n^2 * L)

### Unichromosomal Reversal Sorting with Breakpoints

- 
- Time complexity: O(n^3)

### Multichromosomal Two-break Sorting

- Converts genomes into breakpoint graphs using red (P) and blue (Q) edges.
- Identifies non-trivial cycles (edges not forming isolated cycles) and applies 2-breaks that reduce the number of such cycles.
- Updates the genome after each 2-break and repeats until P = Q, printing the transformation at each step.
- Time complexity: O(n^2)

---

## 🧪 Example Output

- Identity and Signed Permutations:

  Genome 1: [1, 2, 3, 4, 5]
  
  Genome 2: [1, -3, -2, 4, 5]
  
- Reversal Distance and Breakpoints:

  Step 1: [+1 +7 -9 +11 +10 +3 -2 -6 +5 -4 -8] | Breakpoints: 11
  
  Step 2: [+1 +2 -3 -10 -11 +9 -7 -6 +5 -4 -8] | Breakpoints: 9
  
  Step 3: [+1 +2 -3 -10 -11 +9 -7 -6 -5 -4 -8] | Breakpoints: 7
  
  Step 4 [+1 +2 +3 -10 -11 +9 -7 -6 -5 -4 -8] | Breakpoints: 6
  
  Step 5: [+1 +2 +3 +4 +5 +6 +7 -9 +11 +10 -8] | Breakpoints: 5
  
  Step 6 [+1 +2 +3 +4 +5 +6 +7 -11 +9 +10 -8] | Breakpoints: 4
  
  Step 7 [+1 +2 +3 +4 +5 +6 +7 +8 -10 -9 +11] | Breakpoints: 2
  
  Step 8: [+1 +2 +3 +4 +5 +6 +7 +8 +9 +10 +11] | Breakpoints: 0
  
  Reversal distance: 7

---

## 👤 Author

Heitor Gelain do Nascimento
Email: heitorgelain@outlook.com
GitHub: @heitor-sg5

---

## 📚 References

Bioinformatics Algorithms: An Active Learning Approach (Chapter 6) by
Phillip Compeau & Pavel Pevzner
https://bioinformaticsalgorithms.com
