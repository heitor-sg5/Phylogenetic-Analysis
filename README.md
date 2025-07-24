# Phylogenetic Tree Constructor

This repository provides a range of scripts for computing evolutionary distances among sequences and reconstructing phylogenetic trees. It includes tools to build distance matrices from nucleotide sequences using multiple alignment and substitution models, and to infer trees via Unweighted Pair Group Method with Arithmetic Mean (UPGMA), Neighbor‑joining, and Nearest‑neighbour‑interchange (NNI) combined with small‑parsimony algorithms.

---

## 🧬 What is a Phylogenetic Tree?

A phylogenetic tree is a branching diagram that models the evolutionary history of a set of organisms or genes by depicting how they diverged from common ancestors based on shared genetic or phenotypic features. It conveys both the order and relative timing of speciation or duplication events, grouping taxa into clades that reflect common descent. 

Because the number of possible tree topologies grows exponentially with the number of taxa, computational algorithms are essential to efficiently search this vast space and identify the tree that best explains the observed sequence data under explicit evolutionary models.

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

The variables D (distance matrix) and seqs (sequence list) is defined at the bottom of the tree-building scripts, and can be replaced by your own values. To select the type of distance matrix to build, change d to 0 (for Hamming's), 1 (for Kimura's), or 2 (for indel removal).

---

## 🧠 Algorithm Overviews

### MSA and Distance Matrix Builder

- Applies multiple sequence alignment using MUSCLE.
- Calculates the (Kimura's or Hamming's) distance between every pair of sequences to build a 2D distance matrix.
- Removes columns containing indels from aligned sequences, producing indel-free sequences of equal length.
- Time complexity: O(n^2 * L)

### UPGMA Tree Constructor 

- At each step, merges the two clusters with smallest average distance.
- Creates a new internal node with age = half the merge distance.
- Recomputes distances to the new cluster as the arithmetic mean.
- Time complexity: O(n^3)

### Neighbour-joining Tree Constructor

- Computes a transformed distance matrix 𝐷* to correct for unequal rates using; D*{ij} = (n - 2) × D{ij} - Σ_k D_{ik} - Σ_k D_{jk}.
- Joins the pair (i, j) minimizing D*{ij}, attaches them to a new node with calculated limb lengths.
- Reduces matrix size and repeats until two taxa remain.
- Time complexity: O(n^3)

### NNI with Small-parsimony Tree Constructor

- Small‑parsimony computes the minimum number of character changes on a tree for a given assignment of labels to leaves.
- NNI moves swap adjacent subtrees around an internal edge to explore neighboring topologies.
- Iteratively evaluates all single‑move variants and adopts any that reduce the parsimony score, stopping at a local optimum.
- Time complexity: O(I * n^2)

---

## 🧪 Example Output

- Distance Matrix:

  [[0, 295, 306, 497, 1081, 1091, 1003, 956, 954],
  [295, 0, 309, 500, 1084, 1094, 1006, 959, 957],
  [306, 309, 0, 489, 1073, 1083, 995, 948, 946],
  [497, 500, 489, 0, 1092, 1102, 1014, 967, 965],
  [1081, 1084, 1073, 1092, 0, 818, 1056, 1053, 1051],
  [1091, 1094, 1083, 1102, 818, 0, 1066, 1063, 1061],
  [1003, 1006, 995, 1014, 1056, 1066, 0, 975, 973],
  [956, 959, 948, 967, 1053, 1063, 975, 0, 16],
  [954, 957, 946, 965, 1051, 1061, 973, 16, 0]]
  
- Phylogenetic Tree:

  16
  |-- 13 (len=128.50)
    13
      |-- 4 (len=409.00)
        4
      |-- 5 (len=409.00)
        5
  |-- 15 (len=40.33)
    15
      |-- 6 (len=497.17)
        6
      |-- 14 (len=18.92)
        14
          |-- 9 (len=470.25)
            9
              |-- 7 (len=8.00)
                7
              |-- 8 (len=8.00)
                8
          |-- 12 (len=230.58)
            12
              |-- 3 (len=247.67)
                3
              |-- 11 (len=93.92)
                11
                  |-- 2 (len=153.75)
                    2
                  |-- 10 (len=6.25)
                    10
                      |-- 0 (len=147.50)
                        0
                      |-- 1 (len=147.50)
                        1

---

## 👤 Author

Heitor Gelain do Nascimento
Email: heitorgelain@outlook.com
GitHub: @heitor-sg5

---

## 📚 References

Bioinformatics Algorithms: An Active Learning Approach (Chapter 7) by
Phillip Compeau & Pavel Pevzner
https://bioinformaticsalgorithms.com
