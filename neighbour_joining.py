import numpy as np

class Node:
    def __init__(self, name):
        self.name = name
        self.children = []

    def add_child(self, child, length):
        self.children.append((child, length))

def neighbor_joining(D, labels):
    n = len(D)
    if n == 2:
        node = Node(f"({labels[0]},{labels[1]})")
        node.add_child(Node(labels[0]), D[0,1] / 2)
        node.add_child(Node(labels[1]), D[0,1] / 2)
        return node

    total_dist = np.sum(D, axis=1)

    D_star = np.full((n,n), np.inf)
    for i in range(n):
        for j in range(n):
            if i != j:
                D_star[i,j] = (n - 2) * D[i,j] - total_dist[i] - total_dist[j]

    i, j = np.unravel_index(np.argmin(D_star), D_star.shape)
    if j < i:
        i, j = j, i

    delta = (total_dist[i] - total_dist[j]) / (n - 2)
    limb_i = 0.5 * (D[i,j] + delta)
    limb_j = 0.5 * (D[i,j] - delta)

    new_label = f"({labels[i]},{labels[j]})"

    new_row = []
    for k in range(n):
        if k != i and k != j:
            dist = 0.5 * (D[i,k] + D[j,k] - D[i,j])
            new_row.append(dist)

    indices = [x for x in range(n) if x != i and x != j]
    new_D = np.zeros((n-1, n-1))
    for a, idx_a in enumerate(indices):
        for b, idx_b in enumerate(indices):
            new_D[a,b] = D[idx_a, idx_b]
    for a in range(n-1 - 1):
        new_D[a, -1] = new_D[-1, a] = new_row[a]
    new_D[-1, -1] = 0.0

    new_labels = [labels[x] for x in indices] + [new_label]
    subtree = neighbor_joining(new_D, new_labels)

    def find_node(node, target_label):
        if node.name == target_label:
            return node
        for child, _ in node.children:
            found = find_node(child, target_label)
            if found:
                return found
        return None

    new_node = find_node(subtree, new_label)
    new_node.add_child(Node(labels[i]), limb_i)
    new_node.add_child(Node(labels[j]), limb_j)

    return subtree

def print_tree(node, indent=0):
    print("  " * indent + str(node.name))
    for child, length in node.children:
        print("  " * (indent + 1) + f"|-- {child.name} (len={length:.2f})")
        print_tree(child, indent + 2)

D = [
    [0, 295, 306, 497, 1081, 1091, 1003, 956, 954],
    [295, 0, 309, 500, 1084, 1094, 1006, 959, 957],
    [306, 309, 0, 489, 1073, 1083, 995, 948, 946],
    [497, 500, 489, 0, 1092, 1102, 1014, 967, 965],
    [1081, 1084, 1073, 1092, 0, 818, 1056, 1053, 1051],
    [1091, 1094, 1083, 1102, 818, 0, 1066, 1063, 1061],
    [1003, 1006, 995, 1014, 1056, 1066, 0, 975, 973],
    [956, 959, 948, 967, 1053, 1063, 975, 0, 16],
    [954, 957, 946, 965, 1051, 1061, 973, 16, 0]
]

tree = neighbor_joining(np.array(D), list(range(len(D))))
print_tree(tree)
