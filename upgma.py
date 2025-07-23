import numpy as np

class Node:
    def __init__(self, name):
        self.name = name
        self.children = []
        self.age = 0.0

    def add_child(self, child):
        self.children.append(child)

def upgma(D, n):
    D = np.array(D, dtype=float)
    clusters = {i: [i] for i in range(n)}
    nodes = {i: Node(i) for i in range(n)}

    ages = {i: 0.0 for i in range(n)}

    current_id = n

    while len(clusters) > 1:
        min_dist = float('inf')
        to_merge = None
        cluster_keys = list(clusters.keys())
        for i in range(len(cluster_keys)):
            for j in range(i + 1, len(cluster_keys)):
                ci, cj = cluster_keys[i], cluster_keys[j]
                dist = 0.0
                for a in clusters[ci]:
                    for b in clusters[cj]:
                        dist += D[a][b]
                dist /= (len(clusters[ci]) * len(clusters[cj]))
                if dist < min_dist:
                    min_dist = dist
                    to_merge = (ci, cj)

        ci, cj = to_merge
        new_cluster = clusters[ci] + clusters[cj]

        new_node = Node(current_id)
        new_node.add_child(nodes[ci])
        new_node.add_child(nodes[cj])
        new_node.age = min_dist / 2

        nodes[current_id] = new_node
        ages[current_id] = new_node.age

        new_dist_row = []
        for ck in cluster_keys:
            if ck != ci and ck != cj:
                dist = 0.0
                for a in new_cluster:
                    for b in clusters[ck]:
                        dist += D[a][b]
                dist /= (len(new_cluster) * len(clusters[ck]))
                new_dist_row.append((ck, dist))

        clusters.pop(ci)
        clusters.pop(cj)
        clusters[current_id] = new_cluster

        D = np.vstack([D, np.zeros((1, D.shape[1]))])
        D = np.hstack([D, np.zeros((D.shape[0], 1))])
        for ck, dist in new_dist_row:
            D[current_id][ck] = dist
            D[ck][current_id] = dist

        current_id += 1

    root_id = list(clusters.keys())[0]
    root = nodes[root_id]

    edges = []
    def collect_edges(node):
        for child in node.children:
            length = node.age - child.age
            edges.append((node.name, child.name, length))
            collect_edges(child)
    collect_edges(root)

    return root, edges

def print_tree(node, indent=0):
    print("  " * indent + str(node.name))
    for child in node.children:
        length = node.age - child.age
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
n = len(D)

root, edges = upgma(D, n)

print_tree(root)
