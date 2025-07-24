import random

class Node:
    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.label = None
        self.score = {}
        self.tagged = False

    def add_child(self, child):
        self.children.append((child))

    def is_leaf(self):
        return len(self.children) == 0

def small_parsimony(root, alphabet):
    def postorder(node):
        if node.tagged:
            return
        for child in node.children:
            postorder(child)

        if node.is_leaf():
            node.tagged = True
            for k in alphabet:
                node.score[k] = 0 if node.label == k else float('inf')
        else:
            if len(node.children) != 2:
                raise ValueError(f"Node '{node.name}' does not have exactly 2 children (has {len(node.children)})")
            left, right = node.children
            node.tagged = True
            for k in alphabet:
                min_left = min(left.score[i] + (0 if i == k else 1) for i in alphabet)
                min_right = min(right.score[j] + (0 if j == k else 1) for j in alphabet)
                node.score[k] = min_left + min_right

    def assign_labels(node, parent_label=None):
        if node.is_leaf():
            return
        if parent_label is None:
            node.label = min(node.score, key=node.score.get)
        else:
            options = sorted(alphabet, key=lambda k: (node.score[k] + (0 if k == parent_label else 1)))
            node.label = options[0]
        for child in node.children:
            assign_labels(child, node.label)

    def reset_tags(node):
        node.tagged = False
        node.score = {}
        for child in node.children:
            reset_tags(child)

    reset_tags(root)
    postorder(root)
    assign_labels(root)
    return min(root.score.values())

def generate_random_tree(strings):
    leaves = [Node(name=i) for i in range(len(strings))]
    for node, label in zip(leaves, strings):
        node.label = label

    nodes = leaves[:]
    idx = len(strings)
    while len(nodes) > 1:
        a = nodes.pop(random.randint(0, len(nodes) - 1))
        b = nodes.pop(random.randint(0, len(nodes) - 1))
        parent = Node(name=f"internal_{idx}")
        idx += 1
        parent.add_child(a)
        parent.add_child(b)
        nodes.append(parent)
    return nodes[0]

def get_internal_edges(node, edges=None):
    if edges is None:
        edges = []
    for child in node.children:
        if not node.is_leaf() and not child.is_leaf():
            edges.append((node.name, child.name))
        get_internal_edges(child, edges)
    return edges

def deep_copy_tree(node):
    new_node = Node(name=node.name)
    new_node.label = node.label
    for child in node.children:
        new_node.add_child(deep_copy_tree(child))
    return new_node

def find_node(node, target_name):
    if node.name == target_name:
        return node
    for child in node.children:
        found = find_node(child, target_name)
        if found:
            return found
    return None

def generate_nni_variants(tree):
    variants = []
    edges = get_internal_edges(tree)

    for u_name, v_name in edges:
        tree_copy = deep_copy_tree(tree)
        u = find_node(tree_copy, u_name)
        v = find_node(tree_copy, v_name)

        if not u or not v:
            continue
        if v not in u.children or len(v.children) != 2:
            continue

        a, b = v.children
        if len(u.children) != 2:
            continue
        other = [child for child in u.children if child != v]
        if not other:
            continue

        for swap_child in (a, b):
            new_tree = deep_copy_tree(tree)
            u_new = find_node(new_tree, u_name)
            v_new = find_node(new_tree, v_name)

            if not u_new or not v_new:
                continue

            c_new = [child for child in u_new.children if child.name != v_new.name][0]
            a_new, b_new = v_new.children

            if swap_child.name == a_new.name:
                v_new.children = [b_new, c_new]
            else:
                v_new.children = [a_new, c_new]
            u_new.children = [child for child in u_new.children if child.name != v_new.name]
            u_new.children.append(swap_child)

            variants.append(new_tree)

    return variants

def nearest_neighbour_interchange(seqs):
    alphabet = sorted(set("".join(seqs)))
    tree = generate_random_tree(seqs)
    current_score = small_parsimony(tree, alphabet)
    best_tree = tree

    improved = True
    while improved:
        improved = False
        variants = generate_nni_variants(best_tree)
        for variant in variants:
            score = small_parsimony(variant, alphabet)
            if score < current_score:
                best_tree = variant
                current_score = score
                improved = True
                break
    return best_tree

def print_tree(node, indent=0):
    label_str = f" [{node.label}]" if node.label is not None else ""
    print("  " * indent + str(node.name) + label_str)
    for child in node.children:
        print("  " * (indent + 1) + f"|-- {child.name}")
        print_tree(child, indent + 2)

seqs = ['CCCATGTCAGGGCACGAGCGAGA', 
        'GCGTTTTCAGCAGTTATGTTGAT', 
        'GCCAGTGCAGGGGGCACGGAATA', 
        'GCTATGTCAGGGGGCACGAGACA', 
        'GCTATGTCAGGGGGCACGAGCAT']

final_tree = nearest_neighbour_interchange(seqs)
print_tree(final_tree)
