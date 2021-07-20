from nj import *
from collections import defaultdict
import treeswift as tsw
d = [
    [0, 5, 9, 9, 8],
    [5, 0, 10, 10, 9],
    [9, 10, 0, 8, 7],
    [9, 10, 8, 0, 3],
    [8, 9, 7, 3, 0],
]

def parse_d(ts, d):
    D = defaultdict(dict)
    for i, coll in zip(ts, d):
        for j, v in zip(ts, coll):
            D[i][j] = v
    return dict(D)

tree = tsw.read_tree_newick("(a,b,c,d,e);")
nodes = {}
for n in tree.root.child_nodes():
    nodes[n.label] = n
D = parse_d([nodes[l] for l in "abcde"], d)
S = NState(D)
print([n.label for n in S.find_closest()])