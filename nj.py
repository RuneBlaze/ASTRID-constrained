from collections import defaultdict
from icecream import ic
from math import fsum

class NState:
    def __init__(self, D):
        self.D = D
    def find_closest(self):
        Q = defaultdict(dict)
        R = {}
        for i in self.D:
            for j in self.D:
                if i in Q[j]:
                    Q[i][j] = Q[j][i]
                    continue
                if i.get_parent() != j.get_parent():
                    continue
                if i not in R:
                    R[i] = fsum(self.D[i][k] for k in self.D)
                if j not in R:
                    R[j] = fsum(self.D[j][k] for k in self.D)
                # if i.label == 'a' and j.label == 'b':
                    # print("AABB")
                    # print(len(self.D),self.D[i][j])
                    #ic([self.D[i][k] for k in self.D])
                    #ic([self.D[j][k] for k in self.D])
                    #ic((len(self.D) - 2) * self.D[i][j]
                # - sum(self.D[i][k] for k in self.D)
                # - sum(self.D[j][k] for k in self.D))
                #ic(self.D[i][j])
                #ic((len(self.D) - 2) * self.D[i][j]
                # - sum(self.D[i][k] for k in self.D)
                # - sum(self.D[j][k] for k in self.D))
                #ic(self.D[i][j])
                # val = (len(self.D) - 2) * self.D[i][j] - sum(self.D[i][k] for k in self.D) - sum(self.D[j][k] for k in self.D)
                #ic(self.D[i][j])
                #ic(val)
                #ic(len(self.D))
                Q[i][j] = (len(self.D) - 2) * self.D[i][j] - R[i] - R[j]
                # print(f"{(i.label, j.label)}: {Q[i][j]}")
        # for u in Q:
        #     for v in Q:
        #         print(f"{(u.label, v.label)}: {Q[u][v]}")
        mindis = 121231234
        minpair = None

        for a, i in enumerate(self.D):
            for b, j in enumerate(self.D):
                if i == j:
                    continue
                if i.get_parent() != j.get_parent():
                    continue
                if Q[i][j] < mindis:
                    mindis = Q[i][j]
                    minpair = (i, j)
        return minpair

    def join(self, u, v, n):
        # print([u, v, n])
        # print(f"joining {u.newick()} with {v.newick()} to {n.newick()}")
        # print(u.newick(), v.newick())
        # d = lambda i, j: self.D[i][j]
        def d(i, j):
            if i not in self.D:
                print(i.newick())
            if j not in self.D:
                print(j.newick())
            return self.D[i][j]
        for k in self.D:
            self.D[k][n] = 0.5 * (d(u, k) + d(v, k) - d(u, v))
        self.D[n] = {}
        for k in self.D:
            if k == n:
                self.D[n][k] = 0
            self.D[n][k] = self.D[k][n]
        del self.D[u]
        del self.D[v]

import treeswift as ts
import treeswift as tsf

def treeresolve(tree, ts, D):
    # print(D.str())
    dis = defaultdict(dict)
    for i in tree.traverse_postorder(True, False):
        for j in tree.traverse_postorder(True, False):
            dis[i][j] = D[ts[i.label],ts[j.label]]
    state = NState(dis)
    tree.suppress_unifurcations()
    while len(state.D) > 1:
        i, j = state.find_closest()
        if i.get_parent().num_children() == 2:
            state.join(i, j, i.get_parent())
            continue
        u = i.get_parent()
        u.remove_child(i)
        u.remove_child(j)
        nn = tsf.Node(edge_length=0)
        u.add_child(nn)
        nn.add_child(i)
        nn.add_child(j)
        state.join(i, j, nn)
        # print("")
        # print(tree.newick())
        # print([k.newick() for k in state.D])
        # print([k.get_parent().newick() for k in state.D])
        # print(state.D)
        
    # for u in tree.traverse_postorder():
    #     if u.is_leaf():
    #         continue
    #     else:
    #         while len(u.child_nodes()) > 2:
    #             i, j = state.find_closest(u.child_nodes())
    #             c2 = u.children.pop(j)
    #             c1 = u.children.pop(i)
    #             nn = tsf.Node(edge_length=0)
    #             u.add_child(nn)
    #             nn.add_child(c1)
    #             nn.add_child(c2)
    #             state.join(c1, c2, nn)
    #             # print(nn.newick())
    #         if len(u.child_nodes()) == 2: # binary
    #             # print("binary tree!")
    #             state.join(*u.child_nodes(), u)
    return tree


# def whitelist_tree(tree):
#     whitelist = {}
#     for u in tree.traverse_postorder():
#         if u.is_leaf():
#             s = frozenset([u.label])
#             u.cladeset = s
#             whitelist[s] = u
#         else:
#             s = frozenset().union(*[c.cladeset for c in u.child_nodes()])
#             whitelist[s] = u
#             u.cladeset = s
#             u.allowedset = set(c.cladeset for c in u.child_nodes())
#     return whitelist

constraint = "((1,2,3),(4,5));"
constrainttree = ts.read_tree_newick(constraint)


if __name__ == "__main__":
    import asterid as ad
    trees = ["((1,2),(3,(4,5)));"]
    ts = ad.get_ts(trees)
    D = ad.mk_distance_matrix(ts, trees)
    res = treeresolve(constrainttree, ts, D)
    print(res.newick())