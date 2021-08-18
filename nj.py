from collections import defaultdict
from icecream import ic
from math import fsum, inf
from itertools import combinations
import os
# import numpy as np
# import numpy.ma as ma
luasrc = os.path.join(os.path.dirname(os.path.realpath(__file__)), "njext.lua")
with open(luasrc, "r") as fh:
    LUASRC = fh.read()

class LuaNState:
    def __init__(self, ts, AD, T):
        import lupa
        lua = lupa.LuaRuntime(unpack_returned_tuples=True)
        self.lua = lua
        self.lua.execute(LUASRC)
        self.id2node = {}
        # T.suppress_unifurcations()
        for n in T.traverse_postorder(True, True):
            self.id2node[id(n)] = n
        parent = {}
        dis = {}
        leaves = list(T.traverse_postorder(True, False))
        cache = {}
        for i in leaves:
            ii = id(i)
            dis[ii] = {}
            dis[ii][ii] = 0
            if i.is_root():
                parent[ii] = None
            else:
                parent[ii] = id(i.get_parent())
        for i in ts:
            cache[ts[i]] = i
        for i, j in combinations(leaves, 2):
            ii, ij = id(i), id(j)
            dis[ii][ij] = AD[cache[i.label],cache[j.label]]
            dis[ij][ii] = dis[ii][ij]
        lua.globals().D = lua.table_from(dis)
        lua.globals().parent = lua.table_from(parent)
        self.parent = lua.globals().parent
        self.join_node_lua = self.lua.eval("join_node")
    def find_closest(self):
        ii, ij = self.lua.eval("find_closest()")
        return self.id2node[ii], self.id2node[ij]
    def join(self, u, v, n):
        self.id2node[id(n)] = n
        r = self.join_node_lua(id(u), id(v), id(n))
        # p = self.lua.globals().parent
        self.parent[id(n)] = id(n.get_parent())
        return r

class NState:
    def __init__(self, D, earlystopping=False):
        self.D = D
        self.earlystopping = earlystopping
    def find_closest(self):
        mindis = 1231231234
        minpair = None
        R = {}
        N = len(self.D)
        for i, j in combinations(self.D, 2):
            ip = i.get_parent()
            jp = j.get_parent()
            if ip != jp:
                continue
            if i not in R:
                R[i] = fsum(self.D[i][k] for k in self.D)
            if j not in R:
                R[j] = fsum(self.D[j][k] for k in self.D)
            qij = (N - 2) * self.D[i][j] - R[i] - R[j]
            if qij < mindis:
                mindis = qij
                minpair = (i, j)
        return minpair

        # for i in self.D:
        #     for j in self.D:
        #         if i in Q[j]:
        #             Q[i][j] = Q[j][i]
        #             continue
        #         ip = i.get_parent()
        #         jp = j.get_parent()
        #         if ip != jp:
        #             continue
        #         if self.earlystopping and ip.num_children() == 2:
        #             return (i, j)
        #         if i not in R:
        #             R[i] = fsum(self.D[i][k] for k in self.D)
        #         if j not in R:
        #             R[j] = fsum(self.D[j][k] for k in self.D)
        #         Q[i][j] = (len(self.D) - 2) * self.D[i][j] - R[i] - R[j]
        #         if Q[i][j] < mindis:
        #             mindis = Q[i][j]
        #             minpair = (i, j)
        # mindis = 121231234
        # minpair = None

        # for a, i in enumerate(self.D):
        #     for b, j in enumerate(self.D):
        #         if i == j:
        #             continue
        #         if i.get_parent() != j.get_parent():
        #             continue
        #         if Q[i][j] < mindis:
        #             mindis = Q[i][j]
        #             minpair = (i, j)
        # return minpair

    def join(self, u, v, n):
        def d(i, j):
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

def degree_of_resolution(tre):
    tt = 0
    resolved = 0
    for n in tre.traverse_preorder(False, True):
        tt += 1
        resolved += n.num_children() == 2
    return resolved / tt

def treeresolve_lua(tree, ts, D):
    state = LuaNState(ts, D, tree)
    while True:
        # ic(degree_of_resolution(tree))
        i, j = state.find_closest()
        # ic(i.newick(), i.get_parent().num_children())
        # ic(j.newick())
        if i.get_parent().num_children() == 2:
            cnt = state.join(i, j, i.get_parent())
            if cnt <= 2:
                break
            continue
        u = i.get_parent()
        u.remove_child(i)
        u.remove_child(j)
        nn = tsf.Node(edge_length=0)
        u.add_child(nn)
        nn.add_child(i)
        nn.add_child(j)
        r = state.join(i, j, nn)
        if r <= 2:
            break

    return tree

def treeresolve(tree, ts, D):
    dis = defaultdict(dict)
    leaves = list(tree.traverse_postorder(True, False))
    cache = {}
    for i in leaves:
        # cache[i.label] = ts[i.label]
        dis[i][i] = 0
    for i in ts:
        cache[ts[i]] = i
    for i, j in combinations(leaves, 2):
        dis[i][j] = D[cache[i.label],cache[j.label]]
        dis[j][i] = dis[i][j]
    # for i in leaves:
    #     for j in leaves:
    #         dis[i][j] = D[ts[i.label],ts[j.label]]
    state = NState(dis)
    # tree.suppress_unifurcations()
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
    # dis = defaultdict(dict)
    # label2id = {}
    # leaves = tree.traverse_postorder(True, False)
    # for i in ts:
    #     label2id[ts[i]] = i
    # for l in leaves:
    #     label2id[l] = label2id[l.label]
    # for i, j in combinations(leaves, 2):
    #     dis[i][j] = D[label2id[i], label2id[j]]
    #     dis[j][i] = dis[i][j]
    #     dis[i][i] = 0
    # state = NState(dis)
    # tree.suppress_unifurcations()
    # while len(state.D) > 1:
    #     i, j = state.find_closest()
    #     if i.get_parent().num_children() == 2:
    #         state.join(i, j, i.get_parent())
    #         continue
    #     u = i.get_parent()
    #     u.remove_child(i)
    #     u.remove_child(j)
    #     nn = tsf.Node(edge_length=0)
    #     u.add_child(nn)
    #     nn.add_child(i)
    #     nn.add_child(j)
    #     state.join(i, j, nn)
    return tree

if __name__ == "__main__":
    print(LUASRC)
    # import asterid as ad
    # trees = ["((1,2),(3,(4,5)));"]
    # ts = ad.get_ts(trees)
    # D = ad.mk_distance_matrix(ts, trees)
    # LuaNState(ts, D, )
    # res = treeresolve_fast(constrainttree, ts, D)
    # print(res.newick())%