from collections import defaultdict
from icecream import ic
from math import fsum, inf
from itertools import combinations
import os
import treeswift as ts
import treeswift as tsf

luasrc = os.path.join(os.path.dirname(os.path.realpath(__file__)), "njext.lua")
with open(luasrc, "r") as fh:
    LUASRC = fh.read()

class LuaNState:
    def __init__(self, ts, AD, T, rogue):
        import lupa
        lua = lupa.LuaRuntime(unpack_returned_tuples=True)
        self.lua = lua
        self.lua.execute(LUASRC)
        self.rogue = rogue
        self.id2node = {}
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
        lua.globals().rogue = lua.table_from({id(i): True for i in rogue})
        self.parent = lua.globals().parent
        self.join_node_lua = self.lua.eval("join_node")
    def find_closest(self):
        ii, ij = self.lua.eval("find_closest()")
        return self.id2node[ii], self.id2node[ij]
    def join(self, u, v, n):
        self.id2node[id(n)] = n
        r = self.join_node_lua(id(u), id(v), id(n))
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

def degree_of_resolution(tre):
    tt = 0
    resolved = 0
    for n in tre.traverse_preorder(False, True):
        tt += 1
        resolved += n.num_children() == 2
    return resolved / tt

def treeresolve_lua(tree, ts, D, rogue):
    state = LuaNState(ts, D, tree, rogue)
    rogueset = set(r.label for r in rogue)
    while True:
        i, j = state.find_closest()
        ir = i.label in rogueset
        jr = j.label in rogueset
        if ir == jr:
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
        else:
            def unbalanced_merge(l, r): # merge rogue r to non-rogue l
                u = l.get_parent()
                u.remove_child(l)
                nn = tsf.Node(edge_length=0)
                u.add_child(nn)
                r.get_parent().remove_child(r)
                nn.add_child(l)
                nn.add_child(r)
                return state.join(l, r, nn)
            if jr: # j is rogue
                r = unbalanced_merge(i, j)
            else:
                r = unbalanced_merge(j, i)
        if r <= 2:
            break
    return tree

def treeresolve(tree, ts, D):
    dis = defaultdict(dict)
    leaves = list(tree.traverse_postorder(True, False))
    cache = {}
    for i in leaves:
        dis[i][i] = 0
    for i in ts:
        cache[ts[i]] = i
    for i, j in combinations(leaves, 2):
        dis[i][j] = D[cache[i.label],cache[j.label]]
        dis[j][i] = dis[i][j]
    state = NState(dis)
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
    return tree