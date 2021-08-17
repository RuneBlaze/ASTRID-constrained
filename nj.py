from collections import defaultdict
from icecream import ic
from math import fsum, inf
from itertools import combinations
# import numpy as np
# import numpy.ma as ma

LUASRC = """
D = {}
parent = {}

function set_distance(d)
    D = d
end

function set_parent(p)
    parent = p
end

function find_closest()
    local mindis = 1231231234
    local minpair = nil
    local R = {}
    local N = 0
    for _, _ in pairs(D) do
        N = N + 1
    end
    local active = {}
    for i, _ in pairs(D) do
        table.insert(active, i)
    end
    for _i = 1,#active do
        for _j = _i+1,#active do
            local i = active[_i]
            local j = active[_j]
            local ip = parent[i]
            local jp = parent[j]
            if ip == jp then
                if R[i] == nil then
                    R[i] = 0
                    for k, _ in pairs(D) do
                        R[i] = R[i] + D[i][k]
                    end
                end

                if R[j] == nil then
                    R[j] = 0
                    for k, _ in pairs(D) do
                        R[j] = R[j] + D[j][k]
                    end
                end
                local qij = (N - 2) * D[i][j] - R[i] - R[j]
                if qij < mindis then
                    mindis = qij
                    minpair = {i, j}
                end
            end
        end
    end
    return minpair[1], minpair[2]
end

function join_node(u, v, n)
    for k, _ in pairs(D) do
        D[k][n] = 0.5 * (D[u][k] + D[v][k] - D[u][v])
    end
    D[n] = {}
    for k, _ in pairs(D) do
        if k == n then
            D[n][k] = 0
        end
        D[n][k] = D[k][n]
    end
    -- update parent
    -- if not joined against its direct parent,
    -- then the parent of the new node is n
    if parent[u] ~= n then
        parent[n] = parent[u]
    end
    -- if joined against its direct parent
    -- then no need to update
    D[u] = nil
    D[v] = nil
    local cnt = 0
    for i, _ in pairs(D) do
        cnt = cnt + 1
    end
    return cnt
end
"""

class NMState:
    def __init__(self, D, l):
        self.D = D
        self.len = l
        self.rest = l

    def find_closest(self, id2node):
        Q = (self.len - 2.0) * self.D
        for i in range(min(Q.shape[0], self.len)):
            for j in range(min(Q.shape[1], self.len)):
                if id2node[i].get_parent() != id2node[j].get_parent():
                    Q[i, j] = ma.masked
                    continue
                if i == j:
                    Q[i, j] = ma.masked
                    continue
                if i != j:
                    Q[i, j] -= ma.sum(self.D[i]) + ma.sum(self.D[j])
        x, y = np.unravel_index(Q.argmin(fill_value = 1231231234), Q.shape)
        ic(x, y, Q[x, y])
        return x, y

    def join(self, u, v, n):
        D = self.D
        for slices in ma.notmasked_contiguous(self.D, axis=0):
            for k in slices:
                D[k,n] = 0.5 * (D[u, k] + D[v, k] - D[u, v])
                D[n, k] = D[k, n]
        D[n, n] = 0
        # D[u].mask = True
        # D[v].mask = True
        # D[:, u].mask = True
        # D[:, v].mask = True
        D[u, v] = ma.masked
        ma.mask_rowcols(D)
        self.rest -= 1

def treeresolve_fast(tree, ts, D):
    # dis = defaultdict(dict)
    # for i in tree.traverse_postorder(True, False):
    #     for j in tree.traverse_postorder(True, False):
    label2id = {}
    id2node = {}
    for i in ts:
        label2id[ts[i]] = i
    for n in tree.traverse_postorder(True, False):
        id2node[label2id[n.label]] = n
    # ic(id2node)
    nd = ma.zeros((len(ts) * 2, len(ts) * 2))
    nd.harden_mask()
    for i in ts:
        for j in ts:
            nd[i,j] = D[i, j]
    tree.suppress_unifurcations()
    state = NMState(nd, len(ts))
    while state.rest > 0:
        i, j = ic(state.find_closest(id2node))
        ix, jx = id2node[i], id2node[j]
        u = ix.get_parent()
        u.remove_child(ix)
        u.remove_child(jx)
        nn = tsf.Node(edge_length=0)
        u.add_child(nn)
        nn.add_child(ix)
        nn.add_child(jx)
        id2node[state.len] = nn
        state.join(i, j, state.len)
        state.len += 1
    return tree

import lupa
class LuaNState:
    def __init__(self, ts, AD, T):
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
        # ic(r)
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
    import asterid as ad
    trees = ["((1,2),(3,(4,5)));"]
    ts = ad.get_ts(trees)
    D = ad.mk_distance_matrix(ts, trees)
    # LuaNState(ts, D, )
    # res = treeresolve_fast(constrainttree, ts, D)
    # print(res.newick())%
λ ~/scratch/atsume/methods/ASTRID-constrained/ main* git status
On branch main
Your branch is up to date with 'origin/main'.

Changes to be committed:
  (use "git reset HEAD <file>..." to unstage)

	modified:   astrid_constrained.py
	modified:   nj.py

Changes not staged for commit:
  (use "git add/rm <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	deleted:    README.md

Untracked files:
  (use "git add <file>..." to include in what will be committed)

	dev/
	luajit/
	nj_.py
	out.tre
	out.tre.fn
	out1000.tre
	out2.tre
	scratch/
	test/c1000.tre
	test/t1000.tre
	tests/

λ ~/scratch/atsume/methods/ASTRID-constrained/ main* git commit -m "luajit"
error: invalid object 100644 e5bee14827ea9b361defe485a5abb235a7b1852b for 'README.md'
error: invalid object 100644 e5bee14827ea9b361defe485a5abb235a7b1852b for 'README.md'
error: Error building trees
λ ~/scratch/atsume/methods/ASTRID-constrained/ main* git checkout HEAD README.md
error: unable to read sha1 file of README.md (e5bee14827ea9b361defe485a5abb235a7b1852b)
λ ~/scratch/atsume/methods/ASTRID-constrained/ main* git checkout HEAD README.md
error: unable to read sha1 file of README.md (e5bee14827ea9b361defe485a5abb235a7b1852b)
λ ~/scratch/atsume/methods/ASTRID-constrained/ main* ls
__pycache__            dev     nj.py   out.tre     out1000.tre  scratch  testnj.py      tests
astrid_constrained.py  luajit  nj_.py  out.tre.fn  out2.tre     test     testrunner.py  treecmp.py
λ ~/scratch/atsume/methods/ASTRID-constrained/ main* cat nj.py
from collections import defaultdict
from icecream import ic
from math import fsum, inf
from itertools import combinations
# import numpy as np
# import numpy.ma as ma

LUASRC = """
D = {}
parent = {}

function set_distance(d)
    D = d
end

function set_parent(p)
    parent = p
end

function find_closest()
    local mindis = 1231231234
    local minpair = nil
    local R = {}
    local N = 0
    for _, _ in pairs(D) do
        N = N + 1
    end
    local active = {}
    for i, _ in pairs(D) do
        table.insert(active, i)
    end
    for _i = 1,#active do
        for _j = _i+1,#active do
            local i = active[_i]
            local j = active[_j]
            local ip = parent[i]
            local jp = parent[j]
            if ip == jp then
                if R[i] == nil then
                    R[i] = 0
                    for k, _ in pairs(D) do
                        R[i] = R[i] + D[i][k]
                    end
                end

                if R[j] == nil then
                    R[j] = 0
                    for k, _ in pairs(D) do
                        R[j] = R[j] + D[j][k]
                    end
                end
                local qij = (N - 2) * D[i][j] - R[i] - R[j]
                if qij < mindis then
                    mindis = qij
                    minpair = {i, j}
                end
            end
        end
    end
    return minpair[1], minpair[2]
end

function join_node(u, v, n)
    for k, _ in pairs(D) do
        D[k][n] = 0.5 * (D[u][k] + D[v][k] - D[u][v])
    end
    D[n] = {}
    for k, _ in pairs(D) do
        if k == n then
            D[n][k] = 0
        end
        D[n][k] = D[k][n]
    end
    -- update parent
    -- if not joined against its direct parent,
    -- then the parent of the new node is n
    if parent[u] ~= n then
        parent[n] = parent[u]
    end
    -- if joined against its direct parent
    -- then no need to update
    D[u] = nil
    D[v] = nil
    local cnt = 0
    for i, _ in pairs(D) do
        cnt = cnt + 1
    end
    return cnt
end
"""

class NMState:
    def __init__(self, D, l):
        self.D = D
        self.len = l
        self.rest = l

    def find_closest(self, id2node):
        Q = (self.len - 2.0) * self.D
        for i in range(min(Q.shape[0], self.len)):
            for j in range(min(Q.shape[1], self.len)):
                if id2node[i].get_parent() != id2node[j].get_parent():
                    Q[i, j] = ma.masked
                    continue
                if i == j:
                    Q[i, j] = ma.masked
                    continue
                if i != j:
                    Q[i, j] -= ma.sum(self.D[i]) + ma.sum(self.D[j])
        x, y = np.unravel_index(Q.argmin(fill_value = 1231231234), Q.shape)
        ic(x, y, Q[x, y])
        return x, y

    def join(self, u, v, n):
        D = self.D
        for slices in ma.notmasked_contiguous(self.D, axis=0):
            for k in slices:
                D[k,n] = 0.5 * (D[u, k] + D[v, k] - D[u, v])
                D[n, k] = D[k, n]
        D[n, n] = 0
        # D[u].mask = True
        # D[v].mask = True
        # D[:, u].mask = True
        # D[:, v].mask = True
        D[u, v] = ma.masked
        ma.mask_rowcols(D)
        self.rest -= 1

def treeresolve_fast(tree, ts, D):
    # dis = defaultdict(dict)
    # for i in tree.traverse_postorder(True, False):
    #     for j in tree.traverse_postorder(True, False):
    label2id = {}
    id2node = {}
    for i in ts:
        label2id[ts[i]] = i
    for n in tree.traverse_postorder(True, False):
        id2node[label2id[n.label]] = n
    # ic(id2node)
    nd = ma.zeros((len(ts) * 2, len(ts) * 2))
    nd.harden_mask()
    for i in ts:
        for j in ts:
            nd[i,j] = D[i, j]
    tree.suppress_unifurcations()
    state = NMState(nd, len(ts))
    while state.rest > 0:
        i, j = ic(state.find_closest(id2node))
        ix, jx = id2node[i], id2node[j]
        u = ix.get_parent()
        u.remove_child(ix)
        u.remove_child(jx)
        nn = tsf.Node(edge_length=0)
        u.add_child(nn)
        nn.add_child(ix)
        nn.add_child(jx)
        id2node[state.len] = nn
        state.join(i, j, state.len)
        state.len += 1
    return tree

import lupa
class LuaNState:
    def __init__(self, ts, AD, T):
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
        # ic(r)
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
    import asterid as ad
    trees = ["((1,2),(3,(4,5)));"]
    ts = ad.get_ts(trees)
    D = ad.mk_distance_matrix(ts, trees)
    # LuaNState(ts, D, )
    # res = treeresolve_fast(constrainttree, ts, D)
    # print(res.newick())