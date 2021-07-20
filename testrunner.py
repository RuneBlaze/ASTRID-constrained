# this is for testing the correctness of the neighbor-joining code
import asterid as ad
import dendropy
import dendropy as dp
import nj
import argparse
import treeswift as trs
from treecmp import compareTreesFromNewick

def consensus(trees, minfreq=0.5):
    import dendropy
    res = dendropy.TreeList()
    for treenewick in trees:
        res.read(data=treenewick, schema="newick", rooting='force-unrooted')
    con = res.consensus(min_freq = minfreq)
    con.is_rooted = False
    return con.as_string(schema="newick")

def test_nj(gtreepath):
    import math
    genes = open(gtreepath, "r").readlines()
    ts = ad.get_ts(genes)
    D = ad.mk_distance_matrix(ts, genes)
    tree = trs.read_tree_newick(consensus(genes, 2.0))
    tree.contract_low_support(1)
    print(tree)
    merged_tree = nj.treeresolve(tree, ts, D)
    res = merged_tree.newick()
    res2 = ad.fastme_nj(ts, D, 0, 0)
    print(res)
    print("")
    print(res2)
    nl, ei1, ei2, fp, fn, rf = compareTreesFromNewick(res, res2)
    # assert rf == 0, f"{gtreepath}: {nl, ei1, ei2, fp, fn, rf}"
    return rf

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gene', type=str, required=True)
args = parser.parse_args()

rf = test_nj(args.gene)
print(f"testing {args.gene}: {rf}")