import asterid as ad
import dendropy
import dendropy as dp
import nj
import treeswift as ts
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run ASTRID')
    parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output tree file", required=True)
    parser.add_argument("-t", "--tree", type=str, help="constraint trees", required=True)

    args = parser.parse_args()

    with open(args.tree) as fh:
        tree = ts.read_tree_newick(fh.read())
    
    genes = open(args.input, "r").readlines()
    ts = ad.get_ts(genes)
    D = ad.mk_distance_matrix(ts, genes)
    merged_tree = nj.treeresolve(tree, ts, D)
    res = merged_tree.newick()
    if args.output == "-":
        print(res)
    else:
        with open(args.output, "w+") as outf:
            outf.write(res)