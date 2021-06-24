import asterid as ad
import dendropy
import dendropy as dp
import nj
import treeswift as ts
import argparse

def adm_to_dendropy_pdm(dmat, ts):
    pdm = dendropy.PhylogeneticDistanceMatrix()
    pdm.taxon_namespace = dendropy.TaxonNamespace()
    pdm._mapped_taxa = set()

    for i in range(len(ts)):
        for j in range(len(ts)):
            si = ts[i]
            sj = ts[j]
            dij = dmat[i, j]
            xi = pdm.taxon_namespace.get_taxon(si)
            if not xi:
                xi = dendropy.Taxon(si)
                pdm.taxon_namespace.add_taxon(xi)
                pdm._mapped_taxa.add(xi)
                pdm._taxon_phylogenetic_distances[xi] = {}

            xj = pdm.taxon_namespace.get_taxon(sj)
            if not xj:
                xj = dendropy.Taxon(sj)
                pdm.taxon_namespace.add_taxon(xj)
                pdm._mapped_taxa.add(xj)
                pdm._taxon_phylogenetic_distances[xj] = {}
            dij = float(dij)
            pdm._taxon_phylogenetic_distances[xi][xj] = dij
    return pdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run ASTRID')
    parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output tree file", required=True)
    parser.add_argument("-t", "--tree", type=str, help="constraint trees", required=True)

    args = parser.parse_args()

    # trees = []
    # trees.append(dendropy.Tree.get(path=args.tree, schema="newick"))
    with open(args.tree) as fh:
        tree = ts.read_tree_newick(fh.read())
    
    
    genes = open(args.input, "r").readlines()
    # methods = args.methods
    # print(genes)
    ts = ad.get_ts(genes)
    D = ad.mk_distance_matrix(ts, genes)
    merged_tree = nj.treeresolve(tree, ts, D)
    res = merged_tree.newick()
    if args.output == "-":
        print(res)
    else:
        with open(args.output, "w+") as outf:
            outf.write()