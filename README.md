Constrained NJst
===================

NJst-constrained (NJst-J) is a summary method based on [NJst](https://academic.oup.com/sysbio/article/60/5/661/1644054)
that can honor user constraints (analogous to [ASTRAL with user constraints](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6607-z), although this work uses a different approach). This uploaded version works only for x86_64 Linux machines.

This branch contains non-widely tested code for dealing with incomplete constraints
(constraint trees with some taxa missing).

## Command

```
python3 njst_constrained.py -i $genes -j $constraint -o $output
```

 - `$genes`: path to newline-separated Newick formatted single-label gene trees
 - `$constraint`: the user-specified unresolved constraint tree
 - `$output`: output species-tree path

NJst-J currently outputs a tree with only the tree topology. It is likely that ASTRAL should be used for estimating branch lengths and support values.

## Dependencies

Again, NJst-J currently only works on Linux (x86_64 architecture) due
to its dependency on `asterid` (https://pypi.org/project/asterid/), a Python3 binding for [ASTRID](https://github.com/pranjalv123/ASTRID).

NJst constrained requires the following packages from pip:

```python3
treeswift # we used v1.1.19, although latest stable version likely also works
asterid # currently for Linux x86_64 only
icecream
lupa # Lua bindings for Python, it comes with a default Lua prepackaged
```

Additionally, NJst-J (through `lupa`) uses Lua for faster distance computations,
and in the original experiments LuaJIT 2.1.0-beta3 was used to speed up computation.
The `lupa` package thus should be built against LuaJIT for the full speed.