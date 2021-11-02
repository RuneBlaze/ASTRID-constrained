Constrained NJst
===================

Constrained NJst (NJst-J) similar to that of ASTRAL https://doi.org/10.1186/s12864-020-6607-z

## Command

```
python3 njst_constrained.py -i $genes -j $constraint -o $output
```

## Dependencies

Again, NJst-J currently only works on Linux (x86_64 architecture) due
to its dependency on `asterid` (https://pypi.org/project/asterid/), a Python3 binding for [ASTRID](https://github.com/pranjalv123/ASTRID).

NJst constrained requires the following packages from pip:

```python3
asterid # currently for Linux x86_64 only
icecream
lupa
```

Additionally, NJst-J (through `lupa`) uses Lua for faster distance computations,
and in the original experiments LuaJIT 2.1.0-beta3 was used to speed up computation.
The `lupa` package thus should be built against LuaJIT for the full speed.