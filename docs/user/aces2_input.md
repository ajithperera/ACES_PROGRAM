# ACES II input

```{note}
This chapter is being written. It will be based on chapters 6 and 7 of the
ACES II 2.5.0 user manual (2006), checked against the current code.
```

Planned sections:

- The ZMAT file: title line, geometry, variables, the `*ACES2` namelist, and
  per-atom basis lines
- Internal (Z-matrix) and Cartesian coordinates, dummy and ghost atoms,
  molecular orientation, and the Z-matrix analyzer
- GENBAS and ZMAT.BAS basis set libraries; ECPDATA
- GUESS, OLDMOS and OLDAOMOS orbital guesses
- Output and restart files

(zmat-namelist)=
## The `*ACES2` namelist

Keywords are given after a line starting with `*ACES2`, either in
parentheses on that line or on the lines that follow, separated by commas or
new lines. A line starting with `#` is a comment. Every keyword is listed in
{ref}`keyword-reference`.
