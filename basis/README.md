# Basis sets and ECPs

`GENBAS` (basis sets) and `ECPDATA` (effective core potentials) are the
libraries that every test suite in ACES_PROGRAM uses: ACES II, ACES III and
pyaces. See the manual (ACES II input, "GENBAS and ZMAT.BAS" and
"ECPDATA") for the file formats.

- `GENBAS`: the former `aces2/tests/fast/GENBAS`, plus four entries that
  only `aces3/tests/fast/GENBAS` had (`S:AUG-CC-PV(Q+D)Z`,
  `S:AUG-CC-PV(5+D)Z`, `NE:D-AUG-CC-PCVDZ`, `O:6-31G(d,p)`). The two files
  agreed on every entry they shared.
- `ECPDATA`: from `KEEP/` of the old ACES II test archive.

Add new basis sets and ECPs here, not to copies in test directories.
