# Duplicate symbols inside one archive

Some aces2 archives (`aces2/lib/lib*.a`) have two or more members that
define the same routine. The link does not complain: it takes one copy and
ignores the others. Which copy it takes depends on the order of the members
in the archive and on which object first needed the routine. When the
copies differ, a binary can silently call the wrong one.

This page lists every such case found on 2026-10-28 (after the upstream
merge), says which copy each binary actually gets today, and suggests a fix
for each. Use it as a checklist: fix an item, rerun the check, and update
the table.

## Why it matters

GNU ld reads an archive in one pass from front to back. It pulls in a
member only when that member defines a symbol that is undefined at that
moment. If the symbol first becomes undefined while the pass is already
past the first definition (because the member that calls it was itself
pulled in later in the pass), ld takes the next member that defines it.

That is how `pccd_oo_rfo_ccd_h2o` and `pccd_oo_rfo_lccd_h2o` failed on the
local workstation. `pccd/pccd_read_htau.F` defined a second
`Pccd_form_htau` with three arguments, unused, next to the real one in
`pccd/pccd_form_htau.F`. The call in `pccd_form_htau_2d.F` was bound to the
three-argument copy, which read array data as `Maxcor` and stopped with a
bogus "Insufficient memory". On HiPerGator the member order happened to pick
the right copy. Fixed in `0a27e717` by renaming the unused copy to
`Pccd_read_htau`.

## How to check

```
source build_config.sh; aces_program_load_env
export ACES_PROGRAM_MACHNAME
python3 aces2/Scripts/dup_symbols.py            # list duplicates
python3 aces2/Scripts/dup_symbols.py --trace    # also: which copy each binary links
python3 aces2/Scripts/dup_symbols.py --json     # machine-readable
```

It reads the archives in `aces2/lib`, so build aces2 first. For each
duplicate it prints every member with its position in the archive, the
routine header (arguments), whether the routine bodies are identical, and
the files that call it. `--trace` relinks every aces2 binary into a scratch
file with `-Wl,--trace-symbol` (it does not touch `aces2/bin`). On the
shared workstation, run it on the allocated cores
(`taskset -c 20-23 flock -s /tmp/qtoy-bench.lock ...`).

When adding a file to a module, check that its routine names are not
already defined elsewhere in the same module: `nm -A --defined-only
aces2/lib/lib<module>.a | grep ' T <name>_$'`.

## Summary

Risk: **high**, the copies have different arguments or meaning, so the
wrong one would corrupt memory or results. **Medium**, the copies have the
same arguments but different code. **Low**, the copies are identical, or
nothing links the symbol.

| # | Archive | Symbol | Copies | Risk | Linked today | Status |
|---|---|---|---|---|---|---|
| 1 | `libpccd.a` | `pccd_form_htau_` | 2, different arguments | high | wrong copy (local workstation) | **fixed** `0a27e717` |
| 2 | `liba3ecp.a` | `compute_2dernai_integrals_` | 2, different arguments | high | `compute_2dernai_integrals.o` (correct) | open |
| 3 | `libcopex.a` | `drccof_` | 5, 16 or 17 arguments | high | `drccof.o` | open |
| 4 | `libvee.a` | `tablee_` | 2, same arguments, different code | medium | `tablee.o` (newer) | open |
| 5 | `libsbtools.a` | `sb_mem_fin_` | 2, different COMMON blocks | medium | `ks_mem_fin.o` (matches the callers) | open |
| 6 | `libsbtools.a` | `sb_io_init_` | 2, different COMMON blocks | low | not linked | open |
| 7 | `libscrnc.a` | `MAIN__` | 2 main programs | low | `scrnc_main.o` (main object) | open |
| 8 | `libintgrt.a` | 7 symbols (`corlsd_`, `corpbe_`, `corpw91_`, `exchpbe_`, `exchpw91_`, `gcor_`, `gcor2_`) | 2, identical | low | the single-routine files | open |
| 9 | `liboeptools.a` | 7 symbols (`aainer_`, `civpt_`, `evalreplp_`, `evlnucatt_`, `fmc_`, `rhftce_`, `setrhf_`) | 2, identical | low | `nucatt.o` | open |
| 10 | `liboombpt2.a` | `show_rot_` | 2, identical | low | `oombpt_formdiis.o` | open |
| 11 | `libtools.a` | `to_combat_compiler_stupidity_` | 5, identical empty stub | low | not linked | leave |

"Linked today" comes from `--trace` on the local workstation (ifx 2024.0,
sorted member order). pyaces links the same archives into one `.so`; for
items 4 and 5 it gets the same copy as the aces2 binaries (checked by symbol
size). aces3 was not checked.

## Details and suggested fixes

### 2. `compute_2dernai_integrals_` in `liba3ecp.a` (high)

| Pos. | Member | Arguments end with |
|---|---|---|
| 7 | `a3ecp/compute_2dernai_integrals.F` | `out_nai, Hess, …, ERD_index, ERD_scale` |
| 8 | `a3ecp/compute_2dernai_pvp_integrals.F` | `out_pvp, …, ERD_index, ERD_scale, ACES_PVP` (one more) |

The only caller, `a3ecp/main.F:479`, uses the first form, and `xa3ecp` links
the first copy, because `main.o` needs the routine before the archive pass
starts. The `_pvp_` file looks like a copy that was meant to be renamed and
never was.

**Fix:** rename the routine in `compute_2dernai_pvp_integrals.F` to
`compute_2dernai_pvp_integrals` (it has no callers), or delete the file if
it is abandoned.

### 3. `drccof_` in `libcopex.a` (high)

| Pos. | Member | Arguments |
|---|---|---|
| 22 | `copex/drccof.f` | 16 (`…, LWORKEA, NEA, NEAN, NEASH, NIP`), comment "HE AUG-CC-PVTZ" |
| 23 | `copex/drccofh2pbs.f` | 16 (`…, NEA, NVIR, NEASH, NIP`) |
| 24 | `copex/drccofhe1.f` | 17 (extra `CPIN`) |
| 25 | `copex/drccofhe2.f` | 17 (extra `CPIN`) |
| 26 | `copex/drccofne.f` | 17 (extra `CPIN`) |

These are system-specific variants (He, Ne, H2) of the same routine.
`copex/copex.F:249,293` calls it with 16 arguments, and `xcopex` links
`drccof.o`. Which variant is used is decided by archive order, not by the
input.

**Fix:** decide with the copex author which variant is current. Either
keep that one as `DRCCOF` and move the others out of the build (for
example into `copex/archive/`), or give each variant its own name and
choose between them explicitly in `copex.F`.

### 4. `tablee_` in `libvee.a` (medium)

| Pos. | Member | Note |
|---|---|---|
| 351 | `vee/tablee.F` | newer: extra oscillator-strength arrays and `CORE_SEARCH`/`DOUBLE_EES` printing |
| 352 | `vee/tablee_bd.F` | older version of the same routine |

Callers: `vee/dobwpt2.F`, `doeomee_blockdavid.F`, `doeomee_david.F`,
`doeomee_lanczos.F`, `doparteom.F`. `xvee` and pyaces link `tablee.o`.
Only alphabetical order keeps it that way.

**Fix:** delete `tablee_bd.F` (and `.f`), or rename its routine if it is
still wanted for comparison.

### 5. `sb_mem_fin_` in `libsbtools.a` (medium)

| Pos. | Member | Uses |
|---|---|---|
| 7 | `sbtools/ks_mem_fin.F` | `kscore.com`, `ks_mem.com` |
| 45 | `sbtools/sb_mem_fin.F` | `sbcore.com`, `sb_mem.com` |

Callers: `intgrt/intgrt_run.F`, `intpack/intpack.F`, `vksdint/vksdint.F`.
All three allocate with `ks_mem_init`, so `ks_mem_fin.o`'s copy is the
right one, and it is the one linked. Both copies only print memory use
under `_DEBUG`, so a swap would print the wrong pool's numbers rather than
break anything.

**Fix:** rename the routine in `ks_mem_fin.F` to `ks_mem_fin` (matching
`ks_mem_init`) and update the three callers.

### 6. `sb_io_init_` in `libsbtools.a` (low)

`sbtools/ks_io_init.F` and `sbtools/sb_io_init.F` both define `sb_io_init`,
with different COMMON blocks. Nothing calls it (the bodies are mostly
`cTODO` lines), and no binary links it.

**Fix:** rename the copy in `ks_io_init.F` to `ks_io_init`, as for item 5.

### 7. `MAIN__` in `libscrnc.a` (low)

`scrnc/scrnc_main.F` (`program scrnc`) and `scrnc/scrnc_sigma_driver.F`
(`program Scrnc_sigma_driver`) are both main programs in the same archive.
`xscrnc` links `scrnc_main.o` explicitly as its main object, so the second
program is never used and has no executable of its own.

**Fix:** if `Scrnc_sigma_driver` is meant to run, build it as a separate
executable; otherwise keep it out of `libscrnc.a`.

### 8. Seven XC routines in `libintgrt.a` (low)

`intgrt/perdew_corr_lib.F` (pos. 52) and `intgrt/perdew_exch_lib.F` (pos.
53) contain the same routines as the single-routine files `corlsd.f`,
`corpbe.f`, `corpw91.f`, `exchpbe.f`, `exchpw91.f`, `gcor.f` and `gcor2.f`
(pos. 3–28). The bodies are identical, so whichever is linked is
harmless. `xintgrt`, `xvksdint` and `xvscf` link `gcor2.o`.

**Fix:** remove one set (the two `perdew_*_lib.F` files, or the seven
single files) so that a later edit to one copy cannot be silently ignored.

### 9. Seven integral routines in `liboeptools.a` (low)

`oeptools/gridnucatt.F` (pos. 12) and `oeptools/nucatt.F` (pos. 18) share
`AAINER`, `CIVPT`, `EVALREPLP`, `EVLNUCATT`, `FMC`, `RHFTCE` and `SETRHF`,
identical in both. `xvscf` links all of them from `nucatt.o`.

**Fix:** move the shared routines into one file and remove them from the
other.

### 10. `show_rot_` in `liboombpt2.a` (low)

Identical in `oombpt2/oombpt_formdiis.f` (pos. 29) and
`oombpt2/oombptdiis.f` (pos. 32). No caller found by name.

**Fix:** keep one copy.

### 11. `to_combat_compiler_stupidity_` in `libtools.a` (low)

An empty subroutine from `tools/return.f`, which `achar.F`, `blas_dp2sp.F`,
`iachar.F` and `ssum.F` `#include` so that the file always has
some compilable code. Identical everywhere and never called. Leave it.

## Not covered

- Symbols duplicated across *different* archives. pyaces handles those for
  its single `.so` with `pyaces/tools/localize.txt` and
  `collision_audit.py`; the aces2 binaries link only a handful of archives
  each.
- aces3 and the pyaces glue code.
- `xfsip` (main object `fs.o`) and `xrspdens1` (a second binary in
  `rspdens/`): `--trace` cannot get their link lines and skips them. By
  `ACES_LINK_LISTS` they link only `libr2`, `crust` and the default
  libraries, so of the archives above only `libtools.a` (item 11).
