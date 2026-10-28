# Adding to ACES II

## A new keyword

Keywords are defined in two places, which must agree:

1. **`aces2/include/flags.h`**: a handle for the new entry,
   `h_IOPPAR_<name>` (and `h_IFLAGS_<name>` or `h_IFLAGS2_<name>`), equal
   to its position in the table below. The first 100 entries are the
   `IFLAGS` block, the next ones `IFLAGS2`.
2. **`aces2/joda/asv_ctl.c`**: one line in the `ASV_nl[]` table,

   ```c
   {/*283*/ "INEWKEY", "NEW_K#EYWORD", h_ICHAR_handle, 0, ""},
   ```

   with an internal name, the keyword (`#` marks the end of the minimum
   abbreviation), the type (`h_ICHAR_handle`, `h_ICHAR_f_int`,
   `h_ICHAR_string`, ...), the default, and units. Reuse an empty slot
   (`{"","",0,0,""}`) where one exists, and keep the comment number equal to
   the position.
3. **For a handle** (a keyword with named values), add a `case` to
   `asv_handle_proc` in the same file with its list of names. **End the
   list with an empty string `""`**: the lookup scans until it finds the
   value or an empty string, so a missing terminator makes unknown values
   read past the end of the array. Mind the comma after the last name;
   `"LAST" ""` without a comma is one string. Keywords without their own
   `case` are switches (`ON`/`OFF`).

Programs read the value from `IFLAGS`/`IFLAGS2` with the handle, for
example `iflags2(h_IFLAGS2_newkey)`, after `#include "flags.h"`.

Rebuild and install `xjoda`, then check the parsing with a small ZMAT: an
invalid value must stop with "The value token is not recognized". Rebuild
the manual: the keyword then appears in the {ref}`keyword-reference`
automatically. Add a description under its name in
`docs/reference/keyword_notes.md`.

## A new member program

1. Make a directory under `aces2/` with your sources and a `GNUmakefile`
   link to `../Makefiles/GNUmakefile`, as the other program directories
   have.
2. Start the program with `aces_init` (or `crapsi`) and end it with
   `aces_fin`; read inputs with `GETREC`/`GETLST` and leave results in
   `JOBARC` (see {doc}`io`). Write the final energy to `TOTENERG` if the
   program changes it.
3. Add the call to the driver in `aces2/aces2/` (usually
   `job_control.F` or `post_scf_eneg.F`) with `call runit('xnewprog')`,
   under the flag that selects it.
4. Build with `gmake && gmake install`, and add a test case (see
   {doc}`testing`).

## Documentation

When you add a method or keyword, update the manual in the same commit:
the keyword notes, an entry in {doc}`../user/methods` with the test case as
its example, and a description of any new `JOBARC` record in
`docs/programmer/jobarc_notes.md`.
