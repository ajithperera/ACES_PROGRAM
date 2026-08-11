      Subroutine Runee(Iuhf)

      Implicit None
      Integer Iuhf

#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "jobarc.com"

      Call Joda_()
      Call Aces_ja_fin()
      Call aces_init(icore,i0,icrsiz, iuhf, .true.)

      Call Ints_(Icore(i0),Icrsiz,Iuhf)
      Call V2ja(icore(i0),Icrsiz,Iuhf)

C Classic driver's energy.f runs xvprops here (after xvmol2ja, before
C xvscf) whenever IFLAGS(87) (the EXCITE flag) is set -- Runee is only
C ever invoked for EXCITE= cases, so this is unconditional here. vprops's
C own PROPS(S,MAXCOR,IS,IMAXCOR) computes/stores the one-electron property
C integrals (DIPOLE_X/Y/Z etc, see vprops/dintfm.f) that vee's own EOM
C transition-property code (ESTATE_PROP=EXPECTATION) later expects to
C already exist via GETREC. vprops/props.F's own SUBROUTINE PROPS collides
C with props/props.F's differently-shaped SUBROUTINE PROPS(ICORE,MAXCOR,
C IUHF) (needed by Runprops) -- build.sh renames vprops's copy (and this
C call site's own reference to it) via objcopy --redefine-sym, same
C mechanism as the FINISH/GSCHMIDT/MODF/EXPDEN collisions. (Previously
C crashed here due to a separate ICRSIZ-corruption bug in aces_init.F,
C fixed 2026-08-11 -- see acescore/aces_init.F's iMemActual comment.)
      Call Props(Icore(i0),Icrsiz/Iintfp,Icore(i0),Icrsiz)

      Call Vscf(Icore(i0),Icrsiz,Iuhf)
      Call aces_init_chemsys()
      Call Vtran(icore(i0),Icrsiz,Iuhf)
      Call aces_init_chemsys()
      Call Init_legacies()
      Call Intproc(icore(i0),Icrsiz,Iuhf)
      Call Vcc(icore(i0),Icrsiz,Iuhf)
      Call Vlambda(icore(i0),Icrsiz,Iuhf)
      Call Vee(icore(i0),Icrsiz,Iuhf)

C Do not call aces_fin (use might be intending to run something
C else after external manupilations).

CSSS      Call aces_fin
CSSS      Call C_free(Icore)
 
      Return
      End

