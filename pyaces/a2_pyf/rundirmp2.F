      Subroutine Rundirmp2(Iuhf)
C
C Direct MBPT(2). The classic driver reaches this by running xdirmp2 as a
C separate program after xvscf (job_control.F), so there was no way to ask
C for it from here and suite cases 117a/b/c could not run at all.
C
C dirmp2.F's program body has been split into dirmp2_run (same
C wrapper-split as joda/vscf/vtran/vprops), which is everything between
C runtime setup and teardown: it opens the OUT./DUMP. units its own
C routines write to, calls CRAPSI for the core and the JOBARC session, and
C dispatches on the reference. Rank 0 is passed because in process there is
C only ever the one rank.
C
C Iuhf is accepted for signature consistency with the other Run* entry
C points; dirmp2 reads the reference from JOBARC through CRAPSI itself, so
C it is deliberately not forwarded.
C
      Implicit None
      Integer*8 Iuhf
      Integer Length
      Double Precision Tmp
      Logical Scf_done
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "jobarc.com"
C
C Direct MP2 reads the SCF vectors, so run the SCF first if this process
C has not already done one -- same test Runcc uses. dirmp2 works straight
C off the AO integrals and the SCF vectors; there is no Vtran/Intproc
C stage, which is the whole point of the direct method.
C
      Call Getrec(0,"JOBARC","SCFEVCA0",Length,Tmp)
      Scf_done = .False.
      If (Length .Gt. 0) Scf_done = .True.
      If (.Not. Scf_done) Then
         Call Ints_(Icore(i0),Icrsiz,Iuhf)
         Call V2ja(icore(i0),Icrsiz,Iuhf)
         Call Vscf(Icore(i0),Icrsiz,Iuhf)
      End If
      Call Dirmp2_run(0)
      Return
      End
