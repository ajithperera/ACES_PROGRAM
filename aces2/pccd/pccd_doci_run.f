










      Subroutine Pccd_doci_run(Wr,Lwr,Ec,Eref,Enuc)
C
C With pccd_doci.in present: DOCI at the orbitals of the first pair-form
C pCCD point of the run (Pccd_pair_drv, after Pccd_pair_slv), once per
C run (JOBARC PCCDDOCI). For DOCI at converged OO-pCCD orbitals, restart
C from them (oorot.in LOAD + c0.dat). Scratch: the Davidson area of
C /Pccdrsp/ (Wr(Ldv), Ldvl words), free at this point. Prints and
C stores DOCIENER (total, with Enuc); Ec, Eref: pCCD correlation and
C electronic reference energy at the same orbitals.
C
      Implicit Double Precision(A-H,O-Z)
      Logical Exist
      Dimension Wr(Lwr)
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



      Inquire(File="pccd_doci.in",Exist=Exist)
      If (.Not. Exist) Return
      Call Getrec(0,"JOBARC","PCCDDOCI",Length,Junk)
      If (Length .Gt. 0) Return
      Call Putrec(20,"JOBARC","PCCDDOCI",1,1)
      Call Pccd_doci(Wr(Le),Wr(Lw),Wr(Lk),Nb,No,Wr(Ldv),Ldvl,Edoci,
     +               Ndim,Nit,Rn,Iconv)
      If (Iconv .Lt. 0) Then
         Write(6,"(2a)") " @-DOCI not done: determinant space too",
     +                   " large for the memory (or Nb/No limits)."
         Return
      Endif
      Etot  = Edoci + Enuc
      Epccd = Eref + Ec + Enuc
      Write(6,"(a,i12,a,i4,a,1p,E10.3)") " @-DOCI dimension",Ndim,
     +     ",",Nit," Davidson iterations, |r| =",Rn
      If (Iconv .Ne. 1) Write(6,"(a)")
     +     " @-DOCI WARNING: Davidson not converged."
      Write(6,"(a,F20.12)") " @-DOCI energy        ",Etot
      Write(6,"(a,F20.12)") " @-DOCI pCCD energy   ",Epccd
      Write(6,"(a,F16.9,a)") " @-DOCI pCCD - DOCI   ",
     +     (Epccd-Etot)*1.0D3," mEh"
      Call Putrec(20,"JOBARC","DOCIENER",Iintfp,Etot)
      Return
      End
