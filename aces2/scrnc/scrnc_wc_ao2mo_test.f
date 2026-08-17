










      Subroutine Scrnc_wc_ao2mo_test(Work,Maxcor,Iuhf,Jj_lrespn,
     +                               Nao_pairs)
C
C Standing diagnostic (guarded by 1, permanent, not
C exploratory cruft -- same convention as Scrnc_bcktran_test.F),
C called from Scrnc_driver.F right after the real Jj_lrespn (W(omega)
C in the AO-pair basis) is available. Transforms it into the MO
C particle-hole basis via Scrnc_wc_ao2mo (new, this session -- see
C SCRNC_DEBUG_LOG.md 2026-08-17 "self-energy architecture" entry) using
C the REAL alpha-spin SCF eigenvectors for whatever test case is
C running, then independently re-derives every Wc_mo element via a
C brute-force quadruple-nested-loop sum sharing nothing with
C Scrnc_wc_ao2mo's block-Xgemm code path except the raw AO-pair index
C convention -- same methodology as Scrnc_bcktran_left.F's own
C brute-force cross-check.
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Maxcor)
      Double Precision Jj_lrespn(Nao_pairs,Nao_pairs)


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



c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
C
      Character*8 Labelc
      Data Labelc /"SCFEVECA"/
C
      Nao = Nint(Sqrt(Dble(Nao_pairs)))
      Nmo = Nocco(1) + Nvrto(1)
      Nmo_pairs = Nmo * Nmo
C
      I_c      = 1
      I_wcmo   = I_c    + Nao*Nmo
      I_scr    = I_wcmo + Nmo_pairs*Nmo_pairs
      Iend     = I_scr  + Nmo_pairs*Nao_pairs + Nao*Nao + Nao*Nmo +
     +                     Nmo*Nmo
C
      If (Iend .GT. Maxcor) Call Insmem("Scrnc_wc_ao2mo_test",
     +                                  Iend,Maxcor)
C
      Call Getrec(20,"JOBARC",Labelc,Iintfp*Nao*Nmo,Work(I_c))
C
      Call Scrnc_wc_ao2mo(Work(I_scr),Maxcor-I_scr+1,
     +                    Jj_lrespn,Nao,Nao_pairs,
     +                    Work(I_c),Nmo,Work(I_wcmo),Nmo_pairs)
C
      Call Scrnc_wc_ao2mo_bfcheck(Jj_lrespn,Nao,Nao_pairs,
     +                            Work(I_c),Nmo,Work(I_wcmo),
     +                            Nmo_pairs)
C
      Return
      End
C
C-----------------------------------------------------------------
C
      Subroutine Scrnc_wc_ao2mo_bfcheck(Jj_lrespn,Nao,Nao_pairs,
     +                                  C,Nmo,Wc_mo,Nmo_pairs)
C
C Brute-force reference: Wc_mo_bf(Ip,Iq,Ir,Is) = Sum(Mu,Nu,Lam,Nsig)
C   C(Mu,Ip)*C(Nu,Iq)*Jj_lrespn(Iao_pair,Jao_pair)*C(Lam,Ir)*C(Nsig,Is),
C direct quadruple loop, no blocking/Xgemm -- independent of
C Scrnc_wc_ao2mo_do's algorithm. All loop variables deliberately kept
C in the I-N implicit-integer range -- 2026-08-17: using P,Q,R,S,Sig
C (O-Z range, implicit DOUBLE PRECISION here) plus an explicit
C "Integer" override was tried and broke under this module's "-i8"
C build flag (bare explicit Integer declarations don't reliably
C inherit the 8-byte default the way implicit I-N typing does); see
C SCRNC_DEBUG_LOG.md and Scrnc_wc_ao2mo.F's own header comment.
C Jj_lrespn is explicitly Double Precision below -- its name starts
C with 'J', which IS in the implicit I-N range this file's own
C Implicit Double Precision (A-H,O-Z) statement does NOT cover, so a
C bare Dimension here would silently make it Integer*8 (this was the
C actual root cause of the whole corruption bug, see SCRNC_DEBUG_LOG.md).
C
      Implicit Double Precision (A-H,O-Z)
      Double Precision Jj_lrespn(Nao_pairs,Nao_pairs)
      Dimension C(Nao,Nmo)
      Dimension Wc_mo(Nmo_pairs,Nmo_pairs)
C
      Xmaxdiff = 0.0D0
      Ip_worst=0
      Iq_worst=0
      Ir_worst=0
      Is_worst=0
C
      Do Ip = 1, Nmo
       Do Iq = 1, Nmo
        Ipq = (Ip-1)*Nmo + Iq
        Do Ir = 1, Nmo
         Do Is = 1, Nmo
          Irs = (Ir-1)*Nmo + Is
C
          Sum = 0.0D0
          Do Mu = 1, Nao
           Do Nu = 1, Nao
            Iao_pair = (Mu-1)*Nao + Nu
            Cmp = C(Mu,Ip)*C(Nu,Iq)
            If (Cmp .NE. 0.0D0) Then
             Do Lam = 1, Nao
              Do Nsig = 1, Nao
               Jao_pair = (Lam-1)*Nao + Nsig
               Term = Cmp*Jj_lrespn(Iao_pair,Jao_pair)*
     +                     C(Lam,Ir)*C(Nsig,Is)
               Sum = Sum + Term
              Enddo
             Enddo
            Endif
           Enddo
          Enddo
C
          Diff = Abs(Sum - Wc_mo(Ipq,Irs))
          If (Diff .GT. Xmaxdiff) Then
             Xmaxdiff = Diff
             Ip_worst = Ip
             Iq_worst = Iq
             Ir_worst = Ir
             Is_worst = Is
             Sum_worst = Sum
             Fast_worst = Wc_mo(Ipq,Irs)
          Endif
C
         Enddo
        Enddo
       Enddo
      Enddo
C
      Write(6,"(a)") " @-Scrnc_wc_ao2mo_bfcheck"
      Write(6,"(a,1PD15.6)")
     +   "   max|brute force - Scrnc_wc_ao2mo| = ", Xmaxdiff
      Write(6,"(a,4(1x,I3))")
     +   "   worst (p,q,r,s) = ", Ip_worst,Iq_worst,Ir_worst,Is_worst
      Write(6,"(a,1PD15.6)") "   brute value  = ", Sum_worst
      Write(6,"(a,1PD15.6)") "   fast  value  = ", Fast_worst
C
      Return
      End
