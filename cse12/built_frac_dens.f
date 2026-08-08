










      Subroutine Built_frac_dens(Dens_aa,Dens_bb,Dens_tt,Evec,
     +                          Docc,Scr,Nbasis)

      Implicit Double Precision(A-H,O-Z)

      Dimension Dens_aa(Nbasis,Nbasis) 
      Dimension Dens_bb(Nbasis,Nbasis) 
      Dimension Dens_tt(Nbasis,Nbasis)

      Dimension Docc(Nbasis,Nbasis) 
      Dimension Evec(Nbasis*Nbasis)
      Dimension Scr(Nbasis*Nbasis)

      Integer P,Q,R,S

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end

      Integer cc_maxcyc
      Integer Act_min_a,Act_min_b,Act_max_a,Act_max_b
      Logical Ring_cc,Brueck,Active_space
      Double Precision ocn_oa,Ocn_ob,Ocn_va,Ocn_vb
      Double Precision Denom_tol,Brueck_tol
      Dimension E_corr(0:500)

      Common /ccsdlight_vars/Ring_cc,Brueck,cc_conv,cc_maxcyc,
     +                       ocn_oa(Maxbasfn),ocn_ob(Maxbasfn),
     +                       ocn_va(Maxbasfn),ocn_vb(Maxbasfn),
     +                       E_corr,Denom_tol,Brueck_tol,
     +                       Act_min_a,Act_min_b,Act_max_a,
     +                       Act_max_b,Active_space





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




      Nbasis2 = Nbasis*Nbasis 
      Call Dzero(Docc,Nbasis2)

C built AA AO density matrix

      Do P = 1, Nbasis
         Docc(P,P) = Ocn_oa(P)
      Enddo

      Call Getrec(20,"JOBARC","SCFEVCA0",Nbasis*Nbasis*Iintfp,
     +            Evec)
      Call Form_dens(Dens_aa,Evec,Scr,Docc,Nbasis)

C built BB AO density matrix 

      Do P = 1, Nbasis 
         Docc(P,P) = Ocn_ob(P)
      Enddo

      Call Getrec(20,"JOBARC","SCFEVCB0",Nbasis*Nbasis*Iintfp,
     +               Evec)
      Call Form_dens(Dens_bb,Evec,Scr,Docc,Nbasis)
     
      Call Dcopy(Nbasis2,Dens_aa,1,Dens_tt,1)
      Call Daxpy(Nbasis2,1.0D0,Dens_bb,1,Dens_tt,1)
      Write(6,"(a)") "Alpha/beta and total density"
      call output(Dens_aa,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Dens_bb,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Dens_tt,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)

      Return
      End 


