










      Subroutine  a3_symadapt_scfvecs(Scfvec_a, Scfveqc_b, Scfevl_a,
     &                                Scfevl_b, Tmp1, Tmp2,
     &                                Oed2AScale, Ioed2Aorder,
     &                                Tmp2_a, Tmp2_b,
     &                                Nbfns,
     &                                Naobfns, Nbfirr, Nirrep, Iuhf,
     &                                Spherical, Work, Imemleft)

      Implicit Double Precision (A-H, O-Z)
      Logical Spherical
      Parameter (Tol = 1.0D-09)



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



C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)


c ***NOTE*** This is a genuine (though not serious) limit on what Aces3 can do.
c     12 => s,p,d,f,g,h,i,j,k,l,m,n
      integer maxangshell
      parameter (maxangshell=12)



      Dimension Scfvec_a(Naobfns*Naobfns), Scfvec_b(Naobfns*Naobfns), 
     &          Scfevl_a(Nbfns), Scfevl_b(Nbfns), 
     &          Tmp1(Naobfns*Naobfns), Nbfirr(8),
     &          Tmp2(Naobfns*Naobfns),
     &          Tmp2_a(Naobfns*Naobfns), Tmp2_b(Naobfns*Naobfns), 
     &          Oed2AScale(Naobfns), Ioed2Aorder(Naobfns)

      Dimension Work(Imemleft)
C
      Dimension Nocc(8,2)
      Dimension Nprim_shell(Maxangshell*Mxatms)
      Dimension Orig_nprim_shell(Maxangshell*Mxatms)
      Integer   Reorder_Shell(Maxangshell*Mxatms)
C
      Call Getrec(20, "JOBARC", "SCFEVCA0", Nbfns*Nbfns*Iintfp,
     &            Scfvec_a)
      Call Filter(Scfvec_a, Nbfns*Nbfns, Tol)
      Call Getrec(20, "JOBARC", "SCFEVLA0", Nbfns*Iintfp, Scfevl_a)

      If (Iuhf .EQ. 1) then
         Call Getrec(20, "JOBARC", "SCFEVCB0", Nbfns*Nbfns*Iintfp,
     &               Scfvec_b)
         Call Filter(Scfvec_b, Nbfns*Nbfns, Tol)
         Call Getrec(20, "JOBARC", "SCFEVLB0", Nbfns*Iintfp, Scfevl_b)

      Endif 

C
C This skipped block seems to be unnessary. The SCF vectors comming
C out of ACESIII is in the ACES II order (s,p1x,p1y,p1z,p2x,p2y,p2z;..
C and so on. Neverthelles I am not deleting since things can change. 
C
C
C 
C Write the the scf vectors in ACES II order so molden
C interface that handle excited state densities will function
C correctly. Note that I did not overwrite the records SCFEVCA0
C and SCFEVCB0 (which I could have done). That would have prevented
C me from runing xa2proc over and over without getting wrong 
c results!

      Call Putrec(20, "JOBARC", "SCFVECA0", Nbfns*Nbfns*Iintfp,
     &            Scfvec_a)
      If (Iuhf .EQ. 1) Call Putrec(20, "JOBARC", "SCFVECB0",
     &                             Nbfns*Nbfns*Iintfp, Scfvec_b)
C
C Generate the occupation numbers for each irrep based on eigen
C values and the number of basis functions per irrep.

      Call Occupy(Nirrep, Nbfirr, Nbfns, Scfevl_a, Work, Nocc(1,1),
     &            1)
      If (Iuhf .EQ. 1) Then
         Call Occupy(Nirrep, Nbfirr, Nbfns, Scfevl_b, Work, Nocc(1,1),
     &               2)
      Else
         Call Icopy(8, Nocc(1, 1), 1, Nocc(1,2), 1)
      Endif
C
C First convert from Spherical to Cartesian (if the calculation is
C in Cartesian this should do nothing).

      Call Getrec(20, "JOBARC", "CMP2CART", Nbfns*Naobfns*Iintfp,
     &            Tmp1)

      Call Xgemm("N", "N", Naobfns, Nbfns, Nbfns, 1.0D0, Tmp1,
     &            Naobfns, Scfvec_a, Nbfns, 0.0D0, Tmp2, Naobfns)

      Call Dcopy(Naobfns*Naobfns, Tmp2, 1, Scfvec_a, 1)
      
      Call Putrec(20, "JOBARC", "SCFVECA3", Naobfns*Naobfns*Iintfp,
     &            Scfvec_a)

      If (Iuhf .EQ. 1) Then
         Call Xgemm("N", "N", Naobfns, Nbfns, Nbfns, 1.0D0, Tmp1, 
     &               Naobfns, Scfvec_b, Nbfns, 0.0D0, Tmp2, Naobfns)

         Call Dcopy(Naobfns*Naobfns, Tmp2, 1, Scfvec_b, 1)

         Call Putrec(20, "JOBARC", "SCFVECB3", Naobfns*Naobfns*Iintfp,
     &               Scfvec_b)
      Endif
C Note that seventh argument to both calls is the same. This is
C because alpha and beta eigenvectors are kept in two different
C arrays instead of a one in which beta vectors comes latter.
C
      Call Get_irreps_gs(Scfvec_a, Scfevl_a, Work, Imemleft*Iintfp, 
     &                Nbfns, Naobfns, 1, Nocc, Iuhf, 1, "GROUND ")
      If (Iuhf .EQ. 1) Call Get_irreps_gs(Scfvec_b, Scfevl_b, Work, 
     &                                    Imemleft*Iintfp, Nbfns, 
     &                                    Naobfns, 1, Nocc, Iuhf,
     &                                    2, "GROUND ")
C
      Return
      End

