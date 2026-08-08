






































































































































































































      Subroutine Fix_qrhf_evecs(EvecA,EvecB,Dcore,Maxdcor,Nbas,Iuhf)

      Implicit Double Precision(A-H, O-Z)

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end


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



c symm2.com : begin

c This is initialized in vscf/symsiz.

      integer nirrep,      nbfirr(8),   irpsz1(36),  irpsz2(28),
     &        irpds1(36),  irpds2(56),  irpoff(9),   ireps(9),
     &        dirprd(8,8), iwoff1(37),  iwoff2(29),
     &        inewvc(maxbasfn),         idxvec(maxbasfn),
     &        itriln(9),   itriof(8),   isqrln(9),   isqrof(8),
     &        mxirr2
      common /SYMM2/ nirrep, nbfirr, irpsz1, irpsz2, irpds1, irpds2,
     &               irpoff, ireps,  dirprd, iwoff1, iwoff2, inewvc,
     &               idxvec, itriln, itriof, isqrln, isqrof, mxirr2
c symm2.com : end
C
      Dimension EvecA(Nbas*Nbas),EvecB(Nbas*Nbas),Dcore(Maxdcor)
      Dimension Idummy(Maxbasfn),Idummy2(Maxbasfn)
      Dimension Idummy3(Maxbasfn),Ilocate(Maxbasfn)
      Dimension Irem(8,2),Iadd(8,2), Nocc_save(16)
C
      Common /Popul/ Nocc(8,2)
C    
      Ione = 1
      Call Izero(Irem,16)
      Call Icopy(16, Nocc, 1, Nocc_save, 1)
    
      Do I=1, Nbas
         Ilocate(I) = I
      Enddo

      CALL GETREC(-1,'JOBARC','QRHFTOT ',Ione,NMODIFY)
      CALL GETREC(-1,'JOBARC','QRHFIRR ',NMODIFY,IDUMMY)
      CALL GETREC(-1,'JOBARC','QRHFLOC ',NMODIFY,IDUMMY2)
      CALL GETREC(-1,'JOBARC','QRHFSPN ',NMODIFY,IDUMMY3)
      CALL GETREC(-1,'JOBARC','OSCALC  ',Ione,IOS)
      
      Call Getrec(20, "JOBARC", "SCFEVCA0", Nbas*Nbas*Iintfp,
     &            EvecA)
      If (Iuhf .Ne. 0) Call Getrec(20, "JOBARC", "SCFEVCB0", 
     &                      Nbas*Nbas*Iintfp, EvecB)

      Do I = 1, Nmodify

         Irrp = Idummy(I)
         ILoc = Idummy2(I)
         Ispn = Idummy3(I)
         
         If (Irrp .lt. 0) Then
             If (Ispn .ne. 1) Ispn = 2
                Irrp = -Irrp
                Iposorg = Irpoff(Irrp) + Nocc(Irrp,Ispn) + 1 -
     &                    Max(Iloc,1)  + Irem(Irrp,Ispn) - 
     &                    Iadd(Irrp,Ispn)
                Iposabs = Ilocate(Iposorg)
                Irem(Irrp,Ispn) = Irem(Irrp,Ispn) + 1
                Ichk = 1 + Nocc(Irrp,Ispn) - Inewvc(Iposabs) 
                Ibgn = Nocc(Irrp,Ispn)
                Iend = Nocc(Irrp,Ispn) + 1  - Ichk
C                If (Ichk .gt. 1) Then
                    If (Ios .eq. 0) then
                        Ioff_bgn  = Nbas*(Ibgn-1) + 1
                        Ioff_end  = Nbas*(Iend-1) + 1
                        Call Dswap(Nbas, EvecA(Ioff_end), 1, 
     &                             EvecA(Ioff_bgn), 1)
                        Call Dswap(Nbas, EvecB(Ioff_end), 1, 
     &                             EvecB(Ioff_bgn), 1)

                        Itmp = Ilocate(Ibgn)
                        Ilocate(Ibgn) = Ilocate(Iend) 
                        Ilocate(Iend) = Itmp
                    Else
                        if (ispn .eq. 1) then
                           Ioff_bgn  = Nbas*(Ibgn-1) + 1
                           Ioff_end  = Nbas*(Iend-1) + 1
                           Call Dswap(Nbas, EvecA(Ioff_end), 1,
     &                                EvecA(Ioff_bgn), 1)

                           Itmp = Ilocate(Ibgn)
                           Ilocate(Ibgn) = Ilocate(Iend)
                           Ilocate(Iend) = Itmp
                        elseif (ispn .eq. 2) then
                           Ioff_bgn  = Nbas*(Ibgn-1) + 1
                           Ioff_end  = Nbas*(Iend-1) + 1
                           Call Dswap(Nbas, EvecA(Ioff_end), 1,
     &                                EvecB(Ioff_bgn), 1)

                           Itmp = Ilocate(Ibgn)
                           Ilocate(Ibgn) = Ilocate(Iend)
                           Ilocate(Iend) = Itmp
                        endif 
                    Endif 
C                Endif 

                 Nocc(Irrp, Ispn) = Nocc(Irrp,Ispn) - 1
         Else
                if (Ispn .ne. 2) Ispn = 1
                   Iposorg = Irpoff(Irrp) + Nocc(Irrp,Ispn) +
     &                       Max(Iloc,1)  + Irem(Irrp,Ispn) -
     &                       Iadd(Irrp,Ispn)
                Iposabs = Ilocate(Iposorg)
                Iadd(Irrp,Ispn) = Iadd(Irrp,Ispn) + 1
                Ichk = Inewvc(Iposabs) - Nocc(Irrp,Ispn)
                Ibgn = Nocc(Irrp,Ispn) + 1
                Iend = Nocc(Irrp,Ispn) + ichk

                If (Ichk .gt. 1) then
                    If (Ios .eq. 0) then
                        Ioff_bgn  = Nbas*(Ibgn-1) + 1
                        Ioff_end  = Nbas*(Iend-1) + 1
                        Call Dswap(Nbas, EvecA(Ioff_end), 1, 
     &                             EvecA(Ioff_bgn), 1)
                        Call Dswap(Nbas, EvecB(Ioff_end), 1, 
     &                             EvecB(Ioff_bgn), 1)
                   
                       Itmp = Ilocate(Ibgn)
                       Ilocate(Ibgn) = Ilocate(Iend) 
                       Ilocate(Iend) = Itmp
                     Else
                       If (ispn .eq. 1) then
                          Ioff_bgn  = Nbas*(Ibgn-1) + 1
                          Ioff_end  = Nbas*(Iend-1) + 1
                          Call Dswap(Nbas, EvecA(Ioff_end), 1, 
     &                               EvecA(Ioff_bgn), 1)
                   
                          Itmp = Ilocate(Ibgn)
                          Ilocate(Ibgn) = Ilocate(Iend) 
                          Ilocate(Iend) = Itmp
                        Elseif (ispn .eq. 2) then 
                          Ioff_bgn  = Nbas*(Ibgn-1) + 1
                          Ioff_end  = Nbas*(Iend-1) + 1
                          Call Dswap(Nbas, EvecA(Ioff_end), 1,
     &                               EvecB(Ioff_bgn), 1)
                          Itmp = Ilocate(Ibgn)
                          Ilocate(Ibgn) = Ilocate(Iend)
                          Ilocate(Iend) = Itmp
                        Endif 

                     Endif 

                 Endif 
                 Nocc(Irrp, Ispn) = Nocc(Irrp,Ispn) + 1
         Endif 
           
      Enddo 

      Call Putrec(20, "JOBARC", "SCFEVCA0", Nbas*Nbas*Iintfp,
     &            EvecA)
      If (Iuhf .Ne. 0) Call Putrec(20, "JOBARC", "SCFEVCB0", 
     &                      Nbas*Nbas*Iintfp, EvecB)
      Call Icopy(16, Nocc_save, 1, Nocc, 1)
 
      Return
      End

