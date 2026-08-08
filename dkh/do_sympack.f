










      Subroutine Do_sympack(Nbfns, Ovlp, Ao2so, Tmp1, Tmp2,
     +                      Itriang_length)

      Implicit None



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



      Integer Nbfns, Naobfns, Nbfirr(8), Irrep, Nirrep, Idex, I, J
      Integer Itriang_off(8)
      Integer Isquare_off(8)
      Integer Isquar2_off(8)
      Integer Ioff, Itriang_length

      Double Precision Ovlp(Nbfns,Nbfns), Tmp1(Nbfns,Nbfns),
     &                 Tmp2(Nbfns*Nbfns), Ao2So(Nbfns*Nbfns)

      Call Getrec(20, "JOBARC", "NIRREP  ", 1, Nirrep)
      Call Getrec(20, "JOBARC", "NUMBASIR", Nirrep, Nbfirr)
C
      Idex = 0
      Do Irrep = 1, Nirrep
         Do I = 1 + Idex, Nbfirr(Irrep) + idex
              Do J = 1+ Idex, Nbfirr(Irrep) + idex
              Ovlp(j, i) = Tmp1(j,i)
           Enddo
        Enddo
           Idex = Idex + Nbfirr(Irrep)
      Enddo

      Itriang_off(1) = 1
      Isquare_off(1) = 1
      Isquar2_off(1) = 1
      Do Irrep = 1, Nirrep-1
         Itriang_off(Irrep+1) = Itriang_off(Irrep) + (Nbfirr(Irrep)+1)*
     &                          Nbfirr(Irrep)/2
         Isquare_off(Irrep+1) = Isquare_off(Irrep) + Nbfirr(Irrep)
         Isquar2_off(Irrep+1) = Isquar2_off(Irrep) + Nbfirr(Irrep) *
     &                          Nbfirr(Irrep)
      Enddo

      Write(6,*)
      Write(6,"(a)")"Sym. Adapted PVP integrals sym. block basis"
      Do Irrep = 1, Nirrep
      If (Nbfirr(Irrep) .ne. 0) then
      Call output(Ovlp(Isquare_off(Irrep),Isquare_off(Irrep)),
     &            1, Nbfirr(Irrep), 1, Nbfirr(Irrep),
     &            Nbfns, Nbfns, 1)
      endif
      Enddo

      Ioff = 1
      Idex = 0
      Do Irrep = 1, Nirrep
         Ioff = Isquar2_off(Irrep)
         Do I =1 +Idex, Nbfirr(Irrep) + Idex
            Call Dcopy(Nbfirr(Irrep),Ovlp(1+idex,i),1,Ao2so(Ioff),1)
            Ioff = Ioff + Nbfirr(Irrep)
         Enddo
         Idex = Idex + Nbfirr(Irrep)
      Enddo

      Do Irrep = 1, Nirrep
         Call Squez2(AO2so(Isquar2_off(Irrep)),
     +               Tmp2(Itriang_off(Irrep)),Nbfirr(Irrep))
      Enddo

      Itriang_length = Itriang_off(Nirrep) + (Nbfirr(Nirrep)+1)*
     +                                        Nbfirr(Nirrep)/2
      Itriang_length = Itriang_length      -  1

      Return
      End

