










      Subroutine Scale_ints_driver(Work,Maxcor,IUhf)
   
      Implicit Integer (A-Z)

      Double Precision Work(Maxcor), Delta, Scale_factor,Scale 
      Double Precision Buf(600)
      Dimension Ibuf(600)
      Logical UHF
      Character*80 FNAME

c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end


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



c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
C
C This is a routine written for debugging purposes. What this 
C is to add delat to type of MO two electon integrals.

CSSS      Return
C
      Delta        = 0.001D0
      Scale        = 0.0D0
      Scale_factor = 1.0D0 + Delta*Scale 
      Zlich        = 0.0D0

      UHF    = .False.
      UHF    = (Iuhf .EQ. 1)
      Irrepx = 1 
      Imode  = 0 
      Write(6,*) 

C Lets scale all the <pq||rs> integrals 
C
C Remove the dependency of t1int2, dwmbej, and t1int1 on integral
C list 17 and 18. These terms also dependent on 23,24,25 and 26
C but those do not interfer with debugging G(IJ,AB) in the denisty
C code. Dependecy on 17 and 18 cause G(IA,JB) to intefere with
C G(IJ,AB). Make copies of 17 and 18 and store them in 121 and
C 122. Also, make duplicates copies of 14,15 and 16 in 147,148
C and 149
     
      Write(6,*) 
      Write(6,"(a)") "-----Warning! Only do this for debugging-----" 
      Write(6,"(a,a)") " Duplicate copies of W list 17 and 18 are",
     +               " formed on 121 and 122"
      Write(6,"(a,a)") " Duplicate copies of W list 14-16 are",
     +               " formed on 147-149"
      Write(6,*)

      If (Iuhf .Ne. 0) Then
         Call Inipck(1,10,9,121,Imode,0,1)
         Length_121 = IDSYMSZ(IRREPX,ISYTYP(1,17),ISYTYP(2,17))
         Call Getall(Work, Length_121, Irrepx, 17)
         Call Putall(Work, Length_121, Irrepx, 121)

C         Call Inipck(1,2,4, 148,Imode,0,1)
C         Call Inipck(1,13,14,149,Imode,0,1)
C         Length_147 = IDSYMSZ(IRREPX,ISYTYP(1,147),ISYTYP(2,147))
C         Length_148 = IDSYMSZ(IRREPX,ISYTYP(1,148),ISYTYP(2,148))
C         Length_149 = IDSYMSZ(IRREPX,ISYTYP(1,149),ISYTYP(2,149))
C         Call Getall(Work, Length_147, Irrepx, 14)
C         Call Putall(Work, Length_147, Irrepx, 147)
C         Call Getall(Work, Length_148, Irrepx, 15)
C         Call Putall(Work, Length_148, Irrepx, 148)
C         Call Getall(Work, Length_149, Irrepx, 16)
C         Call Putall(Work, Length_149, Irrepx, 149)
      Endif 

      Call Inipck(1,10,10,147,Imode,0,1)
      Length_147 = IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
      Call Getall(Work, Length_147, Irrepx, 25)
      Call Putall(Work, Length_147, Irrepx, 147)
      Call Inipck(1,9,10,122,Imode,0,1)
      Length_122 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
      Call Getall(Work, Length_122, Irrepx, 18)
      Call Putall(Work, Length_122, Irrepx, 122)
      Call Inipck(1,9,10,148,Imode,0,1)

C      Call Inipck(1,13,14,149,Imode,0,1)
C      Length_147 = IDSYMSZ(IRREPX,ISYTYP(1,147),ISYTYP(2,147))
C      Length_149 = IDSYMSZ(IRREPX,ISYTYP(1,149),ISYTYP(2,149))
C      Call Getall(Work, Length_147, Irrepx, 14)
C      Call Putall(Work, Length_147, Irrepx, 147)
C      Call Getall(Work, Length_149, Irrepx, 16)
C      Call Putall(Work, Length_149, Irrepx, 149)


C#ifdef _SCALE_IBJA

      Write(6,"(a,F7.5)") "<ph||ph> integrals are incremented by ", 
     &                    delta*scale 
      If (Uhf) Then
         Length_23  = IDSYMSZ(IRREPX,ISYTYP(1,23),ISYTYP(2,23))
         Length_24  = IDSYMSZ(IRREPX,ISYTYP(1,24),ISYTYP(2,24))
         Length_25  = IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
         Length_26  = IDSYMSZ(IRREPX,ISYTYP(1,26),ISYTYP(2,26))
         Length_121 = IDSYMSZ(IRREPX,ISYTYP(1,17),ISYTYP(2,17))
         Length_122 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))

C         Call Getall(Work, Length_23, Irrepx, 23)
C         Call Dscal(Length_23, Scale_factor, Work, 1)
C         Call Putall(Work, Length_23, Irrepx, 23)
C
C         Call Getall(Work, Length_24, Irrepx, 24)
C         Call Dscal(Length_24, Scale_factor, Work, 1)
C         Call Putall(Work, Length_24, Irrepx, 24)
C          
C         Call Getall(Work, Length_25, Irrepx, 25)
C         Call Dscal(Length_25, Scale_factor, Work, 1)
C         Call Putall(Work, Length_25, Irrepx, 25)
C
C         Call Getall(Work, Length_26, Irrepx, 26)
C         Call Dscal(Length_26, Scale_factor, Work, 1)
C         Call Putall(Work, Length_26, Irrepx, 26)

         Call Getall(Work, Length_121, Irrepx, 121)
         Call Dscal(Length_121, Scale_factor, Work, 1)
         Call Putall(Work, Length_121, Irrepx, 121)

         Call Getall(Work, Length_122, Irrepx, 122)
         Call Dscal(Length_122, Scale_factor, Work, 1)
         Call Putall(Work, Length_122, Irrepx, 122)

      Else

C RHF needs all list 122,147 (instead of 25 (used in G(ij,ab)
C 23 and 25 scalled. This yield the sum of three G(ia,jb) 
C contributions.

        Length_122 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
        Length_147 = IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
        Length_23  = IDSYMSZ(IRREPX,ISYTYP(1,23),ISYTYP(2,23))
        Length_25  = IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))

         Call Getall(Work, Length_122, Irrepx, 122)
         Call Dscal(Length_122, Scale_factor, Work, 1)
         Call Putall(Work, Length_122, Irrepx, 122)

         Call Getall(Work, Length_147, Irrepx, 147)
         Call Dscal(Length_147, Scale_factor, Work, 1)
         Call Putall(Work, Length_147, Irrepx, 147)

         Call Getall(Work, Length_23, Irrepx, 23)
         Call Dscal(Length_23, Scale_factor, Work, 1)
         Call Putall(Work, Length_23, Irrepx, 23)

         Call Getall(Work, Length_25, Irrepx, 25)
         Call Dscal(Length_25, Scale_factor, Work, 1)
         Call Putall(Work, Length_25, Irrepx, 25)

      Endif 
           
C#endif 


      Return 
C
C One electron terms 
C
      Call Getrec(20,"JOBARC","NBASTOT ", 1, Norbs)
      LUNIT   = 10
      ILNBUF  = 600
      LDIM = NORBS*(NORBS+1)/2
      I000 = 1
      I010 = I000 + Ldim
      I020 = I010 + Norbs*Norbs

      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(LUNIT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &           ACCESS='SEQUENTIAL')
      CALL LOCATE(LUNIT,'ONEHAMIL')
      CALL ZERO(Work(I010),LDIM)
      NUT = ILNBUF
      DO WHILE (NUT.EQ.ILNBUF)
         READ(LUNIT) BUF, IBUF, NUT
         DO INT = 1, NUT
            Work((I000-1)+(IBUF(INT))) = BUF(INT)
         END DO
      END DO
      CALL EXPND2(Work(I000),Work(I010),Norbs)

CSSS      Write(6,*) "The one electron ints in AO basis"
CSSS      call output(Work(I010),1,norbs,1,norbs,norbs,norbs,1)
    
      I030 = I020 +  Norbs*Norbs
      I040 = I030 +  Norbs*Norbs
      I050 = I040 +  2*Norbs*Norbs
      If(I050 .GT. Maxcor) Call Insmem("@-Scale_ints",I050,Maxcor)

      Call Ao2mo2(Work(I010),Work(I020),Work(I030),Work(I040),
     &                Norbs,Norbs,1)
      If (UHF) Call Ao2mo2(Work(I010),Work(I020),Work(I030),
     &                     Work(I040),Norbs,Norbs,2)

CSSS      Write(6,*) "The one electron ints in MO basis"
CSSS      call output(Work(I020),1,norbs,1,norbs,norbs,norbs,1)



      Return
      End
     
      Subroutine OO_Block(Hfull,Hblock,Norbs,Pop,Vrt,Nocco,Nvrto,
     &                    Nirrep, Ispin)

      Implicit Integer (A-Z)

      Double Precision Hfull(Norbs,Norbs), Hblock(Norbs,Norbs)
      Integer Pop(8,2),Vrt(8,2),Nocco(2),Nvrto(2)

      Ioff  = 0
      Call Dzero(Hblock,Norbs*Norbs)
      Do Irrep = 1, Nirrep
         Do J = 1, Pop(Irrep,Ispin)
            Do I = 1, Pop(Irrep,Ispin)
               Hblock(I+Ioff,J+Joff) = Hfull(I+Ioff,J+Ioff) 
            Enddo 
         Enddo 
         Ioff = Ioff + Pop(Irrep,Ispin)
      Enddo 
      Return
      End

      Subroutine VV_Block(Hfull,Hblock,Norbs,Pop,Vrt,Nocco,Nvrto,
     &                    Nirrep,Ispin)

      Implicit Integer (A-Z)

      Double Precision Hfull(Norbs,Norbs), Hblock(Norbs,Norbs)
      Integer Pop(8,2),Vrt(8,2),Nocco(2),Nvrto(2)

      Ioff   = Nocco(Ispin)
      Call Dzero(Hblock,Norbs*Norbs)
      Do Irrep = 1, Nirrep
         Do J = 1, Vrt(Irrep,Ispin)
            Do I = 1, Vrt(Irrep,Ispin)
               Hblock(I+Ioff,J+Ioff) = Hfull(I+Ioff,J+Ioff)
            Enddo
         Enddo
         Ioff = Ioff + Vrt(Irrep,Ispin)
      Enddo

      Return
      End

      Subroutine OV_Block(Hfull,Hblock,Norbs,Pop,Vrt,Nocco,Nvrto,
     &                    Nirrep,Ispin)

      Implicit Integer (A-Z)

      Double Precision Hfull(Norbs,Norbs), Hblock(Norbs,Norbs)
      Integer Pop(8,2),Vrt(8,2),Nocco(2),Nvrto(2)

      Ioff  = 0
      Aoff  = Nocco(Ispin)
      Call Dzero(Hblock,Norbs*Norbs)
      Do Irrep = 1, Nirrep
         Do I = 1, Pop(Irrep,Ispin)
            Do A = 1, Vrt(Irrep,Ispin)
               Hblock(A+Aoff,I+Ioff) = Hfull(A+Aoff,I+Ioff)
               Hblock(I+Ioff,A+Aoff) = Hfull(A+Aoff,I+Ioff)
            Enddo
         Enddo
         Ioff = Ioff + Pop(Irrep,Ispin)
         Aoff = Ioff + Vrt(Irrep,Ispin)
      Enddo

      Return
      End

