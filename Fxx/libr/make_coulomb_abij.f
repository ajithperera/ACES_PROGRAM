      Subroutine Make_coulomb_abij(W,Maxcor,ListC,ISPIN,Type)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION W(Maxcor)
      CHARACTER*4 TYPE



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
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 


      LEN_C = IDSYMSZ(1,ISYTYP(1,LISTC),ISYTYP(2,LISTC))
C
C Form <AB|IJ> from <AI|BJ> or retrive <AB|IJ> from the disk.
C
      If (Type .EQ. "AIBJ") Then

         NSCRSZ=NVRTO(ISPIN)*NVRTO(ISPIN)+NOCCO(ISPIN)*NOCCO(ISPIN)+
     &          NVRTO(ISPIN)*NOCCO(ISPIN)

         I000 = 1
         I010 = I000 + Len_c
         I020 = I010 + Len_c
         I030 = I020 + Nscrsz

         CALL GETALL(W(I000), LEN_C, 1, LISTC) 

         CALL SSTGEN(W(I000),W(I010),LEN_C,VRT(1,ISPIN),
     &               POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     &               W(I020),1,"1324") 
         CALL DCOPY(Len_c, W(I010),W(I000))
      Else

         I000 = 1
         I010 = I000 + Len_c
         CALL GETALL(W(I000), LEN_C, 1, LISTC) 

      Endif 

      call checksum("W<AB|IJ>",W(I000),LEN_C)

      I020 = I010 + Len_C

C Take the transpose; <AB|IJ>-><IJ|AB>; WIN(LEN1,LEN2)->WOUT(LEN2,LEN1)

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)
         LEN1   = IRPDPD(IRREPR,20+ISPIN)
         LEN2   = IRPDPD(IRREPL,18+ISPIN)
         CALL TRANSP(W(I000+Incrm1),W(I010+Incrm),LEN2,LEN1)
         INCRM = IINCRM + LEN1 * LEN2
      ENDDO 

      call checksum("W<IJ|AB>",W(I020),LEN_C)

      I030 = I020 + Len_C

C Change <IJ|AB> to <I<J|AB>; W(NDSSIZ,NUMDISC) -> W(NDSSIZW,NUMDISC) 

      IOFF2 = I020 
      IOFF1 = I010

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)
         NUMDISC=IRPDPD(IRREPR,18+ISPIN)

         NDSSIZC=IRPDPD(IRREPL,20+ISPIN)
         NDSSIZW=IRPDPD(IRREPL,2+ISPIN)

         CALL SQSYM(IRREPL,POP(1,ISPIN),NDSSIZW,NDSSIZC,NUMDISC,
     &              W(IOFF2),W(IOFF1))

         IOFF1 = IOFF1 + NDSSIZC * NUMDISC
         IOFF2 = IOFF2 + NDSSIZW * NUMDISC 
      ENDDO 

      length = idsymsz(1,2+ISPIN,18+ISPIN)
      call checksum("W<I<J|AB>",W(I020),length)

C Take the transpose; <I<J|AB> -> <AB|I<J>; WIN(LEN1,LEN2)->WOUT(LEN2,LEN1)

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)
         LEN1   = IRPDPD(IRREPR,12+ISPIN)
         LEN2   = IRPDPD(IRREPL,18+ISPIN)
         CALL TRANSP(W(I020+Incrm1),W(I010+Incrm),LEN2,LEN1)
         INCRM = IINCRM + LEN1 * LEN2
      ENDDO 

      call checksum("W<AB|I<J>",W(I010),length)

C Change <IJ|AB> to <I<J|AB>; W(NDSSIZ,NUMDISC) -> W(NDSSIZW,NUMDISC) 

      IOFF1 = I010 
      IOFF2 = I020

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)

         NUMDISC=IRPDPD(IRREPR,2+ISPIN)
         NDSSIZC=IRPDPD(IRREPL,18+ISPIN)
         NDSSIZW=IRPDPD(IRREPL,0+ISPIN)

         CALL SQSYM(IRREPL,VRT(1,ISPIN),NDSSIZW,NDSSIZC,NUMDISC,
     &             W(IOFF2),W(IOFF1))

         IOFF1 = IOFF1 + NDSSIZC * NUMDISC
         IOFF2 = IOFF2 + NDSSIZW * NUMDISC 

      ENDDO 

      length = idsymsz(1,ISPIN,2+ISPIN)
      call checksum("W<A<B|I<J>",W(I020),length)

      RETURN
      END
