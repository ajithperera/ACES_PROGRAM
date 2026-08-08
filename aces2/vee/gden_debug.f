










      Subroutine Gden_debug(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end


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



      COMMON /CALCINFO/ NROOT(8) 

      I000 = 1

      MAXLEN=0
      DO IRREPX = 1, NIRREP
         LEN = 0
         DO ISPIN = 1, IUHF+1
            LEN=LEN+IRPDPD(IRREPX,8+ISPIN)
         ENDDO
         LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         IF(IUHF.NE.0)THEN
           LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
           LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         ENDIF
         MAXLEN=MAX(MAXLEN,LEN)
      ENDDO

      Print*, ITOP,MAXLEN,MAXCOR

      IF (ITOP + MAXLEN .GT. MAXCOR) CALL
     +   INSMEM("@-gtden_debug",ITOP+MAXLEN,MAXCOR)

      NTOTAL = 0
      DO IRREP = 1, NIRREP
         NTOTAL = NTOTAL + NROOT(IRREP)
      ENDDO

      IDONE = 0
      DO IRREPX = 1, 1

         LENSZ = 0
         DO ISPIN = 1, IUHF+1
            LENSZ=LENSZ+IRPDPD(IRREPX,8+ISPIN)
         ENDDO
         LENSZ=LENSZ+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         IF(IUHF.NE.0)THEN
           LENSZ=LENSZ+IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
           LENSZ=LENSZ+IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         ENDIF

         Print*, IRREPX, NROOT(IRREPX),LENSZ
         DO IROOT = 1, NROOT(IRREPX)
            IDONE = IDONE + 1
            IF (IDONE .EQ. NTOTAL) THEN
               IGET = 2
               CALL GETLST(WORK(I000),IGET,1,1,1,472)
            ELSE
               IGET = 17 + (2*IDONE-1) + 1
               CALL GETLST(WOrk(I000),IGET,1,1,3,472)

               Write(6,*) "From Tdens_debug"
               write(6,"(2a,3(1x,i2))") " Reading the right vector",
     +                                  " of state: ",
     +                             Idone,iget,iroot
               Call checksum("list-472      :",work(i000),lensz,s)
            ENDIF

            JDONE = 0
            DO JROOT = 1, IROOT
               JDONE = JDONE + 1
               JGET = 17 + (2*JDONE-1)  + 2
               IF (NDONE .EQ. NTOTAL) THEN
                  JGET = 3
                  CALL GETLST(WORK(I000),JGET,1,1,1,472)
               ELSE
                  CALL GETLST(Work(I000),JGET,1,1,3,472)
               ENDIF
               write(6,"(2a,3(1x,i2))") " Reading the left vector",
     +                                  " of state : ",
     +                                   jdone,jget,jroot
               Call checksum("list-472      :",work(i000),lensz,s)
               Write(6,*)
            ENDDO
         ENDDO
      ENDDO

      Return
      End

           
      
