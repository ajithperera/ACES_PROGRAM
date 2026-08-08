










      SUBROUTINE INITNUMR_IP(SCR, MAXCOR, IUHF)
C
C IF THE NUMBER OF ROOTS TO BE DETERMINED IS NOT GIVEN (I.E. ALL ARE ZERO)
C THEN THE NUMBER OF ROOTS TO BE DETERMINED IN THE IP EOM CALCULATION IS 
C DETERMINED FROM INSPECTION OF THE SCF ENERGY EIGENVALUES. POSSIBLE 
C DEGENERACIES BETWEEN BLOCKS IS ACCOUNTED FOR.
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR, IP(8,2), IPMIN, DIFFE
      DIMENSION SCR(MAXCOR),IOFFVRT(8,2), IOFFPOP(8,2)
C
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


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



c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      COMMON/IPCALC/LEFTHAND,SINGONLY, DROPCORE
      COMMON/IpINFO/NUMROOT(8,3)
C
      IROOT = 0
      DO IRREP = 1, NIRREP
         DO ISPIN = 1, 1+IUHF
            IROOT = IROOT + NUMROOT(IRREP, ISPIN)
         ENDDO
      ENDDO
      IF (IROOT .GT. 0) THEN
         RETURN
      ELSE
         WRITE(6,*)
         WRITE(6,"(a,a)") ' NUMBER OF DESIRED ROOTS',
     $                    ' IS ESTIMATED FROM ORBITAL EIGENVALUES'
         WRITE(6,*)
      ENDIF
C
      NBASA = NOCCO(1) + NVRTO(1)
      IF (IUHF.NE.0) THEN
         NBASB = NOCCO(2)  + NVRTO(2)
      ELSE
         NBASB = 0
      ENDIF

      I000 = 1
      I010 = I000 + NBASA
      I020 = I010 + NBASB
      CALL GETREC(20, 'JOBARC', 'SCFEVALA', IINTFP*NBASA, SCR(I000))

         write(6,*) ' Hartree Fock orbital energies : alfa'
         call output(scr(i000), 1, 1, 1, nbasa, 1, nbasa, 1)
      IF (IUHF. NE. 0) THEN
         CALL GETREC(20, 'JOBARC', 'SCFEVALB', IINTFP*NBASB, SCR(I010))

         write(6,*) ' Hartree Fock orbital energies : beta'
         call output(scr(i010), 1, 1, 1, nbasb, 1, nbasb, 1)
      ENDIF
C
C  CALCULATE OFFSETS IN ORBITAL ENERGY ARRAYS
C
      IOFFPOP(1,1) = I000
      IOFFPOP(1,2) = I010

      DO IRREP = 2, NIRREP
         DO ISPIN = 1, 1 + IUHF
            IOFFPOP(IRREP,ISPIN) = IOFFPOP(IRREP-1,ISPIN) +
     $                             POP(IRREP-1, ISPIN)
         ENDDO 
      ENDDO 

      IOFFVRT(1,1) = IOFFPOP(NIRREP,1) + POP(NIRREP,1)

      IF (IUHF .NE. 0) THEN
         IOFFVRT(1,2) = IOFFPOP(NIRREP,2) + POP(NIRREP,2)
      ENDIF

      DO IRREP = 2, NIRREP
         DO ISPIN = 1, 1 + IUHF
            IOFFVRT(IRREP,ISPIN) = IOFFVRT(IRREP-1, ISPIN) +
     $                             VRT(IRREP-1, ISPIN)
         ENDDO 
      ENDDO 
C
C DETERMINE MINIMUM SCF IP PER IRREP
C
      DO IRREP = 1, NIRREP
         DO ISPIN = 1, 1 + IUHF
            IF (POP(IRREP,ISPIN) .NE. 0) THEN
            IP(IRREP,ISPIN) = DABS(SCR(IOFFPOP(IRREP,ISPIN)+
     &                                         POP(IRREP,ISPIN)-1))
            ENDIF 
         ENDDO 
      ENDDO 

         write(6,*)
         WRITE(6,*) '  Minumm SCF Ionization potentials'
         CALL OUTPUT(IP,1,nirrep,1,2,8,2,1)
C
C DETERMINE minimum Ionization Potential.
C
      IPMIN = IP(1,1)
      DO  ISPIN = 1, 1 + IUHF
         DO  IRREP = 1, NIRREP
            IF (POP(IRREP,ISPIN).GT.0) IPMIN=MIN(IPMIN,IP(IRREP,ISPIN))
         ENDDO 
      ENDDO
     
      write(6,*)
      Write(6,"(a,F10.6)") " IPMIN =",IPMIN
C      
      DO  ISPIN =1, 1+IUHF
         DO  IRREP = 1, NIRREP
             NUMROOT(IRREP,ISPIN) = 0
             IF (POP(IRREP,ISPIN).GT.0) THEN
                 IF (ABS(IPMIN-IP(IRREP,ISPIN)).LT.1.0D0) THEN
                     NUMROOT(IRREP,ISPIN) = 1
                 ENDIF
            ENDIF
         ENDDO
      ENDDO
C
C CHECK FOR DEGENERACIES
C
      DO ISPIN = 1, 1+IUHF
         DO IRREP = 1, NIRREP - 1
            IF (NUMROOT(IRREP,ISPIN) .NE. 0) THEN
               DO JRREP = IRREP + 1, NIRREP
                  IF ((NUMROOT(JRREP,ISPIN) .NE. 0) .AND.
     +               (POP(IRREP,ISPIN) .EQ. POP(JRREP,ISPIN)) .AND.
     +               (VRT(IRREP,ISPIN) .EQ. VRT(JRREP,ISPIN))) THEN
                     DIFFE = 0.0D0
                     DO I = 0, POP(IRREP,ISPIN) - 1
                     DIFFE = DIFFE + ABS(SCR(IOFFPOP(IRREP,ISPIN) + I) -
     +                                   SCR(IOFFPOP(JRREP,ISPIN) + I))
                     ENDDO
                     DO I = 0, VRT(IRREP,ISPIN) - 1
                     DIFFE = DIFFE + ABS(SCR(IOFFVRT(IRREP,ISPIN) + I) -
     $                               SCR(IOFFVRT(JRREP,ISPIN) + I))
                     ENDDO
                     IF (ABS(DIFFE) .LT. 0.0001)  THEN
                        NUMROOT(JRREP,ISPIN) = 0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
