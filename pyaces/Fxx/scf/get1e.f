






































































































































































































      SUBROUTINE GET1E(ONEH,OVRLP,BUF,IBUF,IBAS,KINETIC,LDIM,ILNBUF,
     &                 REPULS, NBAS,NIRREP,NBFIRR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION KINETIC
      Logical ECP, Response
C
      PARAMETER(LUINT=10, LUNITSE=25)
C
      DIMENSION ONEH(LDIM),OVRLP(LDIM),KINETIC(LDIM)
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF),IBAS(NBAS)
C
c molcas.com : begin
      logical seward, petite_list
      character*8 fnord
      integer luord
      integer ipmat(8,2), lbbt, lbbs
      common /molcas_com/ seward, petite_list, fnord, luord,
     &                    ipmat, lbbt, lbbs
c molcas.com : end


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



      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)

c   o pick up the Hamiltonian 1-e elements and the overlap matrix
      IF (SEWARD) THEN
         REWIND LUNITSE
         CALL LOCATE(LUNITSE,'Kinetic ')
         CALL ZERO(ONEH,LDIM)
         NUT = 0
         DO WHILE (NUT.NE.-1)
            READ(LUNITSE) BUF, IBUF, NUT
            DO INT = 1, NUT
               ONEH(IBUF(INT)) = BUF(INT)
            END DO
         END DO
         REWIND LUNITSE
         CALL LOCATE(LUNITSE,'Attract ')
         NUT = 0
         DO WHILE (NUT.NE.-1)
            READ(LUNITSE) BUF, IBUF, NUT
            DO INT = 1, NUT
               ONEH(IBUF(INT)) = ONEH(IBUF(INT)) + BUF(INT)
            END DO
         END DO
         REWIND LUNITSE
         CALL LOCATE(LUNITSE,'Mltpl  0')
         CALL ZERO(OVRLP,LDIM)
         NUT = ILNBUF
         DO WHILE (NUT.EQ.ILNBUF)
            READ(LUNITSE) BUF, IBUF, NUT
            DO INT = 1, NUT
               OVRLP(IBUF(INT)) = BUF(INT)
            END DO
         END DO
         CLOSE(LUNITSE,STATUS='KEEP')

      ELSE
c        (in other words: VMOL)
C The kinetic energies are required in wave function analysis to
C compute virial coefficients. These integrlas are use in a2proc
C rotine a2make_wfn.F. Ajith Perera, 08/2010. 
C
         CALL LOCATE(LUINT,'ONEHAMIL')
         CALL ZERO(ONEH,LDIM)
         NUT = ILNBUF
         DO WHILE (NUT.EQ.ILNBUF)
            READ(LUINT) BUF, IBUF, NUT
            DO INT = 1, NUT
               ONEH(IBUF(INT)) = BUF(INT)
            END DO
         END DO
C
C The following calls are experimental. This is added for Johanness work
C on July-Augutst, 2016. Should not interfere with others work since 
C he only can write the SEMIPMCR record. Ajith Perera
C
        CALL Getrec(0,"JOBARC","SEMIPMCR",Length,KINETIC)
        If (Length .Gt. 0) then 
            Write(6,"(a,a)") " Adding a semi-empirical correction",
     &                       " to the Fock matrix"
           CALL Getrec(20,"JOBARC","SEMIPMCR",LDIM,KINETIC)

           Do Int = 1, Ldim
              ONEH(INT) = ONEH(INT) + Kinetic(INT)
           Enddo 
        Endif 

         CALL LOCATE(LUINT,'KINETINT')
         CALL ZERO(KINETIC,LDIM)
         NUT = ILNBUF
         DO WHILE (NUT.EQ.ILNBUF)
            READ(LUINT) BUF, IBUF, NUT
            DO INT = 1, NUT
               KINETIC(IBUF(INT)) = BUF(INT)
            END DO
         END DO
C
C Read and process the ECP contrbutions to one electron integrals
C and add those to ONEH.
C
         ECP = (Iflags(71) .EQ. 1)

         If (ECP) Call Process_ecp(Ovrlp, Oneh, Ldim, Nbas, Nirrep,
     &                             Nbfirr)

         CALL LOCATE(LUINT,'OVERLAP ')
         CALL ZERO(OVRLP,LDIM)
         NUT = ILNBUF
         DO WHILE (NUT.EQ.ILNBUF)
            READ(LUINT) BUF, IBUF, NUT
            DO INT = 1, NUT
               OVRLP(IBUF(INT)) = BUF(INT)
            END DO
         END DO

c     END IF (SEWARD)
      END IF

c   o add finite field, if requested
      Response = Iflags(19) .EQ. 1

      IF (.NOT. response .AND.
     &     (IFLAGS(23).NE.0.OR.IFLAGS(24).NE.0.OR.IFLAGS(25).NE.0)) THEN
         CALL FINFLD(ONEH,BUF,IBUF,IBAS,REPULS,NBAS,NIRREP,NBFIRR,
     &               ILNBUF)
      END IF
c   o calculate RFSCF parameters
      IF (IFLAGS2(118).NE.0) CALL RFIELD1

      RETURN
      END

