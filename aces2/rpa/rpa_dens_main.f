










      SUBROUTINE RPA_DENS_MAIN(IRREPX,ZVEC,NDIM,SCR,MAXCOR,IUHF)
C
C DRIVE THE CALCULATION OF EXCITED STATE RPA/TDA DENSITY MATRICES
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SCR(MAXCOR), ZVEC(NDIM)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)


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
C
C READ SINGLE EXCITATION VECTORS

      CALL UPDMOI(8,NFMI(1),8,91,0,0)
      CALL UPDMOI(8,NFMI(2),9,91,0,0)
      CALL UPDMOI(8,NFEA(1),8,92,0,0)
      CALL UPDMOI(8,NFEA(2),9,92,0,0)

      IOFF = 1
      DO 10 ISPIN=1,1+IUHF

       LEN=IRPDPD(IRREPX,8+ISPIN)
       IRREPD=1
       IDOO=1
       IDVV=IDOO+IRPDPD(IRREPD,20+ISPIN)
       IVEC=IDVV+IRPDPD(IRREPD,18+ISPIN)
CYCP-b
       IDOO2=IVEC
       IDVV2=IDOO2+IRPDPD(IRREPD,20+ISPIN)
       ICMO =IDVV2+IRPDPD(IRREPD,18+ISPIN)
       CALL GETREC(20, 'JOBARC', 'NBASTOT ', 1, NBAS)
       ICMO2=ICMO+NBAS*NBAS*IINTFP
       IVEC =ICMO2+NBAS*NBAS*IINTFP
CYCP-e
       ITOP=IVEC+LEN
       IF (ITOP .GT. MAXCOR) CALL INSMEM("@-RPA_DENS_MAIN",ITOP,
     &                                    MAXCOR)
       CALL ZERO(SCR,ITOP)

       IOFF = IOFF + (ISPIN-1)*IRPDPD(IRREPX,9)
C
C CALCULATE DENSITY MATRIX CONTRIBUTIONS      
C
       CALL CALC_RPA_DENS(ZVEC(IOFF),SCR(IDOO),SCR(IDVV),IRREPX,
     &                    ISPIN,IUHF)
C
C GET NATURAL TRANSITION ORBITALS (NTO)
C
       CALL GET_NTO(SCR(IDOO),SCR(IDVV),SCR(IDOO2),SCR(IDVV2),
     &      SCR(ICMO),SCR(ICMO2),IRREPX,ISPIN,IUHF)
10    CONTINUE
C
      IONE=1
      CALL PUTREC(20,'JOBARC','IOPTSYM ',IONE,IRREPX)
C
      RETURN
      END
