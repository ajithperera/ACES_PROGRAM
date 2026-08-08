







































































































































































































      SUBROUTINE DO_SVD_OFT2(WORK, MAXCOR, IUHF)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C This is suppose to do a singular value decompsiton of T2(ij,ab).
C  

      DIMENSION WORK(MAXCOR)



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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

C
C Lets do T2(AA), T2(BB)
C
      IRREPX = 1
      IJUNK  = 0

      DO ISPIN = 1, (IUHF+ 1)
        
         DO IRREP = 1, NIRREP

            IRREPR = IRREP
            IRREPL = DIRPRD(IRREPX,IRREPR)

            NCOLS     = IRPDPD(IRREPR, 2+ISPIN)
            NROWS     = IRPDPD(IRREPL, ISPIN)
            NCOLS_EXP = IRPDPD(IRREPR, 20+ISPIN)
            Write(*,*) "NCOLS, NROWS=", NCOLS, NROWS 
            NDIM      = NCOLS*NROWS 
            NDIM_EXP  = NCOLS_EXP*NROWS

            NVRT_PAIRS = VRT(IRREPL, ISPIN) * (VRT(IRREPL, ISPIN)-1)/2
            NOCC_PAIRS = POP(IRREPR, ISPIN) * POP(IRREPR, ISPIN)

            I000 = 1
            I010 = I000 + NDIM_EXP*IINTFP
            I020 = I010 + NDIM_EXP*IINTFP
            I030 = I020 + NOCC_PAIRS*IINTFP
            I040 = I030 + NOCC_PAIRS*IINTFP
C
C Read T2(A<B, I<J) for alpha and beta and expand and then transpose.
C
            CALL GETLST(WORK(I000), 1, NCOLS, 1, IRREPR,  60+ISPIN)
            CALL SYMEXP(IRREPR, POP(IRREPR, ISPIN), NROWS, WORK(I000))
            CALL TRANSP(WORK(I000), WORK(I010), NCOLS_EXP, NROWS)

 
C Create a Set of (I,J) matrices for a pair of AB and digaonalize

            IOFF = 0
            DO IPAIR = 1, NVRT_PAIRS

               CALL DCOPY(NOCC_PAIRS, WORK(I010 + IOFF), 1, 
     &                                       WORK(I020), 1)
      n_occ=pop(irrepr,ispin)
      call output(work(I020), 1,n_occ,1,n_occ,n_occ,n_occ,0)
               IOFF = IPAIR*NOCC_PAIRS 

               CALL EIG(WORK(I020), WORK(I030), IJUNK, 
     &                  POP(IRREP,ISPIN), 0)

            ENDDO
        ENDDO
      ENDDO
      STOP

C T2(AB) blocks
       
      DO IRREP = 1, NIRREP

         IRREPR = IRREP
         IRREPL = DIRPRD(IRREPX,IRREPR)

         NCOLS     = IRPDPD(IRREPR, 14)
         NROWS     = IRPDPD(IRREPL, 15)
      
         NDIM      = NCOLS*NROWS

         NVRT_PAIRS = VRT(IRREP, ISPIN) * VRT(IRREP, ISPIN)
         NOCC_PAIRS = POP(IRREP, ISPIN) * POP(IRRPE, ISPIN)
      
         I000 = 1
         I010 = I000 + NDIM*IINTFP
         I020 = I010 + NDIM*IINTFP
         I030 = I020 + NOCC_PAIRS*IINTFP
         I040 = I030 + NOCC_PAIRS*IINTFP
C     
C Read T2(Ab, Ij) and transpose it,
C
         CALL GETLST(WORK(I000), 1, NCOLS, 1, IRREPR,  63)
         CALL TRANSP(WORK(I000), WORK(I010), NCOLS, NROWS)
C     
C Create a Set of (I,J) matrices for a pair of AB and digaonalize
C     
         IOFF = 0
         DO IPAIR = 1, NVRT_PAIRS
         
            CALL DCOPY(NOCC_PAIRS, WORK(I010 + IOFF), 1,
     &                                WORK(I020), 1)
            IOFF = IPAIR*NOCC_PAIRS 
        
            CALL EIG(WORK(I020), WORK(I030), IJUNK,
     &               POP(IRREP,ISPIN), 0)

      write(*,*)
      write(6,"(a)") "The egienvlaues of the OO block per VV pair"
      write(6,"(6(1x,F15.7))") (Work(i020+(i-1)*pop(irrepr,ispin)+
     &                         (i-1)),i=1,pop(irrepr,ispin))
         ENDDO

      ENDDO 

      RETURN
      END
      






