




















      SUBROUTINE GET_CORR_ORBS(DENS, DENS_STG1, DENS_STG2, SCR1, SCR2,
     &                         IPEA_ORBS_R, IPEA_ORBS_L, SCRA, SCRB, 
     &                         SCRC, G_OCC, G_VRT, EVALS_R, EVALS_I,   
     &                         IP_ORBS_R, IP_ORBS_L, EA_ORBS_R, 
     &                         EA_ORBS_L, WORK, MEMLEFT, LDIM1, 
     &                         LDIM2, NOCC, NVRT, NBAS, NBASX, IUHF)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C


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

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
C
      COMMON /FILES/LUOUT, MOINTS
      DOUBLE PRECISION DENS((IUHF+1)*LDIM1), G_OCC(NOCC,NOCC),
     &                 G_VRT(NVRT,NVRT), IPEA_ORBS_L(NBAS,NBAS), 
     &                 IPEA_ORBS_R(NBAS,NBAS), WORK(MEMLEFT),
     &                 SCR1(NBAS*NBAS), SCR2(NBAS*NBAS), 
     &                 SCRA(LDIM2), SCRB(LDIM2), SCRC(LDIM2), 
     &                 DENS_STG1(NBAS*NBAS), DENS_STG2(NBAS*NBAS), 
     &                 EVALS_R(NBAS), EVALS_I(NBAS), 
     &                 IP_ORBS_R(NOCC, NOCC), IP_ORBS_L(NOCC, NOCC),
     &                 EA_ORBS_R(NVRT, NVRT), EA_ORBS_L(NVRT, NVRT)
     
C
      DATA ZILCH, ONE, IONE / 0.0D0, 1.0D0, 1/
C
      ZTEST = 0.0D+00
C
C Obtain the list numbers for IP/EA STEOM of FOCK-SPACE matrix
C elements.
C
      CALL GETREC(20, 'JOBARC', 'IP_LIST', IONE, IUII_LIST)
      CALL GETREC(20, 'JOBARC', 'EA_LIST', IONE, IUAA_LIST)
       
C
C The dimensions of the IP and EA Matrices.
C       
      IRREP_LIST = 1
      NDIM_II = IRPDPD(IRREP_LIST, 21)
      NDIM_AA = IRPDPD(IRREP_LIST, 19)
C
C Read the IP/EA STEOM or FOCK space lists 
C
      CALL GETLST(SCR1, 1, 1, 1, IRREP_LIST, IUII_LIST)
      CALL GETLST(SCR2, 1, 1, 1, IRREP_LIST, IUAA_LIST)
C
      IOFF1_O_4IRREP = IONE
      IOFF2_O_4IRREP = IONE
      IOFF1_V_4IRREP = IONE
      IOFF2_V_4IRREP = IONE

      DO IRREP_R = 1, NIRREP
         NOCC_4IRREP = POP(IRREP_R, 1)
         NVRT_4IRREP = VRT(IRREP_R, 1)
         JOOF_4IRREP = IOFF1_O_4IRREP 
         JVVF_4IRREP = IOFF1_V_4IRREP
         DO JFNS = 1,  NOCC_4IRREP
             CALL DCOPY(NOCC_4IRREP, SCR1(IOFF2_O_4IRREP), 1, 
     &                 G_OCC(IOFF1_O_4IRREP, JOOF_4IRREP), 1)
             IOFF2_O_4IRREP = IOFF2_O_4IRREP + NOCC_4IRREP
              JOOF_4IRREP = JOOF_4IRREP    + 1
         ENDDO
C
         DO JFNS = 1, NVRT_4IRREP
             CALL DCOPY(NVRT_4IRREP, SCR2(IOFF2_V_4IRREP), 1, 
     &                 G_VRT(IOFF1_V_4IRREP, JVVF_4IRREP), 1)
C
          IOFF2_V_4IRREP = IOFF2_V_4IRREP + NVRT_4IRREP
             JVVF_4IRREP = JVVF_4IRREP    + 1
         END DO
C
         IOFF1_O_4IRREP = IOFF1_O_4IRREP + NOCC_4IRREP
         IOFF1_V_4IRREP = IOFF1_V_4IRREP + NVRT_4IRREP
      END DO
C
      Write(6,*)
      Write(6,*) "The G(occ,occ)" 
      call output(G_occ, 1, Nocc, 1, Nocc, Nocc, Nocc, 1)
      Write(6,*) "The G(vrt,vrt)" 
      call output(G_vrt, 1, Nvrt, 1, Nvrt, Nvrt, Nvrt, 1)
C
C Diagonalize the IP/EA matrices assuiming that they are not
C symmetric.
C     
      CALL DCOPY(NOCC*NOCC, G_OCC, 1, SCR1, 1)
      CALL XGEEV("N", "V", NOCC, G_OCC, NOCC, EVALS_R, EVALS_I, 
     &            IP_ORBS_R, NOCC, IP_ORBS_L, NOCC, WORK(1), 
     &            3*NOCC*NOCC, IERR)
      CALL SIMPLE_SORT(EVALS_R, IP_ORBS_R, NOCC)
      CALL DCOPY(NOCC*NOCC, SCR1, 1, G_OCC, 1)
      CALL XGEEV("V", "N", NOCC, G_OCC, NOCC, EVALS_R, EVALS_I,
     &            IP_ORBS_R, NOCC, IP_ORBS_L, NOCC, WORK(1), 
     &            4*NOCC*NOCC, IERR)
      CALL SIMPLE_SORT(EVALS_R, IP_ORBS_L, NOCC)
C
      Write(6,*) 
      Write(6,*) "The Principal IPs"
      Write(6, "(4F10.5)") (EVALS_R(I), I=1, NOCC)
      Write(6,*) "The right and left IP orbitals "
      CALL OUTPUT(IP_ORBS_R, 1, NOCC, 1, NOCC, NOCC, NOCC, 1)
      CALL OUTPUT(IP_ORBS_L, 1, NOCC, 1, NOCC, NOCC, NOCC, 1)
C
      IOFFE = NOCC + IONE
      CALL DCOPY(NVRT*NVRT, G_VRT, 1, SCR1, 1)
      CALL XGEEV("N", "V", NVRT, G_VRT, NVRT, EVALS_R(IOFFE), 
     &            EVALS_I(IOFFE), EA_ORBS_R, NVRT, EA_ORBS_L, 
     &            NVRT, WORK(1), 3*NVRT*NVRT, IERR)
      CALL SIMPLE_SORT(EVALS_R(IOFFE), EA_ORBS_R, NVRT)
      CALL DCOPY(NVRT*NVRT, SCR1, 1, G_VRT, 1)
      CALL XGEEV("V", "N", NVRT, G_VRT, NVRT, EVALS_R(IOFFE), 
     &            EVALS_I(IOFFE), EA_ORBS_R, NVRT, EA_ORBS_L, 
     &            NVRT, WORK(1), 4*NVRT*NVRT, IERR)
      CALL SIMPLE_SORT(EVALS_R(IOFFE), EA_ORBS_L, NVRT)
C
      Write(6,*)
      Write(6,*) "The Principal EAs"
      Write(6, "(4F10.5)") (EVALS_R(IOFFE+I-1), I=1, NVRT)
      Write(6,*) "The right and left EA orbitals "
      CALL OUTPUT(EA_ORBS_R, 1, NVRT, 1, NVRT, NVRT, NVRT, 1)
      CALL OUTPUT(EA_ORBS_L, 1, NVRT, 1, NVRT, NVRT, NVRT, 1)
      CALL ZERO(IPEA_ORBS_R, NBAS*NBAS)
      CALL ZERO(IPEA_ORBS_L, NBAS*NBAS)
      DO JFNS = 1, NOCC 
         CALL DCOPY(NOCC, IP_ORBS_R(1,JFNS), 1, IPEA_ORBS_R(1,JFNS),
     &              1)
         CALL DCOPY(NOCC, IP_ORBS_L(1,JFNS), 1, IPEA_ORBS_L(1,JFNS),
     &              1)
      ENDDO
      IOFF_V = NOCC + 1 
      JOFF_V = IOFF_V
C
      DO JFNS = 1, NVRT
         CALL DCOPY(NVRT, EA_ORBS_R(1, JFNS), 1, IPEA_ORBS_R(IOFF_V, 
     &              JOFF_V), 1)
         CALL DCOPY(NVRT, EA_ORBS_L(1, JFNS), 1, IPEA_ORBS_L(IOFF_V,
     &              JOFF_V), 1)
         JOFF_V = JOFF_V + 1
         
      END DO
C
      Write(6,*)
      Write(6,*) "The right and left IP/EA orbitals "
      CALL OUTPUT(IPEA_ORBS_R, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
      CALL OUTPUT(IPEA_ORBS_L, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
C For plotting purposes, lets put the IP/EA orbitals to the JOBARC
C file. These will then pickup by a2proc and write them in MOLDEN 
C format
C 
      CALL PUTREC(20, "JOBARC", "IPEAORBR", NBAS*NBAS*IINTFP, 
     &            IPEA_ORBS_R)
      CALL PUTREC(20, "JOBARC", "IPEAORBL", NBAS*NBAS*IINTFP, 
     &            IPEA_ORBS_L)
C
      CALL DZERO(SCR2, NBAS*NBAS)
      DO IOCC = 1, NOCC
         SCR2(NBAS*(IOCC -1) + IOCC) = 1.0D0
      ENDDO
C
      CALL MKDDEN_FULL(IPEA_ORBS_R, SCR1, DENS_STG1, SCR2, NBAS, IUHF)
C
      Write(6,*) "The Alpha Density after MKDDEN"
      call output(DENS_STG1, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
      CALL MO2AO2(DENS_STG1, DENS_STG2, SCR1, SCR2, NBAS, 1)
C
      DO I=1,NIRREP
          CALL GETBLK(DENS_STG1, DENS_STG2, NBFIRR(I),NBAS, IREPS(I))
          CALL SQUEZ2(DENS_STG2, DENS(ITRIOF(I)), NBFIRR(I))
      ENDDO
C
      RETURN
      END
