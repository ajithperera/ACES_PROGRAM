      Subroutine Pccd_dens(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)
      LOGICAL PCCD,CCD,LCCD



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
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

      COMMON/CALC/PCCD,CCD,LCCD
      COMMON/NHFREF/NONHF_REF
 
      Data Ione,Onem /1,-1.0D0/

      I0 = Ione
       call pccd_checkintms(work(I0),maxcor,iuhf,0)
C G(ij,kl),G(ab,cd),G(ai,bj) is constructed in gamdrv. The other
C Gamma elements are zero for pCCD.

      Call Pccd_gamdrv(Work(I0),Maxcor,Iuhf)

      IDOO   = I0 + MAXCOR-(Nfmi(1)+IUHF*Nfmi(2))
      MAXCOR = MAXCOR - (Nfmi(1)+IUHF*Nfmi(2))
C
      IDVV   = IDOO-(Nfea(1)+IUHF*Nfea(2))
      MAXCOR = MAXCOR-(Nfea(1)+IUHF*Nfea(2))
C
      IDOV   = IDVV-(Nt(1)+IUHF*Nt(2))
      MAXCOR = MAXCOR-(Nt(1)+IUHF*NT(2))
C
      Call Dzero(Work(IDOO),Nfmi(1)+IUHF*Nfmi(2))
      Call Dzero(Work(IDVV),Nfea(1)+IUHF*Nfea(2))
      Call Dzero(Work(IDOV),Nt(1)+IUHF*Nt(2))
      CALL Pccd_densoo(Work(IDOO),Work(I0),Maxcor,Iuhf)
      CALL Pccd_densvv(Work(IDVV),Work(I0),Maxcor,Iuhf)
      If (Nonhf_ref) Call Pccd_densvo(Work(IDOV),Work(I0),Maxcor,
     +                                Iuhf)

C G(ij,ab) + G(ab,ij)
      write(6,*)
      write(6,"(a)") "Checksums of G(ij,ab)"
      Call Pccd_gamma1(Work(I0),Maxcor,Iuhf)
CSSS      If (Lccd) Call Pccd_setg(Work(I0),Maxcor)

       Write(6,*)
       Write(6,"(a)") "The checksums of the Density blocks"
CSSS       Call Dscal(Nfmi(1)+IUHF*Nfmi(2),Onem, Work(IDOO),1)
CSSS       Call Dscal(Nfea(1)+IUHF*Nfea(2),Onem, Work(IDVV),1)
       Call checksum("DENSOO   :", Work(IDOO),Nfmi(1)+IUHF*Nfmi(2))
       Call checksum("DENSVV   :", Work(IDVV),Nfea(1)+IUHF*Nfea(2))
       Call checksum("DENSVO   :", Work(IDOV),Nt(1)+IUHF*Nt(2))

      Return
      End 

