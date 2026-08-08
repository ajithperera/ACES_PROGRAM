













































































































































































































      Subroutine Psi4dbg_dens(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)
      Logical Non_hf



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
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      Data Ione /1/

      I0 = Ione
      Non_hf = (Iflags(38) .gt.0)
      Call Psi4dbg_init2pdm(Work(I0),Maxcor)
      Call Psi4dbg_gamdrv(work(I0),Maxcor,Iuhf)

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
      CALL Psi4dbg_densoo(Work(IDOO),Work(I0),Maxcor,Iuhf)
      CALL Psi4dbg_densvv(Work(IDVV),Work(I0),Maxcor,Iuhf)
      If (Non_hf) CALL Psi4dbg_densvo(Work(IDOV),Work(I0),
     +                                Maxcor,Iuhf)

C G(ij,ab) + G(ab,ij)
      write(6,*)
      write(6,"(a)") "Checksums of G(ij,ab)"
       Call Pccd_Gamma1(Work(I0),Maxcor,Iuhf)
      
       Call Psi4dbg_setg(Work(I0),Maxcor)

       Write(6,*)
       Write(6,*) "The checksums of the Density blocks"
       Call checksum("DENSOO ", Work(IDOO),Nfmi(1)+IUHF*Nfmi(2))
       Call checksum("DENSVV ", Work(IDVV),Nfea(1)+IUHF*Nfea(2))
       If (non_hf) Call checksum("DENSOV ", Work(IDOV),
     +                            Nt(1)+IUHF*Nt(2))

      Return
      End 

