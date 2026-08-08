










      Subroutine Pccd_gossip(Grd,Grd_stata,Grd_statb,Sgrad_stata,
     +                       Sgrad_statb,Nbas,Ispin)

      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry 

      Dimension Grd_stata(6)
      Dimension Grd_statb(6)
      Dimension SGrd_stata(6)
      Dimension SGrd_statb(6)



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

      Common /Symm/Symmetry

      If (Ispin .Eq. 1) Then
          Call Pccd_vstat(Grd,Sgrd_stata,Nbas*Nbas)
          Call Putrec(20,"JOBARC","SGRSTATA",6*Iintfp,Sgrd_stata)
      Elseif (Ispin .Eq. 2) Then
          Call Pccd_vstat(Grd,Sgrd_statb,Nbas*Nbas)
          Call Putrec(20,"JOBARC","SGRSTATB",6*Iintfp,Sgrd_statb)
      Endif

      Write(6,"(a)") "Scalled gradients (K=-H^-1g)"
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
      Write(6,*)
      If (Ispin .Eq. 1) Write(6,"(2a)") " The statistics of the scaled",
     +                  " alpha orbital rotation gradient matrix"
      If (Ispin .Eq. 2) Write(6,"(2a)") " The statistics of the scaled",
     +                  " beta  orbital rotation gradient matrix"
      Write(6,"(1x,2a)") "--------------------------------------------",
     +                   "-------"
      Write(6,"(5x,a,5xa,5xa)") "Minimum grad.", "Maximum grad.",
     +                           "RMS grad"
      If (Ispin .Eq. 1) Then
          Write(6,*)
          Write(6,"(3(5x,E12.6))") SGrd_stata(3), SGrd_stata(4),
     +                             SGrd_stata(5)
      Elseif (Ispin .Eq. 2) Then
          Write(6,*)
          Write(6,"(3(5x,E12.6))") SGrd_statb(3), SGrd_statb(4),
     +                             SGrd_statb(5)
      Endif
      Write(6,*)
      Write(6,"(1x,2a)") "--------------------------------------------",
     +                   "-------"

      Return
      End 
