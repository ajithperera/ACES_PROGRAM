










      Subroutine Zero_490_lists(Work,Length,Irrepx,Iuhf)
 
      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Length)

      Integer pha_length,phb_length

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
  
      Pha_length = Irpdpd(Irrepx,9)
      Phb_length = Irpdpd(Irrepx,10)

      Max_length = Max(Pha_length,Phb_length)

      Call Dzero(Work,Max_length)
      Call Putlst(Work,1,1,1,3,490)

      If (Iuhf .Ne. 0) Then
         Call Putlst(Work,1,1,1,4,490)
      Endif

      Return
      End
