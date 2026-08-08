










      Subroutine Pccd_scale_grads(Grd_oo,Grd_vv,Grd_vo,Grd_ov,Hes_vo,
     +                            Hes_ov,Hes_oo,Hes_vv,Lenoo,Lenvv,
     +                            Lenvo,Nocc,Nvrt,Icycle)

      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry 
      Logical Energy_exist

      Dimension Grd_oo(Lenoo)
      Dimension Grd_vv(Lenvv)
      Dimension Grd_vo(Lenvo)
      Dimension Grd_ov(Lenvo)
      Dimension Hes_ov(Lenvo)
      Dimension Hes_vo(Lenvo)
      Dimension Hes_oo(Lenoo)
      Dimension Hes_vv(Lenvv)

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




      Common /Symm/Symmetry 

      Data Onem,Izero,Thres,Half,One /-1.0D0,0,1.0D-07,0.50D0,
     +                                 1.0D0/
       Write(6,*)
       Write(6,"(a)") "Printing from Pccd_scale grads:@-entry"
       Write(6,"(a)") "The OO-MO gradients"
       call output(Grd_oo,1,Nocc,1,Nocc,Nocc,Nocc,1)
       Write(6,"(a)") "The VV-MO gradients"
       call output(Grd_vv,1,Nvrt,1,Nvrt,Nvrt,Nocc,1)
       Write(6,"(a)") "The VO-MO gradients"
       call output(Grd_vo,1,Nvrt,1,Nocc,Nvrt,Nocc,1)
       Write(6,"(a)") "The OV-MO gradients"
       call output(Grd_ov,1,Nocc,1,Nvrt,Nocc,Nvrt,1)
       Write(6,"(a)") "The VO-MO Hessians"
       call output(Hes_vo,1,Nvrt,1,Nocc,Nvrt,Nocc,1)
       Write(6,"(a)") "The OV-MO Hessians"
       call output(Hes_ov,1,Nocc,1,Nvrt,Nocc,Nvrt,1)
       Call Getrec(20,"JOBARC","TOTENERG",Iintfp,E_new)
       Inquire(File="Energy",Exist=Energy_exist)
       If (Energy_exist) Then
           Open(Unit=15,File="Energy",Form="Formatted",Status="old",
     +          Position="Append")
       Else 
           Open(Unit=15,File="Energy",Form="Formatted",Status="new")
       Endif 
       If (Icycle .Eq. 0) Then
           Direction = Onem
           Write(15,*) E_new, Direction 
       Endif 

       If (Icycle .Gt. 0) Then
           Backspace(Unit=15)
           Read(15,*) E_old,Direction 
       Endif 
       If (Icycle .Gt. 0) Then
          If (E_new .Gt. E_old) Then
              Print*, E_new, E_old, Direction
              Write(6,"(a,a)") "The energy is going up: changing the",
     +                         " direction of search"
              If (Direction .Eq. Onem) Then
                  Direction=One
              Elseif (Direction .Eq. One) Then
                  Direction=Onem
              Endif 
          Endif 
       Endif 
       Write(15,*) E_new, Direction 
       Close(Unit=15)

C First I thought that changing the sign of the gradient if the energy
C goes up is useful, but it turns out not to be the case and that is why
C Direction is -1 (from G(k+1) = -H^(-1)G(k) (G(k) is the current gradient)

       Fact = Half 
       Direction = Onem
       Do I = 1, Lenvo
          If (Hes_vo(I) .Gt. Thres) Hes_vo(I) = Direction/Hes_vo(I)
          If (Hes_ov(I) .Gt. THres) Hes_ov(I) = Direction/Hes_ov(I)
       Enddo
       Do I = 1, Lenvo
          Grd_vo(I) = Grd_vo(I)*Hes_vo(I)*Fact
          Grd_Ov(I) = Grd_ov(I)*Hes_ov(I)*Fact
       Enddo

       Write(6,*)
       Write(6,"(a)") "Printing from Pccd_scale grads:@-exit"
       Write(6,"(a)") "The OO-MO gradients"
       call output(Grd_oo,1,Nocc,1,Nocc,Nocc,Nocc,1)
       Write(6,"(a)") "The VV-MO gradients"
       call output(Grd_vv,1,Nvrt,1,Nvrt,Nvrt,Nocc,1)
       Write(6,"(a)") "The VO-MO gradients"
       call output(Grd_vo,1,Nvrt,1,Nocc,Nvrt,Nocc,1)
       Write(6,"(a)") "The OV-MO gradients"
       call output(Grd_ov,1,Nocc,1,Nvrt,Nocc,Nvrt,1)

      Return
      End 
