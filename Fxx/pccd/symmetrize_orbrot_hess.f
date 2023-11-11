      Subroutine Symmetrize_orbrot_hess(W,WT,WA,Nrow,Ncol,Maxd,Irrepl,
     +                                  Irrepr,Type)
  
      Implicit Double Precision(A-H,O-Z)
      Dimension W(Nrow,Ncol),WT(Ncol,Nrow)
      Dimension Wa(Maxd)
      Character*2 Type

      Data Ione /1/

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

      Ispin = Ione 
      Do Icol = 1, Ncol
         Ioff  = Ione
         Do Irrep_p = 1, Nirrep 
            Irrep_q = Dirprd(Irrep_p,Irrepl)
            If (Type .Eq. "PP") Then
               Np = Vrt(Irrep_p,Ispin)
               Nq = Vrt(Irrep_q,Ispin)
            Else 
               Np = Pop(Irrep_p,Ispin)
               Nq = Pop(Irrep_q,Ispin)
            Endif
            Call pccd_asymm(W(Ioff,Icol),WA,Np,Nq)
            Call Dcopy(Np*Nq,Wa,1,W(Ioff,Icol),1)
            Ioff = Ioff + Np*Nq
         Enddo  
      Enddo 

      Call Transp(W,WT,Ncol,Nrow)

      Do Irow = 1, Nrow
         Ioff  = Ione
         Do Irrep_p = 1, Nirrep 
            Irrep_q = Dirprd(Irrep_p,Irrepr)
            If (Type .Eq. "PP") Then
               Np = Vrt(Irrep_p,Ispin)
               Nq = Vrt(Irrep_q,Ispin)
            Else 
               Np = Pop(Irrep_p,Ispin)
               Nq = Pop(Irrep_q,Ispin)
            Endif
            Call pccd_asymm(WT(Ioff,Irow),Wa,Np,Nq)
            Call Dcopy(Np*Nq,Wa,1,WT(Ioff,Irow),1)
            Ioff = Ioff + Np*Nq
         Enddo  
      Enddo 

      Call Transp(WT,W,Nrow,Ncol)

      Return
      End
