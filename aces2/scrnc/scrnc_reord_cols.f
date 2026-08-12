










      Subroutine Scrnc_reord_cols(G2,G3,Nbfns,Irrep_ao,Iord)

      Implicit Double Precision(A-H,O-Z)

      Dimension G2(Nbfns*Nbfns,Nbfns*Nbfns)
      Dimension G3(Nbfns*Nbfns,Nbfns*Nbfns)
      Dimension Irrep_ao(Nbfns*Nbfns),Iordr_ao(Nbfns*Nbfns)

c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

      Nbfns2 = Nbfns * Nbfns 

      Ithru = 0
      Do Irrepr = 1, Nirrep
         Do I = 1, Nbfns2 
            If (Irrep_ao(I) .eq. Irrepr) Then
                Ithru = Ithru + 1
                Iordr_ao(Ithru) = I
            Endif
         Enddo
      Enddo

      write(6,"(a,10(1x,I4))")"Reorder array      :", (Iordr_ao(i),
     +                       i=1,Nbfns2)
      Do I = 1, Nbfns2 
         Iold = Iordr_ao(I)
         Inew = I
            Do J = 1, Nbfns2
               G3(J,Inew) = G2(J,Iold)
            Enddo
      Enddo 
      Call Dcopy(Nbfns2*Nbfns2,G3,1,G2,1)

      Do I = 1, Nbfns2-1
         Do J = I+1, Nbfns2 
            If (Irrep_ao(I) .Gt. Irrep_ao(J)) Then
               Ikeep = Irrep_ao(J)
               Irrep_ao(J) = Irrep_ao(I)
               Irrep_ao(I) = Ikeep
            Endif 
         Enddo 
      Enddo

      Write(6,*)
      Write(6,"(a)") "The I<bra|ket>"
      Call output(G2,1,Nbfns2,1,Nbfns2,Nbfns2,Nnbfns2,1)
    
      Return
      End
