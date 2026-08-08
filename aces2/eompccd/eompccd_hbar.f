










      Subroutine Eompccd_hbar(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

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
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Data One,Onem,Ione,Halfm /1.0D0,-1.0D0,1,-0.50D0/
      
C Note that in pccd_ldriver calls to pccd_formwl, pccd_wtwtw,
C pccd_formw3 and pccd_formw5 formed Hbar(mb,ej),Hbar(ia,jk),
C Hbar(ab,ci). 

C For ccd and its approximate methods like pCCD:j
C Hbar(ij,ka) is the same as W(ij,ka) and there is no change 
C to the W lists. Similarly there is no change to W(ai,bc) and
C Hbar(ai,bc) is W(ai,bc), Also, the F(ae) and F(mi) needs 
C no changes. 

C Construct the redundent Hbar(mb,ej) AAAA from ABAB and ABBA.
      Irrepx   = Ione
      Nsize_54 = Idsymsz(Irrepx,9,9)
      Nsize_56 = Idsymsz(Irrepx,9,10)
      Nsize_58 = Idsymsz(Irrepx,11,11)

      I000 = Ione
      I010 = I000 + Nsize_56
      Iend = I010 + Nsize_58
      If (Iend .Gt. Maxcor) Call Insmem("eom_pccd_hbar",Iend,Maxcor)
      Call Getall(Work(I000),Nsize_56,Irrepx,56)
      Call Getall(Work(I010),Nsize_58,Irrepx,58)
      Call Daxpy(Nsize_54,Onem,Work(I010),1,Work(I000),1)
      Call Dscal(Nsize_54,Halfm,Work(I000),1)
      Call Putall(Work(I000),Nsize_54,Irrepx,54)

      Write(6,"(a)") "The OO-CCD/pCCD Hbar elements"
      Call checkhbar(Work,Maxcor,Iuhf) 

      Return
      End 
