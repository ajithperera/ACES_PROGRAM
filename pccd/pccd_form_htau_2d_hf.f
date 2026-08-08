










      Subroutine pccd_form_htau_2d_hf(Hvo,Hov,Dpq,Dhf,Dcc,Work,Maxcor,
     +                                Nocc,Nvrt,Nbas,List_g1,List_g2,
     +                                Fact)

      Implicit Double Precision(A-H,O-Z)
      Integer D,C,Dn,Cn
      Integer B,A,Bn,An

      Dimension Work(Maxcor)
      Dimension Hvo(Nocc*Nvrt)
      Dimension Hov(Nocc*Nvrt)
      Dimension Dpq(Nbas,Nbas)
      Dimension Dhf(Nbas,Nbas)
      Dimension Dcc(Nbas,Nbas)
      Dimension Ioffo(8),Ioffv(8)

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

      Data One,Onem,Dnull,Ione,Two,Quart,Half/1.0D0,-1.0D0,0.0D0,
     +                                        1,2.0D0,0.250D0,
     +                                        0.50D0/
      Ioffo(1) = 0
      Ioffv(1) = Nocco(1)

      Ispin  = Ione
      Irrepx = Ione

      Do Irrep =2, Nirrep
         Ioffo(Irrep)=Ioffo(Irrep-1)+Pop(irrep-1,1)
         Ioffv(Irrep)=Ioffv(Irrep-1)+Vrt(irrep-1,1)
      Enddo 

      Call Dzero(Dhf,Nbas*Nbas)
      Do Irrep = 1, Nirrep
         N = Pop(Irrep,1)
         Do I = 1, N
            Dhf(Ioffo(Irrep)+I,Ioffo(Irrep)+I) =  Two
         Enddo
      Enddo
C Let us work Alpha+Beta for non-normal ordered (simply for convenince and
C as result of the fact that two-particle only contributions are also
C Alpha+Beta

      Call Dscal(Nbas*Nbas,Two,Dpq,1)

      write(6,"(2a)") "@-pccd_form_htau_2d_hf: The correlated",
     +                " density matrix"
      call output(dpq,1,nbas,1,nbas,nbas,nbas,1)

C Dco(d,c)*Dhf(k,i) is D(ck,di) 

      Ioff = Ione
      Do Irrep_di = 1, Nirrep
         Irrep_ck = Dirprd(Irrep_di,Irrepx)

         Ndi = Irpdpd(Irrep_di,11)
         Nck = Irpdpd(Irrep_ck,11)

         I010 = Ioff + Ndi*Nck
         Iend = I010 + Ndi*Nck
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph",
     +                                     Iend,Maxcor)
         Call Dzero(Work(Ioff),Ndi*Nck)
         Icount = Inull
         Do Irrep_i = 1, Nirrep
            Irrep_d = Dirprd(Irrep_i,Irrep_di)
            Do I = 1, Pop(Irrep_i,Ispin)
               Do D = 1, Vrt(Irrep_d,Ispin)
                  Do Irrep_k = 1, Nirrep
                     Irrep_c = Dirprd(Irrep_k,Irrep_ck)
                     Do K = 1, Pop(Irrep_k,Ispin)
                        Do C = 1, Vrt(Irrep_c,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Dn = D + Ioffv(Irrep_d)
                           Kn = K + Ioffo(Irrep_k)
                           Cn = C + Ioffv(Irrep_c)
                           Work(Ioff+Icount)=Dpq(Dn,Cn)*Dhf(Kn,In)*Fact 
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          Call Getlst(Work(I010),1,Ndi,2,Irrep_ck,List_g1)
          write(6,*)
          call checksum("G-in    :",Work(i010),Ndi*Nck)
          Call Daxpy(Ndi*Nck,One,Work(I010),1,Work(Ioff),1)
          Call Putlst(Work(Ioff),1,Ndi,2,Irrep_ck,List_g1)
          call checksum("G-out   :",Work(ioff),Ndi*Nck)
          Ioff = Ioff + Ndi*Nck
      Enddo 

C Dco(m,n)*Dhf(k,i) is D(nk,mi)

      Ioff = Ione
      Do Irrep_mi = 1, Nirrep
         Irrep_nk = Dirprd(Irrep_mi,Irrepx)

         Nnk = Irpdpd(Irrep_nk,14)
         Nmi = Irpdpd(Irrep_mi,14)
         I010 = Ioff + Nnk*Nmi
         Iend = I010 + Nnk*Nmi 
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph",
     +                                     Iend,Maxcor)
         Call Dzero(Work(Ioff),Nnk*Nmi)
         Icount = Inull
         Do Irrep_i = 1, Nirrep
            Irrep_m = Dirprd(Irrep_i,Irrep_mi)
            Do I = 1, Pop(Irrep_i,Ispin)
               Do M = 1, Pop(Irrep_m,Ispin)
                  Do Irrep_k = 1, Nirrep
                     Irrep_n = Dirprd(Irrep_k,Irrep_nk)
                     Do K = 1, Pop(Irrep_k,Ispin)
                        Do N  = 1, Pop(Irrep_n,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Mn = M + Ioffo(Irrep_m)
                           Kn = K + Ioffo(Irrep_k)
                           Nn = N + Ioffo(Irrep_n)
                           Work(Ioff+Icount)=Dpq(Mn,Nn)*Dhf(Kn,In)*Fact
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          Call Getlst(Work(I010),1,Nmi,2,Irrep_mi,List_g2)
          call checksum("G-in    :",Work(i010),Nmi*Nnk)
          Call Daxpy(Nmi*Nnk,One,Work(I010),1,Work(Ioff),1)
          Call Putlst(Work(Ioff),1,Nmi,2,Irrep_mi,List_g2)
          call checksum("G-out   :",Work(ioff),Nmi*Nnk)
          write(6,*)
          Ioff = Ioff + Nnk*Nmi
      Enddo

      Return
      End

 

