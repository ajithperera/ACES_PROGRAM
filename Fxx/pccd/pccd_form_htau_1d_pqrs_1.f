      Subroutine Pccd_form_htau_1d_pqrs_1(Htau_pq,Htau_oo,Htau_vv,
     +                                    Htau_ov,Htau_vo,Ioo,Ivv,
     +                                    Work,Ioffo,Ioffv,Maxcor,
     +                                    Nocc,Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Ioo,Ivv
      Integer a,b,c,d,an,bn,cn,dn
      Logical ONEP_ONLY
     
      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_oo(Nocc,Nocc)
      Dimension Htau_vv(Nocc,Nocc)
      Dimension Htau_ov(Nocc,Nocc)
      Dimension Htau_vo(Nocc,Nocc)
      Dimension Ioo(Nocc,Nocc)
      Dimension Ivv(Nvrt,Nvrt)
      Dimension Work(Maxcor)
      Dimension Ioffo(8)
      Dimension Ioffv(8)
 
      COMMON/ORBR_HESS/ONEP_ONLY

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

      Data One,Ione,Inull /1.0D0,1,0/

      Call Dzero(Ioo,Nocc*Nocc)
      Call Dzero(Ivv,Nvrt*Nvrt)

      Ioffo(1) = 0
      Ioffv(1) = Nocco(1)
       
      Do Irrep =2, Nirrep
         Ioffo(Irrep)=Ioffo(Irrep-1)+Pop(irrep-1,1)
         Ioffv(Irrep)=Ioffv(Irrep-1)+Vrt(irrep-1,1)
      Enddo 

      Do I = 1, Nocc
         Ioo(I,I) = One
      Enddo 
      Do I = 1, Nvrt
         Ivv(I,I) = One
      Enddo 
     
      Irrepx = Ione
      Ioff   = Ione
      Ispin  = Ione 

C OOOO block (H(ij,kl)=-delta(ik)(F(jl)+D(l,j)) and all permutations)

      List_h = 213 
      Do Irrep_kl = 1, Nirrep
         Irrep_ij = Dirprd(Irrep_kl,Irrepx)

         Nij = Irpdpd(Irrep_ij,14)
         Nkl = Irpdpd(Irrep_kl,14)
         I000 = Ione
         I010 = I000 + Nij*Nkl
         Iend = I010 + Nij*Nkl
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nij*Nkl)
         Icount = Inull
         Do Irrep_l = 1, Nirrep
            Irrep_k = Dirprd(Irrep_l,Irrep_kl)
            Do l = 1, Pop(Irrep_l,Ispin)
               Do k = 1, Pop(Irrep_k,Ispin)
                  Do Irrep_j = 1, Nirrep
                     Irrep_i = Dirprd(Irrep_j,Irrep_ij)
                     Do j = 1, Pop(Irrep_j,Ispin)
                        Do i = 1, Pop(Irrep_i,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           Kn = K + Ioffo(Irrep_k)
                           Ln = L + Ioffo(Irrep_l)
                           Work(Ioff+Icount) =  
     +                              -Ioo(i,j)*
     +                               (Htau_pq(kn,ln)+Htau_pq(ln,kn)) 
     +                              -Ioo(i,l)*
     +                               (Htau_pq(kn,jn)+Htau_pq(kn,jn))
     +                              +Ioo(j,k)*
     +                               (Htau_pq(in,ln)+Htau_pq(ln,in))
     +                              +Ioo(i,l)*
     +                               (Htau_pq(jn,kn)+Htau_pq(kn,jn))
                          Icount = Icount + 1
                        Enddo 
                      Enddo
                  Enddo                   
               Enddo
            Enddo 
          Enddo
          If (.Not. ONEP_ONLY) Then 
             Call Getlst(Work(I010),1,Nkl,1,Irrep_kl,List_h)
             Call Daxpy(Nij*Nkl,One,Work(I010),1,Work(Ioff),1)
          Endif
          Call Putlst(Work(Ioff),1,Nkl,1,Irrep_kl,List_h)
          Call checksum("OOOO    :",Work(Ioff),Nij*Nkl)
      Enddo 

C VVVV block (H(ab,cd)=-delta(ac)(F(bd)+F(db)) and all permutations)

      Ioff   = Ione
      List_h = 234
      Do Irrep_cd = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_cd,Irrepx)

         Nab = Irpdpd(Irrep_ab,15)
         Ncd = Irpdpd(Irrep_cd,15)
         I000 = Ione
         I010 = I000 + Nab*Ncd
         Iend = I010 + Nab*Ncd
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nab*Ncd)
         Icount = Inull 
         Do Irrep_d = 1, Nirrep
            Irrep_c = Dirprd(Irrep_d,Irrep_cd)
            Do d = 1, Vrt(Irrep_d,Ispin)
               Do c = 1, Vrt(Irrep_c,Ispin)
                  Do Irrep_b = 1, Nirrep
                     Irrep_a = Dirprd(Irrep_b,Irrep_ab)
                     Do b = 1, Vrt(Irrep_b,Ispin)
                        Do a = 1, Vrt(Irrep_a,Ispin)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Cn = C + Ioffv(Irrep_c)
                           Dn = D + Ioffv(Irrep_d)
                           Work(Ioff+Icount) =
     +                              -Ivv(a,b)*
     +                               (Htau_pq(cn,dn)+Htau_pq(dn,cn))
     +                              -Ivv(c,d)*
     +                               (Htau_pq(an,bn)+Htau_pq(bn,an))
     +                              +Ivv(b,c)*
     +                               (Htau_pq(an,dn)+Htau_pq(dn,an))
     +                              +Ivv(a,d)*
     +                               (Htau_pq(bn,cn)+Htau_pq(cn,bn))
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          If (.Not. ONEP_ONLY) Then
             Call Getlst(Work(I010),1,Ncd,1,Irrep_cd,List_h)
             Call Daxpy(Nab*Ncd,One,Work(I010),1,Work(Ioff),1)
          Endif 
          Call Putlst(Work(Ioff),1,Ncd,1,Irrep_cd,List_h)
          Call checksum("VVVV    :",Work(Ioff),Nab*Ncd)
      Enddo

C OOVV block (H(ij,ab)=-delta(ij)(F(ab)+F(ba)))

      Ioff   = Ione
      List_h = 217
      Do Irrep_ab = 1, Nirrep
         Irrep_ij = Dirprd(Irrep_ab,Irrepx)

         Nij = Irpdpd(Irrep_ij,14)
         Nab = Irpdpd(Irrep_ab,15)
         I000 = Ione
         I010 = I000 + Nab*Nij
         Iend = I010 + Nab*Nij
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nab*Nij)
         Icount = Inull 
         Do Irrep_b = 1, Nirrep
            Irrep_a = Dirprd(Irrep_b,Irrep_ab)
            Do b = 1, Vrt(Irrep_b,Ispin)
               Do a = 1, Vrt(Irrep_a,Ispin)
                  Do Irrep_j = 1, Nirrep
                     Irrep_i = Dirprd(Irrep_j,Irrep_ij)
                     Do j = 1, Pop(Irrep_j,Ispin)
                        Do i = 1, Pop(Irrep_i,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Work(Ioff+Icount) =  
     +                              -Ioo(i,j)*
     +                               (Htau_pq(an,bn)+Htau_pq(bn,an))
     +                              -Ivv(a,b)*
     +                               (Htau_pq(in,jn)+Htau_pq(jn,in))
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          If (.Not. ONEP_ONLY) Then
             Call Getlst(Work(I010),1,Nab,1,Irrep_ab,List_h)
             Call Daxpy(Nab*Nij,One,Work(I010),1,Work(Ioff),1)
          Endif 
          Call Putlst(Work(Ioff),1,Nab,1,Irrep_ab,List_h)
          Call checksum("OOVV    :",Work(Ioff),Nab*Nij)
      Enddo

C VVOO block (H(ab,ij)=-delta(ab)(F(ij)+F(ji)))

      Ioff   = Ione
      List_h = 216
      Do Irrep_ij = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)

         Nab = Irpdpd(Irrep_ab,15)
         Nij = Irpdpd(Irrep_ij,14)
         I000 = Ione
         I010 = I000 + Nab*Nij
         Iend = I010 + Nab*Nij
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_pqrs_1",
     +                                     Iend,Maxcor)
         Call Dzero(Work(I000),Nab*Nij)
         Icount = Inull 
         Do Irrep_j = 1, Nirrep
            Irrep_i = Dirprd(Irrep_j,Irrep_ij)
            Do j = 1, Pop(Irrep_j,Ispin)
               Do i = 1, Pop(Irrep_i,Ispin)
                  Do Irrep_b = 1, Nirrep
                     Irrep_a = Dirprd(Irrep_b,Irrep_ab)
                     Do b = 1, Vrt(Irrep_b,Ispin)
                        Do a = 1, Vrt(Irrep_a,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Jn = J + Ioffo(Irrep_j)
                           An = A + Ioffv(Irrep_a)
                           Bn = B + Ioffv(Irrep_b)
                           Work(Ioff+Icount) = 
     +                              -Ivv(a,b)*
     +                               (Htau_pq(in,jn)+Htau_pq(jn,in))
     +                              -Ioo(i,j)*
     +                               (Htau_pq(an,bn)+Htau_pq(bn,an))
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
          If (.Not. ONEP_ONLY) Then
             Call Getlst(Work(I010),1,Nij,1,Irrep_ij,List_h)
             Call Daxpy(Nab*Nij,One,Work(I010),1,Work(Ioff),1)
          Endif 
          Call Putlst(Work(Ioff),1,Nij,1,Irrep_ij,List_h)
          Call checksum("VVOO    :",Work(Ioff),Nab*Nij)
      Enddo

C No OVOV or VOVO like combinations 

C These blocks are not necessary for orbital optimizations 


      Call Pccd_check_hess(Work,Maxcor)

      Return
      End 


