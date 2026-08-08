










      Subroutine Tdcc_form_mutilde_dot_sc(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Logical Rhf
C
      Dimension Work(Memleft)
      Integer Hp_sym,Hh_sym
      Integer Hbar_Ppph_list,A
C
      Data One,Onem,Zero/1.0D0,-1.0D0,0.0D0/



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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
c symloc.com : begin
c The numbering scheme is as follows:
c   a<b    (alpha) [1]
c   a<b    (beta)  [2]
c   i<j    (alpha) [3]
c   i<j    (beta)  [4]
c   a<=b   (alpha) [5]
c   a<=b   (beta)  [6]
c   i<=j   (alpha) [7]
c   i<=j   (beta)  [8]
c   a,i    (alpha) [9]
c   a,i    (beta)  [10]
c   a,i    (AB)    [11]
c   a,i    (BA)    [12]
c   a,b    (AB)    [13]
c   i,j    (AB)    [14]
c   a,b    (AB)    [15]
c   i,a    (alpha) [16]
c   i,a    (beta)  [17]
c   i,a    (AB)    [18]
c   a,b    (alpha) [19]
c   a,b    (beta)  [20]
c   i,j    (alpha) [21]
c   i,j    (beta)  [22]
c   a,b    (BA)    [23]
c   i,j    (BA)    [24]
c   i,a    (BA)    [25]
      integer         isymoff(8,8,25)
      common /symloc/ isymoff
c symloc.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

        Write(6,*)
        Write(6,"(a)") "---Entered  tdcc_form_mutilde_dot_sc---"
        Write(6,*)
      Rhf = .False.
      If (Iuhf .Eq. 0) Rhf = .True.
      Mut_d_list     = 396

C Z(A,I) = Hbar(FE,AM) * Mut(IM,FE) +  Hbar(Fe,Am) * Mut(Im,Fe)
C Z(a,i) = Hbar(fe,am) * Mut(im,fe) +  Hbar(fE,aM) * Mut(iM,fE)
  
      Do Ispin = 1, 1+Iuhf

         Ndim_ai = Irpdpd(Irrepx,8+Ispin)

C Hbar Read as Hbar(Fe,Am)  Ispin=1
C              Hbar(Ef,Ma)  Ispin=2
C Mu_tilde read as Mut(Fe,Im) 

         Ibgn = 1
         I000 = 1
         I010 = I000 + Ndim_ai
         Call Dzero(Work(I000),Ndim_ai)
         
         Do Irrep_im = 1, Nirrep 

            Hbar_ppph_list = 131-Ispin
            Mut_list       = 326

            Irrep_ef = Dirprd(Irrep_im,Irrepx)
            Irrep_am = Irrep_ef 

            Ndim_ef_z = Irpdpd(Irrep_ef,Isytyp(1,Mut_list))
            Ndim_im_z = Irpdpd(Irrep_im,Isytyp(2,Mut_list))
            Ndim_ef_h = Irpdpd(Irrep_ef,Isytyp(1,Hbar_Ppph_list))
            Ndim_am_h = Irpdpd(Irrep_am,Isytyp(2,Hbar_Ppph_list))
            Max_z     = Max(Ndim_ef_z,Ndim_im_z)
            Max_h     = Max(Ndim_ef_h,Ndim_am_h)
            Max_zh    = Max(Max_h,Max_z)
            
           Ibgn = I010
           I020 = I010 + Ndim_ef_z * Ndim_im_z
           Iend = I020 
           Call Getlst(Work(I010),1,Ndim_im_z,1,Irrep_im,Mut_list)

C Mut(Fe,Im) Ispin=1 or 2
            
           If (Ispin .EQ. 1) Then
       
C Mut(Fe,Im) -> Mut(Fe,mI)

              Ibgn  = I020
              Iptr1 = Ibgn
              Iptr2 = Iptr1 + Max_zh
              Iptr3 = Iptr2 + Max_zh
              Iptr4 = Iptr3 + Max_zh
              Iend  = Iptr4

              If (Iend .Gt. Memleft) Call Insmem(
     +                     "@-Tdcc_form_mutilde_dot_sa",
     +                      Memleft,Iend)

              Call Symtr1(Irrep_im,Pop(1,1),Pop(1,2),Ndim_ef_z,
     +                    Work(I010),Work(IPtr1),Work(Iptr2),
     +                    Work(IPtr3))
              If (RHF) Then
                 Call Spinad3(Irrep_ef,Vrt(1,1),Ndim_ef_z,Ndim_im_z,
     +                        Work(I010),Work(IPtr1),Work(Iptr2))
              Endif 
           Endif 

C Hbar(Fe,Am) Ispin=1
C Hbar(Ef,Ma) Ispin=2

           Ibgn = I020
           I030 = I020 + Ndim_ef_h * Ndim_am_h
           I040 = I030 + 3*Max_zh
           Iend = I040
        
           If (Iend .LT. Memleft) Then

               Call Getlst(Work(I020),1,Ndim_am_h,1,Irrep_am,
     +                     Hbar_ppph_list)

C Hbar(Fe,Am) -> Hbar(Fe,mA)  Ispin=1

               If (Ispin .EQ. 1) Then
                  Ibgn = I030 
                  Iptr1 = Ibgn
                  Iptr2 = Iptr1 + Max_zh
                  Iptr3 = Iptr2 + Max_zh
                  Iptr4 = Iptr3 + Max_zh
                  Iend  = Iptr4

                  If (Iend .Gt. Memleft) Call Insmem(
     +                     "@-Tdcc_form_mutilde_dot_sa",
     +                      Memleft,Iend)

                   Call Symtr1(Irrep_am,Vrt(1,1),Pop(1,2),
     +                         Ndim_ef_h,Work(I020),Work(Iptr1),
     +                         Work(Iptr2),Work(Iptr3))
             
               Endif

C Z(A,I) =  Hbar(Fe,mA)^t * Mut(Fe,mI) Ispin=1
C Z(a,i) =  Hbar(Ef,Ma)^t * Mut(Fe,im) Ispin=2
 
              Ioff_mut_d = I000
              If (Ispin .EQ.1) Then
                 Hp_sym=25
                 Hh_sym=24
              Else
                 Hp_sym=18
                 Hh_sym=14
              Endif 
 
              Do Irrep_i = 1, Nirrep 
                 Irrep_m = Dirprd(Irrep_i,Irrep_im)
                 Irrep_a = Dirprd(Irrep_m,Irrep_am)

                 Nrow = Vrt(Irrep_a,Ispin)
                 Ncol = Pop(Irrep_i,Ispin)
                 Nsum = Ndim_ef_h * Pop(Irrep_m,3-Ispin)
                 
                 Ioff_hb    = I020 + Ndim_ef_h*(Isymoff(Irrep_a,
     +                                   Irrep_am,hP_sym)-1)
                 Ioff_mut   = I010 + Ndim_ef_z*(Isymoff(Irrep_i,
     +                                   Irrep_im,Hh_sym)-1)

                 Call Xgemm("T","N",Nrow,Ncol,Nsum,One,
     +                       Work(Ioff_hb),Nsum,
     +                       Work(Ioff_mut),Nsum,One,
     +                       Work(Ioff_mut_d),Nrow)

                 Ioff_mut_d = Ioff_mut_d + Nrow * Ncol
              Enddo 

           Else
c 
              If (Ispin .EQ. 1) Then

C Hbar(Fe,Am) * Mut(Fe,mI)

                 Hh_sym      = 24
                 Ihbar_start = 1
                 Do irrep_m = 1, Nirrep
                    Irrep_a = Dirprd(Irrep_m,Irrep_am)
                    Irrep_i = Dirprd(Irrep_a,Irrepx)
                    
                     Ndim_m = Pop(irrep_m,2)
                     Ndim_a = Vrt(Irrep_a,1)
                     Ndim_i = Pop(Irrep_i,1)

                     Ibgn = I020 
                     I030 = I020 + Ndim_ef_h * Ndim_a 
                     Iend = I030 
                     If (Iend .Gt. Memleft) Call Insmem(
     +                        "@-Tdcc_form_mutilde_dot_sa",
     +                         Memleft,Iend)

                     Ioff_mut = I010 + Ndim_ef_z*(Isymoff(Irrep_i,
     +                                 Irrep_im,Hh_sym)-1)
                     Do m = 1, Ndim_m

                         Call Getlst(Work(I020),Ihbar_start,
     +                               Ndim_a,1,Irrep_am,
     +                               Hbar_Ppph_list)
                         Nrow = Ndim_a
                         Ncol = Ndim_i
                         Nsum = Ndim_ef_h

                         Ioff_mut_d = I000 + (Isymoff(Irrep_i,
     +                                        Irrepx,9)-1)
                         Ioff_hb    = I020
                         Call Xgemm("T","N",Nrow,Ncol,Nsum,One,
     +                               Work(Ioff_hb),Nsum,
     +                               Work(Ioff_mut),Nsum*Ndim_m,One,
     +                               Work(Ioff_mut_d),Nrow)

                         Ioff_mut    = Ioff_mut + Ndim_ef_z
                         Ihbar_start = Ihbar_start + Ndim_a
                     Enddo 
                  Enddo

              Else 

C Hbar(Ef,Ma) * Mut(fE,Mi)

                 Ihbar_start = 1
                 Hh_sym      = 14
                 Do irrep_a = 1, Nirrep
                    Irrep_m = Dirprd(Irrep_a,Irrep_am)
                    Irrep_i = Dirprd(Irrep_a,Irrepx)

                     Ndim_m = Pop(irrep_m,1)
                     Ndim_a = Vrt(Irrep_a,2)
                     Ndim_i = Pop(Irrep_i,2)

                     Ibgn = I020
                     I030 = I020 + Ndim_ef_h * Ndim_m
                     Iend = I030
                     If (Iend .Gt. Memleft) Call Insmem(
     +                        "@-Tdcc_form_mutilde_dot_sa",
     +                         Memleft,Iend)

                     Ioff_mut = I010 + Ndim_ef_z*(Isymoff(Irrep_i,
     +                                 Irrep_im,Hh_sym)-1)
                     Do A = 1, Ndim_a

                         Call Getlst(Work(I020),Ihbar_start,
     +                               Ndim_m,1,Irrep_am,
     +                               Hbar_Ppph_list)
                         Nrow = 1
                         Ncol = Ndim_i
                         Nsum = Ndim_ef_h * Ndim_m

                         Ioff_mut_d = I000 + (Isymoff(Irrep_i,
     +                                        Irrepx,10)-1)+(A-1)
                         Ioff_hb    = I020
                         Call Xgemm("T","N",NroW,Ncol,Nsum,One,
     +                               Work(Ioff_hb),Nsum,
     +                               Work(Ioff_mut),Nsum,One,
 
     +                               Work(Ioff_mut_d),Nrow)

                         Ihbar_start = Ihbar_start + Ndim_m
                     Enddo
                 Enddo

              Endif 

           Endif 

C AAAA and BBBB contributions 
C Z(A,I) = Hbar(FE,AM) * Mut(IM,FE) 
C Z(a,i) = Hbar(fe,am) * Mut(im,fe) 

C Read Hbar(F<E,AM) and Hbar(f<e,am) and Mut(F<E,I<M) and Mut(f<e,i<m)
           If (.NOT. RHF) Then

              Hbar_Ppph_list = 126 + Ispin
              Mut_list       = 323 + Ispin

              Ndim_fe_c_h = Irpdpd(Irrep_ef,
     +                           Isytyp(1,Hbar_ppph_list))
              Ndim_am_f_h = Irpdpd(Irrep_am,
     +                           Isytyp(2,Hbar_ppph_list))

              Ndim_fe_c_z = Irpdpd(Irrep_ef,Isytyp(1,Mut_list))
              Ndim_im_c_z = Irpdpd(Irrep_im,Isytyp(2,Mut_list))

              Ndim_im_f   = Irpdpd(Irrep_im,20+Ispin)
              Ndim_fe_f   = Irpdpd(Irrep_ef,18+Ispin)

              Max_buf = Max(Ndim_fe_c_h,Ndim_im_c_z,Ndim_fe_f)

              Ibgn = I010
              I020 = I010 + Ndim_fe_c_z * Ndim_im_f
              I030 = I020 + Ndim_fe_c_h * Ndim_am_f_h
              I040 = I030 + 3*Max_buf 
              Iend = I040 

C Read Mut(F<E,I<M) Ispin=1 and Mut(f<e,i<m) Ispin=2
              
              Call Getlst(Work(I010),1,Ndim_im_c_z,1,Irrep_im,
     +                       Mut_list)

C Expand ket to Mut(F<E,IM) Ispin=1 and Mut(f<e,im) Ispin=2

              Call Symexp(Irrep_im,Pop(1,Ispin),Ndim_fe_c_z,
     +                       Work(I010))

              If (Iend .LT. Memleft) Then

C Hbar(F<E,AM) ispin=1
C Hbar(f<e,am) ispin=2
               
                 Call Getlst(Work(I020),1,Ndim_am_f_h,1,Irrep_am,
     +                       Hbar_ppph_list)

C Hbar(F<E,AM) -> Hbar(F<E,MA) ispin=1
C Hbar(f<e,am) -> Hbar(f<e,ma) ispin=2

                 Ibgn  = I030 
                 Iptr1 = Ibgn 
                 Iptr2 = Iptr1 + Max_buf
                 IPtr3 = Iptr2 + Max_buf
                 Iptr4 = Iptr3 + Max_buf
                 Iend  = Iptr4

                 If (Iend .Gt. Memleft) Call Insmem(
     +                    "@-Tdcc_form_mutilde_dot_sc",
     +                     Memleft,Iend)

                 Call Symtr1(Irrep_am,Vrt(1,Ispin),Pop(1,Ispin),
     +                      Ndim_fe_c_h,Work(I020),Work(Iptr1),
     +                      Work(Iptr2),Work(Iptr3))

C Z(A,I) = Mut(F<E,IM) * Hbar(F<E,MA)  Ispin=1
C Z(a,i) = Mut(f<E,im) * Hbar(f<e,ma)  Ispin=2 

                 Ioff_mut_d = I000
                 Do Irrep_i = 1, Nirrep
                    Irrep_a = Dirprd(Irrep_i,Irrepx)
                    Irrep_m = Dirprd(Irrep_i,Irrep_im)

                    Nrow = Vrt(Irrep_a,Ispin)
                    Ncol = Pop(Irrep_i,Ispin)
                    Nsum = Ndim_fe_c_h * Pop(Irrep_m,Ispin)
          
                    Ioff_hb    = I020 + Ndim_fe_c_h *
     +                                  (Isymoff(Irrep_a,
     +                                  Irrep_am,15+Ispin)-1)

                    Ioff_mut   = I010 + Ndim_fe_c_z *
     +                                  (Isymoff(Irrep_i,
     +                                  Irrep_im,20+Ispin)-1)

                    Call Xgemm("T","N",Nrow,Ncol,Nsum,Onem,
     +                         Work(Ioff_hb),
     +                         Nsum,Work(Ioff_mut),Nsum,One,
     +                         Work(Ioff_mut_d),Nrow)

                    Ioff_mut_d = Ioff_mut_d + Nrow * Ncol
                 Enddo

              Else

C Z(A,I) = Mut(F<E,IM) * Hbar(F<E,AM)  Ispin=1
C Z(a,i) = Mut(f<E,im) * Hbar(f<e,am)  Ispin=2

                 Ihbar_start= 1
                 Do Irrep_m = 1, Nirrep
                    Irrep_a = Dirprd(Irrep_m,Irrep_am)
                    Irrep_i = Dirprd(Irrep_a,Irrepx)

                    Ndim_m = Pop(Irrep_m,Ispin)
                    Ndim_a = Vrt(Irrep_a,Ispin)
                    Ndim_i = Pop(Irrep_i,Ispin)

                    Ibgn = I020
                    I030 = I020 + Ndim_fe_c_h * Ndim_a
                    Iend = I030 
                    If (Iend .Gt. Memleft) Call Insmem(
     +                    "@-Tdcc_form_mutilde_dot_sa",
     +                     Memleft,Iend)

                    Ioff_mut = I010 + Ndim_fe_c_z * 
     +                                (Isymoff(Irrep_i,
     +                                Irrep_im,20+Ispin)-1)

                    Do M = 1, Ndim_m
                       Call Getlst(Work(I020),Ihbar_start,
     +                               Ndim_a,1,Irrep_am,
     +                               Hbar_Ppph_list)
                       Nrow = Ndim_a
                       Ncol = Ndim_i
                       Nsum = Ndim_fe_c_h

                       Ioff_hb    = I020 
                       Ioff_mut_d = I000 + (Isymoff(Irrep_i,
     +                                      Irrepx,8+Ispin)-1)

                       Call Xgemm("T","N",Nrow,Ncol,Nsum,Onem,
     +                            Work(Ioff_hb),
     +                            Nsum,Work(Ioff_mut),Nsum*Ndim_m,
     +                            One,Work(Ioff_mut_d),Nrow)
                 
                       Ihbar_start = Ihbar_start + Ndim_a
                       Ioff_mut    = Ioff_mut    + Ndim_fe_c_z
                   Enddo
                Enddo

              Endif 

           Endif 

         Enddo

      call Checksum("Tdcc_form_mu_dot_sc,Mutilde_dot(f):",
     +               Work(I000),Ndim_ai,S)
         Call Getlst(Work(I010),1,1,1,Ispin,Mut_d_list)
         Call Daxpy(Ndim_ai,One,Work(I010),1,Work(I000),1)
         Call Putlst(Work(I000),1,1,1,Ispin,Mut_d_list)

      Enddo 

      Return
      End 

