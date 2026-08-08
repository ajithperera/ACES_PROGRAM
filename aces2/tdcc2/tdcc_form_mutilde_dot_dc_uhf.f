










      Subroutine Tdcc_form_mutilde_dot_dc_uhf(Work,Memleft,Irrepx)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Logical Rhf
      Integer Tmp_list,Hb_ph_list
      Character*4 Spcase(2)



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

      Data Spcase /'AAAA','BBBB'/

      Data One /1.0D0/

        Write(6,*)
        Write(6,"(a)") "---Entered tdcc_form_mutilde_dot_dc_uhf---"
        Write(6,*)
C The ABAB block.

      Iz2_td_list = 346
      Hb_ph_list  = 93
      Mut_s_list  = 392
      Tmp_list    = 302

      Isize_pphh = Idsymsz(Irrepx,Isytyp(1,Iz2_td_list),
     +                     Isytyp(2,Iz2_td_list))

C NDeP Eqn. 29, fourth term: Z(Ib,Aj) = Hbar(j,b) * Mut(A,I) and
C Z is built as Z(AI,bj).

      Do Irrep_bj_z = 1, Nirrep
         Irrep_ai_z = Dirprd(Irrep_bj_z,Irrepx)

         Ndim_ai_z   = Irpdpd(Irrep_ai_z,9)
         Ndim_bj_z   = Irpdpd(Irrep_bj_z,10)

         Ndim_aibj_z = Ndim_bj_z * Ndim_ai_z

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_aibj_z
         I020 = I010 + Irpdpd(1,9)
         Iend = I020 + Irpdpd(1,10)

         Ioff_hb     = I010
         Ioff_mut    = I020
         Ioff_mut_d  = I000

         Call Dzero(Work(I000),Ndim_aibj_z)

C Z(AI,bj) = Mut(b,j) Hbar(A,I)

         If (Irrep_bj_z .EQ. Irrepx) Then

            Call Getlst(Work(I010),1,1,1,1,Hb_ph_list)
            Call Getlst(Work(I020),1,1,1,2,Mut_s_list)

            Nrow = Ndim_ai_z
            NCol = Ndim_bj_z
            Nsum = 1
            Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(Ioff_hb),Nrow,
     +                  Work(Ioff_mut),1,Zero,Work(Ioff_mut_d),Nrow)
         Endif

C Z(AI,bj) = Mut(A,I) Hbar(b,j)

         If (Irrep_ai_z .EQ. Irrepx) Then

            Call Getlst(Work(I020),1,1,1,2,Hb_ph_list)
            Call Getlst(Work(I010),1,1,1,1,Mut_s_list)

            Nrow = Ndim_ai_z
            NCol = Ndim_bj_z
            Nsum = 1
            Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(Ioff_mut),Nrow,
     +                  Work(Ioff_hb),1,One,Work(Ioff_mut_d),Nrow)
         Endif

         Call Putlst(Work(I000),1,Ndim_bj_z,1,Irrep_bj_z,
     +               Tmp_list)
      Enddo

      Isize_phph = Idsymsz(Irrepx,Isytyp(1,Tmp_list),
     +                            Isytyp(2,Tmp_list))
      Isize_pphh = Idsymsz(Irrepx,Isytyp(1,Iz2_td_list),
     +                            Isytyp(2,Iz2_td_list))

      Ibgn = 1
      I000 = Ibgn
      I010 = I000  + Max(Isize_phph,Isize_pphh)
      I020 = I010  + Max(Isize_phph,Isize_pphh)
      I030 = I020  + 100
      Iend = I030

      If (Iend .GT. Memleft) Call Insmem
     +   ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

      Call Getall(Work(I000),Isize_phph,Irrepx,Tmp_list)

      Call Sstgen(Work(I000),Work(I010),Isize_phph,Vrt(1,1),Pop(1,1),
     +            Vrt(1,2),Pop(1,2),Work(I020),Irrepx,'1324')

      Call Putall(Work(I000),Isize_pphh,Irrepx,Iz2_td_list)

      Write(6,*)
      call checksum("Tdcc_form_mutilde_dot_dc,mu^t_dot(Ab,Ij)(f)  :",
     +               Work(I010),Isize_pphh,S)

      Do Ispin = 1, 2

         Tmp_list    = 299 + Ispin
         Iz2_td_list = 343 + Ispin

C Target, Z(AI,BJ): Ispin =1
C Target, Z(ai,bj): Ispin =2

         Do Irrep_bj_z = 1, Nirrep
            Irrep_ai_z = Dirprd(Irrep_bj_z,Irrepx)

            Ndim_ai_z = Irpdpd(Irrep_ai_z,8+Ispin)
            Ndim_bj_z = Irpdpd(Irrep_bj_z,8+Ispin)

            Ndim_aibj_z = Ndim_ai_z * Ndim_bj_z

            Ibgn = 1
            I000 = Ibgn
            I010 = I000 + Ndim_aibj_z 
            I020 = I010 + Irpdpd(1,8+Ispin)
            Iend = I020 + Irpdpd(Irrepx,8+Ispin)

            If (Iend .Gt. Memleft) Call Insmem
     +         ("Tdcc_form_mutilde_0_d_uhf",Iend,Memleft)

            Call Dzero(Work(I000),Ndim_aibj_z)

            Call Getlst(Work(I010),1,1,1,Ispin,Hb_ph_list)
            Call Getlst(Work(I020),1,1,1,Ispin,Mut_s_list)

            Ioff_hb     = I010
            Ioff_mut    = I020
            Ioff_mut_d  = I000

C The term Mut(B,J) Hbar(I,A) performed as 

C Z(AI,BJ) = Mut(B,J) * Hbar(A,I) ISPIN=1
C Z(ai,bj) = Mut(b,j) * Hbar(a,i) ISPIN=2

            If (Irrep_bj_z .EQ. Irrepx) Then

               Nrow = Ndim_ai_z
               Ncol = Ndim_bj_z
               Nsum = 1
               Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(Ioff_hb),
     +                     Nrow,Work(Ioff_mut),1,Zero,
     +                     Work(Ioff_mut_d),Nrow)
            Endif

C Z(AI,BJ) = Mut(A,I) * Hbar(B,J) ISPIN=1 
C Z(ai,bj) = Mut(a,i) * Hbar(b,j) ISPIN=2

            If (Irrep_ai_z .EQ. Irrepx) Then

               Nrow = Ndim_ai_z
               Ncol = Ndim_bj_z
               Nsum = 1
               Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(Ioff_hb),
     +                     Nrow,Work(Ioff_mut),1,One,
     +                     Work(Ioff_mut_d),Nrow)
            Endif

            Call Putlst(Work(I000),1,Ndim_bj_z,1,Irrep_bj_z,
     +                  Tmp_list)
         Enddo

        Isize_phph = Idsymsz(Irrepx,Isytyp(1,Tmp_list),
     +                       Isytyp(2,Tmp_list))
        Isize_pphh = Idsymsz(Irrepx,Isytyp(1,Iz2_td_list),
     +                       Isytyp(2,Iz2_td_list))
        Ibgn = 1
        I000 = Ibgn
        I010 = I000 + Max(Isize_phph,Isize_pphh)
        I020 = I010 + Max(Isize_phph,Isize_pphh)
        I030 = I020 + 100
        Iend = I030

        If (Iend .GT. Memleft) Call Insmem
     +     ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

        Call Getall(Work(I000),Isize_phph,Irrepx,Tmp_list)

        Call Dsst03i(Work(I000),Work(I010),Isize_phph,
     +               Isize_Phph,Work(I020),Spcase(Ispin),Irrepx)
    
        Call Putall(Work(I000),Isize_pphh,Irrepx,Iz2_td_list)

      Write(6,*)
      call checksum("Tdcc_form_mutilde_dot_dc,mu^t_dot(A<B,I<J)(f):",
     +               Work(I000),Isize_pphh,S)
      Enddo 

      Return
      End
