













































































































































































































      Subroutine Tdcc_form_mu_dot_sa(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)

      Double Precision Mone,Mtwo 
      Integer Hbr_pp_list,Hbr_hh_list,Hbar_phph_list
      Integer Hbar_phph_list_a,Hbar_phph_list_b
 
      Data Zero, One, Mone, Two, Mtwo /0.0D0,  1.0D0, -1.0D0, 
     +                                2.0D0, -2.0D0/


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
        Write(6,"(a)") "---Entered  tdcc_form_mu_dot_sa---"
        Write(6,*)
      Mu0_ph_list = 390
      Hbr_pp_list = 92
      Hbr_hh_list = 91
      Mud_ph_list = 394

      Do Ispin = 1, 1+Iuhf 
         
         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Irpdpd(Irrepx,8+Ispin)
         I020 = I010 + Irpdpd(1,18+Ispin) 
         I030 = I020 + Irpdpd(1,20+Ispin) 
         I040 = I030 + Irpdpd(Irrepx,8+Ispin)
         Iend = I040

         If (Iend .Gt. Memleft) Call Insmem("@-Tdcc_form_mu_0_S",
     +                                       Memleft,Iend)

         Call Getlst(Work(I000),1,1,1,Ispin,Mu0_ph_list)
         Call Getlst(Work(I010),1,1,1,Ispin,Hbr_pp_list)
         Call Getlst(Work(I020),1,1,1,Ispin,Hbr_hh_list)


         Ioff_m0_bgn = I000 
         Ioff_pp_bgn = I010 
         Ioff_hh_bgn = I020 
         Ioff_zd_bgn = I030 

         Ndim_zd = Irpdpd(Irrepx,8+Ispin)
         Call Dzero(Work(I030),Ndim_zd)

C NDeP Eqn. 25, First term: Mu(e,i) * Hbar(a,e) = Zdot(a,i)
C Stored as Hbar(e,a)

         Do Irrep_i = 1, Nirrep
            Irrep_e = Dirprd(Irrep_i,Irrepx)
            Irrep_a = Irrep_e

            Ndim_i = Pop(Irrep_i,Ispin)
            Ndim_e = Vrt(Irrep_e,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
               
            Ioff_pp = Ioff_pp_bgn + (Isymoff(Irrep_e,1,18+Ispin)-1)
            Ioff_m0 = Ioff_m0_bgn + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)
            Ioff_zd = Ioff_zd_bgn + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)

            Nrow = Ndim_a
            Nsum = Ndim_e
            Ncol = Ndim_i

            Call Xgemm("T","N",Nrow,Ncol,Nsum,One,Work(Ioff_pp),
     +                  Nrow,Work(Ioff_m0),Nsum,One,Work(Ioff_zd),
     +                  Nrow)
         Enddo

C NDeP Eqn. 25, Second term: Mu(a,m) * Hbar(m,i) = Zdot(a,i)

         Do Irrep_m = 1, Nirrep
            Irrep_a = Dirprd(Irrep_m,Irrepx)
            Irrep_i = Irrep_m

            Ndim_m = Pop(Irrep_m,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
            Ndim_i = Pop(Irrep_i,Ispin)

            Ioff_hh = Ioff_hh_bgn + (Isymoff(Irrep_i,1,20+Ispin)-1)
            Ioff_m0 = Ioff_m0_bgn + (Isymoff(Irrep_m,Irrepx,8+Ispin)-1)
            Ioff_zd = Ioff_zd_bgn + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)
            
            Nrow = Ndim_a
            Nsum = Ndim_m
            Ncol = Ndim_i

            Call Xgemm("N","N",Nrow,Ncol,Nsum,Mone,Work(Ioff_m0),
     +                 Nrow,Work(Ioff_hh),Nsum,One,Work(Ioff_zd),
     +                 Nrow)
         Enddo 

      call Checksum("Tdcc_form_mu_dot_sa,Mu_dot(f1) :",
     +               Work(Ioff_zd_bgn),Ndim_zd,S)
         Call Putlst(Work(Ioff_zd_bgn),1,1,1,Ispin,Mud_ph_list)

      Enddo 

C NDep Eqn. 25, Third term:

C UHF:  Zdot(A,I) = Mu(E,M)*Hbar(MA,EI) + Mu(e,m)Hbar(mA,eI)
C       Zdot(a,i) = Mu(E,m)*Hbar(ma,Ei) + Mu(E,M)Hbar(Ma,Ei)
C RHF:  Zdot(A,I) = Mu(E,M){(2 Hbar(Ma,Ei) - Hba(Ma,Ie)}

      Do Ispin = 1, 1+Iuhf

         If (IUHF .EQ. 0) Then
            Hbar_phph_list = 56
         Else
            Hbar_phph_list = 58 - Ispin
         Endif 
         
         Ndim_hb_em = Irpdpd(Irrepx,Isytyp(1,Hbar_phph_list))
         Ndim_hb_ai = Irpdpd(Irrepx,Isytyp(2,Hbar_phph_list))
         Ndim_mu_em = Irpdpd(Irrepx,11-Ispin)
         Ndim_zs_ai = Irpdpd(Irrepx,8+Ispin)
         Ndim_zd_ai = Ndim_zs_ai 

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Max(Ndim_zd_ai,Ndim_mu_em)
         I020 = I010 + Max(Ndim_zs_ai,Ndim_mu_em)
         I030 = I020 + Ndim_hb_em * Ndim_hb_ai
         Iend = I030
         If (Iend .GT. Memleft) Call Insmem("@-Tdcc_form_mu_0_S",
     +                                       Memleft,Iend)
         If (Iuhf .NE. 0) Then

C Form Z(A,I) = Hbar(em,AI)*mu(e,m) (Ispin=1)
C Form Z(a,i) = Hbar(EM,ai)*mu(E,M) (Ispin=2)

            Call Getlst(Work(I000),1,1,1,Ispin,Mud_ph_list)
            Call Getlst(Work(I010),1,1,1,3-Ispin,Mu0_ph_list)
            Call Getlst(Work(I020),1,Ndim_hb_ai,1,Irrepx,
     +                       Hbar_phph_list)
            Nsum = Ndim_mu_em
            Nrow = 1
            Ncol = Ndim_hb_ai

         Else 

C List 54=Hbar(EM,AI), List 58 = (Em,Ai)

            Call Getlst(Work(I000),1,1,1,1,Mud_ph_list)
            Call Getlst(Work(I010),1,1,1,1,Mu0_ph_list)
            Hbar_phph_list_a = 54
            Hbar_phph_list_b = 58

            Ibgn = I030
            I040 = I030 + Ndim_hb_em * Ndim_hb_ai
            Iend = I040
            If (Iend .Gt. Memleft) Call Insmem("@-Tdcc_form_mu_0_S",
     +                                          Memleft,Iend)
            
            Call Getlst(Work(I020),1,Ndim_hb_ai,1,Irrepx,
     +                  Hbar_phph_list_b)
            Call Getlst(Work(I030),1,Ndim_hb_ai,1,Irrepx,
     +                  Hbar_phph_list_a)
            Call Daxpy(Ndim_hb_em * Ndim_hb_ai,Mtwo,Work(I030),1,
     +                 Work(I020),1)
            Nsum = Ndim_mu_em
            Nrow = 1
            Ncol = Ndim_hb_ai

         Endif

         Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(I010),
     +               Nrow,Work(I020),Nsum,One,Work(I000),
     +               Nrow)

         If (Iuhf .NE. 0) Then

C Form Z(A,I) = Hbar(EM,AI)*mu(E,M) (Ispin=1)
C Form Z(a,i) = Hbar(em,ai)*mu(e,m) (Ispin=2)

            Hbar_phph_list = 53 + Ispin

            Ndim_hb_em = Irpdpd(Irrepx,Isytyp(1,Hbar_phph_list))
            Ndim_hb_ai = Irpdpd(Irrepx,Isytyp(2,Hbar_phph_list))
            Ndim_zs_ai = Irpdpd(Irrepx,8+Ispin)
            Ndim_mu_em = Irpdpd(Irrepx,11-Ispin)

            Ibgn = 1
            I000 = Ibgn
            I010 = I000 + Max(Ndim_zs_ai,Ndim_zs_em)
            I020 = I010 + Max(Ndim_mu_ai,Ndim_mu_em)
            I030 = I020 + Ndim_hb_em * Ndim_hb_ai 
            Iend = I030
            If (Iend .Gt. Memleft) Call Insmem("@-Tdcc_form_mu_0_S",
     +                                         Memleft,Iend)
       
            Call Getlst(Work(I010),1,1,1,Ispin,Mu0_ph_list)
            Call Getlst(Work(I020),1,Ndim_hb_ai,1,Irrepx,Hbar_phph_list)

            Nsum = Ndim_hb_em
            Nrow = 1
            Ncol = Ndim_hb_ai

C Double check the minus sign (Stored as Hbar(ME,AI))
         
            Call Xgemm("N","N",Nrow,Ncol,Nsum,Mone,Work(I010),
     +                  Nrow,Work(I020),Nsum,One,Work(I000),
     +                  Nrow)
         Endif

      call Checksum("Tdcc_form_mu_dot_sa,Mu_dot(f2) :",
     +               Work(I000),Ndim_zd_ai,S)


      Enddo 

      Return
      End
     
           
        
