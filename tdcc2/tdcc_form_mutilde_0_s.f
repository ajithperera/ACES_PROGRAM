













































































































































































































      Subroutine Tdcc_form_mutilde_0_s(Work,Memleft,Doo,Dvv,Dvo,
     +                                 Mu_0_t_S,Irrepx,Lenoo,Lenvv,
     +                                 Lenvo,Iuhf,Mu_0)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Doo(Lenoo),Dvv(Lenvv),Dvo(Lenvo)

      Double Precision Mu_0, Mu_0_t_s(Lenvo)
      Double Precision Mone,Zero,One,Two
 
      Data Zero, One, Mone, Two /0.0D0, 1.0D0, -1.0D0, 2.0D0/



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

C Make mu_tilde(i,a) (stored as mu_tilde(a,i))

        Write(6,*)
        Write(6,"(a)") "---Entered tdcc_form_mutilde_0_s---"
        Write(6,*)
      L1_list      = 190 
      Mut_ph_list  = 392
      Ivo_0         = 1

      Do Ispin = 1, 1+Iuhf 
         
         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Irpdpd(1,8+Ispin)
         I020 = I010 + Irpdpd(1,18+Ispin) * Irpdpd(1,20+Ispin)
         Iend = I020

         If (Iend .Gt. Memleft) Call Insmem("-Tdvee_form_mu_0_S",
     +                                       Iend,Memleft)
         Call Getlst(Work(I000),1,1,1,Ispin,L1_list)

         Ndim_vo = Irpdpd(Irrepx,8+Ispin)

         Ivo_l0 = I000
         Ivo_z0 = 1 + (Ispin-1) * Irpdpd(Irrepx,9)
         Ivv_d0 = 1 + (Ispin-1) * Irpdpd(Irrepx,19)
         Ioo_d0 = 1 + (Ispin-1) * Irpdpd(Irrepx,21)

      call Checksum("Tdcc_form_mutilde_0_s,Dvo      :",
     +               Dvo(Ivo_z0),Ndim_vo,S)
         Call Dzero(Mu_0_t_s(Ivo_z0),Irpdpd(Irrepx,8+Ispin))

         Do Irrep_i = 1, Nirrep
            Irrep_a = Dirprd(Irrep_i,Irrepx)

            Ndim_i = Pop(Irrep_i,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
               
            Ivo_z = Ivo_z0 + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)
            Ivo_d = Ivo_z

C NDeP Eqn. 20, First term: z(a,i) = Mu(a,i) (and Z(A,I)=Mu(A,I))

            Ndim_ai = Ndim_i * Ndim_a
            Call Dcopy(Ndim_ai,Dvo(Ivo_d),1,Mu_0_t_s(Ivo_z),1)

C NDep Eqn. 20 Second term: z(A,I) = Mu(E,A) * L(I,E) (and 
C z(a,i) = Mu(e,a) * L(i,e)). The L is stored as L(E,I).
            
            Irrep_e = Irrep_i

            Ndim_e = Vrt(Irrep_e,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
            Ndim_i = Pop(Irrep_i,Ispin)
               
            Nrow = Ndim_a
            Nsum = Ndim_e
            Ncol = Ndim_i

            Ivv_d = Ivv_d0 + (Isymoff(Irrep_a,Irrepx,18+Ispin)-1)
            Ivo_l = Ivo_l0 + (Isymoff(Irrep_i,1,8+Ispin)-1)

            Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Dvv(Ivv_d),Nrow,
     +                  Work(Ivo_l),Nsum,One,Mu_0_t_s(Ivo_z),Nrow)
         Enddo


C NDeP Eqn. 20, Two third terms: 
C a. z(a,i) = -L(m,a)*Mu(i,m)  (and z(A,I) = -L(M,A)*Mu(I,M))

         Ivo_l0 = I000
         Do Irrep_i = 1, Nirrep
            Irrep_a = Dirprd(Irrep_i,Irrepx)
            Irrep_m = Irrep_a

            Ndim_m = Pop(Irrep_m,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
            Ndim_i = Pop(Irrep_i,Ispin)
            
            Ioo_d = Ioo_d0 + (Isymoff(Irrep_m,Irrepx,20+Ispin)-1)
            Ivo_l = Ivo_l0 + (Isymoff(Irrep_m,1,8+Ispin)-1)
            Ivo_z = Ivo_z0 + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)

            Nrow = Ndim_a
            Nsum = Ndim_m
            Ncol = Ndim_i

            Call Xgemm("N","T",Nrow,Ncol,Nsum,Mone,Work(Ivo_l),
     +                 Nrow,Doo(Ioo_d),Ncol,One,Mu_0_t_s(Ivo_z),
     +                 Nrow)
C           call checksum("check1",Mu_0_t_s(Ivo_z),Nrow*Ncol,s)
         Enddo 


C z(a,i) = L(m,a)*Mu(i,m)*delta(i,m) (and z(A,I) = L(M,A) * 
C                                         Mu(I,M)Delta(I,M)
         If (Irrepx .EQ. 1) Then
             Call Daxpy(Lenvo,Mu_0,Work(Ivo_l0),1,Mu_0_t_s(Ivo_z0),1)
         Endif 

      Enddo

C NDep Eqn. 20, The third term Z(A,I) = L(MI,EA)*mu(E,M) 
C UHF  Z(A,I) = L(MI,EA) * mu(E,M) and Z(a,i) = L(mi,ea) * mu(e,m)
c Lists 137,136: L(AI,em) and L(ai,EM)
C RHF  Z(A,I) = {2*L(AI,em) -L(Am,eI)} * mu(e,m)

      Do Ispin = 1, 1+Iuhf 

C UHF Ispin=1, L(AI,em)*Mu(e,m); Ispin=2, L(ai,EM)*Mu(E,M)
C     Ispin=2, L(AI,EM)*Mu(E,M);          L(ai,em)*Mu(e,m)
C
C RHF Ispin=1  {2*L(AI,em)-L(Am,eI)}Mu(E,M)

         Ivo_z  = 1 + (Ispin-1) * Irpdpd(Irrepx,8+Ispin)

         L2_list_ab0 = 138 - Ispin
         L2_list_ab1 = 139
         L2_list_pq0 = 133 + Ispin

C Irrep_em=Irrep_ai=Irrepx

         Irrep_em = Irrepx
         Irrep_ai = Irrepx
         Irrep_ei = Irrepx

         Len_em_pq = Irpdpd(Irrepx,11-Ispin)
         Len_ai_pq = Irpdpd(Irrepx,8+Ispin)
         Len_em_pp = Irpdpd(Irrepx,8+Ispin)
         Len_ai_pp = Irpdpd(Irrepx,8+Ispin)
         Len_am_ab = Irpdpd(Irrepx,11)
         Len_ei_ab = Irpdpd(Irrepx,12)

         Max_pq = Max(Len_em_pq,Len_ai_pq,Len_ai_pp,Len_em_pp,
     +                Len_am_ab,Len_ei_ab)
         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Max_pq*Max_pq 
         I020 = I010 + Max_pq*Max_pq
         Iend = I020 

         If (Iend .GT. Memleft) Call Insmem("-Tdvee_form_mutilde_0_s",
     +                                       Iend,Memleft)

C UHF Ispin=1, Z(A,I) = L(AI,em)*Mu(e,m)
C     Ispin=2  Z(a,i) = L(ai,EM)*Mu(E,M)
C RHF Ispin=1  Z(A,I) = {2*L(AI,em) -L(Am,eI)} Mu(e,m)

         Call Getlst(Work(I000),1,Len_em_pq,1,Irrep_em,L2_list_ab0)

         If (Iuhf .EQ. 0) Then
             Len_aiem_ab = Len_am_ab * Len_ei_ab
             Call Getlst(Work(I010),1,Len_ei_ab,1,Irrep_ei,
     +                    L2_list_ab1)
             Call Dscal(Len_aiem_ab,Two,Work(I000),1)
             Call Daxpy(Len_aiem_ab,Mone,Work(I010),1,Work(I000),1)
         Endif

         If (Ispin .EQ. 1) Then
            Ivo_d0 = 1
            Ivo_d1 = 1 + Iuhf * Irpdpd(Irrepx,9)
            Ivo_z0 = 1
            Ivo_z1 = 1 
            Nrow0  = Len_ai_pp
            Ncol0  = 1
            Nsum0  = Len_em_pp
            If (Iuhf .Eq. 0) Then 
                Nrow1 = Len_am_ab
                Ncol1 = 1
                Nsum1 = Len_ei_ab 
            Else
                Nrow1 = Len_ai_pq
                Ncol1 = 1
                Nsum1 = Len_em_pq
            Endif 
              
         Elseif (Ispin .EQ. 2) Then
            Ivo_d0 = 1 + Iuhf * Irpdpd(Irrepx,9)
            Ivo_d1 = 1 
            Ivo_z0 = 1 + Iuhf * Irpdpd(Irrepx,9)
            Ivo_z1 = 1 + Iuhf * Irpdpd(Irrepx,9)
            Nrow0  = Len_ai_pp
            Ncol0  = 1
            Nsum0  = Len_em_pp
            Nrow1  = Len_ai_pq
            Ncol1  = 1
            Nsum1  = Len_em_pq
         Endif 

            Call Xgemm("N","N",Nrow1,Ncol1,Nsum1,One,Work(I000),Nrow1,
     +                 Dvo(Ivo_d1),Nrow1,One,Mu_0_t_s(Ivo_z1),Nrow1)

            If (Iuhf .NE. 0) Then

C UHF Ispin=1, L(AI,EM)*Mu(E,M)
C     Ispin=2  L(ai,em)*Mu(e,m)

                Call Getlst(Work(I000),1,Len_em_pp,1,Irrep_em,
     +                      L2_list_pq0)
                Call Xgemm("N","N",Nrow0,Ncol0,Nsum0,Mone,Work(I000),
     +                     Nrow0,Dvo(Ivo_d0),Nrow0,One,
     +                     Mu_0_t_s(Ivo_z0),Nrow0)
            Endif 

         Call Putlst(Mu_0_t_s(Ivo_z0),1,1,1,Ispin,Mut_ph_list)

      Write(6,*)
      Ndim_vo = Irpdpd(Irrepx,8+Ispin)
      call Checksum("Tdcc_form_mutilde_0_s,Mu^t(f)  :",
     +               Mu_0_t_s(Ivo_z0),Ndim_vo,S)

      Enddo 
  
      Return
      End
