













































































































































































































      Subroutine Tdcc_form_mutilde_0_0(Work,Memleft,Doo,Dvv,Dvo,
     +                                 Irrepx,Lenoo,Lenvv,Lenvo,
     +                                 Iuhf,Mu_t_0)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Doo(Lenoo),Dvv(Lenvv),Dvo(Lenvo)
      Double Precision Mu_ii_0,Mu_t_l,Mu_t_0



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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

        Write(6,*)
        Write(6,"(a)") "---Entered tdcc_form_mutilde_0_0---"
        Write(6,*)
      I000 = 1
      Iend = I000 + Lenvo 
      If (Iend .Gt. Memleft) Call Insmem("@-Tdvee_form_mu_0_0",
     +                                    Memleft,Iend)

C Both contributions contribute only when irrepx=1.
C NDeP Eqn. 19, first-term, M(I,J)*delta(I,J) + M(i,j)*delta(i,j)

      L1_list = 190 
      Ioff_vo = 1

      Do Ispin = 1, 1+Iuhf 

         Mu_ii_0 = 0.0D0
         Mu_t_l  = 0.0D0
         Mu_t_0  = 0.0D0

         Istart = (Ispin-1)*Irpdpd(Irrepx,21)

         Do Irrep_j = 1, Nirrep
            Irrep_i = Dirprd(Irrep_j,Irrepx)
            If (Irrep_i .EQ. Irrep_j) Then
               Do Jocc = 1, Pop(Irrep_j,Ispin)
                  Ioff = Istart + Jocc +  
     +                   (Jocc-1)* Pop(Irrep_i,Ispin)
                   Mu_ii_0 = Mu_ii_0 + Doo(Ioff) 
               Enddo
            Endif
            Istart = Ioff
         Enddo

C NDeP Eqn. 19, second-term, Mu(A,I)*L(A,I) + Mu(a,i)*L(a,i)

         If (Irrepx .EQ. 1) Then
            Ioff_vo = Ioff_vo + (Ispin-1) * Irpdpd(Irrepx,9)
            Ndim_vo =  Irpdpd(Irrepx,8+ispin)

            Call Getlst(Work(I000),1,1,1,Ispin,L1_list)
            Mu_t_l =  Ddot(Ndim_vo,Work(I000),1,Dvo(ioff_vo),1)
            Mu_t_0 = Mu_ii_0 + Mu_t_l

         Endif 

      If (Ispin .Eq.1) Write(6,"(a,a,3(1x,F12.8))") "The mu_ii_0,",
     +                 "mu_t_l,mu_t_0, alpha", Mu_ii_0,Mu_t_l,
     +                  Mu_t_0
      If (Ispin .Eq.2) Write(6,"(a,a,3(1x,F12.8))") "The mu_ii_0,",
     +                 "mu_t_l,mu_t_0, beta ", Mu_ii_0,Mu_t_l,
     +                  Mu_t_0
      Enddo 

      Return
      End
