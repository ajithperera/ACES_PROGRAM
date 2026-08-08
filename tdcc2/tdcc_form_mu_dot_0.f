













































































































































































































      Subroutine Tdcc_form_mu_dot_0(Work,Memleft,Irrepx,Iuhf,Mu0_d)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Double Precision Mu0_d
      Integer Hbr_hp_list,Mu0_hp_list



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
       
      Hbr_hp_list = 93
      Mu0_hp_list = 390
      
      Max_vo = Max(Irpdpd(Irrepx,9),Irpdpd(Irrepx,10))

      I000 = 1
      I010 = I000 + Max_vo
      Iend = I010 + Max_vo 
      If (Iend .Gt. Memleft) Call Insmem("@-Tdvee_form_mu_0_0",
     +                                    Memleft,Iend)

C NDeP Eqn. 24, first-term, m(a,i)*Hbar(a,i)

      Do Ispin = 1, 1+Iuhf 
       
         Mu0_d = 0.0D0
         Lenvo = Irpdpd(Irrepx,8+Ispin)

         Call Getlst(Work(I000),1,1,1,Ispin,Mu0_hp_list)
         Call Getlst(Work(I010),1,1,1,Ispin,Hbr_hp_list)

         Mu0_d =  Mu0_d + Ddot(Lenvo,Work(I000),1,Work(1010),1)

        Write(6,*) 
CSSS        call checksum("Tdcc_form_mu_dot_0:",Work(I000),Lenvo,s)
CSSS        call checksum("Tdcc_form_mu_dot_0:",Work(I010),Lenvo,s)
        if (Ispin .Eq.1) Write(6,"(a,F12.6)") "The mu_0_dot alpha:",
     +                   Mu_0_d
        if (Ispin .Eq.2) Write(6,"(a,F12.6)") "The mu_0_dot  beta:",
     +                   Mu_0_d

      Enddo 

      Return
      End
