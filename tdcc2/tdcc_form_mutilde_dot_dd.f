










      Subroutine Tdcc_form_mutilde_dot_dd(Work,Memleft,Irrepx,Iuhf,
     +                               Mut_dot)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Double Precision Mut_dot_ab, Mut_dot_pp
      Integer D2_list
      Logical Rhf



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

      Rhf = (Iuhf.Eq.0)
      Mut_dot_ab = 0.0D0

C Mut_dot = Mut(Ab,Ij)*(e_I+e_j-e_A-e_b)

      Mut_list = 326
      D2_list  = 66

      Do Irrep_ij=1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)

c Since Irrepx can only be one; irrep_ab=Irrep_ij

            Ndim_ij_z = Irpdpd(Irrep_ij,14)
            Ndim_ab_z = Irpdpd(Irrep_ab,13)
            Ndim_ij_d = Ndim_ij_z
           
            Call Getlst(Work(I000),1,Ndim_ij_z,1,Irrep_ij,Mut_list)
            Call Getlst(Work(I010),1,Ndim_ij_d,1,Irrep_ij,D2_list)
        
            Mut_dot_ab = Mut_dot_ab + Ddot(Ndim_ij_z*Ndim_ab_z,
     +                                Work(I000),1,Work(I010),1)
  
      Enddo

C Mut_dot = Mut(AB,IJ)*(e_I+e_J-e_A-e_B) ISPIN=1
C Mut_dot = Mut(ab,ij)*(e_i+e_j-e_a-e_b) ISPIN=2
C Read Mut as Mut(A<B,I<J) and D(A<B,I<J)

      If (Ispin .Eq. 1) Write(6,"(a,1x,F12.6)")
     +                   "@-Tdcc_form_mutilde_dot_dd, AbAb:",
     +                       Mut_dot_ab
      If (Rhf) Return 
   
      Mut_dot_aa = 0.0D0
      Mut_dot_bb = 0.0D0

      Do Ispin = 1, 1+Iuhf

         Mut_list = 323 + Ispin 
         D2_list  = 63  + Ispin
 
         Do Irrep_ij=1, Nirrep
            Irrep_ab = Dirprd(Irrep_ij,Irrepx)

c Since Irrepx can only be one; irrep_ab=Irrep_ij

            Ndim_ij_z = Irpdpd(Irrep_ij,2+Ispin)
            Ndim_ab_z = Irpdpd(Irrep_ab,Ispin)
            Ndim_ij_d = Ndim_ij_z
         
            Call Getlst(Work(I000),1,Ndim_ij_z,1,Irrep_ij,Mut_list)
            Call Getlst(Work(I010),1,Ndim_ij_d,1,Irrep_ij,D2_list)

            If (Ispin .EQ. 1) Then
               Mut_dot_aa = Mut_dot_aa + Ddot(Ndim_ij_z*Ndim_ab_z,
     +                                   Work(I000),1,Work(I010),1)
            Else
               Mut_dot_bb = Mut_dot_bb + Ddot(Ndim_ij_z*Ndim_ab_z,
     +                                   Work(I000),1,Work(I010),1)
            Endif 
         Enddo 
         
      Write(6,*)
      If (Ispin .Eq. 1) Write(6,"(a,1x,F12.6)")
     +                   "@-Tdcc_form_mutilde_dot_dd, AAAA:",
     +                    Mut_dot_aa
      If (Ispin .Eq. 1) Write(6,"(a,1x,F12.6)")
     +                   "@-Tdcc_form_mutilde_dot_dd, BBBB:",
     +                    Mut_dot_bb

      Enddo 
      Return
      End
    
