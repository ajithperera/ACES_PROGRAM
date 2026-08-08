










      Subroutine Tdcc_built_ac_func(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Memleft)
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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
     
C Singles left, right dipole functions and derivative lists; 
C Mu : Right
C Mut: Left 

      Mu_s      = 390
      Mu_s_dot  = 394 
      Mut_s     = 392 
      Mut_s_dot = 396

C Corresponding doubles lists; There are no right dipole doubles
C list 

      Mut_d_aa = 324
      Mut_d_bb = 325
      Mut_d_ab = 326
      
      Mu_d_dot_aa = 334
      Mu_d_dot_bb = 335
      Mu_d_dot_ab = 336

      Mut_d_dot_aa = 344
      Mut_d_dot_bb = 345
      Mut_d_dot_ab = 346

      Rhf = (Iuhf .EQ. 0)

      If (Rhf) Then 

         Ndim_ai = Irpdpd(Irrepx,9)
         Ndim_pphh_ab = Idsymsz(Irrepx,13,14)

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_ai
         I020 = I010 + Ndim_ai
         I030 = I020 + Ndim_ai
         I040 = I030 + Ndim_ai
         I050 = I040 + Ndim_pphh_ab
         I060 = I050 + Ndim_pphh_ab
         I070 = I060 + Ndim_pphh_ab
         Iend = I070
         If (Iend .GT. Memleft) Then
            Call Insmem("@-Tdcc_built_ac_func",Iend,Memleft)
         Endif

         Call Getlst(Work(I000),1,1,1,1,Mu_s) 
         Call Getlst(Work(I010),1,1,1,1,Mu_s_dot) 
         Call Getlst(Work(I020),1,1,1,1,Mut_s) 
         Call Getlst(Work(I030),1,1,1,1,Mut_s_dot) 

         Write(6,*)
         call checksum("Mu_s(a,i)    :",Work(I000),Ndim_ai,s)
         call checksum("Mu_s_dot(a,i):",Work(I010),Ndim_ai,s)
         call checksum("Mut_s(a,i)   :",Work(I020),Ndim_ai,s)
         call checksum("Mut_s_dot(a,i):",Work(I030),Ndim_ai,s)
         Call Getall(Work(I040),Ndim_pphh_ab,Irrepx,Mut_d_ab)
         Call Getall(Work(I050),Ndim_pphh_ab,Irrepx,Mu_d_dot_ab)
         Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mut_d_dot_ab)

         Write(6,*)
         call checksum("Mut_d(Ab,Ij)     :",Work(I040),Ndim_pphh_ab,s)
         call checksum("Mu_d_dot(Ab,Ij) :",Work(I050),Ndim_pphh_ab,s)
         call checksum("Mut_d_dot(Ai,Bj):",Work(I060),Ndim_pphh_ab,s)
         Return
      Endif 

      Do Ispin = 1, (Iuhf+1)
         Ndim_ai = Irpdpd(Irrepx,8+Ispin)

         Ndim_pphh_aa = Idsymsz(Irrepx,1,3)
         Ndim_pphh_bb = Idsymsz(Irrepx,2,4)
         Ndim_pphh_ab = Idsymsz(Irrepx,13,14)

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_ai
         I020 = I010 + Ndim_ai
         I030 = I020 + Ndim_ai
         I040 = I030 + Ndim_ai
         I050 = I040 + Ndim_pphh_aa
         I060 = I050 + Ndim_pphh_bb
         I070 = I060 + Ndim_pphh_ab
         I080 = I070 + Ndim_pphh_aa
         I090 = I080 + Ndim_pphh_bb
         I100 = I090 + Ndim_pphh_ab
         Iend = I070
         If (Iend .GT. Memleft) Then
            Call Insmem("@-Tdcc_built_ac_func",Iend,Memleft)
         Endif

         Call Getlst(Work(I000),1,1,1,Ispin,Mu_s) 
         Call Getlst(Work(I010),1,1,1,Ispin,Mu_s_dot)
         Call Getlst(Work(I020),1,1,1,Ispin,Mut_s)
         Call Getlst(Work(I030),1,1,1,Ispin,Mut_s_dot)

         write(6,*)
         call checksum("Mu_s(a,i)    :",Work(I000),Ndim_ai,s)
         call checksum("Mu_s_dot(a,i):",Work(I010),Ndim_ai,s)
         call checksum("Mut_s(a,i)   :",Work(I020),Ndim_ai,s)
         call checksum("Mut_s_dot(a,i):",Work(I030),Ndim_ai,s)
         Call Getall(Work(I040),Ndim_pphh_aa,Irrepx,Mut_d_aa)
         Call Getall(Work(I050),Ndim_pphh_bb,Irrepx,Mut_d_bb)
         Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mut_d_ab)

         write(6,*)
         call checksum("Mut_d(A<B,I<J):",Work(I040),Ndim_pphh_aa,s)
         call checksum("Mut_d(a<b,i<j):",Work(I050),Ndim_pphh_bb,s)
         call checksum("Mut_d(Ab,Ij)  :",Work(I060),Ndim_pphh_ab,s)
         Call Getall(Work(I040),Ndim_pphh_aa,Irrepx,Mut_d_dot_aa)
         Call Getall(Work(I050),Ndim_pphh_bb,Irrepx,Mut_d_dot_bb)
         Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mut_d_dot_ab)

         Write(6,*)
         call checksum("Mut_d_dot(A<B,I<J):",Work(I040),Ndim_pphh_aa,s)
         call checksum("Mut_d_dot(a<b,i<j):",Work(I050),Ndim_pphh_bb,s)
         call checksum("Mut_d_dot(Ab,Ij)  :",Work(I060),Ndim_pphh_ab,s)
         Call Getall(Work(I040),Ndim_pphh_aa,Irrepx,Mu_d_dot_aa)
         Call Getall(Work(I050),Ndim_pphh_bb,Irrepx,Mu_d_dot_bb)
         Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mu_d_dot_ab)

         Write(6,*)
         call checksum("Mu_d_dot(A<B,I<J):",Work(I040),Ndim_pphh_aa,s)
         call checksum("Mu_d_dot(a<b,i<j):",Work(I050),Ndim_pphh_bb,s)
         call checksum("Mu_d_dot(Ab,Ij)  :",Work(I060),Ndim_pphh_ab,s)
      Enddo

      Return 
      End
