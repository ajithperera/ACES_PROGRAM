










      Subroutine Tdcc_load_mu_vo(Work,Memleft,Memleft_modf,Iuhf,
     +                           Ioff_vo,Irrepx)

      Implicit Integer (A-Z)

      Dimension Work(Memleft), Ov_start(2)
      Dimension Ioff_vo(8,2)



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

C The Mu_vo vector reside at the end of the Work array ordered as
C UHF Mu_vo(AA)Mu_vo(BB)
C RHF Mu_vo(AA)
C The vo_off array give offset to the right-hand irrep. 

      List_mu0 = 390
      Memleft_modf = Memleft 
      Do Ispin = 2, 2-Iuhf, -1

         If (Iuhf .EQ. 0) Then
            Ndim_vo = Irpdpd(Irrepx,9) 
         Else
           If (Ispin .EQ. 2) Then
            Ndim_vo = Irpdpd(Irrepx,10) 
           Else
            Ndim_vo = Irpdpd(Irrepx,9) 
           Endif 
         Endif

         Ov_start(Ispin) = Memleft_modf-Ndim_vo+1
         If (Iuhf .EQ. 0) Ov_start(1) = Ov_start(2)
         Memleft_modf = Memleft_modf - Ndim_vo
         Ioff = Ov_start(Ispin)

         Do Irrep_i = 1, Nirrep
            Irrep_a = Dirprd(Irrep_i,Irrepx)
            Ioff_vo(Irrep_i,Ispin) = Ioff 
            If (Iuhf .EQ. 0) Ioff_vo(Irrep_i,1) = Ioff
            Ndim = Vrt(Irrep_a,Ispin) * Pop(Irrep_i,Ispin)
            Ioff = Ioff + Ndim 
         Enddo 

         If (Iuhf .EQ. 0) then
            Call Getlst(Work(Ov_start(Ispin)),1,1,1,1,List_mu0)
         Else 
            Call Getlst(Work(Ov_start(Ispin)),1,1,1,ISpin,List_mu0)
         Endif 

      Enddo 

      Return
      End
