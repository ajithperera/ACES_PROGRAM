













































































































































































































      Subroutine Tdcc_form_mu_0_s(Work,Memleft,Doo,Dvv,Dvo,Mu_0_S,
     +                            Irrepx,Lenoo,Lenvv,Lenvo,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Doo(Lenoo),Dvv(Lenvv),Dvo(Lenvo)

      Double Precision Mu_0_s(Lenvo)
      Double Precision Mone
      Integer dissiz_aa,dissiz_bb,dissiz_ab
 
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
      Mu0_ph_list = 390

      Do Ispin = 1, 1+Iuhf 
         
         Ivo_0 = 1 + (Ispin-1) * Irpdpd(Irrepx,9)

         Do Irrep_i = 1, Nirrep
            Irrep_a = Dirprd(Irrep_i,Irrepx)

            Ndim_i = Pop(Irrep_i,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
               
            Ivo_z = Ivo_0 + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)
            Ivo_d = Ivo_z

C z(A,I) = Mu(A,I) (and Z(a,i) = Mu(a,i))

            Ndim_ai = Ndim_i * Ndim_a
            Call Dcopy(Ndim_ai,Dvo(Ivo_d),1,Mu_0_S(Ivo_z),1)

         Enddo

         Call Putlst(Mu_0_S(Ivo_0),1,1,1,Ispin,Mu0_ph_list)

      Enddo

      do ispin=1, 1+iuhf
      Call Getlst(Mu_0_S,1,1,1,Ispin,Mu0_ph_list)
      ioff = 1 + (Ispin-1)*Irpdpd(Irrepx,9)
      do irrep_i=1,Nirrep
         irrep_a=dirprd(irrep_i,Irrepx)
         ndim_i = Pop(irrep_i,ispin)
         Ndim_a = Vrt(irrep_a,ispin)
         Write(6,*)
         write(6,"(a,a,5(1x,i2))") "irrep_i,Irrep_a,Ndim_i,",
     +                     "Ndim_a,Ioff: ",irrep_i,Irrep_a,
     +                      Ndim_i,Ndim_a,Ioff
         call output(Mu_0_S(ioff),1,ndim_a,1,Ndim_i,Ndim_a,
     +                Ndim_i,1)
         Ioff = Ioff + ndim_i*Ndim_a
      enddo
      enddo

      Return
      End
     
           
        
