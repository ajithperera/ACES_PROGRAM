










      Subroutine Reorient_grd(Grd,R,Work,Maxcor,Nreals)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Grd(Nreals,3),R(3,3)



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




      Data Ione,Done,Dnull /1,1.0D0,0.0D0/

      I000 = Ione 
      Iend = I000 + Ndims*3

      Write(6,*)
      Write(6,"(a)") "The reorientation matrix"
      call output(R,1,3,1,3,3,3,1)

      Call Dgemm("N","N",Nreals,3,3,Done,Grd,Nreals,R,3,Dnull,
     +           Work(I000),Nreals)
      Call Dcopy(Nreals*3,Work(I000),1,Grd,1)

       Write(6,*)
       Write(6,"(a)") "The reoriented gradients"
       Do I = 1, Nreals 
       Write(6,"(6(1x,F15.6))") (Grd(I,J),J=1,3)
       Enddo

      Return
      End
