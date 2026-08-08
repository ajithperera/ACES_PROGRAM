










      Subroutine Reorient_cfc(CFc,R,Work,Maxcor,Ndims,Nreals,Nmodes,
     +                        Imode)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Cfc(Ndims,Ndims,Nmodes),R(3,3)



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




      Data Ione /1/

      I000 = Ione 
      I010 = I000 + Ndims*Ndims 
      Iend = I010 + Ndims*Ndims 

      Write(6,"(a)") "The incoming cubic force constant matrix"
      call output(Cfc(1,1,Imode),1,Ndims,1,Ndims,Ndims,Ndims,1)

      Call Big_reori_mat(Work(I000),R,Nreals)

      Write(6,*)
      Write(6,"(a)") "The reorientation matrix"
      call output(R,1,3,1,3,3,3,1)
      Write(6,"(a)") "The big reorientation matrix"
      call output(Work(I000),1,Ndims,1,Ndims,Ndims,Ndims,1)

      Call Reorder(Cfc(1,1,Imode),Work(I010),Work(I000),Ndims)

      Write(6,"(a)") "The reoriented cubic force constant matrix"
      call output(Cfc(1,1,Imode),1,Ndims,1,Ndims,Ndims,Ndims,1)

      Return
      End
