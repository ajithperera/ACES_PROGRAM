










      Subroutine Tdee_pes_ode(T,Data_in_out,Mubar_dot)

      Implicit Double Precision (A-H,O-Z)
      Double Precision Mubar_dot, Mone
      Logical Real


      Common /Tdee_vars/Nsize_dummy,Irrepx_dummy,Iside_dummy,
     +                  Iuhf_dummy,Memleft_Dummy




c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end


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




      Dimension Data_in_out(2*Nsize_dummy),Mubar_dot(2*Nsize_dummy)
      Data Mone /-1.0D0/

      Irrepx  = Irrepx_dummy
      Iside   = Iside_dummy
      Iuhf    = Iuhf_dummy 
      Nsize   = Nsize_dummy
      Memleft = Memleft_dummy 

      Write(6,*)
      Write(6,"(a)") "Entered - tdee_pes_ode"
      Write(6,"(a,3i2,3i10)") "The global vars: ", Irrepx,Iside,
     +                        Iuhf,Nsize,Memleft
CSSS      call prvecr(Data_in_out,2*Nsize)
CSSS      call checksum("@-tdee_pes_ode;Mubar:",Data_in_out,2*Nsize,s)

C First do the -i Hbar * Mubar_R and put as the imaginary component 
C of the derivative (Real=.True. means we are working on the real 
C component of Mubar(t)).

      Real = .True. 

      Call Tdee_get_derivatives(Irrepx,Iside,Iuhf,Memleft,Nsize,
     +                          Real,Data_in_out,Mubar_dot(Nsize+1))

C Now do the -i*i Hbar * Mbar_I and put as the real component 
C of the derivative.

      Real = .False.

      Call Tdee_get_derivatives(Irrepx,Iside,Iuhf,Memleft,Nsize,
     +                          Real,Data_in_out(Nsize+1),Mubar_dot)

CSSS      call checksum("@-tdee_pes_ode;Mubar_dot:",Mubar_dot,2*Nsize,s)

      Return 
      End

 
