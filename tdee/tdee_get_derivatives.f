










      Subroutine Tdee_get_derivatives(Irrepx,Iside,Iuhf,Memleft,Nsize,
     +                                Type,Dummy1,Dummy2)

      Implicit Double Precision (A-H,O-Z)
      Logical Type 
      Dimension Dummy1(Nsize),Dummy2(Nsize)
      Double Precision Mone



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




      Data Mone /-1.0D0/
    
C At t=0 (the initial state mubar(0) is already on the correct lists.
C Multiplication is entirly stand-alone and works from the data
C from the following lists.
C UHF   Right Left
C List: 390   392 AA(1),BB(2)
C List: 314   324 AAAA
C List: 315   325 BBBB
C List: 316   326 ABAB
C RHF
C List: 390   392 AA(1)
C List: 316   326 ABAB
C
C First do the -i Hbar * Mbar_R and put it in imaginary component

      Mubar_s_pq_t0 = 390
      Mubar_d_pq_t0 = 313

      If (Iside .EQ. 1) Then
         Ioffr1 = 0
         Ioffr2 = 0
         Ioffsp = 0

C This write Mubar(t_n),  iside=1

         Call Tdee_dump_vec(Irrepx,Dummy1,Nsize,Mubar_s_pq_t0,
     +                      Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
     +                      Iuhf,.False.)

      Elseif (Iside .EQ. 2) Then
         Ioffr1 = 2
         Ioffr2 = 10
         Ioffsp = 0

C This write Mubar(t_n),  iside=2

         Call Tdee_dump_vec(Irrepx,Dummy1,Nsize,Mubar_s_pq_t0,
     +                      Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
     +                      Iuhf,.False.)
      Endif

      Write(6,*)
      If (type) Then
         write(6,"(a)") "The real part of the Mbar(t)"
         Call checksum("tdee_pes_ode_in:",Dummy1,Nsize,s)
      else
         write(6,"(a)") "The imginary part of the Mbar(t)"
         Call checksum("tdee_pes_ode_in:",Dummy1,Nsize,s)
      endif 

      Call Tdee_return_mubar_dot(Icore(I0),Memleft,Iuhf,Irrepx,
     +                           Iside)

C The data returned from the multiplication is stored as

C The data returned from the multiplication is stored as
C follows.
C UHF   Right Left
C List: 394   396 AA(1),BB(2)
C List: 334   344 AAAA
C List: 335   345 BBBB
C List: 336   346 ABAB
C RHF
C List: 394   396 AA(1)
C List: 336   346 ABAB

      Mubar_s_pq_tn = 394
      Mubar_d_pq_tn = 333

      If (Iside .EQ. 1) Then

        Ioffr1 = 0
        Ioffr2 = 0

C This Returns Mubar(t_n+1), iside=1

        Call Tdee_Load_vec(Irrepx,Dummy2,Nsize,Mubar_s_pq_tn,
     +                     Ioffr1,Mubar_d_pq_tn,Ioffr2,Iuhf,.False.)


      Elseif (Iside .EQ.2) Then

        Ioffr1 = 2
        Ioffr2 = 10

C This Returns Mubar(t_n+1), iside=2

        Call Tdee_Load_vec(Irrepx,Dummy2,Nsize,Mubar_s_pq_tn,
     +                     Ioffr1,Mubar_d_pq_tn,Ioffr2,Iuhf,.False.)

      Endif

      If (Type) Call Dscal(Nsize,Mone,Dummy2,1)

      Write(6,*)
      If (type) Then
         write(6,"(a)") "The real part of the Mbar_dot(t)"
         Call checksum("tdee_pes_ode_in:",Dummy2,Nsize,s)
      else
         write(6,"(a)") "The imginary part of the Mbar_dot(t)"
         Call checksum("tdee_pes_ode_in:",Dummy2,Nsize,s)
      endif 

      Return
      End

