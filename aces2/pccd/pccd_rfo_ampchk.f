















































































































































































































      Subroutine Pccd_rfo_ampchk(Work,Maxcor,List,Name,Ibad)
C
C OO_METHOD=RFO guard against unphysical amplitude solutions. Away from
C equilibrium the orbital steps can carry the pCCD solution onto an
C unphysical branch: the T2 pair amplitudes grow without bound while
C the energy slides far below the physical minimum (N2 cc-pVDZ 2.5 and
C 3.0 A: max|t| 42 and 48), or the amplitude Jacobian approaches a fold
C and the lambda equations blow up while the energy stays flat (N2
C 1.5 A: analytic gradient 1e-2 -> 5e3 over 60 iterations, max|L2|
C 6.9). At the verified minima max|t| <= 0.9, but intermediate points
C of paths that reach them go up to max|T2| 2.3 and max|L2| 1.7 (N2
C cc-pVDZ 1.5-2.5 A, H2O cc-pVDZ 2Re); Ampcap = 5 lies between. Both
C bounds are from N2 and H2O only. Ibad = 1 when max|amplitude| of
C List (46: T2, 146: L2) exceeds Ampcap; the caller then treats the
C point like unconverged amplitudes. No-op (Ibad = 0) unless
C OO_METHOD=RFO. RHF.
C
      Implicit Double Precision(A-H,O-Z)
      Character*2 Name
      Dimension Work(Maxcor)


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



c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
      Parameter (Ampcap=5.0D0)

      Ibad = 0
      If (Iflags2(178) .Ne. 6) Return
      Nsize = Isymsz(Isytyp(1,List),Isytyp(2,List))
      If (Nsize .Gt. Maxcor) Call Insmem("Pccd_rfo_ampchk",Nsize,Maxcor)
      Call Getall(Work,Nsize,1,List)
      Amax = 0.0D0
      Anrm = 0.0D0
      Do I = 1, Nsize
         Amax = Max(Amax,Dabs(Work(I)))
         Anrm = Anrm + Work(I)*Work(I)
      Enddo
      Anrm = Dsqrt(Anrm)
      Write(6,"(a,a2,a,F12.6,a,F12.6)") " @-RFO ",Name," max|amp| =",
     +      Amax,"  norm =",Anrm
      If (List .Eq. 46) Call Putrec(20,"JOBARC","RFO_T2MX",Iintfp,Amax)
      If (Amax .Gt. Ampcap) Then
         Ibad = 1
         Write(6,"(a,a2,a,F6.2,a)") " @-RFO ",Name," amplitudes exceed",
     +         Ampcap," : unphysical solution, point rejected."
      Endif
      Return
      End
