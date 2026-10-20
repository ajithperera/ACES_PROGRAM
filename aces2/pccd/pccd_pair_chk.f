










      Subroutine Pccd_pair_chk(Work,Maxcor)
C Pair-form check (each RFO point): the orbital gradient of Pccd_ohess
C Imode 2 (from the pair densities n, A, S and the integrals, OHESGRAD)
C against that of the pair solver (Pccd_rsp_zgr, OBRTGRDA), two
C independent formulas (antisymmetric parts).
      Implicit Double Precision(A-H,O-Z)
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



      Call Getrec(0,"JOBARC","OHESGRAD",Length,Junk)
      If (Length .Le. 0) Return
      Nsq  = Length/Iintfp
      Nbas = Nint(Dsqrt(Dble(Nsq)))
      If (2*Nsq .Gt. Maxcor) Return
      Call Getrec(20,"JOBARC","OHESGRAD",Nsq*Iintfp,Work)
      Call Getrec(20,"JOBARC","OBRTGRDA",Nsq*Iintfp,Work(Nsq+1))
      Dif = 0.0D0
      Gmx = 0.0D0
      Do Iq = 1, Nbas
         Do Ip = 1, Nbas
            Ipq = Ip + (Iq-1)*Nbas
            Iqp = Iq + (Ip-1)*Nbas
            Gpq = Work(Nsq+Ipq) - Work(Nsq+Iqp)
            Dif = Max(Dif,Dabs(Work(Ipq) - Gpq))
            Gmx = Max(Gmx,Dabs(Gpq))
         Enddo
      Enddo
      Write(6,"(2a,1p,E10.3,a,E10.3)") " @-PAIR gradient check,",
     +     " densities vs pair form: max|dG| =",Dif,"  max|G| =",Gmx
      Return
      End
