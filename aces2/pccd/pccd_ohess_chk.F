      Subroutine Pccd_ohess_chk(Work,Maxcor)
C
C Debug (file pccd_ohess.in, see Pccd_do_orbrots), after
C Pccd_form_htau_pq: compare the orbital gradient of Pccd_ohess
C (OHESGRAD) with OBRTGRDA (antisymmetric parts) and write the
C analytic frozen Hessian OHESFRZ to pccd_ohess.dat (Npar, then the
C matrix, rows fastest).
C
      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)
#include "machsp.com"

      Call Getrec(0,"JOBARC","OHESGRAD",Length,Junk)
      If (Length .Le. 0) Return
      Nsq  = Length/Iintfp
      Nbas = Nint(Dsqrt(Dble(Nsq)))
      Call Getrec(0,"JOBARC","OHESFRZ ",Length,Junk)
      Nhs  = Length/Iintfp
      Npar = Nint(Dsqrt(Dble(Nhs)))
      If (2*Nsq+Nhs .Gt. Maxcor) Call Insmem("Pccd_ohess_chk",
     +                                       2*Nsq+Nhs,Maxcor)
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
      Write(6,"(2a,1p,E10.3,a,E10.3)") " @-Pccd_ohess_chk:",
     +     " max|G(ohess) - G(OBRTGRDA)| =",Dif,"  max|G| =",Gmx
      Call Getrec(20,"JOBARC","OHESFRZ ",Nhs*Iintfp,Work(2*Nsq+1))
      Open(Unit=43,File="pccd_ohess.dat",Status="Unknown")
      Write(43,*) Npar
      Write(43,"(4(1x,E25.17))") (Work(2*Nsq+I),I=1,Nhs)
      Close(43)
      Return
      End
