#include "flags.h"
      Subroutine Erg_dens(Iuhf,Lcorr,N,B,Ha,Hb,Xa,Xb,Dt,W,Dso,C,Occ,
     &                    Eps,Scr)
C
C Densities over ERD functions, D = B D(SO) B(T) (Erg_bso). Per spin
C
C    H   = C n C(T)        the SCF (reference) density, n = ORB_OCC*
C
C and for an SCF gradient W = sum_spin C n e C(T) (e = SCFEVL*); for
C ROHF with standard (not semicanonical) orbitals, where the occupied
C block of the Fock matrix is not diagonal, W = sum_spin H F H with F
C the SO Fock matrix FOCKA/FOCKB (xvdint's Formfij). For
C a correlated gradient (Lcorr) the relaxed density R (RELDENS*, MO
C basis) and the intermediates I (IINTERM*) give, as in xvdint's
C Onedso and Dsofso,
C
C    X   = C R C(T) - H    the correlation part of the density
C    W   = - C sym(I) C(T)
C
C Once the correlated programs have run (Lcord: a correlated or a
C second derivative calculation), SCFEVC* and SCFEVL* are in the
C correlated MO order, all occupied orbitals first, like R
C and I (xvdint's Symcext and Symd turn all of them back to the SCF
C order instead). ORB_OCC* stays in the SCF order, so in that case
C the occupied orbitals are the first ones, as many as OCCUPY*0 says
C (SYMPOPO* would leave out dropped core orbitals, which are still
C there).
C For RHF, RELDENSA and IINTERMA are the totals (both spins): Xa = Xb
C is half of C R C(T) - 2 H. Ha, Hb, Xa, Xb are per spin, Dt is the
C total one-particle density (H + X over both spins), W the total
C energy weighted density, Dso the total one-particle density in the
C SO basis (for the ECP gradient).
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"
#include "machsp.com"
#include "flags.com"

      Logical Lcorr, Lcord
      Dimension B(N,N), Ha(N,N), Hb(N,N), Xa(N,N), Xb(N,N), Dt(N,N),
     &          W(N,N), Dso(N,N), C(N,N), Occ(N), Eps(N), Scr(N,N,2)
      Integer Iocc(8,2)
      Logical Lrohf
      Character*8 Labc(2), Labe(2), Labo(2), Labd(2), Labi(2), Labf(2)
      Data Labc /'SCFEVCA0','SCFEVCB0'/
      Data Labe /'SCFEVLA0','SCFEVLB0'/
      Data Labo /'ORB_OCCA','ORB_OCCB'/
      Data Labd /'RELDENSA','RELDENSB'/
      Data Labi /'IINTERMA','IINTERMB'/
      Data Labf /'FOCKA   ','FOCKB   '/

      Lcord = Lcorr .Or. Iflags(h_IFLAGS_deriv_lev) .Eq. 2
      Lrohf = Iflags(h_IFLAGS_ref) .Eq. 2 .And.
     &        Iflags(h_IFLAGS_orbitals) .Eq. 0
      If (Lcord) Then
         Call Getrec(20,'JOBARC','OCCUPYA0',Nirr,Iocc(1,1))
         If (Iuhf .Eq. 1) Then
            Call Getrec(20,'JOBARC','OCCUPYB0',Nirr,Iocc(1,2))
         Endif
      Endif
      Call Zero(W,N*N)
      Call Zero(Dso,N*N)
      Call Zero(Xa,N*N)
      Call Zero(Xb,N*N)
      Fac = 1.0D0
      If (Iuhf .Eq. 0) Fac = 2.0D0
      Do Ispin = 1, Iuhf+1
         Call Getrec(20,'JOBARC',Labc(Ispin),N*N*Iintfp,C)
         Call Getrec(20,'JOBARC',Labe(Ispin),N*Iintfp,Eps)
         If (Lcord) Then
            Nocc = 0
            Do Ir = 1, Nirr
               Nocc = Nocc + Iocc(Ir,Ispin)
            Enddo
            Do J = 1, N
               Occ(J) = 0.0D0
               If (J .Le. Nocc) Occ(J) = 1.0D0
            Enddo
         Else
            Call Getrec(20,'JOBARC',Labo(Ispin),N*Iintfp,Occ)
         Endif

C H = B C n C(T) B(T).

         Do J = 1, N
            Do I = 1, N
               Scr(I,J,2) = Occ(J)*C(I,J)
            Enddo
         Enddo
         Call Dgemm('N','T',N,N,N,1.0D0,Scr(1,1,2),N,C,N,0.0D0,
     &              Scr(1,1,1),N)
         If (.Not. Lcorr) Call Daxpy(N*N,Fac,Scr(1,1,1),1,Dso,1)
         If (Ispin .Eq. 1) Then
            Call Erg_btr(N,B,Scr(1,1,1),Ha,Scr(1,1,2),0.0D0)
         Else
            Call Erg_btr(N,B,Scr(1,1,1),Hb,Scr(1,1,2),0.0D0)
         Endif

         If (.Not. Lcorr) Then

C W += B C (fac n e) C(T) B(T), or B H F H B(T) for ROHF.

            If (Lrohf) Then
               Call Getrec(20,'JOBARC',Labf(Ispin),N*N*Iintfp,C)
               Call Dgemm('N','N',N,N,N,1.0D0,Scr(1,1,1),N,C,N,0.0D0,
     &                    Scr(1,1,2),N)
               Call Dgemm('N','N',N,N,N,1.0D0,Scr(1,1,2),N,Scr(1,1,1),
     &                    N,0.0D0,C,N)
               Call Erg_btr(N,B,C,W,Scr(1,1,2),1.0D0)
               Goto 100
            Endif
            Do J = 1, N
               Do I = 1, N
                  Scr(I,J,2) = Fac*Occ(J)*Eps(J)*C(I,J)
               Enddo
            Enddo
            Call Dgemm('N','T',N,N,N,1.0D0,Scr(1,1,2),N,C,N,0.0D0,
     &                 Scr(1,1,1),N)
            Call Erg_btr(N,B,Scr(1,1,1),W,Scr(1,1,2),1.0D0)
  100       Continue
         Else

C X = B C R C(T) B(T) / fac - H, W -= B C sym(I) C(T) B(T).

            Do Ipass = 1, 2
               If (Ipass .Eq. 1) Then
                  Call Getrec(20,'JOBARC',Labd(Ispin),N*N*Iintfp,
     &                        Scr(1,1,1))
               Else
                  Call Getrec(20,'JOBARC',Labi(Ispin),N*N*Iintfp,
     &                        Scr(1,1,1))
                  Do J = 1, N
                     Do I = 1, J-1
                        X = 0.5D0*(Scr(I,J,1) + Scr(J,I,1))
                        Scr(I,J,1) = X
                        Scr(J,I,1) = X
                     Enddo
                  Enddo
               Endif
               Call Dgemm('N','N',N,N,N,1.0D0,C,N,Scr(1,1,1),N,
     &                    0.0D0,Scr(1,1,2),N)
               Call Dgemm('N','T',N,N,N,1.0D0,Scr(1,1,2),N,C,N,
     &                    0.0D0,Scr(1,1,1),N)
               If (Ipass .Eq. 1) Then
                  Call Daxpy(N*N,1.0D0,Scr(1,1,1),1,Dso,1)
                  Call Dscal(N*N,1.0D0/Fac,Scr(1,1,1),1)
                  If (Ispin .Eq. 1) Then
                     Call Erg_btr(N,B,Scr(1,1,1),Xa,Scr(1,1,2),0.0D0)
                     Call Daxpy(N*N,-1.0D0,Ha,1,Xa,1)
                  Else
                     Call Erg_btr(N,B,Scr(1,1,1),Xb,Scr(1,1,2),0.0D0)
                     Call Daxpy(N*N,-1.0D0,Hb,1,Xb,1)
                  Endif
               Else
                  Call Dscal(N*N,-1.0D0,Scr(1,1,1),1)
                  Call Erg_btr(N,B,Scr(1,1,1),W,Scr(1,1,2),1.0D0)
               Endif
            Enddo
         Endif
      Enddo
      If (Iuhf .Eq. 0) Then
         Call Dcopy(N*N,Ha,1,Hb,1)
         Call Dcopy(N*N,Xa,1,Xb,1)
      Endif

      Call Dcopy(N*N,Ha,1,Dt,1)
      Call Daxpy(N*N,1.0D0,Hb,1,Dt,1)
      Call Daxpy(N*N,1.0D0,Xa,1,Dt,1)
      Call Daxpy(N*N,1.0D0,Xb,1,Dt,1)

      Return
      End
