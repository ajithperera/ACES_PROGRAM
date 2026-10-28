      Subroutine Pccd_ohess_lst(Buf,Lbuf,List,Is1,Is2,Is3,Is4,Nocc,
     +                          Mode,W,A,S,Dmax,H,Xj,Xk,Kp,Nbas,Npar)
C
C Visit every element of the totally symmetric list List (RHF, rows pq,
C columns rs, each DPD pair ordered by the irrep of its second index,
C first index fastest; index spaces Is = 1 occupied, 2 virtual) with
C its MOs numbered over all MOs (occupied by irrep, then virtual by
C irrep, as Irp of Pccd_rfo_irp). The columns are read in batches that
C fit Buf(Lbuf) (Lbuf at least one column).
C   Mode 1, 2: density list; W times it goes into the seniority-zero
C              A and S (Pccd_ohess_s0el; Dmax: largest dropped value);
C   Mode 3:    integral list <pq|rs> = (pr|qs), into the frozen Hessian
C              H and the integral classes Xj, Xk (Pccd_ohess_scat);
C   Mode 4:    as 3, Xj and Xk only (H is not touched).
C
      Implicit Double Precision(A-H,O-Z)
      Dimension Buf(Lbuf),Is(4),Iof(8,2),Nsp(8,2),Igo(2)
      Dimension A(*),S(*),H(*),Xj(*),Xk(*),Kp(*)
#include "sym.com"
#include "syminf.com"
#include "sympop.com"

      Is(1) = Is1
      Is(2) = Is2
      Is(3) = Is3
      Is(4) = Is4
      Iob = 0
      Ivb = 0
      Do Irr = 1, Nirrep
         Iof(Irr,1) = Iob
         Iof(Irr,2) = Ivb
         Nsp(Irr,1) = Pop(Irr,1)
         Nsp(Irr,2) = Vrt(Irr,1)
         Iob = Iob + Pop(Irr,1)
         Ivb = Ivb + Vrt(Irr,1)
      Enddo
      Igo(1) = 0
      Igo(2) = Nocc

      Do Ix = 1, Nirrep
         Nrow = Irpdpd(Ix,Isytyp(1,List))
         Ncol = Irpdpd(Ix,Isytyp(2,List))
         If (Nrow*Ncol .Eq. 0) Go to 10
         Ncb   = Min(Ncol,Lbuf/Nrow)
         Icol  = 0
         Icol0 = 1
         Icol1 = 0
         Do I4r = 1, Nirrep
         I3r = Dirprd(I4r,Ix)
         Do I4 = 1, Nsp(I4r,Is(4))
         Is_ = Igo(Is(4)) + Iof(I4r,Is(4)) + I4
         Do I3 = 1, Nsp(I3r,Is(3))
            Ir = Igo(Is(3)) + Iof(I3r,Is(3)) + I3
            Icol = Icol + 1
            If (Icol .Gt. Icol1) Then
               Icol0 = Icol
               Nc    = Min(Ncb,Ncol-Icol+1)
               Call Getlst(Buf,Icol0,Nc,1,Ix,List)
               Icol1 = Icol + Nc - 1
            Endif
            Ioff = (Icol-Icol0)*Nrow
            Irow = 0
            Do I2r = 1, Nirrep
            I1r = Dirprd(I2r,Ix)
            Do I2 = 1, Nsp(I2r,Is(2))
            Iq = Igo(Is(2)) + Iof(I2r,Is(2)) + I2
            Do I1 = 1, Nsp(I1r,Is(1))
               Ip = Igo(Is(1)) + Iof(I1r,Is(1)) + I1
               Irow = Irow + 1
               X = Buf(Ioff+Irow)
               If (X .Ne. 0.0D0) Then
                  If (Mode .Ge. 3) Then
                     Call Pccd_ohess_scat(Ip,Ir,Iq,Is_,X,Is(1),Is(3),
     +                                    Is(2),Is(4),A,S,H,Xj,Xk,Kp,
     +                                    Nbas,Npar,Mode-3)
                  Else
                     Call Pccd_ohess_s0el(Ip,Iq,Ir,Is_,Mode,W*X,A,S,
     +                                    Nbas,Dmax)
                  Endif
               Endif
            Enddo
            Enddo
            Enddo
         Enddo
         Enddo
         Enddo
 10      Continue
      Enddo
      Return
      End
