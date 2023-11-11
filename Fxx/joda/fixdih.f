      Subroutine Fixdih(Iangle,Nx)

      Implicit Double Precision(A-H,O-Z)
      Integer Dep_ind 
      Logical Dep_dih  

C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end



      Data Done,Ione /1.0D0,1/

      Write(6,"(a,6(1x,I3))") "The connectivity pattern: ",(Ncon(i),
     &                         i=1,Nx)
      Icon1 = Ncon(Iangle)
      Icon2 = Ncon(Iangle-1)

      Write(6,"(a,2i2)") " The atom making the angle: ", Icon1,Icon2

      Ioff = 1
      Do Iatom = 2, (NX/3+1)
         Dep_dih = ((Ncon(Ioff+2) .Eq. Icon1 .And. Ncon(Ioff+3) .Eq.
     +               Icon2) .Or. (Ncon(Ioff+2) .Eq. Icon2 .And.
     +               Ncon(Ioff+3) .Eq.Icon1))
         If (Dep_dih) Dep_ind =  Ioff + 3
         Ioff = (Iatom - 1)*3 
      Enddo 

      Write(6,"(a,2i2)") " The dependent diherdral: ", Dep_ind
      R(Dep_ind) = Dacos(-Done) - R(Dep_ind)
      Call Putrec(20,"JOBARC","DIH_FIX ",Ione,Ione)

      Return
      End 
