      Subroutine Form_genrlizd_oscstrength(Work,Maxcor,Iuhf,Eigval)
C
C State-to-state electric-dipole transition moments and oscillator
C strengths from the generalized transition densities written by
C form_genrlizd_tdens ("GTDENSll" = spin-summed AO densities
C <L_k|{p+ q}|R_l> for all k, one record per R state l):
C
C   Tm(c; k<-l) = <L_k| mu_c |R_l>
C   f(k,l)      = 2/3 |E_k - E_l| Sum_c Tm(c; k<-l) Tm(c; l<-k)
C
C (the EOM-CC non-symmetric form). Tm_c(Jroot,Irrepy,Iroot,Irrepx)
C is the moment with L state (Jroot,Irrepy) and R state
C (Iroot,Irrepx). GEN_OSCS holds f in the loop order
C Irrepx, Iroot, Irrepy, Jroot (read by tablee/tablee_bd).
C
C 11/2026 (Ajith Perera): rewritten with the AO-density records of
C form_genrlizd_tdens; dipole integrals read as in edens.F (VPROP
C triangle or symmetry-packed); the 3x1M-element pair-of-pairs
C scratch arrays (osc_init/osc_fin) replaced by the direct product.
C
      Implicit Double Precision(A-H,O-Z)

      Character*8 Label_xyz(3),Label
      Character*2 Suffix(32)
      Logical Vprop

      Dimension Eigval(100,8)
      Dimension Tm_x(4,8,4,8), Tm_y(4,8,4,8), Tm_z(4,8,4,8)
      Dimension Osc(1024)
      Dimension Work(Maxcor)

#include "flags.com"
#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"
#include "machsp.com"

      Common /Calcinfo/Nroot(8)
      Common /Vdintlen/ Lens(8),Lent(8)
      Common /Aosym/Iaopop(8),Ioffao(8),Ioffv(8,2),Ioffo(8,2),
     &              Irpdpdao(8),Irpdpdaos(8),
     &              Istart(8,8),Istartmo(8,3)
      Common /Vdintprt/Ntpert,Npert(8),Kpert(8),Idipsm(3),
     &                 Iyzpert,Ixzpert,Ixypert,Itransx,
     &                 Itransy,Itransz,Nucind
      Common /Intprog/ Vprop

      Data Label_xyz/'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      Data Suffix /"01","02","03","04","05","06","07","08","09","10",
     +             "11","12","13","14","15","16","17","18","19","20",
     +             "21","22","23","24","25","26","27","28","29","30",
     +             "31","32"/
      Data Ione,Fact,Dnull /1,0.66666666666666667D0,0.0D0/

      Call Getrec(20,"JOBARC","NBASTOT ",Ione,Nao)
      Nao2 = Nao*Nao

      Ntotal = 0
      Do Irrep = 1, Nirrep
         Ntotal = Ntotal + Nroot(Irrep)
      Enddo

      Itmp1 = Ione
      Itmpx = Itmp1 + Nao2
      Itmpy = Itmpx + Nao2
      Itmpz = Itmpy + Nao2
      Igt   = Itmpz + Nao2
      Iend  = Igt   + Ntotal*Nao2
      If (Iend .Gt. Maxcor) Call Insmem("form_genrlizd_oscstrength",
     +                                   Iend,Maxcor)
C
C Dipole integrals in the full AO (computational) basis
C
      Do Ixyz = 1, 3
         Ioff = Itmpx + (Ixyz-1)*Nao2
         If (Vprop) Then
            Length = Nao*(Nao+1)/2
            Call Getrec(20,"JOBARC",Label_xyz(Ixyz),Length*Iintfp,
     +                  Work(Itmp1))
            Call Expnd2(Work(Itmp1),Work(Ioff),Nao)
         Else
            Irrprt = Idipsm(Ixyz)
            Length = Lent(Irrprt)
            Call Getrec(20,"JOBARC",Label_xyz(Ixyz),Length*Iintfp,
     +                  Work(Itmp1))
            Call Zero(Work(Ioff),Nao2)
            Call Matexp(Irrprt,Iaopop,Work(Itmp1),Work(Ioff))
            Call Matexp2(Irrprt,Work(Ioff),Work(Itmp1),Nao)
            Call Dcopy(Nao2,Work(Itmp1),1,Work(Ioff),1)
         Endif
      Enddo

      Call Dzero(Tm_x,1024)
      Call Dzero(Tm_y,1024)
      Call Dzero(Tm_z,1024)

      Idone = 0
      Do Irrepx = 1, Nirrep
      Do Iroot  = 1, Nroot(Irrepx)
         Idone = Idone + 1
         Label = "GTDENS"//Suffix(Idone)
         Call Getrec(20,"JOBARC",Label,Ntotal*Nao2*Iintfp,Work(Igt))
         Jdone = 0
         Do Irrepy = 1, Nirrep
         Do Jroot  = 1, Nroot(Irrepy)
            Jdone = Jdone + 1
            Ioff  = Igt + (Jdone-1)*Nao2
            Tm_x(Jroot,Irrepy,Iroot,Irrepx) =
     +                Ddot(Nao2,Work(Ioff),1,Work(Itmpx),1)
            Tm_y(Jroot,Irrepy,Iroot,Irrepx) =
     +                Ddot(Nao2,Work(Ioff),1,Work(Itmpy),1)
            Tm_z(Jroot,Irrepy,Iroot,Irrepx) =
     +                Ddot(Nao2,Work(Ioff),1,Work(Itmpz),1)
         Enddo
         Enddo
      Enddo
      Enddo

      Write(6,*)
      Write(6,"(2x,2a)") "State-to-state transition moments ",
     +                   "<L_k|mu|R_l> (a.u., electronic only)"
      Write(6,"(2x,a,6x,a,9x,a,12x,a,12x,a)") "k: irrep root",
     +      "l: irrep root","x","y","z"
      Do Irrepx = 1, Nirrep
      Do Iroot  = 1, Nroot(Irrepx)
      Do Irrepy = 1, Nirrep
      Do Jroot  = 1, Nroot(Irrepy)
         Write(6,"(5x,i2,3x,i2,10x,i2,3x,i2,3(1x,F12.6))")
     +            Irrepy,Jroot,Irrepx,Iroot,
     +            Tm_x(Jroot,Irrepy,Iroot,Irrepx),
     +            Tm_y(Jroot,Irrepy,Iroot,Irrepx),
     +            Tm_z(Jroot,Irrepy,Iroot,Irrepx)
      Enddo
      Enddo
      Enddo
      Enddo

      Call Dzero(Osc,1024)
      Length = 0
      Do Irrepx = 1, Nirrep
      Do Iroot  = 1, Nroot(Irrepx)
         E_i = Eigval(Iroot,Irrepx)
         Do Irrepy = 1, Nirrep
         Do Jroot  = 1, Nroot(Irrepy)
            E_j = Eigval(Jroot,Irrepy)
            If (Irrepx .Eq. Irrepy .And. Iroot .Eq. Jroot) Then
               E = Dnull
            Else
               E = Dabs(E_i - E_j)
            Endif
            Length = Length + 1
            Osc(Length) = Fact*E*
     +       (Tm_x(Jroot,Irrepy,Iroot,Irrepx)*
     +        Tm_x(Iroot,Irrepx,Jroot,Irrepy) +
     +        Tm_y(Jroot,Irrepy,Iroot,Irrepx)*
     +        Tm_y(Iroot,Irrepx,Jroot,Irrepy) +
     +        Tm_z(Jroot,Irrepy,Iroot,Irrepx)*
     +        Tm_z(Iroot,Irrepx,Jroot,Irrepy))
         Enddo
         Enddo
      Enddo
      Enddo

      Call Putrec(20,"JOBARC","GEN_OSCS",Length*Iintfp,Osc)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(20x,2a)") "root(irrepy),irrepy,root(irrepx)",
     +                     " irrepx", " Osc. strength"
      Write(6,"(20x,a)") "---------------------------------------------"
      Ndone = 0
      Do Irrepx = 1, Nirrep
      Do Iroot  = 1, Nroot(Irrepx)
      Do Irrepy = 1, Nirrep
      Do Jroot  = 1, Nroot(Irrepy)
         Ndone = Ndone + 1
         write(6,"(a,(4(7x,i2)),f10.5)") "Osc. Strength: ",
     +            Jroot,Irrepy,Iroot,Irrepx,Osc(Ndone)
      Enddo
      Enddo
      Enddo
      Enddo
      Write(6,"(20x,a)") "---------------------------------------------"
#endif
      Return
      End
