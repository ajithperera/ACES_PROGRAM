      Subroutine Form_genrlizd_osc_bd(Work,Maxcor,Iuhf,Eigval)
C
C State-to-state oscillator strengths beyond the dipole
C approximation from the generalized transition densities
C ("GTDENSll", spin-summed AO <L_k|{p+ q}|R_l>, see
C form_genrlizd_tdens). For every pair the multipole transition
C moments of both directions, <L_k|O|R_l> and <L_l|O|R_k>, are
C combined by Genrlizd_multipole_init/_fin (unchanged).
C
C 11/2026 (Ajith Perera): rewritten for the per-R-state AO density
C records. The 8-index pair-of-pairs arrays of osc_bd_init/osc_bd_fin
C (about 416 MB, with mutually inconsistent declared shapes) are
C replaced by gathering the two needed moment vectors per pair.
C
      Implicit Double Precision(A-H,O-Z)

      Double Precision Mq_tl,Mq_tr,Mq2_t,M_iso
      Character*8 Label
      Character*2 Suffix(32)

      Dimension Eigval(100,8)
      Dimension Work(Maxcor)
      Dimension Qm_tl(6),Om_tl(10),Am_tl(3),Dm_tl(3),Mq_tl(9)
      Dimension Qm_tr(6),Om_tr(10),Am_tr(3),Dm_tr(3),Mq_tr(9)
      Dimension Qm2_t(9),Om2_t(10),Am2_t(3),Dm2_t(3),Mq2_t(27)

#include "flags.com"
#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"
#include "machsp.com"

      Common /Calcinfo/Nroot(8)

      Data Suffix /"01","02","03","04","05","06","07","08","09","10",
     +             "11","12","13","14","15","16","17","18","19","20",
     +             "21","22","23","24","25","26","27","28","29","30",
     +             "31","32"/
      Data Ione,Dnull /1,0.0D0/
C
C Moment arrays M(Jroot,Irrepy,Iroot,Irrepx,comp) = <L_J|O_comp|R_I>
C
      N1 = 1024
      Call Getrec(20,"JOBARC","NBASTOT ",Ione,Nao)
      Nao2 = Nao*Nao

      Ntotal = 0
      Do Irrep = 1, Nirrep
         Ntotal = Ntotal + Nroot(Irrep)
      Enddo

      Iqm   = Ione
      Iom   = Iqm   + N1*6
      Iam   = Iom   + N1*10
      Idm   = Iam   + N1*3
      Imq   = Idm   + N1*3
      Ioscq = Imq   + N1*9
      Iosco = Ioscq + N1
      Iosca = Iosco + N1
      Ioscd = Iosca + N1
      Ioscm = Ioscd + N1
      Iscr1 = Ioscm + N1
      Iscr2 = Iscr1 + Nao2
      Igt   = Iscr2 + Nao2
      Iend  = Igt   + Ntotal*Nao2
      If (Iend .Ge. Maxcor) Call Insmem("form_genrlizd_osc_bd",
     +                                   Iend,Maxcor)
      Mxcor = Maxcor - Iend + 1

      Call Dzero(Work(Iqm),Iscr1-Iqm)

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
            Call Genrlizd_beyond_dipole(Work(Ioff),Work(Iscr1),
     +                                  Work(Iscr2),Work(Iqm),
     +                                  Work(Iom),Work(Iam),
     +                                  Work(Idm),Work(Imq),
     +                                  Work(Iend),Mxcor,Nao,
     +                                  Jroot,Irrepy,Iroot,Irrepx)
         Enddo
         Enddo
      Enddo
      Enddo

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
C
C "r" = <L_J|O|R_I>, "l" = <L_I|O|R_J>
C
            Ir = Jroot + 4*(Irrepy-1) + 32*(Iroot-1) + 128*(Irrepx-1)
            Il = Iroot + 4*(Irrepx-1) + 32*(Jroot-1) + 128*(Irrepy-1)
            Do Ic = 1, 6
               Qm_tr(Ic) = Work(Iqm-1+Ir+N1*(Ic-1))
               Qm_tl(Ic) = Work(Iqm-1+Il+N1*(Ic-1))
            Enddo
            Do Ic = 1, 10
               Om_tr(Ic) = Work(Iom-1+Ir+N1*(Ic-1))
               Om_tl(Ic) = Work(Iom-1+Il+N1*(Ic-1))
            Enddo
            Do Ic = 1, 3
               Am_tr(Ic) = Work(Iam-1+Ir+N1*(Ic-1))
               Am_tl(Ic) = Work(Iam-1+Il+N1*(Ic-1))
               Dm_tr(Ic) = Work(Idm-1+Ir+N1*(Ic-1))
               Dm_tl(Ic) = Work(Idm-1+Il+N1*(Ic-1))
            Enddo
            Do Ic = 1, 9
               Mq_tr(Ic) = Work(Imq-1+Ir+N1*(Ic-1))
               Mq_tl(Ic) = Work(Imq-1+Il+N1*(Ic-1))
            Enddo

            Call Genrlizd_multipole_init(Qm_tl,Om_tl,Am_tl,Mq_tl,Dm_tl,
     +                                   Qm_tr,Om_tr,Am_tr,Mq_tr,Dm_tr,
     +                                   Qm2_t,Om2_t,Am2_t,Mq2_t,Dm2_t)
            Call Genrlizd_multipole_fin(Qm2_t,Om2_t,Am2_t,Dm2_t,Mq2_t,
     +                                  E,D_iso,A_iso,O_iso,Q_iso,
     +                                  M_iso)
            Length = Length + 1
            Work(Ioscq-1+Length) = Q_iso
            Work(Iosco-1+Length) = O_iso
            Work(Iosca-1+Length) = A_iso
            Work(Ioscd-1+Length) = D_iso
            Work(Ioscm-1+Length) = M_iso
         Enddo
         Enddo
      Enddo
      Enddo

      Call Putrec(20,"JOBARC","QUADMOMC",Length*Iintfp,Work(Ioscq))
      Call Putrec(20,"JOBARC","OCTPMOMC",Length*Iintfp,Work(Iosco))
      Call Putrec(20,"JOBARC","ANGLMOMC",Length*Iintfp,Work(Iosca))
      Call Putrec(20,"JOBARC","MAGQUADC",Length*Iintfp,Work(Ioscm))
      Call Putrec(20,"JOBARC","DIPOLEMC",Length*Iintfp,Work(Ioscd))

      Return
      End
