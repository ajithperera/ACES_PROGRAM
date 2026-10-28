      Subroutine Erd_wrja(Icmpao,Ncrt,U,Ui,Csy,Scr,Ipiv,Lwork)
C
C JOBARC records of the SO basis that vmol2ja (driver.F, v2ja.F)
C writes from VMOL's SYMTRANS and CSYMTRAN: NIRREP, NUMBASIR, NBASTOT,
C NAOBASFN, FAOBASIR, AO2SO (SO -> computational AO, rows are the
C Cartesian AOs; the name is the reverse of what it does), FULLSOAO,
C FULLAOSO, AO2SOINV, the same pair with the AOs in ZMAT order
C (CMP2ZMAT, ZMAT2CMP), and for spherical bases CART3CMP, CART2CMP,
C CMP3CART, CMP2CART.
C
C U = FULLSOAO is square: the Nso SOs, then (spherical basis only) the
C symmetry adapted contaminants r^2j Y_(L-2j) that complete the
C Cartesian space, all made of VMOL's functions (Erd_vharm). Csy (= VMOL's CSYMTRAN) is the purely Cartesian SO
C basis. As in vmol2ja the inverses are taken with LAPACK.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"
#include "machsp.com"

      Parameter (Lmax = 7, Mxyz = (Lmax+1)*(Lmax+2)/2)
      Integer Icmpao(Ncrt), Ipiv(Ncrt), Nfao(8), Irc(8), Isg(8,8),
     &        Nconir(8),
     &        Iparv(3,Mxyz), Ixyz(3,Mxyz), Ipar(3)
      Dimension U(Ncrt,Ncrt), Ui(Ncrt,Ncrt), Csy(Ncrt,Ncrt),
     &          Scr(Lwork), Cv(Mxyz,Mxyz), W(13,13)

      Call Putrec(20,'JOBARC','NIRREP  ',1,Nirr)
      Call Putrec(20,'JOBARC','NUMBASIR',Nirr,Nbfirr)
      Call Putrec(20,'JOBARC','NBASTOT ',1,Nso)
      Call Putrec(20,'JOBARC','NAOBASFN',1,Ncrt)

C SOs -> AOs.

      Call Dzero(U,Ncrt*Ncrt)
      Do Ip = 1, Nso
         Io  = Isoorb(Ip)
         Is0 = Iorbsh(1,Io)
         L   = Lshl(Is0)
         Nk  = (L+1)*(L+2)/2
         Nc  = Nfnsh(Is0)/Nctsh(Is0)
         It  = (Isofn(Ip)-1)/Nc + 1
         Mc  = Mod(Isofn(Ip)-1,Nc) + 1
         If (Spher .And. L .Ge. 2)
     &      Call Erd_vsph(L,Cv,Mxyz,W,13,Iparv)
         Do M = 1, Norbm(Io)
            I0 = Iofcrt(Iorbsh(M,Io)) + (It-1)*Nk
            If (Spher .And. L .Ge. 2) Then
               Do K = 1, Nk
                  U(Icmpao(I0+K),Ip) = Dble(Isosgn(M,Ip))*Cv(K,Mc)
               Enddo
            Else
               U(Icmpao(I0+Mc),Ip) = Dble(Isosgn(M,Ip))
            Endif
         Enddo
      Enddo

C Contaminants (irrep by irrep) and the Cartesian SO basis.

      Do Ig = 1, Nirr
         Nfao(Ig) = 0
      Enddo
      If (Spher) Then
         Ncol = Nso
         Do Ig = 1, Nirr
            Nconir(Ig) = 0
         Enddo
         Do Ig = 1, Nirr
            Do Io = 1, Norb
               Is0 = Iorbsh(1,Io)
               L   = Lshl(Is0)
               Nk  = (L+1)*(L+2)/2
               If (L .Ge. 2) Then
                  Call Erd_vsph(L,Cv,Mxyz,W,13,Iparv)
                  Ns = 2*L + 1
                  Do It = 1, Nctsh(Is0)
                     Do Ic = Ns+1, Nk
                        Call Erd_irrfn(Io,Iparv(1,Ic),Irc,Isg)
                        Do J = 1, Norbm(Io)
                           If (Irc(J) .Eq. Ig) Then
                              Ncol = Ncol + 1
                              Nconir(Ig) = Nconir(Ig) + 1
C VMOL's phase: the member that comes first in the computational
C order gets the + sign (as in Erd_soord).
                              Kmn = Ncrt + 1
                              Do M = 1, Norbm(Io)
                                 I0 = Iofcrt(Iorbsh(M,Io))
     &                              + (It-1)*Nk
                                 If (Icmpao(I0+1) .Lt. Kmn) Then
                                    Kmn = Icmpao(I0+1)
                                    Isg0 = Isg(M,J)
                                 Endif
                              Enddo
                              Do M = 1, Norbm(Io)
                                 I0 = Iofcrt(Iorbsh(M,Io))
     &                              + (It-1)*Nk
                                 Do K = 1, Nk
                                    U(Icmpao(I0+K),Ncol) =
     &                              Dble(Isg(M,J)*Isg0)*Cv(K,Ic)
                                 Enddo
                              Enddo
                           Endif
                        Enddo
                     Enddo
                  Enddo
               Endif
            Enddo
         Enddo
         If (Ncol .Ne. Ncrt) Then
            Write(6,"(a,2i6)") " @Erd_wrja: contaminant count error ",
     &                         Ncol, Ncrt
            Call Errex
         Endif

C VMOL lists the contaminants of an irrep by their first AO.

         Call Erd_consrt(U,Ncrt,Nso,Nconir,Nirr,Ipiv,Scr)

         Call Dzero(Csy,Ncrt*Ncrt)
         Ncol = 0
         Do Ig = 1, Nirr
            Do Io = 1, Norb
               Is0 = Iorbsh(1,Io)
               L   = Lshl(Is0)
               Call Erd_monos(L,Ixyz,Nk)
               Do It = 1, Nctsh(Is0)
                  Do K = 1, Nk
                     Do Kk = 1, 3
                        Ipar(Kk) = Mod(Ixyz(Kk,K),2)
                     Enddo
                     Call Erd_irrfn(Io,Ipar,Irc,Isg)
                     Do J = 1, Norbm(Io)
                        If (Irc(J) .Eq. Ig) Then
                           Ncol = Ncol + 1
                           Nfao(Ig) = Nfao(Ig) + 1
                           Do M = 1, Norbm(Io)
                              I0 = Iofcrt(Iorbsh(M,Io)) + (It-1)*Nk
                              Csy(Icmpao(I0+K),Ncol) =
     &                           Dble(Isg(M,J))
                           Enddo
                        Endif
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo
      Endif
      If (Spher) Call Erd_csyord(Csy,Ncrt,Nfao,Nirr,Ipiv,Scr)
      Call Putrec(20,'JOBARC','FAOBASIR',Nirr,Nfao)

      Call Putrec(20,'JOBARC','AO2SO   ',Ncrt*Nso*Iintfp,U)
      Call Putrec(20,'JOBARC','FULLSOAO',Ncrt*Ncrt*Iintfp,U)
      Call Erd_minv(U,Ui,Ncrt,Ipiv,Scr,Lwork)
      Call Putrec(20,'JOBARC','FULLAOSO',Ncrt*Ncrt*Iintfp,Ui)
      Call Erd_blkcp(Ui,Ncrt,Nso,Ncrt,Scr)
      Call Putrec(20,'JOBARC','AO2SOINV',Nso*Ncrt*Iintfp,Scr)

C Cartesian SOs <-> computational SOs (spherical basis only).

      If (Spher) Then
         Call Xgemm('N','N',Ncrt,Ncrt,Ncrt,1.0D0,Ui,Ncrt,Csy,Ncrt,
     &              0.0D0,Scr,Ncrt)
         Call Putrec(20,'JOBARC','CART3CMP',Ncrt*Ncrt*Iintfp,Scr)
         Call Dcopy(Ncrt*Ncrt,Scr,1,Csy,1)
         Call Erd_blkcp(Csy,Ncrt,Nso,Ncrt,Scr)
         Call Putrec(20,'JOBARC','CART2CMP',Nso*Ncrt*Iintfp,Scr)
         Call Dcopy(Ncrt*Ncrt,Csy,1,Scr,1)
         Call Erd_minv(Scr,Csy,Ncrt,Ipiv,Ui,Ncrt*Ncrt)
         Call Putrec(20,'JOBARC','CMP3CART',Ncrt*Ncrt*Iintfp,Csy)
         Call Putrec(20,'JOBARC','CMP2CART',Ncrt*Nso*Iintfp,Csy)
         Call Erd_minv(U,Ui,Ncrt,Ipiv,Scr,Lwork)
      Endif

C Same transformations with the AOs in ZMAT order.

      Do Ip = 1, Nso
         Do Iao = 1, Ncrt
            Scr(Iao + (Ip-1)*Ncrt) = U(Icmpao(Iao),Ip)
         Enddo
      Enddo
      Call Putrec(20,'JOBARC','CMP2ZMAT',Ncrt*Nso*Iintfp,Scr)
      Do Iao = 1, Ncrt
         Do Ip = 1, Nso
            Scr(Ip + (Iao-1)*Nso) = Ui(Ip,Icmpao(Iao))
         Enddo
      Enddo
      Call Putrec(20,'JOBARC','ZMAT2CMP',Nso*Ncrt*Iintfp,Scr)

      Return
      End
