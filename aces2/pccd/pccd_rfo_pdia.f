










      Subroutine Pccd_rfo_pdia(P,C,F,Scr,Wk,Lwork,Irp,Nbas,Nocc,Npar)
C
C Diagonal preconditioner for the Davidson iterations of the RFO
C stability check (Pccd_rfo_stab), at the current (stationary) point:
C   ov pairs:     4|f(q,q)-f(p,p)|                 (as the model H0)
C   oo, vv pairs: max(2|D(p)-D(q)| |f(q,q)-f(p,p)|, Pfloor)
C with D(p) = 2(n_ref(p) + D_corr(p,p)) from DENSOO/DENSVV (blocks per
C irrep, MOs occupied by irrep then virtual by irrep, as Irp). Used only
C to precondition: as the RFO starting Hessian the oo/vv formula failed
C for stretched N2. With the exact FD Hessian at the N2 cc-pVDZ 1.1016 A
C symmetric saddle, the negative mode (-3.6e-4) is found after 92 Hv
C products, against 117 with the model H0 diagonal. RHF.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension P(Npar),C(Nbas,Nbas),F(Nbas,Nbas),Scr(Nbas,Nbas)
      Dimension Wk(Lwork),Irp(Nbas),Nob(8),Nvb(8),Ioff(8),Ivff(8)
      Dimension Nj(8)


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



      Parameter (Pfloor=1.0D-2)
      Data Zero,One,Two,Four /0.0D0,1.0D0,2.0D0,4.0D0/

      Nirrep = 0
      Do I = 1, Nbas
         Nirrep = Max(Nirrep,Irp(I))
      Enddo
      Do Irr = 1, 8
         Nob(Irr) = 0
         Nvb(Irr) = 0
      Enddo
      Do I = 1, Nbas
         If (I .Le. Nocc) Then
            Nob(Irp(I)) = Nob(Irp(I)) + 1
         Else
            Nvb(Irp(I)) = Nvb(Irp(I)) + 1
         Endif
      Enddo
      Noo = 0
      Nvv = 0
      Do Irr = 1, Nirrep
         Ioff(Irr) = Noo
         Ivff(Irr) = Nvv
         Noo = Noo + Nob(Irr)*Nob(Irr)
         Nvv = Nvv + Nvb(Irr)*Nvb(Irr)
      Enddo
C Wk: DENSOO, DENSVV, then the occupations D(p) and f(p,p).
      Idoo = 1
      Idvv = Idoo + Noo
      Iocc = Idvv + Nvv
      Ifpp = Iocc + Nbas
      If (Ifpp+Nbas-1 .Gt. Lwork) Call Insmem("Pccd_rfo_pdia",
     +                                        Ifpp+Nbas-1,Lwork)
      Call Getrec(20,"JOBARC","DENSOO  ",Noo*Iintfp,Wk(Idoo))
      Call Getrec(20,"JOBARC","DENSVV  ",Nvv*Iintfp,Wk(Idvv))
C D(p): the J-th orbital of its irrep block (dimension Nob or Nvb) has
C its diagonal element at block offset + (J-1)*dimension + J.
      Do Irr = 1, 8
         Nj(Irr) = 0
      Enddo
      Do I = 1, Nocc
         Irr = Irp(I)
         Nj(Irr) = Nj(Irr) + 1
         Ii = Ioff(Irr) + (Nj(Irr)-1)*Nob(Irr) + Nj(Irr)
         Wk(Iocc+I-1) = Two*(One + Wk(Idoo+Ii-1))
      Enddo
      Do Irr = 1, 8
         Nj(Irr) = 0
      Enddo
      Do I = Nocc+1, Nbas
         Irr = Irp(I)
         Nj(Irr) = Nj(Irr) + 1
         Ia = Ivff(Irr) + (Nj(Irr)-1)*Nvb(Irr) + Nj(Irr)
         Wk(Iocc+I-1) = Two*Wk(Idvv+Ia-1)
      Enddo

      Call Getrec(20,"JOBARC","FOCKA   ",Nbas*Nbas*Iintfp,F)
      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,F,Nbas,C,Nbas,Zero,
     +           Scr,Nbas)
      Do I = 1, Nbas
         Wk(Ifpp+I-1) = Ddot(Nbas,C(1,I),1,Scr(1,I),1)
      Enddo

      K = 0
      Do Ip = 1, Nbas
         Do Iq = Ip+1, Nbas
            If (Irp(Ip) .Eq. Irp(Iq)) Then
               K = K + 1
               Df = Dabs(Wk(Ifpp+Iq-1) - Wk(Ifpp+Ip-1))
               If (Ip .Le. Nocc .And. Iq .Gt. Nocc) Then
                  P(K) = Four*Df
               Else
                  P(K) = Max(Two*Dabs(Wk(Iocc+Ip-1)-Wk(Iocc+Iq-1))*Df,
     +                       Pfloor)
               Endif
            Endif
         Enddo
      Enddo
      If (K .Ne. Npar) Then
         Write(6,"(a,2i6)") " @-Pccd_rfo_pdia: pair count mismatch",
     +                      K,Npar
         Call Errex
      Endif
      Return
      End
