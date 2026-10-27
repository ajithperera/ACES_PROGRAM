      Subroutine Erg_gam(Mode,Ipx,N,Nirr,Nbg,Irg,Npr,Ioffg,Lgso,Gso,
     &                   Buf)
C
C The nonseparable two-particle density of a correlated gradient in
C the SO basis, as xbcktrn leaves it on the MOINTS lists (the lists
C xvdint's Gamsrt reads), stored for random access:
C
C    Gam(pq,rs) = Gso(Ioffg(g) + (Ipx(p,q)-1)*Npr(g) + Ipx(r,s))
C
C with g the irrep of the pair pq (= that of rs, other elements are
C zero), Ipx(p,q) = Ipx(q,p) the position of the pair among the
C Npr(g) pairs p >= q of irrep g. Both (pq,rs) and (rs,pq) are kept.
C
C The SO basis of the lists has N functions, Nbg(i) in irrep i,
C numbered irrep by irrep (Irg(p) the irrep of p): the SOs of the SCF,
C or for spherical functions the Cartesian SOs xbcktrn works in.
C
C Mode 1: set up Ipx, Npr, Ioffg and the length Lgso of Gso.
C Mode 2: read the lists into Gso (Buf: N*N scratch).
C
C The lists, in Gamsrt's order (Getlst(Buf,column,1,1,type,ind)):
C  type 1, ind = irrep i      : columns p >= q of irrep i, each a
C                               packed triangle r >= s of irrep i
C  type 2, irreps a < b       : columns p >= q of a, triangles of b
C  type 3, irreps a < b       : columns p (a) x q (b), full r (a) x
C                               s (b) blocks, s fastest
C  type 4, irreps (a1..a4)    : columns p (a1) x q (a2), full r (a3)
C                               x s (a4) blocks, s fastest
C
      Implicit Double Precision (A-H,O-Z)

      Integer Ipx(N,N), Nbg(8), Irg(N), Npr(8), Ioffg(8), Iof(8),
     &        Idid(8), Itp(4)
      Dimension Gso(*), Buf(*)

      Iof(1) = 0
      Do I = 2, Nirr
         Iof(I) = Iof(I-1) + Nbg(I-1)
      Enddo

      If (Mode .Eq. 1) Then
         Do I = 1, 8
            Npr(I) = 0
         Enddo
         Do Ip = 1, N
            Do Iq = 1, Ip
               Ig = Ieor(Irg(Ip)-1,Irg(Iq)-1) + 1
               Npr(Ig) = Npr(Ig) + 1
               Ipx(Ip,Iq) = Npr(Ig)
               Ipx(Iq,Ip) = Npr(Ig)
            Enddo
         Enddo
         Lgso = 0
         Do Ig = 1, Nirr
            Ioffg(Ig) = Lgso
            Lgso = Lgso + Npr(Ig)*Npr(Ig)
         Enddo
         Return
      Endif

      Call Zero(Gso,Lgso)

C Type 1.

      Do Ia = 1, Nirr
         Na = Nbg(Ia)
         Kc = 0
         Do L1 = 1, Na
         Do L2 = 1, L1
            Kc = Kc + 1
            Call Getlst(Buf,Kc,1,1,1,Ia)
            Ind = 0
            Do L3 = 1, Na
            Do L4 = 1, L3
               Ind = Ind + 1
               Call Erg_gset(Iof(Ia)+L1,Iof(Ia)+L2,Iof(Ia)+L3,
     &                       Iof(Ia)+L4,Buf(Ind),Ipx,N,Irg,Npr,Ioffg,
     &                       Gso)
            Enddo
            Enddo
         Enddo
         Enddo
      Enddo
      If (Nirr .Eq. 1) Return

C Type 2.

      Ind2 = 0
      Do Ib = 2, Nirr
      Do Ia = 1, Ib-1
         Ind2 = Ind2 + 1
         Na = Nbg(Ia)
         Nb = Nbg(Ib)
         If (Na .Eq. 0 .Or. Nb .Eq. 0) Goto 200
         Kc = 0
         Do L1 = 1, Na
         Do L2 = 1, L1
            Kc = Kc + 1
            Call Getlst(Buf,Kc,1,1,2,Ind2)
            Ind = 0
            Do L3 = 1, Nb
            Do L4 = 1, L3
               Ind = Ind + 1
               Call Erg_gset(Iof(Ia)+L1,Iof(Ia)+L2,Iof(Ib)+L3,
     &                       Iof(Ib)+L4,Buf(Ind),Ipx,N,Irg,Npr,Ioffg,
     &                       Gso)
            Enddo
            Enddo
         Enddo
         Enddo
  200    Continue
      Enddo
      Enddo

C Type 3.

      Ind3 = 0
      Do Ir = 2, Nirr
      Do Ia = 1, Nirr
         Ib = Ieor(Ir-1,Ia-1) + 1
         If (Ib .Le. Ia) Goto 300
         Ind3 = Ind3 + 1
         Na = Nbg(Ia)
         Nb = Nbg(Ib)
         If (Na .Eq. 0 .Or. Nb .Eq. 0) Goto 300
         Kc = 0
         Do L1 = 1, Na
         Do L2 = 1, Nb
            Kc = Kc + 1
            Call Getlst(Buf,Kc,1,1,3,Ind3)
            Do L3 = 1, Na
            Do L4 = 1, Nb
               Call Erg_gset(Iof(Ia)+L1,Iof(Ib)+L2,Iof(Ia)+L3,
     &                       Iof(Ib)+L4,Buf(L4+(L3-1)*Nb),
     &                       Ipx,N,Irg,Npr,Ioffg,Gso)
            Enddo
            Enddo
         Enddo
         Enddo
  300    Continue
      Enddo
      Enddo
      If (Nirr .Lt. 4) Return

C Type 4.

      Ind4 = 0
      Do Ir = 2, Nirr
      Do I1 = 1, Nirr
         I2 = Ieor(Ir-1,I1-1) + 1
         If (I2 .Lt. I1) Goto 400
         Ibot = Max(I1,I2) + 1
         Do I = 1, 8
            Idid(I) = 0
         Enddo
         Do Itmp = Ibot, Nirr
            I4 = Ieor(Itmp-1,Ir-1) + 1
            I3 = Min(Itmp,I4)
            I4 = Max(Itmp,I4)
            If (Max(Idid(I4),Idid(I3)) .Ne. 0) Goto 410
            Idid(I3) = 1
            Idid(I4) = 1
            Ind4 = Ind4 + 1
            Itp(1) = I3
            Itp(2) = I4
            Itp(3) = I1
            Itp(4) = I2
            N1 = Nbg(Itp(1))
            N2 = Nbg(Itp(2))
            N3 = Nbg(Itp(3))
            N4 = Nbg(Itp(4))
            If (N1*N2*N3*N4 .Eq. 0) Goto 410
            Kc = 0
            Do L1 = 1, N1
            Do L2 = 1, N2
               Kc = Kc + 1
               Call Getlst(Buf,Kc,1,1,4,Ind4)
               Do L3 = 1, N3
               Do L4 = 1, N4
                  Call Erg_gset(Iof(Itp(1))+L1,Iof(Itp(2))+L2,
     &                          Iof(Itp(3))+L3,Iof(Itp(4))+L4,
     &                          Buf(L4+(L3-1)*N4),
     &                          Ipx,N,Irg,Npr,Ioffg,Gso)
               Enddo
               Enddo
            Enddo
            Enddo
  410       Continue
         Enddo
  400    Continue
      Enddo
      Enddo

      Return
      End
