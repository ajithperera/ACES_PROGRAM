










      Subroutine Pccd_ohess_occ(Occ,Dd,A,S,Irp,Nbas,Nocc)
C
C pCCD occupations n(p) = a(p) + 2 b(p), a = 2 (occupied) or 0, b the
C diagonal of DENSOO/DENSVV (blocks per irrep, MOs numbered as Irp),
C and the parts of A(p,q) = Gam(pp|qq) and S(p,q) = Gam(pq|qp) +
C Gam(pq|pq) from the reference and the reference-correlation cross
C terms, Gam_ref(a) + 4 gam2(a,b) (see Pccd_ohess):
C   A(p,q) = a(p)a(q) + 2(a(p)b(q) + b(p)a(q))
C          - d(p,q)(a(p)a(p)/2 + 2a(p)b(p))
C   S(p,q) = -a(p)a(q)/2 - (a(p)b(q) + b(p)a(q)),  p.ne.q
C Dd: scratch for the DENSOO and DENSVV records. RHF.
C
      Implicit Double Precision(A-H,O-Z)
      Dimension Occ(Nbas),Dd(*),A(Nbas,Nbas),S(Nbas,Nbas),Irp(Nbas)
      Dimension Nob(8),Nvb(8),Ioff(8),Ivff(8),Nj(8)


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



      Data Zero,Half,Two /0.0D0,0.5D0,2.0D0/

      Nirr = 0
      Do I = 1, Nbas
         Nirr = Max(Nirr,Irp(I))
      Enddo
      Do Irr = 1, 8
         Nob(Irr) = 0
         Nvb(Irr) = 0
         Nj(Irr)  = 0
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
      Do Irr = 1, Nirr
         Ioff(Irr) = Noo
         Ivff(Irr) = Nvv
         Noo = Noo + Nob(Irr)*Nob(Irr)
         Nvv = Nvv + Nvb(Irr)*Nvb(Irr)
      Enddo
C b(p) goes to S(p,1) for the moment; a(p) to A(p,1).
      Call Getrec(20,"JOBARC","DENSOO  ",Noo*Iintfp,Dd)
      Do I = 1, Nocc
         Irr = Irp(I)
         Nj(Irr) = Nj(Irr) + 1
         S(I,1) = Dd(Ioff(Irr) + (Nj(Irr)-1)*Nob(Irr) + Nj(Irr))
         A(I,1) = Two
      Enddo
      Call Getrec(20,"JOBARC","DENSVV  ",Nvv*Iintfp,Dd)
      Do Irr = 1, 8
         Nj(Irr) = 0
      Enddo
      Do I = Nocc+1, Nbas
         Irr = Irp(I)
         Nj(Irr) = Nj(Irr) + 1
         S(I,1) = Dd(Ivff(Irr) + (Nj(Irr)-1)*Nvb(Irr) + Nj(Irr))
         A(I,1) = Zero
      Enddo
      Do I = 1, Nbas
         Occ(I)   = A(I,1) + Two*S(I,1)
         Dd(I)        = A(I,1)
         Dd(Nbas+I)   = S(I,1)
      Enddo
      Do Iq = 1, Nbas
         Aq = Dd(Iq)
         Bq = Dd(Nbas+Iq)
         Do Ip = 1, Nbas
            Ap = Dd(Ip)
            Bp = Dd(Nbas+Ip)
            A(Ip,Iq) = Ap*Aq + Two*(Ap*Bq + Bp*Aq)
            S(Ip,Iq) = -Half*Ap*Aq - (Ap*Bq + Bp*Aq)
         Enddo
         A(Iq,Iq) = A(Iq,Iq) - Half*Aq*Aq - Two*Aq*Bq
         S(Iq,Iq) = Zero
      Enddo
      Return
      End
