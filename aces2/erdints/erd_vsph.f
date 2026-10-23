










      Subroutine Erd_vsph(L,Cv,Ldc,W,Ldw,Iparv)
C
C VMOL's spherical functions of an L shell (L >= 2) and their relation
C to ERD's:
C
C   Cv(k,j)  : VMOL function j over the Cartesian AOs (Erd_vharm),
C              reordered: j = 1..2L+1 the harmonics, then the
C              contaminants, each in VMOL's table order
C   W(m,j)   : VMOL harmonic j = sum_m W(m,j) ERD spherical function m
C   Iparv    : parity (see Erd_cmppar) of every VMOL function
C
C W solves Cv(:,j) = sum_m C(:,m) W(m,j) with C = ERD's functions in
C the same AOs (Erd_sphcof). Stops unless there are 2L+1 harmonic
C columns, every function has one parity, and the
C fit is exact (VMOL's and ERD's harmonics span the same space).
C (VMOL's table lists the two kinds interleaved; they are told apart
C by their Laplacian.)
C
      Implicit Double Precision (A-H,O-Z)
      Parameter (Lmax = 7, Mxyz = (Lmax+1)*(Lmax+2)/2)
      Dimension Cv(Ldc,*), W(Ldw,*), C(Mxyz,Mxyz), A(2*Lmax+1,
     &          2*Lmax+1), B(2*Lmax+1), Q(Mxyz), Ct(Mxyz,Mxyz)
      Integer Iparv(3,*), Ixyz(3,Mxyz), Ipiv(2*Lmax+1), Iord(Mxyz)

      Ncrt = (L+1)*(L+2)/2
      Nsph = 2*L + 1
      Call Erd_vharm(L,Cv,Ldc)
      Call Erd_sphcof(L,C,Mxyz,Nsph0)
      Call Erd_monos(L,Ixyz,N)

C VMOL's table interleaves harmonics and contaminants (for f, x r^2
C is column 6); the harmonics are the columns with zero Laplacian.
C Keep VMOL's table order within each kind, harmonics first.

      Nh = 0
      Nc = 0
      Do J = 1, Ncrt
         Call Erd_lapl(L,Cv(1,J),Q)
         Qn = 0.0D0
         Do K = 1, (L-1)*L/2
            Qn = Qn + Dabs(Q(K))
         Enddo
         If (Qn .Eq. 0.0D0) Then
            Nh = Nh + 1
            Iord(Nh) = J
         Else
            Nc = Nc + 1
            Iord(Nsph+Nc) = J
         Endif
      Enddo
      If (Nh .Ne. Nsph) Then
         Write(6,"(a,2i3)") " @Erd_vsph: harmonic count in VMOL table ",
     &                      L, Nh
         Call Errex
      Endif
      Do J = 1, Ncrt
         Do K = 1, Ncrt
            Ct(K,J) = Cv(K,Iord(J))
         Enddo
      Enddo
      Do J = 1, Ncrt
         Do K = 1, Ncrt
            Cv(K,J) = Ct(K,J)
         Enddo
      Enddo

C Parities.

      Do J = 1, Ncrt
         Kfirst = 0
         Do K = 1, Ncrt
            If (Cv(K,J) .Ne. 0.0D0 .And. Kfirst .Eq. 0) Kfirst = K
         Enddo
         Do Kk = 1, 3
            Iparv(Kk,J) = Mod(Ixyz(Kk,Kfirst),2)
         Enddo
         Do K = 1, Ncrt
            If (Cv(K,J) .Ne. 0.0D0) Then
               Do Kk = 1, 3
                  If (Mod(Ixyz(Kk,K),2) .Ne. Iparv(Kk,J)) Then
                     Write(6,"(a,2i3)") " @Erd_vsph: mixed parity ",
     &                                  L, J
                     Call Errex
                  Endif
               Enddo
            Endif
         Enddo
      Enddo

C Least squares fit of the VMOL harmonics by ERD's: normal equations
C (C^T C) W = C^T Cv, then the residual must vanish.

      Do M = 1, Nsph
         Do M2 = 1, Nsph
            A(M,M2) = 0.0D0
            Do K = 1, Ncrt
               A(M,M2) = A(M,M2) + C(K,M)*C(K,M2)
            Enddo
         Enddo
      Enddo
      Call Dgetrf(Nsph,Nsph,A,2*Lmax+1,Ipiv,Info)
      If (Info .Ne. 0) Call Errex
      Resmx = 0.0D0
      Do J = 1, Nsph
         Do M = 1, Nsph
            B(M) = 0.0D0
            Do K = 1, Ncrt
               B(M) = B(M) + C(K,M)*Cv(K,J)
            Enddo
         Enddo
         Call Dgetrs('N',Nsph,1,A,2*Lmax+1,Ipiv,B,2*Lmax+1,Info)
         Do M = 1, Nsph
            W(M,J) = B(M)
         Enddo
         Do K = 1, Ncrt
            X = Cv(K,J)
            Do M = 1, Nsph
               X = X - C(K,M)*W(M,J)
            Enddo
            Resmx = Max(Resmx,Dabs(X))
         Enddo
      Enddo
      If (Resmx .Gt. 1.0D-10) Then
         Write(6,"(a,i3,1p,e10.2)") " @Erd_vsph: VMOL harmonics off ",
     &                              L, Resmx
         Call Errex
      Endif

      Return
      End
