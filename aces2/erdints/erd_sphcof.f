










      Subroutine Erd_sphcof(L,C,Ldc,Nsph)
C
C Real spherical functions of an L shell, as ERD makes them, written
C in VMOL's Cartesian AOs (ERD order x^p y^q z^r, p>q>r):
C
C    Y_m = sum_k C(k,m) AO_k,   m = 1..Nsph = 2L+1 (m = +L first)
C
C ERD's Y_m = N_m sum_k T(k,m) x^p y^q z^r R(r), T from
C ERD__XYZ_TO_RY_MATRIX and N_m normalizing it. ERD normalizes every
C Cartesian component separately, and the ACES AO is fac_k times that
C (Erd_cmpfac), which makes all ACES AOs x^p y^q z^r R(r)/sqrt(K) with
C one common K. Over one exponent <x^a y^b z^c|x^a' y^b' z^c'> is
C K I(a+a') I(b+b') I(c+c'), I(n) = (n-1)!! for even n and 0 for odd,
C so C(k,m) = T(k,m)/sqrt(Q_m), Q_m = sum T(k,m) T(k',m) I I I.
C (Checked against OED by Erd_sphchk.) For L <= 1 the functions are
C the Cartesian ones.
C
      Implicit Double Precision (A-H,O-Z)
      Parameter (Lmax = 7, Mxyz = (Lmax+1)*(Lmax+2)/2,
     &           Mrow = (Lmax/2+1)*(Lmax/2+2)/2)
      Dimension C(Ldc,*), Temp(Mxyz), Tmat(Mrow,2*Lmax+1)
      Integer Ixyz(3,Mxyz), Nrow(2*Lmax+1), Irow(Mrow,2*Lmax+1)
      Double Precision Ovmono

      Ncrt = (L+1)*(L+2)/2
      Call Erd_monos(L,Ixyz,N)

      If (L .Le. 1) Then
         Nsph = Ncrt
         Do M = 1, Nsph
            Do K = 1, Ncrt
               C(K,M) = 0.0D0
            Enddo
            C(M,M) = 1.0D0
         Enddo
         Return
      Endif

      Nsph = 2*L + 1
      Call Erd__xyz_to_ry_matrix(Ncrt,Nsph,Mrow,L,Temp,Nrow,Irow,Tmat)
      Do M = 1, Nsph
         Do K = 1, Ncrt
            C(K,M) = 0.0D0
         Enddo
         Do J = 1, Nrow(M)
            C(Irow(J,M),M) = Tmat(J,M)
         Enddo
         Qm = 0.0D0
         Do K = 1, Ncrt
            Do K2 = 1, Ncrt
               Qm = Qm + C(K,M)*C(K2,M)*Ovmono(Ixyz(1,K),Ixyz(1,K2))
            Enddo
         Enddo
         Scl = 1.0D0/Dsqrt(Qm)
         Do K = 1, Ncrt
            C(K,M) = C(K,M)*Scl
         Enddo
      Enddo

      Return
      End
