      SUBROUTINE THR_MOMENT(L,M,N,GA,V,NT,D)
C
C Raw (not traceless) third moments of a Gaussian product about the
C origin shifted by D, no charge factor, in the order
C xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz. Same algebra as the 3rd
C moments inside OPAB3 (opaa0.F) before the octupole transform.
C Ajith Perera, 11/2026 (for r2rmom.F).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FNU(20), DP(20), V(10), D(3)
C
      NT=10
      FNU(1)=OLAP( L,M,N,GA)
      FNU(2)=OLAP(L+1,M,N,GA)
      FNU(3)=OLAP(L,M+1,N,GA)
      FNU(4)=OLAP(L,M,N+1,GA)
      FNU(5)=OLAP(L+2,M,N,GA)
      FNU(6)=OLAP(L,M+2,N,GA)
      FNU(7)=OLAP(L,M,N+2,GA)
      FNU(8)=OLAP(L+1,M+1,N,GA)
      FNU(9)=OLAP(L+1,M,N+1,GA)
      FNU(10)=OLAP(L,M+1,N+1,GA)
      FNU(11)=OLAP(L+3,M,N,GA)
      FNU(12)=OLAP(L,M+3,N,GA)
      FNU(13)=OLAP(L,M,N+3,GA)
      FNU(14)=OLAP(L+2,M+1,N,GA)
      FNU(15)=OLAP(L+2,M,N+1,GA)
      FNU(16)=OLAP(L+1,M+2,N,GA)
      FNU(17)=OLAP(L,M+2,N+1,GA)
      FNU(18)=OLAP(L+1,M,N+2,GA)
      FNU(19)=OLAP(L,M+1,N+2,GA)
      FNU(20)=OLAP(L+1,M+1,N+1,GA)
      DP(2)=D(1)
      DP(3)=D(2)
      DP(4)=D(3)
      DP(5)=D(1)**2
      DP(6)=D(2)**2
      DP(7)=D(3)**2
      DP(8)=DP(2)*DP(3)
      DP(9)=DP(2)*DP(4)
      DP(10)=DP(3)*DP(4)
      DP(11)=DP(5)*DP(2)
      DP(12)=DP(6)*DP(3)
      DP(13)=DP(7)*DP(4)
      DP(14)=DP(5)*DP(3)
      DP(15)=DP(5)*DP(4)
      DP(16)=DP(6)*DP(2)
      DP(17)=DP(6)*DP(4)
      DP(18)=DP(7)*DP(2)
      DP(19)=DP(7)*DP(3)
      DP(20)=DP(8)*DP(4)
      V(1)=FNU(11)+3.*DP(2)*FNU(5)+3.*DP(5)*FNU(2)+DP(11)*FNU(1)
      V(2)=FNU(12)+3.*DP(3)*FNU(6)+3.*DP(6)*FNU(3)+DP(12)*FNU(1)
      V(3)=FNU(13)+3.*DP(4)*FNU(7)+3.*DP(7)*FNU(4)+DP(13)*FNU(1)
      V(4)=FNU(14)+2.*DP(2)*FNU(8)+DP(3)*FNU(5)+DP(5)*FNU(3)+
     1 2.*DP(8)*FNU(2)+DP(14)*FNU(1)
      V(5)=FNU(15)+2.*DP(2)*FNU(9)+DP(4)*FNU(5)+DP(5)*FNU(4)+
     1 2.*DP(9)*FNU(2)+DP(15)*FNU(1)
      V(6)=FNU(16)+2.*DP(3)*FNU(8)+DP(2)*FNU(6)+DP(6)*FNU(2)+
     1 2.*DP(8)*FNU(3)+DP(16)*FNU(1)
      V(7)=FNU(17)+2.*DP(3)*FNU(10)+DP(4)*FNU(6)+DP(6)*FNU(4)+
     1 2.0*DP(10)*FNU(3)+DP(17)*FNU(1)
      V(8)=FNU(18)+2.*DP(4)*FNU(9)+DP(2)*FNU(7)+DP(7)*FNU(2)+
     1 2.*DP(9)*FNU(4)+DP(18)*FNU(1)
      V(9)=FNU(19)+2.*DP(4)*FNU(10)+DP(3)*FNU(7)+DP(7)*FNU(3)+
     1 2.*DP(10)*FNU(4)+DP(19)*FNU(1)
      V(10)=FNU(20)+DP(2)*FNU(10)+DP(3)*FNU(9)+DP(4)*FNU(8)+
     1 DP(10)*FNU(2)+DP(9)*FNU(3)+
     2 DP(8)*FNU(4)+DP(20)*FNU(1)
C
      RETURN
      END
