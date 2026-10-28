










      DOUBLE PRECISION FUNCTION LINDH_RHO(IZA, IZB, R2)
C
C rho_ab = exp(alpha_ab (r_ref,ab^2 - r_ab^2)) of the Lindh model
C Hessian (R. Lindh et al., Chem. Phys. Lett. 241, 423 (1995)). IZA,
C IZB are atomic numbers, R2 the squared distance (bohr^2). alpha and
C r_ref depend on the periodic-table rows (rows beyond the third use
C the third-row values). 10/2026.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ALPHA(3,3), RREF(3,3)
      DATA ALPHA /1.0000D0, 0.3949D0, 0.3949D0,
     &            0.3949D0, 0.2800D0, 0.2800D0,
     &            0.3949D0, 0.2800D0, 0.2800D0/
      DATA RREF  /1.35D0, 2.10D0, 2.53D0,
     &            2.10D0, 2.87D0, 3.40D0,
     &            2.53D0, 3.40D0, 3.40D0/
C
      IRA = 3
      IF (IZA .LE. 10) IRA = 2
      IF (IZA .LE. 2)  IRA = 1
      IRB = 3
      IF (IZB .LE. 10) IRB = 2
      IF (IZB .LE. 2)  IRB = 1
      LINDH_RHO = DEXP(ALPHA(IRA,IRB)*(RREF(IRA,IRB)**2 - R2))
C
      RETURN
      END
