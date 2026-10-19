      Subroutine Rotate_pairs(Scf_mos,Eps,Index1,Index2,Ip1Sym,Ip2sym,
     +                         Nbasis,Type)

      Implicit Double Precision(A-H,O-Z)
      Character*2 Type 

      Dimension Scf_mos(Nbasis,Nbasis)

      Cos_eps = Cos(Eps)
      Sin_eps = Sin(EPs)

C p' = cos(eps) p + sin(eps) q, q' = sin(eps) p - cos(eps) q for every
C pair type (OO, VV, OV, VO); Get_orbsym sets which of p, q is occupied.
C OV used to store the two results in swapped columns, so eps = 0
C exchanged the occupied and virtual orbitals. p and q are orthonormal
C SCF orbitals and the 2x2 transformation is orthogonal, so p' and q'
C stay orthonormal (no normalization needed). The rows run over all
C functions rather than only those of the pair's irrep.

      Do Irow = 1, Nbasis

         P = Scf_mos(Irow,Index1)
         Q = Scf_mos(Irow,Index2)

         Scf_mos(Irow,Index1) = Cos_eps*P +
     +                          Sin_eps*Q 

         Scf_mos(Irow,Index2) = Sin_eps*P -
     +                          Cos_eps*Q 

      Enddo

      Return
      End
