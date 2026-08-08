










C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
C
C  Added 2026-07-17 (debugging session): independent standalone check of
C  oed__gener_ovl_derv_batch, mirroring compute_1dernai_integrals.F, to
C  cross-check ACES III's overlap-derivative integral via a completely
C  different, much simpler wrapper (no SIAL/contract_oed_derv_integrals.F).
C
      subroutine compute_1derovl_integrals(a1,a2,b1,b2,scr,
     &                              maxblk,iscr,coords,coeffs,
     &                              alphas, ccbeg,
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps,
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical,
     &                              out_ovl_x,out_ovl_y,out_ovl_z,
     &                              nsend,
     &                              Nalpha, Npcoef,
     &                              max_centers)

c---------------------------------------------------------------------------
c  Computes the set of overlap-derivative integrals via oed__gener_ovl_derv_batch
c---------------------------------------------------------------------------
      implicit none
      integer a1, a2, b1, b2, max_centers
      integer adim, bdim
      integer m1, m2, n1, n2
      integer i, n, m
      integer a,b
      integer nshells, nbasis, zmax, intmax, ispherical

      integer nsend
      integer nints, maxblk
      integer npcoeff
      integer ncsum, nfirst
      integer me, nalpha, npcoef
      integer matom, natom

      integer ivangmom(*), ixalpha(*), ixpcoef(*), npfps(*)
      integer ncfps(*), indx_cc(*), Atm_4shell(*)

      logical*8 l8true, l8spherical
      logical spherical

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision Fact
      integer ix, icomponent
      integer DER1X, DER1Y,DER1Z,DER2X,DER2Y,DER2Z
      integer Ea1,Ea2,Eb1,Eb2

      double precision coords(3,*), coeffs(*), alphas(*)
      double precision out_ovl_x(a1:a2,b1:b2)
      double precision out_ovl_y(a1:a2,b1:b2)
      double precision out_ovl_z(a1:a2,b1:b2)

      double precision scr(*)
      integer iscr(*)

      integer ccbeg(*), ccend(*), end_nfps(*)

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff),
     *                 pcoeff_pack(max_dim_coeff)
      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack

      spherical = (ispherical .eq. 1)
      l8spherical = spherical
      l8true = .false.

      adim = Nbasis
      bdim = Nbasis
      nsend = adim*bdim
      Fact = 1.0D0

      call lookup_shell(end_nfps, nshells, a1, m1)
      call lookup_shell(end_nfps, nshells, a2, m2)
      call lookup_shell(end_nfps, nshells, b1, n1)
      call lookup_shell(end_nfps, nshells, b2, n2)

c-------------------------------------------------------------------------
c   Calculate the overlap-derivative integrals over the necessary shell
c   blocks. ix = 1 -> derivative w.r.t. center 1 (matom); ix = 2 ->
c   derivative w.r.t. center 2 (natom). No nuclear-center loop needed
c   (overlap only depends on the two basis centers).
c-------------------------------------------------------------------------
         do ix = 1, 2

           do icomponent = 1, 3

            do b = 1, Nbasis
            do a = 1, Nbasis
               Out_ovl_x(a,b) = 0.0d0
               Out_ovl_y(a,b) = 0.0d0
               Out_ovl_z(a,b) = 0.0d0
            end do
            end do

            do m = m1, m2
               matom = Atm_4shell(m)
               x1 = coords(1,matom)
               y1 = coords(2,matom)
               z1 = coords(3,matom)

            do n = n1, n2
               natom = Atm_4shell(n)
               x2 = coords(1,natom)
               y2 = coords(2,natom)
               z2 = coords(3,natom)

            call pack_coeffs_oed(alphas, ixalpha, coeffs, ixpcoef,
     *                           ncfps, npfps, m, n,
     *                           alpha_pack, nalpha,
     *                           pcoeff_pack, npcoeff,
     *                           ccbeg, ccend, indx_cc,
     *                           ccbeg_pack, ccend_pack,max_dim_coeff)

            ncsum = ncfps(m) + ncfps(n)

            der1x = 0
            der1y = 0
            der1z = 0
            der2x = 0
            der2y = 0
            der2z = 0

            if (ix .eq. 1) then
               if (icomponent .eq. 1) der1x = 1
               if (icomponent .eq. 2) der1y = 1
               if (icomponent .eq. 3) der1z = 1
               if (matom .eq. natom) then
                  if (icomponent .eq. 1) der2x = 1
                  if (icomponent .eq. 2) der2y = 1
                  if (icomponent .eq. 3) der2z = 1
               endif
            else
               if (icomponent .eq. 1) der2x = 1
               if (icomponent .eq. 2) der2y = 1
               if (icomponent .eq. 3) der2z = 1
               if (matom .eq. natom) then
                  if (icomponent .eq. 1) der1x = 1
                  if (icomponent .eq. 2) der1y = 1
                  if (icomponent .eq. 3) der1z = 1
               endif
            endif

            call oed__gener_ovl_derv_batch
     *                (intmax, zmax, nalpha, npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)

            if (nints .gt. 0) then
                Ea2 = end_nfps(m)
                if (m .eq. 1) then
                   Ea1 = 1
                else
                   Ea1 = end_nfps(m-1)+1
                endif

                Eb2 = end_nfps(n)
                if (n .eq. 1) then
                   Eb1 = 1
                else
                   Eb1 = end_nfps(n-1)+1
                endif

              If (Icomponent .Eq. 1) Then
              call add_integrals2(out_ovl_x,a1,a2,b1,b2,scr(nfirst),Ea1,
     *                            Ea2,Eb1,Eb2,1.0D0)
              Elseif (Icomponent .Eq. 2) Then
              call add_integrals2(out_ovl_y,a1,a2,b1,b2,scr(nfirst),Ea1,
     *                            Ea2,Eb1,Eb2,1.0D0)
              Elseif (Icomponent .Eq. 3) Then
              call add_integrals2(out_ovl_z,a1,a2,b1,b2,scr(nfirst),Ea1,
     *                            Ea2,Eb1,Eb2,1.0D0)
              Endif

            endif ! nints

         enddo   ! n shells
         enddo   ! m shells

      Write(6,*)
      Write(6,"(a,1x,i3)") " OVLDER ix=", ix
      Write(6,"(a)") " The d/dr(x)<mu|nu> overlap-derivative integral"
      call output(out_ovl_x, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,"(a)") " The d/dr(y)<mu|nu> overlap-derivative integral"
      call output(out_ovl_y, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,"(a)") " The d/dr(z)<mu|nu> overlap-derivative integral"
      call output(out_ovl_z, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)

         enddo   ! icomponent
         enddo   ! ix

      return
      end
