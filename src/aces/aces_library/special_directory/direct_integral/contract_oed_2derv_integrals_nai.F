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
      subroutine contract_oed_2derv_integrals_nai(type, density, z,
     *                       da1, da2, db1, db2, hess)
c---------------------------------------------------------------------------
c  Computes the set of OED second derivative integrals matching "type" and 
c  contracts them with density block (da1:da2, db1:db2). 
c
c   type = 1 --> NAI integrals
c
      include 'int_gen_parms.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      include 'parallel_info.h'
      include 'trace.h'

      integer type
      integer da1, da2, db1, db2
      double precision density(da1:da2, db1:db2)
      double precision z(da1:da2, db1:db2)
      double precision hess(3*ncenters, 3*ncenters)
      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, dbuf)
#else
      common dbuf
      equivalence (ibuf, dbuf)
#endif

      integer i, j, m, n, nfirst, nints, nalpha, ncoeff, iatom, jatom 
      integer icomponent, jcomponent, ix, jx, ncsum  
      integer a, b, a1,a2,b1,b2,lastao, ihess, jhess 
      integer m1, m2, n1, n2
      integer der1x, der1y, der1z
      integer der2x, der2y, der2z
      integer dercx, dercy, dercz
      integer ixderc, jxderc 
      integer ind(2)
      integer flopcount
      logical spherical
      double precision sumxx, sumxy, sumxz
      double precision sumyx, sumyy, sumyz
      double precision sumzx, sumzy, sumzz
   
      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers)
      double precision x1,y1,z1,x2,y2,z2, x, y, tess, fact

      if (type .lt. 1 .or. type .gt. 3) then
         print *,'Error: Invalid type in contract_oed_2derv_integrals'
         print *,'type = ',type
         call abort_job()
      endif

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif
      spherical = (ispherical .eq. 1)
      lastao    = end_nfps(nshells)

      ind(1) = 1
      ind(2) = 2

#ifdef _DEBUG_LVL0
      sumxx=0.0d0
      sumxy=0.0d0
      sumxz=0.0d0
      sumyx=0.0d0
      sumyy=0.0d0
      sumyz=0.0d0
      sumzx=0.0d0
      sumzy=0.0d0
      sumzz=0.0d0
#endif

      if (intpkg .ne. flocke_package) then
         print *,'Error: contract_oed_2derv_integrals requires ERD ',
     *           'integral package.'
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Find the shells we need to loop over for this block of density data.
c--------------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, da1, m1)
      call lookup_shell(end_nfps, nshells, da2, m2)
      call lookup_shell(end_nfps, nshells, db1, n1)
      call lookup_shell(end_nfps, nshells, db2, n2)
      
c-------------------------------------------------------------------------
c   Build an array of coordinates per atom.
c-------------------------------------------------------------------------

      do m = 1, nshells
         iatom = atom(m)
         xn(iatom) = dbuf(icoord+3*(m-1))
         yn(iatom) = dbuf(icoord+3*(m-1)+1)
         zn(iatom) = dbuf(icoord+3*(m-1)+2)
      enddo

c-------------------------------------------------------------------------
c   Calculate the integrals over the necessary shell blocks.
c-------------------------------------------------------------------------
       do ixderc = 1, Ncenters 
       Iatom = ixderc

       do icomponent = 3, 3  
       do jcomponent = 3, 3   

       do m = m1, m2
          x1 = dbuf(icoord+3*(m-1))
          y1 = dbuf(icoord+3*(m-1)+1)
          z1 = dbuf(icoord+3*(m-1)+2)
       do n = n1, n2
          x2 = dbuf(icoord+3*(n-1))
          y2 = dbuf(icoord+3*(n-1)+1)
          z2 = dbuf(icoord+3*(n-1)+2)

          call pack_coeffs_oed(dbuf(ialpha), ixalpha, dbuf(ipcoeff), 
     *                       ixpcoef, ncfps, npfps, m, n,
     *                       alpha_pack, nalpha, pcoeff_pack,
     *                       ncoeff, ibuf(iccbeg), ibuf(iccend), 
     *                       indx_cc,
     *                       ccbeg_pack, ccend_pack,
     *                       max_dim_coeff )
           ncsum = ncfps(m) + ncfps(n)          


C#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-entry to ERD_batch"
      Write(6,"((a,4(1x,i4)))") "m, n :", m, n
      Write(6,"(a,2(1x,I4))") "The Nalpha_pack and Npcoeff_pack:",
     &                          nalpha, npcoeff
      Write(6,"((a,2(1x,i4)))") "Ivanagmom:", Ivangmom(m),
     &      Ivangmom(n)
      Write(6,"((a,2(1x,i4)))") "Ncfps:", Ncfps(m),
     &      Ncfps(n)
      Write(6,"((a,4(1x,i4)))") "Npfps:", Npfps(m),
     &      Npfps(n)
      Write(6,"(4(1x,F15.7))") (Alpha_pack(i),i=1,Nalpha)
      Write(6,*)
      Write(6,"(4(1x,F15.7))")(Pcoeff_pack(i),i=1,npcoeff)
      Write(6,*)
      Write(6,"(6(1x,I4))")(CCbeg_pack(i),i=1,Nalpha)
      Write(6,*)
      Write(6,"(6(1x,I4))")(CCend_pack(i),i=1,Nalpha)
      Write(6,*)
      Write(6,"(a,3(1x,I10))") "The intmax,Zmax,ncsum:",intmax, zmax,
     &                       ncsum
      Write(6,*) "The shell and nuclear coordinates"
      Write(6,*)
      Write(6,"(3(1x,F12.6))") x1, y1, z1
      Write(6,"(3(1x,F12.6))") x2, y2, z2
CSSS      Write(6,"(3(1x,F12.6))") xn, yn, zn
C#endif

c-------------------------------------------------------------------------
c   Clear the integral block.
c-------------------------------------------------------------------------

            do b = db1, db2
            do a = da1, da2
               Z(a,b) = 0.0D0
            enddo
            enddo

            der1x = 0
            der1y = 0
            der1z = 0
            der2x = 0
            der2y = 0
            der2z = 0
            dercx = 0
            dercy = 0
            dercz = 0

C XX-derivatives

            If (icomponent.eq.1.and.jcomponent.eq.1) then

               dercx = 2
               If (ixderc .eq. atom(m) .and.
     &             ixderc .eq. atom(n)) Then
                   der1x = 2
                   der2x = 2
               Else
                   If (ixderc .eq. atom(m)) then
                       der1x = 2
                   Elseif (ixderc .eq. atom(n)) then
                       der2x = 2
                   Endif
               Endif

C XY-derivatives

            Elseif (icomponent.eq.1.and.jcomponent.eq.2) then

               dercx = 1
               dercy = 1
               If (ixderc .eq. atom(m) .and.
     &             ixderc .eq. atom(n)) Then
                   der1x = 1
                   der1y = 1
                   der2x = 1
                   der2y = 1
               Else
                   If (ixderc .eq. atom(m)) then
                       der1x = 1
                       der1y = 1
                   Elseif (ixderc .eq. atom(n)) then
                       der2x = 1
                       der2y = 1
                  Endif
               Endif

C XZ-derivatives

            Elseif (icomponent.eq.1.and.jcomponent.eq.3) then

               dercx = 1
               dercz = 1
               If (ixderc .eq. atom(m) .and.
     &             ixderc .eq. atom(n)) Then
                   der1x = 1
                   der1z = 1
                   der2x = 1
                   der2z = 1
               Else
                   If (ixderc .eq. atom(m)) then
                       der1x = 1
                       der1z = 1
                   Elseif (ixderc .eq. atom(n)) then
                       der2x = 1
                       der2z = 1
                   Endif
               Endif

C YX-derivatives

            Elseif (icomponent.eq.2.and.jcomponent.eq.1) then

               dercx = 1
               dercy = 1
               If (ixderc .eq. atom(m) .and.
     &             ixderc .eq. atom(n)) Then
                   der1x = 1
                   der1y = 1
                   der2x = 1
                   der2y = 1
               Else
                   If (ixderc .eq. atom(m)) then
                       der1x = 1
                       der1y = 1
                   Elseif (ixderc .eq. atom(n)) then
                       der2x = 1
                       der2y = 1
                   Endif
               Endif

C YY-Derivatives

            Elseif (icomponent.eq.2.and.jcomponent.eq.2) then

               dercy = 2
               If (ixderc .eq. atom(m) .and.
     &             ixderc .eq. atom(n)) Then
                   der1y = 2
                   der2y = 2
               Else
                   If (ixderc .eq. atom(m)) then
                       der1y = 2
                   Elseif (ixderc .eq. atom(n)) then
                       der2y = 2
                   Endif
               Endif

C YZ-Derivatives

            Elseif (icomponent.eq.2.and.jcomponent.eq.3) then

               dercy = 1
               dercz = 1
               If (ixderc .eq. atom(m) .and.
     &             ixderc .eq. atom(n)) Then
                   der1y = 1
                   der1z = 1
                   der2y = 1
                   der2z = 1
               Else
                   If (ixderc .eq. atom(m)) then
                       der1y = 1
                       der1z = 1
                   Elseif (ixderc .eq. atom(n)) then
                       der2y = 1
                       der2z = 1
                   Endif
               Endif

C ZX-Derivatives

            Elseif (icomponent.eq.3.and.jcomponent.eq.1) then

               dercx = 1
               dercz = 1
               If (ixderc .eq. atom(m) .and.
     &             ixderc .eq. atom(n)) Then
                   der1x = 1
                   der1z = 1
                   der2x = 1
                   der2z = 1
               Else
                   If (ixderc .eq. atom(m)) then
                       der1x = 1
                       der1z = 1
                    Elseif (ixderc .eq. atom(n)) then
                       der2x = 1
                       der2z = 1
                    Endif
               Endif

C ZY-Derivatives

            Elseif (icomponent.eq.3.and.jcomponent.eq.2) then

               dercy = 1
               dercz = 1
               If (ixderc .eq. atom(m) .and.
     &             ixderc .eq. atom(n)) Then
                   der1y = 1
                   der1z = 1
                   der2y = 1
                   der2z = 1
               Else
                   If (ixderc .eq. atom(m)) then
                      der1y = 1
                      der1z = 1
                   Elseif (ixderc .eq. atom(n)) then
                      der2y = 1
                      der2z = 1
                   Endif
               Endif

C ZZ-Derivatives

            Elseif (icomponent.eq.3.and.jcomponent.eq.3) then

              If (ixderc .eq. atom(m) .and.
     &             ixderc .eq. atom(n)) Then
                   der1z = 2
                   der2z = 2
                   dercz = 2
               Else
                   If (atom(m) .eq. atom(n)) then
                       der1z = 2
                       der2z = 2
                       dercz = 0
                   Else
                       If (ixderc .eq. atom(m)) then
                           der1z = 2
                           dercz = 2
                       Elseif (ixderc .eq. atom(n)) then
                          der2z = 2
                          dercz = 2
                       Endif
                   Endif
               Endif

            Endif

C#ifdef _DEBUG_LVL0
      Write(6,*)
      write(6,"(a,3I2)") "dercx,dercy,dercz : ", dercx,dercy,dercz
      write(6,"(a,3I2)") "der1x,der1y,der1z : ", der1x,der1y,der1z
      write(6,"(a,3I2)") "der2x,der2y,der2z : ", der2x,der2y,der2z
C#endif

            Call oed__gener_nai_derv_batch(intmax,zmax,
     *                 nalpha,ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n),
     *                 x1,y1,z1,x2,y2,z2,
     *                 ncenters, xn, yn, zn, charge, ixderc,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 dercx, dercy, dercz, alpha_pack,
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., ibuf(ipkgscr), nints, nfirst,
     *                 dbuf(dpkgscr))

            if (nints .gt. 0) then
               a2 = end_nfps(m)
               if (m .eq. 1) then
                  a1 = 1
               else
                  a1 = end_nfps(m-1)+1
               endif

               b2 = end_nfps(n)
               if (n .eq. 1) then
                  b1 = 1
               else
                  b1 = end_nfps(n-1)+1
               endif

              call move_integrals2(z, da1,da2,db1,db2,
     *                             dbuf(dpkgscr+nfirst-1),a1,a2,b1,b2)
      call output(z,da1,da2,db1,db2,da2,db2,1)
c--------------------------------------------------------------------------
c   Contract the density block with the block of integrals.
c--------------------------------------------------------------------------

              call cwork0222(y,
     *                       density,da1,da2,db1,db2,
     *                       da1,da2,db1,db2, ind,
     *                       z,da1,da2,db1,db2,
     *                       da1,da2,db1,db2,ind(1), ind(2),
     *                       flopcount)

c---------------------------------------------------------------------------
c   Sum into the hessian.
c---------------------------------------------------------------------------

              ihess = (iatom-1)*3 + icomponent 
              jhess = (jatom-1)*3 + jcomponent 

              hess(jhess,ihess) = hess(jhess,ihess) + y 

             endif

#ifdef _DEBUG_LVL0
      do b = db1, db2
      do a = da1, da2
         if (icomponent .eq. 1 .and. jcomponent .eq. 1)
     &       sumxx = sumxx + z(a,b) *  z(a,b)
         if (icomponent .eq. 1 .and. jcomponent .eq. 2)
     &       sumxy = sumxy + z(a,b) *  z(a,b)
         if (icomponent .eq. 1 .and. jcomponent .eq. 3)
     &       sumxz = sumxz + z(a,b) * z(a,b)
         if (icomponent .eq. 2 .and. jcomponent .eq. 1)
     &       sumyx = sumyx + z(a,b) *  z(a,b)
         if (icomponent .eq. 2 .and. jcomponent .eq. 2)
     &       sumyy = sumyy + z(a,b) *  z(a,b)
         if (icomponent .eq. 2 .and. jcomponent .eq. 3)
     &       sumyz = sumyz + z(a,b) *  z(a,b)
         if (icomponent .eq. 3 .and. jcomponent .eq. 1)
     &       sumzx = sumzx + z(a,b) *  z(a,b)
         if (icomponent .eq. 3 .and. jcomponent .eq. 2)
     &       sumzy = sumzy + z(a,b) *  z(a,b)
         if (icomponent .eq. 3 .and. jcomponent .eq. 3)
     &       sumzz = sumzz + z(a,b) *  z(a,b)
      enddo
      enddo
#endif
         enddo ! jcomponent
         enddo ! icomponent
         enddo ! ncenters 

         enddo   ! n shells
         enddo   ! m shells

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a,F20.10)")" The 1-el (NAI) hes. integral checksumxx =",
     &sumxx
      write(6,"(a,F20.10)")" The 1-el (NAI) hes. integral checksumxy =",
     &sumxy
      write(6,"(a,F20.10)")" The 1-el (NAI) hes. integral checksumxz =",
     &sumxz
      write(6,"(a,F20.10)")" The 1-el (NAI) hes. integral checksumyx =",
     &sumyx
      write(6,"(a,F20.10)")" The 1-el (NAI) hes. integral checksumyy =",
     &sumyy
      write(6,"(a,F20.10)")" The 1-el (NAI) hes. integral checksumxz =",
     &sumyz
      write(6,"(a,F20.10)")" The 1-el (NAI) hes. integral checksumzx =",
     &sumzx
      write(6,"(a,F20.10)")" The 1-el (NAI) hes. integral checksumzy =",
     &sumzy
      write(6,"(a,F20.10)")" The 1-el (NAI) hes. integral checksumzz =",
     &sumzz
      endif
#endif

      return
      end




