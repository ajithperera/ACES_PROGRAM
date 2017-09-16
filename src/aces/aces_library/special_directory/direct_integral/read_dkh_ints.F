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
      Subroutine read_dkh_ints()
C
      Implicit none
C
c-----------------------------------------------------------------------------
C Read dkhint.dat file and prepare a common block that contain
C DKH-1el integrals. Ajith Perera, 12/2016.
c----------------------------------------------------------------------------
C
      include 'mpif.h'
      include 'dropmo.h' 


      integer i, j, k, ii, jj, kk, istart,jstart,iend 
      integer ndim, ncheck, Ierror 

      double precision Dkh_ints,ints(max_basis*max_basis)
      double precision W, Ws

      logical bexist 

      common /DKH/ndim, Dkh_ints(2500,2500)  
C
      INQUIRE(file='dkhint.dat',exist=bExist) 

      if (bexist) then 

      OPEN (Unit = 88, File = 'dkhint.dat', Form = 'Unformatted', 
     *      Status = 'Old', Iostat = Ierror) 

      rewind 88
      read(88) ndim 
      ncheck = (-1 + Int(Dsqrt(8.0D0 * ndim + 1.0D0)))/2

      If (ncheck .gt. 2500) Then
         Write(6,"(a,a,a)") " The space allocated for DKH integrals",
     *                      " is limited to maximum of" "2500", 
     *                      " functions."
         Call abort_job()
      Endif 

      Istart = 0
      Iend   = Ndim

      Do j =1, Ndim
         read(88) (ints(Istart+i), i = 1, Iend)
         Istart = Istart + Iend 
         Iend   = Iend - 1
      Enddo

      close(88) 

      else 
           return 
      endif 

      Iend = ndim
      k    = 0

      do i = 1, ndim
         jstart = i 
      do j = 1, Iend
           k = k + 1

           dkh_ints(jstart,i) = ints(k) 
           jstart = jstart + 1

      enddo
         Iend = Iend - 1
      enddo

      do j = 1, ndim
         do i = j+1, ndim
           dkh_ints(j,i) = dkh_ints(i,j)
         enddo 
      enddo
C
C#ifdef _DEBUG_LVL0
C      write(6,*) 
C      write(6,*) "@-read_dkh_ints: The DKH transformed 1-el matrix" 
C      do j=1,ndim
C         Write(6,"(5(1x,F10.5))") (dkh_ints(i,j),i=1,Ndim) 
C      enddo 

      w =0.0d0
      ws=0.0d0

      write(6,*) "@-read_dkh_ints: The checksums" 
      do j=1,ndim
         do i=1,ndim
            w = w + dkh_ints(i,j)
            ws= ws + dkh_ints(i,j) * dkh_ints(i,j)
         enddo
      enddo 
      Write(6,"(a,1x,2F20.10)") "The linear and square-root checksum",
     +                           w,dsqrt(ws)
C#endif 
      
      return
      end
