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
      subroutine comp_return_1el_dkh(z, da1, da2, db1, db2, nsend)
c---------------------------------------------------------------------------
c  Returns 1el-DKH integrals.
c
c---------------------------------------------------------------------------
      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'parallel_info.h'
      include 'dropmo.h'

      integer nsend
      integer da1, da2, db1, db2
      integer a,b,i,k,j,ierror 
      integer adim, bdim, ndim
      double precision z(da1:da2, db1:db2)
      double precision w,ws
      
      logical DKH,Bexist

      double precision Etemp(max_basis*max_basis)
      double precision  Dkh_ints(max_basis,max_basis)

      integer iflags(100),iflags2(500)
      common /flags/ iflags
      common /flags2/ iflags2

c --------------------------------------------------------------------------- 
c  Set up send data. 
c --------------------------------------------------------------------------- 

      adim = da2-da1+1
      bdim = db2-db1+1

      nsend = adim*bdim
      if (nsend .lt. 0) then
          Write(6, "(1x,a,a)") '@-comp_return_1el_dkh:Error: ',
     +          "Dimension of the requested data is zero."
          Write(6,"(1x,a,2I5)") 'The requested block sizes = ',
     +                        adim,bdim
         call abort_job()
      endif

c----------------------------------------------------------------------------
c   Clear the output array.
c----------------------------------------------------------------------------

      DKH = (iflags2(157) .GT. 0)

      If (.NOT. DKH) Return

#ifdef _DEBUG_LVL0
      If (me .EQ. 0 .and. DKH) 
     *    Write(6, "(a,i1,a)")  " The DKH-1el of order ", 
     *    iflags2(157), " is added to the 1-el Hamiltonian" 
#endif 

      if (intpkg .ne. flocke_package) then
         Write(6, "(1x,a,a)") "@-comp_return_1el_dkh:Error: ",
     +           "comp_return_1el_dkh requires OED 'integral package."
         call abort_job()
      endif

      Inquire(file='dkhint.dat',exist=Bexist)

      if (Bexist) then
      Open (Unit = 71, File = 'dkhint.dat', Form = 'Unformatted',
     *      Status = 'Old', Iostat = Ierror)

      rewind 71
      read(71)  Ndim
      If (Ndim .Gt. max_basis) then
         Write(6,"(a,a,i4)") " The maximum number of basis functions",
     +                    " allowed with the DKH option is",max_basis
         Write(6,"(a,i5)") " The current calculation requested", Ndim,
     +                    "functions."
         call abort_job()
      Endif 
      Read(71) (Etemp(i), i=1, Ndim*(Ndim+1)/2)
      Close(71)

      else
           Write(6,"(a,a)") " DKH option is requested but dkhint.dat",
     +                      " file is not found!"
           Call abort_job()
      endif 

      k = 0
      do i = 1, ndim
      do j = 1, i
           k = k + 1
           Dkh_ints(i,j) = Etemp(k) 
           Dkh_ints(j,i) = Etemp(k) 
      enddo
      enddo

      do b = db1,db2
      do a = da1,da2
         z(a,b) = Dkh_ints(a,b)
      enddo
      enddo

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a)") "@-comp_return_1el_dkh: The DKH 1-el Hamiltonian"
C      write(6,"(6(1x,F10.5))")((z(a,b),a=da1,da2),b=db1,db2)
      call output(z,db1,db2,da1,da2,adim,bdim,1)
      ws = 0.0d0
      w  = 0.0d0
      do b = db1,db2
      do a = da1,da2
      w  = w +   z(a,b)
      ws = ws +  z(a,b) *  z(a,b)
      enddo
      enddo
      write(6,"(a,2F20.10)")"The linear and square-root checksum:",
     +                       w, dsqrt(ws)
#endif

      return
      end

