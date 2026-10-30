      subroutine lrexch(buf,ibuf,dens,ehflx)

c This routine determines the long-range exact exchange energy of the
c SCF density from the long-range integrals erf(omega*r12)/r12 that
c VMOL writes for range-separated hybrids (IIIILR, IJIJLR).  The
c contraction is the exchange part of exchange.F.

      implicit none
#include <aces.h>
#include <iuhf.com>
#include <vmol.com>
#include <machsp.com>
#include <mol.com>
#include <symm2_ks.com>

      integer iiiiio,ijijio
      character*(*) iiiifil,ijijfil
      parameter (iiiiio=10)
      parameter (ijijio=21)
      parameter (iiiifil='IIIILR')
      parameter (ijijfil='IJIJLR')

      integer ibuf(ilnbuf)
      M_REAL
     &    buf(ilnbuf),dens(nbastot,nbastot,iuhf+1),ehflx

      integer
     &    iupki,iupkj,iupkk,iupkl,i,j,k,l,int,nut,iunit,nunit,
     &    alpha,beta,unit(2)
      M_REAL
     &    eexchaa,eexchbb,x,tmp3a,tmp3b
      character*20 files(2)
      logical exist

c Statement functions for unpacking indices.
      iupki(int)=iand(int,ialone)
      iupkj(int)=iand(ishft(int,-ibitwd),ialone)
      iupkk(int)=iand(ishft(int,-2*ibitwd),ialone)
      iupkl(int)=iand(ishft(int,-3*ibitwd),ialone)

      call callstack_push('LREXCH')
      alpha=1
      beta=2
      if (iuhf.eq.0) beta=1

c IIIILR always, IJIJLR when there is symmetry (IIJJLR and IJKLLR give
c no exchange energy)
      nunit=2
      if (c1symmet.eq.1) nunit=1
      unit(1)=iiiiio
      unit(2)=ijijio
      files(1)=iiiifil
      files(2)=ijijfil
      do iunit=1,nunit
        inquire(file=files(iunit),exist=exist)
        if (.not.exist) then
          write (*,9200) files(iunit)
          call errex
        endif
      end do
 9200 format('@LREXCH-F, file ',a12,' not found')

      eexchaa=0.d0
      eexchbb=0.d0

      do 400 iunit=1,nunit
        open(unit(iunit),file=files(iunit),status='old',
     &      form='unformatted',access='sequential')
        call locate(unit(iunit),'TWOELSUP')
    1   read(unit(iunit))buf,ibuf,nut
        do 7 int=1,nut
          i=iupki(ibuf(int))
          j=iupkj(ibuf(int))
          k=iupkk(ibuf(int))
          l=iupkl(ibuf(int))
          x=buf(int)

          tmp3a=x*dens(l,i,alpha)*dens(k,j,alpha)
          tmp3b=x*dens(l,i,beta)*dens(k,j,beta)

          if (i.ne.j .or. k.ne.l) then
            tmp3a=tmp3a+x*dens(k,i,alpha)*dens(l,j,alpha)
            tmp3b=tmp3b+x*dens(k,i,beta)*dens(l,j,beta)
          endif

          if (i.ne.j .and. k.ne.l) then
            tmp3a=2.d0*tmp3a
            tmp3b=2.d0*tmp3b
          endif

c         off diagonal elements
          if (i.ne.k .or. j.ne.l) then
            tmp3a=2.d0*tmp3a
            tmp3b=2.d0*tmp3b
          endif

          eexchaa=eexchaa+tmp3a
          eexchbb=eexchbb+tmp3b
    7   continue
        if (nut.ne.-1) goto 1
        close(unit(iunit),status='keep')
  400 continue

      ehflx=-0.5d0*(eexchaa+eexchbb)

      call callstack_pop
      return
      end
