










      Subroutine blt_hfksFai(C,iuhf,ld3,Fock,scr1,ld2,ld1,
     &        scrtmp,CT,tfia,tfib,
     &        occvrt1,occvrt2,
     &        pop)

      implicit none







c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>




c symm2.com : begin

c This is initialized in vscf/symsiz.

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      integer nirrep,      nbfirr(8),   irpsz1(36),  irpsz2(28),
     &        irpds1(36),  irpds2(56),  irpoff(9),   ireps(9),
     &        dirprd(8,8), iwoff1(37),  iwoff2(29),
     &        inewvc(maxbasfn),         idxvec(maxbasfn),
     &        itriln(9),   itriof(8),   isqrln(9),   isqrof(8),
     &        mxirr2
      common /SYMM2/ nirrep, nbfirr, irpsz1, irpsz2, irpds1, irpds2,
     &               irpoff, ireps,  dirprd, iwoff1, iwoff2, inewvc,
     &               idxvec, itriln, itriof, isqrln, isqrof, mxirr2
c symm2.com : end


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



       integer ld1,ld2,ld3,iuhf
       integer occvrt2,occvrt1
       integer pop(8,2)
       integer indx2

       double precision
     &       C((IUHF+1)*LD3),
     &       fock((iuhf+1)*ld1),
     &       scr1(ld2),
     &       scrtmp(ld2),
     &       CT(ld2),
     &       tfia(occvrt1),
     &       tfib(occvrt2)
      integer i,j,n,a,ind,ispin
      INDX2(I,J,N)=I+(J-1)*N

      do ispin=1,iuhf+1
         ind=1
        do i=1,nirrep
c         if (ispin .eq. 1 ) goto 102
         if(nbfirr(i) .eq.0) goto 101
         call expnd2(fock(((ispin-1)*ld1)+itriof(i)),scr1,
     &              nbfirr(i))
         call transp(C(((ispin-1)*ld3)+isqrof(i)),CT,
     &             nbfirr(i),nbfirr(i))
         call mxm(CT,nbfirr(i),scr1,nbfirr(i),scrtmp,nbfirr(i))
         call mxm(scrtmp,nbfirr(i),C(((ispin-1)*ld3)+isqrof(i)),
     &          nbfirr(i),scr1,nbfirr(i))
c       call kkk(nbfirr(i)*nbfirr(i),scr1)
c 102    continue
         do j=1,pop(i,ispin)
           do a=pop(i,ispin)+1,nbfirr(i)
            if(ispin .eq. 1) then
c                tfia(ind)=scr1(j+(a-1)*nbfirr(i))
               tfia(ind)=scr1(indx2(j,a,nbfirr(i)))
            ind=ind+1
           else
               tfib(ind)=scr1(indx2(j,a,nbfirr(i)))
               ind=ind+1
            end if
           end do 
         end do
101    continue
      end do
      end do
      call putrec(20,'JOBARC','HFKSFIAA',occvrt1*iintfp,
     &             tfia)
      if(iuhf .eq. 1) then
      call putrec(20,'JOBARC','HFKSFIAB',occvrt2*iintfp,
     &             tfib)
      end if
      return
      end   
