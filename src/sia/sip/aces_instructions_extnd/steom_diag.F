c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It diagonalize STEOM-EE CI-singles matrices                  |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine steom_diag(x,nindex, type, bval,eval, bdim, edim,
     *           x2, nindex2, type2, bval2, eval2, bdim2, edim2)


      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)
      integer ncis 
c----------------------------------------------------------------

c        write(6,*) nindex,'>>---->',nindex2
      
        ncis = nalpha_occupied * nalpha_virtual 
        call do_steom_diag( x, bval(1),eval(1),type(1),
     *                         bval(2),eval(2),type(2),
     *                   x2, bval2(1),eval2(1),type2(1),
     *                       bval2(2),eval2(2),type2(2),ncis)
       return
       end 
c-------------------------------------------------------------------
c-------------------------------------------------------------------
       subroutine do_steom_diag( A_matr,a1,a2,type1,b1,b2,type2,EVAL,
     & c1,c2,type4,d1,d2,type5,ncis)

      implicit none 
      include 'parallel_info.h'
      include 'daveom.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,aa,bb,IWRK,LCI,info,ICK,KK
      integer type1,type2,type3,type4,type5,type6,type7
      integer ncis 
      double precision A_matr(a1:a2,b1:b2),EVAL(c1:c2,d1:d2)
      double precision HTEST(ncis,ncis),PPP
      double precision EIG(ncis),EIG1(ncis),
     & VECTOR(ncis,ncis),VECTOR1(ncis,ncis),LWRK(8*(ncis)),WR5(ncis),
     & TV(ncis),WR4(ncis)

      IWRK = 8*ncis
      LCI = ncis

      if(me.eq.0) write(6,*) a1,a2,b1,b2,'===>>',c1,c2,ncis
c-----------Matrix for diagonalization--------------------------------

           DO i=1,LCI
           DO j=1,LCI
           HTEST(i,j) = A_matr(i,j)
c           if(me.eq.0) write(6,*) HTEST(i,j),'!!!!!!!!!!!!!!!!!',i,j 
           ENDDO
           ENDDO
c-------------Diagonalization of A - matrix-----------------------------
       CALL DGEEV('V','V',LCI,HTEST,LCI,EIG,EIG1,VECTOR1,LCI,VECTOR,LCI,
     & LWRK,IWRK,info)
          if(me.eq.0) write(6,*) 'info=',info
c--------Sort the eigenvalues in the ascending order

      DO i=1,LCI
      WR4(i) = EIG(i)
      WR5(i) = WR4(i)
      ENDDO

63    CONTINUE

      ICK=0
      DO KK=1,LCI-1
      IF( (WR4(KK+1)).lt.(WR4(KK)) ) then
      PPP = WR4(KK)
      WR4(KK) = WR4(KK+1)
      WR4(KK+1) = PPP
      ICK=ICK+1
      ENDIF
      ENDDO
      IF(ICK.NE.0) goto 63

c-----------------------------------------------------------------------
          if(me.eq.0) write(6,*) 'STEOM eigenvalues are:'
          DO i=1,LCI
          if(me.eq.0) write(6,*) WR4(i)*27.2114,'~~~~~~~>',i
          ENDDO
          if(me.eq.0) write(6,*) '====================================='
c----------Copy target roots into static array--------------------------

          DO i=c1,c2
          DO j=d1,d2
          EVAL(i,j) = WR4(i)*27.2114
          ENDDO 
          ENDDO 
c-----------------------------------------------------------------------
       return
       end
