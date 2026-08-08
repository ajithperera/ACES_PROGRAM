c--------------------------------------------------------------------|
c   This is a superinstruction for SIAL code                         |
c   It performs diagonalization of complex  EOM-SO-CI matrix and     |  
C      prints out roots to the summary.out file                      |
c   Added by Denis Bokhan, Moscow Lomonosov State University,        |
c                  MOSCOW, RUSSIAN FEDERATION                        |
c--------------------------------------------------------------------|



      subroutine eom_so_diag(x,nindex, type, bval,eval,
     & bdim, edim, x2, nindex2, type2, bval2, eval2, bdim2, edim2)

      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'spin_orb.h'
      include 'int_gen_parms.h'

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*),matr_size
    
      n_si = n_singlets 
      n_tr = n_triplets 
        
      matr_size = 1 + n_si + n_tr*3
      if(ip_roots.ne.0.and.ea_roots.ne.0) matr_size = n_si + n_tr*3
c----------------------------------------------------------------
      call do_eom_so_diag(x,bval(1),eval(1),type(1),
     &                      bval(2),eval(2),type(2),
     &                  x2, bval2(1),eval2(1),type2(1),
     &                      bval2(2),eval2(2),type2(2),matr_size )
c----------------------------------------------------------------
      return
      end
C=====================================================================
C=====================================================================
      subroutine do_eom_so_diag(x1,a1,a2,type1,b1,b2,type2,
     &                          x2,i1,i2,type3,j1,j2,type4,msize)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2, msize,KK,ICK,jj1,jj2,zz1,zz2
      integer type1,type2,type3,type4,type5,type6,type7,type8
      integer n_si,n_tr
      double precision x1(a1:a2,b1:b2),CRE,CIM,PPP,
     & x2(i1:i2,j1:j2), AAA(msize,msize),BBB(msize,msize)
      double complex CMATR(msize,msize),cnum
c-------------Arrays for diagonalization------------------------------
      CHARACTER  JOBVL, JOBVR
      INTEGER    INFO, LDA, LDVL, LDVR, LWORK, N, III(msize)      
      double complex  VL(msize, msize), VR(msize,msize),
     & eval(msize), WORK(4*msize),UUU,VVV
      double precision RWORK(4*msize),WR4(msize),TTT,TTT1,QQQ,QQQ1
      double precision Conv 
c---------------------------------------------------------------------

         n_si = n_singlets 
         n_tr = n_triplets 
c-----------Copy static arrays to the temporal arrays-----------------
         do i = 1,  msize
         do j = 1,  msize
         AAA(i,j) = x1(i,j)
         BBB(i,j) = x2(i,j) 
         enddo
         enddo 
c-----------Create complex matrix for the diagonalization-------------
         do i=1,msize
         do j=1,msize  
         CRE = AAA(i,j)
         CIM = BBB(i,j)
c         if (me.eq.0) write(6,*) '>>>',i,j,'===',CIM
c         if (me.eq.0) write(6,*) '>>>',i,j,'===',CRE 
         cnum = dcmplx(CRE,CIM)
         CMATR(i,j) = cnum
c          if (me.eq.0) write(6,*) '>>>',i,j,'===',CMATR(i,j)
         enddo
         enddo 
      
         do i=1,msize
         do j=1,msize
         CRE = AAA(i,j)
c         if (me.eq.0) write(6,*) '<<<',i,j,'===',CRE
         enddo
         enddo


c-----------Diagonalize complex EOM-SO-CI matrix----------------------- 
         JOBVL='N'
         JOBVR='V' 
         LWORK = 4*msize
         call zgeev(JOBVL,JOBVR,msize,CMATR,msize,eval,VL,msize,VR,
     &  msize, WORK, LWORK, RWORK, INFO)
         if(me.eq.0) then 
C          write(6,*) 'Info=',INFO
          do i = 1, msize
c          write(6,*) eval(i), '%%%', i 
          WR4(i) = REAL(eval(i))
          enddo 
          endif
c-----------Sort eigenvalues in acending order-------------------------
63    CONTINUE

      ICK=0
      DO KK=1,msize-1
      IF( (WR4(KK+1)).lt.(WR4(KK)) ) then
      PPP = WR4(KK)
      WR4(KK) = WR4(KK+1)
      WR4(KK+1) = PPP
      ICK=ICK+1
      ENDIF
      ENDDO
      IF(ICK.NE.0) goto 63

      call izero(III,msize)
      DO i=1,msize
      DO j=1,msize
      IF(REAL(eval(i)).eq.WR4(j)) III(j)=i
      ENDDO
      ENDDO

      DO i=1,msize
      if(me.eq.0) write(6,*) i,'--->>>',III(i)
      ENDDO

C-------------Zero out arrays x1 and x2----------------------------------
      do i=1,msize
      do j=1,msize 
      x1(i,j) = 0.d0
      x2(i,j) = 0.d0
      enddo
      enddo

c-------Calculate module of EOM-SO-CI coefficients----------------------
      do i=1,msize
      do j=1,msize
      x1(i,j) = REAL((VR(i,III(j)))*CONJG(VR(i,III(j))))
      enddo
      enddo  
  
      jj1 = 1 + n_si
      jj2 = 1 + n_si + 2*n_tr   

      do i =1, msize
      do j=1,n_tr

      VVV = VR(j+jj1,III(i))
      UUU = VR(j+jj2,III(i))

      TTT = REAL(VVV)
      QQQ = AIMAG(VVV)  

      TTT1 = REAL(UUU)
      QQQ1 = AIMAG(UUU)    

      zz1 = 0    
      IF(DABS(TTT-TTT1).ge.1d-4.or.DABS(QQQ-QQQ1).ge.1d-4) zz1 = 1 
      IF(DABS(TTT+TTT1).ge.1d-4.or.DABS(QQQ+QQQ1).ge.1d-4) zz1 = 2

      if(zz1.eq.1) x2(j+jj1,i) = 1.d0
      if(zz1.eq.2) x2(j+jj1,i) = 2.d0
      if(zz1.eq.0) x2(j+jj1,i) = 0.d0

c      if(me.eq.0) write(6,*) VR(j+jj1,III(i)),"<===>",VR(j+jj2,III(i)),
c     & zz1, j+jj1,i,"\\\\\\\\\\//////////",j+jj2,III(i)  
      enddo
      enddo  
c-------------------------------------------------------------------------
            if(me.eq.0) then 
C         write(6,*) '||||||||||||||||||||||||||||||||||||||||||||||' 
c
         do i = 1,  msize
C         write(6,*) eval(i) 
C         write(6,*) '______________________________________________'
         do j = 1,  msize
c         UUU = dcmplx(0.d0,0.d0)
c         VVV = dcmplx(0.d0,0.d0)
C          write(6,*) VR(j,i), '---->',i,j  

c         do k = 1, msize   
c         UUU = UUU + VR(k,i)*CONJG(VR(k,i)) 
c         VVV = VVV + VL(i,k)*CONJG(VL(i,k))
c         enddo

c         write(6,*) UUU,'---->',i,i    
c         write(6,*) VVV,'---->',i,i

         enddo
C         write(6,*) '=============================================='  
         enddo

c         write(6,*) '||||||||||||||||||||||||||||||||||||||||||||||'
            endif   
C-----------------------------------------------------------------------
      Conv = 27.21138505
      If (me .eq. 0)then
         If(ip_roots.ne.0.and.ea_roots.ne.0) then
           write(66,10)
           write(66,100)
           write(66,30)
           write(66,100)
           do i = 1, Msize
              write(66,20) i, Wr4(i),  Wr4(i)*conv
           enddo
           write(66,100)
         Else
           write(66,11)
           write(66,100)
           write(66,30)
           write(66,100)
           do i=1,msize
              write(66,20) i, (Wr4(i)-WR4(1)),(Wr4(i)-WR4(1))*conv
           Enddo 
           write(66,100)
        Endif 
      Endif 

100   format(2x,60('-'))
10    format(2x,"SOC split STEOM-EOM-CCSD energies")
11    format(2x,"SOC split EOM-CCSD energies")
20    format(2x,I4,5x,F20.8,x,F20.8)
30    format(2x,'ROOT',15x,'Hartree',15x,'eV')

c----------------------------------------------------------------------
      return
      end 
