c--------------------------------------------------------------------|
c   This is a superinstruction for SIAL code                         |
c   It performs diagonalization of complex  EOM-SO-CI matrix and     |  
C      prints out roots to the summary.out file                      |
c   Added by Denis Bokhan, Moscow Lomonosov State University,        |
c                  MOSCOW, RUSSIAN FEDERATION                        |
c--------------------------------------------------------------------|



      subroutine steom_so_diag(x,nindex, type, bval,eval,
     & bdim, edim, x2, nindex2, type2, bval2, eval2, bdim2, edim2)

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
      integer bdim2(*), edim2(*),matr_size
      integer ncis 
        
       ncis = Nalpha_occupied * Nalpha_virtual 
c      matr_size = 1 + n_si + n_tr*3
       matr_size = 4*ncis
c----------------------------------------------------------------
       call do_steom_so_diag(x,bval(1),eval(1),type(1),
     &                      bval(2),eval(2),type(2),
     &                  x2, bval2(1),eval2(1),type2(1),
     &                      bval2(2),eval2(2),type2(2),matr_size )
c----------------------------------------------------------------
      return
      end
C=====================================================================
C=====================================================================
      subroutine do_steom_so_diag(x1,a1,a2,type1,b1,b2,type2,
     &                          x2,i1,i2,type3,j1,j2,type4,msize)
      implicit none
      include 'parallel_info.h'
      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,a,b,l1,l2, msize,KK,ICK,jj1,jj2,zz1,zz2
      integer type1,type2,type3,type4,type5,type6,type7,type8
      integer n_si,n_tr,n_s_roots 
      double precision x1(a1:a2,b1:b2),CRE,CIM,PPP,
     & x2(i1:i2,j1:j2), AAA(msize,msize),BBB(msize,msize)
      double complex CMATR(msize,msize),cnum
c-------------Arrays for diagonalization------------------------------
      CHARACTER  JOBVL, JOBVR
      INTEGER    INFO, LDA, LDVL, LDVR, LWORK, N, III(msize)      
      double complex  VL(msize, msize), VR(msize,msize),
     & eval(msize), WORK(4*msize),UUU,VVV
      double precision RWORK(4*msize),WR4(msize),TTT,TTT1,QQQ,QQQ1
c---------------------------------------------------------------------
      double precision E_ccsd

         n_si      = n_singlets 
         n_tr      = n_triplets
         n_s_roots = n_singlets 

         E_ccsd = totenerg - scf_energy
         if(me.eq.0) write(6,*) 'CCSD correlation energy = ',E_ccsd
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
          write(6,*) 'Info=',INFO
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
c      do i=1,msize
c      do j=1,msize
c      x1(i,j) = REAL((VR(i,III(j)))*CONJG(VR(i,III(j))))
c      enddo
c      enddo  
  
c      jj1 = 1 + n_si
c      jj2 = 1 + n_si + 2*n_tr   

c      do i =1, msize
c      do j=1,n_tr

c      VVV = VR(j+jj1,III(i))
c      UUU = VR(j+jj2,III(i))

c      TTT = REAL(VVV)
c      QQQ = AIMAG(VVV)  

c      TTT1 = REAL(UUU)
c      QQQ1 = AIMAG(UUU)    

c      zz1 = 0    
c      IF(DABS(TTT-TTT1).ge.1d-4.or.DABS(QQQ-QQQ1).ge.1d-4) zz1 = 1 
c      IF(DABS(TTT+TTT1).ge.1d-4.or.DABS(QQQ+QQQ1).ge.1d-4) zz1 = 2

c      if(zz1.eq.1) x2(j+jj1,i) = 1.d0
c      if(zz1.eq.2) x2(j+jj1,i) = 2.d0
c      if(zz1.eq.0) x2(j+jj1,i) = 0.d0

c      if(me.eq.0) write(6,*) VR(j+jj1,III(i)),"<===>",VR(j+jj2,III(i)),
c     & zz1, j+jj1,i,"\\\\\\\\\\//////////",j+jj2,III(i)  
c      enddo
c      enddo  
c-------------------------------------------------------------------------
            if(me.eq.0) then 
         write(6,*) '||||||||||||||||||||||||||||||||||||||||||||||' 
c
         do i = 1,  msize
         write(6,*) eval(i) 
         write(6,*) '______________________________________________'
         do j = 1,  msize
c         UUU = dcmplx(0.d0,0.d0)
c         VVV = dcmplx(0.d0,0.d0)
          write(6,*) VR(j,i), '---->',i,j  

c         do k = 1, msize   
c         UUU = UUU + VR(k,i)*CONJG(VR(k,i)) 
c         VVV = VVV + VL(i,k)*CONJG(VL(i,k))
c         enddo

c         write(6,*) UUU,'---->',i,i    
c         write(6,*) VVV,'---->',i,i

         enddo
         write(6,*) '=============================================='  
         enddo

c         write(6,*) '||||||||||||||||||||||||||||||||||||||||||||||'
            endif   
c
c-------------------------------------------------------------------------        

        if(me.eq.0) then
        open(66,file='summary.out')

c         open(66,file='summary.out',form='FORMATTED',status='old',
c     &   access='APPEND')
      write(66,*) '----------------------------------------------------'
      write(66,*) ' Eigenvalues of relativistic STEOM-SO matrix '
      do i=1,n_s_roots
      write(66,*) WR4(i),'----------->',i
      enddo
      write(66,*) '-----------------------------------------'
      write(66,*) '-----Spin-orbit splitted excitation energies------'
      do i=1,n_s_roots
      PPP = WR4(i) 
      PPP = PPP*27.211385
      write(66,*) PPP,'----------->',i 
      enddo   
      write(66,*) '----------------------------------------------------'

        close(66,status='keep')
        endif 
c----------------------------------------------------------------------
      return
      end 
