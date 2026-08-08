










      Subroutine Lanczos_main(Vecs,Work,Memleft,Nsize,Irrepx,Iuhf,
     +                        Maxmem)
 
      Implicit Double Precision(A-H,O-Z)


 
C Maxn : Lead dimension of the tridiagonal matrix
C Info : Serve as both input/output. Only Info(1) is 
C        set by the caller. 
C Maxvw: This is set to 1. This allocate one extra vector 
C        of the diemension of Hbar (this space is not used).
C Convrg: Strict or light 

      Integer Maxn
      Parameter(Maxn_int=10000,Maxvw_int=1)
      Character*6 Convrg 

      Common /lanczos_vars/Maxn,Maxvw,Trd(Maxn_int,Maxn_int),
     +                     Convrg 
      




C This  subroutine uses  the Lanczos  algorithm  with look-ahead to
C compute eigenvalue  estimates.  The algorithm was first described
C in  the  RIACS Technical  Report 90.45, `An Implementation of the
C Look-Ahead Lanczos Algorithm for Non-Hermitian Matrices, Part I`,
C by R.W. Freund, M.H. Gutknecht and N.M. Nachtigal, November 1990.
C The routine will first run the look-ahead Lanczos algorithm for a
C number of steps.  It will then prompt the user  for the number of
C eigenvalues to compute.  In addition, the routine will  also find
C the common eigenvalues in two successive runs.

      Integer Revcom,Colx,Colb,Info(4),Retlabel 
      Dimension Vecs(Nsize,2*Maxvw+6)
      Dimension Work(Memleft)
      Dimension orthoInit(NSize*2)
      Write(6,"(a)") "--- Entered Lanczos_main ---"
      Write(6,"(a)") " The input vars,memleft,Nsize,Irrepx,Iuhf"
      Write(6,"(1x,ES8.2E2,1x,ES8.2E2,1x,i2,1x,i2)") Dble(Memleft),
     +          Dble(Nsize),Irrepx,iuhf

      If (Maxn .Gt. Maxn_int) Then
         Write(6,"(a,a)") " The number of requested Lanczos vectors",
     +                    " are more than the currently allowed"
         Write(6,"(a,I6)")" maximum of ", No_max_vecs
         Call Errex
      Endif

C For Lanczos iterations the following lists are used
C R1_in: 490(1,2)  R1_out: 490(3,4)
C R2_in: 444-446   R2_out: 461-463
C L1_in: 493(1,2)  L1_out: 493(3,4)
C L2_in: 448-450   L2_out: 467-469
C Also availble 448-450 for R2 or L2

C To initilize read both left and right starting vectors 

      Mubar_s_pq_t0 = 490
      Mubar_d_pq_t0 = 443
      Mubar_s_pq_tn = 493
      Mubar_d_pq_tn = 447

C Load Mubar_R0

      Ioffr1 = 0
      Ioffr2 = 0
!      Call Lanczos_Load_vec(Irrepx,Vecs(1,1),Nsize,Mubar_s_pq_t0,Ioffr1,
!     +                      Mubar_d_pq_t0,Ioffr2,Iuhf,.False.)


!        Mubar_s_pq_t0=490
!        Mubar_d_pq_t0=460
!         Call Lanczos_dump_vec(Irrepx,Vecs(1,1),Nsize,Mubar_s_pq_t0,
!     +                      Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
!     +                      Iuhf,.False.)
!!       print*, 'calling hbarxc'
!          call HBARXC(Work,Memleft, Iuhf, 1,Irrepx,.True.)
!        ! Call Tdee_hbar_mult(Work,Memleft,Iuhf,Irrepx,1)
!
!C This Returns Mubar_R(n+1), iside=1
!        Ioffr1=0
!        Ioffr2=0
!        Mubar_s_pq_tn=490
!        Mubar_d_pq_tn=460
!!         Call Lanczos_Load_vec(Irrepx,Vecs(1,1),Nsize,Mubar_s_pq_tn,
!!     +                     Ioffr1,Mubar_d_pq_tn,Ioffr2,Iuhf,.False.)
!        call loadvec1(irrepx,Vecs(1,1),Nsize,iuhf,490,2,460,
!     +          NSIZEC,.false.)




C Load Mubar_L0 

      Ioffr1 = 3
      Ioffr2 = 4
      Ioffsp=0
      call getlst(Vecs(1,1),2,1,1,1,472)
      call getlst(Vecs(1,2),3,1,1,1,472)
!      Call Lanczos_Load_vec(Irrepx,Vecs(1,2),Nsize,Mubar_s_pq_t0,Ioffr1,
!     +                      Mubar_d_pq_t0,Ioffr2,Iuhf,.False.)

         Call Lanczos_dump_vec(Irrepx,Vecs(1,2),Nsize,Mubar_s_pq_t0,
     +                      Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
     +                      Iuhf,.False.)

      Ioffr1 = 0
      Ioffr2 = 0
         Call Lanczos_dump_vec(Irrepx,Vecs(1,1),Nsize,Mubar_s_pq_t0,
     +                      Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
     +                      Iuhf,.False.)
!
!!       print*, 'going to hbarxc'
!        call HBARXC(Work,Memleft, Iuhf, 2,Irrepx,.True.)
!       !  Call Tdee_hbar_mult(Work,Memleft,Iuhf,Irrepx,2)
!!       print*, 'out of hbarxc'
!C This Returns Mubar_L(n+1), iside=2
!        call loadvec1(irrepx,Vecs(1,2),Nsize,iuhf,493,2,466,
!     +          NSIZEC,.false.)
!      Call Getlst(Vecs(1,2),2,1,1,1,94)
!#      call checksum("vecs1,2:",Vecs(1,2),nsize,s)
!#      orthoInit=0.0d0
      call dcopy(Nsize,Vecs(1,1),1,orthoInit(1),1)
!#
      call dcopy(NSize,Vecs(1,2),1,orthoInit(1+Nsize),1)
!#      print*,'printting ortho1,2'
!#      call checksum("orthoinit1,2:",orthoInit(1+Nsize),nsize,s)
!      call GramSchmidt(orthoInit,NSize,2)
!      call dcopy(Nsize,orthoInit(1),1,Vecs(1,1),1)
!      call dcopy(NSize,orthoInit(1+Nsize),1,Vecs(1,2),1)
!#      call checksum("vecs1,2:",Vecs(1,2),nsize,s)
!      print*,'overlap',dot_product(Vecs(:,1),Vecs(:,2))


      !Call Dscal(Nsize,0.75d0,Vecs(1,2),1)
C Du-lal stands for Double Precision Lanczos, LAL refer to three
C term look-ahead Lanczos. 
      Info(1) = 1

      Info(2) = 0
      Info(3) = 0
      Ncount  = 0
      Icount  = 0
      Jcount  = 0
      Retlabel= 1

 10   Continue 

      Ncount = Ncount + 1

      Write(6,*)
      Write(6,"(a)") "Input for Dulal or Lanczos_light"
      If (ncount .eq.1) then
         call checksum("Lanczos_main:",Vecs(1,1),nsize,s)
         call checksum("Lanczos_main:",Vecs(1,2),nsize,s)
      else 
         if (revcom.eq.1) call checksum("Lanczos_main:",
     +                     Vecs(1,colb),nsize,s)
         if (revcom.eq.2) call checksum("Lanczos_main:",
     +                     Vecs(1,colb),nsize,s)
      endif 
C Always do one more iteration than requested. This is useful to asses 
C the convegence of Lanczos chains toward convergence. 

      If (Ncount .gt. 2*Maxn+4) Go to  99
      Call Lanczos_light_v1(Nsize,Vecs,Info,Retlabel) 

      Ierror = Info(1)
      Write(6,"(a,i2)") "Monitoring variable ", Ierror 
      If (Ierror .EQ.0) Go to 99
      
      Revcom = Info(2)
      Colx   = Info(3)
      Colb   = Info(4)

      Write(6,"(a)") "The Revcom,Colx,Colb vars returned"
      Write(6,"(5(1x,i2))") Info(2),Info(3),Info(4),Colx,Colb
      if (revcom.eq.1) call checksum("Lanczos_main:",
     +                                Vecs(1,colx),nsize,s)
      if (revcom.eq.2) call checksum("Lanczos_main:",
     +                                Vecs(1,colx),nsize,s)
      If (Revcom .EQ. 1) Then

         If (Convrg .EQ. "Strict") Then 
            Icount = Icount + 1
            Write(6,*) "Icount & Ncount", Icount, Ncount 
C Do the external biorthogonalization of Lanczos vectors. This should
C prevent the deterioration of biorthogonality of Lanczos vectors. 
C First right vector with all the left vectors.

            If (Ncount .Gt. 1) Call Lanczos_gs_ortho_rl(Work,Memleft,
     +                              Vecs(1,Colx),Nsize,Icount,
     +                              Irrepx)

            If (Ncount .eq.1) Call Updmoi(maxn+4,nsize,Irrepx,497,0,0)
            Call Putlst(Vecs(1,Colx),Icount,1,0,Irrepx,497)
  
         Else if (Convrg .EQ. "Light") Then
            Icount = Icount + 1
            Write(6,*) "Icount & Ncount", Icount, Ncount 
            If (Icount .EQ. Maxn+2) Then
               Call Updmoi(Icount-1,nsize,Irrepx,497,0,0)
               Call Putlst(Vecs(1,Colx),Icount-1,1,0,Irrepx,497)
            Endif 
 
         Endif 

         Ioffr1 = 0
         Ioffr2 = 0
         Ioffsp = 2

C Do the right hand multiplication and continue
      Write(6,*)
      Write(6,"(a,i3)") "Dumping Vector Vec(1,Colx)",Colx
      call checksum("value  :",Vecs(1,Colx),Nsize,S)

        Mubar_s_pq_t0=490
        Mubar_d_pq_t0=443
         Call Lanczos_dump_vec(Irrepx,Vecs(1,Colx),Nsize,Mubar_s_pq_t0,
     +                      Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
     +                      Iuhf,.False.)

         Call Hbarxc(Work,Memleft,Iuhf,1,Irrepx) 

C This Returns Mubar_R(n+1), iside=1
        Ioffr1=0
        Ioffr2=0
        Mubar_s_pq_tn=490
        Mubar_d_pq_tn=460
        call loadvec1(irrepx,Vecs(1,Colb),NSize,IUHF,Mubar_s_pq_tn,2,
     +                  Mubar_d_pq_tn,NSIZEC,.FALSE.)
!         Call Lanczos_Load_vec(Irrepx,Vecs(1,Colb),Nsize,Mubar_s_pq_tn,
!     +                     Ioffr1,Mubar_d_pq_tn,Ioffr2,Iuhf,.False.)
         Go to 10 

      Elseif (Revcom .EQ. 2) Then

         If (Convrg .EQ. "Strict") Then 
            Jcount = Jcount + 1

C Do the external biorthogonalization of Lanczos vectors. This should
C prevent the deterioration of biorthogonality of Lanczos vectors.
C Now left vector with all the right vectors.

            If (Ncount .Gt. 1) Call Lanczos_gs_ortho_lr(Work,Memleft,
     +                              Vecs(1,Colx),Vecs(1,3),Nsize,Jcount,
     +                              Irrepx)
            If (Ncount .eq.2) Call updmoi(maxn+4,nsize,Irrepx,498,0,0)
         Write(6,*) "Jcount & Ncount",jcount, Ncount 
           Call Putlst(Vecs(1,Colx),Jcount,1,0,Irrepx,498)

         Else if (Convrg .EQ. "Light") Then
            Jcount = Jcount + 1
            Write(6,*) "Jcount & Ncount", Jcount, Ncount 
            If (Jcount .EQ. Maxn+2) Then
               Call Updmoi(Jcount-1,nsize,Irrepx,498,0,0)
               Call putlst(Vecs(1,Colx),Jcount-1,1,0,Irrepx,498)
            Endif 
  
         Endif 

!         Ioffr1 = 0
!         Ioffr2 = 4
!         Ioffsp = 0


         Ioffr1 = 0
         Ioffr2 = 0
         Ioffsp = 0
        Mubar_s_pq_t0=493
        Mubar_d_pq_t0=447
      Write(6,*)
      Write(6,"(a,i3)") "Dumping left Vector Vec(1,Colx)",Colx
      call checksum("AA  :",Vecs(1,Colx),Nsize,S)
         Call Lanczos_dump_vec(Irrepx,Vecs(1,Colx),Nsize,Mubar_s_pq_t0,
     +                    Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
     +                    Iuhf,.False.)

         Call Hbarxc(Work,Memleft,Iuhf,2,Irrepx) 

C This Returns Mubar_L(n+1), iside=2
        Ioffr1=0
        Ioffr2=0!-13
        Mubar_s_pq_tn=493
        Mubar_d_pq_tn=466
        call loadvec1(irrepx,Vecs(1,Colb),NSize,IUHF,Mubar_s_pq_tn,2,
     +                  Mubar_d_pq_tn,NSIZEC,.FALSE.)
!         Call Lanczos_Load_vec(Irrepx,Vecs(1,Colb),Nsize,Mubar_s_pq_tn,
!     +                      Ioffr1,Mubar_d_pq_tn,Ioffr2,Iuhf,.False.)
         Go to 10 

      Endif 

  99  Continue

      Write(6,"(a)") "Lanczos iterations are completed!" 

C#ifdef 1
      If (Convrg .EQ. "Strict") Then
          Call check_lanczos(Work,Maxmem,Nsize,irrepx,Maxn)
      Endif 
C#endif 
    
C At this point non-symmetric tridiagonal matirx is complete. Digonalize
C (using dgeev) and obtain left and right eigenvectors and eigenvalues.

      Call Lanczos_diag_trd(Work,Maxmem,Ncount,Nsize,Iuhf,Irrepx)

      Return
      End
