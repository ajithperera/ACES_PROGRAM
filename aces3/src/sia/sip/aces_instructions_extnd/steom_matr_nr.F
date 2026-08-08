c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It construct STEOM-SO CI-singles matrix                      |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine steom_matr_nr(x,nindex, type, bval,eval, bdim, edim,
     *           x2, nindex2, type2, bval2, eval2, bdim2, edim2)


      implicit none
      include 'interpreter.h'
      include 'trace.h'
      include 'parallel_info.h'
     

      double precision x(*)
      integer nindex, type(*), bval(*), eval(*)
      integer bdim(*), edim(*)
      double precision x2(*)
      integer nindex2, type2(*), bval2(*), eval2(*)
      integer bdim2(*), edim2(*)
c----------------------------------------------------------------

c        write(6,*) nindex,'>>---->',nindex2
        
        if(nindex.eq.2) then
        call do_NR_2ind( x, bval(1),eval(1),type(1),
     *                     bval(2),eval(2),type(2),
     *                x2,bval2(1),eval2(1),type2(1),
     *                   bval2(2),eval2(2),type2(2))
        endif

        if(nindex.eq.4) then
        call do_NR_4ind( x, bval(1),eval(1),type(1),
     *                     bval(2),eval(2),type(2),
     *                     bval(3),eval(3),type(3),
     *                     bval(4),eval(4),type(4),
     *               x2,bval2(1),eval2(1),type2(1),
     *                  bval2(2),eval2(2),type2(2))
        endif


       return
       end

c-------------------------------------------------------------------
c-------------------------------------------------------------------

       subroutine do_NR_2ind(x1,a1,a2,type1,b1,b2,type2,x2,c1,c2,type3,
     & d1,d2,type4 )
      implicit none
      include 'parallel_info.h'
      include 'daveom.h'
      include 'int_gen_parms.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,addr_r,addr_c,a,b,aa,bb
      integer ncis,nvirta,nocca 
      integer type1,type2,type3,type4,type5,type6,type7
      double precision x1(a1:a2,b1:b2),x2(c1:c2,d1:d2)
      integer shift_c,shift_r


      double precision E_ccsd
      E_ccsd = totenerg - scf_energy
      if(me.eq.0) write(6,*) 'CCSD correlation energy = ',E_ccsd

      ncis   = nalpha_occupied * nalpha_virtual 
      nvirta = nalpha_virtual
      nocca  = nalpha_occupied 
      
c      if(me.eq.0) write(6,*) a1,a2,'==========>',b1,b2,'----->'      
      if(me.eq.0) write(6,*) d1,d2,'^^^^^^^^^^>',c1,c2,'@@@@@'

c      IF (RL.eq.0) then
c      do i = c1,c2
c       do j=d1,d2
c       if(i.eq.j) then
c       x2(i,j) = 0.d0 
c       endif 
c       enddo
c       enddo
c       if(me.eq.0) write(6,*) 'Initialization of A-matrix'
c      endif 

c----------------------------------------------------------------|
c         Block 1                                                |
      IF(chunk_type.eq.1) then 
c----------------------------------------------------------------|     

c-----------G_ij - part----------------------------------|

      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta              
           addr_r = nvirta*(i-1) + a
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b
             IF(a.eq.b) then
              x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(i,j)
c             if (i.eq.j) x2(addr_r,addr_c) = x2(addr_r,addr_c) + E_ccsd
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
 

c-----------G_ab - part----------------------------------|

      IF(RL.eq.2) then

       DO i=1,nocca
        DO a=a1,a2  
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb
             IF(i.eq.j) then
              x2(addr_r,addr_c) = x2(addr_r,addr_c) + x1(a,b)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO 
      endif
c------------------------------------------------------------------|
       endif

C====================================================================|
C====================================================================|
      
c----------------------------------------------------------------|
c         Block 6                                                |
      IF(chunk_type.eq.6) then
c----------------------------------------------------------------|     
      shift_r = ncis
      shift_c = ncis
c-----------G_ij - part----------------------------------|

      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
              x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(i,j)
c              if (i.eq.j) x2(addr_r,addr_c) = x2(addr_r,addr_c) + E_ccsd
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|

      IF(RL.eq.2) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
              x2(addr_r,addr_c) = x2(addr_r,addr_c) + x1(a,b)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif


c------------------------------------------------------------------|
       endif

C====================================================================|
C====================================================================|

c----------------------------------------------------------------|
c         Block 11                                               |
      IF(chunk_type.eq.11) then
c----------------------------------------------------------------|     
      shift_r = 2*ncis
      shift_c = 2*ncis
c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
              x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(i,j)
c              if (i.eq.j) x2(addr_r,addr_c) = x2(addr_r,addr_c) + E_ccsd
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|
      IF(RL.eq.2) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
              x2(addr_r,addr_c) = x2(addr_r,addr_c) + x1(a,b)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif

C====================================================================|
C====================================================================|


c----------------------------------------------------------------|
c         Block 16                                               |
      IF(chunk_type.eq.16) then
c----------------------------------------------------------------|     
      shift_r = 3*ncis
      shift_c = 3*ncis
c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
              x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(i,j)
c              if (i.eq.j) x2(addr_r,addr_c) = x2(addr_r,addr_c) + E_ccsd
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|
      IF(RL.eq.2) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
              x2(addr_r,addr_c) = x2(addr_r,addr_c) + x1(a,b)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif


       return
       end
c--------------------------------------------------------|
c--------------------------------------------------------|

       subroutine do_NR_4ind(x1,a1,a2,type1,b1,b2,type2,g1,g2,type3,
     & f1,f2,type4,x2,c1,c2,type5,d1,d2,type6 )
      implicit none
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'daveom.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2,g1,g2,f1,f2
      integer i,j,k,l,addr_r,addr_c,a,b,aa,bb
      integer nocca,nvirta,ncis  
      integer type1,type2,type3,type4,type5,type6,type7
      double precision x1(a1:a2,b1:b2,g1:g2,f1:f2),x2(c1:c2,d1:d2)
      integer shift_r, shift_c
 
      nocca  = Nalpha_occupied 
      nvirta = Nalpha_virtual 
      ncis   = Nalpha_occupied * Nalpha_virtual

c      if(me.eq.0) write(6,*) a1,a2,'=>',b1,b2,'->',g1,g2,'->',f1,f2 
      
c      if(me.eq.0) write(6,*) d1,d2,'^^^^^^^^^^>',c1,c2,'@@@@@'

c----------------------------------------------------------------|
c         Block 1                                                |
      IF(chunk_type.eq.1) then
c----------------------------------------------------------------|     

c---------------G_ai_bj-----------------------------------------|
      IF(RL.eq.3) then
       DO i=f1,f2
        DO a=g1,g2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa
           DO j=b1,b2
            DO b=a1,a2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb
               x2(addr_r,addr_c) = x2(addr_r,addr_c) + 2.d0*x1(b,j,a,i)
            ENDDO
           ENDDO
        ENDDO
       ENDDO       
      ENDIF
c---------------G_ab_ij-----------------------------------------|
      IF(RL.eq.4) then
       DO i=g1,g2
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa
           DO j=f1,f2
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb
               x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(a,b,i,j)
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      ENDIF
c------------------------------------------------------------------|
       endif

C====================================================================|
C====================================================================|

c----------------------------------------------------------------|
c         Block 6                                                |
      IF(chunk_type.eq.6) then
c----------------------------------------------------------------|     
      shift_r = ncis
      shift_c = ncis
c---------------G_ab_ij-----------------------------------------|
      IF(RL.eq.4) then
       DO i=g1,g2
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=f1,f2
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
               x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(a,b,i,j)
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      ENDIF
c------------------------------------------------------------------|
       endif

C====================================================================|
C====================================================================|

c----------------------------------------------------------------|
c         Block 11                                               |
      IF(chunk_type.eq.11) then
c----------------------------------------------------------------|     
      shift_r = 2*ncis
      shift_c = 2*ncis

c---------------G_ab_ij-----------------------------------------|
      IF(RL.eq.4) then
       DO i=g1,g2
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=f1,f2
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
               x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(a,b,i,j)
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      ENDIF
c------------------------------------------------------------------|
       endif

C====================================================================|
C====================================================================|


c----------------------------------------------------------------|
c         Block 16                                               |
      IF(chunk_type.eq.16) then
c----------------------------------------------------------------|     
      shift_r = 3*ncis
      shift_c = 3*ncis

c---------------G_ab_ij-----------------------------------------|
      IF(RL.eq.4) then
       DO i=g1,g2
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=f1,f2
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
               x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(a,b,i,j)
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      ENDIF
c------------------------------------------------------------------|
       endif

C====================================================================|
C====================================================================|
      return
      end

