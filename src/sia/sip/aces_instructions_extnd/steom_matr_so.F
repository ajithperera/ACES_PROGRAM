c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It construct STEOM-SO CI-singles matrix                      |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine steom_matr_so(x,nindex, type, bval,eval, bdim, edim,
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
        call do_SO_2ind( x, bval(1),eval(1),type(1),
     *                     bval(2),eval(2),type(2),
     *                x2,bval2(1),eval2(1),type2(1),
     *                   bval2(2),eval2(2),type2(2))
        endif



       return
       end

c-------------------------------------------------------------------
c-------------------------------------------------------------------

       subroutine do_SO_2ind(x1,a1,a2,type1,b1,b2,type2,x2,c1,c2,type3,
     & d1,d2,type4 )
      implicit none
      include 'parallel_info.h'
      include 'daveom.h'
      include 'int_gen_parms.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,addr_r,addr_c,a,b,aa,bb
      integer nocca,nvirta,ncis 
      integer type1,type2,type3,type4,type5,type6,type7
      double precision x1(a1:a2,b1:b2),x2(c1:c2,d1:d2),sqq
      integer shift_c,shift_r

      Nocca  = Nalpha_occupied 
      Nvirta = Nalpha_virtual 
      Ncis   = Nalpha_occupied * Nalpha_virtual

      if(me.eq.0) write(6,*) d1,d2,'^^^^^^^^^^>',c1,c2,'@@@@@'

C====================================================================|
C====================================================================|

      sqq = 1.d0/(dsqrt(2.d0))
      

c----------------------------------------------------------------|
c         Block 11                                               |
      IF(chunk_type.eq.11) then
c----------------------------------------------------------------|     

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|


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
              x2(addr_r,addr_c) = x2(addr_r,addr_c) + x1(i,j)*(-1.d0)
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
              x2(addr_r,addr_c) = x2(addr_r,addr_c) + x1(a,b)*(-1.d0)
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

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|

 
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
              x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(i,j)*(-1.d0)
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
              x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif

c----------------------------------------------------------------|
c          Block 2                                               |
      IF(chunk_type.eq.2) then
c----------------------------------------------------------------|     

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|


      shift_r = 0
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
              x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(i,j)*(-1.d0)
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
              x2(addr_r,addr_c) = x2(addr_r,addr_c) + x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif

c----------------------------------------------------------------|
c          Block 5                                               |
      IF(chunk_type.eq.5) then
c----------------------------------------------------------------|     

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|


      shift_r = ncis
      shift_c = 0
c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
              x2(addr_r,addr_c) = x2(addr_r,addr_c) - x1(i,j)*(-1.d0)
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
              x2(addr_r,addr_c) = x2(addr_r,addr_c) + x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif


c----------------------------------------------------------------|
c          Block 3                                               |
      IF(chunk_type.eq.3) then
c----------------------------------------------------------------|     

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|


      shift_r = 0
      shift_c = 2*ncis


c---------------F_so_x part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) + sqq*x1(i,j)*(-1.d0)
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
             x2(addr_r,addr_c)=x2(addr_r,addr_c) - sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif

c---------------F_so_y part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.10) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) - sqq*x1(i,j)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|
      IF(RL.eq.20) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
             x2(addr_r,addr_c)=x2(addr_r,addr_c) + sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif


c------------------------------------------------------------------|
       endif

c----------------------------------------------------------------|
c          Block 4                                               |
      IF(chunk_type.eq.4) then
c----------------------------------------------------------------|     

      shift_r = 0
      shift_c = 3*ncis


c---------------F_so_x part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) - sqq*x1(i,j)*(-1.d0)
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
             x2(addr_r,addr_c)=x2(addr_r,addr_c) + sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif

c---------------F_so_y part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.10) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) - sqq*x1(i,j)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|
      IF(RL.eq.20) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
             x2(addr_r,addr_c)=x2(addr_r,addr_c) + sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif

c----------------------------------------------------------------|
c          Block 7                                               |
      IF(chunk_type.eq.7) then
c----------------------------------------------------------------|     

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|


      shift_r = ncis
      shift_c = 2*ncis


c---------------F_so_x part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) + sqq*x1(i,j)*(-1.d0)
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
             x2(addr_r,addr_c)=x2(addr_r,addr_c) + sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif

c---------------F_so_y part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.10) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) - sqq*x1(i,j)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|
      IF(RL.eq.20) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
             x2(addr_r,addr_c)=x2(addr_r,addr_c) - sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif

c----------------------------------------------------------------|
c          Block 8                                               |
      IF(chunk_type.eq.8) then
c----------------------------------------------------------------|     

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|


      shift_r = ncis
      shift_c = 3*ncis


c---------------F_so_x part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) + sqq*x1(i,j)*(-1.d0)
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
             x2(addr_r,addr_c)=x2(addr_r,addr_c) + sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif

c---------------F_so_y part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.10) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) + sqq*x1(i,j)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|
      IF(RL.eq.20) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
             x2(addr_r,addr_c)=x2(addr_r,addr_c) + sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif

c----------------------------------------------------------------|
c          Block 9                                               |
      IF(chunk_type.eq.9) then
c----------------------------------------------------------------|     

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|


      shift_r = 2*ncis
      shift_c = 0

c---------------F_so_x part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) + sqq*x1(i,j)*(-1.d0)
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
             x2(addr_r,addr_c)=x2(addr_r,addr_c) - sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif

c---------------F_so_y part-----------------------------------|

c-----------G_ij - part----------------------------------|

      IF(RL.eq.10) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) + sqq*x1(i,j)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|
      IF(RL.eq.20) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
             x2(addr_r,addr_c)=x2(addr_r,addr_c) - sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif

c----------------------------------------------------------------|
c          Block 10                                               |
      IF(chunk_type.eq.10) then
c----------------------------------------------------------------|     

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|


      shift_r = 2*ncis
      shift_c = ncis

c---------------F_so_x part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) + sqq*x1(i,j)*(-1.d0)
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
             x2(addr_r,addr_c)=x2(addr_r,addr_c) + sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif

c---------------F_so_y part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.10) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) + sqq*x1(i,j)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|
      IF(RL.eq.20) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
             x2(addr_r,addr_c)=x2(addr_r,addr_c) + sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif

c----------------------------------------------------------------|
c          Block 13                                              |
      IF(chunk_type.eq.13) then
c----------------------------------------------------------------|     

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|


      shift_r = 3*ncis
      shift_c = 0

c---------------F_so_x part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) - sqq*x1(i,j)*(-1.d0)
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
             x2(addr_r,addr_c)=x2(addr_r,addr_c) + sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif

c---------------F_so_y part-----------------------------------|

c-----------G_ij - part----------------------------------|

      IF(RL.eq.10) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) + sqq*x1(i,j)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|
      IF(RL.eq.20) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
             x2(addr_r,addr_c)=x2(addr_r,addr_c) - sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c------------------------------------------------------------------|
       endif

c----------------------------------------------------------------|
c          Block 14                                              |
      IF(chunk_type.eq.14) then
c----------------------------------------------------------------|     

c---------------------------------------------------------------------|
c       p = -ih d/dx, so we need to multiply matrix elements by -1.0  |
c---------------------------------------------------------------------|


      shift_r = 3*ncis
      shift_c = ncis

c---------------F_so_x part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.1) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) + sqq*x1(i,j)*(-1.d0)
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
             x2(addr_r,addr_c)=x2(addr_r,addr_c) + sqq*x1(a,b)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif

c---------------F_so_y part-----------------------------------|

c-----------G_ij - part----------------------------------|
      IF(RL.eq.10) then

       DO i=a1,a2
        DO a=1,nvirta
           addr_r = nvirta*(i-1) + a + shift_r
           DO j=b1,b2
            DO b=1,nvirta
               addr_c = nvirta*(j-1) + b + shift_c
             IF(a.eq.b) then
             x2(addr_r,addr_c) = x2(addr_r,addr_c) - sqq*x1(i,j)*(-1.d0)
             ENDIF
            ENDDO
           ENDDO
        ENDDO
       ENDDO
      endif
c-----------G_ab - part----------------------------------|
      IF(RL.eq.20) then

       DO i=1,nocca
        DO a=a1,a2
        aa = a - nocca
           addr_r = nvirta*(i-1) + aa + shift_r
           DO j=1,nocca
            DO b=b1,b2
            bb = b - nocca
               addr_c = nvirta*(j-1) + bb + shift_c
             IF(i.eq.j) then
             x2(addr_r,addr_c)=x2(addr_r,addr_c) - sqq*x1(a,b)*(-1.d0)
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


