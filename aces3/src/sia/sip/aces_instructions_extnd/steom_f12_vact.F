c---------------------------------------------------------------|
c  This subroutine is a superinstruction for the SIAL code      |
c  It process F12 corrections to STEOM intermediates            |
c  Added by Denis Bokhan, Moscow Lomonosov State University,    |
c              Moscow, RUSSIAN FEDERATION                       |
c---------------------------------------------------------------|

      subroutine steom_f12_vact(x,nindex, type, bval,eval, bdim, edim,
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
        call do_steom_f12_2ind( x, bval(1),eval(1),type(1),
     *                     bval(2),eval(2),type(2),
     *                x2,bval2(1),eval2(1),type2(1),
     *                   bval2(2),eval2(2),type2(2))
        endif

        if(nindex.eq.4) then
        call do_steom_f12_4ind( x, bval(1),eval(1),type(1),
     *                     bval(2),eval(2),type(2),
     *                     bval(3),eval(3),type(3),
     *                     bval(4),eval(4),type(4),
     *                x2,bval2(1),eval2(1),type2(1),
     *                   bval2(2),eval2(2),type2(2),
     *                   bval2(3),eval2(3),type2(3),
     *                   bval2(4),eval2(4),type2(4) )
        endif


       return
       end

c-------------------------------------------------------------------
c-------------------------------------------------------------------

       subroutine do_steom_f12_2ind(x1,a1,a2,type1,b1,b2,type2,
     & x2,c1,c2,type3, d1,d2,type4 )
      implicit none
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'daveom.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2
      integer i,j,k,l,c,a,b,aa,bb,d,dd
      integer nvirta,nocca 
      integer type1,type2,type3,type4,type5,type6,type7
      double precision x1(a1:a2,b1:b2),x2(c1:c2,d1:d2)

      if(me.eq.0) write(6,*) a1,a2,'==========>',b1,b2,'----->'      
      if(me.eq.0) write(6,*) d1,d2,'^^^^^^^^^^>',c1,c2,'@@@@@'

       nvirta = Nalpha_virtual 
       nocca  = Nalpha_occupied 

       IF(RL.eq.1) then

       DO b=b1,b2
       DO a=a1,a2
       DO d=d1,d2
       dd = d + nocca
       IF(a.eq.dd) then
       if(me.eq.0) write(6,*) a,b,d,'<<<<<',dd
       x2(b,d) = x1(a,b)
       ENDIF
       ENDDO
       ENDDO
       ENDDO

       ENDIF

       return
       end
c--------------------------------------------------------|
c--------------------------------------------------------|

       subroutine do_steom_f12_4ind(x1,a1,a2,type1,b1,b2,type2,
     & g1,g2,type3,f1,f2,type4,x2,i1,i2,type5,j1,j2,type6,k1,k2,type7,
     & l1,l2,type8 )
      implicit none
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'daveom.h'

      integer a1,a2,b1,b2,c1,c2,d1,d2,i1,i2,j1,j2,k1,k2,g1,g2,f1,f2
      integer i,j,k,l,c,a,b,aa,bb,l1,l2,f,g,ii,ll
      integer nvirta,nocca 
      integer type1,type2,type3,type4,type5,type6,type7,type8
      double precision x1(a1:a2,b1:b2,g1:g2,f1:f2),
     & x2(i1:i2,j1:j2,k1:k2,l1:l2)

       nvirta = Nalpha_virtual 
       nocca  = Nalpha_occupied 

      IF(RL.eq.2) then

      DO a=a1,a2
      DO b=b1,b2
      DO f=f1,f2
      DO g=g1,g2 
      DO i=i1,i2
      ii = i + nocca
      IF(g.eq.ii) then
      x2(i,f,a,b) = x1(a,b,g,f)  
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO 


      ENDIF

      IF(RL.eq.3) then

      DO a=a1,a2
      DO b=b1,b2
      DO f=f1,f2
      DO g=g1,g2
      DO l=l1,l2
      ll = l + nocca
      IF(g.eq.ll) then
      x2(a,b,f,l) = x1(a,b,g,f)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO      

      ENDIF 


      IF(RL.eq.4) then
      DO b=b1,b2
      DO f=f1,f2
      DO g=g1,g2
      DO a=a1,a2
      DO i=i1,i2
      ii = i + nocca
      IF(a.eq.ii) then
c      if(me.eq.0) write(6,*) a,b,g,f,'<<<<<',ii
      x2(i,b,g,f) = x1(a,b,g,f)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
  


      ENDIF

      return
      end
c----------------------------------------------------------------------
