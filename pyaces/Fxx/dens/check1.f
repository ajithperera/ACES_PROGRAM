     

 
      SUBROUTINE CHECK1(ICORE,A,LENGTH,LISTF,ISPIN,FACT,iflag,
     &                  num)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD
      double precision icore
      DIMENSION ICORE(1),A(length),num(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/ADD/SUM
      common/info/nocco(2),nvrto(2)
      E=0.0D+0
      call zero(icore,length)
c      CALL GETLST(ICORE,1,1,1,Ispin,LISTf)
      if(ispin.eq.1) then
      call getrec(20,'JOBARC','SCFEVALA',iintfp*(nocco(1)+nvrto(1)),
     &             icore(1+length))
      else
      call getrec(20,'JOBARC','SCFEVALB',iintfp*(nocco(1)+nvrto(1)),
     &             icore(1+length))
      endif
      ioff1=0
      ioff2=0
      if(iflag.eq.1) ioff2=nocco(ispin) 
      do 100 irrep=1,nirrep
       n=num(irrep)
       do 10 i=1,n 
        icore(ioff1+(i-1)*N+i) = icore(length+ioff2+i)
        write(*,*) icore(length+ioff2+i)
10     continue
       ioff1=ioff1+N*N
       ioff2=ioff2+N
100   continue
   
      E=E+SDOT(length,ICORE,1,a,1)
1000  CONTINUE
      write(6,*) 'listf',listf
      write(6,*) 'energy contribution',fact*e
      sum=sum+FACT*e
      write(6,*) ' cumulated e',sum
      return
      end
