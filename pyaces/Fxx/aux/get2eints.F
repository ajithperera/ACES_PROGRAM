C*************************************************************
      subroutine Get2EInts(Label,G2,Nao)
C     
C     Get two electron integrals
C         
C*************************************************************
      implicit None
c
      Logical FileExist
      Character*4 filename(4)
      Character*8 Title(24)
      Character*8 Label 

#include "machsp.com"

      Double Precision G2(Nao,Nao,Nao,Nao)
      Double Precision buf(600)
      Double Precision Repuls 
      Double Precision X
      Integer Ibuf(600)
      Integer Nbfirr(8)
      Integer And
      Integer Or
      Integer Nut
      Integer Nao
      Integer Ifile
      Integer Iupki
      Integer Iupkj
      Integer Iupkk
      Integer Iupkl
      Integer Ipack
      Integer Ijunk
      Integer Int
      Integer N
      Integer Indx,Indx2,Indx3
      Integer Nirrep
      Integer Ix,Jx,Kx,Lx
      Integer I,J,K,L
      Integer Ind1,Ind2,Ind3
c
      iupki(int)=and(int,ialone)
      iupkj(int)=and(ishft(int,-ibitwd),ialone)
      iupkk(int)=and(ishft(int,-2*ibitwd),ialone)
      iupkl(int)=and(ishft(int,-3*ibitwd),ialone)
      ipack(i,j,k,l)=or(or(or(i,ishft(j,ibitwd)),ishft(k,2*ibitwd)),
     &                  ishft(l,3*ibitwd))

      indx(i,j)=j+(i*(i-1))/2
      indx3(i,j)=i+(j*(j-1))/2
      indx2(i,j,n)=i+(j-1)*n
c
      filename(1)='IIII'
      filename(2)='IJIJ'
      filename(3)='IIJJ'
      filename(4)='IJKL'
c
      do 10 ifile=1,4
         FileExist=.false.
         inquire(file=filename(ifile),exist=FileExist)

         if (FileExist) then
            open(unit=ifile,file=filename(ifile),form='UNFORMATTED',
     &           access='SEQUENTIAL')

            rewind ifile
            If (ifile .Eq. 1) Then
               Read(ifile) Title,Nirrep,(Nbfirr(i),i=1,Nirrep),Repuls,
     &                     Ijunk
            Endif 
c     
            call locate(ifile,Label)
            nut = 1
            do while (nut.gt.0)
               read(ifile) buf, ibuf, nut
               do int = 1, nut

                  ix=iupki(ibuf(int))
                  jx=iupkj(ibuf(int))
                  kx=iupkk(ibuf(int))
                  lx=iupkl(ibuf(int))
                  i=max(ix,jx)
                  j=min(ix,jx)
                  k=max(kx,lx)
                  l=min(kx,lx)
                  ind1=indx(i,j)
                  ind2=indx(k,l)
                  if(ind1.lt.ind2) then
                     ix=i
                     jx=j
                     i=k
                     j=l
                     k=ix
                     l=jx
                  endif
c     
                  x = buf(int)

                  G2(i,j,k,l)=x
                  G2(i,j,l,k)=x
                  G2(j,i,k,l)=x
                  G2(j,i,l,k)=x
                  G2(k,l,j,i)=x
                  G2(k,l,i,j)=x
                  G2(l,k,i,j)=x
                  G2(l,k,j,i)=x
c     
                  ind3=indx(ind1,ind2)
c     
               end do
            end do
            close(ifile)
         else
         end if   
 10   continue  

#ifdef _DEBUG_LVLM
      Do k = 1, Nao
         Do l = 1, Nao
            Do j = 1, Nao
               Do i = 1, Nao
      write(6,"(a,1x,4i4,1x,F15.9)") "I,J,L,K:", i,j,l,k,G2(i,j,l,k)
               Enddo
            Enddo
         Enddo
      Enddo
#endif
      return
      end
