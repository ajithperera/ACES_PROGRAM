










      Subroutine  Pccd_dropmo_fix(Grad,Grad_t,Nbas,Ndrop)

      Implicit Double Precision(A-H,O-Z)
      Dimension Grad(Nbas,NBas), Idrop(Ndrop), Grad_t(Nbas,Nbas)

c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end


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



c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Data Azero /0.0D0/

      Call Getrec(20,"JOBARC","MODROPA",Ndrop,Idrop)

      Noc_drp = 0
      Nvr_drp = 0

      Do I = 1, Ndrop
         If (Idrop(i) .Le. Nocco(1)) Then
             Noc_drp = Noc_drp + 1
         Else
             Nvr_drp = Nvr_drp +1 
         Endif
      Enddo 

C We always assume that the first Noc_drp orbitals and last Nbas-Nvr_drp
C orbitals are dropped. Zero out the corresponding blocks of the
C gradient matrix

       Do I = 1, Noc_drp
          Do J = 1, Nbas 
             Grad(J,I) = Azero
          Enddo
       Enddo  

       Do I = 1, Nbas
          Do J = 1, Noc_drp 
             Grad(J,I) = Azero
          Enddo
       Enddo 


       Nstart = Nbas-Nvr_drp+1
       Do I = Nstart, Nbas
          Do J = 1, Nbas 
             Grad(J,I) = Azero
          Enddo
       Enddo  

       Do I = 1, Nbas
          Do J = Nstart, Nbas 
             Grad(J,I) = Azero
          Enddo
       Enddo  

       Call Dzero(Grad_t,Nbas*Nbas)

       K = 0
       L = 0
       Nstart = Nbas-Nvr_drp
       Do I = Noc_drp+1, Nstart
             K = K + 1
             L = 0
          Do J = Noc_drp+1, Nstart
             L = L + 1
             Grad_t(L,K) =Grad(J,I)
             Print*, "J and I", J,I,Grad(J,I)
             Print*, "K and L", K,L,Grad_t(L,K)
          Enddo
       ENddo
       Print*,"Test1"
       Call output(Grad,1,Nbas,1,Nbas,Nbas,Nbas,1)
   
       Ndrop_mo = Nbas - Ndrop 
       Call output(Grad_t,1,Ndrop_mo,1,Ndrop_mo,nbas,nbas1)
       Call Dcopy(Ndrop_mo*Ndrop_mo,Grad_t,1,Grad,1)
       Print*,"Test2"
       Call output(Grad,1,Ndrop_mo,1,Ndrop_mo,Nbas,Nbas,1)
      
       Return 
       End 
