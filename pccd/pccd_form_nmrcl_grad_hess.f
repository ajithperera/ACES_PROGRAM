










      Subroutine Pccd_form_nmrcl_grad_hess(Work,Maxcor,Iuhf,Igrad_calc,
     +                                     IHess_calc,Scale)
      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)

c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end


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




      Nbasis = Nocco(1) + Nvrto(1)
      Isqlen = Nbasis*Nbasis
      Ndis = 2
      Iplus = Ione
      Imins = Iplus + Isqlen
      Iderv = Imins + Isqlen
      Iend  = Iderv + Isqlen

      Do Idis = 1, Ndis
C  Here we need to Call orbital rotation driver, This code is not 
C complete. Needs lot od work to get numerical gradient and Hessians. 
         
         If (Idis .Eq. Itwo) Scale = -Scale
         Call Pccd_form_htau_pq(Work,Maxcor,Iuhf,Igrad_calc,
     +                          IHess_calc,Scale)
         Call Getrec(20,"JOBARC","OBRTGRDA",Isqlen*IIntfp,Work(Iplus))
         If (Idis .Eq. Itwo) Call Getrec(20,"JOBARC","OBRTGRDA",
     +                                   Isqlen*IIntfp,Work(Imins))

         Dinv = Half*(One/Scale)
         Do I = 1, Isqlen
            Work(Iderv-1+I) = (Work(Iplus-1+I) - Work(Imins-1+I))*Dinv
         Enddo
      Enddo 

      Return 
      End
