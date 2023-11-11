










      Subroutine Form_new_t1pamps(T1amps,Work,Iwork,Maxcor,Imaxcor,
     &                            List1,T1ln,T1ln_aa,Ipert,Irrepx,Iuhf)
      
      Implicit Double Precision(A-H,O-Z)
      Character*3 Type
      Integer T1ln,T1off,T1ln_aa
      Integer Aend
      Logical CCt1
      
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end


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



c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c active_space.com : begin
      Parameter(Max_xp=100)
      integer Active_oo(8,2),Active_vv(8,2)
      integer Pactive_oo(100,8,2),Pactive_vv(100,8,2)
      integer Ioff_active_oo(8,2),Ioff_active_vv(8,2)
      integer Pioff_active_oo(100,8,2),Pioff_active_vv(100,8,2)
      Double Precision Oo_threshold,Vv_threshold
      common /actvsp_info/Active_oo,Pactive_oo,Active_vv,Pactive_vv,
     +                    Ioff_active_oo,Pioff_active_oo,
     +                    Ioff_active_vv,Pioff_active_vv,
     +                    Oo_threshold,Poo_threshold,
     +                    Vv_threshold,Pvv_threshold,Eta_val(Max_xp),
     +                    E_k(Max_xp),E_ks(Max_xp)
c active_space.com: end


      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension T1amps(T1ln)
      Dimension Work(Maxcor)

      Data Ione,Inull /1,0/

      Do Ispin =1, 1+Iuhf

         Ioff   = Inull
         Ndim   = Irpdpd(Irrepx,8+Ispin)
         T1off  = (Ispin-1)*Ndim + Ione

         I000 = Ione
         I010 = I000 + Ndim 
         Iend = I010 + Ndim 

         If (Iend .Gt. Maxcor) Call Insmem("transform_2index",Iend,
     +                                      Maxcor)

C EOM-CCSD Tx on the disk while PEOM-CCSD is in memory.

         Call Dcopy(Ndim,T1amps(T1off),1,Work(I010),1)
         Call Getlst(Work(I000),1,1,1,Ispin,List1)

         do Irrep = 1, Nirrep
            Irrepr = Irrep
            Irrepl = Dirprd(Irrepr,Irrepx) 
            
             Nrow = Vrt(Irrepl,Ispin)
             Ncol = Pop(Irrepr,Ispin)

             Aend = Pactive_vv(Ipert,Irrepl,Ispin)
             Iend = Pactive_oo(Ipert,Irrepr,Ispin)

             Call Form(Work(I000+Ioff),Work(I010+Ioff),Iwork,
     +                      Imaxcor,Nrow,Ncol,Aend,Iend)

             Ioff = Ioff + Nrow*Ncol
             
         Enddo

         Call Dcopy(Ndim,Work(I010),1,T1amps,1)

         Call Putlst(T1amps,1,1,1,Ispin,List1)
      Enddo 

      Return
      End
            
