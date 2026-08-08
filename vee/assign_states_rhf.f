










      Subroutine Assign_states_rhf(Coefs,Nature,Bgn,End,Bgn_irp,
     &                             End_irp,Iuhf,Irrepx,Iroot)

      Implicit Double Precision (A-H, O-Z) 

      Dimension Coefs(*)
      Integer Bgn(100,8),End(100,8),NBfirr(8),Bgn_irp(100,8),
     &        End_irp(100,8)
      
      Character*1 Nature(100,8)
      Logical Singlet, Triplet



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
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

      Triplet = .False.
      Singlet = .False.

      Thres1  = 1.0D-05
      Thres2  = 0.05D0
      Ncount = 1 
      Call Getlst(Coefs, 1, 1, 1, 1, 490)

      Nature(Iroot,Irrepx) = "S"
      Ncount = 1
      Weight_max = Thres2
      Do Irepr = 1, Nirrep
         irepl = Dirprd(Irepr,Irrepx)
         Do I = 1, Pop(irepr,1)
            Do J = 1, vrt(irepl,1)
               Weight1    = Dabs(Coefs(Ncount))
               Save       = Weight1
               Weight_max = max(Weight1,Weight_max) 

               If (Dabs(Save-Weight_max) .eq. 0.0D0) Then
                   Bgn(Iroot,Irrepx)       = i
                   Bgn_irp(Iroot,Irrepx)   = irepr
                   End(Iroot,Irrepx)       = j
                   End_irp(Iroot,Irrepx)   = irepl
               Endif 
               Ncount = Ncount + 1 
            Enddo
         Enddo
      Enddo

     
      Return
      End
