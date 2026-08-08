






































































































































































































      Subroutine Get_2ndord_t1s(T2amp,Work,Maxcor,Tln,T2ln,T1ln,
     +                          T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,
     +                          T1ln_bb,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Integer T2ln,T2ln_aa,T2ln_bb,T2ln_ab
      Integer T1ln,T1ln_aa,T1ln_bb
      Integer Tln
      Integer T1off
      Dimension Work(Maxcor),T2amp(Tln)

      Data Ione /1/



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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      Ioff_t1 = T2ln + Ione
      Ioff_t2 = Ione

      Call Dzero(T2amp(Ioff_t1),T1ln_aa)
      If (Iuhf .Eq. 0) Call Dzero(T2amp(Ioff_t1+T1ln_aa),T1ln_bb)

      If (Iuhf .Eq. 1) Then 
         Call T2t1_aa1(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                 Pop(1,1),Vrt(1,1),1,T2ln,T1ln,T2ln_aa,
     +                 T2ln_bb,T2ln_ab,T1ln_aa,T1ln_bb)

         Call T2t1_aa1(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                 Pop(1,2),Vrt(1,2),2,T2ln,T1ln,T2ln_aa,
     +                 T2ln_bb,T2ln_ab,T1ln_aa,T1ln_bb)

      Endif 

      Call T2t1_ab1(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +              Pop(1,1),Pop(1,2),Vrt(1,1),Vrt(1,2),1,Iuhf,
     +              T2ln,T1ln,T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,
     +              T1ln_bb)

      If (Iuhf .Eq. 1) Then
          Call T2t1_ab1(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                  Pop(1,2),Pop(1,1),Vrt(1,2),Vrt(1,1),2,Iuhf,
     +                  T2ln,T1ln,T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,
     +                  T1ln_bb)
      Endif 

      If (Iuhf .Eq. 1) Then 

         Call T2t1_aa2(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                 Pop(1,1),Vrt(1,1),1,T2ln,T1ln,T2ln_aa,
     +                 T2ln_bb,T2ln_ab,T1ln_aa,T1ln_bb)

         Call T2t1_aa2(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                 Pop(1,2),Vrt(1,2),2,T2ln,T1ln,T2ln_aa,
     +                 T2ln_bb,T2ln_ab,T1ln_aa,T1ln_bb)

      Endif

      Call T2t1_ab2(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +              Pop(1,1),Pop(1,2),Vrt(1,1),Vrt(1,2),1,Iuhf,
     +              T2ln,T1ln,T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,
     +              T1ln_bb)

      If (Iuhf .Eq. 1) Then
          Call T2t1_ab2(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                  Pop(1,2),Pop(1,1),Vrt(1,2),Vrt(1,1),2,Iuhf,
     +                  T2ln,T1ln,T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,
     +                  T1ln_bb)
      Endif

      Call checksum("T1AA :",T2amp(t2ln+Ione),t1ln_aa)
      If (Iuhf .Ne.0) Call checksum("T1BB :",T2amp(T2ln+T1ln_aa+Ione),
     +                               t1ln_bb)

      Do Ispin = 1, Iuhf+1
         T1off = Ione + T2ln + (Ispin-1)*T1ln_aa 
         Call Form_t1(T2amp(T1off),Work,Maxcor,Pop(1,Ispin),
     +                 Vrt(1,Ispin),Ispin,T1ln,T1ln_aa)
      Enddo 

      Call checksum("T1AA :",T2amp(t2ln+Ione),t1ln_aa)
      If (Iuhf .Ne.0) Call checksum("T1BB :",T2amp(T2ln+T1ln_aa+Ione),
     +                               t1ln_bb)
   
      Return
      End 

      
      
