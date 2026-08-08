













































































































































































































      Subroutine Tdcc_form_mutilde_dot_d(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)

      Double Precision Mone, Mut_dot 
      Dimension Ioff_vo(8,2)
 
      Data Zero, One, Mone, Two /0.0D0, 1.0D0, -1.0D0, 2.0D0/



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
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
c symloc.com : begin
c The numbering scheme is as follows:
c   a<b    (alpha) [1]
c   a<b    (beta)  [2]
c   i<j    (alpha) [3]
c   i<j    (beta)  [4]
c   a<=b   (alpha) [5]
c   a<=b   (beta)  [6]
c   i<=j   (alpha) [7]
c   i<=j   (beta)  [8]
c   a,i    (alpha) [9]
c   a,i    (beta)  [10]
c   a,i    (AB)    [11]
c   a,i    (BA)    [12]
c   a,b    (AB)    [13]
c   i,j    (AB)    [14]
c   a,b    (AB)    [15]
c   i,a    (alpha) [16]
c   i,a    (beta)  [17]
c   i,a    (AB)    [18]
c   a,b    (alpha) [19]
c   a,b    (beta)  [20]
c   i,j    (alpha) [21]
c   i,j    (beta)  [22]
c   a,b    (BA)    [23]
c   i,j    (BA)    [24]
c   i,a    (BA)    [25]
      integer         isymoff(8,8,25)
      common /symloc/ isymoff
c symloc.com : end

      Call Tdcc_load_mutilde_vo(Work,Memleft,Memleft_modf,Iuhf,
     +                          Ioff_vo,Irrepx)

C NDep paper Eqn. 29 first term(^t indicate transpose)
C
C UHF: Z(AB,IJ) = P(AB) {Mu_t(B,M)^t * Hbar(AM,IJ)^t}
C    : Z(ab,ij) = P(ab) {Mu_t(b,m)^t * Hbar(am,ij)^t}
C    : Z(aB,iJ) = P(aB) {Mu_t(B,M)^t * Hbar(aM,iJ)^t}

      Call Tdcc_form_mutilde_dot_da(Work,Memleft_modf,Irrepx,
     +                              Iuhf,Ioff_vo)
      call tdcc_mutilde_dot_da_debug(Work,Memleft,Iuhf,Irrepx)

C NDep paper Eqn. 29 second term
C
C UHF: Z(AB,IJ) = P(IJ) {Mu_t(E,J)^t * Hbar(AB,IE)^t}
C    : Z(ab,ij) = P(ab) {Mu_t(e,j)^t * Hbar(AB,ie)^t}
C    : Z(aB,iJ) = P(aB) {Mu_t(E,J)^t * Hbar(aB,iJ)^t}
   
      Call Tdcc_form_mutilde_dot_db(Work,Memleft_modf,Irrepx,
     +                              Iuhf,Ioff_vo)
      call tdcc_mutilde_dot_db_debug(Work,Memleft,Iuhf,Irrepx)

C NDep paper Eqn. 29 third term.
C Z(AB,IJ) = P(AB)P(IJ) Mu_t(I,A) F(B,J)

      If (Iuhf .Ne. 0) Then
         Call Tdcc_form_mutilde_dot_dc_uhf(Work,Memleft,Irrepx)
     +                                     
      Else
         Call Tdcc_form_mutilde_dot_dc_rhf(Work,Memleft,Irrepx)
      Endif 

C NDep paper Eqn. 29 Forth-term.
C Z(AB,IJ) = P(AB)P(IJ) Mu_t(I,A) F(B,J)

      If (Irrepx .EQ. 1) Then
         Call Tdcc_form_mutilde_dot_dd(Work,Memleft,Irrepx,
     +                                 Iuhf,Mut_dot)
      Else
         Mut_dot = 0.0D0
      Endif 

      Return
      End
     
           
        
