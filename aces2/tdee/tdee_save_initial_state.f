










      Subroutine Tdee_Save_initial_state(Work,Memleft,Irrepx,Iuhf,
     +                                   Iside,Nsize)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Memleft)
      Logical Rhf



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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
     
C Singles and doubles right and dipole functions list are saved. 

        Write(6,*)
        Write(6,"(a)") "---Entered Tdee_Save_initial_state ---"
        Write(6,*)
      I000 = 1
      I010 = I000 + NSIZE 
      I020 = I010 + NSIZE

      Ioffr1 = 0
      Ioffr2 = 0
      Ioffsp = 2
      Mu_s   = 390
      Mu_d_o = 313
      Mu_d_d = 316

      Call Tdee_load_vec(Irrepx,Work(I000),Memleft,Mu_s,Ioffr1,Mu_d_o,
     +                   Ioffr2,Iuhf,.False.)
      Call Tdee_normaliz_init_vec(Work(I000),Nsize)

      Ioffr1  = 0
      Ioffr2  = 0
      Ioffsp  = 2
      Mut_s   = 392
      Mut_d_o = 323
      Mut_d_d = 326

      Call Tdee_load_vec(Irrepx,Work(I010),Memleft,Mut_s,Ioffr1,Mut_d_o,
     +                  Ioffr2,Iuhf,.False.)
   
      RL_NORM = Ddot(Nsize,Work(I000),1,Work(I010),1)
     
      Call Dscal(Nsize,1.0D0/RL_NORM,Work(I010),1)

C#ifdef 1
      RL_NORM = Ddot(Nsize,Work(I000),1,Work(I010),1)
      Write(6,"(a,F15.7)") "@-Tdee_Save_initial_state,RL_NORM: ",
     +                      RL_NORM
C#endif 
      Call Tdee_dump_vec(Irrepx,Work(I000),Memleft,Mu_s,Ioffr1,Ioffsp,
     +                   Mu_d_d,Ioffr,Iuhf,.False.)
      Call Tdee_dump_vec(Irrepx,Work(I010),Memleft,Mut_s,Ioffr1,Ioffsp,
     +                   Mut_d_d,Ioffr,Iuhf,.False.)


        Write(6,*)
        Write(6,"(a)") "---Exit Tdee_Save_initial_state ---"
        Write(6,*)

      Return 
      End
