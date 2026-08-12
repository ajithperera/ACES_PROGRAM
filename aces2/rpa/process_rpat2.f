










      Subroutine Process_rpat2(Work,Length,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)
      Character*4 Spin(2)
      Dimension Ioff_row(9),Ioff_col(9)



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











































































































































































































c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/




c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Data Ione,One /1,1.0D0/
      Data Spin /"AAAA","BBBB"/

C The T2 vectors generated in the RPA code are ordered in ai,bj form 
C regardless of the spin. We need to reorder them into ab,ij form to
C interface with the rest of the code. The RPA vs rCC identity holds
C only when both singlet and triplets are included. Therefore RHF
C runs are flagged.

      If (Iuhf .EQ. 0) Then
         Write(6,"(a)") " Both singlet and triplet roots are needed.",
     +                  " Therfore switch to UHF reference to proceed."
         Call Errex 
      Else

C First AAAA and BBBB spin cases.

         Do Ispin = 1, IUHF+1  

            Nsize_s = 0
            Nsize_t = 0
            Do Irrepx = 1, Nirrep
               Nrow_s = Irpdpd(Irrepx,8+Ispin)
               Ncol_s = Irpdpd(Irrepx,8+Ispin)
               Nrow_t = Irpdpd(Irrepx,18+Ispin)
               Ncol_t = Irpdpd(Irrepx,20+Ispin)
               Nsize_s = Nsize_s + Nrow_s*Ncol_s
               Nsize_t = Nsize_t + Nrow_t*Ncol_t
            Enddo 

            If (Nsize_s .ne. Nsize_t) Then
             Write(6,"(a,a,a)") " The dimensions of the same spin",
     +                          " incoming and outgoing arrays",
     +                          " do not match."
               Call Errex
            Endif 
            It2_s = Ione
            It2_t = It2_s + Nsize_s
            Iend  = It2_t
            If (Iend .Ge. Length) Call Insmem("Process_rpat2",
     +                                         Iend,Maxcor)

C Read in T2(ai,bj) order and obtain T2(ab,ij). Note that this T2
C is not antisymmetric. 

            Ioff = It2_s
            Do Irrepx = 1, Nirrep 
               Nrow = Irpdpd(Irrepx,8+Ispin)
               Ncol = Irpdpd(Irrepx,8+Ispin)
               Call Getlst(Work(Ioff),1,Ncol,1,Irrepx,199+Ispin)
               Ioff = Ioff + Nrow*Ncol
            Enddo 

            Iscr   = It2_t + Nsize_t
            Iend   = Iscr  + Nsize_t
            If (Iend .Ge. Length) Call Insmem("Process_rpat2",
     +                                         Iend,Maxcor)
            Call Sstgen(Work(It2_s),Work(It2_t),Nsize_t,Vrt(1,Ispin),
     +                  Pop(1,Ispin),Vrt(1,Ispin),Pop(1,Ispin),
     +                  Work(Iscr),1,"1324")
            Call Putall(Work(It2_t),Nsize_t,1,43+Ispin)
        Enddo 

C ABAB block

        Nsize_s = 0
        Nsize_t = 0
        Do Irrepx = 1, Nirrep
           Nrow_s = Irpdpd(Irrepx,9)
           Ncol_s = Irpdpd(Irrepx,10)
           Nrow_t = Irpdpd(Irrepx,15)
           Ncol_t = Irpdpd(Irrepx,14)
           Nsize_s = Nsize_s + Nrow_s*Ncol_s
           Nsize_t = Nsize_t + Nrow_t*Ncol_t
        Enddo 
        If (Nsize_s .ne. Nsize_t) Then
           Write(6,"(a,a,a)") " The dimensions of the opposite",
     +                        " spin incoming and outgoing arrays",
     +                        " do not match."
          Call Errex
        Endif 

        It2_s  = Ione
        It2_t  = It2_s + Nsize_s
        Iend   = It2_t
        If (Iend .Ge. Length) Call Insmem("Process_rpat2",
     +                                     Iend,Maxcor)

C Read in T2(AI,bj) order and obtain T2(Ab,Ij)

        Ioff = It2_s 
        Do Irrepx = 1, Nirrep
           Nrow = Irpdpd(Irrepx,9)
           Ncol = Irpdpd(Irrepx,10)
           Call Getlst(Work(Ioff),1,Ncol,1,Irrepx,199+Ispin)
           Ioff = Ioff + Nrow*Ncol
        Enddo 

        Iscr   = It2_t + Nsize_t
        Iend   = Iscr  + Nsize_t
        If (Iend .Ge. Length) Call Insmem("Process_rpat2",
     +                                     Iend,Maxcor)

        Call Sstgen(Work(It2_s),Work(It2_t),Nsize_t,Vrt(1,1),Pop(1,1),
     +              Vrt(1,2),Pop(1,2),Work(Iscr),1,"1324")
        Call Putall(Work(It2_t),Nsize_t,1,46)

      Endif 
C Add a record to JOBARC to flag that the RPA T2 vectors are available
C on the list 44,45 and 46.

      Call Putrec(20,"JOBARC","RPAT2VEC",Ione,One)

      Return
      End 
