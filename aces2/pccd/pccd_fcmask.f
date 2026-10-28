










      Subroutine Pccd_fcmask(Work,Maxcor,List)
C
C Frozen (dropped) MOs with OO_METHOD=RFO (see runpccd/pccd_fc_setup):
C zero every element T2(Ab,Ij) of list List (RHF, A-b/I-j full
C distributions) with a frozen index: I or j among the first NFZO
C occupied MOs of its irrep, A or b among the last NFZV virtual MOs of
C its irrep. Applied to the guess T2 and to each new increment, this
C keeps the frozen pairs out of the pCCD amplitudes (and Lambda, which
C starts from T2), so that the full-space densities, orbital gradient
C and Hessian are those of the active-space pCCD energy with the
C frozen MOs in the reference. No-op without the PCCDNFZO record or
C when nothing is frozen.
C
      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz
      Dimension Work(Maxcor),Nfzo(8),Nfzv(8)


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

      Call Getrec(0,"JOBARC","PCCDNFZO",Length,Junk)
      If (Length .Le. 0) Return
      Call Getrec(20,"JOBARC","PCCDNFZO",8,Nfzo)
      Call Getrec(20,"JOBARC","PCCDNFZV",8,Nfzv)
      Nfz = 0
      Do Irr = 1, Nirrep
         Nfz = Nfz + Nfzo(Irr) + Nfzv(Irr)
      Enddo
      If (Nfz .Eq. 0) Return
C Layout as in Pccd_reset: rows A-b (type 15), columns I-j (type 14).

      Do Irrep = 1, Nirrep
         Numdis = Irpdpd(Irrep,14)
         Dissiz = Irpdpd(Irrep,15)
         If (Numdis*Dissiz .Eq. 0) Go to 10
         Ncol = 0
         Nrow = 0
         Do Irb = 1, Nirrep
            Ncol = Ncol + Pop(Irb,1)*Pop(Dirprd(Irb,Irrep),1)
            Nrow = Nrow + Vrt(Irb,1)*Vrt(Dirprd(Irb,Irrep),1)
         Enddo
         If (Ncol .Ne. Numdis .Or. Nrow .Ne. Dissiz) Then
            Write(6,"(a,i4)") " @-Pccd_fcmask: layout mismatch, list",
     +                        List
            Call Errex
         Endif
         If (Numdis*Dissiz .Gt. Maxcor) Call Insmem("Pccd_fcmask",
     +                                   Numdis*Dissiz,Maxcor)
         Call Getlst(Work,1,Numdis,1,Irrep,List)
C Columns: distributions Ij, I fastest within the (irrep I, irrep j)
C block, blocks ordered by the irrep of j.
         Icol = 0
         Do Irj = 1, Nirrep
            Iri = Dirprd(Irj,Irrep)
            Do J = 1, Pop(Irj,1)
               Do I = 1, Pop(Iri,1)
                  Icol = Icol + 1
                  If (I .Le. Nfzo(Iri) .Or. J .Le. Nfzo(Irj))
     +               Call Dzero(Work((Icol-1)*Dissiz+1),Dissiz)
               Enddo
            Enddo
         Enddo
C Rows: Ab, A fastest, blocks ordered by the irrep of b.
         Irow = 0
         Do Irb = 1, Nirrep
            Ira = Dirprd(Irb,Irrep)
            Do Ib = 1, Vrt(Irb,1)
               Do Ia = 1, Vrt(Ira,1)
                  Irow = Irow + 1
                  If (Ia .Gt. Vrt(Ira,1)-Nfzv(Ira) .Or.
     +                Ib .Gt. Vrt(Irb,1)-Nfzv(Irb)) Then
                     Do K = 1, Numdis
                        Work((K-1)*Dissiz+Irow) = 0.0D0
                     Enddo
                  Endif
               Enddo
            Enddo
         Enddo
         Call Putlst(Work,1,Numdis,1,Irrep,List)
 10      Continue
      Enddo
      Return
      End
