










      Subroutine Shannon_main(Work, Maxcor, Iuhf)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION Work(MAXCOR)
C
      COMMON /FLAGS/ IFLAGS(100)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)

      If (Iuhf.NE.0) Then
C
        It1off = 1
        Call Getlst(Work(It1off), 1, 1, 1, 1, 90)
        Nt1absiz = IRPDPD(1,9)
        It1off = It1off + Nt1absiz
        Call Getlst(Work(It1off), 1, 1, 1, 2, 90)
        It1off = It1off + Nt1absiz
        Nt1siz = IRPDPD(1,9) + IRPDPD(1,10)

        It2off = 1 + Nt1siz 
        Nt2absiz = Isymsz(1,3)
        Call Getall(Work(it2off), Nt2absiz, 1, 44)

        It2off = It2off + Nt2absiz
        Nt2absiz = Isymsz(2,4) 
        Call Getall(Work(it2off), Nt2absiz, 1, 45)
        It2off = It2off + Nt2absiz
        Nt2absiz =  Isymsz(13,14)
        Call Getall(Work(it2off), Nt2absiz, 1, 46)
 
        Nt2siz = Isymsz(1,3) + Isymsz(2,4) + Isymsz(13,14)

        Shannon = Shannon12(Nt1siz, Nt2siz, Work(It1off),
     &                      Work(It2off))

        Write(6,*)
        Write(6,"(a,F13.7)") " Shannon Index =", Shannon
        Write(6,*)
         
      Else

         It1off = 1
         Nt1siz = IRPDPD(1,9)
         Call Getlst(Work(It1off), 1, 1, 1, 1, 90)

         It2off = It1off + Nt1siz

         Nt2siz = Isymsz(13,14)
         Call Getall(Work(it2off), Nt2siz, 1, 46)

         Shannon= Shannon12(Nt1siz, Nt2siz, Work(It1off),
     &                     Work(It2off))
         Write(6,*)
         Write(6,"(a,F13.7)") " Shannon Index =", Shannon
         Write(6,*) 

      Endif

      RETURN
      END

    
