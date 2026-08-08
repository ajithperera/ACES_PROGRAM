










      Subroutine Check_phph(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)

      Double Precision Work(Maxcor)
      Logical UHF



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

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 


      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)

      IRREPX = 1

      If (Ispar .and. Coulomb) Then
      IF (UHF) Then
         AAAA_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,214),ISYTYP(2,214))
         BBBB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,215),ISYTYP(2,215))
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,216),ISYTYP(2,216))

         Call Getall(Work, AAAA_LENGTH_IJAB, IRREPX, 214)
         Call Checksum("AAAA_IJAB-1", Work, AAAA_LENGTH_IJAB)
         Call Getall(Work(AAAA_LENGTH_IJAB+1), AAAA_LENGTH_IJAB,
     &               IRREPX, 114)
         Call Daxpy(AAAA_LENGTH_IJAB, 1.0D0,Work(AAAA_LENGTH_IJAB+1),
     &              1,Work,1)
         Call Checksum("AAAA_IJAB-2", Work, AAAA_LENGTH_IJAB)

         Call Getall(Work, BBBB_LENGTH_IJAB, IRREPX, 215)
         Call Checksum("BBBB_IJAB-1", Work, BBBB_LENGTH_IJAB)
         Call Getall(Work(BBBB_LENGTH_IJAB+1), BBBB_LENGTH_IJAB,
     &               IRREPX, 114)
         Call Daxpy(BBBB_LENGTH_IJAB, 1.0D0,Work(BBBB_LENGTH_IJAB+1),
     &              1,Work,1)
         Call Checksum("BBBB_IJAB-2", Work, BBBB_LENGTH_IJAB)

         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 216)
         Call Checksum("ABAB_IJAB-1", Work, ABAB_LENGTH_IJAB)
         Call Getall(Work(ABAB_LENGTH_IJAB+1), ABAB_LENGTH_IJAB,
     &               IRREPX, 116)
         Call Daxpy(ABAB_LENGTH_IJAB, 1.0D0,Work(ABAB_LENGTH_IJAB+1),
     &              1,Work,1)
         Call Checksum("ABAB_IJAB-2", Work, ABAB_LENGTH_IJAB)
      Else
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,216),ISYTYP(2,216))
         Write(6,*) "ABAB_LENGTH_IJAB", ABAB_LENGTH_IJAB
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 216)
         Call Checksum("ABAB_IJAB", Work, ABAB_LENGTH_IJAB)
      Endif
      Endif 

      Return
      End 
