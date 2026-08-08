      Subroutine Init_legacies()

      Implicit Integer (A-Z)

c   o from the old MCHPRM routine
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /BUFLEN/ ILNBUF
      LOGICAL MINTPRC,MVDINT1,MVDINT2,MANTI
      COMMON /MINTPC/ MINTPRC
      COMMON /MVDINT/ MVDINT1,MVDINT2
      COMMON /MEMANTI/ MANTI
c   o from the old GETSTF routine
      LOGICAL LINCC,CICALC
      COMMON /LINEAR/ LINCC,CICALC

c ----------------------------------------------------------------------

c   o initialize legacy data

c     /BUFLEN/ (only used in intprc)
      ILNBUF = 1536

c     /FILES/
      LUOUT  = 6
      MOINTS = 50

c     /LINEAR/
      LINCC  = .FALSE.
      CICALC = .FALSE.

#if defined(_UNICOS) || ( defined(_AIX) && defined(F_64BIT) )
      MINTPRC = .FALSE.
      MVDINT2 = .FALSE.
      MANTI   = .FALSE.
#else
      MINTPRC = .TRUE.
      MVDINT2 = .TRUE.
      MANTI   = .TRUE.
#endif /* _UNICOS || (_AIX && F_64BIT) */

c ----------------------------------------------------------------------

      RETURN
      END


