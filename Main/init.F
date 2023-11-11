      Subroutine Init()

      Implicit None

      include "machsp.com"
      include "jobarc.com"
      include "icore.com"
      include "istart.com"
      include "flags.com"
      include "flags2.com"
      include "lists.com"
      include "filspc.com"
      include "auxcache.com"

CSSS      Write(6,"(a)") "Entering aces2_main"

c ----------------------------------------------------------------------
CSSS      If (Flag .Eq. 1) Then
      Call Aces2_main()
CSSS      Endif
CSSS      If (Flag .Eq. 0) Call aces_fin
C
c ----------------------------------------------------------------------
CSSS      Write(6,"(a)") "Exiting aces2_main"
CSSS      call test_includes()
CSSS      call test_core(Icore(i0),Icrsiz)

      Return
      End

