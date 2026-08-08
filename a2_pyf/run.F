      Subroutine Run(String)

      Implicit None

      Integer Iuhf 
#include "flags.com"
     
      Character*(*) String

      if (iflags(11).eq.0) then
         iUHF = 0
      else
         iUHF = 1
      end if

      If (String(1:4) .Eq. "ints") Then
         Call Runints(Iuhf)
      Elseif(String(1:3) .Eq. "scf") Then 
         Call Runscf(Iuhf)
      Elseif(String(1:2) .Eq. "cc") Then
         Call Runcc(Iuhf)
      Elseif(String(1:2) .Eq. "ee") Then
         Call Runee(Iuhf)
      Elseif(String(1:4) .Eq. "hbar") Then
         Call Runhbar(Iuhf)
      Endif

      Return
      End


