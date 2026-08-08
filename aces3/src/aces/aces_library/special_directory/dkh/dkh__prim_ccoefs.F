      Subroutine Dkh__prim_ccoefs(Work, Ixpcoef, Npfps, Nshells,
     &                            Npcoefs)
      
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Npcoefs), Ixpcoef(Nshells), Npfps(Nshells)

      Index = 0
      Do ishell=1, Nshells
         ishel_off = Ixpcoef(ishell)
         iprim_off = npfps(ishell)
         index = 0
         Do i=1, iprim_off
            work(ishel_off+index) = 1.0D0
            Index = Index + iprim_off + 1
         Enddo
      Enddo

      Return
      End
