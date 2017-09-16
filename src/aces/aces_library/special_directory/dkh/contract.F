      Subroutine Contract(Dkh_pints,Ccoef_ishell,Ccoef_jshell,Iprimfns,
     +                    Jprimfns,Icntrfns,Jcntrfns,Iprim_4shell,
     +                    Jprim_4shell,Icntr_4shell,Jcntr_4shell,
     +                    Tot_cnt_fns,Tot_prm_fns,Idegen,Jdegen,
     +                                              Dkh_cints)

      Implicit Double Precision (A-H, O-Z)

      Integer Tot_cnt_fns
      Dimension Dkh_pints(Idegen*Jdegen,Iprim_4shell*Jprim_4shell)
      Dimension Ccoef_ishell(Iprim_4shell,Icntr_4shell)
      Dimension Ccoef_jshell(Jprim_4shell,Jcntr_4shell)
      Dimension Dkh_cints(Tot_cnt_fns)

#ifdef _DEBUG_LVL0
      Write(6,*)
      write(6,"(1x,a,a)") "Printing from Contract: The primitve DKH",
     +           " ints per shell block."
      write(6,"(6(1x,E12.5))")((Dkh_pints(i,j),i=1,Idegen*Jdegen),j=1,
     +                         Iprim_4shell*Jprim_4shell)
      write(6,*) "The contraction coefs."
      write(6,"(6(1x,f10.6))")((CCoef_ishell(i,j),i=1,Iprim_4shell),
     +                          j=1,Icntr_4shell)
      write(6,*)
      write(6,"(6(1x,f10.6))")((CCoef_jshell(i,j),i=1,Jprim_4shell),
     +                          j=1,Jcntr_4shell)
      Write(6,*)
#endif 
      Call Dzero(Dkh_cints,Tot_cnt_fns)

      Int = 0

       Do Jcint = 1, Jcntr_4shell
          Do Icint = 1, Icntr_4shell

             Do J = 1, Jdegen
                Do I = 1, Idegen
                   II = (J - 1) * Idegen + I
                   Int = Int + 1

                   Do Jprm = 1, Jprim_4shell
                      Do Iprm = 1, Iprim_4shell
                         JJ = (Jprm - 1) * Iprim_4shell + Iprm 

                      Weight = Ccoef_ishell(Iprm,Icint) * 
     +                         Ccoef_jshell(Jprm,Jcint)

                      Dkh_cints(Int) = 
     +                                 Dkh_cints(Int) +
     +                                 Dkh_pints(II,JJ) *
     +                                 Weight 

#ifdef _DEBUG_LVL0
      write(6,"(1x,5i4,5F10.5)") II,JJ,Icint,jcint, int,
     +     Ccoef_ishell(Iprm,Icint),Ccoef_jshell(Jprm,Jcint),
     +     Weight,Dkh_pints(II,JJ),Dkh_cints(Int)
#endif 
                     Enddo
                  Enddo
C
                Enddo
            Enddo
C
         Enddo
      Enddo 
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,2i4)") "The contracted DKH ints per shell block"
      write(6,"(6(1x,E12.5))") (Dkh_cints(i),i=1,Tot_cnt_fns)
#endif 
      
      Return
      End 
     
