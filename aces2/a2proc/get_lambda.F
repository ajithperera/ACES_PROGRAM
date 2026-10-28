      Subroutine get_lambda(Grad_K1CM, Vec_K1C, Vec_K1C_On_Evecs, 
     &                      Grad_K1C_On_Evecs, Hess_K1C, Evecs, 
     &                      AtmMass, Vec_K1C_Updated, Work, 
     &                      Coords_K1C, Stride, Delta, Nreals,
     &                      Hess_Eval_Tol, Lamsearch_Tol, 
     &                      Binsearch_Tol)

      Implicit Double Precision (A-H, O-Z)
C
      Double Precision L_bound,  Lamsearch_Tol, Mass, Lamda
C
      Logical lamda_found, Bracket_end
      Parameter (Sym_Tol = 1.0D-06)
C
      Dimension Grad_K1CM(3*Nreals), Hess_K1C(3*Nreals, 3*Nreals), 
     &          AtmMass(Nreals), Evecs(3*Nreals, 3*Nreals),
     &          Vec_K1C(3*Nreals), Vec_K1C_On_Evecs(3*Nreals),
     &          Grad_K1C_On_Evecs(3*Nreals), Work(27*Nreals*Nreals),
     &          Vec_K1C_Updated(3*Nreals), Coords_K1C(3*Nreals)
C
C Mass weigh the Hessian and obtain the eigenvectors.
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,*) "The Hessian"
      Call output(Hess_K1C, 1, 3*Nreals, 1, 3*Nreals, 3*Nreals, 
     &            3*Nreals, 1)
#endif
      Call Dcopy(9*Nreals*Nreals, Hess_K1C, 1,
     &                            Work(18*nreals*nreals + 1), 1)
C
      Do Ideg = 1, 3*Nreals
         Do Jdeg = 1, 3*Nreals  
            Mass = DSQRT(AtmMass(1+(Ideg-1)/3)*AtmMass(1+(Jdeg-1)/3))
            IF (Mass .lt. 1.0D-3) Then
               Hess_K1C(Ideg, Jdeg) = 0.0D0
            Else
               Hess_K1C(Ideg, Jdeg) = Hess_K1C(Ideg, Jdeg)/Mass
            Endif
         Enddo
      Enddo
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,*) "The mass weighted Hessian"
      Call output(Hess_K1C, 1, 3*Nreals, 1, 3*Nreals, 3*Nreals, 
     &            3*Nreals, 1)
#endif
C
C Project the mass-weighted Hessian and diagonalize to define normal
C modes.

      Call Projec_FC(Coords_K1C, Hess_K1C, AtmMass, Grad_K1CM, Work, 
     &               WorK(9*Nreals*Nreals+1), 1.0D-8, Nreals, .True., 
     &               .True., .False.)
C
#ifdef _DEBUG_LVLM1
      NX = 3*Nreals
      Write(6,"(a)") "The projected Hessian"
      CALL OUTPUT(Hess_K1C, 1, NX, 1, Nx, Nx, Nx, 1)
#endif

      Call Eig(Hess_K1C, Evecs, 3*Nreals, 3*Nreals, 0)
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      NX = 3*Nreals
      Write(6,"(a)") "The eigen vectors of the projected Hessian"
      CALL OUTPUT(Evecs, 1, NX, 1, Nx, Nx, Nx, 1)
      Write(6,"(a)") "The eigenvalues of the  projected Hessian"
      Write(6, "(4F17.13)") (Hess_K1C(I,I), I=1, NX)
#endif
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,*) "@-get_lambda Grad_K1CM (ZMAT order)"
      Write(6, "(3F17.13)") (Grad_K1CM(i), i=1,3*Nreals)
      Write(6,*)
      Write(6,*) "@-get_lambda VEC_K1C (ZMAT order)"
      Write(6, "(3F17.13)") (Vec_K1C(i), i=1,3*Nreals)
#endif
C
C Gradient along the eigenvectors of the Hessian
C   
      Do Imode = 1, 3*Nreals
         Grad_K1C_on_Evecs(Imode) = Ddot(3*Nreals, Grad_K1CM, 1, 
     &                                  Evecs(1, Imode), 1)
         Vec_K1C_on_Evecs(Imode)  = Ddot(3*Nreals, Vec_K1C, 1,
     &                                  Evecs(1, Imode), 1)
C
         If (Dabs(Hess_K1C(Imode, Imode)) .LT. Hess_Eval_Tol) Then
             Grad_K1C_on_Evecs(Imode) = 0.0D0
             Vec_K1C_on_Evecs(Imode)  = 0.0D0
             Hess_K1C(Imode, Imode)   = 0.0D0
         Endif 
C
C Modes with neither gradient nor displacement are those forbidden by
C the symmetry of the path; leave them out like the zero modes. Near
C a minimum such a mode can have the lowest eigenvalue, lambda then
C sits just below it and 1/(h-lambda) blew round-off up into a real
C out-of-plane step (HCOH IRC, 10/2026).
C
         If (Dabs(Grad_K1C_on_Evecs(Imode)) .LT. Sym_Tol .AND.
     &       Dabs(Vec_K1C_on_Evecs(Imode))  .LT. Sym_Tol) Then
             Grad_K1C_on_Evecs(Imode) = 0.0D0
             Vec_K1C_on_Evecs(Imode)  = 0.0D0
             Hess_K1C(Imode, Imode)   = 0.0D0
         Endif
C
      Enddo
C 
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,*) "The position and gradient vectors along e-vecs"
      Write(6, "(F17.13)") (Grad_K1C_on_EvecS(I), I=1, 3*Nreals)
      Write(6,*)
      Write(6, "(F17.13)") (Vec_K1C_on_EvecS(I), I=1, 3*Nreals)
#endif
C
C The lowest nonzero eigenvalue (the zeroed ones are the translations,
C rotations and modes below Hess_Eval_Tol). The old Do While ran with
C .OR. Imode .GT. 0 and walked off the array: an endless loop. 10/2026.
C
      Eval_min = 0.0D0
      Do Imode = 1, 3*Nreals
         If (Hess_K1C(Imode, Imode) .NE. 0.0D0 .AND.
     &      (Eval_min .EQ. 0.0D0 .OR.
     &       Hess_K1C(Imode, Imode) .LT. Eval_min))
     &       Eval_min = Hess_K1C(Imode, Imode)
      Enddo
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,*) "The lowest eiegnvalue of hess. and the tolerance"
      Write(6, "(a,F17.13,a,F13.10)") "Eval_min =", Eval_min, 
     &          " Hess_Eval_Tol =",Hess_Eval_Tol 
#endif
C
C The solution to Eq. 26 im JPC, vol. 94, 5525, 1990 must be
C lower than the lowst eigenvalue of the Hessian. Another words
C Lambda < Eval_min. 
C 
      halfs      = 0.50D0*Stride
      L_bound    = Eval_min - (2.0D0)*Delta 
      U_bound    = Eval_min - (1/2.0D0)*Delta 
      Halfs2     = halfs*halfs
      niter      = 1
      Lbound_Tol = max((Eval_min)*Lamsearch_Tol,
     &                   10.0D0*Lamsearch_Tol)
      Ubound_Tol = max(1.0D03, Abs(Eval_min))
      Ubound_Tol = 1.0D03*Ubound_Tol
C
#ifdef _DEBUG_LVLM1
        Write(6,*) "Lambda iterations" 
        Write(6,*) "Starting lower and upper limits"
        Write(6,"(2F17.13)") L_bound, U_bound
#endif
C
#ifdef _DEBUG_LVL0
      Write(6,"(a,ES14.6,a,ES12.4,a,ES12.4)")
     &      " @-get_lambda: Eval_min=", Eval_min, " Delta=", Delta,
     &      " Stride=", Stride
      Write(6,"(a)") " @-get_lambda: eigenvalues, grad, vec on e-vecs"
      Write(6,"(3ES16.7)") (Hess_K1C(I,I), Grad_K1C_on_Evecs(I),
     &                      Vec_K1C_on_Evecs(I), I=1,3*Nreals)
#endif
      Bracket_end = .FALSE. 
      Do while (.NOT. Bracket_end) 
C The bracket must be found well before the bounds run away; stop
C with a message instead of looping (it used to hang). 10/2026.
         If (niter .GT. 100) Then
            Write(6,"(2a,ES14.6,2(a,ES12.4))") " @-get_lambda: no",
     &            " bracket for lambda; Eval_min=", Eval_min,
     &            " F_Lowr=", F_Lowr, " F_Uppr=", F_Uppr
            Call Errex
         Endif
         F_Lowr = 0.0D0
         F_Uppr = 0.0D0
C
         Do Imode = 1, 3*Nreals     
            Dnumer = Hess_K1C(Imode, Imode)*Vec_K1C_on_Evecs(Imode) - 
     &               Grad_K1C_on_Evecs(Imode)
            Denomi = Hess_K1C(Imode, Imode) - L_bound
            F_Lowr = Dnumer**2/Denomi**2 + F_Lowr 
            Denomi = Hess_K1C(Imode, Imode) - U_bound 
            F_Uppr = Dnumer**2/Denomi**2 + F_Uppr
         Enddo
#ifdef _DEBUG_LVLM1
        Write(6,*) "The F_Lowr and F_Uppr"
        Write(6,"(2F17.13)") F_Lowr, F_Uppr
#endif
C
         F_Lowr = F_lowr - Halfs2
         F_Uppr = F_Uppr - Halfs2
         If (F_Lowr*F_Uppr .LT. 0.0D0) Then 
             Bracket_end = .TRUE.
             Go to 10
         Endif
C
         niter   = niter  + 1
         L_bound = Eval_min - (2.0D0)**niter*Delta
         U_bound = Eval_min - (1.0D0/(2.0D0)**niter)*Delta
c
#ifdef _DEBUG_LVLM1
        Write(6,*) "Lower and upper limits during its."
        Write(6,"(2F17.13)") L_bound, U_bound
#endif
         
         If (L_bound .LT. -Ubound_Tol .AND. ABS(Eval_min-U_bound) .LT. 
     &       Lbound_Tol)  Then 
             L_bound = -Ubound_Tol
             U_bound = Eval_min - Lbound_Tol
             Write(6, "(a,a)") "@-get_lambda Bracketing for lambda", 
     &             "is failed"

             Call Errex
         Endif
      Enddo
C
   10 Continue
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The upper and lower bound for lambda"
      Write(6, "(a,F17.13,a,F17.13)") " L_bound =",  L_bound,
     &          " U_bound  =",U_bound
#endif
C
C The binary search portion to extract the value in the range 
C L_bound - U_bound.
C
      Tmp_lamda    = 0.0D0
      Lamda_found = .False.
      Niter        = 0
      
      Do while (.Not. lamda_found) 
         F_Lowr = 0.0D0
         F_Uppr = 0.0D0
         F_Midl = 0.0D0
         Lamda  = 0.50D0*(L_bound + U_bound) 
#ifdef _DEBUG_LVLM1
       Write(6,*) "The starting lambda"
       Write(6, "(F17.13)") lamda
#endif
         Do Imode = 1, 3*Nreals
            Dnumer = Hess_k1C(Imode, Imode)*Vec_K1C_on_Evecs(Imode) -
     &               Grad_K1C_on_Evecs(Imode)
            Denomi = Hess_K1C(Imode, Imode) - L_bound
            F_Lowr = Dnumer**2/Denomi**2 + F_Lowr
            Denomi = Hess_K1C(Imode, Imode) - U_bound
            F_Uppr = Dnumer**2/Denomi**2 + F_Uppr
            Denomi = Hess_K1C(Imode, Imode) - Lamda       
            F_Midl = Dnumer**2/Denomi**2 + F_Midl
         Enddo
#ifdef _DEBUG_LVLM1
        Write(6,*) "The FL, FU, FM during iterations"
        Write(6,"(3F17.13)") F_Lowr, F_Uppr, F_Midl
#endif 

            F_Lowr = F_lowr - Halfs2
            F_Uppr = F_Uppr - Halfs2
            F_Midl = F_midl - Halfs2       
C
            If (Abs(Tmp_lamda -lamda) .LT. Binsearch_Tol) 
     &           lamda_found = .True.
            
            Niter = Niter + 1
            If (Niter .gt. 1000) Then
                Write(6, "(a,a)") "@-Get_lambda: The binary search",
     &                            " for lambda failed!"
                Call Errex
            Endif 
C
            Tmp_lamda = Lamda
C
#ifdef _DEBUG_LVLM1
        Write(6,*) "The FL, FU, FM post halfs2"
        Write(6,"(3F17.13)") F_Lowr, F_Uppr, F_Midl
#endif
            if (F_Lowr*F_Midl .LT. 0.0D0) U_bound = Lamda
            If (F_Uppr*F_Midl .LT. 0.0D0) L_bound = Lamda
C
#ifdef _DEBUG_LVLM1
        Write(6,*) "The L_bound and U_bound during iterations"
        Write(6,"(2F17.13)") L_bound, U_bound
#endif
       Enddo
       
       If (Lamda .GT. Eval_min) Then 
           Write(6, "(a,a,2F10.5)")
     &    "@-Get_lambda: Warning, lambda is grater than the lowest",
     &    "eigenvalue", lamda, Eval_min
          Call Errex
       Endif
       
       If (Abs(Lamda - Eval_min) .LT. Binsearch_Tol) Then
          Write(6, "(a,a,F10.5)") "@-Get_lambda: Lamda is too close",
     &                            "to the lowest eigenvalue", Lamda,
     &                             Eval_min 
          Call Errex
       Endif 
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The value of lambda after binary search"
      Write(6, "(a,F17.13)") " Lambda  =", Lamda 
#endif
C
C Make the Newton-Raphson Step
C
       Call Dzero(Vec_K1C_Updated, 3*Nreals)
       Do Imode = 1, 3*Nreals
          If (Hess_K1C(Imode, Imode) .EQ. 0.0D0) Then
             Do Jmode = 1, 3*Nreals 
                Vec_K1C_Updated(Jmode) = Vec_K1C_Updated(Jmode) + 
     &                                   Hess_K1C(Imode, Imode)*
     &                                   Evecs(Jmode,Imode)
             Enddo
          Else
             Denomi = (Lamda*Vec_K1C_on_Evecs(Imode) - 
     &                 Grad_K1C_on_Evecs(Imode))
             Dnumer = Hess_K1C(Imode, Imode) - Lamda
             Do Jmode = 1, 3*Nreals
                Vec_K1C_Updated(Jmode) = Vec_K1C_Updated(Jmode) + 
     &                                   (Denomi/Dnumer)*
     &                                   Evecs(Jmode,Imode)
  
             Enddo 
          Endif
       Enddo
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The M. W. updated vector (Vec_K1C_Updated)"
      Write(6,"(3F17.13)") (Vec_K1C_Updated(i), i=1,3*Nreals) 
#endif
C
C Eleminate mass weighing from the updated vector.
C
      Ioff = 0
      Do Iatom = 1, Nreals
         Do Ixyz = 1, 3
            Ioff = Ioff + 1
            Vec_K1C_updated(Ioff) = Vec_K1C_updated(Ioff)/
     &                              Dsqrt(AtmMass(Iatom))
         Enddo 
      Enddo    
C 
      Call Dcopy(9*Nreals*Nreals, Work(18*Nreals*Nreals+1), 1,
     &           Hess_K1C, 1)
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6,*) "The restored Hessian"
      Call output(Hess_K1C, 1, 3*Nreals, 1, 3*Nreals, 3*Nreals,
     &            3*Nreals, 1)
#endif                          
C
       Return
       End
             
