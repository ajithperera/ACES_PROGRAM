










      Subroutine Vibavg_lvl0_shifts_2der(Prop_derv2,Exponent,Nreals,
     +                                   Nmodes,Ndata,Avg)

      Implicit Double Precision(A-H,O-Z)

      Dimension Prop_derv2(Nreals,Nmodes)
      Dimension Cfc_contr(Nmodes)
      Dimension Exponent(Nmodes)
      Dimension Avg(Nreals,Ndata)

      Data Done,Factor14,Factor18/1.0D0,0.25D0,0.125D0/
      Data Done/1.0D0/

     
      Do K = 1, Ndata
         Do J = 1, Nreals
            Do Imode = 1, Nmodes
               Avg(J,1) = Avg(J,1) + Prop_derv2(J,Imode)*
     +                    Exponent(Imode) 
            Enddo
         Enddo
      Enddo 


      Call Dscal(Ndata*Nreals,Factor14,Avg,1)

      Write(6,*)
      Write(6,"(2a)") " Zeroth-order vibrational and thermal",
     +                " contribution to NMR shielding constant"
      do i=1,Ndata
      write(6,"(5(1x,F15.8))") (Avg(j,i),j=1,Nreals)
      enddo
      Return
      End 
     
