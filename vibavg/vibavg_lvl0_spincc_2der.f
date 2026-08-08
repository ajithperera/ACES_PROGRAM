










      Subroutine Vibavg_lvl0_spincc_2der(Prop_derv2,Exponent,Nreals,
     +                                  Nmodes,Ndata,Avg)

      Implicit Double Precision(A-H,O-Z)

      Dimension Prop_derv2(Nreals,Nreals,Nmodes)
      Dimension Exponent(Nmodes,Ndata)
      Dimension Avg(Nreals,Nreals,Ndata)

      Data Done,Factor14,Factor18/1.0D0,0.25D0,0.125D0/

      Do K = 1, Ndata 
         Do J = 1, Nreals 
            Do I = Nreals, J-1, -1
               Do Imode = 1, Nmodes 
                  Avg(I,J,K) = Avg(I,J,K) + Prop_derv2(I,J,Imode)*
     +                         Exponent(Imode,K)
               Enddo
            Enddo
         Enddo 
      Enddo 


      Call Dscal(Ndata*Nreals*Nreals,Factor14,Avg,1) 

      Write(6,*)
      Write(6,"(2a)") " Zeroth-order vibrational and thermal", 
     +                " contribution to NMR spin-spin coupling"
      do i=1,Ndata
      call output(Avg(1,1,i),1,Nreals,1,Nreals,Nreals,Nreals,1)
      enddo
      Return
      End 

     
