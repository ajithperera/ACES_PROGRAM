










      Subroutine Do_lanczos_main(Vecs,Work,Memleft,Iuhf,Irrepx,Nsize,
     +                           Maxmem)

      Implicit Double Precision (A-H,O-Z)


 
C Maxn : Lead dimension of the tridiagonal matrix
C Info : Serve as both input/output. Only Info(1) is 
C        set by the caller. 
C Maxvw: This is set to 1. This allocate one extra vector 
C        of the diemension of Hbar (this space is not used).
C Convrg: Strict or light 

      Integer Maxn
      Parameter(Maxn_int=10000,Maxvw_int=1)
      Character*6 Convrg 

      Common /lanczos_vars/Maxn,Maxvw,Trd(Maxn_int,Maxn_int),
     +                     Convrg 
      




      Dimension Work(Memleft),Vecs(Nsize,2*Maxvw+6)

        Write(6,*)
        Write(6,"(a)") "---Entered Do_lanczos_main ---"
        Write(6,*)

      Call  Lanczos_main(Vecs,Work,Memleft,Nsize,Irrepx,Iuhf,
     +                   Maxmem)

        Write(6,*)
        Write(6,"(a)") "---Exit Do_lanczos_main ---"
        Write(6,*)

      Return
      End

      

