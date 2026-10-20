














































































































































































































      Subroutine Pccd_finalize(Work,Maxcor,Iuhf,Tol,Grd_max,Grd_rms)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Character*6 Status 
      Logical Pccd,Ccd,Lccd

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


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



      Common /Calc/Pccd,Ccd,Lccd

      Write(6,*) 
      Status ="Finish"  
      Call Pccd_banner(6,Status)
C SCFENEG0/TOTENEG0 are written by do_brueckner after the first macro
C iteration; converged in that iteration (OO_METHOD=RFO restart at a
C minimum, exact-Hessian check) the starting orbitals are the final ones.
      Call Getrec(0,'JOBARC','SCFENEG0',Length,Junk)
      If (Length .Gt. 0) Then
         Call Getrec(20,'JOBARC','SCFENEG0',Iintfp,Eref_I)
         Call Getrec(20,'JOBARC','TOTENEG0',Iintfp,Ecor_I)
      Else
         Call Getrec(20,'JOBARC','SCFENEG ',Iintfp,Eref_I)
         Call Getrec(20,'JOBARC','TOTENERG',Iintfp,Ecor_I)
      Endif
  

      Call Getrec(20,'JOBARC','SCFENEG ',Iintfp,Eref_F)
      Call Getrec(20,'JOBARC','TOTENERG',Iintfp,Ecor_F)

      write(6,"(2a,F15.8)") " The SCF energy of the starting",
     +                     " orbitals : ", Eref_I
      If (Pccd) Then
      write(6,"(2a,F15.8)") " The pCCD energy of the starting",
     +                     " orbitals: ", Ecor_I
      Elseif(Ccd) Then
      write(6,"(2a,F15.8)") " The CCD energy of the starting",
     +                     " orbitals : ", Ecor_I
      Elseif (Lccd) Then
      write(6,"(2a,F15.8)") " The LCCD energy of the starting",
     +                      " orbitals: ", Ecor_I
      Endif 

      Write(6,*) 

      If (Iflags2(178) .Eq. 6) Then
      Write(6,"(2a)") " The orbital rotation algorithm :",
     +                " RFO (trust-region, SR1)"
      Else
      Write(6,"(3a)") " The orbital rotation algorithm :", " Newton",
     +                "-Raphson with Diis" 
      Endif
      Write(6,*) 

      Write(6,"(a,E12.6)")  " The maximum gradient at convergence: ",
     +                      Grd_max
      Write(6,"(a,E12.6)")  " The RMS gradient at convergence    : ",
     +                      Grd_rms
      Write(6,"(a,E12.6)")  " The selected convergence tolarence : ",
     +                      Tol

      Write(6,*) 

      If (Pccd) Then
      write(6,"(2a,1x,F15.8)") " The orbital optimized pCCD energy ",
     +                      "     :",Ecor_f
      Elseif(Ccd) Then
      write(6,"(2a,1x,F15.8)") " The orbital optimized CCD energy ",
     +                      "      :", Ecor_f
      Elseif (Lccd) Then
      write(6,"(2a,1x,F15.8)") " The orbital optimized LCCD energy ",
     +                      "     :", Ecor_f
      Endif 

      Return
      End
