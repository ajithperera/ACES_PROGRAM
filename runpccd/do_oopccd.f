













































































































































































































      Subroutine Do_oopccd(Iuhf)

      integer, intent(in)::Iuhf
      Logical Ca,Fc,Bca,Bfc,Occnums_file,Move
      Character*80 Fname

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      Data Izz /100/

      print*,'inside job control',Iflags2((178))
      print*,'Iflags for oo:', Iflags2(178)
      print*,'true statement:',Iflags2(178).eq.3

      If (Iflags2(178) .EQ. 1 .or.
     +         Iflags2(178) .EQ. 2  .or.
     +         Iflags2(178).EQ.3    .or.
!SteepestDescent
     +         Iflags2(178).EQ.4) then

         print*,'Running steepest descent (SD) algo to optimize...'
         call do_FullNR(Work,Maxcor,Iuhf)
         !Call do_steepestDescent(Work,Maxcor,Iuhf)
!      Else if (Iflags(178) .EQ. 2) then !L-BFGS
!          print*,'Inside L-BFGS'
!         Call do_LBFGS()
!      Else if (Iflags2(178).EQ.3) then !NR step
!          print*,'Running full NR algo to optimize...'
!         call do_FullNR(Work,Maxcor,Iuhf)
      Else
         print*,'zero off-diag T2'
         Call Runit("xvmol")
         Call Runit("xvmol2ja")
         Call Runit("xvscf")
         Call Runit("xvtran")
         Call Runit("xintprc")
         Call Runit("xpccd")
      Endif 

      Return
      End 


     
