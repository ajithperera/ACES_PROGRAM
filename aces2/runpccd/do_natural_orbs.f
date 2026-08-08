











      subroutine do_natural_orbs
      implicit integer (a-z)
      integer ndropa,icalc,iabcd,igamma_abcd
      logical NOS_Exist
     



































































































































































































c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
      
      ndrgeo1=0
      call aces_ja_init
      call getrec(20,'JOBARC','IFLAGS  ',100,iflags)
      icalc=iflags(2)
      iflags(2)=1
      iflags(19)=1
      iflags(16)=0
      iflags(45)=4
      iflags2(138)=1
      iabcd=iflags(93)
      iflags(93)=0
      ivtran=iflags(83)
      iflags(83)=0
      igamma_abcd=iflags(100)
      iflags(100)=0
      call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
      call putrec(20,'JOBARC','IFLAGS2 ',500,iflags2)
      call getrec(20,'JOBARC','NUMDROPA',1,ndropa)
      call putrec(20,'JOBARC','FNOFREEZ',1,ndropa)
      call putrec(20,'JOBARC','NUMDROPA',1,0)
      if (iflags(11).ne.0) then
         call putrec(20,'JOBARC','NUMDROPB',1,0)
      end if
      if ((geom_opt.or.vib_specs).and.analytical_gradient) then
         print *,'@RUNFNO FNO gradients not yet supported.'
         call aces_exit(1)
         call putrec(20,'JOBARC','NDROPGEO',1,0)
      endif
      call aces_ja_fin
      call runit('xvtran')
      call runit('xintprc')
      call runit('xvcc')
      call runit('xdens')

      call aces_ja_init
      call getrec(20,'JOBARC','IFLAGS  ',100,iflags)
      call getrec(20,'JOBARC','IFLAGS2 ',500,iflags2)
      iflags(2)=icalc
      iflags(93)=iabcd
      iflags(83)=ivtran
      iflags(100)=igamma_abcd
      call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
      call putrec(20,'JOBARC','IFLAGS2 ',500,iflags2)
      call aces_ja_fin
      inquire(file='NATORBS',exist=NOS_Exist)
      if (NOS_Exist) Then
         call runit('xvmol')
         call runit('xvmol2ja')
         call runit('mv NATORBS OLDMOS')
         call runit('xvscf')
	 call runit('xvtran')
         call runit('xintprc')
         call runit('xvcc')
      else
         write(6,"(2a)") " The pCCD natural orbital optimization", 
     &                   "requires NATORBS file and it does not exist"
         call errex
      endif 
      return
      end

