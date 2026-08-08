











      subroutine runchg
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      double precision     dTimeStart
      common /aces2_times/ dTimeStart
      save   /aces2_times/

      integer iflags(100),iflags2(500)
      common /flags/ iflags
      common /flags2/ iflags2

      logical
     &        mrcc,
     &        geom_opt, analytical_gradient,numerical_gradient,
     &        raman, vib_specs,Cmpt_props_numrcl,
     &        fno, hf_scf, ks_scf, sewswitch, plain_scf, direct,
     &        nddo_guess, dirmp2, mopac_guess, bexport, 
     &        Single_point, Single_point_energy,	
     &        Single_point_gradient, First_order_props,
     &        Higher_order_props, NMR_SHIFTS, NMR_SPNSPN,
     &        NLO_PROPS, TDHF, bHyper, bMolden,blCCSDT,
     &        BEextrap, BGextrap, BCextrap, Analytical_hessian,
     &        hf_dft,sl_oep,Manual_FD,dkh,contract 
      common /sv_bool/
     &        mrcc,
     &        geom_opt, analytical_gradient,numerical_gradient,
     &        raman, vib_specs,Cmpt_props_numrcl,
     &        fno, hf_scf, ks_scf, sewswitch, plain_scf, direct,
     &        nddo_guess, dirmp2, mopac_guess, bexport,
     &        Single_point, Single_point_energy,	
     &        Single_point_gradient, First_order_props,
     &        Higher_order_props, NMR_SHIFTS, NMR_SPNSPN,
     &        NLO_PROPS, TDHF, bHyper, bMolden,blCCSDT,
     &        BEextrap, BGextrap, BCextrap, Analytical_hessian,
     &        hf_dft,sl_oep,Manual_FD,dkh,contract
      save   /sv_bool/

      character
     &          xjoda*79,
     &          integral_package*79,
     &          der_integral_package*79
      common /sv_char/
     &          xjoda,
     &          integral_package,
     &          der_integral_package
      save   /sv_char/


c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      DIMENSION ITOTREC(5),ITOTWRD(5)
      DIMENSION MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500),
     &          MOIODS(10,500),MOIOFL(10,500)
      DIMENSION NOIO(10,500),NOIOWD(10,500),NOIOSZ(10,500),
     &          NOIODS(10,500),NOIOFL(10,500)
      DIMENSION ISYTYP(2,500),isyty2(2,500)


      call a2getrec(20,'JOBARC','TOTRECMO',5,ITOTREC)
      call a2getrec(20,'JOBARC','TOTWRDMO',5,ITOTWRD)
      call a2getrec(20,'JOBARC','NDROTVRT',1,ndrvrtt)

      if (fno) call resetfno(eval,nbas,'DM')

      call runit('mv JOBARC JOBARC_DM')
      call runit('mv JAINDX JAINDX_DM')
      call runit('mv JOBARC_AM JOBARC')
      call runit('mv JAINDX_AM JAINDX')

      ndrgeo=2

      call a2putrec(20,'JOBARC','TOTRECMO',5,ITOTREC)
      call a2putrec(20,'JOBARC','TOTWRDMO',5,ITOTWRD)
      call a2putrec(20,'JOBARC','NDROPGEO',1,NDRGEO)

      call runit('xvtran')
      call a2putrec(20,'JOBARC','NDROTVRT',1,ndrvrtt)
      call runit('xintprc')
    
      if (fno) call resetfno(eval,nbas,'AM')

      MOIOSIZ=5000

      CALL A2GETREC(20,'JOBARC','MOIOVEC',MOIOSIZ,MOIO)
      CALL A2GETREC(20,'JOBARC','MOIOWRD',MOIOSIZ,MOIOWD)
      CALL A2GETREC(20,'JOBARC','MOIOSIZ',MOIOSIZ,MOIOSZ)
      CALL A2GETREC(20,'JOBARC','MOIODIS',MOIOSIZ,MOIODS)
      CALL A2GETREC(20,'JOBARC','MOIOFIL',MOIOSIZ,MOIOFL)
      CALL A2GETREC(20,'JOBARC','ISYMTYP',1000,ISYTYP)

      call runit('mv JOBARC JOBARC_AM')
      call runit('mv JAINDX JAINDX_AM')
      call runit('mv JOBARC_DM JOBARC')
      call runit('mv JAINDX_DM JAINDX')

      CALL A2GETREC(20,'JOBARC','MOIOVEC',MOIOSIZ,NOIO)
      CALL A2GETREC(20,'JOBARC','MOIOWRD',MOIOSIZ,NOIOWD)
      CALL A2GETREC(20,'JOBARC','MOIOSIZ',MOIOSIZ,NOIOSZ)
      CALL A2GETREC(20,'JOBARC','MOIODIS',MOIOSIZ,NOIODS)
      CALL A2GETREC(20,'JOBARC','MOIOFIL',MOIOSIZ,NOIOFL)
      CALL A2GETREC(20,'JOBARC','ISYMTYP',1000,ISYTY2)

      if (ndrvrtt.eq.0) then
         istart=291
      else
         istart=201
      end if
      do i=istart,400
         isyty2(1,i)=isytyp(1,i)
         isyty2(2,i)=isytyp(2,i)
         do k=1,10
            noio(k,i)  =moio(k,i)
            noiowd(k,i)=moiowd(k,i)
            noiosz(k,i)=moiosz(k,i)
            noiods(k,i)=moiods(k,i)
            noiofl(k,i)=moiofl(k,i)
         end do
      end do
      do i=istart,400
         isytyp(1,i)=0
         isytyp(2,i)=0
         do k=1,10
            moio(k,i)  =0
            moiowd(k,i)=0
            moiosz(k,i)=0
            moiods(k,i)=0
            moiofl(k,i)=0
         end do
      end do

      CALL A2PUTREC(20,'JOBARC','MOIOVEC',MOIOSIZ,NOIO)
      CALL A2PUTREC(20,'JOBARC','MOIOWRD',MOIOSIZ,NOIOWD)
      CALL A2PUTREC(20,'JOBARC','MOIOSIZ',MOIOSIZ,NOIOSZ)
      CALL A2PUTREC(20,'JOBARC','MOIODIS',MOIOSIZ,NOIODS)
      CALL A2PUTREC(20,'JOBARC','MOIOFIL',MOIOSIZ,NOIOFL)
      CALL A2PUTREC(20,'JOBARC','ISYMTYP',1000,ISYTY2)

      call runit('mv JOBARC JOBARC_DM')
      call runit('mv JAINDX JAINDX_DM')
      call runit('mv JOBARC_AM JOBARC')
      call runit('mv JAINDX_AM JAINDX')

      CALL A2PUTREC(20,'JOBARC','MOIOVEC',MOIOSIZ,MOIO)
      CALL A2PUTREC(20,'JOBARC','MOIOWRD',MOIOSIZ,MOIOWD)
      CALL A2PUTREC(20,'JOBARC','MOIOSIZ',MOIOSIZ,MOIOSZ)
      CALL A2PUTREC(20,'JOBARC','MOIODIS',MOIOSIZ,MOIODS)
      CALL A2PUTREC(20,'JOBARC','MOIOFIL',MOIOSIZ,MOIOFL)
      CALL A2PUTREC(20,'JOBARC','ISYMTYP',1000,ISYTYP)

      call runit('mv JOBARC JOBARC_AM')
      call runit('mv JAINDX JAINDX_AM')
      call runit('mv JOBARC_DM JOBARC')
      call runit('mv JAINDX_DM JAINDX')


      return
      end

