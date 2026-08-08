










      subroutine Tdee_reset_vars_4restart(Component,Irrepx,Naobfns,
     +                                    Nbfns,Ipert,Work_Label,
     +                                    D_pole,Q_pole)

      Implicit Integer (A-Z)
      Character*8 Label_D(3),Work_Label 
      Character*8 Label_Q(6)
      Logical D_pole,Q_pole

      Common /LISTDENS/LDENS

c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      Data Label_D /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      Data Label_Q /'QUAD_XX ','QUAD_YY ','QUAD_ZZ ','QUAD_XY ',
     +              'QUAD_XZ ','QUAD_YZ '/
      Data Ione /1/

      Nbfns   = Nocco(1) + Nvrto(1)
      Call Getrec(20, "JOBARC","NBASTOT ", Ione, Naobfns)

      Ldens = 160
      If (D_pole) Then
         Work_Label = Label_D(Component)
      Elseif (Q_pole) Then 
         Work_Label = Label_Q(Component)
      Endif 

      Ipert = Component 
     
      Return
      End
