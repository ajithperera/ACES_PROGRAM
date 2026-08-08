










      Subroutine Pdcc_fixf(Work,Maxcor,Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Maxcor)


      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 


      Call Dcc_hbar_fae(Work,Maxcor,Iuhf)
      Call Dcc_hbar_fmi(Work,Maxcor,Iuhf)
      Call Dcc_fixfbar(Work,Maxcor,Iuhf)

      REturn
      End 
