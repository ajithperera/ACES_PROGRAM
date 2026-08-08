        subroutine parread(iuhf)

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c
c read in PCCSD parameters
c
        Character*1 A, B
	open(unit=1,file='parcc',
     &     form='formatted',access='sequential',status='old',err=999)

        ispar   = .false.
        coulomb = .false.

        read(1, "(a,1x,a)") A, B

        If (A .eq. "T" .and. B .eq. "T") then
	   ispar   = .true.
           coulomb = .true.
        Else if (A .eq. "T" .and. B .eq. "F") then
	   ispar   = .true.
           coulomb = .false.
        Endif 

	read(1,*) paralpha
        read(1,*) parbeta
        read(1,*) pargamma
        read(1,*) pardelta
c Originally this parameter is called pargamma.
        read(1,*) parepsilon
        close(1)
        return
999	ispar= .false.
        return
        end
