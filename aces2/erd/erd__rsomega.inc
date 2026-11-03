C------------------------------------------------------------------------
C  INCLUDE FILE: ERD__RSOMEGA
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  DESCRIPTION : Range-separation parameter of the two-electron
C                operator. RSOMEG = 0 (the default of a common block)
C                gives the ordinary 1/r12 integrals; RSOMEG = omega > 0
C                gives the long-range integrals erf(omega*r12)/r12
C                that range-separated hybrid functionals need. Only
C                xerdints sets it (see erdints/erd_driver.F).
C
C------------------------------------------------------------------------
         DOUBLE PRECISION  RSOMEG

         COMMON  /ERDRSO/  RSOMEG
