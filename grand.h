C---------------------------------------------------------------
      double precision VECCA(2000,100,12)
      double precision VECCB(2000,100,12)
      double precision VALA(20)
      integer n_tr, n_si
      double precision dg_symsin(100,8),dg_symtrp(100,8)
      double precision dg_symtrp1(100,8),dg_symtrp2(100,8)
      integer sym_orb(1000)
      common /SO/n_tr,n_si, 
     & VECCA,VECCB,VALA,dg_symsin,dg_symtrp,dg_symtrp1,
     & dg_symtrp2,sym_orb

C---------------------------------------------------------------
      integer NPT
      double precision PTX(1000000),PTY(1000000)
      double precision PTZ(1000000),PTW(1000000)
      common /F12_gr/ NPT,PTX,PTY,PTZ,PTW
C----------------------------------------------------------------
      integer ncis,n_s_roots,n_t_roots
      common /STEOM/ ncis,n_s_roots,n_t_roots

C----------------------------------------------------------------



