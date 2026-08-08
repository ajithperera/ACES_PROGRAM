      double precision VECCA(2000,100,12)
      double precision VECCB(2000,100,12)
      double precision VALA(20)
      integer n_tr, n_si
      double precision dg_symsin(100,8),dg_symtrp(100,8)
      double precision dg_symtrp1(100,8),dg_symtrp2(100,8)
      integer sym_orb(1000)  
      common /SO/VECCA,VECCB,VALA,dg_symsin,dg_symtrp,dg_symtrp1,
     & dg_symtrp2,sym_orb
