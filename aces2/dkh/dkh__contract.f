










C
      subroutine dkh__contract(zcore,maxcor,coreham,dkh_ints,coeffs,
     +                         nshells,nprims,ncnfns,npfps,npcoef,
     +                         npcoef_uc,ncfps,Ivangmom,Spherical)

      Implicit Double Precision (a-h, o-z)

      parameter (max_centers = 300)
      parameter (max_shells  = 5000)
      parameter (max_prims   = max_shells)
      Parameter (max_cbf     = 1000)
      parameter (max_primcc  = Max_prims*max_cbf)
      parameter (Ndi4 = 550, Ndi9 = max_shells, Ndi10 = max_centers,
     &           Ndi13 = 350, Ndico = 10, ndi14 = 120, ndi27 = 400)
     &
      parameter (Maxang = 6)
      parameter (Maxproj =5)
      parameter (ndilmx = Maxang+1)
      parameter (nh4=4*(ndilmx-1)-3)
      parameter (maxjco = 10)
      Parameter(Max_ecpmem = 50000)
C
C



      Dimension zcore(maxcor),coreham(nprims*nprims),coeffs(npcoef)
      Dimension dkh_ints(ncnfns*ncnfns)
      Integer npfps(nshells),ncfps(nshells)
      Integer ivangmom(Nshells)
      Integer coef_start,zcor_start,prim_off(max_shells)
      Integer cont_off(max_shells),Ccoef_off(max_shells)
      Integer Tot_prm_fns,Tot_cnt_fns,Tot_cnt_ish,Tot_cnt_jsh
      Integer Icnt_block
      
      logical spherical 

      prim_off(1)  = 1
      cont_off(1)  = 1
      ccoef_off(1) = 1

      Do ishell = 2, Nshells
         ncont_shell = ncfps(ishell-1)
         nprim_shell = npfps(ishell-1)

         if (spherical) then
            ndegen = 2*ivangmom(ishell-1) + 1
         else
            ndegen = (ivangmom(ishell-1)+1) * (ivangmom(ishell-1)+2)/2
         endif

         prim_off(ishell) = prim_off(ishell-1) + nprim_shell * ndegen
         cont_off(ishell) = cont_off(ishell-1) + ncont_shell * ndegen   
         Ccoef_off(ishell)= Ccoef_off(ishell-1) + 
     +                      ncont_shell * nprim_shell
      Enddo



      Do Jshell = 1, Nshells
     
         Do Ishell = 1, Nshells 

            Iprim_4shell = Npfps(Ishell)
            Jprim_4shell = Npfps(Jshell)
            Icntr_4shell = Ncfps(Ishell)
            Jcntr_4shell = Ncfps(Jshell)

            If (Spherical) Then
               Idegen = 2*ivangmom(ishell) + 1
               Jdegen = 2*ivangmom(jshell) + 1
            Else 
               Idegen = (ivangmom(ishell)+1) * (ivangmom(ishell)+2)/2
               Jdegen = (ivangmom(jshell)+1) * (ivangmom(jshell)+2)/2
            Endif
           
            Iprimfns = Iprim_4shell * Idegen
            Jprimfns = Jprim_4shell * Jdegen
            
            Icntrfns = Icntr_4shell * Idegen
            Jcntrfns = Jcntr_4shell * Jdegen

            Tot_prm_fns = Iprimfns * Jprimfns 
            Tot_cnt_fns = Icntrfns * Jcntrfns 
            Tot_cnt_ish = Iprim_4shell * Icntr_4shell 
            Tot_cnt_jsh = Jprim_4shell * Jcntr_4shell 
            Tot_bfns    = Ncnfns * Ncnfns
C

            Istart      = 1
            Ipnt_block  = Istart 
            Icn_block   = Ipnt_block + Tot_prm_fns 
            Jcn_block   = Icn_block  + Tot_cnt_ish
            Ict_block   = Jcn_block  + Tot_cnt_jsh
            Inext       = Ict_block  + Tot_cnt_fns 


            Call get_ints_block(coreham,zcore(Ipnt_block),ishell,
     +                          jshell,Nprims,prim_off,cont_off,
     +                          Iprimfns,jprimfns,Nshells)

            Call get_ccoefs_block(Coeffs,zcore(Icn_block),ishell,
     +                            Npcoef,ccoef_off,
     +                            Tot_cnt_ish,Nshells)

            Call get_ccoefs_block(Coeffs,zcore(Jcn_block),jshell,
     +                            Npcoef,Ccoef_off,
     +                            Tot_cnt_jsh,NShells)

            Call contract(zcore(Ipnt_block),zcore(Icn_block),
     +                    zcore(Jcn_block),Iprimfns,jprimfns,
     +                    Icntrfns,Jcntrfns,Iprim_4shell,
     +                    Jprim_4shell,Icntr_4shell,Jcntr_4shell,
     +                    Tot_cnt_fns,Tot_prm_fns,Idegen,Jdegen,
     +                    zcore(Ict_block))

            Call form_dkh_ints(zcore(Ict_block),Dkh_ints,
     +                         Ishell,jshell,Icntrfns,Jcntrfns,
     +                         Ncnfns,prim_off,cont_off,NShells)

         Enddo
      Enddo


      Return 
      End
