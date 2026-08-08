










      Subroutine Scrnc_jj_cc_lrspn(Work,Maxcor,Iuhf,JJ_lrespn,Isympert,
     +                             Ioffpert,Nao_pair)

      Implicit Double Precision (A-H,O-Z)

c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
      
      Double Precision M_expect,Mbar_00
      Double Precision JJ_lrespn(Nao_pair,Nao_pair)
      Logical Incore
      Logical ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd
      Logical Realfreq
      Integer Occ_pair,Occ_pair_4irrep
      Dimension Isympert(Nao_pair),Ioffpert(Nao_pair)
      Dimension Work(Maxcor)
      Dimension D(Maxbasfn,Maxbasfn)
      Character*8 Label(3),Work_label

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      Common/Sym_pqxx_ints/Irrep_aos(Maxbasfn),
     &                     Irrep_ao_pairs(Maxbasfn*Maxbasfn)
      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +       Drccd
      Common /Freqinfo/Freq,Realfreq,Xfreq(200),No_freqs 
      Common /eomprop/Ieomprop

      Data Label /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/



c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
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

      Irreprl     = 1
      Npert_count = 0
      Do Irrep = 1, Nirrep
         Nperts_4irrepr = Irpdpdao(Irrep)

         Iao_pair_4irrep = 0
         Do Iao_pair= 1, Nperts_4irrepr

         Iao_pair_4irrep = Iao_pair_4irrep + 1 
         Npert_count     = Npert_count     + 1
         Irrep_ao_pair = Irrep_ao_pairs(Npert_count)
         Irrepr = Irrep_ao_pair
         Irrepl = Dirprd(Irrepr,Irreprl)
         Irrepx = Irrepl

      Nbfns   = Nocco(1) + Nvrto(1)
      Call scrnc_get_pert_type(Work,Maxcor,Nbfns,Label,Isympert,
     +                        Iuhf)
      Do Occ_pair= 1, 3
      Irrepx     = Isympert(occ_pair)
      Ipert      = occ_pair
      Work_label = Label(occ_pair)
         Nsizec = Irpdpd(Irrepx,9)
         If (Iuhf .EQ. 0) Then
            Nsizec = Nsizec  + Idsymsz(Irrepx,13,14)
         Else
            Nsizec = Nsizec + Irpdpd(Irrepx,10)
            Nsizec = Nsizec + Idsymsz(Irrepx,1,3)
            Nsizec = Nsizec + Idsymsz(Irrepx,2,4)
            Nsizec = Nsizec + Idsymsz(Irrepx,13,14)
         Endif

         Call Scrnc_reset_lrspn_cc_lists(Work,Maxcor,Iuhf,Irrepx)

C List_ljki, List_abki and List_bjki are the OO, VV and OV integral
C list
         Lenoo = Irpdpd(Irrepx,21) + Iuhf*Irpdpd(Irrepx,22)
         Lenvv = Irpdpd(Irrepx,19) + Iuhf*Irpdpd(Irrepx,20)
         Lenvo = Irpdpd(Irrepx,9)  + Iuhf*Irpdpd(Irrepx,10)

         I000 = 1
         I010 = I000 + Lenoo
         I020 = I010 + Lenvv
         I030 = I020 + Lenvo
         I040 = I030 + Lenoo
         I050 = I040 + Lenvv
         I060 = I050 + Lenvo
         I070 = I060 + Lenvo
         Iend = I070
         Mxcor = Maxcor - Iend 

         If (Iend .Gt. Maxcor) Call Insmem("Scrnc_jj_lrspn",Iend,
     +                                      Maxcor)
         Ioff_oo = I000
         Ioff_vv = I010 
         Ioff_vo = I020 


         Call Scrnc_prep_culomb_ints(Work(Ioff_oo),Work(Ioff_vv),
     +                                Work(Ioff_vo),Iuhf,Lenoo,
     +                                Lenvv,Lenvo,Iao_pair_4irrep,
     +                                Irrepx,Irrepr)
      Nbfns   = Nocco(1) + Nvrto(1)
      Call Getrec(20, "JOBARC","NBASTOT ", 1, Naobfns)
      Length   = Naobfns * (Naobfns+1)/2
      Length21 = Nbfns * Nbfns
      Length22 = Nbfns * Naobfns
      I080 = I070 + Length
      I090 = I080 + Length21
      I100 = I090 + Length22
      I110 = I100 + Length22
      Iend = I110
      Mxcor = Maxcor - Iend
      If (Iend .Gt. Maxcor) Call Insmem("Scrnc_lrspn_main",Iend,
     +                                  Maxcor)
      Call Scrnc_prep_dipole_ints(Work(I070),Work(I080),Work(I090),
     +                           Work(I100),Work(I000),Work(I010),
     +                           Work(I020),Work(Iend),
     +                           Memleft,Lenoo,Lenvv,Lenvo,Nbfns,
     +                           Naobfns,Iuhf,Ipert,Work_label,
     +                           Irrepx)
      Call Dcopy(Lenoo,Work(I000),1,Work(I030),1)
      Call Dcopy(Lenvv,Work(I010),1,Work(I040),1)
      Call Dcopy(Lenvo,Work(I020),1,Work(I050),1)
      Call Dcopy(Lenvo,Work(I020),1,Work(I060),1)
         Call Dcopy(Lenoo,Work(I000),1,Work(I030),1)
         Call Dcopy(Lenvv,Work(I010),1,Work(I040),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I050),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I060),1)

C We construct e^-TMue(pq)^T, similarity transformed Mu(p,q).
            
         Iside = 1
         Call Scrnc_form_mubar_00(Work(Ioff_oo),Work(Ioff_vv),
     +                            Work(Ioff_vo),Work(Iend),
     +                            Mxcor,Iuhf,Irrepx,M_expect,
     +                            Mbar_00)

         Call Scrnc_form_mubar_s(Work(Ioff_oo),Work(Ioff_vv),
     +                           Work(Ioff_vo),Work(I030),
     +                           Work(I040),Work(I050),
     +                           Work(I060),Work(Iend),Mxcor, 
     +                           Iuhf,Irrepx,Iside,Mbar_00,
     +                           Ieomprop)

         If (Iuhf .eq. 0) Then 
            Call Scrnc_form_mubar_d_rhf(Work(I030),Work(I040),
     +                                  Work(I050),Work(Iend),
     +                                  Mxcor,Iuhf,Irrepx,
     +                                  Iside,Mbar_00,Ieomprop)
                                          
         Else
            Call Scrnc_form_mubar_d_uhf(Work(I030),Work(I040),
     +                                  Work(I050),Work(Iend),
     +                                  Mxcor,Iuhf,Irrepx,
     +                                  Iside,Mbar_00,Ieomprop)
         Endif 

         Iside = 2
         Ibgn  = I070 
         I080  = I070 + Lenoo
         I090  = I080 + Lenvv
         I100  = I090 + Lenvo
         I110  = I100 + Lenvo
         Iend  = I110
         Mxcor = Maxcor - Iend 
         If (Iend .Gt. Maxcor) Call Insmem("Scrnc_lrspn_main",Iend,
     +                                     Maxcor)

         Call Dcopy(Lenoo,Work(I000),1,Work(I070),1)
         Call Dcopy(Lenvv,Work(I010),1,Work(I080),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I090),1)
         Call Dcopy(Lenvo,Work(I020),1,Work(I100),1)
 
         Call Scrnc_form_mubar_s(Work(I030),Work(I040),
     +                           Work(I050),Work(I070),
     +                           Work(I080),Work(I090),
     +                           Work(I100),Work(Iend),Mxcor,   
     +                           Iuhf,Irrepx,Iside,Mbar_00-
     +                           M_expect,Ieomprop)

         If (Iuhf .eq. 0) Then   
            Call Scrnc_form_mubar_d_rhf(Work(I030),Work(I040),
     +                                  Work(I060),Work(Iend),
     +                                  Mxcor,Iuhf,Irrepx,
     +                                  Iside,Mbar_00-Mexpect,
     +                                  Ieomprop)
         Else
            Call Scrnc_form_mubar_d_uhf(Work(I030),Work(I040),
     +                                  Work(I060),Work(Iend),
     +                                  Mxcor,Iuhf,Irrepx,
     +                                  Iside,Mbar_00-Mexpect,
     +                                  Ieomprop)
         Endif   

          If (Iside .Eq. 2) Then 
             Mu_s    = 482
             Mu_d    = 453
             Ioff_s  = 0
             Ioff_d  = 0
             Ioffset = 1
             Call Scrnc_load_vec(Irrepx,Work(I000),Maxcor,Mu_s,
     +                           Ioff_s,Mu_d,Ioff_d,Iuhf,.False.)
             Call Putlst(Work(I000),1,1,0,Irrepx,374)

      write(6,*)
      call checksum("Jbar_L<(S,D)|pq>:",Work(I000),Nsizec)
          Endif 

C Save the left amplitudeds.

          Call Getlst(Work(I000),Ioffset,1,1,Irrepx,372+Iside)
          Call Putlst(Work(I000),Npert_count,1,1,Irrepx,372+Iside)

         Enddo ! Enddo for Iao_pair 
      Enddo ! Enddo of irrepr

      Call Dzero(JJ_lrespn,Nao_pair*Nao_pair)

      I_pairs = Nao_pair
      J_pairs = Nao_pair

      Write(6,*)
      Write(6,"(a)") " Symmetry of perturbations"
      Write(6,"(5(1x,I5))") (Isympert(i),i=1,Nao_pair)
      Write(6,"(a)") " Offset for perturbations"
      Write(6,"(5(1x,I5))") (Ioffpert(i),i=1,Nao_pair)
      Write(6,*)
      Do Iao_pair = 1, I_pairs 

         Irrepx = Isympert(Iao_pair)
         Call Scrnc_newsize(Work,Maxcor,Iuhf,Irrepx,Nsizec)

         I000 = 1
         I010 = I000 + Nsizec 
         I020 = I010 + Nsizec 
         Iend = I020

         Memleft = Maxcor - Iend 
         If (Iend .Gt. Maxcor)  Call Insmem("Scrnc_jj_lrspn",Iend,
     +                                      Memleft)

         Call Getlst(Work(I000),Iao_pair,1,1,Irrepx,373)
         If (Iuhf .eq. 0) Then
            Call Spntsing(Isympert(Iocc_pair),Work(I000),Work(Iend),
     +                    Memleft) 
         Endif 

         Do Jao_pair = 1, J_pairs 

            If (Isympert(Iao_pair) .eq. Isympert(Jao_pair)) Then
 
               Call Getlst(Work(I010),Jao_pair,1,1,Irrepx,
     +                     374)
               X = Ddot(Nsizec,Work(I010),1,Work(I000),1)

               JJ_lrespn(Iao_pair,Jao_pair) = 
     +                  JJ_lrespn(Iao_pair,Jao_pair) + X
            Endif 
              
         Enddo
      Enddo 
      Write(6,*)
      Write(6,"(a)") " CI-like JJ matrix"
      call output(JJ_lrespn,1,Nao_pair,1,Nao_pair,Nao_pair,
     +            Nao_pair,1)
      Write(6,*)
      
      If (Ieomprop .eq. 2) Then
          Itop = 1
          Call Calcquad1(Work(I000),JJ_lrespn,Maxcor,Isympert,
     +                   Ioffpert,Nao_pair,Iuhf,1.0D0,
     +                   Nao_pair,1,Itop,.False.)
      Endif 
      Write(6,*)
      Write(6,"(a)") " LR JJ matrix"
      call output(JJ_lrespn,1,Nao_pair,1,Nao_pair,Nao_pair,
     +            Nao_pair,1)
      Nbfns   = Nocco(1) + Nvrto(1)
      i000=1
      i010=i000+Nbfns*Nbfns 
      call getrec(20,"JOBARC","HFDENSTY",Nbfns*Nbfns,work(i000))
      ind=0
      do i=1,nbfns 
      do j=1,nbfns
         ind=ind+1
         d(j,i) = Work(i000-1+ind)
      enddo
      enddo
      ind=0
      do l=1,Nbfns
      do k=1,Nbfns
      do j=1,Nbfns 
      do i=1,Nbfns
         ind= ind + 1
         Work(i010-1+ind)= D(k,l)* D(i,j)
      enddo
      enddo
      enddo
      enddo
      nbfns2=Nbfns*Nbfns 
      write(6,*) "The SCF density matrix product"
      call output(work(i010),1,nbfns2,1,nbfns2,nbfns2,nbfsn2,1)
      e2 = Ddot(Nbfns2*Nbfns2,JJ_lrespn,1,Work(i010),1)
      Write(6,"(a,F15.9)") "The pseudo energy: ", e2
    
      Return
      End 
