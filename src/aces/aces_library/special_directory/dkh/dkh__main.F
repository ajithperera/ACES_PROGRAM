#include "flags.h"
      Subroutine dkh__main(Dkh_order)
      
      Implicit Double Precision (A-H, O-Z)

#include "dkhlims.par"
#include "flags.com"
#include "machsp.com"

      Character*4 INP_Fname
      Character*4 ATMNAM(Max_centers), Namat(Max_centers)
      Logical Spherical, Cartesian 
      Integer*8 Istart_Heap, Ialpha, Ipcoef, Isalpha
      Integer*8 Ispcoef, Ipcoef_uc, Dbuf_start,Inext
      Integer*8 Ioed_ovl, Ioed_kin, Ioed_nai, Ioed_pvp
      Integer*8 Icoreham, Idkh, Ireord, Ileft, Ibuf_start

      Integer CCbeg(Maxbfns), CCend(Maxbfns), Zmax, Zhrr, 
     &        End_nfps, Val1(4), Val2(4), Angmom, Atm_4shell,
     &        Dkh_order 

      DIMENSION Work(1), Dbuf(1),Ibuf(1)

      DIMENSION NFCT(max_centers), NAOATM(max_centers),
     &          NUFCT(max_centers), NAOUATM(max_centers),
     &          iqmstr(max_centers),jcostr(max_centers,Maxjco),
     &          nucstr(max_centers,Maxang,Maxjco),
     &          nrcstr(max_centers,Maxang,Maxjco),
     &          Jstrt(max_shells),Jrs(max_shells),
     &          NCFpS(max_shells),Angmom(max_shells),
     &          NPFpS(max_shells), Atomic_label(max_shells),
     &          Charge(max_centers),Nspc(max_centers), 
     &          NPOP(max_centers),
     &          IREORDER(max_centers), 
     &          COORD(3,max_shells),Cord(3,max_centers),
     &          Ivshloff(max_shells),Nfps(max_shells),
     &          End_nfps(max_shells),Npfs(max_shells),
     &          End_nfps_uc(max_shells),
     &          ivAngMom(max_shells),Ixalpha(max_shells),
     &          Ixpcoef(max_shells),Indx_cc(max_shells),
     &          Ixshells(max_shells),Atm_4shell(max_shells),
     &          Temp(max_shells)

      Ispherical = 0
      Spherical  = .False.
      cartesian  = .False.

      If (IFLAGS(h_IFLAGS_spherical) .EQ. 1) Then
          Ispherical = 1
          Spherical  = .True.
      Else
          Cartesian  = .True.
      Endif 
C
C Note that Iecp is used only as a flag to get certain arrays filled
C (not as to indicate a wholesale ecp calculation). 

      Iecp      = 1
      INP_Fname = "MOL"
      Call SIMPLE_INSPECT_MOL_DKH(INP_Fname, max_centers, max_shells,
     &                        ncenters, nshells, nspc, cartesian, 
     &                        ITFCT, LNP1, lnpo, nfct, nufct, 
     &                        nbasis, NAOBASIS, nCFpS, 
     &                        nPFpS, NAOATM, angmom, atomic_label, 
     &                        vnn, Maxang, Maxjco, Iecp,
     &                        NUcstr, Nrcstr, Iqmstr, Jcostr, 
     &                        Jstrt, Jrs, Atmnam, Charge)


CSS#ifdef _DEBUG_LVL0
      Write(6,"(a,a,7(1x,I4))") "natoms,nshells,ITFCT,LNP1,lnpo,",
     &                          "nbasis,NAOBASIS:",ncenters,nshells,
     &                                         ITFCT,LNP1,lnpo,
     &                                        nbasis,NAOBASIS
      Write(6,*) 
      Write(6,"(a)") "The number of con. function per shell"
      Write(6,"(6(1x,I4))") (nCFpS(i), i=1, nshells)
      Write(6,"(a)") "The number of prim. function per shell"
      Write(6,"(6(1x,I4))") (nPFpS(i), i=1, nshells)
CSS#endif 

CSSS#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The Iqmstr"
      Write(6,"(6(1x,I4))") (Iqmstr(i), i=1, ncenters)
      Write(6,"(a)") "The jcostr"
      Do i=1, Ncenters
         Write(6,"(4(1x,I4))") (jcostr(i,j),j=1,4)
      Enddo
CSSS#endif
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The Nucstr"
      Do i=1, Ncenters
         Write(6,"(10(1x,I4))") ((Nucstr(i,j,k),j=1,4),k=1,ITFCT)
      Enddo
      Write(6,"(a)") "The Nrcstr"
      Do i=1, Ncenters
         Write(6,"(10(1x,I4))") ((Nrcstr(i,j,k),j=1,4),k=1,ITFCT)
      Enddo
      Write(6,"(a)") "The jstrt"
      Write(6,"(6(1x,I4))") (Jstrt(i), i=1, Nshells)
      Write(6,"(a)") "The jrs"
      Write(6,"(6(1x,I4))") (Jrs(i), i=1, Nshells)
#endif 

      nalpha    = 0
      npcoef    = 0
      npcoef_uc = 0 
      do i = 1, nshells
         nalpha    = nalpha + npfps(i)
         npcoef    = npcoef + npfps(i) * ncfps(i)
         npcoef_uc = npcoef_uc + npfps(i) * npfps(i) 
         Ivangmom(i) = Angmom(i)
      enddo

      Nprims = 0
      ncnfns = 0
      Do Ishell = 1, Nshells
         If (Spherical) Then
            lvalue = Ivangmom(Ishell)
            Nprims = Nprims + (2*lvalue+1) *
     &               npfps(ishell)
            Ncnfns = ncnfns  + (2*lvalue+1) * ncfps(ishell)
         Else
            lvalue = Ivangmom(Ishell)
            Nprims = Nprims + (lvalue+1)*(lvalue+2)/2 *
     &               npfps(ishell)
            Ncnfns = Ncnfns + (lvalue+1)*(lvalue+2)/2 *
     &               ncfps(ishell)
         Endif
    
      Enddo
C
      NBFNS = Nprims
C#ifdef _DEBUG_LVL0
      Write(6,"(a,4(1x,I4))") "nprims;nalpha;npcoef,Npcoef_uc: ",
     &                         nprims,nalpha,npcoef,Npcoef_uc
C#endif
C
C The base memory needed to handle DKH transformtion: The ovrlap
C kintic, and repulsion integrals and dkh_transfomed integrals
C need 5 * Nprims * Nprims. There are 3 more Nprims * Nprims arrays
C and 4 arrays of Nprims are needed for DKH (see DKH__MAIN_TRANS). 
C The Npcoef + Npcoef_uc space is for contraction and primitive
C coefficients. The Lnp1 + Lnpo are two temprarory place holders. 
C 
      Mem_needed =  4 * Nprims + 8 * Nprims * Nprims + Nalpha 
     &                + Npcoef + Npcoef_uc + Lnp1 + Lnpo
C
C The DKH order dependent memory requirements. The memory Ileft 
C amount is passed down to DKH routine. The same space is being 
C used during the contractions step. 
C
      If (Dkh_order .EQ. 1) Then
          Mem_factor = 6
          Ileft      =  Mem_factor * Nprims * Nprims 
          Mem_needed =  Mem_needed + Ileft 
          Mem_4dkh   =  Ileft +  4 * Nprims + 3 * Nprims * Nprims 
      Else if (Dkh_order .EQ. 2 .OR. Dkh_order .EQ. 3) Then
          Mem_factor = 8
          Ileft      =  Mem_factor * Nprims * Nprims 
          Mem_needed =  Mem_needed + Ileft 
          Mem_4dkh   =  Ileft +  4 * Nprims + 3 * Nprims * Nprims 
      Else if (Dkh_order .EQ. 4) Then
          Mem_factor = 18
          Ileft      =  Mem_factor * Nprims * Nprims 
          Mem_needed =  Mem_needed + Ileft 
          Mem_4dkh   =  Ileft +  4 * Nprims + 3 * Nprims * Nprims 
      Else if (Dkh_order .EQ. 5) Then
          Mem_factor = 34
          Ileft      =  Mem_factor * Nprims * Nprims 
          Mem_needed =  Mem_needed + Ileft 
          Mem_4dkh   =  Ileft +  4 * Nprims + 3 * Nprims * Nprims 
      Endif
C 
      Ierr = 0
      Call MEM_ALLOC(Work, Mem_needed, Ifltln, Istart_Heap, .True.,
     &               Ierr)

      If (Ierr .ne. 0) then
         Write(6, "(a,i10,a)") "@-Dkhmain: Unable to allocate",
     &         Mem_needed, " words of double precision memory"
        Call Errex
      Endif

      IALPHA    = Istart_Heap
      IPCOEF    = IALPHA + nalpha
      ISALPHA   = IPCOEF + npcoef
      ISPCOEF   = ISALPHA + LNP1
      IPCOEF_UC = ISPCOEF + LNPO
      INEXT     = IPCOEF_UC + npcoef_uc

      do i = 1, ncenters
         npop(i)     = 1
         ireorder(i) = i
         Namat(i) = ATMNAM(i)
      enddo
C
      Nshll = 0
      do iat=1,Ncenters
        do ilq=1,iqmstr(iat)
          do ijco=1,jcostr(iat,ilq)
             nshll=nshll+1
             Atm_4shell(nshll)=iat
           Enddo
        Enddo
      Enddo
C
      Call READ_BASIS_INFO_DKH(INP_FNAME, Ncenters, Ncenters, NPOP,
     &                     IREORDER,
     &                     CARTESIAN, NALPHA,NPCOEF, LNP1, LNPO,
     &                     NFCT, NBASIS, WORK(IALPHA), IXALPHA,
     &                     Work(IPCOEF), IXPCOEF, Work(ISALPHA),
     &                     Work(ISPCOEF), MAX_CENTERS,
     &                     ATMNAM, COORD, CORD, NAOATM)
C
CSSS#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The exponents"
      Write(6,"(6(1x,F12.8))") (Work(Ialpha+i), i=0, Nalpha-1)
      Write(6,*)
      Write(6,"(a)") "The contraction coefs"
      Write(6,"(6(1x,F12.8))") (Work(Ipcoef+i), i=0, Npcoef-1)
      Write(6,"(a)") "The jcostr"
      Do i=1, Ncenters
         Write(6,"(4(1x,I4))") (jcostr(i,j),j=1,4)
      Enddo

      Write(6,*)
      Write(6,*) "Atom Coordiantes of each shell"
      do Iatm =1, Ncenters
         Write(6, "(6(1x,F12.8))") (Cord(i,Iatm), i=1,3)
      Enddo
CSSS#endif
C
      Ixpcoef(1) = 1
      do i = 1, Nshells 
          Ixpcoef(i+1) = Npfps(i) * Npfps(i) + Ixpcoef(i)
      enddo 
C
C Here the number of contracted functions are the same as primitive 
C functions. 

      Call Dkh__int_gen_init(Nshells, Ivshloff, Ivangmom, Nfps, Npfs, 
     &                       Npfps,
     &                       Ixshells, Ixalpha, Ixpcoef, End_nfps, 
     &                       Temp, Ispherical)
CSSS#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The Ivangmom array" 
      Write(6,"(8(1x,I4))") (Ivangmom(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Ixshells array" 
      Write(6,"(8(1x,I4))") (Ixshells(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Nfps array" 
      Write(6,"(8(1x,I4))") (Nfps(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The IVshloff array" 
      Write(6,"(8(1x,I4))") (Ivshloff(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Ixpcoef array" 
      Write(6,"(8(1x,I4))") (Ixpcoef(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Ixalpha array" 
      Write(6,"(8(1x,I4))") (Ixalpha(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The End_nfps array" 
      Write(6,"(8(1x,I4))") (End_nfps(i), i=1, Nshells)
      Write(6,*)
      Write(6,*) "The Npfps array" 
      Write(6,"(8(1x,I4))") (Npfps(i), i=1, Nshells)
CSSS#endif 
C
      call dzero(work(ipcoef_uc), npcoef_uc)
      call dkh__prim_ccoefs(Work(ipcoef_uc), Ixpcoef, Npfps, Nshells,
     &                      Npcoef_uc)

CSSS#ifdef _DEBUG_LVL0
      write(6,*)
      Write(6,"(a)") "The primitive contraction coefs"
      Write(6,"(4(1x,F10.5))") (work(ipcoef_uc +i),i=0,npcoef_uc-1)
      Write(6,*)
CSSS#endif 
C
C Here the number of contracted functions are the same as primitive 
C functions. 
C
      Call setup_ccbeg(WORK(Ialpha), IXalpha, Work(Ipcoef_uc),Ixpcoef,
     &                 Npfps, Npfps, Nshells, CCbeg, CCend, Indx_cc)
C
CSSS#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*)"The CCbeg, CCend arrays and indx_cc arrays" 
      Write(6,"(10(1x,I4))") (CCbeg(i), i=1, Itfct)
      Write(6,"(10(1x,I4))") (CCend(i), i=1, Itfct)
      Write(6,"(10(1x,I4))") (Indx_cc(i), i=1, Nshells)
CSSS#endif 
 
      Call Dkh__erd_scratch_mem_calc(nshells, ivangmom, npfps,
     &                               npfps, atom, Coord, 
     &                               Work(ialpha),
     &                               Work(ipcoef_uc), 
     &                               ixalpha, ixpcoef,
     &                               Ccbeg, Ccend, indx_cc,
     &                               spherical, Ncenters, .true., 
     &                               intmax, zmax)

      Call Erd__memory_hrr_correction(ivAngMom, nshells,
     &                                spherical, ihrr, zhrr)
     
      Intmax = Intmax + Ihrr
      Zmax   = Zmax   + Zhrr

      Ierr = 0
      Call Mem_alloc(dbuf, Zmax, Ifltln, dbuf_start, .true., ierr)
      If (Ierr .ne. 0) then
          Write(6, "(a,i10,a)") "@DKH_main: Unable to allocate",
     &              Zmax, " words of double precision memory"
          Call Errex
      Endif

      Ierr = 0
      Call Mem_alloc(ibuf, Intmax, Iintln, ibuf_start,
     *                     .true., ierr)
      If (Ierr .ne. 0) then
          Write(6, "(a,i10,a)") "@DKH_main: Unable to allocate",
     &              Zmax, " words of integer memory"
          Call Errex
      Endif

CSSS#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,I6,1x,I10)") "ERD memmory requirments Ibuf & Dbuf",
     &                            Intmax, Zmax
CSSS#endif 

      Ioed_ovl      = INext
      Ioed_kin      = Ioed_ovl    + Nprims**2
      Ioed_nai      = Ioed_kin    + Nprims**2
      Ioed_pvp      = Ioed_nai    + Nprims**2
      Icoreham      = Ioed_pvp    + Nprims**2
      Idkh          = Icoreham    + Nprims**2 
      IEnd          = Idkh        + Mem_4dkh
C
      If ((Iend-Istart_Heap) .Gt. Mem_needed) then
         Write(6, "(a,i10,a)") "@DKHAIN: Unable to allocate",
     &         Mem_needed, " words of double precision memory"
        Call Errex
      Endif
C
C  There is one segment in serial work (not necessary, but lets work
C  with that in mind).
C
      Val1(1) = 1
      Val1(2) = 1
      Val2(1) = Nprims
      Val2(2) = Nprims 
C
CSSS#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "Atom Coordinates of each shell"
      Do Iatm =1, Ncenters
         Write(6, "(6(1x,F12.8))") (Cord(i,Iatm), i=1,3)
      Enddo
      Write(6,*)
      Write(6,*) "The atom label for each shell in main",Nshells
      Write(6,"(4(1x,I4))") (Atm_4shell(i), i=1, Nshells)
      Write(6,*)
      Write(6,"(a,i10)") "Memory for DKH:", Mem_4dkh
      Write(6,"(a,i10)") "Dkh memory test:", 34*nprims*nprims
     +                   +  4 * Nprims + 3 * Nprims * Nprims 
CSSS#endif 

      Call Dkh__compute_oed_integrals(
     &                       val1(1),val2(1),val1(2),val2(2),
     &                       Dbuf(dbuf_start),Zmax,
     &                       Ibuf(ibuf_start),
     &                       Cord,Work(ipcoef),Work(ipcoef_uc),
     &                       Work(ialpha),
     &                       Ccbeg,Ccend,End_nfps, Nshells, Nbfns,
     &                       Zmax, Intmax, Npfps, Ncfps, 
     &                       Ivangmom,
     &                       Ixalpha, Indx_cc, Atm_4shell, Ixpcoef,
     &                       Ispherical,
     &                       Work(Ioed_ovl),
     &                       Work(Ioed_nai),
     &                       Work(Ioed_kin),
     &                       Work(Ioed_pvp), 
     &                       Work(Icoreham),
     &                       Nsend,
     &                       Nalpha,Nprims,Ncnfns,Npcoef,Npcoef_uc,
     &                       Ncenters,
     &                       charge,
     &                       max_centers,Work(Idkh),Mem_4dkh,Dkh_order)

      Return
      End
