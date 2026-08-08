





















c#define _DEBUG_GETAOINF

      subroutine Scrnc_getaoinf(iUHF,irp_x)
      IMPLICIT INTEGER (A-Z)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end


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



c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

c   o pick up population vector of AOs (orbitals per irrep)

c JDW 7/1/96
c   It is incorrect to read NUMBASIR when orbitals have been
c dropped, since DROPVC (from xvtran) redefines this record
c from its original meaning. I have written original contents
c of NUMBASIR on NUMBASI0 (in DROPVC) and this should be
c read to get IAOPOP when orbitals have been dropped.

      call getrec(1,'JOBARC','NUMDROPA',1,nDrop)
      if (nDrop.eq.0) then
         call getrec(1,'JOBARC','NUMBASIR',nirrep,IAOPOP)
      else
         call getrec(1,'JOBARC','NUMBASI0',nirrep,IAOPOP)
      end if

c   o calculate cumulative AO indices of each irrep block
      IOFFAO(1) = 1
      do irrep = 1, nirrep-1
         IOFFAO(irrep+1) = IOFFAO(irrep) + IAOPOP(irrep)
      end do

c   o calculate cumulative MO indices of each irrep block (vrt imm. follows occ)
      ioffv(1,1) = nocco(1) + 1
      ioffo(1,1) = 1
      do irrep = 1, nirrep-1
         ioffv(irrep+1,1) = ioffv(irrep,1) + vrt(irrep,1)
         ioffo(irrep+1,1) = ioffo(irrep,1) + pop(irrep,1)
      end do
      if (iUHF.eq.0) then
         do irrep = 1, 8
            ioffv(irrep,2) = ioffv(irrep,1)
            ioffo(irrep,2) = ioffo(irrep,1)
         end do
      else
         ioffv(1,2) = nocco(2) + 1
         ioffo(1,2) = 1
         do irrep = 1, nirrep-1
            ioffv(irrep+1,2) = ioffv(irrep,2) + vrt(irrep,2)
            ioffo(irrep+1,2) = ioffo(irrep,2) + pop(irrep,2)
         end do
      end if

c   o create an irpdpd vector for the full (IRPDPDAO) and packed (IRPDPDAOS)
c     AO representations
      IRPDPDAO(1)      = 0
      IRPDPDAOS(1)     = 0

      IRPDPDAOMO_OCCKB(1,1) = 0
      IRPDPDAOMO_OCCKB(1,2) = 0
      IRPDPDAOMO_VRTKB(1,1) = 0
      IRPDPDAOMO_VRTKB(1,2) = 0

      IRPDPDAOMO_OCCBK(1,1) = 0
      IRPDPDAOMO_OCCBK(1,2) = 0
      IRPDPDAOMO_VRTBK(1,1) = 0
      IRPDPDAOMO_VRTBK(1,2) = 0

      do irp_ket = 1, nirrep
         IRPDPDAO(1)  =   IRPDPDAO(1)
     &                  +   IAOPOP(irp_ket)*IAOPOP(irp_ket)

         Do Ispin = 1, Iuhf+1
             IRPDPDAOMO_OCCKB(1,Ispin) =   IRPDPDAOMO_OCCKB(1,Ispin)
     &                      +   POP(irp_ket,Ispin)*IAOPOP(irp_ket)
             IRPDPDAOMO_VRTKB(1,Ispin) =   IRPDPDAOMO_VRTKB(1,Ispin)
     &                      +   VRT(irp_ket,Ispin)*IAOPOP(irp_ket)

             IRPDPDAOMO_OCCBK(1,Ispin) =   IRPDPDAOMO_OCCBK(1,Ispin)
     &                      +   POP(irp_ket,Ispin)*IAOPOP(irp_ket)
             IRPDPDAOMO_VRTBK(1,Ispin) =   IRPDPDAOMO_VRTBK(1,Ispin)
     &                      +   VRT(irp_ket,Ispin)*IAOPOP(irp_ket)
         Enddo

         IRPDPDAOS(1) =   IRPDPDAOS(1)
     &                  + IAOPOP(irp_ket)*(IAOPOP(irp_ket)-1)/2
      end do

      if (nirrep.gt.1) then
      do irrep = 2, nirrep
         IRPDPDAO(irrep)      = 0
         IRPDPDAOS(irrep)     = 0

         IRPDPDAOMO_OCCBK(irrep,1) = 0
         IRPDPDAOMO_OCCBK(irrep,2) = 0
         IRPDPDAOMO_VRTBK(irrep,1) = 0
         IRPDPDAOMO_VRTBK(irrep,2) = 0

         IRPDPDAOMO_OCCKB(irrep,1) = 0
         IRPDPDAOMO_OCCKB(irrep,2) = 0
         IRPDPDAOMO_VRTKB(irrep,1) = 0
         IRPDPDAOMO_VRTKB(irrep,2) = 0

         do irp_ket = 1, nirrep
            irp_bra = dirprd(irp_ket,irrep)
            if (irp_bra.gt.irp_ket) then
               iTmp = IAOPOP(irp_bra)*IAOPOP(irp_ket)
               IRPDPDAO(irrep)  = IRPDPDAO(irrep)  + iTmp

               Do Ispin = 1, Iuhf+1
                  IRPDPDAOMO_OCCBK(irrep,Ispin)  = 
     &            IRPDPDAOMO_OCCBK(irrep,Ispin)  + IAOPOP(irp_ket)*
     &                                           POP(irp_bra,Ispin)
                  IRPDPDAOMO_VRTBK(irrep,Ispin)  = 
     &            IRPDPDAOMO_VRTBK(irrep,Ispin)  + IAOPOP(irp_ket)*
     &                                           VRT(irp_bra,Ispin)
                  IRPDPDAOMO_OCCKB(irrep,Ispin)  = 
     &            IRPDPDAOMO_OCCKB(irrep,Ispin)  + IAOPOP(irp_bra)*
     &                                           POP(irp_ket,Ispin)
                  IRPDPDAOMO_VRTKB(irrep,Ispin)  = 
     &            IRPDPDAOMO_VRTKB(irrep,Ispin)  + IAOPOP(irp_bra)*
     &                                           VRT(irp_ket,Ispin)

               Enddo 
                 
               IRPDPDAOS(irrep) = IRPDPDAOS(irrep) + iTmp
            else
               IRPDPDAO(irrep) =   IRPDPDAO(irrep)
     &                           + IAOPOP(irp_bra)*IAOPOP(irp_ket)
               Do Ispin = 1, Iuhf+1
                  IRPDPDAOMO_OCCBK(irrep,Ispin)  = 
     &            IRPDPDAOMO_OCCBK(irrep,Ispin)  + IAOPOP(irp_ket)*
     &                                           POP(irp_bra,Ispin)
                  IRPDPDAOMO_VRTBK(irrep,Ispin)  = 
     &            IRPDPDAOMO_VRTBK(irrep,Ispin)  + IAOPOP(irp_ket)*
     &                                           VRT(irp_bra,Ispin)

                  IRPDPDAOMO_OCCKB(irrep,Ispin)  = 
     &            IRPDPDAOMO_OCCKB(irrep,Ispin)  + IAOPOP(irp_bra)*
     &                                          POP(irp_ket,Ispin)
                  IRPDPDAOMO_VRTKB(irrep,Ispin)  = 
     &            IRPDPDAOMO_VRTKB(irrep,Ispin)  + IAOPOP(irp_bra)*
     &                                           VRT(irp_ket,Ispin)
               Enddo 
            end if
         end do
      end do
      end if

c   o calculate cumulative offsets of AO ket irreps per total irrep
      do irrep = 1, nirrep
         ISTART(1,irrep) = 0
         do irp_ket = 1, nirrep-1
            irp_bra = dirprd(irp_ket,irrep)
            ISTART(irp_ket+1,irrep) =   ISTART(irp_ket,irrep)
     &                                + IAOPOP(irp_bra)*IAOPOP(irp_ket)
         end do
      end do

c   o calculate cumulative offsets of MO ket irreps per total irrep, in which
c     each column distribution is a full AO distribution of irp_bra (irp_ab).
      if (iUHF.eq.0) then
         ISTARTMO(1,3) = 0
         do irp_ij = 1, nirrep-1
            irp_ab = dirprd(irp_ij,irp_x)
            ISTARTMO(irp_ij+1,3) =   ISTARTMO(irp_ij,3)
     &                             + (   IRPDPDAO(irp_ab)
     &                                 * irpdpd(irp_ij,14) )
         end do
      else
         iOff = 0
         ISTARTMO(1,3) = 0
         do irp_ij = 1, nirrep-1
            irp_ab = dirprd(irp_ij,irp_x)
            iOff = iOff + (   IRPDPDAO(irp_ab)
     &                      * irpdpd(irp_ij,14) )
            ISTARTMO(irp_ij+1,3) = iOff
         end do
         irp_ab = dirprd(nirrep,irp_x)
         iOff = iOff + IRPDPDAO(irp_ab)*irpdpd(nirrep,14)
         ISTARTMO(1,2) = iOff
         do irp_ij = 1, nirrep-1
            irp_ab = dirprd(irp_ij,irp_x)
            iOff = iOff + (   IRPDPDAO(irp_ab)
     &                      * irpdpd(irp_ij,4) )
            ISTARTMO(irp_ij+1,2) = iOff
         end do
         irp_ab = dirprd(nirrep,irp_x)
         iOff = iOff + IRPDPDAO(irp_ab)*irpdpd(nirrep,4)
         ISTARTMO(1,1) = iOff
         do irp_ij = 1, nirrep-1
            irp_ab = dirprd(irp_ij,irp_x)
            iOff = iOff + (   IRPDPDAO(irp_ab)
     &                      * irpdpd(irp_ij,3) )
            ISTARTMO(irp_ij+1,1) = iOff
         end do
      end if

      print *, '@GETAOINF: AOSYM common block:'
      print *, 'iaopop(irrep) = ',(iaopop(irrep),irrep=1,nirrep)
      print *, 'ioffao(irrep) = ',(ioffao(irrep),irrep=1,nirrep)
      do iSpin = 1, 2
      print *, 'spin = ',iSpin
      print *, 'occ(irrep) = ',(pop(irrep,iSpin),irrep=1,nirrep)
      print *, 'vrt(irrep) = ',(vrt(irrep,iSpin),irrep=1,nirrep)
      print *, 'ioffo(irrep) = ',(ioffo(irrep,iSpin),irrep=1,nirrep)
      print *, 'ioffv(irrep) = ',(ioffv(irrep,iSpin),irrep=1,nirrep)
      end do
      print *, 'irpdpdao(irrep) = ',(irpdpdao(irrep),irrep=1,nirrep)
      do ispin=1,2
      print *, 'irpdpdaomo_vrtbk(irrep) = ',
     +                (irpdpdaomo_vrtbk(irrep,ispin),irrep=1,nirrep)
      print *, 'irpdpdaomo_occbk(irrep) = ',
     +                (irpdpdaomo_occbk(irrep,ispin),irrep=1,nirrep)
      print *, 'irpdpdaomo_vrtkb(irrep) = ',
     +                (irpdpdaomo_vrtkb(irrep,ispin),irrep=1,nirrep)
      print *, 'irpdpdaomo_occkb(irrep) = ',
     +                (irpdpdaomo_occkb(irrep,ispin),irrep=1,nirrep)
      enddo
      print *, 'irpdpdaos(irrep) = ',(irpdpdaos(irrep),irrep=1,nirrep)
      do irrep = 1, nirrep
      print *, 'istart(ket) = ',(istart(irp_ket,irrep),irp_ket=1,nirrep)
      end do
      do iSpin = 3, 1
      print *, 'istartmo(Irrep,Ispin) = ',(istartmo(irrep,iSpin),
     +                                     irrep=1,nirrep)
      end do
      print '(/)'

      return
c     end subroutine getaoinf
      end

