











c THIS ROUTINE LOADS AO INTEGRALS FROM THE INTEGRAL FILE IIII AND
c TRANSFORMS THE FIRST INDEX.

      SUBROUTINE LOAD1(CMO,W,W2,dBuf,iBuf,iSymAO,nBas,
     &                 nFirst,nStart,nEnd,iSize,iOffAO,
     &                 iOff_Irrep_4AO,iOffI,
     &                 nSize,iLnBuf,iSpin,iUnit,bLast)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL bLast
c
      DIMENSION iBuf(iLnBuf),dBuf(iLnBuf),W(1),W2(1),CMO(1)
      DIMENSION nBas(8),nFirst(8),iSize(8),iOffAO(8),iOffI(8)
      DIMENSION nStart(8),nEnd(8),iOff_Irrep_4AO(8)
      DIMENSION iSymAO(100)
c
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,iAlone,iBitWd
      COMMON /FLAGS/  iFlags(100)
      COMMON /FLAGS2/ iFlags2(500)
      COMMON /AOOFST/ INDOCC(8,2)
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /VTINFO/ nPass1,nPass2,nPass3,nPass4,
     &                nLoad1,nLoad2,nLoad3,nLoad4,
     &                NWRIT1, NWRIT2, NWRIT3, NWRIT4,
     &                NWRIT1A,NWRIT2A,NWRIT3A,NWRIT4A,
     &                NWRIT1B,NWRIT2B,NWRIT3B,NWRIT4B

c ----------------------------------------------------------------------
      Indx(i, j) = j+(i*(i-1))/2
      NnP1o2(I)  = (i*(i+1))/2

   

c   o initialize W array
      call zero(W,nSize)

c   o increment Load1 counter
      nPass1 = nPass1 + 1

      nAOInt = 0
      NUT    = iLnBuf
      do while (NUT.eq.iLnBuf)

c   o read in integrals
      read(iUnit) dBuf, iBuf, NUT
      do int = 1, NUT

c      o extract X value and IX, JX, KX, LX indices
         X    = dBuf(int)
         iTmp = iBuf(int)
         IX   = iand(      iTmp,           iAlone)
         JX   = iand(ishft(iTmp,-  iBitWd),iAlone)
         KX   = iand(ishft(iTmp,-2*iBitWd),iAlone)
         LX   = iand(ishft(iTmp,-3*iBitWd),iAlone)

c      o get irrep of integral
         IrrepX = iSymAO(IX)
c      o process the occupied orbitals within IrrepX (if any)

C
         if (nStartX.le.nEndX) then

c         o get offsets within basis functions
            iTmp   = iOffAO(IrrepX)
            IX     = IX - iTmp
            JX     = JX - iTmp
            KX     = KX - iTmp
            LX     = LX - iTmp

c         o condition the packed indices
            I = max(IX,JX)
            J = min(IX,JX)
            K = max(KX,LX)
            L = min(KX,LX)
C
c         o scale the integral

c         o determine number of basis functions and offsets within IrrepX
C

c DETERMINE REDUNDANCY FACTOR FOR PLUGGING IN INTEGRALS
c THERE ARE A TOTAL OF EIGHT CONTRIBUTIONSa
c
c   (IJ|KL) (JI|KL) (IJ|LK) (JI|LK)
c   (KL|IJ) (KL|JI) (LK|IJ) (LK|JI)
c
            INDKL = Indx(K, L)

C
C
c
C
c        end if (nStartX.le.nEndX)
         end if

c     end do int = 1, NUT
      end do

      nAOInt = nAOInt + NUT
C
      Write(6,*)
      Print*, "The AO/MO(if old way) integrals"
      Ioff = 0
      Do Irrep = 1, Nirrep
         Nao = NBAS(IrreP)
         Nmo = NFirst(IrreP)
         KLPAIRS = NAO*(NAO+1)/2
         IJPAIRS = NAO*(NAO+1)/2
         Do KL = 1, KLPAIRS
            Print*, "The KL Pair = ", KL
            Write(*,'(4(1X,F12.7)))'), (W(Ioff + I), I=1, IJPAIRS)
CSSS            Call checksum2("AOLAD1", W(ioff), IJPAIRS)
            Ioff = Ioff + IJPAIRS
         End Do
      End Do 

c     end do while (NUT.eq.iLnBuf)
      end do

      if (NUT.eq.-1) then
         nAOInt = nAOInt + 1
         NUT = 0
      end if
      nLoad1 = nAOInt

      if (bLast) then
         if (iFlags(18).ge.3.or.
     &       iFlags(26).eq.1.or.
     &       iFlags(93).eq.2.or.
     &       (iSpin.eq.1.and.iFlags(11).gt.0).or.
     &       iFlags2(103).eq.1
     &      ) then
            close(unit=iUnit,status='KEEP')
         else
            close(unit=iUnit,status='KEEP')
         end if
      else
         rewind(iUnit)
         call locate(iUnit,'TWOELSUP')
      end if

      write(*,'(t3,a,i15,a)')
     &           '@LOAD1: ',nAOInt,' integrals read in from 
     &           file IIII'

      return
c     end subroutine load1
      end

