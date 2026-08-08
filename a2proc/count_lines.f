












      Subroutine Count_lines(Ilines,Wrk)
      Implicit None

c     Maximum string length of terminal lines
      INTEGER LINELEN
      PARAMETER (LINELEN=80)
c io_units.par : begin

      integer    LuOut
      parameter (LuOut = 6)

      integer    LuErr
      parameter (LuErr = 6)

      integer    LuBasL
      parameter (LuBasL = 1)
      character*(*) BasFil
      parameter    (BasFil = 'BASINF')

      integer    LuVMol
      parameter (LuVMol = 3)
      character*(*) MolFil
      parameter    (MolFil = 'MOL')
      integer    LuAbi
      parameter (LuAbi = 3)
      character*(*) AbiFil
      parameter    (AbiFil = 'INP')
      integer    LuCad
      parameter (LuCad = 3)
      character*(*) CadFil
      parameter    (CadFil = 'CAD')

      integer    LuZ
      parameter (LuZ = 4)
      character*(*) ZFil
      parameter    (ZFil = 'ZMAT')

      integer    LuGrd
      parameter (LuGrd = 7)
      character*(*) GrdFil
      parameter    (GrdFil = 'GRD')

      integer    LuHsn
      parameter (LuHsn = 8)
      character*(*) HsnFil
      parameter    (HsnFil = 'FCM')

      integer    LuFrq
      parameter (LuFrq = 78)
      character*(*) FrqFil
      parameter    (FrqFil = 'FRQARC')

      integer    LuDone
      parameter (LuDone = 80)
      character*(*) DonFil
      parameter    (DonFil = 'JODADONE')

      integer    LuNucD
      parameter (LuNucD = 81)
      character*(*) NDFil
      parameter    (NDFil = 'NUCDIP')

      integer LuFiles
      parameter (LuFiles = 90)

c io_units.par : end
c     Maximum string length of absolute file names
      INTEGER FNAMELEN
      PARAMETER (FNAMELEN=80)

      Logical Btmp,Bopened,Bdone,bDelLastChar
      Character*(Linelen) Wrk,Sztmp,Wrk2
      Integer FNBLNK,LINBLNK,f_strpbrk
      Character*(linelen) TRMBLK
      Character*1 achar
      Integer i,j
      Integer Ilines 
      Integer iLastC, iOpenP, iCloseP, ValueSize, icycle, Ion
      Integer iTmpReg2,  iTmpReg3
      Character*1 czPercent, czAsterisk, czHash, czSpace, czFirstNLChar
      character*1 czTab

      czAsterisk = achar(42)
      czPercent  = achar(37)
      czSpace    = achar(32)
      czHash     = achar(35)
      czTab      = achar(9)

      Bdone  = .false.
      do while (.not.bDone)

c      o recondition the test string and point to the first char
c        (both in the namelist, with czFirstNLChar, and in the line, with i)
         iLastC = index(wrk,czHash)
         if (iLastC.eq.0) then
            wrk2 = trmblk(wrk)
            wrk  = wrk2
            i    = fnblnk(wrk)
         else
c         o pure comments will not terminate the namelist
            i = 1
            if (iLastC.eq.1) then
               wrk = ' '
            else
               wrk2 = trmblk(wrk(1:iLastC-1))
               wrk  = wrk2
            end if
         end if
         if (czFirstNLChar.eq.czSpace) then
            if (iLines.eq.1) then
               j = 6+fnblnk(wrk(7:))
               if (j.ne.6) czFirstNLChar = wrk(j:j)
            else
               j = fnblnk(wrk)
               if (j.ne.0) czFirstNLChar = wrk(j:j)
            end if
         end if

c      o test for a (totally) blank line or a new namelist after the first line
c        (EOF automatically jumps)
         if (i.ne.0) then
            bTmp = (wrk(i:i).eq.czAsterisk)
         else
            bTmp = .false.
         end if
         if ((i.eq.0).or.((iLines.ne.1).and.bTmp)) then
c         o read too far -> jump back and quit counting
            backspace(luz)
            iLines = iLines - 1
            bDone  = .true.
c         o Did we accidentally record an asterisk for the first char of
c           an empty namelist?
            if (czFirstNLChar.eq.czAsterisk) czFirstNLChar = czSpace

         else
c         o wrk has no comments and may be empty (denoting a pure comment
c           in the namelist, which won't terminate parsing)
            iCloseP = index(wrk,')')
            if (iCloseP.eq.0) then
c            o no ')' -> keep going
c              (NOTE: we will check for multiple, unmatched '(' during the
c                     key-value tokenization)
               read(luz,'(a)',end=8012) wrk
               iLines = iLines + 1

c           else if (iCloseP.ne.0) then
            else
c            o measure the parenthetical level (j) from the beginning (i)
               i = 0
               if (iLines.eq.1) then
                  j = 0
               else
c               o This may seem strange, but we have to convince the logic
c                 after the parenthetical count that we are still in the
c                 namelist. If we encounter a closing ')', then we will test
c                 czFirstNLChar to see if it is '('.
                  j = 1
               end if

               iTmpReg2 = f_strpbrk(wrk,'()')
               do while (iTmpReg2.ne.0)

                  i = i + iTmpReg2
                  if (wrk(i:i).eq.'(') then
                     j = j + 1
                  else
                     j = j - 1
                  end if

                  if (i.lt.len(wrk)) then
                     iTmpReg2 = f_strpbrk(wrk(i+1:),'()')
                  else
                     iTmpReg2 = 0
                  end if

c              end do while ([have more parentheses])
               end do
c            o no more parentheses
               if (j.eq.0) then
c               o we should be done with the namelist, but is there any data
c                 after the last ')'?
                  if (i.lt.len(wrk)) then
                     iTmpReg3 = fnblnk(wrk(i+1:))
                     if (iTmpReg3.ne.0) then
                        print *, '@-count_lines: ERROR - ',
     &         'data found after the terminal ")" in the ACES namelist:'
                        print *, wrk
                        call errex
                     end if
                  end if
c               o since a final ')' was found, make sure the namelist was
c                 started with a '('
                  if (czFirstNLChar.ne.'(') then
                  print *, '@-count_lines: ERROR - ")" terminates the ',
     &                     'namelist but "(" does not initiate it.'
                  call errex
                  end if
                  bDelLastChar = .true.
                  bDone = .true.

c              else if (j.ne.0) then
               else

                  if (j.eq.1) then
c                  o we should still be somewhere in the namelist
                     read(luz,'(a)',end=8012) wrk
                     iLines = iLines + 1
                  else
                     print *, '@GTFLGS: Error - ',
     &                        'open parentheses in the ACES namelist:'
                     print *, wrk
                     call errex
                  end if

c              end if (j.eq.0)
               end if

c           end if (iCloseP.eq.0)
            end if

c        end if ([blank line].or.[new namelist])
         end if
c     end do while (.not.bDone)
      end do
 8012 Continue 

c   o If the namelist is empty, then stop.
      if ((iLines.eq.1).and.(czFirstNLChar.eq.czSpace)) then
         print *, 'buildzmat: ACES will not run with an empty namelist.'
         call errex
      end if

      Return
      End 

