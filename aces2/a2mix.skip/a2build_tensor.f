










      SUBROUTINE A2BUILD_TENSOR(WORK, MAXCOR, PRDUTINT, TENSOR_ATR,
     &                          NMBR_OF_PERTS, NAOBFNS, NATOMS, 
     &                          IBEGIN_P_DENS, IBEGIN_P_OPRTS,
     &                          IBEGIN_AO_OVRLP, IBEGIN_MO_OVRLP, 
     &                          IBEGIN_MO_VECTS, IBEGIN_MEM_LEFT,
     &                          JPERT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION WORK(MAXCOR), PRDUTINT(NAOBFNS, NAOBFNS),
     &          TENSOR_ATR(NMBR_OF_PERTS)
C
      DATA IZR0, IONE /0, 1/
C
C

      DO IPERT = 1, NMBR_OF_PERTS
C
         I_TOT_PDEN_N = IBEGIN_P_DENS + (IPERT-1)*
     &                  NAOBFNS*NAOBFNS

         TOT_PDEN_N = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                     WORK(I_TOT_PDEN_N), 1)
C
CSSS         TENSOR_ATR(IPERT, JPERT) = TOT_PDEN_N
C
CSSS         DO JPERT = 1, NMBR_OF_PERTS
C
            I_TOT_POPRT_N = IBEGIN_P_OPRTS + (IPERT-1)*
     &                      NAOBFNS*NAOBFNS
C
            TOT_OPRT_N = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                        WORK(I_TOT_POPRT_N), 1)

            TOT_AO_OVRLP = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                        WORK(IBEGIN_AO_OVRLP), 1)

            TOT_MO_OVRLP = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                        WORK(IBEGIN_MO_OVRLP), 1)
C#ifdef 1
      Write(6,*)
      write(6,'(a,a)') "Computing the perturbed total/spin",
     &                    " density and others for props: "
      Write(6,*) 
      Write(*,'((2I5, 1x),F10.5)')ipert,jpert,TOT_PDEN_N
      Write(6,*)
C#endif
C
            TENSOR_ATR(IPERT) = TOT_PDEN_N
C
            write(6,*)
            Write(6,"(2I5,F10.5)") ipert, jpert,TENSOR_ATR(ipert)
CSSS         ENDDO
C
      ENDDO
C
      RETURN
      END 
C  
