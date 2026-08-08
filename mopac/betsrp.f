










C 
C    BETSRP unit (Ivan Rossi - April 94)                               
C
C include in mopcvar.par the two following paramrters:
C
C     MXATSP : MAXIMUM NUMBER OF atomic species in the system
C     MXSRPB : MAXIMUM NUMBER OF SRP special beta allowed
*
*     COMMON BLOCKS FOR SRP (Ivan Rossi - April 94)                    IR0494
*
*     COMMON /SRPI/ IBTPTR(107), NATPTR(MXATSP), NATSP
*    *       /SRPL/ ISSRP
*    *       /SRPR/ BETSS(MAXBET), BETSP(MXATSP,MXATSP), BETPP(MAXBET)
*     LOGICAL ISSRP
c
      SUBROUTINE INIBET
*
*  Initialize arrays and pointers for special SRP BETAs (Ivan Rossi- April '94)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C  SIZES.i
C
C  This file is based on the MOPAC include file SIZES.
C
************************************************************************
*   THIS FILE CONTAINS ALL THE ARRAY SIZES FOR USE IN MOPAC.            
*                                                                       
*     THERE ARE ONLY  PARAMETERS THAT THE PROGRAMMER NEED SET:          
*     MAXHEV = MAXIMUM NUMBER OF HEAVY ATOMS (HEAVY: NON-HYDROGEN ATOMS)
*     MAXLIT = MAXIMUM NUMBER OF HYDROGEN ATOMS.                        
*     MAXTIM = DEFAULT TIME FOR A JOB. (SECONDS)                        
*     MAXDMP = DEFAULT TIME FOR AUTOMATIC RESTART FILE GENERATION (SECS)
*     MXATSP = MAXIMUM NUMBER OF ATOMIC SPECIES IN THE SYSTEM           IR0494
*                                                                       
      PARAMETER (MAXHEV=50,   MAXLIT=50)
      PARAMETER (MAXTIM=3600, MAXDMP=3600)
      PARAMETER (MXATSP=10)                                             IR0494
*                                                                       
************************************************************************
*                                                                       
*   THE FOLLOWING CODE DOES NOT NEED TO BE ALTERED BY THE PROGRAMMER    
*                                                                       
************************************************************************
*                                                                       
*    ALL OTHER PARAMETERS ARE DERIVED FUNCTIONS OF THESE TWO PARAMETERS 
*                                                                       
*      NAME                   DEFINITION                                
*     NUMATM         MAXIMUM NUMBER OF ATOMS ALLOWED.                   
*     MAXORB         MAXIMUM NUMBER OF ORBITALS ALLOWED.                
*     MAXPAR         MAXIMUM NUMBER OF PARAMETERS FOR OPTIMISATION.     
*     N2ELEC         MAXIMUM NUMBER OF TWO ELECTRON INTEGRALS ALLOWED.  
*     MPACK          AREA OF LOWER HALF TRIANGLE OF DENSITY MATRIX.     
*     MORB2          SQUARE OF THE MAXIMUM NUMBER OF ORBITALS ALLOWED.  
*     MAXHES         AREA OF HESSIAN MATRIX                             
*     MAXDMO         MAXIMUM DIAGONALIZABLE MATRIX ORDER                IR0494
*     MXSRPB         MAXIMUM NUMBER OF SRP SPECIAL BETA ALLOWED         IR0494
************************************************************************
      PARAMETER (VERSON=5.07D0)
      PARAMETER (NUMATM=MAXHEV+MAXLIT)
      PARAMETER (MAXORB=4*MAXHEV+MAXLIT)
      PARAMETER (MAXPAR=3*NUMATM)
      PARAMETER (MAXBIG=MAXORB*MAXORB*2)
      PARAMETER (N2ELEC=2*(50*MAXHEV*(MAXHEV-1)+10*MAXHEV*MAXLIT
     +                     +(MAXLIT*(MAXLIT-1))/2))
      PARAMETER (MAXHES=(MAXPAR*(MAXPAR+1))/2,MORB2=MAXORB**2)
      PARAMETER (MPACK=(MAXORB*(MAXORB+1))/2)
      PARAMETER (MAXBET=(MXATSP*(MXATSP+1))/2)                          IR0494
************************************************************************
*   FOR SHORT VERSION USE LINE WITH NMECI=1, FOR LONG VERSION USE LINE  
*   WITH NMECI=10                                                       
************************************************************************
C     PARAMETER (NMECI=10,  NPULAY=MPACK)
C     PARAMETER (NMECI=1,   NPULAY=1)
      PARAMETER (NMECI=10,  NPULAY=MPACK)
C WARNING : MAXDMO have to be >=  max(MAXORB,NMECI**2)                   IR0494
      PARAMETER (MAXDMO=250)                                             IR0494
************************************************************************
C
      COMMON /SRPI/ IBTPTR(107), NATPTR(MXATSP), NATSP
     *       /SRPR/ BETSS(MAXBET), BETSP(MXATSP,MXATSP), BETPP(MAXBET)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
*
*     COMMON BLOCKS FOR STANDARD BETAs
*
      COMMON /BETAS / BETAS(107),BETAP(107),BETAD(107)
c
c     Initialize pointers to RSP betas
c
      NATSP=0
      do 5 i=1,107
   5     IBTPTR(i)=0
      DO 10 i=1,NATOMS
        if (LABELS(i) .eq. 99 .or. LABELS(i) .eq.107) goto 10
        if (IBTPTR(LABELS(i)) .eq. 0) then
           NATSP=NATSP+1
           if (NATSP .gt. MXATSP) then
             Write(6,'("***  The parameter MXATSP is too SMALL ! ***")')    
             Write(6,'("*Increase it in mopcvar.par and recompile *")')      
             STOP
           ENDIF
           NATPTR(NATSP)=LABELS(I)
           IBTPTR(LABELS(i))=NATSP
        endif
 10   continue
        do 50 j=1,NATSP
          do 50 i=j,NATSP
          BSS= 0.5d0*(BETAS(NATPTR(i))+BETAS(NATPTR(j)))
          BPP= 0.5d0*(BETAP(NATPTR(i))+BETAP(NATPTR(j)))
          call SETBET( NATPTR(i), NATPTR(j), 'BETSS', BSS )
          call SETBET( NATPTR(i), NATPTR(j), 'BETPP', BPP )
 50     continue
        do 60 j=1,NATSP
          do 60 i=1,NATSP
          BSP= 0.5d0*(BETAS(NATPTR(i))+BETAP(NATPTR(j)))
          call SETBET( NATPTR(i), NATPTR(j), 'BETSP', BSP )
 60     continue
      RETURN
      END
*
      SUBROUTINE SETBET( NATM1, NATM2, BETTYP, BETVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*5 BETTYP
      DOUBLE PRECISION BETVAL
* 
*  Set the value of the special SRP BETAs ( Ivan Rossi - April '94 )
*      WARNING: BETSP(x,y) is DIFFERENT from BETSP(y,x)
*
*   INPUT:
*     NATM1, NATM2 : Atomic numbers of the atom couple
*     BETTYP : Name of the Beta parameter type (UPPERCASE!)
*     BETVAL : value of the BETA to set
*
C
C  SIZES.i
C
C  This file is based on the MOPAC include file SIZES.
C
************************************************************************
*   THIS FILE CONTAINS ALL THE ARRAY SIZES FOR USE IN MOPAC.            
*                                                                       
*     THERE ARE ONLY  PARAMETERS THAT THE PROGRAMMER NEED SET:          
*     MAXHEV = MAXIMUM NUMBER OF HEAVY ATOMS (HEAVY: NON-HYDROGEN ATOMS)
*     MAXLIT = MAXIMUM NUMBER OF HYDROGEN ATOMS.                        
*     MAXTIM = DEFAULT TIME FOR A JOB. (SECONDS)                        
*     MAXDMP = DEFAULT TIME FOR AUTOMATIC RESTART FILE GENERATION (SECS)
*     MXATSP = MAXIMUM NUMBER OF ATOMIC SPECIES IN THE SYSTEM           IR0494
*                                                                       
      PARAMETER (MAXHEV=50,   MAXLIT=50)
      PARAMETER (MAXTIM=3600, MAXDMP=3600)
      PARAMETER (MXATSP=10)                                             IR0494
*                                                                       
************************************************************************
*                                                                       
*   THE FOLLOWING CODE DOES NOT NEED TO BE ALTERED BY THE PROGRAMMER    
*                                                                       
************************************************************************
*                                                                       
*    ALL OTHER PARAMETERS ARE DERIVED FUNCTIONS OF THESE TWO PARAMETERS 
*                                                                       
*      NAME                   DEFINITION                                
*     NUMATM         MAXIMUM NUMBER OF ATOMS ALLOWED.                   
*     MAXORB         MAXIMUM NUMBER OF ORBITALS ALLOWED.                
*     MAXPAR         MAXIMUM NUMBER OF PARAMETERS FOR OPTIMISATION.     
*     N2ELEC         MAXIMUM NUMBER OF TWO ELECTRON INTEGRALS ALLOWED.  
*     MPACK          AREA OF LOWER HALF TRIANGLE OF DENSITY MATRIX.     
*     MORB2          SQUARE OF THE MAXIMUM NUMBER OF ORBITALS ALLOWED.  
*     MAXHES         AREA OF HESSIAN MATRIX                             
*     MAXDMO         MAXIMUM DIAGONALIZABLE MATRIX ORDER                IR0494
*     MXSRPB         MAXIMUM NUMBER OF SRP SPECIAL BETA ALLOWED         IR0494
************************************************************************
      PARAMETER (VERSON=5.07D0)
      PARAMETER (NUMATM=MAXHEV+MAXLIT)
      PARAMETER (MAXORB=4*MAXHEV+MAXLIT)
      PARAMETER (MAXPAR=3*NUMATM)
      PARAMETER (MAXBIG=MAXORB*MAXORB*2)
      PARAMETER (N2ELEC=2*(50*MAXHEV*(MAXHEV-1)+10*MAXHEV*MAXLIT
     +                     +(MAXLIT*(MAXLIT-1))/2))
      PARAMETER (MAXHES=(MAXPAR*(MAXPAR+1))/2,MORB2=MAXORB**2)
      PARAMETER (MPACK=(MAXORB*(MAXORB+1))/2)
      PARAMETER (MAXBET=(MXATSP*(MXATSP+1))/2)                          IR0494
************************************************************************
*   FOR SHORT VERSION USE LINE WITH NMECI=1, FOR LONG VERSION USE LINE  
*   WITH NMECI=10                                                       
************************************************************************
C     PARAMETER (NMECI=10,  NPULAY=MPACK)
C     PARAMETER (NMECI=1,   NPULAY=1)
      PARAMETER (NMECI=10,  NPULAY=MPACK)
C WARNING : MAXDMO have to be >=  max(MAXORB,NMECI**2)                   IR0494
      PARAMETER (MAXDMO=250)                                             IR0494
************************************************************************
C
      COMMON /SRPI/ IBTPTR(107), NATPTR(MXATSP), NATSP
     *       /SRPR/ BETSS(MAXBET), BETSP(MXATSP,MXATSP), BETPP(MAXBET)
c
      DIMENSION NATM(2)
C
        NATM(1)=NATM1
        NATM(2)=NATM2
C
C   INPUT Error checking
C
      do 10 j=1,2
        if (IBTPTR(NATM(j)) .eq. 0) then 
          WRITE(6,'("SETBET: Atom type",I4," NOT present!")') NATM(j)
cmgc          STOP
          return
        endif
 10   continue
c
c  Calculate address of beta (packed lower simmetric matrix) and set value
c
     
      IROW=IBTPTR(NATM(1))
      JCOL=IBTPTR(NATM(2))
      if( JCOL .gt. IROW ) then
        IROW=IBTPTR(NATM(2))
        JCOL=IBTPTR(NATM(1))
      ENDIF
      j=NATSP*(JCOL-1)+IROW -(JCOL*(JCOL-1))/2
      if( BETTYP .eq. 'BETPP') then
        BETPP(j)=BETVAL
      else if( BETTYP .eq. 'BETSP') then
        BETSP(IBTPTR(NATM1),IBTPTR(NATM2))=BETVAL
      else if( BETTYP .eq. 'BETSS') then
        BETSS(j)=BETVAL
      else
        write(6,'("SETBET: Wrong beta type",A5)') BETTYP
        STOP
      endif
      return
      END
*
      DOUBLE PRECISION FUNCTION GETBET( NATM1, NATM2, BETTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*5 BETTYP
* 
*  Get the value of the special SRP BETAs ( Ivan Rossi - April '94 )
*      WARNING: BETSP(x,y) is DIFFERENT from BETSP(y,x)
*
*   INPUT:
*     NATM1, NATM2 : Atomic numbers of the atom couple
*     BETTYP : Name of the Beta parameter type (UPPERCASE!)
*   OUTPUT:
*     Value of the SRP BETA requested
*
C
C  SIZES.i
C
C  This file is based on the MOPAC include file SIZES.
C
************************************************************************
*   THIS FILE CONTAINS ALL THE ARRAY SIZES FOR USE IN MOPAC.            
*                                                                       
*     THERE ARE ONLY  PARAMETERS THAT THE PROGRAMMER NEED SET:          
*     MAXHEV = MAXIMUM NUMBER OF HEAVY ATOMS (HEAVY: NON-HYDROGEN ATOMS)
*     MAXLIT = MAXIMUM NUMBER OF HYDROGEN ATOMS.                        
*     MAXTIM = DEFAULT TIME FOR A JOB. (SECONDS)                        
*     MAXDMP = DEFAULT TIME FOR AUTOMATIC RESTART FILE GENERATION (SECS)
*     MXATSP = MAXIMUM NUMBER OF ATOMIC SPECIES IN THE SYSTEM           IR0494
*                                                                       
      PARAMETER (MAXHEV=50,   MAXLIT=50)
      PARAMETER (MAXTIM=3600, MAXDMP=3600)
      PARAMETER (MXATSP=10)                                             IR0494
*                                                                       
************************************************************************
*                                                                       
*   THE FOLLOWING CODE DOES NOT NEED TO BE ALTERED BY THE PROGRAMMER    
*                                                                       
************************************************************************
*                                                                       
*    ALL OTHER PARAMETERS ARE DERIVED FUNCTIONS OF THESE TWO PARAMETERS 
*                                                                       
*      NAME                   DEFINITION                                
*     NUMATM         MAXIMUM NUMBER OF ATOMS ALLOWED.                   
*     MAXORB         MAXIMUM NUMBER OF ORBITALS ALLOWED.                
*     MAXPAR         MAXIMUM NUMBER OF PARAMETERS FOR OPTIMISATION.     
*     N2ELEC         MAXIMUM NUMBER OF TWO ELECTRON INTEGRALS ALLOWED.  
*     MPACK          AREA OF LOWER HALF TRIANGLE OF DENSITY MATRIX.     
*     MORB2          SQUARE OF THE MAXIMUM NUMBER OF ORBITALS ALLOWED.  
*     MAXHES         AREA OF HESSIAN MATRIX                             
*     MAXDMO         MAXIMUM DIAGONALIZABLE MATRIX ORDER                IR0494
*     MXSRPB         MAXIMUM NUMBER OF SRP SPECIAL BETA ALLOWED         IR0494
************************************************************************
      PARAMETER (VERSON=5.07D0)
      PARAMETER (NUMATM=MAXHEV+MAXLIT)
      PARAMETER (MAXORB=4*MAXHEV+MAXLIT)
      PARAMETER (MAXPAR=3*NUMATM)
      PARAMETER (MAXBIG=MAXORB*MAXORB*2)
      PARAMETER (N2ELEC=2*(50*MAXHEV*(MAXHEV-1)+10*MAXHEV*MAXLIT
     +                     +(MAXLIT*(MAXLIT-1))/2))
      PARAMETER (MAXHES=(MAXPAR*(MAXPAR+1))/2,MORB2=MAXORB**2)
      PARAMETER (MPACK=(MAXORB*(MAXORB+1))/2)
      PARAMETER (MAXBET=(MXATSP*(MXATSP+1))/2)                          IR0494
************************************************************************
*   FOR SHORT VERSION USE LINE WITH NMECI=1, FOR LONG VERSION USE LINE  
*   WITH NMECI=10                                                       
************************************************************************
C     PARAMETER (NMECI=10,  NPULAY=MPACK)
C     PARAMETER (NMECI=1,   NPULAY=1)
      PARAMETER (NMECI=10,  NPULAY=MPACK)
C WARNING : MAXDMO have to be >=  max(MAXORB,NMECI**2)                   IR0494
      PARAMETER (MAXDMO=250)                                             IR0494
************************************************************************
C
      COMMON /SRPI/ IBTPTR(107), NATPTR(MXATSP), NATSP
     1       /SRPR/ BETSS(MAXBET), BETSP(MXATSP,MXATSP), BETPP(MAXBET)
c
      DIMENSION NATM(2)
C
        NATM(1)=NATM1
        NATM(2)=NATM2
C
C   INPUT Error checking
C
      do 10 j=1,2
        if (IBTPTR(NATM(j)) .eq. 0) then 
          WRITE(6,'("GETBET: Atom type",I4," NOT present!")') NATM(j)
          STOP
        endif
 10   continue
c
c  Calculate address of beta and set value
c
      IROW=IBTPTR(NATM(1))
      JCOL=IBTPTR(NATM(2))
      if( JCOL .gt. IROW ) then
        IROW=IBTPTR(NATM(2))
        jCOL=IBTPTR(NATM(1))
      ENDIF
      j=NATSP*(JCOL-1)+IROW -(JCOL*(JCOL-1))/2
      if( BETTYP .eq. 'BETPP') then
        GETBET=BETPP(j)
      else if( BETTYP .eq. 'BETSP') then
        GETBET=BETSP(IBTPTR(NATM1), IBTPTR(NATM2))
      else if( BETTYP .eq. 'BETSS') then
        GETBET=BETSS(j)
      else
        write(6,'("GETBET: Wrong beta type",A5)') BETTYP
          STOP
      endif
      return
      END
