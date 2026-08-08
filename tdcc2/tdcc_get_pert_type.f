










      Subroutine Tdcc_get_pert_type(Work,Memleft,Nbfns,Label,Isympert,
     +                              Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Iaopop(8),Isympert(3)
      Character*8 Label(3)
      Dimension Iscr(Nbfns)



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

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end

      Parameter (Tol = 1.0D-05) 
      Label(1) = 'DIPOLE_X'
      Label(2) = 'DIPOLE_Y'
      Label(3) = 'DIPOLE_Z'

      Length = Nbfns * (Nbfns+1)/2

      Call Getrec(20,"JOBARC","NUMBASIR",Nirrep,Iaopop)

      Write(6,*)
      Write(6,"(a,8(1x,i2))") "The Iaopop:",
     +                         (Iaopop(i),i=1,Nirrep)
      Do Npert = 1, 3
         Call Getrec(20,"JOBARC",Label(Npert),Length*Iintfp,Work)
      Write(6,*)
      write(6,"(a)") "The dipole integrals"
      write(6,"(6(1x,F12.6))") (work(i),i=1,Length)
      Write(6,*)
         Ioff = 0
         Do Irrep = 1, Nirrep 
            Numirr = Iaopop(Irrep)
            Do I = 1, Numirr
               Ioff = Ioff + 1
               Iscr(Ioff) = Irrep
            Enddo
         Enddo 
             
         Ithru = 0
         Do Index1 = 1, Nbfns 
            Do Index2 = 1, Index1
               Ithru = Ithru + 1
               If (Dabs(Work(Ithru)) .GT. Tol) Then
                  Isym1 = Iscr(Index1)
                  Isym2 = Iscr(Index2)
                  Isym  = DIRPRD(Isym1,Isym2)
               Endif
            Enddo
         Enddo 

         Isympert(Npert) = Isym

      Enddo
      Write(6,*) 
      Write(6,"(a,3(1x,I2))") "The symmetry of dipole perturbations: ",
     +     (isympert(i),i=1,3)

      Return
      End
