










      Subroutine Orb_kint_enegs(Orbdens,Docc,Work,Maxcor,Naobfns,Nbfns,
     +                          NOrbs)

      Implicit Double Precision(A-H,O-Z)
 
      Dimension Orbdens(Naobfns*Naobfns,Norbs),Docc(Norbs,Norbs)
      Dimension Work(Maxcor)
      Dimension Buf(600), Ibuf(600)
      Logical IIII_present,VPOUT_present
      Character*80 Fname 
      Character*8 String
      Dimension Itriln(9)
      Dimension Isqrln(9)
      Dimension Itriof(9)
      Dimension Isqrof(9)
      Dimension Irpoff(9)
      Dimension Nbfirr(8)



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

      Data Ione,Dnull,Half /1,0.0D0,0.50D0/
     
      Indx2(I,J,N) = I+(J-1)*N

      Call Gfname("IIII   ",Fname,Length) 
      Iunit = 10
      Inquire(File=Fname(1:Length),Exist=IIII_present)
      If (IIII_present) Then
         Open (Unit=Iunit,File=Fname(1:Length),Form="Unformatted",
     +         Access="sequential")
      Endif 

      Call Getrec(20,"JOBARC",'NUMBASIR',Nirrep,Nbfirr)

      Irpoff(1) = 0
      Do Irrep = 1, Nirrep
         Itriln(Irrep)   = Nbfirr(Irrep)*(Nbfirr(Irrep)+1)/2
         Isqrln(Irrep)   = Nbfirr(Irrep)*(Nbfirr(Irrep))
         Irpoff(Irrep+1) = Irpoff(Irrep) + Nbfirr(Irrep)
      Enddo

      Itriof(1) = 1
      Isqrof(1) = 1
      Irpoff(1) = 0
      Do Irrep = 1, Nirrep-1
         Itriof(Irrep+1) = Itriof(Irrep) + Itriln(Irrep)
         Isqrof(Irrep+1) = Isqrof(Irrep) + Isqrln(Irrep)
         Irpoff(Irrep+1) = Irpoff(Irrep) + Nbfirr(Irrep)
      Enddo

      Naobfns2 = Naobfns*Naobfns 
      Ldim  = Naobfns*(Naobfns+1)/2
      Ldim2 = Naobfns2
      I000 = Ione
      I010 = I000 + Ldim
      I020 = I010 + Naobfns2
      Iend = I020 + Naobfns 
      If (Iend+Ldim2 .Gt. Maxcor) Call Insmem("orb_kint_enegs",
     +                                         Iend+Ldim2,Maxcor)

      Call Dzero(Work(I000),Ldim)
      Call Dzero(Work(I010),Naobfns2)

      Call Locate(Iunit,"KINETINT")
      Ilnbuf = 600
      Marker = Ilnbuf
      Do While (Marker .Eq. Ilnbuf)
         Read(Iunit), Buf,Ibuf,Marker
         Do Int = 1, Marker
            Work(I000-1+Ibuf(Int)) = Buf(Int)
         Enddo
      Enddo

      Call Mkfull_kints(Work(I010),Work(I000),Work(Iend),Ldim,Ldim2,
     +                  Naobfns,Ione,Nbfirr,Irpoff,Itriof,Nirrep)
      
      Write(6,1)
1     Format(10x,(63("-")))
      Write(6,"(10x,a,5x,a)") "Orbital Kinetic Energy (in a.u.)", 
     +                        "Natural orbital occ. Nums."
      Write(6,1)
      Do Iorbs = 1, Norbs
         Work(I020-1+Iorbs) = Ddot(Naobfns2,Work(I010),1,
     +                        Orbdens(1,Iorbs),1)
         Write(6,"(15x,F15.8,15X,F15.8)") Work(I020-1+Iorbs)*Half,
     +                                    Docc(Iorbs,Iorbs)
      Enddo 
      Write(6,1)
      Close(Iunit)


      Return
      End 
