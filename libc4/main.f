










      Program libc4
 
      Implicit double precision(A-H,O-Z)
      Logical T2_fine



c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end


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




      Data Ione, Dzero /1, 0.0D0/
     
      Call Crapsi_c4(Icore,I0,Icrsiz,Iuhf,0)

      Write(6,"(10x,a)") " -----------Interface to CFOUR------"
      Write(6,*) 

      Call Getrec_c4(20,"JOBARC","NATOMS ", Ione, Natoms)
      Call Getrec_c4(20,"JOBACR","COORD  ", 3*Natoms,Icore(I0))

      Write(6,"(a,1x,I2)") " Hello CFOUR world!"
 
      Write(6,"(a,1x,I3)") " The number of  atoms: ", Natoms 
      Write(6,"(a)")       " Molecular geomery   :" 
      Call output(Icore(I0),1,Natoms,1,3,Natoms,3,1)
      
      Distance=Ddot(3*Natoms,Icore(I0),1,Icore(I0),1)
      If (Natoms .Gt. 0 .And. Distance .Gt. Dzero) Then 
         Write(6,*)
         Write(6,"(a)") " Communication with JOBARC is established!"
      Endif 
     
      Call Getrec_C4(20,"JOBARC",'SCFENEG ',Ione*IINTFP,ESCF) 
      Call Getrec_C4(20,"JOBARC",'TOTENERG',Ione*IINTFP,ETOT) 

      Ediff = Etot - Escf 
   
      If (Ediff .Ne. Dzero) Then
         Call Check_t2(Icore(I0),Icrsiz,Iuhf,T2_Fine)
      Endif 

      If (T2_fine) Then
         Write(6,*)
         Write(6,"(a)") " Communication with MOIO is established!"
      Endif 

      Write(6,*)
      Call Crapso_c4()

CSSS      Call Aces_io_fin_c4
CSSS      Call aces_ja_fin_c4

      Stop
      End
