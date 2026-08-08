










      Subroutine Psi4dbg_form_dpq(Dcpq,Dspq,Doo,Dvv,Dvo,Work,Maxcor,
     +                            Nocc,Nvrt,Nbas,Iuhf,Non_hf)

      Implicit Double Precision (A-H,O-Z)
      Logical pCCD,CCD,LCCD
      Logical Non_hf

      Dimension Dcpq(Nbas,Nbas)
      Dimension Dspq(Nbas,Nbas)
      Dimension Dvv(Nvrt,Nvrt) 
      Dimension Doo(Nocc,Nocc) 
      Dimension Dvo(Nvrt,Nocc) 
      Dimension Work(Maxcor)



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



      Common /Meth/pCCD,CCD,LCCD

      Data Ione,Inull,One,Two /1,0,1.0D0,2.0D0/

      Call Psi4dbg_frmful(Doo,Dvv,Dvo,Dcpq,Work,Maxcor,Nbas,0,
     +                    Non_hf)

      Scfocc = One
CSSS      Call Dscal(Nbas*Nbas,Two,Dcpq,1)

C This density matrix is post-HF ordered (i.e. all occupied for all irreps
C fist then followed by virtuals in all irreps).

      If (pCCD) Write(6,"(a)") "The pCCD density matrix"
      If (CCD)  Write(6,"(a)") "The CCD density matrix"
      If (LCCD) Write(6,"(a)") "The LCCD density matrix"
      call output(Dcpq,1,Nbas,1,Nbas,Nbas,Nbas,1)


      Return
      End
      
