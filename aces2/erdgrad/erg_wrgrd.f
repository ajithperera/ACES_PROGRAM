










      Subroutine Erg_wrgrd(G,Natm,Natoms,Iz2a,Imap,Gc)
C
C Print the gradient and write the GRADIENT record. xvdint writes it
C in the computational center order: computational center i is ZMAT
C center Imap(i) (MAP2ZMAT), which is atom Iz2a(Imap(i)) of
C MOL_ERDOED (0 for a dummy center).
C
      Implicit Double Precision (A-H,O-Z)



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




      Dimension G(3,Natm), Gc(3,Natoms)
      Integer Iz2a(Natoms), Imap(Natoms)

      Write(6,"(/,a)") "  Molecular gradient (OED/ERD), computational"
     &                 //" order"
      Ncen = 0
      Do I = 1, Natoms
         Iz = Imap(I)
         Ia = 0
         If (Iz .Gt. 0) Ia = Iz2a(Iz)
         If (Ia .Gt. 0) Then
            Ncen = Ncen + 1
            Do K = 1, 3
               Gc(K,Ncen) = G(K,Ia)
            Enddo
            Write(6,"(2x,i4,3f20.10)") Ncen, (Gc(K,Ncen),K=1,3)
         Endif
      Enddo
      Call Putrec(20,'JOBARC','GRADIENT',3*Ncen*Iintfp,Gc)

      Return
      End
