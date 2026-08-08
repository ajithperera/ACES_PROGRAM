










      Subroutine Cor_orbdens(Evecs,Docc,Dens,Norbs,Nbfns,Naobfns,Ispin)

      Implicit Double Precision(A-H,O-Z)

      Dimension Evecs(Naobfns,Nbfns), Dens(Naobfns*Naobfns,Norbs),
     +          Docc(Nbfns,Nbfns)

      Data Done,Dnull /1.0D0,0.0D0/

      Call Dzero(Dens,Naobfns*Naobfns*Norbs)

      Do Iorbs = 1, Norbs 
         Call Mk_corden(Evecs(1,Iorbs),Dens(1,Iorbs),Nbfns,Naobfns,
     +                  Docc(Iorbs,Iorbs))

      Enddo 

      Return
      End 


 
