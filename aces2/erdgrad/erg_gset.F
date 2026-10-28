      Subroutine Erg_gset(Ip,Iq,Ir,Is,Val,Ipx,N,Irg,Npr,Ioffg,Gso)
C
C Store Gam(pq,rs) = Gam(rs,pq) = Val (layout: Erg_gam).
C
      Implicit Double Precision (A-H,O-Z)
      Integer Ipx(N,N), Irg(N), Npr(8), Ioffg(8)
      Dimension Gso(*)

      Ig  = Ieor(Irg(Ip)-1,Irg(Iq)-1) + 1
      Ipq = Ipx(Ip,Iq)
      Irs = Ipx(Ir,Is)
      Gso(Ioffg(Ig) + (Ipq-1)*Npr(Ig) + Irs) = Val
      Gso(Ioffg(Ig) + (Irs-1)*Npr(Ig) + Ipq) = Val

      Return
      End
