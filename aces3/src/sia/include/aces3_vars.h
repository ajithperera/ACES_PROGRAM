C----------------------------------------------------------------
C This include file combined with the SIP set_aces3_vars can be
C used to set variables from SIAL code. 
C 
      Integer AAAA, BBBB, ABAB
      Integer RHS, LHS 
      Integer OCC, VRT
      Integer No_core_sites 
      Logical Ee,Ip,Ea,Dip,Dea 
      Character*8 Exc_block 
      common /Spin_types/ AAAA, BBBB, ABAB
      common /Side/ RHS, LHS
      common /Type/ OCC,VRT
      Common /Direct/Exc_block 
      Common /Coresites/No_core_sites 
      Common /eomtype/Ee,Ip,Ea,Dip,Dea 

C----------------------------------------------------------------


