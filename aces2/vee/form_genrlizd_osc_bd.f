










      Subroutine Form_genrlizd_osc_bd(Work,Maxcor,Iuhf,Eigval)

      Implicit Double Precision(A-H,O-Z)

      Character*8 Label_xyz(3),Label
      Character*2 Suffix(32)

      Dimension Tm(3,4,8,4,8), Eigval(100,8)
      Dimension Osc(1024)
      Dimension Work(Maxcor)

c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end


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



       
      Common /Calcinfo/Nroot(8)

      Data Label_xyz/'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
      Data Suffix /"01","02","03","04","05","06","07","08","09","10",
     +             "11","12","13","14","15","16","17","18","19","20",
     +             "21","22","23","24","25","26","27","28","29","30",
     +             "31","32"/
      Data Ione, Ithree, Fact, Done, Dnull /1,3,0.66666666666666D0,
     +                                      1.0D0,0.0D0/

      write(6,*) "Entering form_genrlizd_osc_bd"
      Call Getrec(20,"JOBARC","NBASTOT ",Ione,Nao)
      Nmo  = Nocco(1) + Nvrto(1)
      Nmo2 = Nmo*Nmo
      Nao2 = Nao*Nao

      Iscr1 = Ione
      Iscr2 = Iscr1 + Nao2
      Iscr3 = Iscr2 + Nao2
      Iend  = Iscr3 + Nao2

      Iqm6  = Iend
      Iom10 = Iqm6  + 1024*6
      Iam3  = Iom10 + 1024*10
      Idm3  = Iam3  + 1024*3
      Imq9  = Idm3  + 1024*3
      Iend2 = Imq9  + 1024*9

      Iqm6_2  = Iend2
      Iom10_2 = Iqm6_2  + 1048576*9
      Iam3_2  = Iom10_2 + 1048576*10
      Idm3_2  = Iam3_2  + 1048576*3
      Imq9_2  = Idm3_2  + 1048576*3
      Iend3   = Imq9_2  + 1048576*27

      Iosc_qm = Iend3
      Iosc_om = Iosc_qm + 1024
      Iosc_am = Iosc_om + 1024
      Iosc_dm = Iosc_am + 1024
      Iosc_mq = Iosc_dm + 1024
      Iend    = Iosc_mq + 1024
      If (Iend .Ge. Maxcor) Call Insmem("form_genrlizd_osc_bd",
     +                                   Iend,Maxor)
      Maxcor = Maxcor - Iend 
       
      Len_drec = Nmo*(Nmo+1)/2
      Call Dzero(Work(Iqm6), 1024*6)
      Call Dzero(Work(Iom10),1024*10)
      Call Dzero(Work(Iam3), 1024*3)
      Call Dzero(Work(Idm3), 1024*3)
      Call Dzero(Work(Imq9), 1024*9)

      Ndone = 0
      Do Irrepx = 1, Nirrep
         Do Iroot = 1, Nroot(Irrepx)
            Do Irrepy = 1, Nirrep
               Do Jroot = 1, Nroot(Irrepy)
                   Do Ispin = 1, Iuhf+1
                      Ndone = Ndone + 1
                      Label = "GTDENS"//Suffix(Ndone)
                      Call Getrec(20,"JOBARC",LABEL,Nao2*Iintfp, 
     +                            Work(Iscr3))
      Write(6,"(a,2(1x,i2))") "Irrepx,Irrepy : ",
     +                         Irrepx,Irrepy
      Write(6,"(a,2(1x,i2))") "Rootx,Rooty   : ",
     +                         Iroot,Jroot 
      Call Checksum("TDENS-READ ",Work(Iscr3),Nao2,s)
                      Call Genrlizd_beyond_dipole(Work(Iscr3),
     +                                            Work(Iscr1),
     +                                            Work(Iscr2),
     +                                            Work(Iqm6),
     +                                            Work(Iom10),
     +                                            Work(Iam3),
     +                                            Work(Idm3),
     +                                            Work(Imq9),
     +                                            Work(Iend),
     +                                            Maxcor,Nao,
     +                                            Jroot,Irrepy,Iroot,
     +                                            Irrepx) 
                   Enddo
               Enddo
            Enddo 
         Enddo
      Enddo
      Write(6,*) 
      Write(6,"(a)") "Printing from print_bd"
      Write(6,"(a)") "DM-coontribution"
      call print_bd1(Work(Idm3),3,4,Nirrep,Nroot,3) 
      Write(6,*) 
      Write(6,"(a)") "QM-coontribution"
      call print_bd1(Work(Iqm6),6,4,Nirrep,Nroot,6) 
      Write(6,*) 
      Write(6,"(a)") "OM-coontribution"
      call print_bd1(Work(Iom10),10,4,Nirrep,Nroot,10) 
      Write(6,*) 
      Write(6,"(a)") "AM-coontribution"
      call print_bd1(Work(Iam3),3,4,Nirrep,Nroot,3) 
      Write(6,*) 
      Write(6,"(a)") "MD-coontribution"
      call print_bd1(Work(Imq9),9,4,Nirrep,Nroot,9) 
      Write(6,*) 

      Length = 0 
      Do Irrepx_1 = 1, Nirrep
         Do Iroot_1 = 1, Nroot(Irrepx_1)
            E_i = Eigval(Iroot_1,Irrepx_1)
            Do Irrepy_1 = 1, Nirrep
               Do Jroot_1 = 1, Nroot(Irrepy_1)
                  E_j = Eigval(jroot_1,Irrepy_1)

                  If (E_i .eq. E_j) Then
                      E = Dnull 
                  Else
                      E = Dabs(E_i - E_j)
                  Endif 

                  Do Irrepx_2 = 1, Nirrep
                     Do Iroot_2 = 1, Nroot(Irrepx_2)
                        Do Irrepy_2 = 1, Nirrep
                           Do Jroot_2 = 1, Nroot(Irrepy_2)
                           
                           Call Osc_bd_init(Work(Iqm6),
     +                                      Work(Iom10),Work(Iam3),
     +                                      Work(Idm3),Work(Imq9),
     +                                      Work(Iqm6_2),
     +                                      Work(Iom10_2),  
     +                                      Work(Iam3_2),
     +                                      Work(Idm3_2),  
     +                                      Work(Imq9_2),
     +                                      Work(Iend),Maxcor,
     +                                      Irrepx_1,Iroot_1,Irrepy_1,
     +                                      jroot_1,
     +                                      Irrepx_2,Iroot_2,Irrepy_2,
     +                                      jroot_2)
     +                                      
                           Enddo
                        Enddo
                     Enddo
                  Enddo
                 Length = Length + 1
                 Call Osc_bd_fin(Work(Iosc_qm),Work(Iosc_om),
     +                           Work(Iosc_am),Work(Iosc_dm),
     +                           Work(Iosc_mq),Work(Iqm6_2),
     +                           Work(Iom10_2),Work(Iam3_2),
     +                           Work(Idm3_2), Work(Imq9_2), 
     +                           Jroot_1,Irrepy_1,Iroot_1,
     +                           Irrepx_1,E,Length)
               Enddo
            Enddo
         Enddo
      Enddo

      Call Putrec(20,"JOBARC","QUADMOMC",Length*Iintfp,Work(IOsc_Qm))
      Call Putrec(20,"JOBARC","OCTPMOMC",Length*Iintfp,Work(IOsc_Om))
      Call Putrec(20,"JOBARC","ANGLMOMC",Length*Iintfp,Work(IOsc_Am))
      Call Putrec(20,"JOBARC","MAGQUADC",Length*Iintfp,Work(IOsc_Mq))
      Call Putrec(20,"JOBARC","MAGQUADC",Length*Iintfp,Work(IOsc_Dm))

      Return
      End
