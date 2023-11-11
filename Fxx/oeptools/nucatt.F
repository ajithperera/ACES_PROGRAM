      subroutine nucatt (
     &     natm, nprm, maxshl, maxang, 
     &     nshlatm, nprmshl, angfct, alpha, 
     &     crdnuc, crdatm, attint)
c
c The subroutines evaluates nuclear attraction integrals in 
c primitive functions representation (<mu|1/(r-R)|nu>). It is
c Ajith's TWODENS/repprdint () routine modified for the INTGRT
c module. The repprdint () uses subroutines borrowed from the
c VMOL(?) property integral package.
c
c Igor Schweigert, Oct 2002 
c $Id: nucatt.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c 
      implicit none
c
c     .. input parameters ..
c
      integer
     &     natm, nprm, maxshl, maxang, i,j,
     &     nshlatm (natm), nprmshl (maxshl, natm), angfct (maxang,3)
c
      double precision
     &     crdnuc (3), crdatm (3, natm),
     &     alpha (nprm), attint (nprm, nprm)
c
c     .. local variables ..
c
      integer
     &     iprm, iatm, ishl, jprm, jatm, jshl, k, 
     &     ip (20), ioffst (15)
c
      double precision
     &     crdatmi (3), crdatmj (3), fac (9,9)
c      
c     .. initialize the factorial array ..
c
      call setrhf(fac, ioffst, ip)
c
c     .. clear the nuclear attraction integrals array ..
c
      call dzero(attint, nprm*nprm)
c
c     .. the first loop over atoms and angular momentum shells  ..
c
      iprm = 1
      do iatm = 1, natm
         do ishl = 1, nshlatm (iatm)
c
c           .. the second loop over atoms and angular momentum shells  ..
c
            jprm = 1
            do jatm = 1, natm
               do jshl = 1, nshlatm (jatm)
c                  
c                 .. coordinates of atoms ..
c
                  do k = 1, 3
                     crdatmi(k)  = crdatm(k, iatm)
                     crdatmj(k)  = crdatm(k, jatm)
                  enddo
c
c                 .. evaluate the nuc.att.int. for all the primitives in the shell ..
c                  
                  call evlnucatt (
     &                 iprm, nprmshl (ishl, iatm), 
     &                 jprm, nprmshl (jshl, jatm), nprm, 
     &                 angfct(ishl,1), angfct(ishl,2), angfct(ishl,3),
     &                 angfct(jshl,1), angfct(jshl,2), angfct(jshl,3),  
     &                 crdatmi, crdatmj, crdnuc,
     &                 alpha, fac, attint)
c
                  jprm = jprm + nprmshl (jshl, jatm)
               enddo
            enddo
c           
            iprm = iprm + nprmshl (ishl, iatm)
         enddo
      enddo
c
      return
      end
