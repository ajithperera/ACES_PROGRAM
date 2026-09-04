
      program intgrt

c This program defines a grid of points and determines the value of
c the MO, the electron density, and several energy densities (such as
c the exchange density, the kinetic energy density, etc.) at each point
c (both at the SCF and correlated levels).  It also integrates over
c the gridpoints to determine total values.

c Written by Nevin Oliphant as hfdft, Quantum Theory Project,
c University of Florida, 1993.

c ReWritten by Nevin Oliphant and Sullivan Beck as potdens, 
c Quantum Theory Project, University of Florida, 1994.

c ReWritten by Nevin Oliphant and Sullivan Beck as intgrt, 
c Quantum Theory Project, University of Florida, 1995.

      implicit none
#include <aces.h>
#include <flags.com>
#include <flags2.com>
#include <iuhf.com>
#include <mol.com>
#include <kscore.com>
#include <symm2_ks.com>
#include <timeinfo.com>
#include <int.com>
#include <intgrtflags.com>

      integer intnumradpts,max_angpts,islateron
      integer maxmem, i0tmp,Length 
      logical kshf
      call timer(0)
      call aces_init_rte
      call aces_ja_init

      call intgrt_run

      call aces_ja_fin
      call timer(1)
      write(*,9020) timetot
 9020 format(t3,'@CHECKOUT-I, Total execution time : ',f12.4,
     &       ' seconds.')
      stop
      end


