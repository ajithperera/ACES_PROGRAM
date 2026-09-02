# GNUmakefile.chssi (included below) sets its own site-default FC/CC/CPP and
# FFLAGS/CFLAGS/CPPFLAGS. Every sibling under special_directory/ then clobbers
# them with configure's substituted values (FC=@FC@ etc. in its Makefile.in) --
# but this directory has no Makefile.in and is not in configure.ac's
# AC_CONFIG_FILES, so what is checked in is a FROZEN copy of one old configure
# run's output instead. Capture the toolchain the top-level Makefile exported
# (blanket `export`, so configure's own values arrive here as environment)
# BEFORE the include can overwrite it, and restore it after -- which is what
# the generated siblings effectively do. Left as-is, this directory compiled
# with chssi's own mpif77, which under Intel MPI is gfortran, and died on
# -Zp8/-traceback/-shared-intel while every sibling built fine.
ENV_FC := $(FC)
ENV_CC := $(CC)
ENV_CPP := $(CPP)
ENV_SERIAL_CPP := $(SERIAL_CPP)
ENV_FFLAGS := $(FFLAGS)
ENV_CFLAGS := $(CFLAGS)
ENV_CPPFLAGS := $(CPPFLAGS)

include ../../../../../GNUmakefile.chssi

# Fall back to the old frozen literals only when built standalone, i.e. with
# nothing exported from a parent make. The -march=core-avx2 in them replaces
# the -msse3 -axcore-avx2,core-avx-i this file used to carry: those are the
# Intel-CPU-only dispatch flags 07938684 removed everywhere else, and they
# make the resulting binary refuse to run on AMD.
#Compilers
FC := $(if $(ENV_FC),$(ENV_FC),mpif77)
CC := $(if $(ENV_CC),$(ENV_CC),mpicc)
CPP := $(if $(ENV_CPP),$(ENV_CPP),mpicxx)
SERIAL_CPP := $(if $(ENV_SERIAL_CPP),$(ENV_SERIAL_CPP),icpc)

#Compiler flags
FFLAGS := $(if $(ENV_FFLAGS),$(ENV_FFLAGS),-D__fortran -D__fortran77 -DMPIF2C -Zp8 -zero -traceback -DMPI2 -shared-intel -O2 -march=core-avx2)
CFLAGS := $(if $(ENV_CFLAGS),$(ENV_CFLAGS),-DMPIF2C -DC_SUFFIX  -DMPI2 -DCB_SUFFIX -shared-intel -O2 -march=core-avx2)
CPPFLAGS := $(if $(ENV_CPPFLAGS),$(ENV_CPPFLAGS),-DMPIF2C -DC_SUFFIX -DCB_SUFFIX  -DMPI2 -shared-intel -O2 -march=core-avx2)

#ACESLIBS variable will be replaced by the autoconf script with appropriate libraries
LIBS = -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -Wl,--end-group
LIBS:= -lstdc++ -lsip1 -lsip2 -lsip_shared -lframelib -laces2 -lgeopt -lsymcor -laces2 -lerd -loed -lprop_ints -ldup -lsip1 -lsip2 $(LIBS)

#ACESFLAGS variable will be replaced by the autoconf script with location of libraries
LIB_DIRS= -L.

SIAL_COMPILER_LIBS= -lsial -lsip_shared -laces2  -lifcore


INCLUDE_DIRS := -I../../../../sia/include -I../aces2/include $(INCLUDE_DIRS)
LIB_DIRS:=-L../../../lib $(LIB_DIRS)

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

PROP_INTSLIB=../../../../../lib/libprop_ints.a

all: tabula_rasa $(PROP_INTSLIB)

tabula_rasa: ;
	rm -f make.out make.err

$(PROP_INTSLIB): libprop_intslib.a
	cp libprop_intslib.a $(PROP_INTSLIB)
	rm -f make.out make.err

libprop_intslib.a: $(F_OBJS) $(f_OBJS) $(c_OBJS)
	ar $(ARFLAGS) libprop_intslib.a $(F_OBJS) $(f_OBJS) $(c_OBJS)

libclean: ;
	rm -f libprop_intsdlib.a

ppclean clean distclean: % : tabula_rasa libclean
	rm -f *.o

rebuild: libclean $(FRAMELIB)

binclean relink archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

