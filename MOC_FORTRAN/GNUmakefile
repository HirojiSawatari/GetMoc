#
# File:  GNUmakefile
#
#  The bourne shell is safest for gmake
SHELL = /bin/sh

TOP   = .

SRCDIRS = $(TOP)

TARGETS = pop_moc pop_moc_ei

#  VPATH is the built-in symbol whose value is the path that gmake will 
#  search for dependencies.
VPATH = $(SRCDIRS)

########################################################################
#                  Macro Definitions
########################################################################
#
#on TCS
#INCLDIR = /usr/local/include
#on paleo
INCLDIR = /usr/include

#SUN INCLDIR = /usr/local/netcdf-3.4/include
#
FC=     ifort
#SUN FFLAGS= -O3 -I$(INCLDIR)
#FFLAGS= -O3 -r4 -mips4 -64 -I$(INCLDIR)
#on TCS
FFLAGS = -O2 -q64 -I$(INCLDIR)


#
LDR =   ifort
LDFLAGS= 
#
#SUN LIBDIR = -L/usr/local/lib
#LIBDIR = -L/usr/local/lib64/r4i4

#on TCS
#LIBDIR = -L/usr/local/lib

#on paleo
LIBDIR = -L/usr/lib64

#
LIBS = $(LIBDIR) -lnetcdff
#
########################################################################
#              Non-standard Suffix Rules
########################################################################
#
.SUFFIXES : .f90 .o
.f90.o :
	$(FC) -c $(FFLAGS) $*.f90
#
########################################################################
#            Main Program Build Rules
########################################################################
#
all: $(TARGETS)

pop_moc : pop_moc.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $@.o $(LIBS)
pop_moc_ei : pop_moc_ei.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $@.o $(LIBS)
#
########################################################################
#
###########################################################################
clean:
	rm -f *.o
