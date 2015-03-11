#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# ThermWire
#
# Elena Diaz, Francisco Dominguez-Adame, Rudolf A Roemer
#
# (C) 2015, all rights reserved
#
# Version: :VERSION:
# Date:    :DATE:
# Time:    :TIME:
# Status:  :RLSTATUS:
# Build:   :BUILD:
# Author:  :AUTHOR:
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#  This file is part of felixsim.
#
#  felixsim is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  felixsim is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# $Id: makefile-CSC.GF,v 1.2 2014/03/27 18:13:30 phsht Exp $
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F90 = gfortran
FC = gfortran -ffixed-form -std=legacy
F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none # debug
#F90FLAGS = -pg -fbounds-check -fno-range-check -fimplicit-none # profiling
F77FLAGS = -fbounds-check -fno-range-check # debug

INCLUDEDIR = /warwick/mathlib/gnu/x86_64/ # CSC machines
#INCLUDEDIR = /usr/lib/ # RAR Linux

#LIBFLAGS= -llapack -L../cephes/ -lmd #-lblas #-lfftw3
LIBFLAGS= -llapack

OBJECTFILES = main.f90 bessel.f90

%.o: %.f90
#	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)
	$(F90) -c $< $(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.f
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

all:	thermwire.GF

thermwire.GF: $(OBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(OBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(OBJECTFILES) $(LIBFLAGS)

print:	
	a2ps -E -o print.ps main.f90; convert -density 150 print.ps print.pdf
	     	   	  
clean:	
	rm -f core *.mod *.o

#clean all:	
#	rm -f core *.mod *.o *.exe
