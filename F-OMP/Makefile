# Set the target system
# Can use one of this list: archer, cirrus or linux

SYSTEM=archer

# System dependent definitions

# For ARCHER - default

ifeq ($(SYSTEM), archer)
FC=	ftn
FFLAGS=	-O3 -h omp
LFLAGS=
endif

# For Cirrus using Intel compilers

ifeq ($(SYSTEM), cirrus)
FC=	ifort
FFLAGS=	-qopenmp -O3
LFLAGS=
endif

# For a standard Linux cluster using GNU compilers

ifeq ($(SYSTEM), linux)
FC=	gfortran
FFLAGS=	-fopenmp -O3
LFLAGS=
endif

# System independent definitions

MF=	Makefile

EXE=	cfd

SRC= \
	boundary.f90 \
	cfd.f90 \
	cfdio.f90 \
	jacobi.f90

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

cfd.o:	boundary.o jacobi.o cfdio.o

tar:
	tar cvf cfd.tar $(MF) $(INC) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) *.mod velocity.dat colourmap.dat cfd.plt core
