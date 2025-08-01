
# ldd ./matlab/bdmk_module.mexa64 

EXEC = int2-bdmk
UNAME := $(shell uname)
ARCH := $(shell uname -m)

# brew install hdf5
ifeq ($(UNAME), Darwin)
	HOST=osx
	ifeq ($(ARCH), arm64)
		OPENBLAS_DIR = /opt/homebrew/Cellar/openblas-singlethread/0.3.29/lib
		SRC_DIR = /Users/hzhu/Documents/Github/job/bdmkf
		HDF5_INCLUDE = -I/opt/homebrew/include
		HDF5_DIR = /opt/homebrew/lib
	else
		OPENBLAS_DIR = /Users/hzhu1/git/OpenBLAS
		SRC_DIR = /Users/hzhu1/Documents/Github/job/bdmkf
		HDF5_INCLUDE = -I/usr/local/include
		HDF5_DIR = /usr/local/lib
	endif
endif
ifeq ($(UNAME), Linux)
HOST=linux-gfortran
OPENBLAS_DIR = /home/linuxbrew/.linuxbrew/Cellar/openblas-singlethread/0.3.29/lib
SRC_DIR = /home/hai/Documents/Github/job/bdmkf
HDF5_INCLUDE = -I/usr/include/hdf5/serial
HDF5_DIR = /usr/lib/x86_64-linux-gnu/hdf5/serial
endif
# HOST=linux-ifort

ifeq ($(HOST),osx)
FC = gfortran
#FFLAGS = -O3 -march=native -std=legacy -funroll-loops -c -w
FFLAGS = -fPIC -O3 -march=native -std=legacy -funroll-loops -c -w
FFLAGS_OMP = $(FFLAGS) -fopenmp
FLINK = gfortran -w -o $(EXEC)
FEND = $(OPENBLAS_DIR)/libopenblas.a -lm
LBLAS = -L$(OPENBLAS_DIR)/ -lopenblas
MWRAP = ~/mwrap/mwrap
MEX = /Applications/MATLAB_R2023b.app/bin/mex
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
# FFLAGS = -fPIC -O3 -march=native -std=legacy -fopenmp -funroll-loops -c -w
FFLAGS = -fPIC -O3 -march=native -std=legacy -funroll-loops -c -w
FFLAGS_OMP = $(FFLAGS) -fopenmp
# FLINK = gfortran -fopenmp -w -o $(EXEC)
FLINK = gfortran -w -o $(EXEC)
FEND = $(OPENBLAS_DIR)/libopenblas.a -lm
LBLAS = -L$(OPENBLAS_DIR)/ -lopenblas
MWRAP = ~/mwrap/mwrap
MEX = mex
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
#FFLAGS = -fPIC -O3 -march=native -qopenmp -funroll-loops -c -w 
FFLAGS = -fPIC -O3 -march=x86-64 -qopenmp -funroll-loops -c -w
#FFLAGS = -O3 -fp-model=strict -qopenmp -xHost -c -w
FLINK = ifort -qopenmp -o $(EXEC)
#FEND = -qmkl=sequential
FEND = -qmkl
endif

# SRC = ./src

.PHONY: all clean list bdmk_wrap.o

SOURCES =  $(SRC_DIR)/bdmk_wrap.f90 \
  $(SRC_DIR)/bdmk.f \
  $(SRC_DIR)/boxfgt_md.f \
  $(SRC_DIR)/src/besseljs3d.f \
  $(SRC_DIR)/src/hank103.f \
  $(SRC_DIR)/src/legeexps.f \
  $(SRC_DIR)/src/chebexps.f \
  $(SRC_DIR)/src/prini_new.f \
  $(SRC_DIR)/src/fmmcommon2d.f \
  $(SRC_DIR)/src/lapack_f77.f \
  $(SRC_DIR)/src/cumsum.f \
  $(SRC_DIR)/src/hkrand.f \
  $(SRC_DIR)/src/dlaran.f \
  $(SRC_DIR)/src/voltab2d.f \
  $(SRC_DIR)/src/voltab3d.f \
  $(SRC_DIR)/src/polytens.f \
  $(SRC_DIR)/src/tree_data_routs.f \
  $(SRC_DIR)/src/tensor_prod_routs.f \
  $(SRC_DIR)/src/pts_tree.f \
  $(SRC_DIR)/src/tree_routs.f \
  $(SRC_DIR)/src/tree_vol_coeffs.f \
  $(SRC_DIR)/src/dmk_routs.f \
  $(SRC_DIR)/src/get_sognodes.f \
  $(SRC_DIR)/src/l2dsognodes.f \
  $(SRC_DIR)/src/l3dsognodes.f \
  $(SRC_DIR)/src/sl3dsognodes.f \
  $(SRC_DIR)/src/y2dsognodes.f \
  $(SRC_DIR)/src/y3dsognodes.f \
  $(SRC_DIR)/src/bdmk_local_tables.f \
  $(SRC_DIR)/src/bdmk_local.f \
  $(SRC_DIR)/src/bdmk_pwterms.f \
  $(SRC_DIR)/src/bdmk_pwrouts.f

ifeq ($(WITH_SECOND),1)
SOURCES += $(SRC_DIR)/src/second-r8.f
endif

# OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))
OBJECTS = $(patsubst $(SRC_DIR)/%.f, %.o, $(patsubst $(SRC_DIR)/src/%.f, %.o, \
          $(patsubst $(SRC_DIR)/%.f90, %.o, $(SOURCES))))
# $(info OBJECTS contains: $(OBJECTS))
#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) $< -o $@

%.o : $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $< -o $@

%.o : $(SRC_DIR)/src/%.f
	$(FC) $(FFLAGS) -c $< -o $@

%.mod : $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $< 

all: bdmk_wrap.o helper.o $(OBJECTS) int2-bdmk-mlscf mex matlab

bdmk_wrap.o: $(SRC_DIR)/bdmk_wrap.f90
	$(FC) $(FFLAGS_OMP) -c $< -o $@

helper.o: helper.f90
	$(FC) $(FFLAGS_OMP) -c $< -o $@

int2-bdmk-mlscf: $(OBJECTS) bdmk_mlscf.f90
	$(FC) -O3 -o int2-bdmk-mlscf -fopenmp bdmk_mlscf.f90 $(OBJECTS) \
				-L. $(OPENBLAS_DIR)/libopenblas.a -lm \
				$(HDF5_INCLUDE) -L$(HDF5_DIR) -lhdf5_fortran -lhdf5

mex: helper_module.mw
	$(MWRAP) -c99complex -mex helper_module -mb -list helper_module.mw
	$(MWRAP) -c99complex -mex helper_module -c helper_module.c helper_module.mw
	mv ./*_mex.m ../

matlab: helper_module.c helper.o
	$(MEX) helper_module.c helper.o -largeArrayDims $(MFLAGS) -lgfortran -lgomp -lm -lstdc++ $(LBLAS) -llapack -outdir ../

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)
	rm -f fort*
	rm -f int*
#	rm -f $(MATLAB_DIR)/*.m $(MATLAB_DIR)/*.mex* 

list: $(SOURCES)
	$(warning Requires:  $^)


