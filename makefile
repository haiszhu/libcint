# makefile for cgto... need to see if there are missing source files
# 
# 07/27/24 Hai
CFLAGS = -I./src -I./include -fPIC -O2 -Wall -D_Float128=__float128 -DWITH_CINT2_INTERFACE -fopenmp
LDFLAGS = -shared -lm -fopenmp -lquadmath
LIBS = -lblas -lm -ldl
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
CC = gcc-14
MEX = /Applications/MATLAB_R2023b.app/bin/mex
MWRAP = ~/mwrap/mwrap
LBLAS = -L/opt/homebrew/opt/openblas/lib -lopenblas
CINTDLIBRARY_NAME = libcint.dylib
NPHELPERDLIBRARY_NAME = libnp_helper.dylib
CGTODLIBRARY_NAME = libcgto.dylib
endif
ifeq ($(UNAME), Linux)
CC = gcc
MEX = mex
MWRAP = ~/mwrap/mwrap
LBLAS = -lopenblas
CINTDLIBRARY_NAME = libcint.so
NPHELPERDLIBRARY_NAME = libnp_helper.so
CGTODLIBRARY_NAME = libcgto.so
endif

# Paths
ROOT_DIR = $(CURDIR)
CINT_SOURCE_DIR = $(ROOT_DIR)/src
CINT_AUTOCODE_SOURCE_DIR = $(ROOT_DIR)/src/autocode
NP_HELPER_SOURCE_DIR = $(ROOT_DIR)/np_helper
CGTO_SOURCE_DIR = $(ROOT_DIR)/gto
CGTO_AUTOCODE_SOURCE_DIR = $(ROOT_DIR)/gto/autocode
LIB_DIR = $(ROOT_DIR)
BINARY_DIR = $(ROOT_DIR)/bin

# Object files
LIBCINT_OBJS = $(BINARY_DIR)/c2f.o $(BINARY_DIR)/cart2sph.o \
               $(BINARY_DIR)/cint1e.o $(BINARY_DIR)/cint2e.o \
               $(BINARY_DIR)/cint_bas.o $(BINARY_DIR)/fblas.o \
               $(BINARY_DIR)/g1e.o $(BINARY_DIR)/g2e.o \
               $(BINARY_DIR)/misc.o $(BINARY_DIR)/optimizer.o \
               $(BINARY_DIR)/fmt.o $(BINARY_DIR)/rys_wheeler.o \
               $(BINARY_DIR)/eigh.o $(BINARY_DIR)/rys_roots.o \
               $(BINARY_DIR)/find_roots.o $(BINARY_DIR)/cint2c2e.o \
               $(BINARY_DIR)/g2c2e.o $(BINARY_DIR)/cint3c2e.o \
               $(BINARY_DIR)/g3c2e.o $(BINARY_DIR)/cint3c1e.o \
               $(BINARY_DIR)/g3c1e.o $(BINARY_DIR)/breit.o \
               $(BINARY_DIR)/cint1e_a.o $(BINARY_DIR)/cint3c1e_a.o \
               $(BINARY_DIR)/cint1e_grids.o \
               $(BINARY_DIR)/g1e_grids.o $(BINARY_DIR)/breit1.o \
               $(BINARY_DIR)/dkb.o $(BINARY_DIR)/gaunt1.o \
               $(BINARY_DIR)/grad1.o $(BINARY_DIR)/grad2.o \
               $(BINARY_DIR)/hess.o $(BINARY_DIR)/int3c1e.o \
               $(BINARY_DIR)/int3c2e.o $(BINARY_DIR)/intor1.o \
               $(BINARY_DIR)/intor2.o $(BINARY_DIR)/intor3.o \
               $(BINARY_DIR)/intor4.o $(BINARY_DIR)/deriv3.o \
               $(BINARY_DIR)/int1e_grids1.o $(BINARY_DIR)/deriv4.o \
               $(BINARY_DIR)/lresc.o 
#							 $(BINARY_DIR)/g2e_f12.o \
#               $(BINARY_DIR)/stg_roots.o $(BINARY_DIR)/cint2e_f12.o

NP_HELPER_OBJS = $(BINARY_DIR)/transpose.o $(BINARY_DIR)/pack_tril.o \
                 $(BINARY_DIR)/npdot.o $(BINARY_DIR)/condense.o \
                 $(BINARY_DIR)/omp_reduce.o $(BINARY_DIR)/np_helper.o

CGTO_OBJS = $(BINARY_DIR)/fill_int2c.o $(BINARY_DIR)/fill_nr_3c.o \
            $(BINARY_DIR)/fill_r_3c.o $(BINARY_DIR)/fill_int2e.o \
            $(BINARY_DIR)/fill_r_4c.o $(BINARY_DIR)/ft_ao.o \
            $(BINARY_DIR)/ft_ao_deriv.o $(BINARY_DIR)/fill_grids_int2c.o \
            $(BINARY_DIR)/grid_ao_drv.o $(BINARY_DIR)/deriv1.o \
            $(BINARY_DIR)/deriv2.o $(BINARY_DIR)/nr_ecp.o \
            $(BINARY_DIR)/nr_ecp_deriv.o $(BINARY_DIR)/auto_eval1.o

MY_OBJS = $(BINARY_DIR)/cgto.o

# Targets
.PHONY: all libcint np_helper cgto clean

all: libcint np_helper cgto libcint.a np_helper.a cgto.a mexfile

$(BINARY_DIR)/%.o: $(CINT_SOURCE_DIR)/%.c 
	$(CC) $(CFLAGS) -I$(ROOT_DIR)/include -I$(ROOT_DIR)/src -c $< -o $@

$(BINARY_DIR)/%.o: $(CINT_AUTOCODE_SOURCE_DIR)/%.c
	$(CC) $(CFLAGS) -I$(ROOT_DIR)/include -I$(ROOT_DIR)/src -c $< -o $@

$(BINARY_DIR)/%.o: $(NP_HELPER_SOURCE_DIR)/%.c
	$(CC) $(CFLAGS) -I$(NP_HELPER_SOURCE_DIR) -I$(ROOT_DIR) -c $< -o $@

$(BINARY_DIR)/%.o: $(CGTO_SOURCE_DIR)/%.c
	$(CC) $(CFLAGS) -I$(CGTO_SOURCE_DIR) -I$(ROOT_DIR) -c $< -o $@

$(BINARY_DIR)/%.o: $(CGTO_AUTOCODE_SOURCE_DIR)/%.c
	$(CC) $(CFLAGS) -I$(CGTO_AUTOCODE_SOURCE_DIR) -I$(ROOT_DIR) -c $< -o $@

$(BINARY_DIR)/cgto.o: $(ROOT_DIR)/cgto.c
	$(CC) $(CFLAGS) -I$(ROOT_DIR)/include -I$(ROOT_DIR) -c $< -o $@

libcint: $(LIBCINT_OBJS)
	$(CC) $(LDFLAGS) -o $(LIB_DIR)/$(CINTDLIBRARY_NAME) $(LIBCINT_OBJS) $(LIBS)

np_helper: $(NP_HELPER_OBJS)
	$(CC) $(LDFLAGS) -o $(LIB_DIR)/$(NPHELPERDLIBRARY_NAME) $(NP_HELPER_OBJS) $(LIBS)

cgto: libcint np_helper $(CGTO_OBJS)
	$(CC) $(LDFLAGS) -o $(LIB_DIR)/$(CGTODLIBRARY_NAME) $(CGTO_OBJS) -L$(LIB_DIR) -lcint -lnp_helper $(LIBS)

libcint.a: $(LIBCINT_OBJS)
	ar rcs $(LIB_DIR)/libcint.a $(LIBCINT_OBJS)

np_helper.a: $(NP_HELPER_OBJS)
	ar rcs $(LIB_DIR)/libnp_helper.a $(NP_HELPER_OBJS)

cgto.a: libcint.a libnp_helper.a $(CGTO_OBJS)
	ar rcs $(LIB_DIR)/libcgto.a $(CGTO_OBJS)

# Generate the MEX interface using MWrap
gateway.c: gateway.mw
	$(MWRAP) -c99complex -mex gateway -mb -list gateway.mw
	$(MWRAP) -c99complex -mex gateway -c gateway.c gateway.mw
# not sure all the flags are needed (seems to work on linux)... need to check with -v and other things... 07/27/24 Hai
mexfile: gateway.c $(BINARY_DIR)/cgto.o
	$(MEX) gateway.c $(BINARY_DIR)/cgto.o -largeArrayDims -lm -lgomp -lquadmath -lstdc++ -llapack $(LBLAS) -L$(LIB_DIR) libcgto.a libcint.a -I$(ROOT_DIR)
#	mex -v gateway.c $(BINARY_DIR)/cgto.o -largeArrayDims -lm -lgomp -lquadmath -lstdc++ -llapack -lopenblas -L$(LIB_DIR) libcgto.a libcint.a -I$(ROOT_DIR)

clean:
	rm -f $(BINARY_DIR)/*.o \
	      $(ROOT_DIR)/$(CINTDLIBRARY_NAME) $(ROOT_DIR)/$(NPHELPERDLIBRARY_NAME) $(ROOT_DIR)/$(CGTODLIBRARY_NAME) \
				$(ROOT_DIR)/libcint.a $(ROOT_DIR)/libnp_helper.a $(ROOT_DIR)/libcgto.a \
				$(ROOT_DIR)/*.mex*
