include ../configure.user

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mpi.o     \
              $(COMMON_DIR)/common_mtx.o     \
              $(COMMON_DIR)/netlib.o         \
              $(COMMON_DIR)/common_letkf.o   \
              $(COMMON_DIR)/common_gfs.o

PROGS   = gfsmean

OBJS    = sscal_tools.o

INCS    = -I$(SIGIO_INC) \
          -I$(SFCIO_INC)

LIBS    = -L$(BACIO_LIB) -lbacio_4 \
          -L$(SIGIO_LIB) -lsigio_4 \
          -L$(SFCIO_LIB) -lsfcio_4 \
          -L$(SP_LIB) -lsp_4 \
          $(LIBBLAS)

all:	$(PROGS)

gfsmean:	gfsmean.o $(OBJS) common_modules
	$(FC) $(FOPTS) -o gfsmean$(shell ../common/enssize) gfsmean.o $(COMMON_OBJS) $(OBJS) $(LIBS)

gfsmean.o:	gfsmean.f90 $(OBJS) common_modules
	$(FC) $(FOPTS) -c gfsmean.f90 -I$(COMMON_DIR) $(INCS)

sscal_tools.o:	common_modules
	$(FC) $(FOPTS) -c sscal_tools.f90 -I$(COMMON_DIR) $(INCS)

common_modules:
	( cd $(COMMON_DIR) ; make )

clean:
	rm -f *.o *.mod $(PROGS) gfsmean[0-9]*
