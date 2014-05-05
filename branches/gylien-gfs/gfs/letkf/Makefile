include ../configure.user

PROGS   = letkf mean obsope efso

OBJS    = letkf_obs.o letkf_tools.o efso_tools.o efso_nml.o

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mpi.o     \
              $(COMMON_DIR)/common_mtx.o     \
              $(COMMON_DIR)/netlib.o         \
              $(COMMON_DIR)/common_letkf.o   \
              $(COMMON_DIR)/common_gfs.o     \
              $(COMMON_DIR)/common_obs_gfs.o \
              $(COMMON_DIR)/common_mpi_gfs.o \
              $(COMMON_DIR)/common_precip.o

INCS    = -I$(SIGIO_INC)

LIBS    = -L$(SIGIO_LIB) -lsigio_4 \
          -L$(SP_LIB) -lsp_4 \
          $(LIBBLAS)

all:	$(PROGS)

letkf:	letkf.o $(OBJS) common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o letkf$(shell ../common/enssize) letkf.o \
	$(COMMON_OBJS) $(OBJS) $(LIBS)
mean:	mean.o letkf_obs.o common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o mean$(shell ../common/enssize) mean.o \
	$(COMMON_OBJS) $(OBJS) $(LIBS)
obsope:	obsope.o common_modules
	$(SFC) $(FOPTS) -o obsope obsope.o \
	$(COMMON_DIR)/SFMT.o $(COMMON_DIR)/common.o \
	$(COMMON_DIR)/common_gfs.o $(COMMON_DIR)/common_obs_gfs.o \
	-L$(SP_LIB) -lsp_4
efso:	efso.o $(OBJS) common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o efso$(shell ../common/enssize) efso.o \
	$(COMMON_OBJS) $(OBJS) $(LIBS)

letkf.o:	letkf.f90 $(OBJS) common_modules
	$(FC) $(FOPTS) -c letkf.f90 -I$(COMMON_DIR) $(INCS)
mean.o:	mean.f90 letkf_obs.o common_modules
	$(FC) $(FOPTS) -c mean.f90 -I$(COMMON_DIR) $(INCS)
obsope.o:	obsope.f90 common_modules
	$(SFC) $(FOPTS) -c obsope.f90 -I$(COMMON_DIR) $(INCS)
efso.o:	efso.f90 $(OBJS) common_modules
	$(FC) $(FOPTS) -c efso.f90 -I$(COMMON_DIR) $(INCS)

letkf_obs.o:	letkf_obs.f90 common_modules
	$(FC) $(FOPTS) -c letkf_obs.f90 -I$(COMMON_DIR) $(INCS)
letkf_tools.o:	letkf_tools.f90 letkf_obs.o efso_nml.o efso_tools.o common_modules
	$(FC) $(FOPTS) -c letkf_tools.f90 -I$(COMMON_DIR) $(INCS)
efso_tools.o:	efso_tools.f90 letkf_obs.o efso_nml.o common_modules
	$(FC) $(FOPTS) -c efso_tools.f90 -I$(COMMON_DIR) $(INCS)
efso_nml.o:	efso_nml.f90 common_modules
	$(SFC) $(FOPTS) -c efso_nml.f90 -I$(COMMON_DIR)

common_modules:
	( cd $(COMMON_DIR) ; make )

clean:
	rm -f *.o *.mod $(PROGS) letkf[0-9]* mean[0-9]* efso[0-9]*
