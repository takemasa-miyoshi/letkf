include ../configure.user

PROGS   = letkf #mean efso

OBJS    = letkf_obs.o letkf_tools.o #efso_tools.o efso_nml.o

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mpi.o     \
              $(COMMON_DIR)/common_mtx.o     \
              $(COMMON_DIR)/netlib.o         \
              $(COMMON_DIR)/common_letkf.o   \
              $(COMMON_DIR)/common_ncio.o    \
              $(COMMON_DIR)/common_scale.o     \
              $(COMMON_DIR)/common_obs_scale.o \
              $(COMMON_DIR)/common_mpi_scale.o \
              $(COMMON_DIR)/common_nml.o
#              $(SCALE_RM_OBJDIR)/mod_admin_time.o \
#              $(SCALE_RM_OBJDIR)/mod_atmos_vars.o
#              $(COMMON_DIR)/common_scalelib.o
#              $(COMMON_DIR)/common_precip.o

OBSOPE_DIR = ../obs

OBSOPE_OBJS = $(OBSOPE_DIR)/obsope_tools.o

INCS    = $(SCALE_INC)

LIBS    = $(SCALE_LIB) $(NETCDF_LIB) $(LAPACK_LIB)

ifeq ($(H08),T)
COMMON_OBJS += $(COMMON_DIR)/scale_H08_fwd.o
LIBS    += $(RT_INCLUDE) $(RT_LDFLAGS)
endif

all:	$(PROGS)

letkf:	letkf.o $(OBJS) common_modules obsope_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o letkf letkf.o \
	$(COMMON_OBJS) $(OBSOPE_OBJS) $(OBJS) $(LIBS)
mean:	mean.o letkf_obs.o common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o mean mean.o \
	$(COMMON_OBJS) $(OBJS) $(LIBS)
efso:	efso.o $(OBJS) common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o efso efso.o \
	$(COMMON_OBJS) $(OBJS) $(LIBS)

letkf.o:	letkf.f90 $(OBJS) common_modules obsope_modules
	$(FC) $(FOPTS) -c letkf.f90 -I$(COMMON_DIR) -I$(OBSOPE_DIR) $(INCS)
mean.o:	mean.f90 letkf_obs.o common_modules
	$(FC) $(FOPTS) -c mean.f90 -I$(COMMON_DIR) $(INCS)
efso.o:	efso.f90 $(OBJS) common_modules
	$(FC) $(FOPTS) -c efso.f90 -I$(COMMON_DIR) $(INCS)

letkf_obs.o:	letkf_obs.f90 common_modules
	$(FC) $(FOPTS) -c letkf_obs.f90 -I$(COMMON_DIR) $(INCS)
letkf_tools.o:  letkf_tools.f90 letkf_obs.o common_modules
	$(FC) $(FOPTS) -c letkf_tools.f90 -I$(COMMON_DIR) $(INCS)
#letkf_tools.o:	letkf_tools.f90 letkf_obs.o efso_nml.o efso_tools.o common_modules
#	$(FC) $(FOPTS) -c letkf_tools.f90 -I$(COMMON_DIR) $(INCS)
efso_tools.o:	efso_tools.f90 letkf_obs.o efso_nml.o common_modules
	$(FC) $(FOPTS) -c efso_tools.f90 -I$(COMMON_DIR) $(INCS)
efso_nml.o:	efso_nml.f90 common_modules
	$(SFC) $(FOPTS) -c efso_nml.f90 -I$(COMMON_DIR)

common_modules:
	( cd $(COMMON_DIR) ; make )

obsope_modules:
	( cd $(OBSOPE_DIR) ; make obsope_tools.o )

clean:
	rm -f *.o *.mod *.lst $(PROGS) letkf mean efso
