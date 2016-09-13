include ../configure.user

PROGS   = scale-rm_pp_ens scale-rm_init_ens scale-rm_ens

OBJS    =

COMMON_DIR = ../.common_for_scale

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mpi.o     \
              $(COMMON_DIR)/common_ncio.o    \
              $(COMMON_DIR)/common_scale.o     \
              $(COMMON_DIR)/common_obs_scale.o \
              $(COMMON_DIR)/common_mpi_scale.o \
              $(COMMON_DIR)/common_nml.o     

INCS    = $(SCALE_INC) -I$(SCALE_RM_OBJDIR)

LIBS    = $(SCALE_LIB) $(NETCDF_LIB) $(LAPACK_LIB)

SCALE_RM_OBJS = $(shell ls $(SCALE_RM_OBJDIR)/*.o | grep -v "scale-rm.o\|scale-rm_pp.o\|scale-rm_init.o")

ifeq ($(H08),T)
COMMON_OBJS += $(COMMON_DIR)/scale_H08_fwd.o
LIBS    += $(RT_INCLUDE) $(RT_LDFLAGS)
endif

all:	$(PROGS)

scale-rm_pp_ens:	scale-rm_pp_ens.o $(OBJS) common_modules
	$(FC) $(FOPTS_SCALE) $(FMCMODEL) -o scale-rm_pp_ens scale-rm_pp_ens.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS)

scale-rm_pp_ens.o:	scale-rm_pp_ens.f90 $(OBJS) common_modules
	$(FC) $(FOPTS_SCALE) -c scale-rm_pp_ens.f90 -I$(COMMON_DIR) $(INCS)

scale-rm_init_ens:	scale-rm_init_ens.o $(OBJS) common_modules
	$(FC) $(FOPTS_SCALE) $(FMCMODEL) -o scale-rm_init_ens scale-rm_init_ens.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS)

scale-rm_init_ens.o:	scale-rm_init_ens.f90 $(OBJS) common_modules
	$(FC) $(FOPTS_SCALE) -c scale-rm_init_ens.f90 -I$(COMMON_DIR) $(INCS)

scale-rm_ens:	scale-rm_ens.o $(OBJS) common_modules
	$(FC) $(FOPTS_SCALE) $(FMCMODEL) -o scale-rm_ens scale-rm_ens.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS)

scale-rm_ens.o:	scale-rm_ens.f90 $(OBJS) common_modules
	$(FC) $(FOPTS_SCALE) -c scale-rm_ens.f90 -I$(COMMON_DIR) $(INCS)

common_modules:
	( cd $(COMMON_DIR) ; make )

clean:
	rm -f *.o *.mod *.lst $(PROGS)
