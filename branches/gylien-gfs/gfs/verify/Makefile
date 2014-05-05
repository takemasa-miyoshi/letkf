include ../configure.user

PROGS   = verify

OBJS    =

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o            \
              $(COMMON_DIR)/common.o          \
              $(COMMON_DIR)/common_gfs.o      \
              $(COMMON_DIR)/common_gfs_pres.o \
              $(COMMON_DIR)/common_obs_gfs.o 

all:	$(PROGS)

verify:	verify.o $(OBJS) common_modules
	$(SFC) $(FOPTS) $(FMCMODEL) -o verify verify.o $(COMMON_OBJS) $(OBJS) -L$(SP_LIB) -lsp_4

verify.o:	verify.f90 $(OBJS) common_modules
	$(SFC) $(FOPTS) -c verify.f90 -I$(COMMON_DIR)

common_modules:
	( cd $(COMMON_DIR) ; make )

clean:
	rm -f *.o *.mod $(PROGS)
