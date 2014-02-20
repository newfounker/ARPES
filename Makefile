#written by Xue Bing
#2014-01-10

#fortran complier
#FC     = gfortran
FC     = ifort

LAPACK = -llapack
#LAPACK = -L/opt/OpenBlas/lib/libopenblas.so
FFLAGS = -static -fpe3 -warn $(LAPACK)

ifeq ($(FC), gfortran)
	DBG_FLAGS = -pipe -O0 -static -fno-exceptions -g -ggdb -fbacktrace -Wall -Wextra -mtune=k8 $(LAPACK) -J$(DBG_DIR)
	RLS_FLAGS = -pipe -O3 -static -fno-exceptions -Wall -Wextra -mtune=k8 $(LAPACK) -J$(RLS_DIR)
endif

ifeq ($(FC), ifort)
	DBG_FLAGS = -O0 -debug -g -traceback $(FFLAGS) -module $(DBG_DIR)
	RLS_FLAGS = -O3 -ip -ipo -xHost -mtune=pentium4m $(FFLAGS) -module $(RLS_DIR)
endif

SRC_DIR = src/
BUILD_DIR = build/
DBG_DIR = $(BUILD_DIR)debug/
RLS_DIR = $(BUILD_DIR)release/
DAT_DIR = data/
DBG_DAT_DIR = $(DAT_DIR)debug/
RLS_DAT_DIR = $(DAT_DIR)release/
EXEC    = hubbard
FILES   = time.f90 system_data.f90 basis.f90 hamiltonian.f90 main.f90
OBJECTS = $(FILES:.f90=.o)
vpath %.f90 $(SRC_DIR)

.PHONY : dbg rls rundbg runrls cldbg clrls clean cndbg cnrls ctags cscope tags doc nuke

default : dbg
	
#debug
DBG_TARGET = $(DBG_DIR)$(EXEC)
DBG_OBJS = $(addprefix $(DBG_DIR), $(OBJECTS))
DBG_MODS = $(wildcard $(DBG_DIR)*.mod)
dbg : DFFLAGS = $(DBG_FLAGS)
dbg : $(DBG_TARGET)
$(DBG_TARGET) : $(DBG_OBJS)
	$(FC) $(DFFLAGS) -o $@ $^
$(DBG_OBJS) : $(DBG_DIR)%.o : %.f90
	@if [ ! -d "$(DBG_DIR)" ]; then mkdir -p $(DBG_DIR); fi
	$(FC) $(DFFLAGS) -c $< -o $@

#release
RLS_TARGET = $(RLS_DIR)$(EXEC)
RLS_OBJS = $(addprefix $(RLS_DIR), $(OBJECTS))
RLS_MODS = $(wildcard $(RLS_DIR)*.mod)
rls : RFFLAGS = $(RLS_FLAGS)
rls : $(RLS_TARGET)
$(RLS_TARGET) : $(RLS_OBJS)
	$(FC) $(RFFLAGS) -o $@ $^
$(RLS_OBJS) : $(RLS_DIR)%.o : %.f90
	@if [ ! -d "$(RLS_DIR)" ]; then mkdir -p $(RLS_DIR); fi
	$(FC) $(RFFLAGS) -c $< -o $@

#clean debug files
cldbg :
	rm -f $(DBG_OBJS) $(DBG_MODS) $(DBG_TARGET)

#clean release files
clrls :
	rm -f $(RLS_OBJS) $(RLS_MODS) $(RLS_TARGET)

#run debug
rundbg : $(DBG_TARGET)
	@if [ ! -d "$(DBG_DAT_DIR)" ]; then mkdir -p $(DBG_DAT_DIR); fi
	@cd $(DBG_DAT_DIR) && ../../$<

#run release
runrls : $(RLS_TARGET)
	@if [ ! -d "$(RLS_DAT_DIR)" ]; then mkdir -p $(RLS_DAT_DIR); fi
	@cd $(RLS_DAT_DIR) && ../../$<

#clean && debug
cndbg : cldbg dbg

#clean && release
cnrls : clrls rls

#clean
clean : cldbg clrls

#clean all
nuke :
	rm -rf build/
	rm -rf data/
	rm -rf doc/

#make tags
ctags : $(FILES)
	cd $(SRC_DIR) && rm -f tags
	cd $(SRC_DIR) && ctags -R --extra=+q --fields=+Saim --c++-kinds=+lpx --c-kinds=+lpx --fortran-kinds=+iL *.f90
cscope : $(FILES)
	cd $(SRC_DIR) && rm -f cscope.*
	cd $(SRC_DIR) && cscope -Rbq *.f90
tags : ctags cscope

#documents by doxygen
doc :
	doxygen Doxyfile
