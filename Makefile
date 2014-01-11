#writtn by Xue Bing 
#2014-01-10

#fortran complier
FC     = ifort
FFLAGS = -llapack

DBGDIR = build/debug/
RLSDIR = build/release/
SRCDIR = src/
DBGDATDIR = data/debug/
RLSDATDIR = data/release/
EXEC    = hubbard
FILES   = time.f90 system_data.f90 basis.f90 main.f90
OBJECTS = $(FILES:.f90=.o)
vpath %.f90 $(SRCDIR)

.PHONY : dbg rls rundbg runrls cldbg clrls clean cndbg cnrls ctags cscope tags doc dir dbgdir rlsdir

#debug
DBGTARGET = $(DBGDIR)$(EXEC)
DBGOBJS = $(addprefix $(DBGDIR), $(OBJECTS))
DBGMODS = $(wildcard $(DBGDIR)*.mod)
dbg : DFFLAGS = -O0 -debug -g -traceback -warn $(FFLAGS)
dbg : $(DBGTARGET)
$(DBGTARGET) : $(DBGOBJS)
	$(FC) $(DFFLAGS) -o $@ $^
$(DBGOBJS) : $(DBGDIR)%.o : %.f90
	$(FC) $(DFFLAGS) -c $< -o $@ -module $(DBGDIR)

#release
RLSTARGET = $(RLSDIR)$(EXEC)
RLSOBJS = $(addprefix $(RLSDIR), $(OBJECTS))
RLSMODS = $(wildcard $(RLSDIR)*.mod)
rls : RFFLAGS = -O3 -ip -ipo -xHost -static -mtune=pentium4m $(FFLAGS)
rls : $(RLSTARGET)
$(RLSTARGET) : $(RLSOBJS)
	$(FC) $(RFFLAGS) -o $@ $^
$(RLSOBJS) : $(RLSDIR)%.o : %.f90
	$(FC) $(RFFLAGS) -c $< -o $@ -module $(RLSDIR)


#make debug directory
dbgdir :
	if [ ! -d "$(DBGDIR)" ]; then mkdir -p $(DBGDIR); fi

#make release directory
rlsdir :
	if [ ! -d "$(RLSDIR)" ]; then mkdir -p $(RLSDIR); fi

#clean debug files
cldbg :
	rm -f $(DBGOBJS) $(DBGMODS) $(DBGTARGET)

#clean release files
clrls :
	rm -f $(RLSOBJS) $(RLSMODS) $(RLSTARGET)

#run debug
rundbg : $(DBGTARGET)
	if [ ! -d "$(DBGDATDIR)" ]; then mkdir -p $(DBGDATDIR); fi
	cd $(DBGDATDIR) && ../../$<

#run release
runrls : $(RLSTARGET)
	if [ ! -d "$(RLSDATDIR)" ]; then mkdir -p $(RLSDATDIR); fi
	cd $(RLSDATDIR) && ../../$<

#make directories
dir: dbgdir rlsdir

#clean && debug
cndbg : cldbg dbg

#clean && release
cnrls : clrls rls

#clean
clean : cldbg clrls

#make tags
ctags : $(FILES)
	cd $(SRCDIR) && rm -f tags
	cd $(SRCDIR) && ctags -R --extra=+q --fields=+Saim --c++-kinds=+lpx --c-kinds=+lpx --fortran-kinds=+iL *.f90
cscope : $(FILES)
	cd $(SRCDIR) && rm -f cscope.*
	cd $(SRCDIR) && cscope -Rbq *.f90
tags : ctags cscope

#documents by doxygen
doc :
	doxygen Doxyfile
