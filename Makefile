# 1. Put this file in the same folder as your 'driver' code 
#    (the code containing the 'main' function).

# 2. Edit LIBRARY_DIR to point at the location of your ITensor Library
#    source folder (this is the folder that has options.mk in it).
#    Also, edit IDMRG_DIR to point at the location of the iDMRG source files.
LIBRARY_DIR=$(HOME)/itensor

# 3. If your 'main' function is in a file called 'myappname.cc', then
#    set APP to 'myappname'. Running 'make' will compile the app.
#    Running 'make debug' will make a program called 'myappname-g'
#    which includes debugging symbols and can be used in gdb (Gnu debugger);
APP=quench

# 4. Add any headers your program depends on here. The make program
#    will auto-detect if these headers have changed and recompile your app.
MYDIR=$(HOME)/itensor.utility/

MYFLAGS=-fmax-errors=3 -Wno-unused-variable -Wno-unused-function -Wno-sign-compare

HEADERS=MyObserver.h MixedBasis.h SortBasis.h SpecialFermion.h tdvp.h TDVPObserver.h basisextension.h InitState.h BdGBasis.h OneParticleBasis.h Hamiltonian.h

# 5. For any additional .cc (source) files making up your project,
#    add their full filenames here.
CCFILES=$(APP).cc

#################################################################
#################################################################
#################################################################
#################################################################


include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h


#Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) $(MYFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) $(MYFLAGS) -o $@ $<

#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP).exe $(LIBFLAGS)

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g.exe $(LIBGFLAGS)

clean:
	rm -fr .debug_objs *.o $(APP).exe $(APP)-g.exe

mkdebugdir:
	mkdir -p .debug_objs

