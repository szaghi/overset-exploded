#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = exe/obj/
DMOD    = exe/mod/
DEXE    = exe/
LIBS    =
FC      = gfortran
OPTSC   = -c -std=f2008 -fall-intrinsics -O2 -J exe/mod
OPTSL   = -O2 -J exe/mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)OVERSET-EXPLODED: $(MKDIRS) $(DOBJ)overset-exploded.o
	@rm -f $(filter-out $(DOBJ)overset-exploded.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) OVERSET-EXPLODED

#compiling rules
$(DOBJ)overset-exploded.o: src/app/overset-exploded.F90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
