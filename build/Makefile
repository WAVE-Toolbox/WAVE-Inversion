################################################################################
# WAVE CODE BUILD SYSTEM
#
# v1.0a, Jan 2018, Thomas Hertweck: Initial version.
# v1.0b, Jan 2018, Tilman Steinweg: GTEST_DIR PATH update.
################################################################################

SRCTREE  := ../
EXENAME  := IFOS
LIBNAME  := IFOS
UTSTNAME := Test_unit
ITSTNAME := Test_CompareMisfit

GPI_ROOT_DEFAULT  := /data14/software/SOFI++/build_dev
SCAI_ROOT_DEFAULT := /data14/software/scai_lama/build_dev
GTEST_DIR_DEFAULT := /data14/software/googletest/googletest

################################################################################
# Usually, there is no need to change anything below this line
################################################################################

.SUFFIXES:
.SECONDARY:

MAKEFLAGS += -rR
SHELL     := /bin/bash

wave.dir  := $(realpath $(SRCTREE))
wave.exe  := $(strip $(EXENAME))
wave.lib  := $(strip $(LIBNAME))
wave.utst := $(strip $(UTSTNAME))
wave.itst := $(strip $(ITSTNAME))
gtest     := $(strip $(shell gtester --version 2>/dev/null))

ifeq ($(realpath $(PWD)),$(wave.dir)/src)
  $(error You are trying to build software in the source tree; please use a separate build directory)
endif

ifdef V
  ifeq ("$(origin V)", "command line")
    BUILD_VERBOSE := $(V)
  endif
endif
ifndef BUILD_VERBOSE
  BUILD_VERBOSE := 0
endif
ifeq ($(BUILD_VERBOSE),1)
  quiet    := 
  Q        := 
else
  quiet    := quiet_
  Q        := @
endif
ifneq ($(filter 4.%,$(MAKE_VERSION)),)  # make-4
ifneq ($(filter %s ,$(firstword x$(MAKEFLAGS))),)
  quiet := silent_
endif
else                                    # make-3.8x
ifneq ($(filter s% -s%,$(MAKEFLAGS)),)
  quiet := silent_
endif
endif

ifndef GPI_ROOT
  GPI_ROOT := $(GPI_ROOT_DEFAULT)
  $(warning GPI_ROOT environment variable undefined, now set to $(GPI_ROOT))
endif

ifndef SCAI_ROOT
  SCAI_ROOT := $(SCAI_ROOT_DEFAULT)
  $(warning SCAI_ROOT environment variable undefined, now set to $(SCAI_ROOT))
endif

ifndef GTEST_DIR
  ifeq ($(gtest),)   
    GTEST_DIR := $(GTEST_DIR_DEFAULT)
    $(warning no local installation found, GTEST_DIR environment variable undefined, now set to $(GTEST_DIR))
  endif
endif
wave.libdir := ./lib
wave.bindir := ./bin
wave.objdir := ./obj
#wave.docdir := ./doc
dummy := $(shell mkdir -p $(wave.libdir) $(wave.bindir) $(wave.objdir))
#dummy := $(shell mkdir -p $(wave.libdir) $(wave.bindir) $(wave.objdir) $(wave.docdir))

wave.deps := $(wildcard .dep*)

wave.srcdir := $(addsuffix /src,$(wave.dir))
wave.commit := $(strip $(shell cd $(wave.dir) && git log --pretty=short --format="%H" 2>/dev/null | head -n1))
ifeq ($(wave.commit),)
  wave.commit := -1
endif

wave.srcsubdirs := $(dir $(wildcard $(wave.srcdir)/*/.))
wave.liballsrc  := $(wildcard $(addsuffix */*.cpp,$(wave.srcsubdirs))) $(wildcard $(addsuffix *.cpp,$(wave.srcsubdirs)))
wave.unitsrc    := $(wildcard $(addsuffix */*UnitTest.cpp,$(wave.srcsubdirs))) $(wildcard $(addsuffix *UnitTest.cpp,$(wave.srcsubdirs)))
wave.integsrc   :=  $(wildcard $(addsuffix *Test_CompareMisfit.cpp,$(wave.srcsubdirs)))
wave.libsrc     := $(filter-out $(wave.unitsrc) $(wave.integsrc),$(wave.liballsrc))
wave.mainsrc    := $(wave.srcdir)/IFOS.cpp

wave.libobj   := $(addprefix $(wave.objdir)/,$(notdir $(patsubst %.cpp,%.o,$(wave.libsrc))))
wave.mainobj  := $(addprefix $(wave.objdir)/,$(notdir $(patsubst %.cpp,%.o,$(wave.mainsrc))))
wave.unitobj  := $(addprefix $(wave.objdir)/,$(notdir $(patsubst %.cpp,%.o,$(wave.unitsrc))))
wave.integobj := $(addprefix $(wave.objdir)/,$(notdir $(patsubst %.cpp,%.o,$(wave.integsrc))))

CXX      := g++
CPPFLAGS := -DSCAI_ASSERT_LEVEL_DEBUG -DSCAI_LOG_LEVEL_DEBUG -DSCAI_TRACE_OFF
CXXFLAGS :=  -std=c++11 -fopenmp -Wall -Wextra -Werror
INCDIRS   = -I${GPI_ROOT}/include -I$(SCAI_ROOT)/include $(addprefix -I,$(wave.srcsubdirs))
CDEPEND   = -MMD -MP -MF .dep-$(subst /,-,$@)
LD       := $(CXX)
LDFLAGS  := -L$(wave.libdir)  -l$(wave.lib) -L${GPI_ROOT}/lib -lSOFI \
            -L$(SCAI_ROOT)/lib -lscai_lama -lscai_solver -lscai_dmemo -lscai_sparsekernel \
            -lscai_utilskernel -lscai_blaskernel -lscai_kregistry -lscai_hmemo -lscai_tasking \
            -lscai_tracing -lscai_logging -lscai_common 
AR       := ar
ARFLAGS  := cr
RM       := rm -f --preserve-root

.PHONY: help all clean distclean vinfo utest itest lib prog #doc doxygen theory

help:
	@echo
	@echo "Usage:"
	@echo "   $(MAKE) [V=0|1] <target>"
	@echo
	@echo "Targets:"
	@echo "   all       :: equivalent to \"$(MAKE) lib prog doc utest itest\""
	@echo "   lib       :: build archive lib$(wave.lib).a"
	@echo "   prog      :: build executable $(wave.exe) (and, implicitly, lib)"
	#@echo "   doc       :: build html (make doxygen) and pdf (make theory) documentation"
	@echo "   utest     :: build unit tests"
	@echo "   itest     :: build integration tests"
	@echo "   clean     :: remove object files"
	@echo "   distclean :: clean up entire build directory"
	@echo "   vinfo     :: print version information stored in executable"
	@echo "                (strings <executable> | grep _KITGPI_LINK_)"
	@echo
	@echo "By default, the build is quiet (V=0) and we do not display the full"
	@echo "command being executed, only a brief summary. If you would like to"
	@echo "see the whole command, use \"$(MAKE) V=1 <target>\"."
	@echo
	@echo "Top-level source tree (git) directory set to: $(wave.dir)"
ifeq "$(origin SCAI_ROOT)" 'environment'
	@echo "SCAI_ROOT environment variable set to: $(SCAI_ROOT)"
endif
ifeq "$(origin GTEST_DIR)" 'environment'
	@echo "GTEST_DIR environment variable set to: $(GTEST_DIR)"
endif
	@echo

all: lib prog utest itest #doc 

prog: $(wave.bindir)/$(wave.exe)

lib: $(wave.libdir)/lib$(wave.lib).a

utest: $(wave.bindir)/$(wave.utst)

itest: $(wave.bindir)/$(wave.itst)

$(wave.bindir)/$(wave.exe): $(wave.mainobj) $(wave.libdir)/lib$(wave.lib).a
	$(call cmd,link)

$(wave.libdir)/lib$(wave.lib).a: $(wave.libobj)
	$(call cmd,ar,$(wave.objdir))

$(wave.unitobj): INCDIRS += -I$(GTEST_DIR)/include
$(wave.unitobj): CXXFLAGS += -isystem -pthread
$(wave.bindir)/$(wave.utst): LDFLAGS += -L$(GTEST_DIR)/ -lgtest -lgtest_main
$(wave.bindir)/$(wave.utst): $(wave.unitobj) $(wave.libdir)/lib$(wave.lib).a
	$(call cmd,link)

$(wave.bindir)/$(wave.itst): $(wave.integobj) $(wave.libdir)/lib$(wave.lib).a
	$(call cmd,link)

clean:
	$(call cmd,rmdir,$(wave.objdir))

distclean: clean
	@echo -e "   \e[31m[RM]\e[0m    .dep-* files"
	-$(Q)find . -iname ".dep-*" | xargs $(RM)
	-$(Q)$(RM) *~
	$(call cmd,rmdir,$(wave.libdir))
	$(call cmd,rmdir,$(wave.bindir))
	$(call cmd,rmdir,$(wave.docdir))

vinfo: prog
	$(Q)strings $(wave.bindir)/$(wave.exe) | grep _KITGPI_LINK_

#doc: doxygen theory
#doxygen: $(wave.docdir)/doxygen/html/index.html
#theory: $(wave.docdir)/theory/WAVE_guide_EN.pdf

# $(wave.docdir)/doxygen/html/index.html: $(wave.dir)/doc/doxygen/Doxyfile $(wave.liballsrc)
# 	@echo "   [DOC]   $@"
# 	$(Q)mkdir -p $(wave.docdir)/doxygen && cd  $(wave.docdir)/doxygen && \
# 	sed -e 's#\(^INPUT[[:space:]]*=\).*$$#\1 $(wave.srcdir) $(wave.dir)/README.md#g' \
#             -e 's#\(^USE_MDFILE_AS_MAINPAGE[[:space:]]*=\).*$$#\1 $(wave.dir)/README.md#g' $< > Doxyfile && \
# 	doxygen Doxyfile &> log 
# 	@echo "   [DOC]   doxygen log file: $(wave.docdir)/doxygen/log"
# 
# $(wave.docdir)/theory/WAVE_guide_EN.pdf: $(wave.dir)/doc/theory/WAVE_guide_EN.tex
# 	@echo "   [DOC]   $@"
# 	$(Q)mkdir -p $(wave.docdir)/theory && cd  $(wave.docdir)/theory && \
# 	sed -e '/\\author{}/ a\\\graphicspath{{$(wave.dir)/doc/theory/}}' $< \
# 	    -e 's#\\bibliography{WAVE_theory}#\\bibliography{$(wave.dir)/doc/theory/WAVE_theory}#g' > WAVE_guide_EN.tex && \
# 	pdflatex WAVE_guide_EN > /dev/null && bibtex WAVE_guide_EN > /dev/null && \
# 	pdflatex WAVE_guide_EN > /dev/null && pdflatex WAVE_guide_EN > log
# 	@echo "   [DOC]   pdflatex log file: $(wave.docdir)/theory/log"

$(wave.objdir)/%.o: $(wave.srcdir)/%.cpp
	$(call cmd,compile_cxx)

$(wave.objdir)/%.o: $(wave.srcdir)/*/%.cpp
	$(call cmd,compile_cxx)

$(wave.objdir)/%.o: $(wave.srcdir)/*/*/%.cpp
	$(call cmd,compile_cxx)

quiet_cmd_compile_cxx = \e[32m[CXX]\e[0m   $@
      cmd_compile_cxx = $(CXX) -o$@ $(CDEPEND) $(CXXFLAGS.$@) $(CPPFLAGS) $(INCDIRS) $(CXXFLAGS) -c $<

quiet_cmd_ar = \e[93m[AR]\e[0m    $@
      cmd_ar = cd $(1) && $(AR) $(ARFLAGS) $(abspath $@) $(notdir $^) && ranlib $(abspath $@)

quiet_cmd_rmdir = \e[31m[RM]\e[0m    $(1)
      cmd_rmdir = $(RM) -rf $(1)

lnkbase = $(wave.objdir)/wave_link_$(notdir $@)
lnkh    = $(addsuffix .hpp,$(lnkbase))
lnkcxx  = $(addsuffix .cpp,$(lnkbase))
lnko    = $(addsuffix .o,$(lnkbase))

os.string      := $(shell awk -F\" '/PRETTY_NAME/ {print $$2}' /etc/os-release)
cpp.string     := $(subst \,,$(CPPFLAGS))
ldflags.string  = $(subst \,,$(LDFLAGS))

quiet_cmd_link = \e[94m[LD]\e[0m    $@
      cmd_link = $(RM) $(lnkh) $(lnkcxx) $(lnko) ; \
		 echo "\#ifndef _KITGPI_LINK_H" > $(lnkh) ; \
		 echo "\#define _KITGPI_LINK_H" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_HOST[] = \"_KITGPI_LINK_HOST $(HOSTNAME)\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_OS[] = \"_KITGPI_LINK_OS $(os.string)\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_DATE[] = \"_KITGPI_LINK_DATE `date`\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_USER[] = \"_KITGPI_LINK_USER `whoami`\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_PWD[] = \"_KITGPI_LINK_PWD `pwd`\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_IMAGE[] = \"_KITGPI_LINK_IMAGE $(notdir $@)\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_CPPFLAGS[] = \"_KITGPI_LINK_CPPFLAGS $(cpp.string)\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_INCDIRS[] = \"_KITGPI_LINK_INCDIRS $(INCDIRS)\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_CXXFLAGS[] = \"_KITGPI_LINK_CXXFLAGS $(CXXFLAGS)\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_LD[] = \"_KITGPI_LINK_LD $(LD)\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_LDFLAGS[] = \"_KITGPI_LINK_LDFLAGS $(ldflags.string)\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_SRCTREE[] = \"_KITGPI_LINK_SRCTREE $(wave.dir)\";" >> $(lnkh) ; \
		 echo "char KITGPI_DAT_GIT[] = \"_KITGPI_LINK_GITCOMMIT $(wave.commit)\";" >> $(lnkh) ; \
		 echo "void kitgpi_link(void);" >> $(lnkh) ; \
		 echo "\#endif" >> $(lnkh) ; \
		 echo "\#include \"$(notdir $(lnkh))\"" >> $(lnkcxx) ; \
		 echo "\#include <cstdio>" >> $(lnkcxx) ; \
		 echo "using namespace std;" >> $(lnkcxx) ; \
		 echo "void kitgpi_link(void) { printf(\"%s\",KITGPI_DAT_HOST); return; }" >> $(lnkcxx) ; \
		 $(CXX) -I$(wave.objdir) $(CXXFLAGS) -c -o$(lnko) $(lnkcxx) && \
		 $(LD) $(CXXFLAGS) -o$@ $(filter-out $(abspath $(wave.libdir)/lib$(wave.lib).a),$(abspath $^)) \
		       $(lnko) $(LDFLAGS.$@) $(LDFLAGS) ; ret=$$? ; \
		 $(RM) $(lnkh) $(lnkcxx) $(lnko) ; exit $$ret

cmd = @$(if $($(quiet)cmd_$(1)),echo -e '   $(call $(quiet)cmd_$(1),$(2))' &&) $(call cmd_$(1),$(2))

ifneq ($(wave.deps),)
  -include $(wave.deps)
endif
