include make.inc
#################################################
### Defining Compiling variables              ###
#################################################

CXX           = g++ -Wno-deprecated -Wall
#-m32  -Wno-deprecated
LD            = g++ -Wno-deprecated -Wall
#-m32  -Wno-deprecated 
INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################

SHAREDLIBFLAGS =-L./CommonUtils/lib -lTauNtuple 
SHAREDCXXFLAGS =-I./CommonUtils 
SUBDIRS = CommonUtils/ 
TAUSPINERFLAG = 

override CXXFLAGS += 
ALLCXXFLAGS = $(CXXFLAGS) -I$(ROOTSYS)/include -c -I/include $(SHAREDCXXFLAGS) -c -I./ -c
LDFLAGS       =  
ROOTLIBS      =  -L$(ROOTSYS)/lib -L/usr/lib/ -L/lib/i686/  -lCore -lCint -lHist -lGraf  -lGraf3d -lGpad -lTree -lRint -lReflexDict -lReflex -lPostscript -lMatrix -lPhysics -lGui -LObj -lMinuit2 -lThread -rdynamic -Wl,--rpath $(ROOTSYS)/lib

# -ldl -lpthread

DEFS  	      = -DSTANDALONE=1
LIBS          = $(ROOTLIBS) $(SHAREDLIBFLAGS)

TARGETS       =	Ntuple_Controller \
		NtupleReader \
		Parameters \
		Selection_Base \
		Selection_Factory \
		Selection \
		Tables \
		Plots \
		DoubleEventRemoval \
		PDG_Var \
		Tools \
		SkimConfig \
		HistoConfig \
		Example \
		TauSolver \
		TauSpinExample \
		TauSpinerInterface \
		ReferenceScaleFactors \
		rochcor2012jan22 \
		Objects \
		UncertaintyValue

CINTTARGETS = 

# now include conditional options set by flags in make.inc
ifdef USE_TauSpinner
	SUBDIRS += TauSpiner/
	SHAREDLIBFLAGS +=	-L./TauSpiner/lhapdf/lib -lLHAPDF \
				-L./TauSpiner/tauola++/1.1.3/lib -lTauolaCxxInterface -lTauolaFortran -lTauolaTauSpinner \
				-L./TauSpiner/hepmc/install/lib -lHepMC
	SHAREDCXXFLAGS +=	-I./TauSpiner/hepmc/install/include \
				-I./TauSpiner/tauola++/1.1.3/include/Tauola/ \
				-I./TauSpiner/tauola++/1.1.3/include/TauSpinner/ \
				-I./TauSpiner/tauola++/1.1.3/include/ \
				-I./TauSpiner/lhapdf/include/ \
				-I./TauSpiner/pythia8/176/include/
	DEFS += -DUSE_TauSpinner=1
endif

ifdef USE_cherepanov
	TARGETS += Ztotautau_hadmu_ControlSample Tau_momentum_calculation Validation
	DEFS += -DUSE_cherepanov=1
endif

ifdef USE_inugent
	TARGETS += TauLifeTime TriggerStudy TriggerStudyMC ZDouble3prong Ztotautau_ControlSample Ztomumu_ControlSample ChargedHiggs_dilepontic ChargedHiggs_tauplusjet 
	DEFS += -DUSE_inugent=1
endif

ifdef USE_kargoll
	TARGETS += HToTaumuTauh Category
	TARGETS += OneJetBoost OneJetHigh OneJetLow VBFLoose VBFTight ZeroJetHigh ZeroJetLow Inclusive
	TARGETS += CategoryStandalone VBFLooseStandalone VBFTightStandalone RelaxedVBFTightStandalone RelaxedVBFLooseStandalone
	#TARGETS += MuTauSync MCDecayChain
	TARGETS += ZeroJet3Prong
	DEFS += -DUSE_kargoll=1
endif 

ifdef USE_nehrkorn
ifndef USE_TauSpinner
	TAUSPINERFLAG = LHAPDF
	SUBDIRS += TauSpiner/
	SHAREDLIBFLAGS +=	-L./TauSpiner/lhapdf/lib -lLHAPDF
	SHAREDCXXFLAGS +=	-I./TauSpiner/lhapdf/include/
endif
	TARGETS += PDFweights ZtoEMu Tvariable_Base Tvariable_EE Tvariable_MuMu Tvariable_EMu
	DEFS += -DUSE_nehrkorn=1
endif

ifdef USE_pistone
	DEFS += -DUSE_pistone=1
endif
ifdef USE_zotz
	TARGETS += ZToTaumuTauh
	DEFS += -DUSE_zotz=1
endif
ifdef USE_stahl
	TARGETS += 
	DEFS += -DUSE_stahl=1
endif


HDRS          = $(TARGETS:=.h) 		
SRCS          = $(TARGETS:=.cxx)  Analysis.cxx
OBJS          = $(SRCS:.cxx=.o)

PROGRAM       = Analysis.exe

$(PROGRAM):     $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@echo $(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@$(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@echo "done"

VPATH = utilities:i386_linux
vpath %.cxx inugent
vpath %.h inugent
vpath %.cxx cherepanov
vpath %.h cherepanov
vpath %.cxx kargoll
vpath %.h kargoll
vpath %.cxx nehrkorn
vpath %.h nehrkorn
vpath %.cxx zotz
vpath %.h zotz
vpath %.cxx stahl
vpath %.h stahl
vpath %.cxx CommonFiles
vpath %.h CommonFiles

#	@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o ../$(PROGRAM)

##########################################################
###                                                    ###
##########################################################

$(OBJS): %.o : %.cxx
	$(CXX) $(ALLCXXFLAGS) $(DEFS) $< -o i386_linux/$@ 

.PHONY: clean cleanall all install sharedlib 

install: Analysis.exe


clean:
	@rm i386_linux/*.o
	@rm Analysis.exe

cleanall:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) -f Makefile.standalone clean && cd ../../; \
	done
	@rm i386_linux/*.o
	@rm Analysis.exe

all: sharedlib install

sharedlib:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		if [ $$subdir == "TauSpiner/" ]; then \
			echo with flag; \
			cd $$subdir && $(MAKE) -f Makefile.standalone $(TAUSPINERFLAG) && cd ../ && pwd; \
		else \
			echo without flag ; \
			cd $$subdir echo $(CXXFLAGS) && $(MAKE) -f Makefile.standalone CXXFLAGS+='$(CXXFLAGS)' && cd ../ && pwd; \
		fi \
	done
	$(install)

