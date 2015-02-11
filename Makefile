CC        = g++
#CCFLAGS   = -Wall -g
CCFLAGS   = -Wall -O2
SOURCES   =
ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs --ldflags`           
ROOFITLIBS = -lRooFit -lRooFitCore -lMinuit -lFoam -I${ROOFITSYS}include -L${ROOFITSYS}lib
ROOSTATSLIBS = -lRooStats -I${ROOFITSYS}include/RooStats/
# boost
#BOOSTFLAGS = -I${BOOST_ROOT}include/boost-1_48
#BOOSTLIBS = -L${BOOST_ROOT}lib -lboost_program_options-gcc43-mt-1_48

TMVA = -L${ROOTSYS}lib -lTMVA

all: R2GGBBBiasStudy_2D.exe

R2GGBBBiasStudy_2D.exe: R2GGBBBiasStudy_2D.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) $(ROOSTATSLIBS) $(ROOFITLIBS) R2GGBBBiasStudy_2D.cc -o R2GGBBBiasStudy_2D.exe

clean:
	rm *.exe
