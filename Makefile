# C++ flags 
#-----------------------------------
CXX = g++

# set library flags for compiler
#-----------------------------------
FJFLAGS	:= $(shell fastjet-config --cxxflags)
FJLIBS	:= $(shell fastjet-config --libs )
ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTGLIBS	:= $(shell root-config --glibs)


INCLUDES = -I$(DELPHES) -I$(DELPHES)/external

# set platform specific flags
#-----------------------------------
PLATFORM_TYPE := $(shell uname -s)
ifeq ($(PLATFORM_TYPE), Linux)
	SOFLAGS := -shared
else
	ifeq ($(PLATFORM_TYPE), Darwin)
		SOFLAGS := -dynamiclib -undefined dynamic_lookup
	endif
endif


CXXFLAGS = -ansi -pedantic -Wall -O3 -Wno-long-long $(FJFLAGS) $(ROOTCFLAGS) $(INCLUDES) -std=c++11

# rules
# ----------------------------------------------------------------------------

.SUFFIXES:      .o .cxx .f .exe .C


# instructions for building a .o file from a .cxx file
# ----------------------------------------------------------------------------

#FILES = Qantikt.o QantiktPlugin.o

#all: Delphes_Analysis

#lib/libQantikt.a: $(FILES) $(FILES:.cc=.o)
#	ar cq lib/libQantikt.a $(FILES)

Delphes_Analysis: Delphes_Analysis.o
	@echo "Building $@	 ..."
	$(CXX) $(CXXFLAGS) Delphes_Analysis.o \
	        $(FJLIBS) \
		-L$(DELPHES) -lDelphes \
		-L./lib \
	        $(ROOTGLIBS) -o $@

.o: %.C %.h
	g++ -fPIC  -O3 -c $(ROOTCFLAGS) $(FJFLAGS) $< -o $@ 

clean:
	rm  Delphes_Analysis *.o



