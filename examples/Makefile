
ROOT_INCLUDE=$(shell root-config --cflags)
ROOT_LIBS=$(shell root-config --libs)

PMTCalib_include = -I../src
PMTCalib_lib     = -L../lib -lPMTCalib

INCLUDES += $(ROOT_INCLUDE)
INCLUDES += $(PMTCalib_include)

LIBS += $(ROOT_LIBS)
LIBS += $(PMTCalib_lib)

src_files=$(wildcard *.cxx)

exec=$(src_files:.cxx=.exe)

all: $(exec)

%.exe: %.cxx
	$(CXX) $(INCLUDES) $(LIBS) -o $@ $<
