SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.c=.o)

.PHONY: all clean 

all: Calibration

Calibration: $(OBJECTS)
	g++ $(shell root-config --cflags) -o Calibration $(OBJECTS) $(shell root-config --libs) -lSpectrum

.cc.o:
	g++ $(shell root-config --cflags) -c $<

clean:
	rm -f Calibration *.o dict.cc *.pcm *.rootmap *.dylib

