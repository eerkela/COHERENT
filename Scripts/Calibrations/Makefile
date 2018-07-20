.PHONY: all clean 

all: Calibration

Calibration: Calibration.o
	g++ $(shell root-config --cflags) -o Calibration Calibration.o $(shell root-config --libs)

.cc.o:
	g++ $(shell root-config --cflags) -c $<

clean:
	rm -f Calibration *.o dict.cc *.pcm *.rootmap *.dylib

