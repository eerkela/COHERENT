#ifndef CALSTRUCT_H
#define CALSTRUCT_H

#include <TH1.h>
#include <TGraphErrors.h>

struct ParWindow {
	Double_t low;
	Double_t high;
};

struct PeakInfo {
	Double_t energy;
	Double_t mu;
	Double_t muErr;
	Double_t sigma;
	Double_t sigmaErr;
	Double_t count;
};

struct FitResults {
	Double_t offset;
	Double_t offsetErr;
	Double_t slope;
	Double_t slopeErr;
};

struct FitInfo {
	std::vector<Double_t> peakEnergies;
	std::string fitFunc;
	std::map<Int_t, Double_t> fitPars;
	std::map<Int_t, ParWindow> fitParLimits;
	ParWindow fitWindow;
	Double_t backgroundRange;
};

struct Measurement {
	Double_t val;
	Double_t err;
};

#endif
