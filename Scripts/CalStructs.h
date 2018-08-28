#ifndef CALSTRUCT_H
#define CALSTRUCT_H

#include <map>
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
	bool includeInCal = true;

	bool operator < (const PeakInfo &other) const {
		return energy < other.energy;
	}
};

struct FitResults {
	Double_t offset;
	Double_t offsetErr;
	Double_t slope;
	Double_t slopeErr;

	//Double_t nonlinear;
	//Double_t nonlinearErr;
};

struct FitInfo {
	std::vector<Double_t> peakEnergies;
	std::string fitFunc;
	std::map<Int_t, Double_t> fitPars;
	std::map<Int_t, ParWindow> fitParLimits;
	ParWindow fitWindow;
	Double_t backgroundRange;
	std::vector<Double_t> excludeFromCal;
};

struct Measurement {
	Double_t val;
	Double_t err;
};

#endif
