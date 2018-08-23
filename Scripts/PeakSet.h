#ifndef PEAKSET_H
#define PEAKSET_H

#include "CalStructs.h"

class PeakSet {
private:
	std::vector<PeakInfo> peakPars;
public:
	PeakSet(std::vector<Double_t> energies);
	PeakSet();
	void put(PeakInfo info);
	PeakInfo get(Double_t energy);
	PeakInfo getAtIndex(Int_t index);
	PeakInfo getHighestEnergyPeak();
	Int_t indexOf(Double_t energy);
	bool contains(Double_t energy);
	Int_t size();
};

#endif
