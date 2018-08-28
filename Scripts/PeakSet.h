#ifndef PEAKSET_H
#define PEAKSET_H

#include "CalStructs.h"

class PeakSet {
private:
	std::set<PeakInfo> peaks;
	Double_t firstPeakEnergy;
public:
	PeakSet(std::vector<Double_t> energies);
	PeakSet();
	void put(PeakInfo info);
	PeakInfo get(Double_t energy);
	std::set<PeakInfo> getSet();
	PeakInfo remove(Double_t energy);
	PeakInfo getFirstPeak();
	bool contains(Double_t energy);
	Int_t size();
};

#endif
