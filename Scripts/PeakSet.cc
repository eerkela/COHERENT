#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <stdio.h>
#include <algorithm>
#include <TApplication.h>
#include <TRint.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TChain.h>
#include <TMultiGraph.h>
#include <TH2D.h>
#include <TLine.h>
#include <TSpectrum.h>

#include "CalStructs.h"
#include "PeakSet.h"

// This is a custom data structure emulating the behavior of a set of PeakInfo structs.  The set is
// searchable by peak energy and will never contain two peaks with the same energy.

PeakSet::PeakSet(std::vector<Double_t> energies) {
	// Constructor: makes a new set of PeakInfo structs from a vector of energies.  All
	// paramaters of the contained PeakInfo structs besides energy are initialized to null.
	for(Int_t i = 0; i < energies.size(); i++) {
		PeakInfo curr;
		curr.energy = energies[i];
		this->put(curr);
	}
}

PeakSet::PeakSet() {
	// Standard Constructor: creates a PeakSet with no elements.
}

void PeakSet::put(PeakInfo info) {
	// puts a new PeakInfo struct into the set or updates it if it already exists
	Int_t index = this->indexOf(info.energy);
	if (index == -1) {
		this->peakPars.push_back(info);
		std::cout << std::endl;
		std::cout << "put inserted new peak: peak " << this->peakPars.size() << std::endl;
		std::cout << "energy: " << info.energy << std::endl;
		std::cout << "positions: " << info.mu << std::endl;
		std::cout << std::endl;
	} else {
		this->peakPars[index].mu = info.mu;
		this->peakPars[index].muErr = info.muErr;
		this->peakPars[index].sigma = info.sigma;
		this->peakPars[index].sigmaErr = info.sigmaErr;
		this->peakPars[index].count = info.count;
	}
}

PeakInfo PeakSet::get(Double_t energy) {
	// returns the PeakInfo struct in the set with the specified energy
	Int_t index = this->indexOf(energy);
	/*if (index == -1) {
		throw invalid_argument("peak not in set: " + std::to_string(energy) + " keV");
	}*/
	return this->peakPars[index];
}

PeakInfo PeakSet::getAtIndex(Int_t index) {
	return this->peakPars[index];
}

PeakInfo PeakSet::getHighestEnergyPeak() {
	// returns the PeakInfo for the highest energy peak.  This is useful because the highest
	// energy peak must be pinned to allow for linear extrapolation of other peaks.
	PeakInfo max;
	for (Int_t i = 0; i < this->peakPars.size(); i++) {
		PeakInfo curr = this->peakPars[i];
		if (curr.energy > max.energy) {
			max = curr;
		}
	}
	return max;
}

Int_t PeakSet::indexOf(Double_t energy) {
	// returns index of the PeakInfo struct with the specified energy (returns -1 if not present)
	for (Int_t i = 0; i < this->peakPars.size(); i++) {
		PeakInfo curr = this->peakPars[i];
		if (curr.energy == energy) {
			return i;
		}
	}
	return -1;
}

bool PeakSet::contains(Double_t energy) {
	// returns true if a PeakInfo struct exists in the set with the specified energy
	return this->indexOf(energy) != -1;
}

Int_t PeakSet::size() {
	return this->peakPars.size();
}
