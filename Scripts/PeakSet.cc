#include <set>

#include "CalStructs.h"
#include "PeakSet.h"

/* This class describes a wrapper class for a set of PeakInfo structs, allowing the set to be
searchable by peak energy.  It also allows the user to easily insert/update information stored in
the set using the put method.
*/

PeakSet::PeakSet(std::vector<Double_t> energies) {
/* Constructor: creates a new PeakSet from a vector of peak energies

Accepts:
	vector<Double_t> energies: a vector of peak energies to be included in the PeakSet.

Returns:
	a PeakSet containing PeakInfo structs with the supplied energies.  All parameters besides
		the energy for each peak is initialized to null.

*/
	for(Double_t en : energies) {
		PeakInfo curr;
		curr.energy = en;
		this->peaks.insert(curr);
	}
}

PeakSet::PeakSet() {
/* default constructor: creates a PeakSet with no elements. */
}

void PeakSet::put(PeakInfo info) {
/* puts a new peak into the set or updates it if it is already present

Accepts:
	PeakInfo info: the peak to be added to the set.  PeakInfo defined at line 12 of 
		CalStructs.h

*/
	PeakInfo peak = this->get(info.energy);
	if (peak.energy == -1) {
		this->peaks.insert(info);
	} else {
		peak.mu = info.mu;
		peak.muErr = info.muErr;
		peak.sigma = info.sigma;
		peak.sigmaErr = info.sigmaErr;
		peak.count = info.count;
		peak.includeInCal = info.includeInCal;
		this->peaks.erase(peak);
		this->peaks.insert(peak);
	}
}

PeakInfo PeakSet::get(Double_t energy) {
/* gets the peak stored in this PeakSet with the provided energy

Accepts:
	Double_t energy: the energy (in keV) of the peak to be retrieved

Returns:
	PeakInfo corresponding to the peak that stored in this PeakSet with the desired energy.
		Returns a PeakInfo struct with energy = -1 if the set does not contain the peak.

*/
	for (PeakInfo pk : this->peaks) {
		if (pk.energy == energy) {
			return pk;
		}
	}
	// throw invalid_argument("peak not in set: " + std::to_string(energy) + " keV");
	PeakInfo notFound;
	notFound.energy = -1;
	return notFound;
}

std::set<PeakInfo> PeakSet::getSet() {
/* returns the actual set in which peaks are being stored

TODO:
	make PeakSet itself iterable so the user doesn't need to directly access the raw set.
 
*/
	return this->peaks;
}

PeakInfo PeakSet::remove(Double_t energy) {
/* removes the peak with the provided energy from the set

Accepts:
	Double_t energy: the energy (in keV) of the peak to be removed.

Returns:
	PeakInfo describing the removed set.  If the set does not contain the peak, returns a
		PeakInfo with energy = -1.

*/
	for (PeakInfo pk : this->peaks) {
		if (pk.energy == energy) {
			this->peaks.erase(pk);
			return pk;
		}
	}
	// throw invalid_argument("peak not in set: " + std::to_string(energy) + " keV");
	PeakInfo notFound;
	notFound.energy = -1;
	return notFound;
}

bool PeakSet::contains(Double_t energy) {
/* returns true if a PeakInfo struct exists in the set with the specified energy */
	for (PeakInfo pk : this->peaks) {
		if (pk.energy == energy) {
			return true;
		}
	}
	return false;
}

Int_t PeakSet::size() {
/* returns the current size of the set */
	return this->peaks.size();
}
