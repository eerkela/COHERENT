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
#include "PeakFinder.h"

/*
This class describes the core of the analysis engine itself, which handles all the heavy lifting 
of the characterization program.  This includes peak fitting, background estimation, and 
generating the calibration itself.
*/

Double_t PeakFinder::snapToMax(TH1D *h, Double_t low, Double_t high) {
/* returns the location of the local maximum within the provided window

Accepts:
	TH1D *h: the histogram to scan over.
	Double_t low: the low edge of the snap window (in units of h's x-axis)
	Double_t high: the high edge of the snap window (in units of h's x-axis)

Returns:
	Double_t describing the location (in units of h's x-axis) of the bin with the 
		maximum contents in the window.

*/
	h->GetXaxis()->SetRangeUser(low, high);
	Double_t maxPos = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
	h->GetXaxis()->SetRangeUser(0, this->getOverflowPos());
	return maxPos;
}

bool PeakFinder::isNumber(std::string input) {
/* tests if a string consists only of numbers

Accepts:
	string input: string to be analyzed.

Returns:
	true if string contains nothing but numerical digits.

*/
	for (Int_t i = 0; i < input.length(); i++) {
		char c = input[i];
		if (!isdigit(c)) {
			return false;
		}
	}
	return input.length() != 0;
}

PeakFinder::PeakFinder(Double_t pinnedEnergy, TChain *c, std::string channel, TApplication *app) {
/* Constructor: builds a PeakFinder object

Accepts:
	Double_t pinnedEnergy: the energy of the pinned peak, which will be used to estimate
		the other peak locations.  As implemented, this should always be the energy
		of the 208Tl peak (2614.511 keV)
	TChain *c: pointer to the TChain containing the raw data to be analyzed.
	string channel: the digitizer channel for which data is to be analyzed.
	TApplication *app: a pointer to a ROOT interactive application, to allow for user input
		and manipulation of plots.

Returns:
	A PeakFinder object initialized with the relevant information to begin analysis.

*/
	this->data = c;
	this->channel = channel;
	
	Int_t numBins = 16384; // 2^14
	TCanvas *tempCanvas = new TCanvas("tempCanvas", "tempCanvas");
	gPad->SetLogy();

	Double_t overflowPos = 1.01 * this->data->GetMaximum("energy");
	TH1D *hTemp = new TH1D("hTemp", "Pinning Highest Energy Peak", numBins, 0, overflowPos);
	hTemp->GetXaxis()->SetTitle("Uncalibrated Energy");
	hTemp->GetYaxis()->SetTitle("Count");
	this->data->Draw("energy >> hTemp", this->channel.c_str());

	// must identify the position of the pinned peak, so that other peaks may be estimated.
	TH1D* hSmoothed = (TH1D*) hTemp->Clone();
	hSmoothed->Smooth(1);
	TSpectrum* s = new TSpectrum();
	Int_t nFound = s->Search(hSmoothed, 2, "", 0.0001);
	while (nFound > 7) {   // arbitrary threshold.  7 seems to work well.
		hSmoothed->Rebin(2);
		nFound = s->Search(hSmoothed, 2, "", 0.0001);
	}
	delete hSmoothed;

	Double_t* sPeaks = s->GetPositionX();
	Double_t TlGuess = 0;
	for (Int_t k = 0; k < nFound; k++) {
		Double_t currPeak = sPeaks[k];
		if (currPeak < 0.9 * overflowPos && currPeak > TlGuess) {
			TlGuess = currPeak;
		}
	}
	Double_t pos = TlGuess;
	pos = this->snapToMax(hTemp, 0.95 * (Double_t) pos, 1.05 * (Double_t) pos);
	hTemp->Draw();
	
	hTemp->GetXaxis()->SetRangeUser(0, 2 * pos);
	TLine *line = new TLine(pos, 0, pos, hTemp->GetBinContent(hTemp->GetMaximumBin()));
	line->SetLineColor(kRed);
	line->Draw();

	// user must manually verify the found peak position.
	/*std::cout << "VISUAL CHECK" << std::endl;
	std::cout << "estimated position for highest energy peak: " << pos << std::endl;
	std::cout << "(type .q to continue)" << std::endl;
	app->Run(true);
	std::cout << "does this make sense? (y/n) ";
	std::string response;
	std::cin >> response;
	while (response != "y" && response != "n") {
		std::cout << "error: cannot interpret response \"" + response + "\"" << std::endl;
		std::cout << "estimated position for first peak: " << pos << std::endl;
		std::cout << "does this make sense? (y/n) ";
		std::cin >> response;
	}
	if (response == "n") {
		std::cout << "new peak position: ";
		std::cin >> response;
		while (!this->isNumber(response)) {
			std::cout << "error: response \"" + response + "\" isn't a number" << std::endl;
			std::cout << "new peak position: ";
			std::cin >> response;
		}
		pos = stod(response);
		pos = this->snapToMax(hTemp, 0.95 * pos, 1.05 * pos);
	} */
	delete tempCanvas;

	// histogram is redrawn with a constant 500 bins below first peak position
	// this helps stabilize fits
	Double_t normPos = pos / overflowPos;
	numBins = (Int_t) (500.0 / normPos);
	this->numBins = numBins;
	TH1D *h = new TH1D("h", "Uncalibrated Spectrum", numBins, 0, overflowPos);
	this->data->Draw("energy >> h", channel.c_str(), "goff");
	this->rawPlot = h;
	
	this->pinnedPeak.energy = pinnedEnergy;
	this->pinnedPeak.mu = pos;
	this->pinnedPeak.count = h->GetBinContent(h->FindBin(pos));
	this->peaks.put(this->pinnedPeak);

	delete hTemp;
}

void PeakFinder::addPeakToSet(PeakInfo info) {
/* adds a peak to this analyzer's set

Accepts:
	PeakInfo info: peak to be inserted into this analyzer's PeakSet

*/
	this->peaks.put(info);
}

PeakInfo PeakFinder::findPeak(Double_t energy) {
/* estimates a peak's position by linear extrapolation from pinned peak

Accepts:
	Double_t energy: the energy of the peak to be found

Returns:
	PeakInfo describing the parameters of the found peak

*/
	Double_t pinnedEnergy = this->pinnedPeak.energy;
	Double_t pinnedPosition = this->pinnedPeak.mu;
	Double_t pos = energy * pinnedPosition / pinnedEnergy;
	// automatically snaps to local max within +/- 5% of estimated position
	pos = this->snapToMax(this->rawPlot, 0.95 * pos, 1.05 * pos);

	PeakInfo peak;
	if (this->peaks.contains(energy)) {
		peak = this->peaks.get(energy);
	} else {
		peak.energy = energy;
	}
	peak.mu = pos;
	peak.count = this->rawPlot->GetBinContent(this->rawPlot->FindBin(pos));
	this->peaks.put(peak);

	return peak;
}

FitResults PeakFinder::backEst(ParWindow win, Double_t range, std::string fitFunc) {
/* Estimates the background underneath a peak by counting inwards from the window's edges.

Accepts:
	ParWindow win: the window identifying the lower and upper bounds to be used in the
		background estimation algorithm.  ParWindow defined at line 7 in CalStructs.h
	Double_t range: the proportion of the provided window to be considered in the background 
		estimation (0 < range < 1.0).  This number should be chosen to maximize the 
		number of points considered background, but should not include any data from
		the region of the peak itself.  In Calibration.cc, the "back" option will
		allow for a visual inspection so that this parameter may be tuned.
	string fitFunc: describes the function to be used for background fitting.  As implemented,
		this should always be an exponential ("expo").

Returns:
	FitResults describing the parameters (with errors) of the background curve.  As 
		implemented, this is always of the form: "exp([offset] + [slope] * x)

*/
	TH1D *h = this->rawPlot;
	
	Int_t lowBin = h->FindBin(win.low);
	Int_t highBin = h->FindBin(win.high);
	Int_t overallBinRange = highBin - lowBin;
	Int_t backWindowRange = (Int_t) ((range / 2.0) * (Double_t) overallBinRange);

	Float_t backVals[2 * backWindowRange];
	Float_t backValErrs[2 * backWindowRange];
	Float_t backPos[2 * backWindowRange];

	for (Int_t i = 0; i < backWindowRange; i++) {
		Int_t bin = i + lowBin;
		backVals[i] = h->GetBinContent(bin);
		backValErrs[i] = TMath::Sqrt(h->GetBinContent(bin));
		backPos[i] = h->GetBinCenter(bin);
	}
	for (Int_t i = 0; i < backWindowRange; i++) {
		Int_t bin = (highBin - backWindowRange) + (i + 1);
		backVals[i + backWindowRange] = h->GetBinContent(bin);
		backValErrs[i + backWindowRange] = TMath::Sqrt(h->GetBinContent(bin));
		backPos[i + backWindowRange] = h->GetBinCenter(bin);
	}

	TGraphErrors* backGraph = new TGraphErrors(2 * backWindowRange, backPos, backVals, 0, backValErrs);
	TF1 backFit = TF1("backFit", fitFunc.c_str(), win.low, win.high);
	backGraph->Fit("backFit");

	FitResults pars;
	pars.offset = backFit.GetParameter(0);
	pars.slope = backFit.GetParameter(1);
	
	backGraph->SetTitle(("Background Estimation Graph (" + fitFunc + ")").c_str());
	backGraph->GetYaxis()->SetTitle("Count");
	backGraph->GetXaxis()->SetTitle("Uncalibrated Energy");
	backGraph->SetMarkerStyle(4);
	backGraph->SetMarkerSize(0.5);
	gPad->SetLogy();

	this->backPlots.push_back(backGraph);

	return pars;
}

void PeakFinder::fit(FitInfo info) {
/* performs a fit to a region on this PeakFinder's histogram.  Automatically estimates background 
and can fit to an arbitrary number of peaks.

Accepts:
	FitInfo info: a struct containing all of the necessary information for the desired fit,
		including the fit function, peak energies, parameter guesses, parameter limits, 
		fit window, and background range.  FitInfo defined at line 36 of CalStructs.h

Returns:
	nothing.  Updates values in this PeakFinder's PeakSet based on results of fit.

*/
	TCanvas *c = new TCanvas("c", "c");
	TH1D *h = this->rawPlot;
	Int_t pos = this->findPeak(info.peakEnergies[0]).mu;
	Int_t count = h->GetBinContent(h->FindBin(pos));

	/*fit function must be of form:

	[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[2])^2) + ...  + exp([n-1]+[n]*x);
	
	note that each peak has the same width set by parameter [2], though there may be
	an arbitrary amount of gaussian peaks to consider.
	*/
	
	TF1 *fit = new TF1("fit", info.fitFunc.c_str(), info.fitWindow.low, info.fitWindow.high);
	Int_t highestKey = 0;
	for (std::pair<Int_t, Double_t> parGuess : info.fitPars) {
		fit->SetParameter(parGuess.first, parGuess.second);
		if (parGuess.first > highestKey) {
			highestKey = parGuess.first;
		}
	}
	FitResults backPars = this->backEst(info.fitWindow, info.backgroundRange, "expo");
	fit->SetParameter(highestKey + 1, backPars.offset);
	fit->SetParameter(highestKey + 2, backPars.slope);
	
	for (std::pair<Int_t, ParWindow> lims : info.fitParLimits) {
		fit->SetParLimits(lims.first, lims.second.low, lims.second.high);
	}

	h->Fit(fit, "R+L"); // Q suppresses some output, ME = "better" fits and error estimation
	h->Draw("goff");
	gPad->SetLogy();

	PeakInfo firstPeak;
	firstPeak.energy = info.peakEnergies[0];
	firstPeak.count = fit->GetParameter(0);
	firstPeak.mu = fit->GetParameter(1);
	firstPeak.muErr = fit->GetParError(1);
	firstPeak.sigma = fit->GetParameter(2);
	firstPeak.sigmaErr = fit->GetParError(2);
	for (Double_t en : info.excludeFromCal) {
		if (en == firstPeak.energy) {
			firstPeak.includeInCal = false;
		}
	}
	this->peaks.put(firstPeak);

	for (Int_t i = 1; i < info.peakEnergies.size(); i++) {
		PeakInfo nextPeak;
		nextPeak.energy = info.peakEnergies[i];
		nextPeak.count = fit->GetParameter(2 * i + 1);
		nextPeak.mu = fit->GetParameter(2 * i + 2);
		nextPeak.muErr = fit->GetParError(2 * i + 2);
		nextPeak.sigma = fit->GetParameter(2);
		nextPeak.sigmaErr = fit->GetParError(2);
		for (Double_t en : info.excludeFromCal) {
			if (en == nextPeak.energy) {
				nextPeak.includeInCal = false;
			}
		}
		this->peaks.put(nextPeak);
	}
	delete fit;
	delete c;
}

FitResults PeakFinder::findCalibration() {
/* finds a calibration for this detector

Accepts:
	nothing.  All relevant information is stored in this PeakFinder's PeakSet

Returns:
	FitResults describing the final parameters of the found calibration

*/
	std::vector<Double_t> expEs;
	std::vector<Double_t> fitEs;
	std::vector<Double_t> fitEErrs;
	for (PeakInfo pk : this->peaks.getSet()) {
		if (pk.includeInCal) {
			expEs.push_back(pk.energy);
			fitEs.push_back(pk.mu);
			fitEErrs.push_back(pk.muErr);
		}
	}

	// calibration is linear for now
	TF1 *calFit = new TF1("calFit", "pol1", 0, this->data->GetMaximum("energy"));
	this->calPlot = new TGraphErrors(expEs.size(), &expEs[0], &fitEs[0], 0, &fitEErrs[0]);
	this->calPlot->Fit("calFit", "R+");

	FitResults pars;
	pars.offset = calFit->GetParameter(0);
	pars.offsetErr = calFit->GetParError(0);
	pars.slope = calFit->GetParameter(1);
	pars.slopeErr = calFit->GetParError(1);
	//pars.nonlinear = calFit->GetParameter(2);
	//pars.nonlinear = calFit->GetParError(2);

	this->calibration = pars;
	return pars;
}

Measurement PeakFinder::calibrate(Measurement uncalibrated) {
/* converts a given measurement to a real energy scale based on this PeakFinder's found 
calibration.  Also performs error analysis.

Accepts:
	Measurement uncalibrated: the measurement (value and error) to be calibrated.  
		Measurement defined at line 46 of CalStructs.h

Returns:
	Measurement describing the calibrated measurement, including error.

*/
	FitResults calPars = this->calibration;
	Double_t calEnergy = (uncalibrated.val - calPars.offset) / calPars.slope;
	
	// Generalized Error Formula = Sqrt(term1 + term2 + term3)
	// term1: ((partial of calEnergy wrt uncalibratedEnergy) * uncalibratedEnergyErr)^2
	// term2: ((partial of calEnergy wrt offset) * offsetErr)^2
	// term3: ((partial of calEnergy wrt slope) * slopeErr)^2
	Double_t term1 = TMath::Power(uncalibrated.err / calPars.slope, 2);
	Double_t term2 = TMath::Power(calPars.offsetErr / calPars.slope, 2);
	Double_t term3 = TMath::Power(calEnergy * (calPars.slopeErr / calPars.slope), 2);
	Double_t calEnergyErr = TMath::Sqrt(term1 + term2 + term3);
		
	Measurement calibrated;
	calibrated.val = calEnergy;
	calibrated.err = calEnergyErr;
	return calibrated;
}

std::vector<TGraphErrors*> PeakFinder::getBackgroundPlots() {
/* returns a vector of TGraphErrors* representing the background plots generated by backEst */
	return this->backPlots;
}

FitResults PeakFinder::getCalibration() {
/* returns this PeakFinder's current calibration as generated by findCalibration */
	return this->calibration;
}

TGraphErrors* PeakFinder::getCalPlot() {
/* returns the TGraphErrors* used to find the calibration in findCalibration */
	return this->calPlot;
}

Double_t PeakFinder::getOverflowPos() {
/* returns the uncalibrated energy corresponding to the maximum bin in the raw histogram */
	return 1.01 * this->data->GetMaximum("energy");
}

PeakSet PeakFinder::getPeakSet() {
/* returns the PeakSet being used to store peak information for this PeakFinder */
	return this->peaks;
}

PeakInfo PeakFinder::getPinnedPeak() {
/* returns the PeakInfo corresponding to the pinned peak as found in constructor */
	return this->pinnedPeak;
}

TH1D* PeakFinder::getRawPlot() {
/* returns the histogram containing raw data being analyed by this PeakFinder */
	return this->rawPlot;
}
