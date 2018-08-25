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

// This class describes the analysis engine itself, which handles all the heavy lifting of the 
// characterization program. 

Double_t PeakFinder::snapToMax(TH1D *h, Int_t pos, Double_t low, Double_t high) {
	h->GetXaxis()->SetRangeUser(low, high);
	Double_t maxPos = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
	h->GetXaxis()->SetRangeUser(0, this->getOverflowPos());
	return maxPos;
}

bool PeakFinder::isNumber(std::string input) {
	for (Int_t i = 0; i < input.length(); i++) {
		char c = input[i];
		if (!isdigit(c)) {
			return false;
		}
	}
	return input.length() != 0;
}


PeakFinder::PeakFinder(TCanvas* canvas, TChain *c, PeakSet peaks, Double_t time, std::string channel, TApplication *app) {
// Constructor: takes in a data chain c, a PeakSet containing each of the peaks to be
// analyzed, the total time within the data run, and the channel to be used.  This
// constructor then draws an initial histogram to guess the location of the highest energy 
// peak, asks for user confirmation/override, then redraws the final histogram with a 
// guaranteed 500 bins below the confirmed highest energy peak position.  This is done to 
// enhance fit quality.
	this->canvas = canvas;
	this->data = c;
	this->time = time;
	this->channel = channel;
	this->peaks = peaks;
	
	Int_t numBins = 16384; // 2^14
	TCanvas *tempCanvas = new TCanvas("tempCanvas", "tempCanvas");
	gPad->SetLogy();

	Double_t overflowPos = 1.01 * this->data->GetMaximum("energy");
	TH1D *hTemp = new TH1D("hTemp", "Pinning Highest Energy Peak", numBins, 0, overflowPos);
	hTemp->GetXaxis()->SetTitle("Uncalibrated Energy");
	hTemp->GetYaxis()->SetTitle("Count");
	this->data->Draw("energy >> hTemp", this->channel.c_str());

	TH1D* hSmoothed = (TH1D*) hTemp->Clone();
	hSmoothed->Smooth(1);
	TSpectrum* s = new TSpectrum();
	Int_t nFound = s->Search(hSmoothed, 2, "", 0.0001);
	while (nFound > 7) {
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
	pos = this->snapToMax(hTemp, pos, 0.95 * (Double_t) pos, 1.05 * (Double_t) pos);
	hTemp->Draw();
	
	hTemp->GetXaxis()->SetRangeUser(0, 2 * pos);
	TLine *line = new TLine(pos, 0, pos, hTemp->GetBinContent(hTemp->GetMaximumBin()));
	line->SetLineColor(kRed);
	line->Draw();

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
		pos = this->snapToMax(hTemp, pos, 0.95 * pos, 1.05 * pos);
	}
	tempCanvas->Close();*/

	Double_t normPos = pos / overflowPos;
	
	numBins = (Int_t) (500.0 / normPos);
	this->numBins = numBins;
	TH1D *h = new TH1D("h", "Uncalibrated Spectrum", numBins, 0, overflowPos);
	this->data->Draw("energy >> h", channel.c_str());
	this->rawPlot = h;
	
	PeakInfo firstPeakInfo = this->peaks.getFirstPeak();
	firstPeakInfo.mu = pos;
	firstPeakInfo.count = h->GetBinContent(h->FindBin(pos));
	this->pinnedPeak = firstPeakInfo;
	this->peaks.put(firstPeakInfo);

	delete hTemp;
}

PeakInfo PeakFinder::findPeak(Double_t energy) {
// Estimates a peak's location by linear extrapolation from the confirmed location of the
// highest energy peak as found in the constructor above.  This function then snaps the
// estimate to the local maximum within +-5% of the estimated position and returns the 
// corresponding position as a final estimate.
	Double_t pinnedEnergy = this->pinnedPeak.energy;
	Double_t pinnedPosition = this->pinnedPeak.mu;
	Double_t scaleFactor = pinnedPosition / pinnedEnergy;
	
	Double_t pos = energy * scaleFactor;
	pos = this->snapToMax(this->rawPlot, pos, 0.95 * pos, 1.05 * pos);
	
	PeakInfo peak;
	if (this->peaks.indexOf(energy) != -1) {
		peak = this->peaks.get(energy);
	}
	peak.energy = energy;
	peak.mu = pos;
	peak.count = this->rawPlot->GetBinContent(this->rawPlot->FindBin(pos));
	this->peaks.put(peak);

	return peak;
}

FitResults PeakFinder::backEst(ParWindow win, Double_t range, std::string fitFunc) {
// This function provides an estimate of the background underneath a peak by counting
// inwards from the edges of the provided ParWindow.  It first translates the window's edges
// into bin numbers, then counts inwards from the window edges a number of bins determined
// by the provided range parameter.  It pushes the contents of these counted bins into a
// TGraph, whereupon they may be fit according to the supplied fitFunc, the results of
// which are returned as a FitResults struct.

// The range parameter works as a proportion of the provided window to be considered in
// background estimation.  1.00 means 100% of the window will be included in the background
// estimation and 0.00 means none will.  This number should be chosen to maximize the number
// of points considered in the background fit, but should not include data from the peak 
// itself.
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
// This function provides a fit to a region specified by the provided ParWindow.  It 
// automatically estimates background and can fit multiple peaks at once if necessary.
// The user must supply energies and initial guesses for each peak to be fitted.

// Make it so that every guassian beyond the first has only 2 parameters.  The resolution
// should not change.  It's a feature of the detector, not the peak.

	TCanvas *c = this->canvas;
	TH1D *h = this->rawPlot;
	Int_t pos = this->findPeak(info.peakEnergies[0]).mu;
	Int_t count = h->GetBinContent(h->FindBin(pos));
	
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

	h->Fit(fit, "R+L"); // Q suppresses some output, ME does "better" fits and error estimation
	//h->Draw("SAME");
	h->Draw();
	gPad->SetLogy();

	PeakInfo firstPeak;
	firstPeak.energy = info.peakEnergies[0];
	firstPeak.count = fit->GetParameter(0);
	firstPeak.mu = fit->GetParameter(1);
	firstPeak.muErr = fit->GetParError(1);
	firstPeak.sigma = fit->GetParameter(2);
	firstPeak.sigmaErr = fit->GetParError(2);
	firstPeak.includeInCal = true;
	for (Double_t en : info.excludeFromCal) {
		if (en == firstPeak.energy) {
			firstPeak.includeInCal = false;
		}
	}
	this->peaks.put(firstPeak);

	std::cout << "fitted " << firstPeak.energy << std::endl;
	std::cout << "include in cal? " << firstPeak.includeInCal << std::endl;

	for (Int_t i = 1; i < info.peakEnergies.size(); i++) {
		PeakInfo nextPeak;
		nextPeak.energy = info.peakEnergies[i];
		nextPeak.count = fit->GetParameter(2 * i + 1);
		nextPeak.mu = fit->GetParameter(2 * i + 2);
		nextPeak.muErr = fit->GetParError(2 * i + 2);
		nextPeak.sigma = fit->GetParameter(2);
		nextPeak.sigmaErr = fit->GetParError(2);
		nextPeak.includeInCal = true;
		for (Double_t en : info.excludeFromCal) {
			if (en == nextPeak.energy) {
				nextPeak.includeInCal = false;
			}
		}
		this->peaks.put(nextPeak);

		std::cout << "fitted " << nextPeak.energy << std::endl;
		std::cout << "include in cal? " << nextPeak.includeInCal << std::endl;
	}
	delete fit;
}

FitResults PeakFinder::findCalibration() {

	std::cout << "Peaks included in calibration:" << std::endl;

	Int_t numPeaks = 0;
	std::vector<Double_t> expEs;
	std::vector<Double_t> fitEs;
	std::vector<Double_t> fitEErrs;
	for (Int_t i = 0; i < this->peaks.size(); i++) {
		PeakInfo currPeak = this->peaks.getAtIndex(i);
		if (currPeak.includeInCal) {
			expEs.push_back(currPeak.energy);
			fitEs.push_back(currPeak.mu);
			fitEErrs.push_back(currPeak.muErr);
			numPeaks++;

			std::cout << "peak " << numPeaks << std::endl;
			std::cout << "expected energy: " << currPeak.energy << std::endl;
			std::cout << "fitted position: " << currPeak.mu << std::endl;
			std::cout << std::endl;
		}
	}

	TF1 *calFit = new TF1("calFit", "pol1", 0, this->getOverflowPos());
	this->calPlot = new TGraphErrors(numPeaks, &expEs[0], &fitEs[0], 0, &fitEErrs[0]);
	this->calPlot->Fit("calFit", "", "", 0, 0);

	FitResults pars;
	pars.slope = calFit->GetParameter(1);
	pars.slopeErr = calFit->GetParError(1);
	pars.offset = calFit->GetParameter(0);
	pars.offsetErr = calFit->GetParError(0);

	this->calibration = pars;
	return pars;
}

Measurement PeakFinder::calibrate(Measurement uncalibrated) {
	FitResults calPars = this->calibration;
	Double_t calEnergy = (uncalibrated.val - calPars.offset) / calPars.slope;
	
	// Generalized Error = Sqrt(term1 + term2 + term3)
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
	return this->backPlots;
}

FitResults PeakFinder::getCalibration() {
	return this->calibration;
}

TGraphErrors* PeakFinder::getCalPlot() {
	return this->calPlot;
}

Double_t PeakFinder::getOverflowPos() {
// returns the maximum uncalibrated position contained in the main histogram.
	return 1.01 * this->data->GetMaximum("energy");
}

PeakSet PeakFinder::getPeakSet() {
// returns the PeakSet being used by this PeakFinder.
	return this->peaks;
}

TH1D* PeakFinder::getRawPlot() {
// returns the main histogram drawn for this PeakFinder
	return this->rawPlot;
}
