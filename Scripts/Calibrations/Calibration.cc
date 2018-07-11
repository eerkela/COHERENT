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
#include <stdexcept>

/*
This script (built using the ROOT Data Analysis Framework from CERN) will analyze a 
characterization suite for a Thalium-doped Sodium Iodide (NaI[Tl]) crystal scintillator, 
the data for which has been collected according to the University of Washington COHERENT group's 
official protocol.  Baseline, it will translate root-readable raw digitizer output for a crystal 
and calibrate the collected data to a real energy scale, but more detailed analysis may be done 
using one of the options described below.

Modes:
"pos"	will execute calibration algorithm for position data.
"volt"	will execute calibration algorithm for voltage data.

Options:
"cal" 		will display the graphs used to generate the calibrations for each run.
"gain"		will display a graph of the calculated calibration slope (gain) as a function of
		the dependent variable specified by the mode.
"over"		will calibrate the spectra then overlay the calibrated histograms with real 
		energies for qualitative comparison.
"res"		will display residues from the calibration algorithm for each run. Residues are 
		the % error by which the calibrated energies deviate from expected.
"sig" 		will display fractional sigmas for each peak in each run.
"AE"		will display an amplitude / energy vs energy plot, useful for pulse shape 
		discrimination.
"backFits"	will display the background fits produced by the backEst function for all peaks.
"noise"		will display the noise wall energy vs the dependent variable set by the mode
		parameter."
"muon" (NYI)	will fit the muon peak.


Required Directory structure for Calibration to work:
Crystal Serial #/
	Characterization.cc
	Position/
		.../position_1/NaI_ET_run*.root
		.../position_2/NaI_ET_run*.root
		.../position_3/NaI_ET_run*.root
		.../position_4/NaI_ET_run*.root
		.../position_5/NaI_ET_run*.root
	Voltage/
		.../600_V/NaI_ET_run*.root
		.../700_V/NaI_ET_run*.root
		.../800_V/NaI_ET_run*.root
		.../900_V/NaI_ET_run*.root
		.../1000_V/NaI_ET_run*.root	


Notes and TODO:

	Data randomly deletes itself during execution.  This is almost certainly just a quirk 
	of root itself, but I can't figure out exactly how to get around it atm.
	
	FIR stuff removed because it turns out the detector response to peaking time is mostly
	flat in the region we're sensitive to.

	Continue making the script agnostic to low/high gain data.  This is mostly complete by
	now, but I haven't had the chance to test it yet.  Ideally, we can take this one version
	of the script and use it for both possible data sets since the underlying algorithm is
	the same between them.  However, all the peaks we're actually fitting are different, 
	which makes it more difficult to generalize.

	Fix Cs fit.  There's a few peaks underneath it that kind of screw it up.  Might want to
	fix some parameters.

	Figure out how to do math on TChains.  I want to store the calibrated data in the RESULTS
	struct somewhere, but I need to apply the calibration directly to the data in that case.

	Figure out how to get the visual check functioning correctly.

*/

using namespace std;

/*
   ###################
   #     STRUCTS     #
   ###################
*/

struct FitWindow {
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
	vector<Double_t> peakEnergies;
	string fitFunc;
	FitWindow fitWindow;
	Double_t backgroundRange;
};

struct Plots {
	TH1D *raw;
	TH1D *calibrated;
	TH1D *muPlot;

	TGraphErrors *calPlot;
	TGraphErrors *sigmaPlot;

};

struct Results {
	FitResults calPars;
	TChain calibratedData;
	
	Double_t noiseWall;
	Double_t bestPeakTime;
	Double_t lowERate;
};

/*
   ###################
   #     PEAKSET     #
   ###################
*/

class PeakSet {
private:
	vector<PeakInfo> peakPars;

public:
	PeakSet(vector<Double_t> energies) {
		for(Int_t i = 0; i < energies.size(); i++) {
			PeakInfo curr;
			curr.energy = energies[i];
			this->put(curr);
		}
	}
	
	PeakSet() {
		// do nothing.  PeakPars will be a vector with no elements.
	}

	void put(PeakInfo info) {
		Int_t index = this->indexOf(info.energy);
		if (index == -1) {
			this->peakPars.push_back(info);
		} else {
			this->peakPars[index].mu = info.mu;
			this->peakPars[index].muErr = info.muErr;
			this->peakPars[index].sigma = info.sigma;
			this->peakPars[index].sigmaErr = info.sigmaErr;
			this->peakPars[index].count = info.count;
		}
	}

	PeakInfo get(Double_t energy) {
		Int_t index = this->indexOf(energy);
		if (index == -1) {
			throw invalid_argument("peak not in set: " + to_string(energy) + " keV");
		}
		return this->peakPars[index];
	}

	PeakInfo getHighestEnergyPeak() {
		PeakInfo max;
		for (Int_t i = 0; i < this->peakPars.size(); i++) {
			PeakInfo curr = this->peakPars[i];
			if (curr.energy > max.energy) {
				max = curr;
			}
		}
		return max;
	}

	Int_t indexOf(Double_t energy) {
		for (Int_t i = 0; i < this->peakPars.size(); i++) {
			PeakInfo curr = this->peakPars[i];
			if (curr.energy == energy) {
				return i;
			}
		}
		return -1;
	
	}

	bool contains(Double_t energy) {
		return this->indexOf(energy) != -1;
	}

	Int_t size() {
		return this->peakPars.size();
	}
};


/*
   ######################
   #     PEAKFINDER     #
   ######################
*/

class PeakFinder {
// This class describes a set of PeakInfo structs that is searchable by their energies.
// An object of this type will house all the pertinent info collected for each peak
private:
	TChain *data;
	TH1D *plot;
	Double_t time;
	Int_t numBins;
	string channel;
	PeakSet peaks;
	PeakInfo pinnedPeak;

	//TApplication* app;

	Double_t snapToMax(TH1D *h, Int_t pos, Double_t low, Double_t high) {
		h->GetXaxis()->SetRangeUser(low, high);
		Double_t maxPos = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
		h->GetXaxis()->SetRangeUser(0, this->getOverflowPos());
		return maxPos;
	}
	
public:
	PeakFinder(TChain *c, PeakSet peaks, Double_t time, string channel) {
		// give me a vector of energies sorted in the order they will be searched in.
		//this->app = new TRint("app", 0, NULL);
		this->data = c;
		this->time = time;
		this->channel = channel;
		this->peaks = peaks;
		
		Int_t numBins = 10000;
		TCanvas *tempCanvas = new TCanvas("tempCanvas", "tempCanvas");
		gPad->SetLogy();

		Double_t overflowPos = 1.01 * this->data->GetMaximum("energy");
		TH1D *hTemp = new TH1D("hTemp", "Finding Max E", numBins, 0, overflowPos);
		this->data->Draw("energy >> hTemp", this->channel.c_str());

		Int_t binNum = numBins;
		Int_t count = 0;
		while ((Double_t) count < 9.0 * time && binNum > 0) {
			count += (Int_t) hTemp->GetBinContent(binNum);
			binNum--;
		}
		Double_t pos = hTemp->GetXaxis()->GetBinCenter(binNum);
		pos = this->snapToMax(hTemp, pos, 0.95 * (Double_t) pos, 1.05 * (Double_t) pos);
		hTemp->Draw();

		cout << "VISUAL CHECK" << endl;
		cout << "estimated position for highest energy peak: " << pos << endl;
		//this->app->Run(true);
		cout << "does this make sense? (y/n) ";
		string response;
		cin >> response;
		while (response != "y" && response != "n") {
			cout << "error: cannot interpret response \"" + response + "\"" << endl;
			cout << "estimated position for highest energy peak: ";
			         cout << pos << endl;
			cout << "does this make sense? (y/n) ";
			cin >> response;
		}
		if (response == "n") {
			cout << "new peak position: ";
			cin >> pos;
			pos = this->snapToMax(hTemp, pos, 0.95 * pos, 1.05 * pos);
		}
		delete hTemp;
		delete tempCanvas;

		Double_t normPos = pos / overflowPos;
		
		numBins = (Int_t) (500.0 / normPos);
		this->numBins = numBins;
		TH1D *h = new TH1D("h", "Uncalibrated Spectrum", numBins, 0, overflowPos);
		this->data->Draw("energy >> h", channel.c_str());
		this->plot = h;
		
		PeakInfo firstPeakInfo = this->peaks.getHighestEnergyPeak();
		firstPeakInfo.mu = pos;
		firstPeakInfo.count = h->GetBinContent(h->FindBin(pos));
		this->pinnedPeak = firstPeakInfo;
		this->peaks.put(firstPeakInfo);
	}

	Int_t findPeak(Double_t energy) {
		Double_t pinnedEnergy = this->pinnedPeak.energy;
		Double_t pinnedPosition = this->pinnedPeak.mu;
		Double_t scaleFactor = pinnedPosition / pinnedEnergy;

		//cout << "pinnedEnergy: " << pinnedEnergy << endl;
		//cout << "pinnedPosition: " << pinnedPosition << endl;
		//cout << "Scale Factor: " << scaleFactor << endl;
		
		Double_t pos = energy * scaleFactor;
		cout << "Estimated position: " << pos << endl;
		pos = this->snapToMax(this->plot, pos, 0.95 * pos, 1.05 * pos);
		
		PeakInfo peak = this->peaks.get(energy);
		peak.mu = pos;
		peak.count = this->plot->GetBinContent(this->plot->FindBin(pos));

		return pos;
	}

	FitResults backEst(FitWindow win, Double_t range, string fitFunc) {
		TH1D *h = this->plot;
		
		Int_t lowBin = h->FindBin(win.low);
		Int_t highBin = h->FindBin(win.high);
		Int_t overallBinRange = highBin - lowBin;
		Int_t backWindowRange = (Int_t) ((range / 2.0) * (Double_t) overallBinRange);

		Float_t backVals[2 * backWindowRange];
		Float_t backPos[2 * backWindowRange];

		for (Int_t i = 0; i < backWindowRange; i++) {
			Int_t bin = i + lowBin;
			backVals[i] = h->GetBinContent(bin);
			backPos[i] = h->GetBinCenter(bin);
		}
		for (Int_t i = 0; i < backWindowRange; i++) {
			Int_t bin = (highBin - backWindowRange) + (i + 1);
			backVals[i + backWindowRange] = h->GetBinContent(bin);
			backPos[i + backWindowRange] = h->GetBinCenter(bin);
		}

		TGraph *backGraph = new TGraph(2 * backWindowRange, backPos, backVals);
		TF1 backFit = TF1("backFit", fitFunc.c_str(), win.low, win.high);
		backGraph->Fit("backFit");

		FitResults pars;
		pars.offset = backFit.GetParameter(0);
		pars.slope = backFit.GetParameter(1);
		
		backGraph->GetYaxis()->SetTitle("Count");
		backGraph->GetXaxis()->SetTitle("Uncalibrated Energy");
		backGraph->SetMarkerStyle(4);
		backGraph->SetMarkerSize(0.5);

		backGraph->Draw();
		delete backGraph;
		return pars;
	}

	void fit(vector<Double_t> energies, FitWindow win, string fitFunc, 
	         vector<Double_t> pars, Double_t backRange) {
		TH1D *h = this->plot;
		Int_t pos = this->findPeak(energies[0]);
		Int_t count = h->GetBinContent(h->FindBin(pos));
		
		FitResults backPars = this->backEst(win, backRange, "expo");
		pars.push_back(backPars.offset);
		pars.push_back(backPars.slope);

		Double_t parArray[pars.size()];
		for (Int_t i = 0; i < pars.size(); i++) {
			parArray[i] = pars[i];
		}
		
		TF1 *fit = new TF1("fit", fitFunc.c_str(), win.low, win.high);
		fit->SetParameters(parArray);
		h->Fit("fit", "R+ll");

		for (Int_t i = 0; i < energies.size(); i++) {
			PeakInfo curr;
			curr.energy = energies[i];
			curr.count = fit->GetParameter(3 * i);
			curr.mu = fit->GetParameter(3 * i + 1);
			curr.muErr = fit->GetParameter(3 * i + 1);
			curr.sigma = fit->GetParameter(3 * i + 2);
			curr.sigmaErr = fit->GetParameter(3 * i + 2);
			this->peaks.put(curr);
		}
		delete fit;
	}

	Double_t getOverflowPos() {
		return 1.01 * this->data->GetMaximum("energy");
	}

	PeakSet getPeakSet() {
		return this->peaks;
	}

	TH1D* getPlot() {
		return this->plot;
	}

};

/*
   ################
   #     MAIN     #
   ################
*/

#define NUMFILES 5

void formatGraph(TGraphErrors *graph);
void getData(TChain *data[NUMFILES], string mode);
string getDepVar(string mode);
string getLabels(string mode, Int_t index);
Double_t getRunTime(TChain *c, string mode);

Int_t VOLTAGES[NUMFILES] = {600, 700, 800, 900, 1000}; 
//Int_t PEAKINGTIMES[5] = {50, 100, 200, 400, 800};
Int_t POSITIONS[NUMFILES] = {1, 2, 3, 4, 5};

vector<PeakSet> PEAKS(NUMFILES);
vector<Plots> PLOTS(NUMFILES);
vector<Results> RESULTS(NUMFILES);

void Calibration(string mode, string option) {

	// global variable gdirectory.
	// cd to groot before allocating memory.  groot->cd().

	//gStyle->SetOptFit(1111);
	TCanvas *fitCanvas = new TCanvas("Fits", "Fits", 1200, 800);
	fitCanvas->Divide(2, 3);
	TChain *data[NUMFILES]; 
	getData(data, mode);
	string channel = "channel==4";
	vector<FitInfo> fitPars;

	// Tl peak parameters:
	FitInfo TlPars;
	TlPars.peakEnergies.push_back(2615.0);
	TlPars.fitFunc = "[0] * exp(-0.5 * ((x - [1]) / 0.05 * [2])^2) + exp([3] + [4] * x)";
	TlPars.fitWindow.low = 0.9;
	TlPars.fitWindow.high = 1.1;
	TlPars.backgroundRange = 0.3;
	fitPars.push_back(TlPars);
	
	// K peak parameters:
	FitInfo KPars;
	KPars.peakEnergies.push_back(1460.0);
	KPars.fitFunc = "[0] * exp(-0.5 * ((x - [1]) / 0.05 * [2])^2) + exp([3] + [4] * x)";
	KPars.fitWindow.low = 0.85;
	KPars.fitWindow.high = 1.15;
	KPars.backgroundRange = 0.3;
	fitPars.push_back(KPars);

	// Cs peak parameters:
	FitInfo CsPars;
	CsPars.peakEnergies.push_back(661.6);
	CsPars.fitFunc = "[0] * exp(-0.5 * ((x - [1]) / 0.05 * [2])^2) + exp([3] + [4] * x)";
	CsPars.fitWindow.low = 0.8;
	CsPars.fitWindow.high = 1.2;
	CsPars.backgroundRange = 0.3;
	fitPars.push_back(CsPars);

	for (Int_t i = 0; i < NUMFILES; i++) {
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << endl;

		cout << "Beginning Calibration " << i + 1 << "..." << endl;
		fitCanvas->cd(i+1);
		gPad->SetLogy();

		Double_t time = getRunTime(data[i], mode);
		cout << "Run time in data chain: " << time << " seconds" << endl;

		PeakSet allPeaks;
		vector<Double_t> allPeakEnergies;
		for (Int_t j = 0; j < fitPars.size(); j++) {
			vector<Double_t> thisPeakEnergies = fitPars[j].peakEnergies;
			for (Int_t k = 0; k < thisPeakEnergies.size(); k++) {
				allPeakEnergies.push_back(thisPeakEnergies[k]);
				PeakInfo curr;
				curr.energy = thisPeakEnergies[k];
				allPeaks.put(curr);
			}
		}
		PeakFinder *analyzer = new PeakFinder(data[i], allPeaks, time, channel.c_str());
		
		for (Int_t j = 0; j < fitPars.size(); j++) {
			FitInfo pars = fitPars[j];

			PeakInfo estimate = analyzer->getPeakSet().get(pars.peakEnergies[0]);
			FitWindow actualWin;
			actualWin.low = pars.fitWindow.low * estimate.mu;
			actualWin.high = pars.fitWindow.high * estimate.mu;

			Int_t numPars = 0;  // counting # of pars needed in fit
			for (Int_t k = 0; k < pars.fitFunc.length(); k++) {
				char c = pars.fitFunc[k];
				if (c == '[') {
					numPars++;
				}
			}
			
			vector<Double_t> initialGuesses;
			for (Int_t k = 0; k < numPars / 3; k++) {
				initialGuesses.push_back(estimate.count);
				initialGuesses.push_back(estimate.mu);
				initialGuesses.push_back(estimate.mu);
			}
			
			analyzer->fit(pars.peakEnergies, actualWin, pars.fitFunc,
                                      initialGuesses, pars.backgroundRange);
		}
		PEAKS[i] = analyzer->getPeakSet();

		Float_t range = 1.15 * PEAKS[i].getHighestEnergyPeak().mu;
		PLOTS[i].raw = analyzer->getPlot();
		PLOTS[i].raw->GetXaxis()->SetRangeUser(0, range);

		// Collecting data:
		vector<Double_t> fitE;
		vector<Double_t> fitEErr;
		vector<Double_t> fitSigmaErr;
		vector<Double_t> fitSigma;
		for (Int_t j = 0; j < allPeakEnergies.size(); j++) {
			PeakInfo currPeak = PEAKS[i].get(allPeakEnergies[j]);
			fitE.push_back(currPeak.mu);
			fitEErr.push_back(currPeak.muErr);
			fitSigma.push_back(currPeak.sigma / currPeak.mu);
			fitSigmaErr.push_back(currPeak.sigmaErr / currPeak.mu);
		} 

		PLOTS[i].calPlot = new TGraphErrors(allPeakEnergies.size(), &allPeakEnergies[0], &fitE[0], 0, &fitEErr[0]);
		PLOTS[i].calPlot->SetTitle(getLabels(mode, i).c_str());
		PLOTS[i].calPlot->SetMarkerColor(i+1);
		if (i == 4) {PLOTS[i].calPlot->SetMarkerColor(i+2);}
		PLOTS[i].calPlot->SetMarkerStyle(21);
		PLOTS[i].calPlot->SetLineColor(1);
		PLOTS[i].calPlot->SetLineWidth(2);		
		
		PLOTS[i].sigmaPlot = new TGraphErrors(allPeakEnergies.size(), &allPeakEnergies[0], &fitSigma[0], 0, &fitSigmaErr[0]); 
		PLOTS[i].sigmaPlot->SetTitle(getLabels(mode, i).c_str());
		PLOTS[i].sigmaPlot->SetMarkerColor(i+1);
		PLOTS[i].sigmaPlot->SetLineColor(1);
		if (i == 4) {PLOTS[i].sigmaPlot->SetMarkerColor(i+2);}
		PLOTS[i].sigmaPlot->SetMarkerStyle(21);
		PLOTS[i].sigmaPlot->SetLineWidth(1);
		
		TF1 *calibrationFit = new TF1("calibrationFit", "pol1", 0, analyzer->getOverflowPos());
		PLOTS[i].calPlot->Fit("calibrationFit", "", "", 0, 0);
		
		RESULTS[i].calPars.slope = calibrationFit->GetParameter(1);
		RESULTS[i].calPars.slopeErr = calibrationFit->GetParError(1);
		RESULTS[i].calPars.offset = calibrationFit->GetParameter(0);
		RESULTS[i].calPars.offsetErr = calibrationFit->GetParError(0);
		//RESULTS[i].calibratedData = (data[i] - RESULTS[i].calPars.offset) 
		//                             / RESULTS[i].calPars.slope;

		string calName = "calibrated" + to_string(i);
		string label = getLabels(mode, i);
		Int_t maxCalBin = (analyzer->getOverflowPos() - RESULTS[i].calPars.offset) / 
		                   RESULTS[i].calPars.slope;
		PLOTS[i].calibrated = new TH1D(calName.c_str(), label.c_str(), 20e3, 0, maxCalBin);
		PLOTS[i].calibrated->SetLineColor(i+1);
		if (i == 4) {PLOTS[i].calibrated->SetLineColor(i+2);}
		PLOTS[i].calibrated->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		PLOTS[i].calibrated->GetYaxis()->SetTitle("Count");
		
		string calibration = "(energy - " + to_string(RESULTS[i].calPars.offset);
		calibration += ") / " + to_string(RESULTS[i].calPars.slope);
		calibration += " >> " + calName;
		data[i]->Draw(calibration.c_str(), "channel==0", "goff");

		Int_t rateCount;
		Int_t rateBin = PLOTS[i].calibrated->FindBin(50);
		while (rateBin > 0) {
			rateCount += PLOTS[i].calibrated->GetBinContent(rateBin);
			rateCount--;
		}
		RESULTS[i].lowERate = (Double_t) rateCount / (Double_t) time;
		cout << "Rate below 50 keV: " << RESULTS[i].lowERate << endl;
		
		Int_t noiseCount;
		Int_t noiseBin = PLOTS[i].calibrated->FindBin(50);
		Int_t noiseMin = 2 * PLOTS[i].calibrated->GetBinContent(noiseBin);
		while (noiseCount < 2 * noiseMin) {
			noiseCount = PLOTS[i].calibrated->GetBinContent(noiseBin);
			if (noiseCount < noiseMin) {noiseMin = noiseCount;}
			noiseBin--;
			if (noiseBin == 0) {break;}
		}
		RESULTS[i].noiseWall = PLOTS[i].calibrated->GetXaxis()->GetBinCenter(noiseBin);
		cout << "Noise Wall at: " << RESULTS[i].noiseWall << " keV" << endl;
		
		cout << endl;
		cout << "-------------------------------------------------------------" << endl;
		cout << endl;

	}

/*
   #####################
   #   MODE HANDLING   #
   #####################
*/

/*	if (mode == "fir") {

		TCanvas *FIRCanvas = new TCanvas("FIR", "FIR", 1);
		Double_t peakingTimes[NUMFILES];
		Double_t FIRResolution[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {
			peakingTimes[k] = PEAKINGTIMES[k];
			//peakingTimes[k] = data[i]->GetLeaf("peaktime")->GetValue();
			FIRResolution[k] = pow(PEAKS[k].K.sigma / PEAKS[k].K.mu, 2.0);
		}
		TGraph *FIROptPlot = new TGraph(5, peakingTimes, FIRResolution);
		string FIRFitFunc = "([0]/x)^2 + [1]^2 + ([3]*x)^2";
		TF1 *FIRFit = new TF1("FIRFit", FIRFitFunc.c_str(), 100.0, 300.0);
		//FIRFit->SetParameters(
		//FIROptPlot->Fit("FIRFit", "", "", 0, 0);
		//Double_t OptPeakingTime = FIRFit->GetMinimum();
		
		FIROptPlot->Draw();
		//cout << "OPTIMAL PEAKING TIME: " << OptPeakingTime << endl;

	} else */
	if (mode == "pos") {

		Double_t resolution = PEAKS[3].get(1460.0).sigma / PEAKS[3].get(1460.0).mu;
		cout << "Normalized detector resolution (width of Cs peak" << endl;
		cout << " / mean at 3rd position):" << resolution << endl;
		
		Double_t CsPos[NUMFILES];
		Double_t CsPosErrs[NUMFILES];
		Double_t CsSig[NUMFILES];
		Double_t CsSigErrs[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {
			CsPos[k] = PEAKS[k].get(1460.0).mu;
			CsPosErrs[k] = PEAKS[k].get(1460.0).muErr;
			CsSig[k] = PEAKS[k].get(1460.0).sigma / PEAKS[k].get(1460.0).mu;
			CsSigErrs[k] = PEAKS[k].get(1460.0).sigmaErr / PEAKS[k].get(1460.0).mu;
		}

		Double_t xAxis[NUMFILES];
		for (Int_t i = 0; i < 5; i++) {
			if (mode == "pos") {
				xAxis[i] = POSITIONS[i];
			} else if (mode == "volt") {
				xAxis[i] = VOLTAGES[i];
			} else if (mode == "fir") {
				//xAxis[i] = PEAKINGTIMES[i];
			}
		}
		
		TCanvas *CsPosCanvas = new TCanvas("CsPosCanvas", "CsPosCanvas");
		TGraphErrors *CsPosGraph = new TGraphErrors(5, xAxis, CsPos, 0, CsPosErrs);
		CsPosGraph->SetTitle(("Cs Peak Energy vs " + getDepVar(mode)).c_str());
		CsPosGraph->GetXaxis()->SetTitle(getDepVar(mode).c_str());
		CsPosGraph->GetYaxis()->SetTitle("Uncalibrated Cs Peak Energy");
	
		formatGraph(CsPosGraph);
		/*CsPosGraph->SetMarkerColor(4);
		CsPosGraph->SetMarkerStyle(21);
		CsPosGraph->SetLineColor(1);
		CsPosGraph->SetLineWidth(2);
		CsPosGraph->GetYaxis()->SetTitleOffset(1.5);
		CsPosGraph->GetXaxis()->SetTitleOffset(1.2);*/
		CsPosGraph->Draw();
		CsPosCanvas->Print("CsEvsPos.pdf[");
		CsPosCanvas->Print("CsEvsPos.pdf");
		CsPosCanvas->Print("CsEvsPos.pdf]");
						
		TCanvas *CsResCanvas = new TCanvas("CsResCanvas", "CsResCanvas");
		TGraphErrors *CsResGraph = new TGraphErrors(5, xAxis, CsSig, 0, CsSigErrs);
		CsResGraph->SetTitle(("Cs Peak Width vs " + getDepVar(mode)).c_str());
		CsResGraph->GetXaxis()->SetTitle(getDepVar(mode).c_str());
		CsResGraph->GetYaxis()->SetTitle("Normalized Cs Peak Width (Cs sigma / Cs mean)");
		
		formatGraph(CsResGraph);
		/*CsResGraph->SetMarkerColor(4);
		CsResGraph->SetMarkerStyle(21);
		CsResGraph->SetLineColor(1);
		CsResGraph->SetLineWidth(2);
		CsResGraph->GetYaxis()->SetTitleOffset(1.5);
		CsResGraph->GetXaxis()->SetTitleOffset(1.2);*/
		CsResGraph->Draw();
		CsResCanvas->Print("CsResvsPos.pdf[");
		CsResCanvas->Print("CsResvsPos.pdf");
		CsResCanvas->Print("CsResvsPos.pdf]");
		
	}

/*
   #######################
   #   OPTION HANDLING   #
   #######################
*/

	if (option == "cal") {	
		
		TCanvas *calCanvas = new TCanvas("Calibration Canvas", "Calibration Canvas", 1);
		TMultiGraph *calComp = new TMultiGraph("calComp", "calComp");
		for (Int_t k = 0; k < NUMFILES; k++) {
			calComp->Add((TGraph*) PLOTS[k].calPlot->Clone());
			// Y-axis is ADC energy, X is real energy
		}
		calComp->Draw("ALP");
		calComp->SetTitle(("Calibration Curves for Each " + getDepVar(mode)).c_str());
		calComp->GetYaxis()->SetTitle("ADC Energy");
		calComp->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		calCanvas->BuildLegend(0.15,0.6,0.30,0.85); 	// legend in top left
		
	} else if (option == "sig") {
		
		TCanvas *sigmaCanvas = new TCanvas("Sigma Canvas", "Sigma Canvas", 1);
		TMultiGraph *sigmaComp = new TMultiGraph("SigmaComp", "SigmaComp");
		for (Int_t k = 0; k < NUMFILES; k++) {
			sigmaComp->Add((TGraphErrors*) PLOTS[k].sigmaPlot->Clone());	
		}
		sigmaComp->Draw("ALP");
		sigmaComp->SetTitle(("Resolution vs " + getDepVar(mode)).c_str());
		sigmaComp->GetYaxis()->SetTitle("Fractional Sigma (sigma/mean)");
		sigmaComp->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		sigmaCanvas->BuildLegend(0.7,0.6,0.85,0.85); 	// legend in top right

	} else if (option == "res") {

		TCanvas *resCanvas = new TCanvas("Residue Canvas", "Residue Canvas", 1);
		TMultiGraph *resComp = new TMultiGraph("ResComp", "ResComp");
		Double_t expE[3] = {2615.0, 1460.0, 661.6};
		for (Int_t k = 0; k < NUMFILES; k++) {
			Double_t residues[3];
			
			Double_t energy = (PEAKS[k].get(2615.0).mu - RESULTS[k].calPars.offset)
			                   / RESULTS[k].calPars.slope;
			residues[0] = 100 * (energy - expE[0]) / expE[0];
				
			energy = (PEAKS[k].get(1460.0).mu - RESULTS[k].calPars.offset) / 
			          RESULTS[k].calPars.slope;
			residues[1] = 100 * (energy - expE[1]) / expE[1];			
			
			energy = (PEAKS[k].get(661.6).mu - RESULTS[k].calPars.offset) / 
			          RESULTS[k].calPars.slope;
			residues[2] = 100 * (energy - expE[2]) / expE[2];
			
			TGraphErrors *residuePlot = new TGraphErrors(3, expE, residues, 0, 0);
			residuePlot->SetTitle(getLabels(mode, k).c_str());
			residuePlot->SetLineColor(1);
			residuePlot->SetMarkerColor(k+1);
			if (k == 4) {residuePlot->SetMarkerColor(k+2);}			
			residuePlot->SetMarkerStyle(21);
			residuePlot->SetLineWidth(1);
			resComp->Add((TGraphErrors*) residuePlot->Clone());
		}
		resComp->Draw("ALP");
		resComp->SetTitle(("Residues for " + getDepVar(mode) + " Variation").c_str());
		resComp->GetXaxis()->SetTitle("ADC Energies");
		resComp->GetXaxis()->SetTitleOffset(1.3);
		resComp->GetYaxis()->SetTitle("\% Error in calibrated energy ((Fit - Exp)/Exp)");
		resComp->GetYaxis()->SetTitleOffset(1.2);
		resCanvas->BuildLegend(0.15,0.15,0.30,0.40); 	// legend in bottom left

	} else if (option == "gain") {

		Double_t gain[NUMFILES];
		Double_t gainErr[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {
			gain[k] = TMath::Log(RESULTS[k].calPars.slope);
		}
		for (Int_t k = 0; k < NUMFILES; k++) {
			gainErr[k] = RESULTS[k].calPars.slopeErr / RESULTS[k].calPars.slope;
		}	
		
		Double_t xAxis[NUMFILES];
		for (Int_t i = 0; i < 5; i++) {
			if (mode == "pos") {
				xAxis[i] = POSITIONS[i];
			} else if (mode == "volt") {
				xAxis[i] = VOLTAGES[i];
			} else if (mode == "fir") {
				//xAxis[i] = PEAKINGTIMES[i];
			}
		}
		
		TCanvas *gainCanvas = new TCanvas("Gain Canvas", "Gain Canvas", 1);
		
		TGraphErrors *gainGraph = new TGraphErrors(5, xAxis, gain, 0, gainErr);
		gainGraph->SetTitle(("Detector Gain vs " + getDepVar(mode)).c_str());
		gainGraph->GetXaxis()->SetTitle(getDepVar(mode).c_str());
		gainGraph->GetYaxis()->SetTitle("Log(Calibration Slope)");
		formatGraph(gainGraph);

		/*GainGraph->SetMarkerColor(4);
		GainGraph->SetMarkerStyle(21);
		GainGraph->SetLineColor(1);
		GainGraph->SetLineWidth(2);*/
		
		TF1 *gainFit = new TF1("gainFit", "pol2");
		gainFit->SetParNames("Log(G0)", "Slope", "Curvature");
		gainGraph->Fit("gainFit");
		gainGraph->Draw();
		
		if (mode == "volt") {
			gainCanvas->Print("Gain.pdf[");
			gainCanvas->Print("Gain.pdf");
			gainCanvas->Print("Gain.pdf]");
		}

	} else if (option == "over") {
				
		TCanvas *overlayCanvas = new TCanvas("overlayCanvas", "Overlay Canvas", 1);
		gPad->SetLogy();
		for (Int_t k = 0; k < NUMFILES; k++) {
			PLOTS[k].calibrated->Draw("SAME");
		}
		overlayCanvas->BuildLegend(0.7,0.6,0.85,0.85); 	// legend in top right

	} else if (option == "noise") {
		
		Double_t noiseWalls[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {
			noiseWalls[k] = RESULTS[k].noiseWall;
		}
		
		Double_t xAxis[NUMFILES];
		for (Int_t i = 0; i < 5; i++) {
			if (mode == "pos") {
				xAxis[i] = POSITIONS[i];
			} else if (mode == "volt") {
				xAxis[i] = VOLTAGES[i];
			} else if (mode == "fir") {
				//xAxis[i] = PEAKINGTIMES[i];
			}
		}
		
		TCanvas *noiseCanvas = new TCanvas("noiseCanvas", "NoiseCanvas", 1);
		TGraphErrors *noiseGraph = new TGraphErrors(5, xAxis, noiseWalls, 0, 0);
		noiseGraph->SetTitle(("Noise Wall Energy vs " + getDepVar(mode)).c_str());
		noiseGraph->GetXaxis()->SetTitle(getDepVar(mode).c_str());
		noiseGraph->GetYaxis()->SetTitle("Noise Wall Energy (keV)");

		/*NoiseGraph->SetMarkerColor(4);
		NoiseGraph->SetMarkerStyle(21);
		NoiseGraph->SetLineColor(1);
		NoiseGraph->SetLineWidth(2);
		NoiseGraph->GetYaxis()->SetTitleOffset(1.5);
		NoiseGraph->GetXaxis()->SetTitleOffset(1.2);	*/	

		formatGraph(noiseGraph);
		noiseGraph->Draw();
		if (mode == "volt") {
			noiseCanvas->Print("Noise.pdf[");
			noiseCanvas->Print("Noise.pdf");
			noiseCanvas->Print("Noise.pdf]");
		}
		
	} else if (option == "backFits") {
		/*
		TCanvas *backCanvas = new TCanvas("backCanvas", "Background Fits", 1400, 900);
		backCanvas->Divide(3, 5);
		for (Int_t i = 0; i < 5; i++) {
			backCanvas->cd(3*i+1);
			PLOTS[i].backPlots.Tl->Draw();
			backCanvas->cd(3*i+2);
			PLOTS[i].backPlots.K->Draw();
			backCanvas->cd(3*i+3);
			PLOTS[i].backPlots.Cs->Draw();
		}*/

	} else if (option == "AE") {
		
		TCanvas *AECanvas = new TCanvas("AECanvas", "A/E Plot", 1);
		// Figure out how to do mathematical manipulations on a TChain, then I can sum
		// all our data.

		/*TChain summed = new TChain("st");
		for (Int_t i = 0; i < NUMFILES; i++) {
			summed.Add((data[i] - RESULTS[i].calPars.offset) / );
		}

		Example(){

		  TFile *_file0 = TFile::Open("tree1.root");

		  Float_t px;
		  t1->SetBranchAddress("px",&px);
		  for(int i=0; i<t1->GetEntries(); i++){
		    t1->GetEntry(i);
		    cout<<"px: "<<px<<endl;
		  }

		}*/

		TH2D *AEHist = new TH2D("AEHist", "Amplitude / Energy vs Energy (calibrated)",
		                        1e3, 0, 3e3, 1e3, 0, 10);

		string calE = "(energy - " + to_string(RESULTS[NUMFILES - 1].calPars.offset);
		calE += ") / " + to_string(RESULTS[NUMFILES - 1].calPars.slope);
		string toPlot = "amp / (" + calE + ") : (" + calE + ") >> AEHist";
		data[NUMFILES - 1]->Draw(toPlot.c_str(), "channel==0", "COLZ");
		
		AEHist->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		AEHist->GetYaxis()->SetTitle("Amplitude / Callibrated Energy");
		
		

	} else if (option == "muon") {
		// muons at about 10x Tl peak position.
		// if it fails, get results and set to unphysical estimates;

		/*
		TCanvas *muonCanvas = new TCanvas("muonCanvas", "Muon Canvas", 1);
		TH1D *muonHist = new TH1D("muonHist", "Muon Peak", "

		data[0]->
		*/
	}

}

void getData(TChain *data[NUMFILES], string mode) {
	cout << "Getting Data..." << endl;
	
	string path[NUMFILES];
	if (mode == "pos") {
		cout << "Finding Position Data..." << endl;
		path[0] = "position/position_1/NaI_ET_run*";
		path[1] = "position/position_2/NaI_ET_run*";
		path[2] = "position/position_3/NaI_ET_run*";
		path[3] = "position/position_4/NaI_ET_run*";
		path[4] = "position/position_5/NaI_ET_run*";
	} else if (mode == "volt") {
		cout << "Finding Voltage Data..." << endl;
		path[0] = "voltage/600_V/NaI_ET_run*";
		path[1] = "voltage/700_V/NaI_ET_run*";
		path[2] = "voltage/800_V/NaI_ET_run*";
		path[3] = "voltage/900_V/NaI_ET_run*";
		path[4] = "voltage/1000_V/NaI_ET_run*";
	} else if (mode == "fir") {
		cout << "Finding FIR Data..." << endl;
		path[0] = "fir/fir_1_lowest/NaI_ET_run*";
		path[1] = "fir/fir_2/NaI_ET_run*";
		path[2] = "fir/fir_3/NaI_ET_run*";
		path[3] = "fir/fir_4/NaI_ET_run*";
		path[4] = "fir/fir_5_highest/NaI_ET_run*";
	}
	
	for (Int_t i = 0; i < NUMFILES; i++) {
		data[i] = new TChain("st");
		data[i]->Add(path[i].c_str());
	}
	cout << "Data successfully loaded" << endl;
}

Double_t getRunTime(TChain *data, string mode) {
	Double_t time;
	if (mode == "fir") {
		time = (Double_t) data->GetNtrees() * 300;
	} else {
		time = (Double_t) data->GetNtrees() * 600;
	}
	return time;
}
/*
FitResults backEst(TH1D *h, FitWindow win, Double_t coverage, string func, string peak, Int_t i) {
	cout << "Estimating Background..." << endl;

	Int_t lowBin = h->FindBin(win.low);
	Int_t highBin = h->FindBin(win.high);
	Int_t overallBinRange = highBin - lowBin;
	Int_t backWindowRange = (Int_t) (coverage / 2.0 * overallBinRange);

	Float_t backVals[2 * backWindowRange];
	Float_t backPos[2 * backWindowRange];

	for (Int_t z = 0; z < backWindowRange; z++) {
		backVals[z] = h->GetBinContent(z + lowBin);
		backPos[z] = h->GetBinCenter(z + lowBin);
	}
	for (Int_t z = 0; z < backWindowRange; z++) {
		backVals[z + backWindowRange] = h->GetBinContent((highBin - backWindowRange) + (z + 1));
		backPos[z + backWindowRange] = h->GetBinCenter((highBin - backWindowRange) + (z + 1));
	}

	TGraph *backGraph = new TGraph(2 * backWindowRange, backPos, backVals);
	TF1 backFit = TF1("backFit", func.c_str(), win.low, win.high);
	backGraph->Fit("backFit");

	FitResults pars;
	pars.offset = backFit.GetParameter(0);
	pars.slope = backFit.GetParameter(1);
	
	string title = "Background Estimation " + to_string(i) + " for " + peak + " peak";
	backGraph->SetTitle((title).c_str());
	backGraph->GetYaxis()->SetTitle("Count");
	backGraph->GetXaxis()->SetTitle("Uncalibrated Energy");
	backGraph->SetMarkerStyle(4);
	backGraph->SetMarkerSize(0.5);

	if (peak == "Cs") {
		PLOTS[i].backPlots.Cs = backGraph;
	} else if (peak == "K") {
		PLOTS[i].backPlots.K = backGraph;
	} else if (peak == "Tl") {
		PLOTS[i].backPlots.Tl = backGraph;
	}

	return pars;
}*/

string getDepVar(string mode) {
	if (mode == "pos") {
		return "Position";
	} else if (mode == "volt") {
		return "Voltage";
	} else if (mode == "fir") {
		return "Peaking Time";
	}
	return "error in getDepVar: bad mode";
}

string getLabels(string mode, Int_t index) {
	if (mode == "pos") {
		return "Position " + to_string(POSITIONS[index]);
	} else if (mode == "volt") {
		return to_string(VOLTAGES[index]) + " V";
	} else if (mode == "fir") {
		//return "Peaking Time: " + to_string(PEAKINGTIMES[index]);
	}
	return "error in getLabels: bad mode";
}

void formatGraph(TGraphErrors *graph) {
	graph->SetMarkerColor(4);
	graph->SetMarkerStyle(21);
	graph->SetLineColor(1);
	graph->SetLineWidth(2);
	graph->GetYaxis()->SetTitleOffset(1.4);
	graph->GetXaxis()->SetTitleOffset(1.2);
}
