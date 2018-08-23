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

TApplication* app = new TRint("app", 0, NULL);
vector<PeakFinder*> ANALYZERS;

Int_t Calibration(string path, string mode, string option) {

	//gStyle->SetOptFit(1111);
	cout << "Collecting Data..." << endl;
	vector<string> filepaths;
	vector<TChain*> DATA;
	vector<Int_t> POSITIONS;
	vector<Int_t> VOLTAGES;
	if (mode == "pos") {
		vector<Int_t> testedPositions = {1, 2, 3, 4, 5}; 
		// will segfault if there's no data for any one of these positions
		cout << "Finding position data..." << endl;
		cout << "Calibrating positions ";
		for (Int_t pos : testedPositions) {
			cout << pos << " ";
			string expectedPath = path + "/position/position_" + to_string(pos);
			expectedPath += "/NaI_ET_run*";
			filepaths.push_back(expectedPath);
		}
		POSITIONS = testedPositions;
		cout << endl;
	} else if (mode == "volt") {
		vector<Int_t> testedVoltages = {600, 700, 800, 900, 1000};
		// will segfault if there's no data for any one of these voltages
		cout << "Finding voltage data..." << endl;
		cout << "Calibrating Votlages ";
		for (Int_t volt : testedVoltages) {
			cout << volt << " ";
			string expectedPath = path + "/voltage/" + to_string(volt);
			expectedPath += "_V/NaI_ET_run*";
			filepaths.push_back(expectedPath);
		}
		VOLTAGES = testedVoltages;
		cout << endl;
	}
	for (Int_t i = 0; i < filepaths.size(); i++) {
		DATA.push_back(new TChain("st"));
		DATA[i]->Add(filepaths[i].c_str());
	}
	Int_t NUMFILES = DATA.size();

/*
#######################
#   USER PARAMETERS   #
#######################

*/

	using ParGuess = pair<Int_t, Double_t>;
	using ParLimit = pair<Int_t, ParWindow>;

	vector<FitInfo> peakPars;
	string CHANNEL = "channel==4"; // digitizer channel to use
	
	// 208Tl peak parameters:
	FitInfo TlPars;
	TlPars.peakEnergies.push_back(2614.511);
	TlPars.fitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2) + exp([3] + [4]*x)";

	TlPars.fitPars.insert(ParGuess (0, 1.0));
	TlPars.fitPars.insert(ParGuess (1, 1.0));
	TlPars.fitPars.insert(ParGuess (2, 0.05));
	// pars [3] and [4] will be set by background estimation in PeakFinder
	
	TlPars.fitWindow.low = 0.9; // window limits = proportions of highest energy peak mu
	TlPars.fitWindow.high = 1.1; // i.e. 1.1 times peak mu
	TlPars.backgroundRange = 0.3; // proportion of window to be considered in background est.
	peakPars.push_back(TlPars);
	
	// 40K peak parameters:
	FitInfo KPars;
	KPars.peakEnergies.push_back(1460.820);
	KPars.fitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2) + exp([3] + [4]*x)";

	KPars.fitPars.insert(ParGuess (0, 1.0));
	KPars.fitPars.insert(ParGuess (1, 1.0));
	KPars.fitPars.insert(ParGuess (2, 0.05));
	
	KPars.fitWindow.low = 0.85;
	KPars.fitWindow.high = 1.15;
	KPars.backgroundRange = 0.3;
	peakPars.push_back(KPars);

	// 137Cs peak parameters:
	FitInfo CsPars;
	CsPars.peakEnergies.push_back(661.657);
	CsPars.peakEnergies.push_back(583.187);
	CsPars.fitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[2])^2)";
	CsPars.fitFunc += " + exp([5]+[6]*x)";
	
	CsPars.fitPars.insert(ParGuess (0, 1.0));
	CsPars.fitPars.insert(ParGuess (1, 1.0));
	CsPars.fitPars.insert(ParGuess (2, 0.05));
	CsPars.fitPars.insert(ParGuess (3, 0.1));
	CsPars.fitPars.insert(ParGuess (4, 583.187 / 661.657));

	/*ParWindow Cs661ParWindow;
	Cs661ParWindow.low = 0.97;
	Cs661ParWindow.high = 1.03;
	CsPars.fitParLimits.insert(ParLimit (1, Cs661ParWindow));*/

	ParWindow Cs583ParWindow;
	Cs583ParWindow.low = 583.187 / 661.657 - 0.03;
	Cs583ParWindow.high = 583.187 / 661.657 + 0.03;
	CsPars.fitParLimits.insert(ParLimit (4, Cs583ParWindow));
	
	CsPars.fitWindow.low = 0.75;
	CsPars.fitWindow.high = 1.25;
	CsPars.backgroundRange = 0.3;
	peakPars.push_back(CsPars);

	TCanvas *fitCanvas = new TCanvas("fitCanvas", "fitCanvas");
	fitCanvas->Divide(NUMFILES / 2, NUMFILES / 2 + 1);

	for (Int_t i = 0; i < NUMFILES; i++) {
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << endl;

		cout << "Beginning Calibration " << i + 1 << "..." << endl;

		Double_t time = DATA[i]->GetNtrees() * 600;
		cout << "Run time in data chain: " << time << " seconds" << endl;

		// Need to initialize and populate a PeakSet to bind to the analyzer engine.  
		// Simultaneously generates vector of all peak energies to be considered.
		PeakSet allPeaks;
		vector<Double_t> allEnergies;
		for (Int_t j = 0; j < peakPars.size(); j++) {
			vector<Double_t> thisPeakEnergies = peakPars[j].peakEnergies;
			for (Int_t k = 0; k < thisPeakEnergies.size(); k++) {
				allEnergies.push_back(thisPeakEnergies[k]);
				PeakInfo curr;
				curr.energy = thisPeakEnergies[k];
				allPeaks.put(curr);
			}
		}
		PeakFinder *analyzer = new PeakFinder(fitCanvas, DATA[i], allPeaks, time, CHANNEL);
		
		for (Int_t j = 0; j < peakPars.size(); j++) {

			FitInfo pars;
			pars.peakEnergies = peakPars[j].peakEnergies;
			pars.fitWindow = peakPars[j].fitWindow;
			pars.fitFunc = peakPars[j].fitFunc;
			pars.fitPars = peakPars[j].fitPars;
			pars.fitParLimits = peakPars[j].fitParLimits;
			pars.backgroundRange = peakPars[j].backgroundRange;

			Double_t energy = pars.peakEnergies[0];
			analyzer->findPeak(energy);
			PeakInfo estimate = analyzer->getPeakSet().get(energy);

			// Rescaling parameter guesses and limits:
			for (Int_t k = 0; k < pars.fitPars.size(); k++) {
				if (k == 0) {
					pars.fitPars[k] = pars.fitPars[0] * estimate.count;
				} else if (k == 1) {
					pars.fitPars[k] = pars.fitPars[1] * estimate.mu;
				} else if (k == 2) {
					pars.fitPars[k] = pars.fitPars[2] * estimate.mu;
				} else if (k % 2 == 0) {
					pars.fitPars[k] = pars.fitPars[k] * estimate.mu;
				} else {
					pars.fitPars[k] = pars.fitPars[k] * estimate.count;
				}
			}
			for (pair<Int_t, ParWindow> lims : pars.fitParLimits) {
				ParWindow scaledByMu;
				scaledByMu.low = lims.second.low * estimate.mu;
				scaledByMu.high = lims.second.high * estimate.mu;

				ParWindow scaledByCount;
				scaledByCount.low = lims.second.low * estimate.count;
				scaledByCount.high = lims.second.high * estimate.count;

				if (lims.first == 0) {
					pars.fitParLimits[lims.first] = scaledByCount;
				} else if (lims.first == 1) {
					pars.fitParLimits[lims.first] = scaledByMu;
				} else if (lims.first == 2) {
					pars.fitParLimits[lims.first] = scaledByMu;
				} else if (lims.first % 2 == 0) {
					pars.fitParLimits[lims.first] = scaledByMu;
				} else {
					pars.fitParLimits[lims.first] = scaledByCount;
				}
			}

			for (Double_t E : pars.peakEnergies) {
				if (E == 583.187) {
					PeakInfo TlPeak = analyzer->getPeakSet().get(2614.511);
					PeakInfo KPeak = analyzer->getPeakSet().get(1460.820);

					Double_t rise = TlPeak.mu - KPeak.mu;
					Double_t run = 2614.511 - 1460.820;
					Double_t slope = rise / run;
					Double_t offset = TlPeak.mu - 2614.511 * slope;

					pars.fitPars[4] = slope * 583.187 + offset;
					
					ParWindow parLimits;
					parLimits.low = pars.fitPars[4];
					parLimits.high = pars.fitPars[4];
					pars.fitParLimits[4] = parLimits;
				}
			}

			// Rescaling fit window:
			pars.fitWindow.low = pars.fitWindow.low * estimate.mu;
			pars.fitWindow.high = pars.fitWindow.high * estimate.mu;
			
			analyzer->fit(pars);
		}
		// save fits as .root / .png
		Float_t range = 1.15 * analyzer->getPeakSet().getHighestEnergyPeak().mu;
		analyzer->getRawPlot()->GetXaxis()->SetRangeUser(0, range);
		analyzer->getRawPlot()->GetXaxis()->SetTitle("Uncalibrated Energy");
		analyzer->getRawPlot()->GetYaxis()->SetTitle("Counts");
		
		string rawPlotTitle = "Fits for ";
		if (mode == "pos") {
			rawPlotTitle += "Position " + to_string(POSITIONS[i]);
		} else if (mode == "volt") {
			rawPlotTitle += to_string(VOLTAGES[i]) + " V";
		}
		analyzer->getRawPlot()->SetTitle(rawPlotTitle.c_str());
		
		fitCanvas->cd(i + 1);
		gPad->SetLogy();
		analyzer->getRawPlot()->Draw();

		FitResults calPars = analyzer->findCalibration();		

		ANALYZERS.push_back(analyzer);
		
		cout << endl;
		cout << "-------------------------------------------------------------" << endl;
		cout << endl;

	}

	// At this point, a calibration for the data has been found.  What goes below is just 
	// more detailed analysis, as defined by the options passed into the calibration script.

/*
   #####################
   #   MODE HANDLING   #
   #####################
*/

	if (mode == "pos") {

		//Double_t CsEnergy = 1460.820;
		Double_t CsEnergy = 661.657;
		
		vector<Double_t> calibratedCsEnergies;
		vector<Double_t> calibratedCsEnergyErrs;
		vector<Double_t> calibratedCsSigmas;
		vector<Double_t> calibratedCsSigmaErrs;
		vector<Double_t> positions;
		for (Int_t i = 0; i < NUMFILES; i++) {
			PeakInfo CsPeak = ANALYZERS[i]->getPeakSet().get(CsEnergy);
			FitResults calib = ANALYZERS[i]->getCalibration();

			Measurement CsPos;
			CsPos.val = CsPeak.mu;
			CsPos.err = CsPeak.muErr;
			
			Measurement CsSigma;
			CsSigma.val = CsPeak.sigma;
			CsSigma.err = CsPeak.sigmaErr;

			Measurement calCsEnergy = ANALYZERS[i]->calibrate(CsPos);
			
			// calibrating resolution
			Double_t term1 = TMath::Power(CsSigma.err / CsSigma.val, 2);
			Double_t term2 = TMath::Power(calib.slopeErr / calib.slope, 2);

			Measurement calCsSigma;
			calCsSigma.val = CsSigma.val / calib.slope;
			calCsSigma.err = calCsSigma.val * TMath::Sqrt(term1 + term2);

			calibratedCsEnergies.push_back(calCsEnergy.val);
			calibratedCsEnergyErrs.push_back(calCsEnergy.err);
			calibratedCsSigmas.push_back(calCsSigma.val);
			calibratedCsSigmaErrs.push_back(calCsSigma.err);

			positions.push_back((Double_t) POSITIONS[i]);		
		}

		cout << "############################################" << endl;
		cout << "Cs PEAK PARS VS POSITION: " << endl;
		Double_t resolution = calibratedCsSigmas[3];
		Double_t resolutionErr = calibratedCsSigmaErrs[3];
		cout << "Cs peak resolution at 3rd position (keV): " << resolution;
		cout << " +/- " << resolutionErr << endl;
		cout << "############################################" << endl;
		
		TCanvas *CsPosCanvas = new TCanvas("CsPosCanvas", "CsPosCanvas");
		TGraphErrors *CsPosGraph = new TGraphErrors(NUMFILES, &positions[0], 
		                                            &calibratedCsEnergies[0], 0, 
		                                            &calibratedCsEnergyErrs[0]);
		string titleVar;
		if (mode == "pos") {
			titleVar = "Position";
		} else if (mode == "volt") {
			titleVar = "Voltage";
		}

		CsPosGraph->SetTitle(("Cs Peak Energy vs " + titleVar).c_str());
		CsPosGraph->GetXaxis()->SetTitle(titleVar.c_str());
		CsPosGraph->GetYaxis()->SetTitle("Calibrated Cs Peak Energy (keV)");
		CsPosGraph->SetMarkerColor(4);
		CsPosGraph->SetMarkerStyle(21);
		CsPosGraph->SetLineColor(1);
		CsPosGraph->SetLineWidth(2);
		CsPosGraph->GetYaxis()->SetTitleOffset(1.4);
		CsPosGraph->GetXaxis()->SetTitleOffset(1.2);
	
		CsPosGraph->Draw();
		CsPosCanvas->Print("CsEvsPos.pdf[");
		CsPosCanvas->Print("CsEvsPos.pdf");
		CsPosCanvas->Print("CsEvsPos.pdf]");
						
		TCanvas *CsResCanvas = new TCanvas("CsResCanvas", "CsResCanvas");
		TGraphErrors *CsResGraph = new TGraphErrors(NUMFILES, &positions[0], 
		                                            &calibratedCsSigmas[0], 0, 
		                                            &calibratedCsSigmaErrs[0]);

		CsResGraph->SetTitle(("Cs Peak Resolution vs " + titleVar).c_str());
		CsResGraph->GetXaxis()->SetTitle(titleVar.c_str());
		CsResGraph->GetYaxis()->SetTitle("Width of Cs Peak (keV)");
		CsResGraph->SetMarkerColor(4);
		CsResGraph->SetMarkerStyle(21);
		CsResGraph->SetLineColor(1);
		CsResGraph->SetLineWidth(2);
		CsResGraph->GetYaxis()->SetTitleOffset(1.4);
		CsResGraph->GetXaxis()->SetTitleOffset(1.2);
		
		CsResGraph->Draw();
		CsResCanvas->Print("CsResvsPos.pdf[");
		CsResCanvas->Print("CsResvsPos.pdf");
		CsResCanvas->Print("CsResvsPos.pdf]");

	} else if (mode == "volt") {

		vector<Double_t> gains;
		vector<Double_t> gainErrs;
		vector<Double_t> voltages;
		for (Int_t i = 0; i < NUMFILES; i++) {
			FitResults calib = ANALYZERS[i]->getCalibration();
			Double_t gain = TMath::Log(calib.slope);
			Double_t gainErr = calib.slopeErr / calib.slope;

			gains.push_back(gain);
			gainErrs.push_back(gainErr);
		
			voltages.push_back((Double_t) VOLTAGES[i]);
		}
		
		TCanvas* gainCanvas = new TCanvas("Gain Canvas", "Gain Canvas", 1);
		
		TGraphErrors* gainGraph = new TGraphErrors(NUMFILES, &voltages[0], 
		                                           &gains[0], 0, &gainErrs[0]);
		string titleVar;
		if (mode == "pos") {
			titleVar = "Position";
		} else if (mode == "volt") {
			titleVar = "Voltage";
		}

		gainGraph->SetTitle(("Detector Gain vs " + titleVar).c_str());
		gainGraph->GetXaxis()->SetTitle((titleVar + " (V)").c_str());
		gainGraph->GetYaxis()->SetTitle("Log(Calibration Slope)");
		gainGraph->SetMarkerColor(4);
		gainGraph->SetMarkerStyle(21);
		gainGraph->SetLineColor(1);
		gainGraph->SetLineWidth(2);
		gainGraph->GetYaxis()->SetTitleOffset(1.4);
		gainGraph->GetXaxis()->SetTitleOffset(1.2);
		
		TF1 *gainFit = new TF1("gainFit", "pol2");
		gainFit->SetParNames("Log(G0)", "Slope", "Curvature");
		gainGraph->Fit("gainFit");
		gainGraph->Draw();

		// test for valid barium/muon fits here

		cout << endl;
		cout << "############################################" << endl;
		cout << "GAIN VS VOLTAGE PARAMETERS (pol2):" << endl;

		Double_t gainOffset = gainFit->GetParameter(0);
		Double_t gainOffsetErr = gainFit->GetParError(0);
		Double_t gainSlope = gainFit->GetParameter(1);
		Double_t gainSlopeErr = gainFit->GetParError(1);
		Double_t gainCurvature = gainFit->GetParameter(2);
		Double_t gainCurvatureErr = gainFit->GetParError(2);

		cout << "Gain Offset (LOG(G0)): " << gainOffset << " +/- " << gainOffsetErr << endl;
		cout << "Gain Slope: " << gainSlope << " +/- " << gainSlopeErr << endl;
		cout << "Gain Curvature: " << gainCurvature << " +/- " << gainCurvatureErr << endl;
		cout << "############################################" << endl;
		
		gainCanvas->Print("GainVsVolt.pdf[");
		gainCanvas->Print("GainVsVolt.pdf");
		gainCanvas->Print("GainVsVolt.pdf]");

	}

/*
   #######################
   #   OPTION HANDLING   #
   #######################
*/

	if (option.find("ba") != string::npos) {
		
		for (Int_t i = 0; i < NUMFILES; i++) {
			TH1D* h = ANALYZERS[i]->getRawPlot();
			
			Measurement maxBin;
			maxBin.val = h->GetMaximumBin();

			Measurement energy = ANALYZERS[i]->calibrate(maxBin);
			if (energy.val < 30) {
				FitInfo BaPars;
				using ParGuess = pair<Int_t, Double_t>;
				using ParLimit = pair<Int_t, ParWindow>;
				
				BaPars.peakEnergies.push_back(356.0129);
				BaPars.peakEnergies.push_back(383.8485);
				BaPars.peakEnergies.push_back(302.8508);
				BaPars.peakEnergies.push_back(276.3989);
				for (Double_t en : BaPars.peakEnergies) {
					PeakInfo currPeak;
					currPeak.energy = en;
					ANALYZERS[i]->getPeakSet().put(currPeak);
				}

				BaPars.fitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2)";
				BaPars.fitFunc += " + [3]*exp(-0.5*((x-[4])/[2])^2)";
				BaPars.fitFunc += " + [5]*exp(-0.5*((x-[6])/[2])^2)";
				BaPars.fitFunc += " + [7]*exp(-0.5*((x-[8])/[2])^2)";
				BaPars.fitFunc += " + exp([9]+[10]*x)";

				PeakInfo Ba356 = ANALYZERS[i]->findPeak(356.0129);
				
				BaPars.fitPars.insert(ParGuess (0, Ba356.count)); 
				BaPars.fitPars.insert(ParGuess (1, Ba356.mu));
				BaPars.fitPars.insert(ParGuess (2, 0.05 * Ba356.mu)); 
				BaPars.fitPars.insert(ParGuess (3, 8.94 / 62.05 * Ba356.count));
				BaPars.fitPars.insert(ParGuess (4, 383.8485 / 356.0129 * Ba356.mu)); 
				BaPars.fitPars.insert(ParGuess (5, 18.34 / 62.05 * Ba356.count));
				BaPars.fitPars.insert(ParGuess (6, 302.8508 / 356.0129 * Ba356.mu));
				BaPars.fitPars.insert(ParGuess (7, 7.16 / 62.05 * Ba356.count));
				BaPars.fitPars.insert(ParGuess (8, 276.3989 / 356.0129 * Ba356.mu));

				ParWindow Ba383Window;
				Ba383Window.low = (383.8485 / 356.0129 - 0.03) * Ba356.mu;
				Ba383Window.high = (383.8485 / 356.0129 + 0.03) * Ba356.mu;
				BaPars.fitParLimits.insert(ParLimit (4, Ba383Window));

				ParWindow Ba356Window;
				Ba356Window.low = 0.97 * Ba356.mu;
				Ba356Window.high = 1.03 * Ba356.mu;
				BaPars.fitParLimits.insert(ParLimit (1, Ba356Window));

				ParWindow Ba302Window;
				Ba302Window.low = (302.8508 / 356.0129 - 0.03) * Ba356.mu;
				Ba302Window.high = (302.8508 / 356.0129 + 0.03) * Ba356.mu;
				BaPars.fitParLimits.insert(ParLimit (6, Ba302Window));

				ParWindow Ba276Window;
				Ba276Window.low = (276.3989 / 356.0129 - 0.03) * Ba356.mu;
				Ba276Window.high = (276.3989 / 356.0129 + 0.03) * Ba356.mu;
				BaPars.fitParLimits.insert(ParLimit (8, Ba276Window));
				
				BaPars.fitWindow.low = 0.6 * Ba356.mu;
				BaPars.fitWindow.high = 1.25 * Ba356.mu;
				BaPars.backgroundRange = 0.3;
				
				ANALYZERS[i]->fit(BaPars);
				ANALYZERS[i]->findCalibration();
				
				// redraw hist
				fitCanvas->cd(i + 1);
				gPad->SetLogy();
				Float_t range = 1.15*ANALYZERS[i]->getPeakSet().get(2614.511).mu;
				ANALYZERS[i]->getRawPlot()->GetXaxis()->SetRangeUser(0, range);
				ANALYZERS[i]->getRawPlot()->Draw();
			}
		}
		
	}
	if (option.find("muon") != string::npos) {
		
		for (Int_t i = 0; i < NUMFILES; i++) {
			Measurement maxEnergy;
			maxEnergy.val = DATA[i]->GetMaximum("energy");
			maxEnergy = ANALYZERS[i]->calibrate(maxEnergy);
			
			if (maxEnergy.val > 30000) {				
				PeakSet peaks = ANALYZERS[i]->getPeakSet();
				PeakInfo highEPeak = peaks.getHighestEnergyPeak();
				Double_t pos = 25000 / highEPeak.energy * highEPeak.mu;
				
				TH1D *h = (TH1D*) ANALYZERS[i]->getRawPlot()->Clone();
				h->Rebin(4);
				h->GetXaxis()->SetRangeUser(0.9 * pos, 1.1 * pos);
				pos = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
				h->GetXaxis()->SetRangeUser(0, h->GetNbinsX());

				TF1 *muonFit = new TF1("muonFit", "landau", 0.8*pos, 1.2*pos);
				
				PeakInfo muonInfo;
				muonInfo.energy = 25000;
				muonInfo.mu = muonFit->GetParameter(0);
				muonInfo.muErr = muonFit->GetParError(0);
				muonInfo.sigma = muonFit->GetParameter(1);
				muonInfo.sigmaErr = muonFit->GetParError(1);

				ANALYZERS[i]->getPeakSet().put(muonInfo);
				ANALYZERS[i]->findCalibration();
			}
		}

	}
	if (option.find("cal") != string::npos) {
		
		TCanvas* calCanvas = new TCanvas("calCanvas", "calCanvas", 1);
		TMultiGraph* calComp = new TMultiGraph("calComp", "calComp");
		
		for (Int_t i = 0; i < NUMFILES; i++) {
			TGraphErrors* calPlot = ANALYZERS[i]->getCalPlot();
			string label;
			if (mode == "pos") {
				label = "Position " + to_string(POSITIONS[i]);
			} else if (mode == "volt") {
				label = to_string(VOLTAGES[i]) + " V";
			}
			calPlot->SetTitle(label.c_str());
			calPlot->SetMarkerColor(i+1);
			calPlot->SetLineColor(i+1);
			if (i == 4) {
				calPlot->SetMarkerColor(i+2);
				calPlot->SetLineColor(i+2);
			}
			calPlot->SetMarkerStyle(21);
			calPlot->SetLineWidth(2);
			calComp->Add(calPlot);
		}

		string titleVar;
		if (mode == "pos") {
			titleVar = "Position";
		} else if (mode == "volt") {
			titleVar = "Voltage";
		}

		calComp->SetTitle(("Calibration Curves for Each " + titleVar).c_str());
		calComp->GetYaxis()->SetTitle("ADC Energy");
		calComp->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		calComp->GetXaxis()->SetRangeUser(300, 3000);

		calComp->Draw("ALP");
		calCanvas->BuildLegend(0.15, 0.6, 0.30, 0.85);   // legend in top left

	}
	if (option.find("sig") != string::npos) {
		
		TCanvas* sigmaCanvas = new TCanvas("Sigma Canvas", "Sigma Canvas", 1);
		TMultiGraph* sigmaComp = new TMultiGraph("SigmaComp", "SigmaComp");

		for (Int_t i = 0; i < NUMFILES; i++) {
			// get all fitted peak energies and calibrate them:
			PeakSet peaks = ANALYZERS[i]->getPeakSet();
			FitResults calib = ANALYZERS[i]->getCalibration();
			
			vector<Double_t> energies;			
			vector<Double_t> calibratedSigmas;
			vector<Double_t> calibratedSigmaErrs;
			for (Int_t j = 0; j < peaks.size(); j++) {
				PeakInfo peak = peaks.getAtIndex(j);
				energies.push_back(peak.energy);

				Measurement fitPos;
				fitPos.val = peak.mu;
				fitPos.err = peak.muErr;

				Measurement fitSigma;
				fitSigma.val = peak.sigma;
				fitSigma.err = peak.sigmaErr;

				Double_t term1 = TMath::Power(fitSigma.err / fitSigma.val, 2);
				Double_t term2 = TMath::Power(calib.slopeErr / calib.slope, 2);

				Measurement calEnergy = ANALYZERS[i]->calibrate(fitPos);
				Measurement calSigma;
				calSigma.val = fitSigma.val / calib.slope;
				calSigma.err = calSigma.val * TMath::Sqrt(term1 + term2);

				calibratedSigmas.push_back(calSigma.val);
				calibratedSigmaErrs.push_back(calSigma.err);
			}
			
			TGraphErrors* sigmaPlot = new TGraphErrors(energies.size(), &energies[0],
			              &calibratedSigmas[0], 0, &calibratedSigmaErrs[0]);
			string label;
			if (mode == "pos") {
				label = "Position " + to_string(POSITIONS[i]);
			} else if (mode == "volt") {
				label = to_string(VOLTAGES[i]) + " V";
			}
			sigmaPlot->SetTitle(label.c_str());
			sigmaPlot->SetMarkerColor(i+1);
			sigmaPlot->SetLineColor(i+1);
			if (i == 4) {
				sigmaPlot->SetMarkerColor(i+2);
				sigmaPlot->SetLineColor(i+2);
			}
			sigmaPlot->SetMarkerStyle(21);
			sigmaPlot->SetLineWidth(1);

			sigmaComp->Add(sigmaPlot);
		}

		string titleVar;
		if (mode == "pos") {
			titleVar = "Position";
		} else if (mode == "volt") {
			titleVar = "Voltage";
		}
		
		sigmaComp->SetTitle(("Resolution vs " + titleVar).c_str());
		sigmaComp->GetYaxis()->SetTitle("Peak Width (keV)");
		sigmaComp->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		sigmaComp->GetXaxis()->SetRangeUser(300, 3000);

		sigmaComp->Draw("ALP");
		sigmaCanvas->BuildLegend(0.7,0.6,0.85,0.85);   // legend in top right

	}
	if (option.find("res") != string::npos) {
		
		TCanvas* resCanvas = new TCanvas("Residue Canvas", "Residue Canvas", 1);
		TMultiGraph *resComp = new TMultiGraph("ResComp", "ResComp");
		
		for (Int_t i = 0; i < NUMFILES; i++) {
			// get all fitted peak energies and calibrate them:
			PeakSet peaks = ANALYZERS[i]->getPeakSet();

			vector<Double_t> energies;
			vector<Double_t> energyResidues;
			vector<Double_t> energyResidueErrs;
			for (Int_t j = 0; j < peaks.size(); j++) {
				PeakInfo peak = peaks.getAtIndex(j);
				energies.push_back(peak.energy);
				
				Measurement fittedEnergy;
				fittedEnergy.val = peak.mu;
				fittedEnergy.err = peak.muErr;

				Measurement calibrated = ANALYZERS[i]->calibrate(fittedEnergy);
				energyResidues.push_back(calibrated.val - peak.energy);
				energyResidueErrs.push_back(calibrated.err);
			}
		
			TGraphErrors *residuePlot = new TGraphErrors(energies.size(), 
			                            &energies[0], &energyResidues[0], 0, 
			                            &energyResidueErrs[0]);
			string label;
			if (mode == "pos") {
				label = "Position " + to_string(POSITIONS[i]);
			} else if (mode == "volt") {
				label = to_string(VOLTAGES[i]) + " V";
			}
			residuePlot->SetTitle(label.c_str());
			residuePlot->SetMarkerColor(i+1);
			residuePlot->SetLineColor(i+1);
			if (i == 4) {
				residuePlot->SetMarkerColor(i+2);
				residuePlot->SetLineColor(i+2);
			}			
			residuePlot->SetMarkerStyle(21);
			residuePlot->SetLineWidth(1);
			resComp->Add(residuePlot);
		}

		string titleVar;
		if (mode == "pos") {
			titleVar = "Position";
		} else if (mode == "volt") {
			titleVar = "Voltage";
		}

		resComp->SetTitle(("Residues for " + titleVar + " Variation").c_str());
		resComp->GetXaxis()->SetTitle("ADC Energies");
		resComp->GetXaxis()->SetTitleOffset(1.3);
		resComp->GetYaxis()->SetTitle("Error in calibrated energy (keV)");
		resComp->GetYaxis()->SetTitleOffset(1.2);
		resComp->GetXaxis()->SetRangeUser(300, 3000);

		resComp->Draw("ALP");
		resCanvas->BuildLegend(0.15,0.15,0.30,0.40);   // legend in bottom left

	}
	if (option.find("over") != string::npos) {
		
		TCanvas *overlayCanvas = new TCanvas("overlayCanvas", "Overlay Canvas", 1);
		gPad->SetLogy();
		
		for (Int_t i = 0; i < NUMFILES; i++) {
			// need to generate a calibrated histogram
			string calName = "calibrated" + to_string(i);
			string label;
			if (mode == "pos") {
				label = "Position " + to_string(POSITIONS[i]);
			} else if (mode == "volt") {
				label = to_string(VOLTAGES[i]) + " V";
			}
			
			Measurement maxEnergy;
			maxEnergy.val = DATA[i]->GetMaximum("energy");
			maxEnergy.err = 0;

			Measurement calibratedMaxEnergy = ANALYZERS[i]->calibrate(maxEnergy);
			Int_t maxCalBin = (Int_t) (1.01 * calibratedMaxEnergy.val);
			
			TH1D *calibrated = new TH1D(calName.c_str(), label.c_str(), 20e3, 
			                               0, maxCalBin);
			calibrated->SetLineColor(i+1);
			if (i == 4) {
				calibrated->SetLineColor(i+2);
			}
			calibrated->GetXaxis()->SetTitle("Calibrated Energy (keV)");
			calibrated->GetYaxis()->SetTitle("Count");
			
			FitResults calib = ANALYZERS[i]->getCalibration();
			string calibration = "(energy - " + to_string(calib.offset);
			calibration += ") / " + to_string(calib.slope);
			calibration += " >> " + calName;
			DATA[i]->Draw(calibration.c_str(), CHANNEL.c_str(), "SAME");
		}

		overlayCanvas->BuildLegend(0.7,0.6,0.85,0.85); 	// legend in top right

	}
	if (option.find("rawOver") != string::npos) {

		TCanvas* rawOverCanvas = new TCanvas("rawOverlayCanvas", "rawOverlayCanvas");
		gPad->SetLogy();
		for (Int_t i = 0; i < NUMFILES; i++) {
			TH1D *raw = ANALYZERS[i]->getRawPlot();
			raw->SetLineColor(i+1);
			if (i == 4) {
				raw->SetLineColor(i+2);
			}
			raw->GetYaxis()->SetTitle("Count");
			raw->GetXaxis()->SetTitle("Uncalibrated Energy");
			raw->Draw("SAME");
		}
	
	}
	if (option.find("gain") != string::npos) {

		vector<Double_t> gains;
		vector<Double_t> gainErrs;
		vector<Double_t> xAxis;
		for (Int_t i = 0; i < NUMFILES; i++) {
			FitResults calib = ANALYZERS[i]->getCalibration();
			Double_t gain = TMath::Log(calib.slope);
			Double_t gainErr = calib.slopeErr / calib.slope;
			gains.push_back(gain);
			gainErrs.push_back(gainErr);

			if (mode == "pos") {
				xAxis.push_back(POSITIONS[i]);
			} else if (mode == "volt") {
				xAxis.push_back(VOLTAGES[i]);
			}
		}
		
		TCanvas *gainCanvas = new TCanvas("Gain Canvas", "Gain Canvas", 1);
		TGraphErrors *gainGraph = new TGraphErrors(NUMFILES, &xAxis[0], &gains[0], 
		                                           0, &gainErrs[0]);
		string titleVar;
		if (mode == "pos") {
			titleVar = "Position";
		} else if (mode == "volt") {
			titleVar = "Voltage";
		}

		gainGraph->SetTitle(("Detector Gain vs " + titleVar).c_str());
		gainGraph->GetXaxis()->SetTitle(titleVar.c_str());
		gainGraph->GetYaxis()->SetTitle("Log(Calibration Slope)");
		gainGraph->SetMarkerColor(4);
		gainGraph->SetMarkerStyle(21);
		gainGraph->SetLineColor(1);
		gainGraph->SetLineWidth(2);
		gainGraph->GetYaxis()->SetTitleOffset(1.4);
		gainGraph->GetXaxis()->SetTitleOffset(1.2);
		
		TF1 *gainFit = new TF1("gainFit", "pol2");
		gainFit->SetParNames("Log(G0)", "Slope", "Curvature");
		gainGraph->Fit("gainFit");
		gainGraph->Draw();
		
	}
	if (option.find("noise") != string::npos) {
		
		vector<Double_t> noiseWallEnergies;
		vector<Double_t> noiseWallEnergyErrs;
		vector<Double_t> xAxis;
		for (Int_t i = 0; i < NUMFILES; i++) {
			TH1D* h = ANALYZERS[i]->getRawPlot();
			FitResults calib = ANALYZERS[i]->getCalibration();
			Double_t slope = calib.slope;
			Double_t offset = calib.offset;
						
			// find uncalibrated energy that corresponds to 50 keV:
			Double_t startEnergy = slope * 50.0 + offset;
			Int_t binNum = h->FindBin(startEnergy);
			Int_t count = h->GetBinContent(binNum);
			Int_t minimum = count;
			while (count < 2 * minimum && binNum > 0) {
				count = h->GetBinContent(binNum);
				if (count < minimum) {
					minimum = count;
				}
				binNum--;
			}
			Measurement noiseWall;
			noiseWall.val = h->GetBinCenter(binNum);
			noiseWall.err = h->GetBinWidth(binNum) / 2.0;
			Measurement calibratedNoiseWall = ANALYZERS[i]->calibrate(noiseWall);
			
			noiseWallEnergies.push_back(calibratedNoiseWall.val);
			noiseWallEnergyErrs.push_back(calibratedNoiseWall.err);

			if (mode == "pos") {
				xAxis.push_back(POSITIONS[i]);
			} else if (mode == "volt") {
				xAxis.push_back(VOLTAGES[i]);
			}
		}
		
		TCanvas *noiseCanvas = new TCanvas("noiseCanvas", "NoiseCanvas", 1);
		TGraphErrors *noiseGraph = new TGraphErrors(NUMFILES, &xAxis[0], 
		                           &noiseWallEnergies[0], 0, &noiseWallEnergyErrs[0]);
		string titleVar;
		if (mode == "pos") {
			titleVar = "Position";
		} else if (mode == "volt") {
			titleVar = "Voltage";
		}
		
		noiseGraph->SetTitle(("Noise Wall Energy vs " + titleVar).c_str());
		noiseGraph->GetXaxis()->SetTitle(titleVar.c_str());
		noiseGraph->GetYaxis()->SetTitle("Noise Wall Energy (keV)");
		noiseGraph->SetMarkerColor(4);
		noiseGraph->SetMarkerStyle(21);
		noiseGraph->SetLineColor(1);
		noiseGraph->SetLineWidth(2);
		noiseGraph->GetYaxis()->SetTitleOffset(1.4);
		noiseGraph->GetXaxis()->SetTitleOffset(1.2);	
		noiseGraph->Draw();
		
		if (mode == "volt") {
			noiseCanvas->Print("Noise.pdf[");
			noiseCanvas->Print("Noise.pdf");
			noiseCanvas->Print("Noise.pdf]");
		}
		
	}
	if (option.find("rate") != string::npos) {

		vector<Double_t> rates;
		vector<Double_t> xAxis;
		for (Int_t i = 0; i < NUMFILES; i++) {
			Double_t nEntry = (Double_t) DATA[i]->GetBranch("energy")->GetEntries();
			Double_t time = 600 * DATA[i]->GetNtrees();
			rates.push_back(nEntry / time);

			if (mode == "pos") {
				xAxis.push_back(POSITIONS[i]);
			} else if (mode == "volt") {
				xAxis.push_back(VOLTAGES[i]);
			}
		}

		TCanvas* rateCanvas = new TCanvas("rateCanvas", "Rate Canvas", 1);
		TGraphErrors* rateGraph = new TGraphErrors(rates.size(), &xAxis[0], 
		                                           &rates[0], 0, 0);
		string titleVar;
		if (mode == "pos") {
			titleVar = "Position";
		} else if (mode == "volt") {
			titleVar = "Voltage";
		}

		rateGraph->SetTitle(("Count Rate vs " + titleVar).c_str());
		rateGraph->GetXaxis()->SetTitle(titleVar.c_str());
		rateGraph->GetYaxis()->SetTitle("Rate (1/seconds)");
		rateGraph->SetMarkerColor(4);
		rateGraph->SetMarkerStyle(21);
		rateGraph->SetLineColor(1);
		rateGraph->SetLineWidth(2);
		rateGraph->GetYaxis()->SetTitleOffset(1.4);
		rateGraph->GetXaxis()->SetTitleOffset(1.2);
		rateGraph->Draw();

	}
	if (option.find("back") != string::npos) {
		
		cout << "checkpoint 1" << endl;
		TCanvas* backCanvas = new TCanvas("backCanvas", "Background Fits", 1400, 900);
		backCanvas->Divide(peakPars.size(), NUMFILES); // 3 columns, 5 rows
		cout << "checkpoint 2" << endl;
		for (Int_t i = 0; i < NUMFILES; i++) {
			vector<TGraphErrors*> backPlots = ANALYZERS[i]->getBackgroundPlots();
			cout << "checkpoint 3" << endl;
			for (Int_t j = 1; j < backPlots.size() + 1; j++) {
				backCanvas->cd(3 * i + j);
				cout << "checkpoint 4: cding to " << 3 * i + j << endl;
				backPlots[i]->Draw();
				cout << "checkpoint 5: successfully drew plot" << endl;
				// fails on printing the 10th graph (1st one in 4th row
			}
			cout << "checkpoint 6" << endl;
		}

	}
	if (option.find("AE") != string::npos && mode == "pos") {
		
		TCanvas *AECanvas = new TCanvas("AECanvas", "A/E Canvas", 1);
		TH2D *AEHist = new TH2D("AEHist", "Amplitude / Energy vs calibrated Energy",
		                        1e3, 0, 50e3, 1e3, 0, 10);
		
		TChain* allData = new TChain("st");
		for (TChain* c : DATA) {
			allData->Add(c);
		}

		FitResults calib = ANALYZERS[NUMFILES / 2]->getCalibration();

		string calE = "(energy-"+to_string(calib.offset)+")/"+to_string(calib.slope);
		string toPlot = "amp / (" + calE + ") : (" + calE + ") >> AEHist";
		allData->Draw(toPlot.c_str(), CHANNEL.c_str(), "COLZ");
		
		AEHist->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		AEHist->GetYaxis()->SetTitle("Amplitude / Callibrated Energy");

	}

	app->Run(false);
	return 0;

}

/*
Double_t landauFunc(Double_t *x, Double_t *par) {
	Double_t landau = TMath::Landau(x, par[0], par[1]);
	Double_t expo = TMath::Exp(x);
}
*/

int main(int argc, char** argv) {
	if (argc == 1) {
		Calibration(argv[1], "pos", "");
		Calibration(argv[1], "volt", "");
	} else if (argc == 3) {
		Calibration(argv[1], argv[2], "");
	} else if (argc == 4) {
		Calibration(argv[1], argv[2], argv[3]);
	} else {
		cout << "Invalid arguments. Allowed arguments: <path> <mode> <option>" << endl;
		cout << "See protocol for more info on usage of calibration script." << endl;
		return 1;
	}
	return 0;
}
