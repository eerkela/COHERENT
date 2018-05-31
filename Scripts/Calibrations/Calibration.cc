#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <stdio.h>

/*
This root script will analyze an NaI characterization suite containing raw digitizer output
data and calibrate them to a real energy scale.  More detailed analysis options are described
below.


Modes:
"pos"	will execute calibration algorithm for position data.
"volt"	will execute calibration algorithm for voltage data.
"fir"	will perform FIR filter optimization.


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


Directory structure:
Crystal Serial #/
	Characterization.cc
	FIR/
		.../1/NaI_ET_run*.root
		.../2/NaI_ET_run*.root
		.../3/NaI_ET_run*.root
		.../4/NaI_ET_run*.root
		.../5/NaI_ET_run*.root
	Position/
		.../1/NaI_ET_run*.root
		.../2/NaI_ET_run*.root
		.../3/NaI_ET_run*.root
		.../4/NaI_ET_run*.root
		.../5/NaI_ET_run*.root
	Voltage/
		.../600V/NaI_ET_run*.root
		.../700V/NaI_ET_run*.root
		.../800V/NaI_ET_run*.root
		.../900V/NaI_ET_run*.root
		.../1000V/NaI_ET_run*.root	


Notes and TODO:
	
	Current script version does not pull peaking times directly from the .roots.  Instead they
	are hardcoded in. We've modified getSpectrum to pull them out of the file itself, but
	when I try to access them, I just get zeros.

	Fix Cs fit.  There's a few peaks underneath it that kind of screw it up.  Might want to
	fix some parameters.

	Figure out how to do math on TChains.  I want to store the calibrated data in the RESULTS
	struct somewhere, but I need to apply the calibration directly to the data in that case.

	Estimate Tl position without drawing the hist.  This wouldn't change the functionality,
	but would save on runtime.

*/

using namespace std;

//////// BEGIN INTERNAL STRUCTS ////////

struct PeakInfo {
	Double_t mu;
	Double_t muErr;
	Double_t sigma;
	Double_t sigmaErr;
	Double_t count;
};

struct FitPars {
	Double_t offset;
	Double_t offsetErr;
	Double_t slope;
	Double_t slopeErr;
};

struct Residues {
	Double_t Tl;
	Double_t K;
	Double_t Cs;
};

struct BackgroundPlots {
	TGraph *Tl;
	TGraph *K;
	TGraph *Cs;
};

//////// END INTERNAL STRUCTS ////////

struct FitWindow {
	Double_t low;
	Double_t high;
};

struct Peaks {
	PeakInfo Tl;
	PeakInfo K;
	PeakInfo Cs;
};

struct Plots {
	TH1D *raw;
	TH1D *calibrated;
	TH1D *muPlot;

	TGraphErrors *calPlot;
	TGraphErrors *sigmaPlot;

	BackgroundPlots backPlots;
};

struct Results {
	FitPars calPars;	
	Residues residues;
	TChain calibratedData;
	
	Double_t noiseWall;
	Double_t bestPeakTime;
	Double_t lowERate;
};

Int_t NUMFILES = 5;

FitPars backEst(TH1D *h, FitWindow win, Double_t coverage, string func, string peak, Int_t i);
void formatGraph(TGraphErrors *graph);
void getData(TChain *data[NUMFILES], string mode);
string getDepVar(string mode);
string getLabels(string mode, Int_t index);
Double_t getRunTime(TChain *c, string mode);
Double_t snapToMax(TH1D *h, FitWindow win);

Int_t NUMBINS = 10000;
Int_t VOLTAGES[5] = {750, 810, 870, 930, 990};
Int_t PEAKINGTIMES[5] = {50, 100, 200, 400, 800};
Int_t POSITIONS[5] = {1, 2, 3, 4, 5};

Plots PLOTS[5];
Peaks PEAKS[5];
Results RESULTS[5];


void Calibration(string mode, string option) {

	//gStyle->SetOptFit(1111);
	
	TCanvas *fitCanvas = new TCanvas("Fits", "Fits", 1200, 800);
	fitCanvas->Divide(2, 3);
	TChain *data[NUMFILES]; 
	getData(data, mode);

	for (Int_t i = 0; i < NUMFILES; i++) {
		cout << endl;
		cout << "---------------------------------------------------------------" << endl;
		cout << endl;

		cout << "Beginning Calibration " << i + 1 << "..." << endl;
		fitCanvas->cd(i+1);
		gPad->SetLogy();

		Double_t time = getRunTime(data[i], mode);
		cout << "Run time in data chain: " << time << " seconds" << endl;
		Int_t overflowBin = data[i]->GetMaximum("energy");
		overflowBin += 0.01 * overflowBin;

		TH1D *hTemp = new TH1D("hTemp", "Finding Tl pos", NUMBINS, 0, overflowBin);
		data[i]->Draw("energy >> hTemp", "channel==0");
		
		// Finding Thalium:
		cout << "Finding Thalium Peak..." << endl;
		UInt_t binNum = NUMBINS;
		UInt_t totalCount = 0;
		while (totalCount < 13.6 * time) {
			totalCount += (UInt_t) hTemp->GetBinContent(binNum);
			binNum--;
			if (binNum == 0) { break; }
		}
		Double_t TlPos = hTemp->GetXaxis()->GetBinCenter(binNum);
		Double_t newTlPos;
		cout << "VISUAL CHECK" << endl;
		cout << "Estimated Tl Position: " << TlPos << endl;
		cout << "Does this make sense?" << endl;
		cout << "New Tl position (enter nothing to continue using the above estimation): " << endl;
		cout << " >> ";
		cin >> newTlPos;
		cout << endl;
		if (newTlPos != NULL && newTlPos > 0) {
			TlPos = newTlPos;
		}
		
		FitWindow TlSnapWindow;
		TlSnapWindow.low = TlPos - 0.05 * TlPos;
		TlSnapWindow.high = TlPos + 0.05 * TlPos;		
		TlPos = snapToMax(hTemp, TlSnapWindow);
		Double_t TlNormalizedPos = (Double_t)hTemp->FindBin(TlPos) / (Double_t)NUMBINS;

		cout << "Snapped Tl Position: " << TlPos << endl;
		
		// Redraw hist with constant number of bins below Tl peak (enhances fits).
		delete hTemp;
		NUMBINS = (Int_t) 500/TlNormalizedPos;
		string hName = "h" + to_string(i+1);
		PLOTS[i].raw = new TH1D(hName.c_str(), "Uncalibrated Spectrum", NUMBINS, 0, overflowBin);
		data[i]->Draw(("energy >> " + hName).c_str(), "channel==0");

		TlSnapWindow.low = TlPos - 0.05 * TlPos;
		TlSnapWindow.high = TlPos + 0.05 * TlPos;
		TlPos = snapToMax(PLOTS[i].raw, TlSnapWindow);
		Double_t TlCount = PLOTS[i].raw->GetBinContent(PLOTS[i].raw->FindBin(TlPos));

		FitWindow TlFitWindow;
		TlFitWindow.low = TlPos - 0.1 * TlPos;
		TlFitWindow.high = TlPos + 0.1 * TlPos;

		FitPars TlBackPars = backEst(PLOTS[i].raw, TlFitWindow, 0.4, "pol1", "Tl", i);

		string TlFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2)";
		TlFitFunc += " + [3] + [4] * x";

		Double_t TlParameters[5];
		TlParameters[0] = TlCount;
		TlParameters[1] = TlPos;
		TlParameters[2] = 0.05 * TlPos;
		TlParameters[3] = TlBackPars.offset;
		TlParameters[4] = TlBackPars.slope;

		TF1 *TlFit = new TF1("TlFit", TlFitFunc.c_str(), TlFitWindow.low, TlFitWindow.high);
		TlFit->SetParameters(TlParameters);
		PLOTS[i].raw->Fit("TlFit", "R+ll");
		
		PEAKS[i].Tl.mu = TlFit->GetParameter(1);
		PEAKS[i].Tl.muErr = TlFit->GetParError(1);
		PEAKS[i].Tl.sigma = TlFit->GetParameter(2);
		PEAKS[i].Tl.sigmaErr = TlFit->GetParError(2);
		PEAKS[i].Tl.count = TlFit->GetParameter(0);

		// Finding Potassium:
		cout << "Fitting Potassium peak..." << endl;
		Double_t KGuess = (1460.0/2615.0) * PEAKS[i].Tl.mu;
		
		FitWindow KSnapWindow;
		KSnapWindow.low = KGuess - 0.05 * KGuess;
		KSnapWindow.high = KGuess + 0.05 * KGuess;

		Double_t KPos = snapToMax(PLOTS[i].raw, KSnapWindow);
		Double_t KCount = PLOTS[i].raw->GetBinContent(PLOTS[i].raw->FindBin(KPos));
		
		FitWindow KFitWindow;
		KFitWindow.low = KPos - 0.15 * KPos;
		KFitWindow.high = KPos + 0.15 * KPos;
		
		FitPars KBackPars = backEst(PLOTS[i].raw, KFitWindow, 0.3, "expo", "K", i);
		
		string KFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2)";
		KFitFunc += " + exp([3] + [4] * x)";

		Double_t KParameters[5];
		KParameters[0] = KCount;
		KParameters[1] = KPos;
		KParameters[2] = 0.05 * KPos;
		KParameters[3] = KBackPars.offset;
		KParameters[4] = KBackPars.slope;

		TF1 *KFit = new TF1("KFit", KFitFunc.c_str(), KFitWindow.low, KFitWindow.high);
		KFit->SetParameters(KParameters);
		PLOTS[i].raw->Fit("KFit", "R+ll");
		
		PEAKS[i].K.mu = KFit->GetParameter(1);
		PEAKS[i].K.muErr = KFit->GetParError(1);
		PEAKS[i].K.sigma = KFit->GetParameter(2);
		PEAKS[i].K.sigmaErr = KFit->GetParError(2);
		PEAKS[i].K.count = KFit->GetParameter(0);

		// Finding Cs:
		cout << "Fitting Cesium Peak" << endl;
		Double_t CsGuess = (661.6/2615.0) * PEAKS[i].Tl.mu;
		
		FitWindow CsSnapWindow;
		CsSnapWindow.low = CsGuess - 0.05 * CsGuess;
		CsSnapWindow.high = CsGuess + 0.05 * CsGuess;

		Double_t CsPos = snapToMax(PLOTS[i].raw, CsSnapWindow);
		Double_t CsCount = PLOTS[i].raw->GetBinContent(PLOTS[i].raw->FindBin(CsPos));
		
		FitWindow CsFitWindow;
		CsFitWindow.low = CsPos - 0.2 * CsPos;
		CsFitWindow.high = CsPos + 0.2 * CsPos;
	
		FitPars CsBackPars = backEst(PLOTS[i].raw, CsFitWindow, 0.2, "expo", "Cs", i);

		string CsFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2)";
		//CsFitFunc += " + [3]*exp(-0.5*((x-[4])/[5])^2)";
		CsFitFunc += " + exp([3] + [4] * x)";

		Double_t CsParameters[5];
		CsParameters[0] = CsCount;
		CsParameters[1] = CsPos;
		CsParameters[2] = 0.05 * CsPos;
		CsParameters[3] = CsBackPars.offset;
		CsParameters[4] = CsBackPars.slope;

		TF1 *CsFit = new TF1("CsFit", CsFitFunc.c_str(), CsFitWindow.low, CsFitWindow.high);
		CsFit->SetParameters(CsParameters);
		PLOTS[i].raw->Fit("CsFit", "R+ll");
		
		PEAKS[i].Cs.mu = CsFit->GetParameter(1);
		PEAKS[i].Cs.muErr = CsFit->GetParError(1);
		PEAKS[i].Cs.sigma = CsFit->GetParameter(2);
		PEAKS[i].Cs.sigmaErr = CsFit->GetParError(2);
		PEAKS[i].Cs.count = CsFit->GetParameter(0);

		Float_t range = PEAKS[i].Tl.mu + 0.15 * PEAKS[i].Tl.mu;
		PLOTS[i].raw->GetXaxis()->SetRangeUser(0, range);
		NUMBINS = 10000;	

		// Collecting data:
		Double_t expE[3] = {2615.0, 1460.0, 661.6};
		Double_t fitE[3]; 
			fitE[0] = PEAKS[i].Tl.mu;
			fitE[1] = PEAKS[i].K.mu;
			fitE[2] = PEAKS[i].Cs.mu;
		Double_t fitEErr[3];
			fitEErr[0] = PEAKS[i].Tl.muErr;
			fitEErr[1] = PEAKS[i].K.muErr;
			fitEErr[2] = PEAKS[i].Cs.muErr;
		Double_t fitSigma[3];
			fitSigma[0] = PEAKS[i].Tl.sigma / PEAKS[i].Tl.mu;
			fitSigma[1] = PEAKS[i].K.sigma / PEAKS[i].K.mu;
			fitSigma[2] = PEAKS[i].Cs.sigma / PEAKS[i].Cs.mu;
		Double_t fitSigmaErr[3];
			fitSigmaErr[0] = PEAKS[i].Tl.sigmaErr / PEAKS[i].Tl.mu;
			fitSigmaErr[1] = PEAKS[i].K.sigmaErr / PEAKS[i].K.mu;
			fitSigmaErr[2] = PEAKS[i].Cs.sigmaErr / PEAKS[i].Cs.mu;
		PLOTS[i].calPlot = new TGraphErrors(3, expE, fitE, 0, fitEErr);
		PLOTS[i].calPlot->SetTitle(getLabels(mode, i).c_str());
		PLOTS[i].calPlot->SetMarkerColor(i+1);
		if (i == 4) {PLOTS[i].calPlot->SetMarkerColor(i+2);}
		PLOTS[i].calPlot->SetMarkerStyle(21);
		PLOTS[i].calPlot->SetLineColor(1);
		PLOTS[i].calPlot->SetLineWidth(2);		
		
		PLOTS[i].sigmaPlot = new TGraphErrors(3, expE, fitSigma, 0, fitSigmaErr); 
		PLOTS[i].sigmaPlot->SetTitle(getLabels(mode, i).c_str());
		PLOTS[i].sigmaPlot->SetMarkerColor(i+1);
		PLOTS[i].sigmaPlot->SetLineColor(1);
		if (i == 4) {PLOTS[i].sigmaPlot->SetMarkerColor(i+2);}
		PLOTS[i].sigmaPlot->SetMarkerStyle(21);
		PLOTS[i].sigmaPlot->SetLineWidth(1);
		
		TF1 *calibrationFit = new TF1("calibrationFit", "pol1", 0, overflowBin);
		PLOTS[i].calPlot->Fit("calibrationFit", "", "", 0, 0);
		
		RESULTS[i].calPars.slope = calibrationFit->GetParameter(1);
		RESULTS[i].calPars.slopeErr = calibrationFit->GetParError(1);
		RESULTS[i].calPars.offset = calibrationFit->GetParameter(0);
		RESULTS[i].calPars.offsetErr = calibrationFit->GetParError(0);
		//RESULTS[i].calibratedData = (data[i] - RESULTS[i].calPars.offset) 
		//                             / RESULTS[i].calPars.slope;

		string calName = "calibrated" + to_string(i);
		string label = getLabels(mode, i);
		Int_t maxCalBin = (overflowBin - RESULTS[i].calPars.offset) / 
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
			rateCount += PLOTS[i].calibrated->GetBinCotent(rateBin);
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
		cout << "---------------------------------------------------------------" << endl;
		cout << endl;

	}

	//////// BEGIN MODE HANDLING ////////

	if (mode == "fir") {

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

	} else if (mode == "pos") {

		Double_t resolution = PEAKS[3].K.sigma / PEAKS[3].K.mu;
		cout << "Normalized detector resolution (width of Cs peak" << endl;
		cout << " / mean at 3rd position):" << resolution << endl;
		
		Double_t CsPos[NUMFILES];
		Double_t CsPosErrs[NUMFILES];
		Double_t CsSig[NUMFILES];
		Double_t CsSigErrs[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {
			CsPos[k] = PEAKS[k].K.mu;
			CsPosErrs[k] = PEAKS[k].K.muErr;
			CsSig[k] = PEAKS[k].K.sigma / PEAKS[k].K.mu;
			CsSigErrs[k] = PEAKS[k].K.sigmaErr / PEAKS[k].K.mu;
		}

		Double_t xAxis[NUMFILES];
		for (Int_t i = 0; i < 5; i++) {
			if (mode == "pos") {
				xAxis[i] = POSITIONS[i];
			} else if (mode == "volt") {
				xAxis[i] = VOLTAGES[i];
			} else if (mode == "fir") {
				xAxis[i] = PEAKINGTIMES[i];
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

	//////// END MODE HANDLING ////////

	//////// BEGIN OPTION HANDLING ////////

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
			
			Double_t energy = (PEAKS[k].Tl.mu - RESULTS[k].calPars.offset) / 
			              RESULTS[k].calPars.slope;
			residues[0] = 100 * (energy - expE[0]) / expE[0];
				
			energy = (PEAKS[k].K.mu - RESULTS[k].calPars.offset) / 
			          RESULTS[k].calPars.slope;
			residues[1] = 100 * (energy - expE[1]) / expE[1];			
			
			energy = (PEAKS[k].Cs.mu - RESULTS[k].calPars.offset) / 
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
				xAxis[i] = PEAKINGTIMES[i];
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
				xAxis[i] = PEAKINGTIMES[i];
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
		
		TCanvas *backCanvas = new TCanvas("backCanvas", "Background Fits", 1400, 900);
		backCanvas->Divide(3, 5);
		for (Int_t i = 0; i < 5; i++) {
			backCanvas->cd(3*i+1);
			PLOTS[i].backPlots.Tl->Draw();
			backCanvas->cd(3*i+2);
			PLOTS[i].backPlots.K->Draw();
			backCanvas->cd(3*i+3);
			PLOTS[i].backPlots.Cs->Draw();
		}

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
		path[0] = "Position/1/NaI_ET_run*";
		path[1] = "Position/2/NaI_ET_run*";
		path[2] = "Position/3/NaI_ET_run*";
		path[3] = "Position/4/NaI_ET_run*";
		path[4] = "Position/5/NaI_ET_run*";
	} else if (mode == "volt") {
		cout << "Finding Voltage Data..." << endl;
		path[0] = "Voltage/750V/NaI_ET_run*";
		path[1] = "Voltage/810V/NaI_ET_run*";
		path[2] = "Voltage/870V/NaI_ET_run*";
		path[3] = "Voltage/930V/NaI_ET_run*";
		path[4] = "Voltage/990V/NaI_ET_run*";
	} else if (mode == "fir") {
		cout << "Finding FIR Data..." << endl;
		path[0] = "FIR/1/NaI_ET_run*";
		path[1] = "FIR/2/NaI_ET_run*";
		path[2] = "FIR/3/NaI_ET_run*";
		path[3] = "FIR/4/NaI_ET_run*";
		path[4] = "FIR/5/NaI_ET_run*";
	}
	
	for (Int_t i = 0; i < NUMFILES; i++) {
		data[i] = new TChain("st");
		data[i]->Add(path[i].c_str());
	}
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

Double_t snapToMax(TH1D *h, FitWindow win) {
	h->GetXaxis()->SetRangeUser(win.low, win.high);
	Double_t pos = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
	h->GetXaxis()->SetRange(0, NUMBINS);
	return pos;
}

FitPars backEst(TH1D *h, FitWindow win, Double_t coverage, string func, string peak, Int_t i) {
	cout << "Estimating Background..." << endl;

	// This would be if you wanted to do an asymmetric fit or something
	/*FitWindow lowWindow;
	lowWindow.low = win.low;
	lowWindow.high = win.low + 0.1 * win.low;
	Int_t lowWindowLowBin = h->FindBin(lowWindow.low);
	Int_t lowWindowHighBin = h->FindBin(lowWindow.high);
	Int_t lowWindowBinRange = lowWindowHighBin - lowWindowLowBin;

	FitWindow highWindow;
	highWindow.low = win.high - 0.1 * win.high;
	highWindow.high = win.high;
	Int_t highWindowLowBin = h->FindBin(highWindow.low);
	Int_t highWindowHighBin = h->FindBin(highWindow.high);
	Int_t highWindowBinRange = highWindowHighBin - highWindowLowBin;
	
	Float_t backVals[lowWindowBinRange + highWindowBinRange];
	Float_t backPos[lowWindowBinRange + highWindowBinRange];*/

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

	FitPars pars;
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
}

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
		return "Peaking Time: " + to_string(PEAKINGTIMES[index]);
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
