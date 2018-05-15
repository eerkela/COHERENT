#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <stdio.h>

// This root script will process a series of sequential NaI .root files containing raw low gain
// data and calibrate them to a real energy scale.  Some useful options are available for 
// displaying relevant information.


// Options:
// "cal" 	will display the graphs used to generate the calibrations for each crystal.
// "res"	will display residues from the calibration algorithm, i.e. the % error by which the
//				calibrated energies deviate from expected. (NYI)
// "sig" 	will display fractional sigmas for each peak in order to measure detector
//				resolution.
// (else)	will simply calibrate the spectra then overlay the histograms with real energies.

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
	PeakInfo muon;
};

struct Plots {
	TH1D *raw;
	TH1D *muRaw;
	TH1D *calibrated;
	TH1D *muCalibrated;

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
};

FitPars backEst(TH1D *h, FitWindow win, Double_t coverage, string func, string peak);
Double_t snapToMax(TH1D *h, FitWindow win);

Int_t NUMBINS = 10000;
string CHANNEL = "channel==0";

Plots PLOTS;
Peaks PEAKS;
Results RESULTS;

void SingleLowGainCalibration(string option) {
	
	TCanvas *fitCanvas = new TCanvas("fitCanvas", "Fits", 1);
	gPad->SetLogy();
	TChain *data = new TChain("st");
	data->Add("voltage/750_V/NaI_ET_run*");

	Double_t time = (Double_t) data->GetNtrees() * 600;
	cout << "run time in seconds: " << time << endl;
	Int_t overflowBin = data->GetMaximum("energy");
	overflowBin += 0.01 * overflowBin;

	TH1D *hTemp = new TH1D("hTemp", "Uncalibrated Spectrum", NUMBINS, 0, overflowBin);
	data->Draw("energy >> hTemp", CHANNEL.c_str());

	// Finding Thalium:
	cout << "Finding Thalium Peak..." << endl;
	UInt_t binNum = NUMBINS; 
	UInt_t totalCount = 0;
	while (totalCount < (9.5 * time)) { // 
		totalCount += (UInt_t) hTemp->GetBinContent(binNum);
		binNum--;
		if (binNum == 0) { break; }
	}
	Double_t TlPos = hTemp->GetXaxis()->GetBinCenter(binNum);
	cout << "Estimated Tl Position: " << TlPos << endl;

	FitWindow TlSnapWindow;
	TlSnapWindow.low = TlPos - 0.1 * TlPos;
	TlSnapWindow.high = TlPos + 0.2 * TlPos;		
	TlPos = snapToMax(hTemp, TlSnapWindow);
	Double_t TlNormalizedPos = (Double_t)hTemp->FindBin(TlPos) / (Double_t)NUMBINS;

	cout << "Snapped Tl Position: " << TlPos << endl;
	
	// Redraw hist with constant number of bins below Tl peak (enhances fits).
	delete hTemp;
	NUMBINS = (Int_t) 500/TlNormalizedPos;
	PLOTS.raw = new TH1D("uncalibrated", "Uncalibrated Spectrum", NUMBINS, 0, overflowBin);
	data->Draw("energy >> uncalibrated", CHANNEL.c_str());

	TlSnapWindow.low = TlPos - 0.05 * TlPos;
	TlSnapWindow.high = TlPos + 0.05 * TlPos;
	TlPos = snapToMax(PLOTS.raw, TlSnapWindow);
	Double_t TlCount = PLOTS.raw->GetBinContent(PLOTS.raw->FindBin(TlPos));

	FitWindow TlFitWindow;
	TlFitWindow.low = TlPos - 0.1 * TlPos;
	TlFitWindow.high = TlPos + 0.1 * TlPos;

	FitPars TlBackPars = backEst(PLOTS.raw, TlFitWindow, 0.4, "pol1", "Tl");

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
	PLOTS.raw->Fit("TlFit", "R+ll");
	
	PEAKS.Tl.mu = TlFit->GetParameter(1);
	PEAKS.Tl.muErr = TlFit->GetParError(1);
	PEAKS.Tl.sigma = TlFit->GetParameter(2);
	PEAKS.Tl.sigmaErr = TlFit->GetParError(2);
	PEAKS.Tl.count = TlFit->GetParameter(0);

	// Finding Potassium:
	cout << "Fitting Potassium peak..." << endl;
	Double_t KPos = (1460.0/2615.0) * PEAKS.Tl.mu;
	
	FitWindow KSnapWindow;
	KSnapWindow.low = KPos - 0.05 * KPos;
	KSnapWindow.high = KPos + 0.05 * KPos;

	KPos = snapToMax(PLOTS.raw, KSnapWindow);
	Double_t KCount = PLOTS.raw->GetBinContent(PLOTS.raw->FindBin(KPos));
	
	FitWindow KFitWindow;
	KFitWindow.low = KPos - 0.15 * KPos;
	KFitWindow.high = KPos + 0.15 * KPos;
	
	FitPars KBackPars = backEst(PLOTS.raw, KFitWindow, 0.3, "expo", "K");
	
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
	PLOTS.raw->Fit("KFit", "R+ll");
	
	PEAKS.K.mu = KFit->GetParameter(1);
	PEAKS.K.muErr = KFit->GetParError(1);
	PEAKS.K.sigma = KFit->GetParameter(2);
	PEAKS.K.sigmaErr = KFit->GetParError(2);
	PEAKS.K.count = KFit->GetParameter(0);

	// Finding Cs:
	cout << "Fitting Cesium Peak..." << endl;
	Double_t CsPos = (661.6/2615.0) * PEAKS.Tl.mu;
	
	FitWindow CsSnapWindow;
	CsSnapWindow.low = CsPos - 0.05 * CsPos;
	CsSnapWindow.high = CsPos + 0.05 * CsPos;

	CsPos = snapToMax(PLOTS.raw, CsSnapWindow);
	Double_t CsCount = PLOTS.raw->GetBinContent(PLOTS.raw->FindBin(CsPos));
	
	FitWindow CsFitWindow;
	CsFitWindow.low = CsPos - 0.2 * CsPos;
	CsFitWindow.high = CsPos + 0.2 * CsPos;

	FitPars CsBackPars = backEst(PLOTS.raw, CsFitWindow, 0.2, "expo", "Cs");

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
	PLOTS.raw->Fit("CsFit", "R+ll");
	
	PEAKS.Cs.mu = CsFit->GetParameter(1);
	PEAKS.Cs.muErr = CsFit->GetParError(1);
	PEAKS.Cs.sigma = CsFit->GetParameter(2);
	PEAKS.Cs.sigmaErr = CsFit->GetParError(2);
	PEAKS.Cs.count = CsFit->GetParameter(0);

	// Finding muon peak:
	cout << "Fitting Muon Peak..." << endl;
	TCanvas *muCanvas = new TCanvas("muCanvas", "muCanvas", 1);
	gPad->SetLogy();
	Int_t muPos = 10 * PEAKS.Tl.mu;

	FitWindow muSnapWindow;
	muSnapWindow.low = 0.7 * muPos;
	muSnapWindow.high = 1.3 * muPos;
	
	PLOTS.muRaw = new TH1D("muRaw", "Muon Plot", 500, 0, overflowBin);
	data->Draw("energy >> muRaw", CHANNEL.c_str());
	muPos = snapToMax(PLOTS.muRaw, muSnapWindow);
	Double_t muCount = PLOTS.muRaw->GetBinContent(PLOTS.muRaw->FindBin(muPos));

	FitWindow muFitWindow;
	muFitWindow.low = 0.65 * muPos;
	muFitWindow.high = 1.8 * muPos;

	FitPars muBackPars = backEst(PLOTS.muRaw, muFitWindow, 0.2, "pol1", "");

	string muFitFunc = "[0]*exp(-0.5*((x-[1])/[2] + exp(-(x-[1])/[2])))";
	muFitFunc = " + [3] + [4] * x";

	Double_t muParameters[5];
	muParameters[0] = muCount;
	muParameters[1] = muPos;
	muParameters[2] = 0.1 * muPos;
	muParameters[3] = muBackPars.offset;
	muParameters[4] = muBackPars.slope;

	TF1 *muFit = new TF1("muFit", muFitFunc.c_str(), muFitWindow.low, muFitWindow.high);
	muFit->SetParameters(muParameters);
	PLOTS.muRaw->Fit("muFit", "R+ll");

	PEAKS.muon.mu = muFit->GetParameter(1);
	PEAKS.muon.muErr = muFit->GetParError(1);
	PEAKS.muon.sigma = muFit->GetParameter(2);
	PEAKS.muon.sigmaErr = muFit->GetParError(2);
	PEAKS.muon.count = muFit->GetParameter(0);

	PLOTS.muRaw->SetTitle("Muon Peak (uncalibrated)");
	PLOTS.muRaw->GetXaxis()->SetTitle("ADC energy");
	PLOTS.muRaw->GetYaxis()->SetTitle("Count");
	PLOTS.muRaw->GetXaxis()->SetRangeUser(0.5 * PEAKS.muon.mu, 2.5 * PEAKS.muon.mu);

	// Resetting parameters:
	Float_t range = PEAKS.Tl.mu + 0.15 * PEAKS.Tl.mu;
	PLOTS.raw->GetXaxis()->SetRangeUser(0, range);
	NUMBINS = 10000;

	// Collecting data:
	Double_t expE[3] = {2615.0, 1460.0, 661.6};
	Double_t fitE[3]; 
		fitE[0] = PEAKS.Tl.mu;
		fitE[1] = PEAKS.K.mu;
		fitE[2] = PEAKS.Cs.mu;
	Double_t fitEErr[3];
		fitEErr[0] = PEAKS.Tl.muErr;
		fitEErr[1] = PEAKS.K.muErr;
		fitEErr[2] = PEAKS.Cs.muErr;
	Double_t fitSigma[3];
		fitSigma[0] = PEAKS.Tl.sigma / PEAKS.Tl.mu;
		fitSigma[1] = PEAKS.K.sigma / PEAKS.K.mu;
		fitSigma[2] = PEAKS.Cs.sigma / PEAKS.Cs.mu;
	Double_t fitSigmaErr[3];
		fitSigmaErr[0] = PEAKS.Tl.sigmaErr / PEAKS.Tl.mu;
		fitSigmaErr[1] = PEAKS.K.sigmaErr / PEAKS.K.mu;
		fitSigmaErr[2] = PEAKS.Cs.sigmaErr / PEAKS.Cs.mu;
	
	PLOTS.calPlot = new TGraphErrors(3, expE, fitE, 0, fitEErr);
	PLOTS.calPlot->SetTitle("Calibration Curve");
	PLOTS.calPlot->SetMarkerColor(1);
	PLOTS.calPlot->SetMarkerStyle(21);
	PLOTS.calPlot->SetLineColor(1);
	PLOTS.calPlot->SetLineWidth(2);
		
	PLOTS.sigmaPlot = new TGraphErrors(3, expE, fitSigma, 0, fitSigmaErr); 
	PLOTS.sigmaPlot->SetTitle("Peak Resolution vs Energy");
	PLOTS.sigmaPlot->SetMarkerColor(1);
	PLOTS.sigmaPlot->SetLineColor(1);
	PLOTS.sigmaPlot->SetMarkerStyle(21);
	PLOTS.sigmaPlot->SetLineWidth(1);
	
	TF1 *calibrationFit = new TF1("calibrationFit", "pol1", 0, overflowBin);
	PLOTS.calPlot->Fit("calibrationFit", "", "", 0, 0);
	
	RESULTS.calPars.slope = calibrationFit->GetParameter(1);
	RESULTS.calPars.slopeErr = calibrationFit->GetParError(1);
	RESULTS.calPars.offset = calibrationFit->GetParameter(0);
	RESULTS.calPars.offsetErr = calibrationFit->GetParError(0);

	Int_t maxCalBin = (overflowBin - RESULTS.calPars.offset) / RESULTS.calPars.slope;
	PLOTS.calibrated = new TH1D("calibrated", "Calibrated Spectrum", 20e3, 0, maxCalBin);
	PLOTS.calibrated->SetLineColor(1);
	PLOTS.calibrated->GetXaxis()->SetTitle("Calibrated Energy (keV)");
	PLOTS.calibrated->GetYaxis()->SetTitle("Count");
	
	string calibration = "(energy - " + to_string(RESULTS.calPars.offset);
	calibration += ") / " + to_string(RESULTS.calPars.slope);
	calibration += " >> calibrated";
	data->Draw(calibration.c_str(), CHANNEL.c_str(), "goff");

	Int_t noiseCount;
	Int_t noiseBin = PLOTS.calibrated->FindBin(50);
	Int_t noiseMin = 2 * PLOTS.calibrated->GetBinContent(noiseBin);
	while (noiseCount < 2 * noiseMin) {
		noiseCount = PLOTS.calibrated->GetBinContent(noiseBin);
		if (noiseCount < noiseMin) {noiseMin = noiseCount;}
		noiseBin--;
		if (noiseBin == 0) {break;}
	}
	RESULTS.noiseWall = PLOTS.calibrated->GetXaxis()->GetBinCenter(noiseBin);
	cout << "Noise Wall at: " << RESULTS.noiseWall << " keV" << endl;

	//////// BEGIN OPTION HANDLING ////////

	if (option == "cal") {	
		
		TCanvas *calCanvas = new TCanvas("Calibration Canvas", "Calibration Canvas", 1);
		PLOTS.calPlot->Draw("ALP");
		PLOTS.calPlot->SetTitle("Calibration Curve");
		PLOTS.calPlot->GetYaxis()->SetTitle("ADC Energy");
		PLOTS.calPlot->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		//PLOTS.calPlot->BuildLegend(0.15,0.6,0.30,0.85); 	// legend in top left
		
	} else if (option == "sig") {
		
		TCanvas *sigmaCanvas = new TCanvas("Sigma Canvas", "Sigma Canvas", 1);
		
		PLOTS.sigmaPlot->Draw("ALP");
		PLOTS.sigmaPlot->SetTitle("Resolution vs Energy");
		PLOTS.sigmaPlot->GetYaxis()->SetTitle("Fractional Sigma (sigma/mean)");
		PLOTS.sigmaPlot->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		//PLOTS.sigmaPlot->BuildLegend(0.7,0.6,0.85,0.85); 	// legend in top right

	} else if (option == "res") {

		TCanvas *resCanvas = new TCanvas("Residue Canvas", "Residue Canvas", 1);
		Double_t expE[3] = {2615.0, 1460.0, 661.6};
		Double_t residues[3];
		
		Double_t energy = (PEAKS.Tl.mu - RESULTS.calPars.offset) / RESULTS.calPars.slope;
		residues[0] = 100 * (energy - expE[0]) / expE[0];
		energy = (PEAKS.K.mu - RESULTS.calPars.offset) / RESULTS.calPars.slope;
		residues[1] = 100 * (energy - expE[1]) / expE[1];		
		energy = (PEAKS.Cs.mu - RESULTS.calPars.offset) / RESULTS.calPars.slope;
		residues[2] = 100 * (energy - expE[2]) / expE[2];
		
		TGraphErrors *residuePlot = new TGraphErrors(3, expE, residues, 0, 0);
		residuePlot->SetLineColor(1);
		residuePlot->SetMarkerColor(1);
		residuePlot->SetMarkerStyle(21);
		residuePlot->SetLineWidth(1);
		
		residuePlot->Draw("ALP");
		residuePlot->SetTitle("Residues");
		residuePlot->GetXaxis()->SetTitle("ADC Energies");
		residuePlot->GetXaxis()->SetTitleOffset(1.3);
		residuePlot->GetYaxis()->SetTitle("\% Error in calibrated energy ((Fit - Exp)/Exp)");
		residuePlot->GetYaxis()->SetTitleOffset(1.2);
		//resCanvas->BuildLegend(0.15,0.15,0.30,0.40); 	// legend in bottom left

	} else if (option == "over") {
				
		TCanvas *overlayCanvas = new TCanvas("overlayCanvas", "Overlay Canvas", 1);
		gPad->SetLogy();
		PLOTS.calibrated->Draw();
		//overlayCanvas->BuildLegend(0.7,0.6,0.85,0.85); 	// legend in top right

	} else if (option == "backFits") {
		
		TCanvas *backCanvas = new TCanvas("backCanvas", "Background Fits", 1000, 400);
		backCanvas->Divide(2, 2);
		backCanvas->cd(1);
		PLOTS.backPlots.Tl->Draw();
		backCanvas->cd(2);
		PLOTS.backPlots.K->Draw();
		backCanvas->cd(3);
		PLOTS.backPlots.Cs->Draw();

	} else if (option == "AE") {
		
		TCanvas *AECanvas = new TCanvas("AECanvas", "A/E Plot", 1);
		TH2D *AEHist = new TH2D("AEHist", "Amplitude / Energy vs Energy (calibrated)",
		                        1e3, 0, 3e3, 1e3, 0, 20);

		string calE = "(energy - " + to_string(RESULTS.calPars.offset);
		calE += ") / " + to_string(RESULTS.calPars.slope);
		string toPlot = "amp / (" + calE + ") : (" + calE + ") >> AEHist";
		data->Draw(toPlot.c_str(), CHANNEL.c_str(), "COLZ");
		
		AEHist->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		AEHist->GetYaxis()->SetTitle("Amplitude / Callibrated Energy");

	} else if (option == "muon") {
		// muons at about 10x Tl peak position.
		// if it fails, get results and set to unphysical estimates;
		
		/*TCanvas *muonCanvas = new TCanvas("muonCanvas", "Muon Canvas", 1200, 800);
		muonCanvas->Divide(2,3);
		for (Int_t i = 0; i < NUMFILES; i++) {
			muonCanvas->cd(i+1);
			gPad->SetLogy();
			
			PLOTS[i].muRaw->GetXaxis()->SetRangeUser(0.6 * PEAKS[i].muon.mu, 
			                                         1.6 * PEAKS[i].muon.mu);
			PLOTS[i].muRaw->Draw();
			
			string title = "Muon Peak (uncalibrated) for " + getDepVar(mode);
			title += " " + to_string(i + 1);
			PLOTS[i].muRaw->SetTitle(title.c_str());
		}*/		
	}
	
}

Double_t snapToMax(TH1D *h, FitWindow win) {
	h->GetXaxis()->SetRangeUser(win.low, win.high);
	Double_t pos = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
	h->GetXaxis()->SetRange(0, NUMBINS);
	return pos;
}

FitPars backEst(TH1D *h, FitWindow win, Double_t coverage, string func, string peak) {
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
	
	string title = "Background Estimation for " + peak + " peak";
	backGraph->SetTitle((title).c_str());
	backGraph->GetYaxis()->SetTitle("Count");
	backGraph->GetXaxis()->SetTitle("Uncalibrated Energy");
	backGraph->SetMarkerStyle(4);
	backGraph->SetMarkerSize(0.5);

	if (peak == "Cs") {
		PLOTS.backPlots.Cs = backGraph;
	} else if (peak == "K") {
		PLOTS.backPlots.K = backGraph;
	} else if (peak == "Tl") {
		PLOTS.backPlots.Tl = backGraph;
	}

	return pars;
}
