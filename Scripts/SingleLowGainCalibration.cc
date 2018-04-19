#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <stdio.h>

// This root script will process a series of sequential NaI .root files containing raw low gain
// data and calibrate them to a real energy scale.  Some useful options are available for 
// displaying relevant information.

// Parameters:
// NUMFILES 	indicates the number of files to be processed.  They must be sequential and carry
//				the "NaI_ET_run####.root" naming scheme.

// Options:
// "cal" 	will display the graphs used to generate the calibrations for each crystal.
// "res"	will display residues from the calibration algorithm, i.e. the % error by which the
//				calibrated energies deviate from expected. (NYI)
// "sig" 	will display fractional sigmas for each peak in order to measure detector
//				resolution.
// (else)	will simply calibrate the spectra then overlay the histograms with real energies.

using namespace std;

typedef struct {
	Double_t Offset;
	Double_t Slope;
} FitPars;

typedef struct {
	Double_t Low;
	Double_t High;
} FitWindow;

Int_t NUMBINS = 30000;

Double_t SnapToMax(TH1D *h, Double_t pos);
FitPars BackEst(TH1D *h, Double_t pos, FitWindow win, Int_t bins, string func);

void SingleLowGainCalibration(string option) {
	
	Double_t X[3];
	Double_t XErr[3];
	Double_t Sigma[3];
	Double_t SigmaErr[3];
	Double_t Count[3];
	
	TCanvas *c1 = new TCanvas("c1", "Peak Fitting", 1);
	TChain *chain = new TChain("st");
	chain->Add("Voltage/600V/NaI_ET_run*");

	Int_t NTrees = chain->GetNtrees();
	Double_t time = NTrees * 600;

	cout << endl;
	cout << "run time in seconds: " << time << endl;
	cout << endl;
	
	TH1D *h = new TH1D("h", "High Gain Uncalibrated Spectrum", NUMBINS, 0, 5100e3);
	h->GetXaxis()->SetTitle("ADC Energy");
	h->GetYaxis()->SetTitle("Count");
	chain->Draw("energy >> h", "", "");

	// Finding Thalium (Tl):
	Int_t binNum = NUMBINS; 
	Double_t binCount = 0;
	while (binCount < (13.5 * time)) {
		binCount += (Double_t)h->GetBinContent(binNum);
		binNum--;
	}
	X[0] = h->GetXaxis()->GetBinCenter(binNum);
	X[0] = SnapToMax(h, X[0]);
	Count[0] = h->GetBinContent(h->FindBin(X[0]));

	cout << "Triggered Tl peak location: " << X[0] << ", Count: ";
	cout << ", Count: " << Count[0] << endl;
	
	FitWindow TlWindow;
	TlWindow.Low = X[0] - 0.08 * X[0];
	TlWindow.High = X[0] + 0.08 * X[0];

	FitPars TlBackPars = BackEst(h, X[0], TlWindow, 500, "expo");

	string TlFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2)";
	//TlFitFunc += " + [3] + [4] * x";
	TlFitFunc += " + exp([3] + [4] * x)";

	Double_t TlParameters[5];
	TlParameters[0] = Count[0];
	TlParameters[1] = X[0];
	TlParameters[2] = 0.05 * X[0];
	TlParameters[3] = TlBackPars.Offset;
	TlParameters[4] = TlBackPars.Slope;

	TF1 *TlFit = new TF1("TlFit", TlFitFunc.c_str(), TlWindow.Low, TlWindow.High);
	TlFit->SetParameters(TlParameters);
	h->Fit("TlFit", "R+");
	
	X[0] = TlFit->GetParameter(1);
	XErr[0] = TlFit->GetParError(1);
	Sigma[0] = abs(1.0/X[0] * TlFit->GetParameter(2));
	SigmaErr[0] = 1.0/X[0] * TlFit->GetParError(2);	

	cout << "Fitted Tl peak location: " << X[0] << endl;

	// Finding Potassium (K):
	X[1] = (1460.0/2615.0) * X[0]; // we guess where K is based on Tl location.
	cout << "Interpolated location of K peak: " << X[1] << endl;

	X[1] = SnapToMax(h, X[1]);
	Count[1] = h->GetBinContent(h->FindBin(X[1]));
	cout << "Local max at: " << X[1] << ", peak count: " << Count[1] << endl;

	FitWindow KWindow;
	KWindow.Low = X[1] - 0.1 * X[1];
	KWindow.High = X[1] + 0.1 * X[1];

	FitPars KBackPars = BackEst(h, X[1], KWindow, 500, "expo");

	string KFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2)";
	//KFitFunc += "+ [3] + [4] * x";
	KFitFunc += " + exp([3] + [4] * x)";

	Double_t KParameters[5];
	KParameters[0] = Count[1];
	KParameters[1] = X[1];
	KParameters[2] = 0.05 * X[1];
	KParameters[3] = KBackPars.Offset;
	KParameters[4] = KBackPars.Slope;

	TF1 *KFit = new TF1("KFit", KFitFunc.c_str(), KWindow.Low, KWindow.High);
	KFit->SetParameters(KParameters);
	h->Fit("KFit", "R+");
	X[1] = KFit->GetParameter(1);
	XErr[1] = KFit->GetParError(1);
	Sigma[1] = abs(1.0/X[1] * KFit->GetParameter(2));
	SigmaErr[1] = 1.0/X[1] * KFit->GetParError(2);

	cout << "Fitted K peak location: " << X[1] << endl;
	cout << endl;

	// Finding Cesium (Cs):
	X[2] = (661.6/2615.0) * X[0];
	cout << "Interpolated location of Cs peak: " << X[2] << endl;

	X[2] = SnapToMax(h, X[2]);
	Count[2] = h->GetBinContent(h->FindBin(X[2]));
	cout << "Local max at: " << X[2] << ", peak count: " << Count[2] << endl;

	FitWindow CsWindow;
	CsWindow.Low = X[2] - 0.1 * X[2];
	CsWindow.High = X[2] + 0.12 * X[2];

	FitPars CsBackPars = BackEst(h, X[2], CsWindow, 500, "expo");

	string CsFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2)";
	//CsFitFunc += "+ [3] + [4] * x";
	CsFitFunc += " + exp([3] + [4] * x)";

	Double_t CsParameters[5];
	CsParameters[0] = Count[2];
	CsParameters[1] = X[2];
	CsParameters[2] = 0.05 * X[2];
	CsParameters[3] = CsBackPars.Offset;
	CsParameters[4] = CsBackPars.Slope;

	TF1 *CsFit = new TF1("CsFit", CsFitFunc.c_str(), CsWindow.Low, CsWindow.High);
	CsFit->SetParameters(CsParameters);
	h->Fit("CsFit", "R+");
	X[2] = CsFit->GetParameter(1);
	XErr[2] = CsFit->GetParError(1);
	Sigma[2] = abs(1.0/X[2] * CsFit->GetParameter(2));
	SigmaErr[2] = 1.0/X[2] * CsFit->GetParError(2);

	cout << "Fitted Cs peak location: " << X[2] << endl;
	cout << endl;

	// Getting calibration parameters by fitting to a TGraph:
	Double_t expEnergy[3] = {2615.0, 1460.0, 661.6};

	TGraphErrors *calibrationPlot = new TGraphErrors(3, expEnergy, X, 0, XErr);
	TF1 *calibrationFit = new TF1("calibrationFit", "pol1", 0, 5000e3);	
	calibrationPlot->Fit("calibrationFit", "", "", 0, 0);
		
	// Linear calibration parameters:
	Double_t Slope = calibrationFit->GetParameter(1);
	Double_t SlopeError = calibrationFit->GetParError(1);
	Double_t Offset = calibrationFit->GetParameter(0);
	Double_t OffsetError = calibrationFit->GetParError(0);

	if (option.compare("cal") == 0) {
		
		c1->Clear();
		calibrationPlot->SetTitle("Expected Peak Energies vs. Found values");
		calibrationPlot->GetXaxis()->SetTitle("ADC Energies");
		calibrationPlot->GetYaxis()->SetTitle("Expected Energies (keV)");
		calibrationPlot->SetMarkerColor(4);
		calibrationPlot->SetMarkerStyle(21);
		calibrationPlot->SetLineColor(1);
		calibrationPlot->SetLineWidth(2);
		calibrationPlot->Draw("");
	
	} else if (option.compare("res") == 0) {

		c1->Clear();
		Double_t Residues[3];
		for (int k = 0; k < 3; k++) {
			Double_t fitEnergy = (X[k] - Offset) / Slope;
			Residues[k] = 100 * (fitEnergy - expEnergy[k]) / expEnergy[k];
		}
		TGraphErrors *residuePlot = new TGraphErrors(3, expEnergy, Residues, 0, 0);
		residuePlot->SetTitle("Calibration Residues vs Expected values");
		residuePlot->GetXaxis()->SetTitle("ADC Energies");
		residuePlot->GetXaxis()->SetTitleOffset(1.3);
		residuePlot->GetYaxis()->SetTitle("\% Error in calibrated energy ((Fit - Exp)/Exp)");
		residuePlot->GetYaxis()->SetTitleOffset(1.2);
		residuePlot->SetMarkerColor(4);
		residuePlot->SetMarkerStyle(21);
		residuePlot->SetLineColor(1);
		residuePlot->SetLineWidth(2);
		residuePlot->Draw("");

	} else if (option.compare("fit") == 0) {

	} else {

		TH1D *h2 = new TH1D("h2", "Calibrated Spectrum", 3000, 0, 3000);
		h2->GetYaxis()->SetTitle("Count");
		h2->GetXaxis()->SetTitle("Calibrated Energy (keV)");

		string calibration = "(energy - " + to_string(Offset) + ") / ";
		calibration += to_string(Slope) + " >> h2";
		chain->Draw(calibration.c_str(), "", "");

	}
	
}

// This function accepts a position (pos) and returns the location of the local maximum within
// 10% of that stated position.  This is used to find the actual maximum of an approximate
// peak position, enhancing fit stability.
Double_t SnapToMax(TH1D *h, Double_t pos) {
	h->GetXaxis()->SetRangeUser(pos - 0.1 * pos, pos + 0.1 * pos);
	pos = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
	h->GetXaxis()->SetRange(0, NUMBINS);
	return pos;
}

// This function will estimate the background underneath a peak by pulling a set number of points 
// from the edges of the fit region, pushing them into a TGraph, and fitting to a standard root
// function (defined by func) in order to obtain an estimate for the background at pos.
FitPars BackEst(TH1D *h, Double_t pos, FitWindow win, Int_t bins, string func) {
	Double_t NormPos = ((Double_t) h->GetBinCenter(pos))/((Double_t) NUMBINS);

	cout << "Bin / NUMBINS = " << NormPos << endl;

	Int_t BackBins = bins * NormPos;
	Int_t LoBin = h->FindBin(win.Low);
	Int_t HiBin = h->FindBin(win.High);
	Double_t BackVals[2 * BackBins];
	Double_t BackX[2 * BackBins];

	//
	for (Int_t j = 0; j < BackBins; j++) {
		BackVals[j] = h->GetBinContent(j + LoBin);
		BackX[j] = h->GetBinCenter(j + LoBin);
	}
	for (Int_t j = 0; j < BackBins; j++) {
		BackVals[j + BackBins] = h->GetBinContent(HiBin - j);
		BackX[j + BackBins] = h->GetBinCenter(HiBin - j);
	}

	TGraph *Back = new TGraph(2 * BackBins, BackX, BackVals);
	TF1 *BackFit = new TF1("BackFit", func.c_str(), win.Low, win.High);
	Back->Fit("BackFit");
	
	FitPars pars;
	pars.Offset = BackFit->GetParameter(0);
	pars.Slope = BackFit->GetParameter(1);
	return pars;
}


// Can print canvases to pdf books
// canvas->Print("filename.pdf[")  when you intitalize the canvas
// canvas->Print("filename.pdf")   this will make a new page in the document
// canvas->Print("filename.pdf]")  this will close the file
