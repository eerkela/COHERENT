#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TVirtualFitter.h>
#include <stdio.h>

// comment


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
	Double_t p383;
	Double_t p356;
	Double_t p302;
	Double_t p276;
	Double_t p81;
	Double_t p31;
};

struct BackgroundPlots {
	TGraph *p383;
	TGraph *p356;
	TGraph *p302;
	TGraph *p276;
	TGraph *p81;
	TGraph *p31;
	
};

//////// END INTERNAL STRUCTS ////////

struct FitWindow {
	Double_t low;
	Double_t high;
};

struct Peaks {
	PeakInfo p383;
	PeakInfo p356;
	PeakInfo p302;
	PeakInfo p276;
	PeakInfo p81;
	PeakInfo p31;
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
};

Int_t NUMFILES = 5;

FitPars backEst(TH1D *h, FitWindow win, Double_t coverage, string func, string peak, Int_t i);
void formatGraph(TGraphErrors *graph);
void getData(TChain *data[NUMFILES], string mode);
string getDepVar(string mode);
string getLabels(string mode, Int_t index);
Double_t getRunTime(TChain *c, string mode);
Double_t snapToMax(TH1D *h, FitWindow win);

Int_t NUMBINS = 30000;
Int_t VOLTAGES[5] = {750, 810, 870, 930, 990};
Int_t PEAKINGTIMES[5] = {50, 100, 200, 400, 800};
Int_t POSITIONS[5] = {1, 2, 3, 4, 5};

Plots PLOTS[5];
Peaks PEAKS[5];
Results RESULTS[5];

void HighGainCalibration(string mode, string option) {

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
		
		TH1D *hTemp = new TH1D("hTemp", "Finding 356 pos", NUMBINS, 0, overflowBin);
		data[i]->Draw("energy >> hTemp", "channel==2");
		
		// Fitting to complicated 356 keV region:
		cout << "Finding 356 keV Region..." << endl;
		UInt_t binNum = NUMBINS;
		UInt_t totalCount = 0;
		while (totalCount < (800.0 * time)) {
			totalCount += (UInt_t)hTemp->GetBinContent(binNum);
			binNum--;
			if (binNum == 0) {break;} 
		}
		Double_t pos356 = hTemp->GetXaxis()->GetBinCenter(binNum);
		cout << "Estimated 356 keV peak position: " << pos356 << endl;

		FitWindow snapWindow356;
		snapWindow356.low = 0.95 * pos356;
		snapWindow356.high = 1.05 * pos356;
		pos356 = snapToMax(hTemp, snapWindow356);
		Double_t pos356Normalized = (Double_t) hTemp->FindBin(pos356) / (Double_t) NUMBINS;

		cout << "Snapped 356 keV position: " << pos356 << endl;
		
		// Redraw histogram with constant number of bins below 356 peak to enhance fits.
		delete hTemp;
		NUMBINS = (Int_t) 500 / pos356Normalized;
		string hName = "h" + to_string(i+1);
		PLOTS[i].raw = new TH1D(hName.c_str(), "Uncalibrated Spectrum", NUMBINS, 0, overflowBin);
		data[i]->Draw(("energy >> " + hName).c_str(), "channel==2", "");
		
		snapWindow356.low = 0.95 * pos356;
		snapWindow356.high = 1.05 * pos356;
		pos356 = snapToMax(PLOTS[i].raw, snapWindow356);
		Double_t count356 = PLOTS[i].raw->GetBinContent(PLOTS[i].raw->FindBin(pos356));

		FitWindow fitWindow356;
		fitWindow356.low = 0.96 * pos356;
		fitWindow356.high = 1.22 * pos356;

		FitPars backPars356 = backEst(PLOTS[i].raw, fitWindow356, 0.2, "pol1", "356", i);
		
		string fitFunc356 = "[0]*exp(-0.5*((x-[1])/[2])^2)";
		fitFunc356 += " + [3]*exp(-0.5*((x-[4])/[2])^2)";
		fitFunc356 += " + [5]*exp(-0.5*((x-[6])/[2])^2)";
		fitFunc356 += " + [7]*exp(-0.5*((x-[8])/[2])^2)";
		fitFunc356 += " + [9] + [10]*x";

		Double_t parameters356[14];
		parameters356[0] = 0.2 * count356;
		parameters356[1] = 1.08 * pos356;
		parameters356[2] = 0.05 * pos356;
		parameters356[3] = count356;
		parameters356[4] = pos356;
		parameters356[5] = 0.35 * count356;
		parameters356[6] = 0.85 * pos356;
		parameters356[7] = 0.14 * count356;
		parameters356[8] = 0.8 * pos356;
		parameters356[9] = backPars356.offset;
		parameters356[10] = backPars356.slope;
		
		TF1 *fit356 = new TF1("fit356", fitFunc356.c_str(), fitWindow356.low, fitWindow356.high);
		fit356->SetParameters(parameters356);
		PLOTS[i].raw->Fit("fit356", "R+ll");
		
		PEAKS[i].p383.mu = fit356->GetParameter(1);
		PEAKS[i].p383.muErr = fit356->GetParError(1);
		PEAKS[i].p383.sigma = fit356->GetParameter(2);
		PEAKS[i].p383.sigmaErr = fit356->GetParError(2);
		PEAKS[i].p383.count = fit356->GetParameter(0);
		PEAKS[i].p356.mu = fit356->GetParameter(4);
		PEAKS[i].p356.muErr = fit356->GetParError(4);
		PEAKS[i].p356.sigma = fit356->GetParameter(2);
		PEAKS[i].p356.sigmaErr = fit356->GetParError(2);
		PEAKS[i].p356.count = fit356->GetParameter(3);	
		PEAKS[i].p302.mu = fit356->GetParameter(6);
		PEAKS[i].p302.muErr = fit356->GetParError(6);
		PEAKS[i].p302.sigma = fit356->GetParameter(2);
		PEAKS[i].p302.sigmaErr = fit356->GetParError(2);
		PEAKS[i].p302.count = fit356->GetParameter(5);
		PEAKS[i].p276.mu = fit356->GetParameter(8);
		PEAKS[i].p276.muErr = fit356->GetParError(8);
		PEAKS[i].p276.sigma = fit356->GetParameter(2);
		PEAKS[i].p276.sigmaErr = fit356->GetParError(2);
		PEAKS[i].p276.count = fit356->GetParameter(7);	

		// Fitting 80 keV peaks:
		cout << "Fitting 80 keV Region..." << endl;
		Double_t pos81 = (80.0/356.0) * PEAKS[i].p356.mu;
		
		FitWindow snapWindow81;
		snapWindow81.low = 0.95 * pos81;
		snapWindow81.high = 1.05 * pos81;

		pos81 = snapToMax(PLOTS[i].raw, snapWindow81);
		Double_t count81 = PLOTS[i].raw->GetBinContent(PLOTS[i].raw->FindBin(pos81));
		
		FitWindow fitWindow81;
		fitWindow81.low = 0.725 * pos81;
		fitWindow81.high = 1.25 * pos81;
		
		FitPars backPars81 = backEst(PLOTS[i].raw, fitWindow81, 0.2, "pol1", "81", i);
		
		string fitFunc81 = "[0]*exp(-0.5*((x-[1])/[2])^2)";
		fitFunc81 += " + [3]*exp(-0.5*((x-[4])/[2])^2)";
		fitFunc81 += " + [5] + [6] * x";

		Double_t parameters81[8];
		parameters81[0] = 0.8 * count81;
		parameters81[1] = pos81;
		parameters81[2] = 0.05 * pos81;
		parameters81[3] = 0.2 * count81;
		parameters81[4] = 0.825 * pos81;
		parameters81[5] = backPars81.offset;
		parameters81[6] = backPars81.slope;
		
		TF1 *fit81 = new TF1("fit81", fitFunc81.c_str(), fitWindow81.low, fitWindow81.high);
		fit81->SetParameters(parameters81);
		PLOTS[i].raw->Fit("fit81", "R+ll");
		
		PEAKS[i].p81.mu = fit81->GetParameter(1);
		PEAKS[i].p81.muErr = fit81->GetParError(1);
		PEAKS[i].p81.sigma = fit81->GetParameter(2);
		PEAKS[i].p81.sigmaErr = fit81->GetParError(2);
		PEAKS[i].p81.count = fit81->GetParameter(0);

		// Fitting 30 keV region:
		cout << "Fitting 31 keV Region" << endl;
		Double_t pos31 = (31.0/356.0) * PEAKS[i].p356.mu;
		
		FitWindow snapWindow31;
		snapWindow31.low = 0.5 * pos31;
		snapWindow31.high = 1.4 * pos31;
	
		pos31 = snapToMax(PLOTS[i].raw, snapWindow31);
		Double_t count31 = PLOTS[i].raw->GetBinContent(PLOTS[i].raw->FindBin(pos31));
		
		FitWindow fitWindow31;
		fitWindow31.low = 0.8 * pos31;
		fitWindow31.high = 1.2 * pos31;

		FitPars backPars31 = backEst(PLOTS[i].raw, fitWindow31, 0.2, "pol1", "31", i);

		string fitFunc31 = "[0]*exp(-0.5*((x-[1])/[2])^2)";
		fitFunc31 += " + [3]*exp(-0.5*((x-[4])/[2])^2)";
		fitFunc31 += " + [5] + [6]*x";

		Double_t parameters31[8];
		parameters31[0] = 0.07 * count31;
		parameters31[1] = 1.2 * pos31;
		parameters31[2] = 0.1 * pos31;
		parameters31[3] = 0.84 * count31;
		parameters31[4] = pos31;
		parameters31[5] = backPars31.offset;
		parameters31[6] = backPars31.slope;

		TF1 *fit31 = new TF1("fit31", fitFunc31.c_str(), fitWindow31.low, fitWindow31.high);
		fit31->SetParameters(parameters31);
		PLOTS[i].raw->Fit("fit31", "R+ll");
		
		PEAKS[i].p31.mu = fit31->GetParameter(1);
		PEAKS[i].p31.muErr = fit31->GetParError(1);
		PEAKS[i].p31.sigma = fit31->GetParameter(2);
		PEAKS[i].p31.sigmaErr = fit31->GetParError(2);
		PEAKS[i].p31.count = fit31->GetParameter(0);

		Float_t range = PEAKS[i].p383.mu + 0.15 * PEAKS[i].p383.mu;
		PLOTS[i].raw->GetXaxis()->SetRangeUser(0, range);
		NUMBINS = 30000;

		// Collecting data:
		Double_t expE[6];
			expE[0] = 383.8485;
			expE[1] = 356.0129;
			expE[2] = 302.8508;
			expE[3] = 276.3989;
			expE[4] = 80.9979;
			expE[5] = 30.973;
		Double_t fitE[6]; 
			fitE[0] = PEAKS[i].p383.mu;
			fitE[1] = PEAKS[i].p356.mu;
			fitE[2] = PEAKS[i].p302.mu;
			fitE[3] = PEAKS[i].p276.mu;
			fitE[4] = PEAKS[i].p81.mu;
			fitE[5] = PEAKS[i].p31.mu;
		Double_t fitEErr[6];
			fitEErr[0] = PEAKS[i].p383.muErr;
			fitEErr[1] = PEAKS[i].p356.muErr;
			fitEErr[2] = PEAKS[i].p302.muErr;
			fitEErr[3] = PEAKS[i].p276.muErr;
			fitEErr[4] = PEAKS[i].p81.muErr;
			fitEErr[5] = PEAKS[i].p31.muErr;
		Double_t fitSigma[6];
			fitSigma[0] = PEAKS[i].p383.sigma / PEAKS[i].p383.mu;
			fitSigma[1] = PEAKS[i].p356.sigma / PEAKS[i].p356.mu;
			fitSigma[2] = PEAKS[i].p302.sigma / PEAKS[i].p302.mu;
			fitSigma[3] = PEAKS[i].p276.sigma / PEAKS[i].p276.mu;
			fitSigma[4] = PEAKS[i].p81.sigma / PEAKS[i].p81.mu;
			fitSigma[5] = PEAKS[i].p31.sigma / PEAKS[i].p31.mu;
		Double_t fitSigmaErr[6];
			fitSigmaErr[0] = PEAKS[i].p383.sigmaErr / PEAKS[i].p383.mu;
			fitSigmaErr[1] = PEAKS[i].p356.sigmaErr / PEAKS[i].p356.mu;
			fitSigmaErr[2] = PEAKS[i].p302.sigmaErr / PEAKS[i].p302.mu;
			fitSigmaErr[3] = PEAKS[i].p276.sigmaErr / PEAKS[i].p276.mu;
			fitSigmaErr[4] = PEAKS[i].p81.sigmaErr / PEAKS[i].p81.mu;
			fitSigmaErr[5] = PEAKS[i].p31.sigmaErr / PEAKS[i].p31.mu;
		
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
		data[i]->Draw(calibration.c_str(), "channel==2", "goff");

		Int_t noiseCount;
		Int_t noiseBin = PLOTS[i].calibrated->FindBin(10);
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
			FIRResolution[k] = pow(PEAKS[k].p356.sigma / PEAKS[k].p356.mu, 2.0);
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
		
		Double_t resolution = PEAKS[3].p356.sigma / PEAKS[3].p356.mu;
		cout << "Normalized detector resolution (width of 356 keV peak" << endl;
		cout << " / mean at 3rd position):" << resolution << endl;
		
		Double_t pos356[NUMFILES];
		Double_t pos356Errs[NUMFILES];
		Double_t sig356[NUMFILES];
		Double_t sig356Errs[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {
			pos356[k] = PEAKS[k].p356.mu;
			pos356Errs[k] = PEAKS[k].p356.muErr;
			sig356[k] = PEAKS[k].p356.sigma / PEAKS[k].p356.mu;
			sig356Errs[k] = PEAKS[k].p356.sigmaErr / PEAKS[k].p356.mu;
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
		
		TCanvas *pos356Canvas = new TCanvas("pos356Canvas", "pos356Canvas");
		TGraphErrors *pos356Graph = new TGraphErrors(5, xAxis, pos356, 0, pos356Errs);
		pos356Graph->SetTitle(("356 Peak Energy vs " + getDepVar(mode)).c_str());
		pos356Graph->GetXaxis()->SetTitle(getDepVar(mode).c_str());
		pos356Graph->GetYaxis()->SetTitle("Uncalibrated 356 Peak Energy");
		
		formatGraph(pos356Graph);
		pos356Graph->Draw();
		pos356Canvas->Print("356EvsPos.pdf[");
		pos356Canvas->Print("356EvsPos.pdf");
		pos356Canvas->Print("356EvsPos.pdf]");
						
		TCanvas *res356Canvas = new TCanvas("res356Canvas", "res356Canvas");
		TGraphErrors *res356Graph = new TGraphErrors(5, xAxis, sig356, 0, sig356Errs);
		res356Graph->SetTitle(("356 Peak Width vs " + getDepVar(mode)).c_str());
		res356Graph->GetXaxis()->SetTitle(getDepVar(mode).c_str());
		res356Graph->GetYaxis()->SetTitle("Normalized 356 Peak Width (356 sigma / 356 mean)");
		
		formatGraph(res356Graph);
		res356Graph->Draw();
		res356Canvas->Print("356ResvsPos.pdf[");
		res356Canvas->Print("356ResvsPos.pdf");
		res356Canvas->Print("356ResvsPos.pdf]");
		
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

		TCanvas *resCanvas = new TCanvas("resCanvas", "Residue Canvas", 1);
		TMultiGraph *resComp = new TMultiGraph("resComp", "resComp");
		Double_t expE[6];
			expE[0] = 383.8485;
			expE[1] = 356.0129;
			expE[2] = 302.8508;
			expE[3] = 276.3989;
			expE[4] = 80.9979;
			expE[5] = 30.973;
		for (Int_t k = 0; k < NUMFILES; k++) {
			Double_t residues[6];
			
			Double_t energy = (PEAKS[k].p383.mu - RESULTS[k].calPars.offset)
			                   / RESULTS[k].calPars.slope;
			residues[0] = 100 * (energy - expE[0]) / expE[0];
				
			energy = (PEAKS[k].p356.mu - RESULTS[k].calPars.offset) 
			          / RESULTS[k].calPars.slope;
			residues[1] = 100 * (energy - expE[1]) / expE[1];			
			
			energy = (PEAKS[k].p302.mu - RESULTS[k].calPars.offset) 
			          / RESULTS[k].calPars.slope;
			residues[2] = 100 * (energy - expE[2]) / expE[2];

			energy = (PEAKS[k].p276.mu - RESULTS[k].calPars.offset) 
			          / RESULTS[k].calPars.slope;
			residues[3] = 100 * (energy - expE[3]) / expE[3];

			energy = (PEAKS[k].p81.mu - RESULTS[k].calPars.offset) 
			          / RESULTS[k].calPars.slope;
			residues[4] = 100 * (energy - expE[4]) / expE[4];
			
			energy = (PEAKS[k].p31.mu - RESULTS[k].calPars.offset) 
			          / RESULTS[k].calPars.slope;
			residues[5] = 100 * (energy - expE[5]) / expE[5];
			
			TGraphErrors *residuePlot = new TGraphErrors(6, expE, residues, 0, 0);
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
			PLOTS[i].backPlots.p356->Draw();
			backCanvas->cd(3*i+2);
			PLOTS[i].backPlots.p81->Draw();
			backCanvas->cd(3*i+3);
			PLOTS[i].backPlots.p31->Draw();
		}

	} else if (option == "AE") {

		TCanvas *AECanvas = new TCanvas("AECanvas", "A/E Plot", 1);
		TH2D *AEHist = new TH2D("AEHist", "Amplitude / Energy vs Energy (calibrated)",
		                        1e3, 0, 3e3, 1e3, 0, 10);

		string calE = "(energy - " + to_string(RESULTS[NUMFILES - 1].calPars.offset);
		calE += ") / " + to_string(RESULTS[NUMFILES - 1].calPars.slope);
		string toPlot = "amp / (" + calE + ") : (" + calE + ") >> AEHist";
		data[NUMFILES - 1]->Draw(toPlot.c_str(), "channel==2", "COLZ");
		
		AEHist->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		AEHist->GetYaxis()->SetTitle("Amplitude / Callibrated Energy");

	}
	
}

// This function will populate the chain array with the appropriate kind of data for the specified
// mode.  For instance, if mode = "pos", it will grab position data.
void GetData(TChain *DataChain[NUMFILES], string mode) {
	string path[NUMFILES];
	if (mode.compare("pos") == 0) {
		cout << "Finding Position Data..." << endl;
		path[0] = "Position/1/NaI_ET_run*";
		path[1] = "Position/2/NaI_ET_run*";
		path[2] = "Position/3/NaI_ET_run*";
		path[3] = "Position/4/NaI_ET_run*";
		path[4] = "Position/5/NaI_ET_run*";
	} else if (mode.compare("volt") == 0) {
		cout << "Finding Voltage Data..." << endl;
		path[0] = "Voltage/600V/NaI_ET_run*";
		path[1] = "Voltage/700V/NaI_ET_run*";
		path[2] = "Voltage/800V/NaI_ET_run*";
		path[3] = "Voltage/900V/NaI_ET_run*";
		path[4] = "Voltage/1000V/NaI_ET_run*";
	} else if (mode.compare("fir") == 0) {
		cout << "Finding FIR Data..." << endl;
		path[0] = "FIR/1/NaI_ET_run*";
		path[1] = "FIR/2/NaI_ET_run*";
		path[2] = "FIR/3/NaI_ET_run*";
		path[3] = "FIR/4/NaI_ET_run*";
		path[4] = "FIR/5/NaI_ET_run*";
	}
	cout << "Forming data chains..." << endl;
	for (Int_t i = 0; i < NUMFILES; i++) {
		DataChain[i] = new TChain("st");
		DataChain[i]->Add(path[i].c_str());
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
	
	string title = "Background Estimation " + to_string(i) + " for " + peak + " keV peak";
	backGraph->SetTitle((title).c_str());
	backGraph->GetYaxis()->SetTitle("Count");
	backGraph->GetXaxis()->SetTitle("Uncalibrated Energy");
	backGraph->SetMarkerStyle(4);
	backGraph->SetMarkerSize(0.5);

	if (peak == "356") {
		PLOTS[i].backPlots.p356 = backGraph;
	} else if (peak == "81") {
		PLOTS[i].backPlots.p81 = backGraph;
	} else if (peak == "31") {
		PLOTS[i].backPlots.p31 = backGraph;
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
