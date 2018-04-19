#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <stdio.h>

/*
This root script will process a series of sequential NaI .root files containing raw low gain
data and calibrate them to a real energy scale.  Some useful options are available for 
displaying relevant information.


Modes:
"pos"	will execute calibration algorithm for position data.
"volt"	will execute calibration algorithm for voltage data.
"fir"	will perform FIR filter optimization.


Options:
"cal" 	will display the graphs used to generate the calibrations for each run.
"fit"	will display the fits found using the raw digitizer output for each run.  This is
		a very useful debugging tool.
"gain"	will display a graph of the calculated calibration slope (gain) as a function of the
		dependent variable specified by the mode.
"over"	will calibrate the spectra then overlay the calibrated histograms with real energies 
		for qualitative comparison.
"res"	will display residues from the calibration algorithm for each run. Residues are the % 
		error by which the calibrated energies deviate from expected.
"sig" 	will display fractional sigmas for each peak in each run. 


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
	are hardcoded in. When we get the getSpectrum script to pull the peaking time out of the
	ORCA headers, We can revise this.

	I've kept the muon code snippets in just for safe keeping, though none are implemented.

	I define a Calibrated histogram in the CalibrationPlots structure, but don't implement it
	until the "over" output option. If we stick with this behavior, then the calibrated
	histograms will be inaccessible unless that if statement is called.

	Can print canvases to pdf books:	
		canvas->Print("filename.pdf[")  when you intitalize the canvas
		canvas->Print("filename.pdf")   this will make a new page in the document
		canvas->Print("filename.pdf]")  this will close the file
*/

using namespace std;

typedef struct {
	Double_t Offset;
	Double_t OffsetError;
	Double_t Slope;
	Double_t SlopeError;
} FitPars;

typedef struct {
	Double_t Low;
	Double_t High;
} FitWindow;

typedef struct {
	Double_t X383;
	Double_t X383Err;
	Double_t Sigma383;
	Double_t Sigma383Err;
	Double_t Count383;
	
	Double_t X356;
	Double_t X356Err;
	Double_t Sigma356;
	Double_t Sigma356Err;
	Double_t Count356;
	
	Double_t X302;
	Double_t X302Err;
	Double_t Sigma302;
	Double_t Sigma302Err;
	Double_t Count302;

	Double_t X276;
	Double_t X276Err;
	Double_t Sigma276;
	Double_t Sigma276Err;
	Double_t Count276;

	Double_t X81;
	Double_t X81Err;
	Double_t Sigma81;
	Double_t Sigma81Err;
	Double_t Count81;

	Double_t X30;
	Double_t X30Err;
	Double_t Sigma30;
	Double_t Sigma30Err;
	Double_t Count30;
} PeakStruct;

typedef struct {
	Double_t Offset;
	Double_t OffsetError;
	Double_t Slope;
	Double_t SlopeError;

	Double_t NoiseWall;
} CalStruct;

typedef struct {
	TH1D *Raw;
	TH1D *Calibrated;
	TH1D *MuPlot;
	TGraphErrors *FitPlot;
	TGraphErrors *SigmaPlot;
} PlotStruct;

Int_t NUMBINS = 30000;
Int_t NUMFILES = 5;

FitPars BackEst(TH1D *h, Double_t pos, FitWindow win, Int_t bins, string func);
void FormatGraph(TGraphErrors *Graph);
Double_t SnapToMax(TH1D *h, Double_t pos);

void SingleHighGainCalibration(string option) {

	//gStyle->SetOptFit(1111);

	TCanvas *RawCanvas = new TCanvas("Raw Spectra", "Raw Spectra", 1200, 800);
	TChain *DataChain = new TChain("st");
	DataChain->Add("NaI_ET_run*");
	Double_t Time = (Double_t) DataChain->GetNtrees() * 600;

	//TF1 *muonFit = new TF1("muonFit", "[0]*exp(-0.5*((x-[1])/[2]+exp(-(x-[1])/[2])))+[3]*x+[4]",0,0);
	//Double_t MuX[NUMFILES];
	
	PlotStruct CalPlots;
	PeakStruct PeakInfo;
	CalStruct CalInfo;
	
	cout << "Beginning Calibration..." << endl;
	gPad->SetLogy();
	
	Double_t OverflowBinPos = DataChain->GetMaximum("energy");
	OverflowBinPos = OverflowBinPos + 0.01 * OverflowBinPos;
	
	TH1D *hTest = new TH1D("hTest", "Finding 356 pos", NUMBINS, 0, OverflowBinPos);
	DataChain->Draw("energy >> hTest", "channel==0", "");
	
	//string hMuName = "hMu" + std::to_string();
	//hMu[i] = CalPlots.Raw->Clone(hMuName.c_str());		
	//c[i]->Draw(("energy >> " + hMuName).c_str());
	
	// Fitting to complicated 356 keV region:
	cout << "Fitting 356 keV Region..." << endl;
	Int_t BinNum = NUMBINS;
	Int_t BinCount = 0;
	while (BinCount < (800.0 * Time)) {
		BinCount += (Int_t)hTest->GetBinContent(BinNum);
		BinNum--;
		if (BinNum == 0) {break;} 
	}
	Double_t X356 = hTest->GetXaxis()->GetBinCenter(BinNum);
	X356 = SnapToMax(hTest, X356);
	Double_t NormPos356 = (Double_t)hTest->FindBin(X356)/NUMBINS;
	cout << "NormPos356: " << NormPos356 << endl;

	// Redraw histogram with constant number of bins below 356 peak to enhance fits.
	delete hTest;
	string hName = "h1";
	CalPlots.Raw = new TH1D(hName.c_str(), "Uncalibrated Spectrum", (Int_t)(500.0/NormPos356), 0, OverflowBinPos);
	DataChain->Draw(("energy >> " + hName).c_str(), "channel==0", "");
	
	X356 = SnapToMax(CalPlots.Raw, X356);
	Double_t Count356 = CalPlots.Raw->GetBinContent(CalPlots.Raw->FindBin(X356));

	FitWindow Window356;		
	Window356.Low = X356 - 0.4 * X356;
	Window356.High = X356 + 0.25 * X356;

	FitPars BackPars356 = BackEst(CalPlots.Raw, X356, Window356, 500, "pol1");
	
	string FitFunc356 = "[0]*exp(-0.5*((x-[1])/[2])^2)";
	FitFunc356 += " + [3]*exp(-0.5*((x-[4])/[2])^2)";
	FitFunc356 += " + [5]*exp(-0.5*((x-[6])/[2])^2)";
	FitFunc356 += " + [7]*exp(-0.5*((x-[8])/[2])^2)";
	FitFunc356 += " + [9] + [10]*x";

	Double_t Parameters356[14];
	Parameters356[0] = 0.2 * Count356;
	Parameters356[1] = X356 + 0.08 * X356;
	Parameters356[2] = 0.05 * X356;
	Parameters356[3] = Count356;
	Parameters356[4] = X356;
	Parameters356[5] = 0.35 * Count356;
	Parameters356[6] = X356 - 0.15 * X356;
	Parameters356[7] = 0.14 * Count356;
	Parameters356[8] = X356 - 0.2 * X356;
	Parameters356[9] = BackPars356.Offset;
	Parameters356[10] = BackPars356.Slope;
	
	TF1 *Fit356 = new TF1("Fit356", FitFunc356.c_str(), Window356.Low, Window356.High);
	Fit356->SetParameters(Parameters356);
	CalPlots.Raw->Fit("Fit356", "R+ll");
	
	PeakInfo.X383 = Fit356->GetParameter(1);
	PeakInfo.X383Err = Fit356->GetParError(1);
	PeakInfo.Sigma383 = Fit356->GetParameter(2);
	PeakInfo.Sigma383Err = Fit356->GetParError(2);
	PeakInfo.Count383 = Fit356->GetParameter(0);
	PeakInfo.X356 = Fit356->GetParameter(4);
	PeakInfo.X356Err = Fit356->GetParError(4);
	PeakInfo.Sigma356 = Fit356->GetParameter(2);
	PeakInfo.Sigma356Err = Fit356->GetParError(2);
	PeakInfo.Count356 = Fit356->GetParameter(3);	
	PeakInfo.X302 = Fit356->GetParameter(6);
	PeakInfo.X302Err = Fit356->GetParError(6);
	PeakInfo.Sigma302 = Fit356->GetParameter(2);
	PeakInfo.Sigma302Err = Fit356->GetParError(2);
	PeakInfo.Count302 = Fit356->GetParameter(5);
	PeakInfo.X276 = Fit356->GetParameter(8);
	PeakInfo.X276Err = Fit356->GetParError(8);
	PeakInfo.Sigma276 = Fit356->GetParameter(2);
	PeakInfo.Sigma276Err = Fit356->GetParError(2);
	PeakInfo.Count276 = Fit356->GetParameter(7);	

	cout << endl;
	cout<<"383 keV peak X: "<<PeakInfo.X383<<", Count: "<<PeakInfo.Count383<<endl;
	cout<<"356 keV peak X: "<<PeakInfo.X356<<", Count: "<<PeakInfo.Count356<<endl;
	cout<<"302 keV peak X: "<<PeakInfo.X302<<", Count: "<<PeakInfo.Count302<<endl;
	cout<<"276 keV peak X: "<<PeakInfo.X276<<", Count: "<<PeakInfo.Count276<<endl;
	cout << endl;

	// Fitting 80 keV peaks:
	cout << "Fitting 80 keV Region..." << endl;
	Double_t X80 = (80.0/356.0) * PeakInfo.X356;
	X80 = SnapToMax(CalPlots.Raw, X80);
	Double_t Count80 = CalPlots.Raw->GetBinContent(CalPlots.Raw->FindBin(X80));
	
	FitWindow Window80;
	Window80.Low = X80 - 0.275 * X80;
	Window80.High = X80 + 0.25 * X80;
	
	FitPars BackPars80 = BackEst(CalPlots.Raw, X80, Window80, 500, "pol1");
	
	string FitFunc80 = "[0]*exp(-0.5*((x-[1])/[2])^2)";
	FitFunc80 += " + [3]*exp(-0.5*((x-[4])/[2])^2)";
	FitFunc80 += " + [5] + [6] * x";

	Double_t Parameters80[8];
	Parameters80[0] = 0.8 * Count80;
	Parameters80[1] = X80;
	Parameters80[2] = 0.05 * X80;
	Parameters80[3] = 0.2 * Count80;
	Parameters80[4] = X80 - 0.175 * X80;
	Parameters80[5] = BackPars80.Offset;
	Parameters80[6] = BackPars80.Slope;
	
	TF1 *Fit80 = new TF1("Fit80", FitFunc80.c_str(), Window80.Low, Window80.High);
	Fit80->SetParameters(Parameters80);
	CalPlots.Raw->Fit("Fit80", "R+ll");
	
	PeakInfo.X81 = Fit80->GetParameter(1);
	PeakInfo.X81Err = Fit80->GetParError(1);
	PeakInfo.Sigma81 = Fit80->GetParameter(2);
	PeakInfo.Sigma81Err = Fit80->GetParError(2);
	PeakInfo.Count81 = Fit80->GetParameter(0);

	cout << endl;
	cout<<"81.0 keV peak X: "<<PeakInfo.X81<<", Count: "<<PeakInfo.Count81<<endl;
	//cout<<"79.6 keV peak X: "<<PeakInfo.X79<<", Count: "<<PeakInfo.Count79<<endl;
	cout << endl;

	// Fitting 30 keV region:
	cout << "Fitting 30 keV Region" << endl;
	Double_t X30 = (32.0/356.0) * PeakInfo.X356;
	X30 = SnapToMax(CalPlots.Raw, X30);
	Double_t Count30 = CalPlots.Raw->GetBinContent(CalPlots.Raw->FindBin(X30));
	
	FitWindow Window30;
	Window30.Low = X30 - 0.5 * X30;
	Window30.High = X30 + 0.45 * X30;

	FitPars BackPars30 = BackEst(CalPlots.Raw, X30, Window30, 500, "pol1");

	string FitFunc30 = "[0]*exp(-0.5*((x-[1])/[2])^2)";
	FitFunc30 += " + [3]*exp(-0.5*((x-[4])/[2])^2)";
	FitFunc30 += " + [5] + [6]*x";

	Double_t Parameters30[8];
	Parameters30[0] = 0.07 * Count30;
	Parameters30[1] = X30 + 0.2 * X30;
	Parameters30[2] = 0.1 * X30;
	Parameters30[3] = 0.84 * Count30;
	Parameters30[4] = X30;
	Parameters30[5] = BackPars30.Offset;
	Parameters30[6] = BackPars30.Slope;

	TF1 *Fit30 = new TF1("Fit30", FitFunc30.c_str(), Window30.Low, Window30.High);
	Fit30->SetParameters(Parameters30);
	CalPlots.Raw->Fit("Fit30", "R+ll");
	
	PeakInfo.X30 = Fit30->GetParameter(1);
	PeakInfo.X30Err = Fit30->GetParError(1);
	PeakInfo.Sigma30 = Fit30->GetParameter(2);
	PeakInfo.Sigma30Err = Fit30->GetParError(2);
	PeakInfo.Count30 = Fit30->GetParameter(0);

	cout << endl;
	cout<<"30.97 keV peak X: "<<PeakInfo.X30<<", Count: "<<PeakInfo.Count30<<endl;
	cout << endl;

	// Collecting data:
	Double_t expE[6];
		expE[0] = 383.8485;
		expE[1] = 356.0129;
		expE[2] = 302.8508;
		expE[3] = 276.3989;
		expE[4] = 80.9979;
		expE[5] = 30.973;
	Double_t fitE[6]; 
		fitE[0] = PeakInfo.X383;
		fitE[1] = PeakInfo.X356;
		fitE[2] = PeakInfo.X302;
		fitE[3] = PeakInfo.X276;
		fitE[4] = PeakInfo.X81;
		fitE[5] = PeakInfo.X30;
	Double_t fitEErr[6];
		fitEErr[0] = PeakInfo.X383Err;
		fitEErr[1] = PeakInfo.X356Err;
		fitEErr[2] = PeakInfo.X302Err;
		fitEErr[3] = PeakInfo.X276Err;
		fitEErr[4] = PeakInfo.X81Err;
		fitEErr[5] = PeakInfo.X30Err;
	Double_t fitSigma[6];
		fitSigma[0] = PeakInfo.Sigma383 / PeakInfo.X383;
		fitSigma[1] = PeakInfo.Sigma356 / PeakInfo.X356;
		fitSigma[2] = PeakInfo.Sigma302 / PeakInfo.X302;
		fitSigma[3] = PeakInfo.Sigma276 / PeakInfo.X276;
		fitSigma[4] = PeakInfo.Sigma81 / PeakInfo.X81;
		fitSigma[5] = PeakInfo.Sigma30 / PeakInfo.X30;
	Double_t fitSigmaErr[6];
		fitSigmaErr[0] = PeakInfo.Sigma383Err / PeakInfo.X383;
		fitSigmaErr[1] = PeakInfo.Sigma356Err / PeakInfo.X356;
		fitSigmaErr[2] = PeakInfo.Sigma302Err / PeakInfo.X302;
		fitSigmaErr[3] = PeakInfo.Sigma276Err / PeakInfo.X276;
		fitSigmaErr[4] = PeakInfo.Sigma81Err / PeakInfo.X81;
		fitSigmaErr[5] = PeakInfo.Sigma30Err / PeakInfo.X30;
	
	CalPlots.FitPlot = new TGraphErrors(6, expE, fitE, 0, fitEErr);
	CalPlots.FitPlot->SetTitle("Fit Plot");
	FormatGraph(CalPlots.FitPlot);
	
	CalPlots.SigmaPlot = new TGraphErrors(3, expE, fitSigma, 0, fitSigmaErr); 
	CalPlots.SigmaPlot->SetTitle("Sigma Plot");
	FormatGraph(CalPlots.SigmaPlot);
	
	TF1 *CalibrationFit = new TF1("CalibrationFit", "pol1", 0, 4500e3);
	CalPlots.FitPlot->Fit("CalibrationFit", "", "", 0, 0);
	
	CalInfo.Slope = CalibrationFit->GetParameter(1);
	CalInfo.SlopeError = CalibrationFit->GetParError(1);
	CalInfo.Offset = CalibrationFit->GetParameter(0);
	CalInfo.OffsetError = CalibrationFit->GetParError(0);

	Int_t MaxCalBin = (OverflowBinPos - CalInfo.Offset) / CalInfo.Slope;
	CalPlots.Calibrated = new TH1D("Calibrated", "Calibrated Spectrum", 20e3, 0, MaxCalBin);
	CalPlots.Calibrated->SetLineColor(1);
	CalPlots.Calibrated->GetXaxis()->SetTitle("Calibrated Energy (keV)");
	CalPlots.Calibrated->GetYaxis()->SetTitle("Count");
	
	string calibration = "(energy - " + to_string(CalInfo.Offset);
	calibration += ") / " + to_string(CalInfo.Slope) ;
	calibration += " >> Calibrated";
	DataChain->Draw(calibration.c_str(), "", "goff");

	Int_t NoiseCount;
	Int_t NoiseBin = CalPlots.Calibrated->FindBin(10);
	Int_t NoiseMin = 2 * CalPlots.Calibrated->GetBinContent(NoiseBin);
	while(NoiseCount < 2 * NoiseMin) {
		NoiseCount = CalPlots.Calibrated->GetBinContent(NoiseBin);
		if (NoiseCount < NoiseMin) {NoiseMin = NoiseCount;}
		NoiseBin--;
		if (NoiseBin == 0) {break;}
	}
	CalInfo.NoiseWall = CalPlots.Calibrated->GetXaxis()->GetBinCenter(NoiseBin);
	cout << "Noise Wall at: " << CalInfo.NoiseWall << " keV" << endl;

	if (option.compare("cal") == 0) {	
		
		TCanvas *CalCanvas = new TCanvas("Calibration Canvas", "Calibration Canvas", 1);
		CalPlots.FitPlot->Draw();
		CalPlots.FitPlot->SetTitle("Calibration Curve");
		CalPlots.FitPlot->GetYaxis()->SetTitle("ADC Energy");
		CalPlots.FitPlot->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		CalCanvas->BuildLegend(0.15,0.6,0.30,0.85); 	// legend in top left
		
	} else if (option.compare("sig") == 0) {
		
		TCanvas *SigmaCanvas = new TCanvas("Sigma Canvas", "Sigma Canvas", 1);
		CalPlots.SigmaPlot->Draw();	
		CalPlots.SigmaPlot->SetTitle("Resolution");
		CalPlots.SigmaPlot->GetYaxis()->SetTitle("Fractional Sigma (sigma/mean)");
		CalPlots.SigmaPlot->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		SigmaCanvas->BuildLegend(0.7,0.6,0.85,0.85); 	// legend in top right

	} else if (option.compare("fit") == 0) {

		Double_t Range = PeakInfo.X383 + 0.5 * PeakInfo.X383;
		CalPlots.Raw->GetXaxis()->SetRangeUser(0, Range);
		
	} else if (option.compare("res") == 0) {

		TCanvas *ResCanvas = new TCanvas("Residue Canvas", "Residue Canvas", 1);
		Double_t expE[6];
			expE[0] = 383.8485;
			expE[1] = 356.0129;
			expE[2] = 302.8508;
			expE[3] = 276.3989;
			expE[4] = 80.9979;
			expE[5] = 30.973;
		Double_t Residues[6];
		
		Double_t E = (PeakInfo.X383 - CalInfo.Offset) / CalInfo.Slope;
		Residues[0] = 100 * (E - expE[0]) / expE[0];
		E = (PeakInfo.X356 - CalInfo.Offset) / CalInfo.Slope;
		Residues[1] = 100 * (E - expE[1]) / expE[1];			
		E = (PeakInfo.X302 - CalInfo.Offset) / CalInfo.Slope;
		Residues[2] = 100 * (E - expE[2]) / expE[2];
		E = (PeakInfo.X276 - CalInfo.Offset) / CalInfo.Slope;
		Residues[3] = 100 * (E - expE[3]) / expE[3];
		E = (PeakInfo.X81 - CalInfo.Offset) / CalInfo.Slope;
		Residues[4] = 100 * (E - expE[4]) / expE[4];
		E = (PeakInfo.X30 - CalInfo.Offset) / CalInfo.Slope;
		Residues[5] = 100 * (E - expE[5]) / expE[5];
		
		TGraphErrors *ResiduePlot = new TGraphErrors(6, expE, Residues, 0, 0);
		ResiduePlot->SetTitle("Residues for each peak");
		ResiduePlot->GetXaxis()->SetTitle("ADC Energies");
		ResiduePlot->GetYaxis()->SetTitle("\% Error in calibrated energy ((Fit - Exp)/Exp)");
		FormatGraph(ResiduePlot);
		ResiduePlot->Draw();
		ResCanvas->BuildLegend(0.15,0.15,0.30,0.40); 	// legend in bottom left

	} else if (option.compare("over") == 0) {
		
		TCanvas *OverCanvas = new TCanvas("Overlay Canvas", "Overlay Canvas", 1);
		gPad->SetLogy();
		CalPlots.Calibrated->Draw();
		
	}

}

// This function accepts a position (pos) and returns the location of the local maximum within
// 10% of that stated position.  This is used to find the actual maximum of an approximate
// peak position, enhancing fit stability.
Double_t SnapToMax(TH1D *h, Double_t pos) {
	h->GetXaxis()->SetRangeUser(pos - 0.05 * pos, pos + 0.05 * pos);
	pos = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
	h->GetXaxis()->SetRange(0, NUMBINS);
	return pos;
}

// This function will estimate the background underneath a peak by pulling a set number of points 
// from the edges of the fit region, pushing them into a TGraph, and fitting to a standard root
// function (defined by func) in order to obtain an estimate for the background at pos.
FitPars BackEst(TH1D *h, Double_t pos, FitWindow win, Int_t bins, string func) {
	Float_t NormPos = ((Float_t) h->FindBin(pos))/((Float_t) NUMBINS);

	cout << "Initiating Background Estimation..." << endl;
	//cout << "Peak Bin / NUMBINS = " << NormPos << endl;
	
	Int_t BackBins = bins * NormPos;
	Int_t LoBin = h->FindBin(win.Low);
	Int_t HiBin = h->FindBin(win.High);
	Float_t BackVals[2 * BackBins];
	Float_t BackX[2 * BackBins];
	
	for (Int_t j = 0; j < BackBins; j++) {
		BackVals[j] = h->GetBinContent(j + LoBin);
		BackX[j] = h->GetBinCenter(j + LoBin);
	}
	for (Int_t j = 0; j < BackBins; j++) {
		BackVals[j + BackBins] = h->GetBinContent(HiBin - j);
		BackX[j + BackBins] = h->GetBinCenter(HiBin - j);
	}

	TGraph Back = TGraph(2 * BackBins, BackX, BackVals);
	TF1 BackFit = TF1("BackFit", func.c_str(), win.Low, win.High);
	Back.Fit("BackFit");
	
	FitPars pars;
	pars.Offset = BackFit.GetParameter(0);
	pars.Slope = BackFit.GetParameter(1);
	return pars;
}

// This is a simple function to format a TGraphErrors object in a more readable form.
void FormatGraph(TGraphErrors *Graph) {
	Graph->SetMarkerColor(4);
	Graph->SetMarkerStyle(21);
	Graph->SetLineColor(1);
	Graph->SetLineWidth(2);
	Graph->GetYaxis()->SetTitleOffset(1.4);
	Graph->GetXaxis()->SetTitleOffset(1.2);
}
