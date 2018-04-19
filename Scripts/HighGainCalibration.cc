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
void GetData(TChain *c[NUMFILES], string mode);
string GetDepVar(string mode);
string GetLabels(string mode, Int_t k);
Double_t GetRuntime(TChain *c, string mode);
Double_t SnapToMax(TH1D *h, Double_t pos);

void HighGainCalibration_updated(string mode, string option) {

	//gStyle->SetOptFit(1111);

	TCanvas *RawCanvas = new TCanvas("Raw Spectra", "Raw Spectra", 1200, 800);
	RawCanvas->Divide(2,3);
	TChain *DataChain[NUMFILES];
	GetData(DataChain, mode);

	//TF1 *muonFit = new TF1("muonFit", "[0]*exp(-0.5*((x-[1])/[2]+exp(-(x-[1])/[2])))+[3]*x+[4]",0,0);
	//Double_t MuX[NUMFILES];
	
	PlotStruct CalPlots[NUMFILES];
	PeakStruct PeakInfo[NUMFILES];
	CalStruct CalInfo[NUMFILES];
	
	for (Int_t i = 0; i < NUMFILES; i++) {
		cout << endl;
		cout << "---------------------------------------------------------------" << endl;
		cout << endl;
		
		cout << "Beginning Calibration " << i + 1 << "..." << endl;
		RawCanvas->cd(i+1);
		gPad->SetLogy();
		
		Double_t time = GetRuntime(DataChain[i], mode);
		Double_t OverflowBinPos = DataChain[i]->GetMaximum("energy");
		OverflowBinPos = OverflowBinPos + 0.01 * OverflowBinPos;
		
		TH1D *hTest = new TH1D("hTest", "Finding 356 pos", NUMBINS, 0, OverflowBinPos);
		DataChain[i]->Draw("energy >> hTest", "channel==0", "");
		
		//string hMuName = "hMu" + std::to_string(i+1);
		//hMu[i] = CalPlots[i].Raw->Clone(hMuName.c_str());		
		//c[i]->Draw(("energy >> " + hMuName).c_str());
		
		// Fitting to complicated 356 keV region:
		cout << "Fitting 356 keV Region..." << endl;
		Int_t BinNum = NUMBINS;
		Int_t BinCount = 0;
		while (BinCount < (800.0 * time)) {
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
		string hName = "h" + to_string(i+1);
		CalPlots[i].Raw = new TH1D(hName.c_str(), "Uncalibrated Spectrum", (Int_t)(500.0/NormPos356), 0, OverflowBinPos);
		DataChain[i]->Draw(("energy >> " + hName).c_str(), "channel==0", "");
		
		X356 = SnapToMax(CalPlots[i].Raw, X356);
		Double_t Count356 = CalPlots[i].Raw->GetBinContent(CalPlots[i].Raw->FindBin(X356));

		FitWindow Window356;		
		Window356.Low = X356 - 0.4 * X356;
		Window356.High = X356 + 0.22 * X356;

		FitPars BackPars356 = BackEst(CalPlots[i].Raw, X356, Window356, 500, "pol1");
		
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
		CalPlots[i].Raw->Fit("Fit356", "R+ll");
		
		PeakInfo[i].X383 = Fit356->GetParameter(1);
		PeakInfo[i].X383Err = Fit356->GetParError(1);
		PeakInfo[i].Sigma383 = Fit356->GetParameter(2);
		PeakInfo[i].Sigma383Err = Fit356->GetParError(2);
		PeakInfo[i].Count383 = Fit356->GetParameter(0);
		PeakInfo[i].X356 = Fit356->GetParameter(4);
		PeakInfo[i].X356Err = Fit356->GetParError(4);
		PeakInfo[i].Sigma356 = Fit356->GetParameter(2);
		PeakInfo[i].Sigma356Err = Fit356->GetParError(2);
		PeakInfo[i].Count356 = Fit356->GetParameter(3);	
		PeakInfo[i].X302 = Fit356->GetParameter(6);
		PeakInfo[i].X302Err = Fit356->GetParError(6);
		PeakInfo[i].Sigma302 = Fit356->GetParameter(2);
		PeakInfo[i].Sigma302Err = Fit356->GetParError(2);
		PeakInfo[i].Count302 = Fit356->GetParameter(5);
		PeakInfo[i].X276 = Fit356->GetParameter(8);
		PeakInfo[i].X276Err = Fit356->GetParError(8);
		PeakInfo[i].Sigma276 = Fit356->GetParameter(2);
		PeakInfo[i].Sigma276Err = Fit356->GetParError(2);
		PeakInfo[i].Count276 = Fit356->GetParameter(7);	

		cout << endl;
		cout<<"383 keV peak X: "<<PeakInfo[i].X383<<", Count: "<<PeakInfo[i].Count383<<endl;
		cout<<"356 keV peak X: "<<PeakInfo[i].X356<<", Count: "<<PeakInfo[i].Count356<<endl;
		cout<<"302 keV peak X: "<<PeakInfo[i].X302<<", Count: "<<PeakInfo[i].Count302<<endl;
		cout<<"276 keV peak X: "<<PeakInfo[i].X276<<", Count: "<<PeakInfo[i].Count276<<endl;
		cout << endl;

		// Fitting 80 keV peaks:
		cout << "Fitting 80 keV Region..." << endl;
		Double_t X80 = (80.0/356.0) * PeakInfo[i].X356;
		X80 = SnapToMax(CalPlots[i].Raw, X80);
		Double_t Count80 = CalPlots[i].Raw->GetBinContent(CalPlots[i].Raw->FindBin(X80));
		
		FitWindow Window80;
		Window80.Low = X80 - 0.275 * X80;
		Window80.High = X80 + 0.25 * X80;
		
		FitPars BackPars80 = BackEst(CalPlots[i].Raw, X80, Window80, 500, "pol1");
		
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
		CalPlots[i].Raw->Fit("Fit80", "R+ll");
		
		PeakInfo[i].X81 = Fit80->GetParameter(1);
		PeakInfo[i].X81Err = Fit80->GetParError(1);
		PeakInfo[i].Sigma81 = Fit80->GetParameter(2);
		PeakInfo[i].Sigma81Err = Fit80->GetParError(2);
		PeakInfo[i].Count81 = Fit80->GetParameter(0);

		cout << endl;
		cout<<"81.0 keV peak X: "<<PeakInfo[i].X81<<", Count: "<<PeakInfo[i].Count81<<endl;
		//cout<<"79.6 keV peak X: "<<PeakInfo[i].X79<<", Count: "<<PeakInfo[i].Count79<<endl;
		cout << endl;

		// Fitting 30 keV region:
		cout << "Fitting 30 keV Region" << endl;
		Double_t X30 = (32.0/356.0) * PeakInfo[i].X356;
		X30 = SnapToMax(CalPlots[i].Raw, X30);
		Double_t Count30 = CalPlots[i].Raw->GetBinContent(CalPlots[i].Raw->FindBin(X30));
		
		FitWindow Window30;
		Window30.Low = X30 - 0.5 * X30;
		Window30.High = X30 + 0.4 * X30;
	
		FitPars BackPars30 = BackEst(CalPlots[i].Raw, X30, Window30, 500, "pol1");

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
		CalPlots[i].Raw->Fit("Fit30", "R+ll");
		
		PeakInfo[i].X30 = Fit30->GetParameter(1);
		PeakInfo[i].X30Err = Fit30->GetParError(1);
		PeakInfo[i].Sigma30 = Fit30->GetParameter(2);
		PeakInfo[i].Sigma30Err = Fit30->GetParError(2);
		PeakInfo[i].Count30 = Fit30->GetParameter(0);

		cout << endl;
		cout<<"30.97 keV peak X: "<<PeakInfo[i].X30<<", Count: "<<PeakInfo[i].Count30<<endl;
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
			fitE[0] = PeakInfo[i].X383;
			fitE[1] = PeakInfo[i].X356;
			fitE[2] = PeakInfo[i].X302;
			fitE[3] = PeakInfo[i].X276;
			fitE[4] = PeakInfo[i].X81;
			fitE[5] = PeakInfo[i].X30;
		Double_t fitEErr[6];
			fitEErr[0] = PeakInfo[i].X383Err;
			fitEErr[1] = PeakInfo[i].X356Err;
			fitEErr[2] = PeakInfo[i].X302Err;
			fitEErr[3] = PeakInfo[i].X276Err;
			fitEErr[4] = PeakInfo[i].X81Err;
			fitEErr[5] = PeakInfo[i].X30Err;
		Double_t fitSigma[6];
			fitSigma[0] = PeakInfo[i].Sigma383 / PeakInfo[i].X383;
			fitSigma[1] = PeakInfo[i].Sigma356 / PeakInfo[i].X356;
			fitSigma[2] = PeakInfo[i].Sigma302 / PeakInfo[i].X302;
			fitSigma[3] = PeakInfo[i].Sigma276 / PeakInfo[i].X276;
			fitSigma[4] = PeakInfo[i].Sigma81 / PeakInfo[i].X81;
			fitSigma[5] = PeakInfo[i].Sigma30 / PeakInfo[i].X30;
		Double_t fitSigmaErr[6];
			fitSigmaErr[0] = PeakInfo[i].Sigma383Err / PeakInfo[i].X383;
			fitSigmaErr[1] = PeakInfo[i].Sigma356Err / PeakInfo[i].X356;
			fitSigmaErr[2] = PeakInfo[i].Sigma302Err / PeakInfo[i].X302;
			fitSigmaErr[3] = PeakInfo[i].Sigma276Err / PeakInfo[i].X276;
			fitSigmaErr[4] = PeakInfo[i].Sigma81Err / PeakInfo[i].X81;
			fitSigmaErr[5] = PeakInfo[i].Sigma30Err / PeakInfo[i].X30;
		
		CalPlots[i].FitPlot = new TGraphErrors(3, expE, fitE, 0, fitEErr);
		CalPlots[i].FitPlot->SetTitle(GetLabels(mode, i).c_str());
		CalPlots[i].FitPlot->SetMarkerColor(i+1);
		if (i == 4) {CalPlots[i].FitPlot->SetMarkerColor(i+2);}
		CalPlots[i].FitPlot->SetMarkerStyle(21);
		CalPlots[i].FitPlot->SetLineColor(1);
		CalPlots[i].FitPlot->SetLineWidth(2);		
		
		CalPlots[i].SigmaPlot = new TGraphErrors(3, expE, fitSigma, 0, fitSigmaErr); 
		CalPlots[i].SigmaPlot->SetTitle(GetLabels(mode, i).c_str());
		CalPlots[i].SigmaPlot->SetMarkerColor(i+1);
		CalPlots[i].SigmaPlot->SetLineColor(1);
		if (i == 4) {CalPlots[i].SigmaPlot->SetMarkerColor(i+2);}
		CalPlots[i].SigmaPlot->SetMarkerStyle(21);
		CalPlots[i].SigmaPlot->SetLineWidth(1);
		
		TF1 *CalibrationFit = new TF1("CalibrationFit", "pol1", 0, 4500e3);
		CalPlots[i].FitPlot->Fit("CalibrationFit", "", "", 0, 0);
		
		CalInfo[i].Slope = CalibrationFit->GetParameter(1);
		CalInfo[i].SlopeError = CalibrationFit->GetParError(1);
		CalInfo[i].Offset = CalibrationFit->GetParameter(0);
		CalInfo[i].OffsetError = CalibrationFit->GetParError(0);

		string Name = "Calibrated" + to_string(i);
		string Label = GetLabels(mode, i);
		CalPlots[i].Calibrated = new TH1D(Name.c_str(), Label.c_str(), 20e3, 0, 50e3);
		CalPlots[i].Calibrated->SetLineColor(i+1);
		if (i == 4) {CalPlots[i].Calibrated->SetLineColor(i+2);}
		CalPlots[i].Calibrated->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		CalPlots[i].Calibrated->GetYaxis()->SetTitle("Count");
		
		string calibration = "(energy - " + to_string(CalInfo[i].Offset);
		calibration += ") / " + to_string(CalInfo[i].Slope) ;
		calibration += " >> " + Name;
		DataChain[i]->Draw(calibration.c_str(), "", "goff");

		Int_t NoiseCount;
		Int_t NoiseBin = CalPlots[i].Calibrated->FindBin(10);
		Int_t NoiseMin = 2 * CalPlots[i].Calibrated->GetBinContent(NoiseBin);
		while(NoiseCount < 2 * NoiseMin) {
			NoiseCount = CalPlots[i].Calibrated->GetBinContent(NoiseBin);
			if (NoiseCount < NoiseMin) {NoiseMin = NoiseCount;}
			NoiseBin--;
			if (NoiseBin == 0) {break;}
		}
		CalInfo[i].NoiseWall = CalPlots[i].Calibrated->GetXaxis()->GetBinCenter(NoiseBin);
		cout << "Noise Wall at: " << CalInfo[i].NoiseWall << " keV" << endl;
		
		cout << endl;
		cout << "---------------------------------------------------------------" << endl;
		cout << endl;
	}

	if (mode.compare("fir") == 0) {
		
		TCanvas *FIRCanvas = new TCanvas("FIR", "FIR", 1);
		Double_t PeakingTimes[5] = {100, 150, 200, 250, 300};
		Double_t FIRResolution[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {
			FIRResolution[k] = pow(PeakInfo[k].Sigma356 / PeakInfo[k].X356, 2.0);
		}
		TGraph *FIROptPlot = new TGraph(5, PeakingTimes, FIRResolution);
		string FIRFitFunc = "([0]/x)^2 + [1]^2 + ([3]*x)^2";
		TF1 *FIRFit = new TF1("FIRFit", FIRFitFunc.c_str(), 100.0, 300.0);
		//FIRFit->SetParameters(
		//FIROptPlot->Fit("FIRFit", "", "", 0, 0);
		//Double_t OptPeakingTime = FIRFit->GetMinimum();
		
		FIROptPlot->Draw();
		//cout << "OPTIMAL PEAKING TIME: " << OptPeakingTime << endl;
		
	} else if (mode.compare("pos") == 0) {
		
		Double_t Resolution = PeakInfo[3].Sigma356 / PeakInfo[3].X356;
		cout << "Normalized detector resolution (width of Cs peak" << endl;
		cout << " / mean at 3rd position):" << Resolution << endl;
		
		Double_t Pos356[NUMFILES];
		Double_t Pos356Errs[NUMFILES];
		Double_t Sig356[NUMFILES];
		Double_t Sig356Errs[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {
			Pos356[k] = PeakInfo[k].X356;
			Pos356Errs[k] = PeakInfo[k].X356Err;
			Sig356[k] = PeakInfo[k].Sigma356 / PeakInfo[k].X356;
			Sig356Errs[k] = PeakInfo[k].Sigma356Err / PeakInfo[k].X356;
		}

		Double_t XAxis[5];
		if (mode.compare("pos") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j + 1;}
		} else if (mode.compare("volt") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j * 100 + 600;}
		} else if (mode.compare("fir") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j * 50 + 100;}
		}
		
		TCanvas *Pos356Canvas = new TCanvas("Pos356Canvas", "Pos356Canvas");
		TGraphErrors *Pos356Graph = new TGraphErrors(5, XAxis, Pos356, 0, Pos356Errs);
		Pos356Graph->SetTitle(("356 Peak Energy vs " + GetDepVar(mode)).c_str());
		Pos356Graph->GetXaxis()->SetTitle(GetDepVar(mode).c_str());
		Pos356Graph->GetYaxis()->SetTitle("Uncalibrated 356 Peak Energy");
		
		FormatGraph(Pos356Graph);
		/*Pos356Graph->SetMarkerColor(4);
		Pos356Graph->SetMarkerStyle(21);
		Pos356Graph->SetLineColor(1);
		Pos356Graph->SetLineWidth(2);
		Pos356Graph->GetYaxis()->SetTitleOffset(1.5);
		Pos356Graph->GetXaxis()->SetTitleOffset(1.2);*/
		Pos356Graph->Draw();
		Pos356Canvas->Print("356EvsPos.pdf[");
		Pos356Canvas->Print("356EvsPos.pdf");
		Pos356Canvas->Print("356EvsPos.pdf]");
						
		TCanvas *Res356Canvas = new TCanvas("Res356Canvas", "Res356Canvas");
		TGraphErrors *Res356Graph = new TGraphErrors(5, XAxis, Sig356, 0, Sig356Errs);
		Res356Graph->SetTitle(("356 Peak Width vs " + GetDepVar(mode)).c_str());
		Res356Graph->GetXaxis()->SetTitle(GetDepVar(mode).c_str());
		Res356Graph->GetYaxis()->SetTitle("Normalized 356 Peak Width (356 sigma / 356 mean)");
		
		FormatGraph(Res356Graph);
		/*Res356Graph->SetMarkerColor(4);
		Res356Graph->SetMarkerStyle(21);
		Res356Graph->SetLineColor(1);
		Res356Graph->SetLineWidth(2);
		Res356Graph->GetYaxis()->SetTitleOffset(1.5);
		Res356Graph->GetXaxis()->SetTitleOffset(1.2);*/
		Res356Graph->Draw();
		Res356Canvas->Print("356ResvsPos.pdf[");
		Res356Canvas->Print("356ResvsPos.pdf");
		Res356Canvas->Print("356ResvsPos.pdf]");
		
	}


	if (option.compare("cal") == 0) {	
		
		TCanvas *CalCanvas = new TCanvas("Calibration Canvas", "Calibration Canvas", 1);
		TMultiGraph *CalComp = new TMultiGraph("CalComp", "CalComp");
		for (Int_t k = 0; k < NUMFILES; k++) {
			CalComp->Add((TGraph*)CalPlots[k].FitPlot->Clone());
			// Y-axis is ADC energy, X is real energy
		}
		CalComp->Draw("ALP");
		CalComp->SetTitle(("Calibration Curves for Each " + GetDepVar(mode)).c_str());
		CalComp->GetYaxis()->SetTitle("ADC Energy");
		CalComp->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		CalCanvas->BuildLegend(0.15,0.6,0.30,0.85); 	// legend in top left
		
	} else if (option.compare("sig") == 0) {
		
		TCanvas *SigmaCanvas = new TCanvas("Sigma Canvas", "Sigma Canvas", 1);
		TMultiGraph *SigmaComp = new TMultiGraph("SigmaComp", "SigmaComp");
		for (Int_t k = 0; k < NUMFILES; k++) {
			SigmaComp->Add((TGraphErrors*)CalPlots[k].SigmaPlot->Clone());	
		}
		SigmaComp->Draw("ALP");
		SigmaComp->SetTitle(("Resolution vs " + GetDepVar(mode)).c_str());
		SigmaComp->GetYaxis()->SetTitle("Fractional Sigma (sigma/mean)");
		SigmaComp->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		SigmaCanvas->BuildLegend(0.7,0.6,0.85,0.85); 	// legend in top right

	} else if (option.compare("fit") == 0) {

		for (Int_t k = 0; k < NUMFILES; k++) {
			RawCanvas->cd(k+1);
			Double_t Range = PeakInfo[k].X383 + 0.3 * PeakInfo[k].X383;
			CalPlots[k].Raw->GetXaxis()->SetRangeUser(0, Range);
		}
		
	} else if (option.compare("res") == 0) {

		TCanvas *ResCanvas = new TCanvas("Residue Canvas", "Residue Canvas", 1);
		TMultiGraph *ResComp = new TMultiGraph("ResComp", "ResComp");
		Double_t expE[6];
			expE[0] = 383.8485;
			expE[1] = 356.0129;
			expE[2] = 302.8508;
			expE[3] = 276.3989;
			expE[4] = 80.9979;
			expE[5] = 30.973;
		for (Int_t k = 0; k < NUMFILES; k++) {
			Double_t Residues[6];
			
			Double_t E = (PeakInfo[k].X383 - CalInfo[k].Offset) / CalInfo[k].Slope;
			Residues[0] = 100 * (E - expE[0]) / expE[0];
				
			E = (PeakInfo[k].X356 - CalInfo[k].Offset) / CalInfo[k].Slope;
			Residues[1] = 100 * (E - expE[1]) / expE[1];			
			
			E = (PeakInfo[k].X302 - CalInfo[k].Offset) / CalInfo[k].Slope;
			Residues[2] = 100 * (E - expE[2]) / expE[2];

			E = (PeakInfo[k].X276 - CalInfo[k].Offset) / CalInfo[k].Slope;
			Residues[3] = 100 * (E - expE[3]) / expE[3];

			E = (PeakInfo[k].X81 - CalInfo[k].Offset) / CalInfo[k].Slope;
			Residues[4] = 100 * (E - expE[4]) / expE[4];
			
			E = (PeakInfo[k].X30 - CalInfo[k].Offset) / CalInfo[k].Slope;
			Residues[5] = 100 * (E - expE[5]) / expE[5];
			
			TGraphErrors *ResiduePlot = new TGraphErrors(6, expE, Residues, 0, 0);
			ResiduePlot->SetTitle(GetLabels(mode, k).c_str());
			ResiduePlot->SetLineColor(1);
			ResiduePlot->SetMarkerColor(k+1);
			if (k == 4) {ResiduePlot->SetMarkerColor(k+2);}			
			ResiduePlot->SetMarkerStyle(21);
			ResiduePlot->SetLineWidth(1);
			ResComp->Add((TGraphErrors*) ResiduePlot->Clone());
		}
		ResComp->Draw("ALP");
		ResComp->SetTitle(("Residues for " + GetDepVar(mode) + " Variation").c_str());
		ResComp->GetXaxis()->SetTitle("ADC Energies");
		ResComp->GetXaxis()->SetTitleOffset(1.3);
		ResComp->GetYaxis()->SetTitle("\% Error in calibrated energy ((Fit - Exp)/Exp)");
		ResComp->GetYaxis()->SetTitleOffset(1.2);
		ResCanvas->BuildLegend(0.15,0.15,0.30,0.40); 	// legend in bottom left

	} else if (option.compare("gain") == 0) {

		Double_t Gain[NUMFILES];
		Double_t GainErrors[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {Gain[k] = TMath::Log(CalInfo[k].Slope);}
		for (Int_t k = 0; k < NUMFILES; k++) {GainErrors[k] = CalInfo[k].SlopeError/CalInfo[k].Slope;}	
		
		Double_t XAxis[5];
		if (mode.compare("pos") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j + 1;}
		} else if (mode.compare("volt") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j * 100 + 600;}
		} else if (mode.compare("fir") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j * 50 + 100;}
		}
		
		TCanvas *GainCanvas = new TCanvas("Gain Canvas", "Gain Canvas", 1);
		
		TGraphErrors *GainGraph = new TGraphErrors(5, XAxis, Gain, 0, GainErrors);
		GainGraph->SetTitle(("Detector Gain vs " + GetDepVar(mode)).c_str());
		GainGraph->GetXaxis()->SetTitle(GetDepVar(mode).c_str());
		GainGraph->GetYaxis()->SetTitle("Log(Calibration Slope)");
		FormatGraph(GainGraph);

		/*GainGraph->SetMarkerColor(4);
		GainGraph->SetMarkerStyle(21);
		GainGraph->SetLineColor(1);
		GainGraph->SetLineWidth(2);*/
		
		TF1 *GainFit = new TF1("GainFit", "pol2");
		GainFit->SetParNames("Log(G0)", "Slope", "Curvature");
		GainGraph->Fit("GainFit");
		GainGraph->Draw();
		
		if (mode.compare("volt") == 0) {
			GainCanvas->Print("Gain.pdf[");
			GainCanvas->Print("Gain.pdf");
			GainCanvas->Print("Gain.pdf]");
		}

	} else if (option.compare("over") == 0) {
				
		TCanvas *OverlayCanvas = new TCanvas("OverlayCanvas", " Canvas", 1);
		gPad->SetLogy();
		for (Int_t k = 0; k < NUMFILES; k++) {
			CalPlots[k].Calibrated->Draw("SAME");
		}
		OverlayCanvas->BuildLegend(0.7,0.6,0.85,0.85); 	// legend in top right

	} else if (option.compare("noise") == 0) {
		
		Double_t NoiseWalls[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {
			NoiseWalls[k] = CalInfo[k].NoiseWall;
		}
		
		Double_t XAxis[5];
		if (mode.compare("pos") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j + 1;}
		} else if (mode.compare("volt") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j * 100 + 600;}
		} else if (mode.compare("fir") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j * 50 + 100;}
		}
		
		TCanvas *NoiseCanvas = new TCanvas("NoiseCanvas", "NoiseCanvas", 1);
		TGraphErrors *NoiseGraph = new TGraphErrors(5, XAxis, NoiseWalls, 0, 0);
		NoiseGraph->SetTitle(("Noise Wall Energy vs " + GetDepVar(mode)).c_str());
		NoiseGraph->GetXaxis()->SetTitle(GetDepVar(mode).c_str());
		NoiseGraph->GetYaxis()->SetTitle("Noise Wall Energy (keV)");
		/*NoiseGraph->SetMarkerColor(4);
		NoiseGraph->SetMarkerStyle(21);
		NoiseGraph->SetLineColor(1);
		NoiseGraph->SetLineWidth(2);
		NoiseGraph->GetYaxis()->SetTitleOffset(1.5);
		NoiseGraph->GetXaxis()->SetTitleOffset(1.2);	*/	

		FormatGraph(NoiseGraph);
		NoiseGraph->Draw();
		if (mode.compare("volt") == 0) {
			NoiseCanvas->Print("Noise.pdf[");
			NoiseCanvas->Print("Noise.pdf");
			NoiseCanvas->Print("Noise.pdf]");
		}
		
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

// This function will return the time in seconds for a chain of data (c), assuming that the run
// time of each FIR run is 300 seconds and every other run is 600 seconds.
Double_t GetRuntime(TChain *c, string mode) {
	Double_t time;
	if (mode.compare("fir") == 0) {
		time = (Double_t) c->GetNtrees() * 300;
	} else {
		time = (Double_t) c->GetNtrees() * 600;
	}
	cout << "Runtime in chain: " << time << " seconds" << endl;
	return time;
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

// This function will return the dependent variable for a specified run mode.
string GetDepVar(string mode) {
	if (mode.compare("pos") == 0) {
		return "Position";	
	} else if (mode.compare("volt") == 0) {
		return "Voltage";
	} else if (mode.compare("fir") == 0) {
		return "Peaking Time";
	}
	return "Error: Bad mode";
}

// This function will return the labels used to build legends with the correct parameters, 
// depending on the mode.
string GetLabels(string mode, Int_t k) {
	if (mode.compare("pos") == 0) {
		return "Position " + to_string(k+1);
	} else if (mode.compare("volt") == 0) {
		return to_string((100 * k) + 600) + " V";
	} else if (mode.compare("fir") == 0) {
		return "Peaking Time: " + to_string(k * 50 + 100);
	}
	return "Error: Bad mode";
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
