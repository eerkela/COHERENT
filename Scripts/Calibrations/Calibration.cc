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
	Double_t TlX;
	Double_t TlXErr;
	Double_t TlSigma;
	Double_t TlSigmaErr;
	Double_t TlCount;
	
	Double_t KX;
	Double_t KXErr;
	Double_t KSigma;
	Double_t KSigmaErr;
	Double_t KCount;
	
	Double_t CsX;
	Double_t CsXErr;
	Double_t CsSigma;
	Double_t CsSigmaErr;
	Double_t CsCount;
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

void Calibration(string mode, string option) {

	gStyle->SetOptFit(1111);

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
		
		TH1D *hTest = new TH1D("hTest", "Finding Tl pos", NUMBINS, 0, OverflowBinPos);
		DataChain[i]->Draw("energy >> hTest");
		
		//string hMuName = "hMu" + std::to_string(i+1);
		//hMu[i] = CalPlots[i].Raw->Clone(hMuName.c_str());		
		//c[i]->Draw(("energy >> " + hMuName).c_str());
		
		// Finding Thalium:
		cout << "Fitting Thalium Peak..." << endl;
		Int_t BinNum = NUMBINS;
		Int_t BinCount = 0;
		while (BinCount < (2.5 * time)) { // was 13.5
			BinCount += (Int_t)hTest->GetBinContent(BinNum);
			BinNum--;
			if (BinNum == 0) {break;}
		}
		Double_t TlX = hTest->GetXaxis()->GetBinCenter(BinNum);
		TlX = SnapToMax(hTest, TlX);
		Double_t TlNormPos = (Double_t)hTest->FindBin(TlX)/NUMBINS;
		cout << "TlNormPos: " << TlNormPos << endl;

		// Redraw histogram with constant number of bins below Tl peak to enhance fits.
		delete hTest;
		string hName = "h" + to_string(i+1);
		CalPlots[i].Raw = new TH1D(hName.c_str(), "Uncalibrated Spectrum", (Int_t)(500.0/TlNormPos), 0, OverflowBinPos);
		DataChain[i]->Draw(("energy >> " + hName).c_str());
		
		TlX = SnapToMax(CalPlots[i].Raw, TlX);
		Double_t TlCount = CalPlots[i].Raw->GetBinContent(CalPlots[i].Raw->FindBin(TlX));

		FitWindow TlWindow;		
		TlWindow.Low = TlX - 0.1 * TlX;
		TlWindow.High = TlX + 0.1 * TlX;

		FitPars TlBackPars = BackEst(CalPlots[i].Raw, TlX, TlWindow, 500, "pol1");
		
		string TlFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2)";
		TlFitFunc += " + [3] + [4] * x";

		Double_t TlParameters[5];
		TlParameters[0] = TlCount;
		TlParameters[1] = TlX;
		TlParameters[2] = 0.05 * TlX;
		TlParameters[3] = TlBackPars.Offset;
		TlParameters[4] = TlBackPars.Slope;

		TF1 *TlFit = new TF1("TlFit", TlFitFunc.c_str(), TlWindow.Low, TlWindow.High);
		TlFit->SetParameters(TlParameters);
		CalPlots[i].Raw->Fit("TlFit", "R+ll");
		
		PeakInfo[i].TlX = TlFit->GetParameter(1);
		PeakInfo[i].TlXErr = TlFit->GetParError(1);
		PeakInfo[i].TlSigma = TlFit->GetParameter(2);
		PeakInfo[i].TlSigmaErr = TlFit->GetParError(2);
		PeakInfo[i].TlCount = TlFit->GetParameter(0);	

		// Finding Potassium:
		cout << "Fitting Potassium peak..." << endl;
		Double_t KX = (1460.0/2615.0) * PeakInfo[i].TlX;
		KX = SnapToMax(CalPlots[i].Raw, KX);
		Double_t KCount = CalPlots[i].Raw->GetBinContent(CalPlots[i].Raw->FindBin(KX));
		
		FitWindow KWindow;
		KWindow.Low = KX - 0.1 * KX;
		KWindow.High = KX + 0.1 * KX;
		
		FitPars KBackPars = BackEst(CalPlots[i].Raw, KX, KWindow, 500, "expo");
		
		string KFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2)";
		KFitFunc += " + exp([3] + [4] * x)";

		Double_t KParameters[5];
		KParameters[0] = KCount;
		KParameters[1] = KX;
		KParameters[2] = 0.05 * KX;
		KParameters[3] = KBackPars.Offset;
		KParameters[4] = KBackPars.Slope;

		TF1 *KFit = new TF1("KFit", KFitFunc.c_str(), KWindow.Low, KWindow.High);
		KFit->SetParameters(KParameters);
		CalPlots[i].Raw->Fit("KFit", "R+ll");
		
		PeakInfo[i].KX = KFit->GetParameter(1);
		PeakInfo[i].KXErr = KFit->GetParError(1);
		PeakInfo[i].KSigma = KFit->GetParameter(2);
		PeakInfo[i].KSigmaErr = KFit->GetParError(2);
		PeakInfo[i].KCount = KFit->GetParameter(0);

		// Finding Cs:
		cout << "Fitting Cesium Peak" << endl;
		Double_t CsX = (661.6/2615.0) * PeakInfo[i].TlX;
		CsX = SnapToMax(CalPlots[i].Raw, CsX);
		Double_t CsCount = CalPlots[i].Raw->GetBinContent(CalPlots[i].Raw->FindBin(CsX));
		
		FitWindow CsWindow;
		CsWindow.Low = CsX - 0.15 * CsX;
		CsWindow.High = CsX + 0.15 * CsX;
	
		FitPars CsBackPars = BackEst(CalPlots[i].Raw, CsX, CsWindow, 500, "expo");

		string CsFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2)";
		//CsFitFunc += " + [3]*exp(-0.5*((x-[4])/[5])^2)";
		CsFitFunc += " + exp([3] + [4] * x)";

		Double_t CsParameters[5];
		CsParameters[0] = CsCount;
		CsParameters[1] = CsX;
		CsParameters[2] = 0.05 * CsX;
		CsParameters[3] = CsBackPars.Offset;
		CsParameters[4] = CsBackPars.Slope;

		TF1 *CsFit = new TF1("CsFit", CsFitFunc.c_str(), CsWindow.Low, CsWindow.High);
		CsFit->SetParameters(CsParameters);
		CalPlots[i].Raw->Fit("CsFit", "R+ll");
		
		PeakInfo[i].CsX = CsFit->GetParameter(1);
		PeakInfo[i].CsXErr = CsFit->GetParError(1);
		PeakInfo[i].CsSigma = CsFit->GetParameter(2);
		PeakInfo[i].CsSigmaErr = CsFit->GetParError(2);
		PeakInfo[i].CsCount = CsFit->GetParameter(0);

		// Collecting data:
		Double_t expE[3] = {2615.0, 1460.0, 661.6};
		Double_t fitE[3]; 
			fitE[0] = PeakInfo[i].TlX;
			fitE[1] = PeakInfo[i].KX;
			fitE[2] = PeakInfo[i].CsX;
		Double_t fitEErr[3];
			fitEErr[0] = PeakInfo[i].TlXErr;
			fitEErr[1] = PeakInfo[i].KXErr;
			fitEErr[2] = PeakInfo[i].CsXErr;
		Double_t fitSigma[3];
			fitSigma[0] = PeakInfo[i].TlSigma / PeakInfo[i].TlX;
			fitSigma[1] = PeakInfo[i].KSigma / PeakInfo[i].KX;
			fitSigma[2] = PeakInfo[i].CsSigma / PeakInfo[i].CsX;
		Double_t fitSigmaErr[3];
			fitSigmaErr[0] = PeakInfo[i].TlSigmaErr / PeakInfo[i].TlX;
			fitSigmaErr[1] = PeakInfo[i].KSigmaErr / PeakInfo[i].KX;
			fitSigmaErr[2] = PeakInfo[i].CsSigmaErr / PeakInfo[i].CsX;
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
		Int_t MaxCalBin = (OverflowBinPos - CalInfo[i].Offset) / CalInfo[i].Slope;
		CalPlots[i].Calibrated = new TH1D(Name.c_str(), Label.c_str(), 20e3, 0, MaxCalBin);
		CalPlots[i].Calibrated->SetLineColor(i+1);
		if (i == 4) {CalPlots[i].Calibrated->SetLineColor(i+2);}
		CalPlots[i].Calibrated->GetXaxis()->SetTitle("Calibrated Energy (keV)");
		CalPlots[i].Calibrated->GetYaxis()->SetTitle("Count");
		
		string calibration = "(energy - " + to_string(CalInfo[i].Offset);
		calibration += ") / " + to_string(CalInfo[i].Slope) ;
		calibration += " >> " + Name;
		DataChain[i]->Draw(calibration.c_str(), "", "goff");

		Int_t NoiseCount;
		Int_t NoiseBin = CalPlots[i].Calibrated->FindBin(50);
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
			FIRResolution[k] = pow(PeakInfo[k].CsSigma / PeakInfo[k].CsX, 2.0);
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
		
		Double_t Resolution = PeakInfo[3].CsSigma / PeakInfo[3].CsX;
		cout << "Normalized detector resolution (width of Cs peak" << endl;
		cout << " / mean at 3rd position):" << Resolution << endl;
		
		Double_t CsPos[NUMFILES];
		Double_t CsPosErrs[NUMFILES];
		Double_t CsSig[NUMFILES];
		Double_t CsSigErrs[NUMFILES];
		for (Int_t k = 0; k < NUMFILES; k++) {
			CsPos[k] = PeakInfo[k].CsX;
			CsPosErrs[k] = PeakInfo[k].CsXErr;
			CsSig[k] = PeakInfo[k].CsSigma / PeakInfo[k].CsX;
			CsSigErrs[k] = PeakInfo[k].CsSigmaErr / PeakInfo[k].CsX;
		}

		Double_t XAxis[5];
		if (mode.compare("pos") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j + 1;}
		} else if (mode.compare("volt") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j * 100 + 600;}
		} else if (mode.compare("fir") == 0) {
			for (Int_t j = 0; j < 5; j++) {XAxis[j] = j * 50 + 100;}
		}
		
		TCanvas *CsPosCanvas = new TCanvas("CsPosCanvas", "CsPosCanvas");
		TGraphErrors *CsPosGraph = new TGraphErrors(5, XAxis, CsPos, 0, CsPosErrs);
		CsPosGraph->SetTitle(("Cs Peak Energy vs " + GetDepVar(mode)).c_str());
		CsPosGraph->GetXaxis()->SetTitle(GetDepVar(mode).c_str());
		CsPosGraph->GetYaxis()->SetTitle("Uncalibrated Cs Peak Energy");
		
		FormatGraph(CsPosGraph);
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
		TGraphErrors *CsResGraph = new TGraphErrors(5, XAxis, CsSig, 0, CsSigErrs);
		CsResGraph->SetTitle(("Cs Peak Width vs " + GetDepVar(mode)).c_str());
		CsResGraph->GetXaxis()->SetTitle(GetDepVar(mode).c_str());
		CsResGraph->GetYaxis()->SetTitle("Normalized Cs Peak Width (Cs sigma / Cs mean)");
		
		FormatGraph(CsResGraph);
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
			Double_t Range = PeakInfo[k].TlX + 0.1 * PeakInfo[k].TlX;
			CalPlots[k].Raw->GetXaxis()->SetRangeUser(0, Range);
		}
		
	} else if (option.compare("res") == 0) {

		TCanvas *ResCanvas = new TCanvas("Residue Canvas", "Residue Canvas", 1);
		TMultiGraph *ResComp = new TMultiGraph("ResComp", "ResComp");
		Double_t expE[3] = {2615.0, 1460.0, 661.6};
		for (Int_t k = 0; k < NUMFILES; k++) {
			Double_t Residues[3];
			
			Double_t E = (PeakInfo[k].TlX - CalInfo[k].Offset) / CalInfo[k].Slope;
			Residues[0] = 100 * (E - expE[0]) / expE[0];
				
			E = (PeakInfo[k].KX - CalInfo[k].Offset) / CalInfo[k].Slope;
			Residues[1] = 100 * (E - expE[1]) / expE[1];			
			
			E = (PeakInfo[k].CsX - CalInfo[k].Offset) / CalInfo[k].Slope;
			Residues[2] = 100 * (E - expE[2]) / expE[2];
			
			TGraphErrors *ResiduePlot = new TGraphErrors(3, expE, Residues, 0, 0);
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
