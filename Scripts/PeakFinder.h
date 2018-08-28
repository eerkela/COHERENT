#ifndef PEAKFINDER_H
#define PEAKFINDER_H

#include <TCanvas.h>
#include <TChain.h>
#include <TH1.h>
#include <TGraphErrors.h>

#include "CalStructs.h"
#include "PeakSet.h"

class PeakFinder { 
private:
	TCanvas *canvas;
	TChain *data;
	TH1D *rawPlot;
	std::vector<TGraphErrors*> backPlots;
	TGraphErrors *calPlot;
	Double_t time;
	Int_t numBins;
	std::string channel;
	PeakSet peaks;	
	PeakInfo pinnedPeak;
	FitResults calibration;
	bool isNumber(std::string input);
	
public:
	PeakFinder(Double_t pinnedEnergy, TChain *c, std::string channel, TApplication *app);
	void addPeakToSet(PeakInfo info);
	PeakInfo findPeak(Double_t energy);
	FitResults backEst(ParWindow win, Double_t range, std::string fitFunc);
	void fit(FitInfo info);
	FitResults findCalibration();
	Measurement calibrate(Measurement uncalibrated);
	std::vector<TGraphErrors*> getBackgroundPlots();
	FitResults getCalibration();
	TGraphErrors *getCalPlot();
	Double_t getOverflowPos();
	PeakSet getPeakSet();
	PeakInfo getPeakInfo(Double_t energy);
	PeakInfo getPinnedPeak();
	TH1D *getRawPlot();
	Double_t snapToMax(TH1D *h, Double_t low, Double_t high);
};

#endif
