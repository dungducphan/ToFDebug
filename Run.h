#pragma once

#include <TreeAccess.h>
#include <HitFinderAlg.h>
#include <TimeMatchingAlg.h>

#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TString.h>
#include <TH2.h>

struct WaveformStruct {
    Double_t* waveform;
    Double_t  baseline;
    Double_t  noiseband;
    std::map<size_t, hit_t<double> > hitCol;
};

class Run {
public:
    Run();
    virtual ~Run();

    void Go();

private:
    TreeAccess* trA;
    CFDHitFinder<double>* cfd;
    TimeMatchingAlg<double>* tmA;

    WaveformStruct Channel3;
    WaveformStruct Channel5;

    Double_t* timebase;

private:
    void PlotWaveforms(Long64_t evtIndex);
};
