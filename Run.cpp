#include <Run.h>

Run::Run() {
    trA = new TreeAccess();
    cfd = new CFDHitFinder<double>();
    tmA = new TimeMatchingAlg<double>();

    cfd->SetParam(kADCNBits, 12);
    cfd->SetParam(kADCDynamicRange, 1);
    cfd->SetParam(kCFDThreshold, 0.5);
    cfd->SetParam(kADCOffset, 0);
    cfd->SetParam(kCoincidenceWindowLowerLim, 0);
    cfd->SetParam(kCoincidenceWindowUpperLim, 20 );
    cfd->SetParam(kRawHitFinderThresholdInNoiseSigma, 3);
    cfd->SetParam(kGSFilterWindow, 17);
    cfd->SetParam(kGSFilterDegree, 3);
    cfd->SetParam(kConsecutiveHitSeperationDurationInTicks, 10);
    cfd->SetParam(kShortRawHitIgnoringDurationInTicks, 5);
    cfd->SetParam(kNSamplingPoints, 1024);
    cfd->SetParam(kTimeSamplingInterval, 0.2);

    tmA->SetCoincidenceTimeWindow(0, 20, 0.2);

    Channel3 = {};
    Channel5 = {};

    Channel3.waveform = new Double_t[1024];
    Channel5.waveform = new Double_t[1024];
    timebase  = new Double_t[1024];

    for (size_t j = 0; j < 1024; j++) {
        *(timebase + j) = ((Double_t) j) * 0.2; // hardcode
    }
}

Run::~Run() {
}

void Run::Go() {
    Channel3.baseline  = 0;
    Channel5.baseline  = 0;
    Channel3.noiseband = 0;
    Channel5.noiseband = 0;
    TH1D* tofHist = new TH1D("", "", 200, 0, 10);
//    for (int idx = 0; idx < 100; idx++) {
    for (int idx = 0; idx < trA->GetEntries(); idx++) {
        trA->GetWaveformChannel(idx, 3, Channel3.waveform);
        trA->GetWaveformChannel(idx, 5, Channel5.waveform);
        std::vector<uint16_t> wf;

        wf.clear();
        for (size_t j = 0; j < 1024; j++) {
            wf.push_back((uint16_t)*(Channel3.waveform + j));
        }
        cfd->SetWaveform(wf);
        cfd->Go();
        Channel3.hitCol = cfd->GetHitCollection();
        Channel3.baseline  = (Double_t) cfd->GetPedestal();
        Channel3.noiseband = (Double_t) cfd->GetNoiseSigma();

        wf.clear();
        for (size_t j = 0; j < 1024; j++) {
            wf.push_back(*(Channel5.waveform + j));
        }
        cfd->SetWaveform(wf);
        cfd->Go();
        Channel5.hitCol = cfd->GetHitCollection();
        Channel5.baseline  = (Double_t) cfd->GetPedestal();
        Channel5.noiseband = (Double_t) cfd->GetNoiseSigma();

//        PlotWaveforms(idx);

        std::vector<std::pair<double, double> > TOF = tmA->GetAllTimeOfFlight(Channel5.hitCol, Channel3.hitCol);
        for (size_t k = 0; k < TOF.size(); k++) {
            tofHist->Fill(TOF.at(k).first);
        }
    }

    TCanvas* c = new TCanvas();
    tofHist->Draw();
    tofHist->Fit("gaus[3]", "RQ", 0, 4);
    c->SaveAs("testDS-US.pdf");
}

void Run::PlotWaveforms(Long64_t evtIndex) {
    TCanvas* canvas = new TCanvas();
    canvas->Divide(1, 2);

    canvas->cd(1);
    TGraph* gr3 = new TGraph(1024, timebase, Channel3.waveform);
    gr3->SetTitle("Channel 3 - TOF US");
    gr3->GetXaxis()->SetRangeUser(0, 1024 * 0.2);
    gr3->GetXaxis()->SetTitle("ns");
    gr3->GetXaxis()->CenterTitle();
    gr3->GetYaxis()->SetRangeUser(1500, 4500);
    gr3->GetYaxis()->SetTitle("ADC");
    gr3->GetYaxis()->CenterTitle();
    TLine* line_baseline3 = new TLine(0, Channel3.baseline, 1024 * 0.2, Channel3.baseline);
    line_baseline3->SetLineColor(kRed);
    TLine* hitStartTime = new TLine[Channel3.hitCol.size()];
    gr3->Draw("APL");
    line_baseline3->Draw();

    canvas->cd(2);
    TGraph* gr5 = new TGraph(1024, timebase, Channel5.waveform);
    gr5->SetTitle("Channel 5 - TOF US");
    gr5->GetXaxis()->SetRangeUser(0, 1024 * 0.2);
    gr5->GetXaxis()->SetTitle("ns");
    gr5->GetXaxis()->CenterTitle();
    gr5->GetYaxis()->SetRangeUser(1500, 4500);
    gr5->GetYaxis()->SetTitle("ADC");
    gr5->GetYaxis()->CenterTitle();
    TLine* line_baseline5 = new TLine(0, Channel5.baseline, 1024 * 0.2, Channel5.baseline);
    line_baseline5->SetLineColor(kRed);
    gr5->Draw("APL");
    line_baseline5->Draw();

    canvas->SaveAs(Form("Test_%i.pdf", evtIndex));

    return;
}

