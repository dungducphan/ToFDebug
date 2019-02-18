#pragma once

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>

class TreeAccess {
public :
    TreeAccess(TTree *tree = 0);

    virtual ~TreeAccess();
    virtual Int_t Cut(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init(TTree *tree);
    virtual void Loop();
    virtual Bool_t Notify();
    virtual void Show(Long64_t entry = -1);
    virtual Long64_t GetEntries();
    virtual void GetWaveformChannel(Long64_t entry, size_t chID, Double_t* wfch);

private:
    TTree *fChain;
    Int_t fCurrent;

    Long64_t channel3[1024];
    Long64_t channel5[1024];

    TBranch *b_channel3;
    TBranch *b_channel5;
};

