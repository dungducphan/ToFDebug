#include <TreeAccess.h>
#include "TreeAccess.h"


TreeAccess::TreeAccess(TTree *tree) : fChain(0) {
    if (tree == 0) {
        TFile *f = (TFile *) gROOT->GetListOfFiles()->FindObject("../ToFRecoCFD.root.PaddleRun1");
        if (!f || !f->IsOpen()) {
            f = new TFile("../ToFRecoCFD.root.PaddleRun1");
        }
        TDirectory *dir = (TDirectory *) f->Get("../ToFRecoCFD.root.PaddleRun1:/cfd");
        dir->GetObject("WaveformTree", tree);

    }
    Init(tree);
}

TreeAccess::~TreeAccess() {
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t TreeAccess::GetEntry(Long64_t entry) {
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t TreeAccess::LoadTree(Long64_t entry) {
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void TreeAccess::Init(TTree *tree) {
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("channel3", channel3, &b_channel3);
    fChain->SetBranchAddress("channel5", channel5, &b_channel5);
    Notify();
}

Bool_t TreeAccess::Notify() {
    return kTRUE;
}

void TreeAccess::Show(Long64_t entry) {
    if (!fChain) return;
    fChain->Show(entry);
}

Int_t TreeAccess::Cut(Long64_t entry) {
    return 1;
}

void TreeAccess::GetWaveformChannel(Long64_t entry, size_t chID, Double_t* wfch) {
    if (fChain == 0) return;

    Long64_t nbytes = 0, nb = 0;
    Long64_t ientry = LoadTree(entry);
    if (ientry < 0) return;
    nb = fChain->GetEntry(entry);
    nbytes += nb;

    if (chID == 3) {
        for (size_t idx = 0; idx < 1024; idx++) {
            *(wfch + idx) = (Double_t)channel3[idx];
        }
    } else if (chID == 5) {
        for (size_t idx = 0; idx < 1024; idx++) {
            *(wfch + idx) = (Double_t)channel5[idx];
        }
    } else return;
}

Long64_t TreeAccess::GetEntries() {
    return fChain->GetEntries();
}

void TreeAccess::Loop() {
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
    }
}
