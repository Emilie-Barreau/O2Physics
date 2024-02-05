#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TObject.h"
#include "TString.h"
#include "TH1D.h"

#endif
void createFitResFile(TString inPath, TString inFileName, TString outPath, TString outFileName);

void prepareToFit()
{
    TString inPath = "/Users/emiliebarreau/alice/ServiceTask_Analysis/Analysis_LHC22o_pass5_skimQC1/TM_Reassociation/matchedQualityCuts/";
    TString inFileName = "FlatTablebis_LHC22o_Global_KFPCA.root";
    TString outPath = inPath;
    TString outFileName = "histoToFit.root";

    createFitResFile(inPath, inFileName, outPath, outFileName);
    gSystem->Exit(0, kTRUE);
}

void createFitResFile(TString inPath, TString inFileName, TString outPath, TString outFileName)
{
    TFile *fin = TFile::Open(inPath + inFileName);
    TDirectoryFile *df = (TDirectoryFile *)fin->Get("dilepton-reader");
    auto hMass = (TH1D *)df->Get("Mass");

    std::cout << "Got all histos" << std::endl;

    auto res = TFile::Open(outPath + outFileName, "RECREATE");
    hMass->Write("Mass");

    std::cout << "Histos written" << std::endl;
    res->Close();
}