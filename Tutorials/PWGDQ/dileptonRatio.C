#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>

void dileptonRatio()
{
  TFile* Data_KF = TFile::Open("/Users/emiliebarreau/alice/O2Physics/Tutorials/PWGDQ/Histos_pp1_KF.root");
  TFile* Data_PCA = TFile::Open("/Users/emiliebarreau/alice/O2Physics/Tutorials/PWGDQ/Histos_pp1_MassKFV.root");

  TDirectory* Data_KF_D = (TDirectory*)Data_KF->GetDirectory("dilepton-reader");
  TDirectory* Data_PCA_D = (TDirectory*)Data_PCA->GetDirectory("dilepton-reader");

  TH1F* Data_KF_Tauz = static_cast<TH1F*>(Data_KF_D->Get("Mass"));
  TH1F* Data_PCA_Tauz = static_cast<TH1F*>(Data_PCA_D->Get("Mass"));

  TFile* MC_KF = TFile::Open("/Users/emiliebarreau/alice/O2Physics/Tutorials/PWGDQ/Histos_MC1_KF.root");
  TFile* MC_PCA = TFile::Open("/Users/emiliebarreau/alice/O2Physics/Tutorials/PWGDQ/Histos_MC1_PCA.root");

  TCanvas* c = new TCanvas();
  c->SetWindowSize(1000, 1000);
  c->cd(1);
  Data_KF_Tauz->SetLineColor(2);
  Data_KF_Tauz->Draw("HIST E");
  Data_PCA_Tauz->Draw("SAME HIST E");
  auto rp = new TRatioPlot(Data_PCA_Tauz, Data_KF_Tauz);
  rp->Draw();
  rp->GetLowYaxis()->SetNdivisions(5);
  rp->GetUpperPad()->cd();
  TLegend* legend = new TLegend(0.15, 0.9, 0.35, 0.8);
  legend->AddEntry(Data_KF_Tauz, "KFVertexing", "le");
  legend->AddEntry(Data_PCA_Tauz, "GetMass()", "le");
  legend->Draw();
}