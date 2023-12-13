#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>

void HistosRatio()
{
  // ================== MONTE CARLO FILES ==================
  TFile* MC_KF = TFile::Open("/Users/emiliebarreau/alice/O2Physics/Tutorials/PWGDQ/AnalysisResults_MC1_matchedQC_KFVertexing.root");
  TFile* MC_PCA = TFile::Open("/Users/emiliebarreau/alice/O2Physics/Tutorials/PWGDQ/AnalysisResults_MC1_matchedQC_proptoPCA.root");

  TDirectory* MC_KF_D = (TDirectory*)MC_KF->GetDirectory("analysis-same-event-pairing");
  TDirectory* MC_PCA_D = (TDirectory*)MC_PCA->GetDirectory("analysis-same-event-pairing");

  THashList* MC_KF_H = (THashList*)MC_KF_D->Get("output");
  THashList* MC_PCA_H = (THashList*)MC_PCA_D->Get("output");

  TList* MC_KF_L = (TList*)MC_KF_H->FindObject("PairsMuonSEPM_matchedQualityCuts");
  TList* MC_PCA_L = (TList*)MC_PCA_H->FindObject("PairsMuonSEPM_matchedQualityCuts");

  // Mass
  TH1F* MC_KF_Mass = static_cast<TH1F*>(MC_KF_L->FindObject("Mass"));
  TH1F* MC_PCA_Mass = static_cast<TH1F*>(MC_PCA_L->FindObject("Mass"));
  // MC_KF_Mass->Divide(MC_PCA_Mass);
  // Pz
  TH1F* MC_KF_Pz = static_cast<TH1F*>(MC_KF_L->FindObject("Pz"));
  TH1F* MC_PCA_Pz = static_cast<TH1F*>(MC_PCA_L->FindObject("Pz"));
  // PV
  TH1F* MC_KF_PV = static_cast<TH1F*>(MC_KF_L->FindObject("Primary_Vertexing_z"));
  TH1F* MC_PCA_PV = static_cast<TH1F*>(MC_PCA_L->FindObject("Primary_Vertexing_z"));
  // SV
  TH1F* MC_KF_SV = static_cast<TH1F*>(MC_KF_L->FindObject("Secondary_Vertexing_z"));
  TH1F* MC_PCA_SV = static_cast<TH1F*>(MC_PCA_L->FindObject("Secondary_Vertexing_z"));
  // Tauz
  TH1F* MC_KF_Tauz = static_cast<TH1F*>(MC_KF_L->FindObject("Tauz"));
  TH1F* MC_PCA_Tauz = static_cast<TH1F*>(MC_PCA_L->FindObject("Tauz"));

  // ================== DATA PP FILES ==================
  TFile* Data_KF = TFile::Open("/Users/emiliebarreau/alice/O2Physics/Tutorials/PWGDQ/AnalysisResults_pp1_matchedQC_KFVertexing.root");
  TFile* Data_PCA = TFile::Open("/Users/emiliebarreau/alice/O2Physics/Tutorials/PWGDQ/AnalysisResults_pp1_matchedQC_proptoPCA.root");

  TDirectory* Data_KF_D = (TDirectory*)Data_KF->GetDirectory("analysis-same-event-pairing");
  TDirectory* Data_PCA_D = (TDirectory*)Data_PCA->GetDirectory("analysis-same-event-pairing");

  THashList* Data_KF_H = (THashList*)Data_KF_D->Get("output");
  THashList* Data_PCA_H = (THashList*)Data_PCA_D->Get("output");

  TList* Data_KF_L = (TList*)Data_KF_H->FindObject("PairsMuonSEPM_matchedQualityCuts");
  TList* Data_PCA_L = (TList*)Data_PCA_H->FindObject("PairsMuonSEPM_matchedQualityCuts");

  // Mass
  TH1F* Data_KF_Mass = static_cast<TH1F*>(Data_KF_L->FindObject("Mass"));
  TH1F* Data_PCA_Mass = static_cast<TH1F*>(Data_PCA_L->FindObject("Mass"));
  // Pz
  TH1F* Data_KF_Pz = static_cast<TH1F*>(Data_KF_L->FindObject("Pz"));
  TH1F* Data_PCA_Pz = static_cast<TH1F*>(Data_PCA_L->FindObject("Pz"));
  // PV
  TH1F* Data_KF_PV = static_cast<TH1F*>(Data_KF_L->FindObject("Primary_Vertexing_z"));
  TH1F* Data_PCA_PV = static_cast<TH1F*>(Data_PCA_L->FindObject("Primary_Vertexing_z"));
  // SV
  TH1F* Data_KF_SV = static_cast<TH1F*>(Data_KF_L->FindObject("Secondary_Vertexing_z"));
  TH1F* Data_PCA_SV = static_cast<TH1F*>(Data_PCA_L->FindObject("Secondary_Vertexing_z"));
  // Tauz
  TH1F* Data_KF_Tauz = static_cast<TH1F*>(Data_KF_L->FindObject("Tauz"));
  TH1F* Data_PCA_Tauz = static_cast<TH1F*>(Data_PCA_L->FindObject("Tauz"));

  TCanvas* c = new TCanvas();
  c->SetWindowSize(1000, 1000);

  c->cd(1);
  c->SetTitle("Tauz of dimuons for data");
  // Data_KF_Tauz->GetXaxis()->SetTitle("Tauz");
  Data_KF_PV->SetLineColor(2);
  Data_KF_PV->Draw("HIST E");
  Data_KF_SV->Draw("SAME HIST E");
  auto rp = new TRatioPlot(Data_KF_PV, Data_KF_SV);
  rp->Draw();
  rp->GetLowYaxis()->SetNdivisions(5);
  rp->GetLowYaxis()->SetRange(0, 5);
  rp->GetUpperPad()->cd();
  TLegend* legend = new TLegend(0.15, 0.9, 0.2, 0.8);
  legend->AddEntry(Data_KF_PV, "PV", "le");
  legend->AddEntry(Data_KF_SV, "SV", "le");
  legend->Draw();
}