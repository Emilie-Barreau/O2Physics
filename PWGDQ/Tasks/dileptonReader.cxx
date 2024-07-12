// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include <iostream>
#include <vector>
#include "TLorentzVector.h"
#include "TMath.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedMCEventLabels>;
// using MyMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;

// quick task to read typical dileptonAod.root
// for now it only runs on the aod::DimuonsAll table but it can run on any table if you give it one

struct dileptonReader {

  // histograms can be created with OutputObj<TH1F>
  // OutputObj<TH1F> hTauz{TH1F("Tauz", "Tauz", 100, -0.01, 0.01)};
  // OutputObj<TH1F> hRapidityGlobal{TH1F("RapidityGlobal", "RapidityGlobal", 267, 2., 4.)};

  Configurable<Double_t> ptCut{"ptCut", 0.5, "pT cut, default = 0.5"};        // configurable pT cut
  Configurable<Double_t> chi2Cut{"chi2Cut", 0.0, "chi2T cut, default = 0.0"}; // configurable chi2 cut

  // Create registry for histograms (several ways to create histograms)
  HistogramRegistry registry{
    "registry",
    //{{"Mass", "Dimuon mass distribution", {HistType::kTH1F, {{300, 0, 15, "mass (GeV/c^{2})"}}}}}
  };

  void init(o2::framework::InitContext&)
  {
    // define axis for your histograms
    AxisSpec chi2Axis1 = {200, 0, 200, "#chi^2_{MCH-MFT} for mu 1"};
    AxisSpec chi2Axis2 = {200, 0, 200, "#chi^2_{MCH-MFT} for mu 2"};
    AxisSpec massAxis = {750, 0, 15, "dimuon mass (GeV/c^{2})"}; // 750
    AxisSpec tauzAxis = {500, -0.01, 0.01, "#tau_{z} (ns)"};
    AxisSpec ptAxis = {120, 0.0, 30.0, "p_{T} (GeV/c)"};
    AxisSpec pzAxis = {300, 0.0, 100.0, "p_{z} (GeV/c)"};
    AxisSpec SVAxis = {650, -80.0, 50.0, "Secondary_Vertex (cm)"}; // 900 -80 80
    AxisSpec rapAxis = {200, 2., 4., "y"};
    AxisSpec rabsAxis = {100, 0., 100., "R_{abs}"};
    AxisSpec etaAxis = {350, -8., -1., "#eta"};
    AxisSpec mcMaskAxis = {130, -1., 129., "FakeMatch"};
    AxisSpec IndexMcTrackAxis = {250001, -1., 250000., "IndexReducedMCTracks"};
    AxisSpec PdgCodeAxis = {900, -450., 450., "PdgCode"};

    HistogramConfigSpec ptSpec({HistType::kTH1F, {ptAxis}});
    HistogramConfigSpec massSpec({HistType::kTH1F, {massAxis}});
    HistogramConfigSpec chi2Spec({HistType::kTH1F, {chi2Axis1}});
    // HistogramConfigSpec histospec({HistType::kTH3F, {massAxis, chi2Axis1, chi2Axis2}});
    HistogramConfigSpec pzSpec({HistType::kTH1F, {pzAxis}});
    HistogramConfigSpec tauzSpec({HistType::kTH1F, {tauzAxis}});
    HistogramConfigSpec SVSpec({HistType::kTH1F, {SVAxis}});
    HistogramConfigSpec rapSpec({HistType::kTH1F, {rapAxis}});
    HistogramConfigSpec rabsSpec({HistType::kTH1F, {rabsAxis}});
    HistogramConfigSpec etaSpec({HistType::kTH1F, {etaAxis}});
    HistogramConfigSpec mcMaskSpec({HistType::kTH1F, {mcMaskAxis}});
    HistogramConfigSpec indexMcTrackSpec({HistType::kTH1F, {IndexMcTrackAxis}});
    HistogramConfigSpec PdgCodeSpec({HistType::kTH1F, {PdgCodeAxis}});
    HistogramConfigSpec DeltaTauzSpec({HistType::kTH2F, {tauzAxis, tauzAxis}});
    HistogramConfigSpec SVVzSpec({HistType::kTH2F, {SVAxis, SVAxis}});

    // add some histograms to the registry
    registry.add("Pt", "Pt distribution", ptSpec);
    registry.add("Ptbis", "Pt distribution", ptSpec);

    registry.add("Mass", "Invariant mass distribution", massSpec);
    registry.add("Mass_Sa", "Invariant mass distribution", massSpec);
    registry.add("Mass_Sb", "Invariant mass distribution", massSpec);
    registry.add("Mass_Sc", "Invariant mass distribution", massSpec);
    registry.add("Mass_Ba", "Invariant mass distribution", massSpec);
    registry.add("Mass_Bb", "Invariant mass distribution", massSpec);
    registry.add("Mass_Bc", "Invariant mass distribution", massSpec);
    // registry.add("massChi2TH3", "Mass and MCH-MFT Chi2 TH3 Histogram", histospec);
    registry.add("Pz", "Pz distribution", pzSpec);
    registry.add("Pz_Sa", "Pz distribution", pzSpec);
    registry.add("Pz_Sb", "Pz distribution", pzSpec);
    registry.add("Pz_Sc", "Pz distribution", pzSpec);
    registry.add("Pz_Ba", "Pz distribution", pzSpec);
    registry.add("Pz_Bb", "Pz distribution", pzSpec);
    registry.add("Pz_Bc", "Pz distribution", pzSpec);

    registry.add("Chi2_GM1", "Chi2 GM 1 distribution", chi2Spec);
    registry.add("Chi2_FM1", "Chi2 FM 1 distribution", chi2Spec);
    registry.add("Chi2_GM2", "Chi2 GM 2 distribution", chi2Spec);
    registry.add("Chi2_FM2", "Chi2 FM 2 distribution", chi2Spec);

    registry.add("Tauz", "Tauz distribution", tauzSpec);
    registry.add("Tauz_Sa", "Tauz signal GM Jpsi distribution", tauzSpec);
    registry.add("Tauz_Sb", "Tauz signal GM/FM Jpsi distribution", tauzSpec);
    registry.add("Tauz_Sc", "Tauz signal FM Jpsi distribution", tauzSpec);
    registry.add("Tauz_Ba", "Tauz background GM X distribution", tauzSpec);
    registry.add("Tauz_Bb", "Tauz background GM/FM X distribution", tauzSpec);
    registry.add("Tauz_Bc", "Tauz background FM X distribution", tauzSpec);
    registry.add("Delta_Tauz", "(Tauz - TauzMC) distribution", tauzSpec);
    registry.add("TauzMC", "TauzMC distribution", tauzSpec);
    registry.add("TauzMC1", "TauzMC1 distribution", tauzSpec);
    registry.add("TauzMC2", "TauzMC2 distribution", tauzSpec);
    registry.add("Posz", "Posz distribution", SVSpec);
    registry.add("Vz", "Vz distribution", SVSpec);
    registry.add("DeltaMC", "DeltaMC distribution", SVSpec);
    registry.add("PzMC", "Pz MC distribution", pzSpec);
    registry.add("MassMC", "Mass MC distribution", massSpec);

    registry.add("SVVz", "SV - Vz distribution", SVVzSpec);
    registry.add("TauzTauzMC", "Tauz - TauzMC distribution", DeltaTauzSpec);

    registry.add("Rapidity", "Rapidity distribution", rapSpec);
    registry.add("Rapiditybis", "Rapidity distribution", rapSpec);
    registry.add("DeltaZ", "PosZ - SVertex distribution", SVSpec);
    registry.add("PrimaryVertex", "Primary Vertex distribution", SVSpec);

    registry.add("SecondVertex", "Second Vertex distribution", SVSpec);
    registry.add("Sa_SV", "Second Vertex distribution for GM Jpsi", SVSpec);
    registry.add("Sb_SV", "Second Vertex distribution for GM/FM Jpsi", SVSpec);
    registry.add("Sc_SV", "Second Vertex distribution for FM Jpsi", SVSpec);
    registry.add("Ba_SV", "Second Vertex distribution for GM X", SVSpec);
    registry.add("Bb_SV", "Second Vertex distribution for GM/FM X", SVSpec);
    registry.add("Bc_SV", "Second Vertex distribution for FM X", SVSpec);

    registry.add("FakeMatches1", "Fake matches1", mcMaskSpec);
    registry.add("FakeMatches2", "Fake matches2", mcMaskSpec);

    registry.add("Eta", "Pseudo-rapidity distribution", etaSpec);
    registry.add("Eta1", "Pseudo-rapidity distribution of muon 1", etaSpec);
    registry.add("Eta2", "Pseudo-rapidity distribution of muon 2", etaSpec);

    registry.add("IndexMCTrack_GM", "Index ReducedMCTrack Good", indexMcTrackSpec);
    registry.add("IndexMCTrack_FM", "Index ReducedMCTrack Fake", indexMcTrackSpec);
    registry.add("IndexMCTrack_MCMask", "Good/Fake Matches", mcMaskSpec);

    registry.add("PdgCode_mothers", "PdgCode mothers", PdgCodeSpec);
    registry.add("PdgCode_mothers1", "PdgCode mothers of muon 1", PdgCodeSpec);
    registry.add("PdgCode_mothers2", "PdgCode mothers of muon 2", PdgCodeSpec);

    registry.add("IndexMCTrack_Jpsi_GM", "Muons from Jpsi Good Matches", indexMcTrackSpec);
    registry.add("IndexMCTrack_Jpsi_FM", "Muons from Jpsi Fake Matches", indexMcTrackSpec);
    registry.add("IndexMCTrack_Pion_GM", "Muons from Pion Good Matches", indexMcTrackSpec);
    registry.add("IndexMCTrack_Pion_FM", "Muons from Pion Fake Matches", indexMcTrackSpec);

    registry.add("pT_Jpsi_GM", "pT mu from Jpsi Good Matches", ptSpec);
    registry.add("pT_Jpsi_FM", "pT mu from Jpsi Fake Matches", ptSpec);
    registry.add("pT_Pion_GM", "pT mu from Pion Good Matches", ptSpec);
    registry.add("pT_Pion_FM", "pT mu from Pion Fake Matches", ptSpec);

    registry.add("eta_Jpsi_GM", "eta mu from Jpsi Good Matches", etaSpec);
    registry.add("eta_Jpsi_FM", "eta mu from Jpsi Fake Matches", etaSpec);
    registry.add("eta_Pion_GM", "eta mu from Pion Good Matches", etaSpec);
    registry.add("eta_Pion_FM", "eta mu from Pion Fake Matches", etaSpec);

    registry.add("Rabs_Jpsi_GM", "Rabs mu from Jpsi Good Matches", rabsSpec);
    registry.add("Rabs_Jpsi_FM", "Rabs mu from Jpsi Fake Matches", rabsSpec);
    registry.add("Rabs_Pion_GM", "Rabs mu from Pion Good Matches", rabsSpec);
    registry.add("Rabs_Pion_FM", "Rabs mu from Pion Fake Matches", rabsSpec);
  };

  void process(aod::DimuonsAll const& dimuons)
  {
    // Double_t const PI = ROOT::Math::Pi(); // uncomment if you want to use pi, can be useful

    for (auto& dimuon : dimuons) {
      // calculate rapidity
      auto rap = -ROOT::Math::log((ROOT::Math::sqrt(dimuon.mass() * dimuon.mass() + dimuon.pt() * dimuon.pt() * ROOT::Math::cosh(dimuon.eta()) * ROOT::Math::cosh(dimuon.eta())) + dimuon.pt() * ROOT::Math::sinh(dimuon.eta())) / (ROOT::Math::sqrt(dimuon.mass() * dimuon.mass() + dimuon.pt() * dimuon.pt())));

      // MUON standalone tracks studies
      if ((dimuon.eta1() > -4 && dimuon.eta1() < -2.5) && (dimuon.eta2() > -4 && dimuon.eta2() < -2.5)) {
        if (dimuon.pt1() >= ptCut && dimuon.pt2() >= ptCut) {
          if (dimuon.sign() == 0) {
            if (dimuon.chi2MatchMCHMFT1() <= chi2Cut && dimuon.chi2MatchMCHMFT2() <= chi2Cut) {
              if (rap > 2.5 && rap < 4) {
                // fill any histo for MCH-MID tracks (create them first in the histogram registry OR define them befor Init() with OutputObj)
              } // end rapidity cut
            } // end MCH-MID track selection
          } // end OS selection
        } // end pt cut
      } // end eta cut (MUON standalone acceptance)

      // if (dimuon.mcMask1() < 1. && dimuon.mcMask2() < 1.) {}

      ROOT::Math::PtEtaPhiMVector v1MC(dimuon.ptMC1(), dimuon.etaMC1(), dimuon.phiMC1(), 0.105658);
      ROOT::Math::PtEtaPhiMVector v2MC(dimuon.ptMC2(), dimuon.etaMC2(), dimuon.phiMC2(), 0.105658);
      ROOT::Math::PtEtaPhiMVector v12MC = v1MC + v2MC;
      float Tauz1MC = (dimuon.mcPosZ() - dimuon.vz1()) * v12MC.M() / TMath::Abs(v12MC.Pz());
      float Tauz2MC = (dimuon.mcPosZ() - dimuon.vz2()) * v12MC.M() / TMath::Abs(v12MC.Pz());
      float realTauzMC = (Tauz1MC + Tauz2MC) / 2.;

      /*registry.get<TH1>(HIST("Tauz"))->Fill(dimuon.tauz());
      registry.get<TH1>(HIST("Mass"))->Fill(dimuon.mass());
      registry.get<TH1>(HIST("SecondVertex"))->Fill(dimuon.sVertex());
      registry.get<TH1>(HIST("Pt"))->Fill(dimuon.pt());
      registry.get<TH1>(HIST("Pz"))->Fill(dimuon.vertexPz());
      registry.get<TH1>(HIST("DeltaZ"))->Fill(dimuon.posZ() - dimuon.sVertex());
      registry.get<TH1>(HIST("PrimaryVertex"))->Fill(dimuon.posZ());
      registry.get<TH1>(HIST("Eta"))->Fill(dimuon.eta());
      registry.get<TH1>(HIST("Rapidity"))->Fill(rap);
      registry.get<TH1>(HIST("FakeMatches1"))->Fill(dimuon.mcMask1());
      registry.get<TH1>(HIST("FakeMatches2"))->Fill(dimuon.mcMask2());
      registry.get<TH1>(HIST("Eta1"))->Fill(dimuon.eta1());
      registry.get<TH1>(HIST("Eta2"))->Fill(dimuon.eta2());
      registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(dimuon.pdgCode1());
      registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(dimuon.pdgCode2());*/

      /*if (dimuon.pdgCode1() == 443 && dimuon.pdgCode2() == 443) {
        if (dimuon.mcMask1() < 1. && dimuon.mcMask2() < 1.) { // GM & GM Jpsi
          registry.get<TH1>(HIST("Sa_SV"))->Fill(dimuon.sVertex());
          registry.get<TH1>(HIST("Tauz_Sa"))->Fill(dimuon.tauz());
          registry.get<TH1>(HIST("Mass_Sa"))->Fill(dimuon.mass());
          registry.get<TH1>(HIST("Pz_Sa"))->Fill(dimuon.vertexPz());
          registry.get<TH1>(HIST("Chi2_GM1"))->Fill(dimuon.chi2MatchMCHMFT1());
          registry.get<TH1>(HIST("Chi2_GM2"))->Fill(dimuon.chi2MatchMCHMFT2());
        } else if (dimuon.mcMask1() > 1. && dimuon.mcMask2() > 1.) { // FM & FM Jpsi
          registry.get<TH1>(HIST("Sc_SV"))->Fill(dimuon.sVertex());
          registry.get<TH1>(HIST("Tauz_Sc"))->Fill(dimuon.tauz());
          registry.get<TH1>(HIST("Mass_Sc"))->Fill(dimuon.mass());
          registry.get<TH1>(HIST("Pz_Sc"))->Fill(dimuon.vertexPz());
          registry.get<TH1>(HIST("Chi2_FM1"))->Fill(dimuon.chi2MatchMCHMFT1());
          registry.get<TH1>(HIST("Chi2_FM2"))->Fill(dimuon.chi2MatchMCHMFT2());
        } else if ((dimuon.mcMask1() < 1. && dimuon.mcMask2() > 1.) || (dimuon.mcMask1() > 1. && dimuon.mcMask2() < 1.)) { // GM & FM Jpsi
          registry.get<TH1>(HIST("Sb_SV"))->Fill(dimuon.sVertex());
          registry.get<TH1>(HIST("Tauz_Sb"))->Fill(dimuon.tauz());
          registry.get<TH1>(HIST("Mass_Sb"))->Fill(dimuon.mass());
          registry.get<TH1>(HIST("Pz_Sb"))->Fill(dimuon.vertexPz());
        } else {
          continue;
        }
        registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(dimuon.pdgCode1());
        registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(dimuon.pdgCode2());
        if (dimuon.mcMask1() < 1. && dimuon.mcMask2() < 1.) { // GM & GM X
          registry.get<TH1>(HIST("Ba_SV"))->Fill(dimuon.sVertex());
          registry.get<TH1>(HIST("Tauz_Ba"))->Fill(dimuon.tauz());
          registry.get<TH1>(HIST("Mass_Ba"))->Fill(dimuon.mass());
          registry.get<TH1>(HIST("Pz_Ba"))->Fill(dimuon.vertexPz());
        } else if (dimuon.mcMask1() > 1. && dimuon.mcMask2() > 1.) { // FM & FM X
          registry.get<TH1>(HIST("Bc_SV"))->Fill(dimuon.sVertex());
          registry.get<TH1>(HIST("Tauz_Bc"))->Fill(dimuon.tauz());
          registry.get<TH1>(HIST("Mass_Bc"))->Fill(dimuon.mass());
          registry.get<TH1>(HIST("Pz_Bc"))->Fill(dimuon.vertexPz());
        } else if ((dimuon.mcMask1() < 1. && dimuon.mcMask2() > 1.) || (dimuon.mcMask1() > 1. && dimuon.mcMask2() < 1.)) { // GM & FM X
          registry.get<TH1>(HIST("Bb_SV"))->Fill(dimuon.sVertex());
          registry.get<TH1>(HIST("Tauz_Bb"))->Fill(dimuon.tauz());
          registry.get<TH1>(HIST("Mass_Bb"))->Fill(dimuon.mass());
          registry.get<TH1>(HIST("Pz_Bb"))->Fill(dimuon.vertexPz());
        } else {
          continue;
        }
      }*/

      // registry.get<TH1>(HIST("Mass"))->Fill(v12MC.M());
      //  if (dimuon.mcMask1() < 1. && dimuon.mcMask2() < 1.) {
      /*if (dimuon.pdgCode1() == 443 && dimuon.pdgCode2() == 443) {
        registry.get<TH1>(HIST("Mass_Sa"))->Fill(v12MC.M());
        registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(dimuon.pdgCode1());
        registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(dimuon.pdgCode2());
        registry.get<TH1>(HIST("Mass_Sb"))->Fill(dimuon.mass());
      } else {
        registry.get<TH1>(HIST("Mass_Ba"))->Fill(v12MC.M());
        registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(dimuon.pdgCode1());
        registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(dimuon.pdgCode2());
        registry.get<TH1>(HIST("Mass_Bb"))->Fill(dimuon.mass());
      }*/
      //}

      // Global tracks studies
      if ((dimuon.eta1() > -3.6 && dimuon.eta1() < -2.5) && (dimuon.eta2() > -3.6 && dimuon.eta2() < -2.5)) { // cut on eta in mft acceptance
        if (dimuon.chi2MatchMCHMFT1() <= chi2Cut && dimuon.chi2MatchMCHMFT2() <= chi2Cut) {                   // if (dimuon.chi2MatchMCHMFT1() <= chi2Cut && dimuon.chi2MatchMCHMFT2() <= chi2Cut) {
          if (!(dimuon.isAmbig1()) && !(dimuon.isAmbig2())) {                                                 // remove ambiguous tracks
            // hRapidityGlobal->Fill(rap);
            // registry.get<TH3>(HIST("massChi2TH3"))->Fill(dimuon.mass(), dimuon.chi2MatchMCHMFT1(), dimuon.chi2MatchMCHMFT2());
            if (dimuon.sign() == 0) {
              if (rap > 2.5 && rap < 3.6) {
                if (dimuon.mass() > 1.8) {
                  // if (dimuon.mass() > 2.8 && dimuon.mass() < 3.3) {
                  // if (dimuon.pt() > 2.) {
                  if ((dimuon.pdgCode1() == 443 && dimuon.pdgCode2() == 443)) {
                    registry.get<TH1>(HIST("Tauz"))->Fill(dimuon.tauz());
                    registry.get<TH1>(HIST("Mass"))->Fill(dimuon.mass());
                    registry.get<TH1>(HIST("SecondVertex"))->Fill(dimuon.sVertex());
                    registry.get<TH1>(HIST("Pt"))->Fill(dimuon.pt());
                    registry.get<TH1>(HIST("Pz"))->Fill(dimuon.vertexPz());
                    registry.get<TH1>(HIST("DeltaZ"))->Fill(dimuon.posZ() - dimuon.sVertex());
                    registry.get<TH1>(HIST("PrimaryVertex"))->Fill(dimuon.posZ());
                    registry.get<TH1>(HIST("Eta"))->Fill(dimuon.eta());
                    registry.get<TH1>(HIST("Rapidity"))->Fill(rap);
                    registry.get<TH1>(HIST("FakeMatches1"))->Fill(dimuon.mcMask1());
                    registry.get<TH1>(HIST("FakeMatches2"))->Fill(dimuon.mcMask2());
                    registry.get<TH1>(HIST("Eta1"))->Fill(dimuon.eta1());
                    registry.get<TH1>(HIST("Eta2"))->Fill(dimuon.eta2());
                  } else {
                    continue;
                  }

                  if (dimuon.pdgCode1() == 443 && dimuon.pdgCode2() == 443) {
                    if (dimuon.mcMask1() < 1. && dimuon.mcMask2() < 1.) { // GM & GM Jpsi
                      registry.get<TH1>(HIST("Sa_SV"))->Fill(dimuon.sVertex());
                      registry.get<TH1>(HIST("Tauz_Sa"))->Fill(dimuon.tauz());
                      registry.get<TH1>(HIST("Mass_Sa"))->Fill(dimuon.mass());
                      registry.get<TH1>(HIST("Pz_Sa"))->Fill(dimuon.vertexPz());
                      registry.get<TH1>(HIST("Chi2_GM1"))->Fill(dimuon.chi2MatchMCHMFT1());
                      registry.get<TH1>(HIST("Chi2_GM2"))->Fill(dimuon.chi2MatchMCHMFT2());

                      registry.get<TH1>(HIST("TauzMC"))->Fill(realTauzMC);
                      registry.get<TH1>(HIST("TauzMC1"))->Fill(Tauz1MC);
                      registry.get<TH1>(HIST("TauzMC2"))->Fill(Tauz2MC);
                      registry.get<TH1>(HIST("Posz"))->Fill(dimuon.mcPosZ());
                      registry.get<TH1>(HIST("Vz"))->Fill(dimuon.vz1());
                      registry.get<TH1>(HIST("DeltaMC"))->Fill(dimuon.mcPosZ() - dimuon.vz1());
                      registry.get<TH1>(HIST("PzMC"))->Fill(TMath::Abs(v12MC.Pz()));
                      // registry.get<TH1>(HIST("MassMC"))->Fill(v12MC.M());
                      registry.get<TH1>(HIST("Delta_Tauz"))->Fill(dimuon.tauz() - realTauzMC);

                      registry.get<TH2>(HIST("TauzTauzMC"))->Fill(dimuon.tauz(), realTauzMC);
                      registry.get<TH2>(HIST("SVVz"))->Fill(dimuon.sVertex(), dimuon.vz1());
                    } else if (dimuon.mcMask1() > 1. && dimuon.mcMask2() > 1.) { // FM & FM Jpsi
                      registry.get<TH1>(HIST("Sc_SV"))->Fill(dimuon.sVertex());
                      registry.get<TH1>(HIST("Tauz_Sc"))->Fill(dimuon.tauz());
                      registry.get<TH1>(HIST("Mass_Sc"))->Fill(dimuon.mass());
                      registry.get<TH1>(HIST("Pz_Sc"))->Fill(dimuon.vertexPz());
                      registry.get<TH1>(HIST("Chi2_FM1"))->Fill(dimuon.chi2MatchMCHMFT1());
                      registry.get<TH1>(HIST("Chi2_FM2"))->Fill(dimuon.chi2MatchMCHMFT2());
                    } else if ((dimuon.mcMask1() < 1. && dimuon.mcMask2() > 1.) || (dimuon.mcMask1() > 1. && dimuon.mcMask2() < 1.)) { // GM & FM Jpsi
                      registry.get<TH1>(HIST("Sb_SV"))->Fill(dimuon.sVertex());
                      registry.get<TH1>(HIST("Tauz_Sb"))->Fill(dimuon.tauz());
                      registry.get<TH1>(HIST("Mass_Sb"))->Fill(dimuon.mass());
                      registry.get<TH1>(HIST("Pz_Sb"))->Fill(dimuon.vertexPz());
                    } else {
                      continue;
                    }
                    //} else if (dimuon.pdgCode1() == 0 && dimuon.pdgCode2() == 0) {
                    // registry.get<TH1>(HIST("Tauz_Ba"))->Fill(dimuon.tauz());
                    // registry.get<TH1>(HIST("Mass_Ba"))->Fill(dimuon.mass());
                    // registry.get<TH1>(HIST("MassMC"))->Fill(v12MC.M());
                    // continue;
                  } else {
                    registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(dimuon.pdgCode1());
                    registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(dimuon.pdgCode2());
                    if (dimuon.mcMask1() < 1. && dimuon.mcMask2() < 1.) { // GM & GM X
                      registry.get<TH1>(HIST("Ba_SV"))->Fill(dimuon.sVertex());
                      registry.get<TH1>(HIST("Tauz_Ba"))->Fill(dimuon.tauz());
                      registry.get<TH1>(HIST("Mass_Ba"))->Fill(dimuon.mass());
                      registry.get<TH1>(HIST("Pz_Ba"))->Fill(dimuon.vertexPz());
                    } else if (dimuon.mcMask1() > 1. && dimuon.mcMask2() > 1.) { // FM & FM X
                      registry.get<TH1>(HIST("Bc_SV"))->Fill(dimuon.sVertex());
                      registry.get<TH1>(HIST("Tauz_Bc"))->Fill(dimuon.tauz());
                      registry.get<TH1>(HIST("Mass_Bc"))->Fill(dimuon.mass());
                      registry.get<TH1>(HIST("Pz_Bc"))->Fill(dimuon.vertexPz());
                    } else if ((dimuon.mcMask1() < 1. && dimuon.mcMask2() > 1.) || (dimuon.mcMask1() > 1. && dimuon.mcMask2() < 1.)) { // GM & FM X
                      registry.get<TH1>(HIST("Bb_SV"))->Fill(dimuon.sVertex());
                      registry.get<TH1>(HIST("Tauz_Bb"))->Fill(dimuon.tauz());
                      registry.get<TH1>(HIST("Mass_Bb"))->Fill(dimuon.mass());
                      registry.get<TH1>(HIST("Pz_Bb"))->Fill(dimuon.vertexPz());
                    } else {
                      continue;
                    }
                  }

                  //} // end pt cut
                } // end mass cut
              } // end rapidity cut
            } // end sign selection
          } // ambig cut
        } // end MCH-MFT chi2 selection
      } // end eta cut
    } // end loop over dimuons
  };

  Preslice<soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>> perCollisionMuons = aod::fwdtrack::collisionId;

  void processbis(soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsLabels> const& muons, aod::ReducedMCTracks const& mctracks) // MyEventsSelected::iterator const& event,
  {
    // auto groupedMuons = muons.sliceBy(perCollisionReducedMuons, event.globalIndex());

    for (auto& [muon1, muon2] : combinations(muons, muons)) {
      // for (auto& muon1 : muons) {
      if (!(muon1.has_reducedMCTrack()) || !(muon2.has_reducedMCTrack())) {
        continue;
      }
      auto muontrack1 = mctracks.rawIteratorAt(muon1.reducedMCTrackId()); // check si muon a pas de trace MC
      auto muontrack2 = mctracks.rawIteratorAt(muon2.reducedMCTrackId());
      auto motherPdg1 = 0;
      auto motherPdg2 = 0;
      if (muontrack1.has_mothers()) {
        auto motherId1 = muontrack1.mothersIds()[0];
        auto mother1 = mctracks.rawIteratorAt(motherId1);
        motherPdg1 = mother1.pdgCode();
      }
      if (muontrack2.has_mothers()) {
        auto motherId2 = muontrack2.mothersIds()[0];
        auto mother2 = mctracks.rawIteratorAt(motherId2);
        motherPdg2 = mother2.pdgCode();
      }
      // registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(motherPdg2);

      ROOT::Math::PtEtaPhiMVector v1MC(muon1.pt(), muon1.eta(), muon1.phi(), 0.105658);
      ROOT::Math::PtEtaPhiMVector v2MC(muon2.pt(), muon2.eta(), muon2.phi(), 0.105658);
      ROOT::Math::PtEtaPhiMVector v12MC = v1MC + v2MC;
      v12MC = v12MC;

      if (muon1.eta() < -2.5 && muon1.eta() > -3.6) {
        if (muon1.sign() == muon2.sign()) {
          continue;
        }
        registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg1);
        registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg2);
        if (motherPdg1 == 443 && motherPdg2 == 443) {
          registry.get<TH1>(HIST("Mass_Sa"))->Fill(v12MC.M());
          registry.get<TH1>(HIST("Rapidity"))->Fill(abs(v12MC.Rapidity()));
          registry.get<TH1>(HIST("Pt"))->Fill(abs(v12MC.Pt()));
        } else if (motherPdg1 == 0 && motherPdg2 == 0) {
          registry.get<TH1>(HIST("Mass_Ba"))->Fill(v12MC.M());
          registry.get<TH1>(HIST("Rapiditybis"))->Fill(abs(v12MC.Rapidity()));
          registry.get<TH1>(HIST("Ptbis"))->Fill(abs(v12MC.Pt()));
        }
      }
      // end eta selection
    }
  }; // end process muon

  void processbis(soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels> const& tracks, aod::McParticles_001 const& mcTracks, aod::Collisions const& collisions)
  {
    for (auto& collision : collisions) {
      auto groupedMuons = tracks.sliceBy(perCollisionMuons, collision.globalIndex());

      // for (auto& track : groupedMuons) {

      for (auto& [track, track2] : combinations(groupedMuons, groupedMuons)) {
        if (!track.has_mcParticle() || static_cast<int>(track.trackType()) != 0) {
          continue;
        }

        // for (auto& track2 : groupedMuons) {
        if (!track2.has_mcParticle() || static_cast<int>(track2.trackType()) != 0 || track.globalIndex() == track2.globalIndex() || track.sign() == track2.sign()) {
          continue;
        }

        auto trackbis = track.template mcParticle_as<aod::McParticles_001>();
        auto track2bis = track2.template mcParticle_as<aod::McParticles_001>();
        // auto trackPdg = trackbis.pdgCode();
        auto motherPdg = 0;
        auto motherPdg2 = 0;

        if (trackbis.has_mothers()) {
          auto motherId = trackbis.mothersIds()[0];
          auto mother = mcTracks.rawIteratorAt(motherId);
          motherPdg = mother.pdgCode();
        }
        if (track2bis.has_mothers()) {
          auto motherId2 = track2bis.mothersIds()[0];
          auto mother2 = mcTracks.rawIteratorAt(motherId2);
          motherPdg2 = mother2.pdgCode();
        }

        registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg);
        registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(motherPdg2);

        ROOT::Math::PtEtaPhiMVector v1MC(track.pt(), track.eta(), track.phi(), 0.105658);
        ROOT::Math::PtEtaPhiMVector v2MC(track2.pt(), track2.eta(), track2.phi(), 0.105658);
        ROOT::Math::PtEtaPhiMVector v12MC = v1MC + v2MC;
        v12MC = v12MC;
        // registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg1);
        // registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg2);
        /*if (motherPdg == 443 && motherPdg2 == 443) {
          registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg);
          registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg2);
          registry.get<TH1>(HIST("Mass_Sa"))->Fill(v12MC.M());
          registry.get<TH1>(HIST("Rapidity"))->Fill(abs(v12MC.Rapidity()));
        } else {
          registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(motherPdg);
          registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(motherPdg2);
          registry.get<TH1>(HIST("Mass_Ba"))->Fill(v12MC.M());
          registry.get<TH1>(HIST("Rapiditybis"))->Fill(abs(v12MC.Rapidity()));
        }*/
        //}
      }
    }
    /*for (auto& collision : collisions) {
      auto groupedMuons = tracks.sliceBy(perCollisionMuons, collision.globalIndex());
      for (auto& [track1, track2] : combinations(groupedMuons, groupedMuons)) {
        if (!track1.has_mcParticle() || !track2.has_mcParticle()) {
          continue;
        }
        auto muontrack1 = track1.template mcParticle_as<aod::McParticles_001>();
        auto muontrack2 = track2.template mcParticle_as<aod::McParticles_001>();

        auto muonPdg1 = muontrack1.pdgCode();
        auto muonPdg2 = muontrack2.pdgCode();

        auto motherPdg1 = 0;
        auto motherPdg2 = 0;

        if (muontrack1.has_mothers()) {
          auto motherId1 = muontrack1.mothersIds()[0];
          auto mother1 = mcTracks.rawIteratorAt(motherId1);
          motherPdg1 = mother1.pdgCode();
        }
        if (muontrack2.has_mothers()) {
          auto motherId2 = muontrack2.mothersIds()[0];
          auto mother2 = mcTracks.rawIteratorAt(motherId2);
          motherPdg2 = mother2.pdgCode();
        }

        if ((muonPdg1 == 13 && muonPdg2 == -13) || (muonPdg1 == -13 && muonPdg2 == 13)) {
          ROOT::Math::PtEtaPhiMVector v1MC(track1.pt(), track1.eta(), track1.phi(), 0.105658);
          ROOT::Math::PtEtaPhiMVector v2MC(track2.pt(), track2.eta(), track2.phi(), 0.105658);
          ROOT::Math::PtEtaPhiMVector v12MC = v1MC + v2MC;
          v12MC = v12MC;
          if (static_cast<int>(track1.trackType()) == 0 && static_cast<int>(track2.trackType()) == 0) {
            // registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg1);
            // registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg2);
            if (motherPdg1 == 443 && motherPdg2 == 443) {
              registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg1);
              registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(motherPdg2);
              registry.get<TH1>(HIST("Mass_Sa"))->Fill(v12MC.M());
              registry.get<TH1>(HIST("Rapidity"))->Fill(v12MC.Rapidity());
            }
          } // end of trackType selection
        } // end of pdg selection
      } // end of tracks loop
    }*/
    // end of collision loop
  }; // end process AOD
}; // end init

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dileptonReader>(cfgc)};
};