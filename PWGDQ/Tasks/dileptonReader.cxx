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

using namespace o2;
using namespace o2::framework;

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
    AxisSpec tauzAxis = {200, -0.9, 0.9, "#tau_{z} (ns)"};
    AxisSpec ptAxis = {120, 0.0, 30.0, "p_{T} (GeV/c)"};
    AxisSpec pzAxis = {300, 0.0, 100.0, "p_{z} (GeV/c)"};
    AxisSpec SVAxis = {65, -80.0, 50.0, "Secondary_Vertex (cm)"}; // 900 -80 80
    AxisSpec rapAxis = {200, 2., 4., "y"};
    AxisSpec rabsAxis = {100, 0., 100., "R_{abs}"};
    AxisSpec etaAxis = {350, -8., -1., "#eta"};
    AxisSpec mcMaskAxis = {130, -1., 129., "FakeMatch"};
    AxisSpec IndexMcTrackAxis = {250001, -1., 250000., "IndexReducedMCTracks"};
    AxisSpec PdgCodeAxis = {900, -450., 450., "PdgCode"};

    HistogramConfigSpec ptSpec({HistType::kTH1F, {ptAxis}});
    HistogramConfigSpec massSpec({HistType::kTH1F, {massAxis}});
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

    // add some histograms to the registry
    registry.add("Pt", "Pt distribution", ptSpec);
    registry.add("Mass", "Invariant mass distribution", massSpec);
    // registry.add("massChi2TH3", "Mass and MCH-MFT Chi2 TH3 Histogram", histospec);
    registry.add("Pz", "Pz distribution", pzSpec);
    registry.add("Tauz", "Tauz distribution", tauzSpec);
    registry.add("Rapidity", "Rapidity distribution", rapSpec);
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
      registry.get<TH1>(HIST("Mass"))->Fill(dimuon.mass());
      registry.get<TH1>(HIST("Pt"))->Fill(dimuon.pt());
      registry.get<TH1>(HIST("Pz"))->Fill(dimuon.vertexPz());
      registry.get<TH1>(HIST("SecondVertex"))->Fill(dimuon.sVertex());
      registry.get<TH1>(HIST("Tauz"))->Fill(dimuon.tauz());
      registry.get<TH1>(HIST("DeltaZ"))->Fill(dimuon.posZ() - dimuon.sVertex());
      registry.get<TH1>(HIST("PrimaryVertex"))->Fill(dimuon.posZ());
      registry.get<TH1>(HIST("Eta"))->Fill(dimuon.eta());
      registry.get<TH1>(HIST("Rapidity"))->Fill(rap);
      registry.get<TH1>(HIST("FakeMatches1"))->Fill(dimuon.mcMask1());
      registry.get<TH1>(HIST("FakeMatches2"))->Fill(dimuon.mcMask2());
      registry.get<TH1>(HIST("Eta1"))->Fill(dimuon.eta1());
      registry.get<TH1>(HIST("Eta2"))->Fill(dimuon.eta2());

      if (dimuon.pdgCode1() == 443 && dimuon.pdgCode2() == 443) {
        if (dimuon.mcMask1() < 1. && dimuon.mcMask2() < 1.) { // GM & GM Jpsi
          registry.get<TH1>(HIST("Sa_SV"))->Fill(dimuon.sVertex());
        } else if (dimuon.mcMask1() > 1. && dimuon.mcMask2() > 1.) { // FM & FM Jpsi
          registry.get<TH1>(HIST("Sc_SV"))->Fill(dimuon.sVertex());
        } else if ((dimuon.mcMask1() < 1. && dimuon.mcMask2() > 1.) || (dimuon.mcMask1() > 1. && dimuon.mcMask2() < 1.)) { // GM & FM Jpsi
          registry.get<TH1>(HIST("Sb_SV"))->Fill(dimuon.sVertex());
        } else {
          continue;
        }
      } else {
        registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(dimuon.pdgCode1());
        registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(dimuon.pdgCode2());
        if (dimuon.mcMask1() < 1. && dimuon.mcMask2() < 1.) { // GM & GM X
          registry.get<TH1>(HIST("Ba_SV"))->Fill(dimuon.sVertex());
        } else if (dimuon.mcMask1() > 1. && dimuon.mcMask2() > 1.) { // FM & FM X
          registry.get<TH1>(HIST("Bc_SV"))->Fill(dimuon.sVertex());
        } else if ((dimuon.mcMask1() < 1. && dimuon.mcMask2() > 1.) || (dimuon.mcMask1() > 1. && dimuon.mcMask2() < 1.)) { // GM & FM X
          registry.get<TH1>(HIST("Bb_SV"))->Fill(dimuon.sVertex());
        } else {
          continue;
        }
      }

      // Global tracks studies
      if ((dimuon.eta1() > -3.6 && dimuon.eta1() < -2.5) && (dimuon.eta2() > -3.6 && dimuon.eta2() < -2.5)) { // cut on eta in mft acceptance
                                                                                                              // if (dimuon.chi2MatchMCHMFT1() <= chi2Cut && dimuon.chi2MatchMCHMFT2() <= chi2Cut) {                   // if (dimuon.chi2MatchMCHMFT1() <= chi2Cut && dimuon.chi2MatchMCHMFT2() <= chi2Cut) {
        if (!(dimuon.isAmbig1()) && !(dimuon.isAmbig2())) {                                                   // remove ambiguous tracks
          // hRapidityGlobal->Fill(rap);
          // registry.get<TH3>(HIST("massChi2TH3"))->Fill(dimuon.mass(), dimuon.chi2MatchMCHMFT1(), dimuon.chi2MatchMCHMFT2());
          if (dimuon.sign() == 0) {
            if (rap > 2.5 && rap < 3.6) {
              if (dimuon.mass() > 1.8) {

                /*if (dimuon.pdgCode1() == 443 && dimuon.pdgCode2() == 443) {
                  if (dimuon.mcMask1() < 1. && dimuon.mcMask2() < 1.) { // GM & GM Jpsi
                    registry.get<TH1>(HIST("Sa_SV"))->Fill(dimuon.sVertex());
                  } else if (dimuon.mcMask1() > 1. && dimuon.mcMask2() > 1.) { // FM & FM Jpsi
                    registry.get<TH1>(HIST("Sc_SV"))->Fill(dimuon.sVertex());
                  } else if ((dimuon.mcMask1() < 1. && dimuon.mcMask2() > 1.) || (dimuon.mcMask1() > 1. && dimuon.mcMask2() < 1.)) { // GM & FM Jpsi
                    registry.get<TH1>(HIST("Sb_SV"))->Fill(dimuon.sVertex());
                  } else {
                    continue;
                  }
                } else {
                  registry.get<TH1>(HIST("PdgCode_mothers1"))->Fill(dimuon.pdgCode1());
                  registry.get<TH1>(HIST("PdgCode_mothers2"))->Fill(dimuon.pdgCode2());
                  if (dimuon.mcMask1() < 1. && dimuon.mcMask2() < 1.) { // GM & GM X
                    registry.get<TH1>(HIST("Ba_SV"))->Fill(dimuon.sVertex());
                  } else if (dimuon.mcMask1() > 1. && dimuon.mcMask2() > 1.) { // FM & FM X
                    registry.get<TH1>(HIST("Bc_SV"))->Fill(dimuon.sVertex());
                  } else if ((dimuon.mcMask1() < 1. && dimuon.mcMask2() > 1.) || (dimuon.mcMask1() > 1. && dimuon.mcMask2() < 1.)) { // GM & FM X
                    registry.get<TH1>(HIST("Bb_SV"))->Fill(dimuon.sVertex());
                  } else {
                    continue;
                  }
                }*/
              }
            } // end rapidity cut
          } // end sign selection
        } // ambig cut
        //} // end MCH-MFT chi2 selection
      } // end eta cut
    } // end loop over dimuons
  };

  void processbis(soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsLabels> const& muons, aod::ReducedMCTracks const& mctracks)
  {
    for (auto& muon : muons) {
      auto muontrack = mctracks.rawIteratorAt(muon.reducedMCTrackId());
      if (muon.eta() < -2.5 && muon.eta() > -3.6) {
        if (muon.chi2MatchMCHMFT() <= chi2Cut) {
          if (!(muon.isAmbiguous())) {
            if (muontrack.has_mothers()) {
              auto motherId = muontrack.mothersIds()[0];
              auto mother = mctracks.rawIteratorAt(motherId);
              auto motherPdg = mother.pdgCode();
              registry.get<TH1>(HIST("PdgCode_mothers"))->Fill(motherPdg);
              if (muon.mcMask() <= 1.) { // good matches
                registry.get<TH1>(HIST("IndexMCTrack_GM"))->Fill(muon.reducedMCTrackId());
                if (motherPdg == 443) {
                  registry.get<TH1>(HIST("IndexMCTrack_Jpsi_GM"))->Fill(muon.reducedMCTrackId());
                  registry.get<TH1>(HIST("pT_Jpsi_GM"))->Fill(muon.pt());
                  registry.get<TH1>(HIST("eta_Jpsi_GM"))->Fill(muon.eta());
                  registry.get<TH1>(HIST("Rabs_Jpsi_GM"))->Fill(muon.rAtAbsorberEnd());
                } else if (motherPdg == 211 || motherPdg == -211) {
                  registry.get<TH1>(HIST("IndexMCTrack_Pion_GM"))->Fill(muon.reducedMCTrackId());
                  registry.get<TH1>(HIST("pT_Pion_GM"))->Fill(muon.pt());
                  registry.get<TH1>(HIST("eta_Pion_GM"))->Fill(muon.eta());
                  registry.get<TH1>(HIST("Rabs_Pion_GM"))->Fill(muon.rAtAbsorberEnd());
                }
              } else if (muon.mcMask() >= 128.) { // fake matches
                registry.get<TH1>(HIST("IndexMCTrack_FM"))->Fill(muon.reducedMCTrackId());
                if (motherPdg == 443) {
                  registry.get<TH1>(HIST("IndexMCTrack_Jpsi_FM"))->Fill(muon.reducedMCTrackId());
                  registry.get<TH1>(HIST("pT_Jpsi_FM"))->Fill(muon.pt());
                  registry.get<TH1>(HIST("eta_Jpsi_FM"))->Fill(muon.eta());
                  registry.get<TH1>(HIST("Rabs_Jpsi_FM"))->Fill(muon.rAtAbsorberEnd());
                } else if (motherPdg == 211 || motherPdg == -211) {
                  registry.get<TH1>(HIST("IndexMCTrack_Pion_FM"))->Fill(muon.reducedMCTrackId());
                  registry.get<TH1>(HIST("pT_Pion_FM"))->Fill(muon.pt());
                  registry.get<TH1>(HIST("eta_Pion_FM"))->Fill(muon.eta());
                  registry.get<TH1>(HIST("Rabs_Pion_FM"))->Fill(muon.rAtAbsorberEnd());
                }
              } else {
                continue;
              } // end mcmask selection
            }
          }
        }
      }
    }
    for (auto& mctrack : mctracks) {
      if (mctrack.has_mothers()) {                                 // only mctracks with mothers
        if (mctrack.pdgCode() == 13 || mctrack.pdgCode() == -13) { // only for muons
          /*auto motherId = mctrack.mothersIds()[0];
          auto mother = mctracks.rawIteratorAt(motherId);
          auto motherPdg = mother.pdgCode();
          if (motherPdg == 443) { // separate muons from Jpsi
            registry.get<TH1>(HIST("IndexMCTrack_Jpsi"))->Fill(mctrack.index());
          } else if (motherPdg == -211 || motherPdg == 211) { // and muons from pions
            registry.get<TH1>(HIST("IndexMCTrack_Pion"))->Fill(mctrack.index());
          } else {
            registry.get<TH1>(HIST("IndexMCTrack_Badmuon"))->Fill(mctrack.index());
          }*/
        }
        // fill only muons from J/Psi, do the same with other mothers to separate contributions
        // recuperer l'ID du muon
        // parcourir la table des GM/FM et si ID same, récupérer flag GM ou FM -> TTree ou vecteur ?
      }
    }
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dileptonReader>(cfgc)};
};
