// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskCorrelationDplusDplusReduced.cxx
/// \brief Writer of pairs of D mesons candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug, local optimization of analysis on small samples or ML training.
///        In this file are defined and filled the output tables
///
/// \author Valerio DI BELLA <valerio.di.bella@cern.ch>, IPHC Strasbourg
/// Based on the code of Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/HFC/DataModel/ReducedDMesonPairsTables.h"

#include "Tools/ML/MlResponse.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static constexpr double defaultCutsMl[1][3] = {{0.5, 0.5, 0.5}};

struct HfTaskCorrelationDplusDplusReduced {
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 1, "Selection Flag for Dplus"};

  Configurable<bool> applyML{"applyML", 0, "Apply the ML if true"};
  // ML inference
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{1., 5000.}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{cuts_ml::CutSmaller, cuts_ml::CutNot, cuts_ml::CutNot}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {defaultCutsMl[0], 1, 3, {"pT bin 0"}, {"score prompt", "score non-prompt", "score bkg"}}, "ML selections per pT bin"};
  //Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)3, "Number of classes in ML model"};
  Configurable<int> nClassesMl{"nClassesMl", 3, "Number of classes in ML model"};
  // Model file names
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_DplusToPiKPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  // Bonus: CCDB configuration (needed for ML application on the GRID)
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", true, "Flag to enable or disable the loading of models from CCDB"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTSmearedDplus/"}, "Path on CCDB"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};

  //using SelectedCandidates = soa::Filtered<o2::aod::HfCandDpTinys>;
  using SelectedCandidates = soa::Filtered<o2::aod::HfCandDpFulls>;
  using SelectedMcParticles = o2::aod::HfCandDpMcPs;

  Filter filterSelectCandidates = aod::full::candidateSelFlag >= selectionFlagDplus;

  o2::ccdb::CcdbApi ccdbApi;
  
  HistogramConfigSpec hTH1NCand{HistType::kTH1F, {{7, -0.5, 6.5}}};
  HistogramConfigSpec hTH1NMcRec{HistType::kTH1F, {{7, -0.5, 6.5}}};
  HistogramConfigSpec hTH1NMcGen{HistType::kTH1F, {{7, -0.5, 6.5}}};
  HistogramRegistry registry{
    "registry",
    {{"hNCand", "Number of D candidates per event;N", hTH1NCand},
     {"hNMcRec", "Number of reconstructed Mc D mesons per event;N", hTH1NMcRec},
     {"hNMcGen", "Number of generated Mc D mesons per event;N", hTH1NMcGen}}};

  // Add objects needed for ML inference
  std::vector<float> outputMl = {};
  o2::analysis::MlResponse<float> mlResponse;

  void init(InitContext const&)
  {
    registry.add("hMassDplus", "D+ candidates;inv. mass (#pi#pi K) (GeV/#it{c}^{2}))", {HistType::kTH1F, {{120, 1.5848, 2.1848}}});
    registry.add("hMassDminus", "D- candidates;inv. mass (#pi#pi K) (GeV/#it{c}^{2}))", {HistType::kTH1F, {{120, 1.5848, 2.1848}}});
    registry.add("hMassDplusMatched", "D+ matched candidates;inv. mass (#pi#pi K) (GeV/#it{c}^{2}))", {HistType::kTH1F, {{120, 1.5848, 2.1848}}});
    registry.add("hMassDminusMatched", "D- matched candidates;inv. mass (#pi#pi K) (GeV/#it{c}^{2}))", {HistType::kTH1F, {{120, 1.5848, 2.1848}}});
    registry.add("hMassDplusminusPair", "D plus-minus pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});inv. mass (#pi K) (GeV/#it{c}^{2})", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {120, 1.5848, 2.1848}}});
    registry.add("hMassDplusPair", "D plus pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});inv. mass (#pi K) (GeV/#it{c}^{2})", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {120, 1.5848, 2.1848}}});
    registry.add("hMassDminusPair", "D minus pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});inv. mass (#pi K) (GeV/#it{c}^{2})", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {120, 1.5848, 2.1848}}});
    registry.add("hDltPhiMcGen", "Azimuthal correlation for D mesons; #Delta#phi", {HistType::kTH1F, {{100, -3.141593, 3.141593}}});

    registry.add("hPrompt", "Prompt score; Score", {HistType::kTH1F, {{100, 0, 1}}});
    registry.add("hNonPrompt", "Non-prompt score; Score", {HistType::kTH1F, {{100, 0, 1}}});
    registry.add("hBkg", "Background score; Score", {HistType::kTH1F, {{100, 0, 1}}});
    registry.add("hPromptMatched", "Prompt score matched; Score", {HistType::kTH1F, {{100, 0, 1}}});
    registry.add("hNonPromptMatched", "Non-prompt score matched; Score", {HistType::kTH1F, {{100, 0, 1}}});
    registry.add("hBkgMatched", "Background score matched; Score", {HistType::kTH1F, {{100, 0, 1}}});

    // Configure and initialise the ML class
    mlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);

    if (loadModelsFromCCDB) {
      ccdbApi.init(ccdbUrl);
      mlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB.value, timestampCCDB);
    } else {
      mlResponse.setModelPathsLocal(onnxFileNames);
    }

    mlResponse.init();

  }

  void processLocalData(o2::aod::HfCandDpFullEvs::iterator const&,
                        SelectedCandidates const& localCandidates)
  {
    registry.fill(HIST("hNCand"), localCandidates.size());

    for (const auto& cand1 : localCandidates) {
      auto mass1 = cand1.m();
      auto sign1 = 1;
      if (cand1.pt() < 0) {
        sign1 = -1;
        registry.fill(HIST("hMassDminus"), mass1);
      } else {
        registry.fill(HIST("hMassDplus"), mass1);
      }

      for (auto cand2 = cand1 + 1; cand2 != localCandidates.end(); ++cand2) {
        auto mass2 = cand2.m();
        auto sign2 = 1;
        if (cand2.pt() < 0) {
          sign2 = -1;
        }
        if (sign1 == sign2) {
          if (sign1 == 1) {
            registry.fill(HIST("hMassDplusPair"), mass2, mass1);
          } else {
            registry.fill(HIST("hMassDminusPair"), mass2, mass1);
          }
        } else {
          registry.fill(HIST("hMassDplusminusPair"), mass2, mass1);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusDplusReduced, processLocalData, "Process local data", true);

  void processLocalDataMcRec(o2::aod::HfCandDpFullEvs::iterator const&,
                             SelectedCandidates const& localCandidates)
  {
    registry.fill(HIST("hNMcRec"), localCandidates.size());

    for (const auto& cand1 : localCandidates) {
      std::vector<float> inputFeatures{cand1.ptProng0(), cand1.impactParameter0(), cand1.impactParameterZ0(),
                                       cand1.ptProng1(), cand1.impactParameter1(), cand1.impactParameterZ1(),
                                       cand1.ptProng2(), cand1.impactParameter2(), cand1.impactParameterZ2()};
      auto pt=std::abs(cand1.pt());
      //if (pt<1) continue;

      bool isSelMl = mlResponse.isSelectedMl(inputFeatures, pt, outputMl);
      if (outputMl[2] < 0.04) continue;  // Bkg score

      registry.fill(HIST("hPrompt"), outputMl[0]);
      registry.fill(HIST("hNonPrompt"), outputMl[1]);
      registry.fill(HIST("hBkg"), outputMl[2]);
      
      auto mass1 = cand1.m();
      if (cand1.pt() < 0) {
        registry.fill(HIST("hMassDminus"), mass1);
        if (std::abs(cand1.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)
          registry.fill(HIST("hMassDminusMatched"), mass1);
      } else {
        registry.fill(HIST("hMassDplus"), mass1);
        if (std::abs(cand1.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) {
          registry.fill(HIST("hMassDplusMatched"), mass1);
          registry.fill(HIST("hPromptMatched"), outputMl[0]);
          registry.fill(HIST("hNonPromptMatched"), outputMl[1]);
          registry.fill(HIST("hBkgMatched"), outputMl[2]);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusDplusReduced, processLocalDataMcRec, "Process local MC data", false);

  void processLocalDataMcGen(o2::aod::HfCandDpMcEvs::iterator const&,
                             SelectedMcParticles const& localMcParticles)
  {
    registry.fill(HIST("hNMcGen"), localMcParticles.size());

    for (const auto& part1 : localMcParticles) {
      for (auto part2 = part1 + 1; part2 != localMcParticles.end(); ++part2) {
        registry.fill(HIST("hDltPhiMcGen"), part2.phi() - part1.phi());
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusDplusReduced, processLocalDataMcGen, "Process local MC data at the gen level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskCorrelationDplusDplusReduced>(cfgc)};
}
