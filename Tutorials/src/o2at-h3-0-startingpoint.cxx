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
///
/// \brief this is a starting point for the third session of the tutorial
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point
struct startingpoint {

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {{"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
     {"hK0sMass", "hK0sMass", {HistType::kTH1F, {{nBins, -15., 15.}}}}}};

  // void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const & V0s)
  void process(aod::V0Datas const& V0s)
  {
    /*
    if (!collision.sel8()) {
      return;
    }
    */
    // registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());
    for (auto& v0 : V0s) {
      registry.fill(HIST("hK0sMass"), v0.mK0Short());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<startingpoint>(cfgc)};
}
