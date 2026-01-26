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
/// \brief Example how to lopps over V0s ...
/// \author
/// \since

/* Usage :
   o2-analysis-tracks-extra-v002-converter --configuration json://myconfig.json | \  
   o2-analysistutorial-task-yann --configuration json://myconfig.json -b --aod-file AO2D.root
*/

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/Core/trackUtilities.h"

using namespace o2;
using namespace o2::framework;

using myTracks = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra>;
using myTrack = myTracks::iterator;

struct taskYann {

  Configurable<float> cfgZmax{"zMax", 10., "Restriction on the PV position |PVz|<zMax (cm)"};

  OutputObj<TH1F> hVtx{
    TH1F("hVtx", "Primary vertex position after selection; Z (cm)", 100, -20., 20.)};

  OutputObj<TH2F> hdEdx{
    TH2F("hdEdx", "TPC; Momentum (GeV); dE/dx", 400, -4., 4., 200, 0., 1000)};

  void init(o2::framework::InitContext& ic)
  {
  }

  // Collision selector
  bool isCollisionAccepted(aod::Collision const& collision)
  {
    auto z = collision.posZ();
    if (abs(z) > cfgZmax)
      return false;

    // Some other selections
    
    return true;
  }

  // Track-quality selector
  template <typename TrackInstance>
  bool isTrackAccepted(TrackInstance const& track)
  {
    if (abs(track.tgl()) > 0.9)
      return false;
    
    if (track.itsNCls() < 7)
      return false;
    if (!track.hasTPC())
      return false;

    // Some other selections
    
    return true;
  }

  void processData(aod::Collision const& collision, myTracks const& tracks)
  {
    static int ncol = 0;

    if (ncol % 1000 == 0)
      LOG(info) << "Collision: " << ncol;
    ncol++;

    if (!isCollisionAccepted(collision))
      return;

    // Collision counter...
    hVtx->Fill(collision.posZ());

    for (auto& track : tracks) {
      if (!isTrackAccepted(track))
	continue;
      auto sign = track.sign();
      auto mom = track.tpcInnerParam();
      auto dex = track.tpcSignal();
      hdEdx->Fill(sign*mom, dex);
    }
  }
  PROCESS_SWITCH(taskYann, processData, "Process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskYann>(cfgc),
  };
}
