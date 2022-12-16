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
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;

//STEP 0
//This is an empty analysis skeleton: the starting point! 
struct taskdca {
  Configurable<float> Bz{"Bz", 5., "Bz component of the solenoid magnetic field (kG)"};
  Configurable<float> minContrib{"minContrib", 7., "Min. allowed number of contributor to the PV"};
  Configurable<float> maxContrib{"maxContrib", 33., "Max. allowed number of contributor to the PV"};

  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> ncoHistogram{TH1F("ncoHistogram", "ncoHistogram", 100, 0., maxContrib)};
  OutputObj<TH1F> vtxHistogram{TH1F("vtxHistogram", "vtxHistogram", 100, -30., +30.)};
  OutputObj<TH1F> clsHistogram{TH1F("clsHistogram", "clsHistogram", 10, 0., +10.)};
  OutputObj<TH1F> etaHistogram{TH1F("etaHistogram", "etaHistogram", 200, -1., +1)};
  OutputObj<TH2F> dcaHistogram{TH2F("dcaHistogram", "dcaHistogram", 100, 0., 7., 100, -0.1, +0.1)};

  void init(o2::framework::InitContext& initContext)
  {
    // LOG(info) << "MyMyMy "<<maxContrib<<' '<<minContrib<<'\n';
    ncoHistogram->SetBins(100, 0., maxContrib);
  }

  void process(aod::Collisions const& collisions,
	       soa::Join<aod::FullTracks, aod::TracksDCA> const& tracks)
  {
    for (auto& coll : collisions) {
      auto zv=coll.posZ();
      auto nc=coll.numContrib();
      vtxHistogram->Fill(zv);
      ncoHistogram->Fill(nc);
    }
    for (auto& track : tracks) {
      etaHistogram->Fill(track.eta());
      clsHistogram->Fill(track.itsNCls());
      if (abs(track.eta())>1) continue; 
      if (track.itsNCls()!=7) continue;
      if (!track.has_collision()) continue;
      auto const& coll = track.collision();
      auto zv=coll.posZ();
      auto nc=coll.numContrib();
      if (abs(zv)>10) continue;
      if (abs(nc)<minContrib) continue;
      dcaHistogram->Fill(track.pt(), track.dcaXY());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskdca>(cfgc)};
}
