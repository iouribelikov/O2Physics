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

/*
Usage:
   o2-analysis-timestamp --configuration json://myconfig.json | \
   o2-analysis-track-propagation --configuration json://myconfig.json | \
   o2-analysis-taskqa --configuration json://myconfig.json -b --aod-file AO2D.root
*/

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;

struct taskqa {
  Configurable<float> Bz{"Bz", 5., "Bz component of the solenoid magnetic field (kG)"};
  Configurable<float> minContrib{"minContrib", 7., "Min. allowed number of contributor to the PV"};
  Configurable<float> maxContrib{"maxContrib", 7777., "Max. allowed number of contributor to the PV"};

  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> ncoHistogram{TH1F("ncoHistogram", "ncoHistogram", 100, 0., maxContrib)};
  OutputObj<TH1F> chiHistogram{TH1F("chiHistogram", "chiHistogram", 100, 0., +13.)};
  OutputObj<TH1F> vtxHistogram{TH1F("vtxHistogram", "vtxHistogram", 100, -30., +30.)};
  OutputObj<TH1F> ft0Histogram{TH1F("tf0Histogram", "tf0Histogram", 33, 0., +4000 * 30.)};
  OutputObj<TH1F> etaHistogram{TH1F("etaHistogram", "etaHistogram", 200, -3., +3)};
  // OutputObj<TH2F> dcaHistogram{TH2F("dcaHistogram", "dcaHistogram", 100, 0., 7., 100, -0.1, +0.1)};
  OutputObj<TH2F> trdHistogram{TH2F("trdHistogram", "trdHistogram;PV contributors;PV contrib. with TRD", 100, 0., 7777., 100, 0., 2500.)};
  OutputObj<TH2F> mulHistogram{TH2F("mulHistogram", "mulHistogram;TF0 signal;PV contrib. with TRD", 100, 0., 4000 * 30, 100, 0., 2500)};

  void init(o2::framework::InitContext& initContext)
  {
    // LOG(info) << "MyMyMy "<<maxContrib<<' '<<minContrib<<'\n';
    ncoHistogram->SetBins(100, 0., maxContrib);
  }

  using myTracks = soa::Join<aod::FullTracks, aod::TracksDCA>;
  using myTrack = myTracks::iterator;

  bool isSelected(aod::Collision const& coll)
  {
    auto zv = coll.posZ();
    if (abs(zv) > 10)
      return false;

    auto nc = coll.numContrib();
    if (nc < minContrib)
      return false;

    auto chi2 = coll.chi2() / nc;
    if (chi2 > 7.)
      return false;

    chiHistogram->Fill(chi2);
    vtxHistogram->Fill(zv);
    ncoHistogram->Fill(nc);

    return true;
  }

  bool isSelected(myTrack const& track)
  {
    auto eta = track.eta();
    if (abs(eta) > 1)
      return false;

    auto ncl = track.itsNCls();
    if (ncl != 7)
      return false;

    if (!track.has_collision())
      return false;
    // auto const& coll = track.collision();

    etaHistogram->Fill(eta);

    return true;
  }

  void process(aod::BC const& bc, aod::FT0s const& ft0s, aod::Collisions const& collisions, myTracks const& tracks)
  {
    auto ncolls = collisions.size();
    if (ncolls == 0) {
      return;
    }
    if (ncolls > 1) {
      LOG(info) << "!!! In-bunch pileup: " << ncolls;
      return;
    }

    auto nft0s = ft0s.size();
    if (nft0s == 0) {
      return;
    }
    if (nft0s > 1) {
      LOG(error) << "!!! Strange in-bunch pileup in FT0 : " << nft0s;
    }

    auto coll = collisions.begin();
    if (!isSelected(coll))
      return;

    static int ncol = 0;
    if (ncol % 1000 == 0) {
      LOG(info) << "Collision: " << ncol << ' ' << coll.collisionTime() << ' ' << tracks.size();
    }
    ncol++;

    auto collId = coll.globalIndex();

    auto ft0 = ft0s.begin();
    auto ampA = ft0.amplitudeA();
    auto sigA = std::reduce(ampA.begin(), ampA.end());
    auto ampC = ft0.amplitudeC();
    auto sigC = std::reduce(ampC.begin(), ampC.end());
    auto sig = 0.5 * (sigA + sigC);

    ft0Histogram->Fill(sig);

    int nsel = 0;
    int ntrd = 0;
    for (auto& track : tracks) {
      if (!isSelected(track))
        continue;
      if (track.collisionId() != collId)
        continue;
      // dcaHistogram->Fill(track.pt(), track.dcaXY());
      nsel++;
      if (track.isPVContributor() && track.hasTRD())
        ntrd++;
    }
    trdHistogram->Fill(coll.numContrib(), ntrd);
    // mulHistogram->Fill(sig, nsel);
    mulHistogram->Fill(sig, ntrd);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskqa>(cfgc)};
}
