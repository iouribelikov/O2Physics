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

using myTracks = soa::Join<aod::TracksIU, aod::TracksExtra>;
using myTrack = myTracks::iterator;

struct taskYann {

  Configurable<float> cfgZmax{"zMax", 10., "Restriction on the PV position |PVz|<zMax (cm)"};

  OutputObj<TH1F> hVtx{
    TH1F("hVtx", "Primary vertex position after selection; Z (cm)", 100, -20., 20.)};

  /*
  OutputObj<TH1F> hMass{
    TH1F("hMass", "Invariant mass; Mpp (GeV)", 4000, 1.8, 3.8)};
  */
  OutputObj<TH1F> hMass{
    TH1F("hMass", "Invariant mass; Mppi (GeV)", 50, 1.06, 1.16)};
  OutputObj<TH1F> hMassMatch{
    TH1F("hMassMatch", "Invariant mass; Mppi (GeV)", 50, 1.06, 1.16)};

  OutputObj<TH1F> hPdgCode{
    TH1F("hPdgCode", "PDG code; code", 2*3200, -3200, 3200)};

  OutputObj<TH2F> hdEdx{
    TH2F("hdEdx", "TPC; Momentum (GeV); dE/dx", 400, -4., 4., 200, 0., 1000)};
  OutputObj<TH2F> hdEdxPr{
    TH2F("hdEdxPr", "TPC; Momentum (GeV); dE/dx", 400, -4., 4., 200, 0., 1000)};

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

  template <typename TrackInstance>
  bool isProton(TrackInstance const& track)
  {
    auto p = track.tpcInnerParam();
    auto dedx = track.tpcSignal();
    if (p < 0.2) return false;
    if (dedx < 40) return false;
    if (dedx < 40) return false;
    if (dedx < 600 - 500/0.58*p) return false;
    if (dedx < 250 - 200/1.1*p) return false;
    if (dedx < 120 - 120/3.0*p) return false;
    return true;
  }
  
  template <typename TrackInstance>
  float invariantMass(TrackInstance const& neg, TrackInstance const& pos)
  {
    const float pMass=0.938;
    const float piMass=0.138;
    
    auto pxn = neg.px();
    auto pyn = neg.py();
    auto pzn = neg.pz();

    auto pxp = pos.px();
    auto pyp = pos.py();
    auto pzp = pos.pz();

    auto px = pxp + pxn;
    auto py = pyp + pyn;
    auto pz = pzp + pzn;
    auto p2 = px * px + py * py + pz * pz;

    auto p2p = pxp * pxp + pyp * pyp + pzp * pzp;
    auto p2n = pxn * pxn + pyn * pyn + pzn * pzn;

    auto ep = sqrt(piMass * piMass + p2p);
    auto en = sqrt(pMass * pMass + p2n);
    auto e = ep + en;
    auto mass = sqrt(e * e - p2);

    return mass;
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
      auto dedx = track.tpcSignal();
      hdEdx->Fill(sign*mom, dedx);
      if (!isProton(track))
	continue;
      
      hdEdxPr->Fill(sign*mom, dedx);

      if (sign > 0) continue;

      for (auto track1 = track + 1; track1 != tracks.end(); ++track1) {
        if (!isTrackAccepted(track1))
	  continue;
        auto sign1 = track1.sign();
        if (sign1 < 0) continue;

	auto mass = invariantMass(track, track1);
	hMass->Fill(mass);
      }
    }
  }
  PROCESS_SWITCH(taskYann, processData, "Process data", true);

  void processMcGen(aod::McParticles& particles) {
    for (auto &p : particles) {
      if (abs(p.eta())>0.9) continue;
      if (abs(p.pt()) <1.0) continue;
      auto x=p.vx();
      auto y=p.vy();
      if (x*x+y*y > 2*2) continue;
      auto code=p.pdgCode();
      if (code == -3122)
        LOG(info) << "LambdaBar: " << p.p();
      hPdgCode->Fill(code);
    }
  }
  PROCESS_SWITCH(taskYann, processMcGen, "Process MC at the generator level", true);
  
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskYann>(cfgc),
  };
}
