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
/// \brief Example how to enumerate V0s and cascades. Note ...
///        V0s = Join<TransientV0s, StoredV0s>
///        Cascades = Join<TransientCascades, StoredCascades>
///        TransientV0 and TransientCascades are filled by the helper task weak-decay-indices. Hence use ...
///        o2-analysis-weak-decay-indices --aod-file AO2D.root | o2-analysistutorial-weak-decay-iteration
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

double Det(double a00, double a01, double a10, double a11)
{
  //--------------------------------------------------------------------
  // This function calculates locally a 2x2 determinant
  //--------------------------------------------------------------------
  return a00 * a11 - a01 * a10;
}

double Det(double a00, double a01, double a02,
           double a10, double a11, double a12,
           double a20, double a21, double a22)
{
  //--------------------------------------------------------------------
  // This function calculates locally a 3x3 determinant
  //--------------------------------------------------------------------
  return a00 * Det(a11, a12, a21, a22) - a01 * Det(a10, a12, a20, a22) + a02 * Det(a10, a11, a20, a21);
}

struct LoopV0s {
  Configurable<float> cfgDcaDaugh = {"dcaDaugh", 0.05, "Max. allowed DCA between daughters"};
  Configurable<float> cfgRmin = {"rMin", 0.3, "Min. allowed radius for a V0 decay"};
  Configurable<float> cfgRmax = {"rMax", 100., "Max. allowed radius for a V0 decay"};

  float piMass = 0.1396;
  float pMass = 0.9383;
  float lambdaMass = 1.1157;
  float k0sMass = 0.4937;

  OutputObj<TH1F> hVtxSel{
    TH1F("hVtx", "Primary vertex position after selection; Z (cm)", 100, -20., 20.)};

  OutputObj<TH2F> hdEdx{
    TH2F("hdEdx", "TPC; Momentum (GeV); dE/dx", 200, 0., 4., 200, 0., 1000)};
  OutputObj<TH2F> hdEdxSel{
    TH2F("hdEdxSel", "TPC after selection; Momentum (GeV); dE/dx", 200, 0., 4., 200, 0., 1000)};

  OutputObj<TH1F> hK0sMass{
    TH1F("hK0sMass", "K0s mass; Mass (GeV)", 100, k0sMass - 0.2, k0sMass + 0.2)};
  OutputObj<TH1F> hK0sMassSel{
    TH1F("hK0sMassSel", "K0s mass after selection; Mass (GeV)", 100, k0sMass - 0.2, k0sMass + 0.2)};

  OutputObj<TH1F> hLambdaMass{
    TH1F("hLambdaMass", "Lambda+LambdaBar mass; Mass (GeV)", 120, lambdaMass - 0.2, lambdaMass + 0.2)};
  OutputObj<TH1F> hLambdaMassSel{
    TH1F("hLambdaMassSel", "Lambda+LambdaBar mass after selection; Mass (GeV)", 120, lambdaMass - 0.2, lambdaMass + 0.2)};

  template <typename V0Instance>
  float v0Mass(V0Instance const& v0, float pMass, float nMass)
  {
    auto const& p1 = v0.posTrack();
    auto px1 = p1.px();
    auto py1 = p1.py();
    auto pz1 = p1.pz();

    auto const& p2 = v0.negTrack();
    auto px2 = p2.px();
    auto py2 = p2.py();
    auto pz2 = p2.pz();

    auto e1 = sqrt(pMass * pMass + p1.p() * p1.p());
    auto e2 = sqrt(nMass * nMass + p2.p() * p2.p());
    return sqrt((e1 + e2) * (e1 + e2) - (px1 + px2) * (px1 + px2) - (py1 + py2) * (py1 + py2) - (pz1 + pz2) * (pz1 + pz2));
  }

  // Collision selector
  bool isCollisionAccepted(aod::Collision const& collision)
  {
    auto z = collision.posZ();
    if (abs(z) > 10)
      return false;
    return true;
  }

  // Track-quality selector
  template <typename TrackInstance>
  bool isTrackAccepted(TrackInstance const& track)
  {
    if (abs(track.eta()) > 0.9)
      return false;
    return true;
  }

  // V0 selector
  template <typename V0Instance>
  bool isV0Accepted(V0Instance const& v0)
  {
    auto const& p1 = v0.posTrack();
    if (!isTrackAccepted(p1))
      return false;
    auto const& p2 = v0.negTrack();
    if (!isTrackAccepted(p2))
      return false;

    // Topological selections...
    auto px1 = p1.px();
    auto py1 = p1.py();
    auto pz1 = p1.pz();
    auto x1 = p1.x() * cos(p1.alpha()) - p1.y() * sin(p1.alpha());
    auto y1 = p1.x() * sin(p1.alpha()) + p1.y() * cos(p1.alpha());
    auto z1 = p1.z();

    auto px2 = p2.px();
    auto py2 = p2.py();
    auto pz2 = p2.pz();
    auto x2 = p2.x() * cos(p2.alpha()) - p2.y() * sin(p2.alpha());
    auto y2 = p2.x() * sin(p2.alpha()) + p2.y() * cos(p2.alpha());
    auto z2 = p2.z();

    // DCA between daughter tracks
    auto dd = Det(x2 - x1, y2 - y1, z2 - z1, px1, py1, pz1, px2, py2, pz2);
    auto ax = Det(py1, pz1, py2, pz2);
    auto ay = -Det(px1, pz1, px2, pz2);
    auto az = Det(px1, py1, px2, py2);

    auto dca = TMath::Abs(dd) / TMath::Sqrt(ax * ax + ay * ay + az * az);
    if (dca > cfgDcaDaugh)
      return false;

    // V0 vertex position (fiducial volume)
    double t1 = Det(x2 - x1, y2 - y1, z2 - z1, px2, py2, pz2, ax, ay, az) /
                Det(px1, py1, pz1, px2, py2, pz2, ax, ay, az);

    x1 += px1 * t1;
    y1 += py1 * t1; //z1 += pz1*t1;

    if (x1 * x1 + y1 * y1 > cfgRmax * cfgRmax)
      return false;
    if (x1 * x1 + y1 * y1 < cfgRmin * cfgRmin)
      return false;

    return true;
  }

  // Basic PID selectors
  template <typename TrackInstance>
  bool isPion(TrackInstance const& track)
  {
    auto mom = track.tpcInnerParam();
    auto dedx = track.tpcSignal();
    if (mom > 0.1 && dedx > 250.)
      return false;
    if (mom > 0.4 && dedx > 65.)
      return false;
    return true;
  }
  template <typename TrackInstance>
  bool isProton(TrackInstance const& track)
  {
    return !isPion(track);
  }
  template <typename TrackInstance>
  bool isK0sLikeV0(TrackInstance const& pos, TrackInstance const& neg)
  {
    if (!isPion(pos))
      return false;
    if (!isPion(neg))
      return false;
    return true;
  }
  template <typename TrackInstance>
  bool isLambdaLikeV0(TrackInstance const& pos, TrackInstance const& neg)
  {
    if (pos.tpcInnerParam() < 1.2)
      if (!isProton(pos))
        return false;
    if (!isPion(neg))
      return false;
    return true;
  }

  void process(aod::Collision const& collision, aod::V0s const& v0s, aod::Tracks const& tracks, aod::TracksExtra const& exts)
  {

    if (!isCollisionAccepted(collision))
      return;

    // Basic collision counter...
    hVtxSel->Fill(collision.posZ());

    for (auto& track : exts) {
      auto mom = track.tpcInnerParam();
      auto dex = track.tpcSignal();
      hdEdx->Fill(mom, dex);
    }

    for (auto& v0 : v0s) {
      // No selections yet...
      auto massK0s = v0Mass(v0, piMass, piMass);
      hK0sMass->Fill(massK0s);
      auto massLambda = v0Mass(v0, pMass, piMass);
      hLambdaMass->Fill(massLambda);
      auto massLambdaBar = v0Mass(v0, piMass, pMass);
      hLambdaMass->Fill(massLambdaBar);

      if (!isV0Accepted(v0))
        continue;

      auto pos = *(exts.begin() + v0.posTrackId());
      auto neg = *(exts.begin() + v0.negTrackId());
      if (isK0sLikeV0(pos, neg)) {
        hdEdxSel->Fill(pos.tpcInnerParam(), pos.tpcSignal());
        hdEdxSel->Fill(neg.tpcInnerParam(), neg.tpcSignal());
        hK0sMassSel->Fill(massK0s);
      }
      if (isLambdaLikeV0(pos, neg)) {
        hLambdaMassSel->Fill(massLambda);
      }
      if (isLambdaLikeV0(neg, pos)) {
        hLambdaMassSel->Fill(massLambdaBar);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LoopV0s>(cfgc),
  };
}
