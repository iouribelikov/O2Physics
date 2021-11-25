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
  float k0sMass = 0.4937;

  OutputObj<TH1F> hK0sMass{
    TH1F("hK0sMass", "K0s mass; Mass (GeV)", 100, k0sMass - 0.2, k0sMass + 0.2)};
  OutputObj<TH1F> hK0sMassSel{
    TH1F("hK0sMassSel", "K0s mass after selection; Mass (GeV)", 100, k0sMass - 0.2, k0sMass + 0.2)};

  void process(aod::Collision const& collision, aod::V0s const& v0s, aod::Tracks const& tracks)
  {
    for (auto& v0 : v0s) {
      LOGF(debug, "V0 (%d, %d, %d)", v0.posTrack().collisionId(), v0.negTrack().collisionId(), v0.collisionId());

      auto p1 = v0.posTrack();
      auto e1 = sqrt(piMass * piMass + p1.p() * p1.p());
      auto x1 = p1.x() * cos(p1.alpha()) - p1.y() * sin(p1.alpha());
      auto y1 = p1.x() * sin(p1.alpha()) + p1.y() * cos(p1.alpha());
      auto z1 = p1.z();
      auto px1 = p1.px();
      auto py1 = p1.py();
      auto pz1 = p1.pz();

      auto p2 = v0.negTrack();
      auto e2 = sqrt(piMass * piMass + p2.p() * p2.p());
      auto x2 = p2.x() * cos(p2.alpha()) - p2.y() * sin(p2.alpha());
      auto y2 = p2.x() * sin(p2.alpha()) + p2.y() * cos(p2.alpha());
      auto z2 = p2.z();
      auto px2 = p2.px();
      auto py2 = p2.py();
      auto pz2 = p2.pz();

      auto mass = (e1 + e2) * (e1 + e2) - (px1 + px2) * (px1 + px2) - (py1 + py2) * (py1 + py2) - (pz1 + pz2) * (pz1 + pz2);
      mass = sqrt(mass);

      hK0sMass->Fill(mass);

      // DCA between daughter tracks
      auto dd = Det(x2 - x1, y2 - y1, z2 - z1, px1, py1, pz1, px2, py2, pz2);
      auto ax = Det(py1, pz1, py2, pz2);
      auto ay = -Det(px1, pz1, px2, pz2);
      auto az = Det(px1, py1, px2, py2);

      auto dca = TMath::Abs(dd) / TMath::Sqrt(ax * ax + ay * ay + az * az);
      if (dca > cfgDcaDaugh)
        continue;

      // V0 vertex
      double t1 = Det(x2 - x1, y2 - y1, z2 - z1, px2, py2, pz2, ax, ay, az) /
                  Det(px1, py1, pz1, px2, py2, pz2, ax, ay, az);

      x1 += px1 * t1;
      y1 += py1 * t1; //z1 += pz1*t1;

      if (x1 * x1 + y1 * y1 > cfgRmax * cfgRmax)
        continue;
      if (x1 * x1 + y1 * y1 < cfgRmin * cfgRmin)
        continue;

      hK0sMassSel->Fill(mass);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LoopV0s>(cfgc),
  };
}
