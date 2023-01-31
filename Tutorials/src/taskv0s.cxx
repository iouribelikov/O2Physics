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
// o2-analysis-timestamp --configuration json://myconfig.json | \
//   o2-analysis-track-propagation --configuration json://myconfig.json | \
   o2-analysistutorial-taskv0s --configuration json://myconfig.json -b --aod-file AO2D.root
*/

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/Core/trackUtilities.h"
#include "DetectorsVertexing/DCAFitterN.h"

using namespace o2;
using namespace o2::framework;

using myTracks = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra>;
using myTrack = myTracks::iterator;
using myV0s = aod::V0s;
using myV0 = myV0s::iterator;

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

struct taskv0s {
  Configurable<float> cfgBz{"Bz", 5., "Bz component of the solenoid magnetic field (kG)"};
  Configurable<float> cfgDcaDaugh = {"dcaDaugh", 0.05, "Max. allowed DCA between daughters"};
  Configurable<float> cfgRmin = {"rMin", 0.3, "Min. allowed radius for a V0 decay"};
  Configurable<float> cfgRmax = {"rMax", 100., "Max. allowed radius for a V0 decay"};

  float piMass = 0.1396;
  float pMass = 0.9383;
  float lambdaMass = 1.1157;
  float k0sMass = 0.4937;

  float mChi2;
  float mR;
  float mD;

  float mK0sMass;
  float mLambdaMass;
  float mLambdaBarMass;

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

  void v0Mass(float pxp, float pyp, float pzp, float pxn, float pyn, float pzn)
  {
    auto px = pxp + pxn;
    auto py = pyp + pyn;
    auto pz = pzp + pzn;
    auto p2 = px * px + py * py + pz * pz;
    auto p2p = pxp * pxp + pyp * pyp + pzp * pzp;
    auto p2n = pxn * pxn + pyn * pyn + pzn * pzn;

    auto ep = sqrt(piMass * piMass + p2p);
    auto en = sqrt(piMass * piMass + p2n);
    auto e = ep + en;
    mK0sMass = sqrt(e * e - p2);

    ep = sqrt(pMass * pMass + p2p);
    en = sqrt(piMass * piMass + p2n);
    e = ep + en;
    mLambdaMass = sqrt(e * e - p2);

    ep = sqrt(piMass * piMass + p2p);
    en = sqrt(pMass * pMass + p2n);
    e = ep + en;
    mLambdaBarMass = sqrt(e * e - p2);
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
    // if (track.itsNCls() != 7)
    //   return false;
    if (abs(track.eta()) > 0.9)
      return false;
    return true;
  }

  // V0 selector
  bool isV0Accepted(myV0 const& v0)
  {
    const auto& pos = v0.posTrack_as<myTracks>();
    const auto& neg = v0.negTrack_as<myTracks>();

    auto ptrack = getTrackParCov(pos);
    auto ntrack = getTrackParCov(neg);

    int nc = 0;
    try {
      nc = mFitter.process(ptrack, ntrack);
    } catch (...) {
      return false;
    }
    if (nc == 0)
      return false;

    int ibest = 0;
    mChi2 = 1e7;
    for (int i = 0; i < nc; i++) {
      auto chi2 = mFitter.getChi2AtPCACandidate(i);
      if (chi2 > mChi2)
        continue;
      mChi2 = chi2;
      ibest = i;
    }

    auto vtx = mFitter.getPCACandidate(ibest);
    auto x = vtx[0];
    auto y = vtx[1];
    auto z = vtx[2];
    auto r2 = x * x + y * y;
    if (r2 < cfgRmin * cfgRmin)
      return false;
    if (r2 > cfgRmax * cfgRmax)
      return false;

    const auto& t0 = mFitter.getTrack(0, ibest);
    const auto& t1 = mFitter.getTrack(1, ibest);

    auto r0 = t0.getXYZGlo();
    auto r1 = t1.getXYZGlo();
    x = r0.X() - r1.X();
    y = r0.Y() - r1.Y();
    z = r0.Z() - r1.Z();
    auto d2 = x * x + y * y + z * z;
    if (d2 > cfgDcaDaugh * cfgDcaDaugh)
      return false;

    mR = sqrt(r2);
    mD = sqrt(d2);

    std::array<float, 3> p0;
    t0.getPxPyPzGlo(p0);
    std::array<float, 3> p1;
    t1.getPxPyPzGlo(p1);

    v0Mass(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2]);

    return true;
  }

  // My V0 selector
  bool isV0AcceptedMy(myV0 const& v0)
  {
    auto const& p1 = v0.posTrack_as<myTracks>();
    if (!isTrackAccepted(p1))
      return false;
    auto const& p2 = v0.negTrack_as<myTracks>();
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
    y1 += py1 * t1; // z1 += pz1*t1;

    float r2 = x1 * x1 + y1 * y1;
    if (r2 > cfgRmax * cfgRmax)
      return false;
    if (r2 < cfgRmin * cfgRmin)
      return false;

    mChi2 = 0.;
    mR = sqrt(r2);
    mD = dca;

    v0Mass(px1, py1, pz1, px2, py2, pz2);

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

  o2::vertexing::DCAFitter2 mFitter;

  void init(o2::framework::InitContext& ic)
  {
    mFitter.setBz(cfgBz);
    mFitter.setPropagateToPCA(true);
    mFitter.setMaxR(30);
    mFitter.setMaxDZIni(0.1);
    mFitter.setMaxDXYIni(0.1);
    mFitter.setMinParamChange(1e-3);
    mFitter.setMinRelChi2Change(0.9);
    mFitter.setMaxChi2(10);
  }

  void process(aod::Collision const& collision, myTracks const& tracks, myV0s const& v0s)
  {
    static int ncol = 0;

    LOG(info) << "Collision: " << ncol++;

    if (!isCollisionAccepted(collision))
      return;

    // Basic collision counter...
    hVtxSel->Fill(collision.posZ());

    for (auto& track : tracks) {
      auto mom = track.tpcInnerParam();
      auto dex = track.tpcSignal();
      hdEdx->Fill(mom, dex);
    }

    // LOG(info) << "V0s: " << v0s.size() << " tracks: " << tracks.size();

    for (auto& v0 : v0s) {
      const auto& pos = v0.posTrack_as<myTracks>();
      const auto& neg = v0.negTrack_as<myTracks>();

      // No selections yet...
      v0Mass(pos.px(), pos.py(), pos.pz(), neg.px(), neg.py(), neg.pz());
      hK0sMass->Fill(mK0sMass);
      hLambdaMass->Fill(mLambdaMass);
      hLambdaMass->Fill(mLambdaBarMass);

      /*
      if (!isV0AcceptedMy(v0))  // Masses are re-calculated here
        continue;
      */
      if (!isV0Accepted(v0)) // Masses are re-calculated here
        continue;

      // PID
      if (isK0sLikeV0(pos, neg)) {
        hdEdxSel->Fill(pos.tpcInnerParam(), pos.tpcSignal());
        hdEdxSel->Fill(neg.tpcInnerParam(), neg.tpcSignal());
        hK0sMassSel->Fill(mK0sMass);
      }
      if (isLambdaLikeV0(pos, neg)) {
        hLambdaMassSel->Fill(mLambdaMass);
      }
      if (isLambdaLikeV0(neg, pos)) {
        hLambdaMassSel->Fill(mLambdaBarMass);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskv0s>(cfgc),
  };
}
