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
   o2-analysis-collision-converter --configuration json://myconfig.json | \
   o2-analysis-track-propagation --configuration json://myconfig.json | \
   o2-analysistutorial-taskdca --configuration json://myconfig.json -b --aod-file AO2D.root
*/

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;

using myTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksDCA, aod::TracksExtra>;
using myTrack = myTracks::iterator;

// Linear Vertex
#include <array>
#include <TMath.h>
#include <TMatrixD.h>

class LinearVertex
{
 public:
  LinearVertex() = default;
  LinearVertex(const LinearVertex& t) = default;
  LinearVertex& operator=(const LinearVertex& tr) = default;
  ~LinearVertex() = default;

  auto getX() const { return mX; }
  auto getY() const { return mY; }
  auto getZ() const { return mZ; }
  const auto& getCovariance() const { return mCov; }
  auto getCovX2() const { return mCov[0]; }
  auto getCovXY() const { return mCov[1]; }
  auto getCovY2() const { return mCov[2]; }
  auto getCovXZ() const { return mCov[3]; }
  auto getCovYZ() const { return mCov[4]; }
  auto getCovZ2() const { return mCov[5]; }
  auto getNumberOfProngs() const { return mProngs; }
  auto getChi2() const { return mChi2; }

  Bool_t update(myTrack const&);
  Bool_t update(const std::array<Double_t, 3>& p, const std::array<Double_t, 3>& v, Double_t sy2, Double_t sz2);

 private:
  Double_t mX = 0;                                       ///< vertex position
  Double_t mY = 0;                                       ///< vertex position
  Double_t mZ = 0;                                       ///< vertex position
  std::array<Double_t, 6> mCov{0.04, 0, 0.04, 0, 0, 36}; ///< vertex covariance matrix
  Int_t mProngs = 0;                                     ///< number of prongs
  Double_t mChi2 = 0;                                    ///< chi2
  static constexpr Double_t mMaxChi2 = 4;                ///< maximal accepted chi2 increment
};

Bool_t LinearVertex::update(const std::array<Double_t, 3>& p, const std::array<Double_t, 3>& v, Double_t sy2, Double_t sz2)
{
  //--------------------------------------------------------------------
  // Linear update of the vertex parameters
  //--------------------------------------------------------------------

  // Get the vertex weight matrix
  TMatrixD wv(3, 3);
  wv(0, 0) = mCov[0];
  wv(0, 1) = mCov[1];
  wv(0, 2) = mCov[3];
  wv(1, 0) = mCov[1];
  wv(1, 1) = mCov[2];
  wv(1, 2) = mCov[4];
  wv(2, 0) = mCov[3];
  wv(2, 1) = mCov[4];
  wv(2, 2) = mCov[5];
  wv.Invert();
  if (!wv.IsValid())
    return kFALSE;

  // Get the tracklet weight matrix
  TMatrixD wt(3, 3);
  wt(0, 0) = 0.;
  wt(0, 1) = 0.;
  wt(0, 2) = 0.;
  wt(1, 0) = 0.;
  wt(1, 1) = 1 / sy2;
  wt(1, 2) = 0.;
  wt(2, 0) = 0.;
  wt(2, 1) = 0.;
  wt(2, 2) = 1 / sz2; // FIXME : In "parallel" system

  auto l = v[0];
  auto m = v[1];
  auto n = v[2];
  auto phi = TMath::ATan2(m, l);
  auto sp = TMath::Sin(phi);
  auto cp = TMath::Cos(phi);
  auto tgl = n / TMath::Sqrt(l * l + m * m);
  auto cl = 1 / TMath::Sqrt(1. + tgl * tgl);
  auto sl = tgl * cl;

  TMatrixD p2g(3, 3); //"parallel" --> global transformation
  p2g(0, 0) = cp * cl;
  p2g(1, 0) = sp * cl;
  p2g(2, 0) = sl;
  p2g(0, 1) = -sp;
  p2g(1, 1) = cp;
  p2g(2, 1) = 0.;
  p2g(0, 2) = -sl * cp;
  p2g(1, 2) = -sl * sp;
  p2g(2, 2) = cl;
  wt = p2g * wt * TMatrixD(TMatrixD::kTransposed, p2g); // Now, in global system

  // Check the possible chi2 increment
  TMatrixD cv(wv);
  cv += wt;
  cv.Invert();
  if (!cv.IsValid())
    return kFALSE;

  wv = wv * cv * wt;

  auto lmn = TMath::Sqrt(l * l + m * m + n * n);
  auto cosx = l / lmn, cosy = m / lmn, cosz = n / lmn;
  auto pp = (l * (p[0] - mX) + m * (p[1] - mY) + n * (p[2] - mZ)) / lmn;
  auto tx = p[0] - pp * cosx;
  auto ty = p[1] - pp * cosy;
  auto tz = p[2] - pp * cosz;
  Double_t delta[3]{tx - mX, ty - mY, tz - mZ};
  Double_t chi2 = 0;
  for (Int_t i = 0; i < 3; i++) {
    Double_t s = 0.;
    for (Int_t j = 0; j < 3; j++)
      s += wv(i, j) * delta[j];
    chi2 += s * delta[i];
  }
  if (chi2 > mMaxChi2)
    return kFALSE;

  // Update the vertex covariance
  mCov[0] = cv(0, 0);
  mCov[1] = cv(1, 0);
  mCov[2] = cv(1, 1);
  mCov[3] = cv(2, 0);
  mCov[4] = cv(2, 1);
  mCov[5] = cv(2, 2);

  // Update the vertex position, chi2, and the number of prongs
  TMatrixD cwt(cv, TMatrixD::kMult, wt);
  for (Int_t j = 0; j < 3; j++)
    mX += cwt(0, j) * delta[j];
  for (Int_t j = 0; j < 3; j++)
    mY += cwt(1, j) * delta[j];
  for (Int_t j = 0; j < 3; j++)
    mZ += cwt(2, j) * delta[j];
  mChi2 += chi2;
  mProngs++;

  return kTRUE;
}

Bool_t LinearVertex::update(myTrack const& track)
{
  auto alpha = track.alpha();
  auto xt = track.x();
  auto yt = track.y();
  auto z = track.z();
  auto cs = cos(alpha), sn = sin(alpha);
  auto x = cs * xt - sn * yt, y = sn * xt + cs * yt;
  std::array<Double_t, 3> p{x, y, z};
  std::array<Double_t, 3> v{track.px(), track.py(), track.pz()};

  return update(p, v, track.cYY(), track.cZZ());
}

// The task...
struct taskdca {
  Configurable<float> Bz{"Bz", 5., "Bz component of the solenoid magnetic field (kG)"};
  Configurable<float> minContrib{"minContrib", 7., "Min. allowed number of contributor to the PV"};
  Configurable<float> maxContrib{"maxContrib", 33., "Max. allowed number of contributor to the PV"};

  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> ncoHistogram{TH1F("ncoHistogram", "ncoHistogram", 100, 0., maxContrib)};
  OutputObj<TH1F> chiHistogram{TH1F("chiHistogram", "chiHistogram", 100, 0., +13.)};
  OutputObj<TH1F> vtxHistogram{TH1F("vtxHistogram", "vtxHistogram", 100, -30., +30.)};

  OutputObj<TH1F> clsHistogram{TH1F("clsHistogram", "clsHistogram", 10, 0., +10.)};
  OutputObj<TH1F> etaHistogram{TH1F("etaHistogram", "etaHistogram", 200, -3., +3)};

  OutputObj<TH2F> dcaHistogram{TH2F("dcaHistogram", "dcaHistogram", 100, 0., 7., 100, -0.1, +0.1)};
  OutputObj<TH2F> ip0Histogram{TH2F("ip0Histogram", "ip0Histogram", 100, 0., 7., 100, -0.1, +0.1)};

  void init(o2::framework::InitContext& initContext)
  {
    // LOG(info) << "MyMyMy "<<maxContrib<<' '<<minContrib<<'\n';
    ncoHistogram->SetBins(100, 0., maxContrib);

    Int_t nb = 100;
    Double_t xbins[nb + 1], ptcutl = 0.1, ptcuth = 10.;
    Double_t ybins[nb + 1], ycutl = -0.1, ycuth = +0.1;
    Double_t a = TMath::Log(ptcuth / ptcutl) / nb;
    Double_t d = (ycuth - ycutl) / nb;
    for (Int_t i = 0; i <= nb; i++) {
      xbins[i] = ptcutl * TMath::Exp(i * a);
      ybins[i] = ycutl + d * i;
    }
    dcaHistogram->SetBins(nb, xbins, nb, ybins);
    ip0Histogram->SetBins(nb, xbins, nb, ybins);
  }

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

  bool isSelected(myTrack const& track, bool fill = true)
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

    if (fill) {
      etaHistogram->Fill(eta);
      clsHistogram->Fill(ncl);
    }
    return true;
  }

  void getImpactParams(myTrack const& track, float x, float y, float z, float bz, float ip[2])
  {
    //------------------------------------------------------------------
    // This function calculates the transverse and longitudinal impact parameters
    // with respect to a point with global coordinates (x,y,0)
    // in the magnetic field "bz" (kG)
    //------------------------------------------------------------------
    auto f1 = track.snp();
    auto r1 = sqrt((1. - f1) * (1. + f1));
    auto xt = track.x(), yt = track.y();
    auto sn = sin(track.alpha()), cs = cos(track.alpha());
    float a = x * cs + y * sn;
    y = -x * sn + y * cs;
    x = a;
    xt -= x;
    yt -= y;

    float rp4 = track.signed1Pt() * bz * o2::constants::math::B2C;
    if ((abs(bz) < 1e-33) || (abs(rp4) < 1e-33)) {
      ip[0] = -(xt * f1 - yt * r1);
      ip[1] = track.z() + (ip[0] * f1 - xt) / r1 * track.tgl() - z;
      return;
    }

    sn = rp4 * xt - f1;
    cs = rp4 * yt + r1;
    a = 2 * (xt * f1 - yt * r1) - rp4 * (xt * xt + yt * yt);
    auto rr = sqrt(sn * sn + cs * cs);
    ip[0] = -a / (1 + rr);
    auto f2 = -sn / rr;
    auto r2 = sqrt((1. - f2) * (1. + f2));
    ip[1] = track.z() + track.tgl() / rp4 * asin(f2 * r1 - f1 * r2) - z;
  }

  void process(aod::Collision const& coll, myTracks const& tracks)
  {
    if (!isSelected(coll))
      return;

    LinearVertex vertexer;
    for (auto& track : tracks) {
      if (!isSelected(track))
        continue;
      if (track.pt() < 0.3)
        continue;

      vertexer.update(track);
    }

    if (vertexer.getNumberOfProngs() < minContrib)
      return;
    auto vx = vertexer.getX(), vy = vertexer.getY(), vz = vertexer.getZ();

    for (auto& track : tracks) {
      if (!isSelected(track, false))
        continue;

      float ip[2];
      getImpactParams(track, vx, vy, vz, Bz, ip);
      ip0Histogram->Fill(track.pt(), ip[0]);

      dcaHistogram->Fill(track.pt(), track.dcaXY());

      // LOG(info) << ip[0] << ' ' << track.dcaXY() << " ip_dca";
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskdca>(cfgc)};
}
