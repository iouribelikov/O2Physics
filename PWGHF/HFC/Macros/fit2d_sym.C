#include "TF2.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"

Double_t sig_sig(Double_t* x, Double_t* par)
{
  auto r1 = (x[0] - par[1]);
  auto r2 = (x[1] - par[1]);
  return par[0] * TMath::Exp(-0.5 * (r1 * r1 / par[2] + r2 * r2 / par[2]));
}

Double_t sig_bkg(Double_t* x, Double_t* par)
{
  auto r1 = (x[0] - par[1]);
  auto r2 = (x[1] - par[1]);
  return par[3] * TMath::Exp(-0.5 * r1 * r1 / par[2]) * (1 + r2 * par[4] + r2 * r2 * par[6]);
}

Double_t bkg_sig(Double_t* x, Double_t* par)
{
  auto r1 = (x[0] - par[1]);
  auto r2 = (x[1] - par[1]);
  return par[3] * TMath::Exp(-0.5 * r2 * r2 / par[2]) * (1 + r1 * par[4] + r1 * r1 * par[6]);
}

Double_t bkg_bkg(Double_t* x, Double_t* par)
{
  auto r1 = (x[0] - par[1]);
  auto r2 = (x[1] - par[1]);
  return par[5] * (1 + r1 * par[4] + r1 * r1 * par[6]) * (1 + r2 * par[4] + r2 * r2 * par[6]);
}

Double_t fun2(Double_t* x, Double_t* par)
{
  return sig_sig(x, par) + sig_bkg(x, par) + bkg_sig(x, par) + bkg_bkg(x, par);
}


void fit2d_sym()
{
  Double_t sig_sig = 135./2, mass = 1.871, sigma2 = 0.009 * 0.009;
  Double_t sig_bkg = 600./2, slope = -1.0;
  Double_t bkg_bkg = 750./2;
  Double_t f2params[] = {
    sig_sig, mass, sigma2,
    sig_bkg, slope,
    bkg_bkg,
    0. // Quadratic coefficients
  };
  const Int_t npar = sizeof(f2params) / sizeof(Double_t);

  Float_t mmin = 1.80, mmax = 1.94;
  TF2* f2 = new TF2("f2", fun2, mmin, mmax, mmin, mmax, npar);

  f2->SetParameters(f2params);

  f2->SetParName(0, "\nSigSig");
  f2->SetParName(1, "Mass");
  f2->SetParName(2, "Sigma^2");
  f2->SetParName(3, "\nSigBkg");
  f2->SetParName(4, "Slope");
  f2->SetParName(5, "\nBkgBkg");
  f2->SetParName(6, "\nQuad");

  // f2->FixParameter(0,sigsig);
  // f2->FixParameter(1,mass);
  // f2->FixParameter(2,sigma2);
  // f2->FixParameter(3,sig_bkg);
  // f2->FixParameter(4,slope);
  // f2->FixParameter(5,bkg_bkg);
  // f2->FixParameter(6,0.);

  TFile::Open("AnalysisResults_for_fit.root");
  TH2F* h2 = (TH2F*)gFile->Get("hf-task-correlation-dplus-dplus-reduced/hMassDPair");
  h2->GetXaxis()->SetRangeUser(mmin, mmax);
  h2->GetYaxis()->SetRangeUser(mmin, mmax);

  h2->Fit("f2", "R", "lego2");

  auto ff = (TF2*)h2->GetFunction("f2");
  ff->SetBit(1 << 9); // Do not draw the function associated with the histogram
  ff->SetLineWidth(1);
  ff->SetNpx(50);
  ff->SetNpy(50);

  ff->Draw("surf same");
}
