#include "TF2.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"

Double_t sig1_sig2(Double_t* x, Double_t* par)
{
  auto r1 = (x[0] - par[1]);
  auto r2 = (x[1] - par[1]);
  return par[0] * TMath::Exp(-0.5 * (r1 * r1 / par[2] + r2 * r2 / par[2]));
}

Double_t sig1_bkg2(Double_t* x, Double_t* par)
{
  auto r1 = (x[0] - par[1]);
  auto r2 = (x[1] - par[1]);
  return par[3] * TMath::Exp(-0.5 * r1 * r1 / par[2]) * (1 + r2 * par[4] + r2 * r2 * par[8]);
}

Double_t sig2_bkg1(Double_t* x, Double_t* par)
{
  auto r1 = (x[0] - par[1]);
  auto r2 = (x[1] - par[1]);
  return par[5] * TMath::Exp(-0.5 * r2 * r2 / par[2]) * (1 + r1 * par[6] + r1 * r1 * par[9]);
}

Double_t bkg1_bkg2(Double_t* x, Double_t* par)
{
  auto r1 = (x[0] - par[1]);
  auto r2 = (x[1] - par[1]);
  return par[7] * (1 + r1 * par[6] + r1 * r1 * par[9]) * (1 + r2 * par[4] + r2 * r2 * par[8]);
}

Double_t fun2(Double_t* x, Double_t* par)
{
  return sig1_sig2(x, par) + sig1_bkg2(x, par) + sig2_bkg1(x, par) + bkg1_bkg2(x, par);
}

void fit2d_asym()
{
  Double_t sig1sig2 = 135., mass = 1.871, sigma2 = 0.009 * 0.009;
  Double_t sig1bkg2 = 600., slope2 = -1.0;
  Double_t sig2bkg1 = 900., slope1 = -1.4;
  Double_t bkg1bkg2 = 750.;
  Double_t f2params[] = {
    sig1sig2, mass, sigma2,
    sig1bkg2, slope2,
    sig2bkg1, slope1,
    bkg1bkg2,
    0., 0., // Quadratic coefficients
  };
  const Int_t npar = sizeof(f2params) / sizeof(Double_t);

  Float_t mmin = 1.80, mmax = 1.94;
  TF2* f2 = new TF2("f2", fun2, mmin, mmax, mmin, mmax, npar);

  f2->SetParameters(f2params);

  f2->SetParName(0, "\nSig1Sig2");
  f2->SetParName(1, "Mass");
  f2->SetParName(2, "Sigma^2");
  f2->SetParName(3, "\nSig1Bkg2");
  f2->SetParName(4, "Slope2");
  f2->SetParName(5, "\nSig2Bkg1");
  f2->SetParName(6, "Slope1");
  f2->SetParName(7, "\nBkg1Bkg2");

  // f2->FixParameter(0,sig1sig2);
  // f2->FixParameter(1,mass);
  // f2->FixParameter(2,sigma2);
  // f2->FixParameter(3,sig1bkg2);
  // f2->FixParameter(4,slope2);
  // f2->FixParameter(5,sig2bkg1);
  // f2->FixParameter(6,slope1);
  // f2->FixParameter(7,bkg1bkg2);
  // f2->FixParameter(8,0);
  // f2->FixParameter(9,0);

  TFile::Open("AnalysisResults_for_fit.root");
  TH2F* h2 = (TH2F*)gFile->Get("hf-task-correlation-dplus-dplus-reduced/hMassDplusPair");
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
