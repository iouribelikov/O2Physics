#include "TF2.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"

Double_t sig1_sig2(Double_t* x, Double_t* par)
{
  Double_t r1 = Double_t((x[0] - par[1]));
  Double_t r2 = Double_t((x[1] - par[1]));
  return par[0] * TMath::Exp(-0.5 * (r1 * r1 / par[2] + r2 * r2 / par[2]));
}

Double_t sig1_bkg2(Double_t* x, Double_t* par)
{
  Double_t r1 = Double_t((x[0] - par[1]));
  return par[3] * TMath::Exp(-0.5 * r1 * r1 / par[2]) * (par[4] + x[1]);
}

Double_t sig2_bkg1(Double_t* x, Double_t* par)
{
  Double_t r1 = Double_t((x[1] - par[1]));
  return par[5] * TMath::Exp(-0.5 * r1 * r1 / par[2]) * (par[6] + x[0]);
}

Double_t bkg1_bkg2(Double_t* x, Double_t* par)
{
  return par[7] * (par[6] + x[0]) * (par[4] + x[1]);
}

Double_t fun2(Double_t* x, Double_t* par)
{
  return sig1_sig2(x, par) + sig1_bkg2(x, par) + sig2_bkg1(x, par) + bkg1_bkg2(x, par);
}

void fit2d()
{
  Double_t mass = 1.871, sigma2 = 0.009 * 0.009;
  Double_t offset1 = 1760., slope1 = -694.;
  Double_t offset2 = 1793., slope2 = -632.;
  Double_t f2params[] = {
    260, mass, sigma2,
    1.3 * slope2, offset2 / slope2,
    1.8 * slope1, offset1 / slope1,
    3e-3 * slope1 * slope2};
  const Int_t npar = sizeof(f2params) / sizeof(Double_t);

  // Float_t mmin=1.78, mmax=1.97;
  Float_t mmin = 1.80, mmax = 1.94;
  TF2* f2 = new TF2("f2", fun2, mmin, mmax, mmin, mmax, npar);

  f2->SetParameters(f2params);

  f2->SetParName(0, "Sig1_Sig2");
  f2->SetParName(1, "Mass");
  f2->SetParName(2, "Sigma^2");
  f2->SetParName(3, "Sig1_Bkg2");
  f2->SetParName(4, "Offset2");
  f2->SetParName(5, "Sig2_Bkg1");
  f2->SetParName(6, "Offset1");
  f2->SetParName(7, "Bkg1_Bkg2");

  // f2->FixParameter(1,mass);
  // f2->FixParameter(2,sigma2);
  // f2->FixParameter(4,offset2);
  // f2->FixParameter(5,slope2);
  // f2->FixParameter(6,offset1);

  TFile::Open("AnalysisResults.root");
  TH2F* h2 = (TH2F*)gFile->Get("hf-task-correlation-dplus-dplus-reduced/hMassDMesonPair");
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
