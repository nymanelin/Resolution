
#include <stdio.h>
#include <vector>
#include <TRandom.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TComplex.h>
#include "resolution.h"

typedef unsigned int uint;

#define N_VN 2


//void resolution(){
int main(int argc, char **pargv){
  uint seed = argc > 1?atol(pargv[1]):1000;
  uint evtc = argc > 2?atol(pargv[2]):1000;
  printf("seed:\t%u\nevents:\t%u\n",seed,evtc);

  //Read flow coefficients ----------------------------
  const char *pgn[3] = {"Hist1D_y1","Hist1D_y1_e1","Hist1D_y1_e2"};
  const char *pglabel[4] = {"Hist1D_y%u","Hist1D_y%u_e1","Hist1D_y%u_e2plus", "Hist1D_y%u_e2minus"};
  TFile *pff = new TFile("anizo-flow.root","read");
  TFile *highf = new TFile("276_HighPtFlow.root");

  TGraphErrors *pgr_v[N_VN];
  TGraphErrors *high_v[N_VN+1];

  for(uint i = 0; i < N_VN; ++i){
    TH1D *pgr[N_VN];//[3];
    for(uint j = 0; j < 2; ++j)
      pgr[j] = (TH1D*)pff->Get(Form("Table %u/%s",i+1,pgn[j]));

    uint n = pgr[0]->GetNbinsX();
    pgr_v[i] = new TGraphErrors(n);
    for(uint ic = 0; ic < n; ++ic){
      pgr_v[i]->SetPoint(ic,
          pgr[0]->GetBinCenter(ic+1),pgr[0]->GetBinContent(ic+1));
      pgr_v[i]->SetPointError(ic,0,pgr[1]->GetBinContent(ic+1));
    }
  }

  for (uint i = 0; i < N_VN+1; ++i) {
    high_v[i] =  new TGraphErrors();
    for (uint c = 0; c < NC; ++c) {
      TH1D *temp_hist[5];
      for (int ilbl = 0; ilbl < 3; ilbl++) {
        temp_hist[ilbl] = (TH1D*)highf->Get(Form("Table %u/%s", ilbl+1, Form(pglabel[ilbl], c+1)));
      }  // y1 : v2, y2 : v2{4}, y3 : v3;

      high_v[i]->SetPoint(high_v[i]->GetN(), (CentBins[c] + CentBins[c+1])/2, temp_hist[0]->GetBinContent( temp_hist[0]->FindBin(8.1)) );
      high_v[i]->SetPointError(high_v[i]->GetN() - 1 , 0, temp_hist[1]->FindBin(8.1));
    }
  }

  pff->Close();
  delete pff;
  highf->Close();
  delete highf;

  // PDF for Pure flow, High_pT flow, and Jet
  TF1 *pdf = new TF1("df",
      "[0]*(1+2*[1]*cos(x-[5])+2*[2]*cos(2*(x-[6]))+2*[3]*cos(3*(x-[7]))+2*[4]*cos(4*(x-[8])))",-pi,pi);
  pdf->SetParameter(0,100.0); // v0
  pdf->SetParameter(1,0.0); // v1

  TF1 *pdf_high = new TF1("dfhigh",
      "[0]*(1+2*[1]*cos(x-[5])+2*[2]*cos(2*(x-[6]))+2*[3]*cos(3*(x-[7]))+2*[4]*cos(4*(x-[8])))",-pi,pi);
  pdf_high->SetParameter(0,100.0); // v0
  pdf_high->SetParameter(1,0.0); // v1

  TF2 *jetpdf = new TF2("jetpdf", "[0]*exp(-0.5*((x-[4])/[1])**2)*exp(-0.5*((y-[2])/[3])**2)",-0.8, 0.8, -pi, pi);
  jetpdf->SetParameter(0, 100.0);
  jetpdf->SetParameter(1, 0.1); // eta
  jetpdf->SetParameter(3, 0.1); //phi

  TRandom *prng = new TRandom(seed);

  for (uint j = 0; j < NC; ++j) {
  }

  for(uint i = 0; i < R_COUNT; ++i){
    for(uint j = 0; j < NC; ++j){
      pah[i][j] = new TH1D(Form("h_res_%s_a%02u",presn[i],j),"h_res",1024,-1.5,1.5);    // resolution factor for detector i, i-a
      pbh[i][j] = new TH1D(Form("h_res_%s_b%02u",presn[i],j),"h_res",1024,-1.5,1.5);    // resolution factor for detector i, i-b
      pch[i][j] = new TH1D(Form("h_res_%s_c%02u",presn[i],j),"h_res",1024,-1.5,1.5);    // resolution factor for detector i, i-c
      
      evph[i][j] = new TH1D(Form("h_evp_%s_%02u", presn[i], j), "h_evp", 64, -1.6, 1.6);

      evpdifference[i][j] = new TH1D(Form("h_evpdiff_%s_c%02u", presn[i], j), Form("Event Plane Difference %s, %.0f-%.0fC", presn[i], CentBins[j], CentBins[j+1]), 128, -4, 4); 
      evpdifference[i][j]->GetXaxis()->SetTitle("Event plane Difference");
      evpdifference[i][j]->GetYaxis()->SetTitle("N Events");

    }
  }

  //Create multiplicity(eta) distribution & integrate ---------------
  /*double etadst[ETADST_N];
    for(uint i = 0; i < ETADST_N; ++i)
    etadst[i] = -3.5+(double)i*0.25+0.125;*/

  TGraph *pgr_nch[D_COUNT];
  for(uint i = 0; i < D_COUNT; ++i)
    pgr_nch[i] = new TGraph(NC);

  for(uint i = 0; i < NC; ++i){
    TGraph gr_eta(ETADST_N,etadst,etanch[i]);
    TF1 f("etach",[&](double *px, double *pp)->double{
        return gr_eta.Eval(px[0]);
        },-3.5,5.1,0);

    for(uint j = 0; j < D_COUNT; ++j){
      double nch = f.Integral(cov[j][0],cov[j][1])/(cov[j][1]-cov[j][0]);
      pgr_nch[j]->SetPoint(i,0.5*(CentBins[i]+CentBins[i+1]),nch);
    }
  }

  //#define MAX_N 10000
  //double eta[MAX_N], phi[MAX_N];
  //double phi[MAX_N];

  for(uint evt = 0; evt < evtc; ++evt){
    //Event generation ----------------------------

    double cent = prng->Uniform(0,50.0);
    pdf->SetParameter(2,pgr_v[0]->Eval(cent) * 0.65); //v2
    pdf->SetParameter(3,pgr_v[1]->Eval(cent)); //v3
    pdf->SetParameter(4,0.01);//pgr_v[2]->Eval(cent)); //v4

    pdf->SetParameter(5,prng->Uniform(-pi,pi));  //EP for v1
    pdf->SetParameter(6,prng->Uniform(-pi/2.0,pi/2.0)); //EP for v2
    pdf->SetParameter(7,prng->Uniform(-pi/3.0,pi/3.0)); //EP for v3
    pdf->SetParameter(8,prng->Uniform(-pi/4.0,pi/4.0)); //EP for v4

    pdf_high->SetParameter(2,high_v[0]->Eval(cent) * 0.65); //v2
    pdf_high->SetParameter(3,high_v[2]->Eval(cent)); //v3
    pdf_high->SetParameter(4,0.01);//pgr_v[2]->Eval(cent)); //v4

    pdf_high->SetParameter(5,pdf->GetParameter(5));  //EP for v1
    pdf_high->SetParameter(6,pdf->GetParameter(6)); //EP for v2
    pdf_high->SetParameter(7,pdf->GetParameter(7)); //EP for v3
    pdf_high->SetParameter(8,pdf->GetParameter(8)); //EP for v4


    uint cid = 0;
    for(; cid < NC; ++cid)
      if(cent < CentBins[cid+1])
        break;

    //Q-vectors -----------------------------------
    double trueevp = pdf->GetParameter(6);

    //containers of Pure flow tracks and jet tracks
    std::vector <double> trackphi[D_COUNT];
    std::vector <double> tracketa[D_COUNT];
    std::vector <double> jettrackphi[D_COUNT];
    std::vector <double> jettracketa[D_COUNT];
    TComplex Qsd[D_COUNT];
    for(uint s = 0; s < D_COUNT; ++s){
      uint ntracks = (uint)pgr_nch[s]->Eval(cent) *0.9;

      TComplex Qa2 = TComplex(0,0);
      for(uint i = 0; i < ntracks; ++i){
        double tphi = pdf->GetRandom();
        //				phi = TMath::Floor(8.0*(phi+pi)/(2.0*pi))*(2.0*pi)/8-pi;
        double phi = CheckDetectorPhi(tphi);
        //        double phi = tphi;
        trackphi[s].push_back(phi);
        tracketa[s].push_back(prng->Uniform(-0.8, 0.8));
        Qa2 += TComplex(TMath::Cos(2.0*phi),TMath::Sin(2.0*phi));
      }
      
      //      double thphi = pdf_high->GetRandom();
      //      double hphi = CheckDetectorPhi(thphi);
      //      trackphi[s].push_back(hphi);
      //			Qa2 += TComplex(TMath::Cos(2.0*hphi),TMath::Sin(2.0*hphi));

      //		Qa2 /= (double) (ntracks+1);   // for highpt
      Qa2 /= (double) (ntracks);
      Qsd[s] = Qa2/TComplex::Abs(Qa2);
    }

    //Calculate Evp using Q Vectors

    for (uint s = 0; s < R_COUNT; ++s) {
      double recoevp = TMath::ATan2(Qsd[s].Im(), Qsd[s].Re())/2;
      double evpdiff = trueevp - recoevp;

      evph[s][cid]->Fill(recoevp);
      evpdifference[s][cid]->Fill(evpdiff);
    }

    //    evpcorrvsdet2d[cid]->Fill(TMath::ATan2(Qsd[0].Im(), Qsd[0].Re()), TMath::ATan2(Qsd[1].Im(), Qsd[1].Re()));


    //Calculate the resolution components
    TComplex ab, ac, bc;
    /*ab = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC]);
      ac = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_V0C]);
      bc = Qsd[D_TPC]*TComplex::Conjugate(Qsd[D_V0C]);
      pah[R_V0A][cid]->Fill(ab.Re());
      pbh[R_V0A][cid]->Fill(ac.Re());
      pch[R_V0A][cid]->Fill(bc.Re());

      ab = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_V0A]);
      ac = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_TPC]);
      bc = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC]);
      pah[R_V0C][cid]->Fill(ab.Re());
      pbh[R_V0C][cid]->Fill(ac.Re());
      pch[R_V0C][cid]->Fill(bc.Re());*/
    ab = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
    ac = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
    bc = Qsd[D_TPC_ETAA]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
    pah[R_V0A][cid]->Fill(ab.Re());
    pbh[R_V0A][cid]->Fill(ac.Re());
    pch[R_V0A][cid]->Fill(bc.Re());

    ab = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
    ac = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
    //bc = Qsd[D_TPC_ETAA]*TComplex::Conjugate(Qsd[D_TPC]);
    pah[R_V0C][cid]->Fill(ab.Re());
    pbh[R_V0C][cid]->Fill(ac.Re());
    pch[R_V0C][cid]->Fill(bc.Re());

    ab = Qsd[D_V0P]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
    ac = Qsd[D_V0P]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
    bc = Qsd[D_TPC_ETAA]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
    pah[R_V0P][cid]->Fill(ab.Re());
    pbh[R_V0P][cid]->Fill(ac.Re());
    pch[R_V0P][cid]->Fill(bc.Re());

  }

  for(uint i = 0; i < D_COUNT; ++i)
    delete pgr_nch[i];
  for(uint i = 0; i < N_VN; ++i) {
    delete pgr_v[i];
    delete high_v[i];
  }

  delete prng;
  delete pdf;
  delete pdf_high;

  TFile *pfo = new TFile(argc > 3?pargv[3]:"results.root","recreate");
  pfo->cd();

  for(uint i = 0; i < R_COUNT; ++i){
    for(uint j = 0; j < NC; ++j){
      //double a[2] = {pah[i]->GetMean(),pah[i]->GetMeanError()};
      //double b[2] = {pbh[i]->GetMean(),pah[i]->GetMeanError()};
      //double c[2] = {pch[i]->GetMean(),pah[i]->GetMeanError()};
      //
      //double R = TMath::Sqrt(a[0]*b[0]/c[0]);
      //double e = TMath::Sqrt(b[0]*a[1]*a[1]/(a[0]*c[0])
      //	+a[0]*b[1]*b[1]/(b[0]*c[0])+a[0]*b[0]*c[1]/(c[0]*c[0]*c[0]));
      //printf("R2(%u) = %lf pm %lf\n",i,R,e);

      pah[i][j]->Write(Form("h_%s_a%02u",presn[i],j));
      pbh[i][j]->Write(Form("h_%s_b%02u",presn[i],j));
      pch[i][j]->Write(Form("h_%s_c%02u",presn[i],j));

      evph[i][j]->Write(Form("h_evp_%s_%02u", presn[i], j));
      evpdifference[i][j]->Write(Form("h_evpdiff_%s_%02u", presn[i], j));




      delete pah[i][j];
      delete pbh[i][j];
      delete pch[i][j];
      delete evph[i][j];
      delete evpdifference[i][j];
    }
  }

  pfo->Close();
  delete pfo;

  return 0;
}

