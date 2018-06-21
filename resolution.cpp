#include <stdio.h>
#include <vector>
#include <TRandom.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include "TProfile.h"
#include <TComplex.h>
#include "resolution.h"

typedef unsigned int uint;

#define N_VN 2  // Number of Harmonics n=2, 3


//void resolution(){
int main(int argc, char **pargv){
  uint seed = argc > 1?atol(pargv[1]):1000;
  uint evtc = argc > 2?atol(pargv[2]):1000;
  printf("seed:\t%u\nevents:\t%u\n",seed,evtc);

  //Read flow coefficients ----------------------------
  const char *pgn[3] = {"Hist1D_y1","Hist1D_y1_e1","Hist1D_y1_e2"};
  const char *pglabel[4] = {"Hist1D_y%u","Hist1D_y%u_e1","Hist1D_y%u_e2plus", "Hist1D_y%u_e2minus"};
  TFile *pff = new TFile("anizo-flow.root","read");

  TGraphErrors *pgr_v[N_VN];

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

  pff->Close();
  delete pff;

  // PDF for Pure flow
  TF1 *pdf = new TF1("df",
      "[0]*(1+2*[1]*cos(x-[5])+2*[2]*cos(2*(x-[6]))+2*[3]*cos(3*(x-[7]))+2*[4]*cos(4*(x-[8])))",-pi,pi);
  pdf->SetParameter(0,100.0); // const/norm 
  pdf->SetParameter(1,0.0); // v1

  // Random #
  TRandom *prng = new TRandom(seed);

  // Declare histograms...
  // 
  for(uint i = 0; i < R_COUNT; ++i){ // R_COUNT is the number of detectors
    for(uint j = 0; j < NC; ++j){
      pah[i][j] = new TH1D(Form("h_res_%s_a%02u",presn[i],j),"h_res",1024,-1.5,1.5);    // resolution factor for detector i, i-a
      pbh[i][j] = new TH1D(Form("h_res_%s_b%02u",presn[i],j),"h_res",1024,-1.5,1.5);    // resolution factor for detector i, i-b
      pch[i][j] = new TH1D(Form("h_res_%s_c%02u",presn[i],j),"h_res",1024,-1.5,1.5);    // resolution factor for detector i, i-c
      
      evph[i][j] = new TH1D(Form("h_evp_%s_%02u", presn[i], j), "h_evp", 64, -1.6, 1.6);

      evpdifference[i][j] = new TH1D(Form("h_evpdiff_%s_c%02u", presn[i], j), Form("Event Plane Difference %s, %.0f-%.0fC", presn[i], CentBins[j], CentBins[j+1]), 128, -4, 4); 
      evpdifference[i][j]->GetXaxis()->SetTitle("Event plane Difference");
      evpdifference[i][j]->GetYaxis()->SetTitle("N Events");
    
    }
    
    h_v2[i] = new TProfile(Form("h_v2_%s", presn[i]), "h_v2", NC+2, CentBins);
    h_resolution[i] = new TProfile(Form("h_resolution_%s", presn[i]), "h_resolution", NC+2, CentBins);
  }

  TProfile *true_v2 = new TProfile("true_v2", "true_v2", NC+2, CentBins);

  // How many particles we want to put into the simulation.
  // Create multiplicity (eta) distribution & integrate ---------------
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

  for(uint evt = 0; evt < evtc; ++evt){
    //Event generation ----------------------------

    double cent = prng->Uniform(0,50.0);
    // adding vn..
    pdf->SetParameter(2,pgr_v[0]->Eval(cent) ); //v2 0.65 is for v2 weightening because of the pt distribution in the detectors.
    pdf->SetParameter(3,pgr_v[1]->Eval(cent)); //v3
    pdf->SetParameter(4,0.01);//pgr_v[2]->Eval(cent)); //v4
    // producing symmetric angles
    pdf->SetParameter(5,prng->Uniform(-pi,pi));  //EP for v1
    pdf->SetParameter(6,prng->Uniform(-pi/2.0,pi/2.0)); //EP for v2
    pdf->SetParameter(7,prng->Uniform(-pi/3.0,pi/3.0)); //EP for v3
    pdf->SetParameter(8,prng->Uniform(-pi/4.0,pi/4.0)); //EP for v4

    uint cid = 0;
    for(; cid < NC; ++cid)
      if(cent < CentBins[cid+1])
        break;

    //Q-vectors Only for 2nd order-----------------------------------
    double trueevp = pdf->GetParameter(6);
    
    true_v2->Fill(cent, pdf->GetParameter(2));

    //containers of Pure flow tracks and jet tracks
    std::vector <double> trackphi[D_COUNT];
    TComplex Qsd[D_COUNT]; // QVector per each detector 2nd order only
    // Calculating QVector per each detector
    for(uint s = 0; s < D_COUNT; ++s){
      uint ntracks = (uint)pgr_nch[s]->Eval(cent) *0.9;

      TComplex Qa2 = TComplex(0,0);
      for(uint i = 0; i < ntracks; ++i){
        double tphi = pdf->GetRandom();
        double phi = CheckDetectorPhi(tphi);
        trackphi[s].push_back(phi);
        Qa2 += TComplex(TMath::Cos(2.0*phi),TMath::Sin(2.0*phi));
      }
      Qa2 /= (double) (ntracks);
      Qsd[s] = Qa2/TComplex::Abs(Qa2);
    }

    //Calculate Event plane using Q Vectors

    for (uint s = 0; s < R_COUNT; ++s) {
      double recoevp = TMath::ATan2(Qsd[s].Im(), Qsd[s].Re())/2; 
      double evpdiff = trueevp - recoevp;

      evph[s][cid]->Fill(recoevp); // reconstructed EP
      evpdifference[s][cid]->Fill(evpdiff); // for resolution
      
      int ntracks = trackphi[s].size();
      for (uint i = 0; i < ntracks; ++i) {
	h_v2[s]->Fill(cent, TMath::Cos(2*(trackphi[s][i]-recoevp)));
        h_resolution[s]->Fill(cent, TMath::Cos(2*evpdiff));
      }
    }

    //Calculate the resolution components with different methods.
    TComplex ab, ac, bc;
    ab = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
    ac = Qsd[D_V0A]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
    bc = Qsd[D_TPC_ETAA]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
    pah[R_V0A][cid]->Fill(ab.Re());
    pbh[R_V0A][cid]->Fill(ac.Re());
    pch[R_V0A][cid]->Fill(bc.Re());

    ab = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_TPC_ETAA]);
    ac = Qsd[D_V0C]*TComplex::Conjugate(Qsd[D_TPC_ETAC]);
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

  for (int i=0; i<R_COUNT; i++) {
    p_v2[i] = h_v2[i]->ProjectionX(Form("p_v2_corr_%s", presn[i]));
    p_reso[i] = h_resolution[i]->ProjectionX("p_reso");
    p_v2[i]->Divide(p_reso[i]);
  }

  // cleaning up the variables
  for(uint i = 0; i < D_COUNT; ++i)
    delete pgr_nch[i];
  for(uint i = 0; i < N_VN; ++i) {
    delete pgr_v[i];
  }

  delete prng;
  delete pdf;

  TFile *pfo = new TFile(argc > 3?pargv[3]:"results.root","recreate");
  pfo->cd();

  true_v2->Write();

  for(uint i = 0; i < R_COUNT; ++i){
    h_v2[i]->Write();
    h_resolution[i]->Write();
    p_v2[i]->Write();
    for(uint j = 0; j < NC; ++j){
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

