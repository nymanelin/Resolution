
#define NC 6             //Changed for high pT test
#define NSectors 6
static const char *presn[] = {"V0A","V0C","V0P"};

enum RESOLUTION{
	R_V0A,
	R_V0C,
	R_V0P,
	R_COUNT
};

bool saveflow = kFALSE;
bool savejet = kFALSE;

static double CentBins[NC+3] = {0,5,10,20,30,40,50,60,70};

void SubtractBg(TH1D *h, double bg, double ebg) ;


normalize() {

  double pi = TMath::Pi();
  string name;
  TFile *f1 = new TFile("results.root");                     // Resolution & smearing histos
  TFile *f2 = new TFile("normal.root", "recreate");  // Output after normalization
  TH2D *rawfile[R_COUNT][NC];

  TH1D *pah[R_COUNT][NC], *pbh[R_COUNT][NC], *pch[R_COUNT][NC];

  TH1D *resolution[R_COUNT];


  for (int s = 0; s < R_COUNT; s++) {
    resolution[s] = new TH1D(Form("resolution_%s", presn[s]), Form("resolution_%s", presn[s]), NC, CentBins);
  }
  for (int s = 0; s < R_COUNT; s++) {
    for (int c = 0; c < NC; c++) {
			pah[s][c]=(TH1D*)f1->Get(Form("h_%s_a%02u",presn[s],c));
			pbh[s][c]=(TH1D*)f1->Get(Form("h_%s_b%02u",presn[s],c));
			pch[s][c]=(TH1D*)f1->Get(Form("h_%s_c%02u",presn[s],c));

//      evpdifference[s][c] = (TH1D*) f1->Get(Form("h_evpdiff_%s_%02d", presn[s], c));
//      evpdifference[s][c]->GetXaxis()->SetTitle("Event Plane Difference #Psi_{true}-#Psi{reco}");
//      evpdifference[s][c]->GetYaxis()->SetTitleOffset(1.3);
//      evpdifference[s][c]->SetTitle("");
      
      resolution[s]->SetBinContent(c+1, TMath::Sqrt(pah[s][c]->GetMean()*pbh[s][c]->GetMean()/pch[s][c]->GetMean()));



    }
  }


  f2->cd();
  TCanvas *c1 = new TCanvas();
  TCanvas *c2 = new TCanvas();
 
  /* Drawing and Writing to file */


  for (int s = 0; s < R_COUNT; s++) {
    resolution[s]->Write();

    resolution[s]->Draw();
        if (saveflow)
    c1->SaveAs(Form("plots/Resolution%s.eps", presn[s]));

  }

}

void SubtractBg(TH1D *h, double bg, double ebg){
  int nb =h->GetNbinsX();
  TString hname = h->GetName();

  for(int ib=1; ib<=nb; ib++){
    double val = h->GetBinContent(ib);
    double err = h->GetBinError(ib);
    h->SetBinContent(ib,val-bg);
    h->SetBinError(ib,sqrt(err*err +  ebg * ebg));
  }

}
