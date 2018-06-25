#include "resolution.h"
void drawV_2() {
	
	TFile *results = new TFile("results.root");

	TH1D *h_reso_V0A = new TH1D("h_reso_V0A", "h_reso", NC+2, CentBins);

	for (int i=0; i<NC; i++){
		
		TH1D *pah = (TH1D*)results->Get(Form("h_V0A_a0%i", i));
		TH1D *pbh = (TH1D*)results->Get(Form("h_V0A_b0%i", i));
		TH1D *pch = (TH1D*)results->Get(Form("h_V0A_c0%i", i));
	
		h_reso_V0A->SetBinContent(i+1, pah->GetMean() * pbh->GetMean() / pch->GetMean());
	}

	
	// draw resulution calculated with different methods

	TProfile *h_resolution_V0A = (TProfile*)results->Get("h_resolution_V0A");
	h_resolution_V0A->SetTitle("Resolution");
	TCanvas *c2 = new TCanvas("c2", "Resolution");
	gStyle->SetOptStat(0);
	h_resolution_V0A->Draw();
	
	h_reso_V0A->SetLineColor(kRed);
	h_reso_V0A->Draw("same");

	h_resolution_V0A->GetXaxis()->SetTitle("Centrality");

   	TLegend *leg = new TLegend();
        leg->AddEntry(h_resolution_V0A,"Resolution, event plane difference method");
	leg->AddEntry(h_reso_V0A, "Resolution, subevent method");
	leg->Draw();

	
	// draw observed v2, corrected v2 and true v2

	TH1D *p_v2 = (TH1D*)results->Get("p_v2_corr_V0A");
	p_v2->SetTitle("v2");
	TCanvas *c1 = new TCanvas("c1", "v2");
	gStyle->SetOptStat(0);
	p_v2->SetLineColor(kRed);
	p_v2->Draw();
	
	TProfile *h_v2_V0A = (TProfile*)results->Get("h_v2_V0A");
	h_v2_V0A->Draw("same");
	
	TProfile *true_v2 = (TProfile*)results->Get("true_v2");
	true_v2->SetLineColor(kGreen);
	true_v2->Draw("same");

	h_v2_V0A->GetXaxis()->SetTitle("Centrality");
   	
	TLegend *leg2 = new TLegend();
        leg2->AddEntry(h_v2_V0A,"Observed v2");
	leg2->AddEntry(p_v2, "Corrected v2");
	leg2->AddEntry(true_v2, "True v2 (input)");
	leg2->Draw();
	
	TCanvas *ca= new TCanvas("ca", "a");
	h_v2_V0A->Draw();
	TH1D *h = h_v2_V0A->ProjectionX();
	h->SetLineColor(kGreen);
	h->Draw("same");
}
