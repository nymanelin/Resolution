#include "resolution.h"
void draw() {
	
	TFile *results = new TFile("results.root");

	/* Explanations:
	 * t_ - TProfile
	 * h_ - TH1D
	 * p_ - Projection
	 * reso - resolution
	 * sub - sub-event method
	 * obs - observed
	 */

	TProfile *t_reso_true[R_COUNT];
	TH1D *h_reso_true[R_COUNT];
	TH1D *h_reso_sub[R_COUNT];

	TProfile *t_v2_true;
	TProfile *t_v2_obs[R_COUNT];
	TH1D *p_v2_corrtrue[R_COUNT];
	TH1D *p_v2_corrsub[R_COUNT];
	TH1D *p_reso[R_COUNT];

	// loop over detectors
	for (uint i=0; i<R_COUNT; ++i) {

		// load histograms
		
		t_reso_true[i] = (TProfile*)results->Get(Form("t_reso_true_%s", presn[i]));
		
		t_v2_obs[i] = (TProfile*)results->Get(Form("t_v2_obs_%s", presn[i]));
		
		t_v2_true = (TProfile*)results->Get("t_v2_true");

		h_reso_true[i] = new TH1D(Form("h_reso_true_%s", presn[i]), Form("h_reso_true_%s", presn[i]), NC+2, CentBins);

		h_reso_sub[i] = new TH1D(Form("h_reso_sub_%s", presn[i]), Form("h_reso_sub_%s", presn[i]), NC+2, CentBins);

		for (uint j=0; j<NC; ++j){
			
			TH1D *pah = (TH1D*)results->Get(Form("h_%s_a0%i", presn[i], j));
			TH1D *pbh = (TH1D*)results->Get(Form("h_%s_b0%i", presn[i], j));
			TH1D *pch = (TH1D*)results->Get(Form("h_%s_c0%i", presn[i], j));
		
			h_reso_sub[i]->SetBinContent(j+1, pah->GetMean() * pbh->GetMean() / pch->GetMean());
			h_reso_sub[i]->SetBinError(j+1, pah->GetMeanError() * pbh->GetMeanError() / pch->GetMeanError());
			
			TH1D *hist_reso = (TH1D*)results->Get(Form("h_reso_true_%s_%02u", presn[i], j));

			h_reso_true[i]->SetBinContent(j+1, hist_reso->GetMean());
			h_reso_true[i]->SetBinError(j+1, hist_reso->GetMeanError());
		}

		p_v2_corrtrue[i] = t_v2_obs[i]->ProjectionX(Form("p_v2_corrtrue_%s", presn[i]));
		p_reso[i] = t_reso_true[i]->ProjectionX("p_reso");
		p_v2_corrtrue[i]->Divide(p_reso[i]);

		p_v2_corrsub[i] = t_v2_obs[i]->ProjectionX(Form("p_v2_corrtrue_%s", presn[i]));
		p_v2_corrsub[i]->Divide(h_reso_sub[i]);

	}

	// draw resolution calculated with different methods
	
	TCanvas *c1 = new TCanvas("c1", "Resolution");
	c1->Divide(3,1);

	for (int i=0; i<R_COUNT; i++) {
		
		c1->cd(i+1);
		gStyle->SetOptStat(0);
		t_reso_true[i]->Draw();
		t_reso_true[i]->GetXaxis()->SetTitle("Centrality");
		
		h_reso_sub[i]->SetLineColor(kRed);
		h_reso_sub[i]->Draw("same");
		//h_reso_sub[i]->Draw("e");

		h_reso_true[i]->SetLineColor(kGreen);
		h_reso_true[i]->Draw("same");
		//h_reso_true[i]->Draw("e");


		TLegend *leg = new TLegend();
		leg->AddEntry(t_reso_true[i],"Resolution, event plane difference method (profile)");
		leg->AddEntry(h_reso_sub[i], "Resolution, subevent method");
		leg->AddEntry(h_reso_true[i], "Resolution, event plane difference method (mean of histograms)");
		leg->Draw();
	}
		
		
	// draw  observed v2, corrected v2 and true v2
	
	TCanvas *c2 = new TCanvas("c2", "v2");
	c2->Divide(3,1);

	for (int i=0; i<R_COUNT; i++) {

		c2->cd(i+1);
		gStyle->SetOptStat(0);
		p_v2_corrtrue[i]->Draw();
		p_v2_corrtrue[i]->GetXaxis()->SetTitle("Centrality");
		
		t_v2_obs[i]->SetLineColor(kRed);
		t_v2_obs[i]->Draw("same");
		
		t_v2_true->SetLineColor(kGreen);
		t_v2_true->Draw("same");

		p_v2_corrsub[i]->SetLineColor(kOrange);
		p_v2_corrsub[i]->Draw("same");
		
		TLegend *leg2 = new TLegend();
		leg2->AddEntry(t_v2_obs[i],"Observed v2");
		leg2->AddEntry(p_v2_corrtrue[i], "Corrected v2, true resolution");
		leg2->AddEntry(t_v2_true, "True v2 (input)");
		leg2->AddEntry(p_v2_corrsub[i], "Corrected v2, subevent resolution");
		leg2->Draw();
	}

}
