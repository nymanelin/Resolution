#define NFiles 1

compareresolution() {

  TFile *f1[NFiles];
//  f1[0]  = new TFile("normal_highpt.root");
//  f1[1]  = new TFile("normal_lower05.root");
//  f1[2]  = new TFile("normal_lower10.root");
//  f1[3]  = new TFile("normal_lowerv2_40.root");
//  f1[4]  = new TFile("normal_lowerv2_45.root");
//  f1[5]  = new TFile("normal_lowerv2_50.root");
//  f1[6]  = new TFile("normal_lowerv2_45_true.root"); 
  f1[0]  = new TFile("normal_lowerv2_35.root");

  TFile *reference502 = new TFile("v2_res_V0AC_alice5TeVPbPb_alex.root");

  TH1D *resolution_referenceV0C = (TH1D*)reference502->Get("hResV0C");
  TH1D *resolution_referenceV0A = (TH1D*)reference502->Get("hResV0A");

  string legendentry[NFiles] = { //"No tweaks",
//    "-5% Multiplicity, -5% v_{2}", 
//    "-10% Multiplicity, -10% v_{2}", 
//    "-10% Multiplicity, -40% v_{2}",
//    "-10% Multiplicity, -45% v_{2}",
//    "-10% Multiplicity, -50% v_{2}", 
//    "-10% Multiplicity, -45% v_{2}, Resolution Off", 
      "Realistic Toy MC"       // "-10% Multiplicity, -35% v_{2}"
    };

  TH1D *resolutionV0C[NFiles];
  TH1D *resolutionV0A[NFiles];

  gStyle->SetLegendBorderSize(0);
  TLegend *leg = new TLegend(0.3, 0.25, 0.8, 0.48);


  gStyle->SetLegendBorderSize(0);
  TCanvas *c1 = new TCanvas();
  TCanvas *c2 = new TCanvas();

  for (int ifile = 0; ifile < NFiles; ifile++) {
    c1->cd();
    resolutionV0C[ifile] = (TH1D*)f1[ifile]->Get("resolution_V0C");
    resolutionV0C[ifile]->SetTitle("");
   
    resolutionV0C[ifile]->SetXTitle("Centrality (%)");
    resolutionV0C[ifile]->GetXaxis()->SetTitleSize(0.05);
    resolutionV0C[ifile]->GetXaxis()->SetLabelSize(0.05);
    resolutionV0C[ifile]->GetYaxis()->SetTitleSize(0.05);
    resolutionV0C[ifile]->GetYaxis()->SetLabelSize(0.05);
    resolutionV0C[ifile]->SetYTitle("Resolution");
    resolutionV0C[ifile]->SetStats(0);
    resolutionV0C[ifile]->SetLineColor(ifile+2);
    resolutionV0C[ifile]->GetYaxis()->SetRangeUser(0.0, 0.99);
    if (ifile >= 2) 
      resolutionV0C[ifile]->SetLineColor(ifile+2);
    if (ifile+2 == 5) 
      resolutionV0C[ifile]->SetLineColor(8);
    if (ifile+2 == 6) 
      resolutionV0C[ifile]->SetLineColor(15);
    if (ifile+2 == 8) 
      resolutionV0C[ifile]->SetLineColor(kViolet+3);

    leg->AddEntry(resolutionV0C[ifile], legendentry[ifile].c_str(), "l");
    if (ifile == 0)
      resolutionV0C[ifile]->Draw();
    else 
      resolutionV0C[ifile]->Draw("same");

    c2->cd();
    resolutionV0A[ifile] = (TH1D*)f1[ifile]->Get("resolution_V0A");
    resolutionV0A[ifile]->SetTitle("V0A Resolution Comparison");
    resolutionV0A[ifile]->SetXTitle("Centrality (%)");
    resolutionV0A[ifile]->SetYTitle("Resolution");
    resolutionV0A[ifile]->SetStats(0);
    resolutionV0A[ifile]->SetLineColor(ifile+1);
    resolutionV0A[ifile]->GetYaxis()->SetRangeUser(0.0, 0.99);
    if (ifile >= 2) 
      resolutionV0A[ifile]->SetLineColor(ifile+2);
    if (ifile+2 == 5) 
      resolutionV0A[ifile]->SetLineColor(8);
    if (ifile+2 == 6) 
      resolutionV0A[ifile]->SetLineColor(15);
    if (ifile+2 == 8) 
      resolutionV0A[ifile]->SetLineColor(kViolet+3);

    if (ifile == 0)
      resolutionV0A[ifile]->Draw();
    else 
      resolutionV0A[ifile]->Draw("same");
  }

  resolution_referenceV0A->Draw("same");
  leg->AddEntry(resolution_referenceV0A, "Reference from ALICE data", "pl");
  leg->Draw();
  c2->SaveAs("CompareResolution_V0A.eps");

  c1->cd();
  resolution_referenceV0C->Draw("same");
  TLatex l1;
  l1.SetTextSize(0.04);
  l1.DrawLatex(8, 0.54, "V0C Resolution Comparison from toyMC #sqrt{s_{NN}}=5.02TeV");
  l1.DrawLatex(8, 0.50, "V0C Acceptance -3.7<#eta<-1.7");

  leg->Draw();
  c1->SaveAs("CompareResolution_V0C.eps");

}

    
