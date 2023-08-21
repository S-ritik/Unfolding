double lumi = 59.83;

#include "My_Style.C"
#include "TUnfold_example.C"


TH1* removefirstbin(TH1* hist, int nbins_toremove = 1)
{
  TH1* hist_corrected = (TH1*)hist->Clone();
  const int nbins_corrected = hist->GetNbinsX() - nbins_toremove;
  const double *bins = hist->GetXaxis()->GetXbins()->GetArray();
  double bins_corrected[nbins_corrected+1];
  cout<<endl<<hist->GetName()<<" bins = "<<hist->GetNbinsX()<<"-"<<nbins_toremove<<"="<<nbins_corrected<<endl;
  for(int i=0; i < nbins_corrected;i++)
    {
      cout<<bins[i+nbins_toremove]<<" ";
      bins_corrected[i] = bins[i+nbins_toremove];
    }
  //cout<<"hist nbins earlier = "<<
  //hist_corrected->SetBins(nbins_corrected,bins_corrected);
  hist_corrected->SetBins(nbins_corrected,bins+nbins_toremove);

  for(int ib=0; ib < nbins_corrected; ib++)
    {
      hist_corrected->SetBinContent(ib+1, hist->GetBinContent(ib+1+nbins_toremove));
      hist_corrected->SetBinError(ib+1, hist->GetBinError(ib+1+nbins_toremove));
    }
  return hist_corrected;
}

void scale_bin_by_binwidth(TH1D *hist, double effective_binwidth = 0)
{
  if(effective_binwidth == 1)
    effective_binwidth = 1;
  int xbin = hist->GetNbinsX();
  for(int ix = 1; ix <= xbin; ix++)
    {
      float binwidth_factor;
      binwidth_factor = hist->GetBinWidth(ix)/effective_binwidth;
           
      hist->SetBinContent(ix, hist->GetBinContent(ix) / binwidth_factor);
      hist->SetBinError(ix, hist->GetBinError(ix) / binwidth_factor);
    }
  
  hist->GetYaxis()->SetTitle( TString::Format("Events / %.0f GeV",effective_binwidth) );
}

void plot_1dhists(int nhist, TH1D *hist_obs[nhist], string dataname[nhist],string plotname)
{
  
  TCanvas *cv;
  TLegend *legv;
  
  cv = tdrCanvas("canv_d",hist_obs[0],8,0);
  cv->cd();
  legv = tdrLeg(0.2,0.65,0.875,0.915);

  legv->SetTextFont(42);
  legv->SetTextSize(0.045);
  legv->SetBorderSize(0);
  gStyle->SetPaintTextFormat( "1.2f" );

  double max_bincontent = hist_obs[0]->GetMaximum();
  
  for(int ih = 0; ih < nhist; ih++)
    {
      //hist_obs[ih]->Scale(1.0/max(0.00001,hist_obs[ih]->Integral()));
      hist_obs[ih]->SetFillStyle(0);
      hist_obs[ih]->SetFillColor(0);
      hist_obs[ih]->SetLineColor(ih+2);
      //     hist_obs[ih]->SetLineStyle(ih+2);
      hist_obs[ih]->SetLineWidth(2);
      hist_obs[ih]->GetYaxis()->SetTitleOffset(1.5);

      legv->AddEntry(hist_obs[ih],dataname[ih].c_str(),"lep");

      max_bincontent = max(max_bincontent,double(hist_obs[ih]->GetMaximum()));
      hist_obs[ih]->SetMinimum(0.01);
    }

  hist_obs[0]->SetMaximum(1.5*max_bincontent);
  
  //gPad->SetLogy(1);

  hist_obs[0]->Draw("hist");
  for(int ih = 1; ih < nhist; ih++)
    hist_obs[ih]->Draw("histsame");

  legv->Draw("same");
  CMS_lumi( cv, 8, 0 );

  cv->SaveAs(("plots/" + plotname).c_str());
}

void plot_1dhists_withratio(int nhist, TH1D *hist_obs[nhist], string dataname[nhist],string plotname)
{
  TH1D *rat_obs;
  TH1D *rat_obs_2;
  TGraphAsymmErrors *ratioGraph;  
  TCanvas *cv;
  TLegend *legv;
  char name[200];

  double max_bincontent = hist_obs[0]->GetMaximum();
  legv = tdrLeg(0.22,0.6,0.55,0.865);

  for(int ih = 0; ih < nhist; ih++)
    {
      //hist_obs[ih]->Scale(1.0/max(0.00001,hist_obs[ih]->Integral()));
      hist_obs[ih]->SetFillStyle(0);
      hist_obs[ih]->SetFillColor(0);
      hist_obs[ih]->SetLineColor(ih+2);
      //     hist_obs[ih]->SetLineStyle(ih+2);
      hist_obs[ih]->SetLineWidth(2);
      hist_obs[ih]->GetYaxis()->SetTitleOffset(1.5);
      
      legv->AddEntry(hist_obs[ih],dataname[ih].c_str(),"lep");

      for(int ix = 1; ix <= hist_obs[ih]->GetNbinsX(); ix++)
	hist_obs[ih]->SetBinContent(ix, max(0.00001,hist_obs[ih]->GetBinContent(ix)));
	
      max_bincontent = max(max_bincontent,double(hist_obs[ih]->GetMaximum()));
      hist_obs[ih]->SetMinimum(0.0000001);
    }
  hist_obs[0]->SetMaximum(1.5*max_bincontent);
  
  cv = tdrDiCanvas("canv_d",hist_obs[0],hist_obs[1],8,0);
  cv->cd(1);

  legv->SetTextFont(42);
  legv->SetTextSize(0.045);
  legv->SetBorderSize(0);
  gStyle->SetPaintTextFormat( "1.2f" );
  
  //gPad->SetLogy(1);

  hist_obs[0]->Draw("hist");
  for(int ih = 1; ih < nhist; ih++)
    hist_obs[ih]->Draw("histsame");

  legv->Draw("same");
  CMS_lumi( cv, 8, 0 );

  cv->cd(2);

  rat_obs = (TH1D*)hist_obs[4]->Clone();
  rat_obs->Divide(hist_obs[3]);

  rat_obs_2 = (TH1D*)hist_obs[0]->Clone();
  rat_obs_2->Divide(hist_obs[2]);

  sprintf(name,"%s",hist_obs[0]->GetXaxis()->GetTitle());
  rat_obs->GetXaxis()->SetTitle(name);
  rat_obs->GetXaxis()->SetTitleSize(0.1);
  rat_obs->GetXaxis()->SetLabelSize(0.1);
  //rat_obs->GetXaxis()->CenterTitle();
  rat_obs->GetXaxis()->SetNdivisions(406);

  rat_obs->GetYaxis()->SetTitle("Before / After");
  rat_obs->GetYaxis()->SetTitleSize(0.1);
  rat_obs->GetYaxis()->SetLabelSize(0.1);
  rat_obs->GetYaxis()->SetTitleOffset(0.75);
  rat_obs->GetYaxis()->SetNdivisions(406);
  rat_obs->GetYaxis()->CenterTitle();

  rat_obs->SetMinimum(0);
  rat_obs->SetMaximum(2.);

  rat_obs->SetFillStyle(0);
  rat_obs->SetFillColor(0);
  rat_obs->SetMarkerStyle(kFullCircle);
  rat_obs->SetMarkerColor(kBlack);
  rat_obs->SetMarkerSize(0.7);
  rat_obs->SetLineColor(6);
  rat_obs_2->SetMarkerStyle(kFullCircle);
  rat_obs_2->SetMarkerColor(kBlack);
  rat_obs_2->SetMarkerSize(0.7);
  rat_obs_2->SetLineColor(4);
  rat_obs->Draw("e,p");
  rat_obs_2->Draw("e,psame");

  TLegend *legend_ratio = new TLegend(0.77,0.74,0.93,0.87);
  legend_ratio->SetBorderSize(0);
  legend_ratio->SetTextFont(42);
  legend_ratio->SetTextSize(0.045);
  legend_ratio->AddEntry(rat_obs, "Data","l");
  legend_ratio->AddEntry(rat_obs_2, "MC","l");
  legend_ratio->Draw("same");

  TLine *line = new TLine(rat_obs->GetXaxis()->GetXmin(),1,rat_obs->GetXaxis()->GetXmax(),1);
  line->SetLineColor(1);
  line->Draw("sames");
    
  cv->SaveAs(("plots/" + plotname).c_str());
}

void get_resolution_value(string mcfile, string histname)
{
  TFile* sigFile = TFile::Open(mcfile.c_str(), "READ");
  TH2D* hist2d_res = (TH2D*)sigFile->Get(histname.c_str());

  cout<<"Start calculation"<<endl;

  hist2d_res->Rebin2D(1,600);
    
  for(int igen = 1; igen < hist2d_res->GetNbinsX()+1; igen++)
    {
      if(hist2d_res->Integral(igen, igen, 1, hist2d_res->GetNbinsY()+1) < 1)
	continue;
      cout<<"Check for bin "<<igen<<" bin range "<<hist2d_res->GetXaxis()->GetBinLowEdge(igen)<<"-"<<hist2d_res->GetXaxis()->GetBinLowEdge(igen+1)<<" with inetgral ="<<hist2d_res->Integral(igen, igen, 1, hist2d_res->GetNbinsY()+1)<<endl;
      
      /*TH1D * hist = new TH1D("hist","",hist2d_res->GetNbinsY(),-5,5);
      //cout<<hist2d_res->GetNbinsY()<<endl;
      for(int i=0; i<hist2d_res->GetNbinsY(); i++)
	{
	  //cout<<hist2d_res->GetBinContent(i,igen)<<" ";
	  hist->SetBinContent(i+1,hist2d_res->GetBinContent(igen,i+1));
	  }*/

      TH1D * hist = hist2d_res->ProjectionY();
   
      cout<<" = "<<hist->Integral()<<endl;

      TCanvas *c1 = tdrCanvas("canv_d", hist2d_res->ProjectionX(),8,0);
      c1->cd();
      hist->Draw();
      
      TFitResultPtr r = hist->Fit("gaus","S");
      Double_t sigma   = r->Parameter(2);
      cout<<endl<<"Resolution for "<<igen<<" th reco bin = "<<sigma<<endl;

      gPad->Update();

      string ans;
      cout<<"Check the next hist y/n"<<endl;
      //cin>>ans;
      //if(ans == "n")
	break;
      
      delete hist;
      delete c1;
    }
  cout<<"End calculation"<<endl;
  sigFile->Clone();
}

void plot_response_matrix(TH2D *hist2d_res,string genvar, string gentitle, string recotitle)
{
  int xbin = hist2d_res->GetNbinsX();
  int ybin = hist2d_res->GetNbinsY();
  TMatrixD response_matrix(xbin, ybin);
  TH2D *hist2d_response_full = new TH2D("hist2d_response_full","hist2d_response_full", xbin+2,0,xbin+2, ybin+2,0,ybin+2);

  cout<<xbin<<" "<<ybin<<" "<<endl;

  for(int ib=0;ib<xbin+2 ;ib++)
    {
      TString binlabel;
      if(ib == 0) binlabel = "Inefficiency rate";
      else if(ib == xbin+1) binlabel.Form("%.1f - #infty",hist2d_res->GetXaxis()->GetBinLowEdge(ib));
      else  binlabel.Form("%.1f - %.1f",hist2d_res->GetXaxis()->GetBinLowEdge(ib),hist2d_res->GetXaxis()->GetBinUpEdge(ib));
      hist2d_response_full->GetXaxis()->SetBinLabel(ib+1, binlabel.Data());
    }

  for(int ib=0;ib<ybin+2 ;ib++)
    {
      TString binlabel;
      if(ib == 0) binlabel = "Fake rate";
      else if(ib == ybin+1) binlabel.Form("%.1f - #infty",hist2d_res->GetYaxis()->GetBinLowEdge(ib));
      else  binlabel.Form("%.1f - %.1f",hist2d_res->GetYaxis()->GetBinLowEdge(ib),hist2d_res->GetYaxis()->GetBinUpEdge(ib));
      hist2d_response_full->GetYaxis()->SetBinLabel(ib+1, binlabel.Data());
    }
  
  for(int iy = 0; iy < ybin+2; iy++)
    {
      for(int ix = 0; ix < xbin+2; ix++)
	{
	  double nevts = hist2d_res->GetBinContent(ix,iy);
	  if(nevts > 0.00001)
	    {
	      //if(ix !=0 && iy !=0) response_matrix[ix-1][iy-1] = nevts;
	      hist2d_response_full->SetBinContent(ix+1, iy +1, nevts);
	    }
	  else
	    {
	      //if(ix !=0 && iy !=0) response_matrix[ix-1][iy-1] = 0;
	      hist2d_response_full->SetBinContent(ix+1, iy +1, 0);
	    }
	}
    }

  TH2D * hist2d_smearing_matrix = (TH2D*)hist2d_res->Clone("smearing_matrix");
  for(int iy = 0; iy < ybin; iy++)
    {
      float y_integral = hist2d_res->Integral(1, xbin+1, iy +1, iy +1);
      for(int ix = 0; ix < xbin; ix++)
	{
	  if(y_integral > 0.0001 && hist2d_res->GetBinContent(ix+1, iy +1) > 0.0000001)
	    hist2d_smearing_matrix->SetBinContent(ix+1, iy +1, hist2d_res->GetBinContent(ix+1, iy +1) /y_integral );
	  else
	    hist2d_smearing_matrix->SetBinContent(ix+1, iy +1, 0.0000001 );
	  response_matrix[ix][iy] = hist2d_smearing_matrix->GetBinContent(ix+1, iy +1);
	}
    }

  /*for(int iy = 2; iy < ybin+2; iy++)
    {
      for(int ix = 2; ix < xbin+2; ix++)
	{
	  double nevts = hist2d_response_full->GetBinContent(ix,iy);
	  response_matrix[ix-2][iy-2] = nevts;
	}
	}*/

  for(int iy = 0; iy < ybin; iy++)
    {
      for(int ix = 0; ix < xbin; ix++)
	{
  	  cout<<response_matrix[ix][iy]<<"\t";
	}
      cout<<endl;
    }
  
  TDecompSVD * singular_matrix = new TDecompSVD(response_matrix);
  TVectorD singular_values = singular_matrix->GetSig();
  double max_value = singular_values.Max();
  double min_value = singular_values.Min();
  double condition_no = max_value/max(0.000001,min_value);
  cout<<"Max eigen value is "<<max_value<<" Min eigen value is "<<min_value<<" Condition number is "<<condition_no<<endl;

  gStyle->SetPaintTextFormat( "1.2f" );
  gStyle->SetOptStat(0);
  hist2d_response_full->GetYaxis()->SetTitle(gentitle.c_str());
  hist2d_response_full->GetXaxis()->SetTitle(recotitle.c_str());

  TPaveText *text = new TPaveText(0.1,0.85,0.5,0.9,"NDC");
  text->AddText(TString::Format("Condition number = %.2f",condition_no));
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  //text->SetFillColorAlpha();    

  TCanvas *c1 = tdrCanvas("canv_d", hist2d_response_full->ProjectionX(),8,0);

  c1->cd();
  gStyle->SetOptTitle(0);

  hist2d_response_full->GetXaxis()->SetTickSize(0.03);
  hist2d_response_full->GetYaxis()->SetTickSize(0.03);

  hist2d_response_full->GetXaxis()->SetTitleSize(0.03);
  hist2d_response_full->GetXaxis()->SetTitleOffset(1.25);
  hist2d_response_full->GetXaxis()->SetLabelSize(0.03);
  hist2d_response_full->GetXaxis()->SetLabelOffset(0.001);

  hist2d_response_full->GetYaxis()->SetTitleSize(0.035);
  hist2d_response_full->GetYaxis()->SetTitleOffset(1.7);
  hist2d_response_full->GetYaxis()->SetLabelSize(0.025);

  hist2d_response_full->SetMarkerSize(1.4);
  gStyle->SetPaintTextFormat( "1.0f" );
  
  hist2d_response_full->Draw("coltext45");
  CMS_lumi( c1, 8, 0 );
  //text->Draw("same");
  c1->SaveAs(TString::Format("plots/response_matrix_%s_.png",genvar.c_str()));

  TCanvas *c2 = tdrCanvas("canv_d", hist2d_smearing_matrix->ProjectionX(),8,0);

  c2->cd();
  gStyle->SetOptTitle(0);

  hist2d_smearing_matrix->GetXaxis()->SetTickSize(0.03);
  hist2d_smearing_matrix->GetYaxis()->SetTickSize(0.03);

  hist2d_smearing_matrix->GetXaxis()->SetTitleSize(0.04);
  hist2d_smearing_matrix->GetXaxis()->SetTitleOffset(1.0);
  hist2d_smearing_matrix->GetXaxis()->SetLabelSize(0.035);
  hist2d_smearing_matrix->GetXaxis()->SetLabelOffset(0.001);

  hist2d_smearing_matrix->GetYaxis()->SetTitleSize(0.037);
  hist2d_smearing_matrix->GetYaxis()->SetTitleOffset(1.2);
  hist2d_smearing_matrix->GetYaxis()->SetLabelSize(0.035);

  hist2d_smearing_matrix->GetYaxis()->SetTitle(gentitle.c_str());
  hist2d_smearing_matrix->GetXaxis()->SetTitle(recotitle.c_str());

  gStyle->SetPaintTextFormat( "1.2f" );
  
  c2->SetRightMargin( 0.09 );

  hist2d_smearing_matrix->Draw("colztext45");
  CMS_lumi( c2, 8, 0 );
  text->Draw("same");
  c2->SaveAs(TString::Format("plots/smearing_matrix_%s_.png",genvar.c_str()));
  
}
double* get_jackknife_stat_error_ondata(string genvar, string gentitle, string recovar)
{
  TFile* sigFile = TFile::Open("MuMu_Summer20UL18_unfolding_input.root", "READ");
  TH1D* hGen_allbins[11];
  hGen_allbins[0] = (TH1D*)sigFile->Get(("hist_gen_" + genvar).c_str());
  hGen_allbins[0]->SetDirectory(0);

  for(int ih=1; ih < 11; ih++)
    {
      hGen_allbins[ih] = (TH1D*)sigFile->Get(TString::Format("hist_gen_%s_set%d",genvar.c_str(),ih).Data());
      hGen_allbins[ih]->SetDirectory(0);
    }
  
  TH1D* hReco_allbins[11];
  hReco_allbins[0] = (TH1D*)sigFile->Get(TString::Format("hist_reco_%s",recovar.c_str()).Data());
  hReco_allbins[0]->SetDirectory(0);
  
  for(int ih=1; ih < 11; ih++)
    {
      hReco_allbins[ih] = (TH1D*)sigFile->Get(TString::Format("hist_reco_%s_set%d",recovar.c_str(),ih).Data());
      hReco_allbins[ih]->SetDirectory(0);
    }
  for(int ih=0; ih < 11; ih++)
    {
      cout<<hReco_allbins[ih]->Integral()<<endl;
    }
  
  TH2D* hResp_allbins = (TH2D*)sigFile->Get(("hist_responsematrix_recogen_" + genvar).c_str());
  hResp_allbins->SetDirectory(0);
  sigFile->Close();

  TH1D* hGen[11];
  TH1D* hReco[11];
  for(int ih=0; ih < 11; ih++)
    {
      hReco[ih] = (TH1D*)removefirstbin(hReco_allbins[ih],2);
      hGen[ih] = (TH1D*)removefirstbin(hGen_allbins[ih]);
    }
  
  cout<<"end removing bins"<<endl;
  TH1D *hist_unfold[11];
      
  for(int ih=0; ih < 11; ih++)
    {
      TFile* combine_output;
      if(ih==0)
	combine_output = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar + ".root").c_str(), "READ");
      else
	combine_output = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar + "_set" + to_string(ih) + ".root").c_str(), "READ");
      
      TTree * tree = (TTree*)combine_output->Get("limit");
      const int nbins = int(hGen[ih]->GetNbinsX());
      float r_values[nbins];
      float r_bins[nbins];
      float r_bins_high[nbins];
      float r_bins_low[nbins];
      
      for(int ib = 0; ib < nbins; ib++)
	tree->SetBranchAddress(TString::Format("r_bin%d",ib+2),&r_values[ib]);
      
      for(int i=0; i <nbins; i++)   
	{
	  tree->GetEntry(0);
	  r_bins[i] = r_values[i];
	  tree->GetEntry(2*i+3); 
	  r_bins_low[i] = r_values[i];
	  tree->GetEntry(2*i+4); 
	  r_bins_high[i] = r_values[i];      
	}
      tree->GetEntry(0);   //////// get best fit values
      
      hist_unfold[ih] = (TH1D*)hGen[ih]->Clone();
      TH1D *hist_err = (TH1D*)hGen[ih]->Clone();
      TH1D *rat_obs_combine;
      TH1D *rat_obs_tunfold;
      
      TH1D *hist_unfold_tunfold_allbins;
      if(ih == 0)   hist_unfold_tunfold_allbins = (TH1D*)TUnfold_unfold(genvar,gentitle,recovar);
      else          hist_unfold_tunfold_allbins = (TH1D*)TUnfold_unfold(genvar+"_set"+to_string(ih),gentitle,recovar+"_set"+to_string(ih));
      TH1D *hist_unfold_tunfold = (TH1D*)removefirstbin(hist_unfold_tunfold_allbins);
      
      cout<<"Check "<<ih<<endl;
      for(int ib=0; ib < nbins; ib++)
	{
	  float bincontent = 0;//hResp_allbins->Integral(1,hResp_allbins->GetNbinsX(),ib+1,ib+1)   //hGen->GetBinContent(ib+1)
	  for(int ix = 1; ix < hResp_allbins->GetNbinsX()+1; ix++)
	    {
	      if(hResp_allbins->GetBinContent(ix,ib+2) > 0.00001)
		bincontent += hResp_allbins->GetBinContent(ix,ib+2);
	    }

	  cout<<r_bins[ib]<<" "<<r_bins_low[ib]-r_bins[ib]<<" "<<r_bins_high[ib]-r_bins[ib]<<endl;
	  cout<<"gen bin "<<ib<<" content"<<hGen[ih]->GetBinContent(ib+1)<<" "<<hGen[ih]->GetBinLowEdge(ib+1)<<" == "<<bincontent<<endl;

	  hist_unfold[ih]->SetBinContent(ib+1, bincontent*r_bins[ib]);
	  hist_unfold[ih]->SetBinError(ib+1, bincontent*(r_bins_high[ib] - r_bins_low[ib])/2);
	  
	  hist_err->SetBinContent(ib+1, abs(bincontent*r_bins[ib] - hGen[ih]->GetBinContent(ib+1))/hGen[ih]->GetBinContent(ib+1));
	}

      hist_unfold[ih]->SetFillStyle(0);
      hist_unfold[ih]->SetFillColor(0);
      hist_unfold[ih]->SetMarkerStyle(kFullCircle);
      hist_unfold[ih]->SetMarkerColor(kBlack);
      hist_unfold[ih]->SetMarkerSize(1.2);
      hist_unfold[ih]->SetLineColor(kBlack);
      hist_unfold[ih]->SetLineWidth(2);
      hist_unfold[ih]->SetLineStyle(1);
      hist_unfold[ih]->SetLineColor(0);

      hist_unfold_tunfold->SetMarkerStyle(kFullStar);
      hist_unfold_tunfold->SetMarkerColor(kBlack);
      hist_unfold_tunfold->SetLineStyle(1);
      hist_unfold_tunfold->SetMarkerSize(3);
      hist_unfold_tunfold->SetLineWidth(2);
      hist_unfold_tunfold->SetLineColor(kRed);
      
      hGen[ih]->GetXaxis()->SetTitle(gentitle.c_str());
      hGen[ih]->SetLineWidth(3);
      //hGen[ih]->SetLineStyle(2);
      hGen[ih]->SetLineColor(kBlue);
  
      hReco[ih]->SetLineWidth(3);
      //hReco->SetLineStyle(2);
      hReco[ih]->SetLineColor(kGreen);
      hReco[ih]->Rebin(2);

      scale_bin_by_binwidth(hist_unfold[ih],1);
      scale_bin_by_binwidth(hist_unfold_tunfold,1);
      scale_bin_by_binwidth(hGen[ih],1);
      scale_bin_by_binwidth(hReco[ih],1);
      
  
      //hGen->SetMaximum(1.9*max(max(hGen->GetMaximum(),hReco->GetMaximum()),hist_unfold->GetMaximum()));

      hist_err->SetFillStyle(0);
      //hist_err->SetMarkerColor(0);
      hist_err->SetLineColor(0);
      hist_err->SetMarkerSize(1.2);
      
      TCanvas *c1 = tdrDiCanvas("c1",hGen[ih],hReco[ih],8,0);;
      c1->cd(1);
      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      //gPad->SetLogy();
      gStyle->SetPaintTextFormat( "1.2f" );
      
      TLegend *legend = new TLegend(0.34,0.6,0.74,0.85);
      legend->SetBorderSize(0);
      legend->SetTextFont(42);
      legend->SetTextSize(0.035);
      legend->AddEntry(hist_unfold[ih], "Unfolded  distribution using Combine", "pe2");
      legend->AddEntry(hist_unfold_tunfold, "Unfolded  distribution using Tunfold", "l");
      legend->AddEntry(hGen[ih], "Gen Madgraph+Pythia sample", "l");
      legend->AddEntry(hReco[ih], "Reco Madgraph+Pythia sample", "l");
      
      hGen[ih]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen[ih]->GetBinContent(hGen[ih]->GetMaximumBin()),hReco[ih]->GetBinContent(hReco[ih]->GetMaximumBin())),hist_unfold[ih]->GetBinContent(hist_unfold[ih]->GetMaximumBin())));
      
      hGen[ih]->Draw("e1");
      hReco[ih]->Draw("e1same");
      hist_unfold_tunfold->Draw("histsame");
      hist_unfold[ih]->DrawCopy("Psame");
      hist_unfold[ih]->SetFillColor(kBlack);
      hist_unfold[ih]->SetFillStyle(3013);
      hist_unfold[ih]->Draw("e2same");
      //hist_unfold->Draw("Psame");
      //hist_err->Draw("Psame");
      legend->Draw("same");

      CMS_lumi( c1, 8, 0 );
      gPad->Modified();


      c1->cd(2);
      
      rat_obs_combine = (TH1D*)hist_unfold[ih]->Clone();
      rat_obs_combine->Divide(hGen[ih]);
      
      rat_obs_tunfold = (TH1D*)hist_unfold_tunfold->Clone();
      rat_obs_tunfold->Divide(hGen[ih]);
      rat_obs_tunfold->SetLineColor(kRed);
      
      rat_obs_combine->GetXaxis()->SetTitle(gentitle.c_str());
      rat_obs_combine->GetXaxis()->SetTitleSize(0.09);
      rat_obs_combine->GetXaxis()->SetLabelSize(0.1);
      rat_obs_combine->GetXaxis()->CenterTitle();
      rat_obs_combine->GetXaxis()->SetNdivisions(406);

      rat_obs_combine->GetYaxis()->SetTitle("Unfolded / Gen");
      rat_obs_combine->GetYaxis()->SetTitleSize(0.1);
      rat_obs_combine->GetYaxis()->SetLabelSize(0.1);
      rat_obs_combine->GetYaxis()->SetTitleOffset(0.75);
      rat_obs_combine->GetYaxis()->SetNdivisions(406);
      rat_obs_combine->GetYaxis()->CenterTitle();
  
      rat_obs_combine->SetMinimum(0);
      rat_obs_combine->SetMaximum(2.);
      
      rat_obs_combine->SetFillStyle(0);
      rat_obs_combine->SetFillColor(0);
      rat_obs_combine->SetMarkerStyle(kFullCircle);
      rat_obs_combine->SetMarkerColor(kBlack);
      rat_obs_combine->SetMarkerSize(0.7);
      rat_obs_combine->SetLineColor(kBlack);
      
      rat_obs_tunfold->SetFillStyle(0);
      rat_obs_tunfold->SetFillColor(0);
      rat_obs_tunfold->SetMarkerStyle(kFullCircle);
      rat_obs_tunfold->SetMarkerColor(kRed);
      rat_obs_tunfold->SetMarkerSize(0.7);
      rat_obs_tunfold->SetLineColor(kRed);
      
      rat_obs_combine->Draw("e,p");
      rat_obs_tunfold->Draw("esame");
      rat_obs_combine->Draw("e,psame");
  
      TLegend *legend_ratio = new TLegend(0.77,0.74,0.93,0.87);
      legend_ratio->SetBorderSize(0);
      legend_ratio->SetTextFont(42);
      legend_ratio->SetTextSize(0.045);
      legend_ratio->AddEntry(rat_obs_combine, "Combine","l");
      legend_ratio->AddEntry(rat_obs_tunfold, "Tunfold","l");
      legend_ratio->Draw("same");
      
      TLine *line = new TLine(rat_obs_combine->GetXaxis()->GetXmin(),1,rat_obs_combine->GetXaxis()->GetXmax(),1);
      line->SetLineColor(1);
      line->Draw("same");
      
      if(ih==0) c1->SaveAs(TString::Format("plots/withoutfirstbin_data_unfolded_distribution_%s_higgscombine.png",genvar.c_str()));
      else c1->SaveAs(TString::Format("plots/withoutfirstbin_data_unfolded_distribution_%s_set%d_higgscombine.png",genvar.c_str(),ih));

    }
  for(int ih=0; ih < 11; ih++)
    {	   
      hist_unfold[ih]->SetFillStyle(0);
      hist_unfold[ih]->SetFillColor(0);
      hist_unfold[ih]->SetMarkerStyle(kFullCircle);
      hist_unfold[ih]->SetMarkerColor(kBlack);
      hist_unfold[ih]->SetMarkerSize(1.2);
      hist_unfold[ih]->SetLineColor(kBlack);
      hist_unfold[ih]->SetLineWidth(2);
      hist_unfold[ih]->SetLineStyle(1);
      hist_unfold[ih]->SetLineColor(ih);
      
      hist_unfold[ih]->GetXaxis()->SetTitle(gentitle.c_str());
     
      //scale_bin_by_binwidth(hist_unfold[ih],1);
    }
  
  TCanvas *c1 = tdrCanvas("c1",hist_unfold[0],8,0);;
  c1->cd();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gStyle->SetPaintTextFormat( "1.2f" );

  TLegend *legend = new TLegend(0.3,0.6,0.7,0.85);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
  legend->AddEntry(hist_unfold[0], "Unfolded  distribution", "pe2");
  for(int ih=1; ih < 11; ih++)
    legend->AddEntry(hist_unfold[ih], ("Unfolded  distribution using data set " + to_string(ih)).c_str(), "l");
  
  //hGen->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen->GetBinContent(hGen->GetMaximumBin()),hReco->GetBinContent(hReco->GetMaximumBin())),hist_unfold->GetBinContent(hist_unfold->GetMaximumBin())));

  hist_unfold[0]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*hist_unfold[0]->GetBinContent(hist_unfold[0]->GetMaximumBin()));

  //hGen->Draw("e1");
  //hReco->Draw("e1same");
  //hist_unfold_tunfold->Draw("histsame");
  hist_unfold[0]->DrawCopy("Psame");
  hist_unfold[0]->SetFillColor(kBlack);
  hist_unfold[0]->SetFillStyle(3013);
  hist_unfold[0]->Draw("e2same");
  for(int ih=1; ih < 11; ih++)
    hist_unfold[ih]->Draw("histsame");
  //hist_unfold->Draw("Psame");
  //hist_err->Draw("Psame");
  legend->Draw("same");

  CMS_lumi( c1, 8, 0 );
  gPad->Modified();
  c1->SaveAs(TString::Format("plots/jackknife_data_unfolded_distribution_comparsion_%s_higgscombine.png",genvar.c_str()));

  double jackknife_errors[11];

  TH1D *hist_unfold_jacknife_unc = (TH1D*)hist_unfold[0]->Clone();
  for(int ib =0; ib < hist_unfold_jacknife_unc->GetNbinsX(); ib++)
    {
      std::vector<double> bincontents;
      for(int ih =1; ih<11; ih++)
	bincontents.push_back(hist_unfold[ih]->GetBinContent(ib+1));
      double std_dev = TMath::StdDev(bincontents.begin(), bincontents.end() ); 

      hist_unfold_jacknife_unc->SetBinError(ib+1,3*std_dev);
      cout<<"Bin "<<ib<<": Content = "<<hist_unfold_jacknife_unc->GetBinContent(ib+1)<<" with combine uncertainty = "<<hist_unfold[0]->GetBinError(ib+1)<<" and jackknife uncertainty = "<<hist_unfold_jacknife_unc->GetBinError(ib+1)<<endl;
    }

  TCanvas *c2 = tdrCanvas("c2",hist_unfold[0],8,0);;
  c2->cd();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gStyle->SetPaintTextFormat( "1.2f" );

  TLegend *legend_2 = new TLegend(0.3,0.6,0.7,0.85);
  legend_2->SetBorderSize(0);
  legend_2->SetTextFont(42);
  legend_2->SetTextSize(0.035);
  legend_2->AddEntry(hist_unfold[0], "Unfolded  distribution with combine uncertainity", "pe2");
  legend_2->AddEntry(hist_unfold_jacknife_unc, "Unfolded  distribution with jackknife data uncertainity", "l");
   
  //hGen->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen->GetBinContent(hGen->GetMaximumBin()),hReco->GetBinContent(hReco->GetMaximumBin())),hist_unfold->GetBinContent(hist_unfold->GetMaximumBin())));

  hist_unfold[0]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*hist_unfold[0]->GetBinContent(hist_unfold[0]->GetMaximumBin()));

  //hGen->Draw("e1");
  //hReco->Draw("e1same");
  //hist_unfold_tunfold->Draw("histsame");
  hist_unfold[0]->DrawCopy("Psame");
  hist_unfold[0]->SetFillColor(kBlack);
  hist_unfold[0]->SetFillStyle(3013);
  hist_unfold[0]->Draw("e2same");
  hist_unfold_jacknife_unc->SetLineColor(kRed);
  hist_unfold_jacknife_unc->Draw("e1same");
  legend_2->Draw("same");

  CMS_lumi( c2, 8, 0 );
  gPad->Modified();
  c2->SaveAs(TString::Format("plots/jackknife_compare_datauncertainty_unfolded_distribution_%s_combinevsjackknife.png",genvar.c_str()));
  return jackknife_errors;
}


double* get_jackknife_stat_error_onresponsematrix(string genvar, string gentitle, string recovar)
{
  TFile* sigFile = TFile::Open("MuMu_Summer20UL18_unfolding_input.root", "READ");
  TH1D* hGen_allbins[11];
  hGen_allbins[0] = (TH1D*)sigFile->Get(("hist_gen_" + genvar).c_str());
  hGen_allbins[0]->SetDirectory(0);

  for(int ih=1; ih < 11; ih++)
    {
      hGen_allbins[ih] = (TH1D*)sigFile->Get(TString::Format("hist_gen_%s_set%d",genvar.c_str(),ih+10).Data());
      hGen_allbins[ih]->SetDirectory(0);
    }
  
  TH1D* hReco_allbins[11];
  hReco_allbins[0] = (TH1D*)sigFile->Get(TString::Format("hist_reco_%s",recovar.c_str()).Data());
  hReco_allbins[0]->SetDirectory(0);
  
  for(int ih=1; ih < 11; ih++)
    {
      hReco_allbins[ih] = (TH1D*)sigFile->Get(TString::Format("hist_reco_%s_set%d",recovar.c_str(),ih+10).Data());
      hReco_allbins[ih]->SetDirectory(0);
    }
  for(int ih=0; ih < 11; ih++)
    {
      cout<<hReco_allbins[ih]->Integral()<<endl;
    }
  
  TH2D* hResp_allbins[11];
  hResp_allbins[0]= (TH2D*)sigFile->Get(("hist_responsematrix_recogen_" + genvar).c_str());
  hResp_allbins[0]->SetDirectory(0);

  for(int ih=1; ih < 11; ih++)
    {
      hResp_allbins[ih]= (TH2D*)sigFile->Get(TString::Format("hist_responsematrix_recogen_%s_set%d",genvar.c_str(),ih+10).Data());
      hResp_allbins[ih]->SetDirectory(0);
    }
  
  sigFile->Close();

  TH1D* hGen[11];
  TH1D* hReco[11];
  for(int ih=0; ih < 11; ih++)
    {
      hReco[ih] = (TH1D*)removefirstbin(hReco_allbins[ih],2);
      hGen[ih] = (TH1D*)removefirstbin(hGen_allbins[ih]);
    }
  
  cout<<"end removing bins"<<endl;
  TH1D *hist_unfold[11];
      
  for(int ih=0; ih < 11; ih++)
    {
      TFile* combine_output;
      if(ih==0)
	combine_output = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar + ".root").c_str(), "READ");
      else
	combine_output = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar + "_set" + to_string(ih+10) + ".root").c_str(), "READ");
      
      TTree * tree = (TTree*)combine_output->Get("limit");
      const int nbins = int(hGen[ih]->GetNbinsX());
      float r_values[nbins];
      float r_bins[nbins];
      float r_bins_high[nbins];
      float r_bins_low[nbins];
      
      for(int ib = 0; ib < nbins; ib++)
	tree->SetBranchAddress(TString::Format("r_bin%d",ib+2),&r_values[ib]);
      
      for(int i=0; i <nbins; i++)   
	{
	  tree->GetEntry(0);
	  r_bins[i] = r_values[i];
	  tree->GetEntry(2*i+3); 
	  r_bins_low[i] = r_values[i];
	  tree->GetEntry(2*i+4); 
	  r_bins_high[i] = r_values[i];      
	}
      tree->GetEntry(0);   //////// get best fit values
      
      hist_unfold[ih] = (TH1D*)hGen[ih]->Clone();
      TH1D *hist_err = (TH1D*)hGen[ih]->Clone();
      TH1D *rat_obs_combine;
      TH1D *rat_obs_tunfold;
      
      TH1D *hist_unfold_tunfold_allbins;
      if(ih == 0)   hist_unfold_tunfold_allbins = (TH1D*)TUnfold_unfold(genvar,gentitle,recovar);
      else          hist_unfold_tunfold_allbins = (TH1D*)TUnfold_unfold(genvar+"_set"+to_string(ih+10),gentitle,recovar+"_set"+to_string(ih+10));
      TH1D *hist_unfold_tunfold = (TH1D*)removefirstbin(hist_unfold_tunfold_allbins);
      
      cout<<"Check "<<ih<<endl;
      for(int ib=0; ib < nbins; ib++)
	{
	  float bincontent = 0;//hResp_allbins[ih]->Integral(1,hResp_allbins[ih]->GetNbinsX(),ib+1,ib+1)   //hGen->GetBinContent(ib+1)
	  for(int ix = 1; ix < hResp_allbins[ih]->GetNbinsX()+1; ix++)
	    {
	      if(hResp_allbins[ih]->GetBinContent(ix,ib+2) > 0.00001)
		bincontent += hResp_allbins[ih]->GetBinContent(ix,ib+2);
	    }

	  cout<<r_bins[ib]<<" "<<r_bins_low[ib]-r_bins[ib]<<" "<<r_bins_high[ib]-r_bins[ib]<<endl;
	  cout<<"gen bin "<<ib<<" content"<<hGen[ih]->GetBinContent(ib+1)<<" "<<hGen[ih]->GetBinLowEdge(ib+1)<<" == "<<bincontent<<endl;

	  hist_unfold[ih]->SetBinContent(ib+1, bincontent*r_bins[ib]);
	  hist_unfold[ih]->SetBinError(ib+1, bincontent*(r_bins_high[ib] - r_bins_low[ib])/2);
	  
	  hist_err->SetBinContent(ib+1, abs(bincontent*r_bins[ib] - hGen[ih]->GetBinContent(ib+1))/hGen[ih]->GetBinContent(ib+1));
	}

      hist_unfold[ih]->SetFillStyle(0);
      hist_unfold[ih]->SetFillColor(0);
      hist_unfold[ih]->SetMarkerStyle(kFullCircle);
      hist_unfold[ih]->SetMarkerColor(kBlack);
      hist_unfold[ih]->SetMarkerSize(1.2);
      hist_unfold[ih]->SetLineColor(kBlack);
      hist_unfold[ih]->SetLineWidth(2);
      hist_unfold[ih]->SetLineStyle(1);
      hist_unfold[ih]->SetLineColor(0);

      hist_unfold_tunfold->SetMarkerStyle(kFullStar);
      hist_unfold_tunfold->SetMarkerColor(kBlack);
      hist_unfold_tunfold->SetLineStyle(1);
      hist_unfold_tunfold->SetMarkerSize(3);
      hist_unfold_tunfold->SetLineWidth(2);
      hist_unfold_tunfold->SetLineColor(kRed);
      
      hGen[ih]->GetXaxis()->SetTitle(gentitle.c_str());
      hGen[ih]->SetLineWidth(3);
      //hGen[ih]->SetLineStyle(2);
      hGen[ih]->SetLineColor(kBlue);
  
      hReco[ih]->SetLineWidth(3);
      //hReco->SetLineStyle(2);
      hReco[ih]->SetLineColor(kGreen);
      hReco[ih]->Rebin(2);

      scale_bin_by_binwidth(hist_unfold[ih],1);
      scale_bin_by_binwidth(hist_unfold_tunfold,1);
      scale_bin_by_binwidth(hGen[ih],1);
      scale_bin_by_binwidth(hReco[ih],1);
      
  
      //hGen->SetMaximum(1.9*max(max(hGen->GetMaximum(),hReco->GetMaximum()),hist_unfold->GetMaximum()));

      hist_err->SetFillStyle(0);
      //hist_err->SetMarkerColor(0);
      hist_err->SetLineColor(0);
      hist_err->SetMarkerSize(1.2);
      
      TCanvas *c1 = tdrDiCanvas("c1",hGen[ih],hReco[ih],8,0);;
      c1->cd(1);
      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      //gPad->SetLogy();
      gStyle->SetPaintTextFormat( "1.2f" );
      
      TLegend *legend = new TLegend(0.34,0.6,0.74,0.85);
      legend->SetBorderSize(0);
      legend->SetTextFont(42);
      legend->SetTextSize(0.035);
      legend->AddEntry(hist_unfold[ih], "Unfolded  distribution using Combine", "pe2");
      legend->AddEntry(hist_unfold_tunfold, "Unfolded  distribution using Tunfold", "l");
      legend->AddEntry(hGen[ih], "Gen Madgraph+Pythia sample", "l");
      legend->AddEntry(hReco[ih], "Reco Madgraph+Pythia sample", "l");
      
      hGen[ih]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen[ih]->GetBinContent(hGen[ih]->GetMaximumBin()),hReco[ih]->GetBinContent(hReco[ih]->GetMaximumBin())),hist_unfold[ih]->GetBinContent(hist_unfold[ih]->GetMaximumBin())));
      
      hGen[ih]->Draw("e1");
      hReco[ih]->Draw("e1same");
      hist_unfold_tunfold->Draw("histsame");
      hist_unfold[ih]->DrawCopy("Psame");
      hist_unfold[ih]->SetFillColor(kBlack);
      hist_unfold[ih]->SetFillStyle(3013);
      hist_unfold[ih]->Draw("e2same");
      //hist_unfold->Draw("Psame");
      //hist_err->Draw("Psame");
      legend->Draw("same");

      CMS_lumi( c1, 8, 0 );
      gPad->Modified();


      c1->cd(2);
      
      rat_obs_combine = (TH1D*)hist_unfold[ih]->Clone();
      rat_obs_combine->Divide(hGen[ih]);
      
      rat_obs_tunfold = (TH1D*)hist_unfold_tunfold->Clone();
      rat_obs_tunfold->Divide(hGen[ih]);
      rat_obs_tunfold->SetLineColor(kRed);
      
      rat_obs_combine->GetXaxis()->SetTitle(gentitle.c_str());
      rat_obs_combine->GetXaxis()->SetTitleSize(0.09);
      rat_obs_combine->GetXaxis()->SetLabelSize(0.1);
      rat_obs_combine->GetXaxis()->CenterTitle();
      rat_obs_combine->GetXaxis()->SetNdivisions(406);

      rat_obs_combine->GetYaxis()->SetTitle("Unfolded / Gen");
      rat_obs_combine->GetYaxis()->SetTitleSize(0.1);
      rat_obs_combine->GetYaxis()->SetLabelSize(0.1);
      rat_obs_combine->GetYaxis()->SetTitleOffset(0.75);
      rat_obs_combine->GetYaxis()->SetNdivisions(406);
      rat_obs_combine->GetYaxis()->CenterTitle();
  
      rat_obs_combine->SetMinimum(0);
      rat_obs_combine->SetMaximum(2.);
      
      rat_obs_combine->SetFillStyle(0);
      rat_obs_combine->SetFillColor(0);
      rat_obs_combine->SetMarkerStyle(kFullCircle);
      rat_obs_combine->SetMarkerColor(kBlack);
      rat_obs_combine->SetMarkerSize(0.7);
      rat_obs_combine->SetLineColor(kBlack);
      
      rat_obs_tunfold->SetFillStyle(0);
      rat_obs_tunfold->SetFillColor(0);
      rat_obs_tunfold->SetMarkerStyle(kFullCircle);
      rat_obs_tunfold->SetMarkerColor(kRed);
      rat_obs_tunfold->SetMarkerSize(0.7);
      rat_obs_tunfold->SetLineColor(kRed);
      
      rat_obs_combine->Draw("e,p");
      rat_obs_tunfold->Draw("esame");
      rat_obs_combine->Draw("e,psame");
  
      TLegend *legend_ratio = new TLegend(0.77,0.74,0.93,0.87);
      legend_ratio->SetBorderSize(0);
      legend_ratio->SetTextFont(42);
      legend_ratio->SetTextSize(0.045);
      legend_ratio->AddEntry(rat_obs_combine, "Combine","l");
      legend_ratio->AddEntry(rat_obs_tunfold, "Tunfold","l");
      legend_ratio->Draw("same");
      
      TLine *line = new TLine(rat_obs_combine->GetXaxis()->GetXmin(),1,rat_obs_combine->GetXaxis()->GetXmax(),1);
      line->SetLineColor(1);
      line->Draw("same");
      
      if(ih==0) c1->SaveAs(TString::Format("plots/withoutfirstbin_responsematrix_unfolded_distribution_%s_higgscombine.png",genvar.c_str()));
      else c1->SaveAs(TString::Format("plots/withoutfirstbin_responsematrix_unfolded_distribution_%s_set%d_higgscombine.png",genvar.c_str(),ih));

    }
  for(int ih=0; ih < 11; ih++)
    {	   
      hist_unfold[ih]->SetFillStyle(0);
      hist_unfold[ih]->SetFillColor(0);
      hist_unfold[ih]->SetMarkerStyle(kFullCircle);
      hist_unfold[ih]->SetMarkerColor(kBlack);
      hist_unfold[ih]->SetMarkerSize(1.2);
      hist_unfold[ih]->SetLineColor(kBlack);
      hist_unfold[ih]->SetLineWidth(2);
      hist_unfold[ih]->SetLineStyle(1);
      hist_unfold[ih]->SetLineColor(ih);
      
      hist_unfold[ih]->GetXaxis()->SetTitle(gentitle.c_str());
     
      //scale_bin_by_binwidth(hist_unfold[ih],1);
    }
  
  TCanvas *c1 = tdrCanvas("c1",hist_unfold[0],8,0);;
  c1->cd();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gStyle->SetPaintTextFormat( "1.2f" );

  TLegend *legend = new TLegend(0.3,0.6,0.7,0.85);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
  legend->AddEntry(hist_unfold[0], "Unfolded  distribution", "pe2");
  for(int ih=1; ih < 11; ih++)
    legend->AddEntry(hist_unfold[ih], ("Unfolded  distribution using response matrix set " + to_string(ih)).c_str(), "l");
  
  //hGen->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen->GetBinContent(hGen->GetMaximumBin()),hReco->GetBinContent(hReco->GetMaximumBin())),hist_unfold->GetBinContent(hist_unfold->GetMaximumBin())));

  hist_unfold[0]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*hist_unfold[0]->GetBinContent(hist_unfold[0]->GetMaximumBin()));

  //hGen->Draw("e1");
  //hReco->Draw("e1same");
  //hist_unfold_tunfold->Draw("histsame");
  hist_unfold[0]->DrawCopy("Psame");
  hist_unfold[0]->SetFillColor(kBlack);
  hist_unfold[0]->SetFillStyle(3013);
  hist_unfold[0]->Draw("e2same");
  for(int ih=1; ih < 11; ih++)
    hist_unfold[ih]->Draw("histsame");
  //hist_unfold->Draw("Psame");
  //hist_err->Draw("Psame");
  legend->Draw("same");

  CMS_lumi( c1, 8, 0 );
  gPad->Modified();
  c1->SaveAs(TString::Format("plots/jackknife_responsematrix_unfolded_distribution_comparsion_%s_higgscombine.png",genvar.c_str()));


  TH1D *hist_unfold_jacknife_unc = (TH1D*)hist_unfold[0]->Clone();
  for(int ib =0; ib < hist_unfold_jacknife_unc->GetNbinsX(); ib++)
    {
      std::vector<double> bincontents;
      for(int ih =1; ih<11; ih++)
	bincontents.push_back(hist_unfold[ih]->GetBinContent(ib+1));
      double std_dev = TMath::StdDev(bincontents.begin(), bincontents.end() ); 

      hist_unfold_jacknife_unc->SetBinError(ib+1,sqrt(10/9)*std_dev);
      cout<<"Bin "<<ib<<": Content = "<<hist_unfold_jacknife_unc->GetBinContent(ib+1)<<" with combine uncertainty = "<<hist_unfold[0]->GetBinError(ib+1)<<" and jackknife uncertainty = "<<hist_unfold_jacknife_unc->GetBinError(ib+1)<<endl;
    }

  TCanvas *c2 = tdrCanvas("c2",hist_unfold[0],8,0);;
  c2->cd();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gStyle->SetPaintTextFormat( "1.2f" );

  TLegend *legend_2 = new TLegend(0.3,0.6,0.7,0.85);
  legend_2->SetBorderSize(0);
  legend_2->SetTextFont(42);
  legend_2->SetTextSize(0.03);
  legend_2->AddEntry(hist_unfold[0], "Unfolded  distribution with combine uncertainity", "pe2");
  legend_2->AddEntry(hist_unfold_jacknife_unc, "Unfolded  distribution with jackknife response matrix uncertainity", "l");
   
  //hGen->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen->GetBinContent(hGen->GetMaximumBin()),hReco->GetBinContent(hReco->GetMaximumBin())),hist_unfold->GetBinContent(hist_unfold->GetMaximumBin())));

  hist_unfold[0]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*hist_unfold[0]->GetBinContent(hist_unfold[0]->GetMaximumBin()));

  //hGen->Draw("e1");
  //hReco->Draw("e1same");
  //hist_unfold_tunfold->Draw("histsame");
  hist_unfold[0]->DrawCopy("Psame");
  hist_unfold[0]->SetFillColor(kBlack);
  hist_unfold[0]->SetFillStyle(3013);
  hist_unfold[0]->Draw("e2same");
  hist_unfold_jacknife_unc->SetLineColor(kRed);
  hist_unfold_jacknife_unc->Draw("e1same");
  legend_2->Draw("same");

  CMS_lumi( c2, 8, 0 );
  gPad->Modified();
  c2->SaveAs(TString::Format("plots/jackknife_compare_responsematrixuncertainty_unfolded_distribution_%s_combinevsjackknife.png",genvar.c_str()));

}


void get_bin_migration_plots(TH2D *hResp,string genvar,string gentitle)
{
  int nrebins[2] = {1,2};
  //int nrebins[1] = {6};
  float min_events = 1;
  
  TFile* outFile = TFile::Open("MuMu_Summer20UL18_unfolding_input.root", "UPDATE");

  for(int ib = 0; ib < 2; ib++)
    {
      cout<<ib<<" bins = "<<hResp->GetNbinsX()<<" "<<hResp->GetNbinsY()<<endl;
      TH2D *hist_resp = (TH2D*)hResp->Clone();
      hist_resp->Rebin2D(nrebins[ib], nrebins[ib]);
      cout<<nrebins[ib]<<" bins = "<<hist_resp->GetNbinsX()<<" "<<hist_resp->GetNbinsY()<<endl;
      
      TH1D *hist_stability, *hist_fake_rate, *hist_totalreco;
      TH1D *hist_purity, *hist_efficiency, *hist_totalgen;
  
      hist_stability = (TH1D*)hist_resp->ProjectionX()->Clone();
      hist_fake_rate = (TH1D*)hist_resp->ProjectionX()->Clone();
      hist_totalreco = (TH1D*)hist_resp->ProjectionX()->Clone();

      hist_purity = (TH1D*)hist_resp->ProjectionY()->Clone();
      hist_totalgen = (TH1D*)hist_resp->ProjectionY()->Clone();
      hist_efficiency = (TH1D*)hist_resp->ProjectionY()->Clone();


      int nbinx = hist_resp->GetNbinsX();
      int nbiny = hist_resp->GetNbinsY();

      hist_totalgen->SetName(TString::Format("hist_%s_ngen_nbin%i",genvar.c_str(),int(nbinx)));
      hist_purity->SetName(TString::Format("hist_%s_npurity_nbin%i",genvar.c_str(),int(nbinx)));
      hist_efficiency->SetName(TString::Format("hist_%s_nefficiency_nbin%i",genvar.c_str(),int(nbinx)));
      
      hist_totalreco->SetName(TString::Format("hist_%s_nreco_nbin%i",genvar.c_str(),int(nbinx)));
      hist_stability->SetName(TString::Format("hist_%s_nstability_nbin%i",genvar.c_str(),int(nbinx)));
      hist_fake_rate->SetName(TString::Format("hist_%s_nfake_nbin%i",genvar.c_str(),int(nbinx)));

      if(nbiny != nbinx)
	{
	  cout<<endl<<"Please Check bining of response matrix histogram. The purity, stability calculation might not be straight forward."<<endl;
	  exit(0);
	}
      
      double totalreco[200];
      double stability[200];
      double fakerate[200];
      
      double totalgen[200];
      double effi[200];
      double purity[200];
      
      for(int ib=0; ib< 200; ib++)
	{
	  totalreco[ib] = 0;
	  fakerate[ib] = 0;
	  stability[ib] = 0;
	}   
      
      for(int ib=0; ib< 200; ib++)
	{
	  totalgen[ib] = 0;
	  effi[ib] = 0;
	  purity[ib] = 0;
	}
      
      
      for (int ix=0; ix<nbinx+1; ix++) {
	for (int iy=0; iy<nbiny+1; iy++) {

	  if(ix==0&&iy==0) continue ;
	  totalreco[ix] +=hist_resp->GetBinContent(ix, iy);
	  if (iy==0) fakerate[ix] = hist_resp->GetBinContent(ix,iy);
	  
	  totalgen[iy] +=hist_resp->GetBinContent(ix, iy);
	  if (ix==0) effi[iy] = hist_resp->GetBinContent(ix, iy);
	  if (ix == iy)
 	    {
	      purity[iy] = hist_resp->GetBinContent(iy, iy);
	      stability[ix] = hist_resp->GetBinContent(ix, ix);
	    }	
	}//iy
      }//ix
      
      for(int ix=1; ix <nbinx +1; ix++){
	//cout<<totalreco[ix]<<" "<<stability[ix]<<endl;
	hist_totalreco->SetBinContent(ix, totalreco[ix]);
	hist_totalreco->SetBinError(ix, sqrt(totalreco[ix]));
	
	hist_stability->SetBinContent(ix, stability[ix]);
	hist_stability->SetBinError(ix, sqrt(stability[ix]));
	
	hist_fake_rate->SetBinContent(ix, fakerate[ix]);
	hist_fake_rate->SetBinError(ix, sqrt(fakerate[ix]));
	//cout<<hist_totalreco->GetBinContent(ix)<<" "<<hist_stability->GetBinContent(ix)<<endl;

      }
      
      for(int iy=0; iy<nbiny +1; iy++){
	hist_totalgen->SetBinContent(iy, totalgen[iy]);
	hist_totalgen->SetBinError(iy, sqrt(totalgen[iy]));
	
	hist_efficiency->SetBinContent(iy, totalgen[iy] - effi[iy]);
	hist_efficiency->SetBinError(iy, sqrt(totalgen[iy] - effi[iy]));
	
	hist_purity->SetBinContent(iy, purity[iy]);
	hist_purity->SetBinError(iy, sqrt(purity[iy]));
      }
            
      /*      hist_totalgen->Write(TString::Format("hist_%s_ngen_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      hist_purity->Write(TString::Format("hist_%s_npurity_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      hist_efficiency->Write(TString::Format("hist_%s_nefficiency_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      
      hist_totalreco->Write(TString::Format("hist_%s_nreco_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      hist_stability->Write(TString::Format("hist_%s_nstability_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      hist_fake_rate->Write(TString::Format("hist_%s_nfake_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      */
      //cout<<" ib = "<<ib<<endl;
      TH1D *hist_stability_rebinned = (TH1D*)hist_stability->Clone();
      TH1D *hist_fake_rate_rebinned = (TH1D*)hist_fake_rate->Clone();
      TH1D *hist_totalreco_rebinned = (TH1D*)hist_totalreco->Clone();
      
      TH1D *hist_purity_rebinned = (TH1D*)hist_purity->Clone();
      TH1D *hist_totalgen_rebinned = (TH1D*)hist_totalgen->Clone();
      TH1D *hist_efficiency_rebinned = (TH1D*)hist_efficiency->Clone();

      hist_purity_rebinned->Divide(hist_totalgen_rebinned);
      hist_efficiency_rebinned->Divide(hist_totalgen_rebinned);
      hist_stability_rebinned->Divide(hist_totalreco_rebinned);
      hist_fake_rate_rebinned->Divide(hist_totalreco_rebinned);
      
      /*for(int iy=1; iy < hist_purity_rebinned->GetNbinsX() + 1; iy++){	
	//hist_efficiency_rebinned->SetBinContent(iy, (totalgen[iy] - effi[iy]) / totalgen[iy]);
	hist_purity_rebinned->SetBinContent(iy, purity[iy] / max(0.00000001,(totalgen[iy] - effi[iy])) );
	//hist_purity_rebinned->SetBinError(iy, 0);
	cout<<" ib = "<<iy<<" stability = "<<purity[iy] / max(0.00000001,(totalgen[iy] - effi[iy]))<<" total gen = "<<totalgen[iy]<<" ineffici = "<<effi[iy]<<endl;
      }
      
      for(int ix=1; ix < hist_stability_rebinned->GetNbinsX() +1; ix++){
	hist_stability_rebinned->SetBinContent(ix, stability[ix] / max(0.00000001,(totalreco[ix]-fakerate[ix])) );
	hist_stability_rebinned->SetBinError(ix, 0);
	}
      */
      /*hist_purity_rebinned->;
      hist_efficiency_rebinned->;
      hist_stability_rebinned->;
      hist_fake_rate_rebinned->;
      */
      
      hist_purity_rebinned->SetMaximum(1.3);
      hist_efficiency_rebinned->SetMaximum(1.3);
      hist_stability_rebinned->SetMaximum(1.3);
      hist_fake_rate_rebinned->SetMaximum(1.3);

      hist_purity_rebinned->SetMinimum(0.0);
      hist_efficiency_rebinned->SetMinimum(0.0);
      hist_stability_rebinned->SetMinimum(0.0);
      hist_fake_rate_rebinned->SetMinimum(0.0);
      
      /*hist_purity_rebinned->Write(TString::Format("hist_%s_purity_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      hist_efficiency_rebinned->Write(TString::Format("hist_%s_efficiency_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      hist_stability_rebinned->Write(TString::Format("hist_%s_stability_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      hist_fake_rate_rebinned->Write(TString::Format("hist_%s_fake_rate_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      */
      
      hist_purity_rebinned->SetLineWidth(2);
      hist_efficiency_rebinned->SetLineWidth(2);
      hist_stability_rebinned->SetLineWidth(2);
      hist_fake_rate_rebinned->SetLineWidth(2);
      
      hist_purity_rebinned->SetLineColor(kRed);
      hist_efficiency_rebinned->SetLineColor(kBlue);
      hist_stability_rebinned->SetLineColor(kYellow);
      hist_fake_rate_rebinned->SetLineColor(kMagenta);
      
      hist_purity_rebinned->GetXaxis()->SetTitle(gentitle.c_str());
      hist_purity_rebinned->GetYaxis()->SetTitle("Fraction of events");
      hist_purity_rebinned->GetYaxis()->SetTitleSize(0.045);
      hist_purity_rebinned->GetYaxis()->SetTitleOffset(1.2);
      hist_purity_rebinned->GetYaxis()->SetLabelSize(0.035);

      hist_efficiency_rebinned->GetXaxis()->SetTitle(gentitle.c_str());
      hist_efficiency_rebinned->GetYaxis()->SetTitle("Fraction of events");
      hist_efficiency_rebinned->GetYaxis()->SetTitleSize(0.045);
      hist_efficiency_rebinned->GetYaxis()->SetTitleOffset(1.2);
      hist_efficiency_rebinned->GetYaxis()->SetLabelSize(0.035);

      TCanvas *c1 = tdrCanvas("canv_1",hist_stability_rebinned,8,0);
      c1->cd();
      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      
      TLegend *legend = new TLegend(0.54,0.7,0.90,0.88);
      legend->SetBorderSize(0);
      legend->SetTextFont(42);
      legend->SetTextSize(0.035);
      legend->AddEntry(hist_efficiency_rebinned, "Efficiency");
      legend->AddEntry(hist_fake_rate_rebinned, "Fake rate");     

      /*hist_purity_rebinned->Draw("histsame");
      hist_efficiency_rebinned->Draw("histsame");
      hist_stability_rebinned->Draw("histsame");
      hist_fake_rate_rebinned->Draw("histsame");
      */

      hist_efficiency_rebinned->Draw("hist");
      hist_fake_rate_rebinned->Draw("histsame");
      CMS_lumi( c1, 8, 0 );

      legend->Draw("same");
	
      //c1->Write(TString::Format("plots/all_plots_%s_nbin%i",genvar.c_str(),int(nbinx)),TObject::kOverwrite);
      c1->SaveAs(TString::Format("plots/migrationplot_eff_fakerate_%s_nbin%i.png",genvar.c_str(),int(nbinx)));


      TCanvas *c2 = tdrCanvas("canv_2",hist_stability_rebinned,8,0);
      c2->cd();
      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      
      TLegend *legend2 = new TLegend(0.54,0.7,0.90,0.88);
      legend2->SetBorderSize(0);
      legend2->SetTextFont(42);
      legend2->SetTextSize(0.035);
      legend2->AddEntry(hist_stability_rebinned, "Purity");
      legend2->AddEntry(hist_purity_rebinned, "Stability");
      
      hist_purity_rebinned->Draw("hist");
      hist_stability_rebinned->Draw("histsame");
      CMS_lumi( c2, 8, 0 );
  
      legend2->Draw("same");
      
      c2->SaveAs(TString::Format("plots/migrationplot_stability_purity_%s_nbin%i.png",genvar.c_str(),int(nbinx)));
    }
  outFile->Close();

}


void get_unfold_hist(string genvar, string gentitle, string recovar) 
{
  TFile* sigFile = TFile::Open("MuMu_Summer20UL18_unfolding_input.root", "READ");
  TH1D* hGen = (TH1D*)sigFile->Get(("hist_gen_" + genvar).c_str());
  hGen->SetDirectory(0);
  
  //TH1D* hReco = (TH1D*)sigFile->Get(("hist_reco_" + recovar).c_str());
  TH1D* hReco = (TH1D*)sigFile->Get(("hist_gen_" + genvar +"_data_obs").c_str());
  hReco->SetDirectory(0);
  
  TH2D* hResp = (TH2D*)sigFile->Get(("hist_responsematrix_recogen_" + genvar).c_str());
  hResp->SetDirectory(0);
  sigFile->Close();

  TFile* combine_output = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar +".root").c_str(), "READ");
  TTree * tree = (TTree*)combine_output->Get("limit");
  static const int nbins = int(hGen->GetNbinsX());
  float r_values[nbins];
  float r_bins[nbins];
  float r_bins_high[nbins];
  float r_bins_low[nbins];
  
  for(int ib = 0; ib < nbins; ib++)
    tree->SetBranchAddress(TString::Format("r_bin%d",ib+1),&r_values[ib]);

  for(int i=0; i <nbins; i++)   
    {
      tree->GetEntry(0);
      r_bins[i] = r_values[i];
      tree->GetEntry(2*i+1); 
      r_bins_low[i] = r_values[i];
      tree->GetEntry(2*i+2); 
      r_bins_high[i] = r_values[i];      
    }
  tree->GetEntry(0);   //////// get best fit values

  TFile* combine_output_withstatonly = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar +"_withstatonly.root").c_str(), "READ");
  TTree * tree_withstatonly = (TTree*)combine_output_withstatonly->Get("limit");
  float r_values_withstatonly[nbins];
  float r_bins_high_withstatonly[nbins];
  float r_bins_low_withstatonly[nbins];
  
  for(int ib = 0; ib < nbins; ib++)
    tree_withstatonly->SetBranchAddress(TString::Format("r_bin%d",ib+1),&r_values_withstatonly[ib]);

  for(int i=0; i <nbins; i++)   
    {
      tree_withstatonly->GetEntry(0);
      //r_bins[i] = r_values[i];
      tree_withstatonly->GetEntry(2*i+1); 
      r_bins_low_withstatonly[i] = r_values_withstatonly[i];
      tree_withstatonly->GetEntry(2*i+2); 
      r_bins_high_withstatonly[i] = r_values_withstatonly[i];

      double rlow_statonly = (r_bins_low[i]*r_bins_low[i] - r_bins_low_withstatonly[i]*r_bins_low_withstatonly[i]);
      double rhigh_statonly = (r_bins_high[i]*r_bins_high[i] - r_bins_high_withstatonly[i]*r_bins_high_withstatonly[i]);
      
      if( rlow_statonly < 0.0000001 )
	{cout<<"Check errors: Stat only error is > than total error "<<r_bins_low[i]<<" - "<<r_bins_low_withstatonly[i]<<" stat error = "<<rlow_statonly<<"for hist:"<<hReco->GetName()<<endl;exit(0);}
      r_bins_low[i] = sqrt(rlow_statonly);

      if(rhigh_statonly  < 0.0000001 )
	{cout<<"Check errors: Stat only error is > than total error "<<r_bins_high[i]<<" - "<<r_bins_high_withstatonly[i]<<" stat error = "<<rhigh_statonly<<"for hist:"<<hReco->GetName()<<endl;exit(0);}
      r_bins_high[i] = sqrt(rhigh_statonly);
    }
  //tree_withstatonly->GetEntry(0);   //////// get best fit values

  TH1D *hist_unfold = (TH1D*)hGen->Clone();
  TH1D *hist_err = (TH1D*)hGen->Clone();
  TH1D *rat_obs_combine;
  TH1D *rat_obs_tunfold;
  
  TH1D *hist_unfold_tunfold;
  hist_unfold_tunfold = (TH1D*)TUnfold_unfold(genvar,gentitle,recovar);
  cout<<"Check1"<<endl;
  double* dataerror = get_jackknife_stat_error_ondata(genvar, gentitle, recovar);
  double* resperror = get_jackknife_stat_error_onresponsematrix(genvar, gentitle, recovar);

  for(int ib=0; ib < nbins; ib++)
    {
      cout<<r_bins[ib]<<" "<<r_bins_low[ib]-r_bins[ib]<<" "<<r_bins_high[ib]-r_bins[ib]<<endl;
      float bincontent = 0;//hResp->Integral(1,hResp->GetNbinsX(),ib+1,ib+1)   //hGen->GetBinContent(ib+1)
      for(int ix = 1; ix < hResp->GetNbinsX()+1; ix++)
	{
	  if(hResp->GetBinContent(ix,ib+1) > 0.00001)
	    bincontent += hResp->GetBinContent(ix,ib+1);
	}
      hist_unfold->SetBinContent(ib+1, bincontent*r_bins[ib]);

      double error = bincontent*(r_bins_high[ib] - r_bins_low[ib])/2;
      if(ib > 0)
	error = sqrt(bincontent*(r_bins_high[ib] - r_bins_low[ib])/2*bincontent*(r_bins_high[ib] - r_bins_low[ib])/2 + dataerror[ib-1]*dataerror[ib-1] + resperror[ib-1]*resperror[ib-1]);
      hist_unfold->SetBinError(ib+1,error);
      
      hist_err->SetBinContent(ib+1, abs(bincontent*r_bins[ib] - hGen->GetBinContent(ib+1))/hGen->GetBinContent(ib+1));
    }

  cout<<endl;
  hist_unfold->SetFillStyle(0);
  hist_unfold->SetFillColor(0);
  hist_unfold->SetMarkerStyle(kFullCircle);
  hist_unfold->SetMarkerColor(kBlack);
  hist_unfold->SetMarkerSize(1.2);
  hist_unfold->SetLineColor(kBlack);
  hist_unfold->SetLineWidth(2);
  hist_unfold->SetLineStyle(1);
  hist_unfold->SetLineColor(0);

  cout<<"Check2"<<endl;
  hist_unfold_tunfold->SetMarkerStyle(kFullStar);
  hist_unfold_tunfold->SetMarkerColor(kBlack);
  hist_unfold_tunfold->SetLineStyle(1);
  hist_unfold_tunfold->SetMarkerSize(3);
  hist_unfold_tunfold->SetLineWidth(2);
  hist_unfold_tunfold->SetLineColor(kRed);
  cout<<"Check3"<<endl;

  hGen->GetXaxis()->SetTitle(gentitle.c_str());
  hGen->SetLineWidth(3);
  //hGen->SetLineStyle(2);
  hGen->SetLineColor(kBlue);

  hReco->SetLineWidth(3);
  //hReco->SetLineStyle(2);
  hReco->SetLineColor(kGreen);
  hReco->Rebin(2);

  scale_bin_by_binwidth(hist_unfold,1);
  scale_bin_by_binwidth(hist_unfold_tunfold,1);
  scale_bin_by_binwidth(hGen,1);
  scale_bin_by_binwidth(hReco,1);
  
  //hGen->SetMaximum(1.9*max(max(hGen->GetMaximum(),hReco->GetMaximum()),hist_unfold->GetMaximum()));

  hist_err->SetFillStyle(0);
  //hist_err->SetMarkerColor(0);
  hist_err->SetLineColor(0);
  hist_err->SetMarkerSize(1.2);

  TCanvas *c1 = tdrDiCanvas("c1",hGen,hReco,8,0);;
  c1->cd(1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gStyle->SetPaintTextFormat( "1.2f" );

  TLegend *legend = new TLegend(0.34,0.6,0.74,0.85);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
  legend->AddEntry(hist_unfold, "Unfolded  distribution using Combine", "pe2");
  legend->AddEntry(hist_unfold_tunfold, "Unfolded  distribution using Tunfold", "l");
  legend->AddEntry(hGen, "Gen Madgraph+Pythia sample", "l");
  legend->AddEntry(hReco, "Reco Madgraph+Pythia sample", "l");

  hGen->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen->GetBinContent(hGen->GetMaximumBin()),hReco->GetBinContent(hReco->GetMaximumBin())),hist_unfold->GetBinContent(hist_unfold->GetMaximumBin())));

  hGen->Draw("e1");
  hReco->Draw("e1same");
  hist_unfold_tunfold->Draw("histsame");
  hist_unfold->DrawCopy("Psame");
  hist_unfold->SetFillColor(kBlack);
  hist_unfold->SetFillStyle(3013);
  hist_unfold->Draw("e2same");
  //hist_unfold->Draw("Psame");
  //hist_err->Draw("Psame");
  legend->Draw("same");

  CMS_lumi( c1, 8, 0 );
  gPad->Modified();


  c1->cd(2);
  rat_obs_combine = (TH1D*)hist_unfold->Clone();
  rat_obs_combine->Divide(hGen);

  rat_obs_tunfold = (TH1D*)hist_unfold_tunfold->Clone();
  rat_obs_tunfold->Divide(hGen);
  rat_obs_tunfold->SetLineColor(kRed);

  rat_obs_combine->GetXaxis()->SetTitle(gentitle.c_str());
  rat_obs_combine->GetXaxis()->SetTitleSize(0.09);
  rat_obs_combine->GetXaxis()->SetLabelSize(0.1);
  rat_obs_combine->GetXaxis()->CenterTitle();
  rat_obs_combine->GetXaxis()->SetNdivisions(406);

  rat_obs_combine->GetYaxis()->SetTitle("Unfolded / Gen");
  rat_obs_combine->GetYaxis()->SetTitleSize(0.1);
  rat_obs_combine->GetYaxis()->SetLabelSize(0.1);
  rat_obs_combine->GetYaxis()->SetTitleOffset(0.75);
  rat_obs_combine->GetYaxis()->SetNdivisions(406);
  rat_obs_combine->GetYaxis()->CenterTitle();
  
  rat_obs_combine->SetMinimum(0);
  rat_obs_combine->SetMaximum(2.);
  
  rat_obs_combine->SetFillStyle(0);
  rat_obs_combine->SetFillColor(0);
  rat_obs_combine->SetMarkerStyle(kFullCircle);
  rat_obs_combine->SetMarkerColor(kBlack);
  rat_obs_combine->SetMarkerSize(0.7);
  rat_obs_combine->SetLineColor(kBlack);

  rat_obs_tunfold->SetFillStyle(0);
  rat_obs_tunfold->SetFillColor(0);
  rat_obs_tunfold->SetMarkerStyle(kFullCircle);
  rat_obs_tunfold->SetMarkerColor(kRed);
  rat_obs_tunfold->SetMarkerSize(0.7);
  rat_obs_tunfold->SetLineColor(kRed);

  rat_obs_combine->Draw("e,p");
  rat_obs_tunfold->Draw("esame");
  rat_obs_combine->Draw("e,psame");
  
  TLegend *legend_ratio = new TLegend(0.77,0.74,0.93,0.87);
  legend_ratio->SetBorderSize(0);
  legend_ratio->SetTextFont(42);
  legend_ratio->SetTextSize(0.045);
  legend_ratio->AddEntry(rat_obs_combine, "Combine","l");
  legend_ratio->AddEntry(rat_obs_tunfold, "Tunfold","l");
  legend_ratio->Draw("same");
    
  TLine *line = new TLine(rat_obs_combine->GetXaxis()->GetXmin(),1,rat_obs_combine->GetXaxis()->GetXmax(),1);
  line->SetLineColor(1);
  line->Draw("same");
  
  c1->SaveAs(TString::Format("plots/unfolded_distribution_%s_higgscombine.png",genvar.c_str()));
}

void make_input_histforunfolding(TH1D *hReco, TH1D* hGen, TH2D* hResp, TH1D *hReco_Data = NULL, int nbins = 0,int nsys = 0, TH1D *hSign_bin[6][11] = NULL)
{

  int nbinx = hResp->GetNbinsX();
  int nbiny = hResp->GetNbinsY();
  
  static const int nbinmx = 200; //max(hResp->GetNbinsY(),hResp->GetNbinsX()) + 5;
  double totalgen[nbinmx]={0.};
  double totalreco[nbinmx]={0.};

  for(int ib=0; ib< nbinmx; ib++)
    {
      totalgen[ib] = 0;
      totalreco[ib] = 0;
    }

  double fakerate[nbinmx];
  double effi[nbinmx];
  
  for (int ix=0; ix<nbinx+1; ix++) {
    for (int iy=0; iy<nbiny+1; iy++) {
      
      if(ix==0&&iy==0) continue ;
      if(hResp->GetBinContent(ix,iy) < -0.1)
	hResp->SetBinContent(ix,iy,0);
      
      totalreco[ix] +=hResp->GetBinContent(ix, iy);
      if (iy==0 ) fakerate[ix] = hResp->GetBinContent(ix,iy);
	
      totalgen[iy] +=hResp->GetBinContent(ix, iy);
      if (ix==0) effi[iy] = hResp->GetBinContent(ix, iy);

      if (ix==0 || iy==0) {
	hResp->SetBinContent(ix, iy, 0.0);
	hResp->SetBinError(ix, iy, 0.0);
	}
    }//iy
  }//ix
  
  for (int iy=1; iy<nbiny+1; iy++)
    {
      effi[iy] = (totalgen[iy] - effi[iy])/max(1.e-10, totalgen[iy]);
      if (abs(effi[iy]) > 1) effi[iy] = 1;
      else if ( effi[iy] < 0) effi[iy] = 0.0000001;
    } //iy

  for (int ix=1; ix<nbinx+1; ix++)
    {
      fakerate[ix] = fakerate[ix] / max(1.e-10, totalreco[ix]);
      if(abs(fakerate[ix]) > 1) fakerate[ix] = 0.99999999;
      else if ( fakerate[ix] < 0) fakerate[ix] = 0.0000001;
    }//ix

  if( nbinx != hReco->GetNbinsX()  &&  nbiny != hGen->GetNbinsX() )
  {
      cout<<"Check bining of histograms  Nbins Gen hist("<<hGen->GetNbinsX()<<") should be = NbinsY Response matrix ("<<nbiny<<") and Nbins Reco hist("<<hReco->GetNbinsX()<<") should be = NbinsX Response matrix("<<nbinx<<")"<<endl;
      exit(0);
  }
 
  
  for(int ix=0; ix <((hReco->GetNbinsX())+1); ix++){
    cout<<"reco "<<ix+1<<" = "<<1-fakerate[ix+1]<<" * "<<hReco->GetBinContent(ix+1)<<endl;
    hReco->SetBinContent(ix+1,(1-fakerate[ix+1])*(hReco->GetBinContent(ix+1)));
    hReco->SetBinError(ix+1,(1-fakerate[ix+1])*(hReco->GetBinError(ix+1)));
  }

  for(int ix=0; ix<(hGen->GetNbinsX()+1); ix++){
    cout<<"gen "<<ix+1<<" = "<<effi[ix+1]<<" * "<<hGen->GetBinContent(ix+1)<<endl;
    hGen->SetBinContent(ix+1,(hGen->GetBinContent(ix+1))*effi[ix+1]);
    hGen->SetBinError(ix+1,(hGen->GetBinError(ix+1))*effi[ix+1]);
  }

  cout<<"Check input histograms. They dont have equal yields  Gen hist = "<<hGen->Integral()<<" Reco hist = "<<hReco->Integral()<<" Response matrix = "<<hResp->Integral()<<endl;

  if( abs(hResp->Integral() - hGen->Integral()) > 0.1 &&  abs(hGen->Integral() - hReco->Integral()) > 0.1)
    {
      cout<<"Check input histograms. They dont have equal yields  Gen hist = "<<hGen->Integral()<<" Reco hist = "<<hReco->Integral()<<" Response matrix = "<<hResp->Integral()<<endl;
      //exit(0);
    }
  
  TFile* outFile = TFile::Open("MuMu_Summer20UL18_unfolding_input.root", "UPDATE");

  outFile->cd();
  hGen->Write(hGen->GetName(),TObject::kOverwrite);
  hReco->Write(hReco->GetName(),TObject::kOverwrite);
  hResp->Write(hResp->GetName(),TObject::kOverwrite);
  if(hReco_Data != NULL)
    hReco_Data->Write(hReco_Data->GetName(),TObject::kOverwrite);
  for(int ib =0 ; ib < nbins; ib ++)
    {
      for(int isy =0 ; isy < nsys; isy ++)
	{
	  if(hSign_bin[ib][isy] != NULL)
	    {
	      hSign_bin[ib][isy]->Scale(effi[ib+1]);
	      hSign_bin[ib][isy]->Write(hSign_bin[ib][isy]->GetName(),TObject::kOverwrite);
	    }
	}
    }
  
  outFile->Close();
}

void get_unfold_hist_withoutfirstbin(string genvar, string gentitle, string recovar) 
{
  TFile* sigFile = TFile::Open("MuMu_Summer20UL18_unfolding_input.root", "READ");
  TH1D* hGen_allbins = (TH1D*)sigFile->Get(("hist_gen_" + genvar).c_str());
  hGen_allbins->SetDirectory(0);
  
  TH1D* hReco_allbins = (TH1D*)sigFile->Get(("hist_reco_" + recovar).c_str());
  hReco_allbins->SetDirectory(0);
  
  TH2D* hResp_allbins = (TH2D*)sigFile->Get(("hist_responsematrix_recogen_" + genvar).c_str());
  hResp_allbins->SetDirectory(0);
  sigFile->Close();

  TH1D* hGen = (TH1D*)removefirstbin(hGen_allbins);
  
  TH1D* hReco = (TH1D*)removefirstbin(hReco_allbins,2);
    
  TFile* combine_output = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar +".root").c_str(), "READ");
  TTree * tree = (TTree*)combine_output->Get("limit");
  static const int nbins = int(hGen->GetNbinsX());
  float r_values[nbins];
  float r_bins[nbins];
  float r_bins_high[nbins];
  float r_bins_low[nbins];
  
  for(int ib = 0; ib < nbins; ib++)
    tree->SetBranchAddress(TString::Format("r_bin%d",ib+2),&r_values[ib]);

  for(int i=0; i <nbins; i++)   
    {
      tree->GetEntry(0);
      r_bins[i] = r_values[i];
      tree->GetEntry(2*i+3); 
      r_bins_low[i] = r_values[i];
      tree->GetEntry(2*i+4); 
      r_bins_high[i] = r_values[i];      
    }
  tree->GetEntry(0);   //////// get best fit values

  TH1D *hist_unfold = (TH1D*)hGen->Clone();
  TH1D *hist_err = (TH1D*)hGen->Clone();
  TH1D *rat_obs_combine;
  TH1D *rat_obs_tunfold;
  
  TH1D *hist_unfold_tunfold_allbins;
  hist_unfold_tunfold_allbins = (TH1D*)TUnfold_unfold(genvar,gentitle,recovar);
  TH1D *hist_unfold_tunfold = (TH1D*)removefirstbin(hist_unfold_tunfold_allbins);

  cout<<"Check1"<<endl;
  for(int ib=0; ib < nbins; ib++)
    {
      cout<<r_bins[ib]<<" "<<r_bins_low[ib]-r_bins[ib]<<" "<<r_bins_high[ib]-r_bins[ib]<<endl;
      cout<<"gen bin"<<ib<<" content"<<hGen->GetBinContent(ib+1)<<" "<<hGen->GetBinLowEdge(ib+1)<<endl;
      float bincontent = 0;//hResp_allbins->Integral(1,hResp_allbins->GetNbinsX(),ib+1,ib+1)   //hGen->GetBinContent(ib+1)
      for(int ix = 1; ix < hResp_allbins->GetNbinsX()+1; ix++)
	{
	  if(hResp_allbins->GetBinContent(ix,ib+2) > 0.00001)
	    bincontent += hResp_allbins->GetBinContent(ix,ib+2);
	}
      hist_unfold->SetBinContent(ib+1, bincontent*r_bins[ib]);
      hist_unfold->SetBinError(ib+1, bincontent*(r_bins_high[ib] - r_bins_low[ib])/2);
      
      hist_err->SetBinContent(ib+1, abs(bincontent*r_bins[ib] - hGen->GetBinContent(ib+1))/hGen->GetBinContent(ib+1));
    }

  cout<<endl;
  hist_unfold->SetFillStyle(0);
  hist_unfold->SetFillColor(0);
  hist_unfold->SetMarkerStyle(kFullCircle);
  hist_unfold->SetMarkerColor(kBlack);
  hist_unfold->SetMarkerSize(1.2);
  hist_unfold->SetLineColor(kBlack);
  hist_unfold->SetLineWidth(2);
  hist_unfold->SetLineStyle(1);
  hist_unfold->SetLineColor(0);

  cout<<"Check2"<<endl;
  hist_unfold_tunfold->SetMarkerStyle(kFullStar);
  hist_unfold_tunfold->SetMarkerColor(kBlack);
  hist_unfold_tunfold->SetLineStyle(1);
  hist_unfold_tunfold->SetMarkerSize(3);
  hist_unfold_tunfold->SetLineWidth(2);
  hist_unfold_tunfold->SetLineColor(kRed);
  cout<<"Check3"<<endl;

  hGen->GetXaxis()->SetTitle(gentitle.c_str());
  hGen->SetLineWidth(3);
  //hGen->SetLineStyle(2);
  hGen->SetLineColor(kBlue);

  hReco->SetLineWidth(3);
  //hReco->SetLineStyle(2);
  hReco->SetLineColor(kGreen);
  hReco->Rebin(2);

  scale_bin_by_binwidth(hist_unfold,1);
  scale_bin_by_binwidth(hist_unfold_tunfold,1);
  scale_bin_by_binwidth(hGen,1);
  scale_bin_by_binwidth(hReco,1);
  
  
  //hGen->SetMaximum(1.9*max(max(hGen->GetMaximum(),hReco->GetMaximum()),hist_unfold->GetMaximum()));

  hist_err->SetFillStyle(0);
  //hist_err->SetMarkerColor(0);
  hist_err->SetLineColor(0);
  hist_err->SetMarkerSize(1.2);

  TCanvas *c1 = tdrDiCanvas("c1",hGen,hReco,8,0);;
  c1->cd(1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gStyle->SetPaintTextFormat( "1.2f" );

  TLegend *legend = new TLegend(0.34,0.6,0.74,0.85);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
  legend->AddEntry(hist_unfold, "Unfolded  distribution using Combine", "pe2");
  legend->AddEntry(hist_unfold_tunfold, "Unfolded  distribution using Tunfold", "l");
  legend->AddEntry(hGen, "Gen Madgraph+Pythia sample", "l");
  legend->AddEntry(hReco, "Reco Madgraph+Pythia sample", "l");

  hGen->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen->GetBinContent(hGen->GetMaximumBin()),hReco->GetBinContent(hReco->GetMaximumBin())),hist_unfold->GetBinContent(hist_unfold->GetMaximumBin())));
  
  hGen->Draw("e1");
  hReco->Draw("e1same");
  hist_unfold_tunfold->Draw("histsame");
  hist_unfold->DrawCopy("Psame");
  hist_unfold->SetFillColor(kBlack);
  hist_unfold->SetFillStyle(3013);
  hist_unfold->Draw("e2same");
  //hist_unfold->Draw("Psame");
  //hist_err->Draw("Psame");
  legend->Draw("same");

  CMS_lumi( c1, 8, 0 );
  gPad->Modified();


  c1->cd(2);
  
  rat_obs_combine = (TH1D*)hist_unfold->Clone();
  rat_obs_combine->Divide(hGen);

  rat_obs_tunfold = (TH1D*)hist_unfold_tunfold->Clone();
  rat_obs_tunfold->Divide(hGen);
  rat_obs_tunfold->SetLineColor(kRed);

  rat_obs_combine->GetXaxis()->SetTitle(gentitle.c_str());
  rat_obs_combine->GetXaxis()->SetTitleSize(0.09);
  rat_obs_combine->GetXaxis()->SetLabelSize(0.1);
  rat_obs_combine->GetXaxis()->CenterTitle();
  rat_obs_combine->GetXaxis()->SetNdivisions(406);

  rat_obs_combine->GetYaxis()->SetTitle("Unfolded / Gen");
  rat_obs_combine->GetYaxis()->SetTitleSize(0.1);
  rat_obs_combine->GetYaxis()->SetLabelSize(0.1);
  rat_obs_combine->GetYaxis()->SetTitleOffset(0.75);
  rat_obs_combine->GetYaxis()->SetNdivisions(406);
  rat_obs_combine->GetYaxis()->CenterTitle();
  
  rat_obs_combine->SetMinimum(0);
  rat_obs_combine->SetMaximum(2.);
  
  rat_obs_combine->SetFillStyle(0);
  rat_obs_combine->SetFillColor(0);
  rat_obs_combine->SetMarkerStyle(kFullCircle);
  rat_obs_combine->SetMarkerColor(kBlack);
  rat_obs_combine->SetMarkerSize(0.7);
  rat_obs_combine->SetLineColor(kBlack);

  rat_obs_tunfold->SetFillStyle(0);
  rat_obs_tunfold->SetFillColor(0);
  rat_obs_tunfold->SetMarkerStyle(kFullCircle);
  rat_obs_tunfold->SetMarkerColor(kRed);
  rat_obs_tunfold->SetMarkerSize(0.7);
  rat_obs_tunfold->SetLineColor(kRed);

  rat_obs_combine->Draw("e,p");
  rat_obs_tunfold->Draw("esame");
  rat_obs_combine->Draw("e,psame");
  
  TLegend *legend_ratio = new TLegend(0.77,0.74,0.93,0.87);
  legend_ratio->SetBorderSize(0);
  legend_ratio->SetTextFont(42);
  legend_ratio->SetTextSize(0.045);
  legend_ratio->AddEntry(rat_obs_combine, "Combine","l");
  legend_ratio->AddEntry(rat_obs_tunfold, "Tunfold","l");
  legend_ratio->Draw("same");
    
  TLine *line = new TLine(rat_obs_combine->GetXaxis()->GetXmin(),1,rat_obs_combine->GetXaxis()->GetXmax(),1);
  line->SetLineColor(1);
  line->Draw("same");
  
  c1->SaveAs(TString::Format("plots/withoutfirstbin_unfolded_distribution_%s_higgscombine.png",genvar.c_str()));
}

void all_unfolding_prep(string genvar, string gentitle, string mcfile, float ttbar_signal_wgt, string recovar, string recotitle, string datafile, float data_wgt)
{
  TFile* sigFile = TFile::Open(mcfile.c_str(), "READ");
  TH1D* hGen = (TH1D*)sigFile->Get(("hist_gen_" + genvar).c_str());
  hGen->SetDirectory(0);

  TH1D* hReco = (TH1D*)sigFile->Get(("hist_reco_" + recovar).c_str());
  hReco->SetDirectory(0);

  TH2D* hResp = (TH2D*)sigFile->Get(("hist_responsematrix_recogen_" + genvar).c_str());
  //TH2D* hResp = (TH2D*)sigFile->Get("hist_resolutionmatrix_recogen_ak81pt");
  hResp->SetDirectory(0);

  string sysnames[11] = {"","pileupUp","pileupDown","top_tagger_SFUp","top_tagger_SFDown","trigger_SFUp","trigger_SFDown","JERUp","JERDown","JESUp","JESDown"};
  TH1D* hSign_bin[6][11];
    
  for(int igen = 0; igen < 6; igen++)
    {
      for(int isys= 1; isys < 11; isys++)
	{
	  string newtitle = "hist_gen_" + genvar + "_genbin" + to_string(igen +1)+ "_sys_" + sysnames[isys];
	  hSign_bin[igen][isys] = (TH1D*)sigFile->Get(newtitle.c_str());
	  newtitle = "hist_gen_" + genvar + "_genbin" + to_string(igen +1)+ "_sys_" + sysnames[isys];
	  hSign_bin[igen][isys]->SetTitle(newtitle.c_str());
	  hSign_bin[igen][isys]->SetName(newtitle.c_str());
	  hSign_bin[igen][isys]->SetDirectory(0);
	}
      string newtitle = "hist_gen_" + genvar + "_genbin" + to_string(igen +1)+ "_sys_";
      hSign_bin[igen][0]= (TH1D*)sigFile->Get(newtitle.c_str());
      newtitle = "hist_gen_" + genvar + "_genbin" + to_string(igen +1);
      hSign_bin[igen][0]->SetTitle(newtitle.c_str());
      hSign_bin[igen][0]->SetName(newtitle.c_str());
      hSign_bin[igen][0]->SetDirectory(0);
    }
  
  
  sigFile->Close();

	
  TFile* dataFile = TFile::Open(datafile.c_str(), "READ");
  TH1D* hGen_altmc = (TH1D*)dataFile->Get(("hist_gen_" + genvar).c_str());
  hGen_altmc->SetDirectory(0);

  TH1D* hData = (TH1D*)dataFile->Get(("hist_reco_" + recovar).c_str());
  hData->SetDirectory(0);

  TH1D* hReco_Data = (TH1D*)dataFile->Get(("hist_reco_" + recovar).c_str());
  hReco_Data->SetTitle(("hist_gen_" + genvar + "_data_obs").c_str());
  hReco_Data->SetName(("hist_gen_" + genvar + "_data_obs").c_str());
  hReco_Data->SetDirectory(0);
  
  //TH2D* hResp_finebins = (TH2D*)sigFile->Get("hist_responsematrix_recogen_ak81p_finebinst");
  //hResp_finebins->SetDirectory(0);

  //TH2D* hist_reso = (TH2D*)sigFile->Get("hist_delptbypt_genak8");
  //hist_reso->SetDirectory(0);
 
  dataFile->Close();

  hGen->Scale(lumi * ttbar_signal_wgt);
  hResp->Scale(lumi * ttbar_signal_wgt);
  hReco->Scale(lumi * ttbar_signal_wgt);
  for(int igen = 0; igen < 6; igen++)
    {
      for(int isys= 0; isys < 11; isys++)
	{
	  hSign_bin[igen][isys]->Scale(lumi * ttbar_signal_wgt);
	}
    }
  hGen_altmc->Scale(lumi * data_wgt);
  hData->Scale(lumi * data_wgt);
  hReco_Data->Scale(lumi * data_wgt);
  
  get_bin_migration_plots(hResp,genvar,gentitle);
  
  TH1D *hists_recoandgen[2];
  hists_recoandgen[0] = (TH1D*)hReco->Clone();
  hists_recoandgen[1] = (TH1D*)hGen->Clone();

  hGen->Rebin(2);
  hGen_altmc->Rebin(2);
  hReco->Rebin(1);
  hData->Rebin(1);
  hResp->Rebin2D(1,2);

  TH1D *hists_reco[2] = {hReco,hData};
  string dataname_reco[2] = {"reco MC","reco Data"};
  TH1D *hists_gen[2] = {hGen,hGen_altmc};
  string dataname_gen[2] = {"gen MC","gen Data"};

  hists_recoandgen[0]->Rebin(1);
  hists_recoandgen[1]->Rebin(1);
  scale_bin_by_binwidth(hists_recoandgen[0],1);
  scale_bin_by_binwidth(hists_recoandgen[1],1);
  hists_recoandgen[0]->GetXaxis()->SetTitle(gentitle.c_str());  

  string dataname_recoandgen[2] = {"reco","gen"};

  plot_1dhists(2,hists_recoandgen,dataname_recoandgen,"compare_" + genvar + "_recovsgen.png");
  
  plot_1dhists(2,hists_reco,dataname_reco,"compare_l1pt_recovsdata.png");
  plot_1dhists(2,hists_gen,dataname_gen,"compare_l1pt_genvsdata.png");

  plot_response_matrix(hResp,genvar,gentitle,recotitle);
  make_input_histforunfolding(hData,hGen_altmc,hResp,hReco_Data,6,11,hSign_bin);
  //make_input_histforunfolding(hReco,hGen,hResp);
 
  
  //cout<<"enddl"<<endl;
  //get_resolution_value(hist_reso);
  
  
  /*TFile* sigFile = TFile::Open("MuMu_Summer20UL18_unfolding_input.root", "READ");
  TH1D* hGen = (TH1D*)sigFile->Get(("hist_gen_" + varname).c_str());
  hGen->SetDirectory(0);

  TH1D* hReco = (TH1D*)sigFile->Get(("hist_reco_" + varname).c_str());
  hReco->SetDirectory(0);

  TH2D* hResp = (TH2D*)sigFile->Get("hist_responsematrix_recogen_ak81pt");
  hResp->SetDirectory(0);
  sigFile->Close();
  get_unfold_hist(hGen,hReco,hResp);
  
  */  
}

void prep_datacards_for_jackknife_unc(string genvar, string gentitle, string mcfile, float ttbar_signal_wgt, string recovar, string recotitle, string datafile, float data_wgt)
{
  TFile* sigFile = TFile::Open(mcfile.c_str(), "READ");
  //TH1D* hGen = (TH1D*)sigFile->Get(("hist_gen_" + genvar).c_str());
  //hGen->SetDirectory(0);

  //TH1D* hReco = (TH1D*)sigFile->Get(("hist_reco_" + recovar).c_str());
  //hReco->SetDirectory(0);
  TH2D* hResp[11];
  
  hResp[0] = (TH2D*)sigFile->Get(("hist_responsematrix_recogen_" + genvar).c_str());
  hResp[0]->SetDirectory(0);

  for(int ih=1; ih < 11; ih++)
    {
      hResp[ih] = (TH2D*)sigFile->Get(("hist_responsematrix_recogen_" + genvar + "_set" + to_string(ih)).c_str());
      hResp[ih]->SetDirectory(0);
    }

  sigFile->Close();

  TFile* dataFile = TFile::Open(datafile.c_str(), "READ");
  TH1D* hGen_altmc[11];
  for(int ih=0; ih < 10; ih++)
    {
      hGen_altmc[ih] = (TH1D*)dataFile->Get(TString::Format("hist_gen_%s_set%d",genvar.c_str(),ih+1).Data());
      hGen_altmc[ih]->SetDirectory(0);
    }

  hGen_altmc[10] = (TH1D*)dataFile->Get(TString::Format("hist_gen_%s",genvar.c_str()).Data());
  hGen_altmc[10]->SetDirectory(0);

  TH1D* hData[11];
  for(int ih=0; ih < 10; ih++)
    {
      hData[ih] = (TH1D*)dataFile->Get(TString::Format("hist_reco_%s_set%d",recovar.c_str(),ih+1).Data());
      hData[ih]->SetDirectory(0);
    }

  hData[10] = (TH1D*)dataFile->Get(TString::Format("hist_reco_%s",recovar.c_str()).Data());
  hData[10]->SetDirectory(0);

  dataFile->Close();

  //hGen->Scale(lumi * ttbar_signal_wgt);
  //hResp->Scale(lumi * ttbar_signal_wgt);
  //hReco->Scale(lumi * ttbar_signal_wgt);
  for(int ih=0; ih < 11; ih++)
    {
      hResp[ih]->Scale(lumi * ttbar_signal_wgt);
      hGen_altmc[ih]->Scale(lumi * data_wgt);
      hData[ih]->Scale(lumi * data_wgt);
    }
  
  //hGen->Rebin(2);
  //hReco->Rebin(1);
  for(int ih=0; ih < 11; ih++)
    {
      hGen_altmc[ih]->Rebin(2);
      hData[ih]->Rebin(1);
      hResp[ih]->Rebin2D(1,2);
    }

  //plot_response_matrix(hResp,genvar,gentitle,recotitle);
  for(int ih=0; ih < 10; ih++)
    {
      TH2D *hResp_copy = (TH2D*)hResp[0]->Clone();
      hResp_copy->SetName(("hist_responsematrix_recogen_" + genvar + "_set" + to_string(ih+1)).c_str());
      make_input_histforunfolding(hData[ih],hGen_altmc[ih],hResp_copy);
    }
  
  for(int ih=0; ih < 10; ih++)
    {
      TH1D *hData_copy = (TH1D*)hData[10]->Clone();
      hData_copy->SetName(TString::Format("hist_reco_%s_set%d",recovar.c_str(),ih+11).Data());
      TH1D *hGen_copy = (TH1D*)hGen_altmc[10]->Clone();
      hGen_copy->SetName(TString::Format("hist_gen_%s_set%d",genvar.c_str(),ih+11).Data());
      hResp[ih+1]->SetName(("hist_responsematrix_recogen_" + genvar + "_set" + to_string(ih+11)).c_str());
      cout<<hData_copy->GetName()<<"-"<<hGen_copy->GetName()<<"-"<<hResp[ih+1]->GetName()<<endl;
      make_input_histforunfolding(hData_copy,hGen_copy,hResp[ih+1]);
    }  
}





void addhisttofile(string outfile, string outhist, string infile, string inhist)
{
    TFile *fin = new TFile(infile.c_str(),"read");
    if(!fin->IsOpen())
    {cout<<"Check: file "<<infile<<" not present"<<endl;exit(0);}
    TH1D *hist = (TH1D*)fin->Get(inhist.c_str());
    if(hist == NULL)
    {cout<<"Check: file "<<infile<<" does not have Tobject "<<inhist<<endl;exit(0);}
    TFile *fout = new TFile(outfile.c_str(),"update");
    fout->cd();
    hist->SetTitle(outhist.c_str());
    hist->SetName(outhist.c_str());
    hist->Write(outhist.c_str(),TObject::kOverwrite);
    fout->Close();
    delete hist;
    fin->Close();

}

void makedatacard_rootfile(string varname, string reconame)
{
  string sysnames[11] = {"","pileupUp","pileupDown","top_tagger_SFUp","top_tagger_SFDown","trigger_SFUp","trigger_SFDown","JERUp","JERDown","JESUp","JESDown"};

  	for(int igen = 0; igen < 6; igen++)
	  {
	    for(int isys= 1; isys < 11; isys++)
	      {
		addhisttofile("input_combine_unfolding_" + varname + "_withsys.root","hist_gen_" + varname + "_genbin" + to_string(igen +1) + "_sys_" + sysnames[isys],"MuMu_Summer20UL18_TTBar_DiLeptonic_signal_output.root","hist_gen_" + varname + "_genbin" + to_string(igen +1)+ "_sys_" + sysnames[isys]);
	      }
	    addhisttofile("input_combine_unfolding_" + varname + "_withsys.root","hist_gen_" + varname + "_genbin" + to_string(igen +1),"MuMu_Summer20UL18_TTBar_DiLeptonic_signal_output.root","hist_gen_" + varname + "_genbin" + to_string(igen +1)+ "_sys_");
	    
	  }
	
     	addhisttofile("input_combine_unfolding_" + varname + "_withsys.root","hist_gen_" + varname + "_data_obs","MuMu_Summer20UL18_TTBar_DiLeptonic_signal_output.root","hist_reco_" + reconame);
 
}

void unfolding_function()
{
  string mcfile = "MuMu_Summer20UL18_TTBar_DiLeptonic_signal_output.root";
  double ttbar_signal_wgt = 0.00000842839823203154;

  //string altmcfile = "MuMu_Summer20UL18_TTBar_inclusive_signal_output.root";
  //double alt_ttbar_wgt = 0.00288178878903645 ; /////// TT inclusive
  
  string altmcfile = "MuMu_Summer20UL18_TTJets_DiLept_signal_output.root";
  double alt_ttbar_wgt = 0.00305971804433108;  /////// TTjets dilep  * 52.69/88.29

  string datafile = "mergerd_MuMu_Summer20UL18_Single_Muon_output.root";
  
  string recovarnames[] = {"l1pt","l2pt","ak8jet1pt","ak8jet2pt"};
  string genvarnames[] = {"l1pt","l2pt","blep1pt","blep2pt","top1pt","top2pt","ak81pt", "ak82pt"};

  string gentitle[] = {"gen lepton pt from top quark (in GeV)","gen lepton pt from antitop quark (in GeV)","gen (b quark + lepton) system pt from top quark (in GeV)","gen (b quark + lepton) system pt from antitop quark (in GeV)","top quark pt (in GeV)","antitop quark pt (in GeV)","gen leading AK8 pt (in GeV)","gen subleading AK8 pt (in GeV)" };
  string recotitle[] = {"reco positively charged lepton pt (in GeV)", "reco negatively charged lepton pt (in GeV)","reco AK8 jet (containing positively charged lepton) pt (in GeV)","reco AK8 jet (containing negatively charged lepton) pt (in GeV)"};

  all_unfolding_prep(genvarnames[0],gentitle[0],mcfile,ttbar_signal_wgt,recovarnames[0],recotitle[0],altmcfile,alt_ttbar_wgt);
  all_unfolding_prep(genvarnames[2],gentitle[2],mcfile,ttbar_signal_wgt,recovarnames[2],recotitle[2],altmcfile,alt_ttbar_wgt);
  all_unfolding_prep(genvarnames[4],gentitle[4],mcfile,ttbar_signal_wgt,recovarnames[2],recotitle[2],altmcfile,alt_ttbar_wgt);
  
  /*all_unfolding_prep(genvarnames[0],gentitle[0],mcfile,ttbar_signal_wgt,recovarnames[0],recotitle[0],altmcfile,alt_ttbar_wgt);
  all_unfolding_prep(genvarnames[1],gentitle[1],mcfile,ttbar_signal_wgt,recovarnames[1],recotitle[1],altmcfile,alt_ttbar_wgt);
  
  all_unfolding_prep(genvarnames[2],gentitle[2],mcfile,ttbar_signal_wgt,recovarnames[2],recotitle[2],altmcfile,alt_ttbar_wgt);
  all_unfolding_prep(genvarnames[3],gentitle[3],mcfile,ttbar_signal_wgt,recovarnames[3],recotitle[3],altmcfile,alt_ttbar_wgt);
  
  all_unfolding_prep(genvarnames[4],gentitle[4],mcfile,ttbar_signal_wgt,recovarnames[2],recotitle[2],altmcfile,alt_ttbar_wgt);
  all_unfolding_prep(genvarnames[5],gentitle[5],mcfile,ttbar_signal_wgt,recovarnames[3],recotitle[3],altmcfile,alt_ttbar_wgt);

  all_unfolding_prep(genvarnames[6],gentitle[6],mcfile,ttbar_signal_wgt,recovarnames[2],recotitle[2],altmcfile,alt_ttbar_wgt);
  all_unfolding_prep(genvarnames[7],gentitle[7],mcfile,ttbar_signal_wgt,recovarnames[3],recotitle[3],altmcfile,alt_ttbar_wgt); */
  
  /*all_unfolding_prep(genvarnames[0],gentitle[0],mcfile,ttbar_signal_wgt,recovarnames[0],recotitle[0],mcfile,ttbar_signal_wgt);
  all_unfolding_prep(genvarnames[2],gentitle[2],mcfile,ttbar_signal_wgt,recovarnames[2],recotitle[2],mcfile,ttbar_signal_wgt);
  all_unfolding_prep(genvarnames[4],gentitle[4],mcfile,ttbar_signal_wgt,recovarnames[2],recotitle[2],mcfile,ttbar_signal_wgt);*/
 
  /*all_unfolding_prep(genvarnames[0],gentitle[0],altmcfile,alt_ttbar_wgt,recovarnames[0],recotitle[0],altmcfile,alt_ttbar_wgt);
  all_unfolding_prep(genvarnames[2],gentitle[2],altmcfile,alt_ttbar_wgt,recovarnames[2],recotitle[2],altmcfile,alt_ttbar_wgt);
  all_unfolding_prep(genvarnames[4],gentitle[4],altmcfile,alt_ttbar_wgt,recovarnames[2],recotitle[2],altmcfile,alt_ttbar_wgt);
  */

  /*
  all_unfolding_prep(genvarnames[0],gentitle[0],altmcfile,alt_ttbar_wgt,recovarnames[0],recotitle[0],mcfile,ttbar_signal_wgt);
  all_unfolding_prep(genvarnames[2],gentitle[2],altmcfile,alt_ttbar_wgt,recovarnames[2],recotitle[2],mcfile,ttbar_signal_wgt);
  all_unfolding_prep(genvarnames[4],gentitle[4],altmcfile,alt_ttbar_wgt,recovarnames[2],recotitle[2],mcfile,ttbar_signal_wgt);*/
  
  get_unfold_hist(genvarnames[2],gentitle[2],recovarnames[2]);
  
  /*get_unfold_hist(genvarnames[0],gentitle[0],recovarnames[0]);
  get_unfold_hist(genvarnames[1],gentitle[1],recovarnames[1]);

  get_unfold_hist(genvarnames[2],gentitle[2],recovarnames[2]);
  get_unfold_hist(genvarnames[3],gentitle[3],recovarnames[3]);

  get_unfold_hist(genvarnames[4],gentitle[4],recovarnames[2]);
  get_unfold_hist(genvarnames[5],gentitle[5],recovarnames[3]);

  get_unfold_hist(genvarnames[6],gentitle[6],recovarnames[2]);
  get_unfold_hist(genvarnames[7],gentitle[7],recovarnames[3]);
  */
  
  /*get_unfold_hist_withoutfirstbin(genvarnames[0],gentitle[0],recovarnames[0]);
  get_unfold_hist_withoutfirstbin(genvarnames[1],gentitle[1],recovarnames[1]);
  
  get_unfold_hist_withoutfirstbin(genvarnames[2],gentitle[2],recovarnames[2]);
  get_unfold_hist_withoutfirstbin(genvarnames[3],gentitle[3],recovarnames[3]);

  get_unfold_hist_withoutfirstbin(genvarnames[4],gentitle[4],recovarnames[2]);
  get_unfold_hist_withoutfirstbin(genvarnames[5],gentitle[5],recovarnames[3]);
  
  get_unfold_hist_withoutfirstbin(genvarnames[6],gentitle[6],recovarnames[2]);
  get_unfold_hist_withoutfirstbin(genvarnames[7],gentitle[7],recovarnames[3]);
  */

  
  /*prep_datacards_for_jackknife_unc(genvarnames[0],gentitle[0],mcfile,ttbar_signal_wgt,recovarnames[0],recotitle[0],altmcfile,alt_ttbar_wgt);  
  prep_datacards_for_jackknife_unc(genvarnames[2],gentitle[2],mcfile,ttbar_signal_wgt,recovarnames[2],recotitle[2],altmcfile,alt_ttbar_wgt);
  prep_datacards_for_jackknife_unc(genvarnames[4],gentitle[4],mcfile,ttbar_signal_wgt,recovarnames[2],recotitle[2],altmcfile,alt_ttbar_wgt);
  */

  /* get_jackknife_stat_error_ondata(genvarnames[0],gentitle[0],recovarnames[0]);
  get_jackknife_stat_error_ondata(genvarnames[2],gentitle[2],recovarnames[2]);
  get_jackknife_stat_error_ondata(genvarnames[4],gentitle[4],recovarnames[2]);
  
  get_jackknife_stat_error_onresponsematrix(genvarnames[0],gentitle[0],recovarnames[0]);
  get_jackknife_stat_error_onresponsematrix(genvarnames[2],gentitle[2],recovarnames[2]);
  get_jackknife_stat_error_onresponsematrix(genvarnames[4],gentitle[4],recovarnames[2]);
  */
  
  //get_resolution_value("MuMu_Summer20UL18_TTBar_DiLeptonic_signal_output.root","hist_resolution_gentoppt");

  
  /*TFile* sigFile = TFile::Open(mcfile.c_str(), "READ");
  TH1D* hist_ak81pt_spectrum[5];
  TH1D* hist_ak82pt_spectrum[5];
  const char *ak8pt_spectrum_checknames[10] = {"hist_ak81pt_before_toptaggercut_withoutSF","hist_ak81pt_after_toptaggercut_withoutSF","hist_ak81pt_after_toptaggercut_withSF","hist_ak81pt_before_toptaggercut_withoutSF","hist_ak81pt_after_toptaggercut_withoutSF","hist_ak82pt_before_toptaggercut_withoutSF","hist_ak82pt_after_toptaggercut_withoutSF","hist_ak82pt_after_toptaggercut_withSF","hist_ak82pt_before_toptaggercut_withoutSF","hist_ak82pt_after_toptaggercut_withoutSF"};
  
  for(int ih =0; ih< 3; ih++)
    {
      hist_ak81pt_spectrum[ih] = (TH1D*)sigFile->Get(ak8pt_spectrum_checknames[ih]);
      hist_ak81pt_spectrum[ih]->SetDirectory(0);
      hist_ak81pt_spectrum[ih]->Scale(ttbar_signal_wgt * lumi);
      //hist_ak81pt_spectrum[ih]->Scale(1/hist_ak81pt_spectrum[ih]->Integral());
      hist_ak81pt_spectrum[ih]->GetXaxis()->SetTitle("p_{T} of leading AK 8 jet (in GeV)");
	
      hist_ak82pt_spectrum[ih] = (TH1D*)sigFile->Get(ak8pt_spectrum_checknames[ih+5]);
      hist_ak82pt_spectrum[ih]->SetDirectory(0);
      hist_ak82pt_spectrum[ih]->Scale(ttbar_signal_wgt * lumi);
      //hist_ak82pt_spectrum[ih]->Scale(1/1/hist_ak82pt_spectrum[ih]->Integral());
      hist_ak82pt_spectrum[ih]->GetXaxis()->SetTitle("p_{T} of subleading AK 8 jet (in GeV)");
    }
  sigFile->Close();

  TFile* dataFile = TFile::Open("mergerd_MuMu_Summer20UL18_Single_Muon_output.root", "READ");
  for(int ih =3; ih< 5; ih++)
    {
      hist_ak81pt_spectrum[ih] = (TH1D*)dataFile->Get(ak8pt_spectrum_checknames[ih]);
      hist_ak81pt_spectrum[ih]->SetDirectory(0);
      //hist_ak81pt_spectrum[ih]->Scale(1/hist_ak81pt_spectrum[ih]->Integral());
      hist_ak81pt_spectrum[ih]->GetXaxis()->SetTitle("p_{T} of leading AK 8 jet (in GeV)");
	
      hist_ak82pt_spectrum[ih] = (TH1D*)dataFile->Get(ak8pt_spectrum_checknames[ih+5]);
      hist_ak82pt_spectrum[ih]->SetDirectory(0);
      //hist_ak82pt_spectrum[ih]->Scale(1/1/hist_ak82pt_spectrum[ih]->Integral());
      hist_ak82pt_spectrum[ih]->GetXaxis()->SetTitle("p_{T} of subleading AK 8 jet (in GeV)");
    }
  dataFile->Close();

  string ak81pt_spectrum_legends[5] = {"MC - Without tagger cut and tagger SF","MC - With tagger cut and without tagger SF","MC - With tagger cut and tagger SF","Data - Without tagger cut","Data - With tagger cut"};
  string ak82pt_spectrum_legends[5] = {"MC - Without tagger cut and tagger SF","MC - With tagger cut and without tagger SF","MC - With tagger cut and tagger SF","Data - Without tagger cut","Data - With tagger cut"};
  
  plot_1dhists_withratio(5,hist_ak81pt_spectrum,ak81pt_spectrum_legends,"plot_ak81_pt_spectrum_withtaggerconditions.png");
  plot_1dhists_withratio(5,hist_ak82pt_spectrum,ak82pt_spectrum_legends,"plot_ak82_pt_spectrum_withtaggerconditions.png");*/
  
}








/*
void make_datacards(string genvar, string recovar)
{

  ofstream runscript("allvars_unfold_combine.sh");

    std::cout << genvar << " " << recovar << std::endl;
    
    TFile* sigFile = TFile::Open("MuMu_Summer20UL18_unfolding_input.root", "READ");
    TH1D* hGen = (TH1D*)sigFile->Get(("hist_gen_" + genvar).c_str());
    hGen->SetDirectory(0);
    TH1D* hReco = (TH1D*)sigFile->Get(("hist_reco_" + recovar).c_str());
    hReco->SetDirectory(0);
    TH1D* hResp = (TH1D*)sigFile->Get(("hist_responsematrix_recogen_" + genvar).c_str());
    hResp->SetDirectory(0);
    sigFile->Close();

    ofstream datacard("unfolding_" + genvar + "_datacard.txt");
    datacard << "* imax\n";
    datacard << "* jmax\n";
    datacard << "* kmax\n";

    int reco_minbin = 1;
    int reco_maxbin = hReco->GetNbinsX();
    int gen_minbin = 1;
    int gen_maxbin = hGen->GetNbinsX() + 1;
    double min_events_in_response_matrix = 0.0001;
    
    datacard << "----------------\n";
    datacard << "bin ";
    for (int ireco = reco_minbin; ireco <= reco_maxbin; ireco++) {
      datacard << "Reco_" << ireco << " ";
    }
    datacard << "\n";
    datacard << "observation ";
    for (int ireco = reco_minbin; ireco <= reco_maxbin; ireco++) {
      datacard << hReco->GetBinContent(ireco) << " ";
    }
    datacard << "\n";
    datacard << "----------------\n";
    cout << gen_maxbin << "  " << hResp->GetBinContent(2, 0) << endl;

    bool cleanup = true;
    bool hasbkg = false;
    datacard << "bin ";
    for (int ireco = reco_minbin; ireco <= reco_maxbin; ireco++)
      {
	for (int igen = gen_minbin; igen <= gen_maxbin; igen++)
	  {
	    if (cleanup && hResp->GetBinContent(ireco, igen) < min_events_in_response_matrix) 
	      continue;
	    
	    datacard << "Reco_" << ireco << " ";
	  }
      if (hasbkg)
	{
	  datacard << "Reco_" << ireco << " ";
	}
      }
    datacard << "\n";

    datacard << "process ";
    for (int ireco = reco_minbin; ireco <= reco_maxbin; ireco++) {
      for (int igen = gen_minbin; igen <= gen_maxbin; igen++) {
        if (cleanup && hResp->GetBinContent(ireco,igen) < min_events_in_response_matrix) continue;
        datacard << -1*igen << " ";
        if (hasbkg) datacard << "1 ";
      }
    }
    
    datacard << "process ";
    for (int ireco = reco_minbin; ireco <= reco_maxbin; ireco++) {
      for (int igen = gen_minbin; igen <= gen_maxbin; igen++) {
        if (cleanup && hResp->GetBinContent(ireco,igen) < min_events_in_response_matrix) continue;
        datacard << "Gen_" << igen << " ";
        if (hasbkg) datacard << "Bkg ";
      }
    }
    datacard << "\n";
    datacard << "rate ";
    for (int ireco = reco_minbin; ireco <= reco_maxbin; ireco++) {
      for (int igen = gen_minbin; igen <= gen_maxbin; igen++) {
        if (cleanup && hResp->GetBinContent(ireco,igen) < min_events_in_response_matrix) continue;
        datacard << std::fixed << std::setprecision(5) << hResp->GetBinContent(ireco,igen) << " ";
        if (hasbkg) datacard << std::fixed << std::setprecision(5) << hBkg.GetBinContent(ireco) << " ";
      }
    }
    datacard << "\n";
    datacard << "----------------\n";
    for (int ireco = reco_minbin; ireco <= reco_maxbin; ireco++) {
      datacard << "Reco_" << ireco << " autoMCStats 0\n";
    }
    datacard << "----------------\n";
    std::string po = "";
    for (int igen = gen_minbin; igen < gen_maxbin; ++igen) {
      po += " --PO map='.* /Gen_" + std::to_string(igen) + ":r_bin" + std::to_string(igen) + "[1,0,20]'";
    }
    datacard.write("\n");
    datacard.write("\n");
    runscript.write("text2workspace.py unfolding_" + genvar + "_datacard.txt -o unfolding_" + genvar + "_input.root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel " + po + " --X-allow-no-background\n");
    datacard.write("\n");
    runscript.write("combine -M MultiDimFit --algo singles -d unfolding_" + genvar + "_input.root -t 0 --saveWorkspace --saveFitResult\n");
    datacard.write("\n");
    runscript.write("mv higgsCombineTest.MultiDimFit.mH120.root higgsCombineTest.MultiDimFit." + genvar + ".root\n\n"); 
  }
}

void get_jackknife_stat_error_onresponsematrix(string genvar, string gentitle, string recovar)
{
  TFile* sigFile = TFile::Open("MuMu_Summer20UL18_unfolding_input.root", "READ");
  TH1D* hGen_allbins[11];
  hGen_allbins[0] = (TH1D*)sigFile->Get(("hist_gen_" + genvar).c_str());
  hGen_allbins[0]->SetDirectory(0);

  for(int ih=1; ih < 11; ih++)
    {
      hGen_allbins[ih] = (TH1D*)sigFile->Get(TString::Format("hist_gen_%s_set%d",genvar.c_str(),ih+10).Data());
      hGen_allbins[ih]->SetDirectory(0);
    }
  
  TH1D* hReco_allbins[11];
  hReco_allbins[0] = (TH1D*)sigFile->Get(TString::Format("hist_reco_%s",recovar.c_str()).Data());
  hReco_allbins[0]->SetDirectory(0);
  
  for(int ih=1; ih < 11; ih++)
    {
      hReco_allbins[ih] = (TH1D*)sigFile->Get(TString::Format("hist_reco_%s_set%d",recovar.c_str(),ih+10).Data());
      hReco_allbins[ih]->SetDirectory(0);
    }
  for(int ih=0; ih < 11; ih++)
    {
      cout<<hReco_allbins[ih]->Integral()<<endl;
    }
  
  TH2D* hResp_allbins[11];
  hResp_allbins[0]= (TH2D*)sigFile->Get(("hist_responsematrix_recogen_" + genvar).c_str());
  hResp_allbins[0]->SetDirectory(0);

  for(int ih=1; ih < 11; ih++)
    {
      hResp_allbins[ih]= (TH2D*)sigFile->Get(TString::Format("hist_responsematrix_recogen_%s_set%d",genvar.c_str(),ih+10).Data());
      hResp_allbins[ih]->SetDirectory(0);
    }
  
  sigFile->Close();

  TH1D* hGen[11];
  TH1D* hReco[11];
  for(int ih=0; ih < 11; ih++)
    {
      hReco[ih] = (TH1D*)removefirstbin(hReco_allbins[ih],2);
      hGen[ih] = (TH1D*)removefirstbin(hGen_allbins[ih]);
    }
  
  cout<<"end removing bins"<<endl;
  TH1D *hist_unfold[11];
      
  for(int ih=0; ih < 11; ih++)
    {
      TFile* combine_output;
      if(ih==0)
	combine_output = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar + ".root").c_str(), "READ");
      else
	combine_output = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar + "_set" + to_string(ih+10) + ".root").c_str(), "READ");
      
      TTree * tree = (TTree*)combine_output->Get("limit");
      const int nbins = int(hGen[ih]->GetNbinsX());
      float r_values[nbins];
      float r_bins[nbins];
      float r_bins_high[nbins];
      float r_bins_low[nbins];
      
      for(int ib = 0; ib < nbins; ib++)
	tree->SetBranchAddress(TString::Format("r_bin%d",ib+2),&r_values[ib]);
      
      for(int i=0; i <nbins; i++)   
	{
	  tree->GetEntry(0);
	  r_bins[i] = r_values[i];
	  tree->GetEntry(2*i+3); 
	  r_bins_low[i] = r_values[i];
	  tree->GetEntry(2*i+4); 
	  r_bins_high[i] = r_values[i];      
	}
      tree->GetEntry(0);   //////// get best fit values
      
      hist_unfold[ih] = (TH1D*)hGen[ih]->Clone();
      TH1D *hist_err = (TH1D*)hGen[ih]->Clone();
      TH1D *rat_obs_combine;
      TH1D *rat_obs_tunfold;
      
      TH1D *hist_unfold_tunfold_allbins;
      if(ih == 0)   hist_unfold_tunfold_allbins = (TH1D*)TUnfold_unfold(genvar,gentitle,recovar);
      else          hist_unfold_tunfold_allbins = (TH1D*)TUnfold_unfold(genvar+"_set"+to_string(ih+10),gentitle,recovar+"_set"+to_string(ih+10));
      TH1D *hist_unfold_tunfold = (TH1D*)removefirstbin(hist_unfold_tunfold_allbins);
      
      cout<<"Check "<<ih<<endl;
      for(int ib=0; ib < nbins; ib++)
	{
	  float bincontent = 0;//hResp_allbins[ih]->Integral(1,hResp_allbins[ih]->GetNbinsX(),ib+1,ib+1)   //hGen->GetBinContent(ib+1)
	  for(int ix = 1; ix < hResp_allbins[ih]->GetNbinsX()+1; ix++)
	    {
	      if(hResp_allbins[ih]->GetBinContent(ix,ib+2) > 0.00001)
		bincontent += hResp_allbins[ih]->GetBinContent(ix,ib+2);
	    }

	  cout<<r_bins[ib]<<" "<<r_bins_low[ib]-r_bins[ib]<<" "<<r_bins_high[ib]-r_bins[ib]<<endl;
	  cout<<"gen bin "<<ib<<" content"<<hGen[ih]->GetBinContent(ib+1)<<" "<<hGen[ih]->GetBinLowEdge(ib+1)<<" == "<<bincontent<<endl;

	  hist_unfold[ih]->SetBinContent(ib+1, bincontent*r_bins[ib]);
	  hist_unfold[ih]->SetBinError(ib+1, bincontent*(r_bins_high[ib] - r_bins_low[ib])/2);
	  
	  hist_err->SetBinContent(ib+1, abs(bincontent*r_bins[ib] - hGen[ih]->GetBinContent(ib+1))/hGen[ih]->GetBinContent(ib+1));
	}

      hist_unfold[ih]->SetFillStyle(0);
      hist_unfold[ih]->SetFillColor(0);
      hist_unfold[ih]->SetMarkerStyle(kFullCircle);
      hist_unfold[ih]->SetMarkerColor(kBlack);
      hist_unfold[ih]->SetMarkerSize(1.2);
      hist_unfold[ih]->SetLineColor(kBlack);
      hist_unfold[ih]->SetLineWidth(2);
      hist_unfold[ih]->SetLineStyle(1);
      hist_unfold[ih]->SetLineColor(0);

      hist_unfold_tunfold->SetMarkerStyle(kFullStar);
      hist_unfold_tunfold->SetMarkerColor(kBlack);
      hist_unfold_tunfold->SetLineStyle(1);
      hist_unfold_tunfold->SetMarkerSize(3);
      hist_unfold_tunfold->SetLineWidth(2);
      hist_unfold_tunfold->SetLineColor(kRed);
      
      hGen[ih]->GetXaxis()->SetTitle(gentitle.c_str());
      hGen[ih]->SetLineWidth(3);
      //hGen[ih]->SetLineStyle(2);
      hGen[ih]->SetLineColor(kBlue);
  
      hReco[ih]->SetLineWidth(3);
      //hReco->SetLineStyle(2);
      hReco[ih]->SetLineColor(kGreen);
      hReco[ih]->Rebin(2);

      scale_bin_by_binwidth(hist_unfold[ih],1);
      scale_bin_by_binwidth(hist_unfold_tunfold,1);
      scale_bin_by_binwidth(hGen[ih],1);
      scale_bin_by_binwidth(hReco[ih],1);
      
  
      //hGen->SetMaximum(1.9*max(max(hGen->GetMaximum(),hReco->GetMaximum()),hist_unfold->GetMaximum()));

      hist_err->SetFillStyle(0);
      //hist_err->SetMarkerColor(0);
      hist_err->SetLineColor(0);
      hist_err->SetMarkerSize(1.2);
      
      TCanvas *c1 = tdrDiCanvas("c1",hGen[ih],hReco[ih],8,0);;
      c1->cd(1);
      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      //gPad->SetLogy();
      gStyle->SetPaintTextFormat( "1.2f" );
      
      TLegend *legend = new TLegend(0.34,0.6,0.74,0.85);
      legend->SetBorderSize(0);
      legend->SetTextFont(42);
      legend->SetTextSize(0.035);
      legend->AddEntry(hist_unfold[ih], "Unfolded  distribution using Combine", "pe2");
      legend->AddEntry(hist_unfold_tunfold, "Unfolded  distribution using Tunfold", "l");
      legend->AddEntry(hGen[ih], "Gen Madgraph+Pythia sample", "l");
      legend->AddEntry(hReco[ih], "Reco Madgraph+Pythia sample", "l");
      
      hGen[ih]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen[ih]->GetBinContent(hGen[ih]->GetMaximumBin()),hReco[ih]->GetBinContent(hReco[ih]->GetMaximumBin())),hist_unfold[ih]->GetBinContent(hist_unfold[ih]->GetMaximumBin())));
      
      hGen[ih]->Draw("e1");
      hReco[ih]->Draw("e1same");
      hist_unfold_tunfold->Draw("histsame");
      hist_unfold[ih]->DrawCopy("Psame");
      hist_unfold[ih]->SetFillColor(kBlack);
      hist_unfold[ih]->SetFillStyle(3013);
      hist_unfold[ih]->Draw("e2same");
      //hist_unfold->Draw("Psame");
      //hist_err->Draw("Psame");
      legend->Draw("same");

      CMS_lumi( c1, 8, 0 );
      gPad->Modified();


      c1->cd(2);
      
      rat_obs_combine = (TH1D*)hist_unfold[ih]->Clone();
      rat_obs_combine->Divide(hGen[ih]);
      
      rat_obs_tunfold = (TH1D*)hist_unfold_tunfold->Clone();
      rat_obs_tunfold->Divide(hGen[ih]);
      rat_obs_tunfold->SetLineColor(kRed);
      
      rat_obs_combine->GetXaxis()->SetTitle(gentitle.c_str());
      rat_obs_combine->GetXaxis()->SetTitleSize(0.09);
      rat_obs_combine->GetXaxis()->SetLabelSize(0.1);
      rat_obs_combine->GetXaxis()->CenterTitle();
      rat_obs_combine->GetXaxis()->SetNdivisions(406);

      rat_obs_combine->GetYaxis()->SetTitle("Unfolded / Gen");
      rat_obs_combine->GetYaxis()->SetTitleSize(0.1);
      rat_obs_combine->GetYaxis()->SetLabelSize(0.1);
      rat_obs_combine->GetYaxis()->SetTitleOffset(0.75);
      rat_obs_combine->GetYaxis()->SetNdivisions(406);
      rat_obs_combine->GetYaxis()->CenterTitle();
  
      rat_obs_combine->SetMinimum(0);
      rat_obs_combine->SetMaximum(2.);
      
      rat_obs_combine->SetFillStyle(0);
      rat_obs_combine->SetFillColor(0);
      rat_obs_combine->SetMarkerStyle(kFullCircle);
      rat_obs_combine->SetMarkerColor(kBlack);
      rat_obs_combine->SetMarkerSize(0.7);
      rat_obs_combine->SetLineColor(kBlack);
      
      rat_obs_tunfold->SetFillStyle(0);
      rat_obs_tunfold->SetFillColor(0);
      rat_obs_tunfold->SetMarkerStyle(kFullCircle);
      rat_obs_tunfold->SetMarkerColor(kRed);
      rat_obs_tunfold->SetMarkerSize(0.7);
      rat_obs_tunfold->SetLineColor(kRed);
      
      rat_obs_combine->Draw("e,p");
      rat_obs_tunfold->Draw("esame");
      rat_obs_combine->Draw("e,psame");
  
      TLegend *legend_ratio = new TLegend(0.77,0.74,0.93,0.87);
      legend_ratio->SetBorderSize(0);
      legend_ratio->SetTextFont(42);
      legend_ratio->SetTextSize(0.045);
      legend_ratio->AddEntry(rat_obs_combine, "Combine","l");
      legend_ratio->AddEntry(rat_obs_tunfold, "Tunfold","l");
      legend_ratio->Draw("same");
      
      TLine *line = new TLine(rat_obs_combine->GetXaxis()->GetXmin(),1,rat_obs_combine->GetXaxis()->GetXmax(),1);
      line->SetLineColor(1);
      line->Draw("same");
      
      if(ih==0) c1->SaveAs(TString::Format("plots/withoutfirstbin_responsematrix_unfolded_distribution_%s_higgscombine.png",genvar.c_str()));
      else c1->SaveAs(TString::Format("plots/withoutfirstbin_responsematrix_unfolded_distribution_%s_set%d_higgscombine.png",genvar.c_str(),ih));

    }
  for(int ih=0; ih < 11; ih++)
    {	   
      hist_unfold[ih]->SetFillStyle(0);
      hist_unfold[ih]->SetFillColor(0);
      hist_unfold[ih]->SetMarkerStyle(kFullCircle);
      hist_unfold[ih]->SetMarkerColor(kBlack);
      hist_unfold[ih]->SetMarkerSize(1.2);
      hist_unfold[ih]->SetLineColor(kBlack);
      hist_unfold[ih]->SetLineWidth(2);
      hist_unfold[ih]->SetLineStyle(1);
      hist_unfold[ih]->SetLineColor(ih);
      
      hist_unfold[ih]->GetXaxis()->SetTitle(gentitle.c_str());
     
      //scale_bin_by_binwidth(hist_unfold[ih],1);
    }
  
  TCanvas *c1 = tdrCanvas("c1",hist_unfold[0],8,0);;
  c1->cd();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gStyle->SetPaintTextFormat( "1.2f" );

  TLegend *legend = new TLegend(0.3,0.6,0.7,0.85);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
  legend->AddEntry(hist_unfold[0], "Unfolded  distribution", "pe2");
  for(int ih=1; ih < 11; ih++)
    legend->AddEntry(hist_unfold[ih], ("Unfolded  distribution using response matrix set " + to_string(ih)).c_str(), "l");
  
  //hGen->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen->GetBinContent(hGen->GetMaximumBin()),hReco->GetBinContent(hReco->GetMaximumBin())),hist_unfold->GetBinContent(hist_unfold->GetMaximumBin())));

  hist_unfold[0]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*hist_unfold[0]->GetBinContent(hist_unfold[0]->GetMaximumBin()));

  //hGen->Draw("e1");
  //hReco->Draw("e1same");
  //hist_unfold_tunfold->Draw("histsame");
  hist_unfold[0]->DrawCopy("Psame");
  hist_unfold[0]->SetFillColor(kBlack);
  hist_unfold[0]->SetFillStyle(3013);
  hist_unfold[0]->Draw("e2same");
  for(int ih=1; ih < 11; ih++)
    hist_unfold[ih]->Draw("histsame");
  //hist_unfold->Draw("Psame");
  //hist_err->Draw("Psame");
  legend->Draw("same");

  CMS_lumi( c1, 8, 0 );
  gPad->Modified();
  c1->SaveAs(TString::Format("plots/jackknife_responsematrix_unfolded_distribution_comparsion_%s_higgscombine.png",genvar.c_str()));


  TH1D *hist_unfold_jacknife_unc = (TH1D*)hist_unfold[0]->Clone();
  for(int ib =0; ib < hist_unfold_jacknife_unc->GetNbinsX(); ib++)
    {
      std::vector<double> bincontents;
      for(int ih =1; ih<11; ih++)
	bincontents.push_back(hist_unfold[ih]->GetBinContent(ib+1));
      double std_dev = TMath::StdDev(bincontents.begin(), bincontents.end() ); 

      hist_unfold_jacknife_unc->SetBinError(ib+1,sqrt(10/9)*std_dev);
      cout<<"Bin "<<ib<<": Content = "<<hist_unfold_jacknife_unc->GetBinContent(ib+1)<<" with combine uncertainty = "<<hist_unfold[0]->GetBinError(ib+1)<<" and jackknife uncertainty = "<<hist_unfold_jacknife_unc->GetBinError(ib+1)<<endl;
    }

  TCanvas *c2 = tdrCanvas("c2",hist_unfold[0],8,0);;
  c2->cd();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gStyle->SetPaintTextFormat( "1.2f" );

  TLegend *legend_2 = new TLegend(0.3,0.6,0.7,0.85);
  legend_2->SetBorderSize(0);
  legend_2->SetTextFont(42);
  legend_2->SetTextSize(0.03);
  legend_2->AddEntry(hist_unfold[0], "Unfolded  distribution with combine uncertainity", "pe2");
  legend_2->AddEntry(hist_unfold_jacknife_unc, "Unfolded  distribution with jackknife response matrix uncertainity", "l");
   
  //hGen->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen->GetBinContent(hGen->GetMaximumBin()),hReco->GetBinContent(hReco->GetMaximumBin())),hist_unfold->GetBinContent(hist_unfold->GetMaximumBin())));

  hist_unfold[0]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*hist_unfold[0]->GetBinContent(hist_unfold[0]->GetMaximumBin()));

  //hGen->Draw("e1");
  //hReco->Draw("e1same");
  //hist_unfold_tunfold->Draw("histsame");
  hist_unfold[0]->DrawCopy("Psame");
  hist_unfold[0]->SetFillColor(kBlack);
  hist_unfold[0]->SetFillStyle(3013);
  hist_unfold[0]->Draw("e2same");
  hist_unfold_jacknife_unc->SetLineColor(kRed);
  hist_unfold_jacknife_unc->Draw("e1same");
  legend_2->Draw("same");

  CMS_lumi( c2, 8, 0 );
  gPad->Modified();
  c2->SaveAs(TString::Format("plots/jackknife_compare_responsematrixuncertainty_unfolded_distribution_%s_combinevsjackknife.png",genvar.c_str()));

}

double* get_jackknife_stat_error_ondata(string genvar, string gentitle, string recovar)
{
  TFile* sigFile = TFile::Open("MuMu_Summer20UL18_unfolding_input.root", "READ");
  TH1D* hGen_allbins[11];
  hGen_allbins[0] = (TH1D*)sigFile->Get(("hist_gen_" + genvar).c_str());
  hGen_allbins[0]->SetDirectory(0);

  for(int ih=1; ih < 11; ih++)
    {
      hGen_allbins[ih] = (TH1D*)sigFile->Get(TString::Format("hist_gen_%s_set%d",genvar.c_str(),ih).Data());
      hGen_allbins[ih]->SetDirectory(0);
    }
  
  TH1D* hReco_allbins[11];
  hReco_allbins[0] = (TH1D*)sigFile->Get(TString::Format("hist_reco_%s",recovar.c_str()).Data());
  hReco_allbins[0]->SetDirectory(0);
  
  for(int ih=1; ih < 11; ih++)
    {
      hReco_allbins[ih] = (TH1D*)sigFile->Get(TString::Format("hist_reco_%s_set%d",recovar.c_str(),ih).Data());
      hReco_allbins[ih]->SetDirectory(0);
    }
  for(int ih=0; ih < 11; ih++)
    {
      cout<<hReco_allbins[ih]->Integral()<<endl;
    }
  
  TH2D* hResp_allbins = (TH2D*)sigFile->Get(("hist_responsematrix_recogen_" + genvar).c_str());
  hResp_allbins->SetDirectory(0);
  sigFile->Close();

  TH1D* hGen[11];
  TH1D* hReco[11];
  for(int ih=0; ih < 11; ih++)
    {
      hReco[ih] = (TH1D*)removefirstbin(hReco_allbins[ih],2);
      hGen[ih] = (TH1D*)removefirstbin(hGen_allbins[ih]);
    }
  
  cout<<"end removing bins"<<endl;
  TH1D *hist_unfold[11];
      
  for(int ih=0; ih < 11; ih++)
    {
      TFile* combine_output;
      if(ih==0)
	combine_output = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar + ".root").c_str(), "READ");
      else
	combine_output = TFile::Open(("higgsCombineTest.MultiDimFit." + genvar + "_set" + to_string(ih) + ".root").c_str(), "READ");
      
      TTree * tree = (TTree*)combine_output->Get("limit");
      const int nbins = int(hGen[ih]->GetNbinsX());
      float r_values[nbins];
      float r_bins[nbins];
      float r_bins_high[nbins];
      float r_bins_low[nbins];
      
      for(int ib = 0; ib < nbins; ib++)
	tree->SetBranchAddress(TString::Format("r_bin%d",ib+2),&r_values[ib]);
      
      for(int i=0; i <nbins; i++)   
	{
	  tree->GetEntry(0);
	  r_bins[i] = r_values[i];
	  tree->GetEntry(2*i+3); 
	  r_bins_low[i] = r_values[i];
	  tree->GetEntry(2*i+4); 
	  r_bins_high[i] = r_values[i];      
	}
      tree->GetEntry(0);   //////// get best fit values
      
      hist_unfold[ih] = (TH1D*)hGen[ih]->Clone();
      TH1D *hist_err = (TH1D*)hGen[ih]->Clone();
      TH1D *rat_obs_combine;
      TH1D *rat_obs_tunfold;
      
      TH1D *hist_unfold_tunfold_allbins;
      if(ih == 0)   hist_unfold_tunfold_allbins = (TH1D*)TUnfold_unfold(genvar,gentitle,recovar);
      else          hist_unfold_tunfold_allbins = (TH1D*)TUnfold_unfold(genvar+"_set"+to_string(ih),gentitle,recovar+"_set"+to_string(ih));
      TH1D *hist_unfold_tunfold = (TH1D*)removefirstbin(hist_unfold_tunfold_allbins);
      
      cout<<"Check "<<ih<<endl;
      for(int ib=0; ib < nbins; ib++)
	{
	  float bincontent = 0;//hResp_allbins->Integral(1,hResp_allbins->GetNbinsX(),ib+1,ib+1)   //hGen->GetBinContent(ib+1)
	  for(int ix = 1; ix < hResp_allbins->GetNbinsX()+1; ix++)
	    {
	      if(hResp_allbins->GetBinContent(ix,ib+2) > 0.00001)
		bincontent += hResp_allbins->GetBinContent(ix,ib+2);
	    }

	  cout<<r_bins[ib]<<" "<<r_bins_low[ib]-r_bins[ib]<<" "<<r_bins_high[ib]-r_bins[ib]<<endl;
	  cout<<"gen bin "<<ib<<" content"<<hGen[ih]->GetBinContent(ib+1)<<" "<<hGen[ih]->GetBinLowEdge(ib+1)<<" == "<<bincontent<<endl;

	  hist_unfold[ih]->SetBinContent(ib+1, bincontent*r_bins[ib]);
	  hist_unfold[ih]->SetBinError(ib+1, bincontent*(r_bins_high[ib] - r_bins_low[ib])/2);
	  
	  hist_err->SetBinContent(ib+1, abs(bincontent*r_bins[ib] - hGen[ih]->GetBinContent(ib+1))/hGen[ih]->GetBinContent(ib+1));
	}

      hist_unfold[ih]->SetFillStyle(0);
      hist_unfold[ih]->SetFillColor(0);
      hist_unfold[ih]->SetMarkerStyle(kFullCircle);
      hist_unfold[ih]->SetMarkerColor(kBlack);
      hist_unfold[ih]->SetMarkerSize(1.2);
      hist_unfold[ih]->SetLineColor(kBlack);
      hist_unfold[ih]->SetLineWidth(2);
      hist_unfold[ih]->SetLineStyle(1);
      hist_unfold[ih]->SetLineColor(0);

      hist_unfold_tunfold->SetMarkerStyle(kFullStar);
      hist_unfold_tunfold->SetMarkerColor(kBlack);
      hist_unfold_tunfold->SetLineStyle(1);
      hist_unfold_tunfold->SetMarkerSize(3);
      hist_unfold_tunfold->SetLineWidth(2);
      hist_unfold_tunfold->SetLineColor(kRed);
      
      hGen[ih]->GetXaxis()->SetTitle(gentitle.c_str());
      hGen[ih]->SetLineWidth(3);
      //hGen[ih]->SetLineStyle(2);
      hGen[ih]->SetLineColor(kBlue);
  
      hReco[ih]->SetLineWidth(3);
      //hReco->SetLineStyle(2);
      hReco[ih]->SetLineColor(kGreen);
      hReco[ih]->Rebin(2);

      scale_bin_by_binwidth(hist_unfold[ih],1);
      scale_bin_by_binwidth(hist_unfold_tunfold,1);
      scale_bin_by_binwidth(hGen[ih],1);
      scale_bin_by_binwidth(hReco[ih],1);
      
  
      //hGen->SetMaximum(1.9*max(max(hGen->GetMaximum(),hReco->GetMaximum()),hist_unfold->GetMaximum()));

      hist_err->SetFillStyle(0);
      //hist_err->SetMarkerColor(0);
      hist_err->SetLineColor(0);
      hist_err->SetMarkerSize(1.2);
      
      TCanvas *c1 = tdrDiCanvas("c1",hGen[ih],hReco[ih],8,0);;
      c1->cd(1);
      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      //gPad->SetLogy();
      gStyle->SetPaintTextFormat( "1.2f" );
      
      TLegend *legend = new TLegend(0.34,0.6,0.74,0.85);
      legend->SetBorderSize(0);
      legend->SetTextFont(42);
      legend->SetTextSize(0.035);
      legend->AddEntry(hist_unfold[ih], "Unfolded  distribution using Combine", "pe2");
      legend->AddEntry(hist_unfold_tunfold, "Unfolded  distribution using Tunfold", "l");
      legend->AddEntry(hGen[ih], "Gen Madgraph+Pythia sample", "l");
      legend->AddEntry(hReco[ih], "Reco Madgraph+Pythia sample", "l");
      
      hGen[ih]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen[ih]->GetBinContent(hGen[ih]->GetMaximumBin()),hReco[ih]->GetBinContent(hReco[ih]->GetMaximumBin())),hist_unfold[ih]->GetBinContent(hist_unfold[ih]->GetMaximumBin())));
      
      hGen[ih]->Draw("e1");
      hReco[ih]->Draw("e1same");
      hist_unfold_tunfold->Draw("histsame");
      hist_unfold[ih]->DrawCopy("Psame");
      hist_unfold[ih]->SetFillColor(kBlack);
      hist_unfold[ih]->SetFillStyle(3013);
      hist_unfold[ih]->Draw("e2same");
      //hist_unfold->Draw("Psame");
      //hist_err->Draw("Psame");
      legend->Draw("same");

      CMS_lumi( c1, 8, 0 );
      gPad->Modified();


      c1->cd(2);
      
      rat_obs_combine = (TH1D*)hist_unfold[ih]->Clone();
      rat_obs_combine->Divide(hGen[ih]);
      
      rat_obs_tunfold = (TH1D*)hist_unfold_tunfold->Clone();
      rat_obs_tunfold->Divide(hGen[ih]);
      rat_obs_tunfold->SetLineColor(kRed);
      
      rat_obs_combine->GetXaxis()->SetTitle(gentitle.c_str());
      rat_obs_combine->GetXaxis()->SetTitleSize(0.09);
      rat_obs_combine->GetXaxis()->SetLabelSize(0.1);
      rat_obs_combine->GetXaxis()->CenterTitle();
      rat_obs_combine->GetXaxis()->SetNdivisions(406);

      rat_obs_combine->GetYaxis()->SetTitle("Unfolded / Gen");
      rat_obs_combine->GetYaxis()->SetTitleSize(0.1);
      rat_obs_combine->GetYaxis()->SetLabelSize(0.1);
      rat_obs_combine->GetYaxis()->SetTitleOffset(0.75);
      rat_obs_combine->GetYaxis()->SetNdivisions(406);
      rat_obs_combine->GetYaxis()->CenterTitle();
  
      rat_obs_combine->SetMinimum(0);
      rat_obs_combine->SetMaximum(2.);
      
      rat_obs_combine->SetFillStyle(0);
      rat_obs_combine->SetFillColor(0);
      rat_obs_combine->SetMarkerStyle(kFullCircle);
      rat_obs_combine->SetMarkerColor(kBlack);
      rat_obs_combine->SetMarkerSize(0.7);
      rat_obs_combine->SetLineColor(kBlack);
      
      rat_obs_tunfold->SetFillStyle(0);
      rat_obs_tunfold->SetFillColor(0);
      rat_obs_tunfold->SetMarkerStyle(kFullCircle);
      rat_obs_tunfold->SetMarkerColor(kRed);
      rat_obs_tunfold->SetMarkerSize(0.7);
      rat_obs_tunfold->SetLineColor(kRed);
      
      rat_obs_combine->Draw("e,p");
      rat_obs_tunfold->Draw("esame");
      rat_obs_combine->Draw("e,psame");
  
      TLegend *legend_ratio = new TLegend(0.77,0.74,0.93,0.87);
      legend_ratio->SetBorderSize(0);
      legend_ratio->SetTextFont(42);
      legend_ratio->SetTextSize(0.045);
      legend_ratio->AddEntry(rat_obs_combine, "Combine","l");
      legend_ratio->AddEntry(rat_obs_tunfold, "Tunfold","l");
      legend_ratio->Draw("same");
      
      TLine *line = new TLine(rat_obs_combine->GetXaxis()->GetXmin(),1,rat_obs_combine->GetXaxis()->GetXmax(),1);
      line->SetLineColor(1);
      line->Draw("same");
      
      if(ih==0) c1->SaveAs(TString::Format("plots/withoutfirstbin_data_unfolded_distribution_%s_higgscombine.png",genvar.c_str()));
      else c1->SaveAs(TString::Format("plots/withoutfirstbin_data_unfolded_distribution_%s_set%d_higgscombine.png",genvar.c_str(),ih));

    }
  for(int ih=0; ih < 11; ih++)
    {	   
      hist_unfold[ih]->SetFillStyle(0);
      hist_unfold[ih]->SetFillColor(0);
      hist_unfold[ih]->SetMarkerStyle(kFullCircle);
      hist_unfold[ih]->SetMarkerColor(kBlack);
      hist_unfold[ih]->SetMarkerSize(1.2);
      hist_unfold[ih]->SetLineColor(kBlack);
      hist_unfold[ih]->SetLineWidth(2);
      hist_unfold[ih]->SetLineStyle(1);
      hist_unfold[ih]->SetLineColor(ih);
      
      hist_unfold[ih]->GetXaxis()->SetTitle(gentitle.c_str());
     
      //scale_bin_by_binwidth(hist_unfold[ih],1);
    }
  
  TCanvas *c1 = tdrCanvas("c1",hist_unfold[0],8,0);;
  c1->cd();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gStyle->SetPaintTextFormat( "1.2f" );

  TLegend *legend = new TLegend(0.3,0.6,0.7,0.85);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
  legend->AddEntry(hist_unfold[0], "Unfolded  distribution", "pe2");
  for(int ih=1; ih < 11; ih++)
    legend->AddEntry(hist_unfold[ih], ("Unfolded  distribution using data set " + to_string(ih)).c_str(), "l");
  
  //hGen->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen->GetBinContent(hGen->GetMaximumBin()),hReco->GetBinContent(hReco->GetMaximumBin())),hist_unfold->GetBinContent(hist_unfold->GetMaximumBin())));

  hist_unfold[0]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*hist_unfold[0]->GetBinContent(hist_unfold[0]->GetMaximumBin()));

  //hGen->Draw("e1");
  //hReco->Draw("e1same");
  //hist_unfold_tunfold->Draw("histsame");
  hist_unfold[0]->DrawCopy("Psame");
  hist_unfold[0]->SetFillColor(kBlack);
  hist_unfold[0]->SetFillStyle(3013);
  hist_unfold[0]->Draw("e2same");
  for(int ih=1; ih < 11; ih++)
    hist_unfold[ih]->Draw("histsame");
  //hist_unfold->Draw("Psame");
  //hist_err->Draw("Psame");
  legend->Draw("same");

  CMS_lumi( c1, 8, 0 );
  gPad->Modified();
  c1->SaveAs(TString::Format("plots/jackknife_data_unfolded_distribution_comparsion_%s_higgscombine.png",genvar.c_str()));

  double jackknife_errors[11];

  TH1D *hist_unfold_jacknife_unc = (TH1D*)hist_unfold[0]->Clone();
  for(int ib =0; ib < hist_unfold_jacknife_unc->GetNbinsX(); ib++)
    {
      std::vector<double> bincontents;
      for(int ih =1; ih<11; ih++)
	bincontents.push_back(hist_unfold[ih]->GetBinContent(ib+1));
      double std_dev = TMath::StdDev(bincontents.begin(), bincontents.end() ); 

      hist_unfold_jacknife_unc->SetBinError(ib+1,3*std_dev);
      cout<<"Bin "<<ib<<": Content = "<<hist_unfold_jacknife_unc->GetBinContent(ib+1)<<" with combine uncertainty = "<<hist_unfold[0]->GetBinError(ib+1)<<" and jackknife uncertainty = "<<hist_unfold_jacknife_unc->GetBinError(ib+1)<<endl;
    }

  TCanvas *c2 = tdrCanvas("c2",hist_unfold[0],8,0);;
  c2->cd();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gStyle->SetPaintTextFormat( "1.2f" );

  TLegend *legend_2 = new TLegend(0.3,0.6,0.7,0.85);
  legend_2->SetBorderSize(0);
  legend_2->SetTextFont(42);
  legend_2->SetTextSize(0.035);
  legend_2->AddEntry(hist_unfold[0], "Unfolded  distribution with combine uncertainity", "pe2");
  legend_2->AddEntry(hist_unfold_jacknife_unc, "Unfolded  distribution with jackknife data uncertainity", "l");
   
  //hGen->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*max(max(hGen->GetBinContent(hGen->GetMaximumBin()),hReco->GetBinContent(hReco->GetMaximumBin())),hist_unfold->GetBinContent(hist_unfold->GetMaximumBin())));

  hist_unfold[0]->SetMaximum((gPad->GetLogy() == 0 ? 1.9 : 1400)*hist_unfold[0]->GetBinContent(hist_unfold[0]->GetMaximumBin()));

  //hGen->Draw("e1");
  //hReco->Draw("e1same");
  //hist_unfold_tunfold->Draw("histsame");
  hist_unfold[0]->DrawCopy("Psame");
  hist_unfold[0]->SetFillColor(kBlack);
  hist_unfold[0]->SetFillStyle(3013);
  hist_unfold[0]->Draw("e2same");
  hist_unfold_jacknife_unc->SetLineColor(kRed);
  hist_unfold_jacknife_unc->Draw("e1same");
  legend_2->Draw("same");

  CMS_lumi( c2, 8, 0 );
  gPad->Modified();
  c2->SaveAs(TString::Format("plots/jackknife_compare_datauncertainty_unfolded_distribution_%s_combinevsjackknife.png",genvar.c_str()));

  return jackknife_errors;
}

*/
