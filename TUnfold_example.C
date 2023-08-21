#include "TDecompSVD.h"
//#include "CLHEP/Vector/TwoVector.h"
//#include "CLHEP/Vector/ThreeVector.h"
//#include "CLHEP/Vector/LorentzVector.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include "TVector.h"
#include <vector>
#include <TF1.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TStyle.h>
#include "TPostScript.h"
#include <TPad.h>
#include <TLine.h>
#include <TRandom.h>

#include "TLegend.h"

#include "TUnfold.h"
//#include "RooUnfoldTUnfold.h"

using namespace std;
//using namespace CLHEP;


int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -1;
  for (int ix=1; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -1;
}

int LastBin_Counter(TH1D *h1){
int lastbin = 0;
for(int bn=h1->GetNbinsX()-1; bn>=0; bn--){
if((h1->GetBinContent(bn+1)) > 1.e-12){
lastbin = bn+1;
}
else { break;}
}
cout<<h1->GetName()<<" totbins "<<lastbin<<endl;
return lastbin;
}

//int main()                                                                                                                                     {
TH1* TUnfold_unfold(string genvar, string gentitle, string recovar){

const char* rootfiles[3] = {"MuMu_Summer20UL18_unfolding_input.root","MuMu_Summer20UL18_unfolding_input.root","MuMu_Summer20UL18_unfolding_input.root"};

char name[100];
 
string name_histdata = "hist_reco_" + recovar;
string name_histreco = "hist_reco_" + recovar;
string name_histgen = "hist_gen_" + genvar;
string name_histresponse = "hist_responsematrix_recogen_" + genvar;
cout<<name_histreco<<"  "<<name_histgen<<endl;
 
TH1D *hist_data;
TH1D *hist_reco;
TH1D *hist_gen;
TH1D *hist_truth;
TH2D *mat_rm;

int lastbin = 0;

for(int ij=0; ij<3; ij++){
sprintf(name,"%s",rootfiles[ij]);
TFile *file1 = new TFile(name,"read");
file1->cd();

if(ij==0) {
        sprintf(name,"%s",name_histdata.c_str());
        hist_data = (TH1D*)file1->Get(name);
/*
	sprintf(name,"Gen_MC_JetpT_ak4_EtaBin1");
	hist_truth = (TH1D*)file1->Get(name);
	cout<<hist_truth->GetMean()<<endl;
*/
        }
if(ij==1) {
        sprintf(name,"%s",name_histreco.c_str());
        hist_reco = (TH1D*)file1->Get(name);

        sprintf(name,"%s",name_histgen.c_str());
        hist_gen = (TH1D*)file1->Get(name);
        }
if(ij==2) {
        sprintf(name,"%s",name_histresponse.c_str());
        mat_rm = (TH2D*)file1->Get(name);
	for(int bn=0; bn<mat_rm->GetNbinsX(); bn++){
	mat_rm->SetBinContent(bn+1,0,0);
	}
	//mat_rm->Rebin2D(1,2);
        }
 }//ij

//const int noptbins1 = hist_gen->GetNbinsX();
 //const double *ptbins1 = hist_gen->GetXaxis()->GetXbins()->GetArray();
 const int noptbins1 = hist_gen->GetNbinsX();
 const double *ptbins1 = hist_gen->GetXaxis()->GetXbins()->GetArray();
 /* double ptbins1[noptbins1+1] = {0};
 Int_t N = hist_gen->GetXaxis()->GetNbins(); // bin 1 is the first bin, bin (N + 1) is the "overflow bin"
 for (Int_t i = 1; i <= (N + 1); i++)
   {
     ptbins1[i-1] = hist_gen->GetXaxis()->GetBinLowEdge(i);
   }
   //double ptbins1[noptbins1+1] = {200,322,370,435,570,1000};*/
 cout<<noptbins1<<endl;


 const int noptbins = hist_reco->GetNbinsX();
 const double *ptbins = hist_reco->GetXaxis()->GetXbins()->GetArray();
 //double ptbins[noptbins+1] = {300,316,342,376,420,472,534,604,684,770,1000};

 for(int it = 0; it < noptbins1; it++)
   cout<<ptbins1[it]<<",";
 cout<<endl;
 
   TH1D *hist_unf;
   TFile *fileout ;
   fileout = new TFile("TUnfold_example.root","recreate") ;
   sprintf(name,"TUnfold_%s",genvar.c_str());
   hist_unf = new TH1D(name,name,noptbins1,ptbins1);

//   hist_unf = (TH1D*)hist_data->Clone();

   lastbin = LastBin_Counter(hist_data);

//   TUnfold unfoldBbB(mat_rm,TUnfold::kHistMapOutputVert,TUnfold::kRegModeCurvature,TUnfold::kEConstraintNone);
   TUnfold unfoldBbB(mat_rm,TUnfold::kHistMapOutputVert);

//   TUnfold::ERegMode regMode=TUnfold::kRegModeCurvature;
   TUnfold::ERegMode regMode=TUnfold::kRegModeSize;
   unfoldBbB.SetInput(hist_data);

   //////Regularistion
   unfoldBbB.RegularizeBins(1,1,2,regMode);  
   Double_t tauMin=-3.0;
   Double_t tauMax=3.0;
   Int_t nScan=30;
   Int_t iBest;
   TSpline *logTauX,*logTauY;
   TGraph *lCurve;
   
   unfoldBbB.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
   std::cout<<"chi**2="<<unfoldBbB.GetChi2A()<<"+"<<unfoldBbB.GetChi2L()<<" / "<<unfoldBbB.GetNdf()<<"\n";
   

   ////// Without Regularistion
   //unfoldBbB.DoUnfold(0) ;

  Int_t *binMap=new Int_t[noptbins+2];
  for(Int_t i=1;i<=noptbins;i++) binMap[i]=i;
  binMap[0]=-1;
  binMap[noptbins+1]=-1;

  unfoldBbB.GetOutput(hist_unf,binMap);  

  hist_unf->SetMarkerStyle(kFullCircle);
  hist_unf->SetMarkerColor(kBlack);
  hist_unf->SetMarkerSize(0.8);

  hist_data->SetLineColor(kGreen);
  hist_reco->SetLineColor(kBlue);
  hist_gen->SetLineColor(kBlue);
  hist_unf->SetLineColor(kRed);
//  hist_truth->SetLineColor(7);

  hist_data->SetStats(0); hist_data->SetMinimum(0);
  hist_reco->SetStats(0);
  hist_gen->SetStats(0);
  hist_unf->SetStats(0);

  /*  for(int bn=0; bn<(hist_data->GetNbinsX()); bn++){
  hist_data->SetBinContent(bn+1,hist_data->GetBinContent(bn+1)*1./hist_data->GetBinWidth(bn+1));
  hist_reco->SetBinContent(bn+1,hist_reco->GetBinContent(bn+1)*1./hist_reco->GetBinWidth(bn+1));
  hist_gen->SetBinContent(bn+1,hist_gen->GetBinContent(bn+1)*1./hist_gen->GetBinWidth(bn+1));

  hist_data->SetBinError(bn+1,hist_data->GetBinError(bn+1)*1./hist_data->GetBinWidth(bn+1));
  hist_reco->SetBinError(bn+1,hist_reco->GetBinError(bn+1)*1./hist_reco->GetBinWidth(bn+1));
  hist_gen->SetBinError(bn+1,hist_gen->GetBinError(bn+1)*1./hist_gen->GetBinWidth(bn+1));

//  hist_truth->SetBinContent(bn+1,hist_truth->GetBinContent(bn+1)*1./hist_gen->GetBinWidth(bn+1));
  }

  for(int bn=0; bn<(hist_unf->GetNbinsX()); bn++){
  hist_unf->SetBinContent(bn+1,hist_unf->GetBinContent(bn+1)*1./hist_unf->GetBinWidth(bn+1));
  hist_unf->SetBinError(bn+1,hist_unf->GetBinError(bn+1)*1./hist_unf->GetBinWidth(bn+1));
  }
  */
  TCanvas *c1 = new TCanvas("Tunfold_output","Tunfold_output",800,600);
  TLegend *leg1 = new TLegend(0.7,0.7,0.85,0.85);
  leg1->SetBorderSize(0);
  c1->cd(1);
  
  //gPad->SetLogx(1);
  //gPad->SetLogy(1);
  hist_data->Rebin(2);

  hist_data->SetMaximum(1.9*max(max(hist_gen->GetBinContent(hist_gen->GetMaximumBin()),hist_data->GetBinContent(hist_data->GetMaximumBin())),hist_unf->GetBinContent(hist_unf->GetMaximumBin())));

  //hist_gen->Draw();  leg1->AddEntry(hist_gen,"GEN","l");
  hist_data->SetTitle("");
  hist_data->GetXaxis()->SetTitle(gentitle.c_str());
  
  hist_data->Draw("samee"); leg1->AddEntry(hist_data,"DATA","6ple1");
  hist_unf->Draw("textsames"); leg1->AddEntry(hist_unf,"Unfolded","l");
  hist_gen->Draw("sames");  leg1->AddEntry(hist_gen,"GEN","l");
  
//  hist_truth->Draw("sames"); leg1->AddEntry(hist_truth,"Unfolded","l");
  leg1->Draw();

//  TH1D *hist_gen_rat = (TH1D*)hist_gen->Clone(); hist_gen_rat->Divide(hist_data);
  
  sprintf(name,"./plots/Tunfold_output_%s.png",genvar.c_str());
  c1->SaveAs(name);
  
  hist_unf->SetDirectory(0);
    
  fileout->cd();
  fileout->Write();
  fileout->Close();

  return hist_unf;
}


void TUnfold_example(){

 string genvarnames[] = {"l1pt","l2pt","blep1pt","blep2pt","top1pt","top2pt","ak81pt", "ak82pt"};
 string recovarnames[] = {"l1pt","l2pt","ak8jet1pt","ak8jet2pt"};
 string genvar = genvarnames[0];
 string recovar = recovarnames[0];
 string gentitle[] = {"gen lepton pt from top quark (in GeV)","gen lepton pt from antitop quark (in GeV)","gen (b quark + lepton) system pt from top quark (in GeV)","gen (b quark + lepton) system pt from antitop quark (in GeV)","top quark pt (in GeV)","antitop quark pt (in GeV)","gen leading AK8 pt (in GeV)","gen subleading AK8 pt (in GeV)" };

 TUnfold_unfold(genvarnames[0],gentitle[0],recovarnames[0]);
 TUnfold_unfold(genvarnames[1],gentitle[1],recovarnames[1]);

 TUnfold_unfold(genvarnames[2],gentitle[2],recovarnames[2]);
 TUnfold_unfold(genvarnames[3],gentitle[3],recovarnames[3]);

 TUnfold_unfold(genvarnames[4],gentitle[4],recovarnames[2]);
 TUnfold_unfold(genvarnames[5],gentitle[5],recovarnames[3]);

 TUnfold_unfold(genvarnames[6],gentitle[6],recovarnames[2]);
 TUnfold_unfold(genvarnames[7],gentitle[7],recovarnames[3]);

}
