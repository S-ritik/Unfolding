/*
  weight = event_weight; 
  should be in the same ntuple, event by event we are going to read it.
    
  Order lepton according to the matching of leading Ak8 jets. 
  Lepton tagging
  
*/

#define Anal_Leptop_PROOF_cxx
//#include "Anal_Leptop_PROOF.h"
#include "getobjects.h"

#include <TH2.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <fstream>
#include <TProofOutputFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TProofServ.h>

//#define E_MU_TTBar
//#define E_E_TTBar
#define MU_MU_TTBar

void Anal_Leptop_PROOF::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();

}
 
void Anal_Leptop_PROOF::SlaveBegin(TTree * /*tree*/)
{
  //The SlaveBegin() function is called after the Begin() function.
  //When running with PROOF SlaveBegin() is called on each slave server.
  //The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  OutFile = new TProofOutputFile("test_output.root");
  
  fileOut = OutFile->OpenFile("RECREATE");
  if ( !(fileOut = OutFile->OpenFile("RECREATE")) )
    {
      Warning("SlaveBegin", "problems opening file: %s/%s",OutFile->GetDir(), OutFile->GetFileName());
    }
   
  isMC = true;
  isTT = true;
  isST = false;
  isDIB = false;
  isWJ = false;
  isDY = false;
  isbQCD = false;
  isTTX = false;
  
  //GMA
  //Tout = new TTree("leptop","leptop");
  //Tout->Branch("weight",&weight,"weight/F");
  
  Tnewvar = new TTree("newvars","newvars");

  Tnewvar->Branch("irun", &irun, "irun/I");
  Tnewvar->Branch("ilumi", &ilumi, "ilumi/I");
  Tnewvar->Branch("ievt", &ievt, "ievt/i");

  Tnewvar->Branch("event_pt_weight",&event_pt_weight,"event_pt_weight/F");
  Tnewvar->Branch("weight",&weight,"weight/F");
  Tnewvar->Branch("qscale",&qscale,"qscale/F");
  Tnewvar->Branch("M_l1l2",&M_l1l2,"M_l1l2/F");
  Tnewvar->Branch("rat_l2pt_l1pt",&rat_l2pt_l1pt,"rat_l2pt_l1pt/F");
  Tnewvar->Branch("deltaPhi_l1l2",&deltaPhi_l1l2,"deltaPhi_l1l2/F");
  Tnewvar->Branch("l1pt_nearjet",&l1pt_nearjet,"l1pt_nearjet/F");
  Tnewvar->Branch("l2pt_nearjet",&l2pt_nearjet,"l2pt_nearjet/F");
  Tnewvar->Branch("met_pt",&met_pt,"met_pt_chs/F");
  Tnewvar->Branch("met_phi",&met_phi,"met_phi_chs/F");
  Tnewvar->Branch("met_pt_puppi",&met_pt_puppi,"met_pt_puppi/F");
  Tnewvar->Branch("met_phi_puppi",&met_phi_puppi,"met_phi_puppi/F");
  Tnewvar->Branch("M_bl1",&M_bl1,"M_bl1/F"); 
  Tnewvar->Branch("M_bl2",&M_bl2,"M_bl2/F");
  Tnewvar->Branch("delta_phil1_met",&delta_phil1_met,"delta_phil1_met/F");
  Tnewvar->Branch("delta_phil2_met",&delta_phil2_met,"delta_phil2_met/F");
  Tnewvar->Branch("delta_phibl1_met",&delta_phibl1_met,"delta_phibl1_met/F"); 
  Tnewvar->Branch("delta_phibl2_met",&delta_phibl2_met,"delta_phibl2_met/F");
  Tnewvar->Branch("rat_metpt_ak4pt",&rat_metpt_ak4pt,"rat_metpt_ak4pt/F");
  Tnewvar->Branch("rat_metpt_ak8pt",&rat_metpt_ak8pt,"rat_metpt_ak8pt/F");
  Tnewvar->Branch("rat_metpt_eventHT",&rat_metpt_eventHT,"rat_metpt_eventHT/F");
  Tnewvar->Branch("mt_of_l1met",&mt_of_l1met,"mt_of_l1met/F");
  Tnewvar->Branch("mt_of_l2met",&mt_of_l2met,"mt_of_l2met/F");
  Tnewvar->Branch("no_ak4jets",&no_ak4jets,"no_ak4jets/F");
  Tnewvar->Branch("no_ak4bjets",&no_ak4bjets,"no_ak4bjets/F");
  Tnewvar->Branch("no_ak8jets",&no_ak8jets,"no_ak8jets/F");
  Tnewvar->Branch("EventHT",&EventHT,"EventHT/F");
  Tnewvar->Branch("extra_ak4j",&extra_ak4j,"extra_ak4j/F");
  Tnewvar->Branch("ptsum_extra_ak4",&ptsum_extra_ak4,"ptsum_extra_ak4/F");
  Tnewvar->Branch("extra_ak4pt",&extra_ak4pt,"extra_ak4pt/F");
  Tnewvar->Branch("extra_ak4mass",&extra_ak4mass,"extra_ak4mass/F");
  Tnewvar->Branch("extra_ak4jqgl",&extra_ak4jqgl,"extra_ak4jqgl/F");
  Tnewvar->Branch("extra_ak4jdeepb",&extra_ak4jdeepb,"extra_ak4jdeepb/F");
  Tnewvar->Branch("rat_extra_ak4jpt_lpt",&rat_extra_ak4jpt_lpt,"rat_extra_ak4jpt_lpt/F");
  Tnewvar->Branch("ak41pt",&ak41pt,"ak41pt/F");
  Tnewvar->Branch("ak41mass",&ak41mass,"ak41mass/F");
  Tnewvar->Branch("ak42pt",&ak42pt,"ak42pt/F");
  Tnewvar->Branch("ak42mass",&ak42mass,"ak42mass/F");
  Tnewvar->Branch("lep1pt",&lep1pt,"lep1pt/F");
  Tnewvar->Branch("lep2pt",&lep2pt,"lep2pt/F");
  Tnewvar->Branch("ak81pt",&ak81pt,"ak81pt/F"); 
  Tnewvar->Branch("ak81y",&ak81y,"ak81y/F");
  Tnewvar->Branch("ak81mass",&ak81mass,"ak81mass/F");
  Tnewvar->Branch("ak81sdmass",&ak81sdmass,"ak81sdmass/F");
  Tnewvar->Branch("ak81chrad",&ak81chrad,"ak81chrad/F");
  Tnewvar->Branch("ak81NHad",&ak81NHad,"ak81NHad/F");                
  Tnewvar->Branch("ak81neuhad",&ak81neuhad,"ak81neuhad/F");
  Tnewvar->Branch("ak81subhaddiff",&ak81subhaddiff,"ak81subhaddiff/F");
  Tnewvar->Branch("ak81tau21",&ak81tau21,"ak81tau21/F");
  Tnewvar->Branch("ak81tau32",&ak81tau32,"ak81tau32/F");
  Tnewvar->Branch("ak81deep_tvsqcd",&ak81deep_tvsqcd,"ak81deep_tvsqcd/F");
  Tnewvar->Branch("ak81deep_wvsqcd",&ak81deep_wvsqcd,"ak81deep_wvsqcd/F");
  Tnewvar->Branch("ak81hase",&ak81hase,"ak81hase/I"); 
  Tnewvar->Branch("ak81hasmu",&ak81hasmu,"ak81hasmu/I");
  Tnewvar->Branch("ak82pt",&ak82pt,"ak82pt/F");
  Tnewvar->Branch("ak82y",&ak82y,"ak82y/F");
  Tnewvar->Branch("ak82mass",&ak82mass,"ak82mass/F");
  Tnewvar->Branch("ak82sdmass",&ak82sdmass,"ak82sdmass/F");
  Tnewvar->Branch("ak82chrad",&ak82chrad,"ak82chrad/F");
  Tnewvar->Branch("ak82NHad",&ak82NHad,"ak82NHad/F");
  Tnewvar->Branch("ak82neuhad",&ak82neuhad,"ak82neuhad/F");
  Tnewvar->Branch("ak82subhaddiff",&ak82subhaddiff,"ak82subhaddiff/F");
  Tnewvar->Branch("ak82tau21",&ak82tau21,"ak82tau21/F");
  Tnewvar->Branch("ak82tau32",&ak82tau32,"ak82tau32/F");
  Tnewvar->Branch("ak82deep_tvsqcd",&ak82deep_tvsqcd,"ak82deep_tvsqcd/F");
  Tnewvar->Branch("ak82deep_wvsqcd",&ak82deep_wvsqcd,"ak82deep_wvsqcd/F");
  Tnewvar->Branch("ak82hase",&ak82hase,"ak82hase/I"); 
  Tnewvar->Branch("ak82hasmu",&ak82hasmu,"ak82hasmu/I");
  Tnewvar->Branch("delta_phibl1bl2",&delta_phibl1bl2,"delta_phibl1bl2/F"); 
  
  //Tnewvar->Branch("delta_phijl1jl2",&delta_phijl1jl2,"delta_phijl1jl2/F");
  Tnewvar->Branch("deltaR_l1l2",&deltaR_l1l2,"deltaR_l1l2/F");
  
  //Tnewvar->Branch("deltaR_l1b1",&deltaR_l1b1,"deltaR_l1b1/F");
  //Tnewvar->Branch("deltaR_l2b1",&deltaR_l2b1,"deltaR_l2b1/F");
  //Tnewvar->Branch("deltaR_l1b2",&deltaR_l1b2,"deltaR_l1b2/F");
  //Tnewvar->Branch("deltaR_l2b2",&deltaR_l2b2,"deltaR_l2b2/F");
  
  Tnewvar->Branch("deltaR_l1j1",&deltaR_l1j1,"deltaR_l1j1/F");
  Tnewvar->Branch("deltaR_l2j1",&deltaR_l2j1,"deltaR_l2j1/F");
  Tnewvar->Branch("deltaR_l1j2",&deltaR_l1j2,"deltaR_l1j2/F");
  Tnewvar->Branch("deltaR_l2j2",&deltaR_l2j2,"deltaR_l2j2/F");
  Tnewvar->Branch("deltaR_j1j2",&deltaR_j1j2,"deltaR_j1j2/F");
  Tnewvar->Branch("j1_btag_sc",&j1_btag_sc,"j1_btag_sc/F"); 
  Tnewvar->Branch("j2_btag_sc",&j2_btag_sc,"j2_btag_sc/F");

  Tnewvar->Branch("j1_btag_sc_ptsort",&j1_btag_sc_ptsort,"j1_btag_sc_ptsort/F"); 
  Tnewvar->Branch("j2_btag_sc_ptsort",&j2_btag_sc_ptsort,"j2_btag_sc_ptsort/F");

  Tnewvar->Branch("genmatch_sc1",&genmatch_sc1,"genmatch_sc1/I");	
  Tnewvar->Branch("genmatch_sc2",&genmatch_sc2,"genmatch_sc2/I");
  
  Tnewvar->Branch("dirgltrthr",&dirgltrthr,"dirgltrthr/F");
  Tnewvar->Branch("dirglthrmin",&dirglthrmin,"dirglthrmin/F");
  
  Tnewvar->Branch("response_ak81",&response_ak81,"response_ak81/F");
  Tnewvar->Branch("response_ak82",&response_ak82,"response_ak82/F");

  Tnewvar->Branch("tt_decay_mode",&tt_decay_mode,"tt_decay_mode/F");
  Tnewvar->Branch("tau_decay_mode",&tau_decay_mode,"tau_decay_mode/F");
  
  Tnewvar->Branch("response_eventclass",&response_eventclass,"response_eventclass/F");

  Tnewvar->Branch("pt_gen_l1",&pt_gen_l1,"pt_gen_l1/F");
  Tnewvar->Branch("pt_gen_b1",&pt_gen_b1,"pt_gen_b1/F");
  Tnewvar->Branch("pt_gen_lb1",&pt_gen_lb1,"pt_gen_lb1");
  Tnewvar->Branch("deltaR_gen_lb1",&deltaR_gen_lb1,"deltaR_gen_lb1/F");
  Tnewvar->Branch("delptbypt_gen_l1",&delptbypt_gen_l1,"delptbypt_gen_l1/F");
  Tnewvar->Branch("delptbypt_gen_b1",&delptbypt_gen_b1,"delptbypt_gen_b1/F");
  Tnewvar->Branch("delptbypt_gen_lb1",&delptbypt_gen_lb1,"delptbypt_gen_lb1/F");

  Tnewvar->Branch("pt_gen_l2",&pt_gen_l2,"pt_gen_l2/F");
  Tnewvar->Branch("pt_gen_b2",&pt_gen_b2,"pt_gen_b2/F");
  Tnewvar->Branch("pt_gen_lb2",&pt_gen_lb2,"pt_gen_lb2/F");
  Tnewvar->Branch("deltaR_gen_lb2",&deltaR_gen_lb2,"deltaR_gen_lb2/F");
  Tnewvar->Branch("delptbypt_gen_l2",&delptbypt_gen_l2,"delptbypt_gen_l2/F");
  Tnewvar->Branch("delptbypt_gen_b2",&delptbypt_gen_b2,"delptbypt_gen_b2/F");
  Tnewvar->Branch("delptbypt_gen_lb2",&delptbypt_gen_lb2,"delptbypt_gen_lb2/F");

  Tnewvar->Branch("ak81hasgene",&ak81hasgene,"ak81hasgene/I"); 
  Tnewvar->Branch("ak81hasgenmu",&ak81hasgenmu,"ak81hasgenmu/I");
  Tnewvar->Branch("ak81hasgenb",&ak81hasgenb,"ak81hasgenb/I");
  Tnewvar->Branch("ak81hasgenhasalldecay",&ak81hasgenhasalldecay,"ak81hasgenhasalldecay/I");
  
  Tnewvar->Branch("ak82hasgene",&ak82hasgene,"ak82hasgene/I"); 
  Tnewvar->Branch("ak82hasgenmu",&ak82hasgenmu,"ak82hasgenmu/I");
  Tnewvar->Branch("ak82hasgenb",&ak82hasgenb,"ak82hasgenb/I");
  Tnewvar->Branch("ak82hasgenhasalldecay",&ak82hasgenhasalldecay,"ak82hasgenhasalldecay/I");

  Tnewvar->Branch("hlt_Mu37Ele27",&hlt_Mu37Ele27,"hlt_Mu37Ele27/O");
  Tnewvar->Branch("hlt_Mu27Ele37",&hlt_Mu27Ele37,"hlt_Mu27Ele37/O");
  Tnewvar->Branch("hlt_Mu37TkMu27",&hlt_Mu37TkMu27,"hlt_Mu37TkMu27/O");
  Tnewvar->Branch("hlt_DoubleEle25",&hlt_DoubleEle25,"hlt_DoubleEle25/O");
  Tnewvar->Branch("hlt_Ele50_PFJet165",&hlt_Ele50_PFJet165,"hlt_Ele50_PFJet165/O");
  Tnewvar->Branch("hlt_Mu50",&hlt_Mu50,"hlt_Mu50/O");
  
  Tnewvar->Branch("leptonsf_weight",&leptonsf_weight,"leptonsf_weight/F"); 
  Tnewvar->Branch("leptonsf_weight_stat",&leptonsf_weight_stat,"leptonsf_weight_stat/F"); 
  Tnewvar->Branch("leptonsf_weight_syst",&leptonsf_weight_syst,"leptonsf_weight_syst/F");

  Tnewvar->Branch("puWeight",&puWeight,"puWeight/F"); 
  Tnewvar->Branch("puWeightup",&puWeightup,"puWeightup/F"); 
  Tnewvar->Branch("puWeightdown",&puWeightdown,"puWeightdown/F");

  Tnewvar->Branch("btag_SF",&btag_SF,"btag_SF/F");
  Tnewvar->Branch("btag_SF_up",&btag_SF_up,"btag_SF_up/F");
  Tnewvar->Branch("btag_SF_dn",&btag_SF_dn,"btag_SF_dn/F");
  //Tnewvar->Branch("",&,"/F");

  calib_deepflav = BTagCalibration("DeepJet", (dir + "DeepJet_106XUL18SF_WPonly_V1p1.csv").Data());
  reader_deepflav_med = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"}); 
  reader_deepflav_med.load(calib_deepflav, BTagEntry::FLAV_B, "comb");
  reader_deepflav_med.load(calib_deepflav, BTagEntry::FLAV_C, "comb");
  reader_deepflav_med.load(calib_deepflav, BTagEntry::FLAV_UDSG, "incl");

  reader_deepflav_loose = BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"}); 
  reader_deepflav_loose.load(calib_deepflav, BTagEntry::FLAV_B, "comb");
  reader_deepflav_loose.load(calib_deepflav, BTagEntry::FLAV_C, "comb");
  reader_deepflav_loose.load(calib_deepflav, BTagEntry::FLAV_UDSG, "incl");
  
  char name[1000];

  for (int ij=0; ij<nprvar; ij++) {
    hist_prvar[ij] = new TH1D(prvar_name[ij], prvar_name[ij], prvar_bins[ij], prvar_low[ij], prvar_high[ij]);
  }
	
  
  for (int ij=0; ij <ntypes; ij++) {
    for (int jk=0; jk<ntcount; jk++) {
      sprintf (name, "%s_cnt%i", ptprvar_name[ij], jk);
      hist_prptvar[ij][jk] = new TH2D(name, name, 60, -3.0, 3.0, nybin, ybins); 
    }
  }
  
  for (int ij=0; ij<npr_angle; ij++) {
    if (ij<2) {
      hist_prptangle[ij] = new TH2D(pr_angle[ij], pr_angle[ij], 60, 30.0, 510.0, 120, 0.0, 4.8);
    } else if (ij<8) {
      hist_prptangle[ij] = new TH2D(pr_angle[ij], pr_angle[ij], 60, 180.0, 840.0, 120, 0.0, 4.8);
    } else {
      hist_prptangle[ij] = new TH2D(pr_angle[ij], pr_angle[ij], 60, 30.0, 510.0, 120, 0.0, 2.0);
    }
  }

  for(int init=0; init<nhist_in; init++){
    char namein[1000]; //nameinup[1000], nameindn[1000];
    char titlein[1000];
    
    sprintf(namein,"hist_%s",initnames[init]);
    sprintf(titlein,"%s",titlenames[init]);
    hist_init[init] = new TH1D(namein,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
    hist_init[init]->Sumw2();
    
    /*
      sprintf(nameinup,"hist_%s_puup",initnames[init]);
      sprintf(nameindn,"hist_%s_pudn",initnames[init]);
      hist_init_pu_sys[init][0] = new TH1D(nameinup,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
      hist_init_pu_sys[init][0]->Sumw2();
      hist_init_pu_sys[init][1] = new TH1D(nameindn,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
      hist_init_pu_sys[init][1]->Sumw2();
    */
  }
  
  for(int lvar=0; lvar<nobshist; lvar++){
    char lnamein[1000];
    sprintf(lnamein,"Obs_%s",obsnames[lvar]);
    hist_obs[lvar] = new TH1D(lnamein,lnamein,obs_nbins[lvar],obs_low[lvar],obs_up[lvar]);
    hist_obs[lvar]->Sumw2();
  }
  
  for(int init=0; init<6; init++){
    char gen_namein[1000]; //nameinup[1000], nameindn[1000];
    
    if(init <2)
      {
	sprintf(gen_namein,"hist_gen_%s",unfoldvar_initnames[init]);
	hist_genlevl[init][0] = new TH1D(gen_namein,gen_namein,12,unfold_recolepbins_finebins);
	for(int ih=1; ih < 11; ih++)
	  {
	    sprintf(gen_namein,"hist_gen_%s_set%d",unfoldvar_initnames[init],ih);
	    hist_genlevl[init][ih] = new TH1D(gen_namein,gen_namein,12,unfold_recolepbins_finebins);
	    hist_genlevl[init][ih]->Sumw2();
	  }
	for(int igen = 0; igen < 6; igen++)
	  {
	    for(int isys= 0; isys < nsys_names; isys++)
	      {
		sprintf(gen_namein,"hist_gen_%s_genbin%d_sys_%s",unfoldvar_initnames[init],igen+1,sysnames[isys].c_str());
		hist_genlevl_binned[init][igen][isys] = new TH1D(gen_namein,gen_namein,12,unfold_recolepbins_finebins);
		hist_genlevl_binned[init][igen][isys]->Sumw2();
	      }
	  }
      }
    else
      {
	sprintf(gen_namein,"hist_gen_%s",unfoldvar_initnames[init]);
	hist_genlevl[init][0] = new TH1D(gen_namein,gen_namein,12,unfold_recoak8bins_finebins);
	for(int ih=1; ih < 11; ih++)
	  {
	     sprintf(gen_namein,"hist_gen_%s_set%d",unfoldvar_initnames[init],ih);
	     hist_genlevl[init][ih] = new TH1D(gen_namein,gen_namein,12,unfold_recoak8bins_finebins);
	     hist_genlevl[init][ih]->Sumw2();	   
	  }
	for(int igen = 0; igen < 6; igen++)
	  {
	    for(int isys= 0; isys < nsys_names; isys++)
	      {
		sprintf(gen_namein,"hist_gen_%s_genbin%d_sys_%s",unfoldvar_initnames[init],igen+1,sysnames[isys].c_str());
		hist_genlevl_binned[init][igen][isys] = new TH1D(gen_namein,gen_namein,12,unfold_recoak8bins_finebins);
		hist_genlevl_binned[init][igen][isys]->Sumw2();
	      }
	  }
      }
    //hist_genlevl[init]->Rebin(6);
    hist_genlevl[init][0]->Sumw2();
  }
  
  for(int init=0; init<4; init++){
    char rec_namein[1000];
    if(init <2)
      {
	sprintf(rec_namein,"hist_reco_%s",unfold_recovar_initnames[init]);
	hist_recolevl[init][0] = new TH1D(rec_namein,rec_namein,12,unfold_recolepbins_finebins);
	for(int ih=1; ih < 11; ih++)
	  {
	    sprintf(rec_namein,"hist_reco_%s_set%d",unfold_recovar_initnames[init],ih);
	    hist_recolevl[init][ih] = new TH1D(rec_namein,rec_namein,12,unfold_recolepbins_finebins);
	    hist_recolevl[init][ih]->Sumw2();
	  }
      //hist_recolevl[init] = new TH1D(rec_namein,rec_namein,80,0,800);
      }
    else
      {
	sprintf(rec_namein,"hist_reco_%s",unfold_recovar_initnames[init]);
	hist_recolevl[init][0] = new TH1D(rec_namein,rec_namein,12,unfold_recoak8bins_finebins);
	for(int ih=1; ih < 11; ih++)
	  {
	    sprintf(rec_namein,"hist_reco_%s_set%d",unfold_recovar_initnames[init],ih);
	    hist_recolevl[init][ih] = new TH1D(rec_namein,rec_namein,12,unfold_recoak8bins_finebins);
	    hist_recolevl[init][ih]->Sumw2();
	  }   
	//hist_recolevl[init] = new TH1D(rec_namein,rec_namein,80,200,1000);
	//hist_recolevl[init]->Rebin(3);
      }
      hist_recolevl[init][0]->Sumw2();
  }


  hist_response_matrix[0][0] = new TH2D("hist_responsematrix_recogen_l1pt","hist_responsematrix_recogen_lep1pt",12,unfold_recolepbins_finebins,12,unfold_recolepbins_finebins);
  hist_response_matrix[1][0] = new TH2D("hist_responsematrix_recogen_l2pt","hist_responsematrix_recogen_lep2pt",12,unfold_recolepbins_finebins,12,unfold_recolepbins_finebins);
  hist_response_matrix[2][0] = new TH2D("hist_responsematrix_recogen_blep1pt","hist_responsematrix_recogen_blep1pt",12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);
  hist_response_matrix[3][0] = new TH2D("hist_responsematrix_recogen_blep2pt","hist_responsematrix_recogen_blep2pt",12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);
  hist_response_matrix[4][0] = new TH2D("hist_responsematrix_recogen_top1pt","hist_responsematrix_recogen_top1pt",12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);
  hist_response_matrix[5][0] = new TH2D("hist_responsematrix_recogen_top2pt","hist_responsematrix_recogen_top2pt",12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);

  hist_response_matrix[0][0]->Sumw2();
  hist_response_matrix[1][0]->Sumw2();
  hist_response_matrix[2][0]->Sumw2();
  hist_response_matrix[3][0]->Sumw2();
  hist_response_matrix[4][0]->Sumw2();
  hist_response_matrix[5][0]->Sumw2();

    for(int ih=1; ih < 11; ih++)
    {
      hist_response_matrix[0][ih] = new TH2D(TString::Format("hist_responsematrix_recogen_l1pt_set%d",ih).Data(),TString::Format("hist_responsematrix_recogen_lep1pt_set%d",ih).Data(),12,unfold_recolepbins_finebins,12,unfold_recolepbins_finebins);
      hist_response_matrix[1][ih] = new TH2D(TString::Format("hist_responsematrix_recogen_l2pt_set%d",ih).Data(),TString::Format("hist_responsematrix_recogen_lep2pt_set%d",ih).Data(),12,unfold_recolepbins_finebins,12,unfold_recolepbins_finebins);
      hist_response_matrix[2][ih] = new TH2D(TString::Format("hist_responsematrix_recogen_blep1pt_set%d",ih).Data(),TString::Format("hist_responsematrix_recogen_blep1pt_set%d",ih).Data(),12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);
      hist_response_matrix[3][ih] = new TH2D(TString::Format("hist_responsematrix_recogen_blep2pt_set%d",ih).Data(),TString::Format("hist_responsematrix_recogen_blep2pt_set%d",ih).Data(),12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);
      hist_response_matrix[4][ih] = new TH2D(TString::Format("hist_responsematrix_recogen_top1pt_set%d",ih).Data(),TString::Format("hist_responsematrix_recogen_top1pt_set%d",ih).Data(),12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);
      hist_response_matrix[5][ih] = new TH2D(TString::Format("hist_responsematrix_recogen_top2pt_set%d",ih).Data(),TString::Format("hist_responsematrix_recogen_top2pt_set%d",ih).Data(),12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);

      hist_response_matrix[0][ih]->Sumw2();
      hist_response_matrix[1][ih]->Sumw2();
      hist_response_matrix[2][ih]->Sumw2();
      hist_response_matrix[3][ih]->Sumw2();
      hist_response_matrix[4][ih]->Sumw2();
      hist_response_matrix[5][ih]->Sumw2();     
    }
    
  /*hist_response_matrix[4]->Rebin2D(3,6);
  hist_response_matrix[5]->Rebin2D(3,6);
  hist_response_matrix[4]->Rebin2D(3,6);
  hist_response_matrix[5]->Rebin2D(3,6);
  hist_response_matrix[4]->Rebin2D(3,6);
  hist_response_matrix[5]->Rebin2D(3,6);
  */
  
  hist_resolution_unfolding[0] = new TH2D("hist_resolution_genleppt","hist_resolution_genleppt",30,0,900,300,-0.4,0.4);
  hist_resolution_unfolding[1] = new TH2D("hist_resolution_genlepphi","hist_resolution_genlepphi",10,-3.14,3.14,20,-0.2,0.2);
  hist_resolution_unfolding[2] = new TH2D("hist_resolution_genlepeta","hist_resolution_genlepeta",10,-5,5,20,-0.2,0.2 );

  hist_resolution_unfolding[3] = new TH2D("hist_resolution_genbleppt","hist_resolution_genbleppt",30,0,1500,300,-0.4,0.4);
  hist_resolution_unfolding[4] = new TH2D("hist_resolution_genblepphi","hist_resolution_genblepphi",10,-3.14,3.14,20,-0.2,0.2);
  hist_resolution_unfolding[5] = new TH2D("hist_resolution_genblepeta","hist_resolution_genblepeta",10,-5,5,20,-0.2,0.2);
  
  hist_resolution_unfolding[6] = new TH2D("hist_resolution_gentoppt","hist_resolution_gentoppt",30,0,1500,300,-0.4,0.4);
  hist_resolution_unfolding[7] = new TH2D("hist_resolution_gentopphi","hist_resolution_gentopphi",10,-3.14,3.14,20,-0.2,0.2);
  hist_resolution_unfolding[8] = new TH2D("hist_resolution_gentopeta","hist_resolution_gentopeta",10,-5,5,20,-0.2,0.2);

  hist_resolution_unfolding[1]->Sumw2();
  hist_resolution_unfolding[2]->Sumw2();
  hist_resolution_unfolding[3]->Sumw2();
  hist_resolution_unfolding[4]->Sumw2();
  hist_resolution_unfolding[5]->Sumw2();
  hist_resolution_unfolding[6]->Sumw2();
  hist_resolution_unfolding[7]->Sumw2();
  hist_resolution_unfolding[8]->Sumw2();
    
  hist_delr_genvsreco[0] = new TH1D("hist_delr_genlep_recolep","hist_delr_genlep_recolep",300,0,1);
  hist_delr_genvsreco[1] = new TH1D("hist_delr_genblep_recoak8","hist_delr_genblep_recoak8",300,0,3);
  hist_delr_genvsreco[2] = new TH1D("hist_delr_gentop_recoak8","hist_delr_gentop_recoak8",300,0,3);

  hist_delr_genvsreco[0]->Sumw2();
  hist_delr_genvsreco[1]->Sumw2();
  hist_delr_genvsreco[2]->Sumw2();
    
  hist_particlelevl[0] = new TH1D("hist_gen_ak81pt","hist_gen_ak81pt",12,unfold_recoak8bins_finebins);
  //hist_particlelevl[0]->Rebin(6);
  hist_particlelevl[0]->Sumw2();

  hist_particlelevl[1] = new TH1D("hist_gen_ak82pt","hist_gen_ak82pt",12,unfold_recoak8bins_finebins);
  //hist_particlelevl[1]->Rebin(6);
  hist_particlelevl[0]->Sumw2();

  hist_response_matrix_particlelvl[0] = new TH2D("hist_responsematrix_recogen_ak81pt","hist_responsematrix_recogen_ak81pt",12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);
  //hist_response_matrix_particlelvl[0]->Rebin2D(3,6);

  hist_response_matrix_particlelvl[1] = new TH2D("hist_responsematrix_recogen_ak82pt","hist_responsematrix_recogen_ak82pt",12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);
  //hist_response_matrix_particlelvl[1]->Rebin2D(3,6);

  hist_response_matrix_finebins[0] = new TH2D("hist_responsematrix_recogen_blep1pt_finebins","hist_responsematrix_recogen_blep1pt",12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);
  hist_response_matrix_finebins[1] = new TH2D("hist_responsematrix_recogen_top1pt_finebins","hist_responsematrix_recogen_top1pt",12,unfold_recoak8bins_finebins,12,unfold_recoak8bins_finebins);
  
  for(int lvar=0; lvar<ntreevarhist; lvar++){
    char lnamein[1000];
    sprintf(lnamein,"Var_%s",treevarnames[lvar]);
    hist_treevar[lvar] = new TH1D(lnamein,lnamein,treevar_nbins[lvar],treevar_low[lvar],treevar_up[lvar]);
    hist_treevar[lvar]->Sumw2();
  }

  /*  for(int lvar=0; lvar<nhistbtag; lvar++){
      char lnamein[1000];
      sprintf(lnamein,"Obs_%s",obsnames_btag[lvar]);
      hist_obs_btag[lvar] = new TH1D(lnamein,lnamein,obs_nbins_btag[lvar],obs_low_btag[lvar],obs_up_btag[lvar]);
      hist_obs_btag[lvar]->Sumw2();
      }*/

  for(int lvar=0; lvar<ncutflow; lvar++){
    char lnamein[1000],ltitlein[1000];
    sprintf(lnamein,"hist_inv_mass_cut_%s",((string)to_string(lvar)).c_str());
    sprintf(ltitlein,"pt of leading AK8 jet after cut no %s",((string)to_string(lvar+1)).c_str());
    hist_met_cutflow[lvar] = new TH1D(lnamein,ltitlein,60,0,3000);
    hist_met_cutflow[lvar]->Sumw2();		
  }

  for(int lvar=0; lvar<nhistdeltaR; lvar++){
    char lnamein[1000],ltitlein[1000];
    sprintf(lnamein,"hist_%s",names_genmatch_deltaR[lvar]);
    sprintf(ltitlein,"%s",titles_genmatch_deltaR[lvar]);
    hist_genmatch_deltaR[lvar] = new TH1D(lnamein,ltitlein,120,0,0.01+0.1*(lvar%4));
    hist_genmatch_deltaR[lvar]->Sumw2();		
  }

  for(int lvar=0; lvar<nhistdrsc; lvar++){
    char lnamein[1000],ltitlein[1000];
    sprintf(lnamein,"hist_%s",names_genmatch_deltaR_score[lvar]);
    sprintf(ltitlein,"%s",titles_genmatch_deltaR_score[lvar]);
    hist_genmatch_deltaR_score[lvar] = new TH1D(lnamein,ltitlein,223,-0.5,222.5);
    hist_genmatch_deltaR_score[lvar]->Sumw2();		
  }
  hist_2d_deltaR_vsbtagsc[0]=new TH2D("hist_2D_deltaR_vsbtagsc_1","deltaR(leading AK4 jet, leading AK8 jet) vs b tag score of leading AK4 jet",60,0,1,60,0,6);
  hist_2d_deltaR_vsbtagsc[0]->Sumw2();
  hist_2d_deltaR_vsbtagsc[1]=new TH2D("hist_2D_deltaR_vsbtagsc_2","deltaR(sub-leading AK4 jet, sub-leading AK8 jet) vs b tag score of sub-leading AK4 jet",60,0,1,60,0,6);
  hist_2d_deltaR_vsbtagsc[1]->Sumw2();
  
  hist_2d_pt_vsbtagsc[0]=new TH2D("hist_2D_ak4pt_vsbtagsc_1","pt of leading AK4 jet vs b tag score of leading AK4 jet",120,0,1000,60,0,1);
  hist_2d_pt_vsbtagsc[0]->Sumw2();
  hist_2d_pt_vsbtagsc[1]=new TH2D("hist_2D_ak4pt_vsbtagsc_2","pt of sub leading AK4 jet vs b tag score of sub-leading AK4 jet",120,0,1000,60,0,1);
  hist_2d_pt_vsbtagsc[1]->Sumw2();
  hist_2d_pt_vsbtagsc[2]=new TH2D("hist_2D_ak8pt_vsbtagsc_1","pt of leading AK8 jet vs b tag score of leading AK4 jet",120,0,2000,60,0,1);
  hist_2d_pt_vsbtagsc[2]->Sumw2();
  hist_2d_pt_vsbtagsc[3]=new TH2D("hist_2D_ak8pt_vsbtagsc_2","pt of sub leading AK8 jet vs b tag score of sub-leading AK4 jet",120,0,2000,60,0,1);
  hist_2d_pt_vsbtagsc[3]->Sumw2();
  
  hist_2d_deltaR_vspt[0]=new TH2D("hist_2d_deltaR_vspt_1","deltaR(leading AK4 jet, leading lep) vs pt of leading AK8 jet",120,200,2500,120,0,1.5);
  hist_2d_deltaR_vspt[0]->Sumw2();  
  hist_2d_deltaR_vspt[1]=new TH2D("hist_2d_deltaR_vspt_2","deltaR(sub-leading AK4 jet, sub-leading lep) vs pt of sub-leading AK8 jet",120,200,2500,120,0,1.5);
  hist_2d_deltaR_vspt[1]->Sumw2();  
  
  /*
    for(int nvar=0; nvar<15; nvar++){
    char namein_nvar[1000], 
    char titlein_nvar[1000];
    
    sprintf(namein_nvar,"hist_%s",new_var_names[nvar]);
    sprintf(titlein_nvar,"%s",new_var_title[nvar]);
    hist_new_var[nvar] = new TH1D(namein_nvar,titlein_nvar,new_var_nbins[nvar],new_var_low[nvar],new_var_up[nvar]);
    hist_new_var[nvar]->Sumw2();
    }
  */

  hist_npu_vert_true = new TH1D("hist_npu_vert_true","npu_vert_true",100,-0.5,99.5);//80,-0.1,79.9);
  hist_npu_vert_true->Sumw2();

  hist_npu_vert = new TH1D("hist_npu_vert","npu_vert",100,-0.5,99.5);//80,-0.1,79.9);
  hist_npu_vert->Sumw2();

  char title[1000];
  sprintf(name,"N_PV");
  sprintf(title,"# of Primary Vertices");
  hist_npv = new TH1D(name,title,100,-0.5,99.5);//80,-0.1,79.9);
  hist_npv->Sumw2();
    
  sprintf(name,"N_PV_nopuwt");
  sprintf(title,"# of Primary Vertices");
  hist_npv_nopuwt = new TH1D(name,title,100,-0.5,99.5);
  hist_npv_nopuwt->Sumw2();
  
  hist_2D_msd_deepak8 = new TH2D("hist_2D_msd_deepak8","hist_2D_msd_deepak8",25,0,300,25,0,1);
  hist_2D_msd_deepak8->Sumw2();
  hist_2D_bpass_flavb = new TH2D("h2d_btagpass_flavb","h2d_btagpass_flavb",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_bpass_flavb->Sumw2();
  hist_2D_bpass_flavc = new TH2D("h2d_btagpass_flavc","h2d_btagpass_flavc",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_bpass_flavc->Sumw2();
  hist_2D_bpass_flavq = new TH2D("h2d_btagpass_flavq","h2d_btagpass_flavq",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_bpass_flavq->Sumw2();
  hist_2D_ball_flavb = new TH2D("h2d_flavb","h2d_flavb",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_ball_flavb->Sumw2();
  hist_2D_ball_flavc = new TH2D("h2d_flavc","h2d_flavc",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_ball_flavc->Sumw2();
  hist_2D_ball_flavq = new TH2D("h2d_flavq","h2d_flavq",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_ball_flavq->Sumw2();  

  hist_count = new TH1D("Counter","Counter",22,-0.5,21.5);
  hist_count->Sumw2();  
  
  hist_count_lep[1] = new TH1D("hist_tt_cutflow","Event with both top pt > 300",22,-0.5,21.5);
  hist_count_lep[1]->Sumw2();
  hist_count_top    = new TH1D("hist_cut_eff_tt","Total no of ttbar events and surviving no. of ttbar events",5,-0.5,4.5);
  hist_count_top->Sumw2();  
  
  hist_2d_pt_genlepvsb1 = new TH2D("hist_2d_ptgenlep_ptgenb","pt of gen lep vs pt of gen b ",120,20,1000,120,20,1000);
  hist_2d_pt_genlepvsb1->Sumw2();

  hist_2d_pt_genlepvsb2 = new TH2D("hist_2d_ptgenlep_ptgenblep","pt of gen lep vs pt of gen (b + lep) ",120,20,2000,120,20,1000);
  hist_2d_pt_genlepvsb2->Sumw2();

  hist_2d_pt_genlepvsb3 = new TH2D("hist_2d_ptgenlep_ptgentop","pt of gen lep vs pt of gen top ",120,20,2000,120,20,1000);
  hist_2d_pt_genlepvsb3->Sumw2();
  
  hist_2d_pt_genlepvsb4 = new TH2D("hist_2d_ptgenb_ptgenblep","pt of gen b vs pt of gen (b + lep)",120,20,2000,120,20,1000);
  hist_2d_pt_genlepvsb4->Sumw2();

  hist_2d_pt_genlepvsb5 = new TH2D("hist_2d_ptgenb_ptgentop","pt of gen b vs pt of gen top ",120,20,2000,120,20,1000);
  hist_2d_pt_genlepvsb5->Sumw2();

  hist_2d_pt_gentopvsgentop = new TH2D("hist_2d_ptgentop1_ptgentop2","pt of subleading gen top vs pt of leading gen top ",120,20,2000,120,20,2000);
  hist_2d_pt_gentopvsgentop->Sumw2();

  /* hist_delptbypt[0]= new TH2D("hist_delptbypt_genlep","del pt by pt between gen and reco lep",60,-50,50,60,30,800);
  hist_delptbypt[1] = new TH2D("hist_delptbypt_gen4","del pt by pt between gen b quark and reco AK4",60,-50,50,60,30,800);
  hist_delptbypt[2] = new TH2D("hist_delptbypt_genak8","del pt by pt between gen AK8 sys and reco AK8",60,-200,200,60,300,1000);  
  */
  hist_ortho_trig_2d[0] = new TH2D("hist_orthotrig_lepptvsleppt","hist_orthotrig_lepptvsleppt",nbins_ptlep,ptbins_lep,nbins_ptlep,ptbins_lep);
  hist_ortho_trig_2d[1] = new TH2D("hist_orthotrig_lepetavslepeta","hist_orthotrig_lepetavslepeta",5,-2.5,2.5,5,-2.5,2.5);
  hist_ortho_trig_2d[2] = new TH2D("hist_orthotrig_lepphivslepphi","hist_orthotrig_lepphivslepphi",5,-1*TMath::Pi(),TMath::Pi(),5,-1*TMath::Pi(),TMath::Pi());
  hist_ortho_trig_2d[3] = new TH2D("hist_orthotrig_ak8ptvsak8pt","hist_orthotrig_ak8ptvsak8pt",nbins_ptlep,ptbins_ak8,nbins_ptlep,ptbins_ak8);
  hist_ortho_trig[0] = new TH1D("hist_orthotrig_leadak8leppt","hist_orthotrig_leadak8leppt",30,0,900);
  hist_ortho_trig[1] = new TH1D("hist_orthotrig_dilepmass","hist_orthotrig_dilepmass",8,50,130);
  hist_ortho_trig[2] = new TH1D("hist_orthotrig_leadak8lepphi","hist_orthotrig_leadak8lepphi",5,-1*TMath::Pi(),TMath::Pi());
  hist_ortho_trig[3] = new TH1D("hist_orthotrig_leadak8lepeta","hist_orthotrig_leadak8lepeta",5,-2.5,2.5);
  hist_ortho_trig[4] = new TH1D("hist_orthotrig_subleadak8leppt","hist_orthotrig_subleadak8leppt",30,0,900);
  hist_ortho_trig[5] = new TH1D("hist_orthotrig_subleadak8lepphi","hist_orthotrig_subleadak8lepphi",5,-1*TMath::Pi(),TMath::Pi());
  hist_ortho_trig[6] = new TH1D("hist_orthotrig_subleadak8lepeta","hist_orthotrig_subleadak8lepeta",5,-2.5,2.5);
  hist_ortho_trig[7] = new TH1D("hist_orthotrig_leadak8pt","hist_orthotrig_leadak8pt",30,0,1500);
  hist_ortho_trig[8] = new TH1D("hist_orthotrig_leadak8phi","hist_orthotrig_leadak8phi",5,-1*TMath::Pi(),TMath::Pi());
  hist_ortho_trig[9] = new TH1D("hist_orthotrig_leadak8eta","hist_orthotrig_leadak8eta",5,-2.5,2.5);
  hist_ortho_trig[10] = new TH1D("hist_orthotrig_subleadak8pt","hist_orthotrig_subleadak8pt",30,0,1500);
  hist_ortho_trig[11] = new TH1D("hist_orthotrig_subleadak8phi","hist_orthotrig_subleadak8phi",5,-1*TMath::Pi(),TMath::Pi());
  hist_ortho_trig[12] = new TH1D("hist_orthotrig_subleadak8eta","hist_orthotrig_subleadak8eta",5,-2.5,2.5);
  
  hist_trig_eff_2d[0] = new TH2D("hist_trig_eff_lepptvsleppt","hist_trig_eff_lepptvsleppt",nbins_ptlep,ptbins_lep,nbins_ptlep,ptbins_lep);
  hist_trig_eff_2d[1] = new TH2D("hist_trig_eff_lepetavslepeta","hist_trig_eff_lepetavslepeta",5,-2.5,2.5,5,-2.5,2.5);
  hist_trig_eff_2d[2] = new TH2D("hist_trig_eff_lepphivslepphi","hist_trig_eff_lepphivslepphi",5,-1*TMath::Pi(),TMath::Pi(),5,-1*TMath::Pi(),TMath::Pi());
  hist_trig_eff_2d[3] = new TH2D("hist_trig_eff_ak8ptvsak8pt","hist_trig_eff_ak8ptvsak8pt",nbins_ptlep,ptbins_ak8,nbins_ptlep,ptbins_ak8);
  hist_trig_eff[0] = new TH1D("hist_trig_eff_leadleppt","hist_trig_eff_leadleppt",30,0,900);
  hist_trig_eff[1] = new TH1D("hist_trig_eff_dilepmass","hist_trig_eff_dilepmass",8,50,130);
  hist_trig_eff[2] = new TH1D("hist_trig_eff_leadlepphi","hist_trig_eff_leadlepphi",5,-1*TMath::Pi(),TMath::Pi());
  hist_trig_eff[3] = new TH1D("hist_trig_eff_leadlepeta","hist_trig_eff_leadlepeta",5,-2.5,2.5);
  hist_trig_eff[4] = new TH1D("hist_trig_eff_subleadak8leppt","hist_trig_eff_subleadak8leppt",30,0,900);
  hist_trig_eff[5] = new TH1D("hist_trig_eff_subleadak8lepphi","hist_trig_eff_subleadak8lepphi",5,-1*TMath::Pi(),TMath::Pi());
  hist_trig_eff[6] = new TH1D("hist_trig_eff_subleadak8eta","hist_trig_eff_subleadak8eta",5,-2.5,2.5);
  hist_trig_eff[7] = new TH1D("hist_trig_eff_leadak8pt","hist_trig_eff_leadak8pt",30,0,1500);
  hist_trig_eff[8] = new TH1D("hist_trig_eff_leadak8phi","hist_trig_eff_leadak8phi",5,-1*TMath::Pi(),TMath::Pi());
  hist_trig_eff[9] = new TH1D("hist_trig_eff_leadak8eta","hist_trig_eff_leadak8eta",5,-2.5,2.5);
  hist_trig_eff[10] = new TH1D("hist_trig_eff_subleadak8pt","hist_trig_eff_subleadak8pt",30,0,1500);
  hist_trig_eff[11] = new TH1D("hist_trig_eff_subleadak8phi","hist_trig_eff_subleadak8phi",5,-1*TMath::Pi(),TMath::Pi());
  hist_trig_eff[12] = new TH1D("hist_trig_eff_subleadak8eta","hist_trig_eff_subleadak8eta",5,-2.5,2.5);

  for(int ich=0; ich<3; ich++)
    {
      char lnamein[1000];
 
      for(int lvar=0;lvar<ncutflow_sf;lvar++)
	{
	  sprintf(lnamein,"hist_cutvar_%s_%s",prvar_name_sf[lvar][ich],channelname[ich].c_str());
	  hist_cutvar_sf[lvar][ich] = new TH1D(lnamein,lnamein,prvar_bins_sf[lvar],prvar_low_sf[lvar],prvar_high_sf[lvar]); // plot variable before putting cut on it to give an idea of cut value
	  hist_cutvar_sf[lvar][ich]->Sumw2();
	}
	    
      for(int lvar=0;lvar<ncutflow_sf+1;lvar++)
	{
	  sprintf(lnamein,"hist_cutflow_after_%s_cut_%s",to_string(lvar).c_str(),channelname[ich].c_str()); // plot same variable after each cut to understand how Data vs MC changes after each cut
	  hist_met_cutflow_sf[lvar][ich] = new TH1D(lnamein,lnamein,30,0,900);
	  hist_met_cutflow_sf[lvar][ich]->Sumw2();
	}
    }


  /*    hist_compare1[0] = new TH1D("no_ofelectrons_withoutcut","no_ofelectrons_withoutcut",30,0,900);
	hist_compare1[0]->Sumw2();
    
	hist_compare1[1] = new TH1D("no_ofmuons_withoutcut","no_ofmuons_withoutcut",30,0,900);
	hist_compare1[1]->Sumw2();*/
  /*  hist_compare1[0] = new TH1D("ak81_toptaggerscore_withoutcut","ak81toptaggerscore_withoutcut",30,-1,1);
      hist_compare1[0]->Sumw2();
   
      hist_compare1[1] = new TH1D("ak82_toptaggerscore","ak82_toptaggerscore",30,-1,1);
      hist_compare1[1]->Sumw2();
   
      hist_compare1[2] = new TH1D("min_ak8_toptaggerscore","min_ak8_toptaggerscore",30,-1,1);
      hist_compare1[2]->Sumw2();
   
      hist_compare1[3] = new TH1D("electron1pt","electron1pt",30,0,900);
      hist_compare1[3]->Sumw2();
   
      hist_compare1[4] = new TH1D("lepton1pt","lepton1pt",30,0,900);
      hist_compare1[4]->Sumw2();
   
      hist_compare1[5] = new TH1D("electron2pt","electron2pt",30,0,900);
      hist_compare1[5]->Sumw2();
   
      hist_compare1[6] = new TH1D("lepton2pt","lepton2pt",30,0,900);
      hist_compare1[6]->Sumw2();*/

  /* hist_compare1[0] = new TH1D("hist_elecvar_0_checktmp","hist_elecvar_0_checktmp",300,-2.5,2.5);
     hist_compare1[0]->Sumw2();

     hist_compare1[1] = new TH1D("hist_elecvar_1_checktmp","hist_elecvar_1_checktmp",300,0,900);
     hist_compare1[1]->Sumw2();

     hist_compare1[2] = new TH1D("hist_elecvar_2_checktmp","hist_elecvar_2_checktmp",300,0,0.05);
     hist_compare1[2]->Sumw2();

     hist_compare1[3] = new TH1D("hist_elecvar_3_checktmp","hist_elecvar_3_checktmp",300,0,0.05);
     hist_compare1[3]->Sumw2();

     hist_compare1[4] = new TH1D("hist_elecvar_4_checktmp","hist_elecvar_4_checktmp",300,0,1.5);
     hist_compare1[4]->Sumw2();

     hist_compare1[5] = new TH1D("hist_elecvar_5_checktmp","hist_elecvar_5_checktmp",300,0,0.04);
     hist_compare1[5]->Sumw2();

     hist_compare1[6] = new TH1D("hist_elecvar_6_checktmp","hist_elecvar_6_checktmp",300,0,0.2);
     hist_compare1[6]->Sumw2();

     hist_compare1[7] = new TH1D("hist_elecvar_7_checktmp","hist_elecvar_7_checktmp",300,0,4);
     hist_compare1[7]->Sumw2();

     hist_compare1[8] = new TH1D("hist_elecvar_8_checktmp","hist_elecvar_8_checktmp",300,0,30);
     hist_compare1[8]->Sumw2();

     hist_compare1[9] = new TH1D("hist_elecvar_9_checktmp","hist_elecvar_9_checktmp",300,0,6);
     hist_compare1[9]->Sumw2();

     hist_compare1[10] = new TH1D("hist_elecvar_10_checktmp","hist_elecvar_10_checktmp",300,0,1);
     hist_compare1[10]->Sumw2();

     hist_compare1[11] = new TH1D("hist_elecvar_11_checktmp","hist_elecvar_11_checktmp",300,0,10);
     hist_compare1[11]->Sumw2();

     hist_compare1[12] = new TH1D("hist_elecvar_12_checktmp","hist_elecvar_12_checktmp",300,0,5);
     hist_compare1[12]->Sumw2();

     hist_compare1[13] = new TH1D("hist_elecvar_13_checktmp","hist_elecvar_13_checktmp",300,0,30);
     hist_compare1[13]->Sumw2();

     hist_compare1[14] = new TH1D("hist_elecvar_14_checktmp","hist_elecvar_14_checktmp",300,0,10);
     hist_compare1[14]->Sumw2();

     hist_compare1[15] = new TH1D("hist_elecvar_15_checktmp","hist_elecvar_15_checktmp",300,-0.2,0.1);
     hist_compare1[15]->Sumw2();

     hist_compare1[16] = new TH1D("hist_elecvar_16_checktmp","hist_elecvar_16_checktmp",300,-0.1,0.1);
     hist_compare1[16]->Sumw2();

     hist_compare1[17] = new TH1D("hist_elecvar_17_checktmp","hist_elecvar_17_checktmp",300,0,0.2);
     hist_compare1[17]->Sumw2();

     hist_compare1[18] = new TH1D("hist_elecvar_18_checktmp","hist_elecvar_18_checktmp",300,0,100);
     hist_compare1[18]->Sumw2();

     hist_compare1[19] = new TH1D("hist_elecvar_19_checktmp","hist_elecvar_19_checktmp",300,0,100);
     hist_compare1[19]->Sumw2();

     hist_compare1[20] = new TH1D("hist_elecvar_20_checktmp","hist_elecvar_20_checktmp",300,0,50);
     hist_compare1[20]->Sumw2();

     hist_compare1[21] = new TH1D("hist_elecvar_21_checktmp","hist_elecvar_21_checktmp",300,-1,1);
     hist_compare1[21]->Sumw2();

     hist_compare1[22] = new TH1D("hist_elecvar_22_checktmp","hist_elecvar_22_checktmp",300,-0.05,0.05);
     hist_compare1[22]->Sumw2();

     hist_compare1[23] = new TH1D("hist_elecvar_23_checktmp","hist_elecvar_23_checktmp",300,-2,1);
     hist_compare1[23]->Sumw2();

     hist_compare1[24] = new TH1D("hist_elecvar_24_checktmp","hist_elecvar_24_checktmp",300,0,10);
     hist_compare1[24]->Sumw2();

     hist_compare1[25] = new TH1D("hist_elecvar_25_checktmp","hist_elecvar_25_checktmp",300,0,2);
     hist_compare1[25]->Sumw2();

     hist_compare1[26] = new TH1D("hist_elecvar_26_checktmp","hist_elecvar_26_checktmp",300,-1,1);
     hist_compare1[26]->Sumw2();

     hist_compare1[27] = new TH1D("hist_elecvar_27_checktmp","hist_elecvar_27_checktmp",300,0,1);
     hist_compare1[27]->Sumw2();

     hist_compare1[28] = new TH1D("hist_elecvar_28_checktmp","hist_elecvar_28_checktmp",300,0,1);
     hist_compare1[28]->Sumw2();

     hist_compare1[29] = new TH1D("hist_elecvar_29_checktmp","hist_elecvar_29_checktmp",300,0,1);
     hist_compare1[29]->Sumw2();

     hist_compare1[30] = new TH1D("hist_elecvar_30_checktmp","hist_elecvar_30_checktmp",300,0,1);
     hist_compare1[30]->Sumw2();

     hist_compare1[31] = new TH1D("hist_elecvar_31_checktmp","hist_elecvar_31_checktmp",300,-2.5,2.5);
     hist_compare1[31]->Sumw2();

     hist_compare1[32] = new TH1D("hist_elecvar_32_checktmp","hist_elecvar_32_checktmp",300,0,300);
     hist_compare1[32]->Sumw2();

     hist_compare1[33] = new TH1D("hist_elecvar_33_checktmp","hist_elecvar_33_checktmp",300,-1,1);
     hist_compare1[33]->Sumw2();

     hist_compare1[34] = new TH1D("hist_elecvar_34_checktmp","hist_elecvar_34_checktmp",300,-1,1);
     hist_compare1[34]->Sumw2();
  */

  for(int ih=0; ih < 6; ih++)
    {
      hist_ak8pt_spectrum[ih] = new TH1D(ak8pt_spectrum_checknames[ih],ak8pt_spectrum_checknames[ih],obs_nbins[0],obs_low[0],obs_up[0]);
      hist_ak8pt_spectrum[ih]->Sumw2();
    }
  
  reader1 = new TMVA::Reader( "BDTG_Re" );
  reader1->AddVariable("selpfjetAK8NHadF", &in_pfjetAK8NHadF);
  reader1->AddVariable("selpfjetAK8neunhadfrac", &in_pfjetAK8neunhadfrac);
  reader1->AddVariable("selpfjetAK8subhaddiff", &in_pfjetAK8subhaddiff);
  reader1->AddVariable("selpfjetAK8tau21", &in_pfjetAK8tau21);
  reader1->AddVariable("selpfjetAK8chrad", &in_pfjetAK8chrad);
  reader1->AddVariable("selpfjetAK8sdmass", &in_pfjetAK8sdmass);
  reader1->AddVariable("selpfjetAK8matchedeldxy_sv", &in_pfjetAK8matchedeldxy_sv);
  reader1->AddVariable("selpfjetAK8matchedelcleta", &in_pfjetAK8matchedelcleta);
  reader1->AddVariable("selpfjetAK8matchedelpt", &in_pfjetAK8matchedelpt);
  reader1->AddVariable("selpfjetAK8matchedelsigmaieta", &in_pfjetAK8matchedelsigmaieta);
  reader1->AddVariable("selpfjetAK8matchedelsigmaiphi", &in_pfjetAK8matchedelsigmaiphi);
  reader1->AddVariable("selpfjetAK8matchedelr9full", &in_pfjetAK8matchedelr9full);
  reader1->AddVariable("selpfjetAK8matchedelsupcl_etaw", &in_pfjetAK8matchedelsupcl_etaw);
  reader1->AddVariable("selpfjetAK8matchedelsupcl_phiw", &in_pfjetAK8matchedelsupcl_phiw);
  reader1->AddVariable("selpfjetAK8matchedelhcaloverecal", &in_pfjetAK8matchedelhcaloverecal);
  reader1->AddVariable("selpfjetAK8matchedelcloctftrkn", &in_pfjetAK8matchedelcloctftrkn);
  reader1->AddVariable("selpfjetAK8matchedelcloctftrkchi2", &in_pfjetAK8matchedelcloctftrkchi2);
  reader1->AddVariable("selpfjetAK8matchedele1x5bye5x5", &in_pfjetAK8matchedele1x5bye5x5);
  reader1->AddVariable("selpfjetAK8matchedelnormchi2", &in_pfjetAK8matchedelnormchi2);
  reader1->AddVariable("selpfjetAK8matchedelhitsmiss", &in_pfjetAK8matchedelhitsmiss);
  reader1->AddVariable("selpfjetAK8matchedeltrkmeasure", &in_pfjetAK8matchedeltrkmeasure);
  reader1->AddVariable("selpfjetAK8matchedelecloverpout", &in_pfjetAK8matchedelecloverpout);
  reader1->AddVariable("selpfjetAK8matchedelecaletrkmomentum", &in_pfjetAK8matchedelecaletrkmomentum);
  reader1->AddVariable("selpfjetAK8matchedeldeltaetacltrkcalo", &in_pfjetAK8matchedeldeltaetacltrkcalo);
  reader1->AddVariable("selpfjetAK8matchedelsupcl_preshvsrawe", &in_pfjetAK8matchedelsupcl_preshvsrawe);
  reader1->AddVariable("selpfjetAK8matchedelpfisolsumphet", &in_pfjetAK8matchedelpfisolsumphet);
  reader1->AddVariable("selpfjetAK8matchedelpfisolsumchhadpt", &in_pfjetAK8matchedelpfisolsumchhadpt);
  reader1->AddVariable("selpfjetAK8matchedelpfisolsumneuhadet", &in_pfjetAK8matchedelpfisolsumneuhadet);
  reader1->AddVariable("selpfjetAK8matchedeletain", &in_pfjetAK8matchedeletain);
  reader1->AddVariable("selpfjetAK8matchedelphiin", &in_pfjetAK8matchedelphiin);
  reader1->AddVariable("selpfjetAK8matchedelfbrem", &in_pfjetAK8matchedelfbrem);
  reader1->AddVariable("selpfjetAK8matchedeleoverp", &in_pfjetAK8matchedeleoverp);
  reader1->AddVariable("selpfjetAK8matchedelhovere", &in_pfjetAK8matchedelhovere);
  reader1->AddVariable("selpfjetAK8matchedelRho", &in_pfjetAK8matchedelRho);
  reader1->BookMVA("BDTG method", weightfile1);
    
  reader4 = new TMVA::Reader( "BDTG_Rmu" );
  reader4->AddVariable( "selpfjetAK8NHadF", &in_mupfjetAK8NHadF);
  reader4->AddVariable( "selpfjetAK8neunhadfrac", &in_mupfjetAK8neunhadfrac);
  reader4->AddVariable( "selpfjetAK8subhaddiff", &in_mupfjetAK8subhaddiff);
  reader4->AddVariable( "selpfjetAK8tau21", &in_mupfjetAK8tau21);
  reader4->AddVariable( "selpfjetAK8chrad", &in_mupfjetAK8chrad);
  reader4->AddVariable( "selpfjetAK8sdmass", &in_mupfjetAK8sdmass);

  
  reader4->AddVariable("selpfjetAK8matchedmuonchi", &in_pfjetAK8matchedmuonchi);
  reader4->AddVariable("selpfjetAK8matchedmuonposmatch", &in_pfjetAK8matchedmuonposmatch);
  reader4->AddVariable("selpfjetAK8matchedmuontrkink", &in_pfjetAK8matchedmuontrkink);
  reader4->AddVariable("selpfjetAK8matchedmuonsegcom", &in_pfjetAK8matchedmuonsegcom);
  reader4->AddVariable("selpfjetAK8matchedmuonhit", &in_pfjetAK8matchedmuonhit);
  reader4->AddVariable("selpfjetAK8matchedmuonmst", &in_pfjetAK8matchedmuonmst);
  reader4->AddVariable("selpfjetAK8matchedmuontrkvtx", &in_pfjetAK8matchedmuontrkvtx);
  reader4->AddVariable("selpfjetAK8matchedmuondz", &in_pfjetAK8matchedmuondz);
  reader4->AddVariable("selpfjetAK8matchedmuonpixhit", &in_pfjetAK8matchedmuonpixhit);
  reader4->AddVariable("selpfjetAK8matchedmuontrklay", &in_pfjetAK8matchedmuontrklay);
  reader4->AddVariable("selpfjetAK8matchedmuonvalfrac", &in_pfjetAK8matchedmuonvalfrac);
  reader4->AddVariable("selpfjetAK8muinsubptrat", &in_pfjetAK8muinsubptrat);
  reader4->AddVariable("selpfjetAK8muinsubmassrat", &in_pfjetAK8muinsubmassrat);
  reader4->AddVariable("selpfjetAK8muinsubinvmass", &in_pfjetAK8muinsubinvmass);
  reader4->AddVariable("selpfjetAK8muinsubIfarbyI0", &in_pfjetAK8muinsubIfarbyI0);
  reader4->AddVariable("selpfjetAK8muinsubInearbyI0", &in_pfjetAK8muinsubInearbyI0);
  reader4->BookMVA("BDTG method", weightfile4);

  reader2 = new TMVA::Reader( "BDTG_signalevents" );
  /*reader2->AddVariable("l1pt_nearjet", &l1pt_nearjet);
    reader2->AddVariable("deltaR_l1j1", &deltaR_l1j1);
    reader2->AddVariable("ak81deep_tvsqcd", &ak81deep_tvsqcd);
    reader2->AddVariable("no_ak4bjets", &no_ak4bjets);
    reader2->AddVariable("M_bl1", &M_bl1);
    reader2->AddVariable("j1_btag_sc", &j1_btag_sc);
    reader2->AddVariable("ak82deep_tvsqcd", &ak82deep_tvsqcd); 
    reader2->AddVariable("l2pt_nearjet", &l2pt_nearjet);
    reader2->AddVariable("j2_btag_sc", &j2_btag_sc);
    reader2->AddVariable("deltaR_l2j2", &deltaR_l2j2);
    reader2->AddVariable("M_bl2", &M_bl2);
    reader2->AddVariable("rat_metpt_eventHT", &rat_metpt_eventHT);  
    reader2->AddVariable("deltaR_l1l2", &deltaR_l1l2);
    reader2->AddVariable("no_ak4jets", &no_ak4jets);
    reader2->AddVariable("rat_l2pt_l1pt", &rat_l2pt_l1pt);
  */
  reader2->AddVariable("deltaR_l1j1",&deltaR_l1j1);
  reader2->AddVariable("j1_btag_sc",&j1_btag_sc);
  reader2->AddVariable("M_bl1",&M_bl1);
  reader2->AddVariable("deltaR_l2j2",&deltaR_l2j2);
  reader2->AddVariable("no_ak4jets",&no_ak4jets);
  reader2->AddVariable("j2_btag_sc",&j2_btag_sc);
  reader2->AddVariable("M_bl2",&M_bl2);
  reader2->AddVariable("rat_l2pt_l1pt",&rat_l2pt_l1pt);
  reader2->AddVariable("ak81deep_tvsqcd",&ak81deep_tvsqcd);
  reader2->AddVariable("l1pt_nearjet",&l1pt_nearjet); 
  reader2->AddVariable("no_ak4bjets",&no_ak4bjets);
  reader2->AddVariable("ak82deep_tvsqcd",&ak82deep_tvsqcd); 
  reader2->AddVariable("l2pt_nearjet",&l2pt_nearjet);
  reader2->AddVariable("rat_metpt_eventHT",&rat_metpt_eventHT);
  reader2->AddVariable("deltaR_l1l2",&deltaR_l1l2);
  
  /*  reader2->AddVariable("M_l1l2", &M_l1l2);
      reader2->AddVariable("mt_of_l2met", &mt_of_l2met);
      reader2->AddVariable("mt_of_l1met", &mt_of_l1met);
      reader2->AddVariable("rat_metpt_ak4pt", &rat_metpt_ak4pt);  
      reader2->AddVariable("dirglthrmin", &dirglthrmin);
      reader2->AddVariable("ak82deep_wvsqcd", &ak82deep_wvsqcd);
      reader2->AddVariable("ak81deep_wvsqcd", &ak81deep_wvsqcd);
      reader2->AddVariable("EventHT", &EventHT);
      reader2->AddVariable("deltaPhi_l1l2", &deltaPhi_l1l2);
      reader2->AddVariable("dirgltrthr", &dirgltrthr);
      reader2->AddVariable("rat_metpt_ak8pt", &rat_metpt_ak8pt);
      reader2->AddVariable("met_pt", &met_pt);
      reader2->AddVariable("deltaR_j1j2", &deltaR_j1j2);
      reader2->AddVariable("deltaR_l2j1", &deltaR_l2j1);
      reader2->AddVariable("delta_phil2_met", &delta_phil2_met);
      reader2->AddVariable("delta_phibl2_met", &delta_phibl2_met);
      reader2->AddVariable("delta_phil1_met", &delta_phil1_met);
      reader2->AddVariable("delta_phibl1_met", &delta_phibl1_met);
      reader2->AddVariable("rat_extra_ak4jpt_lpt", &rat_extra_ak4jpt_lpt);
      reader2->AddVariable("met_phi", &met_phi);
      reader2->AddVariable("delta_phibl1bl2", &delta_phibl1bl2);
      reader2->AddVariable("extra_ak4jdeepb", &extra_ak4jdeepb);
      reader2->AddVariable("ptsum_extra_ak4", &ptsum_extra_ak4);
      reader2->AddVariable("extra_ak4j", &extra_ak4j);
      reader2->AddVariable("extra_ak4pt", &extra_ak4pt);
      reader2->AddVariable("extra_ak4mass", &extra_ak4mass);
      reader2->AddVariable("extra_ak4jqgl", &extra_ak4jqgl);
      reader2->AddVariable("no_ak8jets", &no_ak8jets);
      reader2->AddVariable("response_ak81", &response_ak81);
      reader2->AddVariable("response_ak82", &response_ak82);
  */  
#ifdef E_MU_TTBar
  weightfile2 = dir + TString("TMVAClassification_BDTG_for_EMu_channel_weights.xml");
#elif defined(E_E_TTBar)
  weightfile2 = dir + TString("TMVAClassification_BDTG_for_EE_channel_weights.xml");
#elif defined(MU_MU_TTBar)
  weightfile2 = dir + TString("TMVAClassification_BDTG_for_MuMu_channel_weights.xml");
#endif
  
  //reader2->BookMVA("BDTG method", weightfile2);
}

Bool_t Anal_Leptop_PROOF::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either Anal_Leptop_PROOF::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  int icheck = GetEntry(entry);
  if(isMC){
    weight = event_weight;
  }else{
    weight = 1;
  }
  //Tout->Fill();
  int unfolding_set_nottoinclude = (int)gRandom->Integer(10) +1;
  
  TString str;
  //str = TString::Format("check 3 for evt %d ",ievt);          
  //if(gProofServ) gProofServ->SendAsynMessage(str);
  hist_prvar[0]->Fill(weight,weight);
 
  vector<GenParton> genpartons;
  vector<GenParton> LHEtops;
  vector<TopQuark> gentops;
  vector<LHEparticle> lheparticles;
  int cntt = 0; int nel=0; int nmu=0;
  tau_decay_mode = 0;
  
  if( isMC ){
    getPartons(genpartons);
    getLHEParticles(lheparticles);
    // Get GEN-level top quarks
    tt_decay_mode = -1;
    if(isTT || isST || isTTX) {
      getGENTops(gentops,genpartons); // after shower (get top quarks from its daughters) --> will tell details about the signature of ttbar events at GEN level
      getLHETops(LHEtops,genpartons); // before shower (get original top quarks which have decayed) --> will be usedto derive top pt reweighting
      tt_decay_mode = 0;
      for(int ip=0; ip <(int)lheparticles.size();ip++)
	{
	  if(abs(lheparticles[ip].pdgId) == 11)
	    { tt_decay_mode +=10; cntt++; nel++;}
	  else if(abs(lheparticles[ip].pdgId) == 13)
	    { tt_decay_mode +=100; cntt++; nmu++;}
	  else if(abs(lheparticles[ip].pdgId) == 15)
	    { tt_decay_mode +=1000; cntt++;}
	}
      tt_decay_mode += 2 - cntt;
      
      // top pt reweighting //
      float toppt_wt = 1;
      if(LHEtops.size()==2 && isTT){
	toppt_wt = SF_TOP(0.0615,0.0005,TMath::Min(float(500),float(LHEtops[0].pt)),TMath::Min(float(500),float(LHEtops[1].pt)));
	weight *= toppt_wt;
      }
    }
  }    
  // event selection starts
  /*  str = TString::Format("event start evt %d ",ievt);          
      if(gProofServ) gProofServ->SendAsynMessage(str);  */
  // Pile up rewighing //
  if(isMC && npu_vert_true>=0 && npu_vert_true<100){
    float *puweights = Get_PU_Weights(npu_vert_true);
    puWeight = puweights[0];
    puWeightup = puweights[1];
    puWeightdown = puweights[2];
    weight *= puWeight;
    weight *= prefiringweight;
  }
  else if(isMC)
    {
      weight = 0;
      str = TString::Format("npu_vert_true = %d weight = %f nprimi = %d is out of range for evt %d  ; ttree = %d",npu_vert_true,weight,nprimi,ievt,icheck);          
      if(gProofServ) gProofServ->SendAsynMessage(str);
      return kFALSE;
    }
  hist_npu_vert_true->Fill(npu_vert_true,weight);
  hist_npu_vert->Fill(npu_vert,weight);
  hist_prvar[1]->Fill(nprimi,weight);
  hist_count->Fill(0.0,weight);

  if(isTT && isMC){
#ifdef E_MU_TTBar
    if(nel == 1 && nmu == 1)  hist_count_top->Fill(0.0,weight);
    if(nel == 1 && nmu == 1 && gentops.size()==2 && (gentops[0].daughter[0].p4+gentops[0].daughter[2].p4).Pt() > 300 && delta2R(gentops[0].daughter[0].p4,gentops[0].daughter[2].p4)<1.6 && (gentops[1].daughter[0].p4+gentops[1].daughter[2].p4).Pt() > 300 && delta2R(gentops[1].daughter[0].p4,gentops[1].daughter[2].p4)<1.6)     hist_count_top->Fill(1.0,weight);
   
#elif defined(E_E_TTBar)
    if(nel == 2 && nmu == 0)     hist_count_top->Fill(0.0,weight);

    if(nel == 2 && nmu == 0 && gentops.size()==2 && (gentops[0].daughter[0].p4+gentops[0].daughter[2].p4).Pt() > 300 && delta2R(gentops[0].daughter[0].p4,gentops[0].daughter[2].p4)<1.6 && (gentops[1].daughter[0].p4+gentops[1].daughter[2].p4).Pt() > 300 && delta2R(gentops[1].daughter[0].p4,gentops[1].daughter[2].p4)<1.6)     hist_count_top->Fill(1.0,weight);
    
#elif defined(MU_MU_TTBar)
    if(nel == 0 && nmu == 2)  hist_count_top->Fill(0.0,weight);
    if(nel == 0 && nmu == 2 && gentops.size()==2 && (gentops[0].daughter[0].p4+gentops[0].daughter[2].p4).Pt() > 300 && delta2R(gentops[0].daughter[0].p4,gentops[0].daughter[2].p4)<1.6 && (gentops[1].daughter[0].p4+gentops[1].daughter[2].p4).Pt() > 300 && delta2R(gentops[1].daughter[0].p4,gentops[1].daughter[2].p4)<1.6)     hist_count_top->Fill(1.0,weight);
    
#endif
  }
  
  if(isTT && gentops.size()==2 && (gentops[0].daughter[0].p4+gentops[0].daughter[2].p4).Pt() > 300 && delta2R(gentops[0].daughter[0].p4,gentops[0].daughter[2].p4)<1.6 && (gentops[1].daughter[0].p4+gentops[1].daughter[2].p4).Pt() > 300 && delta2R(gentops[1].daughter[0].p4,gentops[1].daughter[2].p4)<1.6) hist_count_lep[1]->Fill(0.0,weight);  

  // hist_met_cutflow[0]->Fill(miset_PUPPI,weight);
  hist_count->Fill(1,weight);
  if(isTT && gentops.size()==2 && (gentops[0].daughter[0].p4+gentops[0].daughter[2].p4).Pt() > 300 && delta2R(gentops[0].daughter[0].p4,gentops[0].daughter[2].p4)<1.6 && (gentops[1].daughter[0].p4+gentops[1].daughter[2].p4).Pt() > 300 && delta2R(gentops[1].daughter[0].p4,gentops[1].daughter[2].p4)<1.6) hist_count_lep[1]->Fill(1,weight);
  hist_npv_nopuwt->Fill(nprimi,weight);
  
  //pileup weight factor needs to be updated with UL/
  
  //if(isMC){
  //if(npu_vert>=0 && npu_vert<100){
  //puWeight = pu_rat18[npu_vert];
  //puWeightUp = pu_rat18_up[npu_vert];
  //puWeightDown = pu_rat18_dn[npu_vert];
  //}
  //if(!isnan(puWeightUp) || fabs(puWeightUp)<1.e+6){
  //weight_puwup = weight*puWeightUp;                                           
               
  //}
  //if(!isnan(puWeightDown) || fabs(puWeightDown)<1.e+6){
  //weight_puwdown = weight*puWeightDown;                                        
  
  //}
  //if(!isnan(puWeight) || fabs(puWeight)<1.e+6){
  // weight *= puWeight;                                                         
             
  //}
  //}
  //hist_npv->Fill(nprimi,weight);

  // First demand that the event should at least one of the triggers considered in various topologies            
  // !CAUTION: Change/Remove this while deriving trigger efficiency 
    
  bool itrig_pass = false;
  bool itrig_onlysinglee = false;
  hist_prvar[2]->Fill(int(hlt_Mu37Ele27)*32+int(hlt_Mu27Ele37)*16+int(hlt_Mu37TkMu27)*8+ int(hlt_DoubleEle25)*4+int(hlt_Mu50)*2+int(hlt_Ele50_PFJet165), weight);

  itrig_pass = (hlt_Mu37Ele27||hlt_Mu27Ele37||hlt_Mu37TkMu27||hlt_DoubleEle25||hlt_Mu50||hlt_Ele50_PFJet165);
  itrig_onlysinglee = (!hlt_Mu37Ele27 && !hlt_Mu27Ele37 && !hlt_Mu37TkMu27 && !hlt_DoubleEle25 && !hlt_Mu50 && hlt_Ele50_PFJet165);
		
  //if(!itrig_pass) return kFALSE; //event should at least fire a dileptonic/single lepton trigger            
  // end of basic trigger criterion //
  if(isTT && gentops.size()==2 && (gentops[0].daughter[0].p4+gentops[0].daughter[2].p4).Pt() > 300 && delta2R(gentops[0].daughter[0].p4,gentops[0].daughter[2].p4)<1.6 && (gentops[1].daughter[0].p4+gentops[1].daughter[2].p4).Pt() > 300 && delta2R(gentops[1].daughter[0].p4,gentops[1].daughter[2].p4)<1.6) hist_count_lep[1]->Fill(2,weight);
  hist_count->Fill(2,weight);
  // hist_met_cutflow[1]->Fill(miset_PUPPI,weight);
  // Get generator-level particles //
       
  
  int topblep[2][3]={{-1,-1,-1},{-1,-1,-1}};  /// [0][0,1,2] = top,b,antilep  [1][0,1,2] = antitop,antitop,lep  

  if(isMC){
    for(int igen=0; igen<ngenparticles; igen++){
      if(abs(genpartmompdg[igen])==15  ) {tau_decay_mode = genpartpdg[igen];}
      if(genpartpdg[igen]==6  ) {topblep[0][0]=igen;}
      if(genpartpdg[igen]==-6 ) {topblep[1][0]=igen;}
      if( (genpartpdg[igen]==11||genpartpdg[igen]==13) && genpartmompdg[igen]==-24 ) { topblep[1][2]=igen; }  //lep index at 1 and b quark index at 0
      if(genpartpdg[igen]==-5 && genpartmompdg[igen]==-6) { topblep[1][1]=igen; }
      if( (genpartpdg[igen]==-11||genpartpdg[igen]==-13) && genpartmompdg[igen]==24 ) {  topblep[0][2]=igen; }
      if(genpartpdg[igen]==5 && genpartmompdg[igen]==6) {  topblep[0][1]=igen; }
    }

  }
  if(isMC && isTT){
    nleptop = nhadtop = 0;
    int leptop_id_daught[2];
    
    for(auto & top: gentops){
      if(abs(top.daughter[0].pdgId)==11 || abs(top.daughter[0].pdgId)==13 || abs(top.daughter[0].pdgId)==15){

	leptop_id_daught[nleptop] = abs(top.daughter[0].pdgId);
	nleptop++;
      }
      if(abs(top.daughter[0].pdgId)==1 || abs(top.daughter[0].pdgId)==3 || abs(top.daughter[0].pdgId)==5){
	nhadtop++;
      }
    }

    // tagging top structure of events //
    DiLeptt = SemiLeptt = Hadtt = EE = MUMU = EMU = EJets = MUJets = TauTau = ETau = MuTau = false;
    if(nleptop==2 && nhadtop==0) { DiLeptt = true; }
    if(nleptop==1 && nhadtop==1) { SemiLeptt = true; }
    if(nleptop==0 && nhadtop==2) { Hadtt = true; }
    
    if(DiLeptt && abs(leptop_id_daught[0])==11 && abs(leptop_id_daught[1])==11) { EE = true; }
    if(DiLeptt && abs(leptop_id_daught[0])==13 && abs(leptop_id_daught[1])==13) { MUMU = true; }
    if(DiLeptt && abs(leptop_id_daught[0])==15 && abs(leptop_id_daught[1])==15) { TauTau = true; }
    if(DiLeptt && ((abs(leptop_id_daught[0])==11 && abs(leptop_id_daught[1])==13) || (abs(leptop_id_daught[0])==13 && abs(leptop_id_daught[1])==11)) ) { EMU = true; }
    if(DiLeptt && ((abs(leptop_id_daught[0])==11 && abs(leptop_id_daught[1])==15) || (abs(leptop_id_daught[0])==15 && abs(leptop_id_daught[1])==11)) ) { ETau = true; }
    if(DiLeptt && ((abs(leptop_id_daught[0])==13 && abs(leptop_id_daught[1])==15) || (abs(leptop_id_daught[0])==15 && abs(leptop_id_daught[1])==13)) ) { MuTau = true; }
    
    if(SemiLeptt && abs(leptop_id_daught[0])==11) { EJets = true; }
    if(SemiLeptt && abs(leptop_id_daught[0])==13) { MUJets = true; }
    if(SemiLeptt && abs(leptop_id_daught[0])==15) { TAUJets = true; }
    
    bool boosted = false;
    //    boosted = (LHEtops.size()>1 && LHEtops[0].pt>300. && LHEtops[1].pt>300)    
    
#ifdef E_MU_TTBar
    
    if(!(DiLeptt && EMU/* && boosted*/))   return kFALSE;      //for signal EMU
    //if((DiLeptt && EMU/* && boosted*/)) return kFALSE; //for non-signal EMU TTbar
    
#elif defined(E_E_TTBar)
    
    if(!(DiLeptt && EE/* && boosted*/)) return kFALSE; //for signal EE
    //if((DiLeptt && EE/* && boosted*/)) return kFALSE; //for non-signal EE TTbar

#elif defined(MU_MU_TTBar)
    
    if(!(DiLeptt && MUMU/* && boosted*/)) return kFALSE; //for signal MUMU
    //if((DiLeptt && MUMU/* && boosted*/)) return kFALSE; //for non-signal MUMU TTbar
  
#endif
    
    //str = TString::Format("check event passed with ttbar decay mode = %f for evt %d ",tt_decay_mode,ievt);          
    //if(gProofServ) gProofServ->SendAsynMessage(str);   

    hist_count_top->Fill(2.0,weight);
   
    if(isTT && gentops.size()==2 && (gentops[0].daughter[0].p4+gentops[0].daughter[2].p4).Pt() > 300 && delta2R(gentops[0].daughter[0].p4,gentops[0].daughter[2].p4)<1.6 && (gentops[1].daughter[0].p4+gentops[1].daughter[2].p4).Pt() > 300 && delta2R(gentops[1].daughter[0].p4,gentops[1].daughter[2].p4)<1.6)    hist_count_top->Fill(3.0,weight);

    if(isTT){
      for(int it=0;it<2;it++)
	{
	  if(topblep[it][0]>=0 && topblep[it][1]>=0 && topblep[it][2]>=0 && genpartpt[topblep[it][0]]>300)
	    {
	      TLorentzVector gentopvec,genbvec,genlepvec,genblepsys;
	      gentopvec.SetPtEtaPhiM(genpartpt[topblep[it][0]],genparteta[topblep[it][0]],genpartphi[topblep[it][0]],genpartm[topblep[it][0]]);
	      genbvec.SetPtEtaPhiM(genpartpt[topblep[it][1]],genparteta[topblep[it][1]],genpartphi[topblep[it][1]],genpartm[topblep[it][1]]);
	      genlepvec.SetPtEtaPhiM(genpartpt[topblep[it][2]],genparteta[topblep[it][2]],genpartphi[topblep[it][2]],genpartm[topblep[it][2]]);
	      genblepsys = genbvec + genlepvec;
	      hist_2d_pt_genlepvsb1->Fill(genbvec.Pt(),genlepvec.Pt(),weight);
	      hist_2d_pt_genlepvsb3->Fill(gentopvec.Pt(),genlepvec.Pt(),weight);
	      hist_2d_pt_genlepvsb2->Fill(genblepsys.Pt(),genlepvec.Pt(),weight);
	      hist_2d_pt_genlepvsb5->Fill(gentopvec.Pt(),genbvec.Pt(),weight);
	      hist_2d_pt_genlepvsb4->Fill(genblepsys.Pt(),genbvec.Pt(),weight);
	    }
	}
    }
    
  } // if(isMC && isTT)
  //Get RECO-level objects //
  //First sorted out electron and muon

  //Some histograms before getting electrons with analysis criteria
  for(int ie=0; ie<nelecs; ie++) {
    hist_prptvar[3][min(ntcount-1,ie)]->Fill(min(float(2.99), max(float(-2.99), eleta[ie])), min(ybins[nybin]-1.0, max(ybins[0]+0.1, double(elpt[ie]))),weight);
  }
  //Here you get electrons with your criteria
  vector <Electron> velectrons;
  getelectrons(velectrons,30,absetacut);    // electron_pt_cut & absetacut defined in Proof.h
  
  
  //Some histograms before getting muons with analysis criteria
  for(int mu=0; mu<nmuons; mu++){
    hist_prptvar[4][min(ntcount-1,mu)]->Fill(min(float(2.99), max(float(-2.99), muoneta[mu])), min(ybins[nybin]-1.0, max(ybins[0]+0.1, double(muonpt[mu]))),weight);
  }

  //Here you get muons with your criteria
  vector <Muon> vmuons;
  getmuons(vmuons,30,absetacut);        

  
  //Make lepton collection from electrons & muons (using only common variables)
  vector <Lepton> vleptons;
  getLeptons(vleptons,vmuons,velectrons,30);
    
  //int nbjet_cut = -1; Need to use when btagwt will be considered in mc event 
  //float btagwt = 1.;
  //float btagwtup = 1.;
  //float btagwtdown = 1.;
  //float btag_eff = 1;
  //int nbjetAK4 = 0;
  //int nbjetAK4_lead = 0;

  //Some histograms before getting jets with analysis criteria
  for(int ijet=0; ijet<npfjetAK4; ijet++){
    for(int ie=0; ie<(int)velectrons.size(); ie++) {                         
      double delr = delta2R(pfjetAK4y[ijet],pfjetAK4phi[ijet], velectrons[ie].y, velectrons[ie].phi);               
      hist_prptangle[0]->Fill(pfjetAK4pt[ijet], delr,weight);
    } 
    for(int imu=0; imu<(int)vmuons.size(); imu++)                                
      {                                                                         
	double delr = delta2R(pfjetAK4y[ijet], pfjetAK4phi[ijet], vmuons[imu].y, vmuons[imu].phi);                  
	hist_prptangle[1]->Fill(pfjetAK4pt[ijet], delr,weight);
      }
  }

  //Here you get AK4 jets with your criteria                                     
  vector <AK4Jet> Jets;
  getAK4jets(Jets,AK4jet_pt_cut,absetacut,isMC);
  
  // Add b tag SF (if not in ntuple)//
  for(auto & jet: Jets){
    
    BTagEntry::JetFlavor btv_flav;
    if(abs(jet.hadronFlavour)==5){ btv_flav = BTagEntry::FLAV_B; }
    else if (abs(jet.hadronFlavour)==4){ btv_flav = BTagEntry::FLAV_C; }
    else { btv_flav = BTagEntry::FLAV_UDSG; }
		
    jet.btag_DeepFlav_med_SF = reader_deepflav_med.eval_auto_bounds("central",btv_flav,fabs(jet.eta),jet.pt); 
    jet.btag_DeepFlav_med_SF_up = reader_deepflav_med.eval_auto_bounds("up",btv_flav,fabs(jet.eta),jet.pt);
    jet.btag_DeepFlav_med_SF_dn = reader_deepflav_med.eval_auto_bounds("down",btv_flav,fabs(jet.eta),jet.pt);

    jet.btag_DeepFlav_loose_SF = reader_deepflav_loose.eval_auto_bounds("central",btv_flav,fabs(jet.eta),jet.pt); 
    jet.btag_DeepFlav_loose_SF_up = reader_deepflav_loose.eval_auto_bounds("up",btv_flav,fabs(jet.eta),jet.pt);
    jet.btag_DeepFlav_loose_SF_dn = reader_deepflav_loose.eval_auto_bounds("down",btv_flav,fabs(jet.eta),jet.pt);
  }

  //Get b-tagged jets from AK4 jets                                            
  vector <AK4Jet> BJets;
  for(auto & jet: Jets){
    if(isBJet(jet,deep_btag_cut)){
      //See similar histograms as did early for other objects
      //hist_prptvar[2][min(ntcount-1,ij)]->Fill(min(2.99, max(-2.99, jet.y)), min(ybins[nybin]-1.0, max(ybins[0]+0.1, jet.pt)),weight);
      BJets.push_back(jet);
    }
  }
  
  //Fill similar histograms for largejet as did early for other objects
  /*for(int ijet=0; ijet<npfjetAK8; ijet++){                                       
    if(!pfjetAK8jetID[ijet]) continue;                                           
    pfjetAK8pt[ijet] *= pfjetAK8JEC[ijet] ;                                      
    pfjetAK8mass[ijet] *= pfjetAK8JEC[ijet];
    if(isMC){                                                                    
    pfjetAK8pt[ijet] *= max(float(0.0),(1+pfjetAK8reso[ijet])) ;                               
    pfjetAK8mass[ijet] *= max(float(0.0),(1+pfjetAK8reso[ijet])) ;
    }                                                                            
    hist_prptvar[0][min(ntcount-1,ijet)]->Fill(min(float(2.99), max(float(-2.99), pfjetAK8y[ijet])), min(ybins[nybin]-1.0, max(ybins[0]+0.1, double(pfjetAK8pt[ijet]))),weight);
    }*/
  
  //Here you get AK8 jets with your criteria 
  vector <AK8Jet> LJets;
  getAK8jets(LJets,300,absetacut,isMC);

  if(isMC){
    // assign information from GEN-matching
    AssignGen(LJets,genpartons);
    if(isTT || isTTX || isST) TopAssignment_toJet(LJets,LHEtops,gentops);
  }
  //Get index of AK8 jet nearest to each lepton                                  
  if (LJets.size()>0) {
    /*    for(auto & lep: vleptons){
	  lep.AK8_neighbor_index = get_nearest_AK8Jet(LJets,lep.p4);
	  }*/
    
    for(unsigned ijet=0; ijet<Jets.size(); ijet++){
      double dr = delta2R(LJets[0].p4,Jets[ijet].p4);
      hist_prptangle[2]->Fill(LJets[0].pt,dr,weight);
    }
    if (LJets.size()>1) {
      for(unsigned ijet=0; ijet<Jets.size(); ijet++){
	double dr2 = delta2R(LJets[1].p4,Jets[ijet].p4);
	hist_prptangle[3]->Fill(LJets[1].pt,dr2,weight);
      }
    }
    //Get indices of nearest lepton, AK4 jet for each AK8 jet                    
    //for(auto & jet: LJets){
    

    for(int ijet=0;ijet<min(2,int(LJets.size()));ijet++){
      //Match lepton with AK8 jets
      LJets[ijet].match_lepton_index = get_nearest_lepton(vleptons,LJets[ijet].p4);
      if (LJets[ijet].match_lepton_index >=0 && ijet!=LJets[ijet].match_lepton_index && ijet<int(vleptons.size())) {
	Lepton tmplep = vleptons[LJets[ijet].match_lepton_index];
	vleptons.erase(vleptons.begin() + LJets[ijet].match_lepton_index);
	vleptons.insert(vleptons.begin()+ijet, tmplep);
      }

      // Matched highest b tag AK4 jets  with AK8 jets //                                         
      LJets[ijet].match_AK4_index = get_nearest_AK4(Jets,LJets[ijet].p4);
      if(LJets[ijet].match_AK4_index>=0 && LJets[ijet].match_AK4_index<int(Jets.size())){
        LJets[ijet].matchAK4deepb = Jets[LJets[ijet].match_AK4_index].btag_DeepFlav;	
      }
      if (LJets[ijet].match_AK4_index >=0 && ijet!=LJets[ijet].match_AK4_index && ijet<int(Jets.size())) {
	AK4Jet tmpjet = Jets[LJets[ijet].match_AK4_index];
	Jets.erase(Jets.begin() + LJets[ijet].match_AK4_index);
	Jets.insert(Jets.begin() + ijet, tmpjet);
      }
    }
  
    //Debarati : Keeping histograms of GMA to see the distributions before selection
    if (LJets.size()>0) {
      for (int jk=0; jk<(int)vleptons.size(); jk++) {
	int ityp = (abs(vleptons[jk].pdgId)==11) ? 4 : 5;
	double dr = delta2R(LJets[0].p4,vleptons[jk].p4);
	hist_prptangle[ityp]->Fill(LJets[0].pt, dr,weight);
	if (LJets.size()>1) {
	  int ityp = (abs(vleptons[jk].pdgId)==11) ? 6 : 7;
	  double dr = delta2R(LJets[1].p4, vleptons[jk].p4);
	  hist_prptangle[ityp]->Fill(LJets[1].pt, dr,weight);
	}
      }
      
      if( Jets.size()>0 && LJets.size()>0){
	j1_btag_sc_ptsort=Jets[0].btag_DeepFlav;
	hist_2d_deltaR_vsbtagsc[0]->Fill(j1_btag_sc_ptsort,delta2R(LJets[0].p4,Jets[0].p4),weight);
	hist_2d_pt_vsbtagsc[0]->Fill(Jets[0].pt,j1_btag_sc_ptsort,weight);
      }
      if( Jets.size()>1 && LJets.size()>1){
	j2_btag_sc_ptsort=Jets[1].btag_DeepFlav;
	hist_2d_deltaR_vsbtagsc[1]->Fill(j2_btag_sc_ptsort,delta2R(LJets[1].p4, Jets[1].p4),weight);
	hist_2d_pt_vsbtagsc[1]->Fill(Jets[1].pt,j2_btag_sc_ptsort,weight);
      }
    }
  }

  //if(isnan(weight) || weight>1.e+12) { weight = 0; }
  
  /******trigger object along with pdgid*****/
  std::vector<std::pair<int,TLorentzVector> > TrigRefObj;
  
  for (int tr=0; tr<ntrigobjs; tr++) {
    TLorentzVector trigobj;
    trigobj.SetPtEtaPhiM(trigobjpt[tr],trigobjeta[tr],trigobjphi[tr],trigobjmass[tr]);
    TrigRefObj.push_back(std::make_pair(trigobjpdgId[tr],trigobj));
  }
  
  // end of object selection //
  //// Event selection ////

  leptonsf_weight = leptonsf_weight_stat = leptonsf_weight_syst = 1;
  for(int ie=0; ie<(int)velectrons.size(); ie++) {
    float *sfvalues;
    if( LJets.size() > 0 && (delta2R(velectrons[ie].p4,LJets[min(0,int(LJets.size()-1))].p4)<0.7 || delta2R(velectrons[ie].p4,LJets[min(1,int(LJets.size()-1))].p4)<0.7))
      sfvalues = Electron_SF(velectrons[ie].pt,velectrons[ie].eta,"reco");
    else if(!velectrons[ie].looseid)
      sfvalues = Electron_SF(velectrons[ie].pt,velectrons[ie].eta,"reco");
    else
      sfvalues = Electron_SF(velectrons[ie].pt,velectrons[ie].eta,"Tight");
    leptonsf_weight *= sfvalues[0];
    leptonsf_weight_stat *= (sfvalues[0] + sqrt(sfvalues[1]*sfvalues[1] + sfvalues[2]*sfvalues[2]));  // like this for time being 
    leptonsf_weight_syst *= (sfvalues[0] + sqrt(sfvalues[3]*sfvalues[3] + sfvalues[4]*sfvalues[4] + sfvalues[5]*sfvalues[5] + sfvalues[6]*sfvalues[6]));  // like this for time being 
  }
  
  for(int im=0; im<(int)vmuons.size(); im++) {
    float *sfvalues;
    if( LJets.size() > 0 && (delta2R(vmuons[im].p4,LJets[min(0,int(LJets.size()-1))].p4)<0.7 || delta2R(vmuons[im].p4,LJets[min(1,int(LJets.size()-1))].p4)<0.7))
      sfvalues = Muon_SF(vmuons[im].pt,vmuons[im].eta,"reco");
    else if(!vmuons[im].looseid)
      sfvalues = Muon_SF(vmuons[im].pt,vmuons[im].eta,"reco");
    else
      sfvalues = Muon_SF(vmuons[im].pt,vmuons[im].eta,"Tight");
    leptonsf_weight *= sfvalues[0];
    leptonsf_weight_stat *= sfvalues[1];
    leptonsf_weight_syst *= sfvalues[2];
  }
  
  weight *= leptonsf_weight;
  hist_prvar[3]->Fill(vleptons.size(), weight);

  if (vleptons.size()<2)  return kFALSE;
  /*  str = TString::Format("check start evt %d ",ievt);          
      if(gProofServ) gProofServ->SendAsynMessage(str);*/
  
  //at least two leptons with pT > 30 GeV at this stage
  hist_count->Fill(3,weight);
    
  if(isTT && gentops.size()==2 && (gentops[0].daughter[0].p4+gentops[0].daughter[2].p4).Pt() > 300 && delta2R(gentops[0].daughter[0].p4,gentops[0].daughter[2].p4)<1.6 && (gentops[1].daughter[0].p4+gentops[1].daughter[2].p4).Pt() > 300 && delta2R(gentops[1].daughter[0].p4,gentops[1].daughter[2].p4)<1.6) hist_count_lep[1]->Fill(3,weight);
  hist_met_cutflow[2]->Fill(LJets.size()>0?LJets[0].pt:-100,weight);

  //// Condition on number of leptons is put here only. We will need at least two leptons for trigger matching
  
  bool emu_ch = false;
  bool mumu_ch = false;
  bool ee_ch = false;

      
  if (((abs(vleptons[0].pdgId)==11 && abs(vleptons[1].pdgId)==13) || (abs(vleptons[0].pdgId)==13 && abs(vleptons[1].pdgId)==11)) && vleptons[0].pdgId*vleptons[1].pdgId <0) emu_ch = true;
  else if (abs(vleptons[0].pdgId) == 13 && abs(vleptons[1].pdgId) == 13 && vleptons[0].pdgId*vleptons[1].pdgId <0) mumu_ch = true;
  else if (abs(vleptons[0].pdgId) == 11 && abs(vleptons[1].pdgId) == 11 && vleptons[0].pdgId*vleptons[1].pdgId <0) ee_ch =true;

 
  // Composition of two leading leptons should be correct in respective event categories //    
  bool correct_lepton_combo = false;
#ifdef E_MU_TTBar
  if(emu_ch) {correct_lepton_combo = true;}
#elif defined(E_E_TTBar)
  if(ee_ch) {correct_lepton_combo = true;}
#elif defined(MU_MU_TTBar)
  if(mumu_ch) {correct_lepton_combo = true;}
#endif
	
  // hist_count->Fill(4,weight);

  //Now comes individual channel the triggers consideration//
  
  vector<Single_Trigger> vsinglelep_trig;
  vector<Double_Trigger> vdoublelep_trig;
  // _pt_cuts <-- offline pt thresholds (should be finalized from the trigger efficiency curves                                                                                  
  // _pids <-- corresponding pdgId of the objects used in triggering
  
#ifdef E_MU_TTBar

  Double_Trigger dtrig;
  dtrig.double_hlts = hlt_Mu37Ele27;
  dtrig.double_pt_cuts = {37+3,27+3};
  dtrig.double_pids = {13,11};
  vdoublelep_trig.push_back(dtrig);
  
  dtrig.double_hlts = hlt_Mu27Ele37;
  dtrig.double_pt_cuts  = {37+3,27+3};
  dtrig.double_pids = {11,13};
  vdoublelep_trig.push_back(dtrig);

  Single_Trigger strig;
  strig.single_hlts = hlt_Mu50;
  strig.single_pt_cuts = 50+3;
  strig.single_pids = 13;
  strig.single_other_pt_cuts = -99;
  strig.single_other_pids = 0;
  vsinglelep_trig.push_back(strig);

  strig.single_hlts = hlt_Ele50_PFJet165;
  strig.single_pt_cuts = 50+3;
  strig.single_pids = 11;
  strig.single_other_pt_cuts = 165+15;
  strig.single_other_pids = 0;
  vsinglelep_trig.push_back(strig);

#elif defined(E_E_TTBar)

  Double_Trigger dtrig;
  dtrig.double_hlts = hlt_DoubleEle25;
  dtrig.double_pt_cuts = {25+15,25+5};
  dtrig.double_pids = {11,11};
  vdoublelep_trig.push_back(dtrig);

  Single_Trigger strig;
  strig.single_hlts = hlt_Ele50_PFJet165;
  strig.single_pt_cuts = (50+3);
  strig.single_pids = 11;
  strig.single_other_pt_cuts = 165+15;
  strig.single_other_pids = 0;
  vsinglelep_trig.push_back(strig);
  
#elif defined(MU_MU_TTBar)

  Double_Trigger dtrig;
  dtrig.double_hlts = hlt_Mu37TkMu27;
  dtrig.double_pt_cuts = {37+3,27+3};
  dtrig.double_pids = {13,13};
  vdoublelep_trig.push_back(dtrig);

  Single_Trigger strig;
  strig.single_hlts = hlt_Mu50;
  strig.single_pt_cuts = (50+3);
  strig.single_pids = 13;
  strig.single_other_pt_cuts = -99;
  strig.single_other_pids = 0;
  vsinglelep_trig.push_back(strig);
    
#endif

  // remember that the other object in single lepton triggers is jet by default

  bool anytrig_pass(false);
  if(vdoublelep_trig.size()>0){
    for(unsigned ihlt=0; ihlt<vdoublelep_trig.size(); ihlt++){
      if(vdoublelep_trig[ihlt].double_hlts){ anytrig_pass = true; break; }
    }
  }
  
  if(vsinglelep_trig.size()>0){
    for(unsigned ihlt=0; ihlt<vsinglelep_trig.size(); ihlt++){
      if(vsinglelep_trig[ihlt].single_hlts){ anytrig_pass = true; break; }
    }
  }
  
	
  bool trig_threshold_pass(false), trig_matching_pass(false);

  vector<TH1D*> hists;
  hists.push_back(hist_init[0]);
  hists.push_back(hist_init[1]);

  vector<TH2D*> hist2d_prptangle;
  hist2d_prptangle.push_back(hist_prptangle[8]);
  hist2d_prptangle.push_back(hist_prptangle[9]);
	 
  Match_trigger(vsinglelep_trig,vdoublelep_trig,		
                TrigRefObj,
                vleptons[0],vleptons[1],Jets,
                trig_threshold_pass,
                trig_matching_pass,
                hists,hist2d_prptangle
                );	   
  // end of trigger stuffs //
    
  
  bool is_additional_lepton = false; float additional_muonpt=-99; float additional_elecpt=-99;
  int extra_el=0;
  int extra_mu=0;

  for(int ie=0; ie<(int)velectrons.size(); ie++) {
    if(velectrons[ie].looseid && (int)LJets.size()>=2 && delta2R(velectrons[ie].p4,LJets[0].p4)>0.7 && delta2R(velectrons[ie].p4,LJets[1].p4)>0.7)
      extra_el++;
  }
  for(int im=0; im<(int)vmuons.size(); im++) {
    if(vmuons[im].looseid && (int)LJets.size()>=2 && delta2R(vmuons[im].p4,LJets[0].p4)>0.7 && delta2R(vmuons[im].p4,LJets[1].p4)>0.7)
      extra_mu++;
  }

#ifdef E_MU_TTBar
  is_additional_lepton = ( extra_mu >0 || extra_el >0);
  if(int(vmuons.size())>1)
    {
      additional_muonpt = vmuons[1].pt;
    }
  if(int(velectrons.size())>1)
    {
      additional_elecpt = velectrons[1].pt;
    }
#elif defined(E_E_TTBar)
  is_additional_lepton = ( extra_mu >0 || extra_el >0);
  if(int(vmuons.size())>0)
    {
      additional_muonpt = vmuons[0].pt;
    }
  if(int(velectrons.size())>2)
    {
      additional_elecpt = velectrons[2].pt;
    }
#elif defined(MU_MU_TTBar)
  is_additional_lepton = (extra_mu >0 || extra_el > 0);
  if(int(vmuons.size())>2)
    {
      additional_muonpt = vmuons[2].pt;
    }
  if(int(velectrons.size())>0)
    {
      additional_elecpt = velectrons[0].pt;
    }
#endif

  //some variable distributions before using cuts on them//                        
  //if (LJets.size()<2) return kFALSE;
  //hist_count->Fill(4,weight);

  vector <AK4GenJet> GenJets;
  if(isMC){
    getAK4genjets(GenJets,AK4jet_pt_cut,absetacut);
  }
  
  vector <AK8GenJet> GenLJets;
  if(isMC){
    getAK8genjets(GenLJets,270,absetacut);
  }

  // Assign electronic top tagger score //                                          
  
  if (LJets.size()>1 && vleptons.size()>1 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7 && delta2R(LJets[1].p4,vleptons[1].p4) < 0.7)
    {
      ReadTagger(LJets,vleptons,vmuons,velectrons,reader1,reader4);
      response_ak81 = LJets[0].re_tvsb;
      response_ak82 = LJets[1].re_tvsb;

      if (Jets.size()>0){
	ak41pt=Jets[0].pt;
	ak41mass=Jets[0].mass;
      }

      if (Jets.size()>1){
    
	ak42pt=Jets[1].pt;
	ak42mass=Jets[1].mass;
      }
  
#ifdef E_MU_TTBar
      if(abs(vleptons[0].pdgId) == 11){
	response_ak81 = LJets[0].re_tvsb;
	response_ak82 = LJets[1].re_tvsb;
      }
      else{
	response_ak81 = LJets[1].re_tvsb;
	response_ak82 = LJets[0].re_tvsb;
      }
#endif
    }
  //////////////////// Calculation for trigger SF  ////////////////////////////////////////////////////////////

  bool triggersf_selection = false;
  //bool jet_trigger = hlt_AK8PFHT800_TrimMass50 || hlt_AK8PFHT900_TrimMass50 || hlt_AK8PFJet400_TrimMass30 || hlt_AK8PFJet420_TrimMass30 || hlt_AK8PFJet550 || hlt_PFHT1050 || hlt_PFJet500;
 
  bool top_tagger_cut = false;
#ifdef E_MU_TTBar
  top_tagger_cut = response_ak81 > 0.8 && response_ak82 > 0.5;
#elif defined(E_E_TTBar)
  top_tagger_cut = min(response_ak81,response_ak82) > 0.748 && max(response_ak81,response_ak82) > 0.888;
#elif defined(MU_MU_TTBar)
  top_tagger_cut = min(response_ak81,response_ak82) > 0.33 && max(response_ak81,response_ak82) > 0.9;
#endif
  /*  
#ifdef E_MU_TTBar
  triggersf_selection = (int)velectrons.size() >0 && (int)vmuons.size() > 0 && (int)LJets.size() > 1 && LJets[0].pt > 300 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -143 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7 && delta2R(LJets[1].p4,vleptons[1].p4) < 0.7 && jet_trigger && (vleptons[0].p4+vleptons[1].p4).M()>100  && top_tagger_cut;
#elif defined(E_E_TTBar)
  triggersf_selection = (int)velectrons.size() > 1 && (int)LJets.size() > 1 && LJets[0].pt > 300 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -121 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7 && delta2R(LJets[1].p4,vleptons[1].p4) < 0.7 && jet_trigger && (vleptons[0].p4+vleptons[1].p4).M()>100  && top_tagger_cut;
#elif defined(MU_MU_TTBar)
  triggersf_selection = (int)vmuons.size() > 1 && (int)LJets.size() > 1 && LJets[0].pt > 300 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -169 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7 && delta2R(LJets[1].p4,vleptons[1].p4) < 0.7 && jet_trigger && (vleptons[0].p4+vleptons[1].p4).M()>100  && top_tagger_cut;
#endif
  
  
  if((int)vleptons.size()>1 && (int)velectrons.size()>0 && (int)LJets.size()>0)
    {

      hist_cutvar_sf[0][0]->Fill(jet_trigger,weight);

      if(jet_trigger)
	{
	  hist_met_cutflow_sf[1][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[1][0]->Fill((int)velectrons.size(),weight);

	}
      if(jet_trigger && (int)velectrons.size() > 1)
	{
	  hist_met_cutflow_sf[2][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[2][0]->Fill((int)LJets.size(),weight);	  
	}
      if(jet_trigger && (int)velectrons.size() > 1 && (int)LJets.size() > 1)
	{
	  hist_met_cutflow_sf[3][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[3][0]->Fill(LJets[1].pt,weight);
	}

      if(jet_trigger && (int)velectrons.size() > 1 && (int)LJets.size() > 1 && LJets[1].pt > 300)
	{
	  hist_met_cutflow_sf[4][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[4][0]->Fill(vleptons[0].pdgId*vleptons[1].pdgId,weight);
	}
      if(jet_trigger && (int)velectrons.size() > 1 && (int)LJets.size() > 1 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -121)
	{
	  hist_met_cutflow_sf[5][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[5][0]->Fill(delta2R(LJets[0].p4,vleptons[0].p4),weight);
	}
      if(jet_trigger && (int)velectrons.size() > 1 && (int)LJets.size() > 1 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -121 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7)
	{
	  hist_met_cutflow_sf[6][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[6][0]->Fill(delta2R(LJets[1].p4,vleptons[1].p4),weight);
	}

      if(jet_trigger && (int)velectrons.size() > 1 && (int)LJets.size() > 1 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -121 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7 && delta2R(LJets[1].p4,vleptons[1].p4) < 0.7)
	{
	  hist_met_cutflow_sf[7][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[7][0]->Fill((vleptons[0].p4+vleptons[1].p4).M(),weight);
	}
      if(jet_trigger && (int)velectrons.size() > 1 && (int)LJets.size() > 1 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -121 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7 && delta2R(LJets[1].p4,vleptons[1].p4) < 0.7 && (vleptons[0].p4+vleptons[1].p4).M()>100)
	{	  
	  hist_met_cutflow_sf[8][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[8][0]->Fill(min(response_ak81,response_ak82),weight);
	}

      if(jet_trigger && (int)velectrons.size() > 1 && (int)LJets.size() > 1 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -121 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7 && delta2R(LJets[1].p4,vleptons[1].p4) < 0.7 && (vleptons[0].p4+vleptons[1].p4).M()>100 && min(response_ak81,response_ak82) > 0.748)
	{
	  hist_met_cutflow_sf[9][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[9][0]->Fill(max(response_ak81,response_ak82),weight);
	}
      if(jet_trigger && (int)velectrons.size() > 1 && (int)LJets.size() > 1 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -121 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7 && delta2R(LJets[1].p4,vleptons[1].p4) < 0.7 && (vleptons[0].p4+vleptons[1].p4).M()>100 && min(response_ak81,response_ak82) > 0.748 && max(response_ak81,response_ak82) > 0.888)
	{
	  hist_met_cutflow_sf[10][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[10][0]->Fill(hlt_DoubleEle25,weight);
	}
      if(jet_trigger && (int)velectrons.size() > 1 && (int)LJets.size() > 1 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -121 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7 && delta2R(LJets[1].p4,vleptons[1].p4) < 0.7 && (vleptons[0].p4+vleptons[1].p4).M()>100 && min(response_ak81,response_ak82) > 0.748 && max(response_ak81,response_ak82) > 0.888 && hlt_DoubleEle25)
	{
	  hist_met_cutflow_sf[11][0]->Fill(velectrons[0].pt,weight);
	  hist_cutvar_sf[11][0]->Fill(trig_threshold_pass && trig_matching_pass,weight);
	}

      if(jet_trigger && (int)velectrons.size() > 1 && (int)LJets.size() > 1 && LJets[1].pt > 300 && vleptons[0].pdgId*vleptons[1].pdgId == -121 && delta2R(LJets[0].p4,vleptons[0].p4) < 0.7 && delta2R(LJets[1].p4,vleptons[1].p4) < 0.7 && (vleptons[0].p4+vleptons[1].p4).M()>100 && min(response_ak81,response_ak82) > 0.748 && max(response_ak81,response_ak82) > 0.888 && hlt_DoubleEle25 && trig_threshold_pass && trig_matching_pass)
	{
	  hist_met_cutflow_sf[12][0]->Fill(velectrons[0].pt,weight);
	}
    }
  
  if(triggersf_selection)
    {
      hist_ortho_trig_2d[0]->Fill(min(float(499.0),(float)vleptons[0].pt),min(float(499.0),(float)vleptons[1].pt), weight);
      hist_ortho_trig_2d[1]->Fill(vleptons[0].eta,vleptons[1].eta, weight);
      hist_ortho_trig_2d[2]->Fill(vleptons[0].phi,vleptons[1].phi, weight);
      hist_ortho_trig_2d[3]->Fill(LJets[0].pt,LJets[1].pt, weight);
      hist_ortho_trig[0]->Fill(vleptons[0].pt, weight);
      hist_ortho_trig[1]->Fill((vleptons[0].p4 + vleptons[1].p4).M(), weight);
      hist_ortho_trig[2]->Fill(vleptons[0].phi, weight);
      hist_ortho_trig[3]->Fill(vleptons[0].eta, weight);
      hist_ortho_trig[4]->Fill(vleptons[1].pt, weight);
      hist_ortho_trig[5]->Fill(vleptons[1].phi, weight);
      hist_ortho_trig[6]->Fill(vleptons[1].eta, weight);
      hist_ortho_trig[7]->Fill(LJets[0].pt, weight);
      hist_ortho_trig[8]->Fill(LJets[0].phi, weight);
      hist_ortho_trig[9]->Fill(LJets[0].eta, weight);
      hist_ortho_trig[10]->Fill(LJets[1].pt, weight);
      hist_ortho_trig[11]->Fill(LJets[1].phi, weight);
      hist_ortho_trig[12]->Fill(LJets[1].eta, weight);
    }
  if(triggersf_selection && anytrig_pass && trig_threshold_pass && trig_matching_pass)
    {
      hist_trig_eff_2d[0]->Fill(min(float(499.0),(float)vleptons[0].pt),min(float(499.0),(float)vleptons[1].pt), weight);
      hist_trig_eff_2d[1]->Fill(vleptons[0].eta,vleptons[1].eta, weight);
      hist_trig_eff_2d[2]->Fill(vleptons[0].phi,vleptons[1].phi, weight);
      hist_trig_eff_2d[3]->Fill(LJets[0].pt,LJets[1].pt, weight);
      hist_trig_eff[0]->Fill(vleptons[0].pt, weight);
      hist_trig_eff[1]->Fill((vleptons[0].p4 + vleptons[1].p4).M(), weight);
      hist_trig_eff[2]->Fill(vleptons[0].phi, weight);
      hist_trig_eff[3]->Fill(vleptons[0].eta, weight);
      hist_trig_eff[4]->Fill(vleptons[1].pt, weight);
      hist_trig_eff[5]->Fill(vleptons[1].phi, weight);
      hist_trig_eff[6]->Fill(vleptons[1].eta, weight);
      hist_trig_eff[7]->Fill(LJets[0].pt, weight);
      hist_trig_eff[8]->Fill(LJets[0].phi, weight);
      hist_trig_eff[9]->Fill(LJets[0].eta, weight);
      hist_trig_eff[10]->Fill(LJets[1].pt, weight);
      hist_trig_eff[11]->Fill(LJets[1].phi, weight);
      hist_trig_eff[12]->Fill(LJets[1].eta, weight);
    }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  */

  int genak8_index[2]={-10,-10};
  int genak4_index[2]={-10,-10};
  int genlep_index[2]={-10,-10};

  bool isgen_phasespace = false;
  bool isgen_phasespace_genak8 = false;

  TLorentzVector gentop_vec,gentop_bvec,gentop_lepvec,gentop_blep;
  TLorentzVector genantitop_vec,genantitop_bvec,genantitop_lepvec,genantitop_blep;

  if(isMC && isTT) {

    int ngenmatched_tops = 0;
    for(int itop = 0; itop<(int)GenLJets.size() ;itop++)
      {
	for(int ptop = 0; ptop < (int)gentops.size(); ptop++)
	  {
	    if(delta2R(GenLJets[itop].p4,gentops[ptop].p4) < 0.8 && delta2R(GenLJets[itop].p4,gentops[ptop].daughter[0].p4) < 0.8 && delta2R(GenLJets[itop].p4,gentops[ptop].daughter[2].p4) < 0.8)
	      {
		ngenmatched_tops++;
		GenLJets[itop].gentopmatched = true;
		GenLJets[itop].index_gentopmatched = ptop;
		break;
	      }
	  }
      }

    if((int)gentops.size() > 1)
      {
	int itop = -1; int iantitop = -1;
	if(gentops[0].daughter[2].pdgId == 5 && gentops[1].daughter[2].pdgId == -5 )
	  {itop = 0; iantitop = 1;}
	else if(gentops[0].daughter[2].pdgId ==  -5 && gentops[1].daughter[2].pdgId == 5)
	  {itop = 1; iantitop = 0;}
	else
	  {
	    str = TString::Format("Check Event no %d : Both top quarks have same pdg id of b quark",ievt);          
	    if(gProofServ) gProofServ->SendAsynMessage(str);
		return kFALSE;
	  }
#ifdef E_MU_TTBar
	if(abs(gentops[0].daughter[0].pdgId) == 11)
	  {itop = 0; iantitop = 1;}
	else if(abs(gentops[0].daughter[0].pdgId ==  -5))
	  {itop = 1; iantitop = 0;}
	else
	  {
	    str = TString::Format("Check Event no %d : Both top quarks have same pdg id of lepton",ievt);          
	    if(gProofServ) gProofServ->SendAsynMessage(str);
		return kFALSE;
	  }
#endif
	
	gentop_vec = gentops[itop].p4;
	gentop_bvec = gentops[itop].daughter[2].p4;
	gentop_lepvec = gentops[itop].daughter[0].p4;
	gentop_blep = gentop_bvec + gentop_lepvec;
	genantitop_vec = gentops[iantitop].p4;
	genantitop_bvec = gentops[iantitop].daughter[2].p4;
	genantitop_lepvec = gentops[iantitop].daughter[0].p4;
	genantitop_blep = genantitop_bvec + genantitop_lepvec;

	for(int il=0; il < (int)vleptons.size(); il++)
	  {
	    if(delta2R(vleptons[il].p4,gentop_lepvec) < 0.1)
	      hist_resolution_unfolding[0]->Fill(min(double(899),gentop_lepvec.Pt()), (vleptons[il].p4.Pt() - gentop_lepvec.Pt()) / gentop_lepvec.Pt(), weight);
	    if(abs(vleptons[il].p4.Eta() - gentop_lepvec.Eta()) < 0.01)
	      hist_resolution_unfolding[1]->Fill(gentop_lepvec.Phi(), (vleptons[il].p4.Phi() - gentop_lepvec.Phi()), weight);
	    if(abs(vleptons[il].p4.Phi() - gentop_lepvec.Phi()) < 0.01)
	      hist_resolution_unfolding[2]->Fill(gentop_lepvec.Eta(), (vleptons[il].p4.Eta() - gentop_lepvec.Eta()), weight);
	    
	    if(delta2R(vleptons[il].p4,genantitop_vec) < 0.1)
	      hist_resolution_unfolding[0]->Fill(min(double(899),genantitop_vec.Pt()), (vleptons[il].p4.Pt() - genantitop_vec.Pt()) / genantitop_vec.Pt(), weight);
	    if(abs(vleptons[il].p4.Eta() - genantitop_vec.Eta()) < 0.01)
	      hist_resolution_unfolding[1]->Fill(genantitop_vec.Phi(), (vleptons[il].p4.Phi() - genantitop_vec.Phi()), weight);
	    if(abs(vleptons[il].p4.Phi() - genantitop_vec.Phi()) < 0.01)
	      hist_resolution_unfolding[2]->Fill(genantitop_vec.Eta(), (vleptons[il].p4.Eta() - genantitop_vec.Eta()), weight);
	    hist_delr_genvsreco[0]->Fill(delta2R(vleptons[il].p4,gentop_lepvec),weight);
	    hist_delr_genvsreco[0]->Fill(delta2R(vleptons[il].p4,genantitop_lepvec),weight);
	  }
	
	for(int ij = 0; ij < (int)LJets.size(); ij++)
	  {
	    if(delta2R(LJets[ij].p4,gentop_blep) < 0.8)
	      hist_resolution_unfolding[3]->Fill(min(double(899),gentop_blep.Pt()), (LJets[ij].p4.Pt() - gentop_blep.Pt()) / gentop_blep.Pt(), weight);
	    if(abs(LJets[ij].p4.Eta() - gentop_blep.Eta()) < 0.1)
	      hist_resolution_unfolding[4]->Fill(gentop_blep.Phi(), (LJets[ij].p4.Phi() - gentop_blep.Phi()), weight);
	    if(abs(LJets[ij].p4.Phi() - gentop_blep.Phi()) < 0.1)
	      hist_resolution_unfolding[5]->Fill(gentop_blep.Eta(), (LJets[ij].p4.Eta() - gentop_blep.Eta()), weight);

	    if(delta2R(LJets[ij].p4,genantitop_blep) < 0.8)
	      hist_resolution_unfolding[3]->Fill(min(double(899),genantitop_blep.Pt()), (LJets[ij].p4.Pt() - genantitop_blep.Pt()) / genantitop_blep.Pt(), weight);
	    if(abs(LJets[ij].p4.Eta() - genantitop_blep.Eta()) < 0.1)
	      hist_resolution_unfolding[4]->Fill(genantitop_blep.Phi(), (LJets[ij].p4.Phi() - genantitop_blep.Phi()), weight);
	    if(abs(LJets[ij].p4.Phi() - genantitop_blep.Phi()) < 0.1)
	      hist_resolution_unfolding[5]->Fill(genantitop_blep.Eta(), (LJets[ij].p4.Eta() - genantitop_blep.Eta()), weight);
	      
	    hist_delr_genvsreco[1]->Fill(delta2R(LJets[ij].p4,gentop_blep),weight);
	    hist_delr_genvsreco[1]->Fill(delta2R(LJets[ij].p4,genantitop_blep),weight);
	    
	    /*if(delta2R(LJets[ij].p4,genantitop_vec) < 0.8)
	      {
	      hist_resolution_unfolding[6]->Fill(min(double(899),genantitop_vec.Pt()), (LJets[ij].p4.Pt() - genantitop_vec.Pt()) / genantitop_vec.Pt(), weight);
		hist_resolution_unfolding[7]->Fill(genantitop_vec.Phi(), (LJets[ij].p4.Phi() - genantitop_vec.Phi()) / genantitop_vec.Phi(), weight);
		hist_resolution_unfolding[8]->Fill(genantitop_vec.Eta(), (LJets[ij].p4.Eta() - genantitop_vec.Eta()) / genantitop_vec.Eta(), weight);
	      }
	    if(delta2R(LJets[ij].p4,gentop_vec) < 0.8)
	      {
		hist_resolution_unfolding[6]->Fill(min(double(899),gentop_vec.Pt()), (LJets[ij].p4.Pt() - gentop_vec.Pt()) / gentop_vec.Pt(), weight);
		hist_resolution_unfolding[7]->Fill(gentop_vec.Phi(), (LJets[ij].p4.Phi() - gentop_vec.Phi()) / gentop_vec.Phi(), weight);
		hist_resolution_unfolding[8]->Fill(gentop_vec.Eta(), (LJets[ij].p4.Eta() - gentop_vec.Eta()) / gentop_vec.Eta(), weight);
		}*/

	    for(int itop = 0; itop<(int)GenLJets.size() ;itop++)
	      {
		if(delta2R(GenLJets[itop].p4,LJets[ij].p4) < 0.8)
		  hist_resolution_unfolding[6]->Fill(min(double(899),GenLJets[itop].p4.Pt()), (LJets[ij].p4.Pt() - GenLJets[itop].p4.Pt()) / GenLJets[itop].p4.Pt(), weight);
		if(abs(GenLJets[itop].p4.Eta() - LJets[ij].p4.Eta()) < 0.1)
		  hist_resolution_unfolding[7]->Fill(GenLJets[itop].p4.Phi(), (LJets[ij].p4.Phi() - GenLJets[itop].p4.Phi()), weight);
		if(abs(GenLJets[itop].p4.Phi() - LJets[ij].p4.Phi()) < 0.1)
		  hist_resolution_unfolding[8]->Fill(GenLJets[itop].p4.Eta(), (LJets[ij].p4.Eta() - GenLJets[itop].p4.Eta()), weight);
		  
	      }
	    
	    hist_delr_genvsreco[2]->Fill(delta2R(LJets[ij].p4,gentop_vec),weight);
	    hist_delr_genvsreco[2]->Fill(delta2R(LJets[ij].p4,genantitop_vec),weight);

	  }

      }
    
    isgen_phasespace_genak8 = (int)GenLJets.size() > 1 ;
    isgen_phasespace = (int)gentops.size() > 1 && gentop_blep.Pt() > 270 && gentop_bvec.Pt() > 25 && gentop_lepvec.Pt() > 25 && genantitop_blep.Pt() > 270 && genantitop_bvec.Pt() > 25 && genantitop_lepvec.Pt() > 25 && delta2R(gentop_bvec,gentop_blep) < 0.8 && delta2R(gentop_lepvec,gentop_blep) < 0.8 && delta2R(genantitop_bvec,genantitop_blep) < 0.8 && delta2R(genantitop_lepvec,genantitop_blep) < 0.8 && abs(gentop_blep.Eta()) < 2.5 && abs(gentop_bvec.Eta()) < 2.5 && abs(gentop_lepvec.Eta()) < 2.5 && abs(genantitop_blep.Eta()) < 2.5 && abs(genantitop_bvec.Eta()) < 2.5 && abs(genantitop_lepvec.Eta()) < 2.5 && gentop_vec.Pt() > 270 && abs(gentop_vec.Eta()) < 2.5 && genantitop_vec.Pt() > 270 && abs(genantitop_vec.Eta()) < 2.5;
    //if(delta2R(GenLJets[itop].p4,LJets[ijet].p4) < dr)
     
    if(isgen_phasespace)
      {	
	bool is_reco = true;
	/*is_reco = (int)LJets.size() > 0 && LJets[0].re_tvsb > 0.9 && (int)Jets.size() > 0 && delta2R(Jets[0].p4,LJets[0].p4) < 0.8;
	if( (int)LJets.size() > 0 && delta2R(GenLJets[0].p4,LJets[0].p4) < 0.8  && is_reco)
	  {
	    hist_response_matrix[4]->Fill(min(float(999),LJets[0].pt),min(float(999),GenLJets[0].pt),weight);
	    hist_response_matrix_finebins[0]->Fill(LJets[0].pt,GenLJets[0].pt,weight);
	    }*/
	    
	if( !(correct_lepton_combo && anytrig_pass && trig_threshold_pass && trig_matching_pass && !is_additional_lepton && ((int)LJets.size()>=2 && LJets[1].pt>300.) && (int)Jets.size()>=2 && delta2R(LJets[0].p4, Jets[0].p4) < 0.8 && delta2R(LJets[0].p4,vleptons[0].p4)<0.7 && delta2R(LJets[1].p4,Jets[1].p4)<0.8 && delta2R(LJets[1].p4,vleptons[1].p4)<0.7 && (((vleptons[0].p4+vleptons[1].p4).M())>100) && top_tagger_cut) )
	  {
	    fill_unfolding_variable_withsets(hist_genlevl[0],min(double(799),gentop_lepvec.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_variable_withsets(hist_genlevl[1],min(double(799),genantitop_lepvec.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_variable_withsets(hist_genlevl[2],min(double(999),gentop_blep.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_variable_withsets(hist_genlevl[3],min(double(999),genantitop_blep.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_variable_withsets(hist_genlevl[4],min(double(999),gentop_vec.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_variable_withsets(hist_genlevl[5],min(double(999),genantitop_vec.Pt()),weight,unfolding_set_nottoinclude);

	    fill_unfolding_responsematrix_withsets(hist_response_matrix[0],-1,min(double(799),gentop_lepvec.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[1],-1,min(double(799),genantitop_lepvec.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[2],-1,min(double(999),gentop_blep.Pt()),weight,unfolding_set_nottoinclude);	    
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[3],-1,min(double(999),genantitop_blep.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[4],-1,min(double(999),gentop_vec.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[5],-1,min(double(999),genantitop_vec.Pt()),weight,unfolding_set_nottoinclude);

	    hist_response_matrix_finebins[0]->Fill(-1,min(double(999),(gentops[0].daughter[0].p4 + gentops[0].daughter[2].p4).Pt()),weight);
	    hist_response_matrix_finebins[1]->Fill(-1,min(double(999),gentop_vec.Pt()),weight);

	    if(isgen_phasespace_genak8)
	      {
		if((int)GenLJets.size() > 1  && GenLJets[0].gentopmatched && GenLJets[1].gentopmatched && GenLJets[0].index_gentopmatched < 0){
		  str = TString::Format("check evt %d top matched index for GenJet 0 is %d",ievt,GenLJets[0].index_gentopmatched);
		  if(gProofServ) gProofServ->SendAsynMessage(str);
		}

		hist_particlelevl[0]->Fill(min(double(999),GenLJets[0].p4.Pt()),weight);
		hist_particlelevl[1]->Fill(min(double(999),GenLJets[1].p4.Pt()),weight);

		hist_response_matrix_particlelvl[0]->Fill(-1,min(double(999),GenLJets[0].p4.Pt()),weight);
		hist_response_matrix_particlelvl[1]->Fill(-1,min(double(999),GenLJets[1].p4.Pt()),weight);
	      }
	  }
	
	/*is_reco = (int)LJets.size() > 1 && LJets[1].re_tvsb > 0.33 && (int)Jets.size() > 1 && delta2R(Jets[1].p4,LJets[1].p4) < 0.8;
	hist_genlevl[5]->Fill(min(float(999),GenLJets[1].pt),weight);
	if( (int)LJets.size() >1 && delta2R(GenLJets[1].p4,LJets[1].p4) < 0.8 && is_reco)
	  {
	    hist_response_matrix[5]->Fill(min(float(999),LJets[1].pt),min(float(999),GenLJets[1].pt),weight);
	    hist_response_matrix_finebins[1]->Fill(LJets[1].pt,GenLJets[1].pt,weight);
	    }
	if((int)LJets.size() < 2)
	  {
	    hist_response_matrix[5]->Fill(-1,min(float(999),GenLJets[1].pt),weight);
	    hist_response_matrix_finebins[1]->Fill(-1,GenLJets[1].pt,weight);
	    }*/
      }
    
    /* for(int itop = 0; itop<(int)gentops.size() ;itop++)
      {
	float dr = 0.4;
	for(int ij=0; ij < (int)GenJets.size(); ij++)
	  {
	    if(delta2R(gentops[itop].daughter[2].p4,GenJets[ij].p4) < dr)
	      {
		genak4_index[itop] = ij;
		dr = delta2R(gentops[itop].daughter[2].p4,GenJets[ij].p4);
	      }
	  }
	  }*/
    
    //////////////// For resolution   ///////////////////////////////

    /*for(int itop = 0; itop < (int)GenLJets.size(); itop++)
      {
	for(int ij = 0 ; ij < (int)LJets.size(); ij++)
	  {
	    if(delta2R(GenLJets[itop].p4,LJets[ij].p4) < 0.8)
	      {
		float genpt = GenLJets[itop].p4.Pt();
		hist_delptbypt[2]->Fill( ( genpt - LJets[ij].pt), LJets[ij].pt ,weight);
	      }
	  }
      }

    for(int ilep = 0; ilep < (int)genpartons.size(); ilep++)
      {
	for(int ij = 0 ; ij < (int)vleptons.size(); ij++)
	  {
	    if(vleptons[ij].pdgId == genpartons[ilep].pdgId && delta2R(vleptons[ij].p4, genpartons[ilep].p4) < 0.4)
	      hist_delptbypt[0]->Fill( (genpartons[ilep].pt - vleptons[ij].pt), vleptons[ij].pt ,weight);
	  }
      }
    
    for(int ib = 0; ib < (int)GenJets.size(); ib++)
      {
	for(int ij = 0 ; ij < (int)Jets.size(); ij++)
	  {
	    if(delta2R(GenJets[ib].p4, Jets[ij].p4) < 0.4)
	      hist_delptbypt[1]->Fill( (GenJets[ib].pt - Jets[ij].pt), Jets[ij].pt ,weight);
	  }
	  }*/
  }
        
  vector<Cuts> event_cuts;
  //Before this stage 3 event selection cuts applied offline => at least 1 primary vertex, at least one listed trigger fired, at least two leptons with pt > 30 GeV and they are closest (< 0.7) to leading AK8 jets, at least two large radius jets.
  Cuts sel_cut;
  
  //Here follows other event cuts offline//
  bool  histfill =true;
  
  if(histfill) {
    hist_prvar[4]->Fill(int(emu_ch)*4+int(mumu_ch)*2+int(ee_ch));
    hist_prvar[7]->Fill(vleptons[0].pdgId*vleptons[1].pdgId, weight);
  }
  // Composition of two leading leptons should be correct in respective event categories // \4.
  sel_cut.cut_pass = correct_lepton_combo;
  sel_cut.name = "Leptons must be ee/emu/mumu and should carry opposite charges";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;
    
  if(histfill) hist_prvar[5]->Fill(int(anytrig_pass), weight);
  // Did the event fire any trigger in the topology considered? \5.                     
  sel_cut.cut_pass = anytrig_pass;
  sel_cut.name = "Event should fire trigger depending upon working channel";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;
  
  // offline objects should pass kinematic thresholds in triggers \6.
  if ((int)vleptons.size()>0 && histfill)  hist_prvar[6]->Fill(min(float(899.0),vleptons[0].pt), weight);
  sel_cut.cut_pass = trig_threshold_pass;
  sel_cut.name = "Lepton in the event should have threshold pt";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;

  // offline objects should be matched to trigger objects   \7.
  if (histfill)  hist_prvar[6]->Fill(int(trig_matching_pass), weight);
  sel_cut.cut_pass = trig_matching_pass;
  sel_cut.name = "Offline lepton should match with trigger object";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;

  if(histfill){
    hist_prvar[8]->Fill(extra_el);
    hist_prvar[9]->Fill(extra_mu);
    hist_prvar[23]->Fill(min(float(299.0),additional_muonpt),weight);
    hist_prvar[24]->Fill(min(float(299.0),additional_elecpt),weight);
  }
  //Cut on the number of additional muons and electrons (pt>30 GeV) //   \8.                          
  sel_cut.cut_pass = !is_additional_lepton;
  sel_cut.name = "No extra e/mu in the event";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;
  
  if(histfill) {
    hist_prvar[10]->Fill(LJets.size(), weight);
    if(LJets.size()>=2) hist_prvar[25]->Fill(min(float(999.0),LJets[1].pt), weight);
  }
  //Event should have at least AK8 jet with pt > 300GeV \9.
  sel_cut.cut_pass = ((int)LJets.size()>=2 && LJets[1].pt>300.);
  sel_cut.name = "pt of both AK8 jet > 300 and no. of AK8 > 1";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;

  
  if(histfill) hist_prvar[11]->Fill(Jets.size(), weight);
  //Cut on at least two AK4 jets satisfying object selection criteria// \10.
  sel_cut.cut_pass = (int)Jets.size()>=2; 
  sel_cut.name = "No. of AK4 jets >= 2";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;
  
  double opangle = (LJets.size()>=1 && Jets.size()>=1) ? delta2R(LJets[0].p4, Jets[0].p4) : 10.0;
  if(histfill) hist_prvar[12]->Fill(opangle,weight);
  //Cut on AK8 and AK4 jet matching// \11.
  sel_cut.cut_pass = (opangle<0.8);
  sel_cut.name = "deltaR(1 AK4, 1 AK8) > 0.8";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;
  
  
  if(histfill)  hist_prvar[13]->Fill(delta2R(LJets[0].p4,vleptons[0].p4),weight);
  //Cut on AK8 and lepton matching// \12.
  sel_cut.cut_pass = (LJets.size()>=1 && vleptons.size()>=1 && delta2R(LJets[0].p4,vleptons[0].p4)<0.7);
  sel_cut.name = "deltaR(1 lep, 1 AK8) > 0.7";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;

  
  
  if(histfill) hist_prvar[14]->Fill(delta2R(LJets[1].p4,Jets[1].p4), weight);
  //Cut on 2nd AK8 and AK4 jet matching// \13.
  sel_cut.cut_pass = LJets.size()>=2 && Jets.size()>=2 && delta2R(LJets[1].p4,Jets[1].p4)<0.8;
  sel_cut.name = "deltaR(2 AK4, 2 AK8) > 0.8";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;


  if(histfill) hist_prvar[15]->Fill(delta2R(LJets[1].p4, vleptons[1].p4), weight);
  //Cut on 2nd AK8 and lepton matching// \14.
  sel_cut.cut_pass = LJets.size()>=2 && vleptons.size()>=2 && delta2R(LJets[1].p4,vleptons[1].p4)<0.7;
  sel_cut.name = "deltaR(2 lep, 2 AK8) > 0.7";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;
  

  /*  if(histfill) hist_prvar[16]->Fill(BJets.size(), weight);
  //Cut on btag score of AK4 jets though (as loose cut as possible)// \17.
  sel_cut.cut_pass = ;
  sel_cut.name = "";
  event_cuts.push_back(sel_cut);
  event_cuts.push_back(Jets.size()>=2 && Jets[0].btag_DeepFlav>0 && Jets[1].btag_DeepFlav>0);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;*/

  
  if(histfill){
    int ixtyp=-1;
    if (abs(vleptons[0].pdgId)==11 && abs(vleptons[1].pdgId)==11) {ixtyp=0;
    } else if (abs(vleptons[0].pdgId)==13 && abs(vleptons[1].pdgId)==13) {ixtyp=2;
    } else { ixtyp=1;}
    //str = TString::Format("event no %d ixtyp %i",ievt,ixtyp);          
    //if(gProofServ) gProofServ->SendAsynMessage(str);       
    hist_prvar[18+ixtyp]->Fill(min(float(899.0),float((vleptons[0].p4+vleptons[1].p4).M())), weight);
  }
  // Inv mass of selected leptons should be > 20 GeV \15.
  sel_cut.cut_pass = vleptons.size()>=2 && (((vleptons[0].p4+vleptons[1].p4).M())>100);
  sel_cut.name = "Inv mass of leptons > 100";
  event_cuts.push_back(sel_cut);
  histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;

  bool checknan = LJets.size()>=2 && !isnan(LJets[0].subhaddiff) && !isnan(LJets[1].subhaddiff);
  if(histfill) hist_prvar[26]->Fill(int(checknan), weight);
    
  // MET > 50 GeV \19.
  /*  if(histfill)   hist_prvar[21]->Fill(miset_PUPPI, weight);
      sel_cut.cut_pass = miset_PUPPI>50;
      sel_cut.name = "MET > 50";
      event_cuts.push_back(sel_cut);
      histfill *= event_cuts[int(event_cuts.size()-1)].cut_pass;*/
  //event_cuts array should be useful to derive (n-1) cut efficiency

  // To estimate no of events rejected by nan in large jet variables \19.
  
  bool event_pass = true;
  for(unsigned icut=0; icut<event_cuts.size(); icut++){
    event_pass *= event_cuts[icut].cut_pass;
    //  str = TString::Format("cut %u pass %o",icut+1,bool(event_cuts[icut]));          
    //  if(gProofServ) gProofServ->SendAsynMessage(str);                                
    if(!event_pass) break;
    if(event_pass){
      hist_count->Fill(4+icut,weight);
      if(isTT && gentops.size()==2 && (gentops[0].daughter[0].p4+gentops[0].daughter[2].p4).Pt() > 300 && delta2R(gentops[0].daughter[0].p4,gentops[0].daughter[2].p4)<1.6 && (gentops[1].daughter[0].p4+gentops[1].daughter[2].p4).Pt() > 300 && delta2R(gentops[1].daughter[0].p4,gentops[1].daughter[2].p4)<1.6) hist_count_lep[1]->Fill(4+icut,weight);
      hist_met_cutflow[3+icut]->Fill(LJets.size()>0?LJets[0].pt:-100,weight);
    }
  }
  //    str = TString::Format("check end evt %d ",ievt);          
  //	if(gProofServ) gProofServ->SendAsynMessage(str);
  if(!event_pass) return kFALSE;
  // end of event selection // 

  if(isTT && gentops.size()==2 && (gentops[0].daughter[0].p4+gentops[0].daughter[2].p4).Pt() > 300 && delta2R(gentops[0].daughter[0].p4,gentops[0].daughter[2].p4)<1.6 && (gentops[1].daughter[0].p4+gentops[1].daughter[2].p4).Pt() > 300 && delta2R(gentops[1].daughter[0].p4,gentops[1].daughter[2].p4)<1.6)    hist_count_top->Fill(4.0,weight);
  
  hist_prvar[22]->Fill(int(hlt_Mu37Ele27)*32+int(hlt_Mu27Ele37)*16+int(hlt_Mu37TkMu27)*8+ int(hlt_DoubleEle25)*4+int(hlt_Mu50)*2+int(hlt_Ele50_PFJet165), weight);

  hist_obs[23]->Fill(LJets[1].pt,weight);
  hist_obs[24]->Fill(LJets[1].y,weight);

  float trigger_sf_up=1,trigger_sf_down=1,top_tagger_sf_up=1,top_tagger_sf_down=1,trigger_sf=1,top_tagger_sf=1;
if(isMC){

  float *triggersfvalues;
#ifdef E_MU_TTBar
  triggersfvalues = Trigger_SF(LJets[0].pt,LJets[1].pt,"emu");
#elif defined(E_E_TTBar)
  triggersfvalues = Trigger_SF(LJets[0].pt,LJets[1].pt,"ee");
#elif defined(MU_MU_TTBar)
  triggersfvalues = Trigger_SF(LJets[0].pt,LJets[1].pt,"mumu");
#endif

  weight *= triggersfvalues[0];
  trigger_sf = triggersfvalues[0];
  trigger_sf_up = triggersfvalues[1];
  trigger_sf_down = triggersfvalues[2];
  }
  
  hist_obs[0]->Fill(LJets[0].pt,weight);
  hist_obs[1]->Fill(LJets[0].y,weight);

  hist_ak8pt_spectrum[0]->Fill(LJets[0].pt,weight);
  hist_ak8pt_spectrum[3]->Fill(LJets[1].pt,weight);

  if(top_tagger_cut)
    {
      hist_ak8pt_spectrum[1]->Fill(LJets[0].pt,weight);
      hist_ak8pt_spectrum[4]->Fill(LJets[1].pt,weight);
    }
  
  if(top_tagger_cut && isMC && (isTT || isTTX || isST) )
    {
#ifdef E_MU_TTBar
      int  ieltop = -1, imutop = -1;
      if(abs(vleptons[0].pdgId) == 11){
	ieltop = 0;
	imutop = 1;
      }
      else{
	ieltop = 1;
	imutop = 0;
      }
      float* factor1, *factor2;
      
      if(LJets[ieltop].hasleptop_alldecay)
	factor1 = ElectronToptagger_SF(LJets[ieltop].pt,"0p8");
      if(LJets[imutop].hasleptop_alldecay)
	factor2 = MuonToptagger_SF(LJets[imutop].pt,"0p5");

      if(LJets[imutop].hasleptop_alldecay && LJets[ieltop].hasleptop_alldecay)
	{
	  weight *= factor1[0]*factor2[0];
	  top_tagger_sf = factor1[0]*factor2[0];
	  top_tagger_sf_up = factor1[1]*factor2[1];
	  top_tagger_sf_down = factor1[2]*factor2[2];
	}
      else if(LJets[imutop].hasleptop_alldecay && !LJets[ieltop].hasleptop_alldecay)
	{
	  weight *= factor2[0];
	  top_tagger_sf = factor2[0];
	  top_tagger_sf_up = factor2[1];
	  top_tagger_sf_down = factor2[2];
	}
      else if(!LJets[imutop].hasleptop_alldecay && LJets[ieltop].hasleptop_alldecay)
	{
	  weight *= factor1[0];
	  top_tagger_sf = factor1[0];
	  top_tagger_sf_up = factor1[1];
	  top_tagger_sf_down = factor1[2];
	}
      
#elif defined(E_E_TTBar)
      int  itop1 = -1, itop2 = -1;
      if(LJets[0].re_tvsb > LJets[1].re_tvsb){
	itop1 = 0;
	itop2 = 1;
      }
      else{
	itop1 = 1;
	itop2 = 0;
      }      
      float *factor1, *factor2;
      if(LJets[itop1].hasleptop_alldecay)
	factor1 = ElectronToptagger_SF(LJets[itop1].pt,"0p888");
      if(LJets[itop2].hasleptop_alldecay)
	factor2 = ElectronToptagger_SF(LJets[itop2].pt,"0p748");
      //weight *= factor1*factor2;

      if(LJets[itop1].hasleptop_alldecay && LJets[itop2].hasleptop_alldecay)
	{
	  weight *= factor1[0]*factor2[0];
	  top_tagger_sf = factor1[0]*factor2[0];
	  top_tagger_sf_up = factor1[1]*factor2[1];
	  top_tagger_sf_down = factor1[2]*factor2[2];	  
	}
      else if(!LJets[itop1].hasleptop_alldecay && LJets[itop2].hasleptop_alldecay)
	{
	  weight *= factor2[0];
	  top_tagger_sf = factor2[0];
	  top_tagger_sf_up = factor2[1];
	  top_tagger_sf_down = factor2[2];	  
	}
      else if(LJets[itop1].hasleptop_alldecay && !LJets[itop2].hasleptop_alldecay)
	{
	  weight *= factor1[0];
	  top_tagger_sf = factor1[0];
	  top_tagger_sf_up = factor1[1];
	  top_tagger_sf_down = factor1[2];	  
	}
	    
#elif defined(MU_MU_TTBar)
      int  itop1 = -1, itop2 = -1;
      if(LJets[0].re_tvsb > LJets[1].re_tvsb){
	itop1 = 0;
	itop2 = 1;
      }
      else{
	itop1 = 1;
	itop2 = 0;
      }

      float *factor1, *factor2;

      if(LJets[itop1].hasleptop_alldecay)
	factor1 = MuonToptagger_SF(LJets[itop1].pt,"0p9");
      if(LJets[itop2].hasleptop_alldecay)
	factor2 = MuonToptagger_SF(LJets[itop2].pt,"0p33");

      if(LJets[itop1].hasleptop_alldecay && LJets[itop2].hasleptop_alldecay)
	{
	  weight *= factor1[0]*factor2[0];
	  top_tagger_sf = factor1[0]*factor2[0];
	  top_tagger_sf_up = factor1[1]*factor2[1];
	  top_tagger_sf_down = factor1[2]*factor2[2];	  
	}
      else if(!LJets[itop1].hasleptop_alldecay && LJets[itop2].hasleptop_alldecay)
	{
	  weight *= factor2[0];
	  top_tagger_sf = factor2[0];
	  top_tagger_sf_up = factor2[1];
	  top_tagger_sf_down = factor2[2];	  
	}
      else if(LJets[itop1].hasleptop_alldecay && !LJets[itop2].hasleptop_alldecay)
	{
	  weight *= factor1[0];
	  top_tagger_sf = factor1[0];
	  top_tagger_sf_up = factor1[1];
	  top_tagger_sf_down = factor1[2];	  
	}
     
      //weight *= factor1*factor2;
      //str = TString::Format("weight top tagger factor 1 = %f and factor 2 = %f",factor1,factor2);          
      //if(gProofServ) gProofServ->SendAsynMessage(str);                               
#endif
    }

  // fill basic distributions here: 
  if(LJets.size()>0){
	for(unsigned ijet=0; ijet<LJets.size(); ijet++){
	  hist_2D_msd_deepak8->Fill(LJets[ijet].sdmass,LJets[ijet].DeepTag_TvsQCD,weight);
	}//ijet
	if(Jets.size()>0){
	  for(unsigned ijet=0; ijet<Jets.size(); ijet++){
	    if(isBJet(Jets[ijet],deep_btag_cut)) {
	      if(abs(Jets[ijet].hadronFlavour)==5){hist_2D_bpass_flavb->Fill(Jets[ijet].pt,fabs(Jets[ijet].y),weight); }
	      else if(abs(Jets[ijet].hadronFlavour)==4){hist_2D_bpass_flavc->Fill(Jets[ijet].pt,fabs(Jets[ijet].y),weight); }
	      else{ hist_2D_bpass_flavq->Fill(Jets[ijet].pt,fabs(Jets[ijet].y),weight); }
	    }
	
	    if(abs(Jets[ijet].hadronFlavour)==5){ hist_2D_ball_flavb->Fill(Jets[ijet].pt,fabs(Jets[ijet].y),weight); }
	    else if(abs(Jets[ijet].hadronFlavour)==4){ hist_2D_ball_flavc->Fill(Jets[ijet].pt,fabs(Jets[ijet].y),weight); }
	    else{  hist_2D_ball_flavq->Fill(Jets[ijet].pt,fabs(Jets[ijet].y),weight); }
	  }
	}
      }
  // end of basic histograms 


  if(top_tagger_cut)
    {
      hist_ak8pt_spectrum[2]->Fill(LJets[0].pt,weight);
      hist_ak8pt_spectrum[5]->Fill(LJets[1].pt,weight);
      
      hist_init[2]->Fill(nprimi,weight);
      hist_init[3]->Fill(Jets.size(),weight);
      hist_init[4]->Fill(BJets.size(),weight);
      hist_init[5]->Fill(LJets.size(),weight);

      hist_obs[2]->Fill(LJets[0].mass,weight);
      hist_obs[3]->Fill(LJets[0].NHadF,weight);
      hist_obs[4]->Fill(LJets[0].neunhadfrac,weight);
      hist_obs[5]->Fill(LJets[0].sdmass,weight);
      hist_obs[6]->Fill(LJets[0].chrad,weight);
      hist_obs[7]->Fill(LJets[0].subhaddiff,weight);
      hist_obs[8]->Fill(LJets[0].tau21,weight);

      hist_obs[45]->Fill(LJets[0].tau32,weight);
      
  
      hist_2d_deltaR_vspt[0]->Fill(LJets[0].pt,delta2R(Jets[0].p4,vleptons[0].p4),weight);
      //TString str;                                                                 
      //  str = TString::Format("NHadF %f neunhadfrac %f subhaddiff %f subhaddiff %f",LJets[0].NHadF,LJets[0].neunhadfrac,LJets[0].subhaddiff,diff_func(LJets[0].sub1hadfrac,LJets[0].sub2hadfrac));
      //if(gProofServ) gProofServ->SendAsynMessage(str);

      hist_obs[9]->Fill(LJets[0].DeepTag_TvsQCD,weight);
      hist_obs[10]->Fill(LJets[0].DeepTag_WvsQCD,weight);
      hist_obs[11]->Fill(LJets[0].DeepTag_ZvsQCD,weight);
      hist_obs[12]->Fill(LJets[0].re_tvsb,weight);
      hist_obs[13]->Fill(LJets[0].re_tvsb,weight);
      hist_obs[14]->Fill(LJets[0].haspfelectron,weight);
      hist_obs[15]->Fill(LJets[0].haspfmuon,weight);
      hist_obs[16]->Fill(LJets[0].hasmatche,weight);
      hist_obs[17]->Fill(LJets[0].hasmatchmu,weight);
      hist_obs[18]->Fill(delta2R(LJets[0].p4,vleptons[0].p4),weight);
      hist_obs[19]->Fill(delta2R(LJets[0].p4,vleptons[1].p4),weight);


      if (Jets.size()>0) hist_obs[20]->Fill(delta2R(LJets[0].p4, Jets[0].p4), weight);
      if (Jets.size()>1)
	{
	  hist_obs[21]->Fill(delta2R(LJets[0].y, LJets[0].phi, Jets[1].y, Jets[1].phi), weight);
      
	  if (LJets.size()>1 ) {
	    hist_2d_deltaR_vspt[1]->Fill(LJets[1].pt,delta2R(Jets[1].y,Jets[1].phi,vleptons[1].y,vleptons[1].phi),weight);
	    hist_obs[22]->Fill(delta2R(LJets[0].y, LJets[0].phi, LJets[1].y, LJets[1].phi), weight);
	    hist_obs[25]->Fill(LJets[1].mass,weight);
	    hist_obs[26]->Fill(LJets[1].NHadF,weight);
	    hist_obs[27]->Fill(LJets[1].neunhadfrac,weight);
	    hist_obs[28]->Fill(LJets[1].sdmass,weight);
	    hist_obs[29]->Fill(LJets[1].chrad,weight);
	    hist_obs[30]->Fill(LJets[1].subhaddiff,weight);
	    hist_obs[31]->Fill(LJets[1].tau21,weight);
	    hist_obs[46]->Fill(LJets[1].tau32,weight);
     	        
	    hist_obs[32]->Fill(LJets[1].DeepTag_TvsQCD,weight);
	    hist_obs[33]->Fill(LJets[1].DeepTag_WvsQCD,weight);
	    hist_obs[34]->Fill(LJets[1].DeepTag_ZvsQCD,weight);
	    hist_obs[35]->Fill(LJets[1].re_tvsb,weight);
	    hist_obs[36]->Fill(LJets[1].re_tvsb,weight);
	 
	    hist_obs[37]->Fill(LJets[1].haspfelectron,weight);
	    hist_obs[38]->Fill(LJets[1].haspfmuon,weight);
	    hist_obs[39]->Fill(LJets[1].hasmatche,weight);
	    hist_obs[40]->Fill(LJets[1].hasmatchmu,weight);
	    hist_obs[41]->Fill(delta2R(LJets[1].p4,vleptons[0].p4),weight);
	    hist_obs[42]->Fill(delta2R(LJets[1].p4,vleptons[1].p4),weight);
	    hist_obs[43]->Fill(delta2R(LJets[1].p4, Jets[0].p4), weight);
	    hist_obs[44]->Fill(delta2R(LJets[1].p4, Jets[1].p4), weight);
	  }
	}
      if((int)vleptons.size()>1)
	{
	  	  /*hist_recolevl[0]->Fill(vleptons[0].pt,weight);
	  hist_recolevl[1]->Fill(vleptons[1].pt,weight);
	  if(genlep_index[0] > -0.5 )
	    {
	      if(delta2R(vleptons[0].p4,gentops[genlep_index[0]].daughter[0].p4) < 0.4){
		hist_genlevl[0]->Fill(gentops[genlep_index[0]].daughter[0].pt,weight);
		hist_response_matrix[0]->Fill(vleptons[0].p4.Pt(),gentops[genlep_index[0]].daughter[0].pt,weight);
	      }
	      else
		{
		  hist_genlevl[0]->Fill(gentops[genlep_index[0]].daughter[0].pt,weight);
		  hist_response_matrix[0]->Fill(-100,gentops[genlep_index[0]].daughter[0].pt,weight);	  
		}
	    }
	  else
	    {
	      hist_genlevl[0]->Fill(-100,weight);
	      hist_response_matrix[0]->Fill(vleptons[0].p4.Pt(),-100,weight);
	    }
	  
	  if(genlep_index[1] > -0.5)
	    {
	      if(delta2R(vleptons[1].p4,gentops[genlep_index[1]].daughter[0].p4) < 0.4){
		hist_genlevl[1]->Fill(gentops[genlep_index[1]].daughter[0].pt,weight);
		hist_response_matrix[1]->Fill(vleptons[1].p4.Pt(),gentops[genlep_index[1]].daughter[0].pt,weight);
	      }
	      else
		{
		  hist_genlevl[1]->Fill(gentops[genlep_index[1]].daughter[0].pt,weight);
		  hist_response_matrix[1]->Fill(-100,gentops[genlep_index[1]].daughter[0].pt,weight);	  
		}
	    }
	  else
	    {
	      hist_genlevl[1]->Fill(-100,weight);
	      hist_response_matrix[1]->Fill(vleptons[1].p4.Pt(),-100,weight);
	      }*/

	  hist_init[6]->Fill((vleptons[0].p4+vleptons[1].p4).M(),weight);
	  hist_init[7]->Fill(vleptons[0].p4.Pt(),weight);
	  hist_init[8]->Fill(vleptons[0].p4.Eta(),weight);
	  hist_init[9]->Fill(vleptons[0].p4.Phi(),weight);
	  
	  hist_init[10]->Fill(vleptons[1].p4.Pt(),weight);
	  hist_init[11]->Fill(vleptons[1].p4.Eta(),weight);
	  hist_init[12]->Fill(vleptons[1].p4.Phi(),weight);
	}
      if((int)Jets.size()>0) {
	hist_init[13]->Fill(Jets[0].pt, weight);
	hist_init[14]->Fill(Jets[0].eta, weight);
	hist_init[15]->Fill(Jets[0].phi, weight);
	/*hist_recolevl[2]->Fill(Jets[0].pt,weight);
	if(genak4_index[0] > -0.5)
	  {
	    if( delta2R(Jets[0].p4,GenJets[genak4_index[0]].p4) < 0.4)
	      {
		hist_genlevl[2]->Fill(GenJets[genak4_index[0]].pt,weight);
		hist_response_matrix[2]->Fill(Jets[0].pt,GenJets[genak4_index[0]].pt,weight);
	      }
	    else
	      {
		hist_genlevl[2]->Fill(GenJets[genak4_index[1]].pt,weight);
		hist_response_matrix[2]->Fill(-100,GenJets[genak4_index[1]].pt,weight);		
	      }
	  }
	else
	  {
	    hist_genlevl[2]->Fill(-100,weight);
	    hist_response_matrix[2]->Fill(Jets[0].pt,-100,weight);
	    }*/

      }
      if((int)Jets.size()>1) {
	hist_init[16]->Fill(Jets[1].pt, weight);
	hist_init[17]->Fill(Jets[1].eta, weight);
	hist_init[18]->Fill(Jets[1].phi, weight);

	/*hist_recolevl[3]->Fill(Jets[1].pt,weight);
	if(genak4_index[1] > -0.5)
	  {
	    if(delta2R(Jets[1].p4,GenJets[genak4_index[1]].p4) < 0.4)
	      {
		hist_genlevl[3]->Fill(GenJets[genak4_index[1]].pt,weight);
		hist_response_matrix[3]->Fill(Jets[1].pt,GenJets[genak4_index[1]].pt,weight);
	      }
	    else
	      {
		hist_genlevl[3]->Fill(GenJets[genak4_index[1]].pt,weight);
		hist_response_matrix[3]->Fill(-100,GenJets[genak4_index[1]].pt,weight);
	      }
	  }
	else
	  {
	    hist_genlevl[3]->Fill(-100,weight);
	    hist_response_matrix[3]->Fill(Jets[1].pt,-100,weight);
	    }*/

      }
      if((int)LJets.size()>0) {
	hist_init[19]->Fill(LJets[0].pt, weight);
	hist_init[20]->Fill(LJets[0].eta, weight);
	hist_init[21]->Fill(LJets[0].phi, weight);

	int itop_reco = -1;
	int iantitop_reco = -1;
	
	if(vleptons[0].pdgId < 0)
	  {
	    itop_reco = 0;
	    iantitop_reco = 1;
	  }
	else
	  {
	    itop_reco = 1;
	    iantitop_reco = 0;
	  }
#ifdef E_MU_TTBar
	if(abs(vleptons[0].pdgId) == 11)
	  {
	    itop_reco = 0;
	    iantitop_reco = 1;
	  }
	else
	  {
	    itop_reco = 1;
	    iantitop_reco = 0;
	  }
#endif
	fill_unfolding_variable_withsets(hist_recolevl[0],min(float(799),vleptons[itop_reco].pt),weight,unfolding_set_nottoinclude);
	fill_unfolding_variable_withsets(hist_recolevl[1],min(float(799),vleptons[iantitop_reco].pt),weight,unfolding_set_nottoinclude);
	fill_unfolding_variable_withsets(hist_recolevl[2],min(float(999),LJets[itop_reco].pt),weight,unfolding_set_nottoinclude);
	fill_unfolding_variable_withsets(hist_recolevl[3],min(float(999),LJets[iantitop_reco].pt),weight,unfolding_set_nottoinclude);

	float weightsysarray [7] = {1.0,puWeightup/puWeight,puWeightdown/puWeight,trigger_sf_up/trigger_sf,trigger_sf_down/trigger_sf,top_tagger_sf_up/top_tagger_sf,top_tagger_sf_down/top_tagger_sf};
	
	for(int igen = 0; igen<6; igen++)
	  {
	    for(int iwegsys = 0; iwegsys < 7; iwegsys++)
	      {
		//genantitop_lepvec.Pt()    gentop_blep.Pt()  genantitop_blep.Pt()  gentop_vec.Pt() genantitop_vec.Pt()
								 
		if(min(double(799),gentop_lepvec.Pt()) > unfold_recolepbins_finebins[2*igen] && min(double(799),gentop_lepvec.Pt()) < unfold_recolepbins_finebins[2*igen+2] && delta2R(gentop_lepvec,vleptons[itop_reco].p4) < 0.4)
		  hist_genlevl_binned[0][igen][iwegsys]->Fill(min(float(799),vleptons[itop_reco].pt),weight*weightsysarray[iwegsys]);
		if(min(double(799),genantitop_lepvec.Pt()) > unfold_recolepbins_finebins[2*igen] && min(double(799),genantitop_lepvec.Pt()) < unfold_recolepbins_finebins[2*igen+2] && delta2R(genantitop_lepvec,vleptons[iantitop_reco].p4) < 0.4)
		  hist_genlevl_binned[1][igen][iwegsys]->Fill(min(float(799),vleptons[iantitop_reco].pt),weight*weightsysarray[iwegsys]);

		if(min(double(999),gentop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
		  hist_genlevl_binned[2][igen][iwegsys]->Fill(min(float(999),LJets[itop_reco].pt),weight*weightsysarray[iwegsys]);
		if(min(double(999),genantitop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
		  hist_genlevl_binned[3][igen][iwegsys]->Fill(min(float(999),LJets[iantitop_reco].pt),weight);

		if(min(double(999),gentop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
		  hist_genlevl_binned[4][igen][iwegsys]->Fill(min(float(999),LJets[itop_reco].pt),weight*weightsysarray[iwegsys]);
		if(min(double(999),genantitop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
		  hist_genlevl_binned[5][igen][iwegsys]->Fill(min(float(999),LJets[iantitop_reco].pt),weight);

	      }
	    
		if(min(double(999),gentop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
		  hist_genlevl_binned[2][igen][7]->Fill(min(float(999), LJets[itop_reco].pt_resoup),weight);
		if(min(double(999),genantitop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
		  hist_genlevl_binned[3][igen][7]->Fill(min(float(999), LJets[iantitop_reco].pt_resoup),weight);

		if(min(double(999),gentop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
		  hist_genlevl_binned[4][igen][7]->Fill(min(float(999), LJets[itop_reco].pt_resoup),weight);
		if(min(double(999),genantitop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
		  hist_genlevl_binned[5][igen][7]->Fill(min(float(999), LJets[iantitop_reco].pt_resoup),weight);

		
		if(min(double(999),gentop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
		  hist_genlevl_binned[2][igen][8]->Fill(min(float(999), LJets[itop_reco].pt_resodn),weight);
		if(min(double(999),genantitop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2]  && delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
		  hist_genlevl_binned[3][igen][8]->Fill(min(float(999), LJets[iantitop_reco].pt_resodn),weight);

		if(min(double(999),gentop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
		  hist_genlevl_binned[4][igen][8]->Fill(min(float(999), LJets[itop_reco].pt_resodn),weight);
		if(min(double(999),genantitop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2]  && delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
		  hist_genlevl_binned[5][igen][8]->Fill(min(float(999), LJets[iantitop_reco].pt_resodn),weight);

		
		if(min(double(999),gentop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
		  hist_genlevl_binned[2][igen][9]->Fill(min(float(999), LJets[itop_reco].pt*LJets[itop_reco].jesup_total),weight);
		if(min(double(999),genantitop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2]  && delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
		  hist_genlevl_binned[3][igen][9]->Fill(min(float(999), LJets[iantitop_reco].pt*LJets[iantitop_reco].jesup_total),weight);

		if(min(double(999),gentop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
		  hist_genlevl_binned[4][igen][9]->Fill(min(float(999), LJets[itop_reco].pt*LJets[itop_reco].jesup_total),weight);
		if(min(double(999),genantitop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2]  && delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
		  hist_genlevl_binned[5][igen][9]->Fill(min(float(999), LJets[iantitop_reco].pt*LJets[iantitop_reco].jesup_total),weight);

		if(min(double(999),gentop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
		  hist_genlevl_binned[2][igen][10]->Fill(min(float(999), LJets[itop_reco].pt*LJets[itop_reco].jesdn_total),weight);
		if(min(double(999),genantitop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2]  && delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
		  hist_genlevl_binned[3][igen][10]->Fill(min(float(999), LJets[iantitop_reco].pt*LJets[iantitop_reco].jesdn_total),weight);

		if(min(double(999),gentop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2] && delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
		  hist_genlevl_binned[4][igen][10]->Fill(min(float(999), LJets[itop_reco].pt*LJets[itop_reco].jesdn_total),weight);
		if(min(double(999),genantitop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2]  && delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
		  hist_genlevl_binned[5][igen][10]->Fill(min(float(999), LJets[iantitop_reco].pt*LJets[iantitop_reco].jesdn_total),weight);
		
		/*if(min(double(999),gentop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2])  hist_genlevl_binned[2][igen][9]->Fill(min(float(999), LJets[itop_reco].pt*LJets[itop_reco].jesup_total),weight);
		if(min(double(999),genantitop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2])  hist_genlevl_binned[3][igen][9]->Fill(min(float(999), LJets[iantitop_reco].pt*LJets[iantitop_reco].jesup_total),weight);

		if(min(double(999),gentop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2])  hist_genlevl_binned[4][igen][9]->Fill(min(float(999),LJets[itop_reco].jesup_total),weight);
		if(min(double(999),genantitop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2])  hist_genlevl_binned[5][igen][9]->Fill(min(float(999), LJets[iantitop_reco].jesup_total),weight);

		if(min(double(999),gentop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2])  hist_genlevl_binned[2][igen][10]->Fill(min(float(999), LJets[itop_reco].jesdn_total),weight);
		if(min(double(999),genantitop_blep.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_blep.Pt()) < unfold_recoak8bins_finebins[2*igen+2])  hist_genlevl_binned[3][igen][10]->Fill(min(float(999), LJets[iantitop_reco].jesdn_total),weight);

		if(min(double(999),gentop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),gentop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2])  hist_genlevl_binned[4][igen][10]->Fill(min(float(999), LJets[itop_reco].jesdn_total),weight);
		if(min(double(999),genantitop_vec.Pt()) > unfold_recoak8bins_finebins[2*igen] && min(double(999),genantitop_vec.Pt()) < unfold_recoak8bins_finebins[2*igen+2])  hist_genlevl_binned[5][igen][10]->Fill(min(float(999),LJets[iantitop_reco].jesdn_total),weight);
		*/
	  }
	
	if(isgen_phasespace)
	  {
	    fill_unfolding_variable_withsets(hist_genlevl[0],min(double(799),gentop_lepvec.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_variable_withsets(hist_genlevl[1],min(double(799),genantitop_lepvec.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_variable_withsets(hist_genlevl[2],min(double(999),gentop_blep.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_variable_withsets(hist_genlevl[3],min(double(999),genantitop_blep.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_variable_withsets(hist_genlevl[4],min(double(999),gentop_vec.Pt()),weight,unfolding_set_nottoinclude);
	    fill_unfolding_variable_withsets(hist_genlevl[5],min(double(999),genantitop_vec.Pt()),weight,unfolding_set_nottoinclude);
	  
	    if(delta2R(gentop_lepvec,vleptons[itop_reco].p4) < 0.4)
	      fill_unfolding_responsematrix_withsets(hist_response_matrix[0],min(float(799),vleptons[itop_reco].pt),min(double(799),gentop_lepvec.Pt()),weight,unfolding_set_nottoinclude);
	    else
	      {
		fill_unfolding_responsematrix_withsets(hist_response_matrix[0],-1,min(double(799),gentop_lepvec.Pt()),weight,unfolding_set_nottoinclude);
		fill_unfolding_responsematrix_withsets(hist_response_matrix[0],min(float(799),vleptons[itop_reco].pt),-1,weight,unfolding_set_nottoinclude);
	      }
	    
	    if(delta2R(genantitop_lepvec,vleptons[iantitop_reco].p4) < 0.4)
	      fill_unfolding_responsematrix_withsets(hist_response_matrix[1],min(float(799),vleptons[iantitop_reco].pt),min(double(799),genantitop_lepvec.Pt()),weight,unfolding_set_nottoinclude);
	    else
	      {
		fill_unfolding_responsematrix_withsets(hist_response_matrix[1],-1,min(double(799),genantitop_lepvec.Pt()),weight,unfolding_set_nottoinclude);
		fill_unfolding_responsematrix_withsets(hist_response_matrix[1],min(float(799),vleptons[iantitop_reco].pt),-1,weight,unfolding_set_nottoinclude);		
	      }
	  
	    if(delta2R(gentop_blep,LJets[itop_reco].p4) < 0.8)
	      {
		fill_unfolding_responsematrix_withsets(hist_response_matrix[2],min(float(999),LJets[itop_reco].pt),min(double(999),gentop_blep.Pt()),weight,unfolding_set_nottoinclude);
		//fill_unfolding_responsematrix_withsets(hist_response_matrix_finebins[0],min(float(999),LJets[itop_reco].pt),min(double(999),gentop_blep.Pt()),weight,unfolding_set_nottoinclude);
	      }
	    else
	      {
		fill_unfolding_responsematrix_withsets(hist_response_matrix[2],-1,min(double(999),gentop_blep.Pt()),weight,unfolding_set_nottoinclude);
		fill_unfolding_responsematrix_withsets(hist_response_matrix[2],min(float(999),LJets[itop_reco].pt),-1,weight,unfolding_set_nottoinclude);

		hist_response_matrix_finebins[itop_reco]->Fill(-1,min(double(999),gentop_blep.Pt()),weight);
		hist_response_matrix_finebins[itop_reco]->Fill(min(float(999),LJets[itop_reco].pt),-1,weight);
	      }
	  
	    if(delta2R(genantitop_blep,LJets[iantitop_reco].p4) < 0.8)
	      fill_unfolding_responsematrix_withsets(hist_response_matrix[3],min(float(999),LJets[iantitop_reco].pt),min(double(999),genantitop_blep.Pt()),weight,unfolding_set_nottoinclude);
	    else
	      {
		fill_unfolding_responsematrix_withsets(hist_response_matrix[3],-1,min(double(999),genantitop_blep.Pt()),weight,unfolding_set_nottoinclude);
		fill_unfolding_responsematrix_withsets(hist_response_matrix[3],min(float(999),LJets[iantitop_reco].pt),-1,weight,unfolding_set_nottoinclude);
	      }
	  	  
	    if(delta2R(gentop_vec,LJets[itop_reco].p4) < 0.8)
	      {
		fill_unfolding_responsematrix_withsets(hist_response_matrix[4],min(float(999),LJets[itop_reco].pt),min(double(999),gentop_vec.Pt()),weight,unfolding_set_nottoinclude);
		//fill_unfolding_responsematrix_withsets(hist_response_matrix_finebins[iantitop_reco],min(float(999),LJets[itop_reco].pt),min(double(999),gentop_vec.Pt()),weight,unfolding_set_nottoinclude);
	      }
	    else
	      {
		fill_unfolding_responsematrix_withsets(hist_response_matrix[4],-1,min(double(999),gentop_vec.Pt()),weight,unfolding_set_nottoinclude);
		fill_unfolding_responsematrix_withsets(hist_response_matrix[4],min(float(999),LJets[itop_reco].pt),-1,weight,unfolding_set_nottoinclude);

		hist_response_matrix_finebins[1]->Fill(min(float(999),LJets[itop_reco].pt),-1,weight);
		hist_response_matrix_finebins[1]->Fill(-1,min(double(999),gentop_vec.Pt()),weight);
	      } 
	  	
	    if(delta2R(genantitop_vec,LJets[iantitop_reco].p4) < 0.8)
	      fill_unfolding_responsematrix_withsets(hist_response_matrix[5],min(float(999),LJets[iantitop_reco].pt),min(double(999),genantitop_vec.Pt()),weight,unfolding_set_nottoinclude);
	    else
	      {
		fill_unfolding_responsematrix_withsets(hist_response_matrix[5],-1,min(double(999),genantitop_vec.Pt()),weight,unfolding_set_nottoinclude);
		fill_unfolding_responsematrix_withsets(hist_response_matrix[5],min(float(999),LJets[iantitop_reco].pt),-1,weight,unfolding_set_nottoinclude);
	      }  
	  }
	else
	  {
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[0],min(float(799),vleptons[itop_reco].pt),-1,weight,unfolding_set_nottoinclude);
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[1],min(float(799),vleptons[iantitop_reco].pt),-1,weight,unfolding_set_nottoinclude);		
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[2],min(float(999),LJets[itop_reco].pt),-1,weight,unfolding_set_nottoinclude);
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[3],min(float(999),LJets[iantitop_reco].pt),-1,weight,unfolding_set_nottoinclude);
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[4],min(float(999),LJets[itop_reco].pt),-1,weight,unfolding_set_nottoinclude);
	    fill_unfolding_responsematrix_withsets(hist_response_matrix[5],min(float(999),LJets[iantitop_reco].pt),-1,weight,unfolding_set_nottoinclude);
	    
	    hist_response_matrix_finebins[0]->Fill(min(float(999),LJets[itop_reco].pt),-1,weight);
	    hist_response_matrix_finebins[1]->Fill(min(float(999),LJets[itop_reco].pt),-1,weight);

	  }

	if(isgen_phasespace_genak8 && isgen_phasespace)
	  {
	    hist_particlelevl[0]->Fill(min(double(999),GenLJets[0].p4.Pt()),weight);		
	    hist_particlelevl[1]->Fill(min(double(999),GenLJets[1].p4.Pt()),weight);
	    
	    if(delta2R(GenLJets[0].p4,LJets[0].p4) < 0.8) hist_response_matrix_particlelvl[0]->Fill(min(float(999),LJets[0].pt),min(double(999),GenLJets[0].p4.Pt()),weight);
	    if(delta2R(GenLJets[1].p4,LJets[1].p4) < 0.8) hist_response_matrix_particlelvl[1]->Fill(min(float(999),LJets[1].pt),min(double(999),GenLJets[1].p4.Pt()),weight);
	    
	    if(delta2R(GenLJets[0].p4,LJets[0].p4) >= 0.8) hist_response_matrix_particlelvl[0]->Fill(-1,min(double(999),GenLJets[0].p4.Pt()),weight);
	    if(delta2R(GenLJets[1].p4,LJets[1].p4) >= 0.8) hist_response_matrix_particlelvl[1]->Fill(-1,min(double(999),GenLJets[1].p4.Pt()),weight);
	    
	    if(delta2R(GenLJets[0].p4,LJets[0].p4) >= 0.8) hist_response_matrix_particlelvl[0]->Fill(min(float(999),LJets[0].pt),-1,weight);
	    if(delta2R(GenLJets[1].p4,LJets[1].p4) >= 0.8) hist_response_matrix_particlelvl[1]->Fill(min(float(999),LJets[1].pt),-1,weight);	    
	  }
	else
	  {
	    hist_response_matrix_particlelvl[0]->Fill(min(float(999),LJets[0].pt),-1,weight);
	    hist_response_matrix_particlelvl[1]->Fill(min(float(999),LJets[1].pt),-1,weight);
	  }
      
	/*if((int)GenLJets.size() > 0 && delta2R(GenLJets[0].p4,LJets[0].p4) < 0.8 && isgen_phasespace)
	  {
	    hist_response_matrix[4]->Fill(min(float(999),LJets[0].pt),min(float(999),GenLJets[0].pt),weight);
	    hist_response_matrix_finebins[0]->Fill(LJets[0].pt,GenLJets[0].pt,weight);

	    hist_response_matrix_partonlvl[0]->Fill(min(float(999),LJets[0].pt),min(double(999),gentops[GenLJets[0].index_gentopmatched].p4.Pt()),weight);

	  }
	else //if( !((int)GenLJets.size() > 0 && delta2R(GenLJets[0].p4,LJets[0].p4) < 0.8) )
	  {
	    hist_response_matrix[4]->Fill(min(float(999),LJets[0].pt),-1,weight);
	    hist_response_matrix_finebins[0]->Fill(LJets[0].pt,-1,weight);

	    hist_response_matrix_partonlvl[0]->Fill(min(float(999),LJets[0].pt),-1,weight);
	    
	    if(isgen_phasespace)
	      {
		hist_response_matrix[4]->Fill(-1,min(float(999),GenLJets[0].pt),weight);
		hist_response_matrix_finebins[0]->Fill(-1,GenLJets[0].pt,weight);

		hist_response_matrix_partonlvl[0]->Fill(-1,min(double(999),gentops[GenLJets[0].index_gentopmatched].p4.Pt()),weight);
	      }
	      }*/
	  /*if(genak8_index[0] >= 0 && LJets[0].pt > 200 && LJets[0].pt < 700 && GenLJets[genak8_index[0]].pt > 200 && GenLJets[genak8_index[0]].pt < 700)
	  {
	    hist_recolevl[4]->Fill(LJets[0].pt,weight);
	    hist_genlevl[4]->Fill(GenLJets[genak8_index[0]].pt,weight);
	    hist_response_matrix[4]->Fill(LJets[0].pt,GenLJets[genak8_index[0]].pt,weight);
	    }*/
	
      } 
      if((int)Jets.size()>1) {
	hist_init[22]->Fill(LJets[1].pt, weight);
	hist_init[23]->Fill(LJets[1].eta, weight);
	hist_init[24]->Fill(LJets[1].phi, weight);
	
	/*hist_recolevl[5]->Fill(min(float(999),LJets[1].pt),weight);
	if((int)GenLJets.size() > 1 && delta2R(GenLJets[1].p4,LJets[1].p4) < 0.8 && !isgen_phasespace)
	  {
	    hist_response_matrix[5]->Fill(min(float(999),LJets[1].pt),min(float(999),GenLJets[1].pt),weight);
	    hist_response_matrix_finebins[1]->Fill(LJets[1].pt,GenLJets[1].pt,weight);
	  }
	else if( !((int)GenLJets.size() > 0 && delta2R(GenLJets[1].p4,LJets[1].p4) < 0.8) )
	  {
	    hist_response_matrix[5]->Fill(min(float(999),LJets[1].pt),-1,weight);
	    hist_response_matrix_finebins[1]->Fill(LJets[1].pt,-1,weight);
	    }*/
	/*if(genak8_index[1] >= 0 &&  LJets[1].pt > 200 && LJets[1].pt < 700 && GenLJets[genak8_index[1]].pt > 200 && GenLJets[genak8_index[1]].pt < 700)
	  {
	    hist_recolevl[5]->Fill(LJets[1].pt,weight);
	    hist_genlevl[5]->Fill(GenLJets[genak8_index[1]].pt,weight);
	    hist_response_matrix[5]->Fill(LJets[1].pt,GenLJets[genak8_index[1]].pt,weight);
	    }*/

      }

      
    }
  
    
  /*  //Computation of lepton related Suman's variables//
  
  M_l1l2= rat_l2pt_l1pt= deltaPhi_l1l2= l1pt_nearjet= l2pt_nearjet= met_pt= met_phi= met_pt_puppi= met_phi_puppi= delta_phil1_met= delta_phil2_met= delta_phibl1_met= delta_phibl2_met= rat_metpt_ak4pt= rat_metpt_ak8pt= rat_metpt_eventHT= mt_of_l1met= mt_of_l2met= no_ak4jets= no_ak4bjets= no_ak8jets= EventHT= extra_ak4j= ptsum_extra_ak4= extra_ak4jqgl= extra_ak4jdeepb= rat_extra_ak4jpt_lpt= ak81pt= ak81y= ak81mass= ak81sdmass= ak81deep_tvsqcd= ak81deep_wvsqcd= ak82pt= ak82y= ak82mass= ak82sdmass= ak82deep_tvsqcd= ak82deep_wvsqcd= M_bl1= M_bl2 = delta_phibl1bl2= delta_phijl1jl2= deltaR_l1l2= deltaR_l1j1= deltaR_l2j1= deltaR_l1j2= deltaR_l2j2 = deltaR_j1j2 = j1_btag_sc= j2_btag_sc = ak81NHad = ak81chrad = ak81neuhad = ak81tau21 = ak81subhaddiff = ak81tau32 = ak82NHad = ak82chrad = ak82neuhad = ak82tau21 = ak82subhaddiff = ak82tau32 = ak41pt = ak42pt = ak41mass = ak42mass = lep2pt = lep1pt = extra_ak4mass = extra_ak4pt = response_eventclass = -99;
  
  genmatch_sc1 = genmatch_sc2 = 0;
  ak81hase = ak81hasmu = ak82hase = ak82hasmu = -99;
  //hist_new_var[0]->Fill(M_l1l2,weight);
  M_l1l2 = (vleptons[0].p4 + vleptons[1].p4).M();
  
  //Computation of selected lepton related variables (Suman's proposal)
  rat_l2pt_l1pt = min(vleptons[1].p4.Pt(),vleptons[0].p4.Pt())/max(vleptons[1].p4.Pt(),vleptons[0].p4.Pt());
  
  deltaPhi_l1l2 = PhiInRange(vleptons[0].p4.Phi() - vleptons[1].p4.Phi());
 	
  //2d iso variables for l1 and l2//                                             
  float dRl1_min(0.8), dRl2_min(0.8);
  int nearjet_l1(-1), nearjet_l2(-1);
	
  if(Jets.size()>0){
    l1pt_nearjet = ((vleptons[0].p4.Vect()).Perp(Jets[0].p4.Vect()));
    if (Jets.size()>1) { 
      //hist_new_var[3]->Fill(l1pt_nearjet,weight);
      l2pt_nearjet = ((vleptons[1].p4.Vect()).Perp(Jets[1].p4.Vect()));
    }
  }
  //hist_new_var[4]->Fill(l2pt_nearjet,weight);
  
  //Computation of MET related Suman's proposed variables//
  if (miset_PUPPI > -999) { //  && miset_PUPPIPhi != -1000) { 
    TLorentzVector metvector (miset_PUPPI*std::cos(misphi_PUPPI),miset_PUPPI*std::sin(misphi_PUPPI),0,miset_PUPPI);  // To calculate transvere mass

    met_pt = chs_met;
    met_pt_puppi = miset_PUPPI;
    
    met_phi = chs_misphi;
    met_phi_puppi = misphi_PUPPI;
		
    delta_phil1_met = PhiInRange(vleptons[0].p4.Phi() - metvector.Phi());
  		
    delta_phil2_met = PhiInRange(vleptons[1].p4.Phi() - metvector.Phi());
		
    TLorentzVector bl1_syst;
    TLorentzVector nearbl1;
    if (Jets.size()>0) { 
      nearbl1 = Jets[0].p4;
      bl1_syst = nearbl1 + vleptons[0].p4;
      M_bl1 = bl1_syst.M();
      delta_phibl1_met = PhiInRange(bl1_syst.Phi() - metvector.Phi());
      if (Jets.size()>1) { 
	TLorentzVector nearbl2;
	TLorentzVector bl2_syst;
	nearbl2 = Jets[1].p4;
	bl2_syst = nearbl2 + vleptons[1].p4;
	M_bl2 = bl2_syst.M();
	delta_phibl2_met = PhiInRange(bl2_syst.Phi() - metvector.Phi());
	delta_phibl1bl2 = PhiInRange(bl1_syst.Phi() - bl2_syst.Phi());
      }
    }
    if (Jets.size()>0) { 
      rat_metpt_ak4pt = metvector.Pt()/max(float(1.0),Jets[0].pt);
      rat_metpt_ak8pt = metvector.Pt()/max(float(1.0),LJets[0].pt);
    }
    //GMA defined later on    rat_metpt_eventHT = metvector.Pt()/Event_HT;
    
    mt_of_l1met = (metvector+vleptons[0].p4).Mt();
    mt_of_l2met = (metvector+vleptons[1].p4).Mt();
  }
  
  no_ak4jets = Jets.size();
  no_ak4bjets = BJets.size();
  no_ak8jets = LJets.size();
 
  int nAK4inAK8=0;
  bool found=false;  /// whether we have found the leading extra ak4 jet or not
  int extra_leadak4_index=-1;
  ptsum_extra_ak4 = 0;
  
  for(int ijet=0;ijet<(int)Jets.size();ijet++)
    {
      if( delta2R(LJets[0].p4,Jets[ijet].p4)<0.8  || delta2R(LJets[min(1,(int)LJets.size()-1)].p4,Jets[ijet].p4)<0.8) 
	//if( delta2R(LJets[fjet].y,LJets[fjet].phi,Jets[ijet].y,Jets[ijet].phi)<0.8 )
	{
	  nAK4inAK8++;
	}
      else{
	ptsum_extra_ak4 = ptsum_extra_ak4 + Jets[ijet].pt;
	if(!found)
	  {
	    extra_leadak4_index=ijet;
	    found=true;
	  }
      }

    }
  extra_ak4j = Jets.size() - nAK4inAK8;
  
  if(extra_ak4j>=1 && extra_leadak4_index >=0){
    extra_ak4jqgl = Jets[extra_leadak4_index].qgl;
    extra_ak4jdeepb = Jets[extra_leadak4_index].btag_DeepFlav;
    extra_ak4pt = Jets[extra_leadak4_index].pt;
    extra_ak4mass = Jets[extra_leadak4_index].mass;
    if (delta2R(Jets[extra_leadak4_index].p4,vleptons[0].p4) < delta2R(Jets[extra_leadak4_index].p4,vleptons[1].p4)) {
      rat_extra_ak4jpt_lpt = Jets[extra_leadak4_index].pt/max(1.0,vleptons[0].p4.Pt());
    } else {
      rat_extra_ak4jpt_lpt = Jets[extra_leadak4_index].pt/max(1.0,vleptons[1].p4.Pt());
    }
  }
  //  EventHT = Event_HT;

  //GMA add lepton also there;
  //calculate the directly global transverse thrust i.e. dirgltrthr
  dirgltrthr = 0;
  dirglthrmin = 0;

  if (Jets.size()>=1) {
    std::vector<TLorentzVector> allsjets_4v;
    double Pt_sum(0.);
    for(unsigned ijet=0; ijet< Jets.size(); ijet++){
      TLorentzVector sjv = Jets[ijet].p4;
      allsjets_4v.push_back(sjv);
      Pt_sum += Jets[ijet].pt;
    }
		
    for(unsigned ilep=0; ilep<velectrons.size(); ilep++){
      TLorentzVector sjv = velectrons[ilep].p4;
      allsjets_4v.push_back(sjv);
      Pt_sum += velectrons[ilep].p4.Pt();
    }

    for(unsigned ilep=0; ilep<vmuons.size(); ilep++){
      TLorentzVector sjv = vmuons[ilep].p4;
      allsjets_4v.push_back(sjv);
      Pt_sum += vmuons[ilep].p4.Pt();
    }
    
    EventHT = Pt_sum;	
    rat_metpt_eventHT = miset_PUPPI/max(float(1.0), EventHT);
		
    std::vector<double> ThrustAxis;
    std::vector<double> Thrust;
    
    for(unsigned int j =0;j<4;j++){                                                            
      Thrust.push_back(0.);                                                                    
    }
    Thrust_calculate(allsjets_4v,Thrust);
    dirgltrthr = Thrust[3];
    
    //now comes directly global thrust minor dirglthrmin                         
    //rotate the coordinate system around the beam axis such that 
    //the thrust axis is the new y'-Axis - the projections are                   
    //simply the new y-values then                                               

    double alpha=atan2(Thrust[1],Thrust[0]);
    for(unsigned int i=0; i<allsjets_4v.size(); i++){
      dirglthrmin += fabs(-sin(alpha)*allsjets_4v[i].Px()+cos(alpha)*allsjets_4v[i].Py());
    }
    dirglthrmin = dirglthrmin/max(1.0, Pt_sum);
  }

  // ********thrust calculation ended******
  ak81pt = LJets[0].pt;
  ak81y = LJets[0].y;
  ak81mass = LJets[0].mass;
  ak81sdmass = LJets[0].sdmass;
  ak81deep_tvsqcd = LJets[0].DeepTag_TvsQCD;
  ak81deep_wvsqcd = LJets[0].DeepTag_WvsQCD;
  ak81hasmu = (int)LJets[0].hasmatchmu;
  ak81hase = (int)LJets[0].hasmatche;
  if (LJets.size()>1 ) {
    ak82pt = LJets[1].pt;
    ak82y = LJets[1].y;
    ak82mass = LJets[1].mass;
    ak82sdmass = LJets[1].sdmass;
    ak82deep_tvsqcd = LJets[1].DeepTag_TvsQCD;
    ak82deep_wvsqcd = LJets[1].DeepTag_WvsQCD;
    ak82hasmu = (int)LJets[1].hasmatchmu;
    ak82hase = (int)LJets[1].hasmatche;
  }

  deltaR_l1l2 = delta2R(vleptons[0].p4, vleptons[1].p4);
  
  //deltaR_l1b1 = delta2R(l1.Rapidity(),l1.Phi(),bjv[0].Rapidity(),bjv[0].Phi());
  //if (bjv.size() >1) deltaR_l1b2 = delta2R(vleptons[0].p4.Rapidity(),vleptons[0].p4.Phi(),bjv[1].Rapidity(),bjv[1].Phi());
  //deltaR_l2b1 = delta2R(l2.Rapidity(),l2.Phi(),bjv[0].Rapidity(),bjv[0].Phi());
  //if (bjv.size() >1) deltaR_l2b2 = delta2R(l2.Rapidity(),l2.Phi(),bjv[1].Rapidity(),bjv[1].Phi());
  //L1 has to attached with large radius jets
	
  deltaR_l1j1 = delta2R(vleptons[0].p4,Jets[0].p4);
  deltaR_l1j2 = delta2R(vleptons[0].p4,Jets[1].p4);
  deltaR_l2j1 = delta2R(vleptons[1].p4,Jets[0].p4);
  deltaR_l2j2 = delta2R(vleptons[1].p4,Jets[1].p4);
  deltaR_j1j2 = delta2R(Jets[0].p4,Jets[1].p4);

  j1_btag_sc = Jets[0].btag_DeepFlav;
  j2_btag_sc = Jets[1].btag_DeepFlav;

  btag_SF = btag_SF_up = btag_SF_dn = 1;
  int ijet_maxbtag;
  if (Jets[1].btag_DeepFlav > Jets[0].btag_DeepFlav)
    ijet_maxbtag = 1;
  else
    ijet_maxbtag = 0;
  
  btag_SF    *= Jets[ijet_maxbtag].btag_DeepFlav_med_SF;
  btag_SF_up *= Jets[ijet_maxbtag].btag_DeepFlav_med_SF_up;
  btag_SF_dn *= Jets[ijet_maxbtag].btag_DeepFlav_med_SF_dn;

  btag_SF    *= Jets[int(max(0,ijet_maxbtag-1))].btag_DeepFlav_loose_SF;
  btag_SF_up *= Jets[int(max(0,ijet_maxbtag-1))].btag_DeepFlav_loose_SF_up;
  btag_SF_dn *= Jets[int(max(0,ijet_maxbtag-1))].btag_DeepFlav_loose_SF_dn;

  
  if(isMC){
    /// for selected AK8jets
    double dr;
    int it=-1;
    float genmatch_sc13=0,genmatch_sc23=0;

    if(topblep[0][0]>-1 && topblep[0][1] >-1 && topblep[0][2]>-1 && topblep[1][0]>-1 && topblep[1][1] >-1 && topblep[1][2]>-1 && (int)vleptons.size()>1 && (int)Jets.size()>1 && (int)LJets.size()>1)
      {
	it= (genpartpdg[topblep[0][2]]*vleptons[0].pdgId >0) ? 0 : 1;
	TLorentzVector lep1,bjet1;
	bjet1.SetPtEtaPhiM(genpartpt[topblep[it][1]],genparteta[topblep[it][1]],genpartphi[topblep[it][1]],genpartm[topblep[it][1]]);
	lep1.SetPtEtaPhiM(genpartpt[topblep[it][2]],genparteta[topblep[it][2]],genpartphi[topblep[it][2]],genpartm[topblep[it][2]]);
	
	dr = delta2R(vleptons[0].p4,lep1);
	hist_genmatch_deltaR[0]->Fill(dr,weight);
	if( dr < 0.02)
	  genmatch_sc1 += 1;
	
	
	dr = delta2R(Jets[0].p4,bjet1);   
	hist_genmatch_deltaR[1]->Fill(dr,weight);
	if(dr  < 0.1)
	  genmatch_sc1 += 10;
	
	dr = delta2R(LJets[0].p4,(bjet1+lep1));
	hist_genmatch_deltaR[2]->Fill(dr,weight);
	if(dr <  0.2)
	  genmatch_sc1 += 100;

	hist_genmatch_deltaR[3]->Fill(delta2R(LJets[0].y,LJets[0].phi,genparteta[topblep[it][0]],genpartphi[topblep[it][0]]),weight);

	
	if(LJets.size()>1 && (int)vleptons.size()>1 && (int)Jets.size()>1){
	  
	  it = (genpartpdg[topblep[0][2]]*vleptons[1].pdgId >0) ? 0 : 1;
	  TLorentzVector lep2,bjet2;
	  bjet2.SetPtEtaPhiM(genpartpt[topblep[it][1]],genparteta[topblep[it][1]],genpartphi[topblep[it][1]],genpartm[topblep[it][1]]);
	  lep2.SetPtEtaPhiM(genpartpt[topblep[it][2]],genparteta[topblep[it][2]],genpartphi[topblep[it][2]],genpartm[topblep[it][2]]);
	  
	  dr = delta2R(vleptons[1].p4,lep2);
	  hist_genmatch_deltaR[4]->Fill(dr,weight);
	  if(dr < 0.02)
	    genmatch_sc2 += 1;
	      
	  dr = delta2R(Jets[1].p4,bjet2);
	  hist_genmatch_deltaR[5]->Fill(dr,weight);
	  if(dr < 0.1)
	    genmatch_sc2 += 10;
	
	  dr = delta2R(LJets[1].p4,bjet2+lep2);
	  hist_genmatch_deltaR[6]->Fill(dr,weight);
	  if(dr < 0.2)
	    genmatch_sc2 += 100;

	  hist_genmatch_deltaR[7]->Fill(delta2R(LJets[1].y,LJets[1].phi,genparteta[topblep[it][0]],genpartphi[topblep[it][0]]),weight);
	  
	}
      
	hist_genmatch_deltaR_score[0]->Fill(genmatch_sc1,weight);
	hist_genmatch_deltaR_score[1]->Fill(genmatch_sc2,weight);
      }
  }

  lep1pt=vleptons[0].p4.Pt();
  lep2pt=vleptons[1].p4.Pt();
     
  ak81NHad = LJets[0].NHadF;
  ak81chrad =LJets[0].chrad;
  ak81neuhad =LJets[0].neunhadfrac;
  ak81tau21 =LJets[0].tau21;
  ak81subhaddiff =LJets[0].subhaddiff;
  ak81tau32 = LJets[0].tau32;
  if(LJets.size() >1 ){
    ak82NHad =LJets[1].NHadF;
    ak82chrad =LJets[1].chrad;
    ak82neuhad =LJets[1].neunhadfrac;
    ak82tau21 =LJets[1].tau21;
    ak82subhaddiff =LJets[1].subhaddiff;
    ak82tau32 =LJets[1].tau32;
  }

  //response_eventclass = reader2->EvaluateMVA("BDTG method");
  
  
  if(isTT && isMC &&  gentops.size() > 1)
    {
      AssignGen(LJets,genpartons);

      int ilep = (vleptons[0].pdgId > 0) ? 0 : 1;   ///  choosing top first ie anti lepton in final state
      int itop = (gentops[0].daughter[0].pdgId > 0) ? 0 : 1;
      pt_gen_l1 = gentops[itop].daughter[0].p4.Pt();
      pt_gen_b1 = gentops[itop].daughter[2].p4.Pt();
      pt_gen_lb1 = (gentops[itop].daughter[0].p4 + gentops[itop].daughter[2].p4).Pt();
      deltaR_gen_lb1 = delta2R(gentops[itop].daughter[0].p4 , gentops[itop].daughter[2].p4);
      delptbypt_gen_l1 = (gentops[itop].daughter[0].p4 - vleptons[ilep].p4).Pt()/vleptons[ilep].pt;
      if(genak4_index[itop] > -0.5) delptbypt_gen_b1 = (GenJets[genak4_index[itop]].p4 - Jets[ilep].p4).Pt() / Jets[ilep].pt;
      else delptbypt_gen_b1 = 100;
      delptbypt_gen_lb1 = (gentops[itop].daughter[0].p4 + gentops[itop].daughter[2].p4 - LJets[ilep].p4).Pt()/ LJets[ilep].pt;

      ak81hasgene = (int)LJets[ilep].haselectron; 
      ak81hasgenmu = (int)LJets[ilep].hasmuon;
      ak81hasgenb = (int)LJets[ilep].hasb;
      ak81hasgenhasalldecay = (int)LJets[ilep].hasleptop_alldecay;

      
      ilep = 1 - ilep;  ///  choosing anti-top next ie lepton in final state
      itop = 1 - itop;
      pt_gen_l2 = gentops[itop].daughter[0].p4.Pt();
      pt_gen_b2 = gentops[itop].daughter[2].p4.Pt();
      pt_gen_lb2 = (gentops[itop].daughter[0].p4 + gentops[itop].daughter[2].p4).Pt();
      deltaR_gen_lb2 = delta2R(gentops[itop].daughter[0].p4 , gentops[itop].daughter[2].p4);
      delptbypt_gen_l2 = (gentops[itop].daughter[0].p4 - vleptons[ilep].p4).Pt()/vleptons[ilep].pt;
      if(genak4_index[itop] > -0.5) delptbypt_gen_b2 = (GenJets[genak4_index[itop]].p4 - Jets[ilep].p4).Pt() / Jets[ilep].pt;
      else delptbypt_gen_b2 = 100;
      delptbypt_gen_lb2 = (gentops[itop].daughter[0].p4 + gentops[itop].daughter[2].p4 - LJets[ilep].p4).Pt()/ LJets[ilep].pt;

      ak82hasgene = (int)LJets[ilep].haselectron; 
      ak82hasgenmu = (int)LJets[ilep].hasmuon;
      ak82hasgenb = (int)LJets[ilep].hasb;
      ak82hasgenhasalldecay = (int)LJets[ilep].hasleptop_alldecay;
    }
  
  
  hist_treevar[0]->Fill(M_l1l2,weight);
  hist_treevar[1]->Fill(rat_l2pt_l1pt,weight);
  hist_treevar[2]->Fill(deltaPhi_l1l2,weight);
  hist_treevar[3]->Fill(l1pt_nearjet,weight);
  hist_treevar[4]->Fill(l2pt_nearjet,weight);
  hist_treevar[5]->Fill(met_pt,weight);
  hist_treevar[6]->Fill(met_phi,weight);
  hist_treevar[7]->Fill(M_bl1,weight);
  hist_treevar[8]->Fill(M_bl2,weight);
  hist_treevar[9]->Fill(delta_phil1_met,weight);
  hist_treevar[10]->Fill(delta_phil2_met,weight);
  hist_treevar[11]->Fill(delta_phibl1_met,weight);
  hist_treevar[12]->Fill(delta_phibl2_met,weight);
  hist_treevar[13]->Fill(rat_metpt_ak4pt,weight);
  hist_treevar[14]->Fill(rat_metpt_ak8pt,weight);
  hist_treevar[15]->Fill(rat_metpt_eventHT,weight);
  hist_treevar[16]->Fill(mt_of_l1met,weight);
  hist_treevar[17]->Fill(mt_of_l2met,weight);
  hist_treevar[18]->Fill(no_ak4jets,weight);
  hist_treevar[19]->Fill(no_ak4bjets,weight);
  hist_treevar[20]->Fill(no_ak8jets,weight);
  hist_treevar[21]->Fill(EventHT,weight);
  hist_treevar[22]->Fill(extra_ak4j,weight);
  hist_treevar[23]->Fill(ptsum_extra_ak4,weight);
  hist_treevar[24]->Fill(extra_ak4pt,weight);
  hist_treevar[25]->Fill(extra_ak4mass,weight);
  hist_treevar[26]->Fill(extra_ak4jqgl,weight);
  hist_treevar[27]->Fill(extra_ak4jdeepb,weight);
  hist_treevar[28]->Fill(rat_extra_ak4jpt_lpt,weight);
  hist_treevar[29]->Fill(ak81deep_tvsqcd,weight);
  hist_treevar[30]->Fill(ak81deep_wvsqcd,weight);
  hist_treevar[31]->Fill(ak82deep_tvsqcd,weight);
  hist_treevar[32]->Fill(ak82deep_wvsqcd,weight);
  hist_treevar[33]->Fill(delta_phibl1bl2,weight);
  hist_treevar[34]->Fill(deltaR_l1l2,weight);
  hist_treevar[35]->Fill(deltaR_l1j1,weight);
  hist_treevar[36]->Fill(deltaR_l2j1,weight);
  hist_treevar[37]->Fill(deltaR_l1j2,weight);
  hist_treevar[38]->Fill(deltaR_l2j2,weight);
  hist_treevar[39]->Fill(j1_btag_sc,weight);
  hist_treevar[40]->Fill(j2_btag_sc,weight);
  hist_treevar[41]->Fill(dirgltrthr,weight);
  hist_treevar[42]->Fill(dirglthrmin,weight);
  hist_treevar[43]->Fill(response_ak81,weight);
  hist_treevar[44]->Fill(response_ak82,weight);
  hist_treevar[45]->Fill(deltaR_l1l2,weight);  
  hist_treevar[46]->Fill(response_eventclass,weight);
    
  if(isTT && isMC) hist_2d_pt_gentopvsgentop->Fill(max((float)LHEtops[0].pt,(float)LHEtops[1].pt),min((float)LHEtops[0].pt,(float)LHEtops[1].pt),weight);

  if( Jets.size()>0 && LJets.size()>0)  hist_2d_pt_vsbtagsc[2]->Fill(LJets[0].pt,Jets[0].btag_DeepFlav,weight);
  if( Jets.size()>1 && LJets.size()>1)  hist_2d_pt_vsbtagsc[3]->Fill(LJets[1].pt,Jets[1].btag_DeepFlav,weight);
  */

  Tnewvar->Fill();

  // end //                                                                            
  //if(gProofServ) {  str = TString::Format("check end evt %d ",ievt);gProofServ->SendAsynMessage(str);	}
  return kTRUE;     
}
void Anal_Leptop_PROOF::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  fileOut->cd();
  fileOut->Write();
  fOutput->Add(OutFile);
  fileOut->Close();
}

void Anal_Leptop_PROOF::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
}
