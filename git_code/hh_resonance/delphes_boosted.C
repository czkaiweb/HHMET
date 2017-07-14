#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

int delphes_boosted(/*char* string*/){

//  TString finname(string);
//  TString foutname(string);
  
//  finname = "/data2/weikai/calchep/delphes_results/HHMET/"+finname+"/delphes/myprocess_atlas13tev.root";
//  foutname ="analysis_results/"+foutname+"_boosted_results.root";
//  TFile *fin=new TFile(finname,"READONLY");
//  TFile *fout=new TFile(foutname,"RECREATE");

//  TFile *fin=new TFile("/data2/weikai/calchep/delphes_results/HHMET/gYqq0.05MY600MF200/delphes/myprocess_cms8tev.root","READONLY");
  TFile *fin=new TFile("/data2/weikai/calchep/delphes_results/HHMET/ttbar/delphes/myprocess_atlas13tev.root","READONLY");
//  TFile *fout = new TFile("test_boosted_signal.root","RECREATE");
  TFile *fout = new TFile("test_bkg.root","RECREATE");

//----------------SETUP HISTOGRAM------------------
  TH1D *info=new TH1D("info","info",14,0,14);
  TH1D *h_leading_pt = new TH1D("h_leading_pt","leading_pt",20,400,1000);
  TH1D *h_subleading_pt = new TH1D("h_subleading_pt","subleading_pt",20,200,800);
  TH1D *h_m_lead = new TH1D("h_m_lead","m_lead",20,0,200);
  TH1D *h_m_subl = new TH1D("h_m_subl","m_subl",20,0,200);
  TH1D *h_m_2h = new TH1D("h_m_2h","m_2h",50,0,1500);

//----------------END OF HISTOGRAM---------------------  
  TTree *t = (TTree*)fin->Get("Delphes");
  Long64_t nevent = t->GetEntries();
  cout << nevent <<endl;
  t->SetMakeClass(1);
//----------------SETUP BRANCHES-----------------------
  Int_t nMET=-1;
  Float_t MET[10];
  Float_t MET_PHI[10];
  TBranch *met=t->GetBranch("MissingET.MET");
  TBranch *met_phi= t->GetBranch("MissingET.Phi");
  t->SetBranchAddress("MissingET",&nMET);
  met->SetAddress(&MET);
  met_phi->SetAddress(&MET_PHI);

  Int_t nJET=-1;
  const Int_t Jet_size=50;
  Float_t JET_PT[Jet_size];
  Float_t JET_ETA[Jet_size];
  Float_t JET_PHI[Jet_size];
  Float_t JET_M[Jet_size];
  UInt_t JET_BTAG[Jet_size];
  TLorentzVector all_vector[Jet_size];
//  TLorentzVector b_vector[4];

  t->SetBranchAddress("Jet",&nJET);

  TBranch *jet_pt  = t->GetBranch("Jet.PT");
  jet_pt->SetAddress(&JET_PT);
  TBranch *jet_eta = t->GetBranch("Jet.Eta");
  jet_eta->SetAddress(&JET_ETA);
  TBranch *jet_phi = t->GetBranch("Jet.Phi");
  jet_phi->SetAddress(&JET_PHI);
  TBranch *jet_m = t->GetBranch("Jet.Mass");
  jet_m->SetAddress(&JET_M);
  TBranch *jet_btag = t->GetBranch("Jet.BTag");
  jet_btag->SetAddress(&JET_BTAG);
  

  Int_t     kMaxJetAK10Trim = 50;
  Int_t           nJetAK10Trim = -1;
  Float_t         JetAK10Trim_PT[kMaxJetAK10Trim];   //[JetAK10Trim_]
  Float_t         JetAK10Trim_Eta[kMaxJetAK10Trim];   //[JetAK10Trim_]
  Float_t         JetAK10Trim_Phi[kMaxJetAK10Trim];   //[JetAK10Trim_]
  Float_t         JetAK10Trim_Mass[kMaxJetAK10Trim];   //[JetAK10Trim_]
  TBranch        *b_JetAK10Trim_PT;   //!
  TBranch        *b_JetAK10Trim_Eta;   //!
  TBranch        *b_JetAK10Trim_Phi;   //!
  TBranch        *b_JetAK10Trim_Mass;   //!
  t->SetBranchAddress("JetAK10Trim", &nJetAK10Trim);
  t->SetBranchAddress("JetAK10Trim.PT", JetAK10Trim_PT, &b_JetAK10Trim_PT);
  t->SetBranchAddress("JetAK10Trim.Eta", JetAK10Trim_Eta, &b_JetAK10Trim_Eta);
  t->SetBranchAddress("JetAK10Trim.Phi", JetAK10Trim_Phi, &b_JetAK10Trim_Phi);
  t->SetBranchAddress("JetAK10Trim.Mass", JetAK10Trim_Mass, &b_JetAK10Trim_Mass);

  const Int_t     kMaxJetAK2Track = 50;
  Int_t           nJetAK2Track = -1;
  Float_t         JetAK2Track_PT[kMaxJetAK2Track];   //[JetAK2Track_]
  Float_t         JetAK2Track_Eta[kMaxJetAK2Track];   //[JetAK2Track_]
  Float_t         JetAK2Track_Phi[kMaxJetAK2Track];   //[JetAK2Track_]
  Float_t         JetAK2Track_T[kMaxJetAK2Track];   //[JetAK2Track_]
  Float_t         JetAK2Track_Mass[kMaxJetAK2Track];   //[JetAK2Track_]
  UInt_t          JetAK2Track_BTag[kMaxJetAK2Track];   //[JetAK2Track_]
  Int_t           JetAK2Track_Flavor[kMaxJetAK2Track];
  TBranch        *b_JetAK2Track_PT;   //!
  TBranch        *b_JetAK2Track_Eta;   //!
  TBranch        *b_JetAK2Track_Phi;   //!
  TBranch        *b_JetAK2Track_Mass;   //!
  TBranch        *b_JetAK2Track_BTag;   //!
  TBranch        *b_JetAK2Track_Flavor;
  t->SetBranchAddress("JetAK2Track", &nJetAK2Track);
  t->SetBranchAddress("JetAK2Track.PT", JetAK2Track_PT, &b_JetAK2Track_PT);
  t->SetBranchAddress("JetAK2Track.Eta", JetAK2Track_Eta, &b_JetAK2Track_Eta);
  t->SetBranchAddress("JetAK2Track.Phi", JetAK2Track_Phi, &b_JetAK2Track_Phi);
  t->SetBranchAddress("JetAK2Track.Mass", JetAK2Track_Mass, &b_JetAK2Track_Mass);
  t->SetBranchAddress("JetAK2Track.BTag", JetAK2Track_BTag, &b_JetAK2Track_BTag);
  t->SetBranchAddress("JetAK2Track.Flavor", JetAK2Track_Flavor, &b_JetAK2Track_Flavor);

  Long64_t cut_flow[8]={0,0,0,0,0,0,0,0};
  cout << "start loop"<<endl;

  for(Long64_t i=0;i<10000/*nevent*/;i++){
  //  if(i%100==0) cout << i <<endl;
    cut_flow[0]++;
//    t->GetEntry(i);

    jet_pt->GetEntry(i);
    jet_phi->GetEntry(i);
    jet_eta->GetEntry(i);
    jet_m->GetEntry(i);
    jet_btag->GetEntry(i);
    b_JetAK10Trim_PT->GetEntry(i);
    b_JetAK10Trim_Eta->GetEntry(i);
    b_JetAK10Trim_Phi->GetEntry(i);
    b_JetAK10Trim_Mass->GetEntry(i);
    b_JetAK2Track_PT->GetEntry(i);
    b_JetAK2Track_Eta->GetEntry(i);
    b_JetAK2Track_Phi->GetEntry(i);
    b_JetAK2Track_Mass->GetEntry(i);
    b_JetAK2Track_BTag->GetEntry(i);


//////////////////////////////////////////////////////////////////////////////////
////                          Boosted analysis                                ////
////                    based on ATLAS-CONF-2016-049                          ////
//////////////////////////////////////////////////////////////////////////////////

//-------------------------large-R jet(R=1.0)>=2---------------------------
    Int_t Nb_T=0;
    vector< int > fat_index;
    if(nJetAK10Trim<2) continue;
    for(int j=0;j<nJetAK10Trim ;j++){
      if(JetAK10Trim_PT[j]>=250 && fabs(JetAK10Trim_Eta[j])<=2.0 && JetAK10Trim_Mass[j]>=50){
        Double_t pt_leading = 0;
	if(JetAK10Trim_PT[j]>pt_leading) pt_leading=JetAK10Trim_PT[j];
        if(pt_leading>=450){
	  fat_index.push_back(j);	
        }
      }
    }

//    if(JetAK10Trim_PT[0]<450 || fabs(JetAK10Trim_Eta[0])>2.0 || JetAK10Trim_Mass[0]<50) continue;
//    if(JetAK10Trim_PT[1]<250 || fabs(JetAK10Trim_Eta[1])>2.0 || JetAK10Trim_Mass[1]<50) continue;
    if(fat_index.size()<2) continue;
    cut_flow[1]++;     // cut: Nb_T>=2

    TLorentzVector h_vector[2];
    h_vector[0].SetPtEtaPhiM(JetAK10Trim_PT[fat_index[0]],JetAK10Trim_Eta[fat_index[0]],JetAK10Trim_Phi[fat_index[0]],JetAK10Trim_Mass[fat_index[0]]);
    h_vector[1].SetPtEtaPhiM(JetAK10Trim_PT[fat_index[1]],JetAK10Trim_Eta[fat_index[1]],JetAK10Trim_Phi[fat_index[1]],JetAK10Trim_Mass[fat_index[1]]);
    TLorentzVector track_vector[kMaxJetAK2Track];
    for(Int_t ntrack=0;ntrack<nJetAK2Track;ntrack++){
        track_vector[ntrack].SetPtEtaPhiM(JetAK2Track_PT[ntrack],JetAK2Track_Eta[ntrack],JetAK2Track_Phi[ntrack],JetAK2Track_Mass[ntrack]);
    }
    

//--------------------------------delta Eta_hh---------------------------

    if( fabs(h_vector[0].Eta()-h_vector[1].Eta())>1.7 ) continue;
    cut_flow[2]++;
     
////--------------------------b tagged track-jet >= 2-------------------------
    Int_t b_flag[2]={0,0};
    Int_t ntrack_bjet=0;
    for(Int_t nfatjet=0;nfatjet<2;nfatjet++){
        for(Int_t ntrackjet=0; ntrackjet<nJetAK2Track; ntrackjet++){
	    if(track_vector[ntrackjet].Pt()>10.0 && fabs(track_vector[ntrackjet].Eta())<2.5){
                if(track_vector[ntrackjet].DeltaR(h_vector[nfatjet])<1.0){
		    if(JetAK2Track_BTag[ntrackjet]>0){
                        b_flag[nfatjet]++;
			ntrack_bjet++;
		    }
                }
	    }
	}
    }
    if(b_flag[0]==0 || b_flag[1]==0 || ntrack_bjet<2) continue;
    cut_flow[3]++;
    h_leading_pt->Fill(h_vector[0].Pt());
    h_subleading_pt->Fill(h_vector[1].Pt());

////////////////////////////////////////////////////////////////////////
////         SR based on leading and sub-leading mass               ////
///////////////////////////////////////////////////////////////////////
    Double_t Xhh=-10;
    Double_t m_lead,m_subl;
    m_lead = h_vector[0].M();
    m_subl = h_vector[1].M();
    TLorentzVector invariant_mass;
    invariant_mass = h_vector[0]+h_vector[1];
    Xhh = TMath::Sqrt( ((m_lead-124)/(0.1*m_lead))*((m_lead-124)/(0.1*m_lead)) + ((m_subl-115)/(0.1*m_subl))*((m_subl-115)/(0.1*m_subl)) ) ;
    if(Xhh>1.6) continue;
    cut_flow[4]++;
    
    h_m_lead->Fill(m_lead);
    h_m_subl->Fill(m_subl);
    h_m_2h->Fill(invariant_mass.M());
    if(ntrack_bjet==2) cut_flow[5]++;
    if(ntrack_bjet==3) cut_flow[6]++;
    if(ntrack_bjet>=4) cut_flow[7]++;
} 
    

info->Fill("No selection",cut_flow[0]);
info->Fill("large-R jets",cut_flow[1]);
info->Fill("|#Delta#eta_{hh}|",cut_flow[2]);
info->Fill("b-tagged track-jet>=2",cut_flow[3]);
info->Fill("Xhh<1.6",cut_flow[4]);
info->Fill("b-tagged track-jet=2",cut_flow[5]);
info->Fill("b-tagged track-jet=3",cut_flow[6]);
info->Fill("b-tagged track-jet=4",cut_flow[7]);

fout->cd();
info->Write();
h_leading_pt->Write();
h_subleading_pt->Write();
h_m_lead->Write();
h_m_subl->Write();
h_m_2h->Write();

  ofstream fileout("ttbar_boosted_cutflow_check.dat");
  fileout << "No Selection :                              "    << " " << cut_flow[0] << endl;
  fileout << "large-R jets :                              "    << " " << cut_flow[1] << endl;
  fileout << "Delta eta_{hh} :                            "    << " " << cut_flow[2] << endl;
  fileout << "b-tagged track-jet>=2 :                     "    << " " << cut_flow[3] << endl;
  fileout << "Xhh<1.6 :                                   "    << " " << cut_flow[4] << endl;
  fileout << "b-tagged track-jet=2 :                      "    << " " << cut_flow[5] << endl;
  fileout << "b-tagged track-jet=3 :                      "    << " " << cut_flow[6] << endl;
  fileout << "b-tagged track-jet=4 :                      "    << " " << cut_flow[7] << endl;
  fileout.close();


return 1;
}
