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

#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;

int Hbb_analysis2(){

  TFile *fin=new TFile("/data2/weikai/calchep/delphes_results/HHMET/gYqq0.05MY600MF200/delphes/myprocess_cms8tev.root","READONLY");
  TFile *fout=new TFile("signal.root","RECREATE");

  //----------------SETUP HISTOGRAM------------------    
  TH1D *info = new TH1D("info","info",14,0,14);
  TH1D *h_leading_pt = new TH1D("h_leading_pt","leading_pt",50,0,500);
  TH1D *h_subleading_pt = new TH1D("h_subleading_pt","subleading_pt",50,0,500);
  TH1D *h_m_average = new TH1D("h_m_average","m_average",200,0,200);
  TH1D *h_MET = new TH1D("h_MET","h_MET",50,0,500);
  TH1D *h_dRmax = new TH1D("h_dRmax","dRmax",24,0,4);
  //----------------END OF HISTOGRAM--------------------- 

  TTree *t = (TTree*)fin->Get("Delphes");
  Long64_t nevent = t->GetEntries();
  cout << nevent << endl;
  t->SetMakeClass(1);
  
  //----------------SETUP BRANCHES----------------------- 
  Int_t nMET = -1;
  Float_t MET[10];
  Float_t MET_PHI[10];
  t->SetBranchAddress("MissingET", &nMET);
  TBranch *met = t->GetBranch("MissingET.MET");
  met->SetAddress(&MET);
  TBranch *met_phi = t->GetBranch("MissingET.Phi");
  met_phi->SetAddress(&MET_PHI);

  Int_t nJET = -1;
  Float_t JET_PT[20];
  Float_t JET_ETA[20];
  Float_t JET_PHI[20];
  Float_t JET_M[20];
  UInt_t JET_BTAG[20];
  TLorentzVector all_vector[20];
  t->SetBranchAddress("Jet", &nJET);
  TBranch *jet_pt = t->GetBranch("Jet.PT");
  jet_pt->SetAddress(&JET_PT);
  TBranch *jet_eta = t->GetBranch("Jet.Eta");
  jet_eta->SetAddress(&JET_ETA);
  TBranch *jet_phi = t->GetBranch("Jet.Phi");
  jet_phi->SetAddress(&JET_PHI);
  TBranch *jet_m = t->GetBranch("Jet.Mass");
  jet_m->SetAddress(&JET_M);
  TBranch *jet_btag = t->GetBranch("Jet.BTag");
  jet_btag->SetAddress(&JET_BTAG);

  Int_t nElec = -1;
  Float_t ELEC_PT[10];
  Float_t ELEC_ETA[10];
  Float_t ELEC_PHI[10];
  Float_t ELEC_SUMPT[10];
  Float_t ELEC_HOE[10];
  t->SetBranchAddress("Electron", &nElec);
  TBranch *e_pt = t->GetBranch("Electron.PT");
  e_pt->SetAddress(&ELEC_PT);
  TBranch *e_eta = t->GetBranch("Electron.Eta");
  e_eta->SetAddress(&ELEC_ETA);
  TBranch *e_phi = t->GetBranch("Electron.Phi");
  e_phi->SetAddress(&ELEC_PHI);
  TBranch *e_sumpt = t->GetBranch("Electron.SumPtCharged");
  e_sumpt->SetAddress(&ELEC_SUMPT);
  TBranch *e_hoe = t->GetBranch("Electron.EhadOverEem");
  e_hoe->SetAddress(&ELEC_HOE);

  Int_t nMUON = -1;
  Float_t MUON_PT[10];
  Float_t MUON_ETA[10];
  Float_t MUON_PHI[10];
  Float_t MUON_SUMPT[10];
  t->SetBranchAddress("MuonTight", &nMUON);
  TBranch *muon_pt = t->GetBranch("MuonTight.PT");
  muon_pt->SetAddress(&MUON_PT);
  TBranch *muon_eta = t->GetBranch("MuonTight.Eta");
  muon_eta->SetAddress(&MUON_ETA);
  TBranch *muon_phi = t->GetBranch("MuonTight.Phi");
  muon_phi->SetAddress(&MUON_PHI);
  TBranch *muon_sumpt = t->GetBranch("MuonTight.SumPtCharged");
  muon_sumpt->SetAddress(&MUON_SUMPT);

  
  ofstream fileout("cutinfo2.dat");
  /////////////////////////////////
  // Start Cut Analysis          //
  /////////////////////////////////
  Long64_t cut_flow[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for(Long64_t i = 0; i < nevent; i++){
    cut_flow[0]++;
    met->GetEntry(i);
    met_phi->GetEntry(i);
    jet_pt->GetEntry(i);
    jet_phi->GetEntry(i);
    jet_eta->GetEntry(i);
    jet_m->GetEntry(i);
    jet_btag->GetEntry(i);
    e_pt->GetEntry(i);
    e_eta->GetEntry(i);
    e_phi->GetEntry(i);
    e_sumpt->GetEntry(i);
    e_hoe->GetEntry(i);
    muon_pt->GetEntry(i);
    muon_eta->GetEntry(i);
    muon_phi->GetEntry(i);
    muon_sumpt->GetEntry(i);

    ///////////////////////
    //(1) Ol + 4-5 jets  //
    ///////////////////////
    Int_t num_Lep = 0;
    Int_t num_JET = 0;
    
    for(int j=0; j<nJET; j++){
      if(JET_PT[j] > 30 && abs(JET_ETA[j]) < 2.4) num_JET++;
    }
    for(int j=0; j<nMUON; j++){
      if(MUON_PT[j] > 10 && abs(MUON_ETA[j]) < 2.4) num_Lep++;
    }
    for(int j=0; j<nElec; j++){
      if(ELEC_PT[j] > 10 && abs(ELEC_ETA[j]) < 2.5) num_Lep++;
    }
    if(num_JET<4 || num_JET>5 || num_Lep>0) continue;
    cut_flow[1]++;

    //////////////////////
    //(2) N(b,T) >= 2   //
    //////////////////////
    Int_t Nb_T = 0;
    Int_t Nb_M = 0;
    Int_t Nb_L = 0;

    for(int j=0; j<nJET; j++){
      if(Bool_t(JET_BTAG[j] &(1<<2)) && JET_PT[j]>30 && abs(JET_ETA[j])<2.4){
	Nb_T++;
	Nb_M++;
	Nb_L++;
      }
      else if(Bool_t(JET_BTAG[j] *(1<<1)) && JET_PT[j]>30 && abs(JET_ETA[j])<2.4){
	Nb_M++;
	Nb_L++;
      }
      else if(Bool_t(JET_BTAG[j] *(1<<0)) && JET_PT[j]>30 && abs(JET_ETA[j])<2.4){
	Nb_L++;
      }
      else continue;
    }
    if(Nb_T<2) continue;
    cut_flow[2]++;
   
    //////////////////////       
    //(3) MET > 150    //                                       
    //////////////////////                                       
    if(MET[0]<150) continue;
    cut_flow[3]++;
    
    //////////////////////   
    //(4) Track Veto    //
    //////////////////////
    Int_t iso_flag = 0;
    for(int j=0; j<nElec; j++){
      if( ELEC_HOE[j]<1 && (ELEC_PT[j]<5. || ELEC_PT[j]/ELEC_SUMPT[j]<5.) ) continue;
      if( ELEC_HOE[j]<1 && (ELEC_PT[j]<10. || ELEC_PT[j]/ELEC_SUMPT[j]<10.) ) continue;
      Double_t deltaphi = MET_PHI[0] - ELEC_PHI[j];
      if(deltaphi<0) deltaphi += 4*TMath::Pi();
      if(deltaphi>2.*TMath::Pi()) deltaphi -= 2.*TMath::Pi();
      if(deltaphi>TMath::Pi()) deltaphi = 2.*TMath::Pi() - deltaphi;
      if(TMath::Sqrt(2.*MET[0]*ELEC_PT[j]*(1.-TMath::Cos(deltaphi))) < 100) iso_flag=1; 
    }
    for(int j=0; j<nMUON; j++){
      if( MUON_PT[j]<5. || MUON_PT[j]/MUON_SUMPT[j]<5. ) continue;
      Double_t deltaphi = MET_PHI[0] - MUON_PHI[j];
      if(deltaphi<0) deltaphi += 4*TMath::Pi();
      if(deltaphi>2*TMath::Pi()) deltaphi -= 2*TMath::Pi();
      if(deltaphi>TMath::Pi()) deltaphi = 2*TMath::Pi() - deltaphi;
      if(TMath::Sqrt(2.*MET[0]*MUON_PT[j]*(1.-TMath::Cos(deltaphi))) < 100) iso_flag=1;
    }
    if(iso_flag==1) continue;
    cut_flow[4]++;
   
    //////////////////////
    //(5) Delta Phi     //
    //////////////////////
    Double_t deltaphi[4];
    
    for(int j=0; j<4; j++){
      deltaphi[j] = MET_PHI[0] - JET_PHI[j];
      if(deltaphi[j]<0) deltaphi[j] += 4.*TMath::Pi();
      if(deltaphi[j]>2.*TMath::Pi()) deltaphi[j] -= 2.*TMath::Pi();
      if(deltaphi[j]>TMath::Pi()) deltaphi[j] = 2.*TMath::Pi() - deltaphi[j];
    }
    if(deltaphi[0]<0.5 || deltaphi[1]<0.5 || deltaphi[2]<0.3 || deltaphi[3]<0.3) continue;
    cut_flow[5]++;

    //////////////////////////////////
    // Reco Higgs pair              //
    //////////////////////////////////
    TLorentzVector b_vector[4];
    Int_t b_vector_index = 0;
    // b-tag indentification
    if(Nb_T>=2 && Nb_M>=3 && Nb_L>=4){
      for(int j=0; j<nJET; j++){
	if(b_vector_index==4) continue;
	if( Bool_t(JET_BTAG[j] & (1<<2)) ){
	  b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
	  b_vector_index++;
	}
      }
      for(int j=0; j<nJET; j++){
	if(b_vector_index==4) continue;
	if( Bool_t(JET_BTAG[j] & (1<<1)) ){
	  b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
	  b_vector_index++;
	}
      }
      for(int j=0; j<nJET; j++){
	if(b_vector_index==4) continue;
	if( Bool_t(JET_BTAG[j] & (1<<0)) ){
	  b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
          b_vector_index++;
	}
      }
    }
    else if(Nb_T>=2){
      for(int j=0; j<nJET; j++){
	if(b_vector_index==4) continue;
	if( Bool_t(JET_BTAG[j] & (1<<2)) ){
	  b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
          b_vector_index++;
	}
      }
      for(int j=0; j<nJET; j++){
	if(b_vector_index==4) continue;
	if( Bool_t(JET_BTAG[j] & (1<<1)) ){
          b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
          b_vector_index++;
        }
      }
      for(int j=0; j<nJET; j++){
	if(b_vector_index==4) continue;
	if( Bool_t(JET_BTAG[j] & (1<<0)) ){
          b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
          b_vector_index++;
        }
      }
    }
    // Higgs Reco
    TLorentzVector aux1, aux2;
    Double_t delta_m[3], mean_m[3], delta_Ra[3], delta_Rb[3];
    Int_t z1=-1, z2=-1, z3=-1;
    for(int j=0; j<3; j++){
      z1 = round(j + 1);
      z2 = (z1+3) % 3 + 1;
      z3 = (z1+1) % 3 + 1;
      
      cout << 0 <<":"<<z1<<":"<<z2<<":"<<z3<<endl;
      aux1=b_vector[0]+b_vector[z1];
      aux2=b_vector[z2]+b_vector[z3];
      delta_m[j]=abs(aux1.M()-aux2.M());
      cout << abs(aux1.M()-aux2.M())<<endl;
      mean_m[j]=(aux1.M()+aux2.M())/2.0;
      delta_Ra[j]=abs(b_vector[0].DeltaR(b_vector[z1]));
      delta_Rb[j]=abs(b_vector[z2].DeltaR(b_vector[z3]));
      fileout << "delta_m "<<  j  << ":  " << delta_m[j] << endl;
      fileout << "mean_m "<<  j  << ":  " << mean_m[j] << endl;
    }
    
    Double_t mini = -1;
    Int_t index_min = -1;
    if((delta_m[0] - delta_m[1])>=0){ 
      mini = delta_m[1]; index_min = 1; 
    }
    else { 
      mini = delta_m[0]; index_min = 0; 
    }
    if((mini - delta_m[2])>=0){
      index_min = 2;
    }
    else {
      index_min = index_min;
    }   

    //////////////////////////////////
    // End Reco                     //
    //////////////////////////////////

    //////////////////////
    //(6) Delta mass    //
    ////////////////////// 
    if(delta_m[index_min] < 40) continue;
    cut_flow[6]++;

    //////////////////////
    //(7) Delta Rmax    //
    //////////////////////
    Double_t dRmax = -10;
    dRmax = TMath::Max(delta_Ra[index_min], delta_Rb[index_min]);
    if(dRmax>2.2) continue;
    cut_flow[7]++;

    //////////////////////
    //(8) Higg inv mass //
    //////////////////////
    h_m_average->Fill(mean_m[index_min]);
    if(mean_m[index_min]<100 || mean_m[index_min]>140) continue;
    cut_flow[8]++;
    
    //////////////////////
    //(9) 3b+4b         //
    //////////////////////
    if(Nb_M<3 || Nb_L<3) continue;
    cut_flow[9]++;
    //////////////////////
    //(10)  4b           //
    //////////////////////
    if(Nb_M<3 || Nb_L<4) continue;
    cut_flow[10]++;
    
    //////////////////////   
    //(11) MET > 200    //    
    //////////////////////
    if(MET[0]<200) continue;
    cut_flow[11]++;

    ////////////////////// 
    //(12) MET > 300    //
    ////////////////////// 
    if(MET[0]<300) continue;
    cut_flow[12]++;

    ////////////////////// 
    //(13) MET > 450    //
    ////////////////////// 
    if(MET[0]<450) continue;
    cut_flow[13]++;
 
  }
  /////////////////////                                                       
  // output cut info //                                                       
  ///////////////////// 
  cout << "store cut info ..." << endl;
  
 // ofstream fileout("cutinfo2.dat");
  fileout << "No Selection : "    << " " << cut_flow[0] << endl;
  fileout << "Ol + (4~5) Jets : " << " " << cut_flow[1] << endl;
  fileout << "N_{b,T} >= 2 : "    << " " << cut_flow[2] << endl;
  fileout << "MET > 150 : "       << " " << cut_flow[3] << endl;
  fileout << "Track Veto : "      << " " << cut_flow[4] << endl;
  fileout << "#Delta#phi_{1,2}>0.5,#Delta#phi_{3,4}>0.3 : " << " " << cut_flow[5] << endl;
  fileout << "|#Deltam|<40GeV : "    << " " << cut_flow[6] << endl;
  fileout << "#DeltaR_{max}<2.2 : "  << " " << cut_flow[7] << endl;
  fileout << "100<#bar{m}<140GeV : " << " " << cut_flow[8] << endl;
  fileout << "3b+4b : " << " " << cut_flow[9] << endl;
  fileout << "4b : "    << " " << cut_flow[10] << endl;
  fileout << "MET > 200 : " << " " << cut_flow[11] << endl;
  fileout << "MET > 300 : " << " " << cut_flow[12] << endl;
  fileout << "MET > 450 : " << " " << cut_flow[13] << endl;
  fileout.close();

  return 1;
}
