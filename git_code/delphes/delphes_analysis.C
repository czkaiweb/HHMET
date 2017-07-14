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

int delphes_analysis(/*char* string*/){

//  TString finname(string);
//  TString foutname(string);
  
//  finname = "/data2/weikai/calchep/delphes_results/HHMET/"+finname+"/delphes/myprocess_cms8tev.root";
//  foutname ="analysis_results/"+foutname+"_results.root";
//  TFile *fin=new TFile(finname,"READONLY");
//  TFile *fout=new TFile(foutname,"RECREATE");

//  TFile *fin=new TFile("/data2/weikai/calchep/delphes_results/HHMET/gYqq0.05MY600MF200/delphes/myprocess_cms8tev.root","READONLY");
  TFile *fin=new TFile("/data2/weikai/calchep/delphes_results/HHMET/ttbar_cms/delphes/myprocess_cms8tev.root","READONLY");
//  TFile *fout = new TFile("test_signal.root","RECREATE");
  TFile *fout = new TFile("test_bkg.root","RECREATE");

//----------------SETUP HISTOGRAM------------------
  TH1D *info=new TH1D("info","info",14,0,14);
  TH1D *h_leading_pt = new TH1D("h_leading_pt","leading_pt",10,0,500);
  TH1D *h_subleading_pt = new TH1D("h_subleading_pt","subleading_pt",10,0,500);
  TH1D *h_m_diff = new TH1D("h_m_difference","m_differnce",12,0,120);
  TH1D *h_m_average = new TH1D("h_m_average","m_average",20,0,200);
  TH1D *h_MET = new TH1D("h_MET","h_MET",12,0,600);
  TH1D *h_dRmax = new TH1D("h_dRmax","dRmax",20,0,4);

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
  Float_t JET_PT[20];
  Float_t JET_ETA[20];
  Float_t JET_PHI[20];
  Float_t JET_M[20];
  UInt_t JET_BTAG[20];
  TLorentzVector all_vector[20];
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
  
  Int_t nElec=-1;
  Float_t ELEC_PT[20];
  Float_t ELEC_ETA[20];
  Float_t ELEC_PHI[20];
  Float_t ELEC_SUMPT[20];
  Float_t ELEC_HOE[20];
  t->SetBranchAddress("Electron",&nElec);
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
  cout << "electron finish"<<endl; 
  
  Int_t nMUON=-1;
  Float_t MUON_PT[10];
  Float_t MUON_ETA[10];
  Float_t MUON_PHI[10];
  Float_t MUON_SUMPT[10];
  Float_t MUON_HOE[10];
  t->SetBranchAddress("MuonTight",&nMUON);
  TBranch *muon_pt = t->GetBranch("MuonTight.PT");
  muon_pt->SetAddress(&MUON_PT);
  TBranch *muon_eta = t->GetBranch("MuonTight.Eta");
  muon_eta->SetAddress(&MUON_ETA);
  TBranch *muon_phi = t->GetBranch("MuonTight.Phi");
  muon_phi->SetAddress(&MUON_PHI);
  TBranch *muon_sumpt = t->GetBranch("MuonTight.SumPtCharged");
  muon_sumpt->SetAddress(&MUON_SUMPT);
  //cout << "muon sumpt"<<endl;
  //TBranch *muon_hoe = t->GetBranch("MuonTight.EhadOverEem");
  //muon_hoe->SetAddress(&MUON_HOE);
  cout << "muon finish"<<endl;

  Int_t nTRACK = -1;
  Int_t MaxTrkNum=1000;
  Float_t TRACK_PT[MaxTrkNum];
  Float_t TRACK_ETA[MaxTrkNum];
  Float_t TRACK_PHI[MaxTrkNum];
  t->SetBranchAddress("EFlowTrack", &nTRACK);
  TBranch *track_pt = t->GetBranch("EFlowTrack.PT");
  track_pt->SetAddress(&TRACK_PT);
  TBranch *track_eta = t->GetBranch("EFlowTrack.Eta");
  track_eta->SetAddress(&TRACK_ETA);
  TBranch *track_phi = t->GetBranch("EFlowTrack.Phi");
  track_phi->SetAddress(&TRACK_PHI);



  Long64_t cut_flow[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  cout << "start loop"<<endl;

  for(Long64_t i=0;i<nevent;i++){
  //  if(i%100==0) cout << i <<endl;
    cut_flow[0]++;
//    t->GetEntry(i);
    met->GetEntry(i);
    met_phi->GetEntry(i);
    jet_pt->GetEntry(i);
    jet_phi->GetEntry(i);
    jet_eta->GetEntry(i);
    jet_m->GetEntry(i);
    jet_btag->GetEntry(i);
    //cout << nElec <<endl;
    e_pt->GetEntry(i);
   // cout<<ELEC_PT[0]<<endl;
    e_eta->GetEntry(i);
   // cout<<ELEC_ETA[0]<<endl;
   // cout << "checkpoint 0" << endl;
   //   cout << i <<endl;
      e_phi->GetEntry(i);
  //    cout << "checkpoint 0'" << endl;
      e_sumpt->GetEntry(i);
      e_hoe->GetEntry(i);
    
    muon_pt->GetEntry(i);
    muon_eta->GetEntry(i);
    muon_phi->GetEntry(i);
    muon_sumpt->GetEntry(i);
//    muon_hoe->GetEntry(i);
    track_pt->GetEntry(i);
    track_eta->GetEntry(i);
    track_phi->GetEntry(i);
//    cout << nMUON <<endl;
//    cout << "nJET:"<<nJET<<endl;
//----------------ol+(4-5)jets-----------------------
   Int_t nL_d=0;
   Int_t nJET_d=0;
    
   for(int l=0;l<nJET;l++){
     if(JET_PT[l]>30 && fabs(JET_ETA[l])<2.4) nJET_d++;
   }
   for(int l=0;l<nElec;l++){
     if(ELEC_PT[l]>10 && fabs(ELEC_ETA[l])<2.5) nL_d++;
   }
   for(int l=0;l<nMUON;l++){
     if(MUON_PT[l]>10 && fabs(MUON_ETA[l])<2.4) nL_d++;
   }
   //cout << nJET_d << "  "<< nJET <<endl;
   //cout << nL_d  <<endl;
   if(nL_d>0 || nJET_d<4 || nJET_d>5) continue;
   cut_flow[1]++;     // cut: 0l 4-5jets
//   cout << "checkpoint 1" << endl;

//---------------Nb,T>2---------------------------
   Int_t Nb_T=0;
   Int_t Nb_M=0;
   Int_t Nb_L=0;

//----------------------------count btagging (old)-----------------------------
/*   for(int l=0;l<nJET;l++){
     if(nJET>20) break;
     if(Bool_t(JET_BTAG[l] &(1<<2)) && JET_PT[l]>30 && fabs(JET_ETA[l])<2.4){
       Nb_T++;
       Nb_M++;
       Nb_L++;
     }
     else if(Bool_t(JET_BTAG[l] &(1<<1)) && JET_PT[l]>30 && fabs(JET_ETA[l])<2.4){
       Nb_M++;
       Nb_L++;
     }
     else if(Bool_t(JET_BTAG[l] &(1<<0)) && JET_PT[l]>30 && fabs(JET_ETA[l])<2.4){
       Nb_L++;
     }
     else continue;
   }*/  //btaging count OLD

//-----------------------count b-tagging-----------
    vector<int> tightBindex,  mediumBindex, looseBindex, otherIndex;
    for(int j=0; j<nJET; j++){
      if(j>=20) break;
      if(JET_PT[j]>30 && fabs(JET_ETA[j])<2.4){
        if( Bool_t(JET_BTAG[j] & (1<<2)) ) tightBindex.push_back(j);
        else if( Bool_t(JET_BTAG[j] & (1<<1)) ) mediumBindex.push_back(j);
        else if( Bool_t(JET_BTAG[j] & (1<<0)) ) looseBindex.push_back(j);
	else otherIndex.push_back(j);
      }
    }
    Nb_T = tightBindex.size();
    Nb_M = Nb_T + mediumBindex.size();
    Nb_L = Nb_M + looseBindex.size();
   //cout << "Nb,T:"<<Nb_T<<endl;
   if(Nb_T<2) continue;
   cut_flow[2]++;     // cut: Nb_T>=2

//   cout << "checkpoint 2" << endl;
    h_leading_pt->Fill(JET_PT[0]);
    h_subleading_pt->Fill(JET_PT[1]);

//---------------Pt missing-------------------------
   if(nMET!=1) continue;
   if(Nb_L>=4) h_MET->Fill(MET[0]);     
   if(MET[0]<150) continue;
   //cout<< MET[0] << endl;
   cut_flow[3]++;

   /* for(int j=0;j<nJET;j++){
      cout << JET_ETA[j] << endl;
      cout << JET_PHI[j] << endl;
      cout << JET_M[j] << endl;
      cout << " "<<endl;
    }*/
//---------------track veto-------------------------
//
    Int_t iso_flag=0;
    Int_t etrack_num=0;
    for(int l=0;l<nElec;l++){
      if(nElec>20) break;
      if( ELEC_PT[l]<5 || ELEC_PT[l]/ELEC_SUMPT[l]<5 ) continue;
      if(ELEC_HOE[l]>1 && ( ELEC_PT[l]<10 || ELEC_PT[l]/ELEC_SUMPT[l]<10)) continue;
      Double_t deltaphi = ELEC_PHI[l]-MET_PHI[0];
      if(deltaphi<0) deltaphi+=4*TMath::Pi();
      if(deltaphi>2*TMath::Pi()) deltaphi-=2*TMath::Pi();
      if(deltaphi>TMath::Pi()) deltaphi=2*TMath::Pi()-deltaphi;
      if(TMath::Sqrt(2.*MET[0]*ELEC_PT[l]*(1.-TMath::Cos(deltaphi))) < 100.) iso_flag=1;
    }
    for(int l=0;l<nMUON;l++){
      if(MUON_PT[l]<5 || MUON_PT[l]/MUON_SUMPT[l]<5) continue;
      //if(MUON_HOE[l]>1 && ( MUON_PT[l]<10 || MUON_PT[l]/MUON_SUMPT[l]>0.1)) continue;
      Double_t deltaphi = MUON_PHI[l]-MET_PHI[0];
      if(deltaphi<0) deltaphi+=4*TMath::Pi();
      if(deltaphi>2*TMath::Pi()) deltaphi-=2*TMath::Pi();
      if(deltaphi>TMath::Pi()) deltaphi=2*TMath::Pi()-deltaphi;
      if(TMath::Sqrt(2.*MET[0]*MUON_PT[l]*(1.-TMath::Cos(deltaphi))) < 100.) iso_flag=1;
    }
    TVector3 trkpj, trkpk;
    for(int j=0; j< nTRACK; j++){
      if( TRACK_PT[j]<10. ) continue;
      Double_t TRACK_SUMPT=0;
      trkpj.SetPtEtaPhi(TRACK_PT[j], TRACK_ETA[j], TRACK_PHI[j]);
      for(int k=0; k<nTRACK; k++){
        if(k==j) continue;
        if(TRACK_PT[k]<0.5 || fabs(TRACK_ETA[k]) >2.5) continue;
        trkpk.SetPtEtaPhi(TRACK_PT[k], TRACK_ETA[k], TRACK_PHI[k]);
        if( trkpj.DeltaR(trkpk) < 0.3)  TRACK_SUMPT+=TRACK_PT[k];
      }
      if(TRACK_SUMPT/TRACK_PT[j] > 0.1 ) continue;
      Double_t deltaphi = MET_PHI[0] - TRACK_PHI[j];
      if(deltaphi<0) deltaphi += 4.*TMath::Pi();
      if(deltaphi>2*TMath::Pi()) deltaphi -= 2.*TMath::Pi();
      if(deltaphi>TMath::Pi()) deltaphi = 2.*TMath::Pi() - deltaphi;
      if(TMath::Sqrt(2.*MET[0]*TRACK_PT[j]*(1.-TMath::Cos(deltaphi))) < 100.) iso_flag=1;
    }



    if(iso_flag==1) continue; 
    cut_flow[4]++;
//   cout << "checkpoint 3" << endl;

//--------------\delta\phi1,2,3,4
    Int_t Index[4]={-1,-1,-1,-1};
    Int_t JET_PT_CLONE[20];
    for(int l=0;l<nJET;l++){
      JET_PT_CLONE[l]=JET_PT[l];
    }
    for(int loop=0;loop<4;loop++){
      Double_t leading_pt=0;
      Int_t leading_index=-1;
      for(int l=0;l<nJET;l++){
  	if(JET_PT_CLONE[l]>leading_pt&& JET_PT_CLONE[l]>30 && fabs(JET_ETA[l])<2.4){
          leading_pt=JET_PT_CLONE[l];
          leading_index=l;
        }
      }
      JET_PT_CLONE[leading_index]=-10;
      Index[loop]=leading_index;
    }// mark the index according to pt, in fact delphes output has already put them in order.
    Double_t deltaphi[4];

    for(int l=0;l<4;l++){
      deltaphi[l]=MET_PHI[0]-JET_PHI[Index[l]];
      if(deltaphi[l]<0) deltaphi[l]+=4*TMath::Pi();
      if(deltaphi[l]>2*TMath::Pi()) deltaphi[l]-=2*TMath::Pi();
      if(deltaphi[l]>TMath::Pi()) deltaphi[l]=2*TMath::Pi()-deltaphi[l];
    }//calcualte deltaphi
    if(deltaphi[0]<0.5 || deltaphi[1]<0.5 || deltaphi[2]<0.3 || deltaphi[3]<0.3) continue;
//    cout <<i<<":"<< deltaphi[0]<<endl;  
    cut_flow[5]++;

//-------------m_difference---------------------------

///////////////////////////////////////////////////////////////////////////////
////locate the index of b-jet in order (tight/medium/loose)     (new)      ////
///////////////////////////////////////////////////////////////////////////////
    TLorentzVector b_vector[4];
    vector<int> Bindex;
    for(Int_t j=0; j< tightBindex.size(); j++) Bindex.push_back( tightBindex.at(j));
    for(Int_t j=0; j< mediumBindex.size(); j++) Bindex.push_back( mediumBindex.at(j));
    for(Int_t j=0; j< looseBindex.size(); j++) Bindex.push_back( looseBindex.at(j));
    for(Int_t j=0; j< otherIndex.size(); j++) Bindex.push_back(otherIndex.at(j));
    for(Int_t k=0;k< Bindex.size();k++){
      if(k>=4) break;
      Int_t j = Bindex.at(k);
      if(j<0 || j>=nJET ){
	cout << "Warning: b-tagging jet out-of range" << endl;
	};
      b_vector[k].SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
    } 

///////////////////////////////////////////////////////////////////////////////
////locate the index of b-jet in order (tight/medium/loose)     (old)      ////
///////////////////////////////////////////////////////////////////////////////


/*    Int_t b_vector_index=0;
    if(Nb_T>=2&&Nb_M>=3&&Nb_L>=4){
      for(int l=0;l<nJET;l++){
        if(b_vector_index==4) continue; //in case more than 5 b-tagged jets
        if(Bool_t(JET_BTAG[l] &(1<<2))){
          b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[l],JET_ETA[l],JET_PHI[l],JET_M[l]);
          b_vector_index++;
        }
      }
      for(int l=0;l<nJET;l++){
        if(b_vector_index==4) continue; //in case more than 5 b-tagged jets
        if(Bool_t(JET_BTAG[l] &(1<<1)) && !Bool_t(JET_BTAG[l] &(1<<2))){
          b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[l],JET_ETA[l],JET_PHI[l],JET_M[l]);
          b_vector_index++;
        }
      }
      for(int l=0;l<nJET;l++){
        if(b_vector_index==4) continue; //in case more than 5 b-tagged jets
        if(Bool_t(JET_BTAG[l] &(1<<0))&&!Bool_t(JET_BTAG[l] &(1<<2))&& !Bool_t(JET_BTAG[l] &(1<<1))){
          b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[l],JET_ETA[l],JET_PHI[l],JET_M[l]);
          b_vector_index++;
        }
      }
    }
    else if(Nb_T>=2){
      for(int l=0;l<nJET;l++){
        if(Bool_t(JET_BTAG[l] &(1<<2))){
          b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[l],JET_ETA[l],JET_PHI[l],JET_M[l]);
          b_vector_index++;
        }
      }
      for(int l=0;l<nJET;l++){
        if(b_vector_index==4) continue;
        if(Bool_t(JET_BTAG[l] &(1<<1))&& !Bool_t(JET_BTAG[l] &(1<<2))){
          b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[l],JET_ETA[l],JET_PHI[l],JET_M[l]);
          b_vector_index++;
        }
      }
      for(int l=0;l<nJET;l++){
	if(b_vector_index==4) continue;	
        if(Bool_t(JET_BTAG[l] &(1<<0))  && !Bool_t(JET_BTAG[l] &(1<<2))&& !Bool_t(JET_BTAG[l] &(1<<1))){
          b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[l],JET_ETA[l],JET_PHI[l],JET_M[l]);
          b_vector_index++;
        }
      }
      for(int l=0;l<nJET;l++){
        if(b_vector_index==4) continue; //in case more than 5 b-tagged jets
        if(Bool_t(JET_BTAG[l] &(1<<0))==0){
          b_vector[b_vector_index].SetPtEtaPhiM(JET_PT[l],JET_ETA[l],JET_PHI[l],JET_M[l]);
          b_vector_index++;
        }
      } 
    }
*/

//   cout << "checkpoint 4" << endl;
////------------------calculate the invariant masses-------------------------
    cout <<"jet0:" <<b_vector[0].Pt()<<endl;
    cout <<"jet1:" <<b_vector[1].Pt()<<endl;
    cout <<"jet2:" <<b_vector[2].Pt()<<endl;
    cout <<"jet3:" <<b_vector[3].Pt()<<endl;

    TLorentzVector aux1,aux2;
    Double_t delta_m[3][2]; 
    Double_t delta_R[3][2]; 
    aux1=b_vector[0]+b_vector[1];
    aux2=b_vector[2]+b_vector[3];
    delta_m[0][0]=fabs(aux1.M()-aux2.M());
    delta_m[0][1]=(aux1.M()+aux2.M())/2.0;
    delta_R[0][0]=fabs(b_vector[0].DeltaR(b_vector[1]));
    delta_R[0][1]=fabs(b_vector[2].DeltaR(b_vector[3]));

    aux1=b_vector[0]+b_vector[2];
    aux2=b_vector[1]+b_vector[3];
    delta_m[1][0]=fabs(aux1.M()-aux2.M());
    delta_m[1][1]=(aux1.M()+aux2.M())/2.0;
    delta_R[1][0]=fabs(b_vector[0].DeltaR(b_vector[2]));
    delta_R[1][1]=fabs(b_vector[1].DeltaR(b_vector[3]));

    aux1=b_vector[0]+b_vector[3];
    aux2=b_vector[2]+b_vector[1];
    delta_m[2][0]=fabs(aux1.M()-aux2.M());
    delta_m[2][1]=(aux1.M()+aux2.M())/2.0;
    delta_R[2][0]=fabs(b_vector[0].DeltaR(b_vector[3]));
    delta_R[2][1]=fabs(b_vector[2].DeltaR(b_vector[1]));
    Double_t deltam=-10;
    Double_t averagem=-10;
    Double_t dRmax=-10;
    if(delta_m[0][0]<delta_m[1][0] && delta_m[0][0]<delta_m[2][0] &&delta_m[0][0]<40){
      deltam=delta_m[0][0];
      averagem=delta_m[0][1];
      dRmax=TMath::Max(delta_R[0][0],delta_R[0][1]);
    }
    if(delta_m[1][0]<delta_m[0][0] && delta_m[1][0]<delta_m[2][0] &&delta_m[1][0]<40){
      deltam=delta_m[1][0];
      averagem=delta_m[1][1];
      dRmax=TMath::Max(delta_R[1][0],delta_R[1][1]);
    }
    if(delta_m[2][0]<delta_m[1][0] && delta_m[2][0]<delta_m[0][0] &&delta_m[2][0]<40){
      deltam=delta_m[2][0];
      averagem=delta_m[2][1];
      dRmax=TMath::Max(delta_R[2][0],delta_R[2][1]);
    }
    if(Nb_L>=4) h_m_diff->Fill(TMath::Min(delta_m[0][0],TMath::Min(delta_m[0][0],delta_m[2][0])));
    if(deltam<0) continue;
    cut_flow[6]++;
//   cout << "checkpoint 5" << endl;
//    cout << "dm:"<<deltam<<"   am:"<<averagem <<"dRmax:"<< dRmax <<endl;
//---------------deltaRmax---------------------------
    if(Nb_L>=4) h_dRmax->Fill(dRmax);    
    if(dRmax>2.2) continue;
    cut_flow[7]++;
    
//------------- baseline---------------------------
    //h_leading_pt->Fill(JET_PT[0]);
    //h_subleading_pt->Fill(JET_PT[1]);
//    h_MET->Fill(MET[0]);    
    if(Nb_L>=4) h_m_average->Fill(averagem);
//--------------averagem-----------------------------
    if(averagem<100 || averagem>140)  continue;
    cut_flow[8]++;  
    cout << "dm:"<<deltam<<"   am:"<<averagem <<"dRmax:"<< dRmax <<endl;

//   cout << "checkpoint 6" << endl;
//-------------3b+4b--------------------------------
    if(Nb_M<3 || Nb_L<3) continue;
    cut_flow[9]++; 
     
//-------------4b----------------------------------
    if(Nb_M<3 || Nb_L<4) continue;
    cut_flow[10]++;
    cout<< cut_flow[10]<<endl;;
//-------------pt missing >200--------------------------
    if(MET[0]<200) continue;
    cut_flow[11]++;
    cout <<"misssing pt > 200:" <<MET[0]<<endl;
//-------------pt missing >300--------------------------
    if(MET[0]<300) continue;
    cut_flow[12]++;
    cout <<"misssing pt > 300:" <<MET[0]<<endl;
//-------------pt missing >200--------------------------
    if(MET[0]<450) continue;
    cut_flow[13]++;
    cout <<"misssing pt > 450:" <<MET[0]<<endl;
}

info->Fill("No selection",cut_flow[0]);
info->Fill("0j+(4,5)jet",cut_flow[1]);
info->Fill("N_{b,T}>=2",cut_flow[2]);
info->Fill("p^{miss}_{T}>150GeV",cut_flow[3]);
info->Fill("Track veto",cut_flow[4]);
info->Fill("#Delta#phi_{1,2}>0.5,#Delta#phi_{3,4}>0.3",cut_flow[5]);
info->Fill("|#Deltam|<40GeV",cut_flow[6]);
info->Fill("#DeltaR_{max}<2.2",cut_flow[7]);
info->Fill("100<#bar{m}<140GeV",cut_flow[8]);
info->Fill("3b+4b",cut_flow[9]);
info->Fill("4b",cut_flow[10]);
info->Fill("p^{miss}_{T}>200GeV",cut_flow[11]);
info->Fill("p^{miss}_{T}>300GeV",cut_flow[12]);
info->Fill("p^{miss}_{T}>450GeV",cut_flow[13]);

fout->cd();
info->Write();
h_leading_pt->Write();
h_subleading_pt->Write();
h_m_average->Write();
h_m_diff->Write();
h_MET->Write();
h_dRmax->Write();

  ofstream fileout("ttbar_cutflow_check.dat");
  fileout << "No Selection :                              "    << " " << cut_flow[0] << endl;
  fileout << "Ol + (4~5) Jets :                           "    << " " << cut_flow[1] << endl;
  fileout << "N_{b,T} >= 2 :                              "    << " " << cut_flow[2] << endl;
  fileout << "MET > 150 :                                 "    << " " << cut_flow[3] << endl;
  fileout << "Track Veto :                                "    << " " << cut_flow[4] << endl;
  fileout << "#Delta#phi_{1,2}>0.5,#Delta#phi_{3,4}>0.3 : "    << " " << cut_flow[5] << endl;
  fileout << "|#Deltam|<40GeV :                           "    << " " << cut_flow[6] << endl;
  fileout << "#DeltaR_{max}<2.2 :                         "    << " " << cut_flow[7] << endl;
  fileout << "100<#bar{m}<140GeV :                        "    << " " << cut_flow[8] << endl;
  fileout << "3b+4b :                                     "    << " " << cut_flow[9] << endl;
  fileout << "4b :                                        "    << " " << cut_flow[10] << endl;
  fileout << "MET > 200 :                                 "    << " " << cut_flow[11] << endl;
  fileout << "MET > 300 :                                 "    << " " << cut_flow[12] << endl;
  fileout << "MET > 450 :                                 "    << " " << cut_flow[13] << endl;
  fileout.close();


return 1;
}
