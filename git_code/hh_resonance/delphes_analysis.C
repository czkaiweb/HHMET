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
  
//  finname = "/data2/weikai/calchep/delphes_results/HHMET/"+finname+"/delphes/myprocess_atlas13tev.root";
//  foutname ="analysis_results/"+foutname+"_resolved_results.root";

//  TFile *fin=new TFile(finname,"READONLY");
//  TFile *fout=new TFile(foutname,"RECREATE");

//  TFile *fin=new TFile("/data2/weikai/calchep/delphes_results/HHMET/gYqq0.05MY600MF200/delphes/myprocess_cms8tev.root","READONLY");
//  TFile *fin=new TFile("/data2/weikai/calchep/delphes_results/HHMET/ttbar/delphes/myprocess_atlas13tev.root","READONLY");
  TFile *fout = new TFile("test_signal_tchain.root","RECREATE");
//  TFile *fout = new TFile("test_bkg.root","RECREATE");

  TChain *t = new TChain("Delphes");
  t->Add("/home/weikai/work/calchep/delphes_results/HHMET/gYqq0.5MY1000MF200/delphes/myprocess_atlas13tev.root");
  t->Add("/home/weikai/work/calchep/delphes_results/HHMET/gYqq0.5MY1000MF250/delphes/myprocess_atlas13tev.root");


//----------------SETUP HISTOGRAM------------------
  TH1D *info=new TH1D("info","info",14,0,14);
  TH1D *h_leading_pt = new TH1D("h_leading_pt","leading_pt",20,0,1000);
  TH1D *h_subleading_pt = new TH1D("h_subleading_pt","subleading_pt",20,0,1000);
  TH1D *h_m_lead = new TH1D("h_m_lead","m_lead",20,0,200);
  TH1D *h_m_subl = new TH1D("h_m_subl","m_subl",20,0,200);
  TH1D *h_m_4j = new TH1D("h_m_4j","m_4j",50,0,1500);
  TH1D *h_m_2h = new TH1D("h_m_2h","m_2h",50,0,1500);

//----------------END OF HISTOGRAM---------------------  
//  TTree *t = (TTree*)fin->Get("Delphes");
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
  Int_t           nJetAK10Trim = nJetAK10Trim;
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
  Int_t           nJetAK2Track = nJetAK2Track;
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

  Long64_t cut_flow[5]={0,0,0,0,0};
  cout << "start loop"<<endl;

  for(Long64_t i=0;i<nevent;i++){
    if(i%100==0) cout << i <<endl;
    cut_flow[0]++;
    //t->GetEntry(i);
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
////                          Resloved analysis                               ////
////                    based on ATLAS-CONF-2016-049                          ////
//////////////////////////////////////////////////////////////////////////////////
//---------------4 b-jet(R=0.4)---------------------------
   Int_t Nb_T=0;
//-----------------------count b-tagging-----------
    vector<int> tightBindex;
    for(int j=0; j<nJET; j++){
      if(j>=Jet_size) break;
      if(JET_PT[j]>30 && fabs(JET_ETA[j])<2.5){
        if( JET_BTAG[j] ) tightBindex.push_back(j);
      }
    }

    Nb_T = tightBindex.size();
   if(Nb_T<4) continue;
   cut_flow[1]++;     // cut: Nb_T>=2



//-------------deltaR_jj---------------------------

///////////////////////////////////////////////////////////////////////////////
////       locate the index of b-jet(tight)  and higgs pair  (new)         ////
///////////////////////////////////////////////////////////////////////////////
    TLorentzVector b_vector[4];
    vector<int> Bindex;
    for(Int_t j=0; j< tightBindex.size(); j++) Bindex.push_back( tightBindex.at(j));
    for(Int_t k=0;k< Bindex.size();k++){
      if(k>=4) break;
      Int_t j = Bindex.at(k);
      if(j<0 || j>=nJET ){
	cout << "Warning: b-tagging jet out-of range" << endl;
	};
      b_vector[k].SetPtEtaPhiM(JET_PT[j],JET_ETA[j],JET_PHI[j],JET_M[j]);
    } 

////------------------calculate the invariant masses-------------------------
//    cout <<"jet0:" <<b_vector[0].Pt()<<endl;
//    cout <<"jet1:" <<b_vector[1].Pt()<<endl;
//    cout <<"jet2:" <<b_vector[2].Pt()<<endl;
//    cout <<"jet3:" <<b_vector[3].Pt()<<endl;

    TLorentzVector invariant_mass,aux1,aux2;
    invariant_mass = b_vector[0]+b_vector[1]+b_vector[2]+b_vector[3];
    Double_t m_4j = invariant_mass.M();

    h_m_4j->Fill(m_4j);
//-----------------initialize the index-----------------------
    Int_t leading_index[3][2],subleading_index[3][2];
    for(Int_t index=0;index<3;index++){
	leading_index[index][0]=-1;
	leading_index[index][1]=-1;
	subleading_index[index][0]=-1;
	subleading_index[index][1]=-1;
    }

    Double_t delta_R_leading, delta_R_subleading; 
   
    Int_t deltaR_flag=0;
    Int_t j0=-1,j1=-1,j2=-1,j3=-1;
    Int_t l_index1=-1,l_index2=-1,sl_index1=-1,sl_index2=-1;
    for(Int_t pair_index=0;pair_index<3;pair_index++){
	j1 = pair_index+1;
	j2 = (j1+3)%3+1;
	j3 = (j1+1)%3+1;


        aux1=b_vector[0]+b_vector[j1];
        aux2=b_vector[j2]+b_vector[j3];
	if((b_vector[0].Pt()+b_vector[j1].Pt()) > (b_vector[j2].Pt()+b_vector[j3].Pt())){
		delta_R_leading=fabs(b_vector[0].DeltaR(b_vector[j1]));
	        delta_R_subleading=fabs(b_vector[j2].DeltaR(b_vector[j3]));
		l_index1=0;
		l_index2=j1;
		sl_index1=j2;
		sl_index2=j3;
	}
	else{
		delta_R_subleading=fabs(b_vector[0].DeltaR(b_vector[j1]));
                delta_R_leading=fabs(b_vector[j2].DeltaR(b_vector[j3]));
		sl_index1=0;
		sl_index2=j1;
		l_index1=j2;
		l_index2=j3;
	}


	if(m_4j<0) continue;
	if(m_4j<1250){
            if(delta_R_leading>(360/m_4j-0.5)){
		if(delta_R_leading< (655/m_4j+0.475)){
		    if(delta_R_subleading>(235/m_4j)){
			if(delta_R_subleading<(875/m_4j+0.35)){
                            leading_index[pair_index][0]=l_index1;			
                            leading_index[pair_index][1]=l_index2;			
                            subleading_index[pair_index][0]=sl_index1;			
                            subleading_index[pair_index][1]=sl_index2;			
			    deltaR_flag=1;
		         }
		    }
	        }
	    }
        }
        if(m_4j>1250){
	    if(delta_R_leading<1 && delta_R_leading>0){
	        if(delta_R_subleading<1 && delta_R_subleading>0){
                    leading_index[pair_index][0]=l_index1;			
                    leading_index[pair_index][1]=l_index2;			
                    subleading_index[pair_index][0]=sl_index1;			
                    subleading_index[pair_index][1]=sl_index2;			
	            deltaR_flag=1;
	        }
	    }
        }

    }
    if(deltaR_flag==0) continue;
    cut_flow[2]++;
//-----------------------choose the correct higgs pair---------------------
// using the distance to line (0,0)->(120,115) to select the higgs pair
    Double_t Dhh=-10,Dhh_min=-10;
    Double_t m_lead=-10,m_subl=-10;
    Int_t higgs_index[4]={-1,-1,-1,-1};
    for(Int_t pair_index=0;pair_index<3;pair_index++){
	j0=leading_index[pair_index][0];
	j1=leading_index[pair_index][1];
	j2=subleading_index[pair_index][0];
	j3=subleading_index[pair_index][1];
        if( j0<0 || j1<0 || j2<0 || j3<0 ) continue;
//        cout << "j0:" << j0 << "j1:" << j1 << " j2:" << j2 << " j3:" << j3 <<endl;
        TLorentzVector aux_leading,aux_subl;
	aux_leading = b_vector[j0]+b_vector[j1];
	aux_subl = b_vector[j2]+b_vector[j3];
	m_lead=aux_leading.M();
	m_subl=aux_subl.M();
	Dhh = TMath::Sqrt(m_lead*m_lead+m_subl*m_subl)*fabs(TMath::Sin(TMath::ATan(m_subl/m_lead)-TMath::ATan(115.0/120.0)));
        if(Dhh<Dhh_min||Dhh_min<0){
	    higgs_index[0]=j0;
	    higgs_index[1]=j1;
	    higgs_index[2]=j2;
	    higgs_index[3]=j3;
	}  
    }
    TLorentzVector h_leading_4v,h_subl_4v,hh_resonance;

    h_leading_4v = b_vector[higgs_index[0]]+ b_vector[higgs_index[1]];  
    h_subl_4v = b_vector[higgs_index[2]]+ b_vector[higgs_index[3]];  
    hh_resonance =h_leading_4v+h_subl_4v;
    h_leading_pt->Fill(h_leading_4v.Pt());
    h_subleading_pt->Fill(h_subl_4v.Pt());
    h_m_lead->Fill(h_leading_4v.M());
    h_m_subl->Fill(h_subl_4v.M());
    h_m_2h->Fill(hh_resonance.M());
////////////////////////////////////////////////////////////////////////
////                mass dependent deltaR(h,h)                      ////
////////////////////////////////////////////////////////////////////////
    if(h_leading_4v.Pt()<0.5*m_4j-90.0) continue;
    if(h_subl_4v.Pt()<0.33*m_4j-70) continue;
    if(fabs(h_leading_4v.Eta()-h_subl_4v.Eta())>1.1 && m_4j<850) continue;
    if(fabs(h_leading_4v.Eta()-h_subl_4v.Eta())>0.002*m_4j-0.6  && m_4j>850) continue;
    if(fabs(h_leading_4v.DeltaR(h_subl_4v))<1.5) continue;
    cut_flow[3]++;



////////////////////////////////////////////////////////////////////////
////         SR based on leading and sub-leading mass               ////
///////////////////////////////////////////////////////////////////////
    Double_t Xhh=-10;
    m_lead = h_leading_4v.M();
    m_subl = h_subl_4v.M();
    Xhh = TMath::Sqrt( ((m_lead-120)/(0.1*m_lead))*((m_lead-120)/(0.1*m_lead)) + ((m_subl-115)/(0.1*m_subl))*((m_subl-115)/(0.1*m_subl)) ) ;
    if(Xhh>1.6) continue;
    cut_flow[4]++;
} 


info->Fill("No selection",cut_flow[0]);
info->Fill("4 b-tagged jets",cut_flow[1]);
info->Fill("Delta R_{jj}",cut_flow[2]);
info->Fill("m_{4j} Denpendent Cuts,#Delta R_{hh}",cut_flow[3]);
info->Fill("Xhh",cut_flow[4]);

fout->cd();
info->Write();
h_leading_pt->Write();
h_subleading_pt->Write();
h_m_lead->Write();
h_m_subl->Write();
h_m_4j->Write();
h_m_2h->Write();

  ofstream fileout("resloved_cutflow_check_tchain_2file.dat");
  fileout << "No Selection :                              "    << " " << cut_flow[0] << endl;
  fileout << "4 b-tagged jets :                           "    << " " << cut_flow[1] << endl;
  fileout << "Delta R_{jj} :                              "    << " " << cut_flow[2] << endl;
  fileout << "m_{4j} Denpendent Cuts,#Delta R_{hh} :      "    << " " << cut_flow[3] << endl;
  fileout << "Xhh :                                       "    << " " << cut_flow[4] << endl;
  fileout.close();


return 1;
}
