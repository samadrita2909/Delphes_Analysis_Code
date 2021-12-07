
#include "TH1.h"
#include "TSystem.h"
#include <fstream>
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "external/ExRootAnalysis/ExRootConfReader.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "classes/DelphesClasses.h"

#include <iostream>

using namespace std;

//------------------------------------------------------------------------------
// USER DEFINED FUNCTIONS
//------------------------------------------------------------------------------

/* Delta_R */

double DR(double eta1,double phi1,double eta2,double phi2)
{
	double DPhi = TMath::Abs(phi1-phi2);
	double DEta = eta1 - eta2;
	if(DPhi > TMath::Pi())
	DPhi = TMath::TwoPi() - DPhi;
	return TMath::Sqrt(DEta*DEta + DPhi*DPhi);
}

/* Delta_Phi */

double DPhi(double phi1, double phi2)
{	
	double delPhi = TMath::Abs(phi1 - phi2);
	if(delPhi > TMath::Pi())
	delPhi = TMath::TwoPi() - delPhi;
	return delPhi;
}

/* transverse mass */

double MT(double phi1, double phi2 , double ETmiss, double PTrans)
{	
	double delPhi = TMath::Abs(phi1 - phi2);
	return TMath::Sqrt(2*PTrans*ETmiss*(1-cos(delPhi)));
}

double Mcol(double Mvis , double PTtauhad , double PTmet)
{	
	double xtau_vis= TMath::Abs(PTtauhad)/(TMath::Abs(PTtauhad) + TMath::Abs(PTmet));
	return Mvis/(TMath::Sqrt(xtau_vis));
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

struct MyPlots
{
  TH1 *fpTmu1;
  TH1 *fpTe1;
  TH1 *fpTj1;
  TH1 *fpTj2;
  TH1 *fpTj3;
  TH1 *fpTj4;
//  TH1 *ftransmom;
//  TH1 *fIMmutau;
  TH1 *fmissET;
  TH1 *fphil1j1;
  TH1 *fphil1j2;
  TH1 *fphiemiss;
  TH1 *fphimumiss;
  TH1 *fphil2miss;
  //TH1 *fetall;
  TH1 *fetal1;
  TH1 *fetaj1;
  TH1 *fetaj2;
  TH1 *fetaj3;
  TH1 *fetaj4;
  TH1 *fDRl1j1;
  TH1 *fDRl1j2;
  TH1 *fMj1j2;
  TH1 *fMTrans;
  TH1 *fMeff;
};

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, MyPlots *plots)
{
  THStack *stack;
  TLegend *legend;
  TPaveText *comment;




plots->fpTmu1 = result->AddHist1D(
    "pTmu1",  "pTmu1",
   "pTmu1", "Events",
    30, 0.0,300.0);


plots->fpTe1 = result->AddHist1D(
    "pTe1",  "pTe1",
   "pTe1", "Events",
    30, 0.0,300.0);

plots->fpTj1 = result->AddHist1D(
    "pTj1",  "pTj1",
   "pTj1", "Events",
    30, 0.0,300.0);

plots->fpTj2 = result->AddHist1D(
    "pTj2",  "pTj2",
   "pTj2", "Events",
    20, 0.0,200.0);


plots->fpTj3 = result->AddHist1D(
    "pTj3",  "pTj3",
   "pTj3", "Events",
    30, 0.0,300.0);

plots->fpTj4 = result->AddHist1D(
    "pTj4",  "pTj4",
   "pTj4", "Events",
    20, 0.0,200.0);


plots->fmissET = result->AddHist1D(
    "missET ",  "missET ",
   "missET", "Events",
    30, 0.0,300.0);


plots->fphil1j1 = result->AddHist1D(
    "phil1j1 ",  "phil1j1 ",
   "phil1j1", "Events",
    20, 0.0,4.0);

plots->fphil1j2 = result->AddHist1D(
    "phil1j2 ",  "phi1lj2 ",
   "phil1j2", "Events",
    20, 0.0,4.0);

plots->fphiemiss = result->AddHist1D(
    "phiemiss ",  "phiemiss ",
   "phiemiss", "Events",
    20, 0.0,4.0);



plots->fphimumiss = result->AddHist1D(
    "phimumiss ",  "phimumiss ",
   "phimumiss", "Events",
    20, 0.0,4.0);



plots->fetal1 = result->AddHist1D(
    "etal1 ",  "etal1 ",
   "etal1", "Events",
    50, -4.0,4.0);


plots->fetaj1 = result->AddHist1D(
    "etaj1 ",  "etaj1 ",
   "etaj1", "Events",
    50, -4.0,4.0);


plots->fetaj2 = result->AddHist1D(
    "etaj2 ",  "etaj2 ",
   "etalj2", "Events",
    50, -4.0,4.0);




plots->fetaj3 = result->AddHist1D(
    "etaj3 ",  "etaj3 ",
   "etaj3", "Events",
    50, -4.0,4.0);


plots->fetaj4 = result->AddHist1D(
    "etaj4 ",  "etaj4 ",
   "etalj4", "Events",
    50, -4.0,4.0);



plots->fDRl1j1 = result->AddHist1D(
    "DRl1j1 ",  "DRl1j1 ",
   "DRl1j1", "Events",
    50, 0.0,5.0);

plots->fDRl1j2 = result->AddHist1D(
    "DRl1j2 ",  "DRl1j2 ",
   "DRl1j2", "Events",
    50, 0.0,5.0);

plots->fMj1j2 = result->AddHist1D(
    "M_j1j2 ",  "M_j1j2 ",
   "M_j1j2", "Events",
    50, 0.0,500.0);


plots->fMeff = result->AddHist1D(
 "Meff", "Meff",
 "Meff", "Events",
 70, 0.0,700.0);

plots->fMTrans = result->AddHist1D(
 "MTrans", "MTrans",
 "MTrans", "Events",
 70, 0.0,700.0);



}




//------------------------------------------------------------------------------


void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
{
    
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");

  Long64_t allEntries = treeReader->GetEntries();  

  Long64_t entry;
  Int_t i, j, k, el_no, mu_no, jt_no, ph_no, b_no , tau_no;
  Int_t counter_bb, counter_bel, counter_bmu, counter_elb, counter_mub, counter_elel, counter_mumu;
//  Double_t  min_delPhi_mutau = 2.7;
  Double_t  max_delPhi_mumiss = 0.5;
  Double_t max_trans_masstauhad = 50.0 ;
  Double_t max_trans_massmu = 65.0 ;
  Double_t pT_l1, xmet, phi_l1met, eta_l1, MTrans, Meff, Mj1j2sq, Mj1j2,phi_j1met, phi_j2met, phi_j1j2, phi_l1j1, phi_l1j2, eta_j1, eta_j2, eta_j3, eta_j4,
 DRl1j1, DRl1j2, DRj1j2, jet1pT,jet2pT, jet3pT, jet4pT;
   
  Double_t met_max = 40.0;

  Int_t NEV_ok = 0;
  Int_t NEV_basic = 0;
  Int_t NEV_C1 = 0;
  Int_t NEV_C2 = 0;
  Int_t NEV_C3 = 0;
  Int_t NEV_C4 = 0;
  Int_t NEV_C5 = 0;
  Int_t NEV_C6 = 0;
  Int_t NEV_C7 = 0;
  Int_t NEV_C8 = 0;
  Int_t NEV_C = 0;

  Int_t tau_count = 0;
    
  Electron  *el[50];
  Muon      *mu[50];
  Photon    *ph[50];
  Jet       *jt[50];
  MissingET *met;
  
  Float_t ener_el[20];
  Float_t px_el[20] ;
  Float_t py_el[20];
  Float_t pz_el[20];
  Float_t pT_el[20];

  Float_t ener_mu[20];
  Float_t px_mu[20] ;
  Float_t py_mu[20];
  Float_t pz_mu[20];
  Float_t pT_mu[20];

  Double_t px_tauhad[50];
  Double_t py_tauhad[50];
  Double_t pz_tauhad[50];
  Double_t pT_tauhad[50];
  Double_t Eta_tauhad[50];
  Double_t Phi_tauhad[50];
  Double_t ener_tauhad[50];
  Double_t  inv_tauhad[50];

  Float_t px_jet[50];
  Float_t py_jet[50];
  Float_t pz_jet[50];
  Float_t pT_jet[50];
  Float_t Eta_jet[50];
  Float_t Phi_jet[50];
  Float_t ener_jet[50];
  Float_t inv_jet[50];
   	
//------------------------------------------------------------------------------
// EDIT CUTS HERE
//------------------------------------------------------------------------------

/* pT cuts */

  Double_t min_pt_el   = 10.0;
  Double_t min_pt_mu   = 10.0;
  Double_t min_pt_jt   = 20.0;
  Double_t min_pt_ph   = 10.0;

/* eta cuts */

  Double_t max_eta_el  = 2.5;
  Double_t max_eta_mu  = 2.5;
  Double_t max_eta_jt  = 4.7;
  Double_t max_eta_ph  = 2.5;


/* Delta_R cuts */

  Double_t min_DR_bb = 0.4;
  Double_t min_DR_elel = 0.4;
  Double_t min_DR_mumu = 0.4;
  Double_t min_DR_bel = 0.4;
  Double_t min_DR_bmu = 0.4;
  Double_t min_DR_mutau = 0.4;

/* Delta_phi cuts */

  Double_t min_DPhi_mutau = 2.7;

/* More cuts... */


    Double_t sc = 121.2;
    

  for(entry = 0; entry <  allEntries; ++entry)// EVENT LOOP BEGINS HERE.
  {
  el_no = 0, mu_no = 0, jt_no = 0, ph_no = 0, b_no = 0, tau_no = 0;

  
    treeReader->ReadEntry(entry);// Load selected branches with data from specified event



// Electron trigger.

   for(i = 0; i < branchElectron->GetEntriesFast(); i++)
   {
     el[i] = (Electron*) branchElectron->At(i);

     if((el[i]->PT > min_pt_el) && (fabs(el[i]->Eta) < max_eta_el))
     {
        el_no++;  
        ener_el[el_no] = (el[i]->P4()).E();
        px_el[el_no] = (el[i]->P4()).Px();
        py_el[el_no] = (el[i]->P4()).Py();
        pz_el[el_no] = (el[i]->P4()).Pz();
        pT_el[el_no] = el[i]->PT;  
     }

   }


// Muon trigger.

   for(i = 0; i < branchMuon->GetEntriesFast(); i++)
   {
     mu[i] = (Muon*) branchMuon->At(i);

     if((mu[i]->PT > min_pt_mu) && (fabs(mu[i]->Eta) < max_eta_mu))
     {
      mu_no++;      
      ener_mu[mu_no] = (mu[i]->P4()).E();
      px_mu[mu_no] = (mu[i]->P4()).Px();
      py_mu[mu_no] = (mu[i]->P4()).Py();
      pz_mu[mu_no] = (mu[i]->P4()).Pz();
      pT_mu[mu_no] = mu[i]->PT;
     }

   }


// Photon trigger.

   for(i = 0; i < branchPhoton->GetEntriesFast(); i++)
   {
     ph[i] = (Photon*) branchPhoton->At(i);

    if((ph[i]->PT > min_pt_ph) && (fabs(ph[i]->Eta) < max_eta_ph))

     {
        ph_no++;    
     }

   }


// Jet trigger.

   for(i = 0; i < branchJet->GetEntriesFast(); i++)
   {
     jt[i] = (Jet*) branchJet->At(i);

     if((jt[i]->PT > min_pt_jt) && (fabs(jt[i]->Eta) < max_eta_jt))
     {
        jt_no++; 
                
        if( jt[i]->BTag == 1)
        {
          b_no++;
        }  
        if((jt[i]->TauTag == 1))   
        {
          tau_no++;
          tau_count = tau_count +1 ;

         /* px_tauhad[tau_no] = (jt[i]->P4()).Px();
          py_tauhad[tau_no] = (jt[i]->P4()).Py();
          pz_tauhad[tau_no] = (jt[i]->P4()).Pz();
          pT_tauhad[tau_no] = jt[i]->PT;
          Eta_tauhad[tau_no] = jt[i]->Eta;
          Phi_tauhad[tau_no] = jt[i]->Phi;
          ener_tauhad[tau_no] = (jt[i]->P4()).E();
          inv_tauhad[tau_no] = (jt[i]->P4()).M();*/

px_jet[jt_no] = (jt[i]->P4()).Px();
py_jet[jt_no] = (jt[i]->P4()).Py();
pz_jet[jt_no] = (jt[i]->P4()).Pz();
pT_jet[jt_no] = jt[i]->PT;
Eta_jet[jt_no] = jt[i]->Eta;
Phi_jet[jt_no] = jt[i]->Phi;
ener_jet[jt_no] = (jt[i]->P4()).E();
inv_jet[jt_no] = (jt[i]->P4()).M();
        } 

     }

   }


// MissingET measurement.

   if(branchMissingET->GetEntriesFast() > 0)
   {
   met = (MissingET*) branchMissingET->At(0);
//   plots->fmissET->Fill(met->MET);
   }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if(((jt_no == 4) && (tau_no == 0) && (mu_no == 0) && (el_no==0) && (b_no == 0) )) // 4j

{
  NEV_C1 =  NEV_C1 + 1;

   pT_l1 = mu[0]->PT;
   xmet = met->MET;
   phi_l1met = DPhi(mu[0]->Phi, (met->P4()).Phi());
   eta_l1 = mu[0]->Eta;
   MTrans = MT((met->P4()).Phi(), mu[0]->Phi , met->MET, pT_l1);
   Meff = pT_l1 + pT_jet[1] + pT_jet[2] + pT_jet[3] + pT_jet[4] + xmet;

   phi_j1met = DPhi(Phi_jet[1], met->Phi);
   phi_j2met = DPhi(Phi_jet[2], met->Phi);
   phi_j1j2 = DPhi(Phi_jet[1],Phi_jet[2]);
   phi_l1j1 = DPhi(mu[0]->Phi,Phi_jet[1]);
   phi_l1j2 = DPhi(mu[0]->Phi,Phi_jet[2]);
   
   eta_j1 = Eta_jet[1];
   eta_j2 = Eta_jet[2];
   eta_j3 = Eta_jet[3];
   eta_j4 = Eta_jet[4];
   DRl1j1 = DR(eta_l1,mu[0]->Phi,eta_j1,Phi_jet[1]);
   DRl1j2 = DR(eta_l1,mu[0]->Phi,eta_j2,Phi_jet[2]);
   DRj1j2 = DR(eta_j1,Phi_jet[1],eta_j2,Phi_jet[2]);
   
   jet1pT = pT_jet[1];
   jet2pT = pT_jet[2];
   jet3pT = pT_jet[3];
   jet4pT = pT_jet[4];
   
  
   Mj1j2sq = pow((ener_jet[1]+ener_jet[2]),2.0) - pow((px_jet[1]+px_jet[2]),2.0) - pow((py_jet[1]+py_jet[2]),2.0) - pow((pz_jet[1]+pz_jet[2]),2.0) ;
   Mj1j2 = sqrt(Mj1j2sq);
 

   plots->fpTmu1->Fill(mu[0]->PT);
   plots->fpTj1->Fill(jt[0]->PT);
   plots->fpTj2->Fill(jt[1]->PT);
   plots->fpTj3->Fill(jt[2]->PT);
   plots->fpTj4->Fill(jt[3]->PT);
   plots->fmissET->Fill(xmet);
   plots->fphil1j1->Fill(DPhi(mu[0]->Phi, jt[0]->Phi)) ;
   plots->fphil1j2->Fill(DPhi(mu[0]->Phi, jt[1]->Phi)) ;
   plots->fetal1->Fill(mu[0]->Eta) ;
   plots->fetaj1->Fill(jt[0]->Eta) ;
   plots->fetaj2->Fill(jt[1]->Eta) ;
   plots->fetaj3->Fill(jt[2]->Eta) ;
   plots->fetaj4->Fill(jt[3]->Eta) ;
   plots->fDRl1j1->Fill(DRl1j1) ;
   plots->fDRl1j2->Fill(DRl1j2) ;
   plots->fMj1j2->Fill(Mj1j2) ;
   plots->fMeff->Fill(Meff) ;
   plots->fMTrans->Fill(MTrans) ;



//if ((jt[0]->PT) > 80 && (jt[1]->PT) > 60 && (jt[2]->PT) > 40 && (jt[3]->PT) > 30)
if((Mj1j2 > 120) || (Mj1j2 < 40))
{
        NEV_C2 =  NEV_C2 + 1;

if( abs(jt[0]->Eta) < 1.0 && abs(jt[1]->Eta) < 1.0 && abs(jt[2]->Eta) < 1.0 && abs(jt[3]->Eta) < 1.0)
    {
       NEV_C3 =  NEV_C3 + 1;

if(met->MET > 30)
{
NEV_C4 = NEV_C4 + 1;


 //if((MTrans < 40) || (MTrans > 90))
   if(Meff > 200)
{
       NEV_C5 =  NEV_C5 + 1;


//if(DRl1j1 < 3.2 && DRl1j2 < 3.2)

  {
   NEV_C6 =  NEV_C6 + 1;

//if( 300 < Meff < 500)
//if(jt[0]->PT > 70)
//if((Mj1j2 < 40) || (Mj1j2 > 100))
{

   NEV_C7 = NEV_C7 +1;

}

}

}
}
        }
}
    
     }



}//EVENT LOOP CLOSES HERE.

cout << allEntries << " " << NEV_C1 << " " << NEV_C2 << " " << NEV_C3 << " "  << NEV_C4 << " " << NEV_C5 << " " << NEV_C6 << " " << NEV_C7 << " " << endl;



 /* plots->fjetmuliplicity->Scale(sc);
  plots->ftransmom->Scale(sc); 
   plots->fIMmutau->Scale(sc); */


}


//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, MyPlots *plots)
{
  result->Print("C");
  result->Print("eps");
}


//------------------------------------------------------------------------------





void test_VLL_4j(const char *sig_2hdm_BP3)

{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");


chain->Add(sig_2hdm_BP3);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  MyPlots *plots = new MyPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

 
  result->Write("results_bkgd_wwz.root");


  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
