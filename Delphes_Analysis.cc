#include <stdio.h>
#include <sstream>  // needed for internal io
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <exception>
#include <cstdlib>
#include <cstdio>
#include <limits>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLorentzVector.h"


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/internal/BasicRandom.hh"
#include <fastjet/tools/MassDropTagger.hh>
#include <fastjet/tools/Filter.hh>


#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

#include "classes/DelphesClasses.h"

using namespace std;


//------------------------------------------------------------------------------
// USER DEFINED FUNCTIONS
//------------------------------------------------------------------------------

// Delta_R 

double DR(double eta1,double phi1,double eta2,double phi2)
{
	double DPhi = TMath::Abs(phi1-phi2);
	double DEta = eta1 - eta2;
	if(DPhi > TMath::Pi())
	DPhi = TMath::TwoPi() - DPhi;
	return TMath::Sqrt(DEta*DEta + DPhi*DPhi);
}

// Delta_Phi 

double DPhi(double phi1, double phi2)
{	
	double delPhi = TMath::Abs(phi1 - phi2);
	if(delPhi > TMath::Pi())
	delPhi = TMath::TwoPi() - delPhi;
	return delPhi;
}

// transverse mass

double MT(double phi1, double phi2 , double ETmiss, double PTrans)
{	
	double delPhi = TMath::Abs(phi1 - phi2);
	return TMath::Sqrt(2*PTrans*ETmiss*(1-cos(delPhi)));
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------


int main(int argc, char *argv[]){
  // Create chain of root trees
    TString in_root, out_dat, n_events;
    int events_to_analyze(10);// Max number of events to be analyzed
    
    if(argc<3){ // argc should be more than 3 for correct execution
	cout<<"usage: "<< argv[0] <<" input_delphes_root output_name [#_events]\n";
		exit(1);
	}
	else{
		in_root = TString ( argv[1]) ;
        out_dat = TString ( argv[2]) ;
        
        if(argc > 3) events_to_analyze= atoi(argv[3]);
	}
    ofstream myfile;
    myfile.open (out_dat);
    
	cout << "will work on maximum \t" << events_to_analyze << "\t events " << endl;
    cout<<"running" << endl;
  gSystem->Load("libDelphes");
  TChain chain("Delphes");
  chain.Add(in_root);
  

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchMuon= treeReader->UseBranch("Muon");


    // Book histograms
   // TH1 *histMuonPT= new TH1F("muon pt", "muon P_{T}", 50, 0.0, 100.0);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    // If event contains at least 1 lepton
    if(branchMuon->GetEntries() > 0){

      // Take first Muon
      Muon *muon= (Muon*) branchMuon->At(0);
      
      // Print electron transverse momentum
      myfile << muon->PT << endl;
    }
  }

    myfile.close();
    cout<<"Exiting Now" << endl;
    return 0;  // end of main
}

