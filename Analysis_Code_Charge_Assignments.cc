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

// charge location finder

vector < int > LocationFinder(vector < int > chargearray, int &chflag)     // input is a vector integer and will return a vector int too
{									// 'chflag' is pass-by-reference
       int ncharges = chargearray.size();
       vector < int > PLusloc,MInusloc;
       for (int i=0; i < ncharges; i++){
       if (chargearray[i]==1){
       PLusloc.push_back(i);    					// storing those i in PLusloc  
       cout << "i = " << i << ", ncharges = " << ncharges << endl;      // i = 0,1,2 and ncharges should be always = 3, chargearray.size()=3 for 3 muon event                           
       }  // if loop ends here

       if (chargearray[i]== -1){
       MInusloc.push_back(i);    
                                        // storing those i in MInusloc  
       }  // if loop ends here
 
       } // ncharges loop ends here
						// PLusloc and MImusloc array size are always "2", hence chargelocation[2] doesnot exist.
       if (PLusloc.size() > MInusloc.size()){
       chflag=1;
       return PLusloc;	
       }
       else{
       chflag=-1;
       return MInusloc;
       }
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

//  **********************************Adding Pointers to the branches***************************************************  


  // Get pointers to branches used in this analysis
  TClonesArray *branchMuon= treeReader->UseBranch("Muon");
  TClonesArray *branchElectron= treeReader->UseBranch("Electron");
  TClonesArray *branchMissingET= treeReader->UseBranch("MissingET");

// *********************************************************************************************************************

    TFile *outf = new TFile("mhpp800.root","recreate");   // "recreate" will overwrite if any such file exists

    // Book histograms
    TH1 *histMuonPT= new TH1F("muon pt", "muon P_{T}", 50, 0.0, 1000.0);   // first number is the bin size and then the range


// ************************************** Event Loop Starts here *******************************************************
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry); // storing events from the branch to memory
 
    MissingET *misssingET= (MissingET*) branchMissingET->At(0);
   
//   cout << "Missing ET = " << misssingET->MET << endl;

    int  nMuon = branchMuon->GetEntries();           // Will give no of muons / event
    int  nElectron = branchElectron->GetEntries();  // no of electron / event


    if (nMuon==0 || nElectron==0) continue;	   // if nMuon or nElectron = 0 then skip the event	    
//    cout << "nMuon = " << nMuon << ", nElectron = " << nElectron << " " << endl; 
    

// ******************************************* Muon Things **************************************************************


    vector < TLorentzVector > allMuon;		  // flexible size array is vector, variable size array is allMuon,
							// each event's TLorentzVector will be filled up in this array "allMuon"
    vector < int > muoncharge; 			// muoncharge is a vector, whose each element belongs to the class int						  

    for(Int_t iMuon = 0; iMuon < nMuon; ++iMuon){
   
     Muon *muon= (Muon*) branchMuon->At(iMuon); // branchMuon is a vector, (muon) is one muon belonging to the class Muon
    

     if(muon->PT ==0 || muon->PT > 4000) continue; // filtering out bad muons, if this is satisfied skip the event

     muoncharge.push_back(muon->Charge);	// charges of muon's is now inside the vector, muoncharge vector is getting entries from the muon->Charge

//   cout << muon->Charge << endl;

//   cout << "Muon Index = " << iMuon << ", pT = " << muon->PT << ", Eta = " << muon->Eta << endl;

     TLorentzVector muon_4V = TLorentzVector();        // This is empty construtor of muon_4V vector
     muon_4V.SetPtEtaPhiM(muon->PT,muon->Eta,muon->Phi,0.105);  // giving the 4 vector entries to muon_4V
     allMuon.push_back(muon_4V);                // We are pushing muon 4 vectors inside 'allMuon' array
				//    to push elements into a vector from the back. The new value is inserted into the vector at the end

     } // iMuon loop ends here

// ************************************* Electron Things ******************************************************************


    vector < TLorentzVector > allElectron;		  // flexible size array is vector, variable size array is allElectron,
						  // each event's TLorentzVector will be filled up in this array "allElectron"

    for(Int_t iElectron = 0; iElectron < nElectron; ++iElectron){
   
     Electron *electron= (Electron*) branchElectron->At(iElectron); // branchElectron is a vector, (muon) is one muon belonging to the class Electron
    

     if(electron->PT ==0 || electron->PT > 4000) continue; // filtering out bad muons, if this is satisfied skip the event

     TLorentzVector electron_4V = TLorentzVector();                                                   // This is empty construtor of elctron_4V vector
     electron_4V.SetPtEtaPhiM(electron->PT,electron->Eta,electron->Phi,0.0005);
     allElectron.push_back(electron_4V);                // We are pushing electron 4 vectors inside 'allElectron' array

     } // iElectron loop ends here
 
// ************************************************************************************************************************

    
     if (allMuon.size() < 3) continue;   // skipping if the event has less than 3 Muons
     if (allElectron.size() < 1) continue;   // skipping if the event has less than 1 Electrons

    int chargeflag=0;
    vector < int > chargelocation=LocationFinder(muoncharge,chargeflag); // muoncharge is now takingthe position of chargearray

    cout << "Charge-Flag = " << chargeflag << endl;
    cout << "Same-sign Charge size = " << chargelocation.size() << endl;

     histMuonPT->Fill(allMuon[0].Pt()); // Filling the Pt histogram

 //    cout << "e-mu invariant mass = " << (allElectron[0]+allMuon[0]).M() << endl;

     if ((allMuon[1]+allMuon[2]).M() > 0.106 && (allMuon[1]+allMuon[2]).M() < 4000){ //skipping infinities in the Dimuon mass

//   cout << "DiMuon Mass = " << (allMuon[0]+allMuon[1]).M() << endl;	// each element of allMuon is a four-vector, adding 0th and 1st here
								// .M is the member function of TLorentzvector to give invariant mass
      
      // Print output in the file
   
 //     myfile << (allMuon[chargelocation[0]]+allMuon[chargelocation[1]]).M() << endl; // chargelocation vector is by default the location finder

      myfile << allMuon[chargelocation[1]].Pt() << endl;

   } //skipping infinities in the Dimuon mass

//   cout << allMuon[0].Pt() << endl;   // Here in allMuon Vector we have to write 'Pt' which is defined function under TLorentzVector
  
    }  // Event loop ends here

    myfile.close();

    outf->cd();
    histMuonPT->Write();	// writing the histogram after each event loop
    outf->Close();                  // closing the outfile which is a pointer	    

    cout<<"Exiting Now" << endl;
    return 0;  // end of main
} // int main loop ends here










