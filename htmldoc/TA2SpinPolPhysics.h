//--Author	PP Martel    10th Oct 2010
//--Update	PP Martel    19th Feb 2011
//
// TA2SpinPolPhysics
//
// Reconstruction of Compton and Pi0 kinematics. 

#ifndef __TA2SpinPolPhysics_h__
#define __TA2SpinPolPhysics_h__

#include "TAcquRoot.h"
#include "TAcquFile.h"
#include "TA2Physics.h"
#include "TA2Analysis.h"
#include "TA2Tagger.h"
#include "TA2CB.h"
#include "TA2Taps.h"
#include "TA2Ladder.h"
#include "TA2PhotoPhysics.h"
#include "TA2Event.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>

class TA2Apparatus;

class TA2SpinPolPhysics : public TA2Physics {
  
 protected:

  //Root File and Tree Outputs
  Bool_t fEventSave;
  TFile* fEventFile;
  Char_t fEventFileName[100];
  TTree* fEventTree;
  TA2Event fEvent;

  Bool_t fCompWSave;             // Comp With Proton
  TFile* fCompWFile;
  Char_t fCompWFileName[100];
  TTree* fCompWTree;

  Bool_t fCompMSave;             // Comp Missing Proton
  TFile* fCompMFile;
  Char_t fCompMFileName[100];
  TTree* fCompMTree;

  Bool_t fPionWSave;             // Pion With Proton
  TFile* fPionWFile;
  Char_t fPionWFileName[100];
  TTree* fPionWTree;

  Bool_t fPionMSave;             // Pion Missing Proton
  TFile* fPionMFile;
  Char_t fPionMFileName[100];
  TTree* fPionMTree;
  
  TA2Tagger* fTAGG;              // Glasgow photon tagger
  TA2CB* fCB;                    // Crystal Ball
  TA2Taps* fTAPS;                // TAPS
  TA2Ladder* fLADD;              // Ladder
  
  TA2Particle* fTAGGpart;        // TA2Particles from Tagger
  TA2Particle* fCBpart;          // TA2Particles from CB
  TA2Particle* fTAPSpart;        // TA2Particles from TAPS
  
  TA2Particle** fPARTtagged;     // TA2Particle tagger photon
  TA2Particle** fPARTphoton;     // TA2Particle photon
  TA2Particle** fPARTproton;     // TA2Particle proton
  TA2Particle** fPARTpiplus;     // TA2Particle piplus
  TA2Particle** fPARTneutron;    // TA2Particle neutron
  TA2Particle** fPARTrootino;    // TA2Particle rootino
  
  TA2Particle** fPARTpi0;	 // TA2Particle pi0
  TA2Particle** fPARTeta;        // TA2Particle eta
  TA2Particle** fPARTgprime;     // TA2Particle gprime
  
  TLorentzVector fP4photonTot;   // total 4-momentum of gammas

  Int_t fTgRefTDC;               // Tagger TDC for synch check
  Int_t fCBRefTDC;               // CB TDC for synch check
  Int_t fSynchDif;               // Synch check

  Int_t fCherADC;                // Cherenkov ADC
  Double_t fCherTDC0;            // Cherenkov TDC M0
  Double_t fCherTDC1;            // Cherenkov TDC M1
  Double_t fCherTDC2;            // Cherenkov TDC M2

  Int_t fBeamPol;                // Beam Helicity

  Int_t fNphoton;                // # photon
  Int_t fNproton;                // # proton
  Int_t fNpiplus;                // # pi+
  Int_t fNneutron;               // # neutron
  Int_t fNrootino;               // # unknowns
  Int_t fNgprime;                // # gprimes
  Int_t fNprompt;                // tagger prompts
  Int_t fNrandom;                // tagger randoms
  Int_t fMaxTagg;                // max # tagger hits
  
  Int_t fNpi0;                   // # pi0
  Int_t fNeta;                   // # eta
  
  Double_t* fMassDpi0;           // for meson ID by inv. mass
  Double_t* fMassDeta;           // for meson ID by inv. mass
  Int_t* fMassIJ;                // combinatorial indices
  Int_t* fMassIpi0;              // ditto
  Int_t* fMassIeta;              // ditto
  Bool_t* fIsMesonIndex;         // photon derived from a meson?
  Double_t fMaxMDpi0;            // mass-diff limit pi0
  Double_t fMaxMDeta;            // mass-diff limit eta

  Double_t fPromptL;             // Prompt-Random windows
  Double_t fPromptH;
  Double_t fRand1Lo;
  Double_t fRand1Hi;
  Double_t fRand2Lo;
  Double_t fRand2Hi;
  
  Double_t fM2g;                 // 2-photon invariant mass
  Double_t fM6g;                 // 6-photon invariant mass

  Bool_t fProtCB;                // Proton in CB?
  Bool_t fPrTAPS;                // Proton in TAPS?
  Bool_t fIsComp;
  Bool_t fIsPion;
  
  Double_t fPhotEk;              // Photon variables
  Double_t fPhotPx;
  Double_t fPhotPy;
  Double_t fPhotPz;
  Double_t fPhotTh;
  Double_t fPhotPh;
  Double_t fPhotTm;
  
  Double_t fPionEk;              // Pion variables
  Double_t fPionPx;
  Double_t fPionPy;
  Double_t fPionPz;
  Double_t fPionMa;
  Double_t fPionTh;
  Double_t fPionPh;
  Double_t fPionTm;

  Double_t fDec1Ek;              // Decay 1 variables
  Double_t fDec1Px;
  Double_t fDec1Py;
  Double_t fDec1Pz;
  Double_t fDec1Th;
  Double_t fDec1Ph;
  Double_t fDec1Tm;

  Double_t fDec2Ek;              // Decay 2 variables
  Double_t fDec2Px;
  Double_t fDec2Py;
  Double_t fDec2Pz;
  Double_t fDec2Th;
  Double_t fDec2Ph;
  Double_t fDec2Tm;

  Double_t fPionOA;
 
  Double_t fProtUn;              // Proton variables
  Double_t fProtEk;
  Double_t fProtPx;
  Double_t fProtPy;
  Double_t fProtPz;
  Double_t fProtTh;
  Double_t fProtPh;
  Double_t fProtTm;

  Int_t fCBCrystN;
  Int_t fTAPSCstN;

  Double_t fMPrToPh;             // Mass for incorrectly IDed protons
  
  Int_t* fTaggCh;                // Tagged variables
  Int_t* fTaggChP;
  Int_t* fTaggChR;
  Double_t* fTaggEk;
  Double_t* fTaggEkP;
  Double_t* fTaggEkR;
  Double_t* fTaggTm;
  Double_t* fTaggTmP;
  Double_t* fTaggTmR;

  Double_t* fDifTime;

  Double_t* fRecoEk;
  Double_t* fRecoEkP;
  Double_t* fRecoEkR;
  Double_t fRecoPx;
  Double_t fRecoPy;
  Double_t* fRecoPz;
  Double_t* fRecoPzP;
  Double_t* fRecoPzR;
  Double_t* fRecoTh;
  Double_t* fRecoThP;
  Double_t* fRecoThR;
  Double_t* fRecoPh;
  Double_t* fRecoPhP;
  Double_t* fRecoPhR;
  Double_t* fRecoMa;
  Double_t* fRecoMaP;
  Double_t* fRecoMaR;

  Double_t* fMissEk;
  Double_t* fMissEkP;
  Double_t* fMissEkR;
  Double_t fMissPx;
  Double_t fMissPy;
  Double_t* fMissPz;
  Double_t* fMissPzP;
  Double_t* fMissPzR;
  Double_t fMissPt;
  Double_t* fMissPr;
  Double_t* fMissPrP;
  Double_t* fMissPrR;

  Double_t* fProtOA;
  Double_t* fProtOAP;
  Double_t* fProtOAR;

  //Variables required for energy loss correction

  TF1* fECorrCB;
  TF1* fECoTAPS;

  TF1** fELossCB;
  TF1** fELoTAPS;

  Bool_t fTableCB;
  Bool_t fTblTAPS;
  
 public:
  
  TA2SpinPolPhysics(const char*, TA2Analysis*);
  virtual ~TA2SpinPolPhysics();
  virtual void LoadVariable();   // variables for display/cuts
  virtual void PostInit( );      // init after parameter input
  virtual void SetConfig(Char_t*, Int_t);
  
  virtual void Reconstruct();    // reconstruct detector info
  virtual TA2DataManager* CreateChild(const char*, Int_t){return NULL;}
  virtual void CloseFile();

  virtual Bool_t MakeTable( TString );
  
  virtual void Sort2Photon( );
  virtual void SortNPhoton( );
  virtual void MarkUndefined(Int_t);
  virtual void MarkEndBuffer( );
  
  Double_t Sqr(Double_t x){return (x*x);}
  
  ClassDef(TA2SpinPolPhysics,1)
    };

//-----------------------------------------------------------------------------

inline Bool_t TA2SpinPolPhysics::MakeTable(TString det)
{
    
  TString filename = "ELookup"+det+".table";
  ifstream file (filename);
  
  if(!file){
    printf ("Lookup table not found\n");
    return false;
  }
  
  Int_t nX, nE;
  Double_t valX, valE, valL;
  Double_t* fitX;
  Double_t* fitY;    
  char buffer[256] ;
  
  TFile *myfile = new TFile((det+"Fits.root"),"RECREATE");
  
  file >> nX >> nE;
  
  if(det=="CB"){
    /*
    fECorrCB = new TF1("fECorrCB","[0]+([1]*x)+([2]*x*x)+([3]*x*x*x)+([4]*x*x*x*x)+([5]*x*x*x*x*x)",0,500);
    fECorrCB->SetParameters(1.11,-0.0384,0.000458,-0.00000263,0.00000000724,-0.00000000000772);
    */
    fECorrCB = new TF1("fECorrCB","expo+[2]",0,500);
    fECorrCB->SetParameters(0.221,-0.0362,-0.130);
    printf ("Creating CB lookup table\n");
    fELossCB = new TF1*[nX];
  }
  else if(det=="TAPS"){
    /*
    fECoTAPS = new TF1("fECoTAPS","[0]+([1]*x)+([2]*x*x)+([3]*x*x*x)+([4]*x*x*x*x)+([5]*x*x*x*x*x)",0,500);
    fECoTAPS->SetParameters(1.11,-0.0384,0.000458,-0.00000263,0.00000000724,-0.00000000000772);
    */
    fECoTAPS = new TF1("fECoTAPS","expo+[2]",0,500);
    fECoTAPS->SetParameters(0.221,-0.0362,-0.130);
    printf ("Creating TAPS lookup table\n");
    fELoTAPS = new TF1*[nX];
  }
  else{
    printf ("Invalid detector for lookup table\n");
    return false;
  }
  
  fitX = new Double_t[nE];
  fitY = new Double_t[nE];
  
  //Ignore second line
  file.getline(buffer,256);
  file.getline(buffer,256);
  
  //Read in the data
  for(Int_t iX=0; iX<nX; iX++){
    
    for(Int_t iE=0; iE<nE; iE++){   //Read in the data   
      
      file >> valX >> valE >> valL;
      
      fitX[iE] = valE;
      fitY[iE] = valL;
      
      //printf("%d\t%f\t%d\t%f\t%f\n", iX, valX, iE, valE, valL);
      
    }
    
    TGraph *loss = new TGraph(nE, fitX, fitY);
    
    TF1 *fit1 = new TF1("fit1","expo+[2]",0,500);
    
    loss->Fit("fit1","Q","",20,180);
    
    if(det=="CB"){
      fELossCB[iX] = fit1;
    }
    else if(det=="TAPS"){
      fELoTAPS[iX] = fit1;
    }
    
    myfile->cd();
    loss->Write();
    fit1->Write();
    
  }
  
  myfile->Close();
  delete myfile;
  
  file.close();
  
  return true;
      
}

//-----------------------------------------------------------------------------

inline void TA2SpinPolPhysics::Sort2Photon()
{
  // Test if 2 gamma 4-momenta combine to give a  possible pi0 or eta
  // by finding the invariant mass. This is a fast version of SortNPhoton
  // for the frequent pi0 or eta -> 2 gamma situation
  
  TA2Particle phot1 = *fPARTphoton[0];
  TA2Particle phot2 = *fPARTphoton[1];
  
  Double_t time = (phot1.GetTime() + phot2.GetTime())/2;
  
  TLorentzVector p4 = phot1.GetP4() + phot2.GetP4();   // sum 4 momenta
  
  fM2g = p4.M();				       // inv. mass
  
  fMassDpi0[0] = TMath::Abs( p4.M() - fParticleID->GetMassMeV( kPi0));
  fMassDeta[0] = TMath::Abs( p4.M() - fParticleID->GetMassMeV( kEta));
  
  // is it a pi0
  if ( fMassDpi0[0] < fMaxMDpi0 ) {
  //printf("%6.1f  %6.1f  %6.1f  %6.1f\n", fMaxMDpi0, fMassDpi0[0], fM2g, time);
    (*fPARTpi0[0]).SetP4( p4);
    (*fPARTpi0[0]).SetTime( time);
    fNpi0 = 1;
  }
  // wasn't a pi0 so test for eta
  else if ( fMassDeta[0] < fMaxMDeta ) {
    (*fPARTeta[0]).SetP4( p4);
    (*fPARTeta[0]).SetTime( time);
    fNeta = 1;
  }
  // not an eta so assume both photons are gamma-prime
  else {
    fPARTgprime[0] = fPARTphoton[0];
    fPARTgprime[1] = fPARTphoton[1];
    fNgprime = 2;
  }
  return;
}

//-----------------------------------------------------------------------------

inline void TA2SpinPolPhysics::SortNPhoton()
{
  // Take sample of n gamma 4-momenta and sort into possible pi0 or eta
  // by finding the invariant mass of the possible 2-photon combinations.
  // The differences of the 2-photon invariant masses and pi0/eta mass are
  // sorted in ascending order.
  
  TLorentzVector p4temp, q4;
  Double_t mPi0 = fParticleID->GetMassMeV( kPi0);
  Double_t mEta = fParticleID->GetMassMeV( kEta);
  Double_t* pi0diff = fMassDpi0;
  Double_t* etadiff = fMassDeta;
  Int_t* ij = fMassIJ;
  Int_t i,j,k,jk;
  Int_t n = fNphoton;
  fP4photonTot.SetXYZT(0.0,0.0,0.0,0.0);         // zero total photon 4-mom
  
  Double_t time;
  
  // Loop over possible 2-photon combinations i != j, ij = ji
  for( i = 0, k = 0; i < n; i++ ){
    
    fP4photonTot += (*fPARTphoton[i]).GetP4();
    fIsMesonIndex[i] = EFalse;                   // initialise not meson
    
    for( j = i + 1; j < n; j++ ){
      
      // add the 4 momenta
      p4temp =  (*fPARTphoton[i]).GetP4() + (*fPARTphoton[j]).GetP4();
      *pi0diff++ = TMath::Abs( p4temp.M() - mPi0 ); // check inv mass pi0
      *etadiff++ = TMath::Abs( p4temp.M() - mEta ); // check inv. mass eta
      *ij++ = i | (j<<16);                       // store i,j indices
      k++;                                       // permutation counter
      
    }
  }
  
  if ( n == 6 ) fM6g = fP4photonTot.M();
  TMath::Sort(k, fMassDeta, fMassIeta, EFalse);  // sort mass diffs ascending
  TMath::Sort(k, fMassDpi0, fMassIpi0, EFalse);  // sort mass diffs ascending
  Int_t nMeson = k;                              // max # possible mesons
  
  // Check for eta combinations
  for ( i = 0; i < nMeson; i++) {
    
    // exit loop when mass-diff exceeds maximum bound
    if( fMassDeta[fMassIeta[i]] > fMaxMDeta ) break;
    
    // get indices
    jk = fMassIJ[fMassIeta[i]];
    j = jk & 0xffff;
    k = (jk>>16) & 0xffff;
    
    // photon already used to construct meson ?
    if( fIsMesonIndex[j] || fIsMesonIndex[k] ) continue;
    
    // set 4-momentum
    q4 = (*fPARTphoton[j]).GetP4() + (*fPARTphoton[k]).GetP4();
    (*fPARTeta[fNeta]).SetP4( q4);
    
    // set time
    time = ((*fPARTphoton[j]).GetTime() + (*fPARTphoton[k]).GetTime())/2;
    (*fPARTeta[fNeta++]).SetTime( time);
    
    // mark photons as used
    fIsMesonIndex[j] = fIsMesonIndex[k] = ETrue;
  }
  
  // Check for pi0 combinations
  for ( i = 0; i < nMeson; i++) {
    
    // exit loop when mass-diff exceeds maximum bound
    if ( fMassDpi0[fMassIpi0[i]] > fMaxMDpi0 ) break;
    
    // get indices
    jk = fMassIJ[fMassIpi0[i]];
    j = jk & 0xffff;
    k = (jk>>16) & 0xffff;
    
    // photon already used to construct meson?
    if ( fIsMesonIndex[j] || fIsMesonIndex[k] ) continue;
    
    // set 4-momentum
    q4 = (*fPARTphoton[j]).GetP4() + (*fPARTphoton[k]).GetP4();
    (*fPARTpi0[fNpi0]).SetP4( q4);
    
    // set time
    time = ((*fPARTphoton[j]).GetTime() + (*fPARTphoton[k]).GetTime())/2;
    (*fPARTpi0[fNpi0++]).SetTime( time);
    
    // mark photons as used
    fIsMesonIndex[j] = fIsMesonIndex[k] = ETrue;
  }
  
  // Put any photons not combined into a meson into the gamma-primed list
  for ( j = 0; j < n; j++ ) {
    if ( !fIsMesonIndex[j] ) fPARTgprime[fNgprime++] = fPARTphoton[j];
  }
  
}

//-----------------------------------------------------------------------------

inline void TA2SpinPolPhysics::MarkUndefined( Int_t jtagg )
{
  // Initialise undefined
  fTaggChP[jtagg] = EPhotoUndefined;
  fTaggEkP[jtagg] = EPhotoUndefined;
  fTaggTmP[jtagg] = EPhotoUndefined;

  fRecoEkP[jtagg] = EPhotoUndefined;
  fRecoPzP[jtagg] = EPhotoUndefined;
  fRecoThP[jtagg] = EPhotoUndefined;
  fRecoPhP[jtagg] = EPhotoUndefined;
  fRecoMaP[jtagg] = EPhotoUndefined;

  fMissEkP[jtagg] = EPhotoUndefined;
  fMissPzP[jtagg] = EPhotoUndefined;
  fMissPrP[jtagg] = EPhotoUndefined;

  fProtOAP[jtagg] = EPhotoUndefined;
}

//-----------------------------------------------------------------------------

inline void TA2SpinPolPhysics::MarkEndBuffer()
{
  // Ensure the multi-data buffers are marked as ended
  fTaggChP[fNprompt] = EBufferEnd;
  fTaggChR[fNrandom] = EBufferEnd;
  fTaggEkP[fNprompt] = EBufferEnd;
  fTaggEkR[fNrandom] = EBufferEnd;
  fTaggTmP[fNprompt] = EBufferEnd;
  fTaggTmR[fNrandom] = EBufferEnd;

  fRecoEkP[fNprompt] = EBufferEnd;
  fRecoEkR[fNrandom] = EBufferEnd;
  fRecoPzP[fNprompt] = EBufferEnd;
  fRecoPzR[fNrandom] = EBufferEnd;
  fRecoThP[fNprompt] = EBufferEnd;
  fRecoThR[fNrandom] = EBufferEnd;
  fRecoPhP[fNprompt] = EBufferEnd;
  fRecoPhR[fNrandom] = EBufferEnd;
  fRecoMaP[fNprompt] = EBufferEnd;
  fRecoMaR[fNrandom] = EBufferEnd;
  
  fMissEkP[fNprompt*fNproton] = EBufferEnd;
  fMissEkR[fNrandom*fNproton] = EBufferEnd;
  fMissPzP[fNprompt*fNproton] = EBufferEnd;
  fMissPzR[fNrandom*fNproton] = EBufferEnd;
  fMissPrP[fNprompt*fNproton] = EBufferEnd;
  fMissPrR[fNrandom*fNproton] = EBufferEnd;

  fProtOAP[fNprompt*fNproton] = EBufferEnd;
  fProtOAR[fNrandom*fNproton] = EBufferEnd;
}

#endif
