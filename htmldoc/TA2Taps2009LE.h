////////////////////////////////////////////////////////////////////////////
//                                                                        //
//     --Author		H. Berghaeuser June 2008                          //
//                      HenningBerghaeuser@web.de                         //
//     --Update         H.Berghaeuser Okt 2008                            //
//     --Update         H.Berghaeuser Dez 2008 -> TAPS2009LE Version      //
//                                                                        //
//     -- Use the TAPS2009LE.dat !                                        //
//                                                                        //
//      TA2Taps2009LE is a new TAPS Class for Acqu 4v2+                   //
//      Read BaF2, Veto and  use Veto-dE                                  //
//                                                                        //
//      LE stands for light edition                                       //
//      this class is based on JRM Annands TA2Calorimeter class           //
//                                                                        //
// -----------------------------------------------------------------------//
//                                                                        //
// compatible to Acqu 4v2 and 4v3  (4v1 and 4v0 should also work)         //
// use the TA2Taps2009LE in combination with J.Annands standard           //
// AcquUser-Code ()
//                                                                        //
//            TA2Taps2009LE.h            Version 0.95 (alpha)             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef __TA2Taps2009LE_h__
#define __TA2Taps2009LE_h__

#include "TA2TAPS_BaF2.h"
#include "TA2PlasticPID.h"
#include "TA2Calorimeter.h"
#include "TCutG.h"           // important for dEvE BaF2-Veto-Cuts
#include "TH1F.h" 
#include "TH2F.h" 
#include "TNtuple.h"


class TA2Taps2009LE : public TA2Apparatus {
protected:

  //concerning class-management
  UInt_t fTapsRunStep;                  // step of running
  UInt_t fCalibVetoStep;                // step of Veto  calibration
  UInt_t fCalibBaF2Step;                // step of BaF2 of calibration
  //Detectors
  TA2TAPS_BaF2* fBaF2;                  // BaF2-calorimeeter
  TA2Detector* fVeto;                   // Veto-Detectors
  // Concerning BaF2
  Int_t fBaF2_Ncluster;                 // Nr of BaF2 clusters in event
  UInt_t* id_clBaF2;                    // indices of hit clusters
  HitCluster_t** clBaF2;                // -> cluster structs
  HitCluster_t* cl;                     // cluster struct
  UInt_t clBaF2_Nhits;                  // crystal-hits inside a cluster
  UInt_t* clBaF2_elements;              // crystal-elements inside a cluster
  // Concerning Vetos
  Bool_t* fIsVCharged;                  // to save if hit is charged or not
  UInt_t fVeto_NHits;                   // # veto-hits in this event
  Int_t* fVeto_Hits;                    // all veto-elements that fired
  UInt_t fThisVetoFired;                // element-nr of veto that fired
  Double_t fVeto_dE;                    // to save Veto-dE
  UInt_t fMultipleVetoHit;              // will be used in v1
  // Variables concerning reconstruction
  TLorentzVector* fP4_Nphoton;          // Total invariant mass N photons
  Double_t* fM_Nphoton;                 // Inv mass for N cluster totals
  TVector3 fp3;                         // used to make p4
  // reconstruction: make ID of particle:
  Bool_t* fdEvE_IsProton;               // array important for CheckParticleID
  Bool_t* fdEvE_IsChPion;               // array important for CheckParticleID
  Bool_t* fdEvE_IsElectron;             // array important for CheckParticleID
  // Variables concerning SetConfig .... use TAPS2009LE.dat
  Bool_t fSimpleReconstruct;
  Double_t fTapsFudge;
  Double_t fTimeShift;
  Double_t* fDeltaE;
  Double_t* fEcharged;

  TFile* fTAPS_dEvE_CutFile;
  char   fTAPS_dEvE_Cuts[250];
  Bool_t fTAPS_dEvE_Proton;
  TCutG* fTAPS_dEvE_ProtonCut;
  char   fTAPS_dEvE_Proton_CutName[50];
  Bool_t fTAPS_dEvE_ChPion;
  TCutG* fTAPS_dEvE_ChPionCut;
  char   fTAPS_dEvE_ChPion_CutName[50];
  Bool_t fTAPS_dEvE_Electron;
  TCutG* fTAPS_dEvE_ElectronCut;
  char   fTAPS_dEvE_Electron_CutName[50];

  // Calibration outout;
  TH1F* calibTAPS_m1g_Single;
  TH2F* calibTAPS_m1g_AllCh;
  TH2F* calibTAPS_VetoCorr;
  TNtuple* VetoData;
  TH2F* calib_TAPSdEvEcl;
  TH2F* calib_TAPSdEvE;


public:

  TA2Taps2009LE( const char*, TA2System* ); 
  virtual ~TA2Taps2009LE();
  virtual void PostInit();
  virtual void SetConfig( char* line, int key );
  virtual TA2DataManager* CreateChild( const char*, Int_t );
  virtual void MakeAllRootinos();
  virtual void LoadVariable();
  virtual void SimplePhotonReconstruct();
  virtual void Reconstruct();
  virtual void MainReconstruct();
  virtual TLorentzVector SetLorentzVector( Double_t, HitCluster_t* );
  virtual void Cleanup();

  TA2ClusterDetector* GetCal(){ return fBaF2; }
  TA2Detector* GetVeto(){ return fVeto; }
  TLorentzVector* GetP4_Nphoton(){ return fP4_Nphoton; }
  Double_t* GetM_Nphoton(){ return fM_Nphoton; }
  Int_t CheckParticleID(UInt_t i);


  // ----- Calibration Functions ----------------------------
  virtual void CalibrateVetoEnergy();
  virtual void CalibrateBaF2Energy(TLorentzVector *photonCB);
  Float_t TAPS_GetCLInfo(UInt_t module, UInt_t TimeOrEnergy);


  ClassDef(TA2Taps2009LE,1)
};





inline void TA2Taps2009LE::Reconstruct(){

 if(fSimpleReconstruct == kTRUE) 
 {SimplePhotonReconstruct();} 
  else {MainReconstruct();}
 
 fTapsRunStep = 2; 
 } 



inline void TA2Taps2009LE::Cleanup( )
{
  // Clear any arrays holding variables
 TA2DataManager::Cleanup();

 if( !fNparticle ) return;
 Int_t i;
 for( i=0; i<fNparticle; i++ ){
   fP4[i].SetXYZT(0.,0.,0.,0.);
   fMinv[i] = 0.0;
 }
 fP4_Nphoton[i-1].SetXYZT(0.,0.,0.,0.);
 fM_Nphoton[i-1] = 0.0;
 fNparticle = 0;
}

//-----------------------------------------------------------------------------
inline void TA2Taps2009LE::SimplePhotonReconstruct( )
{
 // All in TAPS detected particles are reconstructed as PHOTONS
 // does absolutely the same as TA2Calorimeter does!

 // now VETO is used .... no dEvE is used ! JUST BaF2 readout

 if(fTapsRunStep == 1)
 {
 printf("\n\nINFO: TAPS only reads BaF2s and 'makes' only photons !!!\n\n"); 
 }

  fNparticle = fBaF2->GetNCluster(); 
  if( !fNparticle ) return;
  id_clBaF2 = fBaF2->GetClustHit();
  clBaF2 = fBaF2->GetCluster();


  fP4tot.SetXYZT(0,0,0,0);
  Int_t i;
  for (i=0; i<fNparticle; i++){
    cl = clBaF2[id_clBaF2[i]];
    fPDG_ID[i] = kGamma;
    fParticleID->SetP4( fP4+i,fPDG_ID[i],cl->GetEnergy(),cl->GetMeanPosition() );
    fP4tot += fP4[i]; 
    fMinv[i] = fP4[i].M();

  } 
  fP4_Nphoton[i - 1]= fP4tot;
  fM_Nphoton[i - 1]= fP4tot.M();
  fMtot = fP4tot.M();

}


//-----------------------------------------------------------------------------
inline TLorentzVector TA2Taps2009LE::SetLorentzVector(Double_t mass, HitCluster_t* cl)
{
  // Get cluster kinetic energy and particle direction
  // Construct 4 momentum using input assumed mass

  Double_t E = cl->GetEnergy() + mass;
  Double_t P = TMath::Sqrt(E*E - mass*mass);
  TLorentzVector p4;
  p4.SetE(E);
  p4.SetVect( (cl->GetMeanPosition())->Unit() * P );
  return p4;
}


#endif

