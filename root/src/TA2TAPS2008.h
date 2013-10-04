////////////////////////////////////////////////////////////////////////////
//                                                                        //
//     --Author		H. Berghaeuser June 2008                          //
//                      HenningBerghaeuser@web.de                         //
//     --Update         H.Berghaeuser Okt 2008                            //
//                                                                        //
//     -- Use the TAPS2008v0.XY.dat !                                     //
//                                                                        //
//      TA2TAPS2008 is a new TAPS Class for Acqu 4vXX                     //
//      Read BaF2, Veto and PbWO4  &  use TOF, PSA, Veto-dE               //
//                                                                        //
//            TA2TAPS2008.h            Version 0.1                        //
//                                                                        //
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    //
// !! IMPORTANT:                                                    !!    //
// !! you need to use TA2CB2008 as class for the CrystalBall        !!    //
// !! if you don't want to use  TA2CB2008 scroll down to the bottom !!    //
// !! of this file an read !                                        !!    //
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    //
//                                                                        //
//                                                                        //
//     this class is based on JRM Annands TA2Calorimeter class            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef __TA2TAPS2008_h__
#define __TA2TAPS2008_h__

#include "TA2TAPS_BaF2.h"
#include "TA2PlasticPID.h"
#include "TA2TAPS_PbWO4.h"
#include "TA2CrystalBall.h"
#include "TA2Calorimeter.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TCutG.h"
#include "TA2Particle.h"

class TA2TAPS2008 : public TA2Apparatus
{
 private:
 protected:
  //concerning output
  TFile   *fTAPS_OUTPUT_FILE;           // OUTPUTFILE set in TAPS2008vX.dat
  // Char_t  fTAPS_OUTPUT_FILE_Name[256];  // OUTPUTFILE set in TAPS2008vX.dat
  // Bool_t  fTAPS_OUTPUT_FILE_Exists;     // if false default-file will be created
  Bool_t  fPrintSetUpInfo;              // if kTRUE ... info/values from dat file are printed
  //  TNtuple *fTAPS_TOF;                   // NTUPLE containing TOF-information
  //  TNtuple *fTAPS_dEvE;                  // NTUPLE containing dEvE-information
  //  TNtuple *fTAPS_PSA;                   // NTUPLE containing PSA-information
  //concerning class-management
  UInt_t fTapsRunStep;                  // step of running
  //Detectors
  TA2TAPS_BaF2* fBaF2;                  // BaF2-calorimeeter
  TA2Detector* fVeto;                   // Veto-Detectors
  TA2Detector* fPbWO4;                  // PBwO4-Detectors
  TA2CrystalBall* fCB;                       // Crystal Ball Apparatus
  // Concerning BaF2
  Int_t fBaF2_Ncluster;                 // Nr of BaF2 clusters in event
  UInt_t* id_clBaF2;                    // indices of hit clusters
  HitCluster_t** clBaF2;                // -> cluster structs
  HitCluster_t* cl;                     // cluster struct
  UInt_t clBaF2_Nhits;                  // crystal-hits inside a cluster
  UInt_t* clBaF2_elements;              // crystal-elements inside a cluster
  Double_t* BaF2Time;                   // pure BaF2-cl-time
  Double_t  fTOF;                       // to save TOF.time
  UInt_t  fNoCBTimeCounter;             // count up if CB doesnt provide TIME
  Double_t  fShortGateValue;            // to save short-gate for PSA
  Double_t  fLongGateValue;             // to save long-gate for PSA
  // Concerning Vetos
  Bool_t* fIsVCharged;                  // to save if hit is charged or not
  UInt_t fVeto_NHits;                   // # veto-hits in this event
  Int_t* fVeto_Hits;                    // all veto-elements that fired
  UInt_t fThisVetoFired;                // element-nr of veto that fired
  Int_t fVeto_index;                    // to save Veto-element-number
  Double_t fVeto_dE;                    // to save Veto-dE
  Double_t fVeto_T;                     // to save Veto-time
  UInt_t fMultipleVetoHit;              // will be used in v1
  // Variables concerning reconstruction
  // TLorentzVector* fP4_Nphoton;          // Total invariant mass N photons
  // Double_t* fM_Nphoton;                 // Inv mass for N cluster totals
  TVector3 fp3;                         // used to make p4
  TA2Particle* particles;               // use this to store particle-information
  UInt_t  NParticles;                   // # of reconstructed particles
  Int_t  *fPDG_ID_sec;                  // to store particles second (maybe) id
  Bool_t *fPDG_ID_unclear;              // if fPDG_ID_sec is set then this should be true
  // reconstruction: make ID of particle:
  Bool_t* fdEvE_IsProton;               // array important for CheckParticleID
  Bool_t* fdEvE_IsChPion;               // array important for CheckParticleID
  Bool_t* fdEvE_IsElectron;             // array important for CheckParticleID
  Bool_t* fTOF_IsNeutron;               // array important for CheckParticleID
  Bool_t* fTOF_IsProton;                // array important for CheckParticleID
  Bool_t* fTOF_IsChPion;                // array important for CheckParticleID
  Bool_t* fTOF_IsElectron;              // array important for CheckParticleID
  Bool_t* fTOF_IsGamma;                 // array important for CheckParticleID
  Bool_t* fPSA_IsNucleon;               // array important for CheckParticleID
  Bool_t* fPSA_IsNoNucleon;             // array important for CheckParticleID
  // Variables concerning SetConfig .... use TAPS2008vX.dat
  Bool_t fSimpleReconstruct;
  Double_t fTapsFudge;
  Double_t fTimeShift;
  TFile* fTAPS_TOF_CutFile;
  TFile* fTAPS_dEvE_CutFile;
  TFile* fTAPS_PSA_CutFile;
  Bool_t fTAPS_TOF_USEAGE;
  Char_t fTAPS_TOF_Cuts[256];
  Bool_t fTAPS_TOF_Neutron;
  TCutG* fTAPS_TOF_NeutronCut;
  Char_t fTAPS_TOF_Neutron_CutName[32];
  Bool_t fTAPS_TOF_Proton;
  TCutG* fTAPS_TOF_ProtonCut;
  Char_t fTAPS_TOF_Proton_CutName[32];
  Bool_t fTAPS_TOF_ChPion;
  TCutG* fTAPS_TOF_ChPionCut;
  Char_t fTAPS_TOF_ChPion_CutName[32];
  Bool_t fTAPS_TOF_Electron;
  TCutG* fTAPS_TOF_ElectronCut;
  Char_t fTAPS_TOF_Electron_CutName[32];
  Bool_t fTAPS_TOF_Gamma;
  TCutG* fTAPS_TOF_GammaCut;
  Char_t fTAPS_TOF_Gamma_CutName[32];
  Bool_t fTAPS_dEvE_USEAGE;
  Char_t fTAPS_dEvE_Cuts[256];
  Bool_t fTAPS_dEvE_Proton;
  TCutG* fTAPS_dEvE_ProtonCut;
  Char_t fTAPS_dEvE_Proton_CutName[32];
  Bool_t fTAPS_dEvE_ChPion;
  TCutG* fTAPS_dEvE_ChPionCut;
  Char_t fTAPS_dEvE_ChPion_CutName[32];
  Bool_t fTAPS_dEvE_Electron;
  TCutG* fTAPS_dEvE_ElectronCut;
  Char_t fTAPS_dEvE_Electron_CutName[32];
  Bool_t fTAPS_PSA_USEAGE;
  Char_t fTAPS_PSA_Cuts[256];
  Bool_t fTAPS_PSA_Nucleon;
  TCutG* fTAPS_PSA_NucleonCut;
  Char_t fTAPS_PSA_Nucleon_CutName[32];
  Bool_t fTAPS_PSA_noNucleon;
  TCutG* fTAPS_PSA_noNucleonCut;
  Char_t fTAPS_PSA_noNucleon_CutName[32];
  Bool_t fTAPS_ReadPbWO4;               // no general use ... only for my version. testing
  Double_t* fTimeOfFlight;

public:
  TA2Particle* GetParticles(){ return particles; }
  TA2Particle GetParticles(Int_t index){ return particles[index]; }
  UInt_t GetNparticles(){ return NParticles; }
  UInt_t GetNParticles(){ return NParticles; }
  void SetParticleInfo(UInt_t NrParticle);

  TA2TAPS2008( const char*, TA2System* );
  virtual ~TA2TAPS2008();
  virtual void PostInit();
  virtual void SetConfig(char* line, int key);
  virtual TA2DataManager* CreateChild(const char*, Int_t);
  virtual void MakeAllRootinos();
  virtual void LoadVariable();
  virtual void SimplePhotonReconstruct();
  virtual void Reconstruct();
  virtual void MainReconstruct();
  virtual TLorentzVector SetLorentzVector(Double_t, HitCluster_t*);
  virtual void Cleanup();

  TA2ClusterDetector* GetCal(){ return fBaF2; }
  TA2Detector* GetVeto(){ return fVeto; }
  // TLorentzVector* GetP4_Nphoton(){ return fP4_Nphoton; }
  // Double_t* GetM_Nphoton(){ return fM_Nphoton; }
  Int_t CheckParticleID(UInt_t i);

  ClassDef(TA2TAPS2008, 1)
};

//-----------------------------------------------------------------------------

inline void TA2TAPS2008::Reconstruct()
{
  // if((fTapsRunStep==1) && (fTAPS_OUTPUT_FILE_Exists==kFALSE))
  //   fTAPS_OUTPUT_FILE = new TFile("TAPS2008.root", "RECREATE");

  if(fSimpleReconstruct==kTRUE)
    SimplePhotonReconstruct();
  else
    MainReconstruct();

  fTapsRunStep = 2;
}

//-----------------------------------------------------------------------------

inline void TA2TAPS2008::Cleanup()
{
  // Clear any arrays holding variables
  TA2DataManager::Cleanup();

  if(!fNparticle) return;
  Int_t i;
  for(i=0; i<fNparticle; i++)
  {
    fP4[i].SetXYZT(0.,0.,0.,0.);
    fMinv[i] = 0.0;
  }
  // fP4_Nphoton[i-1].SetXYZT(0.,0.,0.,0.);
  // fM_Nphoton[i-1] = 0.0;
  fNparticle = 0;
  fTimeOfFlight[0] = EBufferEnd;
}

//-----------------------------------------------------------------------------

inline void TA2TAPS2008::SimplePhotonReconstruct()
{
 // All in TAPS detected particles are reconstructed as PHOTONS
 // does absolutely the same as TA2Calorimeter does!

  if(fTapsRunStep == 1)
    printf("\n\nINFO: TAPS is doing SimplePhotonReconstruct, only!!!\n\n");

  fNparticle = fBaF2->GetNCluster();
  NParticles = fNparticle;
  if( !fNparticle ) return;
  id_clBaF2 = fBaF2->GetClustHit();
  clBaF2 = fBaF2->GetCluster();

  fP4tot.SetXYZT(0,0,0,0);
  Int_t i;
  for (i=0; i<fNparticle; i++)
  {
    cl = clBaF2[id_clBaF2[i]];
    fPDG_ID[i] = kGamma;
    fPDG_ID_unclear[i] = kTRUE;
    fPDG_ID_sec[i] = kRootino;
    fParticleID->SetP4(fP4+i,fPDG_ID[i],cl->GetEnergy(),cl->GetMeanPosition());
    fP4tot += fP4[i];
    fMinv[i] = fP4[i].M();
    SetParticleInfo(i);
  }
  // fP4_Nphoton[i - 1]= fP4tot;
  // fM_Nphoton[i - 1]= fP4tot.M();
  fMtot = fP4tot.M();
}

//-----------------------------------------------------------------------------

inline TLorentzVector TA2TAPS2008::SetLorentzVector(Double_t mass, HitCluster_t* cl)
{
  // Get cluster kinetic energy and particle direction
  // Construct 4 momentum using input assumed mass

  Double_t E = cl->GetEnergy() + mass;
  Double_t P = TMath::Sqrt(E*E - mass*mass);
  TLorentzVector p4;
  p4.SetE(E);
  p4.SetVect((cl->GetMeanPosition())->Unit() * P);
  return p4;
}

//-----------------------------------------------------------------------------

inline void TA2TAPS2008::SetParticleInfo(UInt_t NrParticle)
{
  particles[NrParticle].SetP4(fP4[NrParticle]);
  particles[NrParticle].SetTime(fBaF2->GetTime(id_clBaF2[NrParticle]));
  particles[NrParticle].SetClusterSize(clBaF2[id_clBaF2[NrParticle]]->GetNhits());
  particles[NrParticle].SetParticleIDA(fPDG_ID[NrParticle]);
  particles[NrParticle].SetSecondID(fPDG_ID_sec[NrParticle]);
  particles[NrParticle].SetUnclear(fPDG_ID_unclear[NrParticle]);
  particles[NrParticle].SetCentralIndex(id_clBaF2[NrParticle]);
  particles[NrParticle].SetApparatus(EAppTAPS);
  particles[NrParticle].SetDetector(EDetBaF2);

  if(fSimpleReconstruct==kFALSE)
  {
    if(fIsVCharged[NrParticle]==kTRUE)
    {
      particles[NrParticle].SetVetoEnergy(fVeto_dE);
      particles[NrParticle].SetVetoTime(fVeto_T);
      particles[NrParticle].SetVetoIndex(fVeto_index);
      particles[NrParticle].SetDetector(EDetVeto);
    }
    particles[NrParticle].SetPSAShort(fShortGateValue);
  }

  particles[NrParticle].SetSigmaE(fBaF2->GetEnergyResolution(clBaF2[id_clBaF2[NrParticle]]->GetEnergy()));
  particles[NrParticle].SetSigmaPhi(fBaF2->GetPhiResolution());
  particles[NrParticle].SetSigmaTheta(fBaF2->GetThetaResolution());
  particles[NrParticle].SetSigmaTime(fBaF2->GetTimeResolution());
}

//-----------------------------------------------------------------------------

#endif
