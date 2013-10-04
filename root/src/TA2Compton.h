//--Author      JRM Annand   18th Feb 2004   Use def physics
//--Rev
//--Rev         JRM Annand   13th May 2004   Start general physics suite
//--Rev         JRM Annand   28th Apr 2004   General photo-meson methods
//--Rev         JRM Annand   13th Jul 2005   SortNPhoton bugs
//--Rev         JRM Annand   25th Jul 2005   ED bug fix fP4tot
//--Rev         JRM Annand   22nd Sep 2006   pi+n analysis add lin pol
//--Rev         JRM Annand   14th Mar 2007   fEmProton for D(g,np) calib.
//--Update      JRM Annand   21st Jul 2008   Compton scattering
//--Update      DL Hornidge   6th Sep 2008   More Compton scattering
//--Update      DL Hornidge  27th May 2009   Threshold Pi0
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2Compton
//
// General reconstruction of reaction kinematics in Mainz tagged-photon
// meson production experiments.
// Use 4-momenta and PDG-index information from apparati to reconstruct
// reaction kinematics. The PDG index (and 4-momentum) assigned by the
// apparatus is not considered binding, e.g. in cases where n/gamma
// discrimination by an apparatus is not possible, in which case it
// defaults to kGamma. The method TA2ParticleID->SetMassP4( *p4, ipdg )
// may be used to reset the rest-mass of an existing 4 momentum *p4 to that
// corresponding to PDG index ipdg.
// This one deals with Compton scattering on the nucleon.
//
// PDG codes of particles generlly observed MAMI-C
// kElectron 11,     kPositron -11
// kMuonMinus 13     kMuonPlus -13      kGamma 22
// kPi0 111          kPiPlus 211        kPiMinus -211       kEta 221
// kProton 2212      kNeutron 2112
// 

#ifndef __TA2Compton_h__
#define __TA2Compton_h__

#include "TAcquRoot.h"
#include "TAcquFile.h"
#include "TA2Physics.h"
#include "TA2Analysis.h"
#include "TA2Tagger.h"
#include "TA2CrystalBall.h"
#include "TA2TAPS2009.h"
#include "TA2Ladder.h"
#include "TA2PhotoPhysics.h"
#include <iostream>

class TA2Tagger;

class TA2Compton : public TA2Physics {

     protected:

          TA2Tagger* fTAGG;                      // Glasgow photon tagger
          TA2CrystalBall* fCB;                   // Crystal Ball
          TA2TAPS2009* fTAPS;                    // TAPS
          TA2Ladder* fLADD;                      // Ladder

//          TA2Particle* fTAGGpart;                // TA2Particles from Tagger
          TA2Particle* fCBpart;                  // TA2Particles from CB
          TA2Particle* fTAPSpart;                // TA2Particles from TAPS

          TLorentzVector* fTAGGp4;               // Tagger 4-mom store

//          TA2Particle** fPARTtaggphot;           // TA2Particle tagger photon
          TA2Particle** fPARTphoton;             // TA2Particle photon
          TA2Particle** fPARTproton;             // TA2Particle proton
          TA2Particle** fPARTpiplus;             // TA2Particle piplus
          TA2Particle** fPARTneutron;            // TA2Particle neutron
          TA2Particle** fPARTrootino;            // TA2Particle rootino

          TA2Particle** fPARTpi0;                // TA2Particle pi0
          TA2Particle** fPARTeta;                // TA2Particle eta
          TA2Particle** fPARTgprime;             // TA2Particle gprime

          TLorentzVector fP4photonTot;           // total 4-momentum of gammas

          Int_t fNphoton;                        // # photon
          Int_t fNproton;                        // # proton
          Int_t fNpiplus;                        // # pi+
          Int_t fNneutron;                       // # neutron
          Int_t fNrootino;                       // # unknowns
          Int_t fNgprime;                        // # unknowns
          Int_t fNprompt;                        // tagger prompts
          Int_t fNrandom;                        // tagger randoms
          Int_t fNparaP;                         // prompts lin-pol para
          Int_t fNparaR;                         // randoms lin-pol para
          Int_t fNperpP;                         // prompts lin-pol perp
          Int_t fNperpR;                         // randoms lin-pol perp
          Int_t fMaxTagg;                        // max # tagger hits

          Int_t fNpi0;                           // # pi0
          Int_t fNeta;                           // # eta

          Double_t* fMassDpi0;                   // for meson ID by inv. mass
          Double_t* fMassDeta;                   // for meson ID by inv. mass
          Int_t* fMassIJ;                        // combinatorial indices
          Int_t* fMassIpi0;                      // ditto
          Int_t* fMassIeta;                      // ditto
          Int_t fMax2gPerm;                      // max # 2-gamma permutations
          Bool_t* fIsMesonIndex;                 // photon derived from a meson?
          Double_t fMaxMDpi0;                    // mass-diff limit pi0
          Double_t fMaxMDeta;                    // mass-diff limit eta

          Double_t fP1, fP2, fRl1, fRl2, fRh1, fRh2;  // Prompt-Random windows
          Double_t fMMMinLo;                     // missing mass limit low E lo
          Double_t fMMMinHi;                     // missing mass limit low E hi
          Double_t fMMMaxLo;                     // missing mass limit high E lo
          Double_t fMMMaxHi;                     // missing mass limit high E hi
          Double_t fOACut;                       // opening-angle cut

          Double_t fM2g;                         // 2-photon invariant mass
          Double_t fM6g;                         // 6-photon invariant mass

          Double_t fProtKE;                      // Proton variables
          Double_t fProtTheta;
          Double_t fProtPhi;
          Double_t fProtTime;

          Double_t fPhotKE;                      // Photon variables
          Double_t fPhotTheta;
          Double_t fPhotPhi;
          Double_t fPhotTime;

          Double_t* fTaggerTime;                 // tagger time
          Double_t* fPhotTaggTime;               // photon-tagger time
          Int_t* fTChanHit;                      // channel hits

          Int_t* fTChanPhot;                     // tagger channel
          Int_t* fTChanPhotP;
          Int_t* fTChanPhotR;
          Double_t* fPhotonMmiss;                // missing mass for the photon
          Double_t* fPhotonMmissP;
          Double_t* fPhotonMmissR;

          Double_t* fPhotonKECM;                 // KE CM of the photon
          Double_t* fPhotonKECMP;
          Double_t* fPhotonKECMR;
          Double_t* fPhotonThetaCM;              // Theta CM of the photon
          Double_t* fPhotonThetaCMP;
          Double_t* fPhotonThetaCMR;
          Double_t* fPhotonPhiCM;                // Phi CM of the photon
          Double_t* fPhotonPhiCMP;
          Double_t* fPhotonPhiCMR;

          Int_t* fTChanPhotProt;                 // tagger channel
          Int_t* fTChanPhotProtP;
          Int_t* fTChanPhotProtR;
          Double_t* fPhotProtOA;
          Double_t* fPhotProtOAP;
          Double_t* fPhotProtOAR;

          Double_t* fPhotonMmissProt;
          Double_t* fPhotonMmissProtP;
          Double_t* fPhotonMmissProtR;
          Double_t* fPhotonEmissProt;
          Double_t* fPhotonEmissProtP;
          Double_t* fPhotonEmissProtR;
          Double_t* fPhotonThetaCMProt;
          Double_t* fPhotonThetaCMProtP;
          Double_t* fPhotonThetaCMProtR;

          Double_t* fPhotonMmissProtOA;
          Double_t* fPhotonMmissProtOAP;
          Double_t* fPhotonMmissProtOAR;
          Double_t* fTChanPhotProtOA;
          Double_t* fTChanPhotProtOAP;
          Double_t* fTChanPhotProtOAR;
          Double_t* fPhotonThetaCMProtOA;
          Double_t* fPhotonThetaCMProtOAP;
          Double_t* fPhotonThetaCMProtOAR;

     public:

          TA2Compton( const char*, TA2Analysis* );
          virtual ~TA2Compton();
          virtual void LoadVariable();           // variables for display/cuts
          virtual void PostInit( );              // init after parameter input
          virtual void SetConfig(Char_t*, Int_t);

          virtual void Reconstruct();            // reconstruct detector info
          virtual TA2DataManager* CreateChild(const char*,Int_t){return NULL;}

          virtual void Sort2Photon( );
          virtual void SortNPhoton( );
          virtual void MarkUndefined( Int_t );
          virtual void MarkEndBuffer( );

          Double_t Sqr( Double_t);

     ClassDef(TA2Compton,1)
};

//-----------------------------------------------------------------------------
inline void TA2Compton::Sort2Photon()
{
     // Test if 2 gamma 4-momenta combine to give a  possible pi0 or eta
     // by finding the invariant mass. This is a fast version of SortNPhoton
     // for the frequent pi0 or eta -> 2 gamma situation

     TA2Particle phot1 = *fPARTphoton[0];
     TA2Particle phot2 = *fPARTphoton[1];

     Double_t time = (phot1.GetTime() + phot2.GetTime())/2;

     TLorentzVector p4 = phot1.GetP4() + phot2.GetP4();  // sum 4 momenta

     fM2g = p4.M();                                      // inv. mass

     fMassDpi0[0] = TMath::Abs( p4.M() - fParticleID->GetMassMeV( kPi0));
     fMassDeta[0] = TMath::Abs( p4.M() - fParticleID->GetMassMeV( kEta));

     if ( fMassDpi0[0] < fMaxMDpi0 ) {                   // is it a pi0
//          printf( "%6.1f  %6.1f  %6.1f  %6.1f\n",fMaxMDpi0,fMassDpi0[0],fM2g,
//                    time);
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

inline void TA2Compton::SortNPhoton()
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
          fIsMesonIndex[i] = EFalse;                // initialise not meson

          for( j = i + 1; j < n; j++ ){

               // add the 4 momenta
               p4temp =  (*fPARTphoton[i]).GetP4()+(*fPARTphoton[j]).GetP4();
               *pi0diff++ = TMath::Abs( p4temp.M()-mPi0 );//check inv mass pi0
               *etadiff++ = TMath::Abs( p4temp.M()-mEta );//check inv. mass eta
               *ij++ = i | (j<<16);                 // store i,j indices
               k++;                                 // permutation counter

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

inline Double_t TA2Compton::Sqr( Double_t x)
{
     return( x*x);
}

//-----------------------------------------------------------------------------

inline void TA2Compton::MarkUndefined( Int_t jtagg )
{
     // Initialise undefined
     fTChanPhotP[jtagg] = EPhotoUndefined;
     fPhotonMmissP[jtagg] = EPhotoUndefined;

     fPhotonKECMP[jtagg] = EPhotoUndefined;
     fPhotonThetaCMP[jtagg] = EPhotoUndefined;
     fPhotonPhiCMP[jtagg] = EPhotoUndefined;

     fTChanPhotProtP[jtagg] = EPhotoUndefined;
     fPhotProtOAP[jtagg] = EPhotoUndefined;

     fPhotonMmissProtP[jtagg] = EPhotoUndefined;
     fPhotonEmissProtP[jtagg] = EPhotoUndefined;
     fPhotonThetaCMProtP[jtagg] = EPhotoUndefined;

     fPhotonMmissProtOAP[jtagg] = EPhotoUndefined;
     fTChanPhotProtOAP[jtagg] = EPhotoUndefined;
     fPhotonThetaCMProtOAP[jtagg] = EPhotoUndefined;
}

//-----------------------------------------------------------------------------

inline void TA2Compton::MarkEndBuffer()
{
     // Ensure the multi-data buffers are marked as ended
     fTChanPhotP[fNprompt] = EBufferEnd;
     fTChanPhotR[fNrandom] = EBufferEnd;
     fPhotonMmissP[fNprompt] = EBufferEnd;
     fPhotonMmissR[fNrandom] = EBufferEnd;  

     fPhotonKECMP[fNprompt] = EBufferEnd;
     fPhotonKECMR[fNrandom] = EBufferEnd;
     fPhotonThetaCMP[fNprompt] = EBufferEnd;
     fPhotonThetaCMR[fNrandom] = EBufferEnd;
     fPhotonPhiCMP[fNprompt] = EBufferEnd;
     fPhotonPhiCMR[fNrandom] = EBufferEnd;

     fTChanPhotProtP[fNprompt] = EBufferEnd;
     fTChanPhotProtR[fNrandom] = EBufferEnd;
     fPhotProtOAP[fNprompt] = EBufferEnd;
     fPhotProtOAR[fNrandom] = EBufferEnd;

     fPhotonMmissProtP[fNprompt] = EBufferEnd;
     fPhotonMmissProtR[fNrandom] = EBufferEnd;
     fPhotonEmissProtP[fNprompt] = EBufferEnd;
     fPhotonEmissProtR[fNrandom] = EBufferEnd;
     fPhotonThetaCMProtP[fNprompt] = EBufferEnd;
     fPhotonThetaCMProtR[fNrandom] = EBufferEnd;

     fPhotonMmissProtOAP[fNprompt] = EBufferEnd;
     fPhotonMmissProtOAR[fNrandom] = EBufferEnd;
     fTChanPhotProtOAP[fNprompt] = EBufferEnd;
     fTChanPhotProtOAR[fNrandom] = EBufferEnd;
     fPhotonThetaCMProtOAP[fNprompt] = EBufferEnd;
     fPhotonThetaCMProtOAR[fNrandom] = EBufferEnd;
}

#endif
