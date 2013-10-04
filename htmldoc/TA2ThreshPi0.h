//--Author     JRM Annand   18th Feb 2004   Use def physics
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
// TA2ThreshPi0
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
// This one deals with pion photoproduction on the nucleon.
//
// PDG codes of particles generlly observed MAMI-B
// kElectron 11,     kPositron -11
// kMuonMinus 13     kMuonPlus -13      kGamma 22
// kPi0 111          kPiPlus 211        kPiMinus -211       kEta 221
// kProton 2212      kNeutron 2112
// 

#ifndef __TA2ThreshPi0_h__
#define __TA2ThreshPi0_h__

#include "TAcquRoot.h"
#include "TAcquFile.h"
#include "TA2Physics.h"
#include "TA2Analysis.h"
#include "TA2CrystalBall.h"
#include "TA2TAPS2009.h"
#include "TA2Tagger.h"
#include "TA2Ladder.h"
#include "TA2PhotoPhysics.h"
#include <iostream>

class TA2Tagger;

class TA2ThreshPi0 : public TA2Physics {

     protected:

          TA2Tagger* fTAGG;                     // Glasgow photon tagger
          TA2CrystalBall* fCB;                  // Crystal Ball
          TA2TAPS2009* fTAPS;                   // TAPS
          TA2Ladder* fLADD;                     // Ladder

//          TA2Particle* fTAGGpart;               // TA2Particles from Tagger
          TA2Particle* fCBpart;                 // TA2Particles from CB
          TA2Particle* fTAPSpart;               // TA2Particles from TAPS

          TLorentzVector* fTAGGp4;              // Tagger 4-mom store

//          TA2Particle** fPARTtaggphot;          // TA2Particle tagger photon
          TA2Particle** fPARTphoton;            // TA2Particle photon
          TA2Particle** fPARTproton;            // TA2Particle proton
          TA2Particle** fPARTpiplus;            // TA2Particle piplus
          TA2Particle** fPARTneutron;           // TA2Particle neutron
          TA2Particle** fPARTrootino;           // TA2Particle rootino

          TA2Particle** fPARTpi0;               // TA2Particle pi0
          TA2Particle** fPARTeta;               // TA2Particle eta
          TA2Particle** fPARTgprime;            // TA2Particle gprime

          TLorentzVector fP4photonTot;          // total 4-momentum of gammas

          Int_t fNphoton;                       // # photon
          Int_t fNproton;                       // # proton
          Int_t fNpiplus;                       // # pi+
          Int_t fNneutron;                      // # neutron
          Int_t fNrootino;                      // # unknowns
          Int_t fNgprime;                       // # unknowns
          Int_t fNprompt;                       // tagger prompts
          Int_t fNrandom;                       // tagger randoms
          Int_t fNparaP;                        // prompts lin-pol para
          Int_t fNparaR;                        // randoms lin-pol para
          Int_t fNperpP;                        // prompts lin-pol perp
          Int_t fNperpR;                        // randoms lin-pol perp
          Int_t fMaxTagg;                       // max # tagger hits

          Int_t fNpi0;                          // # pi0
          Int_t fNeta;                          // # eta

          Double_t* fMassDpi0;                  // for meson ID by inv. mass
          Double_t* fMassDeta;                  // for meson ID by inv. mass
          Int_t* fMassIJ;                       // combinatorial indices
          Int_t* fMassIpi0;                     // ditto
          Int_t* fMassIeta;                     // ditto
          Int_t fMax2gPerm;                     // max # 2-gamma permutations
          Bool_t* fIsMesonIndex;                // photon derived from a meson?
          Double_t fMaxMDpi0;                   // mass-diff limit pi0
          Double_t fMaxMDeta;                   // mass-diff limit eta

          Double_t fP1, fP2, fRl1, fRl2, fRh1, fRh2; // Prompt-Random windows
          Int_t fMassCorr;                      // pi0 invmass correction

          Double_t fMMWlo[2], fMMWhi[2];
          Double_t TggCutOffset;

          Double_t fM2g;                      // 2-photon invariant mass
          Double_t fM6g;                      // 6-photon invariant mass

          Double_t fTimeDiff;                 // time difference between photons
          Double_t fTime1;                    // time of photon1
          Double_t fTime2;                    // time of photon2

          Double_t fPi0TGG;
          Double_t fPi0LabKE;
          Double_t fPi0LabTheta;
          Double_t fPi0LabPhi;
          Double_t fPi0CMKE;
          Double_t fPi0CMTheta;
          Double_t fPi0CMPhi;

          Double_t* fTaggerTime;                // tagger time
          Double_t* fPi0TaggTime;               // pi0 tagger time
          Int_t* fTChanHit;                     // channel hits

          Int_t* fTChan;                        // tagger channel
          Int_t* fTChanP;
          Int_t* fTChanR;
          Double_t* fMMiss;                     // missing mass of the pion
          Double_t* fMMissP;
          Double_t* fMMissR;
          Double_t* fTGG;                       // TGG
          Double_t* fTGGP;
          Double_t* fTGGR;
          Double_t* fM2gTagg;                   // tagger channel
          Double_t* fM2gTaggP;
          Double_t* fM2gTaggR;

          Double_t* fKELab;                     // KE lab of the pion
          Double_t* fKELabP;
          Double_t* fKELabR;
          Double_t* fThetaLab;                  // theta lab of the pion
          Double_t* fThetaLabP;
          Double_t* fThetaLabR;
          Double_t* fPhiLab;                    // phi Lab of the pion
          Double_t* fPhiLabP;
          Double_t* fPhiLabR;

          Double_t* fKECM;                      // KE CM of the pion
          Double_t* fKECMP;
          Double_t* fKECMR;
          Double_t* fThetaCM;                   // theta CM of the pion
          Double_t* fThetaCMP;
          Double_t* fThetaCMR;
          Double_t* fPhiCM;                      // phi CM of the pion
          Double_t* fPhiCMP;
          Double_t* fPhiCMR;

          Int_t* fTChanCut1;                     // tagger channel
          Int_t* fTChanCut1P;
          Int_t* fTChanCut1R;
          Double_t* fKECMCut1;                   // KE CM of the pion
          Double_t* fKECMCut1P;
          Double_t* fKECMCut1R;
          Double_t* fThetaCMCut1;                // theta CM of the pion
          Double_t* fThetaCMCut1P;
          Double_t* fThetaCMCut1R;
          Double_t* fMMissCut1;                  // missing mass of the pion
          Double_t* fMMissCut1P;
          Double_t* fMMissCut1R;

          Int_t* fTChanCut2;                     // tagger channel
          Int_t* fTChanCut2P;
          Int_t* fTChanCut2R;
          Double_t* fKECMCut2;                   // KE CM of the pion
          Double_t* fKECMCut2P;
          Double_t* fKECMCut2R;
          Double_t* fThetaCMCut2;                // theta CM of the pion
          Double_t* fThetaCMCut2P;
          Double_t* fThetaCMCut2R;
          Double_t* fPhiCMCut2;                  // phi CM of the pion
          Double_t* fPhiCMCut2P;
          Double_t* fPhiCMCut2R;

     public:

          TA2ThreshPi0( const char*, TA2Analysis* );
          virtual ~TA2ThreshPi0();
          virtual void LoadVariable();           // variables for display/cuts
          virtual void PostInit( );              // init after parameter input
          virtual void SetConfig(Char_t*, Int_t);

          virtual void Reconstruct();            // reconstruct detector info
          virtual TA2DataManager* CreateChild(const char*, Int_t){ return NULL;}

          virtual void Sort2Photon( );
          virtual void SortNPhoton( );
          virtual void MarkUndefined( Int_t );
          virtual void MarkEndBuffer( );

          Double_t Sqr( Double_t);
          Double_t Tgg_Min( Double_t, Double_t);
          Double_t Energy( Double_t, Double_t);
          Double_t Momentum( Double_t, Double_t);
          Double_t qp_thcm( Double_t, Double_t, Double_t, Double_t);
          Double_t qT_max( Double_t, Double_t, Double_t);
          Double_t Linear( Double_t, Double_t, Double_t, Double_t, Double_t);

     ClassDef(TA2ThreshPi0,1)
};

//-----------------------------------------------------------------------------
inline void TA2ThreshPi0::Sort2Photon()
{
     // Test if 2 gamma 4-momenta combine to give a  possible pi0 or eta
     // by finding the invariant mass. This is a fast version of SortNPhoton
     // for the frequent pi0 or eta -> 2 gamma situation

     TA2Particle phot1 = *fPARTphoton[0];
     TA2Particle phot2 = *fPARTphoton[1];

     Double_t time = (phot1.GetTime() + phot2.GetTime())/2;

     TLorentzVector p4 = phot1.GetP4() + phot2.GetP4();        // sum 4 momenta

     fM2g = p4.M();                                            // inv. mass

     fMassDpi0[0] = TMath::Abs( p4.M() - fParticleID->GetMassMeV( kPi0));
     fMassDeta[0] = TMath::Abs( p4.M() - fParticleID->GetMassMeV( kEta));

     if ( fMassDpi0[0] < fMaxMDpi0 ) {                         // is it a pi0
//        printf( "%6.1f  %6.1f  %6.1f  %6.1f\n", fMaxMDpi0, fMassDpi0[0], fM2g,
//                   time);
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

inline void TA2ThreshPi0::SortNPhoton()
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
             p4temp =  (*fPARTphoton[i]).GetP4() + (*fPARTphoton[j]).GetP4();
             *pi0diff++ = TMath::Abs( p4temp.M() - mPi0 ); //check inv mass pi0
             *etadiff++ = TMath::Abs( p4temp.M() - mEta ); //check inv. mass eta
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

inline Double_t TA2ThreshPi0::Sqr( Double_t x)
{
     return( x*x);
}

// Calculates the minimum opening angle of the decay photons for a pi0 as a
// function of pion KE and mass.
//
// Returns angle in RADIANS!
//
inline Double_t TA2ThreshPi0::Tgg_Min( Double_t T, Double_t mpi)
{
     Double_t term, tgg_min;

     term = mpi/sqrt( T*(T + 2*mpi));
     tgg_min = 2*atan( term);

     return( tgg_min);
}

// Calculates the energy of a particle given the mass and three-mom.
inline Double_t TA2ThreshPi0::Energy( Double_t mom, Double_t m)
{
     return( sqrt(mom*mom + m*m));
}

// Calculates relativistic momentum from mass and energy.
inline Double_t TA2ThreshPi0::Momentum( Double_t en, Double_t m)
{
     if ( en >= m) return( sqrt( en*en - m*m));
     else return( -1);
}

// Calculates the lab momentum of the pion.
// This one takes the pion CM theta, but still calculates the lab momentum!
inline Double_t TA2ThreshPi0::qp_thcm( Double_t ke, Double_t pm,
          Double_t qth_cm, Double_t qm)
{

     TLorentzVector k, p, q, pIn;
     TLorentzVector k_cm, p_cm, q_cm;
     TVector3 cmBoost, labBoost, q3mom;
     
     Double_t qph_cm, mom_cm;
     Double_t qx_cm, qy_cm, qz_cm, qe_cm;
     Double_t S;

     k.SetPxPyPzE( 0, 0, ke, ke);
     p.SetPxPyPzE( 0, 0, 0, pm);

     pIn = k + p;
     labBoost = pIn.BoostVector();
     cmBoost = -pIn.BoostVector();

     k_cm = k;
     k_cm.Boost( cmBoost);
     p_cm = p;
     p_cm.Boost( cmBoost);

     qth_cm *= TMath::DegToRad();
     qph_cm = 0;

     S = 2*ke*pm + pm*pm;
     qe_cm = (S - pm*pm + qm*qm)/2/sqrt(S);

     mom_cm = Momentum( qe_cm, qm);
     qx_cm = mom_cm*sin( qth_cm)*cos(qph_cm);
     qy_cm = mom_cm*sin( qth_cm)*sin(qph_cm);
     qz_cm = mom_cm*cos( qth_cm);

     q_cm.SetPxPyPzE( qx_cm, qy_cm, qz_cm, qe_cm);

     q = q_cm;
     q.Boost( labBoost);

     q3mom = q.Vect();

     return( q3mom.Mag());
}

// Calculates maximum KE of the pion.
inline Double_t TA2ThreshPi0::qT_max( Double_t eg, Double_t pm, Double_t qm)
{
     Double_t W, W_cm, S, beta, gamma;
     Double_t qcm_e, qcm_p, qcm_x, qcm_z;
     Double_t qx, qz, qp, qe, qT_max, th_cm;

     W = eg + pm;
     S = Sqr( pm) + 2*eg*pm;
     W_cm = sqrt( S);

     beta = eg/W;
     gamma = W/W_cm;

     qcm_e = (S + Sqr(qm) - Sqr(pm))/2/W_cm;
     qcm_p = Momentum( qcm_e, qm);

     th_cm = 0;  // Condition for max pion KE

     qcm_x = qcm_p*sin( th_cm);
     qcm_z = qcm_p*cos( th_cm);
     qx = qcm_x;
     qz = gamma*(qcm_z + beta*qcm_e);
     qp = sqrt( Sqr( qx) + Sqr(qz));
     qe = Energy( qp, qm);
     qT_max = qe - qm;

     return( qT_max);
}

inline Double_t TA2ThreshPi0::Linear( Double_t x0, Double_t y0, Double_t x1,
          Double_t y1, Double_t x)
{
     return( y0 + (y1 - y0)/(x1 - x0)*(x - x0));
}


//-----------------------------------------------------------------------------

inline void TA2ThreshPi0::MarkUndefined( Int_t jtagg )
{
     // Initialise undefined
     fTChanP[jtagg] = EPhotoUndefined;
     fMMissP[jtagg] = EPhotoUndefined;
     fTGGP[jtagg] = EPhotoUndefined;
     fM2gTaggP[jtagg] = EPhotoUndefined;

     fKELabP[jtagg] = EPhotoUndefined;
     fThetaLabP[jtagg] = EPhotoUndefined;
     fPhiLabP[jtagg] = EPhotoUndefined;

     fKECMP[jtagg] = EPhotoUndefined;
     fThetaCMP[jtagg] = EPhotoUndefined;
     fPhiCMP[jtagg] = EPhotoUndefined;

     fTChanCut1P[jtagg] = EPhotoUndefined;
     fKECMCut1P[jtagg] = EPhotoUndefined;
     fThetaCMCut1P[jtagg] = EPhotoUndefined;
     fMMissCut1P[jtagg] = EPhotoUndefined;

     fTChanCut2P[jtagg] = EPhotoUndefined;
     fKECMCut2P[jtagg] = EPhotoUndefined;
     fThetaCMCut2P[jtagg] = EPhotoUndefined;
     fPhiCMCut2P[jtagg] = EPhotoUndefined;
}

//-----------------------------------------------------------------------------

inline void TA2ThreshPi0::MarkEndBuffer()
{
     // Ensure the multi-data buffers are marked as ended
     fTChanP[fNprompt] = EBufferEnd;
     fTChanR[fNrandom] = EBufferEnd;
     fMMissP[fNprompt] = EBufferEnd;
     fMMissR[fNrandom] = EBufferEnd;  
     fTGGP[fNprompt] = EBufferEnd;
     fTGGR[fNrandom] = EBufferEnd;
     fM2gTaggP[fNprompt] = EBufferEnd;
     fM2gTaggR[fNrandom] = EBufferEnd;

     fKELabP[fNprompt] = EBufferEnd;
     fKELabR[fNrandom] = EBufferEnd;
     fThetaLabP[fNprompt] = EBufferEnd;
     fThetaLabR[fNrandom] = EBufferEnd;
     fPhiLabP[fNprompt] = EBufferEnd;
     fPhiLabR[fNrandom] = EBufferEnd;

     fKECMP[fNprompt] = EBufferEnd;
     fKECMR[fNrandom] = EBufferEnd;
     fThetaCMP[fNprompt] = EBufferEnd;
     fThetaCMR[fNrandom] = EBufferEnd;
     fPhiCMP[fNprompt] = EBufferEnd;
     fPhiCMR[fNrandom] = EBufferEnd;

     fTChanCut1P[fNprompt] = EBufferEnd;
     fTChanCut1R[fNrandom] = EBufferEnd;
     fKECMCut1P[fNprompt] = EBufferEnd;
     fKECMCut1R[fNrandom] = EBufferEnd;
     fThetaCMCut1P[fNprompt] = EBufferEnd;
     fThetaCMCut1R[fNrandom] = EBufferEnd;
     fMMissCut1P[fNprompt] = EBufferEnd;
     fMMissCut1R[fNrandom] = EBufferEnd;  

     fTChanCut2P[fNprompt] = EBufferEnd;
     fTChanCut2R[fNrandom] = EBufferEnd;
     fKECMCut2P[fNprompt] = EBufferEnd;
     fKECMCut2R[fNrandom] = EBufferEnd;
     fThetaCMCut2P[fNprompt] = EBufferEnd;
     fThetaCMCut2R[fNrandom] = EBufferEnd;
     fPhiCMCut2P[fNprompt] = EBufferEnd;
     fPhiCMCut2R[fNrandom] = EBufferEnd;
}

#endif
