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
//--Update      AT Laffoley 11th June 2009   Converting the DLH way to JRMA way
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

#include "TA2ThreshPi0.h"

// Valid Keywords for command-line setup of ThreshPi0
enum { EThreshMassLimits = 1000, EThreshPRLimits, EThreshMassCorr,
     EThreshMMissWin, EThreshTggCutOffset};
static const Map_t kPhotoKeys[] = {
     {"Mass-Limits:",           EThreshMassLimits},
     {"Prompt-Random-Windows:", EThreshPRLimits},
     {"Mass-Corr:",             EThreshMassCorr},
     {"Missing-Mass-Windows:",  EThreshMMissWin},
     {"TGG-Cut-Offset:",        EThreshTggCutOffset},
     {NULL,            -1}
};

ClassImp(TA2ThreshPi0)

//-----------------------------------------------------------------------------
TA2ThreshPi0::TA2ThreshPi0( const char* name, TA2Analysis* analysis )
     :TA2Physics( name, analysis ) {
     // Initialise ThreshPi0 variables here
     // Default null pointers, zeroed variables

     fTAGG = NULL;
     fCB = NULL;
     fTAPS = NULL;

     fTAGGp4 = NULL;

//     fPARTtaggphot = NULL;
     fPARTphoton = NULL;
     fPARTproton = NULL;
     fPARTpiplus = NULL;
     fPARTneutron = NULL;
     fPARTrootino = NULL;

     fPARTpi0 = NULL;
     fPARTeta = NULL;
     fPARTgprime = NULL;

     fNphoton = 0;
     fNproton = 0;
     fNpiplus = 0;
     fNneutron = 0;
     fNrootino = 0;

     fNpi0 = 0;
     fNeta = 0;
     fNgprime = 0;

     fMassDpi0 = NULL;
     fMassDeta = NULL;
     fMassIJ = NULL;
     fMassIpi0 = NULL;
     fMassIeta = NULL;
     fIsMesonIndex = NULL;

     fMaxMDpi0 = 0.0;
     fMaxMDeta = 0.0;

     fM2g = ENullHit;
     fM6g = ENullHit;

     fPi0TGG = ENullHit;
     fPi0LabKE = ENullHit;
     fPi0LabTheta = ENullHit;
     fPi0LabPhi = ENullHit;
     fPi0CMKE = ENullHit;
     fPi0CMTheta = ENullHit;
     fPi0CMPhi = ENullHit;

     fTimeDiff = ENullHit;
     fTime1 = ENullHit;
     fTime2 = ENullHit;

     fNprompt = 0;
     fNrandom = 0;
     fMaxTagg = 0;
     fMax2gPerm = 0;

     fTaggerTime = NULL;
     fPi0TaggTime = NULL;
     fTChanHit = NULL;

     fTChanP = NULL;
     fTChanR = NULL;
     fMMissP = NULL;
     fMMissR = NULL;
     fTGGP = NULL;
     fTGGR = NULL;
     fM2gTaggP = NULL;
     fM2gTaggR = NULL;

     fKELabP = NULL;
     fKELabR = NULL;
     fThetaLabP = NULL;
     fThetaLabR = NULL;
     fPhiLabP = NULL;
     fPhiLabR = NULL;

     fKECMP = NULL;
     fKECMR = NULL;
     fThetaCMP = NULL;
     fThetaCMR = NULL;
     fPhiCMP = NULL;
     fPhiCMR = NULL;

     fTChanCut1P = NULL;
     fTChanCut1R = NULL;
     fKECMCut1P = NULL;
     fKECMCut1R = NULL;
     fThetaCMCut1P = NULL;
     fThetaCMCut1R = NULL;
     fMMissCut1P = NULL;
     fMMissCut1R = NULL;

     fTChanCut2P = NULL;
     fTChanCut2R = NULL;
     fKECMCut2P = NULL;
     fKECMCut2R = NULL;
     fThetaCMCut2P = NULL;
     fThetaCMCut2R = NULL;
     fPhiCMCut2P = NULL;
     fPhiCMCut2R = NULL;

     AddCmdList( kPhotoKeys );       // command-line recognition for SetConfig()
}


//-----------------------------------------------------------------------------
TA2ThreshPi0::~TA2ThreshPi0()
{
     // Free up allocated memory...after checking its allocated
     // detector and cuts lists
}

//-----------------------------------------------------------------------------
void TA2ThreshPi0::SetConfig(Char_t* line, Int_t key)
{
  // Any special command-line input for Crystal Ball apparatus

  switch (key){
    case EThreshMassLimits:
         //  Invariant mass limits
         if( sscanf( line, "%lf%lf", &fMaxMDpi0, &fMaxMDeta ) != 2 ){
              PrintError( line, "<ThreshPi0 meson invariant mass limits>");
              return;
         }
         break;
    case EThreshPRLimits:
         //  Prompt-Random Windows
         if( sscanf( line, "%lf%lf%lf%lf%lf%lf", &fP1, &fP2, &fRl1, &fRl2,
                        &fRh1, &fRh2) != 6 ){
              PrintError( line, "<ThreshPi0 prompt-random windows>");
              return;
         }
         break;
    case EThreshMassCorr:
         //  Pi0 Invariant Mass Correction
         if( sscanf( line, "%d", &fMassCorr) != 1 ){
              PrintError( line, "<ThreshPi0 pi0 invariant mass correction>");
              return;
         }
         break;
    case EThreshMMissWin:
         //  Target Missing-Mass Windows
         if( sscanf( line, "%lf%lf%lf%lf", &fMMWlo[0], &fMMWhi[0], &fMMWlo[1],
                        &fMMWhi[1]) != 4 ){
              PrintError( line, "<ThreshPi0 missing-mass windows>");
              return;
         }
         break;
    case EThreshTggCutOffset:
         //  TGG Cut Offset
         if( sscanf( line, "%lf", &TggCutOffset) != 1 ){
              PrintError( line, "<ThreshPi0 Tgg cut offset>");
              return;
         }
         break;
    default:
         // default main apparatus SetConfig()
         TA2Physics::SetConfig( line, key );
         break;
  }
}

//---------------------------------------------------------------------------
void TA2ThreshPi0::PostInit()
{
   // Initialise arrays to contain 4 momenta and plotable scaler variables
   // Missing mass, missing energy, cm momentum, energies, angles
   // Initialisation will abort if CB or Tagger not initialised
   // TAPS is optional

     fCB = (TA2CrystalBall*)((TA2Analysis*)fParent)->GetChild("CB");
     if ( !fCB) PrintError("","<No CB class found in annalysis>", EErrFatal);
     else {
        fCBpart = fCB->GetParticles();
   }

//   fTAGG = (TA2KensTagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
//   if(!fTAGG) PrintError("","<No Tagger class found in annalysis>",EErrFatal);
//   else{
//        fTAGGpart = fTAGG->GetParticles();
//   }
     fTAGG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
     if (!fTAGG) PrintError("","<No Tagger class found in annalysis>",EErrFatal);
     else {
          fTAGGp4 = fTAGG->GetP4();
     }
     fLADD = (TA2Ladder*)((TA2Analysis*)fParent)->GetGrandChild( "FPD");
     if ( !fLADD)
          PrintError( "", "<No Ladder class found in analysis>", EErrFatal);

     fTAPS = (TA2TAPS2009*)((TA2Analysis*)fParent)->GetChild("TAPS");
     if ( !fTAPS) PrintError("Warning!!!","<No TAPS class found in annalysis>");
     else {
          fTAPSpart = fTAPS->GetParticles();
     }

     Int_t i;
     TA2Particle* part;

     // Maximum # of tagger hits
     Int_t maxtagg = fTAGG->GetMaxParticle() + 1;
     fMaxTagg = maxtagg;

     fTaggerTime = new Double_t[maxtagg];
     fPi0TaggTime = new Double_t[maxtagg];
     fTChanHit = new Int_t[maxtagg];

     fTChan = new Int_t[maxtagg * 2];
     fTChanP = fTChan; fTChanR = fTChan + maxtagg;
     fMMiss = new Double_t[maxtagg * 2];
     fMMissP = fMMiss; fMMissR = fMMiss + maxtagg;
     fTGG = new Double_t[maxtagg * 2];
     fTGGP = fTGG; fTGGR = fTGG + maxtagg;
     fM2gTagg = new Double_t[maxtagg * 2];
     fM2gTaggP = fM2gTagg; fM2gTaggR = fM2gTagg + maxtagg;

     fKELab = new Double_t[maxtagg * 2];
     fKELabP = fKELab; fKELabR = fKELab + maxtagg;
     fPhiLab = new Double_t[maxtagg * 2];
     fPhiLabP = fPhiLab; fPhiLabR = fPhiLab + maxtagg;
     fThetaLab = new Double_t[maxtagg * 2];
     fThetaLabP = fThetaLab; fThetaLabR = fThetaLab + maxtagg;

     fKECM = new Double_t[maxtagg * 2];
     fKECMP = fKECM; fKECMR = fKECM + maxtagg;
     fThetaCM = new Double_t[maxtagg * 2];
     fThetaCMP = fThetaCM; fThetaCMR = fThetaCM + maxtagg;
     fPhiCM = new Double_t[maxtagg * 2];
     fPhiCMP = fPhiCM; fPhiCMR = fPhiCM + maxtagg;

     fTChanCut1 = new Int_t[maxtagg * 2];
     fTChanCut1P = fTChanCut1; fTChanCut1R = fTChanCut1 + maxtagg;
     fKECMCut1 = new Double_t[maxtagg * 2];
     fKECMCut1P = fKECMCut1; fKECMCut1R = fKECMCut1 + maxtagg;
     fThetaCMCut1 = new Double_t[maxtagg * 2];
     fThetaCMCut1P = fThetaCMCut1; fThetaCMCut1R = fThetaCMCut1 + maxtagg;
     fMMissCut1 = new Double_t[maxtagg * 2];
     fMMissCut1P = fMMissCut1; fMMissCut1R = fMMissCut1 + maxtagg;

     fTChanCut2 = new Int_t[maxtagg * 2];
     fTChanCut2P = fTChanCut2; fTChanCut2R = fTChanCut2 + maxtagg;
     fKECMCut2 = new Double_t[maxtagg * 2];
     fKECMCut2P = fKECMCut2; fKECMCut2R = fKECMCut2 + maxtagg;
     fThetaCMCut2 = new Double_t[maxtagg * 2];
     fThetaCMCut2P = fThetaCMCut2; fThetaCMCut2R = fThetaCMCut2 + maxtagg;
     fPhiCMCut2 = new Double_t[maxtagg * 2];
     fPhiCMCut2P = fPhiCMCut2; fPhiCMCut2R = fPhiCMCut2 + maxtagg;
/*
     fPARTtaggphot = new TA2Particle*[maxtagg];
     part = new TA2Particle[maxtagg];
     for ( i = 0; i < maxtagg; i++) fPARTtaggphot[i] = part + i;
*/

     // Maximum # of reaction particles
     Int_t maxparticle = fCB->GetMaxParticle();
     if( fTAPS ) maxparticle += fTAPS->GetMaxParticle();

     // Particle from detectors
     fPARTphoton = new TA2Particle*[maxparticle];
     fPARTproton = new TA2Particle*[maxparticle];
     fPARTpiplus = new TA2Particle*[maxparticle];
     fPARTneutron = new TA2Particle*[maxparticle];
     fPARTrootino = new TA2Particle*[maxparticle];

     // Pi0
     fPARTpi0 = new TA2Particle*[maxparticle];
     part = new TA2Particle[maxparticle];
     for ( i = 0; i < maxparticle; i++) fPARTpi0[i] = part + i;

     // Eta
     fPARTeta = new TA2Particle*[maxparticle];
     part = new TA2Particle[maxparticle];
     for ( i = 0; i < maxparticle; i++) fPARTeta[i] = part + i;

     // Gamma prime
     fPARTgprime =  new TA2Particle*[maxparticle];
     part = new TA2Particle[maxparticle];
     for ( i = 0; i < maxparticle; i++) fPARTgprime[i] = part + i;

     // Arrays used to combine photons to mesons
     Int_t maxperm = 0;
     for ( i = 1; i <= maxparticle; i++) maxperm += i;
     fMassDpi0 = new Double_t[maxperm];
     fMassDeta = new Double_t[maxperm];
     fMassIJ = new Int_t[maxperm];
     fMassIpi0 = new Int_t[maxperm];
     fMassIeta = new Int_t[maxperm];
     fMax2gPerm = maxperm;
     fIsMesonIndex = new Bool_t[maxparticle];

     // Default physics initialisation
     TA2Physics::PostInit();
}

//-----------------------------------------------------------------------------
void TA2ThreshPi0::LoadVariable( )
{
     // Input name - variable pointer associations for any subsequent
     // cut or histogram setup
     // LoadVariable( "name", pointer-to-variable, type-spec );
     // NB scaler variable pointers need the preceeding &
     //     array variable pointers do not.
     // type-spec ED prefix for a Double_t variable
     //                     EI prefix for an Int_t variable
     // type-spec SingleX for a single-valued variable
     //                     MultiX  for a multi-valued variable

     TA2Physics::LoadVariable();
     TA2DataManager::LoadVariable("Nphoton",    &fNphoton,     EISingleX);
     TA2DataManager::LoadVariable("Nproton",    &fNproton,     EISingleX);
     TA2DataManager::LoadVariable("Npiplus",    &fNpiplus,     EISingleX);
     TA2DataManager::LoadVariable("Nneutron",   &fNneutron,    EISingleX);
     TA2DataManager::LoadVariable("Nrootino",   &fNrootino,    EISingleX);
     TA2DataManager::LoadVariable("Ngprime",    &fNgprime,     EISingleX);
     TA2DataManager::LoadVariable("Npi0",       &fNpi0,        EISingleX);
     TA2DataManager::LoadVariable("Neta",       &fNeta,        EISingleX);
     TA2DataManager::LoadVariable("M2g",        &fM2g,         EDSingleX);
     TA2DataManager::LoadVariable("M6g",        &fM6g,         EDSingleX);

     TA2DataManager::LoadVariable("TimeDiff",   &fTimeDiff,    EDSingleX);
     TA2DataManager::LoadVariable("Time1",      &fTime1,       EDSingleX);
     TA2DataManager::LoadVariable("Time2",      &fTime2,       EDSingleX);

     TA2DataManager::LoadVariable("Pi0TGG",     &fPi0TGG,      EDSingleX);
     TA2DataManager::LoadVariable("Pi0LabKE",   &fPi0LabKE,    EDSingleX);
     TA2DataManager::LoadVariable("Pi0LabTheta",&fPi0LabTheta, EDSingleX);
     TA2DataManager::LoadVariable("Pi0LabPhi",  &fPi0LabPhi,   EDSingleX);
     TA2DataManager::LoadVariable("Pi0CMKE",    &fPi0CMKE,     EDSingleX);
     TA2DataManager::LoadVariable("Pi0CMTheta", &fPi0CMTheta,  EDSingleX);
     TA2DataManager::LoadVariable("Pi0CMPhi",   &fPi0CMPhi,    EDSingleX);

     TA2DataManager::LoadVariable("TaggerTime",  fTaggerTime,  EDMultiX);
     TA2DataManager::LoadVariable("Pi0TaggTime", fPi0TaggTime, EDMultiX);
     TA2DataManager::LoadVariable("TChanHit",    fTChanHit,    EIMultiX);

     TA2DataManager::LoadVariable("TChanP",      fTChanP,      EIMultiX);
     TA2DataManager::LoadVariable("TChanR",      fTChanR,      EIMultiX);
     TA2DataManager::LoadVariable("MMissP",      fMMissP,      EDMultiX);
     TA2DataManager::LoadVariable("MMissR",      fMMissR,      EDMultiX);
     TA2DataManager::LoadVariable("TGGP",        fTGGP,        EDMultiX);
     TA2DataManager::LoadVariable("TGGR",        fTGGR,        EDMultiX);
     TA2DataManager::LoadVariable("M2gTaggP",    fM2gTaggP,    EDMultiX);
     TA2DataManager::LoadVariable("M2gTaggR",    fM2gTaggR,    EDMultiX);

     TA2DataManager::LoadVariable("KELabP",      fKELabP,      EDMultiX);
     TA2DataManager::LoadVariable("KELabR",      fKELabR,      EDMultiX);
     TA2DataManager::LoadVariable("ThetaLabP",   fThetaLabP,   EDMultiX);
     TA2DataManager::LoadVariable("ThetaLabR",   fThetaLabR,   EDMultiX);
     TA2DataManager::LoadVariable("PhiLabP",     fPhiLabP,     EDMultiX);
     TA2DataManager::LoadVariable("PhiLabR",     fPhiLabR,     EDMultiX);

     TA2DataManager::LoadVariable("KECMP",       fKECMP,       EDMultiX);
     TA2DataManager::LoadVariable("KECMR",       fKECMR,       EDMultiX);
     TA2DataManager::LoadVariable("ThetaCMP",    fThetaCMP,    EDMultiX);
     TA2DataManager::LoadVariable("ThetaCMR",    fThetaCMR,    EDMultiX);
     TA2DataManager::LoadVariable("PhiCMP",      fPhiCMP,      EDMultiX);
     TA2DataManager::LoadVariable("PhiCMR",      fPhiCMR,      EDMultiX);

     TA2DataManager::LoadVariable("TChanCut1P",  fTChanCut1P,  EIMultiX);
     TA2DataManager::LoadVariable("TChanCut1R",  fTChanCut1R,  EIMultiX);
     TA2DataManager::LoadVariable("KECMCut1P",   fKECMCut1P,   EDMultiX);
     TA2DataManager::LoadVariable("KECMCut1R",   fKECMCut1R,   EDMultiX);
     TA2DataManager::LoadVariable("ThetaCMCut1P",fThetaCMCut1P,EDMultiX);
     TA2DataManager::LoadVariable("ThetaCMCut1R",fThetaCMCut1R,EDMultiX);
     TA2DataManager::LoadVariable("MMissCut1P",  fMMissCut1P,  EDMultiX);
     TA2DataManager::LoadVariable("MMissCut1R",  fMMissCut1R,  EDMultiX);

     TA2DataManager::LoadVariable("TChanCut2P",  fTChanCut2P,  EIMultiX);
     TA2DataManager::LoadVariable("TChanCut2R",  fTChanCut2R,  EIMultiX);
     TA2DataManager::LoadVariable("KECMCut2P",   fKECMCut2P,   EDMultiX);
     TA2DataManager::LoadVariable("KECMCut2R",   fKECMCut2R,   EDMultiX);
     TA2DataManager::LoadVariable("ThetaCMCut2P",fThetaCMCut2P,EDMultiX);
     TA2DataManager::LoadVariable("ThetaCMCut2R",fThetaCMCut2R,EDMultiX);
     TA2DataManager::LoadVariable("PhiCMCut2P",  fPhiCMCut2P,  EDMultiX);
     TA2DataManager::LoadVariable("PhiCMCut2R",  fPhiCMCut2R,  EDMultiX);

     return;
}

//-----------------------------------------------------------------------------
void TA2ThreshPi0::Reconstruct()
{
     // General reconstruction of reaction kinematics in Mainz tagged-photon
     // meson production experiments.
     // Use 4-momenta and PDG-index information from apparati to reconstruct
     // reaction kinematics. The PDG index (and 4-momentum) assigned by the
     // apparatus is not considered binding, e.g. in cases where n/gamma
     // discrimination by an apparatus is not possible, in which case it
     // defaults to kGamma. The method TA2ParticleID->SetMassP4( *p4, ipdg )
     // may be used to reset the rest-mass of an existing 4 momentum *p4 to that
     // corresponding to PDG index ipdg.
     // This one deals with pion and eta photoproduction on the nucleon.

     Int_t i;
     Int_t ntagg = fTAGG->GetNparticle();            // # particles in Tagger
     Int_t ncb = fCB->GetNparticle();                // # particles in CB
     Int_t ntaps;
     if ( fTAPS ) ntaps = fTAPS->GetNparticle();     // # particles in TAPS
     else ntaps = 0;

     // Temporary Tagger Stuff
     Double_t* TaggTime = fLADD->GetTimeOR();        // tagger time

     // zero particle counters
     fNphoton = 0;
     fNproton = 0;
     fNpiplus = 0;
     fNneutron = 0;
     fNrootino = 0;

     fNpi0 = 0;
     fNeta = 0;
     fNgprime = 0;

     fM2g = ENullHit;                         // zero 2-gamma inv. mass
     fM6g = ENullHit;                         // zero 6-gamma inv. mass

     fPi0TGG = ENullHit;
     fPi0LabKE = ENullHit;
     fPi0LabTheta = ENullHit;
     fPi0LabPhi = ENullHit;
     fPi0CMKE = ENullHit;
     fPi0CMTheta = ENullHit;
     fPi0CMPhi = ENullHit;

     fNparticle = ncb + ntaps;                // total number particles (hits)

     // Sort 4-momenta provided by apparati according to particle type

     // Tagger
//     for ( i = 0; i < ntag; i++) fPARTtaggphot[i] = fTAGGpart+i;

     // CB
     for ( i = 0; i < ncb; i++ ){
          switch( (fCBpart+i)->GetParticleID() ) {    // PDG code
               case kGamma:                           // photon
                   fPARTphoton[fNphoton] = fCBpart+i; // include in photon list
                   fNphoton++;
                   break;
               case kProton:                          // proton
                   fPARTproton[fNproton] = fCBpart+i; // include in proton list
                   fNproton++;
                   break;
               case kPiPlus:                          // pi+
                   fPARTpiplus[fNpiplus] = fCBpart+i; // include in piplus list
                   fNpiplus++;
                   break;
               default:                                // don't know
                   fPARTrootino[fNrootino] = fCBpart+i;// include in rootino lst
                   fNrootino++;                          
          }
     }

     // TAPS
     for ( i = 0; i < ntaps; i++ ) {
          switch( (fTAPSpart+i)->GetParticleID() ) {   // PDG code
               case kGamma:                            // photon
                  fPARTphoton[fNphoton] = fTAPSpart+i; // include in photon list
                  fNphoton++;
                  break;
               case kProton:                           // proton
                  fPARTproton[fNproton] = fTAPSpart+i; // include in proton list
                  fNproton++;
                  break;
               case kPiPlus:                           // pi+
                  fPARTpiplus[fNpiplus] = fTAPSpart+i; // include in piplus list
                  fNpiplus++;
                  break;
               default:                                // don't know
                  fPARTrootino[fNrootino] = fTAPSpart+i;// include in rootino 
                  fNrootino++;                          
          }
     }

     // Check if detected photons combine to give pi0 or eta
     TLorentzVector p4;
     switch( fNphoton ){
          case 1:
               // Just 1 photon....assume it is a gamma-prime
               fPARTgprime[fNgprime++] = fPARTphoton[0];
               break;
          case 2:
               // 2 photons detected, fast check if they make a pi0 or eta
               Sort2Photon();
               break;
          default:
               // More than 2 photons 
               SortNPhoton();
               // Check for 3-pi0 eta decay mode
               if( fNpi0 == 3 ){
                 p4 = (*fPARTpi0[0]).GetP4() + (*fPARTpi0[1]).GetP4()
                      + (*fPARTpi0[2]).GetP4();
                 fMassDpi0[0] = TMath::Abs( p4.M()-fParticleID->GetMassMeV(kEta));
                    if( fMassDpi0[0] < fMaxMDeta ){
                         (*fPARTeta[0]).GetP4() = p4;
                         fNeta = 1;
                         fNpi0 = 0;
                    }
               }
               break;
     }

     fNprompt = 0;
     fNrandom = 0;
     MarkEndBuffer();
     fTaggerTime[0] = EBufferEnd;
     fPi0TaggTime[0] = EBufferEnd;
     fTChanHit[0] = EBufferEnd;
     Int_t jtagg = 0;

     // Two photons forming a pi0
     if ( ( fNphoton == 2) && ( fNpi0 == 1)) {

          //Decay Photons and Opening Angle
          TA2Particle phot1 = *fPARTphoton[0];
          TA2Particle phot2 = *fPARTphoton[1];
          TLorentzVector p4a, p4b;
          p4a = phot1.GetP4();
          p4b = phot2.GetP4();
          fPi0TGG = p4a.Vect().Angle(p4b.Vect())*TMath::RadToDeg();

          // Gamma-gamma time difference
          fTime1 = phot1.GetTime();
          fTime2 = phot2.GetTime();
          fTimeDiff = fTime1 - fTime2;

/*
          This is still wonky....
          DLH 2010.03.01

          // Pi0 Parameters
          TA2Particle pi0;
          pi0 = *fPARTpi0[0];
          Double_t mPi0 = fParticleID->GetMassMeV( kPi0);

          // Pi0 lab parameters corrected for invariant mass
          Double_t p4mom, dx, dy, dz;
          dx = cos( pi0.GetPhi())*sin( pi0.GetTheta());
          dy = sin( pi0.GetPhi())*sin( pi0.GetTheta());
          dz = cos( pi0.GetTheta());
          if ( fMassCorr == 1) {

               pi0.SetE( pi0.GetE()*mPi0/pi0.GetM());
               p4mom = sqrt( Sqr( pi0.GetE()) - Sqr( mPi0));
               pi0.SetPx( p4mom*dx);
               pi0.SetPy( p4mom*dy);
               pi0.SetPz( p4mom*dz);
               pi0.SetM( mPi0);
          }
*/
          // Pi0 Parameters
          TA2Particle pi0;
          pi0 = *fPARTpi0[0];
          Double_t mPi0 = fParticleID->GetMassMeV( kPi0);

          TLorentzVector p4pi0;
          p4pi0 = pi0.GetP4();

          // Pi0 lab parameters corrected for invariant-mass resolution
          Double_t p4mom, dx, dy, dz;
          dx = cos( p4pi0.Phi())*sin( p4pi0.Theta());
          dy = sin( p4pi0.Phi())*sin( p4pi0.Theta());
          dz = cos( p4pi0.Theta());
          if ( fMassCorr == 1) {

               p4pi0.SetE( p4pi0.E()*mPi0/p4pi0.M());
               p4mom = sqrt( Sqr( p4pi0.E()) - Sqr( mPi0));
               p4pi0.SetPx( p4mom*dx);
               p4pi0.SetPy( p4mom*dy);
               p4pi0.SetPz( p4mom*dz);
          }

          fPi0LabKE = p4pi0.E() - p4pi0.M();
          fPi0LabTheta = p4pi0.Theta()*TMath::RadToDeg();
          fPi0LabPhi = p4pi0.Phi()*TMath::RadToDeg();

          // Loop over tagger hits.
          // If photon energy out of range ignore tagger hit.
          // Check if hit is in the prompt or random region and set index
          // accordingly.
          // Ignore tagger hit if time not in set prompt/random region.
     
          // Tagger Loop
          for ( i = 0; i < ntagg; i++ )
          {
               // Tagger channel hit
               Int_t chan;
               chan = (fLADD->GetHits())[i];
               fTChanHit[i] = chan;

               // Tagger Time
               fTaggerTime[i] = TaggTime[i];

               // Pi0-Tagger Time
               fPi0TaggTime[i] = TaggTime[i] - pi0.GetTime();

               //
               // RANDOM SUBTRACTION SECTION
               //

               // Prompt or Monte Carlo
               if  ( ( ( fPi0TaggTime[i] >= fP1) && ( fPi0TaggTime[i] <= fP2))
                  || ( gAR->GetProcessType() == EMCProcess)) jtagg = fNprompt++;

               // Random
               else if (((fPi0TaggTime[i] >= fRl1) && (fPi0TaggTime[i] <= fRl2))
                  || ( (fPi0TaggTime[i] >= fRh1) && (fPi0TaggTime[i] <= fRh2)))
                  jtagg = fMaxTagg + fNrandom++;

               // Other - skip the rest!
               else continue;

               // Tagger channel
               fTChan[jtagg] = chan;

               // 2-gamma invariant mass
               fM2gTagg[jtagg] = fM2g;

               // and 2-gamma opening angle
               fTGG[jtagg] = fPi0TGG;

               // Lab pi0 variables
               fKELab[jtagg] = fPi0LabKE;
               fThetaLab[jtagg] = fPi0LabTheta;
               fPhiLab[jtagg] = fPi0LabPhi;

               // Incoming and missing 4-mom
               TLorentzVector p4In, p4miss;
               p4In = fP4target[0] + fTAGGp4[i];
               p4miss = p4In - p4pi0;
               fMMiss[jtagg] = p4miss.M();

               // Boost to CM frame
               TLorentzVector p4pi0cm;
               TVector3 cmBoost;
               cmBoost = -p4In.BoostVector();
               p4pi0cm = p4pi0;
               p4pi0cm.Boost( cmBoost);         

               // CM pi0 variables
               fPi0CMKE = p4pi0cm.E() - mPi0;
               fPi0CMTheta = p4pi0cm.Theta()*TMath::RadToDeg();
               fPi0CMPhi = p4pi0cm.Phi()*TMath::RadToDeg();
               fKECM[jtagg] = fPi0CMKE;
               fThetaCM[jtagg] = fPi0CMTheta;
               fPhiCM[jtagg] = fPi0CMPhi;

               // Cut parameters
               Double_t eg, q_pi0, T_pi0, Tgg, TggCut;
               Double_t mProt = fParticleID->GetMassMeV( kProton);
               eg = fTAGGp4[i].E();
               q_pi0 = qp_thcm( eg, mProt, fPi0CMTheta, mPi0);
               T_pi0 = Energy( q_pi0, mPi0) - mPi0;
               Tgg = Tgg_Min( T_pi0, mPi0)*TMath::RadToDeg();
               TggCut = Tgg - TggCutOffset;

               // Opening angle cut
               if ( fPi0TGG >= TggCut) {

                    fTChanCut1[jtagg] = chan;
                    fKECMCut1[jtagg] = fPi0CMKE;
                    fThetaCMCut1[jtagg] = fPi0CMTheta;
                    fMMissCut1[jtagg] = p4miss.M();

                    Double_t Mmisslo, Mmisshi;
                    Mmisslo = Linear( 145, fMMWlo[0], 400, fMMWlo[1], eg);
                    Mmisshi = Linear( 145, fMMWhi[0], 400, fMMWhi[1], eg);

                    // Missing mass cut
                    if ( ( p4miss.M() >= Mmisslo) && ( p4miss.M() <= Mmisshi)) {

                         fTChanCut2[jtagg] = chan;
                         fKECMCut2[jtagg] = fPi0CMKE;
                         fThetaCMCut2[jtagg] = fPi0CMTheta;
                         fPhiCMCut2[jtagg] = fPi0CMPhi;
                    }
               }
          }
     }

     MarkEndBuffer();
     fTaggerTime[i] = EBufferEnd;
     fPi0TaggTime[i] = EBufferEnd;
     fTChanHit[i] = EBufferEnd;
}
