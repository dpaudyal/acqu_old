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
// This one deals with pion photoproduction on the nucleon.
//
// PDG codes of particles generlly observed MAMI-B
// kElectron 11,     kPositron -11
// kMuonMinus 13     kMuonPlus -13      kGamma 22
// kPi0 111          kPiPlus 211        kPiMinus -211       kEta 221
// kProton 2212      kNeutron 2112
// 

#include "TA2Compton.h"

// Valid Keywords for command-line setup of ThreshPi0
enum { EPhotoMassLimits = 1000, EPhotoPRLimits, EPhotoMissMassLimits,
     EPhotoOpenCut};
static const Map_t kPhotoKeys[] = {
     {"Mass-Limits:",                  EPhotoMassLimits},
     {"Prompt-Random-Windows:",        EPhotoPRLimits},
     {"Missing-Mass-Limits:",          EPhotoMissMassLimits},
     {"Opening-Angle-Cut:",            EPhotoOpenCut},
     {NULL,            -1}
};

ClassImp(TA2Compton)

//-----------------------------------------------------------------------------
TA2Compton::TA2Compton( const char* name, TA2Analysis* analysis )
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

     fProtKE = ENullHit;
     fProtTheta = ENullHit;
     fProtPhi = ENullHit;
     fProtTime = ENullHit;

     fPhotKE = ENullHit;
     fPhotTheta = ENullHit;
     fPhotPhi = ENullHit;
     fPhotTime = ENullHit;

     fNprompt = 0;
     fNrandom = 0;
     fMaxTagg = 0;
     fMax2gPerm = 0;

     fTaggerTime = NULL;
     fPhotTaggTime = NULL;
     fTChanHit = NULL;

     fTChanPhotP = NULL;
     fTChanPhotR = NULL;

     fPhotonMmissP = NULL;
     fPhotonMmissR = NULL;

     fPhotonKECMP = NULL;
     fPhotonKECMR = NULL;
     fPhotonThetaCMP = NULL;
     fPhotonThetaCMR = NULL;
     fPhotonPhiCMP = NULL;
     fPhotonPhiCMR = NULL;

     fTChanPhotProtP = NULL;
     fTChanPhotProtR = NULL;
     fPhotProtOAP = NULL;
     fPhotProtOAR = NULL;
     fPhotonMmissProtP = NULL;
     fPhotonMmissProtR = NULL;
     fPhotonEmissProtP = NULL;
     fPhotonEmissProtR = NULL;
     fPhotonThetaCMProtP = NULL;
     fPhotonThetaCMProtR = NULL;

     fPhotonMmissProtOAP = NULL;
     fPhotonMmissProtOAR = NULL;
     fTChanPhotProtOAP = NULL;
     fTChanPhotProtOAR = NULL;
     fPhotonThetaCMProtOAP = NULL;
     fPhotonThetaCMProtOAR = NULL;

     AddCmdList( kPhotoKeys );       // command-line recognition for SetConfig()
}


//-----------------------------------------------------------------------------
TA2Compton::~TA2Compton()
{
     // Free up allocated memory...after checking its allocated
     // detector and cuts lists
}

//-----------------------------------------------------------------------------
void TA2Compton::SetConfig(Char_t* line, Int_t key)
{
     // Any special command-line input for Crystal Ball apparatus

     switch (key){
          case EPhotoMassLimits:
               //  Invariant mass limits
               if( sscanf( line, "%lf%lf", &fMaxMDpi0, &fMaxMDeta ) != 2 ){
                  PrintError( line, "<ThreshPi0 meson invariant mass limits>");
                  return;
               }
               break;
          case EPhotoPRLimits:
               //  Prompt-Random Windows
               if( sscanf( line, "%lf%lf%lf%lf%lf%lf", &fP1, &fP2, &fRl1, &fRl2,
                              &fRh1, &fRh2) != 6 ){
                 PrintError( line, "<ThreshPi0 prompt-random windows>");
                 return;
               }
               break;
          case EPhotoMissMassLimits:
               //  Invariant mass limits
               if( sscanf( line, "%lf%lf%lf%lf", &fMMMinLo, &fMMMinHi, &fMMMaxLo,
                              &fMMMaxHi) != 4 ){
                 PrintError( line, "<ThreshPi0 missing mass limits>");
                 return;
               }
               break;
          case EPhotoOpenCut:
               //  Opening Angle Cut
               if( sscanf( line, "%lf", &fOACut) != 1 ){
                 PrintError( line, "<Opening Angle Cut>");
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
void TA2Compton::PostInit()
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

//     fTAGG = (TA2KensTagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
//     if(!fTAGG)PrintError("","<No Tagger class found in annalysis>",EErrFatal);
//     else{
//          fTAGGpart = fTAGG->GetParticles();
//     }
     fTAGG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
     if (!fTAGG)PrintError("","<No Tagger class found in annalysis>",EErrFatal);
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
     fPhotTaggTime = new Double_t[maxtagg];
     fTChanHit = new Int_t[maxtagg];

     fTChanPhot = new Int_t[maxtagg * 2];
     fTChanPhotP = fTChanPhot; fTChanPhotR = fTChanPhot + maxtagg;

     fPhotonMmiss = new Double_t[maxtagg * 2];
     fPhotonMmissP = fPhotonMmiss; fPhotonMmissR = fPhotonMmiss + maxtagg;

     fPhotonKECM = new Double_t[maxtagg * 2];
     fPhotonKECMP = fPhotonKECM; fPhotonKECMR = fPhotonKECM + maxtagg;
     fPhotonThetaCM = new Double_t[maxtagg * 2];
     fPhotonThetaCMP = fPhotonThetaCM; fPhotonThetaCMR=fPhotonThetaCM+maxtagg;
     fPhotonPhiCM = new Double_t[maxtagg * 2];
     fPhotonPhiCMP = fPhotonPhiCM; fPhotonPhiCMR = fPhotonPhiCM + maxtagg;

     fTChanPhotProt = new Int_t[maxtagg * 2];
     fTChanPhotProtP = fTChanPhotProt; fTChanPhotProtR=fTChanPhotProt+maxtagg;
     fPhotProtOA = new Double_t[maxtagg * 2];
     fPhotProtOAP = fPhotProtOA;
     fPhotProtOAR = fPhotProtOA + maxtagg;
     fPhotonMmissProt = new Double_t[maxtagg * 2];
     fPhotonMmissProtP = fPhotonMmissProt;
     fPhotonMmissProtR = fPhotonMmissProt + maxtagg;
     fPhotonEmissProt = new Double_t[maxtagg * 2];
     fPhotonEmissProtP = fPhotonEmissProt;
     fPhotonEmissProtR = fPhotonEmissProt + maxtagg;
     fPhotonThetaCMProt = new Double_t[maxtagg * 2];
     fPhotonThetaCMProtP = fPhotonThetaCMProt;
     fPhotonThetaCMProtR = fPhotonThetaCMProt + maxtagg;

     fPhotonMmissProtOA = new Double_t[maxtagg * 2];
     fPhotonMmissProtOAP = fPhotonMmissProtOA;
     fPhotonMmissProtOAR = fPhotonMmissProtOA + maxtagg;
     fTChanPhotProtOA = new Double_t[maxtagg * 2];
     fTChanPhotProtOAP = fTChanPhotProtOA;
     fTChanPhotProtOAR = fTChanPhotProtOA + maxtagg;
     fPhotonThetaCMProtOA = new Double_t[maxtagg * 2];
     fPhotonThetaCMProtOAP = fPhotonThetaCMProtOA;
     fPhotonThetaCMProtOAR = fPhotonThetaCMProtOA + maxtagg;

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
void TA2Compton::LoadVariable( )
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
     TA2DataManager::LoadVariable("Nphoton",        &fNphoton,       EISingleX);
     TA2DataManager::LoadVariable("Nproton",        &fNproton,       EISingleX);
     TA2DataManager::LoadVariable("Npiplus",        &fNpiplus,       EISingleX);
     TA2DataManager::LoadVariable("Nneutron",       &fNneutron,      EISingleX);
     TA2DataManager::LoadVariable("Nrootino",       &fNrootino,      EISingleX);
     TA2DataManager::LoadVariable("Ngprime",        &fNgprime,       EISingleX);
     TA2DataManager::LoadVariable("Npi0",           &fNpi0,          EISingleX);
     TA2DataManager::LoadVariable("Neta",           &fNeta,          EISingleX);
     TA2DataManager::LoadVariable("M2g",            &fM2g,           EDSingleX);
     TA2DataManager::LoadVariable("M6g",            &fM6g,           EDSingleX);

     TA2DataManager::LoadVariable("ProtKE",         &fProtKE,        EDSingleX);
     TA2DataManager::LoadVariable("ProtTheta",      &fProtTheta,     EDSingleX);
     TA2DataManager::LoadVariable("ProtPhi",        &fProtPhi,       EDSingleX);
     TA2DataManager::LoadVariable("ProtTime",       &fProtTime,      EDSingleX);

     TA2DataManager::LoadVariable("PhotKE",         &fPhotKE,        EDSingleX);
     TA2DataManager::LoadVariable("PhotTheta",      &fPhotTheta,     EDSingleX);
     TA2DataManager::LoadVariable("PhotPhi",        &fPhotPhi,       EDSingleX);
     TA2DataManager::LoadVariable("PhotTime",       &fPhotTime,      EDSingleX);

     TA2DataManager::LoadVariable("TaggerTime",      fTaggerTime,    EDMultiX);
     TA2DataManager::LoadVariable("PhotTaggTime",    fPhotTaggTime,  EDMultiX);
     TA2DataManager::LoadVariable("TChanHit",        fTChanHit,      EIMultiX);

     TA2DataManager::LoadVariable("TChanPhotP",      fTChanPhotP,    EIMultiX);
     TA2DataManager::LoadVariable("TChanPhotR",      fTChanPhotR,    EIMultiX);

     TA2DataManager::LoadVariable("PhotonMmissP",    fPhotonMmissP,  EDMultiX);
     TA2DataManager::LoadVariable("PhotonMmissR",    fPhotonMmissR,  EDMultiX);

     TA2DataManager::LoadVariable("PhotonKECMP",     fPhotonKECMP,   EDMultiX);
     TA2DataManager::LoadVariable("PhotonKECMR",     fPhotonKECMR,   EDMultiX);
     TA2DataManager::LoadVariable("PhotonThetaCMP",  fPhotonThetaCMP,EDMultiX);
     TA2DataManager::LoadVariable("PhotonThetaCMR",  fPhotonThetaCMR,EDMultiX);
     TA2DataManager::LoadVariable("PhotonPhiCMP",    fPhotonPhiCMP,  EDMultiX);
     TA2DataManager::LoadVariable("PhotonPhiCMR",    fPhotonPhiCMR,  EDMultiX);

     TA2DataManager::LoadVariable("TChanPhotProtP",  fTChanPhotProtP,EIMultiX);
     TA2DataManager::LoadVariable("TChanPhotProtR",  fTChanPhotProtR,EIMultiX);
     TA2DataManager::LoadVariable("PhotProtOAP",     fPhotProtOAP,   EDMultiX);
     TA2DataManager::LoadVariable("PhotProtOAR",     fPhotProtOAR,   EDMultiX);
     TA2DataManager::LoadVariable("PhotonMmissProtP",fPhotonMmissProtP,     
               EDMultiX);
     TA2DataManager::LoadVariable("PhotonMmissProtR",fPhotonMmissProtR,     
               EDMultiX);
     TA2DataManager::LoadVariable("PhotonEmissProtP",fPhotonEmissProtP,     
               EDMultiX);
     TA2DataManager::LoadVariable("PhotonEmissProtR",fPhotonEmissProtR,     
               EDMultiX);
     TA2DataManager::LoadVariable("PhotonThetaCMProtP",fPhotonThetaCMProtP,     
               EDMultiX);
     TA2DataManager::LoadVariable("PhotonThetaCMProtR",fPhotonThetaCMProtR,     
               EDMultiX);
     TA2DataManager::LoadVariable("PhotonMmissProtOAP",fPhotonMmissProtOAP,
               EDMultiX);
     TA2DataManager::LoadVariable("PhotonMmissProtOAR",fPhotonMmissProtOAR,
               EDMultiX);
     TA2DataManager::LoadVariable("TChanPhotProtOAP",  fTChanPhotProtOAP,
               EDMultiX);
     TA2DataManager::LoadVariable("TChanPhotProtOAR",  fTChanPhotProtOAR,
               EDMultiX);
     TA2DataManager::LoadVariable("PhotonThetaCMProtOAP",fPhotonThetaCMProtOAP,
               EDMultiX);
     TA2DataManager::LoadVariable("PhotonThetaCMProtOAR",fPhotonThetaCMProtOAR,
               EDMultiX);

     return;
}

//-----------------------------------------------------------------------------
void TA2Compton::Reconstruct()
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

     fProtKE = ENullHit;
     fProtTheta = ENullHit;
     fProtPhi = ENullHit;
     fProtTime = ENullHit;

     fPhotKE = ENullHit;
     fPhotTheta = ENullHit;
     fPhotPhi = ENullHit;
     fPhotTime = ENullHit;

     fNparticle = ncb + ntaps;                // total number particles (hits)

     // Sort 4-momenta provided by apparati according to particle type

     // Tagger
//     for ( i = 0; i < ntag; i++) fPARTtaggphot[i] = fTAGGpart+i;

     // CB
     for ( i = 0; i < ncb; i++ ){
          switch( (fCBpart+i)->GetParticleID() ) {     // PDG code
               case kGamma:                            // photon
                    fPARTphoton[fNphoton] = fCBpart+i; // include in photon list
                    fNphoton++;
                    break;
               case kProton:                           // proton
                    fPARTproton[fNproton] = fCBpart+i; // include in proton list
                    fNproton++;
                    break;
               case kPiPlus:                           // pi+
                    fPARTpiplus[fNpiplus] = fCBpart+i; // include in piplus list
                    fNpiplus++;
                    break;
               default:                                 // don't know
                    fPARTrootino[fNrootino] = fCBpart+i;// include in rootino 
                    fNrootino++;                          
          }
     }

     // TAPS
     for ( i = 0; i < ntaps; i++ ) {
        switch( (fTAPSpart+i)->GetParticleID() ) {    // PDG code
             case kGamma:                             // photon
                  fPARTphoton[fNphoton] = fTAPSpart+i;// include in photon list
                  fNphoton++;
                  break;
             case kProton:                            // proton
                  fPARTproton[fNproton] = fTAPSpart+i;// include in proton list
                  fNproton++;
                  break;
             case kPiPlus:                            // pi+
                  fPARTpiplus[fNpiplus] = fTAPSpart+i;// include in piplus list
                  fNpiplus++;
                  break;
             default:                                 // don't know
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
                  fMassDpi0[0]=TMath::Abs(p4.M()-fParticleID->GetMassMeV(kEta));
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
     fPhotTaggTime[0] = EBufferEnd;
     fTChanHit[0] = EBufferEnd;
     Int_t jtagg = 0;

     // 1 Proton ignoring everything else
     if ( fNproton == 1)
     {
          TA2Particle proton = *fPARTproton[0];
          fProtKE = proton.GetT();
          fProtTheta = proton.GetThetaDg();
          fProtPhi = proton.GetPhiDg();
          fProtTime = proton.GetTime();
     }

     // 1 Photon ignoring everything else
     if ( fNphoton == 1)
     {
       TA2Particle photon = *fPARTphoton[0];
       fPhotKE = photon.GetT();
       fPhotTheta = photon.GetThetaDg();
       fPhotPhi = photon.GetPhiDg();
       fPhotTime = photon.GetTime();

       // Tagger Loop
       for ( i = 0; i < ntagg; i++ )
       {
            // Tagger Channel Hit
            Int_t chan;
            chan = (fLADD->GetHits())[i];
            fTChanHit[i] = chan;

            // Tagger Time
            fTaggerTime[i] = TaggTime[i];

            // Photon-Tagger Time
            fPhotTaggTime[i] = TaggTime[i] - fPhotTime;

            //
            // RANDOM SUBTRACTION SECTION
            //

            // Prompt or Monte Carlo
            if  ( ( ( fPhotTaggTime[i] >= fP1) && ( fPhotTaggTime[i] <= fP2))
               || ( gAR->GetProcessType() == EMCProcess)) jtagg = fNprompt++;

            // Random
            else if ( ( (fPhotTaggTime[i]>=fRl1) && (fPhotTaggTime[i]<=fRl2))
               || ( ( fPhotTaggTime[i]>=fRh1) && ( fPhotTaggTime[i]<=fRh2)))
                 jtagg = fMaxTagg + fNrandom++;

            // Other - skip the rest!
            else continue;

            // Tagger channel
            fTChanPhot[jtagg] = chan;

            // Total incoming and missing 4-mom
            TLorentzVector p4In, p4miss;
            p4In = fP4target[0] + fTAGGp4[i];
            p4miss = p4In - photon.GetP4();
            fPhotonMmiss[jtagg] = p4miss.M();

            // Boost to CM frame
            TLorentzVector p4photcm;
            TVector3 cmBoost;
            cmBoost = -p4In.BoostVector();
            p4photcm = photon.GetP4();
            p4photcm.Boost( cmBoost);                           

            // Photon CM variables
            fPhotonKECM[jtagg] = p4photcm.E();
            fPhotonThetaCM[jtagg] = p4photcm.Theta()*TMath::RadToDeg();
            fPhotonPhiCM[jtagg] = p4photcm.Phi()*TMath::RadToDeg();

            // In coincidence with a proton
            if ( fNproton == 1) {

              TA2Particle proton = *fPARTproton[0];
              TLorentzVector prot = proton.GetP4();

              // Angle between detected proton and direction of missing
              // momentum vector.
              Double_t angle=p4miss.Vect().Angle(prot.Vect())*TMath::RadToDeg();

              fTChanPhotProt[jtagg] = chan;
              fPhotProtOA[jtagg] = angle;
              fPhotonMmissProt[jtagg] = p4miss.M();
              fPhotonEmissProt[jtagg] = p4miss.E();
              fPhotonThetaCMProt[jtagg] = p4photcm.Theta()*TMath::RadToDeg();

              if ( angle <= fOACut) {

                 fPhotonMmissProtOA[jtagg] = p4miss.M();
                 fTChanPhotProtOA[jtagg] = chan;
                 fPhotonThetaCMProtOA[jtagg]=p4photcm.Theta()*TMath::RadToDeg();
              }
           }
       }
     }

     MarkEndBuffer();
     fTaggerTime[i] = EBufferEnd;
     fPhotTaggTime[i] = EBufferEnd;
     fTChanHit[i] = EBufferEnd;
}
