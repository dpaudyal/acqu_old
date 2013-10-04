//--Author	JRM Annand   18th Feb 2004   Use def physics
//--Rev
//--Rev         JRM Annand   13th May 2004   Start general physics suite
//--Rev         JRM Annand   28th Apr 2004   General photo-meson methods
//--Rev         JRM Annand   13th Jul 2005   SortNPhoton bugs
//--Rev         JRM Annand   25th Jul 2005   ED bug fix fP4tot
//--Rev         JRM Annand   22nd Sep 2006   pi+n analysis add lin pol
//--Rev         JRM Annand   14th Dec 2006   fM2gCBTAPS for online diagnostic
//--Update      JRM Annand   14th Mar 2007   fEmProton for D(g,np) calib.
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2PhotoPhysics
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
// This one deals with pion and eta photoproduction on the nucleon.
// PDG codes of particles generlly observed MAMI-B
// kElectron 11,     kPositron -11
// kMuonMinus 13     kMuonPlus -13      kGamma 22
// kPi0 111          kPiPlus 211        kPiMinus -211       kEta 221
// kProton 2212      kNeutron 2112
// 

#include "TA2PhotoPhysics.h"
#include "TAcquRoot.h"
#include "TA2Analysis.h"
#include "TA2Calorimeter.h"
#include "TA2Tagger.h"
#include "TAcquFile.h"

// Valid Keywords for command-line setup of PhotoPhysics
enum { EPhotoMassLimits = 1000, EThreshPRLimits, EPhotoKinLimits };
static const Map_t kPhotoKeys[] = {
  {"Mass-Limits:",         EPhotoMassLimits},
	{"Prompt-Random-Windows:",		EThreshPRLimits},
  {"Kinematic-Limits:",    EPhotoKinLimits},
  {NULL,            -1}
};

ClassImp(TA2PhotoPhysics)

//-----------------------------------------------------------------------------
TA2PhotoPhysics::TA2PhotoPhysics( const char* name, TA2Analysis* analysis )
  :TA2Physics( name, analysis ) {
  // Initialise PhotoPhysics variables here
  // Default null pointers, zeroed variables

  fTAGG = NULL;
//   fLP = NULL;
  fCB = fTAPS = NULL;
  fTAGGp4 = fCBp4 = fTAPSp4 = NULL;
  fTAGGpdg = fCBpdg = fTAPSpdg = NULL;
  fP4photon = fP4proton = fP4piplus = fP4neutron = fP4pi0 = fP4gprime
    = fP4eta = fP4rootino = NULL;
  fNphoton = fNproton = fNpiplus = fNneutron = fNpi0 = fNgprime
    = fNeta = fNrootino = 0;
  fMassDpi0 = fMassDeta = NULL;
  fMassIJ = fMassIpi0 = fMassIeta = NULL;
  fIsMesonIndex = NULL;
  fMaxMDpi0 = fMaxMDeta = 0.0;
  fEmP = fEmR = NULL;
  fEmPi0P = fEmPi0R = fEmPi0pP = fEmPi0pR = fEmPi0gP = fEmPi0gR = 
    fEmPi0gpP = fEmPi0gpR = NULL;
  fPmPi0pP = fPcmPi0pR = NULL;
  fPcmPi0pP = fPcmPi0pR = NULL;
  fPcmPi0gpP = fPcmPi0gpR = NULL;
  fPmPi0gpP = fPcmPi0gpR = NULL;
  fEgPi0gpP = fEgPi0gpR = NULL;
  fWEgPi0gp = 0;
  fEmEtaP = fEmEtaR = fEmEtapP = fEmEtapR = NULL;
  fPcmEtapP = fPcmEtapR = NULL;
  fEmPiplusP = fEmPiplusR = NULL;
  fEmProtonP = fEmProtonR = NULL;
  fAopenP = fAopenR = NULL;
  fTheta = fThetaP = fThetaR = NULL;
  fPhi = fPhiParaP = fPhiParaR = fPhiPerpP = fPhiPerpR = NULL;
  fEmLow = fEgLow = fALow = 0.0;
  fEmHigh = 1.0E6;
  fEgHigh = 1.0E6;
  fAHigh = 180.0;
  fAngle1P = fAngle1R = fAngle2P = fAngle2R = NULL;
  fM2g = ENullHit;
  fNprompt = fNrandom = fMaxTagg = fMax2gPerm = 0;
  fNparaP = fNparaR = fNperpP = fNperpR = 0;

	fTaggerTime = NULL;
	fPi0TaggTime = NULL;

  AddCmdList( kPhotoKeys );       // command-line recognition for SetConfig()
}


//-----------------------------------------------------------------------------
TA2PhotoPhysics::~TA2PhotoPhysics()
{
  // Free up allocated memory...after checking its allocated
  // detector and cuts lists
}

//-----------------------------------------------------------------------------
void TA2PhotoPhysics::SetConfig(Char_t* line, Int_t key)
{
  // Any special command-line input for Crystal Ball apparatus

  switch (key){
  case EPhotoMassLimits:
    //  Invariant mass limits
    if( sscanf(line,"%lf%lf",&fMaxMDpi0,&fMaxMDeta ) != 2 ){
      PrintError(line,"<PhotoPysics meson invariant mass limits>");
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
  case EPhotoKinLimits:
    //  Kinematics limits
    if( sscanf(line,"%lf%lf%lf%lf%lf%lf",
	       &fEgLow,&fEgHigh,&fEmLow,&fEmHigh,&fALow,&fAHigh) != 6 ){
      PrintError(line,"<PhotoPysics kinematic limits>");
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
void TA2PhotoPhysics::PostInit()
{
  // Initialise arrays to contain 4 momenta and plotable scaler variables
  // Missing mass, missing energy, cm momentum, energies, angles
  // Initialisation will abort if CB or Tagger not initialised
  // TAPS is optional
  Int_t i;

  fCB = (TA2Calorimeter*)((TA2Analysis*)fParent)->GetChild("CB");
  if( !fCB ) PrintError("","<No CB class found in annalysis>",EErrFatal);
  else{
    fCBp4 = fCB->GetP4();
    fCBpdg = fCB->GetPDG_ID();
  }
  fTAGG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
  if(!fTAGG) PrintError("","<No Tagger class found in annalysis>",EErrFatal);
  else{
    fTAGGp4 = fTAGG->GetP4();
    fTAGGpdg = fTAGG->GetPDG_ID();
  }
	fLADD = (TA2Ladder*)((TA2Analysis*)fParent)->GetGrandChild( "FPD");
	if ( !fLADD)
		PrintError( "", "<No Ladder class found in analysis>", EErrFatal);
  fTAPS = (TA2Calorimeter*)((TA2Analysis*)fParent)->GetChild("TAPS");
  if( !fTAPS ) PrintError("Warning!!!","<No TAPS class found in annalysis>");
  else{
    fTAPSp4 = fTAPS->GetP4();
    fTAPSpdg = fTAPS->GetPDG_ID();
  }
//   fLP = (TA2LinearPol*)((TA2Analysis*)fParent)->GetChild("LinPol");
//  if( !fLP ) PrintError("Warning!","<No Linear Pol class found in analysis>");

  // containers for sorted 4-momentum pointers
  Int_t maxparticle = fCB->GetMaxParticle();
  TLorentzVector* p4;
  if( fTAPS ) maxparticle += fTAPS->GetMaxParticle();
  fP4photon =  new TLorentzVector*[maxparticle];   
  fP4proton =  new TLorentzVector*[maxparticle];
  fP4piplus =  new TLorentzVector*[maxparticle];
  fP4neutron = new TLorentzVector*[maxparticle];
  fP4pi0 =     new TLorentzVector*[maxparticle];
  p4 = new TLorentzVector[maxparticle];           // space for pi0 4-mom
  for(i=0; i<maxparticle; i++) fP4pi0[i] = p4 + i;
  fP4gprime =  new TLorentzVector*[maxparticle];
  fP4eta =     new TLorentzVector*[maxparticle];
  p4 = new TLorentzVector[maxparticle];           // space for eta 4-mom
  for(i=0; i<maxparticle; i++) fP4eta[i] = p4 + i;
  fP4rootino = new TLorentzVector*[maxparticle];
  p4 = new TLorentzVector[maxparticle];           // space for unknown 4-mom
  for(i=0; i<maxparticle; i++) fP4rootino[i] = p4 + i;

  // Arrays for hold plotable variables...missing mass etc.
  Int_t maxtagg = fTAGG->GetMaxParticle() + 1;
  fMaxTagg = maxtagg;
  fEm = new Double_t[maxtagg * 2];
  fEmP = fEm; fEmR = fEm + maxtagg;
  fEmPi0 = new Double_t[maxtagg * 2];
  fEmPi0P = fEmPi0; fEmPi0R = fEmPi0 + maxtagg;
  fEmPi0p = new Double_t[maxtagg * 2];
  fEmPi0pP = fEmPi0p; fEmPi0pR = fEmPi0p + maxtagg;
  fPmPi0p = new Double_t[maxtagg * 2];
  fPmPi0pP = fPmPi0p; fPmPi0pR = fPmPi0p + maxtagg;
  fPcmPi0p = new Double_t[maxtagg * 2];
  fPcmPi0pP = fPcmPi0p; fPcmPi0pR = fPcmPi0p + maxtagg;
  fEmPi0g = new Double_t[maxtagg * 2];
  fEmPi0gP = fEmPi0g; fEmPi0gR = fEmPi0g + maxtagg;
  fEmPi0gp = new Double_t[maxtagg * 2];
  fEmPi0gpP = fEmPi0gp; fEmPi0gpR = fEmPi0gp + maxtagg;
  fPmPi0gp = new Double_t[maxtagg * 2];
  fPmPi0gpP = fPmPi0gp; fPmPi0gpR = fPmPi0gp + maxtagg;
  fPcmPi0gp = new Double_t[maxtagg * 2];
  fPcmPi0gpP = fPcmPi0gp; fPcmPi0gpR = fPcmPi0gp + maxtagg;
  fEgPi0gp = new Double_t[maxtagg * 2];
  fEgPi0gpP = fEgPi0gp; fEgPi0gpR = fEgPi0gp + maxtagg;
  fEmEta = new Double_t[maxtagg * 2];
  fEmEtaP = fEmEta; fEmEtaR = fEmEta + maxtagg;
  fEmEtap = new Double_t[maxtagg * 2];
  fEmEtapP = fEmEtap; fEmEtapR = fEmEtap + maxtagg;
  fPcmEtap = new Double_t[maxtagg * 2];
  fPcmEtapP = fPcmEtap; fPcmEtapR = fPcmEtap + maxtagg;
  fEmPiplus = new Double_t[maxtagg * 2];
  fEmPiplusP = fEmPiplus; fEmPiplusR = fEmPiplus + maxtagg;
  fEmProton = new Double_t[maxtagg * 2];
  fEmProtonP = fEmProton; fEmProtonR = fEmProton + maxtagg;
  fAopen = new Double_t[maxtagg * 2];
  fAopenP = fAopen; fAopenR = fAopen + maxtagg;
  fTheta = new Double_t[maxtagg * 2];
  fThetaP = fTheta; fThetaR = fTheta + maxtagg;
  fPhi = new Double_t[maxtagg * 4];
  fPhiParaP = fPhi; fPhiParaR = fPhi + maxtagg;
  fPhiPerpP = fPhi + 2*maxtagg; fPhiPerpR = fPhi + 3*maxtagg;
  fAngle1 = new Double_t[maxtagg*4];
  fAngle1P = fAngle1; fAngle1R = fAngle1P + maxtagg;
  fAngle2 = fAngle2P = fAngle1R + maxtagg; fAngle2R = fAngle2P + maxtagg;

	fTaggerTime = new Double_t[maxtagg];
	fPi0TaggTime = new Double_t[maxtagg];

  // Arrays used to combine photons to mesons
  Int_t maxperm = 0;
  for(i=1; i<= maxparticle; i++) maxperm += i;
  fMassDpi0 = new Double_t[maxperm];       // space for pi0 diffs
  fMassDeta = new Double_t[maxperm];       // space for eta diffs
  fMassIJ = new Int_t[maxperm];
  fMassIpi0 = new Int_t[maxperm];
  fMassIeta = new Int_t[maxperm];
  fMax2gPerm = maxperm;
  fIsMesonIndex = new Bool_t[maxparticle];

  // Default physics initialisation
  TA2Physics::PostInit();
}

//-----------------------------------------------------------------------------
void TA2PhotoPhysics::LoadVariable( )
{
  // Input name - variable pointer associations for any subsequent
  // cut or histogram setup
  // LoadVariable( "name", pointer-to-variable, type-spec );
  // NB scaler variable pointers need the preceeding &
  //    array variable pointers do not.
  // type-spec ED prefix for a Double_t variable
  //           EI prefix for an Int_t variable
  // type-spec SingleX for a single-valued variable
  //           MultiX  for a multi-valued variable

  TA2Physics::LoadVariable();
  TA2DataManager::LoadVariable("Nphoton",      &fNphoton,       EISingleX);
  TA2DataManager::LoadVariable("Nproton",      &fNproton,       EISingleX);
  TA2DataManager::LoadVariable("Npiplus",      &fNpiplus,       EISingleX);
  TA2DataManager::LoadVariable("Nneutron",     &fNneutron,      EISingleX);
  TA2DataManager::LoadVariable("Npi0",         &fNpi0,          EISingleX);
  TA2DataManager::LoadVariable("Ngprime",      &fNgprime,       EISingleX);
  TA2DataManager::LoadVariable("Neta",         &fNeta,          EISingleX);
  TA2DataManager::LoadVariable("Nrootino",     &fNrootino,      EISingleX);
  TA2DataManager::LoadVariable("M2g",          &fM2g,           EDSingleX);
  TA2DataManager::LoadVariable("M6g",          &fM6g,           EDSingleX);
  TA2DataManager::LoadVariable("M2gCBTAPS",    &fM2gCBTAPS,     EDSingleX);
  TA2DataManager::LoadVariable("EmPi0P",       fEmPi0P,         EDMultiX);
  TA2DataManager::LoadVariable("EmPi0R",       fEmPi0R,         EDMultiX);
	TA2DataManager::LoadVariable("TaggerTime",	fTaggerTime,		EDMultiX);
	TA2DataManager::LoadVariable("Pi0TaggTime",	fPi0TaggTime,		EDMultiX);
  return;
}

//-----------------------------------------------------------------------------
void TA2PhotoPhysics::Reconstruct()
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

  Int_t ntagg = fTAGG->GetNparticle();          // # tagger hits
  Int_t ncb = fCB->GetNparticle();              // # particles in CB
  Int_t ntaps, nptaps;
  if( fTAPS ) ntaps = fTAPS->GetNparticle();    // # particles in TAPS
  else ntaps = 0;
  nptaps = 0;

	Double_t* TaggTime = fLADD->GetTimeOR();		// tagger time
	Double_t fCBMeanTime = ((TA2CrystalBall*)fCB)->GetCBMeanTime();

  fNphoton = fNproton = fNpiplus = fNneutron = fNpi0 = fNgprime =
    fNeta = fNrootino = 0;                 // zero particle counters
  fM2g = ENullHit;                         // zero 2-gamma inv. mass
  fM2gCBTAPS = ENullHit;                   // zero 1-hit each in CB/TAPS 2-gamma inv. mass
  fM6g = ENullHit;

  fNparticle = ncb + ntaps;                // total number particles (hits)
  fP4tot.SetXYZT(0.0,0.0,0.0,0.0);         // zero total out 4-momentum
  //
  Int_t i,j,k;
  // Sort 4-momenta provided by apparati according to particle typw
  for( i=0; i<ncb; i++ ){                  // loop over CB hits
    fP4tot += fCBp4[i];
    switch( fCBpdg[i] ){                   // PDG code
    case kGamma:                           // photon
      fP4photon[fNphoton] = fCBp4+i;       // incl 4-momentum in photon list
      fNphoton++;
      break;
    case kProton:                          // proton
      fP4proton[fNproton] = fCBp4+i;       // incl 4-momentum in proton list
      fNproton++;
      break;
    case kPiPlus:                          // pi+
      fP4piplus[fNpiplus] = fCBp4+i;       // incl 4-momentum in pi+ list
      fNpiplus++;
      break;
    default:                               // don't know
      fP4rootino[fNrootino] = fCBp4+i;     // incl 4 momentum in unknown list
      fNrootino++;                          
    }
  }
  for( i=0; i<ntaps; i++ ){                // loop over TAPS hits
    fP4tot += fTAPSp4[i];
    switch( fTAPSpdg[i] ){                 // PDG code
    case kGamma:                           // photon
      fP4photon[fNphoton] = fTAPSp4+i;     // incl 4-momentum in photon list
      fNphoton++;
      break;
    case kProton:                          // proton
      fP4proton[fNproton] = fTAPSp4+i;     // incl 4-momentum in proton list
      fNproton++;
      nptaps++;
      break;
    case kPiPlus:                          // pi+
      fP4piplus[fNpiplus] = fTAPSp4+i;     // incl 4-momentum in proton list
      fNpiplus++;
      break;
    default:                               // don't know
      fP4rootino[fNrootino] = fTAPSp4+i;   // incl 4 momentum in unknown list
      fNrootino++;
    }
  }
  // Check if detected photons combine to give pi0 or eta
  TLorentzVector p4;
  switch( fNphoton ){
  case 1:
    // Just 1 photon....assume its a gamma-prime
    fP4gprime[ fNgprime++ ] = fP4photon[0];
    break;
  case 2:
    // 2 photons detected,  fast check if they make a pi0 or eta
    Sort2Photon();
    if( (ncb == 1) && (ntaps == 1) ) fM2gCBTAPS = fM2g;
    break;
  default:
    // > 2 photons 
    SortNPhoton( );
    // Check for 3-pi0 eta decay mode
    if( fNpi0 == 3 ){
      p4 = *fP4pi0[0] + *fP4pi0[1] + *fP4pi0[2];
      fMassDpi0[0] = TMath::Abs( p4.M() - fParticleID->GetMassMeV(kEta) );
      if( fMassDpi0[0] < fMaxMDeta ){
	*fP4eta[0] = p4;
	fNeta = 1;
	fNpi0 = 0;
      }
    }
    break;
  }

  fNprompt = fNrandom = fNparaP = fNparaR = fNperpP = fNperpR = 0;
  MarkEndBuffer();
  Double_t eGamma;
  Int_t jtagg = 0;
  TLorentzVector p4In, pi0;

  // Loop over tagger hits
  // If photon energy out of range ignore tagger hit
  // Check if hit is in the prompt or random region and set index accordingly
  // Ignore tagger hit if time not in set prompt/random region

	// Single pi0 production
	if (( fNphoton == 2) && ( fNpi0 == 1)) {

		// This section does the pi0 correction to the energy/momentum
		// if the missing mass is less than the physical mpi0.
		TLorentzVector p4pi0;
		Double_t p4mom, dx, dy, dz;
		Double_t mPi0 = fParticleID->GetMassMeV( kPi0);
		p4pi0 = *fP4pi0[0];
		dx = cos( p4pi0.Phi())*sin( p4pi0.Theta());
		dy = sin( p4pi0.Phi())*sin( p4pi0.Theta());
		dz = cos( p4pi0.Theta());
		p4pi0.SetE( p4pi0.E()*mPi0/p4pi0.M());
		p4mom = sqrt( pow( p4pi0.E(), 2) - pow( mPi0, 2));
		p4pi0.SetPx( p4mom*dx);
		p4pi0.SetPy( p4mom*dy);
		p4pi0.SetPz( p4mom*dz);

		for( i = 0; i < ntagg; i++)
		{

			fTaggerTime[i] = TaggTime[i];
			fPi0TaggTime[i] = TaggTime[i] + fCBMeanTime;

  			eGamma = fTAGGp4[i].E();

			if( (eGamma <= fEgLow) || (eGamma >= fEgHigh ) ) continue;

			// Prompt or Monte Carlo
			if  ( ( ( fPi0TaggTime[i] >= fP1) && ( fPi0TaggTime[i] <= fP2))
					|| ( gAR->GetProcessType() == EMCProcess)) jtagg = fNprompt++;

			// Random
			else if ( ( ( fPi0TaggTime[i] >= fRl1) && ( fPi0TaggTime[i] <= fRl2))
					|| ( ( fPi0TaggTime[i] >= fRh1) && ( fPi0TaggTime[i] <= fRh2)))
				jtagg = fMaxTagg + fNrandom++;

			// Other - skip the rest!
			else continue;

			p4In = fP4target[0] + fTAGGp4[i];
			p4 = p4In - p4pi0;
			fEmPi0[jtagg] = p4.M();
		}
	}
	fTaggerTime[i] = EBufferEnd;
	fPi0TaggTime[i] = EBufferEnd;
	MarkEndBuffer();
}
