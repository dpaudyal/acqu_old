//--Author	JRM Annand   18th Feb 2004   Use def physics
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
//--Update		 AT Laffoley 11th June 2009	Converting the DLH way to JRMA way
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2Efficiency
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

#include "TA2Efficiency.h"

// Valid Keywords for command-line setup of ThreshPi0
enum { EPhotoMassLimits = 1000, EPhotoPRLimits, EPhotoMissMassLimits,
	EPhotoOpenCut};
static const Map_t kPhotoKeys[] = {
	{"Mass-Limits:",					EPhotoMassLimits},
	{"Prompt-Random-Windows:",		EPhotoPRLimits},
	{"Missing-Mass-Limits:",		EPhotoMissMassLimits},
	{NULL,            -1}
};

ClassImp(TA2Efficiency)

//-----------------------------------------------------------------------------
TA2Efficiency::TA2Efficiency( const char* name, TA2Analysis* analysis )
	:TA2Physics( name, analysis ) {
	// Initialise ThreshPi0 variables here
	// Default null pointers, zeroed variables

	fTAGG = NULL;
	fCB = NULL;
	fTAPS = NULL;

	fTAGGp4 = NULL;

//	fPARTtaggphot = NULL;
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

	fPhotKE = ENullHit;
	fPhotTheta = ENullHit;
	fPhotPhi = ENullHit;

	fProtKE = ENullHit;
	fProtTheta = ENullHit;
	fProtPhi = ENullHit;

	fPi0KE = ENullHit;
	fPi0Theta = ENullHit;
	fPi0Phi = ENullHit;

	fMaxTagg = 0;
	fMax2gPerm = 0;

	AddCmdList( kPhotoKeys );       // command-line recognition for SetConfig()
}


//-----------------------------------------------------------------------------
TA2Efficiency::~TA2Efficiency()
{
	// Free up allocated memory...after checking its allocated
	// detector and cuts lists
}

//-----------------------------------------------------------------------------
void TA2Efficiency::SetConfig(Char_t* line, Int_t key)
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
		default:
			// default main apparatus SetConfig()
			TA2Physics::SetConfig( line, key );
			break;
	}
}

//---------------------------------------------------------------------------
void TA2Efficiency::PostInit()
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

//	fTAGG = (TA2KensTagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
//	if(!fTAGG) PrintError("","<No Tagger class found in annalysis>", EErrFatal);
//	else{
//		fTAGGpart = fTAGG->GetParticles();
//	}
	fTAGG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
	if ( !fTAGG) PrintError("","<No Tagger class found in annalysis>",EErrFatal);
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
void TA2Efficiency::LoadVariable( )
{
	// Input name - variable pointer associations for any subsequent
	// cut or histogram setup
	// LoadVariable( "name", pointer-to-variable, type-spec );
	// NB scaler variable pointers need the preceeding &
	//	array variable pointers do not.
	// type-spec ED prefix for a Double_t variable
	//				 EI prefix for an Int_t variable
	// type-spec SingleX for a single-valued variable
	//				 MultiX  for a multi-valued variable

	TA2Physics::LoadVariable();
	TA2DataManager::LoadVariable("Nphoton",		&fNphoton,			EISingleX);
	TA2DataManager::LoadVariable("Nproton",		&fNproton,			EISingleX);
	TA2DataManager::LoadVariable("Npiplus",		&fNpiplus,			EISingleX);
	TA2DataManager::LoadVariable("Nneutron",		&fNneutron,			EISingleX);
	TA2DataManager::LoadVariable("Nrootino",		&fNrootino,			EISingleX);
	TA2DataManager::LoadVariable("Ngprime",		&fNgprime,			EISingleX);
	TA2DataManager::LoadVariable("Npi0",			&fNpi0,				EISingleX);
	TA2DataManager::LoadVariable("Neta",			&fNeta,				EISingleX);
	TA2DataManager::LoadVariable("M2g",				&fM2g,				EDSingleX);
	TA2DataManager::LoadVariable("M6g",				&fM6g,				EDSingleX);

	TA2DataManager::LoadVariable("PhotKE",			&fPhotKE,			EDSingleX);
	TA2DataManager::LoadVariable("PhotTheta", 	&fPhotTheta,		EDSingleX);
	TA2DataManager::LoadVariable("PhotPhi",		&fPhotPhi,			EDSingleX);

	TA2DataManager::LoadVariable("ProtKE",			&fProtKE,			EDSingleX);
	TA2DataManager::LoadVariable("ProtTheta", 	&fProtTheta,		EDSingleX);
	TA2DataManager::LoadVariable("ProtPhi",		&fProtPhi,			EDSingleX);

	TA2DataManager::LoadVariable("Pi0KE",			&fPi0KE,				EDSingleX);
	TA2DataManager::LoadVariable("Pi0Theta", 		&fPi0Theta,			EDSingleX);
	TA2DataManager::LoadVariable("Pi0Phi",			&fPi0Phi,			EDSingleX);

	return;
}

//-----------------------------------------------------------------------------
void TA2Efficiency::Reconstruct()
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
	Int_t ntagg = fTAGG->GetNparticle();				// # particles in Tagger
	Int_t ncb = fCB->GetNparticle();					// # particles in CB
	Int_t ntaps;
	if ( fTAPS ) ntaps = fTAPS->GetNparticle();	// # particles in TAPS
	else ntaps = 0;

	// Temporary Tagger Stuff
	Double_t* TaggTime = fLADD->GetTimeOR();		// tagger time

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

	fPhotKE = ENullHit;
	fPhotTheta = ENullHit;
	fPhotPhi = ENullHit;

	fProtKE = ENullHit;
	fProtTheta = ENullHit;
	fProtPhi = ENullHit;

	fPi0KE = ENullHit;
	fPi0Theta = ENullHit;
	fPi0Phi = ENullHit;

	fNparticle = ncb + ntaps;                // total number particles (hits)

	// Sort 4-momenta provided by apparati according to particle type

	// Tagger
//	for ( i = 0; i < ntag; i++) fPARTtaggphot[i] = fTAGGpart+i;

	// CB
	for ( i = 0; i < ncb; i++ ){
		switch( (fCBpart+i)->GetParticleID() ) {		// PDG code
			case kGamma:										// photon
				fPARTphoton[fNphoton] = fCBpart+i;	// include in photon list
				fNphoton++;
				break;
			case kProton:										// proton
				fPARTproton[fNproton] = fCBpart+i;		// include in proton list
				fNproton++;
				break;
			case kPiPlus:										// pi+
				fPARTpiplus[fNpiplus] = fCBpart+i;		// include in piplus list
				fNpiplus++;
				break;
			default:												// don't know
				fPARTrootino[fNrootino] = fCBpart+i;	// include in rootino list
				fNrootino++;                          
		}
	}

	// TAPS
	for ( i = 0; i < ntaps; i++ ) {
//		std::cout << (fTAPSpart+i)->GetParticleID() << std::endl;
		switch( (fTAPSpart+i)->GetParticleID() ) {	// PDG code
			case kGamma:										// photon
				fPARTphoton[fNphoton] = fTAPSpart+i;	// include in photon list
				fNphoton++;
				break;
			case kProton:										// proton
				fPARTproton[fNproton] = fTAPSpart+i;	// include in proton list
				fNproton++;
				break;
			case kPiPlus:										// pi+
				fPARTpiplus[fNpiplus] = fTAPSpart+i;	// include in piplus list
				fNpiplus++;
				break;
			default:												// don't know
				fPARTrootino[fNrootino] = fTAPSpart+i;	// include in rootino list
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
/*			if( fNpi0 == 3 ){
				p4 = (*fPARTpi0[0]).GetP4() + (*fPARTpi0[1]).GetP4()
					+ (*fPARTpi0[2]).GetP4();
				fMassDpi0[0] = TMath::Abs( p4.M() - fParticleID->GetMassMeV( kEta));
				if( fMassDpi0[0] < fMaxMDeta ){
					(*fPARTeta[0]).GetP4() = p4;
					fNeta = 1;
					fNpi0 = 0;
				}
			}*/
			break;
	}

//	MarkEndBuffer();

	if ( fNphoton == 1) {

		TA2Particle phot = *fPARTphoton[0];
		fPhotKE = phot.GetT();
		fPhotTheta = phot.GetThetaDg();
		fPhotPhi = phot.GetPhiDg();

	}

	if ( fNproton == 1) {

		TA2Particle prot = *fPARTproton[0];
		fProtKE = prot.GetT();
		fProtTheta = prot.GetThetaDg();
		fProtPhi = prot.GetPhiDg();

	}

	if ( fNpi0 == 1) {

		TA2Particle pi0 = *fPARTpi0[0];
		fPi0KE = pi0.GetT();
		fPi0Theta = pi0.GetThetaDg();
		fPi0Phi = pi0.GetPhiDg();

	}

	MarkEndBuffer();
}
