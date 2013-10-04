//--Author	PP Martel    10th Oct 2010
//--Update	PP Martel    19th Feb 2011
//
// TA2SpinPolPhysics
//
// Reconstruction of Compton and Pi0 kinematics. 

#include "TA2SpinPolPhysics.h"

// Valid Keywords for command-line setup of Compton
enum { EPhotoMassLimits = 1000, EPromRandWindows, ESaveEventTree,
       ESaveCompWTree, ESaveCompMTree, ESavePionWTree, ESavePionMTree};
static const Map_t kPhotoKeys[] = {
  {"Mass-Limits:",			EPhotoMassLimits},
  {"Prompt-Random-Windows:",		EPromRandWindows},
  {"SaveEventTree:",	                ESaveEventTree},
  {"SaveCompWTree:",	                ESaveCompWTree},
  {"SaveCompMTree:",	                ESaveCompMTree},
  {"SavePionWTree:",	                ESavePionWTree},
  {"SavePionMTree:",	                ESavePionMTree},
  {NULL,                                  -1}
};

ClassImp(TA2SpinPolPhysics)

//-----------------------------------------------------------------------------
TA2SpinPolPhysics::TA2SpinPolPhysics( const char* name, TA2Analysis* analysis ):TA2Physics( name, analysis ) {
  // Initialise Compton variables here
  // Default null pointers, zeroed variables

  fEventSave = false;
  fCompWSave = false;
  fCompMSave = false;
  fPionWSave = false;
  fPionMSave = false;

  fTableCB = false;
  fTblTAPS = false;

  fTAGG = NULL;
  fCB = NULL;
  fTAPS = NULL;
  
  fPARTtagged = NULL; 
  fPARTphoton = NULL;
  fPARTproton = NULL;
  fPARTpiplus = NULL;
  fPARTneutron = NULL;
  fPARTrootino = NULL;
  
  fPARTpi0 = NULL;
  fPARTeta = NULL;
  fPARTgprime = NULL;

  fTgRefTDC = ENullHit;
  fCBRefTDC = ENullHit;
  fSynchDif = ENullHit;

  fCherADC = ENullHit;
  fCherTDC0 = ENullHit;
  fCherTDC1 = ENullHit;
  fCherTDC2 = ENullHit;
  
  fBeamPol = 0;

  fNphoton = 0;
  fNproton = 0;
  fNpiplus = 0;
  fNneutron = 0;
  fNrootino = 0;
  fNgprime = 0;
  fNprompt = 0;
  fNrandom = 0;
  fMaxTagg = 0;
  
  fNpi0 = 0;
  fNeta = 0;
  
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

  fProtCB = false;
  fPrTAPS = false;
  fIsComp = false;
  fIsPion = false;

  // Requires at least one photon and one proton
  
  fPhotEk = ENullHit;
  fPhotPx = ENullHit;
  fPhotPy = ENullHit;
  fPhotPz = ENullHit;
  fPhotTh = ENullHit;
  fPhotPh = ENullHit;
  fPhotTm = ENullHit;
  
  fPionEk = ENullHit;
  fPionPx = ENullHit;
  fPionPy = ENullHit;
  fPionPz = ENullHit;
  fPionMa = ENullHit;
  fPionTh = ENullHit;
  fPionPh = ENullHit;
  fPionTm = ENullHit;
  
  fDec1Ek = ENullHit;
  fDec1Px = ENullHit;
  fDec1Py = ENullHit;
  fDec1Pz = ENullHit;
  fDec1Th = ENullHit;
  fDec1Ph = ENullHit;
  fDec1Tm = ENullHit;
  
  fDec2Ek = ENullHit;
  fDec2Px = ENullHit;
  fDec2Py = ENullHit;
  fDec2Pz = ENullHit;
  fDec2Th = ENullHit;
  fDec2Ph = ENullHit;
  fDec2Tm = ENullHit;

  fPionOA = ENullHit;
  
  fProtUn = ENullHit;
  fProtEk = ENullHit;
  fProtPx = ENullHit;
  fProtPy = ENullHit;
  fProtPz = ENullHit;
  fProtTh = ENullHit;
  fProtPh = ENullHit;
  fProtTm = ENullHit;

  fMPrToPh = ENullHit;

  fTaggChP = NULL;
  fTaggChR = NULL;
  fTaggEkP = NULL;
  fTaggEkR = NULL;
  fTaggTmP = NULL;
  fTaggTmR = NULL;

  fDifTime = NULL;

  fRecoEkP = NULL;
  fRecoEkR = NULL;
  fRecoPx = ENullHit;
  fRecoPy = ENullHit;
  fRecoPzP = NULL;
  fRecoPzR = NULL;
  fRecoThP = NULL;
  fRecoThR = NULL;
  fRecoPhP = NULL;
  fRecoPhR = NULL;
  fRecoMaP = NULL;
  fRecoMaR = NULL;

  fMissEkP = NULL;
  fMissEkR = NULL;
  fMissPx = ENullHit;
  fMissPy = ENullHit;
  fMissPzP = NULL;
  fMissPzR = NULL;
  fMissPt = ENullHit;
  fMissPrP = NULL;
  fMissPrR = NULL;

  fProtOAP = NULL;
  fProtOAR = NULL;
    
  AddCmdList( kPhotoKeys );       // command-line recognition for SetConfig()
}


//-----------------------------------------------------------------------------
TA2SpinPolPhysics::~TA2SpinPolPhysics() {
  // Free up allocated memory...after checking its allocated
  // detector and cuts lists

  delete fEventTree;
  delete fEventFile;

  delete fCompWTree;
  delete fCompWFile;

  delete fCompMTree;
  delete fCompMFile;

  delete fPionWTree;
  delete fPionWFile;

  delete fPionMTree;
  delete fPionMFile;
  
}

//---------------------------------------------------------------------------
void TA2SpinPolPhysics::CloseFile()
{
  if( fEventSave )
  {
    fEventFile->cd();
    fEventTree->Write();
    fEventFile->Close();
    printf("Event Tree saved to %s\n",fEventFileName);
  }  
  if( fCompWSave )
  {
    fCompWFile->cd();
    fCompWTree->Write();
    fCompWFile->Close();
    printf("Compton Tree (with proton) saved to %s\n",fCompWFileName);
  }  
  if( fCompMSave )
  {
    fCompMFile->cd();
    fCompMTree->Write();
    fCompMFile->Close();
    printf("Compton Tree (missing proton) saved to %s\n",fCompMFileName);
  }  
  if( fPionWSave )
  {
    fPionWFile->cd();
    fPionWTree->Write();
    fPionWFile->Close();
    printf("Pi0 Tree (with proton) saved to %s\n",fPionWFileName);
  }  
  if( fPionMSave )
  {
    fPionMFile->cd();
    fPionMTree->Write();
    fPionMFile->Close();
    printf("Pi0 Tree (missing proton) saved to %s\n",fPionMFileName);
  }  
}

//-----------------------------------------------------------------------------
void TA2SpinPolPhysics::SetConfig(Char_t* line, Int_t key) {
  // Any special command-line input for Crystal Ball apparatus
  
  switch (key){
  case EPhotoMassLimits:
    //  Invariant mass limits
    if( sscanf( line, "%lf%lf", &fMaxMDpi0, &fMaxMDeta ) != 2 ) {
      PrintError( line, "<Compton meson invariant mass limits>");
      return;
    }
    break;
  case EPromRandWindows:
    //  Prompt-Random Windows
    if( sscanf( line, "%lf%lf%lf%lf%lf%lf", &fPromptL, &fPromptH, &fRand1Lo,
		&fRand1Hi, &fRand2Lo, &fRand2Hi) != 6 ){
      PrintError( line, "<Prompt-Random windows>");
      return;
    }
    break;
  case ESaveEventTree:
    //  Output Event Tree to ROOT file
    if( sscanf( line, "%s", fEventFileName) != 1 ) {
      PrintError( line, "<Output Event Tree to ROOT file>");
      return;
    }
    else fEventSave = true;
    break;
  case ESaveCompWTree:
    //  Output Compton Tree to ROOT file
    if( sscanf( line, "%s", fCompWFileName) != 1 ) {
      PrintError( line, "<Output Compton Tree to ROOT file>");
      return;
    }
    else fCompWSave = true;
    break;
  case ESaveCompMTree:
    //  Output Compton Tree to ROOT file
    if( sscanf( line, "%s", fCompMFileName) != 1 ) {
      PrintError( line, "<Output Compton Tree to ROOT file>");
      return;
    }
    else fCompMSave = true;
    break;
  case ESavePionWTree:
    //  Output Pion Tree to ROOT file
    if( sscanf( line, "%s", fPionWFileName) != 1 ) {
      PrintError( line, "<Output Pion Tree to ROOT file>");
      return;
    }
    else fPionWSave = true;
    break;
  case ESavePionMTree:
    //  Output Pion Tree to ROOT file
    if( sscanf( line, "%s", fPionMFileName) != 1 ) {
      PrintError( line, "<Output Pion Tree to ROOT file>");
      return;
    }
    else fPionMSave = true;
    break;
  default:
    // default main apparatus SetConfig()
    TA2Physics::SetConfig( line, key );
    break;
  }
}

//---------------------------------------------------------------------------
void TA2SpinPolPhysics::PostInit() {
  // Initialise arrays to contain 4 momenta and plotable scaler variables
  // Missing mass, missing energy, cm momentum, energies, angles
  // Initialisation will abort if CB or Tagger not initialised
  // TAPS is optional

  fCB = (TA2CB*)((TA2Analysis*)fParent)->GetChild("CB");
  if( !fCB ) PrintError("","<No CB class found in annalysis>",EErrFatal);
  else fCBpart = fCB->GetParticles();

  fTAGG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
  if( !fTAGG ) PrintError("","<No Tagger class found in analysis>",EErrFatal);
  else fTAGGpart = fTAGG->GetParticles();

  fLADD = (TA2Ladder*)((TA2Analysis*)fParent)->GetGrandChild("FPD");
  if( !fLADD ) PrintError("","<No Ladder class found in analysis>",EErrFatal);
  
  fTAPS = (TA2Taps*)((TA2Analysis*)fParent)->GetChild("TAPS");
  if( !fTAPS ) PrintError("Warning!!!","<No TAPS class found in annalysis>");
  else fTAPSpart = fTAPS->GetParticles();
  
  Int_t i;
  TA2Particle* part;
  
  // Maximum # of reaction particles
  Int_t maxparticle = fCB->GetMaxParticle();
  if( fTAPS ) maxparticle += fTAPS->GetMaxParticle();

  // Maximum # of tagger hits
  Int_t maxtagg = fTAGG->GetMaxParticle() + 1;
  fMaxTagg = maxtagg;

  fTaggCh = new Int_t[maxtagg * 2];
  fTaggChP = fTaggCh; fTaggChR = fTaggCh + maxtagg;
  fTaggEk = new Double_t[maxtagg * 2];
  fTaggEkP = fTaggEk; fTaggEkR = fTaggEk + maxtagg;
  fTaggTm = new Double_t[maxtagg * 2];
  fTaggTmP = fTaggTm; fTaggTmR = fTaggTm + maxtagg;

  fDifTime = new Double_t[maxtagg];

  fRecoEk = new Double_t[maxtagg * 2];
  fRecoEkP = fRecoEk; fRecoEkR = fRecoEk + maxtagg;
  fRecoPz = new Double_t[maxtagg * 2];
  fRecoPzP = fRecoPz; fRecoPzR = fRecoPz + maxtagg;
  fRecoTh = new Double_t[maxtagg * 2];
  fRecoThP = fRecoTh; fRecoThR = fRecoTh + maxtagg;
  fRecoPh = new Double_t[maxtagg * 2];
  fRecoPhP = fRecoPh; fRecoPhR = fRecoPh + maxtagg;
  fRecoMa = new Double_t[maxtagg * 2];
  fRecoMaP = fRecoMa; fRecoMaR = fRecoMa + maxtagg;

  fMissEk = new Double_t[maxtagg * 2];
  fMissEkP = fMissEk; fMissEkR = fMissEk + maxtagg;
  fMissPz = new Double_t[maxtagg * 2];
  fMissPzP = fMissPz; fMissPzR = fMissPz + maxtagg;
  fMissPr = new Double_t[maxtagg * 2];
  fMissPrP = fMissPr; fMissPrR = fMissPr + maxtagg;

  fProtOA = new Double_t[maxtagg * 2];
  fProtOAP = fProtOA; fProtOAR = fProtOA + maxtagg;
  
  // Particle from detectors
  fPARTtagged = new TA2Particle*[maxtagg];
  part = new TA2Particle[maxtagg];
  for( i=0; i<maxtagg; i++) fPARTtagged[i] = part + i;

  fPARTphoton = new TA2Particle*[maxparticle];
  fPARTproton = new TA2Particle*[maxparticle];
  fPARTpiplus = new TA2Particle*[maxparticle];
  fPARTneutron = new TA2Particle*[maxparticle];
  fPARTrootino = new TA2Particle*[maxparticle];
 
  // Pi0
  fPARTpi0 = new TA2Particle*[maxparticle];
  part = new TA2Particle[maxparticle];
  for( i=0; i<maxparticle; i++ ) fPARTpi0[i] = part + i;
  
  // Eta
  fPARTeta = new TA2Particle*[maxparticle];
  part = new TA2Particle[maxparticle];
  for( i = 0; i < maxparticle; i++ ) fPARTeta[i] = part + i;
  
  // Gamma prime
  fPARTgprime =  new TA2Particle*[maxparticle];
  part = new TA2Particle[maxparticle];
  for( i=0; i<maxparticle; i++ ) fPARTgprime[i] = part + i;

  // Arrays used to combine photons to mesons
  Int_t maxperm = 0;
  for( i=1; i<=maxparticle; i++ ) maxperm += i;
  fMassDpi0 = new Double_t[maxperm];
  fMassDeta = new Double_t[maxperm];
  fMassIJ = new Int_t[maxperm];
  fMassIpi0 = new Int_t[maxperm];
  fMassIeta = new Int_t[maxperm];
  fIsMesonIndex = new Bool_t[maxparticle];

  if( fEventSave ) {
    fEventFile = new TFile(fEventFileName, "RECREATE", "EventFile", 3);
    fEventTree = new TTree("EventTree", "Compton Kinematics");
    fEventTree->Branch("TgRefTDC", &fTgRefTDC,"TgRefTDC/I");
    fEventTree->Branch("CBRefTDC", &fCBRefTDC,"CBRefTDC/I");
    fEventTree->Branch("Event", &fEvent);
  }

  if( fCompWSave ) {
    fCompWFile = new TFile(fCompWFileName, "RECREATE", "CompWFile", 3);
    fCompWTree = new TTree("CompWTree", "Compton Kinematics");

    fCompWTree->Branch("NParts", &fNparticle, "NParts/I");

    fCompWTree->Branch("PhotEk", &fPhotEk, "PhotEk/D");
    fCompWTree->Branch("PhotPx", &fPhotPx, "PhotPx/D");
    fCompWTree->Branch("PhotPy", &fPhotPy, "PhotPy/D");
    fCompWTree->Branch("PhotPz", &fPhotPz, "PhotPz/D");
    fCompWTree->Branch("PhotTh", &fPhotTh, "PhotTh/D");
    fCompWTree->Branch("PhotPh", &fPhotPh, "PhotPh/D");
    fCompWTree->Branch("PhotTm", &fPhotTm, "PhotTm/D");

    fCompWTree->Branch("ProtUn", &fProtUn, "ProtUn/D");
    fCompWTree->Branch("ProtEk", &fProtEk, "ProtEk/D");
    fCompWTree->Branch("ProtPx", &fProtPx, "ProtPx/D");
    fCompWTree->Branch("ProtPy", &fProtPy, "ProtPy/D");
    fCompWTree->Branch("ProtPz", &fProtPz, "ProtPz/D");
    fCompWTree->Branch("ProtTh", &fProtTh, "ProtTh/D");
    fCompWTree->Branch("ProtPh", &fProtPh, "ProtPh/D");
    fCompWTree->Branch("ProtTm", &fProtTm, "ProtTm/D");

    fCompWTree->Branch("MPrToPh", &fMPrToPh, "MPrToPh/D");

    fCompWTree->Branch("RecoPx", &fRecoPx, "RecoPx/D");
    fCompWTree->Branch("RecoPy", &fRecoPy, "RecoPy/D");

    fCompWTree->Branch("MissPx", &fMissPx, "MissPx/D");
    fCompWTree->Branch("MissPy", &fMissPy, "MissPy/D");
    fCompWTree->Branch("MissPt", &fMissPt, "MissPt/D");

    fCompWTree->Branch("NPrompt", &fNprompt, "NPrompt/I");

    fCompWTree->Branch("TaggChP", fTaggChP, "TaggChP[NPrompt]/I");
    fCompWTree->Branch("TaggEkP", fTaggEkP, "TaggEkP[NPrompt]/D");
    fCompWTree->Branch("TaggTmP", fTaggTmP, "TaggTmP[NPrompt]/D");

    fCompWTree->Branch("RecoEkP", fRecoEkP, "RecoEkP[NPrompt]/D");
    fCompWTree->Branch("RecoPzP", fRecoPzP, "RecoPzP[NPrompt]/D");
    fCompWTree->Branch("RecoThP", fRecoThP, "RecoThP[NPrompt]/D");
    fCompWTree->Branch("RecoPhP", fRecoPhP, "RecoPhP[NPrompt]/D");
    fCompWTree->Branch("RecoMaP", fRecoMaP, "RecoMaP[NPrompt]/D");

    fCompWTree->Branch("MissEkP", fMissEkP, "MissEkP[NPrompt]/D");
    fCompWTree->Branch("MissPzP", fMissPzP, "MissPzP[NPrompt]/D");
    fCompWTree->Branch("MissPrP", fMissPrP, "MissPrP[NPrompt]/D");

    fCompWTree->Branch("ProtOAP", fProtOAP, "ProtOAP[NPrompt]/D");

    fCompWTree->Branch("NRandom", &fNrandom, "NRandom/I");

    fCompWTree->Branch("TaggChR", fTaggChR, "TaggChR[NRandom]/I");
    fCompWTree->Branch("TaggEkR", fTaggEkR, "TaggEkR[NRandom]/D");
    fCompWTree->Branch("TaggTmR", fTaggTmR, "TaggTmR[NRandom]/D");

    fCompWTree->Branch("RecoEkR", fRecoEkR, "RecoEkR[NRandom]/D");
    fCompWTree->Branch("RecoPzR", fRecoPzR, "RecoPzR[NRandom]/D");
    fCompWTree->Branch("RecoThR", fRecoThR, "RecoThR[NRandom]/D");
    fCompWTree->Branch("RecoPhR", fRecoPhR, "RecoPhR[NRandom]/D");
    fCompWTree->Branch("RecoMaR", fRecoMaR, "RecoMaR[NRandom]/D");

    fCompWTree->Branch("MissEkR", fMissEkR, "MissEkR[NRandom]/D");
    fCompWTree->Branch("MissPzR", fMissPzR, "MissPzR[NRandom]/D");
    fCompWTree->Branch("MissPrR", fMissPrR, "MissPrR[NRandom]/D");

    fCompWTree->Branch("ProtOAR", fProtOAR, "ProtOAR[NRandom]/D");

    fCompWTree->Branch("TgRefTDC", &fTgRefTDC, "TgRefTDC/I");
    fCompWTree->Branch("CBRefTDC", &fCBRefTDC, "CBRefTDC/I");

    fCompWTree->Branch("CherADC", &fCherADC, "CherADC/I");
    fCompWTree->Branch("CherTDC0", &fCherTDC0, "CherTDC0/I");
    fCompWTree->Branch("CherTDC1", &fCherTDC1, "CherTDC1/I");
    fCompWTree->Branch("CherTDC2", &fCherTDC2, "CherTDC2/I");

    fCompWTree->Branch("BeamPol", &fBeamPol, "BeamPol/I");

    fCompWTree->Branch("ProtCB", &fProtCB, "ProtCB/O");
    fCompWTree->Branch("PrTAPS", &fPrTAPS, "PrTAPS/O");

  }

  if( fCompMSave ) {
    fCompMFile = new TFile(fCompMFileName, "RECREATE", "CompMFile", 3);
    fCompMTree = new TTree("CompMTree", "Compton Kinematics");

    fCompMTree->Branch("NParts", &fNparticle, "NParts/I");

    fCompMTree->Branch("PhotEk", &fPhotEk, "PhotEk/D");
    fCompMTree->Branch("PhotPx", &fPhotPx, "PhotPx/D");
    fCompMTree->Branch("PhotPy", &fPhotPy, "PhotPy/D");
    fCompMTree->Branch("PhotPz", &fPhotPz, "PhotPz/D");
    fCompMTree->Branch("PhotTh", &fPhotTh, "PhotTh/D");
    fCompMTree->Branch("PhotPh", &fPhotPh, "PhotPh/D");
    fCompMTree->Branch("PhotTm", &fPhotTm, "PhotTm/D");

    fCompMTree->Branch("MPrToPh", &fMPrToPh, "MPrToPh/D");

    fCompMTree->Branch("RecoPx", &fRecoPx, "RecoPx/D");
    fCompMTree->Branch("RecoPy", &fRecoPy, "RecoPy/D");

    fCompMTree->Branch("NPrompt", &fNprompt, "NPrompt/I");

    fCompMTree->Branch("TaggChP", fTaggChP, "TaggChP[NPrompt]/I");
    fCompMTree->Branch("TaggEkP", fTaggEkP, "TaggEkP[NPrompt]/D");
    fCompMTree->Branch("TaggTmP", fTaggTmP, "TaggTmP[NPrompt]/D");

    fCompMTree->Branch("RecoEkP", fRecoEkP, "RecoEkP[NPrompt]/D");
    fCompMTree->Branch("RecoPzP", fRecoPzP, "RecoPzP[NPrompt]/D");
    fCompMTree->Branch("RecoThP", fRecoThP, "RecoThP[NPrompt]/D");
    fCompMTree->Branch("RecoPhP", fRecoPhP, "RecoPhP[NPrompt]/D");
    fCompMTree->Branch("RecoMaP", fRecoMaP, "RecoMaP[NPrompt]/D");

    fCompMTree->Branch("NRandom", &fNrandom, "NRandom/I");

    fCompMTree->Branch("TaggChR", fTaggChR, "TaggChR[NRandom]/I");
    fCompMTree->Branch("TaggEkR", fTaggEkR, "TaggEkR[NRandom]/D");
    fCompMTree->Branch("TaggTmR", fTaggTmR, "TaggTmR[NRandom]/D");

    fCompMTree->Branch("RecoEkR", fRecoEkR, "RecoEkR[NRandom]/D");
    fCompMTree->Branch("RecoPzR", fRecoPzR, "RecoPzR[NRandom]/D");
    fCompMTree->Branch("RecoThR", fRecoThR, "RecoThR[NRandom]/D");
    fCompMTree->Branch("RecoPhR", fRecoPhR, "RecoPhR[NRandom]/D");
    fCompMTree->Branch("RecoMaR", fRecoMaR, "RecoMaR[NRandom]/D");

    fCompMTree->Branch("TgRefTDC", &fTgRefTDC, "TgRefTDC/I");
    fCompMTree->Branch("CBRefTDC", &fCBRefTDC, "CBRefTDC/I");

    fCompMTree->Branch("BeamPol", &fBeamPol, "BeamPol/I");

  }

  if( fPionWSave ) {
    fPionWFile = new TFile(fPionWFileName, "RECREATE", "PionWFile", 3);
    fPionWTree = new TTree("PionWTree", "Pion Kinematics");

    fPionWTree->Branch("NParts", &fNparticle, "NParts/I");

    fPionWTree->Branch("PionEk", &fPionEk, "PionEk/D");
    fPionWTree->Branch("PionPx", &fPionPx, "PionPx/D");
    fPionWTree->Branch("PionPy", &fPionPy, "PionPy/D");
    fPionWTree->Branch("PionPz", &fPionPz, "PionPz/D");
    fPionWTree->Branch("PionMa", &fPionMa, "PionMa/D");
    fPionWTree->Branch("PionTh", &fPionTh, "PionTh/D");
    fPionWTree->Branch("PionPh", &fPionPh, "PionPh/D");
    fPionWTree->Branch("PionTm", &fPionTm, "PionTm/D");

    fPionWTree->Branch("Dec1Ek", &fDec1Ek, "Dec1Ek/D");
    fPionWTree->Branch("Dec1Px", &fDec1Px, "Dec1Px/D");
    fPionWTree->Branch("Dec1Py", &fDec1Py, "Dec1Py/D");
    fPionWTree->Branch("Dec1Pz", &fDec1Pz, "Dec1Pz/D");
    fPionWTree->Branch("Dec1Th", &fDec1Th, "Dec1Th/D");
    fPionWTree->Branch("Dec1Ph", &fDec1Ph, "Dec1Ph/D");
    fPionWTree->Branch("Dec1Tm", &fDec1Tm, "Dec1Tm/D");

    fPionWTree->Branch("Dec2Ek", &fDec2Ek, "Dec2Ek/D");
    fPionWTree->Branch("Dec2Px", &fDec2Px, "Dec2Px/D");
    fPionWTree->Branch("Dec2Py", &fDec2Py, "Dec2Py/D");
    fPionWTree->Branch("Dec2Pz", &fDec2Pz, "Dec2Pz/D");
    fPionWTree->Branch("Dec2Th", &fDec2Th, "Dec2Th/D");
    fPionWTree->Branch("Dec2Ph", &fDec2Ph, "Dec2Ph/D");
    fPionWTree->Branch("Dec2Tm", &fDec2Tm, "Dec2Tm/D");

    fPionWTree->Branch("PionOA", &fPionOA, "PionOA/D");

    fPionWTree->Branch("ProtUn", &fProtUn, "ProtUn/D");
    fPionWTree->Branch("ProtEk", &fProtEk, "ProtEk/D");
    fPionWTree->Branch("ProtPx", &fProtPx, "ProtPx/D");
    fPionWTree->Branch("ProtPy", &fProtPy, "ProtPy/D");
    fPionWTree->Branch("ProtPz", &fProtPz, "ProtPz/D");
    fPionWTree->Branch("ProtTh", &fProtTh, "ProtTh/D");
    fPionWTree->Branch("ProtPh", &fProtPh, "ProtPh/D");
    fPionWTree->Branch("ProtTm", &fProtTm, "ProtTm/D");

    fPionWTree->Branch("RecoPx", &fRecoPx, "RecoPx/D");
    fPionWTree->Branch("RecoPy", &fRecoPy, "RecoPy/D");

    fPionWTree->Branch("MissPx", &fMissPx, "MissPx/D");
    fPionWTree->Branch("MissPy", &fMissPy, "MissPy/D");
    fPionWTree->Branch("MissPt", &fMissPt, "MissPt/D");

    fPionWTree->Branch("NPrompt", &fNprompt, "NPrompt/I");

    fPionWTree->Branch("TaggChP", fTaggChP, "TaggChP[NPrompt]/I");
    fPionWTree->Branch("TaggEkP", fTaggEkP, "TaggEkP[NPrompt]/D");
    fPionWTree->Branch("TaggTmP", fTaggTmP, "TaggTmP[NPrompt]/D");

    fPionWTree->Branch("RecoEkP", fRecoEkP, "RecoEkP[NPrompt]/D");
    fPionWTree->Branch("RecoPzP", fRecoPzP, "RecoPzP[NPrompt]/D");
    fPionWTree->Branch("RecoThP", fRecoThP, "RecoThP[NPrompt]/D");
    fPionWTree->Branch("RecoPhP", fRecoPhP, "RecoPhP[NPrompt]/D");
    fPionWTree->Branch("RecoMaP", fRecoMaP, "RecoMaP[NPrompt]/D");

    fPionWTree->Branch("MissEkP", fMissEkP, "MissEkP[NPrompt]/D");
    fPionWTree->Branch("MissPzP", fMissPzP, "MissPzP[NPrompt]/D");
    fPionWTree->Branch("MissPrP", fMissPrP, "MissPrP[NPrompt]/D");

    fPionWTree->Branch("ProtOAP", fProtOAP, "ProtOAP[NPrompt]/D");

    fPionWTree->Branch("NRandom", &fNrandom, "NRandom/I");

    fPionWTree->Branch("TaggChR", fTaggChR, "TaggChR[NRandom]/I");
    fPionWTree->Branch("TaggEkR", fTaggEkR, "TaggEkR[NRandom]/D");
    fPionWTree->Branch("TaggTmR", fTaggTmR, "TaggTmR[NRandom]/D");

    fPionWTree->Branch("RecoEkR", fRecoEkR, "RecoEkR[NRandom]/D");
    fPionWTree->Branch("RecoPzR", fRecoPzR, "RecoPzR[NRandom]/D");
    fPionWTree->Branch("RecoThR", fRecoThR, "RecoThR[NRandom]/D");
    fPionWTree->Branch("RecoPhR", fRecoPhR, "RecoPhR[NRandom]/D");
    fPionWTree->Branch("RecoMaR", fRecoMaR, "RecoMaR[NRandom]/D");

    fPionWTree->Branch("MissEkR", fMissEkR, "MissEkR[NRandom]/D");
    fPionWTree->Branch("MissPzR", fMissPzR, "MissPzR[NRandom]/D");
    fPionWTree->Branch("MissPrR", fMissPrR, "MissPrR[NRandom]/D");

    fPionWTree->Branch("ProtOAR", fProtOAR, "ProtOAR[NRandom]/D");

    fPionWTree->Branch("TgRefTDC", &fTgRefTDC, "TgRefTDC/I");
    fPionWTree->Branch("CBRefTDC", &fCBRefTDC, "CBRefTDC/I");

    fPionWTree->Branch("CherADC", &fCherADC, "CherADC/I");
    fPionWTree->Branch("CherTDC0", &fCherTDC0, "CherTDC0/I");
    fPionWTree->Branch("CherTDC1", &fCherTDC1, "CherTDC1/I");
    fPionWTree->Branch("CherTDC2", &fCherTDC2, "CherTDC2/I");

    fPionWTree->Branch("BeamPol", &fBeamPol, "BeamPol/I");

    fPionWTree->Branch("ProtCB", &fProtCB, "ProtCB/O");
    fPionWTree->Branch("PrTAPS", &fPrTAPS, "PrTAPS/O");

  }

  if( fPionMSave ) {
    fPionMFile = new TFile(fPionMFileName, "RECREATE", "PionMFile", 3);
    fPionMTree = new TTree("PionMTree", "Pion Kinematics");

    fPionMTree->Branch("NParts", &fNparticle, "NParts/I");

    fPionMTree->Branch("PionEk", &fPionEk, "PionEk/D");
    fPionMTree->Branch("PionPx", &fPionPx, "PionPx/D");
    fPionMTree->Branch("PionPy", &fPionPy, "PionPy/D");
    fPionMTree->Branch("PionPz", &fPionPz, "PionPz/D");
    fPionMTree->Branch("PionMa", &fPionMa, "PionMa/D");
    fPionMTree->Branch("PionTh", &fPionTh, "PionTh/D");
    fPionMTree->Branch("PionPh", &fPionPh, "PionPh/D");
    fPionMTree->Branch("PionTm", &fPionTm, "PionTm/D");

    fPionMTree->Branch("Dec1Ek", &fDec1Ek, "Dec1Ek/D");
    fPionMTree->Branch("Dec1Px", &fDec1Px, "Dec1Px/D");
    fPionMTree->Branch("Dec1Py", &fDec1Py, "Dec1Py/D");
    fPionMTree->Branch("Dec1Pz", &fDec1Pz, "Dec1Pz/D");
    fPionMTree->Branch("Dec1Th", &fDec1Th, "Dec1Th/D");
    fPionMTree->Branch("Dec1Ph", &fDec1Ph, "Dec1Ph/D");
    fPionMTree->Branch("Dec1Tm", &fDec1Tm, "Dec1Tm/D");

    fPionMTree->Branch("Dec2Ek", &fDec2Ek, "Dec2Ek/D");
    fPionMTree->Branch("Dec2Px", &fDec2Px, "Dec2Px/D");
    fPionMTree->Branch("Dec2Py", &fDec2Py, "Dec2Py/D");
    fPionMTree->Branch("Dec2Pz", &fDec2Pz, "Dec2Pz/D");
    fPionMTree->Branch("Dec2Th", &fDec2Th, "Dec2Th/D");
    fPionMTree->Branch("Dec2Ph", &fDec2Ph, "Dec2Ph/D");
    fPionMTree->Branch("Dec2Tm", &fDec2Tm, "Dec2Tm/D");

    fPionMTree->Branch("PionOA", &fPionOA, "PionOA/D");

    fPionMTree->Branch("RecoPx", &fRecoPx, "RecoPx/D");
    fPionMTree->Branch("RecoPy", &fRecoPy, "RecoPy/D");

    fPionMTree->Branch("NPrompt", &fNprompt, "NPrompt/I");

    fPionMTree->Branch("TaggChP", fTaggChP, "TaggChP[NPrompt]/I");
    fPionMTree->Branch("TaggEkP", fTaggEkP, "TaggEkP[NPrompt]/D");
    fPionMTree->Branch("TaggTmP", fTaggTmP, "TaggTmP[NPrompt]/D");

    fPionMTree->Branch("RecoEkP", fRecoEkP, "RecoEkP[NPrompt]/D");
    fPionMTree->Branch("RecoPzP", fRecoPzP, "RecoPzP[NPrompt]/D");
    fPionMTree->Branch("RecoThP", fRecoThP, "RecoThP[NPrompt]/D");
    fPionMTree->Branch("RecoPhP", fRecoPhP, "RecoPhP[NPrompt]/D");
    fPionMTree->Branch("RecoMaP", fRecoMaP, "RecoMaP[NPrompt]/D");

    fPionMTree->Branch("NRandom", &fNrandom, "NRandom/I");

    fPionMTree->Branch("TaggChR", fTaggChR, "TaggChR[NRandom]/I");
    fPionMTree->Branch("TaggEkR", fTaggEkR, "TaggEkR[NRandom]/D");
    fPionMTree->Branch("TaggTmR", fTaggTmR, "TaggTmR[NRandom]/D");

    fPionMTree->Branch("RecoEkR", fRecoEkR, "RecoEkR[NRandom]/D");
    fPionMTree->Branch("RecoPzR", fRecoPzR, "RecoPzR[NRandom]/D");
    fPionMTree->Branch("RecoThR", fRecoThR, "RecoThR[NRandom]/D");
    fPionMTree->Branch("RecoPhR", fRecoPhR, "RecoPhR[NRandom]/D");
    fPionMTree->Branch("RecoMaR", fRecoMaR, "RecoMaR[NRandom]/D");

    fPionMTree->Branch("TgRefTDC", &fTgRefTDC, "TgRefTDC/I");
    fPionMTree->Branch("CBRefTDC", &fCBRefTDC, "CBRefTDC/I");

    fPionMTree->Branch("BeamPol", &fBeamPol, "BeamPol/I");

  }

  gROOT->cd();

  fTableCB = MakeTable("CB");
  fTblTAPS = MakeTable("TAPS");

  // Default physics initialisation
  TA2Physics::PostInit();
}

//-----------------------------------------------------------------------------
void TA2SpinPolPhysics::LoadVariable( ) {
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
  TA2DataManager::LoadVariable("TgRefTDC", &fTgRefTDC, EISingleX);
  TA2DataManager::LoadVariable("CBRefTDC", &fCBRefTDC, EISingleX);
  TA2DataManager::LoadVariable("SynchDif", &fSynchDif, EISingleX);
  TA2DataManager::LoadVariable("CherADC", &fCherADC, EISingleX);
  TA2DataManager::LoadVariable("CherTDC0", &fCherTDC0, EDSingleX);
  TA2DataManager::LoadVariable("CherTDC1", &fCherTDC1, EDSingleX);
  TA2DataManager::LoadVariable("CherTDC2", &fCherTDC2, EDSingleX);

  TA2DataManager::LoadVariable("BeamPol", &fBeamPol, EISingleX);

  TA2DataManager::LoadVariable("Nparticle", &fNparticle, EISingleX);
  TA2DataManager::LoadVariable("Nphoton", &fNphoton, EISingleX);
  TA2DataManager::LoadVariable("Nproton", &fNproton, EISingleX);
  TA2DataManager::LoadVariable("Npiplus", &fNpiplus, EISingleX);
  TA2DataManager::LoadVariable("Nneutron", &fNneutron, EISingleX);
  TA2DataManager::LoadVariable("Nrootino", &fNrootino, EISingleX);
  TA2DataManager::LoadVariable("Ngprime", &fNgprime, EISingleX);
  TA2DataManager::LoadVariable("Nprompt", &fNprompt, EISingleX);
  TA2DataManager::LoadVariable("Nrandom", &fNrandom, EISingleX);
  TA2DataManager::LoadVariable("Npi0", &fNpi0, EISingleX);
  TA2DataManager::LoadVariable("Neta", &fNeta, EISingleX);

  TA2DataManager::LoadVariable("M2g", &fM2g, EDSingleX);
  TA2DataManager::LoadVariable("M6g", &fM6g, EDSingleX);
  
  TA2DataManager::LoadVariable("PhotEk", &fPhotEk, EDSingleX);
  TA2DataManager::LoadVariable("PhotPx", &fPhotPx, EDSingleX);
  TA2DataManager::LoadVariable("PhotPy", &fPhotPy, EDSingleX);
  TA2DataManager::LoadVariable("PhotPz", &fPhotPz, EDSingleX);
  TA2DataManager::LoadVariable("PhotTh", &fPhotTh, EDSingleX);
  TA2DataManager::LoadVariable("PhotPh", &fPhotPh, EDSingleX);
  TA2DataManager::LoadVariable("PhotTm", &fPhotTm, EDSingleX);
  
  TA2DataManager::LoadVariable("PionEk", &fPionEk, EDSingleX);
  TA2DataManager::LoadVariable("PionPx", &fPionPx, EDSingleX);
  TA2DataManager::LoadVariable("PionPy", &fPionPy, EDSingleX);
  TA2DataManager::LoadVariable("PionPz", &fPionPz, EDSingleX);
  TA2DataManager::LoadVariable("PionMa", &fPionMa, EDSingleX);
  TA2DataManager::LoadVariable("PionTh", &fPionTh, EDSingleX);
  TA2DataManager::LoadVariable("PionPh", &fPionPh, EDSingleX);
  TA2DataManager::LoadVariable("PionTm", &fPionTm, EDSingleX);

  TA2DataManager::LoadVariable("PionOA", &fPionOA, EDSingleX);
  
  TA2DataManager::LoadVariable("ProtUn", &fProtUn, EDSingleX);
  TA2DataManager::LoadVariable("ProtEk", &fProtEk, EDSingleX);
  TA2DataManager::LoadVariable("ProtPx", &fProtPx, EDSingleX);
  TA2DataManager::LoadVariable("ProtPy", &fProtPy, EDSingleX);
  TA2DataManager::LoadVariable("ProtPz", &fProtPz, EDSingleX);
  TA2DataManager::LoadVariable("ProtTh", &fProtTh, EDSingleX);
  TA2DataManager::LoadVariable("ProtPh", &fProtPh, EDSingleX);
  TA2DataManager::LoadVariable("ProtTm", &fProtTm, EDSingleX);

  TA2DataManager::LoadVariable("MPrToPh", &fMPrToPh, EDSingleX);
  
  TA2DataManager::LoadVariable("TaggChP", fTaggChP, EIMultiX);
  TA2DataManager::LoadVariable("TaggChR", fTaggChR, EIMultiX);
  TA2DataManager::LoadVariable("TaggEkP", fTaggEkP, EDMultiX);
  TA2DataManager::LoadVariable("TaggEkR", fTaggEkR, EDMultiX);
  TA2DataManager::LoadVariable("TaggTmP", fTaggTmP, EDMultiX);
  TA2DataManager::LoadVariable("TaggTmR", fTaggTmR, EDMultiX);

  TA2DataManager::LoadVariable("DifTime", fDifTime, EDMultiX);
  
  TA2DataManager::LoadVariable("RecoEkP", fRecoEkP, EDMultiX);
  TA2DataManager::LoadVariable("RecoEkR", fRecoEkR, EDMultiX);
  TA2DataManager::LoadVariable("RecoPx", &fRecoPx, EDSingleX);
  TA2DataManager::LoadVariable("RecoPy", &fRecoPy, EDSingleX);
  TA2DataManager::LoadVariable("RecoPzP", fRecoPzP, EDMultiX);
  TA2DataManager::LoadVariable("RecoPzR", fRecoPzR, EDMultiX);
  TA2DataManager::LoadVariable("RecoThP", fRecoThP, EDMultiX);
  TA2DataManager::LoadVariable("RecoThR", fRecoThR, EDMultiX);
  TA2DataManager::LoadVariable("RecoPhP", fRecoPhP, EDMultiX);
  TA2DataManager::LoadVariable("RecoPhR", fRecoPhR, EDMultiX);
  TA2DataManager::LoadVariable("RecoMaP", fRecoMaP, EDMultiX);
  TA2DataManager::LoadVariable("RecoMaR", fRecoMaR, EDMultiX);
  
  TA2DataManager::LoadVariable("MissEkP", fMissEkP, EDMultiX);
  TA2DataManager::LoadVariable("MissEkR", fMissEkR, EDMultiX);
  TA2DataManager::LoadVariable("MissPx", &fMissPx, EDSingleX);
  TA2DataManager::LoadVariable("MissPy", &fMissPy, EDSingleX);
  TA2DataManager::LoadVariable("MissPzP", fMissPzP, EDMultiX);
  TA2DataManager::LoadVariable("MissPzR", fMissPzR, EDMultiX);
  TA2DataManager::LoadVariable("MissPt", &fMissPt, EDSingleX);
  TA2DataManager::LoadVariable("MissPrP", fMissPrP, EDMultiX);
  TA2DataManager::LoadVariable("MissPrR", fMissPrR, EDMultiX);

  TA2DataManager::LoadVariable("ProtOAP", fProtOAP, EDMultiX);
  TA2DataManager::LoadVariable("ProtOAR", fProtOAR, EDMultiX);
  
  return;
}

//-----------------------------------------------------------------------------
void TA2SpinPolPhysics::Reconstruct() {
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
  Int_t j;

  Int_t ntagg = fTAGG->GetNparticle();         // # particles in Tagger
  Int_t ncb = fCB->GetNparticle();             // # particles in CB
  Int_t ntaps;
  if ( fTAPS ) ntaps = fTAPS->GetNparticle();  // # particles in TAPS
  else ntaps = 0;

  fNparticle = ncb + ntaps;                    // total number particles (hits)
  
  fTgRefTDC = ENullHit;
  fCBRefTDC = ENullHit;
  fSynchDif = ENullHit;
  
  fCherADC = ENullHit;
  fCherTDC0 = ENullHit;
  fCherTDC1 = ENullHit;
  fCherTDC2 = ENullHit;

  fBeamPol = 0;
  
  // zero particle counters
  fNphoton = 0;
  fNproton = 0;
  fNpiplus = 0;
  fNneutron = 0;
  fNrootino = 0;
  fNgprime = 0;
  fNprompt = 0;
  fNrandom = 0;
  
  fNpi0 = 0;
  fNeta = 0;

  fM2g = ENullHit;                             // zero 2-gamma inv. mass
  fM6g = ENullHit;                             // zero 6-gamma inv. mass

  fProtCB = kFALSE;
  fPrTAPS = kFALSE;

  fCBCrystN = ENullHit;
  fTAPSCstN = ENullHit;

  fPhotEk = ENullHit;
  fPhotPx = ENullHit;
  fPhotPy = ENullHit;
  fPhotPz = ENullHit;
  fPhotTh = ENullHit;
  fPhotPh = ENullHit;
  fPhotTm = ENullHit;

  fPionEk = ENullHit;
  fPionPx = ENullHit;
  fPionPy = ENullHit;
  fPionPz = ENullHit;
  fPionMa = ENullHit;
  fPionTh = ENullHit;
  fPionPh = ENullHit;
  fPionTm = ENullHit;
  
  fDec1Ek = ENullHit;
  fDec1Px = ENullHit;
  fDec1Py = ENullHit;
  fDec1Pz = ENullHit;
  fDec1Th = ENullHit;
  fDec1Ph = ENullHit;
  fDec1Tm = ENullHit;
  
  fDec2Ek = ENullHit;
  fDec2Px = ENullHit;
  fDec2Py = ENullHit;
  fDec2Pz = ENullHit;
  fDec2Th = ENullHit;
  fDec2Ph = ENullHit;
  fDec2Tm = ENullHit;

  fPionOA = ENullHit;
  
  fProtUn = ENullHit;
  fProtEk = ENullHit;
  fProtPx = ENullHit;
  fProtPy = ENullHit;
  fProtPz = ENullHit;
  fProtTh = ENullHit;
  fProtPh = ENullHit;
  fProtTm = ENullHit;

  fMPrToPh = ENullHit;
    
  fRecoPx = ENullHit;
  fRecoPy = ENullHit;

  fMissPx = ENullHit;
  fMissPy = ENullHit;
  fMissPt = ENullHit;
  
  fTgRefTDC = fADC[2000];
  fCBRefTDC = fADC[1400];
  fSynchDif = (fTgRefTDC - fCBRefTDC);

  fCherADC = fADC[124];
  fCherTDC0 = fMulti[2007]->GetHit(0);
  fCherTDC1 = fMulti[2007]->GetHit(1);
  fCherTDC2 = fMulti[2007]->GetHit(2);

  fBeamPol = fADC[6];

  fDifTime[0] = EBufferEnd;

  MarkEndBuffer();

  // Sort 4-momenta provided by apparati according to particle type

  // Tagger
  for ( i=0; i<ntagg; i++ ) {
    fPARTtagged[i] = fTAGGpart+i;
  }

  // CB
  for ( i=0; i<ncb; i++ ) {
    switch( (fCBpart+i)->GetParticleID() ) {   // PDG code
    case kGamma:                               // photon
      fPARTphoton[fNphoton] = fCBpart+i;       // include in photon list
      fNphoton++;
      break;
    case kProton:                              // proton
      fPARTproton[fNproton] = fCBpart+i;       // include in proton list
      fProtCB = kTRUE;
      fNproton++;
      break;
    case kPiPlus:                              // pi+
      fPARTpiplus[fNpiplus] = fCBpart+i;       // include in piplus list
      fNpiplus++;
      break;
    default:                                   // don't know
      fPARTrootino[fNrootino] = fCBpart+i;     // include in rootino list
      fNrootino++;                          
    }
  }
  
  // TAPS
  for ( i=0; i<ntaps; i++ ) {
    switch( (fTAPSpart+i)->GetParticleID() ) { // PDG code
    case kGamma:                               // photon
      fPARTphoton[fNphoton] = fTAPSpart+i;     // include in photon list
      fNphoton++;
      break;
    case kProton:                              // proton
      fPARTproton[fNproton] = fTAPSpart+i;     // include in proton list
      fPrTAPS = kTRUE;
      fNproton++;
      break;
    case kPiPlus:                              // pi+
      fPARTpiplus[fNpiplus] = fTAPSpart+i;     // include in piplus list
      fNpiplus++;
      break;
    default:                                   // don't know
      fPARTrootino[fNrootino] = fTAPSpart+i;   // include in rootino list
      fNrootino++;                          
    }
  }

  // Check if detected photons combine to give pi0 or eta
  TLorentzVector p4;
  switch( fNphoton ) {
  case 1:
    // Just 1 photon....assume it is a gamma-prime
    fPARTgprime[0] = fPARTphoton[0];
    fNgprime = 1;
    break;
  case 2:
    // 2 photons detected, fast check if they make a pi0 or eta
    Sort2Photon();
    break;
  default:
    // More than 2 photons 
    SortNPhoton();
    // Check for 3-pi0 eta decay mode
    if( fNpi0==3 ) {
      p4 = (*fPARTpi0[0]).GetP4() + (*fPARTpi0[1]).GetP4()
	+ (*fPARTpi0[2]).GetP4();
      fMassDpi0[0] = TMath::Abs( p4.M() - fParticleID->GetMassMeV( kEta));
      if( fMassDpi0[0]<fMaxMDeta ) {
	(*fPARTeta[0]).GetP4() = p4;
	fNeta = 1;
	fNpi0 = 0;
      }
    }
    break;
  }

  if(fEventSave){
    fEvent.Clear();
    for( i=0; i<ntagg; i++ ) {
      fEvent.AddBeam(*fPARTtagged[i]);
    }
    for( i=0; i<fNphoton; i++ ) {
      fEvent.AddParticle(*fPARTphoton[i]);
    }
    for( i=0; i<fNproton; i++ ) {
      fEvent.AddParticle(*fPARTproton[i]);
    }
    for( i=0; i<fNpiplus; i++ ) {
      fEvent.AddParticle(*fPARTpiplus[i]);
    }
    fEventTree->Fill();
    fEvent.SetEventNumber();
  }

  Int_t jtagg = 0;

  Double_t TaggTime;
  Double_t ScatTime;

  TA2Particle fTagged;
  TA2Particle fPhoton;
  TA2Particle fPion;
  TA2Particle fDecay1;
  TA2Particle fDecay2;
  TA2Particle fProton;

  TLorentzVector fTaggP4;
  TLorentzVector fPhotP4;
  TLorentzVector fPionP4;
  TLorentzVector fProtP4;
  TLorentzVector fTargP4 = fP4target[0];
  TLorentzVector fScatP4;
  TLorentzVector fRecoP4;
  TLorentzVector fMissP4;
  TLorentzVector fPrPhP4;     // photon mistakenly id'ed as proton

  if ( ( fNphoton == 1 ) || ( fNpi0 == 1 ) ) {

    if ( fNphoton == 1 ) {

      fIsComp = kTRUE;
      fIsPion = kFALSE;

      fPhoton = *fPARTphoton[0];
      
      fPhotEk = fPhoton.GetT();
      fPhotPx = fPhoton.GetPx();
      fPhotPy = fPhoton.GetPy();
      fPhotPz = fPhoton.GetPz();
      fPhotTh = fPhoton.GetThetaDg();
      fPhotPh = fPhoton.GetPhiDg();
      fPhotTm = fPhoton.GetTime();
      
      fPhotP4 = fPhoton.GetP4();
      fScatP4 = fPhoton.GetP4();

      ScatTime = fPhotTm;

      fRecoPx = -fPhotPx;
      fRecoPy = -fPhotPy;
      
    }
    
    else if ( fNpi0 == 1 ) {

      fIsComp = kFALSE;
      fIsPion = kTRUE;
      
      fPion = *fPARTpi0[0];

      if( fNphoton == 2 ) {
	fDecay1 = *fPARTphoton[0];
	fDecay2 = *fPARTphoton[1];
      }
      else {
	j = fNphoton;
	for( i=0; i<fNphoton; i++ ) {
	  if( fIsMesonIndex[i] && ( i<j ) ) {
	    fDecay1 = *fPARTphoton[i];
	    j = i;
	  }
	  else if( fIsMesonIndex[i] && ( i>j ) ) {
	    fDecay2 = *fPARTphoton[i];
	  }
	}
      }


      fPionEk = fPion.GetT();
      fPionPx = fPion.GetPx();
      fPionPy = fPion.GetPy();
      fPionPz = fPion.GetPz();
      fPionMa = fPion.GetM();
      fPionTh = fPion.GetThetaDg();
      fPionPh = fPion.GetPhiDg();
      fPionTm = fPion.GetTime();
      
      fPionP4 = fPion.GetP4();
      fScatP4 = fPion.GetP4();

      ScatTime = fPionTm;
      
      fDec1Ek = fDecay1.GetT();
      fDec1Px = fDecay1.GetPx();
      fDec1Py = fDecay1.GetPy();
      fDec1Pz = fDecay1.GetPz();
      fDec1Th = fDecay1.GetThetaDg();
      fDec1Ph = fDecay1.GetPhiDg();
      fDec1Tm = fDecay1.GetTime();
      
      fDec2Ek = fDecay2.GetT();
      fDec2Px = fDecay2.GetPx();
      fDec2Py = fDecay2.GetPy();
      fDec2Pz = fDecay2.GetPz();
      fDec2Th = fDecay2.GetThetaDg();
      fDec2Ph = fDecay2.GetPhiDg();
      fDec2Tm = fDecay2.GetTime();

      fPionOA = (fDecay1.GetVect()).Angle(fDecay2.GetVect())*TMath::RadToDeg();
      
      fRecoPx = -fPionPx;
      fRecoPy = -fPionPy;
      
    }

    if ( fNproton == 1 ) {
    
      fProton = *fPARTproton[0];

      if ( fIsComp ) {
      
	fPrPhP4.SetPxPyPzE(1,1,1,fProton.GetT());
	fPrPhP4.SetRho(fProton.GetT());
	fPrPhP4.SetTheta(fProton.GetTheta());
	fPrPhP4.SetPhi(fProton.GetPhi());
	
	fPionP4 = (fPrPhP4 + fPhotP4);
	fMPrToPh = fPionP4.M();

      }

      Double_t TempT, CorrT, LossT;
      
      fProtUn = fProton.GetT();
      TempT = fProtUn;
      
      if( fProtCB ) {
	fCBCrystN = fProton.GetCentralIndex();
      }
      else if( fPrTAPS ) {
	fTAPSCstN = fProton.GetCentralIndex();
      }
      
      if( fTableCB && fProtCB ) {
	if( gAR->GetProcessType() == EMCProcess ) CorrT = 0;
	else if( TempT<30 ) CorrT = 0.3;
	else if( TempT>300 ) CorrT = -0.125;
	else CorrT = fECorrCB->Eval(TempT);
	//std::cout << "CB   : "<< TempT << "\t" << CorrT << "\t";
	TempT = (TempT*(1+CorrT));
	LossT = fELossCB[fCBCrystN]->Eval(TempT);
	//std::cout << TempT << "\t" << LossT << "\t";
	TempT = (TempT+LossT);
	//std::cout << TempT << std::endl;
	fProton.SetKinetic(TempT);
      }
      else if( fTblTAPS && fPrTAPS ) {
	if( gAR->GetProcessType() == EMCProcess ) CorrT = 0;
	else if( TempT<30 ) CorrT = 0.3;
	else if( TempT>300 ) CorrT = -0.125;
	else CorrT = fECoTAPS->Eval(TempT);
	//std::cout << "TAPS   : "<< TempT << "\t" << CorrT << "\t";
	TempT = (TempT*(1+CorrT));
	LossT = fELoTAPS[fTAPSCstN]->Eval(TempT);
	//std::cout << TempT << "\t" << LossT << "\t";
	TempT = (TempT+LossT);
	//std::cout << TempT << std::endl;
	fProton.SetKinetic(TempT);
      }
      
      fProtEk = fProton.GetT();
      fProtPx = fProton.GetPx();
      fProtPy = fProton.GetPy();
      fProtPz = fProton.GetPz();
      fProtTh = fProton.GetThetaDg();
      fProtPh = fProton.GetPhiDg();
      fProtTm = fProton.GetTime();
      
      fProtP4 = fProton.GetP4();

      if ( fIsComp ) {

	fMissPx = -(fPhotPx + fProtPx);
	fMissPy = -(fPhotPy + fProtPy);

      }

      if ( fIsPion ) {

	fMissPx = -(fPionPx + fProtPx);
	fMissPy = -(fPionPy + fProtPy);

      }

      fMissPt = sqrt(Sqr(fMissPx) + Sqr(fMissPy));

    }
  
    // Tagger Loop
    for ( i=0; i<ntagg; i++ ) {
      
      fTagged = *fPARTtagged[i];

      TaggTime = fTagged.GetTime();
      fDifTime[i] = TaggTime - ScatTime;
      
      // Prompt or Monte Carlo
      if( ( ( fDifTime[i]>=fPromptL ) && ( fDifTime[i]<=fPromptH ) ) ||
	  ( gAR->GetProcessType() == EMCProcess ) ) {
	
	jtagg = fNprompt++;
	
      }
      
      // Random
      else if ( ( ( fDifTime[i]>=fRand1Lo ) && ( fDifTime[i]<=fRand1Hi ) ) ||
		( ( fDifTime[i]>=fRand2Lo ) && ( fDifTime[i]<=fRand2Hi ) ) ) {
	
	jtagg = fMaxTagg + fNrandom++;
	
      }
      
      // Other - skip the rest!
      else continue;

      fTaggP4 = fTagged.GetP4();
    
      fTaggCh[jtagg] = fTagged.GetCentralIndex();
      fTaggEk[jtagg] = fTagged.GetT();
      fTaggTm[jtagg] = fTagged.GetTime();
    
      fRecoP4 = fTaggP4+fTargP4-fScatP4;
      fRecoEk[jtagg] = fRecoP4.E()-fRecoP4.M();
      fRecoPz[jtagg] = fRecoP4.Pz();
      fRecoMa[jtagg] = fRecoP4.M();

      if ( fNproton == 1 ) {
	fMissP4 = fTaggP4+fTargP4-fScatP4-fProtP4;
	fMissEk[jtagg] = fMissP4.E();
	fMissPz[jtagg] = fMissP4.Pz();
	fMissPr[jtagg] = fMissP4.P();
	
	fProtOA[jtagg] = fRecoP4.Vect().Angle(fProtP4.Vect())*TMath::RadToDeg();

      }

    }
    
    if( fCompWSave && fIsComp && ( fNproton == 1 ) ) fCompWTree->Fill();
    if( fCompMSave && fIsComp && ( fNproton == 0 ) ) fCompMTree->Fill();
    if( fPionWSave && fIsPion && ( fNproton == 1 ) ) fPionWTree->Fill();
    if( fPionMSave && fIsPion && ( fNproton == 0 ) ) fPionMTree->Fill();

  }

  fDifTime[i] = EBufferEnd;
  
  MarkEndBuffer();

}
