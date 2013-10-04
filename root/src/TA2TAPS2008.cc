////////////////////////////////////////////////////////////////////////////
//                                                                        //
//     --Author		H. Berghaeuser June 2008                             //
//                      HenningBerghaeuser@web.de                         //
//     -- Update        H. Berghaeuser Okt 2008                           //
//                                                                        //
//     -- Use the TAPS2008v0.YX.dat !                                     //
//                                                                        //
//            TA2TAPS2008.cc       version 0.1                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TA2TAPS2008.h"
#include "TA2Analysis.h"
#include <iostream>

enum {ETAPS_BaF2=5000, ETAPS_Veto, ETAPS_PbWO4, ETAPS_fixone};

enum {ETAPSFactor=5004, ETAPS_SimpleReconstruct};

enum {ETAPS_TOF_USEAGE=5006, ETAPS_TOF_Cuts, ETAPS_TOF_Neutron,
	ETAPS_TOF_Proton};

enum {ETAPS_TOF_Neutron_CutName=5010, ETAPS_TOF_Proton_CutName,
	ETAPS_TOF_ChPion, ETAPS_TOF_ChPion_CutName, ETAPS_TOF_Electron,
	ETAPS_TOF_Electron_CutName, ETAPS_TOF_Gamma, ETAPS_TOF_Gamma_CutName,
	ETAPS_dEvE_USEAGE};

enum {ETAPS_dEvE_Cuts=5020, ETAPS_dEvE_Proton, ETAPS_dEvE_Proton_CutName,
	ETAPS_dEvE_ChPion, ETAPS_dEvE_ChPion_CutName, ETAPS_dEvE_Electron,
	ETAPS_dEvE_Electron_CutName, ETAPS_PSA_USEAGE, ETAPS_PSA_Cuts,
	ETAPS_PSA_Nucleon, ETAPS_PSA_Nucleon_CutName, ETAPS_PSA_noNucleon,
	ETAPS_PSA_noNucleon_CutName, ETAPS_fReadPbWO4};

static const Map_t kTAPS2008Keys[] = {
//
	{"TAPS-Factor:",						ETAPSFactor},
//
	{"TAPS_SimpleReconstruct:",		ETAPS_SimpleReconstruct},
//
	{"TAPS_TOF_USEAGE:",					ETAPS_TOF_USEAGE},
	{"TAPS_TOF_Cuts:",					ETAPS_TOF_Cuts},
	{"TAPS_TOF_Neutron:",				ETAPS_TOF_Neutron},
	{"TAPS_TOF_Neutron_CutName:",		ETAPS_TOF_Neutron_CutName},
	{"TAPS_TOF_Proton:",					ETAPS_TOF_Proton},
	{"TAPS_TOF_Proton_CutName:",		ETAPS_TOF_Proton_CutName},
	{"TAPS_TOF_ChPion:",					ETAPS_TOF_ChPion},
	{"TAPS_TOF_ChPion_CutName:",		ETAPS_TOF_ChPion_CutName},
	{"TAPS_TOF_Electron:",				ETAPS_TOF_Electron},
	{"TAPS_TOF_Electron_CutName:",	ETAPS_TOF_Electron_CutName},
	{"TAPS_TOF_Gamma:",					ETAPS_TOF_Gamma},
	{"TAPS_TOF_Gamma_CutName:",		ETAPS_TOF_Gamma_CutName},
//
	{"TAPS_dEvE_USEAGE:",				ETAPS_dEvE_USEAGE},
	{"TAPS_dEvE_Cuts:",					ETAPS_dEvE_Cuts},
	{"TAPS_dEvE_Proton:",				ETAPS_dEvE_Proton},
	{"TAPS_dEvE_Proton_CutName:",		ETAPS_dEvE_Proton_CutName},
	{"TAPS_dEvE_ChPion:",				ETAPS_dEvE_ChPion},
	{"TAPS_dEvE_ChPion_CutName:",		ETAPS_dEvE_ChPion_CutName},
	{"TAPS_dEvE_Electron:",				ETAPS_dEvE_Electron},
	{"TAPS_dEvE_Electron_CutName:",	ETAPS_dEvE_Electron_CutName},
//
	{"TAPS_PSA_USEAGE:",					ETAPS_PSA_USEAGE},
	{"TAPS_PSA_Cuts:",					ETAPS_PSA_Cuts},
	{"TAPS_PSA_Nucleon:" ,				ETAPS_PSA_Nucleon},
	{"TAPS_PSA_Nucleon_CutName:" ,	ETAPS_PSA_Nucleon_CutName},
	{"TAPS_PSA_noNucleon:" ,			ETAPS_PSA_noNucleon },
	{"TAPS_PSA_noNucleon_CutName:",	ETAPS_PSA_noNucleon_CutName},
//
	{"TAPS_fReadPbWO4:",					ETAPS_fReadPbWO4},
//
	{NULL,									-1}
};

static Map_t kValidTAPS2008Detectors[] = {
	{"TA2TAPS_BaF2",		ETAPS_BaF2},
	{"TA2TAPS_PbWO4",		ETAPS_PbWO4},
	{"TA2PlasticPID",		ETAPS_Veto},
	{NULL,					-1}
};

ClassImp(TA2TAPS2008)

//-----------------------------------------------------------------------------

TA2TAPS2008::TA2TAPS2008(const char* name, TA2System* analysis)
		        :TA2Apparatus(name, analysis, kValidTAPS2008Detectors)
{
	fTapsRunStep = 0;
	fBaF2 = NULL;
	fVeto = NULL;

	fSimpleReconstruct  = kTRUE;

	fTAPS_TOF_USEAGE    = kFALSE;
	fTAPS_TOF_ChPion    = kFALSE;
	fTAPS_TOF_Proton    = kFALSE;
	fTAPS_TOF_Neutron   = kFALSE;
	fTAPS_TOF_Electron  = kFALSE;
	fTAPS_TOF_Gamma     = kFALSE;

	fTAPS_dEvE_USEAGE   = kFALSE;
	fTAPS_dEvE_Proton   = kFALSE;
	fTAPS_dEvE_ChPion   = kFALSE;
	fTAPS_dEvE_Electron = kFALSE;

	fTAPS_PSA_USEAGE    = kFALSE;
	fTAPS_PSA_Nucleon   = kFALSE;
	fTAPS_PSA_noNucleon = kFALSE;

	fTAPS_ReadPbWO4          = kFALSE;
//	fTAPS_OUTPUT_FILE_Exists = kFALSE;
//	fPrintSetUpInfo          = kTRUE;
	fPrintSetUpInfo          = kFALSE;

	fTimeShift          = 0.0; // maybe will later use SvensTimeShift
	fTapsFudge          = 1.0;

	fTimeOfFlight = NULL; 

	AddCmdList(kTAPS2008Keys);

	if ( fPrintSetUpInfo == kTRUE)
	{
		printf( "\n-- Calling TA2TAPS2008 ........... done\n");
		printf( " ... reading key words from TAPS-config file ...");
	}
}

//-----------------------------------------------------------------------------

TA2TAPS2008::~TA2TAPS2008()
{
	// Free up allocated memory
	if (fParticleID) delete fParticleID;
	if (fTimeOfFlight) delete fTimeOfFlight;
}

//-----------------------------------------------------------------------------

TA2DataManager* TA2TAPS2008::CreateChild(const char* name, Int_t dclass)
{

	if( !name) name = Map2String( dclass);
	switch( dclass)
	{
	 case ETAPS_BaF2:
		if ( fPrintSetUpInfo == kTRUE) {
		  printf("\n ... reading Detector-CHILD ETAPS_BaF2: ");
		  printf(name);
		}
		fBaF2 = new TA2TAPS_BaF2(name, this);
		return fBaF2;
	 case ETAPS_PbWO4:
		if( fPrintSetUpInfo == kTRUE) {
			printf("\n ... reading Detector-CHILD ETAPS_PbWO4: ");
			printf(name);
		}
		fPbWO4 = new TA2TAPS_PbWO4(name, this);
		return fPbWO4;
	 case ETAPS_Veto:
		if ( fPrintSetUpInfo == kTRUE) {
			printf("\n ... reading Detector-CHILD ETAPS_Veto: ");
			printf(name);
		}
		fVeto = new TA2PlasticPID(name, this);
		return fVeto;
	 default:
		return NULL;
	}
}

//-----------------------------------------------------------------------------

void TA2TAPS2008::LoadVariable()
{
	//                            name        pointer          type-spec
	//TA2DataManager::LoadVariable("NphotMinv", fM_Nphoton,      EDSingleX);
	//
	TA2DataManager::LoadVariable( "TimeOfFlight", fTimeOfFlight, EDMultiX);
	TA2Apparatus::LoadVariable();
}

//-----------------------------------------------------------------------------

void TA2TAPS2008::PostInit()
{
	if ( !fParticleID)
		PrintError("", "<Configuration aborted: ParticleID class MUST be "
				"specified>", EErrFatal);

	fCB = (TA2CrystalBall*)((TA2Analysis*)fParent)->GetChild("CB");
	if(!fCB)
		printf("\nError: not getting TA2CrystalBall in TA2TAPS2008 .... "
				"check class name PostInit()");

	fMaxParticle     = fBaF2->GetMaxCluster();
	fTimeOfFlight    = new Double_t[fMaxParticle+1];
	fTapsRunStep     = 1;
	fNoCBTimeCounter = 0;
	// fP4_Nphoton      = new TLorentzVector[fMaxParticle+1];
	// fM_Nphoton       = new Double_t[fMaxParticle+1];
	particles        = new TA2Particle[fMaxParticle+1];
	fIsVCharged      = new Bool_t[fMaxParticle+1];
	fdEvE_IsProton   = new Bool_t[fMaxParticle+1];
	fdEvE_IsChPion   = new Bool_t[fMaxParticle+1];
	fdEvE_IsElectron = new Bool_t[fMaxParticle+1];
	fTOF_IsProton    = new Bool_t[fMaxParticle+1];
	fTOF_IsNeutron   = new Bool_t[fMaxParticle+1];
	fTOF_IsChPion    = new Bool_t[fMaxParticle+1];
	fTOF_IsElectron  = new Bool_t[fMaxParticle+1];
	fTOF_IsGamma     = new Bool_t[fMaxParticle+1];
	fPSA_IsNucleon   = new Bool_t[fMaxParticle+1];
	fPSA_IsNoNucleon = new Bool_t[fMaxParticle+1];
	fPDG_ID_sec      = new Int_t[fMaxParticle+1];
	fPDG_ID_unclear  = new Bool_t[fMaxParticle+1];

	TA2Apparatus::PostInit();
}

//-----------------------------------------------------------------------------

void TA2TAPS2008::MakeAllRootinos()
{
	fMultipleVetoHit = 0;
	fVeto_index      = 0;
	fVeto_dE         = 0.0;
	fVeto_T          = 0.0;
	fShortGateValue  = 0.0;
	fLongGateValue   = 0.0;
	fTOF             = 0.0;

	// set all array-elements to 0, 0.0 or kFALSE
	for ( Int_t j = 0; j < fMaxParticle; j++)
	{
		fIsVCharged[j]	= kFALSE;
		fdEvE_IsProton[j]	= kFALSE;
		fdEvE_IsChPion[j]	= kFALSE;
		fdEvE_IsElectron[j]	= kFALSE;
		fTOF_IsNeutron[j]	= kFALSE;
		fTOF_IsProton[j]	= kFALSE;
		fTOF_IsChPion[j]	= kFALSE;
		fTOF_IsElectron[j]	= kFALSE;
		fTOF_IsGamma[j]	= kFALSE;
		fPSA_IsNucleon[j]	= kFALSE;
		fPSA_IsNoNucleon[j]	= kFALSE;
		// list will be completed soon
	}

	// first all particles are rootinos!
	for ( Int_t i = 0; i < fNparticle; i++)
	{
		cl = clBaF2[id_clBaF2[i]];
		fp3 = *((clBaF2[id_clBaF2[i]])->GetMeanPosition());
		fp3.SetX(fp3.X()); // + PosShift[1]);
		fp3.SetY(fp3.Y()); // + PosShift[2]);
		fp3.SetZ(fp3.Z()); // + PosShift[3]);
		fParticleID->SetP4(&fP4[i], kRootino,
				(fTapsFudge*(clBaF2[id_clBaF2[i]])->GetEnergy()), &fp3);
		fPDG_ID[i] = kRootino;
		SetParticleInfo(i);
	}
}

//-----------------------------------------------------------------------------

void TA2TAPS2008::MainReconstruct()
{

	fBaF2_Ncluster = fBaF2->GetNCluster();
	fNparticle     = fBaF2_Ncluster;
	NParticles     = fNparticle;
	if(!fNparticle) return;
	id_clBaF2 = fBaF2->GetClustHit();
	clBaF2 = fBaF2->GetCluster();

	if ( fVeto) fVeto_NHits = fVeto->GetNhits();
	else fVeto_NHits = 0;
	if ( fVeto) fVeto_Hits = fVeto->GetHits();
	UInt_t Veto_counter;

	MakeAllRootinos(); // first: all detected particles are made rootinos

	// loop over all BaF2-Clusters
	for ( Int_t i = 0; i < fBaF2_Ncluster; i++)
	{
		fTimeOfFlight[i] = EBufferEnd;
		Veto_counter = 0; fThisVetoFired = 0;
		cl = clBaF2[id_clBaF2[i]];       // i-th cluster
		clBaF2_Nhits = cl->GetNhits();   // # hits of i-th cluster
		clBaF2_elements= cl->GetHits();  // elements (nr. of each crystal)
													//  of i-th cluster
		fMultipleVetoHit = 0;
		fVeto_dE         = 0.0;
		fShortGateValue  = 0.0;
		fLongGateValue   = 0.0;
		fTOF             = 0.0;

		// check if there is a Veto-Hit infront of at least one cluster-crystal
		for ( UInt_t n = 0; n < clBaF2_Nhits; n++)
		{
		  for ( UInt_t j = 0; j < fVeto_NHits; j++)
		  {
		    if ( ( fIsVCharged[i] == kTRUE)
					 && ( clBaF2_elements[n] == (UInt_t) fVeto_Hits[j]))
				 fMultipleVetoHit++;
		    if ( ( fIsVCharged[i] == kFALSE)
					 && ( clBaF2_elements[n] == (UInt_t) fVeto_Hits[j]))
		    {
		      fIsVCharged[i] = kTRUE;
		      Veto_counter++;
		      fThisVetoFired = j;
		      fVeto_dE = fVeto->GetEnergy(fVeto_Hits[j]);
		      fVeto_T = fVeto->GetTime(fVeto_Hits[j]);
		      fVeto_index = fVeto_Hits[j];
		    }
		  }
		}

		// check Veto-dE versus BaF2-E cuts (very preliminary)
		// note: mutiple-veto-hits infront of one baf2-cluster not yet supported
		if(fIsVCharged[i]==kTRUE)
		{
		   //  fTAPS_dEvE->Fill(cl->GetEnergy(), fVeto_dE);
		   if((fTAPS_dEvE_Proton==kTRUE) && fTAPS_dEvE_ProtonCut)
		     if(fTAPS_dEvE_ProtonCut->IsInside(cl->GetEnergy(), fVeto_dE))
		       fdEvE_IsProton[i] = kTRUE;
		   if((fTAPS_dEvE_ChPion==kTRUE) && fTAPS_dEvE_ChPionCut)
		     if(fTAPS_dEvE_ChPionCut->IsInside(cl->GetEnergy(), fVeto_dE))
		       fdEvE_IsChPion[i] = kTRUE;
		   if((fTAPS_dEvE_Electron==kTRUE) && fTAPS_dEvE_ElectronCut)
		     if(fTAPS_dEvE_ElectronCut->IsInside(cl->GetEnergy(), fVeto_dE))
		       fdEvE_IsElectron[i] = kTRUE;
		}

		// check "time of flight" (TAPS-BaF2ClustersTime & CB->MeanTime of Photons)
		// before using this we need time-calibration of TAPS and CB against tagger
		// fTOF = fCB->GetCBMeanTime()- cl->GetTime();
		// don't forget: CB-TimeWalk-Corr is a MUST !!!
		if(gAR->GetProcessType() != EMCProcess)fTOF = fCB->GetCBMeanTime() + cl->GetTime();
		else fTOF = cl->GetTime() - fCB->GetCBMeanTime() ;
		fTimeOfFlight[i] = fTOF;
		//  Double_t tof2 = fCB->GetCBMeanTime()- cl->GetTime(); // just to check

		//  if(fIsVCharged[i]==kTRUE)
		  //  fTAPS_TOF->Fill(fTOF, cl->GetEnergy(), tof2);

		if((fTAPS_TOF_Proton==kTRUE) && fTAPS_TOF_ProtonCut){
		  if(fTAPS_TOF_ProtonCut->IsInside(cl->GetEnergy(),fTOF)){
		    fTOF_IsProton[i] = kTRUE;
		  }
		}
		if((fTAPS_TOF_Neutron==kTRUE) && fTAPS_TOF_NeutronCut)
		  if(fTAPS_TOF_NeutronCut->IsInside(fTOF,cl->GetEnergy()))
		    fTOF_IsNeutron[i] = kTRUE;

		if((fTAPS_TOF_ChPion==kTRUE) && fTAPS_TOF_ChPionCut)
		  if(fTAPS_TOF_ChPionCut->IsInside(fTOF,cl->GetEnergy()))
		    fTOF_IsChPion[i] = kTRUE;

		if((fTAPS_TOF_Electron==kTRUE) && fTAPS_TOF_ElectronCut)
		  if(fTAPS_TOF_ElectronCut->IsInside(fTOF,cl->GetEnergy()))
		    fTOF_IsElectron[i] = kTRUE;

		if((fTAPS_TOF_Gamma==kTRUE) && fTAPS_TOF_GammaCut)
		 if(fTAPS_TOF_GammaCut->IsInside(fTOF,cl->GetEnergy()))
		    fTOF_IsGamma[i] = kTRUE;
		// fNoCBTimeCounter++; not important..now

		if(gAR->GetProcessType() != EMCProcess){
		// check puls shape analysis
		// by now, we will use banana-cuts ... later will be changed (then will test sigmas and so on).
		fShortGateValue = fBaF2->GetSGEnergy(id_clBaF2[i]);
		fLongGateValue  = fBaF2->GetLGEnergy(id_clBaF2[i]);
		}
		//  fTAPS_PSA->Fill(fShortGateValue, fLongGateValue);
		if((fTAPS_PSA_Nucleon==kTRUE) && fTAPS_PSA_NucleonCut)
		  if(fTAPS_PSA_NucleonCut->IsInside(fLongGateValue, fShortGateValue))
		    fPSA_IsNucleon[i] = kTRUE;
		if((fTAPS_PSA_noNucleon==kTRUE) && fTAPS_PSA_noNucleonCut)
		  if(fTAPS_PSA_noNucleonCut->IsInside(fLongGateValue, fShortGateValue))
		    fPSA_IsNoNucleon[i] = kTRUE;
		// determine particle ID based on gathered information
		fPDG_ID[i] = CheckParticleID(i);
		// Save information about this particle
		fParticleID->SetP4(fP4+i, fPDG_ID[i], cl->GetEnergy(), cl->GetMeanPosition());
		fP4tot+=fP4[i];
		fMinv[i]=fP4[i].M();
		// for usage of TA2Particle class -> now write this particle-into partilce-object-array
		SetParticleInfo(i);
	}
	fTimeOfFlight[fBaF2_Ncluster] = EBufferEnd;
}

//-----------------------------------------------------------------------------

Int_t TA2TAPS2008::CheckParticleID(UInt_t i)
{
	fPDG_ID_unclear[i] = kFALSE;

 // 0. Veto - charged/uncharged-Check
 // 1, TOF  - expected to be most reliable -> so its most important
 // 2. dEvE - will see in future how well it works
 // 3. PSA  - not available for inner Rings... PbWO4 don't have ShortGateComponent

 // note: the idea is to support overlapping banana-cuts
 // if TOF is not active (not being used) then (only)dEvE is used for charged particles
 // and (only)PSA is used for uncharged particles

 // if wether TOF, dEvE nor PSA is used , then charged particles will be protons
 // and uncharged particles will be photons

 //particle is charged:

	if(gAR->GetProcessType() == EMCProcess){
		if(fTOF_IsProton[i]==kTRUE) return kProton;
		else return kGamma;
	}

	if(fIsVCharged[i]==kTRUE)
	{
		if(fTAPS_TOF_USEAGE==kTRUE)     // 1. check for proton! most favoured!
		{
		  if(fTOF_IsProton[i]==kTRUE)
		  {
		    fPDG_ID_sec[i] = kProton;
		    if(fdEvE_IsProton[i]==kTRUE) return kProton; // if kTRUE->is clearly proton, else TOF dEvE and PSA sub-checks for secondary particleID
		    if(fTOF_IsElectron[i]== kTRUE)
		    {
		      fPDG_ID_sec[i] = kElectron; fPDG_ID_unclear[i] = kTRUE;
		    }
		    if(fTOF_IsChPion[i]==kTRUE)
		    {
		      fPDG_ID_sec[i] = kPiPlus; fPDG_ID_unclear[i] = kTRUE;
		    }
		    else
		    {
		      if(fdEvE_IsProton[i]==kFALSE && fdEvE_IsChPion[i]==kTRUE)
		      {
		        fPDG_ID_sec[i] = kPiPlus; fPDG_ID_unclear[i] = kTRUE;
		      }
		      if(fdEvE_IsProton[i]==kFALSE && fdEvE_IsElectron[i]==kTRUE)
		      {
		        fPDG_ID_sec[i] = kElectron; fPDG_ID_unclear[i] = kTRUE;
		      }
		      if(fPSA_IsNucleon[i]==kFALSE && fPSA_IsNoNucleon[i]==kTRUE)
		      {
		        fPDG_ID_sec[i] = kPiPlus;fPDG_ID_unclear[i] = kTRUE;
		      }
		    }
		    return kProton;
		  }

		  if(fTOF_IsChPion[i]==kTRUE)  // 2. check for charged pions
		  {
		    fPDG_ID_sec[i] = kPiPlus;
		    if(fdEvE_IsChPion[i]== kTRUE) return kPiPlus;  // if kTRUE->is clearly ch.pion else ..TOF dEvE and PSA sub-checks for secondary particleID
		    if(fTOF_IsElectron[i]== kTRUE)
		    {
		      fPDG_ID_sec[i] = kElectron; fPDG_ID_unclear[i]=kTRUE;
		    }
		    else
		    {
		      if(fdEvE_IsChPion[i]==kFALSE && fdEvE_IsProton[i]==kTRUE)
		      {
		        fPDG_ID_sec[i] = kProton; fPDG_ID_unclear[i] = kTRUE;
		      }
		      if(fdEvE_IsChPion[i]==kFALSE && fdEvE_IsElectron[i]==kTRUE)
		      {
		        fPDG_ID_sec[i] = kElectron; fPDG_ID_unclear[i] = kTRUE;
		      }
		      if(fPSA_IsNoNucleon[i]==kFALSE && fPSA_IsNucleon[i]==kTRUE)
		      {
		        fPDG_ID_sec[i] = kProton;fPDG_ID_unclear[i] = kTRUE;
		      }
		    }
		    return kPiPlus;
		  }
		  if(fTOF_IsElectron[i]== kTRUE)  // @last: check if it's a electron
		  {
		    fPDG_ID_sec[i] = kElectron;
		    if(fdEvE_IsElectron[i]==kTRUE) return kElectron; // if kTRUE->is clearly electron else ... dEvE and PSA sub-checks for secondary particleID:
		    if(fdEvE_IsElectron[i]==kFALSE && fdEvE_IsProton[i]==kTRUE)
		    {
		      fPDG_ID_sec[i]=kProton; fPDG_ID_unclear[i] = kTRUE;
		    }
		    if(fdEvE_IsChPion[i]==kTRUE && fdEvE_IsElectron[i]==kFALSE)
		    {
		      fPDG_ID_sec[i] = kPiPlus; fPDG_ID_unclear[i] = kTRUE;
		    }
		    if(fPSA_IsNoNucleon[i]==kFALSE && fPSA_IsNucleon[i]==kTRUE)
		    {
		      fPDG_ID_sec[i] = kProton; fPDG_ID_unclear[i] = kTRUE;
		    }
		    return kElectron;
		  }
		}

		if(fTAPS_dEvE_USEAGE==kTRUE) // if no TOF usage   but dEvE is used
		{
		  if(fdEvE_IsProton[i]==kTRUE) // 1. check for proton! most favoured!
		  {                        // no PSA sub-check
		    fPDG_ID_sec[i] = kProton;  //dEvE sub-checks for secondary particleID:
		    if(fdEvE_IsElectron[i]==kTRUE)
		    {
		      fPDG_ID_sec[i] = kElectron; fPDG_ID_unclear[i] = kTRUE;
		    }
		    if(fdEvE_IsChPion[i]==kTRUE)
		    {
		      fPDG_ID_sec[i] = kPiPlus; fPDG_ID_unclear[i] = kTRUE;
		    }
		    return kProton;
		  }
		  if(fdEvE_IsChPion[i]==kTRUE) // 2.. check for charged pion
		  {                        // no PSA sub-check
		    fPDG_ID_sec[i] = kPiPlus;  // dEvE sub-check for secondary particleID:
		    if(fdEvE_IsElectron[i]==kTRUE)
		    {
		      fPDG_ID_sec[i] = kElectron; fPDG_ID_unclear[i] = kTRUE;
		    }
		    return kPiPlus;
		  }
		  if(fdEvE_IsElectron[i]==kTRUE)// 3.. check for electron
		  {                        // no PSA sub-check
		    fPDG_ID_sec[i] = kElectron;
		    return kElectron;
		  }
		}

		fPDG_ID_sec[i] = kRootino;
		fPDG_ID_unclear[i] = kTRUE;
		return kProton; // default TAPS-charged.particle return: kProton
	}
	else   //particle is NOT charged
	{
		if(fTAPS_TOF_USEAGE==kTRUE)
		{
		  if(fTOF_IsNeutron[i]==kTRUE) // 1. check for neutron! most favoured!
		  {
		    fPDG_ID_sec[i]= kNeutron;
		    if(fPSA_IsNucleon[i]==kTRUE) return kNeutron; //is clearly neutron else ... check TOF PSA for secondary particlyID
		    if(fTOF_IsChPion[i]==kTRUE || fTOF_IsElectron[i]==kTRUE)
		    {
		      fPDG_ID_sec[i] = kGamma; fPDG_ID_unclear[i] = kTRUE;
		    }
		    if(fPSA_IsNoNucleon[i]==kTRUE)
		    {
		      fPDG_ID_sec[i] = kGamma; fPDG_ID_unclear[i] = kTRUE;
		    }
		    return kNeutron;
		  }
		  if(fTOF_IsGamma[i]==kTRUE)
		  {
		    fPDG_ID_sec[i]=kGamma;
		    //check TOF PSA for secondary particleID
		    if(fTOF_IsChPion[i]==kTRUE || fTOF_IsElectron[i]==kTRUE)
		    {
		      fPDG_ID_sec[i] = kNeutron;fPDG_ID_unclear[i] = kTRUE;
		    }
		    if(fPSA_IsNucleon[i]==kTRUE)
		    {
		      fPDG_ID_sec[i] = kNeutron; fPDG_ID_unclear[i] = kTRUE;
		    }
		  return kGamma;
		  }
		  if(fTOF_IsNeutron[i]==kFALSE && fTOF_IsGamma[i]==kFALSE && fTOF_IsChPion[i]==kFALSE && fTOF_IsElectron[i]==kFALSE)
		  {
		    fPDG_ID_sec[i]= kGamma;
		    fPDG_ID_unclear[i]=kTRUE;
		    return kRootino;
		  }
		}

		if(fTAPS_PSA_USEAGE==kTRUE) // if no TOF usage use PSA if possible:
		{
		  if(fPSA_IsNucleon[i]==kTRUE) // 1. check for neutron! most favoured!
		  {
		    fPDG_ID_sec[i] = kNeutron;
		    fPDG_ID_unclear[i] = kTRUE;// too less information...so unclear is kTRUE
		    return kNeutron;
		  }
		  if(fPSA_IsNoNucleon[i]==kTRUE)// 2. check for NotNeutron->then it is set a gamma
		  {
		    fPDG_ID_sec[i]= kGamma;
		    fPDG_ID_unclear[i] = kTRUE;// too less information...so unclear is kTRUE
		    return kGamma;
		  }
		}
	}

	// if no if-statement is fullfilled (above), then....
	fPDG_ID_sec[i] = kRootino; fPDG_ID_unclear[i] = kTRUE;
	return kGamma;
}

//-----------------------------------------------------------------------------

void TA2TAPS2008::SetConfig(char* line, int key)
{
	TFile *f;
	UInt_t check = 0;
	//UInt_t count = 0;
	fTapsRunStep = 1;

	switch(key)
	{
	 case ETAPSFactor:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPSFactor");
		if(sscanf(line, "%lf", &fTapsFudge) < 1) goto error;
		if(fPrintSetUpInfo==kTRUE) printf("                 .... value is %f", fTapsFudge);
		break;

	 case ETAPS_SimpleReconstruct:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_SimpleReconstruct");
		if(sscanf(line, "%d", &check) < 1) printf("\nERROR\n");
		if(check==0)
		  fSimpleReconstruct=kFALSE;
		else
		  fSimpleReconstruct=kTRUE;
		if(fPrintSetUpInfo==kTRUE) printf("     .... value is %i", fSimpleReconstruct);
		break;

	 case ETAPS_TOF_USEAGE:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_TOF_USEAGE");
		sscanf(line, "%d", &check);
		if(check==0)
		  fTAPS_TOF_USEAGE = kFALSE;
		else
		  fTAPS_TOF_USEAGE = kTRUE;
		if(fPrintSetUpInfo==kTRUE) printf("            .... value is %d ",fTAPS_TOF_USEAGE );
		break;

	 case ETAPS_TOF_Cuts:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_TOF_Cuts ");
		if(fTAPS_TOF_USEAGE==kTRUE)
		{
		  sscanf( line, "%s", fTAPS_TOF_Cuts );
		  //strcat(fTAPS_TOF_Cuts, ".root");
		  //fTestfile[0]->Open(fTAPS_TOF_Cuts,"rb");
		  f = new TFile(fTAPS_TOF_Cuts,"read" );
		  //if(fTestfile[0]){fTestfile[0]->Close();}
		  if(f)
		  {
		   f->Close(); delete f;
		  }
		  else
		  {
		    printf("\n!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		    printf("!       File containing TOF-CUTs is missing.        !\n");
		    printf("!            Turning TAPS TOF off                   !\n");
		    printf("!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
		    fTAPS_TOF_USEAGE = kFALSE; delete f;
		  }
		}
		else
		  strcpy(fTAPS_TOF_Cuts, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("             .... value is ");printf(fTAPS_TOF_Cuts); }
	 break;

	 case ETAPS_TOF_Neutron:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_TOF_Neutron");
		if(fTAPS_TOF_USEAGE==kTRUE)
		{
		  sscanf(line, "%d", &check);
		  if(check==0)
		   fTAPS_TOF_Neutron = kFALSE;
		  else
		   fTAPS_TOF_Neutron=kTRUE;
		}
		if(fPrintSetUpInfo==kTRUE) printf("           .... value is %d", fTAPS_TOF_Neutron);
		break;

	 case ETAPS_TOF_Neutron_CutName:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_TOF_Neutron_CutName");
		if((fTAPS_TOF_USEAGE==kTRUE) && (fTAPS_TOF_Neutron==kTRUE))
		{
		  sscanf(line, "%s", fTAPS_TOF_Neutron_CutName);
		  fTAPS_TOF_CutFile = new TFile(fTAPS_TOF_Cuts, "READ");
		  fTAPS_TOF_NeutronCut = (TCutG*)fTAPS_TOF_CutFile->Get(fTAPS_TOF_Neutron_CutName);
		  fTAPS_TOF_CutFile->Close();
		}
		else
		  strcpy(fTAPS_TOF_Neutron_CutName, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("   .... value is "); printf(fTAPS_TOF_Neutron_CutName); }
		break;

	 case ETAPS_TOF_Proton:
		if(fPrintSetUpInfo==kTRUE){ printf("\n ... reading ETAPS_TOF_Proton"); }
		if(fTAPS_TOF_USEAGE==kTRUE)
		{
		  sscanf( line, "%d", &check );
		  if(check==0)
		    fTAPS_TOF_Proton=kFALSE;
		  else
		    fTAPS_TOF_Proton=kTRUE;
		}
		if(fPrintSetUpInfo==kTRUE) printf("            .... value is %d", fTAPS_TOF_Proton);
		break;

	 case ETAPS_TOF_Proton_CutName:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_TOF_Proton_CutName");
		if((fTAPS_TOF_USEAGE==kTRUE) && (fTAPS_TOF_Proton==kTRUE))
		{
		  sscanf(line, "%s", fTAPS_TOF_Proton_CutName);
		  fTAPS_TOF_CutFile = new TFile(fTAPS_TOF_Cuts,"READ");
		  fTAPS_TOF_ProtonCut = (TCutG*)fTAPS_TOF_CutFile->Get(fTAPS_TOF_Proton_CutName);
		  fTAPS_TOF_CutFile->Close();
		}
		else
		  strcpy(fTAPS_TOF_Proton_CutName, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("    .... value is "); printf(fTAPS_TOF_Proton_CutName); }
		break;

	 case ETAPS_TOF_ChPion:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_TOF_ChPion");
		if(fTAPS_TOF_USEAGE==kTRUE)
		{
		  sscanf(line, "%d", &check);
		  if(check==0)
		    fTAPS_TOF_ChPion = kFALSE;
		  else
		    fTAPS_TOF_ChPion = kTRUE;
		}
		if(fPrintSetUpInfo==kTRUE) printf("            .... value is %d", fTAPS_TOF_ChPion);
		break;

	 case ETAPS_TOF_ChPion_CutName:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_TOF_ChPion_CutName");
		if((fTAPS_TOF_USEAGE)==kTRUE && (fTAPS_TOF_ChPion==kTRUE))
		{
		  sscanf(line, "%s", fTAPS_TOF_ChPion_CutName);
		  fTAPS_TOF_CutFile = new TFile(fTAPS_TOF_Cuts,"READ");
		  fTAPS_TOF_ChPionCut = (TCutG*)fTAPS_TOF_CutFile->Get(fTAPS_TOF_ChPion_CutName);
		  fTAPS_TOF_CutFile->Close();
		}
		else
		  strcpy(fTAPS_TOF_ChPion_CutName, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("    .... value is "); printf(fTAPS_TOF_ChPion_CutName); }
	 break;

	 case ETAPS_TOF_Electron:
		if(fPrintSetUpInfo==kTRUE){ printf("\n ... reading ETAPS_TOF_Electron"); }
		if(fTAPS_TOF_USEAGE==kTRUE)
		{
		  sscanf(line, "%d", &check);
		  if(check==0)
		    fTAPS_TOF_Electron=kFALSE;
		  else
		   fTAPS_TOF_Electron=kTRUE;
		}
		if(fPrintSetUpInfo==kTRUE){ printf("          .... value is %d", fTAPS_TOF_Electron); }
	 break;

	 case ETAPS_TOF_Electron_CutName:
		if(fPrintSetUpInfo==kTRUE){printf("\n ... reading ETAPS_TOF_Electron_CutName");}
		if((fTAPS_TOF_USEAGE==kTRUE) && (fTAPS_TOF_Electron==kTRUE))
		{
		  sscanf(line, "%s", fTAPS_TOF_Electron_CutName);
		  fTAPS_TOF_CutFile = new TFile(fTAPS_TOF_Cuts,"READ");
		  fTAPS_TOF_ElectronCut = (TCutG*)fTAPS_TOF_CutFile->Get(fTAPS_TOF_Electron_CutName);
		  fTAPS_TOF_CutFile->Close();
		}
		else
		  strcpy(fTAPS_TOF_Electron_CutName, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("  .... value is "); printf(fTAPS_TOF_Electron_CutName); }
	 break;

	 case ETAPS_TOF_Gamma:
		if(fPrintSetUpInfo==kTRUE){printf("\n ... reading ETAPS_TOF_Gamma");}
		if(fTAPS_TOF_USEAGE==kTRUE)
		{
		  sscanf(line, "%d", &check);
		  if(check==0)
		    fTAPS_TOF_Gamma=kFALSE;
		  else
		    fTAPS_TOF_Gamma=kTRUE;
		}
		if(fPrintSetUpInfo==kTRUE) printf("             .... value is %d", fTAPS_TOF_Gamma);
	 break;

	 case ETAPS_TOF_Gamma_CutName:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_TOF_Gamma_CutName");
		if((fTAPS_TOF_USEAGE==kTRUE) && (fTAPS_TOF_Gamma==kTRUE))
		{
		  sscanf(line, "%s", fTAPS_TOF_Gamma_CutName);
		  fTAPS_TOF_CutFile = new TFile(fTAPS_TOF_Cuts, "READ");
		  fTAPS_TOF_GammaCut = (TCutG*)fTAPS_TOF_CutFile->Get(fTAPS_TOF_Gamma_CutName);
		  fTAPS_TOF_CutFile->Close();
		}
		else
		  strcpy(fTAPS_TOF_Gamma_CutName, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("     .... value is "); printf(fTAPS_TOF_Gamma_CutName); }
		break;

	 case ETAPS_dEvE_USEAGE:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_dEvE_USEAGE");
		sscanf(line, "%d", &check);
		if(check==0)
		  fTAPS_dEvE_USEAGE = kFALSE;
		else
		  fTAPS_dEvE_USEAGE = kTRUE;
		if(fPrintSetUpInfo==kTRUE) printf("           .... value is %d", fTAPS_dEvE_USEAGE);
		break;

	 case ETAPS_dEvE_Cuts:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_dEvE_Cuts");
		if(fTAPS_dEvE_USEAGE==kTRUE)
		{
		  sscanf(line, "%s", fTAPS_dEvE_Cuts);
		  //strcat(fTAPS_dEvE_Cuts, ".root");
		  f = new TFile(fTAPS_dEvE_Cuts,"READ");
		  if(f)
		  {
		    f->Close(); delete f;
		  }
		  else
		  {
		    printf("\n!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		    printf("!       File containing dEvE-CUTs is missing.       !\n");
		    printf("!            Turning TAPS dEvEoff                   !\n");
		    printf("!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
		    fTAPS_dEvE_USEAGE = kFALSE; delete f;
		  }
		}
		else
		  strcpy(fTAPS_dEvE_Cuts, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("             .... value is "); printf(fTAPS_dEvE_Cuts); }
		break;

	 case ETAPS_dEvE_Proton:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_dEvE_Proton");
		if(fTAPS_dEvE_USEAGE==kTRUE)
		{
		  sscanf(line, "%d", &check);
		  if(check==0)
		    fTAPS_dEvE_Proton = kFALSE;
		  else
		    fTAPS_dEvE_Proton = kTRUE;
		}
		if(fPrintSetUpInfo==kTRUE) printf("           .... value is %d", fTAPS_dEvE_Proton);
		break;

	 case ETAPS_dEvE_Proton_CutName:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_dEvE_Proton_CutName");
		if(fTAPS_dEvE_USEAGE==kTRUE && fTAPS_dEvE_Proton == kTRUE )
		{
		  sscanf(line, "%s", fTAPS_dEvE_Proton_CutName);
		  fTAPS_dEvE_CutFile = new TFile(fTAPS_dEvE_Cuts, "READ");
		  fTAPS_dEvE_ProtonCut  = (TCutG*)fTAPS_dEvE_CutFile->Get(fTAPS_dEvE_Proton_CutName);
		  fTAPS_dEvE_CutFile->Close();
		}
		else
		  strcpy(fTAPS_dEvE_Proton_CutName, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("   .... value is "); printf(fTAPS_dEvE_Proton_CutName); }
	 break;

	 case ETAPS_dEvE_ChPion:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_dEvE_ChPion");
		if(fTAPS_dEvE_USEAGE==kTRUE)
		{
		  sscanf( line, "%d", &check );
		  if(check==0)
		    fTAPS_dEvE_ChPion = kFALSE;
		  else
		    fTAPS_dEvE_ChPion = kTRUE;
		}
		if(fPrintSetUpInfo==kTRUE) printf("           .... value is %d", fTAPS_dEvE_ChPion);
		break;

	 case ETAPS_dEvE_ChPion_CutName:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_dEvE_ChPion_CutName");
		if(fTAPS_dEvE_USEAGE==kTRUE && fTAPS_dEvE_ChPion == kTRUE)
		{
		  sscanf(line, "%s", fTAPS_dEvE_ChPion_CutName);
		  fTAPS_dEvE_CutFile = new TFile(fTAPS_dEvE_Cuts, "READ");
		  fTAPS_dEvE_ChPionCut  = (TCutG*)fTAPS_dEvE_CutFile->Get(fTAPS_dEvE_ChPion_CutName);
		  fTAPS_dEvE_CutFile->Close();
		}
		else
		  strcpy(fTAPS_dEvE_ChPion_CutName, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("   .... value is "); printf(fTAPS_dEvE_ChPion_CutName); }
		break;

	 case ETAPS_dEvE_Electron:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_dEvE_Electron");
		if(fTAPS_dEvE_USEAGE==kTRUE)
		{
		sscanf(line, "%d", &check);
		if(check==0)
		  fTAPS_dEvE_Electron = kFALSE;
		else
		  fTAPS_dEvE_Electron = kTRUE;
		}
		if(fPrintSetUpInfo==kTRUE) printf("         .... value is %d", fTAPS_dEvE_Electron );
		break;

	 case ETAPS_dEvE_Electron_CutName:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_dEvE_Electron_CutName");
		if(fTAPS_dEvE_USEAGE==kTRUE && fTAPS_dEvE_Electron == kTRUE)
		{
		  sscanf(line, "%s", fTAPS_dEvE_Electron_CutName);
		  fTAPS_dEvE_CutFile = new TFile(fTAPS_dEvE_Cuts, "READ");
		  fTAPS_dEvE_ElectronCut = (TCutG*)fTAPS_dEvE_CutFile->Get(fTAPS_dEvE_Electron_CutName);
		  fTAPS_dEvE_CutFile->Close();
		}
		else
		  strcpy(fTAPS_dEvE_Electron_CutName, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf(" .... value is "); printf(fTAPS_dEvE_Electron_CutName); }
		break;

	 case ETAPS_PSA_USEAGE:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_PSA_USEAGE");
		sscanf(line, "%d", &check);
		if(check==0)
		  fTAPS_PSA_USEAGE = kFALSE;
		else
		  fTAPS_PSA_USEAGE = kTRUE;
		if(fPrintSetUpInfo==kTRUE) printf("            .... value is %d", fTAPS_PSA_USEAGE);
		break;

	 case ETAPS_PSA_Cuts:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_PSA_Cuts");
		if(fTAPS_PSA_USEAGE==kTRUE)
		{
		  sscanf(line, "%s", fTAPS_PSA_Cuts);
		  //strcat(fTAPS_PSA_Cuts, ".root");
		  f = new TFile(fTAPS_PSA_Cuts,"READ");
		  if(f)
		  {
		    f->Close(); delete f;
		  }
		  else
		  {
		    printf("\n!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		    printf("!       File containing PSA-CUTs is missing.        !\n");
		    printf("!            Turning TAPS PSA off                   !\n");
		    printf("!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
		    fTAPS_PSA_USEAGE = kFALSE; delete f;
		  }
		}
		else
		  strcpy(fTAPS_PSA_Cuts, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("              .... value is "); printf(fTAPS_PSA_Cuts); }
	 break;

	 case ETAPS_PSA_Nucleon:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_PSA_Nucleon");
		if(fTAPS_PSA_USEAGE==kTRUE)
		{
		  sscanf(line, "%d", &check);
		  if(check==0)
		    fTAPS_PSA_Nucleon = kFALSE;
		  else
		    fTAPS_PSA_Nucleon= kTRUE;
		}
		if(fPrintSetUpInfo==kTRUE) printf("           .... value is %d", fTAPS_PSA_Nucleon);
		break;

	 case ETAPS_PSA_Nucleon_CutName:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_PSA_Nucleon_CutName");
		if((fTAPS_PSA_USEAGE==kTRUE) && (fTAPS_PSA_Nucleon== kTRUE))
		{
		  sscanf(line, "%s", fTAPS_PSA_Nucleon_CutName);
		  fTAPS_PSA_CutFile = new TFile(fTAPS_PSA_Cuts,"READ");
		  fTAPS_PSA_NucleonCut = (TCutG*)fTAPS_PSA_CutFile->Get(fTAPS_PSA_Nucleon_CutName);
		  fTAPS_PSA_CutFile->Close();
		}
		else
		  strcpy(fTAPS_PSA_Nucleon_CutName, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE){ printf("   .... value is "); printf(fTAPS_PSA_Nucleon_CutName); }
		break;

	 case ETAPS_PSA_noNucleon:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_PSA_noNucleon");
		if(fTAPS_PSA_USEAGE==kTRUE)
		{
		  sscanf(line, "%d", &check);
		  if(check == 0)
		    fTAPS_PSA_noNucleon = kFALSE;
		  else
		    fTAPS_PSA_noNucleon= kTRUE;
		}
		if(fPrintSetUpInfo==kTRUE) printf("         .... value is %d", fTAPS_PSA_noNucleon);
		break;

	 case ETAPS_PSA_noNucleon_CutName:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_PSA_Nucleon_CutName");
		if((fTAPS_PSA_USEAGE==kTRUE) && (fTAPS_PSA_noNucleon==kTRUE))
		{
		  sscanf(line, "%s", fTAPS_PSA_noNucleon_CutName);
		  fTAPS_PSA_CutFile = new TFile(fTAPS_PSA_Cuts,"READ");
		  fTAPS_PSA_noNucleonCut = (TCutG*)fTAPS_PSA_CutFile->Get(fTAPS_PSA_noNucleon_CutName);
		  fTAPS_PSA_CutFile->Close();
		}
		else
		  strcpy(fTAPS_PSA_noNucleon_CutName, " not 'activated'");
		if(fPrintSetUpInfo==kTRUE) { printf("   .... value is "); printf(fTAPS_PSA_noNucleon_CutName); }
		break;
	 // case ETAPS_OUTPUT_FILE:
	 //  if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_OUTPUT_FILE ");
	 //  sscanf(line, "%s", fTAPS_OUTPUT_FILE_Name);
	 //  //strcat(fTAPS_OUTPUT_FILE_Name, ".root");
	 //  f = new TFile(fTAPS_OUTPUT_FILE_Name, "RECREATE");
	 //  if(f)
	 //  {
	 //    f->Close(); delete f;
	 //    fTAPS_OUTPUT_FILE = new TFile(fTAPS_OUTPUT_FILE_Name, "RECREATE");
	 //    fTAPS_OUTPUT_FILE_Exists = kTRUE;
	 //    if(fPrintSetUpInfo==kTRUE){ printf("          .... value is "); printf(fTAPS_OUTPUT_FILE_Name); }
	 //  }
	 //  else
	 //  {
	 //    strcpy(fTAPS_OUTPUT_FILE_Name, " not 'activated'"); delete f;
	 //    if(fPrintSetUpInfo==kTRUE){ printf("          .... value is "); printf(fTAPS_OUTPUT_FILE_Name); }
	 //    printf("\n!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	 //    printf("!          TAPS_OUTPUT_File is missing.             !\n");
	 //    printf("!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
	 //  }
	 // break;
	 case ETAPS_fReadPbWO4:
		if(fPrintSetUpInfo==kTRUE) printf("\n ... reading ETAPS_fReadPbWO4");
		sscanf(line, "%d", &check);
		if(check==0)
		  fTAPS_ReadPbWO4 = kFALSE;
		else
		  fTAPS_ReadPbWO4 = kTRUE;
		if(fPrintSetUpInfo==kTRUE) printf("            .... value is %d", fTAPS_ReadPbWO4);
		break;

	 default:
		 // Command not found...possible pass to next config
		 TA2Apparatus::SetConfig(line, key);
		 break;
	}
	return;
 error: PrintError(line);
	return;
}

//-----------------------------------------------------------------------------

