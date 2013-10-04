////////////////////////////////////////////////////////////////////////////
//                                                                        //
//     --Author		H. Berghaeuser June 2008                             //
//                      HenningBerghaeuser@web.de                         //
//     -- Update        H. Berghaeuser Okt 2008                           //
//     -- Update        H. Berghaeuser Dez 2008                           //
//     -- Update        H. Berghaeuser FEB 2009                           //
//                                                                        //
//     -- Use the TAPS2009.dat !                                          //
//                                                                        //
//            TA2Taps2009.cc       version 1.0                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TA2TAPS2009.h"
#include "TA2Analysis.h"

enum { ETAPS_BaF2=6000, ETAPS_Veto };
enum { ETAPSFactor=6002, ETAPS_SimpleReconstruct, ETAPS_dEvE_Cuts,
       ETAPS_dEvE_Proton, ETAPS_dEvE_Proton_CutName,
       ETAPS_dEvE_ChPion, ETAPS_dEvE_ChPion_CutName,
       ETAPS_dEvE_Electron,ETAPS_dEvE_Electron_CutName };

static const Map_t kTAPS2009Keys[] = {
  {"TAPS-Factor:",                ETAPSFactor},
  {"TAPS_SimpleReconstruct:",     ETAPS_SimpleReconstruct},
  {"TAPS_dEvE_Cuts:",             ETAPS_dEvE_Cuts},
  {"TAPS_dEvE_Proton:",           ETAPS_dEvE_Proton},
  {"TAPS_dEvE_Proton_CutName:",   ETAPS_dEvE_Proton_CutName},
  {"TAPS_dEvE_ChPion:",           ETAPS_dEvE_ChPion},
  {"TAPS_dEvE_ChPion_CutName:",   ETAPS_dEvE_ChPion_CutName},
  {"TAPS_dEvE_Electron:",         ETAPS_dEvE_Electron},
  {"TAPS_dEvE_Electron_CutName:", ETAPS_dEvE_Electron_CutName},
  {NULL,            -1}
};

static Map_t kValidTAPS2009Detectors[] = {
  {"TA2TAPS_BaF2",      ETAPS_BaF2},
  {"TA2PlasticPID",     ETAPS_Veto},
  {NULL, 		-1}
};

ClassImp(TA2TAPS2009)

//-----------------------------------------------------------------------------

TA2TAPS2009::TA2TAPS2009(const char* name, TA2System* analysis)
            :TA2Apparatus(name, analysis, kValidTAPS2009Detectors)
{
  fTapsRunStep        = 0;
  fCalibVetoStep      = 0;
  fCalibBaF2Step      = 0;
  fBaF2               = NULL;
  fVeto               = NULL;
  fDeltaE             = NULL;
  fEcharged           = NULL;
  fSimpleReconstruct  = kTRUE;
  fTAPS_dEvE_Proton   = kFALSE;
  fTAPS_dEvE_ChPion   = kFALSE;
  fTAPS_dEvE_Electron = kFALSE;
  fTimeShift          = 0.0; // maybe will later use SvensTimeShift
  fTapsFudge          = 1.0;
  AddCmdList( kTAPS2009Keys );

//  printf("\n-- TA2Taps2009 .......... done\n");
}

//-----------------------------------------------------------------------------

TA2TAPS2009::~TA2TAPS2009()
{
  if (fParticleID) delete fParticleID;
  if( fDeltaE ) delete fDeltaE;
  if( fDeltaE ) delete fEcharged;
  // Free up allocated memory
}

//-----------------------------------------------------------------------------
TA2DataManager*  TA2TAPS2009::CreateChild(const char* name, Int_t dclass)
{

  if(!name) name = Map2String(dclass);
  switch (dclass)
  {
  case ETAPS_BaF2:
    fBaF2 = new TA2TAPS_BaF2(name, this);
    return fBaF2;
  case ETAPS_Veto:
    fVeto = new TA2PlasticPID(name, this);
    return fVeto;
  default:
    return NULL;
  }
}

//-----------------------------------------------------------------------------

void TA2TAPS2009::LoadVariable()
{

  //                            name             pointer             type-spec
  TA2DataManager::LoadVariable("DeltaE",         fDeltaE,           EDSingleX);
  TA2DataManager::LoadVariable("Echarged",       fEcharged,         EDSingleX);
  TA2DataManager::LoadVariable("Echargedcl",     fEchargedcl,       EDSingleX);
  TA2DataManager::LoadVariable("EchargedIndex",  fEchargedIndex ,   EDSingleX);
  //
  TA2Apparatus::LoadVariable();
}

//-----------------------------------------------------------------------------

void TA2TAPS2009::PostInit()
{
  if(!fParticleID)
    PrintError("", "<Configuration aborted: ParticleID class MUST be specified>", EErrFatal);

  fMaxParticle          = fBaF2->GetMaxCluster();
  fDeltaE               = new Double_t[fBaF2->GetNelement()+1];
  fEcharged             = new Double_t[fBaF2->GetNelement()+1];
  fEchargedcl           = new Double_t[fBaF2->GetNelement()+1];
  fEchargedIndex        = new Double_t[fBaF2->GetNelement()+1];
  fP4_Nphoton           = new TLorentzVector[fMaxParticle+1];
  fM_Nphoton            = new Double_t[fMaxParticle+1];
  particles             = new TA2Particle[fMaxParticle+1];
  fIsVCharged           = new Bool_t[fMaxParticle+1];
  fdEvE_IsProton        = new Bool_t[fMaxParticle+1];
  fdEvE_IsChPion        = new Bool_t[fMaxParticle+1];
  fdEvE_IsElectron      = new Bool_t[fMaxParticle+1];
  fTapsRunStep          = 1;

  TA2Apparatus::PostInit();
}

//-----------------------------------------------------------------------------

void TA2TAPS2009::MakeAllRootinos()
{
  fMultipleVetoHit      = 0;
  fVeto_dE              = 0.0;
  fVeto_T               = 0.0;
  fVeto_index           = 0;
  fEnergy_BaF2          = 0.0;
  fEnergy_Index         = 0.0;

  // set all array-elements to 0, 0.0 or kFALSE
  for(Int_t j=0; j<fMaxParticle; j++)
  {
    fIsVCharged[j]      = kFALSE;
    fdEvE_IsProton[j]   = kFALSE;
    fdEvE_IsChPion[j]   = kFALSE;
    fdEvE_IsElectron[j] = kFALSE;
  }

  // first all particles are rootinos!
  for(Int_t i=0; i<fNparticle; i++)
  {
    cl = clBaF2[id_clBaF2[i]];
    fp3 = *((clBaF2[id_clBaF2[i]])->GetMeanPosition());
    fp3.SetX(fp3.X()); // + PosShift[1]);
    fp3.SetY(fp3.Y()); // + PosShift[2]);
    fp3.SetZ(fp3.Z()); // + PosShift[3]);
    fParticleID->SetP4(&fP4[i], kRootino, (fTapsFudge*(clBaF2[id_clBaF2[i]])->GetEnergy()), &fp3);
    fPDG_ID[i] = kRootino;

    fDeltaE[i]   = (Double_t)ENullHit;
    fEcharged[i] = (Double_t)ENullHit;
  }
}

//-----------------------------------------------------------------------------

void TA2TAPS2009::MainReconstruct()
{
  fBaF2_Ncluster = fBaF2->GetNCluster();
  fNparticle     = fBaF2_Ncluster;
  if(!fNparticle) return;

  id_clBaF2 = fBaF2->GetClustHit();
  clBaF2 = fBaF2->GetCluster();

  fVeto_NHits = fVeto->GetNhits();
  fVeto_Hits  = fVeto->GetHits();

  MakeAllRootinos();  // first: all detected particles are made rootinos
  Int_t clindex = 0;

  // loop over all BaF2-Clusters
  for(UInt_t i=0; i<(UInt_t)fBaF2_Ncluster; i++)
  {
    cl = clBaF2[id_clBaF2[i]];         // i-th cluster
    clBaF2_Nhits = cl->GetNhits();     // # hits of i-th cluster
    clBaF2_elements= cl->GetHits();    // elements (nr. of each crystal) of i-th cluster
    clindex = (UInt_t) cl->GetIndex(); // get index-crystal of actual cluster
    fVeto_dE    = 0.0;
    fVeto_T     = 0.0;
    fVeto_index = 0;

    // check if there is a Veto-Hit infront of at least one cluster-crystal
    for(UInt_t j=0; j<fVeto_NHits; j++)
    {
      if(clindex==fVeto_Hits[j]) // first check veto infront of BaF2-cluster-index-Crystal
      {
        fIsVCharged[i] = kTRUE;
        fVeto_dE       = fVeto->GetEnergy(fVeto_Hits[j]);
        fVeto_T        = fVeto->GetTime(fVeto_Hits[j]);
        fVeto_index    = fVeto_Hits[j];
        fEnergy_BaF2   = fBaF2->GetEnergy(clindex);
        fEnergy_Index  = fEnergy_BaF2;
      }
      else // .. if not, then check other vetos infront of the complete BaF2cluster
      {
        for(UInt_t n=0; n<clBaF2_Nhits; n++)
        {
          if(fIsVCharged[i]==kFALSE && (Int_t)clBaF2_elements[n]==fVeto_Hits[j])
          {
            fIsVCharged[i] = kTRUE;
            fVeto_dE       = fVeto->GetEnergy(fVeto_Hits[j]);
            fVeto_T        = fVeto->GetTime(fVeto_Hits[j]);
            fVeto_index    = fVeto_Hits[j];
            fEnergy_BaF2   = cl->GetEnergy();
          }
        }
      }
    }

    // check Veto-dE versus BaF2-E Banana-cuts -> dEvE@TAPS
    // note: mutiple-veto-hits infront of one baf2-cluster not yet supported
    if(fIsVCharged[i]==kTRUE)
    {
      fDeltaE[i] = fVeto_dE;
      fEchargedcl[i] = cl->GetEnergy();
      fEcharged[i] = fEnergy_BaF2;
      fEchargedIndex[i] = fEnergy_Index;
      if ( (fTAPS_dEvE_Proton==kTRUE) && (fTAPS_dEvE_ProtonCut->IsInside(fEnergy_BaF2, fVeto_dE)))
        fdEvE_IsProton[i] = kTRUE;
      if ( (fTAPS_dEvE_ChPion==kTRUE) && (fTAPS_dEvE_ChPionCut->IsInside(fEnergy_BaF2, fVeto_dE)))
        fdEvE_IsChPion[i] = kTRUE;
      if ( (fTAPS_dEvE_Electron==kTRUE) && (fTAPS_dEvE_ElectronCut->IsInside(fEnergy_BaF2, fVeto_dE)))
        fdEvE_IsElectron[i] = kTRUE;
     }

    // determine particle ID based on gathered information
    fPDG_ID[i] = CheckParticleID(i);

    SetParticleInfo(i);  // Set TA2Particle

    // Save information about this particle
    fParticleID->SetP4( fP4+i,fPDG_ID[i],cl->GetEnergy(),cl->GetMeanPosition());
    fP4tot+=fP4[i];
    fMinv[i] = fP4[i].M();
  }
}

//-----------------------------------------------------------------------------

Int_t TA2TAPS2009::CheckParticleID(UInt_t i)
{
  // 1. Veto:  no veto -> photon
  // 2. dEvE:  test (Banana)TCuts

  if(fIsVCharged[i]== kFALSE)
    return kGamma;
  else
  {
    if(fdEvE_IsProton[i]==kTRUE) return kProton;
    if(fdEvE_IsChPion[i]== kTRUE) return kPiPlus;
    if(fdEvE_IsElectron[i]==kTRUE) return kElectron;
    return kRootino;
  }
}

//--------------------------------------------------------------------------------

void TA2TAPS2009::CalibrateBaF2Energy(TLorentzVector *photonCB)
{
  if(fCalibBaF2Step==0)
  {
    calibTAPS_m1g_AllCh  = new TH2F("Calib_TAPS_1g_IM_AllCh", "Calib_TAPS_1g_IM_AllCh", 2000, 0, 1000, 385, 0, 385);
    calibTAPS_m1g_Single = new TH1F("Calib_TAPS_1g_IMS", "Calib_TAPS_1g_IMS", 700 ,0, 700);
  }

  Int_t* VetoHits= fVeto->GetHits();
  HitCluster_t** clBaF2; clBaF2 = fBaF2->GetCluster(); // get cluster structs
  HitCluster_t* cl;                                    // cluster struct
  UInt_t* id_clBaF2; id_clBaF2 = fBaF2->GetClustHit();
  UInt_t fBaF2_Ncluster = fBaF2->GetNCluster();
  UInt_t clBaF2_Nhits;                                 // crystal-hits inside a cluster
  UInt_t* clBaF2_elements;
  Bool_t IsVCharged[fBaF2_Ncluster];
  for(UInt_t i=0; i<fBaF2_Ncluster; i++) IsVCharged[i]=kFALSE;
  UInt_t CountUnchargedClusters = 0;
  UInt_t UnchargedClusters[fBaF2_Ncluster];

  // first: find neutral hits !!!  require: veto must NOT have fired in front of any cluster-crystal!
  // loop over all BaF2-clusters
  if(fVeto->GetNhits() < fBaF2_Ncluster)
  {
    for(UInt_t i=0; i< fBaF2_Ncluster; i++)
    {
      Bool_t IsVCharged = kFALSE;
      cl = clBaF2[id_clBaF2[i]];        // i-th cluster
      clBaF2_Nhits = cl->GetNhits();    // # hits of i-th cluster
      clBaF2_elements= cl->GetHits();   // elements (nr. of each crystal) of i-th cluster

      // check if there is a Veto-Hit infront of at least one cluster-crystal
      for(UInt_t n=0; n<clBaF2_Nhits; n++)
      {
        for(UInt_t j= 0; j<fVeto->GetNhits(); j++)
          if((Int_t)clBaF2_elements[n]==VetoHits[j])
            IsVCharged = kTRUE;
      }
      if(IsVCharged==kFALSE)
      {
        UnchargedClusters[CountUnchargedClusters] = i;
        CountUnchargedClusters++;
      }
    }
  }

  // BaF2 - Energy-Calib (pi0-method)
  if(CountUnchargedClusters==1)
  {
    HitCluster_t* cl0 = clBaF2[id_clBaF2[UnchargedClusters[0]]];
    UInt_t index = (UInt_t)cl0->GetIndex();
    TLorentzVector photonTAPS, pion;
    TVector3  photonTAPS_TV3;
    photonTAPS_TV3 = (cl0->GetMeanPosition())->Unit() * (cl0->GetEnergy()) ;
    photonTAPS.SetE( cl0->GetEnergy() );
    photonTAPS.SetVect( photonTAPS_TV3 );
    pion = *photonCB + photonTAPS;
    //printf("\nTAPS_E = %f %f   CB_E = %f  Pion.M= %f",cl0->GetEnergy(),photonTAPS.E(), photon[0].E(),pion.M());
    calibTAPS_m1g_Single->Fill(pion.M());
    calibTAPS_m1g_AllCh->Fill(pion.M(), index);
  }

  fCalibBaF2Step  = 1;
}

// -------------------------------------------------------------------------------

void TA2TAPS2009::CalibrateVetoEnergy()
{
  if(fCalibVetoStep==0)
  {
    VetoData = new TNtuple("TapsData",
                           "Ntuple containing correlated BaF2 and Veto Hits",
                           "module:vetoRawEnergy:baf2Energy:baf2Time:vetoEnergy:vetoTime:baf2ClEnergy:baf2ClTime");
    calibTAPS_VetoCorr = new TH2F("Calib_TAPS_VetoCorrel",  "Calib_TAPS_VetoCorrelXXX" ,384, 0 ,384, 384 ,0 ,384);
    calib_TAPSdEvEcl = new TH2F("Calib_TAPS_dEvEcl","Calib_TAPS_dEvEcl", 1000, 0, 1000, 100, 0, 10);
    calib_TAPSdEvE   = new TH2F("Calib_TAPS_dEvE","Calib_TAPS_dEvE",500, 0, 500, 100, 0, 10);
  }
  // ------  TAPS BaF2 Veto correlation --- Start
  Int_t* VetoHits= fVeto->GetHits();
  Int_t* BaF2Hits = fBaF2->GetHits();

  for (UInt_t i=0; i<fBaF2->GetNhits(); i++)
  {
    // loop over veto hits
    for (UInt_t j=0; j<fVeto->GetNhits(); j++)
         calibTAPS_VetoCorr->Fill(BaF2Hits[i], VetoHits[j]);
  }
  // ------  TAPS BaF2 Veto correlation --- End

  Int_t nElem = fBaF2->GetNelement();
  Bool_t baf2Hit[nElem];
  Bool_t vetoHit[nElem];
  for(Int_t t=0; t<nElem; t++)
  {
    baf2Hit[nElem] = kFALSE;
    vetoHit[nElem] = kFALSE;
  }

  for (UInt_t i=0; i<fBaF2->GetNhits(); ++i)
    baf2Hit[fBaF2->GetHits(i)] = 1;

  for (UInt_t i=0; i<fVeto->GetNhits(); ++i)
    vetoHit[fVeto->GetHits(i)] = 1;

  for (UInt_t module=0; module<fBaF2->GetNelement(); ++module)
  {
    if (baf2Hit[module] && vetoHit[module])
    {
       //baf2RawEnergy = fBaF2->GetElement(module)->GetRawADCValue();
       //baf2RawTime   = fBaF2->GetElement(module)->GetRawTDCValue();
       Float_t vetoRawEnergy = (Float_t)fVeto->GetElement(module)->GetRawADCValue();
       //vetoRawTime   = fVeto->GetElement(module)->GetRawTDCValue();

       Float_t baf2Energy   = (Float_t)fBaF2->GetEnergy(module);
       Float_t baf2ClEnergy = TAPS_GetCLInfo(module, 0);
       Float_t baf2ClTime   = TAPS_GetCLInfo(module, 1);
       Float_t baf2Time     = (Float_t) fBaF2->GetTime(module);
       Float_t vetoEnergy   = (Float_t) fVeto->GetEnergy(module);
       Float_t vetoTime     = (Float_t) fVeto->GetTime(module);


      calib_TAPSdEvEcl->Fill(baf2ClEnergy,vetoEnergy);
      calib_TAPSdEvE->Fill(baf2Energy,vetoEnergy);

      VetoData->Fill((Float_t)module, vetoRawEnergy, baf2Energy, baf2Time, vetoEnergy, vetoTime, baf2ClEnergy,baf2ClTime);
    }
  }

  fCalibVetoStep = 1;
}

// -------------------------------------------------------------------------------

Float_t TA2TAPS2009::TAPS_GetCLInfo(UInt_t module, UInt_t TimeOrEnergy)
{
  HitCluster_t** clBaF2; clBaF2 = fBaF2->GetCluster(); // -> cluster structs
  HitCluster_t* cl;                                    // cluster struct
  UInt_t* id_clBaF2; id_clBaF2 = fBaF2->GetClustHit();
  Int_t fBaF2_Ncluster = fBaF2->GetNCluster();

  for (Int_t i=0; i<fBaF2_Ncluster; i++)
  {
    cl = clBaF2[id_clBaF2[i]];        // i-th cluster
    if(cl->GetIndex()==module && TimeOrEnergy==0) return (Float_t)cl->GetEnergy();
    if(cl->GetIndex()==module && TimeOrEnergy==1) return (Float_t)cl->GetTime();
   }

  return (Float_t)0.0 ;
}

//--------------------------------------------------------------------------------

void TA2TAPS2009::SetConfig( char* line, int key )
{
  //UInt_t count = 0;
  fTapsRunStep = 1;
  UInt_t check = 0;

  switch(key)
  {
   case ETAPSFactor:
   {
    if(sscanf(line, "%lf", &fTapsFudge) < 1 ) goto error;
    break;
   }
   case ETAPS_SimpleReconstruct:
   {
    fSimpleReconstruct = kFALSE;
    if(sscanf(line, "%d", &check) < 1) printf("\nERROR\n");
    if(check) fSimpleReconstruct = kTRUE;
    break;
   }
   case ETAPS_dEvE_Cuts:
   {
    sscanf(line, "%s", fTAPS_dEvE_Cuts);
    TFile *f = new TFile(fTAPS_dEvE_Cuts, "READ");
    if(f)
    {
      f->Close();
      delete f;
    }
    else
    {
      printf("\n!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("!       File containing dEvE-CUTs is missing.       !\n");
      printf("!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
      delete f;
    }
    break;
   }
   case ETAPS_dEvE_Proton:
   {
    fTAPS_dEvE_Proton = kFALSE;
    sscanf(line, "%d", &check);
    if(check) fTAPS_dEvE_Proton = kTRUE; // dEvE-proton-cut active!
    break;
   }
   case ETAPS_dEvE_Proton_CutName:
   {
    if(fTAPS_dEvE_Proton)
    {
      sscanf(line, "%s", fTAPS_dEvE_Proton_CutName);
      fTAPS_dEvE_CutFile = new TFile(fTAPS_dEvE_Cuts, "READ");
      fTAPS_dEvE_ProtonCut = (TCutG*)fTAPS_dEvE_CutFile->Get(fTAPS_dEvE_Proton_CutName);
      fTAPS_dEvE_CutFile->Close();
    }
    break;
   }
   case ETAPS_dEvE_ChPion:
   {
    fTAPS_dEvE_ChPion = kFALSE;
    sscanf(line, "%d", &check);
    if(check) fTAPS_dEvE_ChPion = kTRUE;    // dEvE-pi+ cut active!
    break;
   }
   case ETAPS_dEvE_ChPion_CutName:
   {
    if(fTAPS_dEvE_ChPion == kTRUE)
    {
      sscanf(line, "%s", fTAPS_dEvE_ChPion_CutName);
      fTAPS_dEvE_CutFile = new TFile(fTAPS_dEvE_Cuts, "READ");
      fTAPS_dEvE_ChPionCut = (TCutG*)fTAPS_dEvE_CutFile->Get(fTAPS_dEvE_ChPion_CutName);
      fTAPS_dEvE_CutFile->Close();
    }
    break;
   }
   case ETAPS_dEvE_Electron:
   {
    fTAPS_dEvE_Electron = kFALSE;
    sscanf(line, "%d", &check);
    if(check) fTAPS_dEvE_Electron = kTRUE;     // dEvE-electron  cut  active!
    break;
   }
   case ETAPS_dEvE_Electron_CutName:
   {
    if(fTAPS_dEvE_Electron == kTRUE)
    {
      sscanf(line, "%s", fTAPS_dEvE_Electron_CutName);
      fTAPS_dEvE_CutFile = new TFile(fTAPS_dEvE_Cuts, "READ");
      fTAPS_dEvE_ElectronCut = (TCutG*)fTAPS_dEvE_CutFile->Get(fTAPS_dEvE_Electron_CutName);
      fTAPS_dEvE_CutFile->Close();
    }
    break;
   }
   default:
   {
    // Command not found...possible pass to next config
    TA2Apparatus::SetConfig( line, key );
    break;;
   }
  }
  return;

error:
  PrintError(line);
  return;
}

//-----------------------Thank you for reading ---------------------------

