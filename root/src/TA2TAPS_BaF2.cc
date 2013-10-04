//--Author      JRM Annand   30th Sep 2003
//--Rev         JRM Annand...15th Oct 2003  ReadDecoded...MC data
//--Update      JRM Annand    5th Feb 2004  3v8 compatibility
//--Update        R Gregor    3th Jul 2005  Added veto and pulsshape reconstruction
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2TAPS_BaF2
//
// Decoding and calibration methods for TAPS array of BaF2
// Configured as forward wall to work in conjunction with the CB
//

#include "TA2TAPS_BaF2.h"
#include "HitClusterTAPS_t.h"

// constants for command-line maps below
enum {
  ETAPSSGMaxDet=510, ETAPSSG, EClustDetMaxTAPSCluster, EClustDetTAPSNeighbour,
  ETAPSEnergyResolution, ETAPSTimeResolution, ETAPSThetaResolution, ETAPSPhiResolution,
  ETAPSVetoEnergy, ETAPSVetoEfficiency, ETAPSVetoThreshold,
};

// Command-line key words which determine what to read in
static const Map_t kTAPSClustDetKeys[] = {
  {"TAPSSGMaxDet:",        ETAPSSGMaxDet},
  {"TAPSSG:",              ETAPSSG},
  {"Max-TAPSCluster:",     EClustDetMaxTAPSCluster},
  {"Next-TAPS:",           EClustDetTAPSNeighbour},
  {"Energy-Resolution:",   ETAPSEnergyResolution},
  {"Time-Resolution:",     ETAPSTimeResolution},
  {"Theta-Resolution:",    ETAPSThetaResolution},
  {"Phi-Resolution:",      ETAPSPhiResolution},
  {"Veto-Energy:",         ETAPSVetoEnergy},
  {"Veto-Efficiency:",     ETAPSVetoEfficiency},
  {"Veto-Threshold:",      ETAPSVetoThreshold},
  {NULL,          -1}
};

ClassImp(TA2TAPS_BaF2)

//---------------------------------------------------------------------------

TA2TAPS_BaF2::TA2TAPS_BaF2(const char* name, TA2System* apparatus)
             :TA2ClusterDetector(name, apparatus)
{
  fUseEnergyResolution    = 0;
  fUseTimeResolution      = 0;
  fEnergyResolutionFactor = -1.0;
  fEnergyResolutionConst  = -1.0;
  fTimeResolution         = -1.0;
  fThetaResolution        = -1.0;
  fPhiResolution          = -1.0;

  fVetoThreshold = 0.0;
  for(Int_t t=0; t<5; t++)
  {
    fVetoEfficiency[t] = 1.0;
    fVetoEnergy[t] = 100.0*t;
  }

  fSGEnergy = NULL;
  fLGEnergy = NULL;
  fMaxSGElements = 0;
  // Do not allocate any "new" memory here...Root will wipe
  // Set private variables to zero/false/undefined state
  AddCmdList( kTAPSClustDetKeys );
}

//---------------------------------------------------------------------------

TA2TAPS_BaF2::~TA2TAPS_BaF2()
{
  // Free up all allocated memory
  // ...arrays created at the initialisation stage
  for(UInt_t i=0; i<fMaxSGElements; i++)
    if(fSGEnergy[i]) delete fSGEnergy[i];

  delete fSGEnergy;
  DeleteClusterArrays();
}

//---------------------------------------------------------------------------

void TA2TAPS_BaF2::SaveDecoded()
{
  // Save decoded info to Root Tree file
}

//---------------------------------------------------------------------------

void TA2TAPS_BaF2::SetConfig(Char_t* line, Int_t key)
{
  // Default detector configuration
  // If this flags it has done nothing carry on with local setup
  //TA2ClusterDetector::SetConfig( line, key );
  //if( !fIsConfigPass ) return;
  //fIsConfigPass = EFalse;           // reset for any inherited SetConfig

  // Own special cluster setup
  // Cluster specific configuration
  switch(key)
  {
   case ETAPSEnergyResolution:
    // Energy Resolution Read-in Line
    if(sscanf(line, "%lf%lf%d", &fEnergyResolutionFactor, &fEnergyResolutionConst, &fUseEnergyResolution) < 3)
      PrintError(line,"<TA2TAPS_BaF2 Energy Resolution>");
    break;
   case ETAPSTimeResolution:
    // Time resolution read-in line
    if(sscanf(line, "%lf%d", &fTimeResolution, &fUseTimeResolution) < 2)
      PrintError(line,"<TA2TAPS_BaF2 Time Resolution>");
    break;
   case ETAPSThetaResolution:
    // Time resolution read-in line
    if(sscanf(line, "%lf", &fThetaResolution) < 1)
      PrintError(line,"<TA2TAPS_BaF2 Theta Resolution>");
    break;
   case ETAPSPhiResolution:
    // Time resolution read-in line
    if(sscanf(line, "%lf", &fPhiResolution) < 1)
      PrintError(line,"<TA2TAPS_BaF2 Phi Resolution>");
    break;

   case ETAPSVetoThreshold:
    if(sscanf(line, "%lf", &fVetoThreshold) < 1)
      PrintError(line,"<TA2TAPS_BaF2 Veto Threshold>");
   case ETAPSVetoEfficiency:
    if(sscanf(line, "%lf%lf%lf%lf%lf", &fVetoEfficiency[0], &fVetoEfficiency[1], &fVetoEfficiency[2], &fVetoEfficiency[3], &fVetoEfficiency[4]) < 5)
      PrintError(line,"<TA2TAPS_BaF2 Veto Efficiency>");
    break;
   case ETAPSVetoEnergy:
    if(sscanf(line, "%lf%lf%lf%lf", &fVetoEnergy[0], &fVetoEnergy[1], &fVetoEnergy[2], &fVetoEnergy[3]) < 4)
      PrintError(line,"<TA2TAPS_BaF2 Veto Energy>");
    break;

   case ETAPSSGMaxDet:
    //Max number of BaF2-SG-clusters
    if(sscanf(line, "%d", &fMaxSGElements) < 1) goto error;
    fSGEnergy = new HitD2A_t*[fMaxSGElements+1];
    fSGEnergyValue = new Double_t[fMaxSGElements+1];
    fNSG = 0;
    break;
   case ETAPSSG:
    //SG-Pedestal and Gain of BaF2-Crystals
    // parameters for each element...need to input # of elements 1st
    if( fNSG>=fMaxSGElements )
    {
      PrintError(line,"<Too many detector elements input>");
      break;
    }
    fSGEnergy[fNSG] = new HitD2A_t( line, fNSG, this );
    break;
   case EClustDetMaxTAPSCluster:
    // Max number of clusters
    if(sscanf(line, "%d%lf", &fMaxCluster, &fClEthresh) < 1) goto error;
    fEthresh = fClEthresh;
    fCluster = new HitCluster_t*[fNelement+1];
    fprintf(fLogStream, " HitClusterTAPS_t class used for shower reconstruction\n");
    fClustHit = new UInt_t[fMaxCluster+1];
    fTempHits = new UInt_t[fNelement+1];
    fNClustHitOR = new UInt_t[fNelement+1];
    fTheta = new Double_t[fNelement+1];
    fPhi      = new Double_t[fNelement+1];
    fClEnergyOR  = new Double_t[fNelement+1];
    fClTimeOR  = new Double_t[fNelement+1];
    fClCentFracOR  = new Double_t[fNelement+1];
    fClRadiusOR  = new Double_t[fNelement+1];
    fNCluster = 0;
    break;
 case EClustDetTAPSNeighbour:
    // Nearest neighbout input
    if(fNCluster < fNelement) fCluster[fNCluster] = new HitClusterTAPS_t(line, fNCluster);
    fNCluster++;
    break;
   default:
    // Command not found...possible pass to next config
    TA2ClusterDetector::SetConfig(line, key);
    fIsConfigPass = ETrue;
    break;
  }
  return;
 error: PrintError(line);
  return;
}

//---------------------------------------------------------------------------

void TA2TAPS_BaF2::PostInit()
{
  // Some further initialisation after all setup parameters read in
  // Start with alignment offsets
  // Create space for various output arrays
  fRandom=new TRandom();
  pRandoms = new TRandom3;
  pRandoms->SetSeed(0); //'Randomize Timer'

  EnergyAll = new Double_t[fNelem];
  theta_c = new Double_t[fNelem];
  phi_c = new Double_t[fNelem];
  x_c = new Double_t[fNelem];
  y_c = new Double_t[fNelem];
  z_c = new Double_t[fNelem];
  Z_c = new Double_t[fNelem];

  TVector3 direction;
  for(UInt_t t=0; t<fNelem; t++)
  {
    direction = fPosition[t]->Unit();
    x_c[t] = -direction.X(); //X position patch in G3 simulation
    y_c[t] = direction.Y();
    z_c[t] = direction.Z();
    Z_c[t] = fPosition[t]->Z();
    theta_c[t] = acos(z_c[t]);
    phi_c[t] = atan2(y_c[t], x_c[t]);
  }

  bSimul = false;

  TA2ClusterDetector::PostInit();
  fLGEnergy = TA2Detector::GetElement();
}

//---------------------------------------------------------------------------

// inline void TA2TAPS_BaF2::ReadDecoded()
// {
//   // Read back "decoded" data for the BaF2 array
//   // In this case produced by the GEANT3 CB simulation

//   UInt_t j;
//   Int_t nVetos = 0;
//   Int_t nBaF2s = 0;
//   Double_t total = 0.0;
//   fNhits = *(Int_t*)(fEvent[EI_ntaps]);
//   Float_t* energy = (Float_t*)(fEvent[EI_ectapsl]);
//   Float_t* time = (Float_t*)(fEvent[EI_tctaps]);
//   UInt_t* index = (UInt_t*)(fEvent[EI_ictaps]);
//   Float_t* venergy = (Float_t*)(fEvent[EI_evtaps]);
//   Float_t* dircos = (Float_t*)(fEvent[EI_dircos]);
//   Int_t npart = *(Int_t*)(fEvent[EI_npart]);
//   Int_t* idpart = (Int_t*)(fEvent[EI_idpart]);
//   Double_t E, T;                                    //Energy, time
//   Double_t Lo, Hi;
//   Int_t proton;
//   Double_t dphi, dtheta;
//   Double_t domega2, temp2;
//   UInt_t crystal;
//   Double_t Efficiency = 0.0;
//   Double_t GammaToF;

//   //Find proton index in simulated data
//   for(proton=0; proton<npart; proton++)
//     if(idpart[proton]==14) break; //14 = GEANT proton
//   //Move pointer to correct positions for proton data
//   for(Int_t t=0; t<proton; t++) //For all particles before proton:
//     for(Int_t i=0;i<3; i++) dircos++; //Jump over x,y,z components
//   //Read x,y,z components for proton
//   Double_t x_p = *dircos++;
//   Double_t y_p = *dircos++;
//   Double_t z_p = *dircos;
//   //Calculate theta and phi angles out of cartesian components
//   Double_t theta_p = acos(z_p);
//   Double_t phi_p = atan2(y_p, x_p);

//   //For each crystal: evaluate deviation from proton direction in theta-phi plane
//   domega2 = 64800.0; //180^2 + 180^2
//   crystal = 65535;
//   if(theta_p < (TMath::DegToRad()*20.0))
//     for(UInt_t t=0; t<fNelem; t++)
//     {
//       dtheta = theta_p - theta_c[t];
//       dphi = phi_p - phi_c[t];
//       if(dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
//       if(dphi < -TMath::Pi()) dphi = -TMath::TwoPi() - dphi;
//        //Find crystal which fits best for proton direction in theta-phi plane
//       temp2 = dphi*dphi + dtheta*dtheta;
//       if(temp2 < domega2)
//       {
//         domega2 = temp2;
//         crystal = t;
//       }
//     }

//   //I want this empty...
//   if(fMaxSGElements)
//     for(UInt_t t=0; t<fMaxSGElements; t++)
//       fSGEnergyValue[t] = EBufferEnd;
//   for(UInt_t t=0; t<fNelem; t++)
//   {
//     //    fPatternHits[0][t] = EBufferEnd;
//     EnergyAll[t] = 0.0;
//   }

//   //Evaluate all simulated hits
//   for(UInt_t i=0; i<fNhits; i++)
//   {
//     j = index[i];
//     if(!j) break;
//     //Use following line for cbsim since 2006-06-22!
//     j--; //For reasons Stefan won't tell...
//     E = energy[i] * fEnergyScale;                                   //G3/4 output in GeV
//     if(fUseEnergyResolution) E+=fRandom->Gaus(0.0, GetEnergyResolutionGeV(E));
//     E*=1000.0;                                                      //G3/4 output in MeV
//     EnergyAll[j] = E;
//     GammaToF = (Z_c[j] * TMath::Cos(theta_c[j]))/30.0;
//     T = -time[i] + GammaToF + pRandoms->Gaus(0.0, fTimeResolution);

//     //Fill veto stuff
//     /*if((venergy[i] * 1000.0) > 0.1) //If energy deposited: Veto Hit*/
//     if((venergy[i] * 1000.0) > fVetoThreshold) //If energy deposited: Veto Hit
//     {
//       //Apply Victor's energy depended (but overall) veto efficiencies
//       /*if(E < 50.0) Efficiency = 0.81;
//       else if((E > 50.0) && (E < 100.0)) Efficiency = 0.78;
//       else if((E > 100.0) && (E < 200.0)) Efficiency =  0.74;
//       else if((E > 200.0) && (E < 300.0)) Efficiency = 0.68;
//       else if(E > 300.0) Efficiency = 0.64;*/
//       if(E < fVetoEnergy[0]) Efficiency = fVetoEfficiency[0];
//       else if((E > fVetoEnergy[0]) && (E < fVetoEnergy[1])) Efficiency = fVetoEfficiency[1];
//       else if((E > fVetoEnergy[1]) && (E < fVetoEnergy[2])) Efficiency = fVetoEfficiency[2];
//       else if((E > fVetoEnergy[2]) && (E < fVetoEnergy[3])) Efficiency = fVetoEfficiency[3];
//       else if(E > fVetoEnergy[3]) Efficiency = fVetoEfficiency[4];
//       if(pRandoms->Rndm() < Efficiency)
//       {
// 	//        fPatternHits[0][nVetos] = j;
//         nVetos++;
//       }
//     }

//     //Fill BaF2 stuff
//     //Lo = fElement[j]->GetEnergyLowThr();
//     //Hi = fElement[j]->GetEnergyHighThr();
//     Lo = fElement[j]->GetEnergyLowThr();
//     Hi = fElement[j]->GetEnergyHighThr();
//     if((E > Lo) && (E < Hi))
//     {
//       fEnergy[j] = E;
//       fTime[j] = T;
//       fEnergyOR[nBaF2s] = E;
//       fTimeOR[nBaF2s] = T;
//       fHits[nBaF2s] = j;
//       total+=E;
//       nBaF2s++;

//       //PSA - reasonable?
//       if(fMaxSGElements)
//       {
//         if((j==crystal) || fCluster[j]->IsNeighbour(crystal)) //Hit in or near crystal best fitting for proton
//           //fSGEnergyValue[j] = E * TMath::Power(E, 0.05) * 0.65; //Sven's function
//           fSGEnergyValue[j] = -0.9691614 + 0.7108238*E + 0.1819194e-2*E*E - 0.9019634e-5*E*E*E + 0.1610750e-7*E*E*E*E; //Marc's function
//         else
//           fSGEnergyValue[j] = 0.96 * E;
//         fSGEnergyValue[j]+=pRandoms->Gaus(0.0, 0.25 * TMath::Sqrt(fSGEnergyValue[j]));
//       }
//     }
//   }
//   fHits[nBaF2s] = EBufferEnd;
//   fTimeOR[nBaF2s] = EBufferEnd;
//   fEnergyOR[nBaF2s] = EBufferEnd;
//   fNhits = nBaF2s;
//   fTotalEnergy = total;

//   //  fPatternHits[0][nVetos] = EBufferEnd;
//   //  fNPatternHits[0] = nVetos;

// //   fRawEnergyHits[0] = EBufferEnd;
// //   fRawTimeHits[0] = EBufferEnd;

//   bSimul = true;
// }

// inline void TA2TAPS_BaF2::ReadDecoded( )
// {
//   // Read back "decoded" data for the BaF2 array
//   // In this case produced by the GEANT3 CB simulation

//   UInt_t i,j;
//   Double_t total = 0.0;
//   fNhits = *(Int_t*)(fEvent[EI_ntaps]);
//   Float_t* energy = (Float_t*)(fEvent[EI_ectapsl]);
//   Float_t* time = (Float_t*)(fEvent[EI_tctaps]);
//   Int_t* index = (Int_t*)(fEvent[EI_ictaps]);
//   Float_t* venergy = (Float_t*)(fEvent[EI_evtaps]);

//   for( i=0; i<fNhits; i++ ){
//     j = (*index++);
//     //j--;
//     fHits[i] = j;
//     fEnergy[j] = *energy++ * 1000.;
//     fEnergyOR[i] = fEnergy[j];
//     fTimeOR[i] = *time;
//     fTime[j] = *time++;
//     total += fEnergy[j];
//   }
//   fHits[i] = EBufferEnd;
//   fTimeOR[i] = EBufferEnd;
//   fEnergyOR[i] = EBufferEnd;
//   fTotalEnergy = total;
// }

inline void TA2TAPS_BaF2::ReadDecoded( )
{
  // Read back "decoded" data for the BaF2 array
  // In this case produced by the GEANT3 CB simulation
  // Remove j-- in indexing 2/11/05
  // Add energy thresholds 2/11/05
  // Note in Geant4 A2 simulation the numbering runs from 1-384...include j-- line

  UInt_t i,j,k;
  Double_t total = 0.0;                               // total energy counter
  fNhits = *(Int_t*)(fEvent[EI_ntaps]);               // detector hits
  Float_t* energy = (Float_t*)(fEvent[EI_ectapsl]);   // energy array
  Float_t* time = (Float_t*)(fEvent[EI_tctaps]);      // time array
  Int_t* index = (Int_t*)(fEvent[EI_ictaps]);         // element indices
  Double_t E;                                         // energy
  Double_t EthrLow,EthrHigh;                          // energy thresholds
  Float_t en;
  Float_t resolution=0.03;

  for( i=0,k=0; i<fNhits; i++ ){                      // loop over hits
    j = (*index++);
    j--;                                              // Need this line for A2 G4 simulation
    EthrLow = fElement[j]->GetEnergyLowThr();         // get low thresh
    EthrHigh = fElement[j]->GetEnergyHighThr();       // and high thresh
    en=(*energy++);
    E = 1000 * fEnergyScale* (en + resolution*fRandom->Gaus()*TMath::Power(en,0.75));  
    if( (E < EthrLow) || (E > EthrHigh) ) {*time++;continue;} // within limits ?
    fHits[k] = j;                                     // store index
    fEnergy[j] = E;                                   // store energy
    fEnergyOR[k] = E;                                 
    fTimeOR[k] = *time;                               // store time
    fTime[j] = *time++;
    total += E;                                       // energy total counter
    k++;;                                             // hits inside thresh
  }
  fHits[k] = EBufferEnd;                              // end markers
  fTimeOR[k] = EBufferEnd;
  fEnergyOR[k] = EBufferEnd;
  fNhits = k;                                         // # hits inside thresh
  fTotalEnergy = total;                               // total energy
}

//-----------------------------------------------------------------------------
