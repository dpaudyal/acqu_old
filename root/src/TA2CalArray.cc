//--Author	JRM Annand   27th Apr 2003
//--Rev 	JRM Annand...30th Sep 2003  New Detector bases
//--Rev 	JRM Annand...15th Oct 2003  ReadDecoded...MC data
//--Rev 	JRM Annand... 5th Feb 2004  3v8 compatible
//--Rev 	JRM Annand... 7th Jun 2005  ReadDecoded...use fEnergyScale
//--Rev 	JRM Annand...25th Oct 2005  ReadDecoded...energy thresh
//--Rev 	D.Glazier ...24th Aug 2007  ReadDecoded...include detector time
//--Update	JRM Annand ..17th Nov 2007  ReadDecoded total energy fix
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2CalArray
//
// Decoding and calibration methods for the Crystall Ball array of NaI(Tl)
// Configured as forward wall to work in conjunction with the CB
// This can use standard TA2ClusterDetector methods or ones redifined here
//

#include "TA2CalArray.h"
#include <iostream>

enum { ECAEnergyResolution=1900, ECATimeResolution, ECAThetaResolution, ECAPhiResolution, ECAClusterUCLA };

static const Map_t kCalArrayKeys[] = {
  {"Energy-Resolution:",   ECAEnergyResolution},
  {"Time-Resolution:",     ECATimeResolution},
  {"Theta-Resolution:",    ECAThetaResolution},
  {"Phi-Resolution:",      ECAPhiResolution},
  {"Max-Cluster:",         EClustDetMaxCluster},
  {"Next-Neighbour:",      EClustDetNeighbour},
  {"Cluster-UCLA:",        ECAClusterUCLA},
  {NULL,            -1}
};
ClassImp(TA2CalArray)

//---------------------------------------------------------------------------
TA2CalArray::TA2CalArray(const char* name, TA2System* apparatus)
            :TA2ClusterDetector(name, apparatus)
{
  // Do not allocate any "new" memory here...Root will wipe
  // Set private variables to zero/false/undefined state
  // Pass kLaddKeys (command-control keywords) and
  // kLaddHist (names of histograms) to progenitor classes

  fTryHits = NULL;
  fTempHits2 = NULL;

  fUseSigmaEnergy       = 0;
  fUseSigmaTime         = 0;
  fUseClusterDecodeUCLA = 0;
  fSigmaEnergyFactor    = -1.0;
  fSigmaEnergyPower     = -1.0;
  fSigmaTime            = -1.0;
  fSigmaTheta           = -1.0;
  fSigmaPhi             = -1.0;
  fEthresh              = 0.0;

  fRandom = new TRandom();

  AddCmdList(kCalArrayKeys);                  // for SetConfig()
}

//---------------------------------------------------------------------------

TA2CalArray::~TA2CalArray()
{
  // Free up all allocated memory
  // ...arrays created at the initialisation stage
  // Start with common arrays of TA2Detector class
  DeleteClusterArrays();
}

//---------------------------------------------------------------------------

void TA2CalArray::SetConfig(char* line, int key)
{
  // Load CalArray parameters from file or command line
  // CalArray specific configuration
  switch(key)
  {
  case ECAClusterUCLA:
    // Set Clustering Algorithm to UCLA version 
    fUseClusterDecodeUCLA = 1;
    break;
  case ECAEnergyResolution:
    // Energy Resolution Read-in Line
    if(sscanf(line, "%lf%lf%d", &fSigmaEnergyFactor, &fSigmaEnergyPower, &fUseSigmaEnergy) < 3)
      PrintError(line,"<TA2CalArray Energy Resolution>");
    break;
  case ECATimeResolution:
    // Time resolution read-in line
    if(sscanf(line, "%lf%d", &fSigmaTime, &fUseSigmaTime) < 2)
      PrintError(line,"<TA2CalArray Time Resolution>");
    break;
  case ECAThetaResolution:
    // Time resolution read-in line
    if(sscanf(line, "%lf", &fSigmaTheta) < 1)
      PrintError(line,"<TA2CalArray Theta Resolution>");
    break;
  case ECAPhiResolution:
    // Time resolution read-in line
    if(sscanf(line, "%lf", &fSigmaPhi) < 1)
      PrintError(line,"<TA2CalArray Phi Resolution>");
    break;
  case EClustDetMaxCluster:
    // Max number of clusters
    if(sscanf( line, "%d%lf", &fMaxCluster, &fClEthresh ) < 1)
    {
      PrintError(line,"<Parse UCLA # clusters>");
      break;
    }
    if(fUseClusterDecodeUCLA)
      {
	fUseClusterDecodeUCLA = 1;
	fEthresh = fClEthresh;
	fClusterUCLA = new HitClusterUCLA_t*[fNelement+1];
	fCluster = (HitCluster_t**)fClusterUCLA;
	fClustHit = new UInt_t[fMaxCluster+1];
	fTryHits = new UInt_t[fNelement+1];
	fTempHits = new UInt_t[fNelement+1];
	fTempHits2 = new UInt_t[fNelement+1];
	fNClustHitOR = new UInt_t[fNelement+1];
	fTheta = new Double_t[fNelement+1];
	fPhi = new Double_t[fNelement+1];
	fClEnergyOR = new Double_t[fNelement+1];
	fClTimeOR = new Double_t[fNelement+1];
	fClCentFracOR = new Double_t[fNelement+1];
	fClRadiusOR = new Double_t[fNelement+1];
	fNCluster = 0;
      } else TA2ClusterDetector::SetConfig(line, key);
    break;
  case EClustDetNeighbour:
    // Nearest neighbout input
    if(fUseClusterDecodeUCLA)
      {
	if( fNCluster < fNelement )
	  fClusterUCLA[fNCluster] = new HitClusterUCLA_t(line,fNCluster,fClustSizeFactor);
	fNCluster++;
      }
    else TA2ClusterDetector::SetConfig(line, key);
    break;
  default:
    // Command not found...possible pass to next config
    TA2ClusterDetector::SetConfig(line, key);
    break;
  }
  return;
}
//-----------------------------------------------------------------------------
void TA2CalArray::ParseDisplay(char* line)
{
  // Input private histogram spec to general purpose parse
  // and setup routine

  //  const Name2Variable_t hist[] = {
  // Name          ->variable  single/mult    Fill-condition
  //    {"Nphoton_Minv", fM_Nphoton, EHistSingleX},
  //    {NULL,          0,         0}
  //  };
  // Do not remove the final NULL line

  TA2ClusterDetector::ParseDisplay(line);
  return;
}

//---------------------------------------------------------------------------
void TA2CalArray::PostInit()
{
  // Some further initialisation after all setup parameters read in
  // Start with alignment offsets
  // Create space for various output arrays

  if(fUseClusterDecodeUCLA){
    std::cout<<""<<std::endl;
    std::cout<<"!!!!!! WARNING - Using UCLA ClusterDecode !!!!!!"<<std::endl;
    std::cout<<"!!   This code is not yet fully standardised  !!"<<std::endl;
    std::cout<<"!! Please report any observed faults to cbana !!"<<std::endl;
    std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
  }
  fEnergyAll = new Double_t[fNelem+1];

  TA2ClusterDetector::PostInit();
}

//---------------------------------------------------------------------------

void TA2CalArray::SaveDecoded()
{
  // Save decoded info to Root Tree file
}

//---------------------------------------------------------------------------

