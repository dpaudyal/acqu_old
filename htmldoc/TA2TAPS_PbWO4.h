////////////////////////////////////////////////////////////////////////////
//
//--Author	JRM Annand   30th Sep 2003
//--Update      H.Berghaeuser June 2008
//
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2TAPS_PbWO4.h
//
//
//
///////////////////////////////////////////////////////////////////////////

#ifndef __TA2TAPS_PbWO4_h__
#define __TA2TAPS_PbWO4_h__

#include "MCBranchID.h"
#include "TA2Detector.h"

class TA2TAPS_PbWO4 : public TA2Detector
{
 private:
 public:
  TA2TAPS_PbWO4(const char*, TA2System*); // Normal use
  virtual ~TA2TAPS_PbWO4();
  //virtual void SetConfig(char*, int); // decode and load setup info
  virtual void LoadVariable();            // display/cut setup
  //virtual void PostInit();            // initialise using setup info
  virtual void Decode();              // hits -> energy procedure
  virtual void SaveDecoded();             // save local analysis
  virtual void ReadDecoded();             // read back previous analysis

  ClassDef(TA2TAPS_PbWO4,1)
};

//---------------------------------------------------------------------------

inline void TA2TAPS_PbWO4::ReadDecoded()
{
  // Read back...
  //   either previously analysed data from Root Tree file
  //   or MC simulation results, assumed to have the same data structure
  // Bug fix remove fNhits-- line, (DG 27/5/05, implemented JRMA 2/6/05)
  // D.Glazier addition of time decoding from A2 Geant4 model 24/08/07

  UInt_t i,j;
  Double_t total = 0.0;
  fNhits = *(Int_t*)(fEvent[EI_vhits]);
  Float_t* energy = (Float_t*)(fEvent[EI_eveto]);
  Float_t* time = NULL;
  if(fIsTime) time = (Float_t*)(fEvent[EI_tveto]);
  Int_t* index = (Int_t*)(fEvent[EI_iveto]);

  for( i=0; i<fNhits; i++ )
  {
    j = *index++;
    if(!j) // is it s real hit
    {
      //fNhits--;
      i--;
      energy++;
      if(fIsTime)time++;
      continue;
    }
    j--;
    fHits[i] = j;
    fEnergy[j] = (*energy++) * 1000.;
    fEnergyOR[i] = fEnergy[j];
    if(fIsTime)
    {
      fTime[j] = (*time++);
      fTimeOR[i] = fTime[j];
    }
    total += fEnergy[j];
  }
  fHits[i] = EBufferEnd;
  fEnergyOR[i] = EBufferEnd;
  if(fIsTime)fTimeOR[i] = EBufferEnd;
  fTotalEnergy = total;
}

//---------------------------------------------------------------------------

inline void TA2TAPS_PbWO4::Decode()
{
  // Run basic TA2Detector decode, then anything else

  DecodeBasic();
  // Anything else here
}

//---------------------------------------------------------------------------

#endif
