//////////////////////////////////////////////////////////////////////////////
//
//--Author	JRM Annand   30th Sep 2003
//--Update      H.Berghaeuser June 2008
//
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2TAPS_PbWO4.cc
//
//
//
///////////////////////////////////////////////////////////////////////////

#include "TA2TAPS_PbWO4.h"

ClassImp(TA2TAPS_PbWO4)

//---------------------------------------------------------------------------

TA2TAPS_PbWO4::TA2TAPS_PbWO4( const char* name, TA2System* apparatus )
              :TA2Detector(name, apparatus)
{
  // Default detector initialisation
}

//---------------------------------------------------------------------------

TA2TAPS_PbWO4::~TA2TAPS_PbWO4()
{
  // Free up all allocated memory
  // ...arrays created at the initialisation stage
  DeleteArrays();
}

//-----------------------------------------------------------------------------

void TA2TAPS_PbWO4::LoadVariable()
{
  // Input name - variable pointer associations for any subsequent
  // cut or histogram setup.
  // LoadVariable( "name", pointer-to-variable, type-spec );
  // NB scaler variable pointers need the preceeding &
  //    array variable pointers do not.
  // type-spec ED prefix for a Double_t variable
  //           EI prefix for an Int_t variable
  // type-spec SingleX for a single-valued variable
  //           MultiX  for a multi-valued variable

  //                            name     pointer          type-spec
  //  TA2DataManager::LoadVariable("Energy", &fTotalEnergy,   EDSingleX);

  // Call the standard detector
  // name-variable load
  TA2Detector::LoadVariable();
}

//---------------------------------------------------------------------------

void TA2TAPS_PbWO4::SaveDecoded()
{
  // Save decoded info to Root Tree file
}

//---------------------------------------------------------------------------
