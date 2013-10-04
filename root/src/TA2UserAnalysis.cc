//--Author	JRM Annand   27th Jan 2004...Apapt A.Starostin code
//--Rev         JRM Annand   12th May 2004...TA2Calorimeter etc. added
//--Rev         JRM Annand   21st Nov 2004...TA2CosmicCal added
//--Rev         JRM Annand   15th Jul 2005...TA2CrystalBall, TA2TAPS
//--Update      JRM Annand   15th Jul 2005...Use TA2Analysis::LoadVariable
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2UserAnalysis
//
// User analysis recognises user-written apparatus and physics classes
// New apparati and physics should be entered in Map_t kKnownChild
// and also in the switch statement of CreateChild() where the apparatus
// constructers are called

#include "TA2UserAnalysis.h"
#include "TA2Calorimeter.h"
#include "TA2Tagger.h"
#include "TA2CosmicCal.h"
#include "TA2CrystalBall.h"
#include "TA2CB.h"
#include "TA2TAPS.h"
#include "TA2TAPS2008.h"
#include "TA2TAPS2009.h"
#include "TA2Taps2009LE.h"
#include "TA2Taps.h"
#include "TA2Compton.h"
#include "TA2ThreshPi0.h"
#include "TA2SpinPolPhysics.h"
#include "TA2Efficiency.h"
#include "TA2GenericApp.h"
#include "TA2UserPhysics.h"
#include "TA2PhotoPhysics.h"

// Recognised apparatus classes.
// The "standard" set is held in TA2Analysis
enum { EA2Calorimeter, EA2Tagger, EA2LinearPol, EA2CosmicCal,
       EA2CrystalBall, EA2CB, EA2TAPS, EA2TAPS2008, EA2TAPS2009, 
       EA2Taps2009LE, EA2Taps, 
       EA2GenericApp, EA2Physics, EA2UserPhysics, EA2PhotoPhysics, 
       EA2Compton, EA2ThreshPi0, EA2SpinPolPhysics, EA2Efficiency, };
static const Map_t kKnownChild[] = {
  {"TA2Calorimeter", EA2Calorimeter},
  {"TA2Tagger",      EA2Tagger},
  {"TA2CosmicCal",   EA2CosmicCal},
  {"TA2CrystalBall", EA2CrystalBall},
  {"TA2CB",          EA2CB},
  {"TA2TAPS",        EA2TAPS},
  {"TA2TAPS2008",    EA2TAPS2008},
  {"TA2TAPS2009",    EA2TAPS2009},
  {"TA2Taps2009LE",  EA2Taps2009LE},
  {"TA2Taps",        EA2Taps},
  {"TA2GenericApp",  EA2GenericApp},
  // Physics
  {"TA2Physics",     EA2Physics},
  {"TA2UserPhysics", EA2UserPhysics},
  {"TA2PhotoPhysics",EA2PhotoPhysics},
  {"TA2Compton",     EA2Compton},
  {"TA2ThreshPi0",   EA2ThreshPi0},
  {"TA2SpinPolPhysics",   EA2SpinPolPhysics},
  {"TA2Efficiency",  EA2Efficiency},
  {NULL,             -1}
};

//-----------------------------------------------------------------------------

TA2UserAnalysis::TA2UserAnalysis( const char* name )
  :TA2Analysis( name )
{
  kValidChild = kKnownChild;              // valid apparatus list
}

//---------------------------------------------------------------------------

TA2UserAnalysis::~TA2UserAnalysis()
{
  // Free up allocated memory
}

//---------------------------------------------------------------------------

TA2DataManager* TA2UserAnalysis::CreateChild( const char* name, Int_t a )
{
  // Add creation of user-defined apparati here

  switch( a ){
  case EA2Tagger:
    // Standard tagger
    return new TA2Tagger( name, this );
  case EA2Calorimeter:
    // General Calorimeter...TAPS, CB etc.
    return new TA2Calorimeter( name, this );
  case EA2CosmicCal:
    // Cosmic calibration array of long bars
    return new TA2CosmicCal( name, this );
  case EA2CrystalBall:
    // Moded CB stuff
    return new TA2CrystalBall( name, this );
  case EA2CB:
    // Moded CB stuff
    return new TA2CB( name, this );
  case EA2TAPS:
    // Moded TAPS stuff
    return new TA2TAPS( name, this );
  case EA2TAPS2008:
    // Moded TAPS stuff (2008)
    return new TA2TAPS2008( name, this );
  case EA2TAPS2009:
    // Moded TAPS stuff (2009)
    return new TA2TAPS2009( name, this );
  case EA2Taps2009LE:
    // Moded TAPS stuff (2009)
    return new TA2Taps2009LE( name, this );
  case EA2Taps:
    // Henning's TAPS class
    return new TA2Taps( name, this );
   case EA2GenericApp:
    // Generic
    return new TA2GenericApp( name, this );
    //
    // Physics stuff
  case EA2Physics:
    // Default (dummy physics)
    return new TA2Physics( name, this );
  case EA2UserPhysics:
    // User defined Physics
    return new TA2UserPhysics( name, this );
  case EA2PhotoPhysics:
    // General photonuclear stuff
    return new TA2PhotoPhysics( name, this );
  case EA2Compton:
    // Compton
    return new TA2Compton( name, this );
  case EA2ThreshPi0:
    // Pi0 threshold
    return new TA2ThreshPi0( name, this );
  case EA2SpinPolPhysics:
    // Compton Spin POl
    return new TA2SpinPolPhysics( name, this );
  case EA2Efficiency:
    // Pi0 threshold
    return new TA2Efficiency( name, this );
  default:
    PrintError((char*)name,
	       "<Unknown apparatus..cannot continue>", EErrFatal);
  }
  return NULL;
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::LoadVariable( )
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
  //  TA2DataManager::LoadVariable("Mmiss",  &fMmiss,         EDSingleX);
  TA2Analysis::LoadVariable();
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::ParseDisplay( char* line )
{
  TA2Analysis::ParseDisplay( line );            // standard parse
  //  if( !fIsConfigPass ) return;              // was it a standard command?
  //  TA2HistManager::ParseDisplay( line );
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::InitSaveTree(Char_t* fname)
{
  if(fPhysics)fPhysics->InitSaveTree();
  TA2Analysis::InitSaveTree(fname);
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::ChangeTreeFile(Char_t* fname)
{
  CloseEvent();        //Close tree & file
  InitSaveTree(fname); //New tree & file
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::CloseEvent()
{
  if(fPhysics) fPhysics->CloseEvent();
  TA2Analysis::CloseEvent();
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::Periodic()
{
  if(fPhysics) fPhysics->Periodic();
  TA2Analysis::Periodic();
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::EndFile()
{
  if(fPhysics) fPhysics->EndFile();
  TA2Analysis::EndFile();
  fFileDone = true;
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::Finish()
{
  if(IsEndOfFile()) EndFile();
  if(fPhysics) fPhysics->Finish();
  TA2Analysis::Finish();
  fFileDone = true;
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::Decode()
{
  //Set correct file name. This code tries to get the new file name as long as
  //it differs from the previous one, in order to get along with possible
  //race conditions in file name handling when changing files
  if(gAR->IsOnline() && fFileDone) //Previous file was finished and file name has really changed
  {
    Char_t* TempName = strrchr(gAR->GetFileName(), '/') + 1;
    //If new file name is different from previous one...
    if(strcmp(fFileName, TempName))
    {
      delete[] fFileName;
      fFileName = new Char_t[strlen(TempName) + 1];
      strcpy(fFileName, TempName); //...get new name...
      fFileDone = false;           //...and don't try again
    }
  }

  //Call standard decode
  TA2Analysis::Decode();
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::ZeroHist(Char_t* name)
{
  //Reset histogram, find histogram by name

  //Check there is something to display
  if(!fIsDisplay)
  {
    printf(" Sorry no histograms enabled for this %s\n", this->ClassName());
    return;
  }

  TH1* h = (TH1*)(gROOT->FindObject(name)); //All histograms inherit from TH1
  if(h)
    h->Reset();
  else
    //Histogram not found...print message
    printf("Histogram %s not found for class %s\n", name, this->ClassName());
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::ZeroAll()
{
  if(fPhysics) fPhysics->ZeroAll();
  TA2Analysis::ZeroAll();
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::SaveAll(Char_t* name)
{
  TFile AcquHist(name, "RECREATE");
  gROOT->GetList()->Write();
  AcquHist.Close();
}

//-----------------------------------------------------------------------------

void TA2UserAnalysis::SaveAll()
{
  TFile AcquHist("ARHistograms.root", "RECREATE");
  gROOT->GetList()->Write();
  AcquHist.Close();
}

//-----------------------------------------------------------------------------
