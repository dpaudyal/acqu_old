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

#ifndef __TA2UserAnalysis_h__
#define __TA2UserAnalysis_h__


#include "TA2Analysis.h"            // Acqu-Root histogrammer

class TA2UserAnalysis : public TA2Analysis {
 protected:
  Char_t* fFileName; //Name of the file currently analysed
  Bool_t  fFileDone; //Flag, if analysed file is finished
 public:
  TA2UserAnalysis( const char* );        // normal use, pass ptr to TAcquRoot
  virtual ~TA2UserAnalysis();
  virtual TA2DataManager* CreateChild( const char*, Int_t );
  virtual void ParseDisplay( char* );   // display setup
  virtual void LoadVariable( );         // cut/display setup
  virtual void InitSaveTree(Char_t*);
  virtual void ChangeTreeFile(Char_t*);
  virtual void CloseEvent();
  virtual void Periodic();
  virtual void EndFile();
  virtual void Finish();
  virtual void Decode();
  virtual Char_t* GetFileName() const { return fFileName; };
  void ZeroAll();
  void ZeroHist(Char_t*);
  void SaveAll();
  void SaveAll(Char_t*);
  
  ClassDef(TA2UserAnalysis,1)
};


//inline void TA2UserAnalysis::Cleanup()
//{
//}

#endif


