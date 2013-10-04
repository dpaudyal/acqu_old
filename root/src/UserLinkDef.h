//--Author	JRM Annand    9th Sep 2004
//--Rev
//--Update
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// ROOT dictionary specification
// All classes linked in libUserRoot.so must included in the ROOT
// dictionary so that they are recognised by the CINT command processor

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// Global pointers to core TAcquRoot and TA2Analysis classes
#pragma link C++ global gAR;
#pragma link C++ global gAN;

// Any new user written class must be included in this list
// Main control and analysis
#pragma link C++ class TA2UserControl+;
#pragma link C++ class TA2UserAnalysis+;
// Physics classes
#pragma link C++ class TA2UserPhysics+;
#pragma link C++ class TA2PhotoPhysics+;
#pragma link C++ class TA2Compton+;
#pragma link C++ class TA2ThreshPi0+;
#pragma link C++ class TA2SpinPolPhysics+;
#pragma link C++ class TA2Efficiency+;
// Apparatus classes
#pragma link C++ class TA2GenericApp+;
#pragma link C++ class TA2Calorimeter+;
//#pragma link C++ class TA2Tagger+;
#pragma link C++ class TA2CosmicCal+;
#pragma link C++ class TA2CrystalBall+;
#pragma link C++ class TA2CB+;
#pragma link C++ class TA2TAPS+;
#pragma link C++ class TA2TAPS2008+;
#pragma link C++ class TA2TAPS2009+;
#pragma link C++ class TA2Taps2009LE+;
#pragma link C++ class TA2Taps+;
// Detector classes
#pragma link C++ class TA2PlasticPID+;
#pragma link C++ class TA2CalArray+;
//#pragma link C++ class TA2Ladder+;
#pragma link C++ class TA2TAPS_BaF2+;
#pragma link C++ class TA2TAPS_PbWO4+;
#pragma link C++ class TA2TAPS_Veto+;
#pragma link C++ class TA2CylMWPC+;
#pragma link C++ class TA2FPMicro;
// Utility classes...components of detectors etc
#pragma link C++ class TA2CylStripSven+;
#pragma link C++ class TA2CylWireSven+;
#pragma link C++ class TA2WCLayerSven+;
#pragma link C++ class TA2Event+;
//#pragma link C++ class TA2Particle+;

#endif
