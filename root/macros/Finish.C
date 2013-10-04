// Stuff to do at the end of an analysis run.

#include "includes/physics.h"
#include "includes/functions.h"

void FinishData( TString name)
{
	TString path;

	cout << "End-of-Run macro executing";
	cout << endl;

	path =  "histograms/";
	path.Append( name);
	path.Append( ".root");
	TFile f( path, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	path.Prepend( "All histograms saved to \"");
	path.Append( "\"");
	cout << path;
	cout << endl;
}

void FinishMC( Char_t tgt[2], Char_t cond[32], Int_t eg)
{
	Char_t dir[128];

	sprintf( dir, "histograms/MonteCarlo/%s_%s_%d.root", tgt, cond, eg);

	printf( "End-of-Run macro executing\n");

	TFile f( dir, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	printf( "All histograms saved to %s\n", dir);

	gROOT->ProcessLine( ".L ThetaEff.C");

	ThetaEffBatch( eg, tgt, cond);

	printf( "Theta efficiencies calculated for %s %s %d\n", tgt, cond, eg);
}

void FinishMC2( TString proc, Int_t eg, TString cond)
{
	TString name;

	name = Form( "histograms/MonteCarlo/%s_%d_%s.root", (const char*) proc, eg,
			(const char*) cond);

	printf( "End-of-Run macro executing\n");

	TFile f( name, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	name.Prepend( "All histograms saved to ");
	cout << name << endl;

//	gROOT->ProcessLine( ".L ThetaEff.C");
//	ThetaEffBatch2( tgt, eg, cond);
//	printf( "Theta efficiencies calculated for %s %d %s\n", tgt, eg, cond);
}

void FinishMC3( TString savedir, TString encl, TString process, TString tgt,
		Int_t eg)
{
	TString name;

	name = "End-of-Run macro executing";
	cout << name << endl;

	name = Form( "histograms/MonteCarlo/%s/%s/%s_%s_%d.root",
			(const char*) savedir, (const char*) encl, (const char*) process,
			(const char*) tgt, eg);

	TFile f( name, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	name.Prepend( "All histograms save to ");
	cout << name << endl;
}

void FinishMC4( TString savedir, TString encl, TString process, TString tgt,
		Int_t chan)
{
	TString dir;

	dir = Form( "histograms/MonteCarlo/%s/%s/%s_%s_chan%d.root",
			(const char*) savedir, (const char*) encl, (const char*) process,
			(const char*) tgt, chan);

	printf( "End-of-Run macro executing\n");

	TFile f( dir, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	printf( "All histograms saved to %s\n", (const char*) dir);

//	gROOT->ProcessLine( ".L ThetaEff.C");
//
//	ThetaEff( kFALSE, process, tgt, chan);
//
//	printf( "Efficiencies for process %s, enclosure %s, target %s,"
//				" and tagger channel %d\n", (const char*) process,
//				(const char*) encl, (const char*) tgt, chan);
}

void FinishMC5( TString process, Int_t energy, Int_t chan)
{
	TString name;

	name = Form( "histograms/MonteCarlo/%s_e%d_t%d.root",
			(const char*) process, energy, chan);

	printf( "End-of-Run macro executing\n");

	TFile f( name, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	name.Prepend( "All histograms saved to ");
	cout << name << endl;
}

void FinishTE()
{
	Int_t run;
	TString name;

	run = GetRunDave();

	name = Form( "histograms/TaggEff/TaggEff_%d.root", run);

	cout << "End-of-Run macro executing";
	cout << endl;

	TFile f( name, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	name.Prepend( "All histograms saved to \"");
	name.Append( "\"");
	cout << name;
	cout << endl;
}
