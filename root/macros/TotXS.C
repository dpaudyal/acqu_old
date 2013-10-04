// Macro TotXS.C
//
// This macro calculates the total cross section from the proton missing mass
// as a function of tagger channel converted to incident photon energy.
//
// Adapted from my old version from the June/July 2004 data with CB-TAPS round
// #1.
//
// NOTE:  The tagger channels run from 0 to 351.
//
// DLH		10.12.2009		First "New" Version for December 2008 data.
//

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

typedef struct {

	Int_t egamma;
	Double_t energy;
	Double_t denergy;
	Double_t etag;
	Double_t d_etag;
	Double_t edet[10];
	Double_t d_edet[10];
	Double_t pa_ratio;
	Double_t xs[10];
	Double_t dxs[10];

} TData;
TData tcd[352];

typedef struct {
	Int_t th;
	Double_t theta;
	Double_t dtheta;
	Double_t dom;
} TBins;
TBins tbin[10];

Double_t f_dead_f, f_dead_e;
Int_t broken_ladd_ch[] = { -1};

TString prompt, random, scalers, dtscalers;
prompt = "PhiCMCut2P_v_ThetaCMCut2P_v_TChanCut2P";
random = "PhiCMCut2R_v_ThetaCMCut2R_v_TChanCut2R";
scalers = "SumScalers152to503";
dtscalers = "SumScalers504to535";

// Target Full
TFile full( "histograms/FullPi0.root");
TH3D *hf3dp = (TH3D*)full.Get( prompt);
TH3D *hf3dr = (TH3D*)full.Get( random);
TH1D *hfsc = (TH1D*)full.Get( scalers);
TH1D *hdeadf = (TH1D*)full.Get( dtscalers);

// Target Empty
TFile empty( "histograms/EmptyPi0.root");
TH3D *he3dp = (TH3D*)empty.Get( prompt);
TH3D *he3dr = (TH3D*)empty.Get( random);
TH1D *hesc = (TH1D*)empty.Get( scalers);
TH1D *hdeade = (TH1D*)empty.Get( dtscalers);

void InitXS()
{
	UInt_t i, eg;
	Double_t eff, deff;
	Double_t f_tagg, f_tot;
	TString name;

	gROOT->ProcessLine( ".L ReadParams.C");

	ReadTagEng( "xs/tageng855.dat");

//	ReadSubt( "xs/chan_subt/chan_subt_full.out");
	for ( i = 0; i <= 351; tcd[i++].pa_ratio = 0.0833);

	ReadDetEff();
	ReadTagEff( "xs/eff/tag_eff.out");

	// Full target scaler deadtime correction.
	f_tagg = hdeadf->GetBinContent( 32)/hdeadf->GetBinContent( 31);
	f_tot = hdeadf->GetBinContent( 25)/hdeadf->GetBinContent( 26);
	f_dead_f = f_tagg/f_tot;
	cout << "f_tagg = " << f_tagg;
	cout << "  f_tot = " << f_tot << endl;
	cout << "full tagger scalers deadtime correction = " << f_dead_f << endl;

	// Empty target scaler deadtime correction.
	f_tagg = hdeade->GetBinContent( 32)/hdeade->GetBinContent( 31);
	f_tot = hdeade->GetBinContent( 25)/hdeade->GetBinContent( 26);
	f_dead_e = f_tagg/f_tot;
	cout << "f_tagg = " << f_tagg;
	cout << "  f_tot = " << f_tot << endl;
	cout << "empty tagger scalers deadtime correction = " << f_dead_e << endl;

	// Zero scalers for broken ladder channels.
	Int_t* lch = broken_ladd_ch;
	while ( *lch != -1) {
		hfsc->SetBinContent( *lch, 0);
		hesc->SetBinContent( *lch, 0);
		lch++;
	}

	tbin[0].th = 0;
	tbin[0].theta = 0;
	tbin[0].dtheta = 0;
	tbin[0].dom = 1;
	for ( i = 1; i < 10; i++) {
		Double_t th, dth;
		th = 10*(2*i-1);
		dth = 20;
		tbin[i].th = th;
		tbin[i].theta = th;
		tbin[i].dtheta = dth;
		tbin[i].dom = (2*kPI)*(cos((th-dth/2)*kD2R)-cos((th+dth/2)*kD2R));
	}
}

void TotXS( TString subt_str)
{
	Int_t i, count, chan;
	Double_t eng, x, dx;
	Double_t enn[352], denn[352], xsn[352], dxsn[352];
	Double_t eng1[50], deng1[50], xs1[50], dxs1[50];
	Double_t xx[10], yy[10];
	TString name;

	c1 = new TCanvas( "c1", "Total Cross Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Calculate the channel cross sections
	// For total cross sections, the theta bin is 0 for 0-180 deg.
	for ( i = 0; i <= 351; ChanXS( subt_str, i++, 0));

	name = "xs/tot_xs/tot_xs.out";
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	count = 0;
	for ( chan = 0; chan <= 351; chan++) {
		enn[count] = tcd[chan].energy;
		denn[count] = tcd[chan].denergy/2;
		xsn[count] = tcd[chan].xs[0];
		dxsn[count] = tcd[chan].dxs[0];
		outFile << enn[count] << "  " << xsn[count] << "  " << dxsn[count++]
			<< endl ;
	}
	outFile.close();

	// Plot Results
	gr = new TGraphErrors( count-1, enn, xsn, denn, dxsn);
	gr->SetTitle( "Preliminary #gammap#rightarrowp#pi^{0} Total Cross Section");
//	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 20);
	gr->SetMarkerSize( 1.0);
	gr->SetLineWidth( 3);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle("#sigma (#mub)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser(100,855);
	gr->SetMaximum(400);
	gr->SetMinimum(-5);
	gr->Draw( "AP");

	pt = new TLegend(0.6,0.2,0.85,0.4);
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	pt->SetTextSize(0.04);

	pt->AddEntry( gr, "This Work", "p");

	l2 = new TLine(130,-0.5,130,15);
	l2->SetLineStyle(1);
	l2->SetLineWidth(2);
	l2->Draw();
	l3 = new TLine(130,15,180,15);
	l3->SetLineStyle(1);
	l3->SetLineWidth(2);
	l3->Draw();
	l4 = new TLine(180,-0.5,180,15);
	l4->SetLineStyle(1);
	l4->SetLineWidth(2);
	l4->Draw();
	l5 = new TLine(130,-0.5,180,-0.5);
	l5->SetLineStyle(1);
	l5->SetLineWidth(2);
	l5->Draw();

	l6 = new TLine(130,0,117,185);
	l6->SetLineStyle(2);
	l6->SetLineWidth(2);
	l6->Draw();
	l7 = new TLine(180,0,250,185);
	l7->SetLineStyle(2);
	l7->SetLineWidth(2);
	l7->Draw();

	// Roman's Data
	i = 0;
	name = "xs/previous/tot_xs_rl.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> eng >> x >> dx;
		eng1[i] = eng;
		deng1[i] = 0;
		xs1[i] = x;
		dxs1[i++] = dx;
	}
	inFile.close();

	rl = new TGraphErrors( i-1, eng1, xs1, deng1, dxs1);
	rl->SetMarkerColor( 2);
	rl->SetMarkerStyle( 21);
	rl->SetMarkerSize( 1.0);
	rl->SetLineWidth( 3);
	rl->SetLineColor( 2);
	rl->Draw("P");

	pt->AddEntry( rl, "Leukel - TAPS", "p");
	pt->Draw();

	sp1 = new TPad( "sp1", "Threshold ", 0.13, 0.475, 0.51, 0.875);
	sp1->Draw();
//	sp1->SetFillColor( 40);
	sp1->SetFillColor( 29);
	sp1->SetGrid();
	sp1->GetFrame()->SetFillColor( 21);
	sp1->GetFrame()->SetBorderSize( 12);
	sp1->SetRightMargin(0.01);
	sp1->SetLeftMargin(0.04);
	sp1->SetBottomMargin(0.05);
	sp1->cd();
	gr2 = new TGraphErrors( count-1, enn, xsn, denn, dxsn);
//	gr2->SetTitle( "Threshold Region");
	gr2->SetTitle();
//	gr2->SetMarkerColor( 4);
	gr2->SetMarkerStyle( 20);
	gr2->SetMarkerSize( 1.0);
	gr2->SetLineWidth( 3);
	gr2->GetXaxis()->SetRangeUser(130,170);
	gr2->SetMaximum(5.5);
	gr2->SetMinimum(-0.5);
	gr2->Draw( "AP");
	l1 = new TLine(145,-0.5,145,5.5);
	l1->SetLineStyle(2);
	l1->SetLineColor(2);
	l1->SetLineWidth(3);
	l1->Draw();

	// Axel's Data
	i = 0;
	name = "xs/previous/tot_xs_as.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> eng >> x >> dx;
		eng1[i] = eng;
		deng1[i] = 0;
		xs1[i] = x;
		dxs1[i++] = dx;
	}
	inFile.close();
	as = new TGraphErrors( i-1, eng1, xs1, deng1, dxs1);
	as->SetMarkerColor( 2);
	as->SetMarkerStyle( 22);
	as->SetMarkerSize( 1.0);
	as->SetLineWidth(3);
	as->SetLineColor(2);
	as->Draw("P");

	pt->AddEntry( as, "Schmidt - TAPS", "p");

	TPaveText *pt1 = new TPaveText(145.5,4.2,162.5,4.9);
	pt1->AddText(" p#pi^{0} Threshold");
	pt1->SetBorderSize(0);
	pt1->SetFillColor(29);
	pt1->Draw();

	name = "plots/xstot";
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

void ChanXS( TString subt_str, UInt_t i, UInt_t j) 
{
	Int_t xmin, xmax, ymin, ymax, zmin, zmax;
	Double_t t, fact, emin, emax, sum, n_yld, dn_yld;
	Double_t theta_lo, theta_hi, chan_lo, chan_hi;
	Double_t yield_p, yield_r;
	Double_t fyield, dfyield, fscalers, dfscalers;
	Double_t eyield, deyield, escalers, descalers;
	Double_t r_fe, yield, dyield;

	t = 4.242e23;		// This is in nuclei/cm^2 for the 10-cm cell.

	// The 1e30 converts cm^2 to microbarn.
	fact = 1e30/t;

	// Limits over which to integrate 3D histogram.
	//		- One tagger channel.
	// 	- MMiss region.
	// 	- All CM theta.
	//  NOTE: Bins start at 1 and not 0!!!!

	// Tagger Channel (e.g., channel 0 is bin 1)
	xmin = i+1;
	xmax = i+1;

	// CM Theta (this assumes 9 theta bins)
	if ( j == 0) {
		ymin = 1;
		ymax = 9;
	}
	else {
		ymin = j;
		ymax = j;
	}

	// CM Phi (this assumes 18 phi bins)
	zmin = 1;
	zmax = 18;

	// Full target
	yield_p = hf3dp->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield_r = hf3dr->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	fyield = yield_p - tcd[i].pa_ratio*yield_r;
	dfyield = sqrt( yield_p - Sqr( tcd[i].pa_ratio)*yield_r);
	fscalers = hfsc->GetBinContent( xmin)*f_dead_f;
	dfscalers = sqrt( fscalers);

	// Empty
	yield_p = he3dp->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield_r = he3dr->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	eyield = yield_p - tcd[i].pa_ratio*yield_r;
	deyield = sqrt( yield_p - Sqr( tcd[i].pa_ratio)*yield_r);
//	escalers = hesc->GetBinContent( xmin)*f_dead_e;
	escalers = hesc->GetBinContent( xmin)*f_dead_f;
	descalers = sqrt( escalers);

	if ( ( fscalers != 0) && ( escalers != 0)) {

		if ( subt_str == "full") {
			n_yld = fyield/fscalers;
			dn_yld = 1/fscalers*sqrt( Sqr( dfyield)
					+ Sqr( fyield*dfscalers/fscalers));
			yield = fyield;
			dyield = dfyield;
		}
		else if ( subt_str == "empty") {
			n_yld = eyield/escalers;
			dn_yld = 1/escalers*sqrt( Sqr( deyield)
					+ Sqr( eyield*descalers/escalers));
			yield = eyield;
			dyield = deyield;
		}
		else if ( subt_str == "subt") {
			n_yld = eyield/escalers;
			r_fe = fscalers/escalers;
			yield = fyield - r_fe*eyield;
			dyield = sqrt( Sqr( dfyield) + Sqr( r_fe*deyield));

			n_yld = yield/fscalers;
			dn_yld = 1/fscalers*sqrt( Sqr( dyield)
					+ Sqr( yield*dfscalers/fscalers));
		}
		else {
			cout << "ERROR: must use known subtraction string";
			cout << endl;
			break;
		}

		tcd[i].xs[j] = n_yld*fact/tcd[i].edet[j]/tcd[i].etag/tbin[j].dom;
		tcd[i].dxs[j] = dn_yld*fact/tcd[i].edet[j]/tcd[i].etag/tbin[j].dom;

		cout << "Chan = " << i;
		cout << "  " << yield_p;
		cout << "  " << yield_r;
		cout << "  " << tcd[i].pa_ratio;
		cout << "  Eg = " << tcd[i].energy;
		cout << "  edet = " << tcd[i].edet[j];
		cout << "  etag = " << tcd[i].etag;
		cout << "  dom = " << tbin[j].dom;
		cout << "  Y = " << yield;
		cout << " +/- " << dyield;
		cout << "  S = " << fscalers;
		cout << " +/- " << dfscalers;
		cout << "  xs = " << tcd[i].xs[j];
		cout << " +/- " << tcd[i].dxs[j];
		cout << endl;
	}
}

void ThreshXS( TString subt_str)
{
	Int_t i, count, min, max, chan;
	Double_t eng, x, dx;
	Double_t enn[352], denn[352], xsn[352], dxsn[352];
	Double_t eng1[50], deng1[50], xs1[50], dxs1[50];
	Double_t xx[10], yy[10];
	TString name;

	c1 = new TCanvas( "c1", "Total Cross Sections", 200, 10, 700, 500);
	c1->SetFillColor( 42);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	min = 280;
	max = 320;

	count = 0;
	ofstream outFile( "xs/tot_xs/thresh_xs.out");
	for ( chan = min; chan <= max; chan++) {
		ChanXS( subt_str, chan, 0);
		enn[count] = tcd[chan].energy;
		denn[count] = tcd[chan].denergy/2;
		xsn[count] = tcd[chan].xs[0];
		dxsn[count] = tcd[chan].dxs[0];
		outFile << enn[count] << "  " << xsn[count] << "  " << dxsn[count++]
			<< endl ;
	}
	outFile.close();

	gr = new TGraphErrors( count-1, enn, xsn, denn, dxsn);
	gr->SetTitle( "Preliminary #gammap#rightarrowp#pi^{0} Total Cross Section");
	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 0.5);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle("#sigma (#mub)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser(140,210);
//	gr->SetMaximum(6);
//	gr->SetMinimum(-0.1);
	gr->Draw( "AP");

	// Label for my data
//	xx[0] = 300;
//	yy[0] = 100;
//	TPolyMarker *pm = new TPolyMarker( 1, xx, yy);
//	pm->SetMarkerColor( 4);
//	pm->SetMarkerStyle( 21);
//	pm->SetMarkerSize( 1.0);
//	pm->Draw();
//	TPaveText *pt2 = new TPaveText(305,90,345,110);
//	pt2->AddText("This Work");
//	pt2->SetBorderSize(0);
//	pt2->SetFillColor(42);
//	pt2->Draw();

	// Label for Axel's data
//	xx[0] = 300;
//	yy[0] = 60;
//	TPolyMarker *pm = new TPolyMarker( 1, xx, yy);
//	pm->SetMarkerColor( 6);
//	pm->SetMarkerStyle( 21);
//	pm->SetMarkerSize( 1.0);
//	pm->Draw();
//	TPaveText *pt2 = new TPaveText(305,50,370,70);
//	pt2->AddText("Schmidt - TAPS");
//	pt2->SetBorderSize(0);
//	pt2->SetFillColor(42);
//	pt2->Draw();

	// Axel's Data
	i = 0;
	ifstream inFile( "xs/previous/tot_xs_as.dat");
	while( !inFile.eof()) {
		inFile >> eng >> x >> dx;
		eng1[i] = eng;
		deng1[i] = 0;
		xs1[i] = x;
		dxs1[i++] = dx;
	}
	inFile.close();
	as = new TGraphErrors( i-1, eng1, xs1, deng1, dxs1);
	as->SetMarkerColor( 6);
	as->SetMarkerStyle( 21);
	as->SetMarkerSize( 0.7);
	as->Draw("P");

//	TPaveText *pt1 = new TPaveText(147,11.5,165,13.5);
//	pt1->AddText(" p#pi^{0} Threshold");
//	pt1->SetBorderSize(0);
//	pt1->SetFillColor(40);
//	pt1->SetFillStyle(4000);
//	pt1->Draw();

	name = "plots/xstot_thresh";
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

void DeltaXS( TString subt_str)
{
	Int_t i, count, min, max, chan;
	Double_t eng, x, dx;
	Double_t enn[200], denn[200], xsn[200], dxsn[200];
	Double_t eng1[50], deng1[50], xs1[50], dxs1[50];
	Double_t xx[10], yy[10];
	TString name;

	c1 = new TCanvas( "c1", "Total Cross Sections", 200, 10, 700, 500);
	c1->SetFillColor( 42);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	min = 50;
	max = 150;

	count = 0;
	ofstream outFile( "xs/tot_xs/delta_xs.out");
	for ( chan = min; chan <= max; chan++) {
		ChanXS( subt_str, chan);
		enn[count] = en[chan];
		denn[count] = 0;
		xsn[count] = xs[chan];
		dxsn[count] = dxs[chan];
		outFile << enn[count] << "  " << xsn[count] << "  " << dxsn[count++]
			<< endl ;
	}
	outFile.close();

	gr = new TGraphErrors( count-1, enn, xsn, denn, dxsn);
	gr->SetTitle( "Preliminary #gammap#rightarrowp#pi^{0} Total Cross Section");
	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 0.5);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle("#sigma (#mub)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser(265,375);
	gr->SetMaximum(350);
	gr->SetMinimum(150);
	gr->Draw( "AP");

	// Roman's Data
	i = 0;
	ifstream inFile( "xs/previous/tot_xs_rl.dat");
	while( !inFile.eof()) {
		inFile >> eng >> x >> dx;
		eng1[i] = eng;
		deng1[i] = 0;
		xs1[i] = x;
		dxs1[i++] = dx;
	}
	inFile.close();

	rl = new TGraphErrors( i-1, eng1, xs1, deng1, dxs1);
	rl->SetMarkerColor( 2);
	rl->SetMarkerStyle( 21);
	rl->SetMarkerSize( 0.5);
	rl->Draw("P");

	// Label for my data
//	xx[0] = 300;
//	yy[0] = 100;
//	TPolyMarker *pm = new TPolyMarker( 1, xx, yy);
//	pm->SetMarkerColor( 4);
//	pm->SetMarkerStyle( 21);
//	pm->SetMarkerSize( 1.0);
//	pm->Draw();
//	TPaveText *pt2 = new TPaveText(305,90,345,110);
//	pt2->AddText("This Work");
//	pt2->SetBorderSize(0);
//	pt2->SetFillColor(42);
//	pt2->Draw();

	name = "plots/xstot_delta";
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

void DiffXS( TString subt_str, UInt_t i)
{
	UInt_t j, k;
	Double_t tth[9], dtth[9], xxs[9], dxxs[9];
	TString name;

	cout << "Energy = " << tcd[i].energy << endl;

	// Output results to a file
	name = Form( "xs/diff_xs/diff_xs_%d.out", tcd[i].egamma);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	k = 0;
	for (  j = 1; j <= 9; j++) {
		ChanXS( subt_str, i, j);
		name = Form( "%3d  %6.4f  %5.4f", tbin[j].th, tcd[i].xs[j],
				tcd[i].dxs[j]);
		outFile << name << endl;
		tth[k] = tbin[j].theta;
		dtth[k] = 0;
		xxs[k] = tcd[i].xs[j];
		dxxs[k++] = tcd[i].dxs[j];
	}
	outFile.close();

	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Plot Results
	gr = new TGraphErrors( 9, tth, xxs, dtth, dxxs);
	name = Form( "#gammap#rightarrowp#pi^{0} Differential Cross Section"
			" E_{#gamma} = %5.1f MeV", tcd[i].energy);
	gr->SetTitle( name);
	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 21);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle("#theta^{*}_{#pi^{0}} (MeV)");
	gr->GetYaxis()->SetTitle("d#sigma/d#Omega (#mub/sr)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
//	gr->GetXaxis()->SetRangeUser(100,402);
//	gr->SetMaximum(350);
//	gr->SetMinimum(-20);
	gr->Draw( "AP");

//	c1->cd();
//	TPad *npad = new TPad( "npad", "Transparent Pad", 0, 0, 1, 1);
//	npad->SetFillStyle( 4000);
//	npad->Draw();
//	npad->cd();
//
//	// "Preliminary"
////	TPaveLabel *pl = new TPaveLabel( 100, 50, 380, 380, "PRELIMINARY");
//	TPaveLabel *pl = new TPaveLabel( 0, 0, 1, 1, "PRELIMINARY");
//	pl->SetTextAngle(30);
//	pl->SetTextColor(14);
////	pl->SetTextFont(82);
//	pl->SetBorderSize(0);
//	pl->SetFillStyle(4000);
//	pl->Draw();

	name = Form( "plots/diff_xs_%s_%d", (const char*) subt_str, tcd[i].egamma);
//	name.Append( ".eps");
	name.Append( ".pdf");
//	c1->Print( name);

}

void DiffComp( Int_t eg)
{
	Int_t i, ct;
	Double_t th, x, dx, tth[30], dtth[30], xxs[30], dxxs[30];
	Double_t xx[1], dxx[1], yy[1], dyy[1], max;
	Double_t x1, x2, y1, y2;
	TString name;

	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
//	c1->SetFillColor( 38);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// My results
	i = 0;
	name = Form( "xs/diff_xs/diff_xs_%d.out", eg);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	ifstream inFile( name);
	max = 0;
	while( !inFile.eof()) {
		inFile >> th >> x >> dx;
		tth[i] = th;
		dtth[i] = 0;
		xxs[i] = x;
		if ( max < xxs[i]) max = xxs[i];
		dxxs[i++] = dx;
	}
	inFile.close();
	ct = i-1;

	// Plot Results
	gr = new TGraphErrors( ct, tth, xxs, dtth, dxxs);
	name = Form( "Preliminary #gammap#rightarrowp#pi^{0} "
			"Differential Cross Section for E_{#gamma} = %d MeV", eg);
	gr->SetTitle( name);
	gr->SetMarkerColor( 4);
	gr->SetMarkerSize( 1.2);
	gr->SetMarkerStyle( 21);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta^{cm}_{#pi^{0}} (deg)");
	gr->GetYaxis()->SetTitle( "d#sigma/d#Omega (#mub/sr)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	max *= 1.4;
	gr->GetYaxis()->SetRangeUser( 0, max);
	gr->Draw( "AP");

	// Previous results
	if ( eg == 156) name = Form( "xs/previous/dxs_%d.dat", eg-1);
	else if ( eg == 166) name = Form( "xs/previous/dxs_%d.dat", eg+1);
	else name = Form( "xs/previous/dxs_%d.dat", eg);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while( !inFile.eof()) {
		inFile >> th >> x >> dx;
		tth[i] = th;
		dtth[i] = 0;
		xxs[i] = x;
		dxxs[i++] = dx;
	}
	inFile.close();
	ct = i-1;

	// Plot Results
	gr1 = new TGraphErrors( ct, tth, xxs, dtth, dxxs);
	gr1->SetMarkerColor( 2);
	gr1->SetMarkerSize( 1.2);
	gr1->SetLineWidth( 2);
	gr1->SetLineColor( 2);
	gr1->SetMarkerStyle( 20);
	gr1->Draw( "Psame");

	pt = new TLegend(0.6,0.15,0.8,0.30);
	pt->SetTextSize(0.04);
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	pt->AddEntry( gr, "This Work", "p");
	pt->AddEntry( gr1, "Schmidt - TAPS", "p");
	pt->Draw();

	name = Form( "plots/dxs_comp_%d", eg);
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}
