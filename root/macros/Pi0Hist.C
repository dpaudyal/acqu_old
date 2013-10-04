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
	Double_t xs;
	Double_t dxs;

} TChanData;
TChanData tcd[352];

typedef struct {
	Int_t mid;
	Double_t lo;
	Double_t hi;
} ThetaBins;
ThetaBins theta[10];

Double_t pa = 0.0833;

//TFile full( "histograms/FullPi0.root");
TFile full( "histograms/Full/ARH.root");
TH1D *hfsc = (TH1D*)full.Get( "SumScalers152to503");
TH1D *hdf = (TH1D*)full.Get( "SumScalers504to535");

TFile empty( "histograms/EmptyPi0.root");
TH1D *hesc = (TH1D*)empty.Get( "SumScalers152to503");
TH1D *hde = (TH1D*)empty.Get( "SumScalers504to535");

TH1D* fhist;
TH1D* ehist;
TH1D* shist;

void Init()
{
	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");

	theta[0].mid = 0;
	theta[0].lo = 0;
	theta[0].hi = 180;

	for ( UInt_t i = 1; i <= 9; i++) {
		theta[i].mid = 20*i - 10;
		theta[i].lo = 20*(i-1);
		theta[i].hi = 20*i;
	}
}

void Pi0HistChan( TString hist = "MMiss", UInt_t i = 308, UInt_t j = 0,
		Bool_t save = kFALSE)
{
	UInt_t rebin;
	Int_t th, dth;
	Double_t eg, deg;
	Double_t x1, x2, x3, left, right;
	Double_t x[2], y[2];
	Double_t q_pi, T_pi;
	Double_t TggCutOffset;

	TggCutOffset = 3;

	TString name;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) {
		delete c1;
	}
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 300, 20, 1200, 800);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);

	eg = tcd[i].energy;
	deg = tcd[i].denergy;

	th = theta[j].mid;
	dth = 20;

	if ( hist == "MMiss")
	{
		x1 = Linear( 145, 933, 400, 880, eg);
		x2 = Linear( 145, 945, 400, 1000, eg);
		x3 = Linear( 145, 918, 400, 880, eg);

		left = x1 - 30;
		right = x2 + 20;

		rebin = 1;
	}
	else if ( hist == "TGG")
	{
		// Proton Opening Angle
		q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		x1 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R - TggCutOffset;

		// 12-C Opening Angle
		q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		x2 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

		x3 = -1;

		left = 0;
		right = 180;

		rebin = 4;
	}
	else if ( hist == "KECM")
	{
		q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		
		x1 = T_pi;

		left = 0;
		right = x1 + 20;

		rebin = 1;
	}

	l = new TLine( left, 0, right, 0);
	l->SetLineStyle( 2);

	Pi0Hist3D( i, j, hist);

	c1->cd( 1);
	name = Form( "full%d", i);
	TH1D *h1 = (TH1D*) fhist->Clone( name);
	h1->SetMarkerStyle( 21);
	h1->SetMarkerColor( 4);
	h1->SetLineWidth( 3);
	h1->SetLineColor( 4);
	h1->SetTitle();
	h1->GetXaxis()->SetTitleOffset( 1.1);
	h1->GetXaxis()->SetLabelSize( 0.03);
	h1->GetXaxis()->CenterTitle();
	if ( hist == "MMiss") name = "M_{miss} (MeV)";
	else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
	else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
	h1->GetXaxis()->SetTitle( name);
	h1->Rebin( rebin);
	h1->GetXaxis()->SetRangeUser( left, right);
	Double_t max = h1->GetMaximum();
	max *= 1.4;
	h1->SetMaximum( max);
	h1->Draw();

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineStyle( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineStyle( 2);
	l2->Draw();
	l3 = new TLine( x3, 0, x3, max);
	l3->SetLineStyle( 2);
	if ( hist == "MMiss") l3->Draw();

	l->Draw();

	name = Form( "empty%d", i);
	TH1D *h2 = (TH1D*) ehist->Clone( name);
	h2->SetMarkerStyle( 22);
	h2->SetMarkerColor( 2);
	h2->SetLineWidth( 3);
	h2->SetLineColor( 2);
	h2->Rebin( rebin);
	h2->Draw( "same");

	x[0] = 0.15;
	x[1] = 0.40;
	y[0] = 0.6;
	y[1] = 0.75;
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.06);
	tl->SetTextFont( 132);
	tl->AddEntry( h1, "Full");
	tl->AddEntry( h2, "Empty");
	tl->Draw();

	name = Form( "E_{#gamma} = %5.1f #pm %3.1f MeV", eg, deg/2);
	pl = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
	pl->SetFillStyle( 0);
	pl->SetBorderSize( 0);
	pl->SetTextSize( 0.5);
	pl->SetTextFont( 132);
	pl->Draw();

	c1->cd( 2);
	name = Form( "subt%d", i);
	TH1D *h3 = (TH1D*) shist->Clone( name);
	h3->SetMarkerStyle( 20);
	h3->SetLineWidth( 3);
	h3->SetTitle();
	h3->GetXaxis()->SetTitleOffset( 1.1);
	h3->GetXaxis()->SetLabelSize( 0.03);
	h3->GetXaxis()->CenterTitle();
	if ( hist == "MMiss") name = "M_{miss} (MeV)";
	else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
	else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
	h3->GetXaxis()->SetTitle( name);
	h3->Rebin( rebin);
	h3->GetXaxis()->SetRangeUser( left, right);
	Double_t max = h3->GetMaximum();
	max *= 1.4;
	h3->SetMaximum( max);
	h3->Draw();

	if ( j == 0) name = "#theta_{CM} = 0 - 180 deg";
	else name = Form( "#theta_{CM} = %d #pm %d deg", th, dth/2);
	pl2 = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
	pl2->SetFillStyle( 0);
	pl2->SetBorderSize( 0);
	pl2->SetTextSize( 0.5);
	pl2->SetTextFont( 132);
	pl2->Draw();

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineStyle( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineStyle( 2);
	l2->Draw();
	l3 = new TLine( x3, 0, x3, max);
	l3->SetLineStyle( 2);
	if ( hist == "MMiss") l3->Draw();

	l->Draw();
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.06);
	tl->SetTextFont( 132);
	tl->AddEntry( h3, "Subtracted");
	tl->Draw();

	name = Form( "plots/Pi0_%s_chan%d_t%d.pdf", (const char*) hist, i, th);
	if ( save == kTRUE) c1->Print( name);
}

void Pi0Hist4( TString hist = "MMiss", UInt_t lo = 308, UInt_t tbin = 0,
		Bool_t save = kFALSE)
{
	UInt_t i, j, rebin;
	Int_t th, dth;
	Double_t eg, deg;
	Double_t x1, x2, x3, left, right;
	Double_t x[2], y[2];
	Double_t q_pi, T_pi;
	Double_t TggCutOffset;
	TString name;

	TggCutOffset = 3;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) {
		delete c1;
	}
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 300, 20, 1200, 600);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 5, 2);

	th = theta[tbin].mid;
	dth = 20;

	j = 1;
	for ( i = lo; i > lo-5; i--)
	{
		eg = tcd[i].energy;
		deg = tcd[i].denergy;

		if ( hist == "MMiss")
		{
			x1 = Linear( 145, 933, 400, 880, eg);
			x2 = Linear( 145, 945, 400, 1000, eg);
			x3 = Linear( 145, 918, 400, 880, eg);

			left = x1 - 30;
			right = x2 + 20;

			rebin = 1;
		}
		else if ( hist == "TGG")
		{
			// Proton Opening Angle
			q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x1 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R - TggCutOffset;

			// 12-C Opening Angle
			q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x2 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

			x3 = -1;

			left = 0;
			right = 180;

			rebin = 4;
		}
		else if ( hist == "KECM")
		{
			// Proton
			q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x1 = T_pi;

			// 12-C
			q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x2 = T_pi;

			x3 = -1;

			left = 0;
			right = x2 + 20;

			rebin = 1;
		}

		l = new TLine( left, 0, right, 0);
		l->SetLineStyle( 2);

		Pi0Hist3D( i, tbin, hist);

		c1->cd( j);
		name = Form( "full%d", i);
		TH1D *h1 = (TH1D*) fhist->Clone( name);
		h1->SetMarkerStyle( 21);
		h1->SetMarkerColor( 2);
		h1->SetLineWidth( 3);
		h1->SetLineColor( 2);
		h1->SetTitle();
		h1->GetXaxis()->SetTitleOffset( 1.1);
		h1->GetXaxis()->SetLabelSize( 0.03);
		h1->GetXaxis()->CenterTitle();
		if ( hist == "MMiss") name = "M_{miss} (MeV)";
		else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
		else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
		h1->GetXaxis()->SetTitle( name);
		h1->GetXaxis()->SetRangeUser( left, right);
		h1->Rebin( rebin);
		Double_t max = h1->GetMaximum();
		max *= 1.4;
		h1->SetMaximum( max);
		h1->Draw();

		l1 = new TLine( x1, 0, x1, max);
		l1->SetLineStyle( 2);
		l1->Draw();
		l2 = new TLine( x2, 0, x2, max);
		l2->SetLineStyle( 2);
		l2->Draw();
		l3 = new TLine( x3, 0, x3, max);
		l3->SetLineStyle( 2);
//		l3->Draw();

		l->Draw();
		name = Form( "E_{#gamma} = %5.1f #pm %3.1f MeV", eg, deg/2);
		pl = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
		pl->SetFillStyle( 0);
		pl->SetBorderSize( 0);
		pl->SetTextSize( 0.5);
		pl->SetTextFont( 132);
		pl->Draw();
		if ( tbin == 0) name = "#theta_{CM} = 0 - 180 deg";
		else name = Form( "#theta_{CM} = %d #pm %d deg", th, dth/2);
		pl2 = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
		pl2->SetFillStyle( 0);
		pl2->SetBorderSize( 0);
		pl2->SetTextSize( 0.5);
		pl2->SetTextFont( 132);
		if ( j == 6) pl2->Draw();

		name = Form( "empty%d", i);
		TH1D *h2 = (TH1D*) ehist->Clone( name);
		h2->SetMarkerStyle( 22);
		h2->SetMarkerColor( 4);
		h2->SetLineWidth( 3);
		h2->SetLineColor( 4);
		h2->Rebin( rebin);
		h2->Draw( "same");

		x[0] = 0.15;
		x[1] = 0.40;
		y[0] = 0.6;
		y[1] = 0.75;
		tl = new TLegend( x[0], y[0], x[1], y[1]);
		tl->SetFillStyle( 0);
		tl->SetBorderSize( 0);
		tl->SetTextSize( 0.06);
		tl->SetTextFont( 132);
		tl->AddEntry( h1, "Full");
		tl->AddEntry( h2, "Empty");
		if ( j == 1) tl->Draw();

		c1->cd( j+5);
		name = Form( "subt%d", i);
		TH1D *h3 = (TH1D*) shist->Clone( name);
		h3->SetMarkerStyle( 20);
		h3->SetLineWidth( 3);
		h3->SetTitle();
		h3->GetXaxis()->SetTitleOffset( 1.1);
		h3->GetXaxis()->SetLabelSize( 0.03);
		h3->GetXaxis()->CenterTitle();
		if ( hist == "MMiss") name = "M_{miss} (MeV)";
		else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
		else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
		h3->GetXaxis()->SetTitle( name);
		h3->GetXaxis()->SetRangeUser( left, right);
		h3->Rebin( rebin);
		Double_t max = h3->GetMaximum();
		max *= 1.4;
		h3->SetMaximum( max);
		h3->Draw();

		l1 = new TLine( x1, 0, x1, max);
		l1->SetLineStyle( 2);
		l1->Draw();
		l2 = new TLine( x2, 0, x2, max);
		l2->SetLineStyle( 2);
		l2->Draw();
		l3 = new TLine( x3, 0, x3, max);
		l3->SetLineStyle( 2);
//		l3->Draw();

		l->Draw();
		tl = new TLegend( x[0], y[0], x[1], y[1]);
		tl->SetFillStyle( 0);
		tl->SetBorderSize( 0);
		tl->SetTextSize( 0.06);
		tl->SetTextFont( 132);
		tl->AddEntry( h3, "Subtracted");
		if ( j == 1) tl->Draw();

		j++;
	}

	name = Form( "plots/Pi0_%s_chan%d-%d_t%d.pdf", (const char*) hsit, lo,
			lo-4, theta[tbin].mid);
	if ( save == kTRUE) c1->Print( name);
}

void Pi0HistSimComp( TString tgt = "subt", TString hist = "MMiss",
		UInt_t i = 308, UInt_t j = 0, Bool_t save = kFALSE)
{
	UInt_t tbin_lo, tbin_hi, rebin;
	Int_t th, dth;
	Double_t eg, deg;
	Double_t x1, x2, x3, left, right;
	Double_t x[2], y[2];
	Double_t q_pi, T_pi;
	Double_t TggCutOffset;
	Char_t tchar;
	TString name;

	TggCutOffset = 3;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) {
		delete c1;
	}
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 300, 20, 700, 500);
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	eg = tcd[i].energy;
	deg = tcd[i].denergy;

	th = theta[j].mid;
	dth = 20;

	if ( hist == "MMiss")
	{
		x1 = Linear( 145, 933, 400, 880, eg);
		x2 = Linear( 145, 945, 400, 1000, eg);
		x3 = Linear( 145, 918, 400, 880, eg);

		left = x1 - 30;
		right = x2 + 20;

		rebin = 1;
	}
	else if ( hist == "TGG")
	{
		// Proton Opening Angle
		q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		x1 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R - TggCutOffset;

		// 12-C Opening Angle
		q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		x3 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

		x2 = 180;

		left = 0;
		right = 180;

		rebin = 4;
	}
	else if ( hist == "KECM")
	{
		if ( tgt == "subt")
			q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
		else if ( tgt == "empty")
			q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);

		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;

		x1 = 0;
		x2 = T_pi + 20;
		x3 = -1;

		left = 0;
		right = x2 + 20;

		rebin = 1;
	}

	l = new TLine( left, 0, right, 0);
	l->SetLineStyle( 2);

	Pi0Hist3D( i, j, hist);

	l->Draw();

	name = Form( "%s%d", (const char*) tgt, i);
	if ( tgt == "subt") TH1D *h3 = (TH1D*) shist->Clone( name);
	else if ( tgt == "empty") TH1D *h3 = (TH1D*) ehist->Clone( name);
	h3->SetMarkerStyle( 20);
	h3->SetMarkerSize( 1.5);
	h3->SetLineWidth( 4);
	h3->SetTitle();
	h3->GetXaxis()->SetTitleOffset( 1.2);
	h3->GetXaxis()->SetLabelSize( 0.03);
	h3->GetXaxis()->CenterTitle();
	if ( hist == "MMiss") name = "M_{miss} (MeV)";
	else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
	else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
	h3->GetXaxis()->SetTitle( name);
	h3->Rebin( rebin);
	h3->GetXaxis()->SetRangeUser( left, right);
	Double_t max = h3->GetMaximum();
	max *= 1.4;
	h3->SetMaximum( max);
	h3->Draw();

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineStyle( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineStyle( 2);
	if ( hist == "MMiss") l2->Draw();
	l3 = new TLine( x3, 0, x3, max);
	l3->SetLineStyle( 2);
	l3->Draw();

	l->Draw();

	x[0] = 0.15;
	x[1] = 0.40;
	y[0] = 0.6;
	y[1] = 0.75;
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.06);
	tl->SetTextFont( 132);
	tl->AddEntry( h3, "Data");

	name = Form( "E_{#gamma} = %5.1f #pm %3.1f MeV", eg, deg/2);
	pl = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
	pl->SetFillStyle( 0);
	pl->SetBorderSize( 0);
	pl->SetTextSize( 0.5);
	pl->SetTextFont( 132);
	pl->Draw();

	if ( j == 0) name = "#theta_{CM} = 0 - 180 deg";
	else name = Form( "#theta_{CM} = %d #pm %d deg", th, dth/2);
	pl2 = new TPaveLabel( 0.55, 0.90, 0.85, 1.00, name, "NDC");
	pl2->SetFillStyle( 0);
	pl2->SetBorderSize( 0);
	pl2->SetTextSize( 0.5);
	pl2->SetTextFont( 132);
	pl2->Draw();

	Int_t bin_lo = h3->GetXaxis()->FindBin( x1);
	Int_t bin_hi = h3->GetXaxis()->FindBin( x2);
	Double_t cts1 = h3->Integral( bin_lo, bin_hi);

	if ( tgt == "subt") tchar = 'p';
	else if ( tgt == "empty") tchar = 'w';
	name = Form( "histograms/MonteCarlo/pi0/pi0_%c_chan%d.root", tchar, i);
	sim = new TFile( name);
	name = "THR_" + hist + "P_v_ThetaCMP";
	TH2D *h2sim = (TH2D*)sim->Get( name);
	tbin_lo = h2sim->GetXaxis()->FindBin( theta[j].lo);
	tbin_hi = h2sim->GetXaxis()->FindBin( theta[j].hi)-1;
	h2sim->GetXaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *h1sim = h2sim->ProjectionY( "simhist");
	h1sim->Rebin( rebin);

	Double_t cts2 = h1sim->Integral( bin_lo, bin_hi);

	Double_t fact = cts1/cts2;

	h1sim->Scale( fact);
	h1sim->SetFillStyle( 3004);
	h1sim->SetFillColor( 16);
	h1sim->SetLineWidth( 2);
	h1sim->Draw( "same");

	tl->AddEntry( h1sim, "Monte Carlo");
	tl->Draw();

	name = Form( "plots/Pi0_%c_%s_SimComp_chan%d_t%d.pdf", tchar,
			(const char*) hist, i, th);
	if ( save == kTRUE) c1->Print( name);
}

void Pi0Hist3D( UInt_t chan, UInt_t tbin, TString var_y)
{
	UInt_t cbin, tbin_lo, tbin_hi;
	TString name;
	TString var_x, var_z;
	TString hprompt, hrand;

	var_x = "TChan";
	var_z = "ThetaCM";

	hprompt = var_z + "P_v_" + var_y + "P_v_" + var_x + "P";
	hrand = var_z + "R_v_" + var_y + "R_v_" + var_x + "R";

	// Target FULL
	// Prompt
	TH3D *hP = (TH3D*)full.Get( hprompt);
	cbin = hP->GetXaxis()->FindBin( chan);
	hP->GetXaxis()->SetRange( cbin, cbin);
	tbin_lo = hP->GetZaxis()->FindBin( theta[tbin].lo);
	tbin_hi = hP->GetZaxis()->FindBin( theta[tbin].hi)-1;
	hP->GetZaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *hPf_y = hP->Project3D( "y");

	// Random
	TH3D *hR = (TH3D*)full.Get( hrand);
	hR->GetXaxis()->SetRange( cbin, cbin);
	hR->GetZaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *hRf_y = hR->Project3D( "y2");
	hRf_y->Scale( pa);

	// Subtracted
	hSf = (TH1D*) hPf_y->Clone( "FullSubt");
	hSf->Sumw2();
	hSf->Add( hRf_y, -1.0);

	// Target EMPTY
	// Prompt
	TH3D *hP = (TH3D*)empty.Get( hprompt);
	hP->GetXaxis()->SetRange( cbin, cbin);
	hP->GetZaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *hPe_y = hP->Project3D( "y");

	// Random
	TH3D *hR = (TH3D*)empty.Get( hrand);
	hR->GetXaxis()->SetRange( cbin, cbin);
	hR->GetZaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *hRe_y = hR->Project3D( "y2");
	hRe_y->Scale( pa);

	// Subtracted
	hSe = (TH1D*) hPe_y->Clone( "EmptySubt");
	hSe->Sumw2();
	hSe->Add( hRe_y, -1.0);

	// COMPLETELY SUBTRACTED
	Double_t f_tagg_f, f_tot_f, f_dead_f, fscal;
	Double_t f_tagg_e, f_tot_e, f_dead_e, escal;
	Double_t factor;

	f_tagg_f = hdf->GetBinContent( 32)/hdf->GetBinContent( 31);
	f_tot_f = hdf->GetBinContent( 25)/hdf->GetBinContent( 26);
	f_dead_f = f_tagg_f/f_tot_f;
//	fscal = hfsc->GetBinContent( chan)*f_dead_f;
	fscal = hfsc->GetBinContent( chan);

	f_tagg_e = hde->GetBinContent( 32)/hde->GetBinContent( 31);
	f_tot_e = hde->GetBinContent( 25)/hde->GetBinContent( 26);
	f_dead_e = f_tagg_e/f_tot_e;
//	escal = hesc->GetBinContent( chan)*f_dead_e;
	escal = hesc->GetBinContent( chan);

	factor = escal/fscal;
	hSf->Scale( factor);
	hS = (TH1D*) hSf->Clone( "Subtracted");
	hS->Add( hSe, -1.0);

	fhist = (TH1D*) hSf->Clone( "Full");
	ehist = (TH1D*) hSe->Clone( "Empty");
	shist = (TH1D*) hS->Clone( "Subt");
}

void M2gChan( UInt_t chan = 250, TString tgt = "full")
{
	UInt_t left, right, bin_lo, bin_hi;
	Double_t eg, deg;
	TString name;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;

	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 200, 20, 700, 400);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);

	eg = tcd[chan].energy;
	deg = tcd[chan].denergy;

	name = Form( "Egamma = %5.1f +/- %3.1f MeV", eg, deg/2);
	cout << name << endl;

	left = 0;
	right = 200;

	// Prompt
//	name = "THR_M2g_v_TAGG_FPD_HitsPrompt";
	name = "THR_M2gTaggP_v_TChanP";
	if ( tgt == "full") TH2D *hp2d = (TH2D*)full.Get( name);
	else if ( tgt == "empty") TH2D *hp2d = (TH2D*)empty.Get( name);
	bin_lo = hp2d->GetXaxis()->FindBin( chan);
	bin_hi = hp2d->GetXaxis()->FindBin( chan);
	hp2d->GetXaxis()->SetRange( bin_lo, bin_hi);
	TH1D *hp = hp2d->ProjectionY( "M2gP");

	// Random
//	name = "THR_M2g_v_TAGG_FPD_HitsRand";
	name = "THR_M2gTaggR_v_TChanR";
	if ( tgt == "full") TH2D *hr2d = (TH2D*)full.Get( name);
	else if ( tgt == "empty") TH2D *hr2d = (TH2D*)empty.Get( name);
	hr2d->GetXaxis()->SetRange( bin_lo, bin_hi);
	TH1D *hr = hr2d->ProjectionY( "M2gR");

	c1->cd( 1);
	// Prompt
	hp->SetLineWidth( 3);
	hp->SetLineColor( 4);
	name = Form( "2#gamma Invariant Mass for E_{#gamma} = %5.1f #pm %3.1f MeV",
			eg, deg/2);
	hp->SetTitle( name);
	hp->GetXaxis()->SetTitleOffset( 1.1);
	hp->GetXaxis()->SetLabelSize( 0.03);
	hp->GetXaxis()->CenterTitle();
	hp->GetXaxis()->SetTitle( "M_{#gamma#gamma} (MeV)");
	hp->GetXaxis()->SetRangeUser( left, right);
	hp->DrawCopy();

	// Random
	hr->Scale( pa);
	hr->SetLineWidth( 3);
	hr->SetLineColor( 2);
	hr->Draw( "same");

	// Subtracted
	hp->Sumw2();
	hp->Add( hr, -1.0);
	hp->SetMarkerStyle( 21);
	hp->SetMarkerColor( 1);
	hp->SetLineWidth( 3);
	hp->SetLineColor( 1);
	c1->cd(2);
	hp->SetTitle();
	hp->Draw( "E");
	hp->GetXaxis()->SetRangeUser( left, right);

	l = new TLine( left, 0, right, 0);
	l->SetLineStyle( 2);
	l->Draw();
}
