/*
 *		Asymmetry.C
 *
 *		New asymmetry macro using more intelligent organization.
 *
 *		2010.03.04		DLH		First Version
 *
 */

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

typedef struct {
	Int_t egamma;
	Double_t energy;
	Double_t denergy;
} TData;
TData tcd[352];

// Full
TFile fperpFile( "histograms/Full/FullPerp.root");
TFile fparaFile( "histograms/Full/FullPara.root");

// Empty
TFile eperpFile( "histograms/Empty/EmptyPerp.root");
TFile eparaFile( "histograms/Empty/EmptyPara.root");

TH3D* hS;
TH1D* hsc;
TH2D* hA;

UInt_t tb_lo[] = { 0, 0, 20, 40, 60, 80, 100, 120, 140, 160};
UInt_t tb_hi[] = { 180, 20, 40, 60, 80, 100, 120, 140, 160, 180};

void InitAsym()
{
	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");
}

void AsymBin( UInt_t tbin = 0, UInt_t chan_lo = 300, UInt_t chan_hi = 300,
		TString subt = "subt")
{
	UInt_t tbin_lo, tbin_hi;
	UInt_t theta;
	Double_t eg, deg;
	Double_t par[6], red_chisq;
	TString name;

	TH1D* asym;

	theta = 20*tbin - 10;

	Get2DHist( chan_lo, chan_hi, subt);
	name = Form( "asym_tbin%d", tbin);
	asym = (TH1D*) hA->ProjectionX( name, tb_lo[tbin], tb_hi[tbin]);

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "Data", 400, 10, 700, 500);
	asym->Draw();

	eg = WeightedEg( chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy/2) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy/2))/2;
	if ( tbin == 0)
		name = Form( "E_{#gamma} = %5.1f #pm %4.1f MeV  #theta = 0-180 deg",
				eg, deg);
	else
		name = Form( "E_{#gamma} = %5.1f #pm %4.1f MeV  #theta = %d #pm 10 deg",
				eg, deg, theta);
	asym->SetTitle( name);
	asym->GetXaxis()->SetTitle( "#phi (deg)");
	asym->GetXaxis()->CenterTitle();
	asym->GetYaxis()->SetTitle( "p_{#gamma}#Sigma cos(2#phi + #phi_{0})");
	asym->GetYaxis()->CenterTitle();

	TF1 *f1 = new TF1( "f1", "[0]*cos(2*x*0.01745+[1])", -180, 180);
	asym->Fit( "f1", "REMQ");
	par[0] = f1->GetParameter(0);
	par[1] = f1->GetParError(0);
	par[2] = f1->GetParameter(1);
	par[3] = f1->GetParError(1);
	red_chisq = f1->GetChisquare()/17;
	if ( par[2] < 0) {
		par[0] *= -1;
		par[2] *= -1;
	}
	cout << par[0];
	cout << " +/- " << par[1];
	cout << "   " << par[2];
	cout << " +/- " << par[3];
	cout << "   " << red_chisq;
	cout << endl;

	pt = new TPaveText( 0.5, 0.60, 0.8, 0.85, "NDC");
	pt->SetTextAlign( 12);
	pt->SetTextSize( 0.05);
	pt->SetTextFont( 132);
	pt->SetBorderSize( 0);
	pt->SetFillColor( 0);
	name = Form( "p_{#gamma}#Sigma = %6.4f #pm %6.4f", par[0], par[1]);
	pt->AddText( name);
	name = Form( "#phi_{0} = %5.1f #pm %4.1f deg", par[2]/kD2R, par[3]/kD2R);
	pt->AddText( name);
	name = Form( "#chi^{2}/n_{d.o.f.} = %4.2f", red_chisq);
	pt->AddText( name);
	pt->Draw();

	name = Form( "plots/AsymBin_%s_%d-%d_t%d.pdf", subt, chan_lo, chan_hi,
			tbin);
	c1->Print( name);
}

void Asymmetry( UInt_t chan_lo = 300, UInt_t chan_hi = 300,
		TString subt = "subt")
{
	UInt_t i, j;
	Double_t eg, deg;
	Double_t par[6], red_chisq;
	Double_t theta[9], dtheta[9], as[9], das[9];
	TString name;

	TH1D* asym;
	TH2D* hA2;

	Get2DHist( chan_lo, chan_hi, subt);
	hA2 = (TH2D*) hA->Clone( "asymmetry2d");

	eg = WeightedEg( chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy/2) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy/2))/2;

	name = Form( "xs/asym/asymmetry_%d-%d.out", chan_lo, chan_hi);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( i = 0; i <= 8; i++)
	{
		j = i + 1;
		theta[i] = 10 + 20*i;
		dtheta[i] = 0;

		name = Form( "asym_tbin%d", j);
		asym = (TH1D*) hA2->ProjectionX( name, j, j);

		TF1 *f1 = new TF1( "f1", "[0]*cos(2*x*0.01745+[1])", -180, 180);
		asym->Fit( "f1", "REMQ0");
		par[0] = f1->GetParameter(0);
		par[1] = f1->GetParError(0);
		par[2] = f1->GetParameter(1);
		par[3] = f1->GetParError(1);
		red_chisq = f1->GetChisquare()/17;
		if ( par[2] < 0) {
			par[0] *= -1;
			par[2] *= -1;
		}

		as[i] = par[0];
		das[i] = par[1];

		name = Form( "%3d  %6.3f %4.3f  %4.2f %3.1f  %5.3f", (int) theta[i],
				as[i], das[i], par[2], par[3], red_chisq);
		cout << name << endl ;
		outFile << name << endl ;
	}
	outFile.close();

	c1 = new TCanvas( "c1", "Asymmetry", 100, 10, 700, 500);

	gr = new TGraphErrors( 9, theta, as, dtheta, das);
	name = Form( "Asymmetry for E_{#gamma} = %5.1f #pm %4.1f MeV", eg, deg);
	gr->SetTitle( name);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.2);
	gr->SetLineWidth(2);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta_{cm} (deg)");
	gr->GetYaxis()->SetTitle( "p_{#gamma}#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->SetRangeUser( 0, 180);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->Draw( "AP");

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

	name = Form( "plots/Asymmetry_%s_%d-%d.pdf", subt, chan_lo, chan_hi);
	c1->Print( name);
}

// Sets the pointer to a single 2D histogram of thetaCM vs. phiCM.
// for either "full", "empty", or "subt", and everything is
// weighted to Full Para scalers.
void Get2DHist( UInt_t chan_lo, UInt_t chan_hi, TString subt)
{
	UInt_t bin_lo, bin_hi;
	Double_t factor;
	TFile* file;

	TH2D* h2_perp_f;
	TH1D* hsc_perp_f;
	TH2D* h2_para_f;
	TH1D* hsc_para_f;
	TH2D* h2_perp_e;
	TH1D* hsc_perp_e;
	TH2D* h2_para_e;
	TH1D* hsc_para_e;

	// Tagger bins are offset by one...
	// Channels are 0-351 but bins are 1-352.
	bin_lo = chan_lo + 1;
	bin_hi = chan_hi + 1;

	// Full Perp
	file = fperpFile;
	GetHistFile( file);
	hS->GetXaxis()->SetRange( bin_lo, bin_hi);
	h2_perp_f = (TH2D*) hS->Project3D( "yz");
	hsc_perp_f = (TH1D*) hsc->Clone( "sc_perp_f");

	// Empty Perp
	file = eperpFile;
	GetHistFile( file);
	hS->GetXaxis()->SetRange( bin_lo, bin_hi);
	h2_perp_e = (TH2D*) hS->Project3D( "yz2");
	hsc_perp_e = (TH1D*) hsc->Clone( "sc_perp_e");

	// Corrected Perp
	if ( subt == "subt") {
		h2_perp_f->Sumw2();
		factor = ScaleFactor( hsc_perp_f, hsc_perp_e, bin_lo, bin_hi);
		h2_perp_f->Add( h2_perp_e, -factor);
	}

	// Full Para
	file = fparaFile;
	GetHistFile( file);
	hS->GetXaxis()->SetRange( bin_lo, bin_hi);
	h2_para_f = (TH2D*) hS->Project3D( "yz3");
	hsc_para_f = (TH1D*) hsc->Clone( "sc_para_f");

	// Empty Para
	file = eparaFile;
	GetHistFile( file);
	hS->GetXaxis()->SetRange( bin_lo, bin_hi);
	h2_para_e = (TH2D*) hS->Project3D( "yz4");
	hsc_para_e = (TH1D*) hsc->Clone( "sc_para_e");

	// Corrected Para
	if ( subt == "subt") {
		h2_para_f->Sumw2();
		factor = ScaleFactor( hsc_para_f, hsc_para_e, bin_lo, bin_hi);
		h2_para_f->Add( h2_para_e, -factor);
	}

	// Generate Difference, Sum, and Difference/Sum
	// Scaling perp to para
	if ( subt != "empty") {
		factor = ScaleFactor( hsc_para_f, hsc_perp_f, bin_lo, bin_hi);
		TH2D *diff2d = (TH2D*) h2_para_f->Clone( "diff2d");
		diff2d->Sumw2();
		diff2d->Add( h2_perp_f, -factor);
		TH2D *sum2d = (TH2D*) h2_para_f->Clone( "sum2d");
		sum2d->Sumw2();
		sum2d->Add( h2_perp_f, factor);
	}
	else {
		factor = ScaleFactor( hsc_para_e, hsc_perp_e, bin_lo, bin_hi);
		TH2D *diff2d = (TH2D*) h2_para_e->Clone( "diff2d");
		diff2d->Sumw2();
		diff2d->Add( h2_perp_e, -factor);
		TH2D *sum2d = (TH2D*) h2_para_e->Clone( "sum2d");
		sum2d->Sumw2();
		sum2d->Add( h2_perp_e, factor);
	}
	hA = (TH2D*) diff2d->Clone( "asym2d");
	hA->Sumw2();
	hA->Divide( sum2d);
}
 
// Sets pointers to the scaler and random-subtracted asymmetry 3D histograms
// for a certain input file (full/empty and perp/para).
void GetHistFile( TFile* file)
{
	Double_t pa;
	TString prompt, random, scalers;

	TH3D* hP;
	TH3D* hR;

	pa = 0.0833;

	prompt = "PhiCMCut2P_v_ThetaCMCut2P_v_TChanCut2P";
	random = "PhiCMCut2R_v_ThetaCMCut2R_v_TChanCut2R";
	scalers = "SumScalers152to503";

	hP = (TH3D*) file->Get( prompt);
	hR = (TH3D*) file->Get( random);
	hsc = (TH1D*) file->Get( scalers);

	hS = (TH3D*) hP->Clone( "sub");
	hS->Add( hR, -pa);
}

// Calculates bremsstrahlung-weighted photon energy for Full-Target-Para
// scalers.
Double_t WeightedEg( UInt_t chan_lo, UInt_t chan_hi)
{
	UInt_t i;
	Double_t eg, sum;

	TH1D *hscal = (TH1D*) fparaFile.Get( "SumScalers152to503");

	sum = 0;
	for ( i = chan_lo + 1; i <= chan_hi + 1; i++)
	{
		eg += hscal->GetBinContent( i)*tcd[i].energy;
		sum += hscal->GetBinContent( i);
	}
	eg /= sum;

	return( eg);
}

// This factor is for scaling the histogram belonging to the second set of
// scalers to that belonging to the first.
//
// At the moment, no deadtime correction!
Double_t ScaleFactor( TH1D* hsc1, TH1D* hsc2, UInt_t bin_lo, UInt_t bin_hi)
{
	Double_t c1, c2, factor;

	c1 = hsc1->Integral( bin_lo, bin_hi);
	c2 = hsc2->Integral( bin_lo, bin_hi);
	factor = c1/c2;

	return( factor);
}
