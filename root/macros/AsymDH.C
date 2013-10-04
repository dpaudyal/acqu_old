/*
 *		AsymDH.C
 *
 *		This is a rip-off of Alex's Alex_pol.C macro and calculates
 *		preliminary asymmetries.
 *
 */

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

// Full
TFile fperpFile( "histograms/Full/FullPerp.root");
TFile fparaFile( "histograms/Full/FullPara.root");

// Empty
TFile eperpFile( "histograms/Empty/EmptyPerp.root");
TFile eparaFile( "histograms/Empty/EmptyPara.root");

void Asym( Int_t bin = 0, TString subt = "full")
{
	Int_t tag_lo, tag_hi, deg, eg;
	Double_t sc_perp, sc_para, factor;
	Double_t par[6], red_chisq;
	TString hprompt, hrand, name;

	if ( bin == 0) {
		tag_lo = 229;
		tag_hi = 245;
		hprompt = "THR_PhiCM0P_v_ThetaCM0P";
		hrand = "THR_PhiCM0R_v_ThetaCM0R";
		eg = 320;
		deg = 20;
	}
	else if ( bin == 1) {
		tag_lo = 303;
		tag_hi = 306;
		hprompt = "THR_PhiCM1P_v_ThetaCM1P";
		hrand = "THR_PhiCM1R_v_ThetaCM1R";
		eg = 155;
		deg = 5;
	}
	else if ( bin == 2) {
		tag_lo = 299;
		tag_hi = 302;
		hprompt = "THR_PhiCM2P_v_ThetaCM2P";
		hrand = "THR_PhiCM2R_v_ThetaCM2R";
		eg = 165;
		deg = 5;
	}
	else if ( bin == 3) {
		tag_lo = 295;
		tag_hi = 298;
		hprompt = "THR_PhiCM3P_v_ThetaCM3P";
		hrand = "THR_PhiCM3R_v_ThetaCM3R";
		eg = 175;
		deg = 5;
	}

	if ( ( subt == "full") || ( subt == "empty")) {
		if ( subt == "full") {
			// Perp histograms
			TH2D *perp2D_P = (TH2D*)fperpFile.Get( hprompt);
			TH2D *perp2D_R = (TH2D*)fperpFile.Get( hrand);
			TH1D *perp_scal = (TH1D*)fperpFile.Get( "SumScalers152to503");

			// Para histograms
			TH2D *para2D_P = (TH2D*)fparaFile.Get( hprompt);
			TH2D *para2D_R = (TH2D*)fparaFile.Get( hrand);
			TH1D *para_scal = (TH1D*)fparaFile.Get( "SumScalers152to503");
		}
		else {
			// Perp histograms
			TH2D *perp2D_P = (TH2D*)eperpFile.Get( hprompt);
			TH2D *perp2D_R = (TH2D*)eperpFile.Get( hrand);
			TH1D *perp_scal = (TH1D*)eperpFile.Get( "SumScalers152to503");

			// Para histograms
			TH2D *para2D_P = (TH2D*)eparaFile.Get( hprompt);
			TH2D *para2D_R = (TH2D*)eparaFile.Get( hrand);
			TH1D *para_scal = (TH1D*)eparaFile.Get( "SumScalers152to503");
		}

		// Perp subtracted
		TH2D *perp2D_sub = (TH2D*)perp2D_P->Clone( "perpsub");
		perp2D_sub->Add( perp2D_R, -0.0833);

		// Para subtracted
		TH2D *para2D_sub = (TH2D*)para2D_P->Clone( "parasub");
		para2D_sub->Add( para2D_R, -0.0833);

		// Para scaling factor
		sc_perp = perp_scal->Integral( tag_lo, tag_hi);
		sc_para = para_scal->Integral( tag_lo, tag_hi);
		factor = sc_perp/sc_para;

		// Scale para
		TH2D *para2D_scal_sub = (TH2D*)para2D_sub->Clone();
		para2D_scal_sub->Scale( factor);
	}
	else if ( subt == "subt") {

		Double_t sc_perp_full, sc_para_full, sc_perp_empty, sc_para_empty;
		Double_t perpfactor, parafactor;

		// Full Perp histograms
		TH2D *fperp2D_P = (TH2D*)fperpFile.Get( hprompt);
		TH2D *fperp2D_R = (TH2D*)fperpFile.Get( hrand);
		TH1D *fperp_scal = (TH1D*)fperpFile.Get( "SumScalers152to503");

		// Full Perp subtracted
		TH2D *fperp2D_sub = (TH2D*)fperp2D_P->Clone( "fperpsub");
		fperp2D_sub->Add( fperp2D_R, -0.0833);

		// Full Para histograms
		TH2D *fpara2D_P = (TH2D*)fparaFile.Get( hprompt);
		TH2D *fpara2D_R = (TH2D*)fparaFile.Get( hrand);
		TH1D *fpara_scal = (TH1D*)fparaFile.Get( "SumScalers152to503");

		// Full Para subtracted
		TH2D *fpara2D_sub = (TH2D*)fpara2D_P->Clone( "fparasub");
		fpara2D_sub->Add( fpara2D_R, -0.0833);

		// Empty Perp histograms
		TH2D *eperp2D_P = (TH2D*)eperpFile.Get( hprompt);
		TH2D *eperp2D_R = (TH2D*)eperpFile.Get( hrand);
		TH1D *eperp_scal = (TH1D*)eperpFile.Get( "SumScalers152to503");

		// Empty Perp subtracted
		TH2D *eperp2D_sub = (TH2D*)eperp2D_P->Clone( "eperpsub");
		eperp2D_sub->Add( eperp2D_R, -0.0833);

		// Empty Para histograms
		TH2D *epara2D_P = (TH2D*)eparaFile.Get( hprompt);
		TH2D *epara2D_R = (TH2D*)eparaFile.Get( hrand);
		TH1D *epara_scal = (TH1D*)eparaFile.Get( "SumScalers152to503");

		// Empty Para subtracted
		TH2D *epara2D_sub = (TH2D*)epara2D_P->Clone( "eparasub");
		epara2D_sub->Add( epara2D_R, -0.0833);

		// Perp scaling factor
		sc_perp_full = fperp_scal->Integral( tag_lo, tag_hi);
		sc_perp_empty = eperp_scal->Integral( tag_lo, tag_hi);
		perpfactor = sc_perp_full/sc_perp_empty;
		// Perp Full-Empty subtracted
		TH2D *perp2D_sub = (TH2D*)fpara2D_sub->Clone( "perpsub");
		perp2D_sub->Add( eperp2D_sub, perpfactor);

		// Para scaling factor
		sc_para_full = fpara_scal->Integral( tag_lo, tag_hi);
		sc_para_empty = epara_scal->Integral( tag_lo, tag_hi);
		parafactor = sc_para_full/sc_para_empty;
		// Para Full-Empty subtracted
		TH2D *para2D_sub = (TH2D*)fpara2D_sub->Clone( "parasub");
		para2D_sub->Add( epara2D_sub, parafactor);

		// Perp-Para scaling factor; everthing is scaled to full
		sc_perp = fperp_scal->Integral( tag_lo, tag_hi);
		sc_para = fpara_scal->Integral( tag_lo, tag_hi);
		factor = sc_perp/sc_para;
		// Scale para
		TH2D *para2D_scal_sub = (TH2D*)para2D_sub->Clone();
		para2D_scal_sub->Scale( factor);
	}

	// Projections of Perp/Para onto phi axis.
	TH1D *perp = perp2D_sub->ProjectionY( "perp",1,9,"e");
	TH1D *para = para2D_scal_sub->ProjectionY( "para",1,9,"e");

	// Calculate sum and difference
	TH1D *diff = (TH1D*)para->Clone( "diff");
	diff->Add( perp, -1);
	TH1D *sum = (TH1D*)para->Clone( "sum");
	sum->Add( perp);

	TH1D *asym = (TH1D*)diff->Clone( "asym");
	asym->Divide( sum);

//	c1 = new TCanvas( "c1", "Data", 100, 10, 700, 700);
//	c1->Divide(1,2);
//	c1->cd(1);
//	perp->Draw();
//	c1->cd(2);
//	para->Draw();

	gStyle->SetOptStat(0);
	c2 = new TCanvas( "c2", "Data", 400, 10, 700, 500);
	c2->cd(1);
	asym->Draw();

	name = Form( "E_{#gamma} = %d #pm %d MeV  #theta = 0-180 deg", eg, deg);
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

	name = Form( "plots/Int_Asym_%d.pdf", eg);
	c2->Print( name);
}

void AsymBin( Int_t bin = 0, TString subt = "full")
{
	Int_t i, tag_lo, tag_hi, eg, deg;
	Double_t sc_perp, sc_para, factor;
	Double_t par[6], red_chisq;
	Double_t theta[9], dtheta[9], as[9], das[9];
	TString hprompt, hrand, name;

//	c1 = new TCanvas( "c1", "Data", 100, 10, 700, 1000);
//	c1->Divide(1,3);
//	c1 = new TCanvas( "c1", "Data", 100, 10, 700, 500);

	if ( bin == 0) {
		tag_lo = 229;
		tag_hi = 245;
		hprompt = "THR_PhiCM0P_v_ThetaCM0P";
		hrand = "THR_PhiCM0R_v_ThetaCM0R";
		eg = 320;
		deg = 20;
	}
	else if ( bin == 1) {
		tag_lo = 303;
		tag_hi = 306;
		hprompt = "THR_PhiCM1P_v_ThetaCM1P";
		hrand = "THR_PhiCM1R_v_ThetaCM1R";
		eg = 155;
		deg = 5;
	}
	else if ( bin == 2) {
		tag_lo = 299;
		tag_hi = 302;
		hprompt = "THR_PhiCM2P_v_ThetaCM2P";
		hrand = "THR_PhiCM2R_v_ThetaCM2R";
		eg = 165;
		deg = 5;
	}
	else if ( bin == 3) {
		tag_lo = 295;
		tag_hi = 298;
		hprompt = "THR_PhiCM3P_v_ThetaCM3P";
		hrand = "THR_PhiCM3R_v_ThetaCM3R";
		eg = 175;
		deg = 5;
	}

	if ( subt == "full") {
		// Perp histograms
		TH2D *perp2D_P = (TH2D*)fperpFile.Get( hprompt);
		TH2D *perp2D_R = (TH2D*)fperpFile.Get( hrand);
		TH1D *perp_scal = (TH1D*)fperpFile.Get( "SumScalers152to503");

		// Para histograms
		TH2D *para2D_P = (TH2D*)fparaFile.Get( hprompt);
		TH2D *para2D_R = (TH2D*)fparaFile.Get( hrand);
		TH1D *para_scal = (TH1D*)fparaFile.Get( "SumScalers152to503");
	}
	else if ( subt == "empty") {
		// Perp histograms
		TH2D *perp2D_P = (TH2D*)eperpFile.Get( hprompt);
		TH2D *perp2D_R = (TH2D*)eperpFile.Get( hrand);
		TH1D *perp_scal = (TH1D*)eperpFile.Get( "SumScalers152to503");

		// Para histograms
		TH2D *para2D_P = (TH2D*)eparaFile.Get( hprompt);
		TH2D *para2D_R = (TH2D*)eparaFile.Get( hrand);
		TH1D *para_scal = (TH1D*)eparaFile.Get( "SumScalers152to503");
	}
	else if ( subt == "subt") {
	}

	// Perp subtracted
	TH2D *perp2D_sub = (TH2D*)perp2D_P->Clone( "perpsub");
	perp2D_sub->Add( perp2D_R, -0.0833);

	// Para subtracted
	TH2D *para2D_sub = (TH2D*)para2D_P->Clone( "parasub");
	para2D_sub->Add( para2D_R, -0.0833);

	// Para scaling factor
	sc_perp = perp_scal->Integral( tag_lo, tag_hi);
	sc_para = para_scal->Integral( tag_lo, tag_hi);
	factor = sc_perp/sc_para;

	// Scale para
	TH2D *para2D_scal_sub = (TH2D*)para2D_sub->Clone();
	para2D_scal_sub->Scale( factor);

	name = Form( "xs/asym/asymmetry_%d.out", eg);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( i = 0; i <= 8; i++)
	{
		theta[i] = 10 + 20*i;
		dtheta[i] = 0;

//		c1->WaitPrimitive();

		// Projections of Perp/Para onto phi axis.
		TH1D *perp = perp2D_sub->ProjectionY( "perp", i+1, i+1, "e");
		TH1D *para = para2D_scal_sub->ProjectionY( "para", i+1, i+1, "e");

		// Calculate sum and difference
		TH1D *diff = (TH1D*)para->Clone( "diff");
		diff->Add( perp, -1);
		TH1D *sum = (TH1D*)para->Clone( "sum");
		sum->Add( perp);

		TH1D *asym = (TH1D*)diff->Clone( "asym");
		asym->Divide( sum);

//		c1->cd(1);
//		perp->Draw();
//		c1->cd(2);
//		para->Draw();
//		c1->cd(3);
//		asym->Draw();

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
		cout << par[0];
		cout << " +/- " << par[1];
		cout << "   " << par[2];
		cout << " +/- " << par[3];
		cout << "   " << red_chisq;
		cout << endl;

		as[i] = par[0];
		das[i] = par[1];

		name = Form( "%3d  %6.3f %4.3f  %4.2f %3.1f  %5.3f", (int) theta[i],
				as[i], das[i], par[2], par[3], red_chisq);
		outFile << name << endl ;
	}
	outFile.close();

	c1 = new TCanvas( "c1", "Asymmetry", 100, 10, 700, 500);

	gr = new TGraphErrors( 9, theta, as, dtheta, das);
	name = Form( "Asymmetry for E_{#gamma} = %3d #pm %2d MeV", eg, deg);
	gr->SetTitle( name);
//	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.2);
//	gr->SetLineColor(4);
	gr->SetLineWidth(2);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta_{cm} (deg)");
	gr->GetYaxis()->SetTitle( "p_{#gamma}#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->SetRangeUser(0,180);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->Draw( "AP");

	l1 = new TLine(0,0,180,0);
	l1->SetLineStyle(1);
	l1->SetLineWidth(1);
	l1->Draw();

	name = Form( "plots/Asymmetry_%d.pdf", eg);
	c1->Print( name);

}
