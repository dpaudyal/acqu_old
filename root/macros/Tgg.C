
#include "includes/physics.h"
#include "includes/functions.h"

void PlotTggCut( Double_t Eg = 150)
{
	Double_t min, offset;
	TString tgt, name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 200, 20, 800, 400);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);

	c1->cd( 1);
	TF1 *f1 = new TF1( "f1", TggOut, 0, 180, 3);
	offset = 0;
	f1->SetParameters( Eg, kM_C12_MEV, offset);
	TF1 *f2 = new TF1( "f2", TggOut, 0, 180, 3);
	offset = 3;
	f2->SetParameters( Eg, kMP_MEV, offset);
	f2->SetLineStyle( 2);

	name = Form( "histograms/MonteCarlo/pi0/pi0_p_%d.root", (int) Eg);
	simp = new TFile( name);
	name = "THR_TGGP_v_ThetaCMP";
	TH2D *h2sim = (TH2D*)simp->Get( name);

	min = f2->Eval( 0) - 10;
	h2sim->GetYaxis()->SetRangeUser( min, 180);
	h2sim->Draw( "zcol");

	f1->Draw( "same");
	f2->Draw( "same");

	// Cut Efficiency - Proton
	Int_t cts, cts_cut1;
	Double_t eff;
	name = "THR_ThetaCMP";
	TH1D *h1 = (TH1D*)simp->Get( name);
	cts = h1->Integral();
	name = "THR_ThetaCMCut1P";
	TH1D *h1_cut1 = (TH1D*)simp->Get( name);
	cts_cut1 = h1_cut1->Integral();

	eff = (double) cts_cut1/cts;
	cout << cts_cut1;
	cout << " " << cts;
	cout << " " << eff;
	cout << endl;

	c1->cd( 2);

	name = Form( "histograms/MonteCarlo/pi0/pi0_c_%d.root", (int) Eg);
	simc = new TFile( name);
	name = "THR_TGGP_v_ThetaCMP";
	TH2D *h2sim = (TH2D*)simc->Get( name);

	min = f2->Eval( 0) - 10;
	h2sim->GetYaxis()->SetRangeUser( min, 180);
	h2sim->Draw( "zcol");

	f1->Draw( "same");
	f2->Draw( "same");

	// Cut Efficiency - C12
	name = "THR_ThetaCMP";
	TH1D *h1c = (TH1D*)simc->Get( name);
	cts = h1c->Integral();
	name = "THR_ThetaCMCut1P";
	TH1D *h1c_cut1 = (TH1D*)simc->Get( name);
	cts_cut1 = h1c_cut1->Integral();

	eff = (double) cts_cut1/cts;
	cout << cts_cut1;
	cout << " " << cts;
	cout << " " << eff;
	cout << endl;

}

Double_t TggOut( Double_t *x, Double_t *par)
{
	Double_t thcm, Eg, mtgt, Tgg, offset;
	Double_t q_pi, T_pi;

	thcm = x[0];

	Eg = par[0];
	mtgt = par[1];
	offset = par[2];

	q_pi = qp_thcm( Eg, mtgt, thcm, kMPI0_MEV);
	T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
	Tgg = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

	return( Tgg - offset);
}
