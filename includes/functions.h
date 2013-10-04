/* 
 * Functions.h
 *
 * Some useful kinematics functions for coherent two-body
 * photoproduction of a neutral pion (or any other particle).
 *
 * 22.10.2006
 *
 * DLH 
 *
 */

// Calculates minimum gamma-gamma opening angle.
Double_t Tgg_Min( Double_t T, Double_t mpi)
{
	Double_t term, tgg_min;

	term = mpi/sqrt(T*(T+2*mpi));
	tgg_min = 2*atan(term);

	return( tgg_min);
}

// Calculates maximum KE of the pion.
Double_t qT_max( Double_t eg, Double_t pm, Double_t qm)
{
	Double_t W, W_cm, S, beta, gamma;
	Double_t qcm_e, qcm_p, qcm_x, qcm_z;
	Double_t qx, qz, qp, qe, qT_max, th_cm;

	W = eg + pm;
	S = Sqr( pm) + 2*eg*pm;
	W_cm = sqrt( S);

	beta = eg/W;
	gamma = W/W_cm;

	qcm_e = (S + Sqr(qm) - Sqr(pm))/2/W_cm;
	qcm_p = Momentum( qcm_e, qm);

	th_cm = 0;  // Condition for max pion KE

	qcm_x = qcm_p*sin( th_cm);
	qcm_z = qcm_p*cos( th_cm);
	qx = qcm_x;
	qz = gamma*(qcm_z + beta*qcm_e);
	qp = sqrt( Sqr( qx) + Sqr(qz));
	qe = Energy( qp, qm);
	qT_max = qe - qm;

	return( qT_max);
}

// Calculates relativistic energy from mass and momentum.
Double_t Energy( Double_t mom, Double_t m)
{
	return( sqrt(mom*mom + m*m));
}

// Calculates relativistic momentum from mass and energy.
Double_t Momentum( Double_t en, Double_t m)
{
	if ( en >= m) return( sqrt( en*en - m*m));
	else return( -1);
}

// Calculates the square.
Double_t Sqr( Double_t x)
{
	return( x*x);
}

// Calculates the threshold energy.
Double_t E_thr( Double_t pm, Double_t qm)
{
	return( (Sqr( qm) + 2*qm*pm)/2/pm);
}

// Does a simple linear interpolation.
Double_t Linear( Double_t x0, Double_t y0, Double_t x1, Double_t y1, Double_t x)
{
	return( y0 + (y1 - y0)/(x1 - x0)*(x - x0));
}

// Calculates the lab momentum of the pion.
// This one takes the pion LAB theta!
Double_t qp( Double_t ke, Double_t pm, Double_t qth, Double_t qm)
{
	Double_t qmom;
	Double_t W;
	Double_t A, B, C, a, b, c;

	W = ke + pm;

	A = 2*ke*pm + Sqr( qm);
	B = 2*ke*cos( qth);
	C = 2*W;

	a = Sqr( C) - Sqr( B);
	b = -2*A*B;
	c = Sqr( C*qm) - Sqr( A);

	qmom = (-b + sqrt( Sqr( b) - 4*a*c))/2/a;

	if ( qmom != qmom) return( -1);
	else return( qmom);
}

// Calculates the lab momentum of the pion.
// This one takes the pion CM theta, but still calculates the 
// lab momentum!
Double_t qp_thcm( Double_t ke, Double_t pm, Double_t qth_cm, Double_t qm)
{

	TLorentzVector k, p, q, pIn;
	TLorentzVector k_cm, p_cm, q_cm;
	TVector3 cmBoost, labBoost, q3mom;
	
	Double_t qm, qph_cm, mom_cm;
	Double_t qx_cm, qy_cm, qz_cm, qe_cm;
	Double_t pm, S;

	k.SetPxPyPzE( 0, 0, ke, ke);
	p.SetPxPyPzE( 0, 0, 0, pm);

	pIn = k + p;
	labBoost = pIn.BoostVector();
	cmBoost = -pIn.BoostVector();

	k_cm = k;
	k_cm.Boost( cmBoost);
	p_cm = p;
	p_cm.Boost( cmBoost);

	qth_cm *= kD2R;
	qph_cm = 0;

	S = 2*ke*pm + pm*pm;
	qe_cm = (S - pm*pm + qm*qm)/2/sqrt(S);

	mom_cm = Momentum( qe_cm, qm);
	qx_cm = mom_cm*sin( qth_cm)*cos(qph_cm);
	qy_cm = mom_cm*sin( qth_cm)*sin(qph_cm);
	qz_cm = mom_cm*cos( qth_cm);

	q_cm.SetPxPyPzE( qx_cm, qy_cm, qz_cm, qe_cm);

	q = q_cm;
	q.Boost( labBoost);

	q3mom = q.Vect();

	return( q3mom.Mag());
}

// Calculates the recoil momentum of the target.  It uses the above
// function qp.
Double_t ppr( Double_t ke, Double_t pm, Double_t qth, Double_t qm)
{
	Double_t ppx, ppz, qx, qz, q, pp;

	q = qp( ke, pm, qth, qm);
	qz = q*cos( qth);
	qx = q*sin( qth);

	ppz = ke - qz;
	ppx = -qx;

	pp = sqrt( Sqr( ppx) + Sqr( ppz));

	return( pp);
}

// Returns the file number for the current run.  Note that it gives you the
// NEXT run if you use it with EndFile!!
Int_t GetRunDave()
{
	Int_t i1, i2, run;
	TString name;

	if ( gAR->IsOnline())
	{
		name = gAR->GetFileName();

		i1 = name.Index("_");
		name.Remove(0,i1+1);
		i2 = name.Index(".");
		name.Remove( i2);
		run = name.Atoi();
	}
	else {

		Int_t n = 0;
		while ( gAR->GetTreeFileList( n) != NULL)
			name = gAR->GetTreeFileList( n++);

		i1 = name.Index("_");
		name.Remove(0,i1+1);
		i2 = name.Index(".");
		name.Remove( i2);
		run = name.Atoi();
	}

	return( run);
}
//
// WeightedPar
//
// Calculates bremsstrahlung-weighted parameter for scalers in file.
// Parameter must exist in the tcd structure...
//
// This is a function of tagger channel, and NOT bin number!
//
Double_t WeightedPar( TH1D* hscal, TString par, UInt_t chan_lo, UInt_t chan_hi,
		UInt_t j = 0)
{
	UInt_t i;
	Double_t quant, numer, denom, factor;
	TString scalers;

	numer = 0;
	for ( i = chan_lo; i <= chan_hi; i++)
	{
		if ( par == "energy") quant = tcd[i].energy;
		else if ( par == "tageff") quant = tcd[i].etag;
		else if ( par == "deteff") quant = tcd[i].edet[j];

		numer += hscal->GetBinContent( i+1)*quant;
	}
	denom = hscal->Integral( chan_lo+1, chan_hi+1);
	factor = numer/denom;

	return( factor);
}

//
// HistSF
//
// This factor is for scaling the histogram belonging to the second set of
// scalers to that belonging to the first.
//
// No deadtime correction!
//
Double_t HistSF( TH1D* hsc1, TH1D* hsc2, UInt_t chan_lo, UInt_t chan_hi)
{
	Double_t c1, c2, factor;

	c1 = hsc1->Integral( chan_lo+1, chan_hi+1);
	c2 = hsc2->Integral( chan_lo+1, chan_hi+1);
	factor = c1/c2;

	return( factor);
}

//
// DeadTimeSF
//
// This scale factor is for correcting the tagger scalers.  You must DIVIDE
// the tagger scalers by this amount, or you can multiply the yield by this
// factor; it is equivalent!
//
// See Section x.x on page xx of the Tagger Manual.
//
Double_t DeadTimeSF( TH1D* hdtsc)
{
	Double_t f_tagg, f_tot, f_dead;
	TString dtscalers;

	f_tagg = hdtsc->GetBinContent( 32)/hdtsc->GetBinContent( 31);
	f_tot = hdtsc->GetBinContent( 25)/hdtsc->GetBinContent( 26);
	f_dead = f_tagg/f_tot;

	return( f_dead);
}

//
// WeightedPol
//
// This calculates the weighted polarization for a number of tagger channels.
//
// NOTE: Here we are working in tagger channels and NOT bins!!!
//
// The variable 'jpol' has the following meaning:
// 	'0' is full para
// 	'1' is full perp
// 	'2' is empty para
// 	'3' is empty perp
//
Double_t WeightedPol( TH1D* hscal, UInt_t chan_lo, UInt_t chan_hi, UInt_t jpol)
{
	UInt_t i;
	Double_t numer, denom, factor;
	TString name, scalers;

	numer = 0;
	for ( i = chan_lo; i <= chan_hi; i++) {
		if ( tcd[i].pol[jpol] != 0)
			numer += hscal->GetBinContent( i+1)*tcd[i].pol[jpol];
		else {
			name = Form( "Error: polarization for tagger channel %d is zero!", i);
			cout << name << endl;
			exit( -1);
		}
	}
	denom = hscal->Integral( chan_lo+1, chan_hi+1);
	factor = numer/denom;

	return( factor);
}

//
// TaggEffCorr
//
// This calculates the weighted tagging efficiency correction factor for a
// number of tagger channels.
//
// NOTE: Here we are working in tagger channels and NOT bins!!!
//
// The variable 'jpol' has the following meaning:
// 	'0' is full para
// 	'1' is full perp
// 	'2' is empty para
// 	'3' is empty perp
//
Double_t TaggEffCorr( TH1D* hscal, UInt_t chan_lo, UInt_t chan_hi, UInt_t jpol)
{
	UInt_t i;
	Double_t numer, denom, factor;
	TString name, scalers;

	numer = 0;
	for ( i = chan_lo; i <= chan_hi; i++) {
		if ( tcd[i].tecorr[jpol] != 0)
			numer += hscal->GetBinContent( i+1)*tcd[i].tecorr[jpol];
		else {
			name = Form( "Error: tag eff for tagger channel %d is zero!", i);
			cout << name << endl;
			exit( -1);
		}
	}
	denom = hscal->Integral( chan_lo+1, chan_hi+1);
	factor = numer/denom;

	return( factor);
}
