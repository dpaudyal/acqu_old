// Note that you must define the following structure with the proper
// variables in it before you can use this:
//
//	typedef struct {
//
//		Int_t egamma;
//		Double_t energy;
//		Double_t denergy;
//		Double_t etag;
//		Double_t d_etag;
//		Double_t edet[10];
//		Double_t d_edet[10];
//		Double_t pa_ratio;
//		Double_t xs;
//		Double_t dxs;
//	
//	} TData;
//	TData tcd[352];
//

//
// Read in channel tagger energies.
//
void ReadTagEng( TString file)
{
	UInt_t i;
	Double_t eff, deff;

	cout << "Reading in parameters from " << file << endl;

	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << file;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> eff >> deff;
		tcd[i].energy = eff;
		tcd[i].denergy = deff/2;
		tcd[i].egamma = (int)(eff + 0.5);
	}
	inFile.close();
}

//
// Read in channel random subtraction factors.
//
void ReadSubt( TString file)
{
	UInt_t i;
	Double_t eff, deff;

	cout << "Reading in parameters from " << file << endl;

	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << file;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> eff;
		tcd[i].pa_ratio = eff;
	}
	inFile.close();
}

//
// Read in channel detection efficiencies.
//
void ReadDetEff()
{
	UInt_t i, j;
	Double_t eff, deff, junk;
	TString file;

	cout << "Reading in detection efficiency parameters" << endl;

	for ( i = 0; i <= 351; i++) {
		for ( j = 0; j < 10; j++) {
			tcd[i].edet[j] = 0.8;
			tcd[i].d_edet[j] = 0.01;
		}
		file = Form( "Eff_pi0_p_chan_%d.out", i);
		file.Prepend( "xs/eff/det_eff/");
		ifstream inFile( file);
		if ( !inFile.is_open()) {
			cout << "Efficiency for channel ";
			cout << i;
			cout << " does not exist.  Using 80%.";
			cout << endl;
		}
		else {
			for ( j = 0; j < 10; j++) {
				inFile >> junk >> eff >> deff;
				tcd[i].edet[j] = eff;
				tcd[i].d_edet[j] = deff;
			}
		}
		inFile.close();
	}
}

//
// Read in channel tagging efficiencies.
//
void ReadTagEff( TString file)
{
	UInt_t i;
	Double_t eff, deff;

	cout << "Reading in parameters from " << file << endl;

	ifstream inFile( file);
	while( !inFile.eof()) {
		inFile >> i >> eff >> deff;
		if ( eff != 0) {
			tcd[i].etag = eff;
			tcd[i].d_etag = deff;
		}
		else {
			tcd[i].etag = -1;
			tcd[i].d_etag = -1;
		}
	}
	inFile.close();
}
