TGraphErrors* TaggEffPlot;
TH1F* histEelec;             // electron energy
TH1F* histHits;              // hits histogram
TH1D* histSc;                // scaler histogram
TH1D* histScB;               // background scaler histogram
TH1D* histS;                 // trigger scalers histogram
TH1D* histSB;                // background trigger scalers histogram

Double_t* Eelec;             // electron energy
Double_t* Hits;              // hits array
Double_t* Scalers;           // accumulated scalers
Double_t* ScalersB;          // accumulated scalers background run
Double_t* TaggEff;           // tagging efficiency
Double_t* dTaggEff;          // uncertainty
Double_t* TaggEffB;          // background corrected tagging efficiency
Double_t* dTaggEffB;         // uncertainty

Double_t meanEff, meanEffB, dmeanEff, dmeanEffB;

FILE* TaggEffLog;

void TaggEff(Char_t* ffore, Char_t* fback,
	     Int_t nchan=352, Int_t jclk=5, Int_t jinh=6, Int_t liveCorr=0)
{
  // calculate tagging efficiency
  Eelec = new Double_t[nchan];
  Hits = new Double_t[nchan];
  Scalers = new Double_t[nchan];
  ScalersB = new Double_t[nchan];
  TaggEff = new Double_t[nchan];
  dTaggEff = new Double_t[nchan];
  TaggEffB = new Double_t[nchan];
  dTaggEffB = new Double_t[nchan];

  TFile* fore = new TFile(ffore);
  if( !fore ){ printf("File %s not found\n", ffore) return; }
  TFile* back = new TFile(fback);
  if( !fore ){ printf("File %s not found\n", fback) return; }
  Char_t logfile[128];
  strcpy(logfile,"./");
  strcat(logfile,ffore);
  Char_t* suff = strstr(logfile,".root");
  strcpy(suff,".log");
  printf("\n Text output saved to file: %s\n\n",logfile);
  if( !(TaggEffLog = fopen(logfile,"w+")) ){
    printf("Unable to open TaggEff.log\n");
    return;
  }

  histEelec =  (TH1F*)fore->Get("FPD_Eelect");
  histHits = (TH1F*)fore->Get("FPD_Hits");
  histSc = (TH1D*)fore->Get("SumScalers12to363");
  histScB = (TH1D*)back->Get("SumScalers12to363");
  //  histSc = (TH1D*)fore->Get("FPD_ScalerAcc");
  //  histScB = (TH1D*)back->Get("FPD_ScalerAcc");
  histS = (TH1D*)fore->Get("SumScalers0to11");
  histSB = (TH1D*)back->Get("SumScalers0to11");

  if( !histEelec ) printf("Electron energy hist not found\n");
  if( !histHits ) printf("Hits hist not found\n");
  if( !histSc ) printf("Scalers hist not found\n");
  if( !histS ) printf("Trig-Scalers hist not found\n");
  if( !histScB ) printf("Background Scalers hist not found\n");
  if( !histSB ) printf("Background Trig-Scalers hist not found\n");


  Double_t clk = (Double_t)histS->GetBinContent(jclk);
  Double_t inh = (Double_t)histS->GetBinContent(jinh);
  Double_t live = clk/inh;
  Double_t clkB = (Double_t)histSB->GetBinContent(jclk);
  Double_t inhB = (Double_t)histSB->GetBinContent(jinh);
  Double_t liveB = clkB/inhB;
  printf("%f %f %f %f %f %f\n",clk,inh,live,clkB,inhB,liveB);

  Double_t bscale = inh/inhB;
  Double_t liveC;
  if( liveCorr ) liveC = live;
  else liveC = 1.0;
  Double_t ScB, dScB;

  fprintf(TaggEffLog,
	 " Tagging efficiency calculation\n J.R.M.Annand 28/01/08\n\n");

  fprintf(TaggEffLog,"Clock = %10g   Inhib. Clock = %10g   Live = %8g\n\
Backgr.clock = %10g    Backgr.Inhib.clock = %10g  Live = %8g\n\
Background scale factor = %8g\n\n",
	  clk,inh,1/live,clkB,inhB,1/liveB,bscale);

  fprintf(TaggEffLog,
"   Gi  Hits   Scaler  Background    TaggEff           TaggEff(back-sub)\n");

  for(UInt_t i = 0; i< nchan; i++){
    //    Eelec[i] = (Double_t)histEelec->GetBinContent(i);
    Eelec[i] = (Double_t)i;
    Hits[i]  = (Double_t)histHits->GetBinContent(i);
    Scalers[i] = (Double_t)histSc->GetBinContent(i);
    ScalersB[i] = (Double_t)histScB->GetBinContent(i);

    if( (Scalers[i] == 0) || (Hits[i] == 0) ){
      TaggEff[i] = TaggEffB[i] = dTaggEff[i] = dTaggEffB[i] = 0.0;
    }
    else{
      ScB = Scalers[i] - bscale*ScalersB[i];
      dScB = TMath::Sqrt(ScalersB[i]) * bscale;
      dScB = TMath::Sqrt( Scalers[i] + dScB*dScB )/ScB;
      
      TaggEff[i] = liveC * Hits[i]/Scalers[i];
      TaggEffB[i] = liveC * Hits[i]/( Scalers[i] - bscale*ScalersB[i] );
      dTaggEff[i] = TaggEff[i] * TMath::Sqrt( 1./Hits[i] + 1./Scalers[i] );
      dTaggEffB[i] = TaggEffB[i] * TMath::Sqrt( 1./Hits[i] + dScB*dScB );
    }
    fprintf(TaggEffLog,
	    "%5g %5g %8g %8g   %7.3g +-%6.2g  %7.3g +-%6.2g\n",
	   Eelec[i], Hits[i], Scalers[i], ScalersB[i],
	   TaggEff[i], dTaggEff[i], TaggEffB[i], dTaggEffB[i]);
  }
  TaggEffDraw(nchan,ffore);
  fprintf(TaggEffLog,
   "\nMean tagging efficiency %7.3g +- %6.2g\n\
                        %7.3g +- %6.2g without background subtraction\n",
	  meanEffB, dmeanEffB, meanEff, dmeanEff );
  fclose(TaggEffLog);
  fore->Close();
  back->Close();
}

void TaggEffDraw(Int_t nchan, Char_t* run = ""){
  // refills the TeffData into gr_data avoiding "bad" channels and
  // plots the graph.
  TCanvas* canv;
  if( !( canv = (TCanvas*)(gROOT->FindObject("Teff")) ) ){
    TCanvas* canv = new TCanvas("Teff","Tagging Efficiency",240,180,700,900);
    canv->SetFillStyle(4000);
    canv->Divide(1,2,0.01,0.01);
    canv->cd(1);
  }
  else canv->cd(1);
  char title[255];
  UInt_t j = 0;
  for(UInt_t i=0; i<nchan; i++){ 
    if( TaggEff[i] > 0 ) j++;
  }

  Double_t* ee = new Double_t[j];
  Double_t* tf = new Double_t[j];
  Double_t* dtf = new Double_t[j];
  j = 0;
  meanEff = dmeanEff = 0.0;
  for(UInt_t i=0; i<nchan; i++){
    if( TaggEff[i] > 0 ){
      ee[j] = Eelec[i];
      tf[j] = TaggEff[i];
      meanEff = meanEff + tf[j];
      dtf[j] = dTaggEff[i];
      dmeanEff = dmeanEff + dtf[j]*dtf[j];
      j++;
    }
  }
  meanEff = meanEff/j;
  dmeanEff = TMath::Sqrt(dmeanEff)/j;

  TaggEffPlot = new TGraphErrors(j, ee, tf, NULL, dtf);
  strcpy(title,"Tagging-Efficiency");
  if (strlen(run)){ 
    strcat(title, " of ");
    strcat(title, run);
  }
  strcat(title, " (Livetime Correction: 0)");
  TaggEffPlot->SetTitle(title);
  TaggEffPlot->SetMaximum(0.8);
  TaggEffPlot->SetMinimum(0.0);
  TaggEffPlot->GetXaxis()->SetTitle("Tagger Channel (Goosey)");
  TaggEffPlot->GetYaxis()->SetTitle("#epsilon_{tagg}");
  TaggEffPlot->GetXaxis()->CenterTitle();
  TaggEffPlot->GetYaxis()->CenterTitle();
  TaggEffPlot->SetMarkerColor(kBlue);
  TaggEffPlot->SetMarkerStyle(21);
  TaggEffPlot->SetMarkerSize(0.7);
  TaggEffPlot->Draw("AP");

  //
  canv->cd(2);
  j = 0;
  meanEffB = 0.0;
  for(UInt_t i=0; i<nchan; i++){
    if( TaggEff[i] > 0 ){
      ee[j] = Eelec[i];
      tf[j] = TaggEffB[i];
      meanEffB = meanEffB + tf[j];
      dtf[j] = dTaggEffB[i];
      dmeanEffB = dmeanEffB + dtf[j]*dtf[j];
      j++;
    }
  }
  meanEffB = meanEffB/j;
  dmeanEffB = TMath::Sqrt(dmeanEffB)/j;

  TaggEffPlotB = new TGraphErrors(j, ee, tf, NULL, dtf);
  strcpy(title,"Background-Corrected Tagging-Efficiency");
  if (strlen(run)){ 
    strcat(title, " of ");
    strcat(title, run);
  }
  strcat(title, " (Livetime Correction: 0)");
  TaggEffPlotB->SetTitle(title);
  TaggEffPlotB->SetMaximum(0.8);
  TaggEffPlotB->SetMinimum(0.0);
  TaggEffPlotB->GetXaxis()->SetTitle("Tagger Channel (Goosey)");
  TaggEffPlotB->GetYaxis()->SetTitle("#epsilon_{tagg}");
  TaggEffPlotB->GetXaxis()->CenterTitle();
  TaggEffPlotB->GetYaxis()->CenterTitle();
  TaggEffPlotB->SetMarkerColor(kRed);
  TaggEffPlotB->SetMarkerStyle(20);
  TaggEffPlotB->SetMarkerSize(0.7);
  TaggEffPlotB->Draw("AP");
}
