void SaveHistos()
{
  // Stuff to do at the end of an analysis run
  // Here all spectra are saved to disk

  printf("End-of-Run macro executing\n");

  TFile f("ARHistograms.root","recreate");

  gROOT->GetList()->Write();
  f.Close();

  printf("All histograms saved to ARHistograms.root\n\n");
}
