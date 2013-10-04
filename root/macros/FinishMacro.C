void FinishMacro(Char_t* file = NULL)
{
  // Stuff to do at the end of an analysis run
  // Here all spectra are saved to disk
  printf("End-of-Run macro executing\n");
  if( !file ) file = "ARHistograms.root";
  TFile f(file,"recreate");
  gROOT->GetList()->Write();
  f.Close();
  printf("All histograms saved to %s\n\n",file);
}
