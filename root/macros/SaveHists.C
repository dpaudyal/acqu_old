void SaveHists(char *dir){
  Int_t RunNumber=0;
  Char_t runNumString[20]; 						//temp string for getting run number
  Char_t CurrentRunFileName[200];
  Char_t RootFileName[200];

  if(gAR->IsOnline()){							//if online
    strcpy(CurrentRunFileName,gAR->GetFileName()); 		    	//save curent file name
  }
  else{									//if offline
    int n=0;								//save last name on list
    while(gAR->GetTreeFileList(n)!=NULL){
      strcpy(CurrentRunFileName,gAR->GetTreeFileList(n++));
    }
  }
  strcpy(runNumString,strrchr(CurrentRunFileName,'_')+1);		//copy file name after last "_"
  strcpy(runNumString,strrchr(CurrentRunFileName,'_')+1);		//copy file name after last "_"
  strcpy(strrchr(runNumString,'.'),"     ");				//replace .suffix with "      "
  sscanf(runNumString,"%d",&RunNumber);					//scan the run number
  sprintf(RootFileName,"%s_%d.root",dir,RunNumber);			//construct rootfile name
  TFile f(RootFileName,"recreate");
  gROOT->GetList()->Write();
  f.Close();
  fprintf(stdout,"All histograms saved to %s\n\n",RootFileName);
  if(gAN!=NULL) fflush(gAN->GetLogStream());
  if(gAR!=NULL) fflush(gAR->GetLogStream());
  if(gDS!=NULL) fflush(gDS->GetLogStream());
}
