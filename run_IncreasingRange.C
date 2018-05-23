//#include <stdio>
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "EfficiencyMonitor.h"
#include "DataSetting.h"

int run_IncreasingRange(string refName ="", string storingName ="", string storedName = "", string run0 = "", string run1 = "", string run2 = "", string run3 ="", string run4 ="", string run5 ="", string run6 ="") {

  gROOT->LoadMacro("EfficiencyMonitor.cc++");

  TChain * chain = new TChain("DTTree");

  if(run0 != "" && run1 == ""){
    chain->Add(run0.c_str());
  }

  else if(run0 != "" && run1 != "" && run2 ==""){
    chain->Add(run0.c_str());
    chain->Add(run1.c_str());
  }

  else if(run0 != "" && run1 != "" && run2 ==""){  
    chain->Add(run0.c_str());
    chain->Add(run1.c_str());
    chain->Add(run2.c_str());
  }

  else if(run0 != "" && run1 != "" && run2 !="" && run3 ==""){  
    chain->Add(run0.c_str());
    chain->Add(run1.c_str());
    chain->Add(run2.c_str());
    chain->Add(run3.c_str());
  }
  
 
  else {cout<<"Problems with files "<<endl;    exit(1);}
 
  context runCont;
 
  runCont.name = "run";
  runCont.nVar=2;

  runCont.varTitle  = varTitle_run;
  runCont.varName   = varName_run;
  runCont.varLabel  = varLabel_run;
  runCont.webFolder = WebFolder;

  runCont.slices.push_back(runslice);
  runCont.slices.push_back(bunchXingslice);
  runCont.slices.push_back(bkgslice);

  EfficiencyMonitor *eff = new EfficiencyMonitor(runCont,chain,refName.c_str(),storingName,storedName);

  if(stat(("data/DeadList/DeadList_"+refName).c_str(),&st) != 0) { 
   cout<<"dead cell list named DeadList_"+refName+" doesn't exist in data/DeadList/\nStarting pre-loop to produce it"<<endl;
   eff->PreLoop();
  }
  else {
    cout<<"dead cell list named DeadList_"+refName+" found in data/DeadList/\nIt will be used instead of producing it"<<endl;
  }

  cout<<"Loop"<<endl;
  eff->Loop();
  eff->write();
  eff->plot();

  return 1;
}
