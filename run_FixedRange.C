#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "EfficiencyMonitor.h"
#include "Utilities.h"
#include "Types.h"

int run_FixedRange(string refName ="", string storingName ="", string storedName = "", string run0 = "", string run1 = "", string run2 = "", string run3 ="", string run4 ="", string run5 ="", string run6 ="") {
  
  gROOT->LoadMacro("Hist.cxx++");
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

  context FixedCont;
  
  FixedCont.name = "Fixed";

  Var InsLumi;
  Var Pileup;
  Var BckGr;

  variables::initVar(InsLumi,"InsLumi","variablesSetting.json");
  variables::initVar(Pileup,"Pileup","variablesSetting.json");
  variables::initVar(BckGr,"Bck","variablesSetting.json");

  FixedCont.var   = {InsLumi,Pileup,BckGr};
  FixedCont.nVar = FixedCont.var.size();
  FixedCont.webFolder = "~/www/DT";
  
  EfficiencyMonitor *eff = new EfficiencyMonitor(FixedCont,chain,refName.c_str(),storingName,storedName);

  if(stat(("data/DeadList/DeadList_"+refName).c_str(),&st) != 0) { 
   cout<<"dead cell list named DeadList_"+refName+" doesn't exist in data/DeadList/\nStarting pre-loop to produce it"<<endl;
   eff->PreLoop();
  }
  else {
    cout<<"dead cell list named DeadList_"+refName+" found in data/DeadList/\nIt will be used instead of producing it"<<endl;
  }
  eff->Loop();
  eff->write();
  eff->plot();
  
  return 1;
}
