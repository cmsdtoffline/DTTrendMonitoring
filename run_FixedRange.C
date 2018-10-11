#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "EfficiencyMonitor.h"
#include "Utilities.h"
#include "Types.h"

void run_FixedRange(string refName ="", string storingName ="", string storedName = "", string run0 = "", string run1 = "", string run2 = "", string run3 ="", string run4 ="", string run5 ="") {
  

  gROOT->LoadMacro("DistTrend.cxx++");
  gROOT->LoadMacro("EffTrend.cxx++");
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

  else if(run0 != "" && run1 != "" && run2 !="" && run3 != "" && run4 ==""){  
    chain->Add(run0.c_str());
    chain->Add(run1.c_str());
    chain->Add(run2.c_str());
    chain->Add(run3.c_str());
    chain->Add(run4.c_str());
  }

  else if(run0 != "" && run1 != "" && run2 !="" && run3 != "" && run4 != "" && run5 ==""){  
    chain->Add(run0.c_str());
    chain->Add(run1.c_str());
    chain->Add(run2.c_str());
    chain->Add(run3.c_str());
    chain->Add(run4.c_str());
    chain->Add(run5.c_str());
  }
  
  else {cout<<"Problems with files "<<endl;    exit(1);}

  context FixedCont;
  
  FixedCont.name = "Fixed";

  Var InsLumi;
  Var Pileup;
  Var BkGr;
  Var Empty;

  variables::initVar(Pileup,"Pileup","variablesSetting.json");
  variables::initVar(InsLumi,"InsLumi","variablesSetting.json");
  variables::initVar(BkGr,"Bkg","variablesSetting.json",1,0);
  variables::initVar(Empty,"Empty","variablesSetting.json",0,0);

  FixedCont.var   = { {InsLumi.name,InsLumi},{Pileup.name,Pileup},{BkGr.name,BkGr},{Empty.name,Empty}};
  FixedCont.nVar  = FixedCont.var.size();
  FixedCont.webFolder = "~/www/DT";
  
  EfficiencyMonitor *eff = new EfficiencyMonitor(FixedCont,chain,refName.c_str(),storingName,storedName);

  if(stat(("data/DeadList/DeadList_"+refName).c_str(),&st) != 0) { 
   cout<<"Dead cell list named DeadList_"+refName+" doesn't exist in data/DeadList/\nStarting pre-loop to produce it"<<endl;
   eff->PreLoop();
  }
  else {
    cout<<"Dead cell list named DeadList_"+refName+" found in data/DeadList/\nIt will be used instead of producing it"<<endl;
  }

  eff->Loop();
  cout<<"Write"<<endl;
  eff->write();
  cout<<"Plot"<<endl;
  eff->plot();
  eff->close();  

}
