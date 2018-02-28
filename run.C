//#include <stdio>
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "EfficiencyMonitor.h"


int run(string name ="", string dateName ="", string oldDate = "", string run0 = "", string run1 = "", string run2 = "", string run3 ="", string run4 ="", string run5 ="", string run6 ="") {

  gROOT->LoadMacro("EfficiencyMonitorLumiPU.C++");


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

	EfficiencyMonitor *eff = new EfficiencyMonitor(chain,name.c_str(),dateName,oldDate);
	//eff->PreLoop();
	eff->Loop();
	//eff->PostLoop();

	return 1;
}
