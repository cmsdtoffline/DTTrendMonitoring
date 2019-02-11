#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "EfficiencyMonitor.h"
#include "Utilities.h"
#include "Types.h"

void run_IncreasingRange(string refName ="", string storingName ="", string storedName = "", string run0 = "", string run1 = "", string run2 = "", string run3 ="", string run4 ="", string run5 ="") {

  gROOT->LoadMacro("DistTrend.cxx++");
  gROOT->LoadMacro("EffTrend.cxx++");
  gROOT->LoadMacro("EfficiencyMonitor.cc++");

  bool doOnlyPlot = kFALSE;
  TChain * chain = new TChain("DTTree");

  if(run0 != "" && run1 == ""){
    chain->Add(run0.c_str());
  }
  else if(run0 != "" && run1 != "" && run2 ==""){
    chain->Add(run0.c_str());
    chain->Add(run1.c_str());
  }
  else if(run0 != "" && run1 != "" && run2 !="" && run3 ==""){  
    chain->Add(run0.c_str());
    chain->Add(run1.c_str());
    chain->Add(run2.c_str());
  }
  else if(run0 != "" && run1 != "" && run2 !="" && run3 != "" && run4 ==""){  
    chain->Add(run0.c_str());
    chain->Add(run1.c_str());
    chain->Add(run2.c_str());
    chain->Add(run3.c_str());
  }
  else if(run0 != "" && run1 != "" && run2 !="" && run3 != "" && run4 != "" && run5 ==""){  
    chain->Add(run0.c_str());
    chain->Add(run1.c_str());
    chain->Add(run2.c_str());
    chain->Add(run3.c_str());
    chain->Add(run4.c_str());
  }
  else if(run0 != "" && run1 != "" && run2 !="" && run3 != "" && run4 != "" && run5 !=""){  
    chain->Add(run0.c_str());
    chain->Add(run1.c_str());
    chain->Add(run2.c_str());
    chain->Add(run3.c_str());
    chain->Add(run4.c_str());
    chain->Add(run5.c_str());
  }
    else if(run0=="" && storedName != "") {
    doOnlyPlot=kTRUE;
  }

  else {cout<<"Problems with files "<<endl;    exit(1);}
 
  context incrCont;
 
  Var InsLumi;
  Var IntLumi;
  Var BunchX;
  Var BckGr;
  Var Empty;
  
  
  //initialize variable struct  Var object, object name, json file, do eff plots, do bkg plots, do projections, variable for projections, is external variable

  variables::initVar(Empty,"Empty","variablesSetting.json",0,0); 
  variables::initVar(BckGr,"Bkg","variablesSetting.json",0,0);  
  variables::initVar(InsLumi,"InsLumi","variablesSetting.json",0,0);
  variables::initVar(IntLumi,"IntLumi","variablesSetting.json",1,0,1);

  incrCont.name      = "Increasing";
  incrCont.var       = { {IntLumi.name,IntLumi},{BckGr.name,BckGr},{Empty.name,Empty},{InsLumi.name,InsLumi}};
  incrCont.webFolder = "~/www/DT";

  incrCont.wwCanvas  = 1200;
  incrCont.whCanvas  = 400;
 
  incrCont.legx1 = 0.79;
  incrCont.legy1 = 0.72;
  incrCont.legx2 = 0.89;
  incrCont.legy2 = 0.89;

  incrCont.titleX = 0.4;
  incrCont.titleY = 0.09;

  
  EfficiencyMonitor *eff = new EfficiencyMonitor(incrCont,chain,refName.c_str(),storingName,storedName,doOnlyPlot);
  
  std::ifstream ifs(("data/DeadList/DeadList_"+refName).c_str(), std::ios::ate); // std::ios::ate means open at end
  
  if(!doOnlyPlot){
    if(ifs.tellg() <= 0){
      cout<<"Dead cell list named DeadList_"+refName+" doesn't exist in data/DeadList or is empty/\nStarting pre-loop to produce it"<<endl;
      eff->PreLoop();
    }
    else {
      cout<<"Dead cell list named DeadList_"+refName+" found in data/DeadList/\nIt will be used instead of producing it"<<endl;
    }
    eff->Loop();
    cout<<"Write"<<endl;
    eff->write();
  }
  
  cout<<"Plot"<<endl;
  eff->plot();
  eff->close();
}
