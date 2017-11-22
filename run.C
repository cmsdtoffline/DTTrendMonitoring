//#include <stdio>
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"

#include "EfficiencyMonitor.h"


int run() {

	gROOT->LoadMacro("EfficiencyMonitorRun.C");
	EfficiencyMonitor *eff = new EfficiencyMonitor();

	//eff->PreLoop();
	eff->Loop();
	//eff->PostLoop();

return 1;
}
