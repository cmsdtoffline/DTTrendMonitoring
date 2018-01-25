//#include <stdio>
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "EfficiencyMonitor.h"


int run(string name ="", string run0 = "", string run1 = "", string run2 = "", string run3 ="", string run4 ="", string run5 ="", string run6 ="") {

	gROOT->LoadMacro("EfficiencyMonitorLumiPU.C++");

	TFile * fout = new TFile("temp.root","recreate");
	TTree * DTTree =  new TTree();

	if(run0 != "" && run1 == ""){
	  TFile *_file0 = TFile::Open(run0.c_str());
	  DTTree =  (TTree*)_file0->Get("DTTree");
	  
	}
	else if(run0 != "" && run1 != "" && run2 ==""){
	  TFile *_file0 = TFile::Open(run0.c_str());
	  TFile *_file1 = TFile::Open(run1.c_str());
	  TTree *DTTree0 =  (TTree*)_file0->Get("DTTree");
	  TTree *DTTree1 =  (TTree*)_file1->Get("DTTree");

	  fout->cd();
	  TList *list = new TList; 
	  list->Add(DTTree0); 
	  list->Add(DTTree1); 
  
	  DTTree = TTree::MergeTrees(list); 
	  DTTree->SetName("DTTree"); 
	  // DTTree->Write();
	}

	else if(run0 != "" && run1 != "" && run2 !="" && run3 == ""){
	 
	  TFile *_file0 = TFile::Open(run0.c_str());
	  TFile *_file1 = TFile::Open(run1.c_str());
	  TFile *_file2 = TFile::Open(run2.c_str());

	  TTree *DTTree0 =  (TTree*)_file0->Get("DTTree");
	  TTree *DTTree1 =  (TTree*)_file1->Get("DTTree");
	  TTree *DTTree2 =  (TTree*)_file2->Get("DTTree");

	  fout->cd();
	  TList *list = new TList; 
	  list->Add(DTTree0); 
	  list->Add(DTTree1); 
	  list->Add(DTTree2); 
  
	  list->Print();

	  DTTree = TTree::MergeTrees(list); 
	  DTTree->SetName("DTTree"); 
	  //DTTree->Write();
	}
	
	else if(run0 != "" && run1 != "" && run2 !="" && run3 != "" && run4==""){
	 
	  TFile *_file0 = TFile::Open(run0.c_str());
	  TFile *_file1 = TFile::Open(run1.c_str());
	  TFile *_file2 = TFile::Open(run2.c_str());
	  TFile *_file3 = TFile::Open(run3.c_str());

	  TTree *DTTree0 =  (TTree*)_file0->Get("DTTree");
	  TTree *DTTree1 =  (TTree*)_file1->Get("DTTree");
	  TTree *DTTree2 =  (TTree*)_file2->Get("DTTree");
	  TTree *DTTree3 =  (TTree*)_file3->Get("DTTree");

	  fout->cd();
	  TList *list = new TList; 
	  list->Add(DTTree0); 
	  list->Add(DTTree1); 
	  list->Add(DTTree2); 
	  list->Add(DTTree3); 
  
	  list->Print();

	  DTTree = TTree::MergeTrees(list); 
	  DTTree->SetName("DTTree"); 
	  //DTTree->Write();
	}
	
	else{
	  std::cout<<"Error: something wrong with the file names ";
	  abort();
	}
	EfficiencyMonitor *eff = new EfficiencyMonitor(DTTree,name.c_str());
	
	//eff->PreLoop();
	eff->Loop();
	//eff->PostLoop();
	system("rm temp.root");	
	return 1;
}
