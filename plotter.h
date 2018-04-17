#ifndef plotter_h
#define plotter_h

#include <TROOT.h>
#include <TFile.h>
#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include "TEfficiency.h"
#include <sys/stat.h>
#include <dirent.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/lexical_cast.hpp>
struct context{

  Int_t nVar;
  std::string name;
  std::vector<vector<Double_t > >  slices;
  std::vector<std::string >  varTitle; 
  std::vector<std::string >  varName; 
  std::vector<std::string >  varLabel; 

  string webFolder;
  
};

struct stat st;

class plotter {
public :
  
  TFile          *fIn;
  TFile          *fOut;

  context dataCont; //struct with the setting
  bool isSliceChanged = 0;
  string LumiFileName;

  Float_t TotLumi = 0;
  vector<vector<vector<TEfficiency* > > > Eff_phiMBWh;
  vector<vector<vector<TEfficiency* > > > EffA_phiMBWh;
  
  vector<vector<vector<TEfficiency* > > > Eff_theMBWh;
  vector<vector<vector<TEfficiency* > > > EffA_theMBWh;

  vector<vector<TEfficiency* > >          Eff_phiMB4Top;
  vector<vector<TEfficiency* > >          EffA_phiMB4Top;
  
  vector<TEfficiency* >                   Eff_phiMB4Bot;
  vector<TEfficiency* >                   EffA_phiMB4Bot;
    
  vector<vector<vector<TH2F* > > > Hist_MBWh;
  vector<vector<TH2F* > >          Hist_MB4Top;
  vector<TH2F* >                   Hist_MB4Bot;
  
  vector<vector<vector<TProfile* > > > Gr_MBWh;
  vector<vector<TProfile* > >          Gr_MB4Top;
  vector<TProfile* >                   Gr_MB4Bot;
  
  plotter(context extDataCont, std::string inFileName = "", std::string outFileName="", string LumiFileName_ = "");

  void  write();  
  void  plot(string dateName = "test");  
  void  setPlots();
  void  setAddBins();
  int   getdir(string dir, vector<string> &files);   // Take a vector of list and pusch inside the files inside a directory. Used for lumi per run files.                 

  TEfficiency * setEffBin(     TEfficiency *effIn, Float_t MaxErr = 0.1); // Change bin size to adjuste the uncertanty below MaxErr
  TEfficiency * addEff(        TEfficiency *eff1,  TEfficiency *eff2 );   // Concatenate two teff. Not used so far
  TEfficiency * setEffRun(     TEfficiency *eff1 );                       // From variable bin size to identical bin size. Used for run number.
  TEfficiency * addBins(       TEfficiency *effIn, int ivar);             // Method to add bins above range to a TEfficiency. Used to values variables such run number.
  TEfficiency * getIntLumiEff( TEfficiency *effIn, Int_t ivar);           // Get integrated lumi from run number. It updates the value from the ones stored.
  TH2F *        getIntLumiHisto(  TH2F *hIn, Int_t ivar);
  TH1F *        getPaintHisto(    TEfficiency *effIn, int ivar, bool doLumi);
  TH2F *        set2DHistoBin(    TH2F *hIn, Float_t MaxErr = 0.1);
  TH2F *        set2DHistoBinRun( TH2F *hIn, int ivar);
  TH2F *        add2DHistoBinRun( TH2F *hIn, int ivar);

  float getLumiRun(string Run);

  //  TH1F * SetHistoLimits( const TH1 *hIn);

  ~plotter();  
};


//#ifdef plotter_cxx
plotter::plotter(context extDataCont, std::string inFileName, std::string outFileName, std::string LumiFileName_){
  
  if(inFileName!=""){
    fIn  = new TFile (("data/results/"+inFileName+".root").c_str());
    if (!fIn){ cout<<"File In doesn't exist or can't be open"<<endl;    exit(1);}
    cout<<"Take objects from file "<<inFileName+".root"<<endl;
  }
  else  cout<<"Create new objects and new file "<<(outFileName+".root")<<endl;
 
  dataCont = extDataCont;
  fOut     = new TFile (("data/results/"+outFileName+".root").c_str(),"RECREATE"); //FIXME. RECREATWE->CREATE


  LumiFileName = "DT_"+LumiFileName_;

  for (int ivar=0; ivar<dataCont.nVar; ivar++){
    
    Eff_phiMBWh.push_back(    vector<vector<TEfficiency*> > ());  
    EffA_phiMBWh.push_back(   vector<vector<TEfficiency*> > ());  
    Eff_theMBWh.push_back(    vector<vector<TEfficiency*> > ());  
    EffA_theMBWh.push_back(   vector<vector<TEfficiency*> > ());  
    Eff_phiMB4Top.push_back(  vector<TEfficiency*>  ());  
    EffA_phiMB4Top.push_back( vector<TEfficiency*>  ());  
    Hist_MBWh.push_back(      vector<vector<TH2F*> > ());  
    Hist_MB4Top.push_back(    vector<TH2F*>  ());  
    Gr_MBWh.push_back(        vector<vector<TProfile*> > ());  
    Gr_MB4Top.push_back(      vector<TProfile*>  ());  

    int nPoints = dataCont.slices[ivar].size();
   
    for (int iwh=0; iwh<5; iwh++){
          
      Eff_phiMBWh[ivar].push_back(vector<TEfficiency*> ());  
      EffA_phiMBWh[ivar].push_back(vector<TEfficiency*> ());  
      
      Eff_theMBWh[ivar].push_back(vector<TEfficiency*> ());         
      EffA_theMBWh[ivar].push_back(vector<TEfficiency*> ());         
      
      Hist_MBWh[ivar].push_back(vector<TH2F*> ());         
      Gr_MBWh[ivar].push_back(vector<TProfile*> ());         
      
      for (int ist=0; ist<4; ist++){
	
	if(inFileName==""){
	  Eff_phiMBWh[ivar][iwh].push_back( new TEfficiency(("Eff"+dataCont.varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
							     "St"+std::to_string(ist)).c_str(),"",nPoints-1,dataCont.slices[ivar].data()));

	  EffA_phiMBWh[ivar][iwh].push_back( new TEfficiency(("EffA"+dataCont.varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
							      "St"+std::to_string(ist)).c_str(),"",nPoints-1,dataCont.slices[ivar].data()));

	  Hist_MBWh[ivar][iwh].push_back(new TH2F(("Hist"+dataCont.varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
						   "St"+std::to_string(ist)).c_str(),"",dataCont.slices[ivar].size()-1,
						  dataCont.slices[ivar].data(),dataCont.slices[2].size()-1,dataCont.slices[2].data()));	   
	}
	else{
	  Eff_phiMBWh[ivar][iwh].push_back( (TEfficiency*) fIn->Get(("Eff"+dataCont.varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
								     "St"+std::to_string(ist)).c_str()));	 	 
	  EffA_phiMBWh[ivar][iwh].push_back( (TEfficiency*) fIn->Get(("EffA"+dataCont.varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
								      "St"+std::to_string(ist)).c_str()));	 	 
	  Hist_MBWh[ivar][iwh].push_back( (TH2F*) fIn->Get(("Hist"+dataCont.varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
							    "St"+std::to_string(ist)).c_str()));	 	 
	}

	Gr_MBWh[ivar][iwh].push_back( new TProfile());
	if (ist!=3){
	  if(inFileName==""){ 
	    Eff_theMBWh[ivar][iwh].push_back( new TEfficiency(("EffThe"+dataCont.varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
							       "St"+std::to_string(ist)).c_str(),"",nPoints-1,dataCont.slices[ivar].data()));	 
	    EffA_theMBWh[ivar][iwh].push_back( new TEfficiency(("EffAThe"+dataCont.varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
								"St"+std::to_string(ist)).c_str(),"",nPoints-1,dataCont.slices[ivar].data()));	 
	  }
	  else{
	    Eff_theMBWh[ivar][iwh].push_back(  (TEfficiency*) fIn->Get(("EffThe"+dataCont.varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
									"St"+std::to_string(ist)).c_str()));	 
	    EffA_theMBWh[ivar][iwh].push_back(  (TEfficiency*) fIn->Get(("EffAThe"+dataCont.varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
									 "St"+std::to_string(ist)).c_str()));
	  }
	}
      }
      if(inFileName==""){
	Eff_phiMB4Top[ivar].push_back( new TEfficiency(("Eff"+dataCont.varName[ivar]+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
						       "",nPoints-1,dataCont.slices[ivar].data()));	          
	EffA_phiMB4Top[ivar].push_back( new TEfficiency(("EffA"+dataCont.varName[ivar]+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
							"",nPoints-1,dataCont.slices[ivar].data()));	          
	Hist_MB4Top[ivar].push_back( new TH2F(("Hist"+dataCont.varName[ivar]+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
					      "",dataCont.slices[ivar].size()-1,dataCont.slices[ivar].data(),dataCont.slices[2].size()-1,dataCont.slices[2].data()));
      }
      
      else{
	Eff_phiMB4Top[ivar].push_back(  (TEfficiency*) fIn->Get(("Eff"+dataCont.varName[ivar]+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
	EffA_phiMB4Top[ivar].push_back(  (TEfficiency*) fIn->Get(("EffA"+dataCont.varName[ivar]+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
	Hist_MB4Top[ivar].push_back( (TH2F*) fIn->Get(("Hist"+dataCont.varName[ivar]+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
      }
      Gr_MB4Top[ivar].push_back( new TProfile());
    }
    
    if(inFileName==""){
      Eff_phiMB4Bot.push_back( new TEfficiency(("Eff"+dataCont.varName[ivar]+"_MBBot").c_str(),"",nPoints-1,dataCont.slices[ivar].data()));	      
      EffA_phiMB4Bot.push_back( new TEfficiency(("EffA"+dataCont.varName[ivar]+"_MBBot").c_str(),"",nPoints-1,dataCont.slices[ivar].data()));	      
      Hist_MB4Bot.push_back( new TH2F(("Hist"+dataCont.varName[ivar]+"_MBBot").c_str(),"",dataCont.slices[ivar].size()-1,
				      dataCont.slices[ivar].data(),dataCont.slices[2].size()-1,dataCont.slices[2].data()));	      
    }
    else{
      Eff_phiMB4Bot.push_back(  (TEfficiency*) fIn->Get(("Eff"+dataCont.varName[ivar]+"_MBBot").c_str()));	      
      EffA_phiMB4Bot.push_back(  (TEfficiency*) fIn->Get(("EffA"+dataCont.varName[ivar]+"_MBBot").c_str()));	      
      Hist_MB4Bot.push_back( (TH2F*) fIn->Get(("Hist"+dataCont.varName[ivar]+"_MBBot").c_str()));
    }
    Gr_MB4Bot.push_back( new TProfile());       
  } 

  if(inFileName!=""){
    setAddBins();
  }  
}

plotter::~plotter(){

}  

void plotter::write(){

   fOut->cd();

   for (int ivar=0; ivar<dataCont.nVar; ivar++) {
     for (int iwh=0; iwh<5; iwh++){
       for (int ist=0; ist<4; ist++){
	 Eff_phiMBWh[ivar][iwh][ist]->Write();
	 EffA_phiMBWh[ivar][iwh][ist]->Write();
	 Hist_MBWh[ivar][iwh][ist]->Write();
	 if(ist!=3){
	   Eff_theMBWh[ivar][iwh][ist]->Write();
	   EffA_theMBWh[ivar][iwh][ist]->Write();
	 }
       }
       Eff_phiMB4Top[ivar][iwh]->Write();
       EffA_phiMB4Top[ivar][iwh]->Write();
       Hist_MB4Top[ivar][iwh]->Write();
     }

     Eff_phiMB4Bot[ivar]->Write();
     EffA_phiMB4Bot[ivar]->Write();
     Hist_MB4Bot[ivar]->Write();
   }
}

void plotter::plot(string dateName){
  
  setPlots();
  gStyle->SetOptStat(0);

  if(stat("plot/",&st) != 0)  system("mkdir plot/");
  if(stat(("plot/"+dateName).c_str(),&st) != 0)   system(("mkdir plot/"+dateName).c_str());

  //To check the metod stat works with final dir. It doesn't work with ~ .To be deleted 
  /* cout<<"stat((dataCont.webFolder+'/'+dateName+'/').c_str(),&st)"<<stat((dataCont.webFolder+"/"+dateName+"/").c_str(),&st)<<endl; */ 
  /* cout<<"stat((dataCont.webFolder'/'+dateName).c_str(),&st)"<<stat((dataCont.webFolder+"/"+dateName).c_str(),&st)<<endl; */
  /* cout<<"stat(('/'+dataCont.webFolder'/'+dateName).c_str(),&st)"<<stat(("/"+dataCont.webFolder+"/"+dateName).c_str(),&st)<<endl; */

  if(stat((dataCont.webFolder+"/"+dateName+"/").c_str(),&st) != 0){
    cout<<(dataCont.webFolder+"/"+dateName+"/").c_str()<<endl;
    system(("mkdir "+dataCont.webFolder+"/"+dateName).c_str());
    system(("cp "+dataCont.webFolder+"/index.php " +dataCont.webFolder+"/"+dateName).c_str());
  }

  if(stat((dataCont.webFolder+"/"+dateName+"/Efficiency").c_str(),&st) != 0){
  system(("mkdir "+dataCont.webFolder+"/"+dateName+"/Efficiency").c_str());
  system(("cp "+dataCont.webFolder+"/index.php " +dataCont.webFolder+"/"+dateName+"/Efficiency").c_str());
  }
  
  if(stat((dataCont.webFolder+"/"+dateName+"/Background").c_str(),&st) != 0){
  system(("mkdir "+dataCont.webFolder+"/"+dateName+"/Background").c_str());
  system(("cp "+dataCont.webFolder+"/index.php " +dataCont.webFolder+"/"+dateName+"/Background").c_str());
  }

  fOut->cd();

  for (int ivar=0; ivar<dataCont.nVar; ivar++) {
    for (int iwh=0; iwh<5; iwh++){
      for (int ist=0; ist<4; ist++){
	if(iwh!=4){
	  Eff_phiMBWh[ivar][iwh][ist]->SetLineColor(iwh+1);
	  Eff_phiMBWh[ivar][iwh][ist]->SetMarkerColor(iwh+1);

	  Hist_MBWh[ivar][iwh][ist]->SetLineColor(iwh+1);
	  Hist_MBWh[ivar][iwh][ist]->SetMarkerColor(iwh+1);
	  
	  if(ist!=3){
	    Eff_theMBWh[ivar][iwh][ist]->SetLineColor(iwh+1);
	    Eff_theMBWh[ivar][iwh][ist]->SetMarkerColor(iwh+1);
	  }
	}
	else{ 
	  Eff_phiMBWh[ivar][iwh][ist]->SetLineColor(6);
	  Eff_phiMBWh[ivar][iwh][ist]->SetMarkerColor(6);
	  
	  Hist_MBWh[ivar][iwh][ist]->SetLineColor(6);
	  Hist_MBWh[ivar][iwh][ist]->SetMarkerColor(6);
	  
	  if(ist!=3){
	    Eff_theMBWh[ivar][iwh][ist]->SetLineColor(6);
	    Eff_theMBWh[ivar][iwh][ist]->SetMarkerColor(6);
	  }
	}
	Eff_phiMBWh[ivar][iwh][ist]->SetMarkerStyle(20);
	Hist_MBWh[ivar][iwh][ist]->SetMarkerStyle(20);
	//profile of 2D histograms
	Gr_MBWh[ivar][iwh][ist] = Hist_MBWh[ivar][iwh][ist]->ProfileX();

	//	for(int bin = 1; bin <=Hist_MBWh[ivar][iwh][ist]->GetNbinsX(); bin++)  Gr_MBWh[ivar][iwh][ist]->GetXaxis()->SetBinLabel(bin,(to_string(static_cast<int>(dataCont.slices[ivar][bin-1]))).c_str());  //To add run number labels
  
	
	if(ist!=3){
	  Eff_theMBWh[ivar][iwh][ist]->SetMarkerStyle(20);
	  
	}
      }
      Eff_phiMB4Top[ivar][iwh]->SetMarkerStyle(20);
      Hist_MB4Top[ivar][iwh]->SetMarkerStyle(20);
      
      if(iwh!=4){

	Eff_phiMB4Top[ivar][iwh]->SetLineColor(iwh+1);
	Eff_phiMB4Top[ivar][iwh]->SetMarkerColor(iwh+1);
	Hist_MB4Top[ivar][iwh]->SetLineColor(iwh+1);
	Hist_MB4Top[ivar][iwh]->SetMarkerColor(iwh+1);
      }
      else {
	Eff_phiMB4Top[ivar][iwh]->SetLineColor(6);
	Eff_phiMB4Top[ivar][iwh]->SetMarkerColor(6);
	
	Hist_MB4Top[ivar][iwh]->SetLineColor(6);
	Hist_MB4Top[ivar][iwh]->SetMarkerColor(6);
      }
      Gr_MB4Top[ivar][iwh] = Hist_MB4Top[ivar][iwh]->ProfileX();
      //      for(int bin = 1; bin <=Hist_MB4Top[ivar][iwh]->GetNbinsX(); bin++)  Gr_MB4Top[ivar][iwh]->GetXaxis()->SetBinLabel(bin,(to_string(static_cast<int>(dataCont.slices[ivar][bin-1]))).c_str()); 
    }
    Eff_phiMB4Bot[ivar]->SetMarkerStyle(20);
    Eff_phiMB4Bot[ivar]->SetLineColor(1);
    Eff_phiMB4Bot[ivar]->SetMarkerColor(1);
    Hist_MB4Bot[ivar]->SetMarkerStyle(20);
    Hist_MB4Bot[ivar]->SetLineColor(1);
    Hist_MB4Bot[ivar]->SetMarkerColor(1);
    
    Gr_MB4Bot[ivar] = Hist_MB4Bot[ivar]->ProfileX();
    // for(int bin = 1; bin <=Hist_MB4Bot[ivar]->GetNbinsX(); bin++)  Gr_MB4Bot[ivar]->GetXaxis()->SetBinLabel(bin,(to_string(static_cast<int>(dataCont.slices[ivar][bin-1]))).c_str()); 
  }
  

  for (int ivar=0; ivar<dataCont.nVar; ivar++){
    
    for (int ist=0; ist<4; ist++){

      //phi
      TCanvas *cPhiMB = new TCanvas(("cPhiMB"+(std::to_string(ist+1))+dataCont.varName[ivar]).c_str());
      
      TH1F * hPaint = new TH1F();
      
      Eff_phiMBWh[ivar][0][ist]->SetTitle(("MB"+(std::to_string(ist+1))+" eff vs "+dataCont.varTitle[ivar]+";"+dataCont.varLabel[ivar]+";Eff").c_str());
      
      if(dataCont.varName[ivar]=="Run"){
	hPaint  =  getPaintHisto(Eff_phiMBWh[ivar][0][ist], ivar, false);
	hPaint->SetMinimum(0.89);
	hPaint->SetMaximum(1.02); 
	hPaint->Draw();
	Eff_phiMBWh[ivar][0][ist]->Draw("samep");
      }
      else Eff_phiMBWh[ivar][0][ist]->Draw("ap");
      
      gPad->Update(); 
      auto graph = Eff_phiMBWh[ivar][0][ist]->GetPaintedGraph(); 
      
      cPhiMB->Update();
      graph->SetMinimum(0.89);
      graph->SetMaximum(1.02);

      TLegend * legPhiMB = new TLegend(0.75,0.75,0.9,0.9);
      for (int iwh=0; iwh<5; iwh++){
	Eff_phiMBWh[ivar][iwh][ist]->Draw("samep");
	legPhiMB->AddEntry(Eff_phiMBWh[ivar][iwh][ist],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
      }
      
      legPhiMB->Draw("same");
      cPhiMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+dataCont.varName[ivar]+".png").c_str());
      
      
      if(dataCont.varName[ivar]=="Run"){

	Eff_phiMBWh[ivar][0][ist] = getIntLumiEff(Eff_phiMBWh[ivar][0][ist],ivar);
	Eff_phiMBWh[ivar][0][ist]->Draw("ap");
	for (int iwh=1; iwh<5; iwh++){
	  Eff_phiMBWh[ivar][iwh][ist] = getIntLumiEff(Eff_phiMBWh[ivar][iwh][ist],ivar);
	  Eff_phiMBWh[ivar][iwh][ist]->Draw("samep");
	}

	gPad->Update(); 
	auto graph = Eff_phiMBWh[ivar][0][ist]->GetPaintedGraph(); 
	
	graph->SetMinimum(0.89);
	graph->SetMaximum(1.02);
       
	legPhiMB->Draw("same");
	cPhiMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+dataCont.varName[ivar]+"_var.png").c_str());
      }



      if(dataCont.varName[ivar]=="Run"){


      }

      delete   cPhiMB;
      
      //theta
      
      if(ist!=3){ 
	TCanvas *cTheMB = new TCanvas(("cTheMB"+(std::to_string(ist+1))+dataCont.varName[ivar]).c_str());
	
	
	TH1F * hPaint = new TH1F();
	
	Eff_theMBWh[ivar][0][ist]->SetTitle(("MB"+(std::to_string(ist+1))+" theta eff vs "+dataCont.varTitle[ivar]+";"+dataCont.varLabel[ivar]+";Eff").c_str());
	
	if(dataCont.varName[ivar]=="Run"){
	  hPaint  =  getPaintHisto(Eff_theMBWh[ivar][0][ist], ivar, false);
	  hPaint->SetMinimum(0.89);
	  hPaint->SetMaximum(1.02); 
	  hPaint->Draw();
	  Eff_theMBWh[ivar][0][ist]->Draw("samep");
	}
	else{
	  Eff_theMBWh[ivar][0][ist]->Draw("ap");
	}
	///	
	gPad->Update(); 
	auto graph = Eff_theMBWh[ivar][0][ist]->GetPaintedGraph(); 
	
	graph->SetMinimum(0.89);
	graph->SetMaximum(1.02); 
	gPad->Update(); 
	
	
	TLegend * legTheMB = new TLegend(0.75,0.75,0.9,0.9);
	for (int iwh=0; iwh<5; iwh++){
	  Eff_theMBWh[ivar][iwh][ist]->Draw("samep");
	  legTheMB->AddEntry(Eff_theMBWh[ivar][iwh][ist],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
	}
	legTheMB->Draw("same");
	cTheMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"TheEffVs"+dataCont.varName[ivar]+".png").c_str());
	
	
	
	if(dataCont.varName[ivar]=="Run"){
	  
	  Eff_theMBWh[ivar][0][ist] = getIntLumiEff(Eff_theMBWh[ivar][0][ist],ivar);
	  Eff_theMBWh[ivar][0][ist]->Draw("ap");
	  for (int iwh=1; iwh<5; iwh++){
	    Eff_theMBWh[ivar][iwh][ist] = getIntLumiEff(Eff_theMBWh[ivar][iwh][ist],ivar);
	    Eff_theMBWh[ivar][iwh][ist]->Draw("samep");
	  }
	  
	  gPad->Update(); 
	  auto graph = Eff_theMBWh[ivar][0][ist]->GetPaintedGraph(); 
	  
	  graph->SetMinimum(0.89);
	  graph->SetMaximum(1.02);
	  
	  legTheMB->Draw("same");
	  cTheMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"TheEffVs"+dataCont.varName[ivar]+"_var.png").c_str());
	}
	delete 	cTheMB;
      }
    }
    
    //phi MB4 Top
    
    TCanvas *cPhiMB4Top = new TCanvas(("cPhiMB4Top"+dataCont.varName[ivar]).c_str());
    
    TH1F * hPaint = new TH1F();    
    Eff_phiMB4Top[ivar][0]->SetTitle(("MB4Top eff vs "+dataCont.varTitle[ivar]+" eff vs "+dataCont.varTitle[ivar]+";"+dataCont.varLabel[ivar]+";Eff").c_str());
    
    if(dataCont.varName[ivar]=="Run"){
      hPaint  =  getPaintHisto(Eff_phiMB4Top[ivar][0], ivar, false);
      hPaint->SetMinimum(0.89);
      hPaint->SetMaximum(1.02); 
      hPaint->Draw();
      Eff_phiMB4Top[ivar][0]->Draw("samep");
    }
    else{
      Eff_phiMB4Top[ivar][0]->Draw("ap");
    }
    gPad->Update(); 
    
    auto graph = Eff_phiMB4Top[ivar][0]->GetPaintedGraph(); 
    graph->SetMinimum(0.89);
    graph->SetMaximum(1.02); 
    
    gPad->Update(); 
    
    TLegend * legPhiMB4Top = new TLegend(0.75,0.75,0.9,0.9);
    for (int iwh=0; iwh<5; iwh++){
      Eff_phiMB4Top[ivar][iwh]->Draw("samep");
      legPhiMB4Top->AddEntry(Eff_phiMB4Top[ivar][iwh],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
    }
    legPhiMB4Top->Draw("same");
    cPhiMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4TopPhiEffVs"+dataCont.varName[ivar]+".png").c_str());
    
    
    if(dataCont.varName[ivar]=="Run"){

	Eff_phiMB4Top[ivar][0] = getIntLumiEff(Eff_phiMB4Top[ivar][0],ivar);
	Eff_phiMB4Top[ivar][0]->Draw("ap");
	for (int iwh=1; iwh<5; iwh++){
	  Eff_phiMB4Top[ivar][iwh] = getIntLumiEff(Eff_phiMB4Top[ivar][iwh],ivar);
	  Eff_phiMB4Top[ivar][iwh]->Draw("samep");
	}

	gPad->Update(); 
	auto graph = Eff_phiMB4Top[ivar][0]->GetPaintedGraph(); 
	
	graph->SetMinimum(0.89);
	graph->SetMaximum(1.02);
       
	legPhiMB4Top->Draw("same");
	cPhiMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4TopPhiEffVs"+dataCont.varName[ivar]+"_var.png").c_str());
      }
      
    delete   cPhiMB4Top ;

    //phi MB4 Bot
    TCanvas *cPhiMB4Bot = new TCanvas(("cPhiMB4Bot"+dataCont.varName[ivar]).c_str());
       
    Eff_phiMB4Bot[ivar]->SetTitle(("MB4Bot eff vs "+dataCont.varTitle[ivar]+" eff vs "+dataCont.varTitle[ivar]+";"+dataCont.varLabel[ivar]+";Eff").c_str());

    if(dataCont.varName[ivar]=="Run"){
      hPaint  =  getPaintHisto(Eff_phiMB4Bot[ivar], ivar, false);
      hPaint->SetMinimum(0.89);
      hPaint->SetMaximum(1.02); 
      hPaint->Draw();
      Eff_phiMB4Bot[ivar]->Draw("samep");
    }
    else  Eff_phiMB4Bot[ivar]->Draw("ap"); 
    
    gPad->Update(); 
    graph = Eff_phiMB4Bot[ivar]->GetPaintedGraph(); 
    graph->SetMinimum(0.89);
    graph->SetMaximum(1.02); 
    
    gPad->Update(); 
    
    cPhiMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4BotPhiEffVs"+dataCont.varName[ivar]+".png").c_str());
    
    
    if(dataCont.varName[ivar]=="Run"){

	Eff_phiMB4Bot[ivar] = getIntLumiEff(Eff_phiMB4Bot[ivar],ivar);
	Eff_phiMB4Bot[ivar]->Draw("ap");

	gPad->Update(); 
	auto graph = Eff_phiMB4Bot[ivar]->GetPaintedGraph(); 
	
	graph->SetMinimum(0.89);
	graph->SetMaximum(1.02);
       
	cPhiMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4BotPhiEffVs"+dataCont.varName[ivar]+"_var.png").c_str());
      }
    delete cPhiMB4Bot ;  

  }
  
  //Graph bkg
  
  for (int ivar=0; ivar<2; ivar++){
    for (int ist=0; ist<4; ist++){
      //phi
      TCanvas *cMB = new TCanvas(("cMB"+(std::to_string(ist+1))+dataCont.varName[ivar]).c_str());
      Gr_MBWh[ivar][0][ist]->SetMaximum(35);
      Gr_MBWh[ivar][0][ist]->SetTitle(("MB"+(std::to_string(ist+1))+" bkg vs "+dataCont.varTitle[ivar]+";"+dataCont.varLabel[ivar]+";bkg").c_str());
      Gr_MBWh[ivar][0][ist]->Draw("E1");
      

      TLegend * legMB = new TLegend(0.75,0.75,0.9,0.9);
      for (int iwh=0; iwh<5; iwh++){
	Gr_MBWh[ivar][iwh][ist]->Draw("sameE1");
	legMB->AddEntry(Gr_MBWh[ivar][iwh][ist],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
      }
      legMB->Draw("same");
      cMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB"+(std::to_string(ist+1))+"BkgVs"+dataCont.varName[ivar]+".png").c_str());
      delete  cMB ;
     }


    // MB4 Top
    TCanvas *cMB4Top = new TCanvas(("cMB4Top"+dataCont.varName[ivar]).c_str());
    Gr_MB4Top[ivar][0]->SetMaximum(35);
    Gr_MB4Top[ivar][0]->SetTitle(("MB4Top bkg vs "+dataCont.varTitle[ivar]+";"+dataCont.varLabel[ivar]+";bkg").c_str());
    Gr_MB4Top[ivar][0]->Draw("E1");
    
    TLegend * legMB4Top = new TLegend(0.75,0.75,0.9,0.9);
    for (int iwh=0; iwh<5; iwh++){
      Gr_MB4Top[ivar][iwh]->Draw("sameE1");
      legMB4Top->AddEntry(Gr_MB4Top[ivar][iwh],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
    }
    legMB4Top->Draw("same");
    cMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4TopBkgVs"+dataCont.varName[ivar]+".png").c_str());
    delete cMB4Top ;
    
    // MB4 Bot
    TCanvas *cMB4Bot = new TCanvas(("cMB4Bot"+dataCont.varName[ivar]).c_str());
    
    Gr_MB4Bot[ivar]->SetTitle(("MB4Bot bkg vs "+dataCont.varTitle[ivar]+dataCont.varTitle[ivar]+";"+dataCont.varLabel[ivar]+";bkg").c_str());
    Gr_MB4Bot[ivar]->Draw("E1");
    
    cMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4BotBkgVs"+dataCont.varName[ivar]+".png").c_str());
    delete cMB4Bot; 
  }      
}

// Reorganize bins on order to have a un uncertainty per bin less than MaxErr. To be improved in order to check merged bins uncertainties.
TEfficiency * plotter::setEffBin( TEfficiency *effIn, Float_t MaxErr){

  Int_t nBins  = effIn->GetTotalHistogram()->GetNbinsX();
  const  TArrayD *arr = effIn->GetTotalHistogram()->GetXaxis()->GetXbins(); 

  std::vector<Double_t> xBins;
  std::vector<Double_t> totCont;
  std::vector<Double_t> passCont;

  bool isFirst = true;
  bool isEmpty = true;

  Int_t maxBin = -1;
  bool  foundMax = false;

  for(int bin = 1; bin<=nBins; bin++){

    //loop backward to find last non empty bin
    if(!foundMax && effIn->GetTotalHistogram()->GetBinContent(nBins - bin)>0.){
      maxBin = nBins - bin;
      foundMax = true;
      break;
    }
  }

  for(int bin = 1; bin<=maxBin; bin++){
    if(isEmpty){
      if(effIn->GetTotalHistogram()->GetBinContent(bin)) isEmpty = false; 
      else continue;
    }
    if(isFirst || ( (effIn->GetTotalHistogram()->GetBinError(bin)/effIn->GetTotalHistogram()->GetBinContent(bin) < MaxErr) || bin ==maxBin)){
      isFirst = false;     
      xBins.push_back(arr->GetArray()[bin-1]); 
      totCont.push_back( effIn->GetTotalHistogram()->GetBinContent(bin));
      passCont.push_back(effIn->GetPassedHistogram()->GetBinContent(bin));
    }
    else{
      // if uncertainty is greater than MaxErr sum content in previous bin
      totCont.back()+=effIn->GetTotalHistogram()->GetBinContent(bin);
      passCont.back()+=effIn->GetPassedHistogram()->GetBinContent(bin);
    }
  }


  //  xBins.insert(xBins.begin(),arrX->GetArray()[firstBin-1]);
  //Cont.insert(Cont.begin(),0);

  //insert last  value to have the right edge of the last  bin
  xBins.push_back(arr->GetArray()[maxBin]);

  //insert extra value to have one last empty bin
  xBins.push_back(arr->GetArray()[maxBin+1]);

  totCont.push_back(0.);
  passCont.push_back(0.);
  
  TEfficiency * effNew = new TEfficiency(effIn->GetName(),effIn->GetTitle(),xBins.size()-1,xBins.data());

  for(uint bin = 1; bin<xBins.size(); bin++){
    effNew->SetTotalEvents(bin, totCont[bin-1]);
    effNew->SetPassedEvents(bin, passCont[bin-1]);
  }
  return effNew;
}


/* TH1F *  plotter::SetHistoLimits( const TH1 *hIn){ */
  
/*   Int_t nBins  = hIn->GetNbinsX(); */
  
/*   static bool foundMin; */
/*   static bool foundMax; */
  
/*   Int_t minBin = 0; */
/*   Int_t maxBin = -1; */
  
/*   TH1F * hNew = (TH1F*)hIn->Clone(); */

/*   for(int bin = 1; bin<=nBins; bin++){ */
/*     if (!foundMin && hIn->GetBinContent(bin)>0. ){ */
/*       minBin = bin; */
/*       foundMin = true; */
/*     } */
    
/*     if(!foundMax && hIn->GetBinContent(nBins - bin)>0.){ */
/*       maxBin = nBins - bin; */
/*       foundMax = true; */
/*     } */
/*   } */
/*   hNew->GetXaxis()->SetRange(minBin,maxBin); */
/*   return hNew; */
/* } */


// Reorganize bins on order to have a un uncertainty per bin less than MaxErr. To be improved in order to check merged bins uncertainties.
TH2F * plotter::set2DHistoBin( TH2F *hIn, Float_t MaxErr){
  
  Int_t nXBins  = hIn->GetNbinsX();       
  Int_t nYBins  = hIn->GetNbinsY();       

  const  TArrayD *arrX = hIn->GetXaxis()->GetXbins();
  const  TArrayD *arrY = hIn->GetYaxis()->GetXbins();

  std::vector<Double_t> xBins;
  std::vector<vector<Double_t>> Cont;
  
  Double_t  Err = 0;
  Double_t IntX = 0;
  

  Int_t maxBin = -1;
  Int_t firstBin = -1;

  bool  foundMax = false;
  bool  isEmpty  = true;
  bool  isFirst  = true;

  for(int bin = 1; bin<=nXBins; bin++){
    //loop backward to find last non empty column 
    if(!foundMax && hIn->Integral((nXBins - bin),(nXBins - bin),0,-1 ) >0.){
      maxBin = nXBins - bin;
      foundMax = true;
      break;
    }
  }
  //the array has a one item more then the number of bins. 

  for(int bin = 1; bin<=maxBin; bin++){  

    IntX = hIn->IntegralAndError(bin,bin,0,-1,Err);
    if(isEmpty){
      if(IntX) {
	isEmpty =false;
	firstBin = (bin==1) ? 1 : bin - 1;
	}
      else continue;
    }    

    if(isFirst ||( (Err/IntX < MaxErr) || (bin == maxBin) ) ){ 
	isFirst = false;
	Cont.push_back(vector<Double_t>());
	xBins.push_back(arrX->GetArray()[bin-1]);
	for(int biny = 1; biny<=nYBins; biny++){  
	  Cont.back().push_back(hIn->GetBinContent(bin,biny));
	}
    }
    else{
      for(int biny = 1; biny<=nYBins; biny++){  
	(Cont.back())[biny-1]+=hIn->GetBinContent(bin,biny);
      }
    }
  }
  
  xBins.insert(xBins.begin(),arrX->GetArray()[firstBin-1]);
  Cont.insert(Cont.begin(),vector<Double_t>());

  for(int biny = 1; biny<=nYBins; biny++) Cont[0].push_back(0.);
  
  xBins.push_back(arrX->GetArray()[maxBin]);
  xBins.push_back(arrX->GetArray()[maxBin+1]);
  Cont.push_back(vector<Double_t>());
  for(int biny = 1; biny<=nYBins; biny++) Cont.back().push_back(0.);

  TH2F * hNew = new TH2F((string(hIn->GetName())).c_str(),"",xBins.size()-1,xBins.data(),nYBins,arrY->GetArray());
  for(uint binx = 1; binx<=xBins.size()-1; binx++){  
    for(int biny = 1; biny<=nYBins; biny++){  
      hNew->SetBinContent(binx,biny, Cont[binx-1][biny-1]);
    }
  }
  return hNew;
}

void plotter::setPlots(){


   for (int ivar=0; ivar<dataCont.nVar; ivar++) {
     for (int iwh=0; iwh<5; iwh++){
       for (int ist=0; ist<4; ist++){
	 if(dataCont.varName[ivar]=="Run"){
	   Eff_phiMBWh[ivar][iwh][ist]  = setEffRun(Eff_phiMBWh[ivar][iwh][ist]);
	   EffA_phiMBWh[ivar][iwh][ist] = setEffRun(EffA_phiMBWh[ivar][iwh][ist]);
	   //	   Hist_MBWh[ivar][iwh][ist]    = set2DHistoBinRun(Hist_MBWh[ivar][iwh][ist], ivar);
	   Hist_MBWh[ivar][iwh][ist]    = getIntLumiHisto(Hist_MBWh[ivar][iwh][ist], ivar);
	 }
	 else{
	   Eff_phiMBWh[ivar][iwh][ist]  = setEffBin(Eff_phiMBWh[ivar][iwh][ist]);
	   EffA_phiMBWh[ivar][iwh][ist] = setEffBin(EffA_phiMBWh[ivar][iwh][ist]);
	 }
	 if(dataCont.varName[ivar] != "Run" &&  dataCont.varName[ivar] != "bkg") Hist_MBWh[ivar][iwh][ist]    = set2DHistoBin(Hist_MBWh[ivar][iwh][ist]);
	 
	 if(ist!=3){
	   if(dataCont.varName[ivar]=="Run"){
	     Eff_theMBWh[ivar][iwh][ist]  = setEffRun(Eff_theMBWh[ivar][iwh][ist]);
	     EffA_theMBWh[ivar][iwh][ist] = setEffRun(EffA_theMBWh[ivar][iwh][ist]);
	     
	     //Hist_theMBWh[ivar][iwh][ist] = getIntLumiHisto(Hist_theMBWh[ivar][iwh][ist], ivar); //theta hist don't exit yet
	     //Hist_theMBWh[ivar][iwh][ist] = set2DHistoBinRun(Hist_theMBWh[ivar][iwh][ist], ivar);
	   }
	   else{
	     Eff_theMBWh[ivar][iwh][ist]  = setEffBin(Eff_theMBWh[ivar][iwh][ist]);
	     EffA_theMBWh[ivar][iwh][ist] = setEffBin(EffA_theMBWh[ivar][iwh][ist]);
	   }
	 }
       }
       if(dataCont.varName[ivar]=="Run"){
	 Eff_phiMB4Top[ivar][iwh]  = setEffRun(Eff_phiMB4Top[ivar][iwh]);
	 EffA_phiMB4Top[ivar][iwh] = setEffRun(EffA_phiMB4Top[ivar][iwh]);
	 //Hist_MB4Top[ivar][iwh]    = set2DHistoBinRun(Hist_MB4Top[ivar][iwh], ivar);
	 Hist_MB4Top[ivar][iwh]    = getIntLumiHisto(Hist_MB4Top[ivar][iwh], ivar);
       }
       else{
	 Eff_phiMB4Top[ivar][iwh]  = setEffBin(Eff_phiMB4Top[ivar][iwh]);
	 EffA_phiMB4Top[ivar][iwh] = setEffBin(EffA_phiMB4Top[ivar][iwh]);
       }
       if(dataCont.varName[ivar] != "Run" &&  dataCont.varName[ivar] != "bkg")  Hist_MB4Top[ivar][iwh] = set2DHistoBin(Hist_MB4Top[ivar][iwh]);
     
     }
     if(dataCont.varName[ivar]=="Run"){

       Eff_phiMB4Bot[ivar]  = setEffRun(Eff_phiMB4Bot[ivar]);
       EffA_phiMB4Bot[ivar] = setEffRun(EffA_phiMB4Bot[ivar]);
       //Hist_MB4Bot[ivar]    = set2DHistoBinRun(Hist_MB4Bot[ivar], ivar);
       Hist_MB4Bot[ivar]    = getIntLumiHisto(Hist_MB4Bot[ivar], ivar);
     }
     else{
       Eff_phiMB4Bot[ivar]  = setEffBin(Eff_phiMB4Bot[ivar]);
       EffA_phiMB4Bot[ivar] = setEffBin(EffA_phiMB4Bot[ivar]);
     }
     if(dataCont.varName[ivar] != "Run" &&  dataCont.varName[ivar] != "bkg")  Hist_MB4Bot[ivar] = set2DHistoBin(Hist_MB4Bot[ivar]);
   }
}

//not used so farg
TEfficiency * plotter::addEff( TEfficiency *eff1, TEfficiency *eff2 ){

  Int_t nbins1 = eff1->GetTotalHistogram()->GetNbinsX();
  Int_t nbins2 = eff2->GetTotalHistogram()->GetNbinsX();

  TEfficiency * effNew = new TEfficiency(eff1->GetName(),eff1->GetTitle(),nbins1+nbins2,eff1->GetTotalHistogram()->GetBinLowEdge(1),eff2->GetTotalHistogram()->GetBinLowEdge(nbins2+1));

 //  cout<<effNew->GetTotalHistogram()->GetNbinsX()<< " "<<effNew->GetTotalHistogram()->GetBinLowEdge(effNew->GetTotalHistogram()->GetNbinsX()+1)<<endl;
  for(int bin = 1; bin<=nbins1; bin++){
    effNew->SetTotalEvents(bin,  eff1->GetTotalHistogram()->GetBinContent(bin));
    effNew->SetPassedEvents(bin, eff1->GetPassedHistogram()->GetBinContent(bin));
  }

  for(int bin = nbins1; bin<=nbins1+nbins2; bin++){
    effNew->SetTotalEvents(bin+nbins1,  eff2->GetTotalHistogram()->GetBinContent(bin));
    effNew->SetPassedEvents(bin+nbins1, eff2->GetPassedHistogram()->GetBinContent(bin));
  }
  return effNew;
}


TEfficiency * plotter::setEffRun( TEfficiency *eff1 ){
  
  Int_t nbins = eff1->GetTotalHistogram()->GetNbinsX();
  TEfficiency * effNew = new TEfficiency(eff1->GetName(),eff1->GetTitle(),nbins,0,nbins);
  
  auto graph = effNew->GetPaintedGraph(); 
  
  for(int bin = 1; bin<=nbins; bin++){
    effNew->SetTotalEvents(bin,  eff1->GetTotalHistogram()->GetBinContent(bin));
    effNew->SetPassedEvents(bin, eff1->GetPassedHistogram()->GetBinContent(bin));
  }    
  return effNew;
}


void plotter::setAddBins(){

  for (int ivar=0; ivar<dataCont.nVar; ivar++) {
    if(dataCont.varName[ivar]=="Run"){
      for (int iwh=0; iwh<5; iwh++){
	for (int ist=0; ist<4; ist++){
	  Eff_phiMBWh[ivar][iwh][ist]  = addBins(Eff_phiMBWh[ivar][iwh][ist],ivar);
	  EffA_phiMBWh[ivar][iwh][ist] = addBins(EffA_phiMBWh[ivar][iwh][ist],ivar);	  
	  Hist_MBWh[ivar][iwh][ist] = add2DHistoBinRun(Hist_MBWh[ivar][iwh][ist],ivar);

	  if(ist!=3){
	    Eff_theMBWh[ivar][iwh][ist]  = addBins(Eff_theMBWh[ivar][iwh][ist],ivar);
	    EffA_theMBWh[ivar][iwh][ist] = addBins(EffA_theMBWh[ivar][iwh][ist],ivar);
	    // Hist_theMBWh[ivar][iwh][ist] = add2DHistoBinRun(Hist_theMBWh[ivar][iwh][ist],ivar);
	  }
	  Eff_phiMB4Top[ivar][iwh]  = addBins(Eff_phiMB4Top[ivar][iwh],ivar);
	  EffA_phiMB4Top[ivar][iwh] = addBins(EffA_phiMB4Top[ivar][iwh],ivar);
	  Hist_MB4Top[ivar][iwh]    = add2DHistoBinRun(Hist_MB4Top[ivar][iwh],ivar);	
}
	Eff_phiMB4Bot[ivar]  = addBins(Eff_phiMB4Bot[ivar],ivar);
	EffA_phiMB4Bot[ivar] = addBins(EffA_phiMB4Bot[ivar],ivar);
	Hist_MB4Bot[ivar]    = add2DHistoBinRun(Hist_MB4Bot[ivar],ivar);	
      }    
    }
  } 
}


TH2F * plotter::add2DHistoBinRun( TH2F *hIn, int ivar){
  
  Int_t nXBins  = hIn->GetNbinsX();       
  Int_t nYBins  = hIn->GetNbinsY();       
  
  const  TArrayD *arrY = hIn->GetYaxis()->GetXbins();

  std::vector<Double_t> xBins;
  xBins  =  dataCont.slices[ivar];

  //  for(int bin = 1; bin<=nXBins; bin++)  
  //  xBins.push_back(bin);
  
  TH2F * hNew = new TH2F((string(hIn->GetName())).c_str(),"",xBins.size()-1,xBins.data(),nYBins,arrY->GetArray());
  for(int binx = 1; binx<=nXBins; binx++){  
    for(int biny = 1; biny<=nYBins; biny++){  
      hNew->SetBinContent(binx,biny, hIn->GetBinContent(binx,biny));
    }
  }
  
  for(int bin = 1; bin <= nXBins; bin++){
    hNew->GetXaxis()->SetBinLabel(bin,(to_string(static_cast<int>(dataCont.slices[ivar][bin-1]))).c_str()); 
  }  
  return hNew;
}


//hot
TEfficiency * plotter::addBins( TEfficiency *effIn, int ivar){

  std::vector<Double_t> xBins = {};
  Int_t nBins  = effIn->GetTotalHistogram()->GetNbinsX();

  // run slices has to changed just once
  if(!isSliceChanged){
      
    const  TArrayD *arr = effIn->GetTotalHistogram()->GetXaxis()->GetXbins(); 
 
    for(int bin = 1; bin<=nBins+1; bin++) {
      xBins.push_back(arr->GetArray()[bin-1]); 
    }
    
    //Delete last elements that was added just to have an extra item to have the last bin            
    xBins.erase(xBins.end()-1);
   
    //insert in xBins the new runs
    xBins.insert(xBins.end(), dataCont.slices[ivar].begin(), dataCont.slices[ivar].end());
    
    // update slice with old run numbers
    dataCont.slices[ivar]= xBins;
    isSliceChanged= true; 

 }
  else{
    // dataCont.slices[ivar] has already been changed and now Xbins can be simply makes equal to  dataCont.slices[ivar]
    xBins  =  dataCont.slices[ivar];
  }

  TEfficiency * effNew = new TEfficiency(effIn->GetName(),effIn->GetTitle(),xBins.size()-1,xBins.data());
  
  for(int bin = 1; bin<=nBins; bin++){
    effNew->SetTotalEvents( bin, effIn->GetTotalHistogram()->GetBinContent(bin));
    effNew->SetPassedEvents(bin, effIn->GetPassedHistogram()->GetBinContent(bin));
  }
  return effNew;
}


float plotter::getLumiRun(string run){

  vector<string> files = vector<string>();
  string dir = "data/IntLumi/";

  //get list of file in InLumi directory
  getdir(dir,files);

  if(std::find(files.begin(), files.end(),("data/IntLumi/DT_Run"+run+".csv") ) != files.end()) {
    //push in front of files the file with the name of the run to be fater in case it exist.
    files.insert(files.begin(),"data/IntLumi/DT_Run"+run+".csv"); //check
  }

  Float_t intLumi = 0;
  bool isFound = false; 

  for (unsigned int i = 0;i < files.size();i++) {

     if(files[i].find("~")!=string::npos) continue; 
     else if(files[i].find(".csv")==string::npos) continue; 

     ifstream inFile;
     string line;

     inFile.open((dir+files[i]).c_str()); 

     if(!inFile){
       cout << "Unable to open file "<<dir<<files[i]<<endl;      
       exit(1);
     }

     size_t pos;
     vector<string> runLine;

     while(inFile.good())
       {
	 getline(inFile,line); // get line from file
	 if(line[0]=='#' || line.empty()) continue;
	 pos=line.find(run); // search
	 boost::split(runLine,line,boost::is_any_of(","));

	 string test = runLine[5];
	 intLumi +=   stof(test);

	 // string::npos is returned if string is not found 
	 if(pos!=string::npos){
	   isFound =true;
	   break;

	 }
       }
     if(isFound) break;
     else intLumi=0;
  }

  TotLumi+=intLumi;
  return TotLumi;  

}

// take a vector of list and pusch inside the files inside a directory. Used for lumi per run files.
int plotter::getdir(string dir, vector<string> &files){
  DIR *dp;
  struct dirent *dirp;
  if((dp = opendir(dir.c_str())) == NULL) {
    cout << "Error(" << errno << ") opening " << dir << endl;
    return errno;
  }
  
  while ((dirp = readdir(dp)) != NULL) {
    files.push_back(string(dirp->d_name));
  }
  closedir(dp);
  return 0;
}


//Create histogram with equal bins to be used as graphic plot instead of tgraph.
TH1F * plotter::getPaintHisto( TEfficiency *effIn, int ivar, bool doLumi){
  TH1F * hPaint = new TH1F();
  
  hPaint = (TH1F*)effIn->GetTotalHistogram()->Clone();

  hPaint->Reset();
  hPaint->SetTitle(effIn->GetTitle());

  Int_t nBins = hPaint->GetNbinsX();
  
  
  if(doLumi){ 
    for(int bin = 1; bin <= nBins; bin++){
      hPaint->GetXaxis()->SetBinLabel(bin, (boost::str(boost::format("%.2f") % getLumiRun(to_string(static_cast<int>(dataCont.slices[ivar][bin-1]))))).c_str());
    }
    hPaint->GetXaxis()->SetTitle("Int.Lumi. pb^{-1}");
    TotLumi = 0; // reset total lumi
  }
  else{
    for(int bin = 1; bin <= nBins; bin++){
      hPaint->GetXaxis()->SetBinLabel(bin,(to_string(static_cast<int>(dataCont.slices[ivar][bin-1]))).c_str()); 
    }
  }

  return hPaint;
}


TEfficiency * plotter::getIntLumiEff( TEfficiency *effIn, Int_t ivar){

  std::vector<Double_t> xBins = {};
  Int_t nBins  = effIn->GetTotalHistogram()->GetNbinsX();

  for(int bin = 1; bin <= nBins; bin++){
    xBins.push_back(getLumiRun(to_string(static_cast<int>(dataCont.slices[ivar][bin-1]))));
  }

  TotLumi = 0; // reset total lumi

  //To be updated using SetTotalEvents TEfficiency method.

  TH1F * hTot  = new TH1F( (string(effIn->GetName())+"_Tot").c_str(),effIn->GetTitle(),30*xBins.size(),xBins.at(0),xBins.back()+xBins.back()/10.);
  TH1F * hPass = new TH1F( (string(effIn->GetName())+"_Pass").c_str(),effIn->GetTitle(),30*xBins.size(),xBins.at(0),xBins.back()+xBins.back()/10.);

  for(int bin = 1; bin<=nBins; bin++){
    hTot->AddBinContent( hTot->FindBin(xBins[bin-1])  , effIn->GetTotalHistogram()->GetBinContent(bin));
    hPass->AddBinContent(hPass->FindBin(xBins[bin-1]) , effIn->GetPassedHistogram()->GetBinContent(bin));
  }

  hTot->GetXaxis()->SetTitle("Int.Lumi. pb^{-1}");
  TEfficiency * effNew = new TEfficiency(*hPass,*hTot);

  hTot->Delete();
  hPass->Delete();

  effNew->SetLineColor(effIn->GetLineColor());
  effNew->SetMarkerStyle(effIn->GetMarkerStyle());
  effNew->SetMarkerColor(effIn->GetMarkerColor());
  return effNew;

}


TH2F * plotter::getIntLumiHisto( TH2F *hIn, Int_t ivar){

  Int_t nXBins  = hIn->GetNbinsX();       
  Int_t nYBins  = hIn->GetNbinsY();       

  const  TArrayD *arrY = hIn->GetYaxis()->GetXbins();

  std::vector<Double_t> xBins = {};


  for(int bin = 1; bin <= nXBins+1; bin++){
    xBins.push_back(getLumiRun(to_string(static_cast<int>(dataCont.slices[ivar][bin-1]))));
  }

  TotLumi = 0; // reset total lumi

  TH2F * hNew = new TH2F((string(hIn->GetName())).c_str(),"",30*xBins.size(),xBins.at(0)-xBins.at(0)/10.,xBins.back()+xBins.back()/10.,nYBins,arrY->GetArray());

  for(int binx = 1; binx<=nXBins; binx++){  
    for(int biny = 1; biny<=nYBins; biny++){  
      hNew->AddBinContent(hNew->FindBin(xBins[binx-1],biny), hIn->GetBinContent(binx,biny));
    }
  }

  hNew->GetXaxis()->SetTitle("Int.Lumi. pb^{-1}");

  hNew->SetLineColor(hIn->GetLineColor());
  hNew->SetMarkerStyle(hIn->GetMarkerStyle());
  hNew->SetMarkerColor(hIn->GetMarkerColor());
  return hNew;
}


TH2F * plotter::set2DHistoBinRun( TH2F *hIn, int ivar){
  
  Int_t nXBins  = hIn->GetNbinsX();       
  Int_t nYBins  = hIn->GetNbinsY();       

  const  TArrayD *arrY = hIn->GetYaxis()->GetXbins();
  const  TArrayD *arrX = hIn->GetXaxis()->GetXbins();
  std::vector<Double_t> xBins;

  for(int bin = 1; bin<=nXBins; bin++){  
    xBins.push_back(arrX->GetArray()[bin-1]);
  }

  TH2F * hNew = new TH2F((string(hIn->GetName())).c_str(),"",xBins.size(),xBins[0],xBins.back(),nYBins,arrY->GetArray());

  for(uint binx = 1; binx<=xBins.size(); binx++){  
    for(int biny = 1; biny<=nYBins; biny++){  
      hNew->SetBinContent(binx,biny, hIn->GetBinContent(binx,biny));
    }
  }

  for(int bin = 1; bin <=nXBins; bin++){
    hNew->GetXaxis()->SetBinLabel(bin,(to_string(static_cast<int>(dataCont.slices[ivar][bin-1]))).c_str()); //del
    }

  return hNew;
}

#endif
