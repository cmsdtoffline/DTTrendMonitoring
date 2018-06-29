#ifndef plotter_h
#define plotter_h

#include "Hist.h"
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
#include <TGaxis.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include "TEfficiency.h"
#include <sys/stat.h>
#include <dirent.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include "Utilities.h"
#include <TLegendEntry.h>

struct stat st;

class plotter {
public :
  
  TFile          *fIn;
  TFile          *fOut;

  context dataCont; //struct with the setting
  bool    isSliceChanged = 0;
  string  LumiFileName;

  Float_t TotLumi = 0;

  map<string,vector<vector<TEfficiency* > > > Eff_phiMBWh;
  map<string,vector<vector<TEfficiency* > > > EffA_phiMBWh;
  
  map<string,vector<vector<TEfficiency* > > > Eff_theMBWh;
  map<string,vector<vector<TEfficiency* > > > EffA_theMBWh;

  map<string,vector<TEfficiency* > >          Eff_phiMB4Top;
  map<string,vector<TEfficiency* > >          EffA_phiMB4Top;
  
  map<string,TEfficiency* >                   Eff_phiMB4Bot;
  map<string,TEfficiency* >                   EffA_phiMB4Bot;
    
  map<string,vector<vector<Hist* > > >        Hist_MBWh;
  map<string,vector<vector<Hist* > > >        Hist_SegMBWh;

  map<string,vector<Hist* > >                 Hist_MB4Top;
  map<string,vector<Hist* > >                 Hist_SegMB4Top;

  map<string,Hist* >                          Hist_MB4Bot;
  map<string,Hist* >                          Hist_SegMB4Bot;
  
  plotter(context extDataCont, std::string inFileName = "", std::string outFileName="", string LumiFileName_ = "");

  void   write();
  void   plot(string dateName = "test");
  void   setPlots();
  void   setAddBins();
  int    getdir(string dir, vector<string> &files);   // Take a vector of list and push inside the files inside a directory. Used for lumi per run files.

  float  getLumiRun(string Run);    // Get  integrated luminosity from the begining of the plot range.

  TEfficiency * setEffBin(     TEfficiency *effIn, Float_t MaxErr = 0.1); // Change bin size to adjuste the uncertanty below MaxErr
  TEfficiency * addEff(        TEfficiency *eff1,  TEfficiency *eff2 );   // Concatenate two teff. Not used so far
  TEfficiency * setEffRun(     TEfficiency *eff1 );                       // From variable bin size to identical bin size. Used for run number.
  TEfficiency * addBins(       TEfficiency *effIn, vector<Double_t> slice);             // Method to add bins above range to a TEfficiency. Used to values variables such run number.
  void addBinsSlice(const TArrayD * arr, vector<double> & slices); 
  TEfficiency * getIntLumiEff( TEfficiency *effIn, vector<Double_t>);           // Get integrated lumi from run number. It updates the value from the ones stored.
  TH1F *        getPaintHisto(    TEfficiency *effIn, vector<Double_t> , bool doLumi);

  ~plotter();
};


//#ifdef plotter_cxx
plotter::plotter(context extDataCont, std::string inFileName, std::string outFileName, std::string LumiFileName_){
  
  if(inFileName!=""){
    fIn  = new TFile (("data/results/"+inFileName+".root").c_str());
    if (!fIn){ cout<<"File In doesn't exist or can't be open"<<endl;    abort();}
    cout<<"Take objects from file "<<inFileName+".root"<<endl;
  }
  else  cout<<"Create new objects and new file "<<(outFileName+".root")<<endl;
 
  dataCont = extDataCont;
  fOut     = new TFile (("data/results/"+outFileName+".root").c_str(),"RECREATE"); //FIXME. RECREATWE->CREATE


  LumiFileName = "DT_"+LumiFileName_;

  for(auto const& ivar : dataCont.var){ 

    Eff_phiMBWh[ivar.first]    = vector<vector<TEfficiency*> > ();
    EffA_phiMBWh[ivar.first]   = vector<vector<TEfficiency*> > ();
    Eff_theMBWh[ivar.first]    = vector<vector<TEfficiency*> > ();
    EffA_theMBWh[ivar.first]   = vector<vector<TEfficiency*> > ();
    Eff_phiMB4Top[ivar.first]  = vector<TEfficiency*>  ();
    EffA_phiMB4Top[ivar.first] = vector<TEfficiency*>  ();


    Hist_MBWh[ivar.first]      = vector<vector<Hist*> > ();
    Hist_SegMBWh[ivar.first]   = vector<vector<Hist*> > ();

    Hist_MB4Top[ivar.first]    = vector<Hist*>  ();
    Hist_SegMB4Top[ivar.first] = vector<Hist*>  ();
    
    int nPoints = ivar.second.slice.size();

    for (int iwh=0; iwh<5; iwh++){

      Eff_phiMBWh[ivar.first].push_back(vector<TEfficiency*> ());
      EffA_phiMBWh[ivar.first].push_back(vector<TEfficiency*> ());
      
      Eff_theMBWh[ivar.first].push_back(vector<TEfficiency*> ());
      EffA_theMBWh[ivar.first].push_back(vector<TEfficiency*> ());
      
      Hist_MBWh[ivar.first].push_back(vector<Hist*> ());
      Hist_SegMBWh[ivar.first].push_back(vector<Hist*> ());

      for (int ist=0; ist<4; ist++){
	if(inFileName==""){

	  //cout<<"var "<<ivar.second.name<<" "<<ivar.second.doEff<<" "<<ivar.second.doBkg<<" "<<ivar.second.slice.size()<<endl;

	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh].push_back( new TEfficiency(("Eff"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
							     "St"+std::to_string(ist)).c_str(),"",nPoints-1,ivar.second.slice.data()));
	  if(ivar.second.doEff) EffA_phiMBWh[ivar.first][iwh].push_back( new TEfficiency(("EffA"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
							      "St"+std::to_string(ist)).c_str(),"",nPoints-1,ivar.second.slice.data()));
	  if(ivar.second.doBkg) Hist_MBWh[ivar.first][iwh].push_back(new Hist(("Hist"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
						   "St"+std::to_string(ist)).c_str(),"",ivar.second.slice.size()-1,
						  ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data()));
	  if(ivar.second.doBkg) Hist_SegMBWh[ivar.first][iwh].push_back(new Hist(("HistSeg"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
						      "St"+std::to_string(ist)).c_str(),"",ivar.second.slice.size()-1,
						     ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data()));

	}
	else{
	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh].push_back( (TEfficiency*) fIn->Get(("Eff"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
								     "St"+std::to_string(ist)).c_str()));
	  if(ivar.second.doEff) EffA_phiMBWh[ivar.first][iwh].push_back( (TEfficiency*) fIn->Get(("EffA"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
								      "St"+std::to_string(ist)).c_str()));
	  if(ivar.second.doBkg) Hist_MBWh[ivar.first][iwh].push_back( (Hist*) fIn->Get(("Hist"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
							    "St"+std::to_string(ist)).c_str()));
	  if(ivar.second.doBkg) Hist_SegMBWh[ivar.first][iwh].push_back( (Hist*) fIn->Get(("HistSeg"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
							    "St"+std::to_string(ist)).c_str()));	  
	}

	if (ist!=3){
	  if(inFileName==""){
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh].push_back( new TEfficiency(("EffThe"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
							       "St"+std::to_string(ist)).c_str(),"",nPoints-1,ivar.second.slice.data()));
	    if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh].push_back( new TEfficiency(("EffAThe"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
								"St"+std::to_string(ist)).c_str(),"",nPoints-1,ivar.second.slice.data()));
	  }
	  else{
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh].push_back(  (TEfficiency*) fIn->Get(("EffThe"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
													  "St"+std::to_string(ist)).c_str()));
	    if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh].push_back(  (TEfficiency*) fIn->Get(("EffAThe"+ivar.second.name+"_MBWh"+std::to_string(iwh-2)+
									 "St"+std::to_string(ist)).c_str()));
	  }
	}
      }
      if(inFileName==""){

	if(ivar.second.doEff) Eff_phiMB4Top[ivar.first].push_back( new TEfficiency(("Eff"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
						       "",nPoints-1,ivar.second.slice.data()));
	if(ivar.second.doEff) EffA_phiMB4Top[ivar.first].push_back( new TEfficiency(("EffA"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
							"",nPoints-1,ivar.second.slice.data()));
	if(ivar.second.doBkg) Hist_MB4Top[ivar.first].push_back( new Hist(("Hist"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
					      "",ivar.second.slice.size()-1,ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data()));

	if(ivar.second.doBkg) Hist_SegMB4Top[ivar.first].push_back( new Hist(("HistSeg"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
						 "",ivar.second.slice.size()-1,ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data()));   
   }
      
      else{
	if(ivar.second.doEff) Eff_phiMB4Top[ivar.first].push_back(  (TEfficiency*) fIn->Get(("Eff"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
	if(ivar.second.doEff) EffA_phiMB4Top[ivar.first].push_back(  (TEfficiency*) fIn->Get(("EffA"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str()));

	if(ivar.second.doBkg) Hist_MB4Top[ivar.first].push_back( (Hist*) fIn->Get(("Hist"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
	if(ivar.second.doBkg) Hist_SegMB4Top[ivar.first].push_back( (Hist*) fIn->Get(("HistSeg"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
      }

    }
    
    if(inFileName==""){
      if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first] =  new TEfficiency(("Eff"+ivar.second.name+"_MBBot").c_str(),"",nPoints-1,ivar.second.slice.data());
      if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first] =  new TEfficiency(("EffA"+ivar.second.name+"_MBBot").c_str(),"",nPoints-1,ivar.second.slice.data());
      if(ivar.second.doBkg) Hist_MB4Bot[ivar.first] =  new Hist(("Hist"+ivar.second.name+"_MBBot").c_str(),"",ivar.second.slice.size()-1,
				      ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data());

      if(ivar.second.doBkg) Hist_SegMB4Bot[ivar.first] =  new Hist(("HistSeg"+ivar.second.name+"_MBBot").c_str(),"",ivar.second.slice.size()-1,
				      ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data());
    }
    else{
      if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first] =   (TEfficiency*) fIn->Get(("Eff"+ivar.second.name+"_MBBot").c_str());
      if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first] =   (TEfficiency*) fIn->Get(("EffA"+ivar.second.name+"_MBBot").c_str());
    
      if(ivar.second.doBkg) Hist_MB4Bot[ivar.first] =  (Hist*) fIn->Get(("Hist"+ivar.second.name+"_MBBot").c_str());
      if(ivar.second.doBkg) Hist_SegMB4Bot[ivar.first] =  (Hist*) fIn->Get(("HistSeg"+ivar.second.name+"_MBBot").c_str());
    }
  }

  if(inFileName!=""){
    setAddBins();
  }
}

plotter::~plotter(){

}

void plotter::write(){

   fOut->cd();
   for(auto const& ivar : dataCont.var) {
     for (int iwh=0; iwh<5; iwh++){
       for (int ist=0; ist<4; ist++){
	 if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->Write();
	 if(ivar.second.doEff) EffA_phiMBWh[ivar.first][iwh][ist]->Write();
	 if(ivar.second.doBkg) Hist_MBWh[ivar.first][iwh][ist]->Write();
	 if(ivar.second.doBkg) Hist_SegMBWh[ivar.first][iwh][ist]->Write();

	 if(ist!=3){
	   if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->Write();
	   if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh][ist]->Write();
	 }
       }
       if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->Write();
       if(ivar.second.doEff) EffA_phiMB4Top[ivar.first][iwh]->Write();
       if(ivar.second.doBkg) Hist_MB4Top[ivar.first][iwh]->Write();
       if(ivar.second.doBkg) Hist_SegMB4Top[ivar.first][iwh]->Write();
     }
     if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]->Write();
     if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first]->Write();
     if(ivar.second.doBkg) Hist_MB4Bot[ivar.first]->Write();
     if(ivar.second.doBkg) Hist_SegMB4Bot[ivar.first]->Write();
   }
}

void plotter::plot(string dateName){
  
  setPlots();
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(4);

  if(stat("plot/",&st) != 0)  system("mkdir plot/");
  if(stat(("plot/"+dateName).c_str(),&st) != 0)   system(("mkdir plot/"+dateName).c_str());

  //To check the metod stat works with final dir. It doesn't work with ~ .To be deleted
  /* cout<<"stat((dataCont.webFolder+'/'+dateName+'/').c_str(),&st)"<<stat((dataCont.webFolder+"/"+dateName+"/").c_str(),&st)<<endl; */
  /* cout<<"stat((dataCont.webFolder'/'+dateName).c_str(),&st)"<<stat((dataCont.webFolder+"/"+dateName).c_str(),&st)<<endl; */
  /* cout<<"stat(('/'+dataCont.webFolder'/'+dateName).c_str(),&st)"<<stat(("/"+dataCont.webFolder+"/"+dateName).c_str(),&st)<<endl; */

  if(stat((dataCont.webFolder+"/"+dateName+"/").c_str(),&st) != 0){
    //    cout<<(dataCont.webFolder+"/"+dateName+"/").c_str()<<endl;
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

  for(auto const& ivar : dataCont.var) {
    for (int iwh=0; iwh<5; iwh++){
      for (int ist=0; ist<4; ist++){
	if(iwh!=4){
	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->SetLineColor(iwh+1);
	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->SetMarkerColor(iwh+1);

	  if(ivar.second.doBkg) Hist_MBWh[ivar.first][iwh][ist]->SetColor(iwh+1);
	  if(ivar.second.doBkg) Hist_SegMBWh[ivar.first][iwh][ist]->SetColor(iwh+1);
	 
	  if(ist!=3){
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->SetLineColor(iwh+1);
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->SetMarkerColor(iwh+1);
	  }
	}
	else{
	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->SetLineColor(6);
	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->SetMarkerColor(6);
	  
	  if(ivar.second.doBkg) Hist_MBWh[ivar.first][iwh][ist]->SetColor(6);
	  if(ivar.second.doBkg) Hist_SegMBWh[ivar.first][iwh][ist]->SetColor(6);
	  
	  if(ist!=3){
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->SetLineColor(6);
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->SetMarkerColor(6);
	  }
	}

	if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->SetMarkerStyle(20);
	if(ivar.second.doBkg) Hist_MBWh[ivar.first][iwh][ist]->SetMarkerStyle(20);
	if(ivar.second.doBkg) Hist_SegMBWh[ivar.first][iwh][ist]->SetMarkerStyle(20);
	//profile of 2D histograms

	 if(ivar.second.doBkg) Hist_MBWh[ivar.first][iwh][ist]->ProfileX();
	 if(ivar.second.doBkg) Hist_SegMBWh[ivar.first][iwh][ist]->ProfileX();
  
	
	if(ist!=3){
	  if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->SetMarkerStyle(20);	  
	}
      }
      if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->SetMarkerStyle(20);
      if(ivar.second.doBkg) Hist_MB4Top[ivar.first][iwh]->SetMarkerStyle(20);
      if(ivar.second.doBkg) Hist_SegMB4Top[ivar.first][iwh]->SetMarkerStyle(20);
      
      if(iwh!=4){

	if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->SetLineColor(iwh+1);
	if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->SetMarkerColor(iwh+1);
	if(ivar.second.doBkg) Hist_MB4Top[ivar.first][iwh]->SetColor(iwh+1);
	if(ivar.second.doBkg) Hist_SegMB4Top[ivar.first][iwh]->SetColor(iwh+1);

      }
      else {
	if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->SetLineColor(6);
	if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->SetMarkerColor(6);
	
	if(ivar.second.doBkg) Hist_MB4Top[ivar.first][iwh]->SetColor(6);
	if(ivar.second.doBkg) Hist_SegMB4Top[ivar.first][iwh]->SetColor(6);
      }

      if(ivar.second.doBkg) Hist_MB4Top[ivar.first][iwh]->ProfileX();
      if(ivar.second.doBkg) Hist_SegMB4Top[ivar.first][iwh]->ProfileX();
    }

    if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]->SetMarkerStyle(20);
    if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]->SetLineColor(1);
    if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]->SetMarkerColor(1);

    if(ivar.second.doBkg) Hist_MB4Bot[ivar.first]->SetMarkerStyle(20);
    if(ivar.second.doBkg) Hist_MB4Bot[ivar.first]->SetColor(1);
   
    if(ivar.second.doBkg) Hist_SegMB4Bot[ivar.first]->SetMarkerStyle(20);
    if(ivar.second.doBkg) Hist_SegMB4Bot[ivar.first]->SetColor(1);
    
    if(ivar.second.doBkg) Hist_MB4Bot[ivar.first]->ProfileX();
    if(ivar.second.doBkg) Hist_SegMB4Bot[ivar.first]->ProfileX();

  }
  
    for(auto const& ivar : dataCont.var){ 
      if(! ivar.second.doEff) continue;
      for (int ist=0; ist<4; ist++){
      // cout<<"ist "<<ist<<" "<<ivar.second.name<<endl;
      //phi
      TCanvas *cPhiMB = new TCanvas(("cPhiMB"+(std::to_string(ist+1))+ivar.second.name).c_str());
      
      TH1F * hPaint = new TH1F();
      
      Eff_phiMBWh[ivar.first][0][ist]->SetTitle(("MB"+(std::to_string(ist+1))+" eff vs "+ivar.second.Title+";"+ivar.second.Label+";Eff").c_str());
      
      if(ivar.second.name=="Run"){
	hPaint  =  getPaintHisto(Eff_phiMBWh[ivar.first][0][ist], ivar.second.slice, false);
	hPaint->SetMinimum(0.89);
	hPaint->SetMaximum(1.02);
	hPaint->Draw();
	Eff_phiMBWh[ivar.first][0][ist]->Draw("samep");
      }
      else Eff_phiMBWh[ivar.first][0][ist]->Draw("ap");
      
      gPad->Update();
      auto graph = Eff_phiMBWh[ivar.first][0][ist]->GetPaintedGraph();
      
      cPhiMB->Update();
      graph->SetMinimum(0.89);
      graph->SetMaximum(1.02);

      TLegend * legPhiMB = new TLegend(0.75,0.75,0.9,0.9);
      for (int iwh=0; iwh<5; iwh++){
	Eff_phiMBWh[ivar.first][iwh][ist]->Draw("samep");
	legPhiMB->AddEntry(Eff_phiMBWh[ivar.first][iwh][ist],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
      }
      
      legPhiMB->Draw("same");
      cPhiMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+ivar.second.name+".png").c_str());
      
      
      if(ivar.second.name=="Run"){
	Eff_phiMBWh[ivar.first][0][ist] = getIntLumiEff(Eff_phiMBWh[ivar.first][0][ist],ivar.second.slice);
	Eff_phiMBWh[ivar.first][0][ist]->Draw("ap");
	for (int iwh=1; iwh<5; iwh++){
	  Eff_phiMBWh[ivar.first][iwh][ist] = getIntLumiEff(Eff_phiMBWh[ivar.first][iwh][ist],ivar.second.slice);
	  Eff_phiMBWh[ivar.first][iwh][ist]->Draw("samep");
	}

	gPad->Update();
	auto graph = Eff_phiMBWh[ivar.first][0][ist]->GetPaintedGraph();
	
	graph->SetMinimum(0.89);
	graph->SetMaximum(1.02);
       
	legPhiMB->Draw("same");
	cPhiMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+ivar.second.name+"_intLum.png").c_str());
      }

      if(ivar.second.name=="Run"){

      }

      delete   cPhiMB;
      
      //theta
      
      if(ist!=3){
	TCanvas *cTheMB = new TCanvas(("cTheMB"+(std::to_string(ist+1))+ivar.second.name).c_str());
	
	
	TH1F * hPaint = new TH1F();
	
	Eff_theMBWh[ivar.first][0][ist]->SetTitle(("MB"+(std::to_string(ist+1))+" theta eff vs "+ivar.second.Title+";"+ivar.second.Label+";Eff").c_str());
	
	if(ivar.second.name=="Run"){
	  hPaint  =  getPaintHisto(Eff_theMBWh[ivar.first][0][ist], ivar.second.slice, false);
	  hPaint->SetMinimum(0.89);
	  hPaint->SetMaximum(1.02);
	  hPaint->Draw();
	  Eff_theMBWh[ivar.first][0][ist]->Draw("samep");
	}
	else{
	  Eff_theMBWh[ivar.first][0][ist]->Draw("ap");
	}
	

	gPad->Update();
	auto graph = Eff_theMBWh[ivar.first][0][ist]->GetPaintedGraph();
	
	graph->SetMinimum(0.89);
	graph->SetMaximum(1.02);
	gPad->Update();
	
	
	TLegend * legTheMB = new TLegend(0.75,0.75,0.9,0.9);
	for (int iwh=0; iwh<5; iwh++){
	  Eff_theMBWh[ivar.first][iwh][ist]->Draw("samep");
	  legTheMB->AddEntry(Eff_theMBWh[ivar.first][iwh][ist],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
	}
	legTheMB->Draw("same");
	cTheMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"TheEffVs"+ivar.second.name+".png").c_str());
	
	
       
	if(ivar.second.name=="Run"){
	  
	  Eff_theMBWh[ivar.first][0][ist] = getIntLumiEff(Eff_theMBWh[ivar.first][0][ist],ivar.second.slice);
	  Eff_theMBWh[ivar.first][0][ist]->Draw("ap");
	  for (int iwh=1; iwh<5; iwh++){
	    Eff_theMBWh[ivar.first][iwh][ist] = getIntLumiEff(Eff_theMBWh[ivar.first][iwh][ist],ivar.second.slice);
	    Eff_theMBWh[ivar.first][iwh][ist]->Draw("samep");
	  }
	  
	  gPad->Update();
	  auto graph = Eff_theMBWh[ivar.first][0][ist]->GetPaintedGraph();
	  
	  graph->SetMinimum(0.89);
	  graph->SetMaximum(1.02);
	  
	  legTheMB->Draw("same");
	  cTheMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"TheEffVs"+ivar.second.name+"_intLum.png").c_str());
	}
	delete 	cTheMB;
      }
    }
    
    //phi MB4 Top
    
    TCanvas *cPhiMB4Top = new TCanvas(("cPhiMB4Top"+ivar.second.name).c_str());
    
    TH1F * hPaint = new TH1F();
    Eff_phiMB4Top[ivar.first][0]->SetTitle(("MB4Top eff vs "+ivar.second.Title+" eff vs "+ivar.second.Title+";"+ivar.second.Label+";Eff").c_str());
    
    if(ivar.second.name=="Run"){
      hPaint  =  getPaintHisto(Eff_phiMB4Top[ivar.first][0], ivar.second.slice, false);
      hPaint->SetMinimum(0.89);
      hPaint->SetMaximum(1.02);
      hPaint->Draw();
      Eff_phiMB4Top[ivar.first][0]->Draw("samep");
    }
    else{
      Eff_phiMB4Top[ivar.first][0]->Draw("ap");
    }
    gPad->Update();
    
    auto graph = Eff_phiMB4Top[ivar.first][0]->GetPaintedGraph();
    graph->SetMinimum(0.89);
    graph->SetMaximum(1.02);
    
    gPad->Update();
    
    TLegend * legPhiMB4Top = new TLegend(0.75,0.75,0.9,0.9);
    for (int iwh=0; iwh<5; iwh++){
      Eff_phiMB4Top[ivar.first][iwh]->Draw("samep");
      legPhiMB4Top->AddEntry(Eff_phiMB4Top[ivar.first][iwh],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
    }
    legPhiMB4Top->Draw("same");
    cPhiMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4TopPhiEffVs"+ivar.second.name+".png").c_str());
    
    
    if(ivar.second.name=="Run"){

	Eff_phiMB4Top[ivar.first][0] = getIntLumiEff(Eff_phiMB4Top[ivar.first][0],ivar.second.slice);
	Eff_phiMB4Top[ivar.first][0]->Draw("ap");
	for (int iwh=1; iwh<5; iwh++){
	  Eff_phiMB4Top[ivar.first][iwh] = getIntLumiEff(Eff_phiMB4Top[ivar.first][iwh],ivar.second.slice);
	  Eff_phiMB4Top[ivar.first][iwh]->Draw("samep");
	}

	gPad->Update();
	auto graph = Eff_phiMB4Top[ivar.first][0]->GetPaintedGraph();
	
	graph->SetMinimum(0.89);
	graph->SetMaximum(1.02);
       
	legPhiMB4Top->Draw("same");
	cPhiMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4TopPhiEffVs"+ivar.second.name+"_intLum.png").c_str());
      }
      
    delete   cPhiMB4Top ;

    //phi MB4 Bot
    TCanvas *cPhiMB4Bot = new TCanvas(("cPhiMB4Bot"+ivar.second.name).c_str());
       
    Eff_phiMB4Bot[ivar.first]->SetTitle(("MB4Bot eff vs "+ivar.second.Title+" eff vs "+ivar.second.Title+";"+ivar.second.Label+";Eff").c_str());

    if(ivar.second.name=="Run"){
      hPaint  =  getPaintHisto(Eff_phiMB4Bot[ivar.first], ivar.second.slice, false);
      hPaint->SetMinimum(0.89);
      hPaint->SetMaximum(1.02);
      hPaint->Draw();
      Eff_phiMB4Bot[ivar.first]->Draw("samep");
    }
    else  Eff_phiMB4Bot[ivar.first]->Draw("ap");
    
    gPad->Update();
    graph = Eff_phiMB4Bot[ivar.first]->GetPaintedGraph();
    graph->SetMinimum(0.89);
    graph->SetMaximum(1.02);
    
    gPad->Update();
    
    cPhiMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4BotPhiEffVs"+ivar.second.name+".png").c_str());
    
    
    if(ivar.second.name=="Run"){

	Eff_phiMB4Bot[ivar.first] = getIntLumiEff(Eff_phiMB4Bot[ivar.first],ivar.second.slice);
	Eff_phiMB4Bot[ivar.first]->Draw("ap");

	gPad->Update();
	auto graph = Eff_phiMB4Bot[ivar.first]->GetPaintedGraph();
	
	graph->SetMinimum(0.89);
	graph->SetMaximum(1.02);
       
	cPhiMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4BotPhiEffVs"+ivar.second.name+"_intLum.png").c_str());
      }
    delete cPhiMB4Bot ;

  }
  
  //Graph bkg

    for(auto const& ivar : dataCont.var){
      if(!ivar.second.doBkg) continue;
           
    for (int ist=0; ist<4; ist++){
      // find maximum
      float_t Max = 0;
      for (int iwh=0; iwh<5; iwh++){
      	if(Hist_MBWh[ivar.first][iwh][ist]->GetProfMax() > Max) Max = Hist_MBWh[ivar.first][iwh][ist]->GetProfMax();
      }

      //phi
      TCanvas *cMB = new TCanvas(("cMB"+(std::to_string(ist+1))+ivar.second.name).c_str());
      Hist_MBWh[ivar.first][0][ist]->SetMaximum(Max*1.3);
      Hist_MBWh[ivar.first][0][ist]->SetTitle(("MB"+(std::to_string(ist+1))+" bkg vs "+ivar.second.Title+";Int.Lumi. pb^{-1};Rate(Hz/cm^{-2})").c_str());
      
      Hist_MBWh[ivar.first][0][ist]->draw("E1");
      

      TLegend * legMB = new TLegend(0.75,0.75,0.9,0.9);
      for (int iwh=0; iwh<5; iwh++){
	Hist_MBWh[ivar.first][iwh][ist]->draw("sameE1");
	TLegendEntry *le = legMB->AddEntry(Hist_MBWh[ivar.first][iwh][ist],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
	le->SetMarkerColor(Hist_MBWh[ivar.first][iwh][ist]->GetColor());
      }
      legMB->Draw("same");
      cMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB"+(std::to_string(ist+1))+"BkgVs"+ivar.second.name+".png").c_str());
   
      delete  cMB ;
    }//end stations

    //segment
    for (int ist=0; ist<4; ist++){

      //find maximum
      float_t Max = 0;
      for (int iwh=0; iwh<5; iwh++){
      	if(Hist_SegMBWh[ivar.first][iwh][ist]->GetProfMax() > Max) Max = Hist_SegMBWh[ivar.first][iwh][ist]->GetProfMax();
      }

      //phi
      TCanvas *cMB = new TCanvas(("cMB"+(std::to_string(ist+1))+ivar.second.name).c_str());
      Hist_SegMBWh[ivar.first][0][ist]->SetMaximum(Max*1.3);
      Hist_SegMBWh[ivar.first][0][ist]->SetTitle(("MB"+(std::to_string(ist+1))+" segment bkg vs "+ivar.second.Title+";Int.Lumi. pb^{-1};Rate(Hz/cm^{-2})").c_str());
      Hist_SegMBWh[ivar.first][0][ist]->draw("E1");
      
      TLegend * legMB = new TLegend(0.75,0.75,0.9,0.9);
      for (int iwh=0; iwh<5; iwh++){
	Hist_SegMBWh[ivar.first][iwh][ist]->draw("sameE1");
	TLegendEntry *le = legMB->AddEntry(Hist_SegMBWh[ivar.first][iwh][ist],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
	le->SetMarkerColor(Hist_SegMBWh[ivar.first][iwh][ist]->GetColor());
      }
      
      legMB->Draw("same");
      cMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB"+(std::to_string(ist+1))+"SegBkgVs"+ivar.second.name+".png").c_str());
      delete  cMB ;
     }


    // MB4 Top

    // find maximum
    float_t Max = 0;
    for (int iwh=0; iwh<5; iwh++) if(Hist_MB4Top[ivar.first][iwh]->GetProfMax() > Max) Max = Hist_MB4Top[ivar.first][iwh]->GetProfMax();
    
    TCanvas *cMB4Top = new TCanvas(("cMB4Top"+ivar.second.name).c_str());
    Hist_MB4Top[ivar.first][0]->SetMaximum(Max*1.3);
    Hist_MB4Top[ivar.first][0]->SetTitle(("MB4Top bkg vs "+ivar.second.Title+";Int.Lumi. pb^{-1};Rate(Hz/cm^{-2})").c_str());
    Hist_MB4Top[ivar.first][0]->draw("E1");
    
    TLegend * legMB4Top = new TLegend(0.75,0.75,0.9,0.9);

    for (int iwh=0; iwh<5; iwh++){
      Hist_MB4Top[ivar.first][iwh]->draw("sameE1");
      TLegendEntry *le = legMB4Top->AddEntry(Hist_MB4Top[ivar.first][iwh],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
      le->SetMarkerColor(Hist_MB4Top[ivar.first][iwh]->GetColor());  
    }    
    legMB4Top->Draw("same");
    cMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4TopBkgVs"+ivar.second.name+".png").c_str());
    delete cMB4Top ;
    

    //segment

    //find maximum
    Max = 0;
    for (int iwh=0; iwh<5; iwh++) if(Hist_SegMB4Top[ivar.first][iwh]->GetProfMax() > Max) Max = Hist_SegMB4Top[ivar.first][iwh]->GetProfMax();
    
    TCanvas *cSegMB4Top = new TCanvas(("cSegMB4Top"+ivar.second.name).c_str());
    Hist_SegMB4Top[ivar.first][0]->SetMaximum(Max*1.3);
    Hist_SegMB4Top[ivar.first][0]->SetTitle(("MB4Top segment bkg vs "+ivar.second.Title+";Int.Lumi. pb^{-1};Rate(Hz/cm^{-2})").c_str());
    Hist_SegMB4Top[ivar.first][0]->draw("E1");
    
    TLegend * legSegMB4Top = new TLegend(0.75,0.75,0.9,0.9);
    for (int iwh=0; iwh<5; iwh++){
      Hist_SegMB4Top[ivar.first][iwh]->draw("sameE1");
      TLegendEntry *le = legSegMB4Top->AddEntry(Hist_SegMB4Top[ivar.first][iwh],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
      le->SetMarkerColor(Hist_SegMB4Top[ivar.first][iwh]->GetColor());
    }
   
    legSegMB4Top->Draw("same");
    cSegMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4TopSegBkgVs"+ivar.second.name+".png").c_str());
    delete cSegMB4Top ;
    

    // MB4 Bot
    TCanvas *cMB4Bot = new TCanvas(("cMB4Bot"+ivar.second.name).c_str());
    
    //    Hist_MB4Bot[ivar.first]->SetTitle(("MB4Bot bkg vs "+ivar.second.Title+ivar.second.Title+";"+ivar.second.Label+";Rate(Hz/cm^{-2})").c_str());
    Hist_MB4Bot[ivar.first]->SetTitle(("MB4Bot bkg vs "+ivar.second.Title+ivar.second.Title+";Int.Lumi. pb^{-1};Rate(Hz/cm^{-2})").c_str());
    Hist_MB4Bot[ivar.first]->SetMaximum(Max*1.3);
    Hist_MB4Bot[ivar.first]->draw("E1");
    
    cMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4BotBkgVs"+ivar.second.name+".png").c_str());
    delete cMB4Bot;


    // MB4 Bot segment
    TCanvas *cSegMB4Bot = new TCanvas(("cMB4Bot"+ivar.second.name).c_str());
    
    Hist_SegMB4Bot[ivar.first]->SetTitle(("MB4Bot segment bkg vs "+ivar.second.Title+ivar.second.Title+";Int.Lumi. pb^{-1};Rate(Hz/cm^{-2})").c_str());
    Hist_SegMB4Bot[ivar.first]->SetMaximum(Max*1.3);
    Hist_SegMB4Bot[ivar.first]->draw("E1");
    
    cSegMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4BotSegBkgVs"+ivar.second.name+".png").c_str());
    delete cSegMB4Bot;
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


  // xBins.insert(xBins.begin(),arrX->GetArray()[firstBin-1]);
  //Cont.insert(Cont.begin(),0); //del check

  //insert last  value to have the right edge of the last  bin
  xBins.push_back(arr->GetArray()[maxBin]);

  //insert extra value to have one last empty bin
  xBins.push_back(arr->GetArray()[maxBin+1]);

  totCont.push_back(0.);
  passCont.push_back(0.);
  
  string hName = effIn->GetName();
  string hTitle = effIn->GetTitle();
  effIn->Delete();

  TEfficiency * effNew = new TEfficiency(hName.c_str(),hTitle.c_str(),xBins.size()-1,xBins.data());

  for(uint bin = 1; bin<xBins.size(); bin++){
    effNew->SetTotalEvents(bin, totCont[bin-1]);
    effNew->SetPassedEvents(bin, passCont[bin-1]);
  }
  return effNew;
}



void plotter::setPlots(){

     for(auto const& ivar : dataCont.var) {
     for (int iwh=0; iwh<5; iwh++){
       for (int ist=0; ist<4; ist++){
	 if(ivar.second.name=="Run"){

	   if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]  = setEffRun(Eff_phiMBWh[ivar.first][iwh][ist]);
	   if(ivar.second.doEff) EffA_phiMBWh[ivar.first][iwh][ist] = setEffRun(EffA_phiMBWh[ivar.first][iwh][ist]);
	   //	   if(ivar.second.doBkg) Hist_MBWh[ivar.first][iwh][ist]    = set2DHistoBinRun(Hist_MBWh[ivar.first][iwh][ist], ivar.second.slice);

	   if(ivar.second.doBkg) Hist_MBWh[ivar.first][iwh][ist]->getIntLumiHisto();
	   if(ivar.second.doBkg) Hist_SegMBWh[ivar.first][iwh][ist]->getIntLumiHisto();
	 }
	 else{
	   if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]  = setEffBin(Eff_phiMBWh[ivar.first][iwh][ist]);
	   if(ivar.second.doEff) EffA_phiMBWh[ivar.first][iwh][ist] = setEffBin(EffA_phiMBWh[ivar.first][iwh][ist]);
	 }
	 if(ivar.second.name != "Run" &&  ivar.second.name != "Bkg"){
	   // cout<<ivar.second.name<<" "<<Hist_MBWh[ivar.first][iwh][ist]->getName()<<" "<<setprecision(2)<<Hist_MBWh[ivar.first][iwh][ist]->Integral()<<endl;
	   if(ivar.second.doBkg) Hist_MBWh[ivar.first][iwh][ist]->set2DHistoBin();
	   if(ivar.second.doBkg) Hist_SegMBWh[ivar.first][iwh][ist]->set2DHistoBin();
	
	 }
	 if(ist!=3){
	   if(ivar.second.name=="Run"){
	     if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]  = setEffRun(Eff_theMBWh[ivar.first][iwh][ist]);
	     if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh][ist] = setEffRun(EffA_theMBWh[ivar.first][iwh][ist]);	     
	     //if(ivar.second.doBkg) Hist_theMBWh[ivar.first][iwh][ist]->getIntLumiHisto();
	     //if(ivar.second.doBkg) Hist_theMBWh[ivar.first][iwh][ist] = set2DHistoBinRun();
	   }
	   else{
	     if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]  = setEffBin(Eff_theMBWh[ivar.first][iwh][ist]);
	     if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh][ist] = setEffBin(EffA_theMBWh[ivar.first][iwh][ist]);
	   }
	 }
       }
       if(ivar.second.name=="Run"){
	 if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]  = setEffRun(Eff_phiMB4Top[ivar.first][iwh]);
	 if(ivar.second.doEff) EffA_phiMB4Top[ivar.first][iwh] = setEffRun(EffA_phiMB4Top[ivar.first][iwh]);

	 if(ivar.second.doBkg) Hist_MB4Top[ivar.first][iwh]->getIntLumiHisto();
	 if(ivar.second.doBkg) Hist_SegMB4Top[ivar.first][iwh]->getIntLumiHisto();

       }
       else{
	 if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]  = setEffBin(Eff_phiMB4Top[ivar.first][iwh]);
	 if(ivar.second.doEff) EffA_phiMB4Top[ivar.first][iwh] = setEffBin(EffA_phiMB4Top[ivar.first][iwh]);
       }
       if(ivar.second.name != "Run" &&  ivar.second.name != "Bkg"){
	 if(ivar.second.doBkg) Hist_MB4Top[ivar.first][iwh]->set2DHistoBin();
	 if(ivar.second.doBkg) Hist_SegMB4Top[ivar.first][iwh]->set2DHistoBin();
       }
     }
     if(ivar.second.name=="Run"){

       if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]  = setEffRun(Eff_phiMB4Bot[ivar.first]);
       if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first] = setEffRun(EffA_phiMB4Bot[ivar.first]);
      
       if(ivar.second.doBkg) Hist_MB4Bot[ivar.first]->getIntLumiHisto();
       if(ivar.second.doBkg) Hist_SegMB4Bot[ivar.first]->getIntLumiHisto();
     }
     else{
       if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]  = setEffBin(Eff_phiMB4Bot[ivar.first]);
       if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first] = setEffBin(EffA_phiMB4Bot[ivar.first]);
     }
     if(ivar.second.name != "Run" &&  ivar.second.name != "Bkg"){
       if(ivar.second.doBkg) Hist_MB4Bot[ivar.first]->set2DHistoBin();
       if(ivar.second.doBkg) Hist_SegMB4Bot[ivar.first]->set2DHistoBin();
     }
   }
}

//not used so far
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
  
    for(auto & ivar : dataCont.var) {
    if(ivar.second.name=="Run"){

      addBinsSlice(Eff_phiMBWh[ivar.first][0][0]->GetTotalHistogram()->GetXaxis()->GetXbins(), ivar.second.slice);

      for (int iwh=0; iwh<5; iwh++){
	for (int ist=0; ist<4; ist++){
	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]  = addBins(Eff_phiMBWh[ivar.first][iwh][ist],ivar.second.slice);
	  if(ivar.second.doEff) EffA_phiMBWh[ivar.first][iwh][ist] = addBins(EffA_phiMBWh[ivar.first][iwh][ist],ivar.second.slice);	  
	  if(ivar.second.doBkg)	Hist_MBWh[ivar.first][iwh][ist]->addHisto(ivar.second.slice);
	  if(ivar.second.doBkg) Hist_SegMBWh[ivar.first][iwh][ist]->addHisto(ivar.second.slice);

	  if(ist!=3){
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]  = addBins(Eff_theMBWh[ivar.first][iwh][ist],ivar.second.slice);
	    if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh][ist] = addBins(EffA_theMBWh[ivar.first][iwh][ist],ivar.second.slice);
	    // if(ivar.second.doBkg) Hist_theMBWh[ivar.first][iwh][ist]->addHisto(ivar.second.slice);
	  }

	  if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]  = addBins(Eff_phiMB4Top[ivar.first][iwh],ivar.second.slice);
	  if(ivar.second.doEff) EffA_phiMB4Top[ivar.first][iwh] = addBins(EffA_phiMB4Top[ivar.first][iwh],ivar.second.slice);
	  if(ivar.second.doBkg) Hist_MB4Top[ivar.first][iwh]->addHisto(ivar.second.slice);
	  if(ivar.second.doBkg) Hist_SegMB4Top[ivar.first][iwh]->addHisto(ivar.second.slice);
	}

	if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]  = addBins(Eff_phiMB4Bot[ivar.first],ivar.second.slice);
	if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first] = addBins(EffA_phiMB4Bot[ivar.first],ivar.second.slice);
	if(ivar.second.doBkg) Hist_MB4Bot[ivar.first]->addHisto(ivar.second.slice);
	if(ivar.second.doBkg) Hist_SegMB4Bot[ivar.first]->addHisto(ivar.second.slice);
      }    
    }
  } 
}


// Add the run bins to slice 
void plotter::addBinsSlice(const TArrayD * arr, vector<double> & slices){

  std::vector<Double_t> xBins = {}; 
    for(int bin = 1; bin<=arr->GetSize(); bin++) {
      xBins.push_back(arr->GetArray()[bin-1]); 
    }
    
    //Delete last elements that was added just to have an extra item to have the last bin            
    xBins.erase(xBins.end()-1);
   
    //insert in xBins the new runs
    xBins.insert(xBins.end(), slices.begin(), slices.end());
    
    // update slice with old run numbers
    slices= xBins; //check
    isSliceChanged= true; 
}


TEfficiency * plotter::addBins( TEfficiency *effIn, vector<double> slices){ //fix




 /*  std::vector<Double_t> xBins = {}; */
 /*  Int_t nBins  = effIn->GetTotalHistogram()->GetNbinsX(); */

 /*  // run slices has to changed just once */
 /*  if(!isSliceChanged){ */
      
 /*    const  TArrayD *arr = effIn->GetTotalHistogram()->GetXaxis()->GetXbins();  */
 
 /*    for(int bin = 1; bin<=nBins+1; bin++) { */
 /*      xBins.push_back(arr->GetArray()[bin-1]);  */
 /*    } */
    
 /*    //Delete last elements that was added just to have an extra item to have the last bin             */
 /*    xBins.erase(xBins.end()-1); */
   
 /*    //insert in xBins the new runs */
 /*    xBins.insert(xBins.end(), slices.begin(), slices.end()); */
    
 /*    // update slice with old run numbers */
 /*    slices= xBins; //check */
 /*    isSliceChanged= true;  */

 /* } */
 /*  else{ */
 /*    // ivar.second.slice has already been changed and now Xbins can be simply makes equal to  ivar.second.slice */
 /*    xBins  =  slices; */
 /*  } */

 /*  //  TEfficiency * effNew = new TEfficiency(effIn->GetName(),effIn->GetTitle(),xBins.size()-1,xBins.data()); */

  Int_t nBins  = effIn->GetTotalHistogram()->GetNbinsX();

  TEfficiency * effNew = new TEfficiency("",effIn->GetTitle(),slices.size()-1,slices.data());
  string effName = effIn->GetName();
  effNew->SetName(effName.c_str());

  for(int bin = 1; bin<= nBins; bin++){
    effNew->SetTotalEvents( bin, effIn->GetTotalHistogram()->GetBinContent(bin));
    effNew->SetPassedEvents(bin, effIn->GetPassedHistogram()->GetBinContent(bin));
  }
  effIn->Delete();
  return effNew;
}


float plotter::getLumiRun(string run){

  vector<string> files = vector<string>();
  string dir = "data/IntLumi/";

  //get list of file in InLumi directory
  getdir(dir,files);

  if(std::find(files.begin(), files.end(),("data/IntLumi/DT_Run"+run+".csv") ) != files.end()) {
    //push in front of files the file with the name of the run to be faster in case it exists.
    files.insert(files.begin(),"data/IntLumi/DT_Run"+run+".csv"); //check
  }


  Float_t intLumi = 0;
  bool isFound = false; 

  for (unsigned int i = 0;i < files.size();i++) {

    //    cout<<files[i]<<" "<<endl;

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
     if(isFound)  break;     
     else intLumi=0;
  }

  if(!isFound){
    cout<<"ERROR!!!  Run not found in integrated luminosity directory. Please run first createJSONs.py as written in the README"<<endl;
    abort();
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
TH1F * plotter::getPaintHisto( TEfficiency *effIn, vector<double> slices, bool doLumi){

  TH1F * hPaint = new TH1F();
  
  hPaint = (TH1F*)effIn->GetTotalHistogram()->Clone();
  hPaint->Reset();
  hPaint->SetTitle(effIn->GetTitle());

  Int_t nBins = hPaint->GetNbinsX();
  
  if(doLumi){ 
    for(int bin = 1; bin <= nBins; bin++){
      hPaint->GetXaxis()->SetBinLabel(bin, (boost::str(boost::format("%.2f") % lumi::getTotLumiRun(to_string(static_cast<int>(slices[bin-1]))))).c_str());
    }
    hPaint->GetXaxis()->SetTitle("Int.Lumi. pb^{-1}");
    TotLumi = 0; // reset total lumi
  }
  else{
    for(int bin = 1; bin <= nBins; bin++){
      hPaint->GetXaxis()->SetBinLabel(bin,(to_string(static_cast<int>(slices[bin-1]))).c_str()); 
    }
  }

  return hPaint;
}


TEfficiency * plotter::getIntLumiEff( TEfficiency *effIn, vector<double> slices){

  std::vector<Double_t> xBins = {};
  Int_t nBins  = effIn->GetTotalHistogram()->GetNbinsX();

  for(int bin = 1; bin <= nBins; bin++)  xBins.push_back(lumi::getTotLumiRun(to_string(static_cast<int>(slices[bin-1])))); 

  TotLumi = 0; // reset total lumi //del fix

  //To be updated using SetTotalEvents TEfficiency method.

  TH1F * hTot  = new TH1F( (string(effIn->GetName())+"_Tot").c_str(),effIn->GetTitle(),300*xBins.size(),xBins.at(0),xBins.back()+xBins.back()/10.);
  TH1F * hPass = new TH1F( (string(effIn->GetName())+"_Pass").c_str(),effIn->GetTitle(),300*xBins.size(),xBins.at(0),xBins.back()+xBins.back()/10.);

  for(int bin = 1; bin<=nBins; bin++){
    //  cout<<bin<<" "<<hTot->FindBin(xBins[bin-1])<< " "<<hTot->GetBinCenter(bin-1)<<" "<<effIn->GetPassedHistogram()->GetBinContent(bin)/effIn->GetTotalHistogram()->GetBinContent(bin)<<endl;
    hTot->AddBinContent(hTot->FindBin(xBins[bin-1])  , effIn->GetTotalHistogram()->GetBinContent(bin));
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


#endif
