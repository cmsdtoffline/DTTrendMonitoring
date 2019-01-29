#ifndef plotter_h
#define plotter_h

#include "DistTrend.h"
#include "EffTrend.h"
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
#include <TLine.h>
#include <TText.h>
#include <TPaveText.h>
#include <TLatex.h>

struct stat st;

class plotter {
public :
  
  TFile          *fIn;
  TFile          *fOut;

  context dataCont; //struct with the setting
  bool    isSliceChanged = 0; //is it used? fix it or delete it
  string  LumiFileName;

  Float_t TotLumi = 0;

  map<string,vector<vector<EffTrend* > > > Eff_phiMBWh;
  map<string,vector<vector<EffTrend* > > > EffA_phiMBWh;
  
  map<string,vector<vector<EffTrend* > > > Eff_theMBWh;
  map<string,vector<vector<EffTrend* > > > EffA_theMBWh;

  map<string,vector<EffTrend* > >          Eff_phiMB4Top;
  map<string,vector<EffTrend* > >          EffA_phiMB4Top;
  
  map<string,EffTrend* >                   Eff_phiMB4Bot;
  map<string,EffTrend* >                   EffA_phiMB4Bot;
    
  map<string,vector<vector<DistTrend* > > >        Dist_MBWh;
  map<string,vector<vector<DistTrend* > > >        Dist_SegMBWh;

  map<string,vector<DistTrend* > >                 Dist_MB4Top;
  map<string,vector<DistTrend* > >                 Dist_SegMB4Top;

  map<string,DistTrend* >                          Dist_MB4Bot;
  map<string,DistTrend* >                          Dist_SegMB4Bot;
  
  plotter(context extDataCont, std::string inFileName = "", std::string outFileName="", string LumiFileName_ = "",bool doOnlyPlot = kFALSE);

  void   write();
  void   close();
  void   plot(string dateName = "test");
  void   setPlots();
  void   setAddBins();
  int    getdir(string dir, vector<string> &files);   // Take a vector of list and push inside the files inside a directory. Used for lumi per run files.

  string wheelStr(int iwh);

  float  getLumiRun(string Run);    // Get integrated luminosity from the begining of the plot range.
  void addBinsSlice(const TArrayD * arr, vector<double> & slices); 

  ~plotter();
};

plotter::plotter(context extDataCont, std::string inFileName, std::string outFileName, std::string LumiFileName_, bool doOnlyPlot){


  dataCont = extDataCont;  
  if(inFileName!=""){    
    fIn  = new TFile (("data/results/"+dataCont.name+"/"+inFileName+".root").c_str());
    if (!fIn){ cout<<"File In doesn't exist or can't be open"<<endl;    abort();}
   
    if(doOnlyPlot) cout<<("Plotting the stored results of "+inFileName).c_str()<<endl;
    else cout<<"Taking objects from file "<<inFileName+".root"<<endl;
  }
  else  cout<<"Creat new objects and new file "<<(outFileName+".root")<<endl;
  
  fOut  = new TFile (("data/results/"+dataCont.name+"/"+outFileName+".root").c_str(),"RECREATE"); //FIXME. RECREATE->CREATE

  LumiFileName = "DT_"+LumiFileName_;
  for(auto const& ivar : dataCont.var){ 
    Eff_phiMBWh[ivar.first]    = vector<vector<EffTrend*> > ();
    EffA_phiMBWh[ivar.first]   = vector<vector<EffTrend*> > ();
    Eff_theMBWh[ivar.first]    = vector<vector<EffTrend*> > ();
    EffA_theMBWh[ivar.first]   = vector<vector<EffTrend*> > ();
    Eff_phiMB4Top[ivar.first]  = vector<EffTrend*>  ();
    EffA_phiMB4Top[ivar.first] = vector<EffTrend*>  ();

    Dist_MBWh[ivar.first]      = vector<vector<DistTrend*> > ();
    Dist_SegMBWh[ivar.first]   = vector<vector<DistTrend*> > ();

    Dist_MB4Top[ivar.first]    = vector<DistTrend*>  ();
    Dist_SegMB4Top[ivar.first] = vector<DistTrend*>  ();
    
    int nPoints = ivar.second.slice.size();
    for (int iwh=0; iwh<5; iwh++){

      Eff_phiMBWh[ivar.first].push_back(vector<EffTrend*> ());
      EffA_phiMBWh[ivar.first].push_back(vector<EffTrend*> ());
      
      Eff_theMBWh[ivar.first].push_back(vector<EffTrend*> ());
      EffA_theMBWh[ivar.first].push_back(vector<EffTrend*> ());
      
      Dist_MBWh[ivar.first].push_back(vector<DistTrend*> ());
      Dist_SegMBWh[ivar.first].push_back(vector<DistTrend*> ());

      for (int ist=0; ist<4; ist++){
	if(inFileName==""){

	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh].push_back( new EffTrend(("Eff"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
										      "St"+std::to_string(ist)).c_str(),"",nPoints-1,ivar.second.slice.data(),
										      ivar.second.projSlice.size()-1,ivar.second.projSlice.data() ));

	  	  
	  if(ivar.second.doEff) EffA_phiMBWh[ivar.first][iwh].push_back( new EffTrend(("EffA"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
										       "St"+std::to_string(ist)).c_str(),"",nPoints-1,ivar.second.slice.data(),
										      ivar.second.projSlice.size()-1,ivar.second.projSlice.data() ));
	  
	  if(ivar.second.doBkg) Dist_MBWh[ivar.first][iwh].push_back(new DistTrend(("Hist"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
										    "St"+std::to_string(ist)).c_str(),"",
										   //ivar.second.nBins,ivar.second.X0, ivar.second.X1, // equal bins constructor to be fix for increasing option.
										   //dataCont.var["Bkg"].nBins,dataCont.var["Bkg"].X0,dataCont.var["Bkg"].X1)); 
										   ivar.second.slice.size()-1,
										   ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data()));
	  
	  if(ivar.second.doBkg) Dist_SegMBWh[ivar.first][iwh].push_back(new DistTrend(("HistSeg"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
										       "St"+std::to_string(ist)).c_str(),"",		
										      //ivar.second.nBins,ivar.second.X0, ivar.second.X1,						   
										      //dataCont.var["Bkg"].nBins,dataCont.var["Bkg"].X0,dataCont.var["Bkg"].X1)); 
										      ivar.second.slice.size()-1,
										      ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data()));
	  
	}
	else{   

	  //	  wheelStr(int iwh)

	  //cout<<iwh<<" "<<("Eff"+ivar.second.name+"_MBWh"+wheelStr(iwh)+"St"+std::to_string(ist)).c_str()<<endl;  //fixme

	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh].push_back( (EffTrend*)  fIn->Get(("Eff"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
											   "St"+std::to_string(ist)).c_str()));  
	  if(ivar.second.doEff) EffA_phiMBWh[ivar.first][iwh].push_back( (EffTrend*) fIn->Get(("EffA"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
											    "St"+std::to_string(ist)).c_str()));

	  if(ivar.second.doBkg) Dist_MBWh[ivar.first][iwh].push_back( (DistTrend*)     fIn->Get(("Hist"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
											    "St"+std::to_string(ist)).c_str()));
	  if(ivar.second.doBkg) Dist_SegMBWh[ivar.first][iwh].push_back( (DistTrend*)  fIn->Get(("HistSeg"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
											    "St"+std::to_string(ist)).c_str()));	  
	}
	
	if (ist!=3){
	  if(inFileName==""){


	  if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh].push_back( new EffTrend(("EffThe"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
										     "St"+std::to_string(ist)).c_str(),"",nPoints-1,ivar.second.slice.data(),
										      ivar.second.projSlice.size()-1,ivar.second.projSlice.data() ));

	    if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh].push_back( new EffTrend(("EffAThe"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
										      "St"+std::to_string(ist)).c_str(),"",nPoints-1,ivar.second.slice.data(),
										      ivar.second.projSlice.size()-1,ivar.second.projSlice.data() ));
	  }
	  else{   
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh].push_back(  (EffTrend*) fIn->Get(("EffThe"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
													  "St"+std::to_string(ist)).c_str()));
	    if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh].push_back(  (EffTrend*) fIn->Get(("EffAThe"+ivar.second.name+"_MBWh"+wheelStr(iwh)+
											       "St"+std::to_string(ist)).c_str()));
	  }
	}
      }
      if(inFileName==""){

	if(ivar.second.doEff) Eff_phiMB4Top[ivar.first].push_back( new EffTrend(("Eff"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
									     "",nPoints-1,ivar.second.slice.data(),
										ivar.second.projSlice.size()-1,ivar.second.projSlice.data() ));


	if(ivar.second.doEff) EffA_phiMB4Top[ivar.first].push_back( new EffTrend(("EffA"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
									      "",nPoints-1,ivar.second.slice.data(),
										ivar.second.projSlice.size()-1,ivar.second.projSlice.data() ));


	if(ivar.second.doBkg) Dist_MB4Top[ivar.first].push_back( new DistTrend(("Hist"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
					      "",ivar.second.slice.size()-1,ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data()));

	if(ivar.second.doBkg) Dist_SegMB4Top[ivar.first].push_back( new DistTrend(("HistSeg"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str(),
						 "",ivar.second.slice.size()-1,ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data()));   
   }
      
      else{
	if(ivar.second.doEff) Eff_phiMB4Top[ivar.first].push_back(  (EffTrend*) fIn->Get(("Eff"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
	if(ivar.second.doEff) EffA_phiMB4Top[ivar.first].push_back(  (EffTrend*) fIn->Get(("EffA"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
	if(ivar.second.doBkg) Dist_MB4Top[ivar.first].push_back( (DistTrend*) fIn->Get(("Hist"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
	if(ivar.second.doBkg) Dist_SegMB4Top[ivar.first].push_back( (DistTrend*) fIn->Get(("HistSeg"+ivar.second.name+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
      }

    }
    
    if(inFileName==""){

      if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]  =  new EffTrend(("Eff"+ivar.second.name+"_MBBot").c_str(),"",nPoints-1,ivar.second.slice.data(),
										ivar.second.projSlice.size()-1,ivar.second.projSlice.data() );


      if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first] =  new EffTrend(("EffA"+ivar.second.name+"_MBBot").c_str(),"",nPoints-1,ivar.second.slice.data(),
										ivar.second.projSlice.size()-1,ivar.second.projSlice.data() );

      
      if(ivar.second.doBkg) Dist_MB4Bot[ivar.first]    =  new DistTrend(("Hist"+ivar.second.name+"_MBBot").c_str(),"",ivar.second.slice.size()-1,
								ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data());
      
      if(ivar.second.doBkg) Dist_SegMB4Bot[ivar.first] =  new DistTrend(("HistSeg"+ivar.second.name+"_MBBot").c_str(),"",ivar.second.slice.size()-1,
								   ivar.second.slice.data(),dataCont.var["Bkg"].slice.size()-1,dataCont.var["Bkg"].slice.data());
    }
    else{ 
      if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]  =  (EffTrend*) fIn->Get(("Eff"+ivar.second.name+"_MBBot").c_str());
      if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first] =  (EffTrend*) fIn->Get(("EffA"+ivar.second.name+"_MBBot").c_str());
      if(ivar.second.doBkg) Dist_MB4Bot[ivar.first]    =  (DistTrend*) fIn->Get(("Hist"+ivar.second.name+"_MBBot").c_str());
      if(ivar.second.doBkg) Dist_SegMB4Bot[ivar.first] =  (DistTrend*) fIn->Get(("HistSeg"+ivar.second.name+"_MBBot").c_str());
    }
  }

  if(inFileName!=""){
    if(doOnlyPlot && (dataCont.var.find("Run")!=dataCont.var.end()) && (dataCont.var.find("Run"))->second.doEff ){
      const TArrayD * arr = Eff_phiMBWh["Run"][0][0]->GetArrayX();
      std::vector<Double_t> xBins = {}; 
      for(int bin = 1; bin<=arr->GetSize(); bin++) {
	xBins.push_back(arr->GetArray()[bin-1]); 
      }
      dataCont.var["Run"].slice = xBins;
      //uncomment to see the runs
      /* for(int bin = 1; bin<=arr->GetSize(); bin++) { */
      /* 	cout<<dataCont.var["Run"].slice[bin-1]<<endl; */
      /* } */
     }
    else {  setAddBins();}
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
	 if(ivar.second.doBkg) Dist_MBWh[ivar.first][iwh][ist]->Write();
	 if(ivar.second.doBkg) Dist_SegMBWh[ivar.first][iwh][ist]->Write();

	 if(ist!=3){
	   if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->Write();
	   if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh][ist]->Write();
	 }
       }
       if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->Write();
       if(ivar.second.doEff) EffA_phiMB4Top[ivar.first][iwh]->Write();
       if(ivar.second.doBkg) Dist_MB4Top[ivar.first][iwh]->Write();
       if(ivar.second.doBkg) Dist_SegMB4Top[ivar.first][iwh]->Write();
     }
     if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]->Write();
     if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first]->Write();
     if(ivar.second.doBkg) Dist_MB4Bot[ivar.first]->Write();
     if(ivar.second.doBkg) Dist_SegMB4Bot[ivar.first]->Write();
   }
}

void plotter::close(){
   fOut->Close();
  }

void plotter::plot(string dateName){

  //run period lines:

  map<float,string> runPeriods = {};
  
  //runPeriods[]=A "2016-Apr-22";
  //  runPeriods[33.6968206]="Run2016B";
  runPeriods[39.9798206]="2016-Jun-24";  //"Run2016C";
  runPeriods[43.2788206]="2016-Jul-04";  //"Run2016D";
  runPeriods[48.3698206]="2016-Jul-15";  //"Run2016E";
  runPeriods[52.6348206]="2016-Jul-29";  //"Run2016F";  
  runPeriods[56.4198206]="2016-Aug-14";  //"Run2016G";
  runPeriods[64.8798206]="2016-Sep-16";  //Run2016H";
  // runPeriods[75.1350086]="2017-May-19";//"Run2017A";
  runPeriods[76.1190086]="2017-Jun-16";  //"Run2017B";  
  runPeriods[82.6210086]="2017-Jul-18";  //"Run2017C";
  runPeriods[94.9950086]="2017-Aug-30";  //"Run2017D";
  runPeriods[99.9610086]="2017-Sep-20";  //"Run2017E";
  runPeriods[110.429008]="2017-Sep-28";  //"Run2017F"; 
  // runPeriods[126.139008]="Run2017G";   
  // runPeriods[126.494008]="Run2017H";
  runPeriods[127.588278]="2018-Apr-26";  //"Run2018A";   
  runPeriods[143.196278]="2018-May-28";  //"Run2018B"; 
  runPeriods[151.214278]="2018-Jul-07";  //"Run2018C";
  runPeriods[158.597278]="2018-Jul-27";  //"Run2018D";   
  //  runPeriods[194.942278]="Run2018E"; 

  setPlots();

  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(4);
  gStyle->SetCanvasColor(0);

  //  if(stat("plot/",&st) != 0)  system("mkdir plot/"); 
  //  if(stat(("plot/"+dateName).c_str(),&st) != 0)   system(("mkdir plot/"+dateName).c_str());

  if(stat((dataCont.webFolder+"/"+dateName+"/").c_str(),&st) != 0){
    system(("mkdir "+dataCont.webFolder+"/"+dateName).c_str());
    system(("cp "+dataCont.webFolder+"/index.php " +dataCont.webFolder+"/"+dateName).c_str());
  }

  if(stat((dataCont.webFolder+"/"+dateName+"/Efficiency").c_str(),&st) != 0){
  system(("mkdir "+dataCont.webFolder+"/"+dateName+"/Efficiency").c_str());
  system(("cp "+dataCont.webFolder+"/index.php " +dataCont.webFolder+"/"+dateName+"/Efficiency").c_str());
  }
  
  for(auto const& ivar : dataCont.var) {

    //check if some variables has projections and create folder
    if( ivar.second.doProj  && stat((dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin").c_str(),&st) != 0){
      system(("mkdir "+dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin").c_str());
      system(("cp "+dataCont.webFolder+"/index.php " +dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin").c_str());
    }
  }

  if(stat((dataCont.webFolder+"/"+dateName+"/Background").c_str(),&st) != 0){
  system(("mkdir "+dataCont.webFolder+"/"+dateName+"/Background").c_str());
  system(("cp "+dataCont.webFolder+"/index.php " +dataCont.webFolder+"/"+dateName+"/Background").c_str());
  }

  fOut->cd();

  int markerStyle = 20;
  int wwCanvas = dataCont.wwCanvas;
  int whCanvas = dataCont.whCanvas;

  float legx1 =  dataCont.legx1;
  float legy1 =  dataCont.legy1;
  float legx2 =  dataCont.legx2;
  float legy2 =  dataCont.legy2;

  for(auto const& ivar : dataCont.var) {
    for (int iwh=0; iwh<5; iwh++){
      for (int ist=0; ist<4; ist++){
	if(iwh!=4){
	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->SetColor(iwh+1);
	  if(ivar.second.doBkg) Dist_MBWh[ivar.first][iwh][ist]->SetColor(iwh+1);
	  if(ivar.second.doBkg) Dist_SegMBWh[ivar.first][iwh][ist]->SetColor(iwh+1);
	 
	  if(ist!=3){
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->SetColor(iwh+1);
	  }
	}
	else{
	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->SetColor(6);
	  if(ivar.second.doBkg) Dist_MBWh[ivar.first][iwh][ist]->SetColor(6);
	  if(ivar.second.doBkg) Dist_SegMBWh[ivar.first][iwh][ist]->SetColor(6);
	  
	  if(ist!=3){
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->SetColor(6);
	  }
	}

	if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->SetMarkerStyle(markerStyle);
	if(ivar.second.doBkg) Dist_MBWh[ivar.first][iwh][ist]->SetMarkerStyle(markerStyle);
	if(ivar.second.doBkg) Dist_SegMBWh[ivar.first][iwh][ist]->SetMarkerStyle(markerStyle);
	//profile of 2D histograms

	 if(ivar.second.doBkg) Dist_MBWh[ivar.first][iwh][ist]->ProfileX();
	 if(ivar.second.doBkg) Dist_SegMBWh[ivar.first][iwh][ist]->ProfileX();
  	
	if(ist!=3){
	  if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->SetMarkerStyle(markerStyle);	  
	}
      }

      if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->SetMarkerStyle(markerStyle);
      if(ivar.second.doBkg) Dist_MB4Top[ivar.first][iwh]->SetMarkerStyle(markerStyle);
      if(ivar.second.doBkg) Dist_SegMB4Top[ivar.first][iwh]->SetMarkerStyle(markerStyle);
      
      if(iwh!=4){

	if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->SetColor(iwh+1);
	if(ivar.second.doBkg) Dist_MB4Top[ivar.first][iwh]->SetColor(iwh+1);
	if(ivar.second.doBkg) Dist_SegMB4Top[ivar.first][iwh]->SetColor(iwh+1);

      }
      else {
	if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->SetColor(6);
	if(ivar.second.doBkg) Dist_MB4Top[ivar.first][iwh]->SetColor(6);
	if(ivar.second.doBkg) Dist_SegMB4Top[ivar.first][iwh]->SetColor(6);
      }

      if(ivar.second.doBkg) Dist_MB4Top[ivar.first][iwh]->ProfileX();
      if(ivar.second.doBkg) Dist_SegMB4Top[ivar.first][iwh]->ProfileX();
    }

    if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]->SetMarkerStyle(markerStyle);
    if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]->SetColor(1);


    if(ivar.second.doBkg) Dist_MB4Bot[ivar.first]->SetMarkerStyle(markerStyle);
    if(ivar.second.doBkg) Dist_MB4Bot[ivar.first]->SetColor(1);
   
    if(ivar.second.doBkg) Dist_SegMB4Bot[ivar.first]->SetMarkerStyle(markerStyle);
    if(ivar.second.doBkg) Dist_SegMB4Bot[ivar.first]->SetColor(1);
    
    if(ivar.second.doBkg) Dist_MB4Bot[ivar.first]->ProfileX();
    if(ivar.second.doBkg) Dist_SegMB4Bot[ivar.first]->ProfileX();
    
  }

  for(auto const& ivar : dataCont.var){ 
    if(! ivar.second.doEff) continue;
    for (int ist=0; ist<4; ist++){
      for(int bin = -1; bin<(int)ivar.second.projSlice.size(); bin++){ 

	cout<<ivar.second.name<< " "<<bin<<endl;
	if(bin<0 && ivar.second.name!="Run"){
	  continue;
	}
	if(bin==1 && !ivar.second.doProj) break;
	if(bin!=0 && Eff_phiMBWh[ivar.first][0][ist]->checkProj(ivar.second.projSlice,bin,bin)) continue;

	TLegend * legPhiMB = new TLegend(legx1, legy1, legx2,legy2);
					 //0.79,0.72,0.89,0.89);


	TCanvas * cPhiMB = new TCanvas(("cPhiMB"+(std::to_string(ist+1))+ivar.second.name).c_str(),"",wwCanvas,whCanvas);
	
	cPhiMB->SetGrid();
	Eff_phiMBWh[ivar.first][0][ist]->setTitle((";"+ivar.second.Label+"; Hit Efficiency of #Phi Layers").c_str()); 
	
	if(ivar.second.name=="Run"){
	  // Bin -1 amnd 0 ared used to draw the inclusive plots
	  if(bin <= 0)
	    Eff_phiMBWh[ivar.first][0][ist]->drawWithLumi(ivar.second.slice,0,-1,"ap",abs(bin)); 
	  else       Eff_phiMBWh[ivar.first][0][ist]->drawWithLumi(ivar.second.slice,bin,bin,"ap");
	}
	else{
	  if(bin <= 0)  Eff_phiMBWh[ivar.first][0][ist]->draw(0,-1,"ap"); 
	  else 	  Eff_phiMBWh[ivar.first][0][ist]->draw(bin,bin,"ap");			    
	}
	
	for (int iwh=0; iwh<5; iwh++){
	  if(ivar.second.name=="Run"){
	    if(bin == 0) Eff_phiMBWh[ivar.first][iwh][ist]->drawWithLumi(ivar.second.slice,0,-1,"samep");
	    else Eff_phiMBWh[ivar.first][iwh][ist]->drawWithLumi(ivar.second.slice,bin,bin,"samep");
	  }
	  else{
	    if(bin <= 0) Eff_phiMBWh[ivar.first][iwh][ist]->draw(0,-1,"samep");
	    else Eff_phiMBWh[ivar.first][iwh][ist]->draw(bin,bin,"samep");
	  }
	  
	  TLegendEntry *le  = legPhiMB->AddEntry(Eff_phiMBWh[ivar.first][iwh][ist],("MB"+std::to_string(ist+1)+"  Wheel "+std::to_string(iwh-2)).c_str(),"lpe");
	  le->SetMarkerColor(Eff_phiMBWh[ivar.first][iwh][ist]->GetColor() ); 
	  le->SetMarkerStyle(Eff_phiMBWh[ivar.first][iwh][ist]->GetMarkerStyle() ); 

	}

	if(ivar.second.name=="Run"){
	  for( const auto& run_pair : runPeriods ){
	    TLine *lineRun = new TLine(run_pair.first,0.93,run_pair.first,1.0);
	    lineRun->SetLineStyle(2);
	    lineRun->Draw("same");
	    TPaveText *rText = new TPaveText(run_pair.first+0.2,0.94,run_pair.first+2,0.95);
	    rText->SetFillColor(0);
	    TText *text = rText->AddText((run_pair.second).c_str());
	    text->SetTextAngle(90);
	    text->SetTextAlign(22);
	    text->SetTextSize(.03);
	    rText->Draw("same");
	  }
	}

        TPad  *CMSPad  = new TPad("cmspad","cmspad",0.1,0.93,0.9,0.99);
	CMSPad->SetTopMargin (0.0);
	CMSPad->SetRightMargin (0.045);
	CMSPad->SetLeftMargin (0.17);
	CMSPad->SetBottomMargin(0.50);
	CMSPad->Draw();
	CMSPad->cd();

	TLatex *CMSTitle = new TLatex(0.28,0.09,"CMS Preliminary, pp collisions (13 TeV)");
	TLatex *Period   = new TLatex(0.99,0.09,"Run II Data");

	CMSTitle->SetNDC();
	CMSTitle->SetTextColor(1);
	CMSTitle->SetTextSize(0.71);
	CMSTitle->SetTextFont(42);
	CMSTitle->SetTextAlign(32);
	CMSTitle->SetTextAngle(0);
	CMSTitle->Draw();

	Period->SetNDC();
	Period->SetTextColor(1);
	Period->SetTextSize(0.71);
	Period->SetTextFont(42);
	Period->SetTextAlign(32);
	Period->SetTextAngle(0);
	Period->Draw();
	cPhiMB->cd();

	legPhiMB->Draw("same");

	if(bin == 0)      {
	  cPhiMB->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+ivar.second.name+".png").c_str());
	  cPhiMB->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+ivar.second.name+".root").c_str());
	}
	else if(bin==-1 && ivar.second.name=="Run"){
	  cPhiMB->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+ivar.second.name+"_runs.png").c_str());
	  cPhiMB->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+ivar.second.name+"_runs.root").c_str());
	}
	else{
	  cPhiMB->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+ivar.second.name+"_"+ivar.second.projVar+"_"+(boost::str(boost::format("%1%-%2%") % ivar.second.projSlice.at(bin-1) % ivar.second.projSlice.at(bin) )+".png")).c_str() );
	  cPhiMB->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+ivar.second.name+"_"+ivar.second.projVar+"_"+(boost::str(boost::format("%1%-%2%") % ivar.second.projSlice.at(bin-1) % ivar.second.projSlice.at(bin) )+".root")).c_str() );
	}
	delete   cPhiMB;
      }      
    }

    for (int ist=0; ist<4; ist++){      
      //theta
      
	if(ist!=3){
	  
	  for(uint bin = 0; bin<ivar.second.projSlice.size(); bin++){ 
	    
	    if(bin==1 && !ivar.second.doProj) break;
	    if(bin!=0 && Eff_theMBWh[ivar.first][0][ist]->checkProj(ivar.second.slice,bin,bin)) continue;

 	    TLegend * legTheMB = new TLegend(legx1, legy1, legx2,legy2);//new TLegend(0.79,0.31,0.89,0.49); 	   
	    TCanvas * cTheMB = new TCanvas(("cTheMB"+(std::to_string(ist+1))+ivar.second.name).c_str(),"",wwCanvas,whCanvas);
	    
	    cTheMB->SetGrid();	    


	    Eff_theMBWh[ivar.first][0][ist]->setTitle((";"+ivar.second.Label+"; Hit Efficiency of #theta Layers").c_str()); 
	    
	    
	    if(ivar.second.name=="Run"){
	      if(bin == 0)  Eff_theMBWh[ivar.first][0][ist]->drawWithLumi(ivar.second.slice,0,-1,"ap"); 
	      else 	  Eff_theMBWh[ivar.first][0][ist]->drawWithLumi(ivar.second.slice,bin,bin,"ap");
	    }
	    else{
	      if(bin == 0)  Eff_theMBWh[ivar.first][0][ist]->draw(0,-1,"ap"); 
	      else 	  Eff_theMBWh[ivar.first][0][ist]->draw(bin,bin,"ap");			    
	    }
	    
	    for (int iwh=0; iwh<5; iwh++){
	      
	      if(ivar.second.name=="Run"){
		if(bin == 0) Eff_theMBWh[ivar.first][iwh][ist]->drawWithLumi(ivar.second.slice,0,-1,"samep");
		else Eff_theMBWh[ivar.first][iwh][ist]->drawWithLumi(ivar.second.slice,bin,bin,"samep");
	      }
	      else{
		if(bin == 0) Eff_theMBWh[ivar.first][iwh][ist]->draw(0,-1,"samep");
		else Eff_theMBWh[ivar.first][iwh][ist]->draw(bin,bin,"samep");
	      }
	      
	      TLegendEntry *le  = legTheMB->AddEntry(Eff_theMBWh[ivar.first][iwh][ist],("MB"+std::to_string(ist+1)+"  Wheel "+std::to_string(iwh-2)).c_str(),"lpe");
	      le->SetMarkerColor(Eff_theMBWh[ivar.first][iwh][ist]->GetColor() ); 
	      le->SetMarkerStyle(Eff_phiMBWh[ivar.first][iwh][ist]->GetMarkerStyle() ); 
	    }
	    
	    legTheMB->Draw("same");

	    
	    if(ivar.second.name=="Run"){
	      for( const auto& run_pair : runPeriods ){
		TLine *lineRun = new TLine(run_pair.first,0.93,run_pair.first,1.0);
		lineRun->SetLineStyle(2);
		lineRun->Draw("same");
		TPaveText *rText = new TPaveText(run_pair.first+0.2,0.94,run_pair.first+2,0.95);
		rText->SetFillColor(0);
		TText *text = rText->AddText((run_pair.second).c_str());
		text->SetTextAngle(90);
		text->SetTextAlign(22);
		text->SetTextSize(.03);
		rText->Draw("same");
	      }
	    }

        TPad  *CMSPad  = new TPad("cmspad","cmspad",0.1,0.93,0.9,0.99);

	CMSPad->SetTopMargin (0.0);
	CMSPad->SetRightMargin (0.045);
	CMSPad->SetLeftMargin (0.17);
	CMSPad->SetBottomMargin(0.50);
	CMSPad->Draw();
	CMSPad->cd();

	TLatex *CMSTitle = new TLatex(0.28,0.09,"CMS Preliminary, pp collisions (13 TeV)");
	TLatex *Period   = new TLatex(0.99,0.09,"Run II Data");

	CMSTitle->SetNDC();
	CMSTitle->SetTextColor(1);
	CMSTitle->SetTextSize(0.71);
	CMSTitle->SetTextFont(42);
	CMSTitle->SetTextAlign(32);
	CMSTitle->SetTextAngle(0);
	CMSTitle->Draw();

	Period->SetNDC();
	Period->SetTextColor(1);
	Period->SetTextSize(0.71);
	Period->SetTextFont(42);
	Period->SetTextAlign(32);
	Period->SetTextAngle(0);
	Period->Draw();
	cTheMB->cd();

	    
	    
	    if(bin == 0){
	      cTheMB->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"TheEffVs"+ivar.second.name+".png").c_str());
	      cTheMB->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"TheEffVs"+ivar.second.name+".root").c_str());
	    }
	    else{
	      cTheMB->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin/"+"MB"+(std::to_string(ist+1))+"TheEffVs"+ivar.second.name+"_"+ivar.second.projVar+"_"+(boost::str(boost::format("%1%-%2%") % ivar.second.projSlice.at(bin-1) % ivar.second.projSlice.at(bin) )+".png")).c_str());
	      cTheMB->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin/"+"MB"+(std::to_string(ist+1))+"TheEffVs"+ivar.second.name+"_"+ivar.second.projVar+"_"+(boost::str(boost::format("%1%-%2%") % ivar.second.projSlice.at(bin-1) % ivar.second.projSlice.at(bin) )+".root")).c_str());
	    }
	    delete  cTheMB;
	  }
	}
      }
      //phi MB4 Top
      
      for(uint bin = 0; bin<ivar.second.projSlice.size(); bin++){ 
	
	if(bin==1 && !ivar.second.doProj) break;
	if(bin!=0 && Eff_phiMB4Top[ivar.first][0]->checkProj(ivar.second.slice,bin,bin)) continue;	  
	
	TLegend * legPhiMB4Top = new TLegend(legx1, legy1, legx2,legy2);//new TLegend(0.79,0.72,0.89,0.89); /fixme

	TCanvas * cPhiMB4Top = new TCanvas(("cPhiMB4Top"+ivar.second.name).c_str(),"",wwCanvas,whCanvas);
	cPhiMB4Top->SetGrid();	

	Eff_phiMB4Top[ivar.first][0]->setTitle((";"+ivar.second.Label+"; Hit Efficiency of #Phi Layers").c_str()); 
	//Eff_phiMB4Top[ivar.first][0]->setTitle(("MB4Top  eff vs "+ivar.second.Title+";"+ivar.second.Label+";Eff").c_str()); 
	
	if(ivar.second.name=="Run"){
	  if(bin == 0)  Eff_phiMB4Top[ivar.first][0]->drawWithLumi(ivar.second.slice,0,-1,"ap"); 
	  else 	  Eff_phiMB4Top[ivar.first][0]->drawWithLumi(ivar.second.slice,bin,bin,"ap");
	}
	else{
	  if(bin == 0)  Eff_phiMB4Top[ivar.first][0]->draw(0,-1,"ap"); 
	  else 	  Eff_phiMB4Top[ivar.first][0]->draw(bin,bin,"ap");			    
	}
	
	for (int iwh=0; iwh<5; iwh++){
	  
	  if(ivar.second.name=="Run"){
	    if(bin == 0) Eff_phiMB4Top[ivar.first][iwh]->drawWithLumi(ivar.second.slice,0,-1,"samep");
	    else Eff_phiMB4Top[ivar.first][iwh]->drawWithLumi(ivar.second.slice,bin,bin,"samep");
	  }
	  else{
	    if(bin == 0) Eff_phiMB4Top[ivar.first][iwh]->draw(0,-1,"samep");
	    else Eff_phiMB4Top[ivar.first][iwh]->draw(bin,bin,"samep");
	  }
	  
	  TLegendEntry *le  = legPhiMB4Top->AddEntry(Eff_phiMB4Top[ivar.first][iwh],("MB4Top  Wheel "+std::to_string(iwh-2)).c_str(),"lpe");
	  le->SetMarkerColor(Eff_phiMB4Top[ivar.first][iwh]->GetColor() ); 
	  le->SetMarkerStyle(Eff_phiMB4Top[ivar.first][iwh]->GetMarkerStyle() ); 
	}
	
	legPhiMB4Top->Draw("same");
	

	if(ivar.second.name=="Run"){
	  for( const auto& run_pair : runPeriods ){
	    TLine *lineRun = new TLine(run_pair.first,0.93,run_pair.first,1.0);
	    lineRun->SetLineStyle(2);
	    lineRun->Draw("same");
	    TPaveText *rText = new TPaveText(run_pair.first+0.2,0.94,run_pair.first+2,0.95);
	    rText->SetFillColor(0);
	    TText *text = rText->AddText((run_pair.second).c_str());
	    text->SetTextAngle(90);
	    text->SetTextAlign(22);
	    text->SetTextSize(.03);
	    rText->Draw("same");
	  }
	}

        TPad  *CMSPad  = new TPad("cmspad","cmspad",0.1,0.93,0.9,0.99);

	CMSPad->SetTopMargin (0.0);
	CMSPad->SetRightMargin (0.045);
	CMSPad->SetLeftMargin (0.17);
	CMSPad->SetBottomMargin(0.50);
	CMSPad->Draw();
	CMSPad->cd();

	TLatex *CMSTitle = new TLatex(0.28,0.09,"CMS Preliminary, pp collisions (13 TeV)");
	TLatex *Period   = new TLatex(0.99,0.09,"Run II Data");

	CMSTitle->SetNDC();
	CMSTitle->SetTextColor(1);
	CMSTitle->SetTextSize(0.71);
	CMSTitle->SetTextFont(42);
	CMSTitle->SetTextAlign(32);
	CMSTitle->SetTextAngle(0);
	CMSTitle->Draw();

	Period->SetNDC();
	Period->SetTextColor(1);
	Period->SetTextSize(0.71);
	Period->SetTextFont(42);
	Period->SetTextAlign(32);
	Period->SetTextAngle(0);
	Period->Draw();

	cPhiMB4Top->cd();

	   

	if(bin == 0) {
	  cPhiMB4Top->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4Top"+"PhiEffVs"+ivar.second.name+".png").c_str());
	  cPhiMB4Top->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4Top"+"PhiEffVs"+ivar.second.name+".root").c_str());
	}
	else {
	  cPhiMB4Top->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin/"+"MB4Top"+"PhiEffVs"+ivar.second.name+"_"+ivar.second.projVar+"_"+(boost::str(boost::format("%1%-%2%") % ivar.second.projSlice.at(bin-1) % ivar.second.projSlice.at(bin) )+".png")).c_str());
	  cPhiMB4Top->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin/"+"MB4Top"+"PhiEffVs"+ivar.second.name+"_"+ivar.second.projVar+"_"+(boost::str(boost::format("%1%-%2%") % ivar.second.projSlice.at(bin-1) % ivar.second.projSlice.at(bin) )+".root")).c_str());
	}
	delete   cPhiMB4Top;
      }
      
      //phi MB4 Bot

      for(uint bin = 0; bin<ivar.second.projSlice.size(); bin++){ 
	
	if(bin==1 && !ivar.second.doProj) break;
	if(bin!=0 && Eff_phiMB4Bot[ivar.first]->checkProj(ivar.second.slice,bin,bin)) continue;
	
	TCanvas * cPhiMB4Bot   = new TCanvas(("cPhiMB4Bot"+ivar.second.name).c_str(),"",wwCanvas,whCanvas);

	TLegend * legPhiMB4Bot = new TLegend(legx1, legy1, legx2,legy2);//new TLegend(0.79,0.72,0.89,0.89); 

	cPhiMB4Bot->SetGrid();	

	Eff_phiMB4Bot[ivar.first]->setTitle((";"+ivar.second.Label+"; Hit Efficiency of #Phi Layer").c_str()); 
	//Eff_phiMB4Bot[ivar.first]->setTitle(("MB4Bot  eff vs "+ivar.second.Title+";"+ivar.second.Label+";Eff").c_str()); 
	
	if(ivar.second.name=="Run"){
	  if(bin == 0)  Eff_phiMB4Bot[ivar.first]->drawWithLumi(ivar.second.slice,0,-1,"ap"); 
	  else 	  Eff_phiMB4Bot[ivar.first]->drawWithLumi(ivar.second.slice,bin,bin,"ap");
	}
	else{
	  if(bin == 0)  Eff_phiMB4Bot[ivar.first]->draw(0,-1,"ap"); 
	  else 	  Eff_phiMB4Bot[ivar.first]->draw(bin,bin,"ap");			    
	}

	TLegendEntry *le  = legPhiMB4Bot->AddEntry(Eff_phiMB4Bot[ivar.first],"MB4Bot","lpe");
	le->SetMarkerColor(Eff_phiMB4Bot[ivar.first]->GetColor() ); 
	le->SetMarkerStyle(Eff_phiMB4Bot[ivar.first]->GetMarkerStyle() ); 	
	legPhiMB4Bot->Draw("same");

  
	if(ivar.second.name=="Run"){
	  for( const auto& run_pair : runPeriods ){
	    TLine *lineRun = new TLine(run_pair.first,0.93,run_pair.first,1.0);
	    lineRun->SetLineStyle(2);
	    lineRun->Draw("same");
	    TPaveText *rText = new TPaveText(run_pair.first+0.2,0.94,run_pair.first+2,0.95);
	    rText->SetFillColor(0);
	    TText *text = rText->AddText((run_pair.second).c_str());
	    text->SetTextAngle(90);
	    text->SetTextAlign(22);
	    text->SetTextSize(.03);
	    rText->Draw("same");
	  }
	}


        TPad  *CMSPad  = new TPad("cmspad","cmspad",0.1,0.93,0.9,0.99);

	CMSPad->SetTopMargin (0.0);
	CMSPad->SetRightMargin (0.045);
	CMSPad->SetLeftMargin (0.17);
	CMSPad->SetBottomMargin(0.50);
	CMSPad->Draw();
	CMSPad->cd();

	TLatex *CMSTitle = new TLatex(0.28,0.09,"CMS Preliminary, pp collisions (13 TeV)");
	TLatex *Period   = new TLatex(0.99,0.09,"Run II Data");

	CMSTitle->SetNDC();
	CMSTitle->SetTextColor(1);
	CMSTitle->SetTextSize(0.71);
	CMSTitle->SetTextFont(42);
	CMSTitle->SetTextAlign(32);
	CMSTitle->SetTextAngle(0);
	CMSTitle->Draw();

	Period->SetNDC();
	Period->SetTextColor(1);
	Period->SetTextSize(0.71);
	Period->SetTextFont(42);
	Period->SetTextAlign(32);
	Period->SetTextAngle(0);
	Period->Draw();

	cPhiMB4Bot->cd();

	if(bin == 0){
	  cPhiMB4Bot->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4Bot"+"PhiEffVs"+ivar.second.name+".png").c_str());
	  cPhiMB4Bot->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+"MB4Bot"+"PhiEffVs"+ivar.second.name+".root").c_str());
	}
	else{
	  cPhiMB4Bot->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin/"+"MB4Bot"+"PhiEffVs"+ivar.second.name+"_"+ivar.second.projVar+"_"+(boost::str(boost::format("%1%-%2%") % ivar.second.projSlice.at(bin-1) % ivar.second.projSlice.at(bin) )+".png")).c_str());
	  cPhiMB4Bot->SaveAs( (dataCont.webFolder+"/"+dateName+"/Efficiency/"+ivar.second.projVar+"Bin/"+"MB4Bot"+"PhiEffVs"+ivar.second.name+"_"+ivar.second.projVar+"_"+(boost::str(boost::format("%1%-%2%") % ivar.second.projSlice.at(bin-1) % ivar.second.projSlice.at(bin) )+".root")).c_str());
	}
	delete   cPhiMB4Bot;
      }
    }
  
  //Graph bkg

    for(auto const& ivar : dataCont.var){
      if(!ivar.second.doBkg) continue;
           
    for (int ist=0; ist<4; ist++){
      // find maximum
      float_t Max = 0;
      for (int iwh=0; iwh<5; iwh++){
      	if(Dist_MBWh[ivar.first][iwh][ist]->GetProfMax() > Max) Max = Dist_MBWh[ivar.first][iwh][ist]->GetProfMax();
      }

      //phi
      TCanvas *cMB = new TCanvas(("cMB"+(std::to_string(ist+1))+ivar.second.name).c_str(),"",wwCanvas,whCanvas);
      cMB->SetGrid();

      Dist_MBWh[ivar.first][0][ist]->SetMaximum(Max*1.45);
      Dist_MBWh[ivar.first][0][ist]->setTitle(("MB"+(std::to_string(ist+1))+"bkg vs "+ivar.second.Title+";"+ivar.second.Label+";Rate(Hz/cm^{2})").c_str());
      Dist_MBWh[ivar.first][0][ist]->draw("E1");
      

      TLegend * legMB = new TLegend(legx1, legy1, legx2,legy2);//new TLegend(0.75,0.75,0.9,0.9);
      for (int iwh=0; iwh<5; iwh++){
	Dist_MBWh[ivar.first][iwh][ist]->draw("sameE1");
	TLegendEntry *le = legMB->AddEntry(Dist_MBWh[ivar.first][iwh][ist],("MB"+std::to_string(ist+1)+"Wh"+std::to_string(iwh-2)).c_str(),"lpe");
	le->SetMarkerColor(Dist_MBWh[ivar.first][iwh][ist]->GetColor());
      }

      legMB->Draw("same");
      cMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB"+(std::to_string(ist+1))+"BkgVs"+ivar.second.name+".png").c_str());
      cMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB"+(std::to_string(ist+1))+"BkgVs"+ivar.second.name+".root").c_str());
   
      delete  cMB ;
    }//end stations

    //segment
    for (int ist=0; ist<4; ist++){

      //find maximum
      float_t Max = 0;
      for (int iwh=0; iwh<5; iwh++){
      	if(Dist_SegMBWh[ivar.first][iwh][ist]->GetProfMax() > Max) Max = Dist_SegMBWh[ivar.first][iwh][ist]->GetProfMax();
      }

      //phi
      TCanvas *cMB = new TCanvas(("cMB"+(std::to_string(ist+1))+ivar.second.name).c_str(),"",wwCanvas,whCanvas);
      cMB->SetGrid();
      Dist_SegMBWh[ivar.first][0][ist]->SetMaximum(Max*1.45);
      Dist_SegMBWh[ivar.first][0][ist]->setTitle(("MB"+(std::to_string(ist+1))+" segment bkg vs "+ivar.second.Title+";"+ivar.second.Label+";Rate(Hz/cm^{2})").c_str());
      Dist_SegMBWh[ivar.first][0][ist]->draw("E1");
      
      TLegend * legMB = new TLegend(legx1, legy1, legx2,legy2);//new TLegend(0.75,0.75,0.9,0.9);
      for (int iwh=0; iwh<5; iwh++){
	Dist_SegMBWh[ivar.first][iwh][ist]->draw("sameE1");
	TLegendEntry *le = legMB->AddEntry(Dist_SegMBWh[ivar.first][iwh][ist],("MB"+std::to_string(ist+1)+"Wh"+std::to_string(iwh-2)).c_str(),"lpe");
	le->SetMarkerColor(Dist_SegMBWh[ivar.first][iwh][ist]->GetColor());
      }
      
      legMB->Draw("same");
      cMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB"+(std::to_string(ist+1))+"SegBkgVs"+ivar.second.name+".png").c_str());
      cMB->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB"+(std::to_string(ist+1))+"SegBkgVs"+ivar.second.name+".root").c_str());
      delete  cMB ;
     }


    // MB4 Top

    // find maximum
    float_t Max = 0;
    for (int iwh=0; iwh<5; iwh++) if(Dist_MB4Top[ivar.first][iwh]->GetProfMax() > Max) Max = Dist_MB4Top[ivar.first][iwh]->GetProfMax();
    
    TCanvas *cMB4Top = new TCanvas(("cMB4Top"+ivar.second.name).c_str(),"",wwCanvas,whCanvas);
    cMB4Top->SetGrid();
    Dist_MB4Top[ivar.first][0]->SetMaximum(Max*1.4);
    Dist_MB4Top[ivar.first][0]->setTitle(("MBTop segment bkg vs "+ivar.second.Title+";"+ivar.second.Label+";Rate(Hz/cm^{2})").c_str());
    Dist_MB4Top[ivar.first][0]->draw("E1");
    
    TLegend * legMB4Top = new TLegend(legx1, legy1, legx2,legy2);//new TLegend(0.75,0.75,0.9,0.9);

    for (int iwh=0; iwh<5; iwh++){
      Dist_MB4Top[ivar.first][iwh]->draw("sameE1");
      TLegendEntry *le = legMB4Top->AddEntry(Dist_MB4Top[ivar.first][iwh],("MB4Top Wh"+std::to_string(iwh-2)).c_str(),"lpe");
      le->SetMarkerColor(Dist_MB4Top[ivar.first][iwh]->GetColor());  
    }    
    legMB4Top->Draw("same");
    cMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4TopBkgVs"+ivar.second.name+".png").c_str());
    cMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4TopBkgVs"+ivar.second.name+".root").c_str());
    delete cMB4Top ;
    

    //segment

    //find maximum
    Max = 0;
    for (int iwh=0; iwh<5; iwh++) if(Dist_SegMB4Top[ivar.first][iwh]->GetProfMax() > Max) Max = Dist_SegMB4Top[ivar.first][iwh]->GetProfMax();
    
    TCanvas *cSegMB4Top = new TCanvas(("cSegMB4Top"+ivar.second.name).c_str(),"",wwCanvas,whCanvas);
    cSegMB4Top->SetGrid();

    Dist_SegMB4Top[ivar.first][0]->SetMaximum(Max*1.45);
    Dist_SegMB4Top[ivar.first][0]->setTitle(("MB4Top segment bkg vs "+ivar.second.Title+";"+ivar.second.Label+";Rate(Hz/cm^{2})").c_str());
    Dist_SegMB4Top[ivar.first][0]->draw("E1");
    
    TLegend * legSegMB4Top = new TLegend(legx1, legy1, legx2,legy2);//new TLegend(0.75,0.75,0.9,0.9);
    for (int iwh=0; iwh<5; iwh++){
      Dist_SegMB4Top[ivar.first][iwh]->draw("sameE1");
      TLegendEntry *le = legSegMB4Top->AddEntry(Dist_SegMB4Top[ivar.first][iwh],("MB4Top Wh"+std::to_string(iwh-2)).c_str(),"lpe");
      le->SetMarkerColor(Dist_SegMB4Top[ivar.first][iwh]->GetColor());
    }
   
    legSegMB4Top->Draw("same");
    cSegMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4TopSegBkgVs"+ivar.second.name+".png").c_str());
    cSegMB4Top->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4TopSegBkgVs"+ivar.second.name+".root").c_str());
    delete cSegMB4Top ;
    
    // MB4 Bot
    TCanvas *cMB4Bot = new TCanvas(("cMB4Bot"+ivar.second.name).c_str(),"",wwCanvas,whCanvas);
    cMB4Bot->SetGrid();

    Dist_MB4Bot[ivar.first]->setTitle(("MB4Bot bkg vs "+ivar.second.Title+";"+ivar.second.Label+";Rate(Hz/cm^{2})").c_str());
    Dist_MB4Bot[ivar.first]->SetMaximum( Dist_MB4Bot[ivar.first]->GetProfMax()*1.45);
    Dist_MB4Bot[ivar.first]->draw("E1");
    
    cMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4BotBkgVs"+ivar.second.name+".png").c_str());
    cMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4BotBkgVs"+ivar.second.name+".root").c_str());
    delete cMB4Bot;

    // MB4 Bot segment
    TCanvas *cSegMB4Bot = new TCanvas(("cMB4Bot"+ivar.second.name).c_str(),"",wwCanvas,whCanvas);
    cSegMB4Bot->SetGrid();    

    Dist_SegMB4Bot[ivar.first]->setTitle(("MB4Bot segment bkg vs "+ivar.second.Title+";"+ivar.second.Label+";Rate(Hz/cm^{2})").c_str());
    Dist_SegMB4Bot[ivar.first]->SetMaximum( Dist_SegMB4Bot[ivar.first]->GetProfMax()*1.45);
    Dist_SegMB4Bot[ivar.first]->draw("E1");
    
    cSegMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4BotSegBkgVs"+ivar.second.name+".png").c_str());
    cSegMB4Bot->SaveAs((dataCont.webFolder+"/"+dateName+"/Background/"+"MB4BotSegBkgVs"+ivar.second.name+".root").c_str());
    delete cSegMB4Bot;
  }
}

void plotter::setPlots(){

     for(auto const& ivar : dataCont.var) {
     for (int iwh=0; iwh<5; iwh++){
       for (int ist=0; ist<4; ist++){
	 if(ivar.second.name=="Run"){
	   if(ivar.second.doBkg) Dist_MBWh[ivar.first][iwh][ist]->getIntLumiHisto();
	   if(ivar.second.doBkg) Dist_SegMBWh[ivar.first][iwh][ist]->getIntLumiHisto();
	 }
	 else{
	   if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->setEffBin();
	   if(ivar.second.doEff) EffA_phiMBWh[ivar.first][iwh][ist]->setEffBin();
	 }
	 if(ivar.second.name != "Run" &&  ivar.second.name != "Bkg"){ 
	   if(ivar.second.doBkg) Dist_MBWh[ivar.first][iwh][ist]->set2DHistoBin();
	   if(ivar.second.doBkg) Dist_SegMBWh[ivar.first][iwh][ist]->set2DHistoBin();
	 }
	 if(ist!=3){
	   if(ivar.second.name!="Run"){
	     if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->setEffBin();
	     if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh][ist]->setEffBin();
	   }
	 }
       }
       if(ivar.second.name=="Run"){
	 if(ivar.second.doBkg) Dist_MB4Top[ivar.first][iwh]->getIntLumiHisto();
	 if(ivar.second.doBkg) Dist_SegMB4Top[ivar.first][iwh]->getIntLumiHisto();
       }
       else{
	 if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->setEffBin();
	 if(ivar.second.doEff) EffA_phiMB4Top[ivar.first][iwh]->setEffBin();
       }
       if(ivar.second.name != "Run" &&  ivar.second.name != "Bkg"){
	 if(ivar.second.doBkg) Dist_MB4Top[ivar.first][iwh]->set2DHistoBin();
	 if(ivar.second.doBkg) Dist_SegMB4Top[ivar.first][iwh]->set2DHistoBin();
       }
     }
     if(ivar.second.name=="Run"){
       if(ivar.second.doBkg) Dist_MB4Bot[ivar.first]->getIntLumiHisto();
       if(ivar.second.doBkg) Dist_SegMB4Bot[ivar.first]->getIntLumiHisto();
     }
     else{
       if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]->setEffBin();
       if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first]->setEffBin();
     }
     if(ivar.second.name != "Run" &&  ivar.second.name != "Bkg"){
       if(ivar.second.doBkg) Dist_MB4Bot[ivar.first]->set2DHistoBin();
       if(ivar.second.doBkg) Dist_SegMB4Bot[ivar.first]->set2DHistoBin();
     }
   }

}


void plotter::setAddBins(){
    for(auto & ivar : dataCont.var) {
    if(ivar.second.name=="Run" && ivar.second.doEff )  {

      addBinsSlice(Eff_phiMBWh[ivar.first][0][0]->GetArrayX(), ivar.second.slice);

      for (int iwh=0; iwh<5; iwh++){
	for (int ist=0; ist<4; ist++){

	  if(ivar.second.doEff) Eff_phiMBWh[ivar.first][iwh][ist]->addBins(ivar.second.slice);
	  if(ivar.second.doEff) EffA_phiMBWh[ivar.first][iwh][ist]->addBins(ivar.second.slice);	  
	  if(ivar.second.doBkg)	Dist_MBWh[ivar.first][iwh][ist]->addHisto(ivar.second.slice);
	  if(ivar.second.doBkg) Dist_SegMBWh[ivar.first][iwh][ist]->addHisto(ivar.second.slice);

	  if(ist!=3){
	    if(ivar.second.doEff) Eff_theMBWh[ivar.first][iwh][ist]->addBins(ivar.second.slice);
	    if(ivar.second.doEff) EffA_theMBWh[ivar.first][iwh][ist]->addBins(ivar.second.slice);
	  }
	}
	  if(ivar.second.doEff) Eff_phiMB4Top[ivar.first][iwh]->addBins(ivar.second.slice);
	  if(ivar.second.doEff) EffA_phiMB4Top[ivar.first][iwh]->addBins(ivar.second.slice);
	  if(ivar.second.doBkg) Dist_MB4Top[ivar.first][iwh]->addHisto(ivar.second.slice);
	  if(ivar.second.doBkg) Dist_SegMB4Top[ivar.first][iwh]->addHisto(ivar.second.slice);
      }	

	if(ivar.second.doEff) Eff_phiMB4Bot[ivar.first]->addBins(ivar.second.slice);
	if(ivar.second.doEff) EffA_phiMB4Bot[ivar.first]->addBins(ivar.second.slice);
	if(ivar.second.doBkg) Dist_MB4Bot[ivar.first]->addHisto(ivar.second.slice);
	if(ivar.second.doBkg) Dist_SegMB4Bot[ivar.first]->addHisto(ivar.second.slice);
         
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
  slices = xBins; //check
  isSliceChanged= true; 
}



float plotter::getLumiRun(string run){

  vector<string> files = vector<string>();
  string dir = "data/IntLumi/";

  //get list of file in InLumi directory
  getdir(dir,files);

  if(std::find(files.begin(), files.end(),("data/IntLumi/DT_Run"+run+".csv") ) != files.end()) {
    //push in front of files the file with the name of the run to be faster in case it exists.
    files.insert(files.begin(),"data/IntLumi/DT_Run"+run+".csv");
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

string plotter::wheelStr(int iwh){
  if(iwh-2<0)   return "M"+std::to_string(abs(iwh-2));
  else return std::to_string(iwh-2);
}


#endif
