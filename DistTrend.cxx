#include <iostream>
#include "Hist.h"
#include <vector>
#include <utility>
#include <TROOT.h> 
#include <TH2.h>
#include <TProfile.h>
#include <TColor.h>
#include "Utilities.h"

ClassImp(Hist)

Hist::Hist(const char* name,const char* title,Int_t nbinsx,
	     const Double_t* xbins,Int_t nbinsy,const Double_t* ybins) :TNamed(name,title){
  theName = name;
  fHistogram = new TH2F((theName).c_str(),title,nbinsx,xbins,nbinsy,ybins);
  fProf = new TProfile();
}


Hist::Hist(const char* name,const char* title,Int_t nbinsx,
	   Float_t X0, Float_t X1, Int_t nbinsy,Float_t Y0, Float_t Y1) :TNamed(name,title){
  theName = name;
  fHistogram = new TH2F((theName).c_str(),title,nbinsx,X0,X1,nbinsy,Y0,Y1);
  fProf = new TProfile();
}



Hist::Hist() :TNamed("","")
{
  theName = "";
  fHistogram = new TH2F();
  fProf = new TProfile();
}

Hist::~Hist() 
{
}

void Hist::Fill(Double_t x, Double_t y){ 
  fHistogram->Fill(x,y);
}

float Hist::GetProfMax(){
  return fProf->GetMaximum();
}


void Hist::SetMaximum(float Max){fProf->SetMaximum(Max);}


void Hist::ProfileX(){
  fProf = fHistogram->ProfileX((theName+"_prx").c_str(),1,-1);
}

void Hist::draw(string  option){
  fProf->Draw(option.c_str());
}

void Hist::drawHisto(string  option){
  fHistogram->Draw(option.c_str());
}


void Hist::addHisto(vector<double> slices){
  
  Int_t nXBins  = fHistogram->GetNbinsX();       
  Int_t nYBins  = fHistogram->GetNbinsY();       
  
  const  TArrayD *arrY = fHistogram->GetYaxis()->GetXbins();
  
  TH2F * hNew = new TH2F("","",slices.size()-1,slices.data(),nYBins,arrY->GetArray()); 

  for(int binx = 1; binx<=nXBins; binx++){  
    for(int biny = 1; biny<=nYBins; biny++){  
      hNew->SetBinContent(binx,biny, fHistogram->GetBinContent(binx,biny));
    }
  }
  fHistogram = hNew;
  fHistogram->SetName(theName.c_str());
}


TH2F * Hist::GetHisto(){ 
  return fHistogram;
}

void Hist::setTitle(const char *title){
  fProf->SetTitle(title);
}


// To be used for equal bins constructor. To be fix for increasing range case.
vector<double> Hist::GetHistoArrayX() {
  

  if ( fHistogram->GetXaxis()->GetXbins()->GetSize()!=0){
  
    const    double * arr = fHistogram->GetXaxis()->GetXbins()->GetArray();
    std::vector<double> vec;
    vec.assign(arr, arr + fHistogram->GetNbinsX()+1);
    return vec;
  }
  else{
    
    Int_t nXBins  = fHistogram->GetNbinsX();
    float X0      = fHistogram->GetXaxis()->GetXmin();
    float X1      = fHistogram->GetXaxis()->GetXmax();
    
    double inc = (X1-X0)/nXBins;
    vector<double> vec;
    for(float n = X0; n <= X1; n += inc)  vec.push_back(n);
    return vec;
  }
}


vector<double> Hist::GetHistoArrayY() {
  

  if ( fHistogram->GetYaxis()->GetXbins()->GetSize()!=0){
  
    const    double * arr = fHistogram->GetYaxis()->GetXbins()->GetArray();
    std::vector<double> vec;
    vec.assign(arr, arr + fHistogram->GetNbinsY()+1);
    return vec;
  }
  else{
    
    Int_t nYBins  = fHistogram->GetNbinsY();
    float Y0      = fHistogram->GetYaxis()->GetXmin();
    float Y1      = fHistogram->GetYaxis()->GetXmax();
    
    double inc = (Y1-Y0)/nYBins;
    vector<double> vec;
    for(float n = Y0; n <= Y1; n += inc)  vec.push_back(n);
    return vec;
  }
}



void Hist::set2DHistoBin(Float_t MaxErr){
  
  Int_t nXBins  = fHistogram->GetNbinsX();
  Int_t nYBins  = fHistogram->GetNbinsY();

  vector<double> arrX =  GetHistoArrayX();
  vector<double> arrY =  GetHistoArrayY();

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
    if(!foundMax && fHistogram->Integral((nXBins - bin),(nXBins - bin),0,-1 ) >0.){
      maxBin = nXBins - bin;
      foundMax = true;
      break;
    }
  }

  //the array has a one item more then the number of bins.

  for(int bin = 1; bin<=maxBin; bin++){
    IntX = fHistogram->IntegralAndError(bin,bin,0,-1,Err);
    if(isEmpty){
      if(IntX){
	isEmpty = false;
	firstBin = (bin==1) ? 1 : bin - 1;
      }
      else continue;
    }
    if(isFirst ||( (Err/IntX < MaxErr) || (bin == maxBin) ) ){
	isFirst = false;
      Cont.push_back(vector<Double_t>());
      xBins.push_back(arrX[bin-1]);      

     
	for(int biny = 1; biny<=nYBins; biny++){
	  Cont.back().push_back(fHistogram->GetBinContent(bin,biny));
	}
    }
    else{
      for(int biny = 1; biny<=nYBins; biny++){
	(Cont.back())[biny-1]+=fHistogram->GetBinContent(bin,biny);
      }
    }
  }
  
  xBins.insert(xBins.begin(),arrX[firstBin-1]);
  Cont.insert(Cont.begin(),vector<Double_t>());

  for(int biny = 1; biny<=nYBins; biny++) Cont[0].push_back(0.);
  
  xBins.push_back(arrX[maxBin]);

   xBins.push_back(arrX[maxBin+1]);
  Cont.push_back(vector<Double_t>());
  for(int biny = 1; biny<=nYBins; biny++) Cont.back().push_back(0.);

  TH2F * hNew = new TH2F("","",xBins.size()-1,xBins.data(),nYBins,arrY.data());

  fHistogram = hNew;
  fHistogram->SetName(theName.c_str());

  for(uint binx = 1; binx<=xBins.size()-1; binx++){
    for(int biny = 1; biny<=nYBins; biny++){
      fHistogram->SetBinContent(binx,biny, Cont[binx-1][biny-1]);
    }
  }
}


void Hist::getIntLumiHisto(){

  Int_t nXBins  = fHistogram->GetNbinsX();       
  Int_t nYBins  = fHistogram->GetNbinsY();       

  const  TArrayD *arrY = fHistogram->GetYaxis()->GetXbins();
  const  TArrayD *arrX = fHistogram->GetXaxis()->GetXbins();

  std::vector<Double_t> xBins = {};

  for(int bin = 1; bin <= nXBins; bin++)  xBins.push_back(lumi::getTotLumiRun(to_string(static_cast<int>(arrX->GetArray()[bin-1]))));

  TH2F * hNew = new TH2F("","",100*xBins.size(), xBins.at(0)*(1+1/15.) - xBins.back()/15. ,  xBins.back()*(1+1/12.) - xBins.at(0)/12. ,nYBins,arrY->GetArray());
  
  string hName = string(fHistogram->GetName());

  for(int binx = 1; binx<=nXBins; binx++){  
    for(int biny = 1; biny<=nYBins; biny++){  
      hNew->AddBinContent(hNew->FindBin(xBins[binx-1],arrY->GetArray()[biny-1]), fHistogram->GetBinContent(binx,biny));
    }
  }

  hNew->GetXaxis()->SetTitle("Int.Lumi. pb^{-1}");
  hNew->GetYaxis()->SetTitle("Rate(Hz/cm^{-2})");
  hNew->SetName(hName.c_str());
  hNew->SetLineColor(fHistogram->GetLineColor());
  hNew->SetMarkerStyle(fHistogram->GetMarkerStyle());
  hNew->SetMarkerColor(fHistogram->GetMarkerColor());

  fHistogram = hNew;
}


void Hist::setEqualBin(vector<double> slices){

  Int_t nXBins  = fHistogram->GetNbinsX();       
  Int_t nYBins  = fHistogram->GetNbinsY();       

  const  TArrayD *arrY = fHistogram->GetYaxis()->GetXbins();
  const  TArrayD *arrX = fHistogram->GetXaxis()->GetXbins();

  std::vector<Double_t> xBins;

  for(int bin = 1; bin<=nXBins; bin++){  
    xBins.push_back(arrX->GetArray()[bin-1]);
  }
 
  TH2F * hNew = new TH2F("","",xBins.size(),xBins[0],xBins.back(),nYBins,arrY->GetArray());
  
  for(uint binx = 1; binx<=xBins.size(); binx++){  
    for(int biny = 1; biny<=nYBins; biny++){  
      hNew->SetBinContent(binx,biny, fHistogram->GetBinContent(binx,biny));
    }
  }

  for(int bin = 1; bin <=nXBins; bin++){
    hNew->GetXaxis()->SetBinLabel(bin,(to_string(static_cast<int>(xBins[bin-1]) ).c_str() ));
    }
  
  fHistogram = hNew; 
  fHistogram->SetName(theName.c_str());
}

void Hist::SetColor(Color_t mcolor){
  fHistogram->SetMarkerColor(mcolor);
  fHistogram->SetLineColor(mcolor);
}

Color_t Hist::GetColor(){
  return fHistogram->GetMarkerColor();
}

void Hist::SetMarkerStyle(Style_t mstyle){
  fHistogram->SetMarkerStyle(mstyle);
}

