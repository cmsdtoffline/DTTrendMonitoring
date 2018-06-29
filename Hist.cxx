//standard header

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
  fProf = fHistogram->ProfileX();
}

void Hist::draw(string  option){
  fProf->Draw(option.c_str());
  //fHistogram->Draw(option.c_str());
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



void Hist::set2DHistoBin(Float_t MaxErr){
  
  Int_t nXBins  = fHistogram->GetNbinsX();
  Int_t nYBins  = fHistogram->GetNbinsY();

  const  TArrayD *arrX = fHistogram->GetXaxis()->GetXbins();
  const  TArrayD *arrY = fHistogram->GetYaxis()->GetXbins();

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
	xBins.push_back(arrX->GetArray()[bin-1]);
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
  
  xBins.insert(xBins.begin(),arrX->GetArray()[firstBin-1]);
  Cont.insert(Cont.begin(),vector<Double_t>());


  for(int biny = 1; biny<=nYBins; biny++) Cont[0].push_back(0.);
  
  xBins.push_back(arrX->GetArray()[maxBin]);
  xBins.push_back(arrX->GetArray()[maxBin+1]);
  Cont.push_back(vector<Double_t>());
  for(int biny = 1; biny<=nYBins; biny++) Cont.back().push_back(0.);


  TH2F * hNew = new TH2F("","",xBins.size()-1,xBins.data(),nYBins,arrY->GetArray());
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
      hNew->AddBinContent(hNew->FindBin(xBins[binx-1],biny), fHistogram->GetBinContent(binx,biny));
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
  //Check
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
    cout<<xBins[bin-1]<<endl;
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


//#endif
