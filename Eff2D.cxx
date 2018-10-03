#include <iostream>
#include "Eff2D.h"
#include <vector>
#include <utility>
#include <TROOT.h> 
#include <TH2.h>
#include <TColor.h>
#include "TEfficiency.h"
#include "Utilities.h"
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/algorithm/string.hpp>
#include "TPad.h"

ClassImp(Eff2D)

Eff2D::Eff2D(const char* name,const char* title,Int_t nbinsx,
	     const Double_t* xbins,Int_t nbinsy,const Double_t* ybins) :TNamed(name,title){
  theName = name;

  fPassedHistogram = new TH2F((theName+"_pass").c_str(),title,nbinsx,xbins,nbinsy,ybins);
  fTotalHistogram  = new TH2F((theName+"_tot").c_str(),title,nbinsx,xbins,nbinsy,ybins);
  xLabel = "";
  yLabel = "";
  Title  = title;
  eff = new TGraphAsymmErrors();
}

Eff2D::Eff2D() :TNamed("","")
{
  theName = "none";
  fPassedHistogram = new TH2F();
  fTotalHistogram  = new TH2F();
  eff = new TGraphAsymmErrors();
}

Eff2D::~Eff2D() 
{
}


void Eff2D::Fill(Bool_t bPassed, Double_t x, string test = ""){ 
  fTotalHistogram->Fill(x,1);
  if(bPassed)  fPassedHistogram->Fill(x,1);
}

void Eff2D::Fill(Bool_t bPassed, Double_t x, Double_t y){ 

  fTotalHistogram->Fill(x,y);
  if(bPassed)  fPassedHistogram->Fill(x,y);
}



void Eff2D::setTitle(const char *title){

  vector<string> strs;
  boost::split(strs,title,boost::is_any_of(";"));
  
  if(strs.size()>=3){
    eff->SetTitle(strs[0].c_str());
    xLabel=strs[1];
    yLabel=strs[2];
  }  
  else if(strs.size()==2){
    eff->SetTitle(strs[0].c_str());
    xLabel=strs[1];
  }  
  else if(strs.size()==1)
    eff->SetTitle(strs[0].c_str());
}



void Eff2D:: draw(Int_t  firstxbin,
		  Int_t  lastxbin, string option){
    
  const  TH1D *hPass =  fPassedHistogram->ProjectionX("_passpx",firstxbin,lastxbin);
  const  TH1D *hTot  =  fTotalHistogram->ProjectionX("_totpx",firstxbin,lastxbin);

  TGraphAsymmErrors * effNew = new TGraphAsymmErrors(hPass,hTot,"cp");

  effNew->SetMarkerColor(eff->GetMarkerColor());
  effNew->SetLineColor(eff->GetLineColor());
  effNew->SetMarkerStyle(20);
  effNew->SetTitle(eff->GetTitle());

  effNew->GetXaxis()->SetTitle(xLabel.c_str());
  effNew->GetYaxis()->SetTitle(yLabel.c_str());
  effNew->SetMinimum(0.89);
  effNew->SetMaximum(1.02);

  effNew->SetName((hPass->GetName()+theName).c_str());
  eff = effNew;
  eff->Draw(option.c_str());

}




void Eff2D:: drawWithLumi(vector<double> slices,Int_t  firstxbin,
			  Int_t  lastxbin, string option){
  
  std::vector<Double_t> xBins = {};
  Int_t nXBins  = fPassedHistogram->GetNbinsX();
  Int_t nYBins  = fPassedHistogram->GetNbinsY();
  
  const  TArrayD *arrY = fPassedHistogram->GetYaxis()->GetXbins();
  for(int bin = 1; bin <= nXBins; bin++){
    xBins.push_back(lumi::getTotLumiRun(to_string(static_cast<int>(slices[bin-1])))); 
  }
  
  float redCon = 400;
  
  TH2F * hTot = new TH2F("","",100*xBins.size(), xBins.at(0)*(1-1/redCon) ,  xBins.back()*(1+1/redCon),nYBins,arrY->GetArray());
  TH2F * hPass = new TH2F("","",100*xBins.size(), xBins.at(0)*(1-1/redCon),  xBins.back()*(1+1/redCon),nYBins,arrY->GetArray());
  
  for(int binx = 1; binx<=nXBins; binx++){  
    for(int biny = 1; biny<=nYBins; biny++){  
      hTot->SetBinContent(hTot->FindBin(xBins[binx-1],arrY->GetArray()[biny-1]), fTotalHistogram->GetBinContent(binx,biny));
      hPass->SetBinContent(hPass->FindBin(xBins[binx-1],arrY->GetArray()[biny-1]), fPassedHistogram->GetBinContent(binx,biny));    
    }
  }  
  
  const  TH1D *hPrPass =  hPass->ProjectionX("_passpx",firstxbin,lastxbin);
  const  TH1D *hPrTot  =  hTot->ProjectionX("_totpx",firstxbin,lastxbin);
  
  TGraphAsymmErrors * effNew = new TGraphAsymmErrors(hPrPass,hPrTot,"cp");
  
  effNew->SetMarkerColor(eff->GetMarkerColor());
  effNew->SetLineColor(eff->GetLineColor());
  effNew->SetMarkerStyle(20);

  effNew->GetYaxis()->SetTitle(yLabel.c_str()); 
  effNew->SetTitle(eff->GetTitle());

  eff = effNew;

  eff->GetXaxis()->SetTitle(xLabel.c_str());
  eff->SetMinimum(0.89);
  eff->SetMaximum(1.02);

  eff->SetName((hPass->GetName()+theName).c_str());
  eff->Draw(option.c_str());

}

void Eff2D::setEffRun(){
                 

  Int_t nXbins = fTotalHistogram->GetNbinsX();
  Int_t nYbins = fTotalHistogram->GetNbinsY();

  const  TArrayD *arrY = fPassedHistogram->GetYaxis()->GetXbins();

  TH2F * hNewPass = new TH2F("","",nXbins,0,nXbins,nYbins,arrY->GetArray());
  TH2F * hNewTot  = new TH2F("","",nXbins,0,nXbins,nYbins,arrY->GetArray());
  
  for(int binx = 1; binx<=nXbins; binx++){  
    for(int biny = 1; biny<=nYbins; biny++){  
      hNewPass->SetBinContent(binx,biny, fPassedHistogram->GetBinContent(binx,biny));
      hNewTot->SetBinContent(binx,biny, fTotalHistogram->GetBinContent(binx,biny));
    }
  }

  fPassedHistogram = dynamic_cast<TH2F*>(hNewPass->Clone());
  fTotalHistogram  = dynamic_cast<TH2F*>(hNewTot->Clone());

}

const TArrayD * Eff2D::GetArrayX(){
  return fTotalHistogram->GetXaxis()->GetXbins();
}

void Eff2D::addBins(vector<double> slices){
  
  Int_t nXBins  = fPassedHistogram->GetNbinsX();       
  Int_t nYBins  = fPassedHistogram->GetNbinsY();       
  
  const  TArrayD *arrY = fPassedHistogram->GetYaxis()->GetXbins();
  TH2F * hNewPass = new TH2F(("pass_"+theName).c_str(),"",slices.size()-1,slices.data(),nYBins,arrY->GetArray());
  TH2F * hNewTot  = new TH2F(("tot_"+theName).c_str(),"",slices.size()-1,slices.data(),nYBins,arrY->GetArray());
 
  for(int binx = 1; binx<=nXBins; binx++){  
    for(int biny = 1; biny<=nYBins; biny++){  
      hNewPass->SetBinContent(binx,biny, fPassedHistogram->GetBinContent(binx,biny));
      hNewTot->SetBinContent(binx,biny, fTotalHistogram->GetBinContent(binx,biny));
    }
  }

  fPassedHistogram = hNewPass;
  fTotalHistogram  = hNewTot;
}

void Eff2D::SetColor(Color_t mcolor){
  eff->SetLineColor(mcolor);
  eff->SetMarkerColor(mcolor);

}


Color_t Eff2D::GetColor(){
  return eff->GetMarkerColor();
}

void Eff2D::SetMarkerStyle(Style_t mstyle){
  eff->SetMarkerStyle(mstyle);
}

void Eff2D::getIntLumiEff( vector<double> slices){

  std::vector<Double_t> xBins = {};
  Int_t nXBins  = fPassedHistogram->GetNbinsX();
  Int_t nYBins  = fPassedHistogram->GetNbinsY();

  const  TArrayD *arrY = fPassedHistogram->GetYaxis()->GetXbins();

  for(int bin = 1; bin <= nXBins; bin++){

    xBins.push_back(lumi::getTotLumiRun(to_string(static_cast<int>(slices[bin-1])))); 
  }
  
  float redCon = 400;

  TH2F * hTot = new TH2F("","",100*xBins.size(), xBins.at(0)*(1-1/redCon) ,  xBins.back()*(1+1/redCon),nYBins,arrY->GetArray());
  TH2F * hPass = new TH2F("","",100*xBins.size(), xBins.at(0)*(1-1/redCon),  xBins.back()*(1+1/redCon),nYBins,arrY->GetArray());
 												

  for(int binx = 1; binx<=nXBins; binx++){  
    for(int biny = 1; biny<=nYBins; biny++){  
      hTot->SetBinContent(hTot->FindBin(xBins[binx-1],arrY->GetArray()[biny-1]), fTotalHistogram->GetBinContent(binx,biny));
      hPass->SetBinContent(hPass->FindBin(xBins[binx-1],arrY->GetArray()[biny-1]), fPassedHistogram->GetBinContent(binx,biny));
    }
  }  

 fPassedHistogram =  dynamic_cast<TH2F*>(hPass->Clone());
 fTotalHistogram  =  dynamic_cast<TH2F*>(hTot->Clone());
}



void Eff2D::setEffBin(Float_t MaxErr){
  
  Int_t nXBins  = fTotalHistogram->GetNbinsX();
  Int_t nYBins  = fTotalHistogram->GetNbinsY();

  const  TArrayD *arrX = fTotalHistogram->GetXaxis()->GetXbins();
  const  TArrayD *arrY = fTotalHistogram->GetYaxis()->GetXbins();

  std::vector<Double_t> xBins;
  std::vector<vector<Double_t>> ContTot;
  std::vector<vector<Double_t>> ContPass;

  Double_t  Err = 0;
  Double_t IntX = 0;
  
  Int_t maxBin = -1;
  Int_t firstBin = -1;

  bool  foundMax = false;
  bool  isEmpty  = true;
  bool  isFirst  = true;

  for(int bin = 1; bin<=nXBins; bin++){
    //loop backward to find last non empty column

    if(!foundMax && fTotalHistogram->Integral((nXBins - bin),(nXBins - bin),0,-1 ) >0.){
      maxBin = nXBins - bin;
      foundMax = true;
      break;
    }
  }

  //the array has a one item more then the number of bins.
  for(int bin = 1; bin<=maxBin; bin++){

    IntX = fTotalHistogram->IntegralAndError(bin,bin,0,-1,Err);
    if(isEmpty){
      if(IntX){
	isEmpty = false;
	firstBin = (bin==1) ? 1 : bin - 1;
      }
      else continue;
    }


   if(isFirst || (Err/IntX < MaxErr) || bin ==maxBin){

	isFirst = false;
	ContTot.push_back(vector<Double_t>());
	ContPass.push_back(vector<Double_t>());

	xBins.push_back(arrX->GetArray()[bin-1]);
	for(int biny = 1; biny<=nYBins; biny++){
	  ContTot.back().push_back(fTotalHistogram->GetBinContent(bin,biny));
	  ContPass.back().push_back(fPassedHistogram->GetBinContent(bin,biny));
	}
    }
    else{
      for(int biny = 1; biny<=nYBins; biny++){
	(ContTot.back())[biny-1]+=fTotalHistogram->GetBinContent(bin,biny);
	(ContPass.back())[biny-1]+=fPassedHistogram->GetBinContent(bin,biny);
      }
    }
  }
  
  xBins.insert(xBins.begin(),arrX->GetArray()[firstBin-1]);
  ContTot.insert(ContTot.begin(),vector<Double_t>());
  ContPass.insert(ContPass.begin(),vector<Double_t>());

  for(int biny = 1; biny<=nYBins; biny++) ContTot[0].push_back(0.);
  for(int biny = 1; biny<=nYBins; biny++) ContPass[0].push_back(0.);
  
  xBins.push_back(arrX->GetArray()[maxBin]);
  xBins.push_back(arrX->GetArray()[maxBin+1]);

  ContTot.push_back(vector<Double_t>());
  ContPass.push_back(vector<Double_t>());

  for(int biny = 1; biny<=nYBins; biny++){
    ContTot.back().push_back(0.);
    ContPass.back().push_back(0.);
  }

  TH2F * hTot = new TH2F("","",xBins.size()-1,xBins.data(),nYBins,arrY->GetArray());
  TH2F * hPass = new TH2F("","",xBins.size()-1,xBins.data(),nYBins,arrY->GetArray());


  fTotalHistogram = hTot;
  fPassedHistogram = hPass;

  fTotalHistogram->SetName(theName.c_str());

  for(uint binx = 1; binx<=xBins.size()-1; binx++){
    for(int biny = 1; biny<=nYBins; biny++){
      fTotalHistogram->SetBinContent(binx,biny, ContTot[binx-1][biny-1]);
      fPassedHistogram->SetBinContent(binx,biny, ContPass[binx-1][biny-1]);
    }
  }
}

bool Eff2D::checkProj(vector<double> slices, int bin,int bin2){
  return fTotalHistogram->ProjectionX("_px",bin,bin2)->Integral()==0;
}

//Create histogram with equal bins to be used as graphic plot instead of tgraph.
TH1D * Eff2D::getPaintHisto(vector<double> slices, bool doLumi, bool isRun){

  TH1D * hPaint = new TH1D();
  
  hPaint = fPassedHistogram->ProjectionX("_px",0,-1);
  hPaint->Reset();

  Int_t nBins = hPaint->GetNbinsX();
  
  if(doLumi) xLabel =  "Int.Lumi. pb^{-1}";
  else if(isRun){
    for(int bin = 1; bin <= nBins; bin++){
      hPaint->GetXaxis()->SetBinLabel(bin,(to_string(static_cast<int>(slices[bin-1]))).c_str()); 
    }
  }

  hPaint->SetMinimum(0.89);
  hPaint->SetMaximum(1.02);
  
  hPaint->GetXaxis()->SetTitle(xLabel.c_str());
  hPaint->GetYaxis()->SetTitle(yLabel.c_str());
  
  return hPaint;
}
