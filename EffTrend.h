#ifndef EFFTREND
#define EFFTREND

/* class EffTrend                                                                                                                                                                                                                               
 Jun. 18, 2018                                                                                                                                                                       Gian Luca Pinna Angioni                                                                                                                                                                                                                    
*/

#include <vector>
#include <utility>
#include <TROOT.h> 
#include <TH2.h>
#include "TEfficiency.h"
#include <TNamed.h>
#include <TColor.h>


class EffTrend: public TNamed 
{


 public:
  
  
  EffTrend(const char* name,const char* title,Int_t nbinsx,
	const Double_t* xbins,Int_t nbinsy,const Double_t* ybins);

  EffTrend();

  ~EffTrend();

  void Fill(Bool_t bPassed, Double_t x, string test);

  void Fill(Bool_t bPassed,Double_t x,Double_t y);
  
  void addBins(vector<double> slices);

  void draw(Int_t  firstxbin = 0,
		  Int_t  lastxbin = -1, string option= "");

  void drawWithLumi(vector<double> slices,
		    Int_t  firstxbin,
		    Int_t  lastxbin, string option, bool plotRuns = false);

  bool checkProj(vector<double> slices, int bin, int bin2);

  Int_t GetNbinsX() {return fTotalHistogram->GetNbinsX();};

  void SetColor(Color_t mcolor = 1);  

  void SetMarkerStyle(Style_t mstyle = 1);

  void getIntLumiEff( vector<double> slices);

  void setEffBin(Float_t MaxErr = .1);

  void setEffRun();

  float Integral() {return fTotalHistogram->Integral();};

  void setTitle(const char *title);
  
  const TArrayD * GetArrayX();

  Color_t GetColor();

  TH1D * getPaintHisto(vector<double> slices, bool doLumi, bool isRun = false);

 private:

  TH2* fPassedHistogram; 
  TH2* fTotalHistogram; 

  TGraphAsymmErrors * eff;

  string theName;
  string xLabel;
  string yLabel;
  string Title;

  ClassDef(EffTrend,1)

};

#endif
