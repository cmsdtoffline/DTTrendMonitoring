#ifndef EFF2D
#define EFF2D
 
//standard header
#include <vector>
#include <utility>
#include <TROOT.h> 
#include <TH2.h>
//#include <TObject.h>
#include "TEfficiency.h"
#include <TNamed.h>
#include <TColor.h>


class Eff2D: public TNamed 
{

 public:
  
  
  Eff2D(const char* name,const char* title,Int_t nbinsx,
	const Double_t* xbins,Int_t nbinsy,const Double_t* ybins);

  Eff2D();

  ~Eff2D();

  void Fill(Bool_t bPassed, Double_t x, string test);

  void Fill(Bool_t bPassed,Double_t x,Double_t y);
  
  void addBins(vector<double> slices);

  void draw(Int_t  firstxbin = 0,
		  Int_t  lastxbin = -1, string option= "");

  void drawWithLumi(vector<double> slices,
		    Int_t  firstxbin,
		    Int_t  lastxbin, string option);

  //  void draw(string option = "");

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

  //  TGraphAsymmErrors * GetPaintedGraph();
//  TEfficiency * getEffProj(  Int_t  firstxbin = 0,
  //Int_t  lastxbin = -1 );

  TH1D * getPaintHisto(vector<double> slices, bool doLumi, bool isRun = false);
  //  TEfficiency * get2DEff(vector<double> slices);


 private:

  TH2* fPassedHistogram; 
  TH2* fTotalHistogram; 

  TGraphAsymmErrors * eff;

  string theName;
  string xLabel;
  string yLabel;
  string Title;

  ClassDef(Eff2D,1)

};

#endif
