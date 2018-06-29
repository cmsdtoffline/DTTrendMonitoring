#ifndef HIST
#define HIST
 
//standard header
#include <vector>
#include <utility>
#include <TROOT.h> 
#include <TH2.h>
#include <TGraphAsymmErrors.h>
#include <TNamed.h>
#include <TColor.h>
//#include <RtypesCore.h>

class Hist: public TNamed 
{

 public:
  
  
  Hist(const char* name,const char* title,Int_t nbinsx,
	const Double_t* xbins,Int_t nbinsy,const Double_t* ybins);

  Hist();

  ~Hist();

  void Fill(Double_t x,Double_t y);
  
  void addHisto(vector<double> slices);

  Int_t GetNbinsX() {return fHistogram->GetNbinsX();};

  float Integral() {return fHistogram->Integral();};

  void ProfileX();

  void draw(string option);

  float GetProfMax();

  void getIntLumiHisto();

  void SetColor(Color_t mcolor = 1);

  Color_t GetColor();
  
  void SetMarkerStyle(Style_t mstyle = 1);

  void setEqualBin(vector<double> slices);

  void set2DHistoBin(Float_t MaxErr= 0.1);

  void SetMaximum(float Max);

  const char *getName() {return theName.c_str();};

 private:

  TH2* fHistogram; 
  TProfile *fProf;

  string theName;
  ClassDef(Hist,1)

};

#endif
