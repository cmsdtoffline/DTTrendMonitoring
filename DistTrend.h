#ifndef DISTTREND
#define DISTTREND


/* class DistTrend
 Jun. 19, 2018
 Gian Luca Pinna Angioni 
*/


#include <vector>
#include <utility>
#include <TROOT.h> 
#include <TH2.h>
#include <TGraphAsymmErrors.h>
#include <TNamed.h>
#include <TColor.h>
//#include <RtypesCore.h>

class DistTrend: public TNamed 
{

 public:
    
  DistTrend(const char* name,const char* title,Int_t nbinsx,
	const Double_t* xbins,Int_t nbinsy,const Double_t* ybins);

  DistTrend(const char* name,const char* title,Int_t nbinsx,
      Float_t X0, Float_t X1, Int_t nbinsy,Float_t Y0, Float_t Y1);
    
  DistTrend();

  ~DistTrend();

  void Fill(Double_t x,Double_t y);
  
  void addHisto(vector<double> slices);

  Int_t GetNbinsX() {return fHistogram->GetNbinsX();};

  float Integral() {return fHistogram->Integral();};

  vector<double>  GetHistoArrayX();
  vector<double>  GetHistoArrayY();
  
  void ProfileX();

  void draw(string option);

  void drawHisto(string  option);

  float GetProfMax();

  void setTitle(const char *title);

  void getIntLumiHisto();

  void SetColor(Color_t mcolor = 1);

  Color_t GetColor();

  TH2F *GetHisto();
  
  void SetMarkerStyle(Style_t mstyle = 1);

  void setEqualBin(vector<double> slices);

  void set2DHistoBin(Float_t MaxErr= 0.1);

  void SetMaximum(float Max);

  const char *getName() {return theName.c_str();};

 private:

  TH2F* fHistogram; 
  TProfile *fProf;

  string theName;
  ClassDef(DistTrend,1)

};

#endif
