#ifndef Types_h
#define Types_h

#include <iostream> 
#include <map>

struct Var{

  Int_t nVar;
  std::string name;
  std::vector<Double_t >  slice;
  std::vector<Double_t >  projSlice;
  std::string Title;
  std::string Label;
  bool doEff;
  bool doBkg;
  float value;
};


struct context{

  Int_t nVar;
  std::string name;
  std::map<string,Var> var;
  string webFolder;
  
};


#endif



