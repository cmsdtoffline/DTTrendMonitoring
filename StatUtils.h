//#include "Math/Math.h"
//#include <Math.h>
#include "Math/SpecFunc.h"
#include "Math/DistFunc.h"

Float_t GetError(int n, int k, bool isUp, float level=0.68540158589942957){
  float alpha = (1.0 - level)/2;
  float xhigh = 0;
  float xlow  = 0;
  if(k<=0) xlow = 0.;
  else xlow = ROOT::Math::beta_quantile(alpha,k,n-k+1.);
  if(k==n) xhigh = 1.;
  else  xhigh = ROOT::Math::beta_quantile(1. - alpha,k+1.,n-k);
  if(isUp){
    //std::cout<<"high err "<<xhigh<<std::endl;
    return xhigh - k/float(n);
  }
  else{
    //    std::cout<<"low err "<<xlow<<std::endl;
    return k/float(n) - xlow;
  }
}
