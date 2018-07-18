#ifndef Utilities_h
#define Utilities_h

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "Types.h"

using boost::property_tree::ptree;
using boost::property_tree::read_json;
using boost::property_tree::write_json;

namespace pt = boost::property_tree;

template <typename T>
std::vector<T> as_vector(ptree const& pt, ptree::key_type const& key)

{
  std::vector<T> r;
  for (auto& item : pt.get_child(key))
    r.push_back(item.second.get_value<T>());
  return r;
}


namespace lumi{
  
  float getTotLumiRun(string run){
    
    string year = "2018";
    string path = "data/IntLumi/";
    string inFileTotal = path + "Total"+year+".json";
    
    pt::ptree pTree;
    pt::read_json(inFileTotal.c_str(), pTree);
    vector<float> runInfos;
    runInfos =  as_vector<float>(pTree,run.c_str());
    cout.precision(100000);
    return runInfos.at(4);  
    
  }
}


namespace variables{
  
  bool initVar(Var & var, string varname, string filename, bool DoEff = 1, bool DoBkg =1, bool doProj = 0, string projVar = "", bool isExtProj = 1){    

    pt::ptree pTree;
    pt::read_json(filename.c_str(), pTree);
    
    var.name  =  pTree.get<std::string>( (varname+".name").c_str());
    var.Title =  pTree.get<std::string>( (varname+".Title").c_str()); 
    var.Label =  pTree.get<std::string>( (varname+".Label").c_str()); 
    var.slice =  as_vector<double>(pTree,(varname+".slice").c_str());

    var.doEff = DoEff;
    var.doBkg = DoBkg;
    var.doProj = doProj;

    if(doProj){
      var.projVar = projVar;
      if(isExtProj){
	var.projSlice = as_vector<double>(pTree,(projVar+".slice").c_str());
      }
      else{
	var.projSlice = as_vector<double>(pTree,(varname+".projSlice").c_str());
      }
    } 

    else{ //fake variable.
      var.projVar  = "Empty";
      var.projSlice = {0}; 
    }

    
    return 1;
    
  }
}


#endif



