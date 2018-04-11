#include <iostream>



std::vector<Double_t> lumislice = { 0., 1000., 2000., 3000., 4000., 
				    5000., 6000., 7000., 8000., 9000.,
				    10000.,11000.,12000.,13000.,14000.,
				    15000.,16000.,17000.,18000.,19000.,
				    20000.,21000.,22000.,23000.,24000.,
				    25000};

std::vector<Double_t> PUslice   = { 0., 2., 4., 6.,8.,
				    10.,12.,14.,16.,18.
				    ,20.,22.,24.,26.,28.
				    ,30.,32.,34.,36.,38.,
				    40.,44.,48.,50.,58.,
				    62.,72.,82.,94.,100.,110.,120.};

std::vector<Double_t>  bkgslice   = { 0., 2., 4., 6.,8.,
				    10.,12.,14.,16.,18.
				    ,20.,22.,24.,26.,28.
				    ,30.,32.,34.,36.,38.,
				    40.,44.,48.,50.,58.,
				    62.,72.,82.,94.};



std::vector<std::string> varName_inst  = {"InsLumi",                
					  "PU",     
					  "bkg"};
std::vector<std::string> varTitle_inst = {"Instantaneus Luminosity",
					  "Pile-up",
					  "background"};
std::vector<std::string> varLabel_inst = {"Inst. Luminosity (cm^{-2}s^{-1}10^{30})",
					  "Pu",
					  "bkg yield"};


Int_t nVar = 3;

std::string WebFolder = "~/www/DT";


//RUN

std::vector<Double_t>  runslice         = {};
std::vector<string>  runsliceLabels     = {};
std::vector<Double_t>  bunchXingslice   = {0,500.,1000.,1500.,2000.,3000.,3500.,4000.};

std::vector<std::string> varName_run    = {"Run","bunchXing","bkg"};             

std::vector<std::string> varTitle_run   = {"Run","bunch crossing","bkg"};

std::vector<std::string> varLabel_run   = {"Run number","bx","bkg"};




