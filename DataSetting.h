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



std::vector<Double_t>  bkgslice  = {0., 1e-3, 2e-3, 4e-3, 6e-3, 8e-3,
				    1e-2, 2e-2, 4e-2, 6e-2, 8e-2,
				    1e-1, 2e-1, 4e-1, 6e-1, 8e-1,
				    1, 2, 4, 6, 8,
				    1e1, 2e1, 4e1, 6e1, 8e1,
				    1e2, 2e2, 4e2, 6e2, 8e2,
				    1e3, 2e3, 4e3, 6e3, 8e3,
				    1e4, 2e4, 4e4, 6e4, 8e4,
				    1e5, 2e5, 4e5, 6e5, 8e5};



std::vector<std::string> varName_inst  = {"InsLumi",                
					  "PU",     
					  "bkg"};
std::vector<std::string> varTitle_inst = {"Instantaneus Luminosity",
					  "Pile-up",
					  "background"};
std::vector<std::string> varLabel_inst = {"Inst. Luminosity (cm^{-2}s^{-1}10^{30})",
					  "Pu",
					  "Rate(Hz/cm^{-2})"};


Int_t nVar = 3;

std::string WebFolder = "~/www/DT";


//RUN

std::vector<Double_t>  runslice         = {};
std::vector<string>    runsliceLabels     = {};
std::vector<Double_t>  bunchXingslice   = {0,500.,1000.,1500.,2000.,3000.,3500.,4000.};

std::vector<std::string> varName_run    = {"Run","bunchXing","bkg"};             

std::vector<std::string> varTitle_run   = {"Run","bunch crossing","bkg"};

std::vector<std::string> varLabel_run   = {"Run number","bx","bkg"};




