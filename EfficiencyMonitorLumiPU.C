#define EfficiencyMonitor_cxx
#include "EfficiencyMonitor.h"
#include <TVectorF.h>
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include "StatUtils.h"
#include "TEfficiency.h"
#include "EfficiencyMonitorSet.h"
int dead[10000][6];
int Ndead=0;
int nrequiredhit=7;
int MaxDead=0;

/* 
note: MaxDead = maximum number of dead channels crossed by a selected segment
MaxDead = 0 allows an umbiased selection.
If MaxDead > 0, dead wires will have efficiency zero while alive wires in the same segment won't contribute!
*/

int lessentries=1000000;


void EfficiencyMonitor::PreLoop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   //nentries=lessentries;

   Long64_t nbytes = 0, nb = 0;

   ofstream DeadList;
   std::string deadname;

   if(fileName!="")  deadname.append("DeadList_"+fileName);
   else{
     deadname.append("DeadList_Run2016");
     deadname = deadname + dataset + ".txt";;
   }   cout<<deadname<<endl;
   DeadList.open (deadname.c_str());


   //Create a Histogram per layer 
   char go;
   char hname[50];
   TH1F* occupancy[5][14][4][3][4];


   for (int iwh=0; iwh<5; iwh++) {
     for (int ise=0; ise<14; ise++) {
       for (int ist=0; ist<4; ist++) {
         if (ist!=3 && ise>11) continue;
         for (int isl=0; isl<3; isl++) {
           if (isl==1 && ist==3) continue;
           for (int ilay=0; ilay<4; ilay++) {
             sprintf (hname,"Wheel%u_Sect%u_MB%u_SL%u_Lay%u",iwh,ise+1,ist+1,isl+1,ilay+1);
             occupancy[iwh][ise][ist][isl][ilay] = new TH1F (hname,"",92,1.,93.);
	   }
	 }
       }
     }
   }

   for (Long64_t jentry=2; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      std::cout<<"Entries "<<ientry<<std::endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%5000 == 0) cout<<" Pre-Loop evento "<<jentry<<endl;
  
      for (int idigi=0; idigi<Ndigis; idigi++) {
   
	occupancy[digi_wheel->at(idigi)+2][digi_sector->at(idigi)-1][digi_station->at(idigi)-1]
	         [digi_sl->at(idigi)-1][digi_layer->at(idigi)-1]->Fill(float(digi_wire->at(idigi)));
      }
   }
    
   // analyze occupancy histos and fill dead channel table

   int nwire=0; int NwireTot=0;
   for (int iw=0; iw<5000; iw++) for (int geo=0; geo<6; geo++) dead[iw][geo]=0;
   for (int iwh=0; iwh<5; iwh++) {
     for (int ise=0; ise<14; ise++) {
       for (int ist=0; ist<4; ist++) {
         if (ist!=3 && ise>11) continue;
         for (int isl=0; isl<3; isl++) {
           if (isl==1 && ist==3) continue;
           if (isl==1){
	     //chimney chambers
	     if (( iwh==1 && ise==2 ) || ( iwh==3 && ise==3 ) ) nwire = 48;
	     else  nwire=57; 
	   }
           else if (ist==0) nwire = 49;
           else if (ist==1) nwire = 60;
           else if (ist==2) nwire = 72;
           else if (ist == 3) {
	     // Following conditions are put in order of decrising average of number of entries per chamber to increase velocity.
	     if (ise == 3 || ise == 12 )    nwire = 72;
             else if ( ise==7 || ise ==11 ) nwire = 92;
             else if ( ise==9 || ise ==13 ) nwire = 60;
             else if ( ise==8 || ise ==10 ) nwire = 48;
             else nwire = 96;
	   }

           NwireTot+=(nwire*4);

           for (int ilay=0; ilay<4; ilay++) {
             for (int iw=1; iw<nwire+1; iw++) {
               if (occupancy[iwh][ise][ist][isl][ilay]->GetBinContent(iw)==0){
		 dead[Ndead][0]=iwh-2;
		 dead[Ndead][1]=ise+1;
 		 dead[Ndead][2]=ist+1;
		 dead[Ndead][3]=isl+1;
		 dead[Ndead][4]=ilay+1;
		 dead[Ndead][5]=iw;  
                 DeadList <<Ndead+1<<" ";
                 for (int ip=0; ip<6; ip++) DeadList<< dead[Ndead][ip]<<" "; DeadList<<endl;
                 Ndead++; 
	       } 
	     }
	   }
	 }
       }
     }
   }

   cout<<Ndead<<" dead wires out of "<<NwireTot<<endl;

   DeadList<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "
          <<endl<<Ndead<<" dead wires of out "<<NwireTot<<endl;

   for (int idead=0; idead<Ndead; idead++) {
     cout<<"YB"<<dead[idead][0]<<"/Sec"<<dead[idead][1]<<"/MB"<<dead[idead][2]<<"_SL"<<dead[idead][3]
         <<"_layer"<<dead[idead][4]<<"_cell"<<dead[idead][5]<<endl;
   }
   
   DeadList.close();
}


void EfficiencyMonitor::Loop()
{

  //   TH1F *hextra = new TH1F("hextra","hextra",2000,0,2000); //del
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   char go;

   //  nentries=lessentries;
   std::cout<<"Entries "<<nentries<<std::endl;

   Long64_t nbytes = 0, nb = 0;

   slices.push_back(lumislice);
   slices.push_back(PUslice);
   slices.push_back(bkgslice);


   // COUNTERS:

   const int nLumiPoints = slices[0].size();
   const int nPUpoints   = slices[1].size();


   int nPoints = -1;
   std::cout<<"nPUpoints "<<nPUpoints<<" nLumiPoints "<<nLumiPoints<<std::endl;


   vector<vector<vector<TEfficiency* > > > Eff_phiMBWh;
   vector<vector<vector<TEfficiency* > > > EffA_phiMBWh;

   vector<vector<vector<TEfficiency* > > > Eff_theMBWh;
   vector<vector<vector<TEfficiency* > > > EffA_theMBWh;

   vector<vector<TEfficiency* > >          Eff_phiMB4Top;
   vector<vector<TEfficiency* > >          EffA_phiMB4Top;

   vector<TEfficiency* >                   Eff_phiMB4Bot;
   vector<TEfficiency* >                   EffA_phiMB4Bot;


   vector<vector<vector<TH2F* > > > Hist_MBWh;
   vector<vector<TH2F* > >          Hist_MB4Top;
   vector<TH2F* >                   Hist_MB4Bot;

   vector<vector<vector<TProfile* > > > Gr_MBWh;
   vector<vector<TProfile* > >          Gr_MB4Top;
   vector<TProfile* >                   Gr_MB4Bot;


   
   for (int ivar=0; ivar<nVar; ivar++){
      
     Eff_phiMBWh.push_back(  vector<vector<TEfficiency*> > ());  
     EffA_phiMBWh.push_back( vector<vector<TEfficiency*> > ());  
     
     Eff_theMBWh.push_back(  vector<vector<TEfficiency*> > ());  
     EffA_theMBWh.push_back( vector<vector<TEfficiency*> > ());  
     
     Eff_phiMB4Top.push_back(  vector<TEfficiency*>  ());  
     EffA_phiMB4Top.push_back( vector<TEfficiency*>  ());  

     Hist_MBWh.push_back(  vector<vector<TH2F*> > ());  
     Hist_MB4Top.push_back(  vector<TH2F*>  ());  

     Gr_MBWh.push_back(  vector<vector<TProfile*> > ());  
     Gr_MB4Top.push_back(  vector<TProfile*>  ());  
    
     nPoints = slices[ivar].size();

     for (int iwh=0; iwh<5; iwh++){

       Eff_phiMBWh[ivar].push_back(vector<TEfficiency*> ());  
       EffA_phiMBWh[ivar].push_back(vector<TEfficiency*> ());  
       
       Eff_theMBWh[ivar].push_back(vector<TEfficiency*> ());         
       EffA_theMBWh[ivar].push_back(vector<TEfficiency*> ());         

       Hist_MBWh[ivar].push_back(vector<TH2F*> ());         
       Gr_MBWh[ivar].push_back(vector<TProfile*> ());         


       for (int ist=0; ist<4; ist++){

	 Eff_phiMBWh[ivar][iwh].push_back( new TEfficiency("","",nPoints-1,slices[ivar].data()));	 	 
	 EffA_phiMBWh[ivar][iwh].push_back( new TEfficiency("","",nPoints-1,slices[ivar].data()));	 	 
	 Hist_MBWh[ivar][iwh].push_back(new TH2F(("Hist"+varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
						  "St"+std::to_string(ist)).c_str(),"",slices[ivar].size()-1,
						 slices[ivar].data(),slices[2].size()-1,slices[2].data()));	 	 
	 Gr_MBWh[ivar][iwh].push_back( new TProfile());//"Hist"+varName[ivar]+"_MBWh"+std::to_string(iwh-2)+
	 //"St"+std::to_string(ist)).c_str()));
	 if (ist!=3){ 
	   Eff_theMBWh[ivar][iwh].push_back( new TEfficiency("","",nPoints-1,slices[ivar].data()));	 
	   EffA_theMBWh[ivar][iwh].push_back( new TEfficiency("","",nPoints-1,slices[ivar].data()));	 
	 }
       }
       Eff_phiMB4Top[ivar].push_back( new TEfficiency("","",nPoints-1,slices[ivar].data()));	          
       EffA_phiMB4Top[ivar].push_back( new TEfficiency("","",nPoints-1,slices[ivar].data()));	          
       Hist_MB4Top[ivar].push_back( new TH2F(("Hist"+varName[ivar]+"_MBTopWh"+std::to_string(iwh-2)).c_str(),"",slices[ivar].size()-1,slices[ivar].data(),slices[2].size()-1,slices[2].data()));
       Gr_MB4Top[ivar].push_back( new TProfile());//"Hist"+varName[ivar]+"_MBTopWh"+std::to_string(iwh-2)).c_str()));
     }
     Eff_phiMB4Bot.push_back( new TEfficiency("","",nPoints-1,slices[ivar].data()));	      
     EffA_phiMB4Bot.push_back( new TEfficiency("","",nPoints-1,slices[ivar].data()));	      
     Hist_MB4Bot.push_back( new TH2F(("Hist"+varName[ivar]+"_MBBot").c_str(),"",slices[ivar].size()-1,slices[ivar].data(),slices[2].size()-1,slices[2].data()));	      
     Gr_MB4Bot.push_back( new TProfile());//"Hist"+varName[ivar]+"_MBBot").c_str()));
   }

   
   // DEAD CHANNELS (to skip):
   
   cout<<" within Loop: Ndead "<<Ndead<<endl;
   
   if (Ndead==0) {
     cout<<" READING FROM FILE !! "<<endl;
     ifstream txtin;
     std::string deadname;
     
     
     if(fileName!="")    deadname.append("DeadList_"+fileName);
     else{
       deadname.append("DeadList_Run2016");
       deadname = deadname + dataset + ".txt";;
     }
     
     cout<<" reading from "<<deadname<<endl;
     txtin.open(deadname.c_str(),std::ifstream::in);
     
     int idead    = 0;
     int prevlay  = 0; int prevSL = 0;
     int ndeadlay = 0;
     
     do {
       txtin>>idead>>dead[Ndead][0]>>dead[Ndead][1]>>dead[Ndead][2]>>dead[Ndead][3]>>dead[Ndead][4]>>dead[Ndead][5];
       
       if (dead[Ndead][3] == prevSL && dead[Ndead][4] == prevlay ) {
	 ndeadlay++;
       }
       else {
         if (ndeadlay>20) cout<<"YB"<<dead[Ndead][0]<<"/Sect"<<dead[Ndead][1]<<"/MB"<<dead[Ndead][2]<<" SL"
			      <<dead[Ndead][3]<<" L"
			      <<dead[Ndead][4]<<" "<<ndeadlay+1<<" deads!!"<<endl;
         prevSL=dead[Ndead][3]; 
         prevlay=dead[Ndead][4]; 
         ndeadlay=0;
       }
       Ndead++;
     }
     while (idead!=0);
   }
   //   cout<<"MAX "<<fChain->GetMaximum("lumiblock")<<std::endl;; //del

   std::vector<Float_t>  varVal  = {0,0}; 

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     
     //     std::cout<<"entry "<<jentry<<std::endl;

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;
     
     if (jentry%10000 == 0) cout<<"Evento "<<jentry<<endl;
     
     varVal[0] = lumiperblock;
     varVal[1] = PV_Nvtx;


     if (lumiperblock > slices[0].back()) { cout<<" luminosity out of range!! "<< lumiperblock<<endl; continue; }
     if (PV_Nvtx > slices[1].back())  { cout<<" PU out of range!! "        << PV_Nvtx<<endl;      continue; }
     
     // First search for Phi segments

     //cout<<"Ndtsegments "<<Ndtsegments<<std::endl;     
     Float_t bkg = 0;
     for (int iseg=0; iseg<Ndtsegments; iseg++) {
       bkg=0;
       //selection
       if (!dtsegm4D_hasPhi->at(iseg)) continue;
       // In chambers 1,2,3 select only segments with also Z (theta) contribution.
       if (dtsegm4D_station->at(iseg)!=4 && !dtsegm4D_hasZed->at(iseg)) continue;
       
       int seg_phinhits = dtsegm4D_phinhits->at(iseg);
       
       if (fabs(dtsegm4D_x_dir_loc->at(iseg))>0.7) continue; // angle
       
       TVectorF *expWire=(TVectorF*)dtsegm4D_hitsExpWire->At(iseg);
       //	expWire->Print(); //del
       
       // If a hit is missing, let us check that the extrapolation doesn't fall beyond layer or cross a dead cell!
       int NexpDead=0; bool OutOfLayer=false;
       

       if (seg_phinhits < 8 ) {
         for (int iex=0; iex<12; iex++) {

	   int expSL = 1;
           int expLay = iex+1;
	   //associate layer with right super layer
           if (dtsegm4D_station->at(iseg) != 4){
	     if (iex > 3 && iex<8) {expSL=2; expLay-=4;}
             else if (iex>7) {expSL=3; expLay-=8;}
	   }
           else {
             if (iex > 3 && iex<8) continue;
             else if (iex > 7) {expSL=3; expLay-=8;}
	   }
	   int nwire=0;
	   if (expSL==2){
	     //chimney chambers
	     if (( dtsegm4D_wheel->at(iseg)==-1 && dtsegm4D_station->at(iseg)==3 ) || ( dtsegm4D_wheel->at(iseg)==1 && dtsegm4D_station->at(iseg)==4 ) ) nwire = 48;
	     else nwire=57;
	   }
	   else if (dtsegm4D_station->at(iseg)==1) nwire = 49;
           else if (dtsegm4D_station->at(iseg)==2) nwire = 60;
           else if (dtsegm4D_station->at(iseg)==3) nwire = 72;
           else if (dtsegm4D_station->at(iseg)==4) {
	     if (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)        nwire = 72;
             else if( dtsegm4D_sector->at(iseg)==8  || dtsegm4D_sector->at(iseg)==12 ) nwire = 92;
             else if( dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14 ) nwire = 60;
	     else if( dtsegm4D_sector->at(iseg)==9  || dtsegm4D_sector->at(iseg)==11 ) nwire = 48;
             else nwire = 96;
	   }
	   
           Float_t expW = (*expWire)(iex);
	   //	   std::cout<<"iex "<<iex<<" explay "<<expLay<<std::endl; //del
           if (expW>nwire) {
	     OutOfLayer=true;
	     //	     std::cout<<"Out of layer "<<expW<<" nwire "<<nwire<<" "<<dtsegm4D_wheel->at(iseg)<<" "<<dtsegm4D_sector->at(iseg)<<" "<<dtsegm4D_station->at(iseg)<<" "<<expSL<<" "<<expLay<<" "<<expW<<std::endl; //del
             break;
	   }

           for (int idead=0; idead<Ndead; idead++) {
	     if (dead[idead][0] != dtsegm4D_wheel->at(iseg))   continue;
	     if (dead[idead][1] != dtsegm4D_sector->at(iseg))  continue;
	     if (dead[idead][2] != dtsegm4D_station->at(iseg)) continue;
	     if (dead[idead][3] != expSL)  continue;  
	     if (dead[idead][4] != expLay) continue;
	     if (dead[idead][5] != expW)   continue; 
	     NexpDead++;
	     //	    std::cout<<"dead "<<dtsegm4D_wheel->at(iseg)<<" "<<dtsegm4D_sector->at(iseg)<<" "<<dtsegm4D_station->at(iseg)<<" "<<expSL<<" "<<expLay<<" "<<expW<<std::endl; //del
	     break;
	   }
           if (NexpDead>MaxDead) break;
	 }
         if (OutOfLayer)       continue; // this segment goes out of layer boundary
         if (NexpDead>MaxDead) continue; // this segment crosses dead cell(s): drop it!
       }
       
       int NHits=0; int missingLayer[3][2]; for (int imi=0; imi<3;imi++) {missingLayer[imi][0]=0;missingLayer[imi][1]=0;}
       int nmissing=0;
       
       TVectorF *hitSuperLayerPhi =(TVectorF*)dtsegm4D_phi_hitsSuperLayer->At(iseg);
       TVectorF *hitLayerPhi      =(TVectorF*)dtsegm4D_phi_hitsLayer->At(iseg);
       TVectorF *hitWirePhi       =(TVectorF*)dtsegm4D_phi_hitsWire->At(iseg);
       
       // std::cout<<"SuperLayerphi"<<std::endl;
       // hitSuperLayerPhi->Print();//del
       // std::cout<<"Layerphi"<<std::endl;
       // hitLayerPhi->Print(); //del
       // std::cout<<"Wirephi"<<std::endl;
       // hitWirePhi->Print(); //del
       
       for (int ilay=1; ilay<9; ilay++) {
	 // Search for associated hits
	 bool foundh = false;
	 for (int kk=0; kk<seg_phinhits; kk++) {
           int sl1  = (*hitSuperLayerPhi)(kk);
           int lay1 = (sl1==1) ? (*hitLayerPhi)(kk) : (*hitLayerPhi)(kk)+4; // hit layer 1-8
	   
           if (lay1==ilay) {
	     NHits++;
             foundh=true;
             break;
           }
	 }

	 if (!foundh) {
	   if (nmissing<3) missingLayer[nmissing][0]=ilay;
	   nmissing++;
	 }

       }
       if (nmissing != 8-NHits) {cout<<NHits<<" hits and "<<nmissing<<" missing!!"<<endl; return;}

       //bgk

       for (int idigi=0; idigi<Ndigis; idigi++) {                                                                    
	 if ((digi_time->at(idigi)>320 || digi_time->at(idigi)<700) &&
	     (digi_wheel->at(idigi)   == dtsegm4D_wheel->at(iseg)) &&
	     (digi_sector->at(idigi)  == dtsegm4D_sector->at(iseg)) &&
	     (digi_station->at(idigi) == dtsegm4D_station->at(iseg))) bkg++;
	     }
       
       varVal[2] = bkg;
       
       Hist_MBWh[0][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(varVal[0],bkg); 
       Hist_MBWh[1][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(varVal[1],bkg); 
       
       if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
	 Hist_MB4Top[0][dtsegm4D_wheel->at(iseg)+2]->Fill(varVal[0],bkg); 
	 Hist_MB4Top[1][dtsegm4D_wheel->at(iseg)+2]->Fill(varVal[1],bkg); 
       }
	     
       // extra chamber of sector 10 (sector 14) 
       else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
	 Hist_MB4Bot[0]->Fill(varVal[0],bkg); 
	 Hist_MB4Bot[1]->Fill(varVal[1],bkg); 
       }       
       
       if (NHits<nrequiredhit) continue;
       else if (NHits==8) {
	 for (int sl=0; sl<2; sl++) for (int lay=0; lay<4; lay++) {
	     //variables, points, stations, wheels 
	     
	     for (int ivar=0; ivar<nVar; ivar++) { 
	       Eff_phiMBWh[ivar][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(1,varVal[ivar]); 
	       EffA_phiMBWh[ivar][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(1,varVal[ivar]); 
	       
	       // extra chamber of sector 4 (sector 13)
	       if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
		 Eff_phiMB4Top[ivar][dtsegm4D_wheel->at(iseg)+2]->Fill(1,varVal[ivar]); 
	       }
	       
	       // extra chamber of sector 10 (sector 14) 
	       else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
		 Eff_phiMB4Bot[ivar]->Fill(1,varVal[ivar]); 
		 
	       }
	     }
	   }
       }
       else { // let's see how to treat missing layers
	 for (int imiss=0; imiss<nmissing; imiss++) {
           int sl  = missingLayer[imiss][0] < 5 ? 0 : 1;
           int lay = sl==0 ? missingLayer[imiss][0]-1 : missingLayer[imiss][0]-5;
	   
	   // is there a digi within the expected tube?
	   
           float digiW  = -1.;
           float d      = 1000000; //just a very big number
           int iex      = missingLayer[imiss][0] < 5 ? missingLayer[imiss][0]-1 : missingLayer[imiss][0]+3;
           Float_t expW = (*expWire)(iex);  //perche` float?
	   
           for (int idigi=0; idigi<Ndigis; idigi++) {
	     
             if (digi_time->at(idigi)<320 || digi_time->at(idigi)>700)  continue; //require only digis time inside time box
 	     if (digi_wheel->at(idigi)   != dtsegm4D_wheel->at(iseg))   continue;
	     if (digi_sector->at(idigi)  != dtsegm4D_sector->at(iseg))  continue;
	     if (digi_station->at(idigi) != dtsegm4D_station->at(iseg)) continue;
	     
	     if (digi_sl->at(idigi) == 2)            continue;
	     if (digi_sl->at(idigi) == 1 && sl != 0) continue;
	     if (digi_sl->at(idigi) == 3 && sl != 1) continue;
	     if (digi_layer->at(idigi) != lay+1)     continue;
	     //let's loop all over the digis and take the closest digis to the extrapolated wire. 
	     // Think about an extra condition on time.
             if (fabs(expW-digi_wire->at(idigi))<fabs(d)) {
	       // std::cout<<expW<< " "<<digi_wire->at(idigi)<<" "<<expW-digi_wire->at(idigi)<<std::endl;
               digiW=digi_wire->at(idigi);
	       d=expW-digiW; 
	     }
	   }
           if ( fabs(d)< 1.1) {missingLayer[imiss][1]=1; } //non dovrebbe essere sempre un intero d?
	 }
	 
	 if (NHits==nrequiredhit) {
	   for (int imiss=0; imiss<nmissing; imiss++) {
	     int sl = missingLayer[imiss][0] < 5 ? 0 : 1;
	     int lay = sl==0 ? missingLayer[imiss][0]-1 : missingLayer[imiss][0]-5;
	     
	     for (int ivar=0; ivar<nVar; ivar++) { 
	       
	       Eff_phiMBWh[ivar][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(missingLayer[imiss][1],varVal[ivar]); 
	       if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
		 Eff_phiMB4Top[ivar][dtsegm4D_wheel->at(iseg)+2]->Fill(missingLayer[imiss][1],varVal[ivar]); 

	       }
	       else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {		
		 Eff_phiMB4Bot[ivar]->Fill(missingLayer[imiss][1],varVal[ivar]);

	       }
	     }
	   }
	 }
	 
	 else {

	   for (int sl=0; sl<2; sl++) for (int lay=0; lay<4; lay++) {
	       bool missAss=false; bool missDigi=false;
               for (int imiss=0; imiss<nmissing; imiss++) {
                 if (missingLayer[imiss][0]==sl*4+lay+1) {
		   missAss=true;
		   if (!missingLayer[imiss][1]) missDigi=true;
		 }
	       }
	       for (int ivar=0; ivar<nVar; ivar++) { 
		 Eff_phiMBWh[ivar][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(!(missAss&&missDigi),varVal[ivar]); 
		 EffA_phiMBWh[ivar][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(!(missAss),varVal[ivar]); 
	       
	       if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
		 Eff_phiMB4Top[ivar][dtsegm4D_wheel->at(iseg)+2]->Fill(!(missAss&&missDigi),varVal[ivar]); 
		 EffA_phiMB4Top[ivar][dtsegm4D_wheel->at(iseg)+2]->Fill(!(missAss),varVal[ivar]); 
		 
	       }
	       else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {		
		 Eff_phiMB4Bot[ivar]->Fill(!(missAss&&missDigi),varVal[ivar]);
		 EffA_phiMB4Bot[ivar]->Fill(!(missAss),varVal[ivar]);
	       }
	       }
	     }
	 }
       }
     }
     // Then search for Zed segments
     
     for (int iseg=0; iseg<Ndtsegments; iseg++) {
       //selection
       //if (!dtsegm4D_hasZed->at(iseg) || !dtsegm4D_hasPhi->at(iseg)) continue; 
       if (!dtsegm4D_hasZed->at(iseg) || dtsegm4D_phinhits->at(iseg)<nrequiredhit) continue; 
       //10-01-17 capire perchè theta efficiency va giù con l ainst.luminosity: è il denominatore che sale?!

       int seg_znhits = dtsegm4D_znhits->at(iseg);
       
       //if (fabs(dtsegm4D_y_dir_loc->at(iseg))>0.7) continue; // angle WARNING!!! try and disable this for theta layers!
       if (seg_znhits < 3) continue; // piuttosto ovvio!!!  :-)
       
       TVectorF *expWire=(TVectorF*)dtsegm4D_hitsExpWire->At(iseg);
       //	expWire->Print();
       
       // If a hit is missing, let us check that the extrapolation doesn't fall out of layer or cross a dead cell!
       int NexpDead=0; bool OutOfLayer=false;
       
       if (seg_znhits < 4) {
         for (int iex=4; iex<8; iex++) {
	   int expSL = 2;
           int expLay = iex-3;		
	   
           Float_t expW = (*expWire)(iex);
           if (expW>58) {
	     OutOfLayer=true;
             break;
	   }
           for (int idead=0; idead<Ndead; idead++) {
	     if (dead[idead][0] != dtsegm4D_wheel->at(iseg)) continue;
	     if (dead[idead][1] != dtsegm4D_sector->at(iseg)) continue;
	     if (dead[idead][2] != dtsegm4D_station->at(iseg)) continue;
	     if (dead[idead][3] != expSL) continue;  
	     if (dead[idead][4] != expLay) continue;
	     if (dead[idead][5] != expW) continue; 
	     NexpDead++;
	     break;   
	   }
           if (OutOfLayer) break;
           if (NexpDead>MaxDead) break; 
	 }
	 
         if (NexpDead>MaxDead) {
           continue;
	   // this segment crosses at least 1 dead cell: drop it!
	 }
       }
       
       int NHits=0; int missingLayer=-1;
       
       TVectorF *hitLayerZ = (TVectorF*)dtsegm4D_z_hitsLayer->At(iseg);
       TVectorF *hitWireZ  = (TVectorF*)dtsegm4D_z_hitsWire->At(iseg);
       
       for (int ilay=1; ilay<5; ilay++) {
	 // Search for associated hits
	 bool foundh = false;
	 for (int kk=0; kk<seg_znhits; kk++) {
	   
           int lay1 = (*hitLayerZ)(kk);
	   
           if (lay1==ilay) {
	     NHits++;
             foundh=true;
             break;
           }
	 }
	 if (!foundh) missingLayer=ilay;
       }
       if (NHits<3) continue;
       else if (NHits==4) {
	 
	 for (int lay=0; lay<4; lay++) {
	   for (int ivar=0; ivar<nVar; ivar++) {
	   Eff_theMBWh[ivar][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(1,varVal[ivar]); 
	   EffA_theMBWh[ivar][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(1,varVal[ivar]); 

	   }
	 }
       }
       else if (NHits==3) {
	 int lay = missingLayer-1;
	 
	 // is there a digi within the expected tube?
	 
	 float digiW=-1.;
	 float d =1000000;
	 int iex = missingLayer+3;
	 Float_t expW = (*expWire)(iex);
	 
	 for (int idigi=0; idigi<Ndigis; idigi++) {
	   
	   if (digi_time->at(idigi)<320 || digi_time->at(idigi)>700)  continue;
	   if (digi_wheel->at(idigi)   != dtsegm4D_wheel->at(iseg))   continue;
	   if (digi_sector->at(idigi)  != dtsegm4D_sector->at(iseg))  continue;
	   if (digi_station->at(idigi) != dtsegm4D_station->at(iseg)) continue;
	   
	   if (digi_sl->at(idigi) != 2) continue;
	   if (digi_layer->at(idigi) != lay+1) continue;
           
	   if (fabs(expW-digi_wire->at(idigi))<fabs(d)) {
	     digiW=digi_wire->at(idigi);
	     d=expW-digiW;
	   }
	 }
	 for (int ivar=0; ivar<nVar; ivar++) Eff_theMBWh[ivar][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill((fabs(d)< 1.1) ,varVal[ivar]); 
	 }
       else {
	 cout<<" what do you want?? NHits (z) = "<<NHits<<endl;
	 return;
       }
     }
   }
   
   ofstream results;
   std::string resultname;
   
   resultname.append("Results_"+fileName+".txt");
   //   resultname = resultname + dataset + ".txt";;
   
   cout<<resultname<<endl;
   results.open (resultname.c_str());
   
   
   //to be fixed
   for (int ivar=0; ivar<nVar; ivar++) {
     for (int iwh=0; iwh<5; iwh++){
       for (int ist=0; ist<4; ist++){
	 
	 nPoints = slices[ivar].size();
	 for (int ipoint=0; ipoint<nPoints; ipoint++) {
	   
	   results<<slices[ivar][ipoint];
	   results <<" "<<ist+1<<" "<<iwh-2
		   <<" "<<Eff_phiMBWh[ivar][iwh][ist]->GetTotalHistogram()->GetBinContent(ipoint+1)<<" "<<
	     Eff_phiMBWh[ivar][iwh][ist]->GetPassedHistogram()->GetBinContent(ipoint+1)
		   <<" "<<EffA_phiMBWh[ivar][iwh][ist]->GetTotalHistogram()->GetBinContent(ipoint+1)<<endl;
	   
	   if(ist<3){	 
	     results<<slices[ivar][ipoint];
	     results <<" "<<ist+1<<" "<<iwh-2
		     <<" "<<Eff_theMBWh[ivar][iwh][ist]->GetTotalHistogram()->GetBinContent(ipoint+1)<<" "<<
	       Eff_theMBWh[ivar][iwh][ist]->GetPassedHistogram()->GetBinContent(ipoint+1)
		     <<" "<<EffA_theMBWh[ivar][iwh][ist]->GetTotalHistogram()->GetBinContent(ipoint+1)<<endl;
	   }
	 }
       }      
       for (int ipoint=0; ipoint<nPoints; ipoint++) {
	 results<<slices[ivar][ipoint];      
	 results <<" 4T"<<" "<<iwh-2
		 <<" "<<Eff_phiMB4Top[ivar][iwh]->GetTotalHistogram()->GetBinContent(ipoint+1)<<" "<<
	   Eff_phiMB4Top[ivar][iwh]->GetPassedHistogram()->GetBinContent(ipoint+1)
		 <<" "<<EffA_phiMB4Top[ivar][iwh]->GetTotalHistogram()->GetBinContent(ipoint+1)<<endl;
       }
     }
     
     for (int ipoint=0; ipoint<nPoints; ipoint++) {
       
       results<<slices[ivar][ipoint];  
       results <<" 4B "
	       <<" "<<Eff_phiMB4Bot[ivar]->GetTotalHistogram()->GetBinContent(ipoint+1)<<" "<<
	 Eff_phiMB4Bot[ivar]->GetPassedHistogram()->GetBinContent(ipoint+1)
	       <<" "<<EffA_phiMB4Bot[ivar]->GetTotalHistogram()->GetBinContent(ipoint+1)<<endl;
     }
   }
      
   system("mkdir plot/");
   system(("mkdir plot/"+fileName).c_str());
   
   system(("mkdir "+WebFolder+"/"+fileName).c_str());
   system(("cp "+WebFolder+"/index.php " +WebFolder+"/"+fileName).c_str());

   system(("mkdir "+WebFolder+"/"+fileName+"/Efficiency").c_str());
   system(("cp "+WebFolder+"/index.php " +WebFolder+"/"+fileName+"/Efficiency").c_str());

   system(("mkdir "+WebFolder+"/"+fileName+"/Background").c_str());
   system(("cp "+WebFolder+"/index.php " +WebFolder+"/"+fileName+"/Background").c_str());


   //   cout<<"cp "+WebFolder+"/index.php " +WebFolder+"/"+fileName<<endl;
   for (int ivar=0; ivar<nVar; ivar++) {
     for (int iwh=0; iwh<5; iwh++){
       for (int ist=0; ist<4; ist++){

	 if(iwh!=4){
	   Eff_phiMBWh[ivar][iwh][ist]->SetLineColor(iwh+1);
	   Eff_phiMBWh[ivar][iwh][ist]->SetMarkerColor(iwh+1);

	   Hist_MBWh[ivar][iwh][ist]->SetLineColor(iwh+1);
	   Hist_MBWh[ivar][iwh][ist]->SetMarkerColor(iwh+1);

	   if(ist!=3){
	     Eff_theMBWh[ivar][iwh][ist]->SetLineColor(iwh+1);
	     Eff_theMBWh[ivar][iwh][ist]->SetMarkerColor(iwh+1);
	   }
	 }
	 else{ 
	   Eff_phiMBWh[ivar][iwh][ist]->SetLineColor(6);
	   Eff_phiMBWh[ivar][iwh][ist]->SetMarkerColor(6);

	   Hist_MBWh[ivar][iwh][ist]->SetLineColor(6);
	   Hist_MBWh[ivar][iwh][ist]->SetMarkerColor(6);

	   if(ist!=3){
	     Eff_theMBWh[ivar][iwh][ist]->SetLineColor(6);
	     Eff_theMBWh[ivar][iwh][ist]->SetMarkerColor(6);
	   }
	 }
	 Eff_phiMBWh[ivar][iwh][ist]->SetMarkerStyle(20);
	 Hist_MBWh[ivar][iwh][ist]->SetMarkerStyle(20);
	 //profile of 2D histograms
	 Gr_MBWh[ivar][iwh][ist] = Hist_MBWh[ivar][iwh][ist]->ProfileX();

	 if(ist!=3)  Eff_theMBWh[ivar][iwh][ist]->SetMarkerStyle(20);
       }
       Eff_phiMB4Top[ivar][iwh]->SetMarkerStyle(20);
       Hist_MB4Top[ivar][iwh]->SetMarkerStyle(20);
       Gr_MB4Top[ivar][iwh] = Hist_MB4Top[ivar][iwh]->ProfileX();
       if(iwh!=4){
	 Eff_phiMB4Top[ivar][iwh]->SetLineColor(iwh+1);
	 Eff_phiMB4Top[ivar][iwh]->SetMarkerColor(iwh+1);

	 Hist_MB4Top[ivar][iwh]->SetLineColor(iwh+1);
	 Hist_MB4Top[ivar][iwh]->SetMarkerColor(iwh+1);
       }
       else {
	 Eff_phiMB4Top[ivar][iwh]->SetLineColor(6);
	 Eff_phiMB4Top[ivar][iwh]->SetMarkerColor(6);

	 Hist_MB4Top[ivar][iwh]->SetLineColor(6);
	 Hist_MB4Top[ivar][iwh]->SetMarkerColor(6);
       }
     }
     Eff_phiMB4Bot[ivar]->SetMarkerStyle(20);
     Eff_phiMB4Bot[ivar]->SetLineColor(1);
     Eff_phiMB4Bot[ivar]->SetMarkerColor(1);
     Hist_MB4Bot[ivar]->SetMarkerStyle(20);
     Hist_MB4Bot[ivar]->SetLineColor(1);
     Hist_MB4Bot[ivar]->SetMarkerColor(1);
     Gr_MB4Bot[ivar] = Hist_MB4Bot[ivar]->ProfileX();
   }
   
   
   for (int ivar=0; ivar<nVar; ivar++){
     
     for (int ist=0; ist<4; ist++){
       //phi
       TCanvas *cPhiMB = new TCanvas(("cPhiMB"+(std::to_string(ist+1))+varName[ivar]).c_str());
       
       //       double xmin = TMath::MinElement(g->GetN(),g->GetX()); 

       // int binmax =Eff_phiMBWh[ivar][0][ist]->GetTotalHistogram()->GetMaximumBin();
       // double xMax = Eff_phiMBWh[ivar][0][ist]->GetTotalHistogram()->GetXaxis()->GetBinCenter(binmax);
       // xMax+=xMax*0.1;

       // int binmin =Eff_phiMBWh[ivar][0][ist]->GetTotalHistogram()->GetMinimumBin();
       // double xMin = Eff_phiMBWh[ivar][0][ist]->GetTotalHistogram()->GetXaxis()->GetBinCenter(binmin);
       // xMin+=xMin*0.1;

       // cout<<"xmin "<<xMin<<" xmax "<<xMax<<endl;
       Eff_phiMBWh[ivar][0][ist]->SetTitle(("MB"+(std::to_string(ist+1))+" eff vs "+varTitle[ivar]).c_str());
       Eff_phiMBWh[ivar][0][ist]->Draw("ap");
       
       gPad->Update(); 
       auto graph = Eff_phiMBWh[ivar][0][ist]->GetPaintedGraph(); 
       //       graph->GetXaxis()->SetRangeUser(xMin,xMax);
       graph->SetMinimum(0.89);
       graph->SetMaximum(1.02); 
       gPad->Update(); 
       
       TLegend * legPhiMB = new TLegend(0.75,0.75,0.9,0.9);
       for (int iwh=0; iwh<5; iwh++){
	 Eff_phiMBWh[ivar][iwh][0]->Draw("samep");
	 legPhiMB->AddEntry(Eff_phiMBWh[ivar][iwh][0],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
       }
       legPhiMB->Draw("same");
       cPhiMB->SaveAs((WebFolder+"/"+fileName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"PhiEffVs"+varName[ivar]+".png").c_str());
       
       
       //theta
       TCanvas *cTheMB = new TCanvas(("cPhiMB"+(std::to_string(ist+1))+varName[ivar]).c_str());
       
       Eff_phiMBWh[ivar][0][ist]->SetTitle(("MB"+(std::to_string(ist+1))+" eff vs "+varTitle[ivar]).c_str());
       Eff_phiMBWh[ivar][0][ist]->Draw("ap");
       
       gPad->Update(); 
       graph = Eff_phiMBWh[ivar][0][ist]->GetPaintedGraph(); 
       graph->SetMinimum(0.89);
       graph->SetMaximum(1.02); 
       gPad->Update(); 
       
       TLegend * legTheMB = new TLegend(0.75,0.75,0.9,0.9);
       for (int iwh=0; iwh<5; iwh++){
	 Eff_phiMBWh[ivar][iwh][0]->Draw("samep");
	 legTheMB->AddEntry(Eff_phiMBWh[ivar][iwh][0],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
       }
       legTheMB->Draw("same");
       cTheMB->SaveAs((WebFolder+"/"+fileName+"/Efficiency/"+"MB"+(std::to_string(ist+1))+"TheEffVs"+varName[ivar]+".png").c_str());
       
     }
     
     //phi MB4 Top
     TCanvas *cPhiMB4Top = new TCanvas(("cPhiMB4Top"+varName[ivar]).c_str());
     
     Eff_phiMB4Top[ivar][0]->SetTitle(("MB4Top eff vs "+varTitle[ivar]).c_str());
     Eff_phiMB4Top[ivar][0]->Draw("ap");
     
     gPad->Update(); 
     auto graph = Eff_phiMB4Top[ivar][0]->GetPaintedGraph(); 
     graph->SetMinimum(0.89);
     graph->SetMaximum(1.02); 
     gPad->Update(); 
     
     TLegend * legPhiMB4Top = new TLegend(0.75,0.75,0.9,0.9);
     for (int iwh=0; iwh<5; iwh++){
       Eff_phiMB4Top[ivar][iwh]->Draw("samep");
       legPhiMB4Top->AddEntry(Eff_phiMB4Top[ivar][iwh],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
     }
     legPhiMB4Top->Draw("same");
     cPhiMB4Top->SaveAs((WebFolder+"/"+fileName+"/Efficiency/"+"MB4TopPhiEffVs"+varName[ivar]+".png").c_str());
     
     
     //phi MB4 Top
     TCanvas *cPhiMB4Bot = new TCanvas(("cPhiMB4Bot"+varName[ivar]).c_str());
     
     Eff_phiMB4Bot[ivar]->SetTitle(("MB4Bot eff vs "+varTitle[ivar]).c_str());
     Eff_phiMB4Bot[ivar]->Draw("ap");
     
     gPad->Update(); 
     graph = Eff_phiMB4Bot[ivar]->GetPaintedGraph(); 
     graph->SetMinimum(0.89);
     graph->SetMaximum(1.02); 
     gPad->Update(); 
     
     cPhiMB4Bot->SaveAs((WebFolder+"/"+fileName+"/Efficiency/"+"MB4BotPhiEffVs"+varName[ivar]+".png").c_str());
   }
   
   //Graph bkg
   
   for (int ivar=0; ivar<2; ivar++){
     
     for (int ist=0; ist<4; ist++){
       //phi
       TCanvas *cMB = new TCanvas(("cMB"+(std::to_string(ist+1))+varName[ivar]).c_str());

       
       Gr_MBWh[ivar][0][ist]->SetMaximum(35);
       Gr_MBWh[ivar][0][ist]->SetTitle(("MB"+(std::to_string(ist+1))+" bkg vs "+varTitle[ivar]).c_str());
       Gr_MBWh[ivar][0][ist]->Draw("E1");
       
       
       TLegend * legMB = new TLegend(0.75,0.75,0.9,0.9);
       for (int iwh=0; iwh<5; iwh++){
	 Gr_MBWh[ivar][iwh][0]->Draw("sameE1");
	 legMB->AddEntry(Gr_MBWh[ivar][iwh][0],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
       }
       legMB->Draw("same");
       cMB->SaveAs((WebFolder+"/"+fileName+"/Background/"+"MB"+(std::to_string(ist+1))+"BkgVs"+varName[ivar]+".png").c_str());
       
     }
     // MB4 Top
     TCanvas *cMB4Top = new TCanvas(("cMB4Top"+varName[ivar]).c_str());
     Gr_MB4Top[ivar][0]->SetMaximum(35);
     Gr_MB4Top[ivar][0]->SetTitle(("MB4Top bkg vs "+varTitle[ivar]).c_str());
     Gr_MB4Top[ivar][0]->Draw("E1");
          
     TLegend * legMB4Top = new TLegend(0.75,0.75,0.9,0.9);
     for (int iwh=0; iwh<5; iwh++){
       Gr_MB4Top[ivar][iwh]->Draw("sameE1");
       legMB4Top->AddEntry(Gr_MB4Top[ivar][iwh],("Wh"+std::to_string(iwh-2)).c_str(),"lpe");
     }
     legMB4Top->Draw("same");
     cMB4Top->SaveAs((WebFolder+"/"+fileName+"/Background/"+"MB4TopBkgVs"+varName[ivar]+".png").c_str());
     
     
     // MB4 Top
     TCanvas *cMB4Bot = new TCanvas(("cMB4Bot"+varName[ivar]).c_str());
     
     Gr_MB4Bot[ivar]->SetTitle(("MB4Bot bkg vs "+varTitle[ivar]).c_str());
     Gr_MB4Bot[ivar]->Draw("E1");
    
     cMB4Bot->SaveAs((WebFolder+"/"+fileName+"/Background/"+"MB4BotBkgVs"+varName[ivar]+".png").c_str());
   }      
     
}


void EfficiencyMonitor::PostLoop()
{
  // check how many cells were not crossed by selected segments (efficiency undetermined)


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   //   nentries=lessentries;

   Long64_t nbytes = 0, nb = 0;
   char go;

   ofstream UndefList;
   UndefList.open ("UndefList5Hits.txt");

   int undef[10000][6];
   int Nundef=0;

   if (Ndead==0) {
     cout<<" READING FROM FILE !! "<<endl;
     ifstream txtin;
     std::string deadname;
     if(fileName!="") deadname.append("DeadList_"+fileName);
     else{
       deadname.append("DeadList_Run2016");
       deadname = deadname + dataset + ".txt";;
     }
     txtin.open(deadname+".txt",std::ifstream::in);
     int idead=0;
     int prevlay=0; int prevSL=0;
     int ndeadlay=0;
     do {
       txtin>>idead>>dead[Ndead][0]>>dead[Ndead][1]>>dead[Ndead][2]>>dead[Ndead][3]>>dead[Ndead][4]>>dead[Ndead][5];


       if (dead[Ndead][3] == prevSL && dead[Ndead][4] == prevlay ) {
	 ndeadlay++;
       }
       else {
         if (ndeadlay>20) cout<<"YB"<<dead[Ndead][0]<<"/Sect"<<dead[Ndead][1]<<"/MB"<<dead[Ndead][2]<<" SL"
                             <<dead[Ndead][3]<<" L"
                             <<dead[Ndead][4]<<" "<<ndeadlay+1<<" deads!!"<<endl;
         prevSL=dead[Ndead][3]; 
         prevlay=dead[Ndead][4]; 
         ndeadlay=0;
       }
       Ndead++;
      }
     while (idead!=0);
   }

   char hname[50];
   TH1F* extr_occupancy[5][14][4][3][4];

   for (int iwh=0; iwh<5; iwh++){
     for (int ise=0; ise<14; ise++) {
       for (int ist=0; ist<4; ist++) {
         if (ist!=3 && ise>11) continue;
         for (int isl=0; isl<3; isl++) {
           if (isl==1 && ist==3) continue;
           for (int ilay=0; ilay<4; ilay++) {
             sprintf (hname,"Extr_Wheel%u_Sect%u_MB%u_SL%u_Lay%u",iwh,ise+1,ist+1,isl+1,ilay+1);
             extr_occupancy[iwh][ise][ist][isl][ilay] = new TH1F (hname,"",92,1.,93.);
	   }
	 }        
       }
     }
   }

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%1000 == 0) cout<<" Post-Loop evento "<<jentry<<endl;

      for (int iseg=0; iseg<Ndtsegments; iseg++) {

        if (!dtsegm4D_hasPhi->at(iseg)) continue; 
        if (dtsegm4D_station->at(iseg)!=4 && !dtsegm4D_hasZed->at(iseg)) continue; 

        //selection
        
        if (fabs(dtsegm4D_x_dir_loc->at(iseg))>0.7)   continue; // angle
        if (dtsegm4D_phinhits->at(iseg)<nrequiredhit) continue;

        // fill occupancy histos for segment positions

	//Extrapolated vector of hitted wires. It'used to check if a 
        TVectorF *expWire=(TVectorF*)dtsegm4D_hitsExpWire->At(iseg);

        // If a hit is missing, let us check that the extrapolation doesn't cross a dead cell!
        int NexpDead=0;

        if (dtsegm4D_phinhits->at(iseg) < 8 ) {  // <8 means at least one hit is missing
  
         for (int iex=0; iex<12; iex++) {
	   int expSL = 1;
           int expLay = iex+1;
           if (dtsegm4D_station->at(iseg) != 4){
	     if (iex > 3 && iex<8) {expSL=2; expLay-=4;}
             else if (iex>7) {expSL=3; expLay-=8;}
	   }
           else {
             if (iex > 3 && iex<8) continue;
             else if (iex > 7) {expSL=3; expLay-=8;}
	   }

           Float_t expW = (*expWire)(iex);
        
           for (int idead=0; idead<Ndead; idead++) {
    	    if (dead[idead][0] != dtsegm4D_wheel->at(iseg)) continue;
            if (dead[idead][1] != dtsegm4D_sector->at(iseg)) continue;
            if (dead[idead][2] != dtsegm4D_station->at(iseg)) continue;
            if (dead[idead][3] != expSL) continue;  
            if (dead[idead][4] != expLay) continue;
            if (dead[idead][5] != expW) continue; 
            NexpDead++;
            break;
	   }
           if (NexpDead==MaxDead) break;
	 }
         if (NexpDead==MaxDead) continue;
	}

	// this segment doesn't cross any dead cell, let's fill its positions:
        for (int iex=0; iex<12; iex++) {
	  int expSL = 1;
          int expLay = iex+1;
          if (dtsegm4D_station->at(iseg) != 4){
	    if (iex > 3 && iex<8) {expSL=2; expLay-=4;}
            else if (iex>7) {expSL=3; expLay-=8;}
	  }
          else {
            if (iex > 3 && iex<8) continue;
            else if (iex > 7) {expSL=3; expLay-=8;}
	  }		

          Float_t expW = (*expWire)(iex);
    	  extr_occupancy[dtsegm4D_wheel->at(iseg)+2][dtsegm4D_sector->at(iseg)-1][dtsegm4D_station->at(iseg)-1]
	         [expSL-1][expLay-1]->Fill(float(expW));
	}
      }

      // Now Zed segments
      for (int iseg=0; iseg<Ndtsegments; iseg++) {

        //selection
        //if (!dtsegm4D_hasZed->at(iseg) || !dtsegm4D_hasPhi->at(iseg)) continue; 
        if (!dtsegm4D_hasZed->at(iseg) || dtsegm4D_phinhits->at(iseg)<nrequiredhit) continue; 
        //10-01-17 capire perchè theta efficiency va giù con l ainst.luminosity: è il denominatore che sale?!

        // fill occupancy histos for segment positions
        TVectorF *expWire=(TVectorF*)dtsegm4D_hitsExpWire->At(iseg);

        // If a hit a missing, let us check that the extrapolation doesn't cross a dead cell!
        int NexpDead=0;

        if (dtsegm4D_znhits->at(iseg) < 3 ) {

         for (int iex=4; iex<8; iex++) {
 	   int expSL = 2;
           int expLay = iex-3;		

           Float_t expW = (*expWire)(iex);
           for (int idead=0; idead<Ndead; idead++) {
  	    if (dead[idead][0] != dtsegm4D_wheel->at(iseg)) continue;
            if (dead[idead][1] != dtsegm4D_sector->at(iseg)) continue;
            if (dead[idead][2] != dtsegm4D_station->at(iseg)) continue;
            if (dead[idead][3] != expSL) continue;  
            if (dead[idead][4] != expLay) continue;
            if (dead[idead][5] != expW) continue; 
            NexpDead++;
            break;   
	   }
           if (NexpDead == MaxDead) break; 
	 }
         if (NexpDead == MaxDead) continue; // this segment crosses at least 1 dead cell: drop it!
	}

        for (int iex=4; iex<8; iex++) {
	  int expSL = 2;
          int expLay = iex-3;		

          Float_t expW = (*expWire)(iex);

    	  extr_occupancy[dtsegm4D_wheel->at(iseg)+2][dtsegm4D_sector->at(iseg)-1][dtsegm4D_station->at(iseg)-1]
	         [expSL-1][expLay-1]->Fill(float(expW));
	}       
      }
   }


   // analyze extr_occupancy histos and fill table for undef channels

   int nwire=0; int NwireTot=0;
   for (int iw=0; iw<5000; iw++) for (int geo=0; geo<6; geo++) undef[iw][geo]=0;

   for (int iwh=0; iwh<5; iwh++){
     for (int ise=0; ise<14; ise++) {
       for (int ist=0; ist<4; ist++) {
         if (ist!=3 && ise>11) continue;
         for (int isl=0; isl<3; isl++) {
           if (isl==1 && ist==3) continue;
           if (isl==1){
	     //chimney chambers
	     if (( iwh==1 && ise==2 ) || ( iwh==3 && ise==3 ) ) nwire = 48;
	     else nwire=57; 
	   }
           else if (ist==0) nwire = 49;
           else if (ist==1) nwire = 60;
           else if (ist==2) nwire = 72;

           else if (ist == 3) {
	     if (ise == 3 || ise == 12 )    nwire = 72;
             else if ( ise==7 || ise ==11 ) nwire = 92;
             else if ( ise==8 || ise ==10 ) nwire = 48;
             else if ( ise==9 || ise ==13 ) nwire = 60;
             else nwire = 96;
	   }
	   
           NwireTot+=(nwire*4);

           for (int ilay=0; ilay<4; ilay++) {
             for (int iw=1; iw<nwire+1; iw++) {
               if (extr_occupancy[iwh][ise][ist][isl][ilay]->GetBinContent(iw)==0){

	       bool expDead=false;
               for (int idead=0; idead<Ndead; idead++) {
  	          if (dead[idead][0] != iwh-2)  continue;
                  if (dead[idead][1] != ise+1)  continue;
                  if (dead[idead][2] != ist+1)  continue;
                  if (dead[idead][3] != isl+1)  continue;  
                  if (dead[idead][4] != ilay+1) continue;
                  if (dead[idead][5] != iw)     continue; 
                  expDead=true;
                  break;   
	       }
	         if (expDead) continue;

		 undef[Nundef][0]=iwh-2;
		 undef[Nundef][1]=ise+1;
 		 undef[Nundef][2]=ist+1;
		 undef[Nundef][3]=isl+1;
		 undef[Nundef][4]=ilay+1;
		 undef[Nundef][5]=iw;  
            
                 UndefList <<Nundef+1<<" ";
                 for (int ip=0; ip<6; ip++) UndefList<< undef[Nundef][ip]<<" "; UndefList<<endl;
                 Nundef++; 
	       }
	     }
             for (int iw=nwire+1; iw<93; iw++) {
               if (extr_occupancy[iwh][ise][ist][isl][ilay]->GetBinContent(iw)!=0){
		 cout<<" warning extrapolation out of layer! YB"<<iwh-2<<"/Sec"<<ise+1<<"/MB"<<ist+1
                     <<" SL "<<isl+1<<" layer "<<ilay+1<<" cell "<<iw<<endl;
	       }
	     }
	   }
	 }
       }
     }
   }

   cout<<Nundef<<" undef wires out of "<<NwireTot<<endl;

   UndefList<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "
          <<endl<<Nundef<<" undef wires of out "<<NwireTot<<endl;


   for (int iundef=0; iundef<Nundef; iundef++) {
     cout<<"YB"<<undef[iundef][0]<<"/Sec"<<undef[iundef][1]<<"/MB"<<undef[iundef][2]<<"_SL"<<undef[iundef][3]
         <<"_layer"<<undef[iundef][4]<<"_cell"<<undef[iundef][5]<<endl;
   }

   UndefList.close();
}

