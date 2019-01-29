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
#include "plotter.h"
#include "TStopwatch.h"

int dead[10000][6];
int Ndead=0;
int nrequiredhit=7;
int MaxDead=0;

/* 
note: MaxDead = maximum number of dead channels crossed by a selected segment
MaxDead = 0 allows an umbiased selection.
If MaxDead > 0, dead wires will have efficiency zero while alive wires in the same segment won't contribute!
*/

int lessentries=50000;

void EfficiencyMonitor::PreLoop(){

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nRealEntries = fChain->GetEntries();

   // nentries=lessentries;

   Long64_t nbytes = 0, nb = 0;

   ofstream DeadList;
   std::string deadname;

   deadname.append("data/DeadList/DeadList_"+fileName);

   //   cout<<deadname<<endl;
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
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%50000 == 0) cout<<setprecision (2)<<" Pre-Loop event "<<jentry<<" "<<jentry/float(nRealEntries)*100.<<"%"<<endl;
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
	       //occupancy[iwh][ise][ist][isl][ilay]->Delete(); //check. delete or find solution. 
	     }
	   }
	 }
       }
     }
   }


   cout<<Ndead<<" dead wires out of "<<NwireTot<<endl;
   DeadList<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "
          <<endl<<Ndead<<" dead wires of out "<<NwireTot<<endl;

   //Uncomment to see list of dead channels.

   // for (int idead=0; idead<Ndead; idead++) {
   //   cout<<"YB"<<dead[idead][0]<<"/Sec"<<dead[idead][1]<<"/MB"<<dead[idead][2]<<"_SL"<<dead[idead][3]
   //       <<"_layer"<<dead[idead][4]<<"_cell"<<dead[idead][5]<<endl;
   // }
   
   DeadList.close();
}


void EfficiencyMonitor::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries     = fChain->GetEntriesFast();
   Long64_t nRealEntries = fChain->GetEntries();

   char go;

   // nentries=lessentries;

   Long64_t nbytes = 0, nb = 0;


// DEAD CHANNELS (to skip):

cout<<" within Loop: Ndead "<<Ndead<<endl;

if (Ndead==0) {
     cout<<" READING FROM FILE !! "<<endl;
     ifstream txtin;
     std::string deadname;
          
     deadname.append("data/DeadList/DeadList_"+fileName);
     
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

   cout<<"Reading tree"<<endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     Long64_t ientry = LoadTree(jentry);
     
     if (jentry%50000 == 0) cout<<"Event "<<jentry<<" "<<setprecision (2)<<((jentry/float(nRealEntries))*100)<<" %  "<<endl;
     if (ientry < 0) break;

     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

     if( !(std::find(dataContext.var["Run"].slice.begin(), dataContext.var["Run"].slice.end(), runnumber) != dataContext.var["Run"].slice.end()) ) continue; //fixme

     bool doBkg = kFALSE;
     for(auto const& ivar : dataContext.var) { 
       if (ivar.second.name == "InsLumi")     dataContext.var[ivar.first].value = lumiperblock;  //add secons value for slices here? 
       else if (ivar.second.name =="Pileup" ) dataContext.var[ivar.first].value = PV_Nvtx;
       else if (ivar.second.name =="Run")     dataContext.var[ivar.first].value = (float)runnumber; 
       else if (ivar.second.name =="Empty")   dataContext.var[ivar.first].value = 0.; // Default value for variables with no projection option 
       //else if (ivar.second.name =="BunchX" ) dataContext.var[ivar.first].value = bunchXing; not used now

       //Check if there are variables to plotted as a function of the background. If yes compute the background.
       if(ivar.second.doBkg) doBkg= kTRUE;
     }
     if(doBkg)   getBkgDigi(jentry);     

     // if (lumiperblock > dataContext.slices[0].back()) { cout<<" luminosity out of range!! "<< lumiperblock<<endl; continue; }
     // if (PV_Nvtx > dataContext.slices[1].back())      { cout<<" PU out of range!! "        << PV_Nvtx<<endl;      continue; } 
   
     // First search for Phi segments
          
     for (int iseg=0; iseg<Ndtsegments; iseg++) {
       //selection
       if (!dtsegm4D_hasPhi->at(iseg)) continue;
       // In chambers 1,2,3 select only segments with also Z (theta) contribution.
       if (dtsegm4D_station->at(iseg)!=4 && !dtsegm4D_hasZed->at(iseg)) continue;
       
       int seg_phinhits = dtsegm4D_phinhits->at(iseg);
       
       if (fabs(dtsegm4D_x_dir_loc->at(iseg))>0.7) continue; // angle
       
       TVectorF *expWire=(TVectorF*)dtsegm4D_hitsExpWire->At(iseg);
       
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

       // background value add only here because depend on chamber, section ecc...
       if(dataContext.name=="Fixed") dataContext.var["Bkg"].value = bkgCounts.DigiBkgWhMB[dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)]; //check. Only for fixed? maybe an external option would be better


       if (NHits<nrequiredhit) continue;
       else if (NHits==8) {
	 for (int sl=0; sl<2; sl++) for (int lay=0; lay<4; lay++) {
	     //variables, points, stations, wheels 

	     for(auto const& ivar : dataContext.var) {
	       if(!ivar.second.doEff) continue;
	       plots->Eff_phiMBWh[ivar.first][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(1,ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value );    
	       plots->EffA_phiMBWh[ivar.first][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(1,ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
	       // extra chamber of sector 4 (sector 13)

	       if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
		 plots->Eff_phiMB4Top[ivar.first][dtsegm4D_wheel->at(iseg)+2]->Fill(1,ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
		 plots->EffA_phiMB4Top[ivar.first][dtsegm4D_wheel->at(iseg)+2]->Fill(1,ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
	       }
  
	       // extra chamber of sector 10 (sector 14) 
	       else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
		 plots->Eff_phiMB4Bot[ivar.first]->Fill(1,ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
		 plots->EffA_phiMB4Bot[ivar.first]->Fill(1,ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
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

	     // let's loop all over the digis and take the closest digis to the extrapolated wire. 
	     // Think about an extra condition on time.

             if (fabs(expW-digi_wire->at(idigi))<fabs(d)) {
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
	     
	     for(auto const& ivar : dataContext.var) { 
	       if(!ivar.second.doEff) continue;
	       plots->Eff_phiMBWh[ivar.first][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(missingLayer[imiss][1],ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
	       if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
		 plots->Eff_phiMB4Top[ivar.first][dtsegm4D_wheel->at(iseg)+2]->Fill(missingLayer[imiss][1],ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 

	       }
	       else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {		
		 plots->Eff_phiMB4Bot[ivar.first]->Fill(missingLayer[imiss][1],ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value );

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
	       for(auto const& ivar : dataContext.var) { 
		 if(!ivar.second.doEff) continue;
		 plots->Eff_phiMBWh[ivar.first][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(!(missAss&&missDigi),ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
		 plots->EffA_phiMBWh[ivar.first][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(!(missAss),ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
		 
		 if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
		   plots->Eff_phiMB4Top[ivar.first][dtsegm4D_wheel->at(iseg)+2]->Fill(!(missAss&&missDigi),ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
		   plots->EffA_phiMB4Top[ivar.first][dtsegm4D_wheel->at(iseg)+2]->Fill(!(missAss),ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
		 }
		 else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {		
		   plots->Eff_phiMB4Bot[ivar.first]->Fill(!(missAss&&missDigi),ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value );
		   plots->EffA_phiMB4Bot[ivar.first]->Fill(!(missAss),ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value );
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
	   for(auto const& ivar : dataContext.var) {
	     if(!ivar.second.doEff) continue;
	     plots->Eff_theMBWh[ivar.first][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(1,ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
	     plots->EffA_theMBWh[ivar.first][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill(1,ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
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
	 for(auto const& ivar : dataContext.var){ 
	   if(!ivar.second.doEff) continue;
	   plots->Eff_theMBWh[ivar.first][dtsegm4D_wheel->at(iseg)+2][dtsegm4D_station->at(iseg)-1]->Fill((fabs(d)< 1.1) ,ivar.second.value,dataContext.var[(ivar.second.projVar).c_str()].value ); 
	 }
       }
       else {
	 cout<<" what do you want?? NHits (z) = "<<NHits<<endl;
	 return;
       }
     }
    
   }
}

void EfficiencyMonitor::write(){
  plots->write();     
}

void EfficiencyMonitor::plot(){
  plots->plot(outName.c_str());     
}

void EfficiencyMonitor::close(){
  plots->close();
}

void EfficiencyMonitor::SetRunSlices()
{

  if (fChain == 0) return;   
  Long64_t nentries     = fChain->GetEntriesFast();
  Long64_t nRealEntries = fChain->GetEntries();

  /* Loop used to collect the list of runs containend inside the TTree */

  cout<<"Loop over the tree to create run list"<<endl;
  cout<<"Tree composed of "<<nRealEntries<<" entries "<<endl; 
  std::set<Int_t> runNumber_Set; // ordered list of runnumbers
  
  for (Long64_t jentry=2; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
   
    if (ientry < 0)  break;  
    if (jentry%500000 == 0) cout<<"\rLoop event "<<jentry<<" "<<setprecision (2)<<(100*jentry/float(nRealEntries))<<" %    "<<endl;    

    b_runnumber->GetEntry(jentry); 
    runNumber_Set.insert(runnumber);
  }

  for (std::set<int>::iterator it= runNumber_Set.begin(); it!= runNumber_Set.end(); ++it)
    dataContext.var["Run"].slice.push_back((float)*it);

  cout<<"original runs "<<dataContext.var["Run"].slice.size()<<endl;
  
  dataContext.var["Run"].slice.push_back( dataContext.var["Run"].slice.back()+1); //Add another element in order to have another bin for the last lumi point. It will be removed in case of concatantion with another file.
  
  // uncoment to check the run list
  

  cout<<"Run slice "<<endl;
  for(uint i =0; i<dataContext.var["Run"].slice.size(); i++)
  cout<<setprecision(6)<<dataContext.var["Run"].slice.at(i)<<endl;

}

void EfficiencyMonitor::checkPuLumiRatio(){

  cout<<"Loop over the tree to create run list and check lumi PU ratio"<<endl;
  if (fChain == 0) return;   
  Long64_t nentries     = fChain->GetEntriesFast();
  Long64_t nRealEntries = fChain->GetEntries();

  std::set<Int_t> runNumber_Set;

  TH2F * hPULumiVSrun = new TH2F("hPULumiVSrun","hPULumiVSrun", dataContext.var["Run"].slice.size()-1, 0,dataContext.var["Run"].slice.size()-1, 200,0,0.025);

  TH2F * hPUvsLumi = new TH2F("hPUvsLumi","PUvsLumi", dataContext.var["InsLumi"].slice.size()-1, dataContext.var["InsLumi"].slice.data(),dataContext.var["Pileup"].slice.size()-1, dataContext.var["Pileup"].slice.data());

  TH2F * hPUvsLumiCut = new TH2F("hPUvsLumiCut","PUvsLumiCut", dataContext.var["InsLumi"].slice.size()-1, dataContext.var["InsLumi"].slice.data(),dataContext.var["Pileup"].slice.size()-1, dataContext.var["Pileup"].slice.data());


  hPULumiVSrun->SetLabelSize(0.03);
  hPULumiVSrun->SetMarkerSize(0.8);

  map<int,int> runnumberMap;

  int i = 0;
  for (std::vector<double>::iterator it = dataContext.var["Run"].slice.begin() ; it != dataContext.var["Run"].slice.end()-1; ++it){
    runnumberMap[*it]=i;
    i++;
    hPULumiVSrun->GetXaxis()->SetBinLabel(i,to_string((int)*it).c_str());
  }

  TH1F * hRatio = new TH1F("hratio","hratio",200,0,0.025);
  
  dataContext.var["Run"].slice.clear();

  Int_t preRunNumber = 0;


  Float_t puLumiCut = 0.0026;


  for (Long64_t jentry = 0; jentry<nentries; jentry++) {

    if( jentry==nRealEntries-1 ) cout<<jentry<< " "<<nRealEntries-1<<endl;
 
    Long64_t ientry = LoadTree(jentry);

    b_runnumber->GetEntry(jentry);
    b_lumiperblock->GetEntry(jentry);
    b_PV_Nvtx->GetEntry(jentry);

    if (ientry < 0)  break;
    if (jentry%50000 == 0) cout<<"\rLoop event "<<jentry<<" "<<setprecision (2)<<(100*jentry/float(nRealEntries))<<" %    "<<endl;

    if(jentry==0) preRunNumber = runnumber;

    hPULumiVSrun->Fill(runnumberMap[runnumber], PV_Nvtx/lumiperblock);
    hPUvsLumi->Fill(lumiperblock,PV_Nvtx);
  }


  TProfile *pr =  hPULumiVSrun->ProfileX("_pfx",1,-1,"dsame");

  for(int bin = 1; bin<=hPULumiVSrun->GetXaxis()->GetNbins(); bin++){
    //    cout<<bin<<" "<<hPULumiVSrun->GetXaxis()->GetBinLabel(bin)<<" "<<pr->GetBinContent(bin)<<endl;
    if(pr->GetBinContent(bin) < puLumiCut){
      runNumber_Set.insert(atoi(hPULumiVSrun->GetXaxis()->GetBinLabel(bin)));
    }
  }


  for (Long64_t jentry = 0; jentry<nentries; jentry++) {

    if( jentry==nRealEntries-1 ) cout<<jentry<< " "<<nRealEntries-1<<endl;
 
    Long64_t ientry = LoadTree(jentry);

    b_runnumber->GetEntry(jentry);
    b_lumiperblock->GetEntry(jentry);
    b_PV_Nvtx->GetEntry(jentry);

    if (ientry < 0)  break;
    if (jentry%50000 == 0) cout<<"\rLoop event "<<jentry<<" "<<setprecision (2)<<(100*jentry/float(nRealEntries))<<" %    "<<endl;

    if(runNumber_Set.find(runnumber)!=runNumber_Set.end())    hPUvsLumiCut->Fill(lumiperblock,PV_Nvtx);
  }

  for (std::set<int>::iterator it= runNumber_Set.begin(); it!= runNumber_Set.end(); ++it)
    dataContext.var["Run"].slice.push_back((float)*it);
  
  dataContext.var["Run"].slice.push_back( dataContext.var["Run"].slice.back()+1); //Add another element in order to have another bin for the last lumi point. It will be removed in case of concatantion with another file.
  
  cout<<"Run slices"<<endl;
  for(uint i =0; i<dataContext.var["Run"].slice.size(); i++)
    cout<<setprecision(6)<<dataContext.var["Run"].slice.at(i)<<endl;

  TCanvas *cGlob = new TCanvas("cGlob","",1200,400);

  hPULumiVSrun->SetMarkerStyle(20);
  hPULumiVSrun->SetMarkerColor(2);
  hPULumiVSrun->Draw("colz");
  hPULumiVSrun->ProfileX("_pfx",1,-1,"dsame"); 

  hPULumiVSrun->SetTitle("Pileup/Instant lumi vs run ;run;Pileup/Instant lumi"); 

  if(stat((dataContext.webFolder+"/"+outName+"/").c_str(),&st) != 0){
    system(("mkdir "+dataContext.webFolder+"/"+outName).c_str());
    system(("cp "+dataContext.webFolder+"/index.php " +dataContext.webFolder+"/"+outName).c_str());
  }

  if(stat((dataContext.webFolder+"/"+outName+"/Global").c_str(),&st) != 0){
    system(("mkdir "+dataContext.webFolder+"/"+outName+"/Global").c_str());
    system(("cp "+dataContext.webFolder+"/"+outName+"/index.php " +dataContext.webFolder+"/"+outName+"/Global").c_str());
  }

  cGlob->SaveAs((dataContext.webFolder+"/"+outName+"/Global/PileUpVsLumiVsRun.png").c_str());
  cGlob->SaveAs((dataContext.webFolder+"/"+outName+"/Global/PileUpVsLumiVsRun.root").c_str());

  hPUvsLumi->Draw("Colz");
  hPUvsLumi->SetTitle((dataContext.var["Pileup"].Title+" vs "+dataContext.var["InsLumi"].Title+";"+dataContext.var["Pileup"].Label+";"+dataContext.var["InsLumi"].Label).c_str());

  cGlob->SaveAs((dataContext.webFolder+"/"+outName+"/Global/PileUpVsLumi.png").c_str());
  cGlob->SaveAs((dataContext.webFolder+"/"+outName+"/Global/PileUpVsLumi.root").c_str());


  hPUvsLumiCut->Draw("Colz");
  hPUvsLumiCut->SetTitle((dataContext.var["Pileup"].Title+" vs "+dataContext.var["InsLumi"].Title+";"+dataContext.var["Pileup"].Label+";"+dataContext.var["InsLumi"].Label).c_str());

  cGlob->SaveAs((dataContext.webFolder+"/"+outName+"/Global/PileUpVsLumiCut.png").c_str());
  cGlob->SaveAs((dataContext.webFolder+"/"+outName+"/Global/PileUpVsLumiCut.root").c_str());

  if(dataContext.var["Run"].slice.size()==1) {
    cout<<"No runs with pu/lumi < "<<puLumiCut<<endl;
    abort();

  }
}


void EfficiencyMonitor::checkRunStat(){

  cout<<"Loop over the tree to create run list and select runs with enough statistic"<<endl;
  if (fChain == 0) return;   
  Long64_t nentries     = fChain->GetEntriesFast();
  Long64_t nRealEntries = fChain->GetEntries();

  Int_t preRunNumber = 0;
 
  map<int,int> runnumberMap;


  dataContext.var["Run"].slice.clear();

  for (Long64_t jentry = 0; jentry<nentries; jentry++) {
 
    Long64_t ientry = LoadTree(jentry);

    b_runnumber->GetEntry(jentry);

    if (ientry < 0)  break;
    //    if (jentry%50000 == 0) cout<<"\rLoop event "<<jentry<<" "<<setprecision (2)<<(100*jentry/float(nRealEntries))<<" %    "<<endl;

    runnumberMap[runnumber]+=1;
  }

  int minRunStat = 75000;

  for( const auto& run_pair : runnumberMap ){

    //to remove single runs with problems
    //    if(run_pair.first==302635 || run_pair.first==304737 || run_pair.first== 304740 || run_pair.first==304739 || run_pair.first==304738 ) continue;
    if( run_pair.second > minRunStat )  dataContext.var["Run"].slice.push_back((float)run_pair.first); 

    //    if( run_pair.first > 304737 && run_pair.first < 304740 && run_pair.second > 15000)  dataContext.var["Run"].slice.push_back((float)run_pair.first); 
  }
  
  if(dataContext.var["Run"].slice.size()==0) {
    cout<<"No runs with more than "<<minRunStat<<" over the original "<<runnumberMap.size()<<" runs. Change run/file or change threshold "<<endl;
    abort();
  }
  else
    cout<<"selected "<<dataContext.var["Run"].slice.size()<<" runs over the original "<<runnumberMap.size()<<" runs, corresponding to "<<(dataContext.var["Run"].slice.size()/static_cast<double>(runnumberMap.size()))*100<<" \%"<<endl;
  
  dataContext.var["Run"].slice.push_back( dataContext.var["Run"].slice.back()+1); //Add another element in order to have another bin for the last lumi point. It will be removed in case of concatantion with another file.

  cout<<"Run slices"<<endl;
  for(uint i =0; i<dataContext.var["Run"].slice.size(); i++)
    cout<<setprecision(6)<<dataContext.var["Run"].slice.at(i)<<endl;
}


void EfficiencyMonitor::getBkgDigi(Int_t jentry){
 
  double timewindowdigi =  1275;
  double timewindowseg  =  800;
 
 float MBarea[3]; float MB4area[14];

 MBarea[0] = 49 * 4.2 * 58 * 4.2;  //n wires phi X cell size X n wires theta.
 MBarea[1] = 60 * 4.2 * 58 * 4.2;
 MBarea[2] = 72 * 4.2 * 58 * 4.2;
 
 MB4area[1-1]  = 96 * 4.2 * 58 * 4.2;
 MB4area[2-1]  = MB4area[1-1];
 MB4area[3-1]  = MB4area[1-1];
 
 MB4area[4-1]  = 72 * 4.2 * 58 * 4.2;
 MB4area[5-1]  = MB4area[1-1];
 MB4area[6-1]  = MB4area[1-1];
 MB4area[7-1]  = MB4area[1-1];

 MB4area[8-1]  = 92 * 4.2 * 58 * 4.2;
 MB4area[9-1]  = 49 * 4.2 * 58 * 4.2;
 MB4area[10-1] = 60 * 4.2 * 58 * 4.2;
 MB4area[11-1] = MB4area[9-1];
 MB4area[12-1] = MB4area[8-1];
 MB4area[13-1] = MB4area[4-1];
 MB4area[14-1] = MB4area[10-1];


  Long64_t nentries     = fChain->GetEntriesFast();  

   int NdigiSig[5][12][4];        int NdigiBkg[5][12][4]; //wheel sector stations
   int NsegmSig[5][12][4];        int NsegmBkg[5][12][4]; 
   
   //intilize to 0 the elements
   for (int iwh=0; iwh<5; iwh++) for (int ise=0; ise<12; ise++) for (int ist=0; ist<4;ist++) {
	 
	 NdigiSig[iwh][ise][ist]=0;    NdigiBkg[iwh][ise][ist]=0;
	 NsegmSig[iwh][ise][ist]=0;    NsegmBkg[iwh][ise][ist]=0;
       }
   

   int ChambCross[100][3];
   for (int i=0; i<100; i++) for (int geo=0; geo<3; geo++) ChambCross[i][geo]=-999;
   
   int NChambCross=0;

   //   cout<<"Nmuons "<<Nmuons<<endl;;   
   for (int imu=0; imu<Nmuons; imu++) {
     
     int Ncross = Mu_nMatches->at(imu); 

     //     if (Mu_isMuGlobal->at(imu) || Mu_isMuTracker->at(imu)) check->Fill(Mu_eta->at(imu),float(Ncross));                  
     TVectorF *Mu_cross_Wh=(TVectorF*)Mu_matches_Wh->At(imu);
     TVectorF *Mu_cross_Se=(TVectorF*)Mu_matches_Sec->At(imu);
     TVectorF *Mu_cross_St=(TVectorF*)Mu_matches_St->At(imu);
     
     for (int ich=0; ich<Ncross; ich++) {
       
       int wheel    = (*Mu_cross_Wh)(ich);
       int sector   = (*Mu_cross_Se)(ich);
       int station  = (*Mu_cross_St)(ich);
       
       ChambCross[NChambCross][0] = wheel;
       ChambCross[NChambCross][1] = sector;
       ChambCross[NChambCross][2] = station;
       
       if (NChambCross<99) NChambCross++;
       else {
	 cout<<" warning event "<<jentry<<" "<<Nmuons<<" muons, more than 100 chamber crossed! "<<endl;
	 break;
       }
     }

     // add search for associated segments (standalone muons may not have matches)
  
     // if (Ncross>0) continue;
     
     int word = Mu_segmentIndex_sta->at(imu);
     int count=0;
       
     while (word>1) {
       while (word % 2 == 0) {
	 count++;
	 word=word>>1;
       }

       ChambCross[NChambCross][0] = dtsegm4D_wheel->at(count);
       ChambCross[NChambCross][1] = dtsegm4D_sector->at(count);
       ChambCross[NChambCross][2] = dtsegm4D_station->at(count);
       
       if (Mu_isMuGlobal->at(imu) || Mu_isMuTracker->at(imu)) {
	 float qual = dtsegm4D_phinhits->at(count) > 0 ?  dtsegm4D_phinhits->at(count) : 0.;
	 qual += dtsegm4D_znhits->at(count) > 0 ?  dtsegm4D_znhits->at(count) : 0.;
	 // qualsig->Fill(qual); check
	 // if (qual>5) t0sig->Fill(dtsegm4D_t0->at(count));
       }
       if (NChambCross<99) NChambCross++;
       else {
	 cout<<" warning event "<<jentry<<" "<<Nmuons<<" muons, more than 100 chamber crossed! "<<endl;
	 break;
       }
       count++; word=word>>1;
     }
   } // leave muon loop: now I have all chambers crossed + associated segments!
   
   for (int idigi=0; idigi<Ndigis; idigi++) {
     
     bool digibkg = true; 
       
     for (int icross=0; icross<NChambCross; icross++) {
       
       if (digi_wheel->at(idigi)   == ChambCross[icross][0] &&
	   digi_sector->at(idigi)  == ChambCross[icross][1] &&
	   digi_station->at(idigi) == ChambCross[icross][2] ) {
	 digibkg=false;
	 
	 //	 timesig->Fill(digi_time->at(idigi));
	 if (digi_sector->at(idigi)<13)        NdigiSig[digi_wheel->at(idigi)+2][digi_sector->at(idigi)-1][digi_station->at(idigi)-1]++;
	 else  if (digi_sector->at(idigi)==13) NdigiSig[digi_wheel->at(idigi)+2][4-1][digi_station->at(idigi)-1]++;
	 else  if (digi_sector->at(idigi)==14) NdigiSig[digi_wheel->at(idigi)+2][10-1][digi_station->at(idigi)-1]++;
	 
	 break;
       }
     }
     
     if (digibkg) {
       //       timebkg->Fill(digi_time->at(idigi));
       if (digi_sector->at(idigi)<13)        NdigiBkg[digi_wheel->at(idigi)+2][digi_sector->at(idigi)-1][digi_station->at(idigi)-1]++;
       else if (digi_sector->at(idigi)==13)  NdigiBkg[digi_wheel->at(idigi)+2][4-1][digi_station->at(idigi)-1]++;
       else if (digi_sector->at(idigi)==14)  NdigiBkg[digi_wheel->at(idigi)+2][10-1][digi_station->at(idigi)-1]++;
     }
   } // loop sui digi
   
   for (int iseg=0; iseg<Ndtsegments; iseg++) {
     
     bool segmbkg = true; 
     
     float qual = dtsegm4D_phinhits->at(iseg) > 0 ?  dtsegm4D_phinhits->at(iseg) : 0.;
     qual += dtsegm4D_znhits->at(iseg) > 0 ?  dtsegm4D_znhits->at(iseg) : 0.;
     
     for (int icross=0; icross<NChambCross; icross++) {
       
       if (dtsegm4D_wheel->at(iseg)   == ChambCross[icross][0] &&
	   dtsegm4D_sector->at(iseg)  == ChambCross[icross][1] &&
	   dtsegm4D_station->at(iseg) == ChambCross[icross][2] ) {
	 segmbkg=false;
	 
	 if(dtsegm4D_sector->at(iseg)<13)       NsegmSig[dtsegm4D_wheel->at(iseg)+2][dtsegm4D_sector->at(iseg)-1][dtsegm4D_station->at(iseg)-1]++;
	 else if(dtsegm4D_sector->at(iseg)==13) NsegmSig[dtsegm4D_wheel->at(iseg)+2][4-1][dtsegm4D_station->at(iseg)-1]++;
	 else if(dtsegm4D_sector->at(iseg)==14) NsegmSig[dtsegm4D_wheel->at(iseg)+2][10-1][dtsegm4D_station->at(iseg)-1]++;
	 
	 break;
       }
     }
     
     if (segmbkg) {
       //       qualbkg->Fill(qual);
       //       if (qual>5) t0bkg->Fill(dtsegm4D_t0->at(iseg));
       
       if(dtsegm4D_sector->at(iseg)<13)       {
	 NsegmBkg[dtsegm4D_wheel->at(iseg)+2][dtsegm4D_sector->at(iseg)-1][dtsegm4D_station->at(iseg)-1]++;
       }
       else if(dtsegm4D_sector->at(iseg)==13) {
	 NsegmBkg[dtsegm4D_wheel->at(iseg)+2][4-1][dtsegm4D_station->at(iseg)-1]++;
       }
       else if(dtsegm4D_sector->at(iseg)==14) {
	 NsegmBkg[dtsegm4D_wheel->at(iseg)+2][10-1][dtsegm4D_station->at(iseg)-1]++;
       }
     }
   }

   for (int iwh=0; iwh<5; iwh++) {
     for (int ise=0; ise<12; ise++) {
       for (int ist=0; ist<4; ist++) {
	 double normdigi;  double normseg;
	 
	 if (ist+1 < 4) {
	   normdigi = MBarea[ist] * 12; // 12 layers
	   normseg  = MBarea[ist]; 
	 }
	 
	 else {
	   normdigi = MB4area[ise] * 8; // 8 layers
	   normseg  = MB4area[ise]; 
	 }
	 
	 if (ist+1 ==4 && (ise+1 == 4 ||  ise+1 == 10) ) {normdigi *=2; normseg *=2;}
	 normdigi *= timewindowdigi/1e9;// * nentries; //1e9 = conversion from ns to s
	 normseg  *= timewindowseg/1e9;  //* nentries;
	 
	 // if(iwh==0 && ist==0) cout<<"norm digi "<<normdigi <<endl;

	 bkgCounts.DigiSigWhMB[iwh][ist]  = float(NdigiSig[iwh][ise][ist]) /normdigi;
	 bkgCounts.DigiBkgWhMB[iwh][ist]  = float(NdigiBkg[iwh][ise][ist]) /normdigi;
	 bkgCounts.DigiSigSecMB[ise][ist] = float(NdigiSig[iwh][ise][ist]) /normdigi;
	 bkgCounts.DigiBkgSecMB[ise][ist] = float(NdigiBkg[iwh][ise][ist]) /normdigi;
	 
	 bkgCounts.SegmSigWhMB[iwh][ist]  = float(NsegmSig[iwh][ise][ist]) /normseg;
	 bkgCounts.SegmBkgWhMB[iwh][ist]  = float(NsegmBkg[iwh][ise][ist]) /normseg;
	 bkgCounts.SegmSigSecMB[ise][ist] = float(NsegmSig[iwh][ise][ist]) /normseg;
	 bkgCounts.SegmBkgSecMB[ise][ist] = float(NsegmBkg[iwh][ise][ist]) /normseg;

	 for(auto const& ivar : dataContext.var) { 
	   if(!ivar.second.doBkg) continue;

	   plots->Dist_MBWh[ivar.first][iwh][ist]->Fill(ivar.second.value ,bkgCounts.DigiBkgWhMB[iwh][ist]); 
	   plots->Dist_SegMBWh[ivar.first][iwh][ist]->Fill(ivar.second.value ,bkgCounts.SegmBkgWhMB[iwh][ist]); 
	  
	   if (ist+1==4 && (ise+1==4 || ise+1 ==13)) {
	     plots->Dist_MB4Top[ivar.first][iwh]->Fill(ivar.second.value ,bkgCounts.DigiBkgWhMB[iwh][ist]);      
	     plots->Dist_SegMB4Top[ivar.first][iwh]->Fill(ivar.second.value ,bkgCounts.SegmBkgWhMB[iwh][ist]);   
	   }

	   else if (ist+1==4 && (ise+1==10 || ise+1==14)) {
	     plots->Dist_MB4Bot[ivar.first]->Fill(ivar.second.value,bkgCounts.DigiBkgWhMB[iwh][ist]);     
	     plots->Dist_SegMB4Bot[ivar.first]->Fill(ivar.second.value,bkgCounts.SegmBkgWhMB[iwh][ist]);     
	   }
	 } //End variables
       } // End stations
     } // End sector 
   } // End wheel
}
