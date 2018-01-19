#define EfficiencyMonitor_cxx
#include "EfficiencyMonitor.h"
#include <TVectorF.h>
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

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
   }
   cout<<deadname<<endl;
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


   // Lumi and PU bins:

   std::vector<float> lumislice = { 0., 1000., 2000., 3000., 4000., 
				    5000., 6000., 7000., 8000., 9000.,
				    10000.,11000.,12000.,13000.,14000.,
				    15000.,16000.,17000.,18000.,19000.,
				    20000.,21000.,22000.,23000.,24000.,
				    25000};

   std::vector<float> PUslice   = { 0., 2., 4., 6.,8.,
				    10.,12.,14.,16.,18.
				    ,20.,22.,24.,26.,28.
				    ,30.,32.,34.,36.,38.,
				    40.,44.,48.,60.,70.,
				    90.};
				    //				    50.,58.,62.,94.};

   std::vector<float> PUe;  //  = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
   std::vector<float> Lumie; //  = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1,};

   for (std::vector<float>::iterator lumi = lumislice.begin() ; lumi != lumislice.end(); ++lumi) Lumie.push_back(1.);
   for (std::vector<float>::iterator pu = PUslice.begin() ; pu != PUslice.end(); ++pu) PUe.push_back(1.);


   // COUNTERS:

   const int nPUpoints   = PUe.size();
   const int nLumiPoints = Lumie.size();

   std::cout<<"nPUpoints "<<nPUpoints<<" nLumiPoints "<<nLumiPoints<<std::endl;

   vector<vector<vector<vector<int>>>> v4
   for (int ivar=0; ivar<2; ivar++){
     v4.push_back( vector<vector<vector<int>>>())
     for (int ipoint=0; ipoint<nLumiPoints; ipoint++){
       v4[ivar].push_back(vector<vector<int> >())
       for (int iwh=0; iwh<5; iwh++){
	 v4[ivar][ipoint].push_back(vector<int>())
        for (int ist=0; ist<4; ist++){
	  v4[ivar][ipoint][ist].push_back(0);
	}
       }
     }
   }
		 
   int Num_phiMBWh[2][nLumiPoints][4][5];  int Den_phiMBWh[2][nLumiPoints][4][5]; //  Lumi, nLumiPoints punti, 4 stazioni, 5 ruote
   int Num_theMBWh[2][nLumiPoints][3][5];  int Den_theMBWh[2][nLumiPoints][3][5]; //  Lumi, nLumiPoints punti, 3 stazioni, 5 ruote
   int NumA_phiMBWh[2][nLumiPoints][4][5]; int NumA_theMBWh[2][nLumiPoints][3][5];// 'A' stands for 'Associated'

   int Num_phiMB4Top[2][5][nLumiPoints];   int Den_phiMB4Top[2][5][nLumiPoints];  //  PU, nLumiPoints punti
   int Num_phiMB4Bot[2][nLumiPoints];      int Den_phiMB4Bot[2][nLumiPoints];     //  Lumi nLumiPoints punti
   int NumA_phiMB4Top[2][5][nLumiPoints];  int NumA_phiMB4Bot[2][nLumiPoints];    // 'A' stands for 'Associated'



   // Set all to 0
   for (int ivar=0; ivar<2; ivar++){
     for (int ipoint=0; ipoint<nLumiPoints; ipoint++){
       for (int iwh=0; iwh<5; iwh++){
        for (int ist=0; ist<4; ist++){

 	  Num_phiMBWh[ivar][ipoint][ist][iwh]  = 0;
          NumA_phiMBWh[ivar][ipoint][ist][iwh] = 0;
	  Den_phiMBWh[ivar][ipoint][ist][iwh]  = 0;

          if (ist==3) continue;
	  Num_theMBWh[ivar][ipoint][ist][iwh]  = 0;
          NumA_theMBWh[ivar][ipoint][ist][iwh] = 0;
	  Den_theMBWh[ivar][ipoint][ist][iwh]  = 0;
        }
        Num_phiMB4Top[ivar][iwh][ipoint] = 0;
        NumA_phiMB4Top[ivar][iwh][ipoint]= 0;
        Den_phiMB4Top[ivar][iwh][ipoint] = 0;

       }
       Num_phiMB4Bot[ivar][ipoint] = 0;
       NumA_phiMB4Bot[ivar][ipoint]= 0;
       Den_phiMB4Bot[ivar][ipoint] = 0;
     }
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


   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     //  std::cout<<"entry "<<jentry<<std::endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%1000 == 0) cout<<" evento "<<jentry<<endl;

      // identify Lumibin and PUbin
      int Lumibin=-1; int PUbin=-1;


      //Find lumi e PU bin
      for (uint islice=0; islice<lumislice.size(); islice++) 
	if (lumiperblock>=lumislice[islice]&&lumiperblock<lumislice[islice+1]) Lumibin=islice;
      for (uint islice=0; islice<PUslice.size(); islice++) 
	if (PV_Nvtx>=PUslice[islice]&&PV_Nvtx<PUslice[islice+1]) PUbin=islice;
      
      if (Lumibin<0) { cout<<" luminosity out of range!! "<< lumiperblock<<endl; continue; }
      if (PUbin<0)   { cout<<" PU out of range!! "        << PV_Nvtx<<endl;      continue; }
         
      // First search for Phi segments

      for (int iseg=0; iseg<Ndtsegments; iseg++) {

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

        if (NHits<nrequiredhit) continue;
        else if (NHits==8) {

	  for (int sl=0; sl<2; sl++) for (int lay=0; lay<4; lay++) {


	     // 2 variabili (Lumi, PU), 22 punti, 4 stazioni, 5 ruote 
            Num_phiMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            NumA_phiMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            Den_phiMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

            Num_phiMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            NumA_phiMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            Den_phiMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

	    // extra chamber of sector 4 (sector 13)
            if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {

	      // 2 variabili (Lumi, PU), 22 punti
	      Num_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Lumibin]++;  
              NumA_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Lumibin]++; 
              Den_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Lumibin]++; 
              Num_phiMB4Top[1][dtsegm4D_wheel->at(iseg)+2][PUbin]++; 
              NumA_phiMB4Top[1][dtsegm4D_wheel->at(iseg)+2][PUbin]++; 
              Den_phiMB4Top[1][dtsegm4D_wheel->at(iseg)+2][PUbin]++; 
	    }

	    // extra chamber of sector 10 (sector 14) 
            else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
	      // 2 variabili (Lumi, PU), 22 punti
              Num_phiMB4Bot[0][Lumibin]++; 
              NumA_phiMB4Bot[0][Lumibin]++; 
              Den_phiMB4Bot[0][Lumibin]++;
              Num_phiMB4Bot[1][PUbin]++;
              NumA_phiMB4Bot[1][PUbin]++; 
              Den_phiMB4Bot[1][PUbin]++;
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

	       //denominator

               Den_phiMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
               Den_phiMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

               if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
                 Den_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Lumibin]++; 
                 Den_phiMB4Top[1][dtsegm4D_wheel->at(iseg)+2][PUbin]++; 
	       }
               else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
                 Den_phiMB4Bot[0][Lumibin]++;
                 Den_phiMB4Bot[1][PUbin]++;
	       }


               if (missingLayer[imiss][1]) {
		// numerator
                Num_phiMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
                Num_phiMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

                if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
                 Num_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Lumibin]++; 
                 Num_phiMB4Top[1][dtsegm4D_wheel->at(iseg)+2][PUbin]++; 
		}
                else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
                 Num_phiMB4Bot[0][Lumibin]++; 
                 Num_phiMB4Bot[1][PUbin]++; 
		}
	       }
	      }
	  }
	    
          else {

              for (int sl=0; sl<2; sl++) for (int lay=0; lay<4; lay++) {

	       //denominator
               Den_phiMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
               Den_phiMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

               if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
                 Den_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Lumibin]++; 
                 Den_phiMB4Top[1][dtsegm4D_wheel->at(iseg)+2][PUbin]++; 
	       }
               else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
                 Den_phiMB4Bot[0][Lumibin]++;
                 Den_phiMB4Bot[1][PUbin]++;
	       }

               bool missAss=false; bool missDigi=false;
               for (int imiss=0; imiss<nmissing; imiss++) {
                 if (missingLayer[imiss][0]==sl*4+lay+1) {
                      missAss=true;
		      if (!missingLayer[imiss][1]) missDigi=true;
		 }
	       }
               if (!(missAss&&missDigi)) {
		// numerator
                Num_phiMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
                Num_phiMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

                if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
                 Num_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Lumibin]++; 
                 Num_phiMB4Top[1][dtsegm4D_wheel->at(iseg)+2][PUbin]++; 
		}
                else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
                 Num_phiMB4Bot[0][Lumibin]++; 
                 Num_phiMB4Bot[1][PUbin]++; 
		}
	       }
               if (!(missAss)) {
		// numerator
                NumA_phiMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
                NumA_phiMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

                if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
                 NumA_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Lumibin]++; 
                 NumA_phiMB4Top[1][dtsegm4D_wheel->at(iseg)+2][PUbin]++; 
		}
                else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
                 NumA_phiMB4Bot[0][Lumibin]++; 
                 NumA_phiMB4Bot[1][PUbin]++; 
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

	    //denominator
            Den_theMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            Den_theMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

	    // numerator
            Num_theMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            NumA_theMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            Num_theMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            NumA_theMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
	  }
	}
        else if (NHits==3) {

	  //denominator
          Den_theMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
          Den_theMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

          int lay = missingLayer-1;

	  // is there a digi within the expected tube?

          float digiW=-1.;
          float d =1000000;
          int iex = missingLayer+3;
          Float_t expW = (*expWire)(iex);

          for (int idigi=0; idigi<Ndigis; idigi++) {

            if (digi_time->at(idigi)<320 || digi_time->at(idigi)>700) continue;

	    if (digi_wheel->at(idigi) != dtsegm4D_wheel->at(iseg)) continue;
	    if (digi_sector->at(idigi) != dtsegm4D_sector->at(iseg)) continue;
	    if (digi_station->at(idigi) != dtsegm4D_station->at(iseg)) continue;

	    if (digi_sl->at(idigi) != 2) continue;
	    if (digi_layer->at(idigi) != lay+1) continue;
            
            if (fabs(expW-digi_wire->at(idigi))<fabs(d)) {
              digiW=digi_wire->at(idigi);
	      d=expW-digiW;
	    }
	  }
   
          if ( fabs(d)< 1.1) {

	   // numerator
           Num_theMBWh[0][Lumibin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
           Num_theMBWh[1][PUbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
	  }
	}
        else {
          cout<<" what do you want?? NHits (z) = "<<NHits<<endl;
          return;
	}
       }
   }

   // computing efficiencies

   float effPhiMBWh[2][nLumiPoints][4][5]; float errPhiMBWh[2][nLumiPoints][4][5];
   float effTheMBWh[2][nLumiPoints][3][5]; float errTheMBWh[2][nLumiPoints][3][5];
   float effMB4Top[2][5][nLumiPoints];     float errMB4Top[2][5][nLumiPoints];
   float effMB4Bot[2][nLumiPoints];        float errMB4Bot[2][nLumiPoints];

   for (int ivar=0; ivar<2; ivar++) {
     for (int ipoint=0; ipoint<nLumiPoints; ipoint++) {
       for (int iwh=0; iwh<5; iwh++){
         for (int ist=0; ist<4; ist++){

	   if (Den_phiMBWh[ivar][ipoint][ist][iwh]>0.) {
            effPhiMBWh[ivar][ipoint][ist][iwh]=
             double(Num_phiMBWh[ivar][ipoint][ist][iwh])/double(Den_phiMBWh[ivar][ipoint][ist][iwh]);
            errPhiMBWh[ivar][ipoint][ist][iwh]=
	      sqrt(  double(Num_phiMBWh[ivar][ipoint][ist][iwh])
		     *(double(Den_phiMBWh[ivar][ipoint][ist][iwh])-double(Num_phiMBWh[ivar][ipoint][ist][iwh]))
		     /double(Den_phiMBWh[ivar][ipoint][ist][iwh]))
                     /double(Den_phiMBWh[ivar][ipoint][ist][iwh]);
	   }
           else {effPhiMBWh[ivar][ipoint][ist][iwh]=0.; errPhiMBWh[ivar][ipoint][ist][iwh]=0.;}

           if (ist==3) continue;

	   if (Den_theMBWh[ivar][ipoint][ist][iwh]>0.) {
            effTheMBWh[ivar][ipoint][ist][iwh]=
             double(Num_theMBWh[ivar][ipoint][ist][iwh])/double(Den_theMBWh[ivar][ipoint][ist][iwh]);
            errTheMBWh[ivar][ipoint][ist][iwh]=
	      sqrt(  double(Num_theMBWh[ivar][ipoint][ist][iwh])
		     *(double(Den_theMBWh[ivar][ipoint][ist][iwh])-double(Num_theMBWh[ivar][ipoint][ist][iwh]))
		     /double(Den_theMBWh[ivar][ipoint][ist][iwh]))
                     /double(Den_theMBWh[ivar][ipoint][ist][iwh]);
	   }
           else {effTheMBWh[ivar][ipoint][ist][iwh]=0.; errTheMBWh[ivar][ipoint][ist][iwh]=0.;}	 
	 }
         if (Den_phiMB4Top[ivar][iwh][ipoint]>0.) {
            effMB4Top[ivar][iwh][ipoint]=double(Num_phiMB4Top[ivar][iwh][ipoint])/double(Den_phiMB4Top[ivar][iwh][ipoint]);
            errMB4Top[ivar][iwh][ipoint]=
	        sqrt(  double(Num_phiMB4Top[ivar][iwh][ipoint])
	 	       *(double(Den_phiMB4Top[ivar][iwh][ipoint])-double(Num_phiMB4Top[ivar][iwh][ipoint]))
		       /double(Den_phiMB4Top[ivar][iwh][ipoint]))
                       /double(Den_phiMB4Top[ivar][iwh][ipoint]);
         }
         else {effMB4Top[ivar][iwh][ipoint]=0.; errMB4Top[ivar][iwh][ipoint]=0.;}
         if (ivar==1) cout<<" MB4Top (PU="<<PUslice[ipoint]<<") "<<effMB4Top[ivar][iwh][ipoint]
                          <<" +- "<<errMB4Top[ivar][iwh][ipoint]<<endl;
       }



       if (Den_phiMB4Bot[ivar][ipoint]>0.) {
          effMB4Bot[ivar][ipoint]=double(Num_phiMB4Bot[ivar][ipoint])/double(Den_phiMB4Bot[ivar][ipoint]);
          errMB4Bot[ivar][ipoint]=
	      sqrt(  double(Num_phiMB4Bot[ivar][ipoint])
		     *(double(Den_phiMB4Bot[ivar][ipoint])-double(Num_phiMB4Bot[ivar][ipoint]))
		     /double(Den_phiMB4Bot[ivar][ipoint]))
                     /double(Den_phiMB4Bot[ivar][ipoint]);
       }
       else {effMB4Bot[ivar][ipoint]=0.; errMB4Bot[ivar][ipoint]=0.;}
     }
   }

   ofstream results;
   std::string resultname;

   resultname.append("Results_Run2016");
   resultname = resultname + dataset + ".txt";;

   cout<<resultname<<endl;
   results.open (resultname.c_str());


   for (int ivar=0; ivar<2; ivar++) {
     for (int ipoint=0; ipoint<nLumiPoints; ipoint++) {
       for (int ist=0; ist<4; ist++){
         for (int iwh=0; iwh<5; iwh++){

           if (ivar==0) results<<lumislice[ipoint]+500.;
           else results<<PUslice[ipoint]+1.;
           results <<" "<<ist+1<<" "<<iwh-2
                   <<" "<<Den_phiMBWh[ivar][ipoint][ist][iwh]<<" "<<Num_phiMBWh[ivar][ipoint][ist][iwh]
                   <<" "<<NumA_phiMBWh[ivar][ipoint][ist][iwh]<<endl;

	 }	 
       }

       for (int iwh=0; iwh<5; iwh++){
         for (int ist=0; ist<3; ist++){

           if (ivar==0) results<<lumislice[ipoint]+500.;
           else results<<PUslice[ipoint]+1.;
           results <<" "<<ist+1<<" "<<iwh-2
                   <<" "<<Den_theMBWh[ivar][ipoint][ist][iwh]<<" "<<Num_theMBWh[ivar][ipoint][ist][iwh]
                   <<" "<<NumA_theMBWh[ivar][ipoint][ist][iwh]<<endl;
	 }
         if (ivar==0) results<<lumislice[ipoint]+500.;
         else results<<PUslice[ipoint]+1.;
         results <<" 4T"<<iwh-2<<" "<<Den_phiMB4Top[ivar][iwh][ipoint]<<" "<<Num_phiMB4Top[ivar][iwh][ipoint]
                 <<" "<<NumA_phiMB4Top[ivar][iwh][ipoint]<<endl;	 
       }

       if (ivar==0) results<<lumislice[ipoint]+500.;
       else results<<PUslice[ipoint]+1.;
       results <<" 4B   "<<Den_phiMB4Bot[ivar][ipoint]<<" "<<Num_phiMB4Bot[ivar][ipoint]
               <<" "<<NumA_phiMB4Bot[ivar][ipoint]<<endl;
     }
   }

   TGraphErrors* MB4TopYB2PU      = new TGraphErrors (nLumiPoints, &PUslice[0],  effMB4Top[1][4],&PUe[0],   errMB4Top[1][4]);
   TGraphErrors* MB4TopYB2Lumi    = new TGraphErrors (nLumiPoints, &lumislice[0],effMB4Top[0][4],&Lumie[0], errMB4Top[0][4]);

   float eff[nLumiPoints], err[nLumiPoints];

   for (int i=0; i<nLumiPoints; i++) {eff[i]=effPhiMBWh[1][i][0][4]; err[i]=errPhiMBWh[1][i][0][4];}
   TGraphErrors* MB1Wh2PU   = new TGraphErrors (nLumiPoints,&PUslice[0],  eff,         &PUe[0],   err);

   for (int i=0; i<nLumiPoints; i++) {eff[i]=effPhiMBWh[0][i][0][4]; err[i]=errPhiMBWh[0][i][0][4];}
   TGraphErrors* MB1Wh2Lumi   = new TGraphErrors (nLumiPoints,&lumislice[0],  eff,         &Lumie[0],   err);

   for (int i=0; i<nLumiPoints; i++) {eff[i]=effPhiMBWh[1][i][0][2]; err[i]=errPhiMBWh[1][i][0][2];}
   TGraphErrors* MB1Wh0PU   = new TGraphErrors (nLumiPoints,&PUslice[0],  eff,         &PUe[0],   err);

   for (int i=0; i<nLumiPoints; i++) {eff[i]=effPhiMBWh[0][i][0][2]; err[i]=errPhiMBWh[0][i][0][2];}
   TGraphErrors* MB1Wh0Lumi   = new TGraphErrors (nLumiPoints,&lumislice[0],  eff,         &Lumie[0],   err);

   for (int i=0; i<nLumiPoints; i++) {eff[i]=effTheMBWh[0][i][0][2]; err[i]=errTheMBWh[0][i][0][2];}
   TGraphErrors* MB1TheWh0Lumi   = new TGraphErrors (nLumiPoints,&lumislice[0],  eff,         &Lumie[0],   err);


   //   TCanvas *cextra  = new TCanvas(); //del
   //   hextra->Draw();
   //   hextra->SaveAs("extran.png");;


   system("mkdir plot/");
   system(("mkdir plot/"+fileName).c_str());


   TCanvas *c1  = new TCanvas();
   MB4TopYB2PU->SetTitle("MB4TopYB2PU");
   MB4TopYB2PU->SetMarkerStyle(20);
   MB4TopYB2PU->GetXaxis()->SetTitle("PU");
   MB4TopYB2PU->GetYaxis()->SetTitle("Eff.");
   MB4TopYB2PU->Draw("ap");
   c1->SaveAs(("plot/"+fileName+"MB4TopYB2PU.png").c_str());
   
   TCanvas *c2  = new TCanvas();
   MB4TopYB2Lumi->SetTitle("MB4TopYB2Lumi");
   MB4TopYB2Lumi->SetMarkerStyle(20);
   MB4TopYB2Lumi->GetXaxis()->SetTitle("Ins. Luminosity (cm^{-2}s^{-1}10^{30})");
   MB4TopYB2Lumi->GetYaxis()->SetTitle("Eff.");
   MB4TopYB2Lumi->Draw("ap");
   c2->SaveAs(("plot/"+fileName+"MB4TopYB2Lumi.png").c_str());
   
   
   TCanvas *c3  = new TCanvas();
   MB1Wh0PU->SetTitle("MB1Wh0PU");
   MB1Wh0PU->SetMarkerStyle(20);
   MB1Wh0PU->GetXaxis()->SetTitle("PU");
   MB1Wh0PU->GetYaxis()->SetTitle("Eff.");
   MB1Wh0PU->Draw("ap");
   c3->SaveAs(("plot/"+fileName+"MB1Wh0PU.png").c_str());
   
   TCanvas *c4  = new TCanvas();
   MB1Wh2PU->SetTitle("MB1Wh2PU");
   MB1Wh2PU->SetMarkerStyle(20);
   MB1Wh2PU->GetXaxis()->SetTitle("PU");
   MB1Wh2PU->GetYaxis()->SetTitle("Eff.");
   MB1Wh2PU->Draw("ap");
   c3->SaveAs(("plot/"+fileName+"MB1Wh2PU.png").c_str());
   
   TCanvas *c5  = new TCanvas();
   MB1Wh0Lumi->SetTitle("MB1Wh0Lumi");
   MB1Wh0Lumi->SetMarkerStyle(20);
   MB1Wh0Lumi->GetYaxis()->SetTitle("Eff.");
   MB1Wh0Lumi->GetXaxis()->SetTitle("Ins. Luminosity (cm^{-2}s^{-1}10^{30})");
   MB1Wh0Lumi->Draw("ap");
   c5->SaveAs(("plot/"+fileName+"MB1Wh0Lumi.png").c_str());
   
   TCanvas *c6  = new TCanvas();
   MB1Wh2Lumi->SetTitle("MB1Wh2Lumi");
   MB1Wh2Lumi->SetMarkerStyle(20);
   MB1Wh2Lumi->GetXaxis()->SetTitle("Ins. Luminosity (cm^{-2}s^{-1}10^{30})");
   MB1Wh2Lumi->GetYaxis()->SetTitle("Eff.");
   MB1Wh2Lumi->Draw("ap");
   c6->SaveAs(("plot/"+fileName+"MB1Wh2Lumi.png").c_str());

   TCanvas *c7  = new TCanvas();
   MB1TheWh0Lumi->SetTitle("MB1TheWh0Lumi");
   MB1TheWh0Lumi->SetMarkerStyle(20);
   MB1TheWh0Lumi->GetYaxis()->SetTitle("Eff.");
   MB1TheWh0Lumi->GetXaxis()->SetTitle("Ins. Luminosity (cm^{-2}s^{-1}10^{30})");
   MB1TheWh0Lumi->Draw("ap");
   c5->SaveAs(("plot/"+fileName+"MB1TheWh0Lumi.png").c_str());


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
     txtin.open("DeadList.txt",std::ifstream::in);
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
        
        if (fabs(dtsegm4D_x_dir_loc->at(iseg))>0.7) continue; // angle
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

