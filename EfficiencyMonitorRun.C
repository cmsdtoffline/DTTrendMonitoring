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


int lessentries = 100;


void EfficiencyMonitor::PreLoop()
{

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  nentries=lessentries;

  Long64_t nbytes = 0, nb = 0;

  ofstream DeadList;
  std::string deadname;

  deadname.append("DeadList_Run2016");
  //deadname.append("DeadList_run");
  deadname = deadname + dataset + ".txt";;
  //  deadname = deadname +std::to_string(runnumber)+".txt";;

  cout<<deadname<<endl;
  DeadList.open (deadname.c_str());

  char go;
  char hname[50];
  TH1F* occupancy[5][14][4][3][4];
  for (int iwh=0; iwh<5; iwh++){
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

    if (jentry%1000 == 0) cout<<" Pre-Loop evento "<<jentry<<endl;

    for (int idigi=0; idigi<Ndigis; idigi++) {

      occupancy[digi_wheel->at(idigi)+2][digi_sector->at(idigi)-1][digi_station->at(idigi)-1]
      [digi_sl->at(idigi)-1][digi_layer->at(idigi)-1]->Fill(float(digi_wire->at(idigi)));

    }
  }

  // analyze occupancy histos and fill dead channel table

  int nwire=0; int NwireTot=0;
  for (int iw=0; iw<5000; iw++) for (int geo=0; geo<6; geo++) dead[iw][geo]=0;

  for (int iwh=0; iwh<5; iwh++){ // Loop on wheels BEGIN
    for (int ise=0; ise<14; ise++) { // Loop on sectors BEGIN
      for (int ist=0; ist<4; ist++) { // Loop on stations BEGIN
        if (ist!=3 && ise>11) continue;
        for (int isl=0; isl<3; isl++) { // Loop on super layers BEGIN
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

          for (int ilay=0; ilay<4; ilay++) { // Loop on layers BEGIN
            for (int iw=1; iw<nwire+1; iw++) { // Loop on wires BEGIN
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
            } // Loop on wires END
          } // Loop on layers END
        } // Loop on super layers END
      } // Loop on stations
    } // Loop on sectors
  } // Loop on Wheels

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

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  char go;

  nentries=lessentries;

  Long64_t nbytes = 0, nb = 0;

  // COUNTERS:

  int Num_phiMBWh[2][22][4][5];  int Den_phiMBWh[2][22][4][5];  // 2 variabili (run, time), 22 punti, 4 stazioni, 5 ruote
  int Num_theMBWh[2][22][3][5];  int Den_theMBWh[2][22][3][5];  // 2 variabili (run, time), 22 punti, 3 stazioni, 5 ruote
  int NumA_phiMBWh[2][22][4][5]; int NumA_theMBWh[2][22][3][5]; // 'A' stands for 'Associated'

  int Num_phiMB4Top[2][5][22];   int Den_phiMB4Top[2][5][22];   // 2 variabili (run, rime), 22 punti
  int Num_phiMB4Bot[2][22];      int Den_phiMB4Bot[2][22];      // 2 variabili (run, time), 22 punti
  int NumA_phiMB4Top[2][5][22];  int NumA_phiMB4Bot[2][22];     // 'A' stands for 'Associated'

  for (int ivar=0; ivar<2; ivar++){
    for (int ipoint=0; ipoint<22; ipoint++) {
      for (int iwh=0; iwh<5; iwh++){
        for (int ist=0; ist<4; ist++){

          Num_phiMBWh[ivar][ipoint][ist][iwh]=0;
          NumA_phiMBWh[ivar][ipoint][ist][iwh]=0;
          Den_phiMBWh[ivar][ipoint][ist][iwh]=0;

          if (ist==3) continue;
          Num_theMBWh[ivar][ipoint][ist][iwh]=0;
          NumA_theMBWh[ivar][ipoint][ist][iwh]=0;
          Den_theMBWh[ivar][ipoint][ist][iwh]=0;
        }
        Num_phiMB4Top[ivar][iwh][ipoint]=0;
        NumA_phiMB4Top[ivar][iwh][ipoint]=0;
        Den_phiMB4Top[ivar][iwh][ipoint]=0;

      }
      Num_phiMB4Bot[ivar][ipoint]=0;
      NumA_phiMB4Bot[ivar][ipoint]=0;
      Den_phiMB4Bot[ivar][ipoint]=0;
    }
  }

  // DEAD CHANNELS (to skip):

  cout<<" within Loop: Ndead "<<Ndead<<endl;

  if (Ndead==0) {
    cout<<" READING FROM FILE !! "<<endl;
    ifstream txtin;
    std::string deadname;

    deadname.append("DeadList_Run2016");
    deadname = deadname + dataset + ".txt";

    cout<<" reading from "<<deadname.c_str()<<endl;
    //ifstream txtin(deadname.c_str(),std::ifstream::in);
    txtin.open(deadname.c_str(),std::ifstream::in);

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

  // run and time bins:

  /*------------------------
               From     To
  --------------------------
  Run2016A |	271036	271658
  Run2016B |	272007	275376
  Run2016C |	275657	276283
  Run2016D |	276315	276811
  Run2016E |	276831	277420
  Run2016F |	277772	278808
  Run2016G |	278820	280385
  Run2016H |	280919	284044
  ---------------------------*/


  std::set<Int_t> runNumber_Set; // ordered list of runnumbers

  /* Loop used to collect the list of runs containend inside the TTree*/
  for (Long64_t jentry=2; jentry<nentries;jentry++) {

    // Getting only runNumber values
    b_runnumber->GetEntry(jentry);
    // std::set store no duplicates of the values inserted. We will end with an ordered list of runs
    runNumber_Set.insert(runnumber);

  }

  /* - container::RunID
       is a container that stores the runNUmber plus the Numerators
       and Denominator, reletated to that runNumber.
     - std::map<Int_t,container::RunID>
       is a collection of all the container::RunID mapped by the runNumber related*/
  std::map<Int_t,container::RunID> mapRuns;
  for ( auto& run : runNumber_Set ){
    container::RunID tempRun(run);
    mapRuns[run] = tempRun;
  }

  int firstrun= *(runNumber_Set.begin());
  int lastrun= *(runNumber_Set.rbegin());

  float runslice[22]; float rune[22];
  //if (dataset=='G') firstrun=278800;
  int runRange = lastrun - firstrun;
  int numBins = 21;
  float binWidth = (float)runRange / (float)numBins;

  if (dataset == 'H' ) firstrun = 282500;
  for (int it=0; it<numBins+1; it++) {
    runslice[it]= firstrun+it*binWidth;
    rune[it]=28.;
  }


  for (Long64_t jentry=2; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // if (Cut(ientry) < 0) continue;

    if (jentry%1000 == 0) cout<<" evento "<<jentry<< "Progress:  "<< ((float)jentry*100.)/(float)nentries<<" \% "<<endl;

    // identify Runbin and Timebin
    int Runbin=-1;
    for (int islice=0; islice<22; islice++) {
      if ( runnumber >= runslice[islice]   &&
        runnumber  < runslice[islice+1]    ) Runbin=islice;
      }

      if (Runbin<0) {
        cout<<" run number out of range!! "<<runnumber<<endl;
        continue;
      }

      // First search for Phi segments

      for (int iseg=0; iseg<Ndtsegments; iseg++) {

        //selection
        if (!dtsegm4D_hasPhi->at(iseg)) continue;
        if (dtsegm4D_station->at(iseg)!=4 && !dtsegm4D_hasZed->at(iseg) ) continue;

        int seg_phinhits = dtsegm4D_phinhits->at(iseg);
        if (fabs(dtsegm4D_x_dir_loc->at(iseg))>0.7) continue; // angle

        TVectorF *expWire=(TVectorF*)dtsegm4D_hitsExpWire->At(iseg);

        // If a hit a missing, let us check that the extrapolation doesn't fall beyond layer or cross a dead cell!
        int NexpDead=0; bool OutOfLayer=false;

        if (seg_phinhits < 8 ) {
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
            if (NexpDead>MaxDead) break;
          }
          if (OutOfLayer)       continue; // this segment goes out of layer boundary
          if (NexpDead>MaxDead) continue; // this segment crosses dead cell(s): drop it!
        }

        int NHits=0; int missingLayer[3][2]; for (int imi=0; imi<3;imi++) {missingLayer[imi][0]=0;missingLayer[imi][1]=0;}
        int nmissing=0;

        TVectorF *hitLayerPhi=(TVectorF*)dtsegm4D_phi_hitsLayer->At(iseg);
        TVectorF *hitSuperLayerPhi=(TVectorF*)dtsegm4D_phi_hitsSuperLayer->At(iseg);
        TVectorF *hitWirePhi=(TVectorF*)dtsegm4D_phi_hitsWire->At(iseg);

        for (int ilay=1; ilay<9; ilay++) {
          // Search for associated hits
          bool foundh=false;
          for (int kk=0; kk<seg_phinhits; kk++) {

            int sl1 = (*hitSuperLayerPhi)(kk);
            int lay1 = (sl1==1) ?  (*hitLayerPhi)(kk) : (*hitLayerPhi)(kk)+4; // hit layer 1-8

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

            Num_phiMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            NumA_phiMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            Den_phiMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc

            // here you fill the container::RunID for the runnumber considered
            mapRuns[runnumber].Num_phiMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            mapRuns[runnumber].NumA_phiMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            mapRuns[runnumber].Den_phiMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;


            if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
              Num_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Runbin]++;
              NumA_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Runbin]++;
              Den_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Runbin]++;
              // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
            }
            else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
              Num_phiMB4Bot[0][Runbin]++;
              NumA_phiMB4Bot[0][Runbin]++;
              Den_phiMB4Bot[0][Runbin]++;
              // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
            }
          }
        }
        else { // let's see how to treat missing layers

        for (int imiss=0; imiss<nmissing; imiss++) {

          int sl = missingLayer[imiss][0] < 5 ? 0 : 1;
          int lay = sl==0 ? missingLayer[imiss][0]-1 : missingLayer[imiss][0]-5;

          // is there a digi within the expected tube?

          float digiW=-1.;
          float d =1000000;
          int iex = missingLayer[imiss][0] < 5 ? missingLayer[imiss][0]-1 : missingLayer[imiss][0]+3;
          Float_t expW = (*expWire)(iex);

          for (int idigi=0; idigi<Ndigis; idigi++) {

            if (digi_time->at(idigi)<320 || digi_time->at(idigi)>700) continue;

            if (digi_wheel->at(idigi) != dtsegm4D_wheel->at(iseg)) continue;
            if (digi_sector->at(idigi) != dtsegm4D_sector->at(iseg)) continue;
            if (digi_station->at(idigi) != dtsegm4D_station->at(iseg)) continue;

            if (digi_sl->at(idigi) == 2) continue;
            if (digi_sl->at(idigi) == 1 && sl != 0) continue;
            if (digi_sl->at(idigi) == 3 && sl != 1) continue;

            if (digi_layer->at(idigi) != lay+1) continue;

            if (fabs(expW-digi_wire->at(idigi))<fabs(d)) {
              digiW=digi_wire->at(idigi);
              d=expW-digiW;
            }
          }

          if ( fabs(d)< 1.1) {missingLayer[imiss][1]=1; }
        }

        if (NHits==nrequiredhit) {

          for (int imiss=0; imiss<nmissing; imiss++) {

            int sl = missingLayer[imiss][0] < 5 ? 0 : 1;
            int lay = sl==0 ? missingLayer[imiss][0]-1 : missingLayer[imiss][0]-5;

            //denominator

            Den_phiMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc

            // here you fill the container::RunID for the runnumber considered
            mapRuns[runnumber].Den_phiMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

            if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
              Den_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Runbin]++;
              // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
            }
            else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
              Den_phiMB4Bot[0][Runbin]++;
              // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
            }

            if (missingLayer[imiss][1]) {
              // numerator
              Num_phiMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
              // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc

              // here you fill the container::RunID for the runnumber considered
              mapRuns[runnumber].Num_phiMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

              if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
                Num_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Runbin]++;
                // here you mat want to increment counters for a second variable (time??) e.g. using [1]  etc
              }
              else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
                Num_phiMB4Bot[0][Runbin]++;
                // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
              }
            }
          }
        }

        else {

          for (int sl=0; sl<2; sl++) for (int lay=0; lay<4; lay++) {

            //denominator
            Den_phiMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
            // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc

            // here you fill the container::RunID for the runnumber considered
            mapRuns[runnumber].Den_phiMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

            if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
              Den_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Runbin]++;
              // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
            }
            else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
              Den_phiMB4Bot[0][Runbin]++;
              // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
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
              Num_phiMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
              // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc

              // here you fill the container::RunID for the runnumber considered
              mapRuns[runnumber].Num_phiMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

              if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
                Num_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Runbin]++;
                // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
              }
              else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
                Num_phiMB4Bot[0][Runbin]++;
                // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
              }
            }
            if (!(missAss)) {
              // numerator
              NumA_phiMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
              // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc

              // here you fill the container::RunID for the runnumber considered
              mapRuns[runnumber].NumA_phiMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

              if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==4 || dtsegm4D_sector->at(iseg)==13)) {
                NumA_phiMB4Top[0][dtsegm4D_wheel->at(iseg)+2][Runbin]++;
                // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
              }
              else if (dtsegm4D_station->at(iseg)==4 && (dtsegm4D_sector->at(iseg)==10 || dtsegm4D_sector->at(iseg)==14)) {
                NumA_phiMB4Bot[0][Runbin]++;
                // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc
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

      TVectorF *hitLayerZ=(TVectorF*)dtsegm4D_z_hitsLayer->At(iseg);
      TVectorF *hitWireZ=(TVectorF*)dtsegm4D_z_hitsWire->At(iseg);

      for (int ilay=1; ilay<5; ilay++) {
        // Search for associated hits
        bool foundh=false;
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
          Den_theMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
          // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc

          // here you fill the container::RunID for the runnumber considered
          mapRuns[runnumber].Den_theMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

          // numerator
          Num_theMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
          NumA_theMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
          // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc

          // here you fill the container::RunID for the runnumber considered
          mapRuns[runnumber].Num_theMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
          mapRuns[runnumber].NumA_theMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
        }
      }
      else if (NHits==3) {

        //denominator
        Den_theMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
        // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc

        // here you fill the container::RunID for the runnumber considered
        mapRuns[runnumber].Den_theMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

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
          Num_theMBWh[0][Runbin][dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;
          // here you mat want to increment counters for a second variable (time??) e.g. using [1][Timebin]  etc

          // here you fill the container::RunID for the runnumber considered
          mapRuns[runnumber].Num_theMBWh[dtsegm4D_station->at(iseg)-1][dtsegm4D_wheel->at(iseg)+2]++;

        }
      }
      else {
        cout<<" what do you want?? NHits (z) = "<<NHits<<endl;
        return;
      }
    }
  }

  // computing efficiencies

  float effPhiMBWh[2][22][4][5]; float errPhiMBWh[2][22][4][5];
  float effTheMBWh[2][22][3][5]; float errTheMBWh[2][22][3][5];
  float effMB4Top[2][5][22];     float errMB4Top[2][5][22];
  float effMB4Bot[2][22];        float errMB4Bot[2][22];

  for (int ivar=0; ivar<1; ivar++) {
    for (int ipoint=0; ipoint<22; ipoint++) {
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
        if (ivar==0) cout<<" MB4Top (Run="<<runslice[ipoint]<<") "<<effMB4Top[ivar][iwh][ipoint]
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
  ofstream results_csv;
  std::string resultname;
  std::string resultname_csv;

  resultname.append("ResultsRun_Run2016");
  resultname_csv.append("ResultsRun_Run2016");
  resultname = resultname + dataset + ".txt";
  resultname_csv = resultname_csv + dataset + ".csv";

  cout<<resultname<<endl;
  cout<<resultname_csv<<endl;
  results.open (resultname.c_str());
  results_csv.open (resultname_csv.c_str());

  results<<"Run,Station,Wheel,Den,Num,NumA"<<endl;
  results_csv<<"Run,Station,Wheel,Den,Num,NumA"<<endl;

  for ( auto& run: runNumber_Set){
    for (int ist=0; ist<4; ist++){
      for (int iwh=0; iwh<5; iwh++){

        results_csv<<mapRuns[run].runNumber;//+50;
        results_csv <<","<<ist+1<<","<<iwh-2
        <<","<<mapRuns[run].Den_phiMBWh[ist][iwh]<<","<<mapRuns[run].Num_phiMBWh[ist][iwh]
        <<","<<mapRuns[run].NumA_phiMBWh[ist][iwh]<<endl;

      }
    }
  }
  //results<<"\n ============================================================ \n"<<endl;
  results<<"RunBin,Station,Wheel,Den,Num,NumA"<<endl;
  for (int ivar=0; ivar<1; ivar++) {
    for (int ipoint=0; ipoint<22; ipoint++) {
      for (int ist=0; ist<4; ist++){
        for (int iwh=0; iwh<5; iwh++){

          if (ivar==0) results<<runslice[ipoint];//+50;
          else cout<<" no data to write out for second variable, as yet!! "<<endl;
          results <<" "<<ist+1<<" "<<iwh-2
          <<" "<<Den_phiMBWh[ivar][ipoint][ist][iwh]<<" "<<Num_phiMBWh[ivar][ipoint][ist][iwh]
          <<" "<<NumA_phiMBWh[ivar][ipoint][ist][iwh]<<endl;

        }
      }

      for (int iwh=0; iwh<5; iwh++){
        for (int ist=0; ist<3; ist++){

          if (ivar==0) results<<runslice[ipoint];//+50;
          else  cout<<" no data to write out for second variable, as yet!! "<<endl;
          results <<" "<<ist+1<<" "<<iwh-2
          <<" "<<Den_theMBWh[ivar][ipoint][ist][iwh]<<" "<<Num_theMBWh[ivar][ipoint][ist][iwh]
          <<" "<<NumA_theMBWh[ivar][ipoint][ist][iwh]<<endl;
        }
        if (ivar==0) results<<runslice[ipoint];//+50;
        else cout<<" no data to write out for second variable, as yet!! "<<endl;
        results <<" 4T"<<iwh-2<<" "<<Den_phiMB4Top[ivar][iwh][ipoint]<<" "<<Num_phiMB4Top[ivar][iwh][ipoint]
        <<" "<<NumA_phiMB4Top[ivar][iwh][ipoint]<<endl;
      }

      if (ivar==0) results<<runslice[ipoint];//+50.;
      else cout<<" no data to write out for second variable, as yet!! "<<endl;
      results <<" 4B   "<<Den_phiMB4Bot[ivar][ipoint]<<" "<<Num_phiMB4Bot[ivar][ipoint]
      <<" "<<NumA_phiMB4Bot[ivar][ipoint]<<endl;
    }
  }


  TGraphErrors* MB4TopYB2Run    = new TGraphErrors (22,runslice,effMB4Top[0][4],rune, errMB4Top[0][4]);

  float eff[22], err[22];

  for (int i=0; i<22; i++) {eff[i]=effPhiMBWh[0][i][0][4]; err[i]=errPhiMBWh[0][i][0][4];}
  TGraphErrors* MB1Wh2Run   = new TGraphErrors (22,runslice,  eff,         rune,   err);

  for (int i=0; i<22; i++) {eff[i]=effPhiMBWh[0][i][0][2]; err[i]=errPhiMBWh[0][i][0][2];}
  TGraphErrors* MB1Wh0Run   = new TGraphErrors (22,runslice,  eff,         rune,   err);


  TCanvas *c1  = new TCanvas();
  MB4TopYB2Run->SetTitle("MB4TopYB2Run");
  MB4TopYB2Run->SetMarkerStyle(20);
  MB4TopYB2Run->Draw("ap");
  c1->SaveAs("MB4.png");


  TCanvas *c2  = new TCanvas();

  MB1Wh0Run->SetTitle("MB1Wh0Run");
  MB1Wh0Run->SetMarkerStyle(20);
  MB1Wh0Run->Draw("ap");
  c2->SaveAs("MB1Wh0.png");

  TCanvas *c3  = new TCanvas();

  MB1Wh2Run->SetTitle("MB1Wh2Run");
  MB1Wh2Run->SetMarkerStyle(20);
  MB1Wh2Run->Draw("ap");
  c3->SaveAs("MB1Wh2.png");

  results.close();
  results_csv.close();
}

void EfficiencyMonitor::PostLoop()
{
  // check how many cells were not crossed by selected segments (efficiency undetermined)


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  nentries=lessentries;

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
      TVectorF *expWire=(TVectorF*)dtsegm4D_hitsExpWire->At(iseg);

      // If a hit a missing, let us check that the extrapolation doesn't cross a dead cell!
      int NexpDead=0;

      if (dtsegm4D_phinhits->at(iseg) < 8 ) {

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
                  if (dead[idead][0] != iwh-2) continue;
                  if (dead[idead][1] != ise+1) continue;
                  if (dead[idead][2] != ist+1) continue;
                  if (dead[idead][3] != isl+1) continue;
                  if (dead[idead][4] != ilay+1) continue;
                  if (dead[idead][5] != iw) continue;
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
