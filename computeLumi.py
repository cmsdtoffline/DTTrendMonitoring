#! /usr/bin/env python                                                                                                                                                              # *-* coding=utf-8 *-*

import sys
import subprocess
import json
from optparse import OptionParser
from operator import itemgetter
import collections
import os

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option("-e", "--ext", dest="extFile",
                  default="",
help="Use an external file")

parser.add_option("-y", "--year", dest="year",
                  default="2018",
help="change year, deafault is 2018")

parser.add_option("-r", "--run", dest="run",
                  default="",
help="Use it to calculte luminosity of a single run")

parser.add_option("-v", "--verbose", dest="verbose",
                  default=0,
help="Show more print out")

parser.add_option("-l", "--list", dest="runList",
                  default="",
help="file list of runs") # to be added


(options, args) = parser.parse_args()
extFile = options.extFile
year    = options.year
run     = options.run
verbose = options.verbose
runList = options.runList

sys.path.append('./Tools')

import lumi_utils as lu

# 2016 Dataset dictionary
Dataset2016 = {"Run2016B": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2016/Run2016BZMu23Sep2016-v1.root",
               "Run2016C": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2016/Run2016CZMu23Sep2016-v1.root",
               "Run2016D": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2016/Run2016DZMu23Sep2016-v1.root",
               "Run2016E": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2016/Run2016EZMu23Sep2016-v1.root",
               "Run2016F": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2016/Run2016FZMu23Sep2016-v1.root",
               "Run2016G": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2016/Run2016GZMu23Sep2016-v1.root",
               "Run2016H": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2016/Run2016HZMuPromptReco-v2.root"
               }

# 2017 Dataset dictionary
Dataset2017 = {"Run2017A": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017A.root ",
               "Run2017B": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017B_94X.root",
               "Run2017C": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017C_94X.root",
               "Run2017D": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017D_94X.root",
               "Run2017E": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017E_94X.root",
               "Run2017F": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017F_94X.root"
               }


# 2018 Dataset dictionary
Dataset2018 = {"Run2018A": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2018/Prompt/DTTree_SingleMuon_Run2018A_Prompt_v1v2_ZMu_GOLDEN.root",
               "Run2018B": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2018/Prompt/DTTree_SingleMuon_Run2018B_Prompt_v1v2_ZMu_GOLDEN.root",
               "Run2018C": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2018/Prompt/DTTree_SingleMuon_Run2018C_Prompt_v1v2v3_ZMu_GOLDEN.root"
               }


TotLumiRun1      = 29445 *1e-3   # check fix. Find delivered and recorded
TotLumiRun1Pbp   = 0.036 *1e-3

TotLumipp2015    = 4220    *1e-3 # from 03/06/15 to 03/111/15
TotLumippRef2015 = 28.82   *1e-3 # from 19/11/15 to 23/11/15
TotLumiPbPb2015  = 0.0006  *1e-3 # from 25/11/15 to 13/12/15

TotLumipp2016    = 41441   *1e-3 # from 22/04/16 to 27/07/216 //OK
TotLumiPbp2016   = 0.188   *1e-3 # from 18/11/16 to 02/12/16

TotLumipp2017    = 51656   *1e-3 # from 30/05/17 to 26/11/17 //OK
TotLumippRef2017 = 334.27  *1e-3 # from 11/11/17 to 21/11/17 ?

TotLumipp2018    = 68180   *1e-3 # from 17/04/18 to /11/18 # check with brilcalc
TotLumiPbPb2018  = 0.00056 *1e-3 # from 08/11/18 to 20/11/18 


def getLumis(inputFile):

    RunLumis = lu.RunLumis({})
    vals = []
    vals.append(lu.getRunLumis(fname))
    allRunLumis = sum(vals, RunLumis)
    return allRunLumis



def writeLumi(run,outFileTotal):
    subprocess.call(["brilcalc", "lumi", "-u/fb", "-o","temp.csv", "--begin",firstRun,"--end",str(run)])
#    print "brilcalc", "lumi", "-u/fb", "-o","temp.csv", "--begin",firstRun,"--end",str(run)

    fileTemp  = open("temp.csv", "r")
    runDic = {}
    mode   = 'r+' if os.path.exists(outFileTotal) else 'w+'

    with open(outFileTotal, mode ) as outjson:
        line = subprocess.check_output(['tail', '-1', "temp.csv"])

        try:
            runDic = json.load(outjson)
        except ValueError:
            print "Creating new file"

        if str(run) in runDic: 
            print "Run",run,"already in the file"
            return         
#            continue

        outjson.seek(0)  #should reset file position to the beginning.

        runDic[run]    = line.split(",")

        if runDic[run][0]=="#run:fill": sys.exit("Error: Something wrong with the selected run")

        runDic[run][0] = runDic[run][0].replace("#","") 
        runDic[run][5] = runDic[run][5].replace("\n","") 

#        print "Lumi pre",runDic[run][4]
        
        if  year=="2018": 
            runDic[run][4] = str(float(runDic[run][4]) + TotLumiRun1 + TotLumipp2015 +TotLumippRef2015 +TotLumiPbPb2015 + TotLumipp2016 +TotLumiPbp2016 + TotLumipp2017 + TotLumippRef2017)
 #           print TotLumiRun1 + TotLumipp2015 +TotLumippRef2015 +TotLumiPbPb2015 + TotLumipp2016 +TotLumiPbp2016 + TotLumipp2017 + TotLumippRef2017
    
        elif year=="2017": 
            runDic[run][4] = str(float(runDic[run][4]) + TotLumiRun1 + TotLumipp2015 +TotLumippRef2015 +TotLumiPbPb2015 + TotLumipp2016 +TotLumiPbp2016)
#            print TotLumiRun1 + TotLumipp2015 +TotLumippRef2015 +TotLumiPbPb2015 + TotLumipp2016 +TotLumiPbp2016

        elif year=="2016": 
            runDic[run][4] = str(float(runDic[run][4]) + TotLumiRun1 + TotLumipp2015 +TotLumippRef2015 +TotLumiPbPb2015)

        elif year=="2015": 
            runDic[run][4] = str(float(runDic[run][4]) + TotLumiRun1)

        if verbose: print "Run",run,"Lumi = ",runDic[run][4]

        #runDic[run][5] = str(float(runDic[run][5])+TotLumiRun1) To be fixed if recorded is needed

        #json.dumps('"comments":["nfill","nrun","nls","ncms","totdelivered(/pb)","totrecorded(/pb)"]',outjson, indent=4) #fixme
        json.dump(runDic,outjson, indent=4) 
        subprocess.call(["rm", "temp.csv"])
        #outjson.truncate() 


if __name__ == '__main__':


    if extFile=="" and run =="" and runList == "":
        if   year=="2016": Dataset = Dataset2016
        elif year=="2017": Dataset = Dataset2017
        elif year=="2018": Dataset = Dataset2018
        else: sys.exit("Wrong year")

    elif extFile!="": Dataset = {"singlefile":extFile}  


    # print "#################################"


    jsonFilesList = []

    # dictionary are not ordered, so the loop is randomly distributed on the keys
    # for k, fname in Dataset2016.items():
    # The following looping method guaranitee and ordered loop onto the keys

    if   year=="2015": firstRun = "247607" # first run 2015A
    elif year=="2016": firstRun = "271036" # first run 2016A
    elif year=="2017": firstRun = "294645" # first run 2017A
    elif year=="2018": firstRun = "314472" # Collisions2018Commiss
    else: sys.exit("Wrong year")

    path = "data/"
    outFileTotal = path + 'IntLumi/Total'+year+'.json'

    if run=="":
        print " Producing Lumi JSON starting from DT DQM TTree  "
        for key in sorted(Dataset):
            print "running ",key
            fname = Dataset[key]
        # Class define in lumi_utils that handles Runs and Lumi Blocks
            datasetRunLumis = getLumis(fname)
            
            for k, val in datasetRunLumis.getJson().items():
                writeLumi(k,outFileTotal)
                
    elif run!="":
        print "hei"
        writeLumi(run,outFileTotal)

    elif runList!="":
        try:
            fRunList = open(runList,"r")
            try:
                lines = fRunList.readlines()
                for line in line:
                    runNumber = line.split(",")
                    writeLumi(runNumber,outFileTotal)
            finally:
                fRunList.close()        
        except IOError:
            pass
    else:
        print "Something wrong with options. Type computeLumi_2.py -h for help"
