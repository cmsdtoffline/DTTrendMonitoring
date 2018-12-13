#! /usr/bin/env python                                                                                                                                                              # *-* coding=utf-8 *-*

import sys
import subprocess
import json
from optparse import OptionParser
import collections
import os

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option("-e", "--ext", dest="extFile",
                  default="",
help="Use an externale file")

parser.add_option("-y", "--year", dest="year",
                  default="2018",
help="change year, deafault is 2018")

(options, args) = parser.parse_args()
extFile = options.extFile
year    = options.year


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

TotLumipp2016    = 40820   *1e-3 # from 22/04/16 to 27/07/216
TotLumiPbp2016   = 0.188   *1e-3 # from 18/11/16 to 02/12/16

TotLumipp2017    = 49790   *1e-3 # from 30/05/17 to 26/11/17
TotLumippRef2017 = 334.27  *1e-3 # from 11/11/17 to 21/11/17 ?

TotLumipp2018    = 68180   *1e-3 # from 17/04/18 to /11/18 
TotLumiPbPb2018  = 0.00056 *1e-3 # from 08/11/18 to 20/11/18 

def getLumis(inputFile):

    RunLumis = lu.RunLumis({})
    vals = []

    vals.append(lu.getRunLumis(fname))
    allRunLumis = sum(vals, RunLumis)
    return allRunLumis

if __name__ == '__main__':

    if extFile=="":
        if   year=="2016": Dataset = Dataset2016
        elif year=="2017": Dataset = Dataset2017
        elif year=="2018": Dataset = Dataset2018
        else: sys.exit("Wrong year")

    else:  Dataset = {"singlefile":extFile}

    # print "#################################"
    print "   Producing Lumi JSON starting from DT DQM TTree   "

    jsonFilesList = []

    # dictionary are not ordered, so the loop is randomly distributed on the keys
    # for k, fname in Dataset2016.items():
    # The following looping method guaranitee and ordered loop onto the keys

    if   year=="2015": firstRun = "247607" # first run 2015A
    elif year=="2016": firstRun = "271036" # first run 2016A
    elif year=="2017": firstRun = "294645" # first run 2017A
    elif year=="2018": firstRun = "314472" # Collisions2018Commiss
    else: sys.exit("Wrong year")

    for key in sorted(Dataset):
        print "running ",key
        fname = Dataset[key]
        # Class define in lumi_utils that handles Runs and Lumi Blocks
        datasetRunLumis = getLumis(fname)

        path = "data/"
        outFileTotal = path + 'IntLumi/Total'+year+'.json'

        for k, val in datasetRunLumis.getJson().items():

            subprocess.call(["brilcalc", "lumi", "-u/fb", "-o","temp.csv", "--begin",firstRun,"--end",str(k)])

            fileTemp  = open("temp.csv", "r")
            runDic = {}
            mode   = 'r+' if os.path.exists(outFileTotal) else 'w+'

            with open(outFileTotal, mode ) as outjson:
                line = subprocess.check_output(['tail', '-1', "temp.csv"])

                try:
                    runDic = json.load(outjson)
                except ValueError:
                    print "Creating new file"

                if str(k) in runDic: 
                    print "Run",k,"already in the file"
                    continue

                outjson.seek(0)  #should reset file position to the beginning.

                runDic[k]    = line.split(",")

                if runDic[k][0]=="#run:fill": sys.exit("Error: Something wrong with the selected run")

                runDic[k][0] = runDic[k][0].replace("#","") 
                runDic[k][5] = runDic[k][5].replace("\n","") 
                
                if   year=="2018": 
                    runDic[k][4] = str(float(runDic[k][4]) + TotLumiRun1 + TotLumipp2015 +TotLumippRef2015 +TotLumiPbPb2015 + TotLumipp2016 +TotLumiPbp2016 + TotLumipp2017 + TotLumippRef2017)

                elif year=="2017": 
                    runDic[k][4] = str(float(runDic[k][4]) + TotLumiRun1 + TotLumipp2015 +TotLumippRef2015 +TotLumiPbPb2015 + TotLumipp2016 +TotLumiPbp2016)

                elif year=="2016": 
                    runDic[k][4] = str(float(runDic[k][4]) + TotLumiRun1 + TotLumipp2015 +TotLumippRef2015 +TotLumiPbPb2015)

                elif year=="2015": 
                    runDic[k][4] = str(float(runDic[k][4]) + TotLumiRun1)


                #runDic[k][5] = str(float(runDic[k][5])+TotLumiRun1) To be fixed if recorded is needed

                json.dumps('"comments":["nfill","nrun","nls","ncms","totdelivered(/pb)","totrecorded(/pb)"]',outjson, indent=4) #fixme
                json.dump(runDic,outjson, indent=4)
                subprocess.call(["rm", "temp.csv"])
                #outjson.truncate() 


