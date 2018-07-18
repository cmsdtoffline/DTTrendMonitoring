#! /usr/bin/env python                                                                                                                                                                                                                       # *-* coding=utf-8 *-*

# from ROOT import *
# import datetime
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
                  default="2017",
help="change year, deafault is 2017")

parser.add_option("-n", "--name", dest="name",
                  default="",
help="choose name when use single file")


(options, args) = parser.parse_args()
extFile = options.extFile
year    = options.year
name    = options.name
# import json
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

Dataset2017 = {"Run2017A": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017A.root ",
               "Run2017B": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017B_94X.root",
               "Run2017C": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017C_94X.root",
               "Run2017D": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017D_94X.root",
               "Run2017E": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017E_94X.root",
               "Run2017F": "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2017/DTTree_ZMuSkim2017F_94X.root"
               }



def getLumis(inputFile):

    RunLumis = lu.RunLumis({})
    vals = []
   # print "lu",lu.getRunLumis(fname),type(lu.getRunLumis(fname))
    vals.append(lu.getRunLumis(fname))
    allRunLumis = sum(vals, RunLumis)
    print "Vals",vals,type((vals.__getitem__(0)).__str__())
    return allRunLumis

if __name__ == '__main__':

    if extFile=="":
        if   year=="2016": Dataset = Dataset2016
        elif year=="2017": Dataset = Dataset2017
        else: sys.exit("Wrong year")

    else:  Dataset = {name:extFile}

    # print "#################################"
    print "   Producing Lumi JSON starting from DT DQM TTree   "

    jsonFilesList = []
    # dictionary are not ordered, so the loop is randomly distributed on the keys
    # for k, fname in Dataset2016.items():
    # The following looping method guaranitee and ordered loop onto the keys


    if year=="2016":   firstRun = "271036"
    elif year=="2017": firstRun = "294645"
    elif year=="2018": firstRun = "314472" # precomissioning
    else: sys.exit("Wrong year")

    for key in sorted(Dataset):
        fname = Dataset[key]
        # Class define in lumi_utils that handles Runs and Lumi Blocks
        datasetRunLumis = getLumis(fname)

        path = "data/"

        outFileName = "DT_" + key + ".json"
        outFileCSV  = "DT_" + key + ".csv"

        outFileName = path + "RunJson/" + outFileName
        outFileCSV  = path + "IntLumi/" + outFileCSV

        outFileTotal = path + 'IntLumi/Total'+year+'.json'

        jsonFilesList.append(outFileName)
        print "outFileName",outFileName
        datasetRunLumis.writeToFile(outFileName)
        print "\t SAVED!! - {:s}".format(fname)

        ### Run brilcalc
        subprocess.call(["brilcalc", "lumi", "-u/pb", "-o"+outFileCSV, "-i"+outFileName])

        for k, val in datasetRunLumis.getJson().items():
            subprocess.call(["brilcalc", "lumi", "-u/pb", "-o","temp.csv", "--begin",firstRun,"--end",str(k)])

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

                runDic[k] = line.split(",")
                runDic[k][0] = runDic[k][0].replace("#","") 
                runDic[k][5] = runDic[k][5].replace("\n","") 

#                json.dumps('"comments":["nfill","nrun","nls","ncms","totdelivered(/pb)","totrecorded(/pb)"]',outjson, indent=4) //fixme
                json.dump(runDic,outjson, indent=4)

                #outjson.truncate() 


