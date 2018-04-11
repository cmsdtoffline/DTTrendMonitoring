# *-* coding=utf-8 *-*

# from ROOT import *
# import datetime
import sys
import subprocess
from optparse import OptionParser


parser = OptionParser(usage="usage: %prog [options]")

parser.add_option("-e", "--ext", dest="extFile",
                  default="",
help="Use an externale file")

parser.add_option("-y", "--year", dest="year",
                  default="2016",
help="change year, deafault is 2016")

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
    vals.append(lu.getRunLumis(fname))
    allRunLumis = sum(vals, RunLumis)
    return allRunLumis


if __name__ == '__main__':

    print "extFile",extFile

    if extFile=="":
        if   year=="2016": Dataset = Dataset2016
        elif year=="2017": Dataset = Dataset2017
        else: sys.exit("Wrong year")

    else:  Dataset = {name:extFile}

    # print "#################################"
    print "    Producing Lumi JSON starting from DT DQM TTree   "

    jsonFilesList = []
    # dictionary are not ordered, so the loop is randomly distributed on the keys
    # for k, fname in Dataset2016.items():
    #
    # The following looping method guaranitee and ordered loop onto the keys
    for key in sorted(Dataset):
        print "key",key
        fname = Dataset[key]

        # Class define in lumi_utils that handles Runs and Lumi Blocks
        datasetRunLumis = getLumis(fname)

        path = "data/"
        # print allRunLumis
        outFileName = "DT_" + key + ".json"
        outFileCSV  = "DT_" + key + ".csv"

        outFileName = path + "RunJson" + outFileName
        outFileCSV  = path + "IntLumi" + outFileCSV

        jsonFilesList.append(outFileName)
        print "outFileName",outFileName
        datasetRunLumis.writeToFile(outFileName)
        print "fname",fname
        print "\t SAVED!! - {:s}".format(fname)
    ### Run brilcalc
        subprocess.call(["brilcalc", "lumi", "-u/pb", "-o"+outFileCSV, "-i"+outFileName])
