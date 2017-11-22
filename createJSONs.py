# *-* coding=utf-8 *-*

# from ROOT import *
# import datetime
import sys
import subprocess
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


def getLumis(inputFile):
    RunLumis = lu.RunLumis({})

    vals = []
    vals.append(lu.getRunLumis(fname))
    allRunLumis = sum(vals, RunLumis)
    return allRunLumis


if __name__ == '__main__':

    # print "#################################"
    print "    Producing Lumi JSON starting from DT DQM TTree   "

    jsonFilesList = []
    # dictionary are not ordered, so the loop is randomly distributed on the keys
    # for k, fname in Dataset2016.items():
    #
    # The following looping method guaranitee and ordered loop onto the keys
    for key in sorted(Dataset2016):
        print key
        fname = Dataset2016[key]
        # Class define in lumi_utils that handles Runs and Lumi Blocks
        datasetRunLumis = getLumis(fname)

        # print allRunLumis
        outFileName = "DT_" + key + ".json"
        outFileCSV  = "DT_" + key + ".csv"
        jsonFilesList.append(outFileName)
        # print outFileName
        datasetRunLumis.writeToFile(outFileName)
        print "\t SAVED!! - {:s}".format(fname)
    ### Run brilcalc
        subprocess.call(["brilcalc", "lumi", "-u/pb", "-o"+outFileCSV, "-i"+outFileName])
