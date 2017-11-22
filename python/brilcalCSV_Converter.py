# *-* coding=utf-8 *-*

# Converting .csv obtained by BrilCalc in a plain format for Pandas

import sys
sys.path.append('/usr/local/lib/python2.7/site-packages')

import csv
import pandas as pd
import numpy as np

dataset = 'DT_Run2016C.csv'
csvPath = '../data/IntLumi/'

with open(csvPath+dataset, 'rb') as csvfile:
    outFile = open('LumiTable_'+dataset,'w')

    reader = csv.reader(csvfile)
    writer = csv.writer(outFile)
    row1 = next(reader)
    header = next(reader)
    header[0] = 'RunNumber'
    #print header
    #print ( (header[0].split(':'))[0] ).translate(None,\#)

    writer.writerow(header)

    for row in reader:
        firstChar = row[0][0]
        firstEntry = row[0]
        if firstChar == '#':
            continue
        print row
        # print firstEntry.split(':')
        row[0] = (firstEntry.split(':'))[0]
        # print row
        # print "="*10
        writer.writerow(row)
        # newrow = ' '.join(row)
        # splited = newrow.split(':')
        # splited2 = splited[1].split(' ')
        # print splited2
