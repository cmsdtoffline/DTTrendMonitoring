# *-* coding=utf-8 *-*
import sys
sys.path.append('/usr/local/lib/python2.7/site-packages')
import csv
import xlrd
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
%matplotlib inline
import matplotlib.pyplot as plt

def IntegratedLumi(lumiSeries):
    integratedLumi = []
    partialSum = 0
    for item in lumiSeries:
        partialSum+=item
        integratedLumi.append(partialSum)
        # print ' Integrated Luminosity (/ub): {:.2f} \t Delivered: {:.2f}'.format(partialSum,item)
    return integratedLumi

#####
# Function that compute Eff and EffErr, rebinning the recordedLumi or the Num value
#######
def ntupleForEffPlots(dataFrame):
    eff  = []
    effErr  = []
    lumi = []
    tempLumi = 0.
    tempNum  = 0.
    tempDen  = 0.
    tempErr  = 0.
    for row in dataFrame.itertuples(index=False):
        if row.Den == 0:
            continue
        tempLumi += row.recordedLumi
        tempNum  += row.Num
        tempDen  += row.Den
        tempErr  = np.sqrt(tempNum*( 1 - tempNum/tempDen ))/tempDen
        #tempRow[col_list].add(row[col_list])
        if tempLumi > 40 and tempNum > 1000: #or not next(st0wh1.itertuples(index=False)):
            eff.append(tempNum/tempDen)
            effErr.append(tempErr)
            lumi.append(tempLumi)
            tempLumi = 0.
            tempErr  = 0.
            tempNum  = 0.
            tempDen  = 0.
    return eff, effErr, lumi

dataset = 'DT_Run2016C.csv'

# .csv file obtained with EfficiencyMonitorRun.C
# It contains on each lines: Run - Station - Wheel - Den - Num - NumA
results = '../ResultsRun_Run2016C.csv'

# Pandas DataFrame reading .csv obationed with brilcalCSV_Converter.py
# It contains the links btw runNumber and integrated luminosity
DF_Run = pd.read_csv('LumiTable_'+dataset)

#print DF_Run
#DF_Run['recorded(/pb)']


with open(results, 'rb') as csvfile:
    # Pandas DataFrame reading .csv obationed with EfficiencyMonitorRun.C
    DF_Results = pd.read_csv(csvfile)
    # Adding "recordedLumi" column
    DF_Results['recordedLumi'] = 0.

    #print DF_Results.head(5)

    # Looping on runNumber containend in the DF_Run DataFrame
    for run in DF_Run.RunNumber:
        # Getting the integrated luminosity associated to the runNumber
        recordedLumi = DF_Run.loc[(DF_Run.RunNumber == run), 'recorded(/pb)']
        if np.isnan(recordedLumi).bool():
            continue
        # print recordedLumi
        # typecast for int and float added otherwise it does not work
        # The .ix method allows you to add the recordedLumi value to all the
        # DataFrame entries that match the condition: DF_Results.Run == int(run)
        DF_Results.ix[DF_Results.Run == int(run), 'recordedLumi'] = float(recordedLumi)

        #print DF_Results[(DF_Results.Run == int(run))].recordedLumi


#-----# Trying Database
#
# import MySQLdb
# from sqlalchemy import create_engine
#
# #engine = create_engine('mysql+mysqlconnector://[user]:[pass]@[host]:[port]/[schema]', echo=False)
# #data.to_sql(name='sample_table2', con=engine, if_exists = 'append', index=False)
# db = MySQLdb.connect(host="localhost",    # your host, usually localhost
#                      user="root",         # your username
#                      passwd="MySQLpassword",  # your password
#                      db="CMS_DT")
# DF_Results.to_sql(con=db, name='CMS_DT')


######
# Plotting the efficiency Vs integrated luminosity for the 5 DT station (phi layers)
######

# Loop on the DT stations
for station in range(DF_Results.Station.min(),DF_Results.Station.max()+1,1):
    # matplotlib figure
    fig2 = plt.figure(station)
    ax2 = fig2.add_subplot(1,1,1)

    for wheel in range(DF_Results.Wheel.min(), DF_Results.Wheel.max()+1,1):
        #print "Station: {:d} - Wheel: {:d} ".format(station,wheel)

        # Subset of the main DataFrame for a particular DT wheel and station
        df_temp = DF_Results[(DF_Results.Wheel == wheel) & (DF_Results.Station == station)]
        # Function that compute Eff and EffErr, rebinning the recordedLumi or the Num value
        eff, effErr, lumi = ntupleForEffPlots(df_temp)

        # ntuple for the x axis
        xi = np.array(IntegratedLumi(lumi))
        # ntuple for y axis plus errors on the efficiencies
        y  = np.array(eff)
        y_err = np.array(effErr)

        ax2.errorbar(xi,y,yerr=y_err, fmt='.', label='Wheel {:d}'.format(wheel))

        #     fig.savefig('runC_st{:d}_wh{:d}.pdf'.format(station, wheel), dpi=4)
    #
    #     fig = plt.figure()
    #     ax = fig.add_subplot(1,1,1)
    #     ax.errorbar(xi,y,yerr=y_err, fmt='.')
    #     ax.set_xlabel('Integrated Lumi (/pb)')
    #     ax.set_xticks(np.arange(0, 2500, 250))
    #     ax.set_yticks(np.arange(0.94, 1.0, 0.01))
    #     ax.set_ylabel('Eff')
    #     ax.grid(color='gray', linestyle='--', linewidth=0.2)
    #     #fig.set_grid(color='gray', linestyle='--', linewidth=0.2)
    #     ax.set_title('Run C 2016 at 13 TeV',fontsize=12, fontweight='bold')
    #     ax.set_title('Station = {:d}, Wheel = {:d}'.format(station, wheel),fontsize=9)
    #     fig.savefig('runC_st{:d}_wh{:d}.pdf'.format(station, wheel), dpi=400)

    ax2.set_xlabel(r'Integrated Lumi (pb$^{-1}$)')
    ax2.set_xticks(np.arange(0, 2500, 250))
    ax2.set_yticks(np.arange(0.94, 1.0, 0.01))
    ax2.set_ylabel(r'DT Efficincy of $\phi$ layers')
    ax2.grid(color='gray', linestyle='--', linewidth=0.2)
    ax2.set_title('Run C 2016 - MB{:d} hit efficiency'.format(station), fontsize=12, fontweight='bold')
    ax2.legend(loc='best', fontsize='x-small')#'Wheel = {:d}'.format(wheel),fontsize=9)
    fig2.savefig('runC_st{:d}.pdf'.format(station), dpi=400)
