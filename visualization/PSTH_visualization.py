'''
Created on Jan 28, 2013

ongoing activity L2 neuron model

@author: robert
'''

import sys
import os, os.path
import glob
import numpy as np
import matplotlib.pyplot as plt

def ongoing_activity_subthreshold(folder, suffix, contains, outName):
    '''
    Vm histogram of ongoing activity
    '''
    fnames = []
    scan_directory2(folder, fnames, suffix, contains)
    
    VmMin = -80.0
    VmMax = 40.0
    binSize = 1.0
    bins_ = np.arange(VmMin, VmMax, binSize)
    totalHist = np.zeros_like(bins_[:-1])
    tOngoingBegin = 120.0
    tOngoingEnd = 220.0
    dt = 0.025
    tOngoingBeginBin = int(tOngoingBegin/dt + 0.5)
    tOngoingEndBin = int(tOngoingEnd/dt + 0.5)
    for fname in fnames:
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        for i in range(1, len(data)):
            hist, bins = np.histogram(data[i][tOngoingBeginBin:tOngoingEndBin], bins_)
            totalHist += hist
    norm = np.sum(totalHist)
    if norm:
        totalHist = 1.0/norm*totalHist
    
    if not outName.endswith('.csv'):
        outName += '.csv'
    with open(outName, 'w') as outFile:
        header = 'Vm bin\thist\n'
        for i in range(len(totalHist)):
            line = str(i + VmMin)
            line += '\t'
            line += str(totalHist[i])
            line += '\n'
            outFile.write(line)

def ongoing_activity_by_location(folder, suffix):
    '''
    3x3 plot of ongoing activity by postsynaptic location
    within C2 column
    '''
    fnames = []
    scan_directory(folder, fnames, suffix)
    
    print 'Creating spike raster plots from %d files' % len(fnames)
    
    whiskerLocations = {'B1border': 'B1', 'B2border': 'B2', 'B3border': 'B3',\
                        'C1border': 'C1', 'C2center': 'C2', 'C3border': 'C3',\
                        'D1border': 'D1', 'D2border': 'D2', 'D3border': 'D3'}
    BRowWhiskers = ('B3', 'B2', 'B1')                    
    CRowWhiskers = ('C3', 'C2', 'C1')
    DRowWhiskers = ('D3', 'D2', 'D1')
    tOffset = 100.0
    tStop = 350.0
    tStim = 245.0
    tStimWindow = 50.0
    tOngoingBegin = 120.0
    tOngoingWindow = 100.0
    binSize = 5.0
    ongoingStartBin = int((tOngoingBegin - tOffset)/binSize)
    ongoingEndBin = int((tOngoingBegin + tOngoingWindow - tOffset)/binSize)
    
    ongoingRates = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': []}
    
    for n in range(len(fnames)):
        fname = fnames[n]
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        PSTH = data[2]
        cumulativePSTH = np.sum(PSTH[ongoingStartBin:ongoingEndBin])
        ongoingRate = cumulativePSTH*10.0
        foundLocation = False
        for location in whiskerLocations:
            if location in fname:
                ongoingRates[whiskerLocations[location]].append(ongoingRate)
                foundLocation = True
        if not foundLocation:
            errstr = 'Could not find location string in filename %s' % (fname)
            raise RuntimeError(errstr)
    
    BRowRates = []
    for whisker in BRowWhiskers:
        BRowRates.append(np.mean(ongoingRates[whisker]))
        print '%s: %.2f Hz' % (whisker, np.mean(ongoingRates[whisker]))
    CRowRates = []
    for whisker in CRowWhiskers:
        CRowRates.append(np.mean(ongoingRates[whisker]))
        print '%s: %.2f Hz' % (whisker, np.mean(ongoingRates[whisker]))
    DRowRates = []
    for whisker in DRowWhiskers:
        DRowRates.append(np.mean(ongoingRates[whisker]))
        print '%s: %.2f Hz' % (whisker, np.mean(ongoingRates[whisker]))
    
    locationRates = []
    locationRates.append(DRowRates)
    locationRates.append(CRowRates)
    locationRates.append(BRowRates)
    ratesImage = np.array(locationRates)
    plt.figure(1)
    plt.imshow(ratesImage, interpolation='nearest', cmap='Blues')
    plt.colorbar()
    plt.show()

def evoked_activity_all(folder, suffix, outName):
    """docstring for evoked_activity_all"""
    fnames = []
    #scan_directory2(folder, fnames, suffix, 'UpState_INH_PW_1.0_SuW_0.5')
    scan_directory2(folder, fnames, suffix, 'UpState_dynamic_syns_INH_PW_1.0_SuW_0.5')
#    scan_directory2(folder, fnames, suffix, 'UpState_L6cc_timing_like_active_INH_PW_1.0_SuW_0.5')
    #scan_directory2(folder, fnames, suffix, 'UpState_L6ccinact_INH_PW_1.0_SuW_0.5')
#    scan_directory2(folder, fnames, suffix, 'UpState_INH_PW_1.0_SuW_0.5_dynamic_syns')
    #scan_directory2(folder, fnames, suffix, 'UpState_L5tt_evoked_inact_INH_PW_1.0_SuW_0.5')
    #scan_directory2(folder, fnames, suffix, 'UpState_VPM_evoked_inact_INH_PW_1.0_SuW_0.5')
    #scan_directory2(folder, fnames, suffix, 'UpState_L3py-L4sp_evoked_inact_INH_PW_1.0_SuW_0.5')
    #scan_directory2(folder, fnames, suffix, 'UpState_L4ss_evoked_inact_INH_PW_1.0_SuW_0.5')
    
    print 'Creating plots from %d files' % len(fnames)
    
    whiskerLocations = {'B1border': 'B1', 'B2border': 'B2', 'B3border': 'B3',\
                        'C1border': 'C1', 'C2center': 'C2', 'C3border': 'C3',\
                        'D1border': 'D1', 'D2border': 'D2', 'D3border': 'D3'}
    whiskerDeflections = {'B1_evoked': 'B1', 'B2_evoked': 'B2', 'B3_evoked': 'B3',\
                        'C1_evoked': 'C1', 'C2_evoked': 'C2', 'C3_evoked': 'C3',\
                        'D1_evoked': 'D1', 'D2_evoked': 'D2', 'D3_evoked': 'D3', 'E2_evoked': 'E2'}
    BRowWhiskers = ('B3', 'B2', 'B1')                    
    CRowWhiskers = ('C3', 'C2', 'C1')
    DRowWhiskers = ('D3', 'D2', 'D1')
    ERowWhiskers = ('E2',)
    tOffset = 100.0
    tStop = 350.0
    tStim = 245.0
    tStimWindow = 50.0
    tStimWindow25 = 25.0
    tOngoingBegin = 120.0
    tOngoingWindow = 100.0
    binSize = 5.0
    ongoingStartBin = int((tOngoingBegin - tOffset)/binSize)
    ongoingEndBin = int((tOngoingBegin + tOngoingWindow - tOffset)/binSize)
    evokedStartBin = int((tStim - tOffset)/binSize)
    evokedEndBin = int((tStim + tStimWindow - tOffset)/binSize)
    evokedEndBin25 = int((tStim + tStimWindow25 - tOffset)/binSize)
    
    evokedSpikesPerWhisker = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    evokedSpikesPerWhisker_0_10 = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    evokedSpikesPerWhisker_10_20 = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    evokedSpikesPerWhisker_20_30 = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    evokedSpikesPerWhisker_30_40 = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    evokedSpikesPerWhisker_40_50 = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    
    ongoingMeanTotal_ = []
    for n in range(len(fnames)):
        fname = fnames[n]
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        PSTH = data[2]
        ongoingMean = np.mean(PSTH[ongoingStartBin:ongoingEndBin])
        ongoingMeanTotal_.append(ongoingMean)
#        evokedPSTH = PSTH - ongoingMean
        evokedPSTH = PSTH
        #evokedSpikes = np.sum(evokedPSTH[evokedStartBin:evokedEndBin])
        evokedSpikes = np.sum(evokedPSTH[evokedStartBin:evokedEndBin25])
        evokedSpikes_0_10 = np.sum(evokedPSTH[evokedStartBin:evokedStartBin+2])
        evokedSpikes_10_20 = np.sum(evokedPSTH[evokedStartBin+2:evokedStartBin+4])
        evokedSpikes_20_30 = np.sum(evokedPSTH[evokedStartBin+4:evokedStartBin+6])
        evokedSpikes_30_40 = np.sum(evokedPSTH[evokedStartBin+6:evokedStartBin+8])
        evokedSpikes_40_50 = np.sum(evokedPSTH[evokedStartBin+8:evokedStartBin+10])
        foundWhisker = False
        for deflection in whiskerDeflections:
            if deflection in fname:
                evokedSpikesPerWhisker[whiskerDeflections[deflection]].append(evokedSpikes)
                evokedSpikesPerWhisker_0_10[whiskerDeflections[deflection]].append(evokedSpikes_0_10)
                evokedSpikesPerWhisker_10_20[whiskerDeflections[deflection]].append(evokedSpikes_10_20)
                evokedSpikesPerWhisker_20_30[whiskerDeflections[deflection]].append(evokedSpikes_20_30)
                evokedSpikesPerWhisker_30_40[whiskerDeflections[deflection]].append(evokedSpikes_30_40)
                evokedSpikesPerWhisker_40_50[whiskerDeflections[deflection]].append(evokedSpikes_40_50)
                foundWhisker = True
        if not foundWhisker:
            errstr = 'Could not find whisker string in filename %s' % (fname)
            raise RuntimeError(errstr)
    
    ongoingMeanTotal = np.mean(ongoingMeanTotal_)*10.0 #per 50ms
    print 'ongoingMeanTotal: ', ongoingMeanTotal
    
    print 'PSTH 0-10ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'E\t',
    for whisker in ERowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker])),
    print ''
    print 'D\t',
    for whisker in DRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker])),
    print ''
    print 'C\t',
    for whisker in CRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker])),
    print ''
    print 'B\t',
    for whisker in BRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker])),
    print ''
    
    print 'PSTH 10-20ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'E\t',
    for whisker in ERowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker])),
    print ''
    print 'D\t',
    for whisker in DRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker])),
    print ''
    print 'C\t',
    for whisker in CRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker])),
    print ''
    print 'B\t',
    for whisker in BRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker])),
    print ''
    print 'PSTH 20-30ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'E\t',
    for whisker in ERowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker])),
    print ''
    print 'D\t',
    for whisker in DRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker])),
    print ''
    print 'C\t',
    for whisker in CRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker])),
    print ''
    print 'B\t',
    for whisker in BRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker])),
    print ''
    
    print 'PSTH 30-40ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'E\t',
    for whisker in ERowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker])),
    print ''
    print 'D\t',
    for whisker in DRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker])),
    print ''
    print 'C\t',
    for whisker in CRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker])),
    print ''
    print 'B\t',
    for whisker in BRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker])),
    print ''
    
    print 'PSTH 40-50ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'E\t',
    for whisker in ERowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker])),
    print ''
    print 'D\t',
    for whisker in DRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker])),
    print ''
    print 'C\t',
    for whisker in CRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker])),
    print ''
    print 'B\t',
    for whisker in BRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker])),
    print ''
    
    print 'PSTH 0-25ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'E\t',
    ERowSpikes = []
    for whisker in ERowWhiskers:
#        ERowSpikes.append(max(np.mean(evokedSpikesPerWhisker[whisker]), ongoingMeanTotal))
        ERowSpikes.append(np.mean(evokedSpikesPerWhisker[whisker]))
#        DRowSpikes.append(max(np.mean(evokedSpikesPerWhisker[whisker]), 0.0))
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker])),
    print ''
    print 'D\t',
    DRowSpikes = []
    for whisker in DRowWhiskers:
#        DRowSpikes.append(max(np.mean(evokedSpikesPerWhisker[whisker]), ongoingMeanTotal))
        DRowSpikes.append(np.mean(evokedSpikesPerWhisker[whisker]))
#        DRowSpikes.append(max(np.mean(evokedSpikesPerWhisker[whisker]), 0.0))
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker])),
    print ''
    print 'C\t',
    CRowSpikes = []
    for whisker in CRowWhiskers:
#        CRowSpikes.append(max(np.mean(evokedSpikesPerWhisker[whisker]), ongoingMeanTotal))
        CRowSpikes.append(np.mean(evokedSpikesPerWhisker[whisker]))
#        CRowSpikes.append(max(np.mean(evokedSpikesPerWhisker[whisker]), 0.0))
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker])),
    print ''
    print 'B\t',
    BRowSpikes = []
    for whisker in BRowWhiskers:
#        BRowSpikes.append(max(np.mean(evokedSpikesPerWhisker[whisker]), ongoingMeanTotal))
        BRowSpikes.append(np.mean(evokedSpikesPerWhisker[whisker]))
#        BRowSpikes.append(max(np.mean(evokedSpikesPerWhisker[whisker]), 0.0))
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker])),
    print ''
    
    with open(outName+'.csv', 'w') as outFile:
        line =  'PSTH 0-10ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'E\t'
        for whisker in ERowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker]))
        line += '\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker]))
        line += '\n'
        
        line += 'PSTH 10-20ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'E\t'
        for whisker in ERowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker]))
        line += '\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker]))
        line += '\n'
        line += 'PSTH 20-30ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'E\t'
        for whisker in ERowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker]))
        line += '\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker]))
        line += '\n'
        
        line += 'PSTH 30-40ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'E\t'
        for whisker in ERowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker]))
        line += '\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker]))
        line += '\n'
        
        line += 'PSTH 40-50ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'E\t'
        for whisker in ERowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker]))
        line += '\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker]))
        line += '\n'
        
        line += 'PSTH 0-25ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'E\t'
        for whisker in ERowWhiskers:
#            line += '%.3f\t' % (max(np.mean(evokedSpikesPerWhisker[whisker]), ongoingMeanTotal))
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker]))
        line += '\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
#            line += '%.3f\t' % (max(np.mean(evokedSpikesPerWhisker[whisker]), ongoingMeanTotal))
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
#            line += '%.3f\t' % (max(np.mean(evokedSpikesPerWhisker[whisker]), ongoingMeanTotal))
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
#            line += '%.3f\t' % (max(np.mean(evokedSpikesPerWhisker[whisker]), ongoingMeanTotal))
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker]))
        line += '\n'
        outFile.write(line)
    
    whiskerRF = []
    whiskerRF.append(DRowSpikes)
    whiskerRF.append(CRowSpikes)
    whiskerRF.append(BRowSpikes)
    RFImage = np.array(whiskerRF)
    plt.figure(1)
    plt.imshow(RFImage, interpolation='nearest', cmap='Blues')
    plt.colorbar()
    plt.savefig(outName+'.pdf')
    plt.show()

def evoked_activity_by_location(folder, suffix, outName, locationStr):
    """docstring for evoked_activity_by_location"""
    fnames = []
#    scan_directory2(folder, fnames, suffix, 'UpState_INH_PW_1.0_SuW_0.5')
#    scan_directory2(folder, fnames, suffix, 'UpState_L6ccinact_INH_PW_1.0_SuW_0.5')
#    scan_directory2(folder, fnames, suffix, 'UpState_INH_PW_1.0_SuW_0.5_dynamic_syns')
    scan_directory2(folder, fnames, suffix, 'UpState_INH_PW_1.0_SuW_0.5_L6cc_timing_like_active')
    
    print 'Creating plots from %d files' % len(fnames)
    
    whiskerLocations = ('B1border', 'B2border', 'B3border', 'C1border', 'C2center',\
                        'C3border', 'D1border', 'D2border', 'D3border')
    if locationStr not in whiskerLocations:
        errstr = 'Unknown location: %s' % locationStr
        raise RuntimeError(errstr)
    
    whiskerDeflections = {'B1_evoked': 'B1', 'B2_evoked': 'B2', 'B3_evoked': 'B3',\
                        'C1_evoked': 'C1', 'C2_evoked': 'C2', 'C3_evoked': 'C3',\
                        'D1_evoked': 'D1', 'D2_evoked': 'D2', 'D3_evoked': 'D3', 'E2_evoked': 'E2'}
    BRowWhiskers = ('B3', 'B2', 'B1')                    
    CRowWhiskers = ('C3', 'C2', 'C1')
    DRowWhiskers = ('D3', 'D2', 'D1')
    tOffset = 100.0
    tStop = 350.0
    tStim = 245.0
    tStimWindow = 50.0
    tOngoingBegin = 120.0
    tOngoingWindow = 100.0
    binSize = 5.0
    ongoingStartBin = int((tOngoingBegin - tOffset)/binSize)
    ongoingEndBin = int((tOngoingBegin + tOngoingWindow - tOffset)/binSize)
    evokedStartBin = int((tStim - tOffset)/binSize)
    evokedEndBin = int((tStim + tStimWindow - tOffset)/binSize)
    
    evokedSpikesPerWhisker = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    evokedSpikesPerWhisker_0_10 = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    evokedSpikesPerWhisker_10_20 = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    evokedSpikesPerWhisker_20_30 = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    evokedSpikesPerWhisker_30_40 = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    evokedSpikesPerWhisker_40_50 = {'B1': [], 'B2': [], 'B3': [], 'C1': [], 'C2': [], 'C3': [], 'D1': [], 'D2': [], 'D3': [], 'E2': []}
    
    for n in range(len(fnames)):
        fname = fnames[n]
        if fname.find(locationStr) == -1:
            continue
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        PSTH = data[2]
        ongoingMean = np.mean(PSTH[ongoingStartBin:ongoingEndBin])
        evokedPSTH = PSTH - ongoingMean
        evokedSpikes = np.sum(evokedPSTH[evokedStartBin:evokedEndBin])
        evokedSpikes_0_10 = np.sum(evokedPSTH[evokedStartBin:evokedStartBin+2])
        evokedSpikes_10_20 = np.sum(evokedPSTH[evokedStartBin+2:evokedStartBin+4])
        evokedSpikes_20_30 = np.sum(evokedPSTH[evokedStartBin+4:evokedStartBin+6])
        evokedSpikes_30_40 = np.sum(evokedPSTH[evokedStartBin+6:evokedStartBin+8])
        evokedSpikes_40_50 = np.sum(evokedPSTH[evokedStartBin+8:evokedStartBin+10])
        foundWhisker = False
        for deflection in whiskerDeflections:
            if deflection in fname:
                evokedSpikesPerWhisker[whiskerDeflections[deflection]].append(evokedSpikes)
                evokedSpikesPerWhisker_0_10[whiskerDeflections[deflection]].append(evokedSpikes_0_10)
                evokedSpikesPerWhisker_10_20[whiskerDeflections[deflection]].append(evokedSpikes_10_20)
                evokedSpikesPerWhisker_20_30[whiskerDeflections[deflection]].append(evokedSpikes_20_30)
                evokedSpikesPerWhisker_30_40[whiskerDeflections[deflection]].append(evokedSpikes_30_40)
                evokedSpikesPerWhisker_40_50[whiskerDeflections[deflection]].append(evokedSpikes_40_50)
                foundWhisker = True
        if not foundWhisker:
            errstr = 'Could not find whisker string in filename %s' % (fname)
            raise RuntimeError(errstr)
    
    print '***************************'
    print 'Location: %s' % locationStr
    print '***************************'
    print 'PSTH 0-10ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'D\t',
    for whisker in DRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker])),
    print ''
    print 'C\t',
    for whisker in CRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker])),
    print ''
    print 'B\t',
    for whisker in BRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker])),
    print ''
    
    print 'PSTH 10-20ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'D\t',
    for whisker in DRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker])),
    print ''
    print 'C\t',
    for whisker in CRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker])),
    print ''
    print 'B\t',
    for whisker in BRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker])),
    print ''
    print 'PSTH 20-30ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'D\t',
    for whisker in DRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker])),
    print ''
    print 'C\t',
    for whisker in CRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker])),
    print ''
    print 'B\t',
    for whisker in BRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker])),
    print ''
    
    print 'PSTH 30-40ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'D\t',
    for whisker in DRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker])),
    print ''
    print 'C\t',
    for whisker in CRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker])),
    print ''
    print 'B\t',
    for whisker in BRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker])),
    print ''
    
    print 'PSTH 40-50ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'D\t',
    for whisker in DRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker])),
    print ''
    print 'C\t',
    for whisker in CRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker])),
    print ''
    print 'B\t',
    for whisker in BRowWhiskers:
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker])),
    print ''
    
    print 'PSTH 0-50ms'
    print 'Row\\Arc\t3\t2\t1'
    print 'D\t',
    DRowSpikes = []
    for whisker in DRowWhiskers:
        DRowSpikes.append(np.mean(evokedSpikesPerWhisker[whisker]))
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker])),
    print ''
    print 'C\t',
    CRowSpikes = []
    for whisker in CRowWhiskers:
        CRowSpikes.append(np.mean(evokedSpikesPerWhisker[whisker]))
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker])),
    print ''
    print 'B\t',
    BRowSpikes = []
    for whisker in BRowWhiskers:
        BRowSpikes.append(np.mean(evokedSpikesPerWhisker[whisker]))
        print '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker])),
    print ''
    
    with open(outName+locationStr+'.csv', 'w') as outFile:
        line =  'PSTH 0-10ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_0_10[whisker]))
        line += '\n'
        
        line += 'PSTH 10-20ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_10_20[whisker]))
        line += '\n'
        line += 'PSTH 20-30ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_20_30[whisker]))
        line += '\n'
        
        line += 'PSTH 30-40ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_30_40[whisker]))
        line += '\n'
        
        line += 'PSTH 40-50ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker_40_50[whisker]))
        line += '\n'
        
        line += 'PSTH 0-50ms'
        line += 'Row\\Arc\t3\t2\t1\n'
        line += 'D\t'
        for whisker in DRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker]))
        line += '\n'
        line += 'C\t'
        for whisker in CRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker]))
        line += '\n'
        line += 'B\t'
        for whisker in BRowWhiskers:
            line += '%.3f\t' % (np.mean(evokedSpikesPerWhisker[whisker]))
        line += '\n'
        outFile.write(line)
    
    whiskerRF = []
    whiskerRF.append(DRowSpikes)
    whiskerRF.append(CRowSpikes)
    whiskerRF.append(BRowSpikes)
    RFImage = np.array(whiskerRF)
    plt.figure(1)
    plt.imshow(RFImage, interpolation='nearest', cmap='Blues')
    plt.colorbar()
    plt.savefig(outName+locationStr+'.pdf')
    plt.show()

def scan_directory(path, fnames, suffix):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory(fname, fnames, suffix)
        elif fname.endswith(suffix):
            fnames.append(fname)
        else:
            continue

def scan_directory2(path, fnames, suffix, contains):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory2(fname, fnames, suffix, contains)
        elif fname.endswith(suffix) and fname.find(contains) != -1:
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 4:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        outName = sys.argv[3]
        evoked_activity_all(folder, suffix, outName)
    elif len(sys.argv) == 5:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        contains = sys.argv[3]
        outName = sys.argv[4]
#        outName = sys.argv[3]
#        location = sys.argv[4]
#        evoked_activity_by_location(folder, suffix, outName, location)
        ongoing_activity_subthreshold(folder, suffix, contains, outName)
    else:
        print 'Error! Number of arguments is %d; should be 3 or 4' % (len(sys.argv)-1)
    
    
    
    