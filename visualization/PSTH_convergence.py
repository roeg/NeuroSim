'''
Created on Jan 28, 2013

ongoing activity L2 neuron model

@author: robert
'''

import sys
import os, os.path
import glob
import single_cell_analyzer as sca
import single_cell_parser as scp
import numpy as np
import matplotlib.pyplot as plt

def estimate_PSTH_convergence_rate(folder, suffix, contains):
    '''
    load all traces, compute spike times and PSTH
    and determine empirical rate of convergence
    '''
    fnames = []
    scan_directory2(folder, fnames, suffix, contains)
    
    print 'Creating spike raster plots from %d files' % len(fnames)
    
    tOffset = 100.0
    tStop = 345.0
    tOngoing = 120.0
    tOngoingWindow = 100.0
    tStim = 245.0
    tStimWindow = 50.0
    binSize = 5.0
    ongoingBeginBin = int((tOngoing-tOffset)/binSize + 0.5)
    ongoingEndBin = int((tOngoing+tOngoingWindow-tOffset)/binSize + 0.5)
    stimBeginBin = int((tStim-tOffset)/binSize + 0.5)
    stimEndBin = int((tStim+tStimWindow-tOffset)/binSize + 0.5)
    
    allOngoingRMSE = []
    allPSTH_RMSE = []
    
    for n in range(len(fnames)):
        fname = fnames[n]
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        t = data[0]
        allSpikeTimes = []
        cumulativePSTHs = []
        for i in range(1, len(data)):
            v = data[i]
            spikeTimes = sca.simple_spike_detection(t, v)
            allSpikeTimes.append(spikeTimes)
            # PSTH of all traces so far
            cumulativePSTHs.append(sca.PSTH_from_spike_times(allSpikeTimes, binSize, tOffset, tStop))
        hist, bins = sca.PSTH_from_spike_times(allSpikeTimes, binSize, tOffset, tStop)
        refOngoing = np.mean(hist[ongoingBeginBin:ongoingEndBin])
        refPSTH = np.array(hist[stimBeginBin:stimEndBin]) - refOngoing
        
        ongoingRMSE = []
        PSTH_RMSE = []
        for PSTH in cumulativePSTHs:
            hist, bins = PSTH
            currentOngoing = np.mean(hist[ongoingBeginBin:ongoingEndBin])
            currentPSTH = np.array(hist[stimBeginBin:stimEndBin]) - currentOngoing
            ongoingRMSE.append(np.abs(currentOngoing - refOngoing)/refOngoing)
            PSTH_RMSE.append(np.sqrt(np.dot(currentPSTH - refPSTH, currentPSTH - refPSTH))/np.max(refPSTH))
        
        allOngoingRMSE.append(ongoingRMSE)
        allPSTH_RMSE.append(PSTH_RMSE)
    
    avgOngoingRMSE = np.mean(np.array(allOngoingRMSE), axis=0)
    stdOngoingRMSE = np.std(np.array(allOngoingRMSE), axis=0)
    avgPSTH_RMSE = np.mean(np.array(allPSTH_RMSE), axis=0)
    stdPSTH_RMSE = np.mean(np.array(allPSTH_RMSE), axis=0)
    trials = [i+1 for i in range(len(avgOngoingRMSE))]
    fig = plt.figure(1)
    fig.add_subplot(2,1,1)
    plt.plot(trials, avgOngoingRMSE, 'k')
    plt.plot(trials, avgOngoingRMSE+stdOngoingRMSE, 'k--')
    plt.plot(trials, avgOngoingRMSE-stdOngoingRMSE, 'k--')
    plt.ylabel('ongoing RMSE')
    plt.ylim([0, np.max(avgOngoingRMSE+stdOngoingRMSE)])
    fig.add_subplot(2,1,2)
    plt.plot(trials, avgPSTH_RMSE, 'k')
    plt.plot(trials, avgPSTH_RMSE+stdPSTH_RMSE, 'k--')
    plt.plot(trials, avgPSTH_RMSE-stdPSTH_RMSE, 'k--')
    plt.ylim([0, np.max(avgPSTH_RMSE+stdPSTH_RMSE)])
    plt.xlabel('nr. of trials')
    plt.ylabel('PSTH RMSE')
    plt.show()
#    outName = fname[:-4]
#    outName += '_PSTH_total_10ms.pdf'
#    outName += '_PSTH_total_5ms.pdf'
#    plt.savefig(outName)
#    scp.write_PSTH(outName[:-4]+'.csv', hist, bins)

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
        contains = sys.argv[3]
        estimate_PSTH_convergence_rate(folder, suffix, contains)
    else:
        print 'Error! Number of arguments is %d; should be 3' % (len(sys.argv)-1)
    
    
    
    