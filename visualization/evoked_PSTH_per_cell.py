#!/usr/bin/python

import sys
import numpy as np
import os, os.path
import glob
from scipy.stats import poisson, binom

# anatomical PC + surround columns (3x3)
# ranging from (potentially) 1-9, starting at row-1, arc-1,
# then increasing by arc and then by row up to row+1, arc+1
# e.g. for C2: B1=1, B2=2, B3=3, C1=4, C2=5, C3=6, D1=7, D2=8, D3=9
surroundColumns = {'A1': {'Alpha': 4, 'A1': 5, 'A2': 6, 'B1': 8, 'B2': 9},\
                   'A2': {'A1': 4, 'A2': 5, 'A3': 6, 'B1': 7, 'B2': 8, 'B3': 9},\
                   'A3': {'A2': 4, 'A3': 5, 'A4': 6, 'B2': 7, 'B3': 8, 'B4': 9},\
                   'A4': {'A3': 4, 'A4': 5, 'B3': 7, 'B4': 8},\
                   'Alpha': {'Alpha': 5, 'A1': 6, 'Beta': 8, 'B1': 9},\
                   'B1': {'Alpha': 1, 'A1': 2, 'A2': 3, 'Beta': 4, 'B1': 5, 'B2': 6, 'C1': 8, 'C2': 9},\
                   'B2': {'A1': 1, 'A2': 2, 'A3': 3, 'B1': 4, 'B2': 5, 'B3': 6, 'C1': 7, 'C2': 8, 'C3': 9},\
                   'B3': {'A2': 1, 'A3': 2, 'A4': 3, 'B2': 4, 'B3': 5, 'B4': 6, 'C2': 7, 'C3': 8, 'C4': 9},\
                   'B4': {'A3': 1, 'A4': 2, 'B3': 4, 'B4': 5, 'C3': 7, 'C4': 8},\
                   'Beta': {'Alpha': 2, 'Beta': 5, 'B1': 6, 'Gamma': 8, 'C1': 9},\
                   'C1': {'Beta': 1, 'B1': 2, 'B2': 3, 'Gamma': 4, 'C1': 5, 'C2': 6, 'D1': 8, 'D2': 9},\
                   'C2': {'B1': 1, 'B2': 2, 'B3': 3, 'C1': 4, 'C2': 5, 'C3': 6, 'D1': 7, 'D2': 8, 'D3': 9},\
                   'C3': {'B2': 1, 'B3': 2, 'B4': 3, 'C2': 4, 'C3': 5, 'C4': 6, 'D2': 7, 'D3': 8, 'D4': 9},\
                   'C4': {'B3': 1, 'B4': 2, 'C3': 4, 'C4': 5, 'D3': 7, 'D4': 8},\
                   'Gamma': {'Beta': 2, 'Gamma': 5, 'C1': 6, 'Delta': 8, 'D1': 9},\
                   'D1': {'Gamma': 1, 'C1': 2, 'C2': 3, 'Delta': 4, 'D1': 5, 'D2': 6, 'E1': 8, 'E2': 9},\
                   'D2': {'C1': 1, 'C2': 2, 'C3': 3, 'D1': 4, 'D2': 5, 'D3': 6, 'E1': 7, 'E2': 8, 'E3': 9},\
                   'D3': {'C2': 1, 'C3': 2, 'C4': 3, 'D2': 4, 'D3': 5, 'D4': 6, 'E2': 7, 'E3': 8, 'E4': 9},\
                   'D4': {'C3': 1, 'C4': 2, 'D3': 4, 'D4': 5, 'E3': 7, 'E4': 8},\
                   'Delta': {'Gamma': 2, 'Delta': 5, 'D1': 6, 'E1': 9},\
                   'E1': {'Delta': 1, 'D1': 2, 'D2': 3, 'E1': 5, 'E2': 6},\
                   'E2': {'D1': 1, 'D2': 2, 'D3': 3, 'E1': 4, 'E2': 5, 'E3': 6},\
                   'E3': {'D2': 1, 'D3': 2, 'D4': 3, 'E2': 4, 'E3': 5, 'E4': 6},\
                   'E4': {'D3': 1, 'D4': 2, 'E3': 4, 'E4': 5}}
index2WhiskerLUT = {1: 'B1', 2: 'B2', 3: 'B3',\
            4: 'C1', 5: 'C2', 6: 'C3',\
            7: 'D1', 8: 'D2', 9: 'D3'}

def create_average_celltype_PSTH_from_clusters(cellTypeFolder, outFileName):
    '''
    loads cluster files from CDK recordings and automatically
    determines avg ongoing activity and evoked PSTH for cell type
    '''
    stimulusOnset = 145.0
    ongoingBegin = 20.0 # 100ms pre-stimulus
    ongoingDur = 100.0
#    ongoingBegin = 0.0 # 100ms pre-stimulus
    PSTHEnd = 195.0 # 0-50ms here
    # load all spike time files for all whiskers of all recorded cells
    print 'Calculating average evoked PSTH for all cells in folder %s' % cellTypeFolder
    fnames = []
    scan_directory(cellTypeFolder, fnames, '.cluster1')
    cellSpikeTimes = {}
    for fname in fnames:
        whiskerSpikeTimes = load_cluster_trials(fname)
        splitName = fname.split('/')
        trialName = splitName[-1]
        cellName = splitName[-2]
        if not cellSpikeTimes.has_key(cellName):
            cellSpikeTimes[cellName] = {}
        cellSpikeTimes[cellName][trialName] = whiskerSpikeTimes
    # calculate spontaneous activity (all spikes 100ms pre-stim
    # across all cells, whiskers, trials)
    rates = {}
    for cell in cellSpikeTimes:
        ongoingSpikes = 0.0
        ongoingTrials = 0.0
        for whisker in cellSpikeTimes[cell]:
            for trial in cellSpikeTimes[cell][whisker]:
                ongoingTrials += 1
                for t in cellSpikeTimes[cell][whisker][trial]:
                    if ongoingBegin <= t < ongoingBegin + ongoingDur:
                        ongoingSpikes += 1
        spontRate = ongoingSpikes/ongoingTrials*0.01 # per ms
        rates[cell] = spontRate
        print '\tcell name: %s' % cell
        print '\tSpontaneous firing rate = %.2f Hz' % (spontRate*1000.0)
    
    # collect all spike times and repetitions per cell
    # 0-50ms post-stimulus PW-aligned
    trialsPerWhisker = {}
    spikesPerWhisker = {}
    for cell in cellSpikeTimes:
        trialsPerWhisker[cell] = dict([(i,0) for i in range(1,10)])
        spikesPerWhisker[cell] = dict([(i,[]) for i in range(1,10)])
    for cell in cellSpikeTimes:
        PW = ''
        for whisker in cellSpikeTimes[cell]:
            if 'PW' in whisker:
                splitName = whisker.split('_')
                PW = splitName[0]
        for whisker in cellSpikeTimes[cell]:
            splitName = whisker.split('_')
            whiskerName = splitName[0]
            if whiskerName in surroundColumns[PW]:
                tmpSpikes = 0
                tmpTrials = 0
                col = surroundColumns[PW][whiskerName]
                for trial in cellSpikeTimes[cell][whisker]:
                    trialsPerWhisker[cell][col] += 1
                    tmpTrials += 1
                    for t in cellSpikeTimes[cell][whisker][trial]:
                        if stimulusOnset <= t < PSTHEnd:
                            spikesPerWhisker[cell][col].append(t)
                            tmpSpikes += 1
                if col == 5:
                    print cell
                    print 'PW: ', PW
                    print 'APs per stim: ', float(tmpSpikes)/tmpTrials
    # create 1ms resolution PSTH per whisker and subtract
    # spontaneous activity per 1ms bin
    numberOfBins = 50
    PSTHrange = (stimulusOnset,PSTHEnd)
    whiskerPSTH = {}
    for cell in cellSpikeTimes:
        whiskerPSTH[cell] = {}
        for col in spikesPerWhisker[cell]:
            nrOfTrials = trialsPerWhisker[cell][col]
            hist, bins = np.histogram(spikesPerWhisker[cell][col], numberOfBins, PSTHrange)
            if nrOfTrials:
                hist = 1.0/nrOfTrials*hist - rates[cell]
            else:
                hist = hist - rates[cell]
            whiskerPSTH[cell][col] = hist, bins
            print 'cell: ', cell
            print 'whisker: ', index2WhiskerLUT[col]
            print 'sum PSTH: ', np.sum(hist)
    
    cellNames = whiskerPSTH.keys()
    cellNames.sort()
    whiskers = whiskerPSTH[cellNames[0]].keys()
    whiskers.sort()
    if not outFileName.endswith('.csv'):
        outFileName += '.csv'
    with open(outFileName, 'w') as PSTHFile:
        header = 'Cell name'
        for whisker in whiskers:
            header += '\t'
            header += index2WhiskerLUT[whisker]
        header += '\n'
        PSTHFile.write(header)
        for cell in cellNames:
            line = cell
            for whisker in whiskers:
                line += '\t'
                PSTH = whiskerPSTH[cell][whisker]
                totalPSTH = 0.0
                if PSTH is not None:
                    hist, bins = PSTH
                    totalPSTH = np.sum(hist)
                line += str(totalPSTH)
            line += '\n'
            PSTHFile.write(line)

def create_single_cell_PSTH_from_clusters(cellTypeFolder):
    '''
    load cluster files from CDK recordings and automatically
    determines ongoing activity and evoked PSTH for individual cell
    '''
    stimulusOnset = 145.0
    ongoingBegin = 20.0 # 100ms pre-stimulus
    ongoingDur = 100.0
#    ongoingBegin = 0.0 # 100ms pre-stimulus
    PSTHEnd = 195.0 # 0-50ms here
    # load all spike time files for all whiskers of recorded cell
    print 'Calculating average evoked PSTH for all cells in folder %s' % cellTypeFolder
    fnames = []
    scan_directory(cellTypeFolder, fnames, '.cluster1')
    cellSpikeTimes = {}
    for fname in fnames:
        whiskerSpikeTimes = load_cluster_trials(fname)
        splitName = fname.split('/')
        trialName = splitName[-1]
        cellName = splitName[-2]
        whiskerName = trialName.split('_')[0]
        if not cellSpikeTimes.has_key(cellName):
            cellSpikeTimes[cellName] = {}
        cellSpikeTimes[cellName][whiskerName] = whiskerSpikeTimes
    # calculate spontaneous activity (all spikes 100ms pre-stim
    # across all cells, whiskers, trials)
    rates = {}
    for cell in cellSpikeTimes:
        ongoingSpikes = 0.0
        ongoingTrials = 0.0
        for whisker in cellSpikeTimes[cell]:
            for trial in cellSpikeTimes[cell][whisker]:
                ongoingTrials += 1
                for t in cellSpikeTimes[cell][whisker][trial]:
                    if ongoingBegin <= t < ongoingBegin + ongoingDur:
                        ongoingSpikes += 1
        spontRate = ongoingSpikes/ongoingTrials*0.01 # per ms
        rates[cell] = spontRate
        print '\tcell name: %s' % cell
        print '\tSpontaneous firing rate = %.2f Hz' % (spontRate*1000.0)
    
    # collect all spike times and repetitions per cell
    # 0-50ms post-stimulus PW-aligned
    trialsPerWhisker = {}
    spikesPerWhisker = {}
    for cell in cellSpikeTimes:
        trialsPerWhisker[cell] = dict([(whisker,0) for whisker in cellSpikeTimes[cell]])
        spikesPerWhisker[cell] = dict([(whisker,[]) for whisker in cellSpikeTimes[cell]])
    for cell in cellSpikeTimes:
        for whisker in cellSpikeTimes[cell]:
            tmpSpikes = 0
            tmpTrials = 0
            for trial in cellSpikeTimes[cell][whisker]:
                trialsPerWhisker[cell][whisker] += 1
                tmpTrials += 1
                for t in cellSpikeTimes[cell][whisker][trial]:
                    if stimulusOnset <= t < PSTHEnd:
                        spikesPerWhisker[cell][whisker].append(t)
                        tmpSpikes += 1
    
    # create 1ms resolution PSTH per whisker and subtract
    # spontaneous activity per 1ms bin
    numberOfBins = 50
    PSTHrange = (stimulusOnset,PSTHEnd)
    whiskerPSTH = {}
    for cell in cellSpikeTimes:
        whiskerPSTH[cell] = {}
        for whisker in spikesPerWhisker[cell]:
            nrOfTrials = trialsPerWhisker[cell][whisker]
            hist, bins = np.histogram(spikesPerWhisker[cell][whisker], numberOfBins, PSTHrange)
            if nrOfTrials:
                hist = 1.0/nrOfTrials*hist - rates[cell]
            else:
                hist = hist - rates[cell]
            whiskerPSTH[cell][whisker] = hist, bins
            print 'cell: ', cell
            print 'whisker: ', whisker
            print 'nr. trials: ', nrOfTrials
            print 'sum PSTH: ', np.sum(hist)
    
    cellNames = whiskerPSTH.keys()
    cellNames.sort()
    whiskers = whiskerPSTH[cellNames[0]].keys()
    whiskers.sort()
    for cell in cellNames:
        outFileName = os.path.join(cellTypeFolder, cell + '_evoked_50ms.csv')
        with open(outFileName, 'w') as PSTHFile:
            header = 'Cell name'
            header += '\t'
            header += 'Ongoing rate (Hz)'
            for whisker in whiskers:
                header += '\t'
#                header += index2WhiskerLUT[whisker]
                header += whisker
            header += '\n'
            PSTHFile.write(header)
            line = cell
            line += '\t'
            line += str(rates[cell]*1000.0)
            for whisker in whiskers:
                line += '\t'
                PSTH = whiskerPSTH[cell][whisker]
                totalPSTH = 0.0
                if PSTH is not None:
                    hist, bins = PSTH
                    totalPSTH = np.sum(hist)
                line += str(totalPSTH)
            line += '\n'
            PSTHFile.write(line)
    for cell in cellNames:
        outFileName = os.path.join(cellTypeFolder, cell + '_totalPSTH_1ms.csv')
        with open(outFileName, 'w') as PSTHFile:
            header = 'Cell name'
            header += '\t'
            header += 'time (ms)'
            for whisker in whiskers:
                header += '\t'
#                header += index2WhiskerLUT[whisker]
                header += whisker
            header += '\n'
            PSTHFile.write(header)
            bins = whiskerPSTH[cell][whiskers[0]][1]
            for i in range(len(bins)-1):
                tBin = 0.5*(bins[i] + bins[i+1])
                line = cell
                line += '\t'
                line += str(tBin)
                for whisker in whiskers:
                    line += '\t'
                    PSTH = whiskerPSTH[cell][whisker][0] + rates[cell]
                    line += str(PSTH[i])
                line += '\n'
                PSTHFile.write(line)

def create_single_cell_latency_from_clusters(cellTypeFolder):
    '''
    load cluster files from CDK recordings and automatically
    determines ongoing activity and evoked PSTH for individual cell
    '''
    stimulusOnset = 145.0
    ongoingBegin = 20.0 # 100ms pre-stimulus
    ongoingDur = 100.0
#    ongoingBegin = 0.0 # 100ms pre-stimulus
    #PSTHEnd = 195.0 # 0-50ms here
    PSTHEnd = 245.0 # 0-100ms as in CDK 07
    # load all spike time files for all whiskers of recorded cell
    print 'Calculating latencies (a la CDK 2007) for all cells in folder %s' % cellTypeFolder
    fnames = []
    scan_directory(cellTypeFolder, fnames, '.cluster1')
    cellSpikeTimes = {}
    for fname in fnames:
        whiskerSpikeTimes = load_cluster_trials(fname)
        splitName = fname.split('/')
        trialName = splitName[-1]
        cellName = splitName[-2]
        whiskerSplitNames = trialName.split('_')
        if len(whiskerSplitNames) == 2:
            whiskerName = whiskerSplitNames[0]
        else:
            whiskerName = whiskerSplitNames[0]
            for i in range(1, len(whiskerSplitNames) - 1):
                whiskerName += '_'
                whiskerName += whiskerSplitNames[i]
        if not cellSpikeTimes.has_key(cellName):
            cellSpikeTimes[cellName] = {}
        cellSpikeTimes[cellName][whiskerName] = whiskerSpikeTimes
    # calculate spontaneous activity (all spikes 100ms pre-stim
    # across all cells, whiskers, trials)
    rates = {}
    latencyThresholds = {}
    for cell in cellSpikeTimes:
        ongoingSpikes = 0.0
        ongoingTrials = 0.0
        for whisker in cellSpikeTimes[cell]:
            for trial in cellSpikeTimes[cell][whisker]:
                ongoingTrials += 1
                for t in cellSpikeTimes[cell][whisker][trial]:
                    if ongoingBegin <= t < ongoingBegin + ongoingDur:
                        ongoingSpikes += 1
        spontRate = ongoingSpikes/ongoingTrials*0.01 # per ms
        rates[cell] = spontRate
        # CDK 07: latency threshold is 1ms bin after stimulus
        # where PSTH is larger then 99th percentile of 
        # Poisson distribution based on spontaneous rate
        if spontRate > 0:
            percentile = 0.99
            threshold = poisson.ppf(percentile, spontRate*1000.0)
            #cumulativePercent = poisson.cdf(threshold, spontRate*1000.0)
            #print 'threshold = %.2f, cdf = %.2f ' % (threshold, cumulativePercent)
            threshold /= 1000.0
            latencyThresholds[cell] = threshold
        else:
            latencyThresholds[cell] = 0.0
        
        print '\tcell name: %s' % cell
        print '\tSpontaneous firing rate = %.2f Hz' % (spontRate*1000.0)
        print '\tLatency threshold = %.2f Hz' % (latencyThresholds[cell]*1000.0)
    
    # collect all spike times and repetitions per cell
    # 0-100ms post-stimulus PW-aligned
    trialsPerWhisker = {}
    spikesPerWhisker = {}
    for cell in cellSpikeTimes:
        trialsPerWhisker[cell] = dict([(i,0) for i in range(1,10)])
        spikesPerWhisker[cell] = dict([(i,[]) for i in range(1,10)])
    for cell in cellSpikeTimes:
        PW = ''
        for whisker in cellSpikeTimes[cell]:
            if 'PW' in whisker:
                splitName = whisker.split('_')
                PW = splitName[0]
        for whisker in cellSpikeTimes[cell]:
            splitName = whisker.split('_')
            whiskerName = splitName[0]
            if whiskerName in surroundColumns[PW]:
                tmpSpikes = 0
                tmpTrials = 0
                col = surroundColumns[PW][whiskerName]
                for trial in cellSpikeTimes[cell][whisker]:
                    trialsPerWhisker[cell][col] += 1
                    tmpTrials += 1
                    for t in cellSpikeTimes[cell][whisker][trial]:
                        if stimulusOnset <= t < PSTHEnd:
                            spikesPerWhisker[cell][col].append(t)
                            tmpSpikes += 1
                if col == 5:
                    print cell
                    print 'PW: ', PW
                    print 'APs per stim: ', float(tmpSpikes)/tmpTrials
    
    # create 1ms resolution PSTH per whisker
    #numberOfBins = 50
    numberOfBins = 100 # CDK 07
    PSTHrange = (stimulusOnset,PSTHEnd)
    PWLatencies = {}
    SWLatencies = {}
    for cell in cellSpikeTimes:
        SWSpikeTimes = []
        nrSWTrials = 0
        for whisker in spikesPerWhisker[cell]:
            if whisker == 5:
                continue
            nrSWTrials += trialsPerWhisker[cell][whisker]
            for t in spikesPerWhisker[cell][whisker]:
                SWSpikeTimes.append(t)
        
        percentile = 0.99
        poissonProb = poisson.pmf(1, rates[cell])
        nrPWTrials = trialsPerWhisker[cell][5]
        PWThreshold = binom.ppf(percentile, nrPWTrials, poissonProb)
        print 'cell ', cell
        print 'PWThreshold = %d spikes in %d trials' % (PWThreshold, nrPWTrials)
        PWhist, bins = np.histogram(spikesPerWhisker[cell][5], numberOfBins, PSTHrange)
        #if nrPWTrials:
            #PWhist = 1.0/nrPWTrials*PWhist
        tmpLatency = 1000.0
        for i in range(len(PWhist)):
            #if PWhist[i] > latencyThresholds[cell]:
            if PWhist[i] > PWThreshold:
                tmpLatency = i
                break
        PWLatencies[cell] = tmpLatency
        
        SWThreshold = binom.ppf(percentile, nrSWTrials, poissonProb)
        print 'SWThreshold = %d spikes in %d trials' % (SWThreshold, nrSWTrials)
        hist, bins = np.histogram(SWSpikeTimes, numberOfBins, PSTHrange)
        #if nrSWTrials:
            #hist = 1.0/nrSWTrials*hist
        tmpLatency = 1000.0
        for i in range(len(hist)):
            #if hist[i] > latencyThresholds[cell]:
            if hist[i] > SWThreshold:
                tmpLatency = i
                break
        SWLatencies[cell] = tmpLatency\
        
        print 'cell: ', cell
        print 'nr. PW spikes: ', len(spikesPerWhisker[cell][5])
        print 'nr. PW trials: ', nrPWTrials
        print 'nr. SW spikes: ', len(SWSpikeTimes)
        print 'nr. SW trials: ', nrSWTrials
        print 'PW latency: ', PWLatencies[cell]
        print 'SW latency: ', SWLatencies[cell]
        
        #if SWLatencies[cell] > 1000.0:
        if SWLatencies[cell]:
            print 'latencyThreshold = ', latencyThresholds[cell]
            import matplotlib.pyplot as plt
            plt.plot(bins[:-1], hist, 'k')
            plt.plot(bins[:-1], PWhist, 'r')
            plt.show()
    
    #cellNames = PWLatencies.keys()
    #cellNames.sort()
    #cellTypeName = cellTypeFolder.split('/')[-1]
    ##outFileName = os.path.join(cellTypeFolder, 'latency_100ms.csv')
    #outFileName = cellTypeName + '_latency_100ms.csv'
    #with open(outFileName, 'w') as PSTHFile:
        #header = 'Cell name'
        #header += '\t'
        #header += 'PW latency (ms)'
        #header += '\t'
        #header += 'SW latency (ms)'
        #header += '\n'
        #PSTHFile.write(header)
        #for cell in cellNames:
            #line = cell
            #line += '\t'
            #line += str(PWLatencies[cell])
            #line += '\t'
            #line += str(SWLatencies[cell])
            #line += '\n'
            #PSTHFile.write(line)

def load_spike_times(spikeTimesName):
    whiskers = {0: 'B1', 1: 'B2', 2: 'B3',\
                3: 'C1', 4: 'C2', 5: 'C3',\
                6: 'D1', 7: 'D2', 8: 'D3'}
    whiskerDeflectionTrials = {'B1': 0, 'B2': 0, 'B3': 0,\
                                'C1': 0, 'C2': 0, 'C3': 0,\
                                'D1': 0, 'D2': 0, 'D3': 0}
    whiskerSpikeTimes = {'B1': [], 'B2': [], 'B3': [],\
                        'C1': [], 'C2': [], 'C3': [],\
                        'D1': [], 'D2': [], 'D3': []}
    with open(spikeTimesName, 'r') as spikeTimesFile:
        lineCnt = 0
        for line in spikeTimesFile:
            if line:
                lineCnt += 1
            if lineCnt <= 1:
                continue
            splitLine = line.strip().split('\t')
            for i in range(len(splitLine)):
                t = float(splitLine[i])
                whisker = whiskers[i]
                if np.isnan(t):
                    whiskerDeflectionTrials[whisker] += 1
                    continue
                # Use -1 as marker for non-stimulated trials
                if t < 0:
                    continue
                whiskerDeflectionTrials[whisker] += 1
                whiskerSpikeTimes[whisker].append(t)
    
    return whiskerSpikeTimes, whiskerDeflectionTrials

def load_cluster_trials(fname):
    data = np.loadtxt(fname, unpack=True)
    trialNumber = data[1]
    spikeTimes = data[3]
    spikeTimes = 0.1*spikeTimes # CDK files in 0.1ms
    trialsSpikeTimes = {}
    for i in range(len(trialNumber)):
        nr = int(trialNumber[i])
        if not trialsSpikeTimes.has_key(nr):
            trialsSpikeTimes[nr] = []
        t = spikeTimes[i]
        if t > 0.0:
            trialsSpikeTimes[nr].append(t)
    
    return trialsSpikeTimes

def scan_directory(path, fnames, suffix):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory(fname, fnames, suffix)
        elif fname.endswith(suffix):
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 2:
        folderName = sys.argv[1]
        #create_single_cell_PSTH_from_clusters(folderName)
        create_single_cell_latency_from_clusters(folderName)
    elif len(sys.argv) == 3:
        folderName = sys.argv[1]
        outFileName = sys.argv[2]
        create_average_celltype_PSTH_from_clusters(folderName, outFileName)
    else:
        print 'parameters: [folderName] [outFileName]'
    