#!/usr/bin/python

import sys
import numpy as np
import os, os.path
import glob
import matplotlib.pyplot as plt

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

def create_average_celltype_PSTH_from_clusters(cellTypeFolder):
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
    spikeTrialsPerWhisker = {}
    for cell in cellSpikeTimes:
        trialsPerWhisker[cell] = dict([(i,0) for i in range(1,10)])
        spikesPerWhisker[cell] = dict([(i,[]) for i in range(1,10)])
        spikeTrialsPerWhisker[cell] = dict([(i,[]) for i in range(1,10)])
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
                    tmpTrials += 1
                    trialsPerWhisker[cell][col] += 1
                    for t in cellSpikeTimes[cell][whisker][trial]:
                        if stimulusOnset <= t < PSTHEnd:
                        #if ongoingBegin <= t < PSTHEnd:  # change by MO
                            spikesPerWhisker[cell][col].append(t-stimulusOnset)
                            spikeTrialsPerWhisker[cell][col].append(tmpTrials)
                            tmpSpikes += 1
                        #else:
                            #spikesPerWhisker[cell][col].append(-100.0)
                if col == 5:
                    print cell
                    print 'PW: ', PW
                    print 'APs per stim: ', float(tmpSpikes)/tmpTrials
    
    cellNames = spikesPerWhisker.keys()
    for i in range(len(cellNames)):
        cell = cellNames[i]
        plt.figure(i)
        # PW 
        #plt.subplot(211)
        plt.subplot(111)
        plt.title(cell)
        trials = spikeTrialsPerWhisker[cell][5]
        spikeTimes = spikesPerWhisker[cell][5]
        plt.plot(spikeTimes, trials, 'k|')
        plt.ylabel('PW trial')
        plt.xlim([-25, 25])
        plt.ylim([1,trialsPerWhisker[cell][5]])
        plt.minorticks_on()
        ## SuW 
        #plt.subplot(212)
        #SuWTrials = spikeTrialsPerWhisker[cell][1][:]
        #SuWSpikeTimes = spikesPerWhisker[cell][1][:]
        #totalTrials = trialsPerWhisker[cell][1]
        #for j in range(2, 10):
            #if j == 5:
                #continue
            #tmpTrials = np.array(spikeTrialsPerWhisker[cell][j])
            #if j == 6:
                #tmpTrials += totalTrials
            #else:
                #tmpTrials += totalTrials
            #SuWTrials.extend(tmpTrials)
            #SuWSpikeTimes.extend(spikesPerWhisker[cell][j])
            #totalTrials += trialsPerWhisker[cell][j]
        #plt.plot(SuWSpikeTimes, SuWTrials, 'ko')
        #plt.ylabel('SuW trial')
        #plt.xlim([0,50])
        #plt.ylim([1,totalTrials])
        #plt.minorticks_on()
        #plt.xlabel('Time post-stimulus (ms)')
    
    plt.show()
    
    ## create 1ms resolution PSTH per whisker and subtract
    ## spontaneous activity per 1ms bin
    ##numberOfBins = 50 # change by MO
    #numberOfBins = 175
    ##PSTHrange = (stimulusOnset,PSTHEnd) # change by MO
    #PSTHrange = (ongoingBegin,PSTHEnd)
    #whiskerPSTH = {}
    #for cell in cellSpikeTimes:
        #whiskerPSTH[cell] = {}
        #for col in spikesPerWhisker[cell]:
            #nrOfTrials = trialsPerWhisker[cell][col]
            #hist, bins = np.histogram(spikesPerWhisker[cell][col], numberOfBins, PSTHrange)
            #if nrOfTrials:
                #hist = 1.0/nrOfTrials*hist # change by MO (no substraction of spont rate)
                ##hist = 1.0/nrOfTrials*hist - rates[cell]
            #else:
                #hist = hist # change by MO (no substraction of spont rate)
                ##hist = hist - rates[cell]
            #whiskerPSTH[cell][col] = hist, bins
            #print 'cell: ', cell
            #print 'whisker: ', index2WhiskerLUT[col]
            #print 'nr. of trials: ', nrOfTrials
            #print 'sum PSTH: ', np.sum(hist)
    
    #cellNames = whiskerPSTH.keys()
    #cellNames.sort()
    #whiskers = whiskerPSTH[cellNames[0]].keys()
    #whiskers.sort()
    #if not outFileName.endswith('.csv'):
        #outFileName += '.csv'
    #with open(outFileName, 'w') as PSTHFile:
        #header = 'Cell name'
        #header += '\t'
        #header += 'Time'
        #for whisker in whiskers:
            #header += '\t'
            #header += index2WhiskerLUT[whisker]
        #header += '\n'
        #PSTHFile.write(header)
        #for cell in cellNames:
            #for i in range (numberOfBins):
                #line = cell
                #line+= '\t'
                #line+= str(i+ongoingBegin)
                #for whisker in whiskers:
                    #line += '\t'
                    #PSTH = whiskerPSTH[cell][whisker]
                    ##totalPSTH = 0.0
                    ##if PSTH is not None:
                    #hist, bins = PSTH
                        ##totalPSTH = np.sum(hist) # change by MO
                    ##line += str(totalPSTH)
                    #line += str(hist[i])
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
        create_average_celltype_PSTH_from_clusters(folderName)
    else:
        print 'parameters: [folderName] [outFileName]'
    