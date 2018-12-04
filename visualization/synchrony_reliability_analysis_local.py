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

def automated_analysis(folder, contains, vmSuffix, synapseSuffix, averageSynapseName, outName):
    vmNames = []
#    scan_directory(folder, vmNames, vmSuffix)
    scan_directory2(folder, vmNames, vmSuffix, contains)
    synNames = []
#    scan_directory(folder, synNames, synapseSuffix)
    scan_directory2(folder, synNames, synapseSuffix, contains)
    
    #===========================================================================
    # first, determine avg. ongoing nr. of synapses per cell type
    # later subtracted from total nr. to give evoked nr. of synapses per type
    #===========================================================================
    averageSynsPerCellType = load_average_synapse_file(averageSynapseName)
    
    vmData = {}
    for n in range(len(vmNames)):
        fname = vmNames[n]
        print 'Loading traces from file %s' % fname
        data_ = np.loadtxt(fname, skiprows=1, unpack=True)
        vmData[fname] = data_
    
    print 'Loading %d synapse activation files...' % len(synNames)
    synapseData = {}
    for synTrialFile in synNames:
        synapseData[synTrialFile] = read_complete_synapse_activation_file(synTrialFile)
    
#    for synchronyWindow in range(25,0,-1):
##    for synchronyWindow in range(2,0,-2):
##    for synchronyWindow in range(12,0,-2):
##        for synchronyWindowOffset in range(0,22-synchronyWindow,2):
##        for synchronyWindowOffset in range(0,4,2):
##        for synchronyWindowOffset in range(8,22-synchronyWindow,2):
#        for synchronyWindowOffset in range(0,26-synchronyWindow,1):
#            if synchronyWindow%2 == 0 and synchronyWindowOffset%2 == 0:
#                continue
#            tmpOutName = outName + '_window_' + str(synchronyWindow) + '_offset_' + str(synchronyWindowOffset)
#            print '****************************************'
#            print 'Analyzing synaptic inputs'
#            print 'synchrony window: %dms' % synchronyWindow
#            print 'window offset: %dms' % synchronyWindowOffset
#            print 'output base name: %s' % tmpOutName
#            print '****************************************'
#            synaptic_input_window_analysis(vmData, synapseData, averageSynsPerCellType, tmpOutName, synchronyWindow, synchronyWindowOffset)
    
    # PW: optimal window 8-13ms
    synchronyWindow = 50
    synchronyWindowOffset = 0
    tmpOutName = outName + '_window_' + str(synchronyWindow) + '_offset_' + str(synchronyWindowOffset)
    print '****************************************'
    print 'Analyzing synaptic inputs'
    print 'synchrony window: %dms' % synchronyWindow
    print 'window offset: %dms' % synchronyWindowOffset
    print 'output base name: %s' % tmpOutName
    print '****************************************'
    synaptic_input_window_analysis(vmData, synapseData, averageSynsPerCellType, tmpOutName, synchronyWindow, synchronyWindowOffset)
    #synaptic_input_PCA(vmData, synapseData, averageSynsPerCellType, tmpOutName, synchronyWindow, synchronyWindowOffset)

#def synaptic_input_window_analysis(folder, contains, vmSuffix, synapseSuffix, averageSynapseName, outName, synchronyWindow, synchronyWindowOffset):
#def synaptic_input_window_analysis(folder, contains, vmSuffix, synapseSuffix, averageSynapseName, outName):
def synaptic_input_window_analysis(vmData, synData, averageSynsPerCellType, outName, synchronyWindow, synchronyWindowOffset):
    '''
    load all traces, compute spike times
    and create raster plots
    
    Here: analyze total/evoked number and SD (timing variability)
    of full simulation trials for comparison with number/timing
    of generic synapse.
    Can be done for all cell types separately, or by combining
    different excitatory presynaptic cell types and analyzing
    them as a combined excitatory type.
    '''
    # Generic excitatory synapse analysis
    #excTypes = ('Generic',)
    # all separately
    #excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
    #            'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    # only VPM + L6cc
    excTypes = ('L6cc', 'VPM')
    # only (supra-)granular active types
    #excTypes = ('L34', 'L4sp', 'L4ss')
    # no L6cc
#    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6ccinv', 'L6ct', 'VPM')
    # no VPM
#    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6cc', 'L6ccinv', 'L6ct')
    # no L6cc, VPM
#    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6ccinv', 'L6ct')
    # no L34
#    excTypes = ('L2', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    inhTypes = ('L1','L23Trans','L45Sym','L45Peak','L56Trans',\
                'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
    cellTypeColorMap = {'L2': 'dodgerblue', 'L34': 'blue', 'L4py': 'palegreen',\
                    'L4sp': 'green', 'L4ss': 'lime', 'L5st': 'yellow', 'L5tt': 'orange',\
                    'L6cc': 'indigo', 'L6ccinv': 'violet', 'L6ct': 'magenta', 'VPM': 'black',\
                    'INH': 'grey', 'EXC': 'red'}
    # Generic excitatory synapse analysis
    #plotTypes = ('EXC', 'INH')
    # all separately
    #plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
    #            'L6cc', 'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
    # L6cc + VPM only (use types for evoked nr, and EXC for STD
    plotTypes = ('L6cc', 'VPM', 'EXC', 'INH')
    # only (supra-)granular active types (use types for evoked nr, and EXC for STD
    #plotTypes = ('L34', 'L4sp', 'L4ss', 'EXC', 'INH')
    # no L6cc
#    plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
    # no VPM
#    plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6cc', 'L6ccinv', 'L6ct', 'EXC', 'INH')
    # no L6cc, VPM
#    plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6ccinv', 'L6ct', 'EXC', 'INH')
    # no L34
#    plotTypes = ('L2', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6cc', 'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
    
#    vmNames = []
##    scan_directory(folder, vmNames, vmSuffix)
#    scan_directory2(folder, vmNames, vmSuffix, contains)
#    synNames = []
##    scan_directory(folder, synNames, synapseSuffix)
#    scan_directory2(folder, synNames, synapseSuffix, contains)
#    
#    #===========================================================================
#    # first, determine avg. ongoing nr. of synapses per cell type
#    # later subtracted from total nr. to give evoked nr. of synapses per type
#    #===========================================================================
#    averageSynsPerCellType = load_average_synapse_file(averageSynapseName)
    excAvg = 0.0
    for excType in excTypes:
        excAvg += averageSynsPerCellType[excType]
    averageSynsPerCellType['EXC'] = excAvg
#    
#    print 'Computing reverse correlation of synaptic input from %d files' % len(synNames)
    
    vmNames = vmData.keys()
    synNames = synData.keys()
    
    # look at spikes in early phase (0-20ms)
    # and late phase (20-50ms) separately
    spikeTrialSyns = {}
    spikeTrialSynsProx = {}
    spikeTrialSynsDistal = {}
    spikeTrialSynsEarly = {}
    spikeTrialSynsEarlyProx = {}
    spikeTrialSynsEarlyDistal = {}
    spikeTimesEarly = []
    spikeTrialSynsLate = {}
    spikeTrialSynsLateProx = {}
    spikeTrialSynsLateDistal = {}
    noSpikeTrialSyns = {}
    noSpikeTrialSynsProx = {}
    noSpikeTrialSynsDistal = {}
    noSpikeTrialSynsEarly = {}
    noSpikeTrialSynsEarlyProx = {}
    noSpikeTrialSynsEarlyDistal = {}
    noSpikeTrialSynsLate = {}
    noSpikeTrialSynsLateProx = {}
    noSpikeTrialSynsLateDistal = {}
    spikeTrialSynsTimingVarEarlyProx = {}
    noSpikeTrialSynsTimingVarEarlyProx = {}
    for cellType in plotTypes:
        spikeTrialSyns[cellType] = []
        spikeTrialSynsProx[cellType] = []
        spikeTrialSynsDistal[cellType] = []
        spikeTrialSynsEarly[cellType] = []
        spikeTrialSynsEarlyProx[cellType] = []
        spikeTrialSynsEarlyDistal[cellType] = []
        spikeTrialSynsLate[cellType] = []
        spikeTrialSynsLateProx[cellType] = []
        spikeTrialSynsLateDistal[cellType] = []
        noSpikeTrialSyns[cellType] = []
        noSpikeTrialSynsProx[cellType] = []
        noSpikeTrialSynsDistal[cellType] = []
        noSpikeTrialSynsEarly[cellType] = []
        noSpikeTrialSynsEarlyProx[cellType] = []
        noSpikeTrialSynsEarlyDistal[cellType] = []
        noSpikeTrialSynsLate[cellType] = []
        noSpikeTrialSynsLateProx[cellType] = []
        noSpikeTrialSynsLateDistal[cellType] = []
        spikeTrialSynsTimingVarEarlyProx[cellType] = []
        noSpikeTrialSynsTimingVarEarlyProx[cellType] = []
    
    tOffset = 245.0
#    tStim = 245.0
    tStim = 253.0 # after VPM activation only
    #tStim = 259.0 # PW late phase
    #tStim = 261.0 # SuW late phase
    tStimWindow = 50.0
    offsetWindow50 = 42.0 # 50 - 8ms offset
#    earlyWindow = 25.0
    earlyWindow = 17.0 # after VPM activation only: 25-8
    #earlyWindow = 6.0 # after VPM activation early phase PW (14-8)
    #earlyWindow = 11.0 # after VPM activation late phase PW (25-14)
    #earlyWindow = 8.0 # after VPM activation early phase SuW (16-8)
    #earlyWindow = 9.0 # after VPM activation late phase SuW (25-16)
#    correlationWindowBegin = 7.5
#    correlationWindowEnd = 0.0 # i.e. only synapses in [tSpike-correlationWindowBegin,tSpike-correlationWindowEnd]
#    correlationWindowBegin = synchronyWindow + synchronyWindowOffset
#    correlationWindowEnd = synchronyWindowOffset # i.e. only synapses in [tSpike-correlationWindowBegin,tSpike-correlationWindowEnd]
    correlationWindowBegin = synchronyWindowOffset
    correlationWindowEnd = synchronyWindow + synchronyWindowOffset # i.e. only synapses in [correlationWindowBegin,correlationWindowEnd]
    binWidth = 1.0
    data = []
    trialSpikeTimes = [[] for j in range(len(vmNames))]
    trialWithSpikes = {}
    
    for n in range(len(vmNames)):
        fname = vmNames[n]
#        pathName = fname[:fname.rfind('/')]
#        print 'pathName = %s' % pathName
        print 'Analyzing traces in file %s' % fname
#        data_ = np.loadtxt(fname, skiprows=1, unpack=True)
        data_ = vmData[fname]
        data.append(data_)
        t = data_[0]
        trialWithSpikes[n] = []
        for i in range(1, len(data_)):
            trialSpikeTimes[n].append([])
            v = data_[i]
            trialSpikeTimes_ = simple_spike_detection(t, v)
            trialWithSpikes_ = False
            for tSpike in trialSpikeTimes_:
                if tSpike >= tStim and tSpike < tStim + earlyWindow:
                    trialSpikeTimes[n][-1].append(tSpike-tStim)
                    trialWithSpikes_ = True
            if trialWithSpikes_:
                tSpikeMin = np.min(trialSpikeTimes[n][-1])
                spikeTimesEarly.append(tSpikeMin)
            trialWithSpikes[n].append(trialWithSpikes_)
    
    allTrialSpikeTimes = [tSpike for traceFile in trialSpikeTimes for trace in traceFile for tSpike in trace]
    totalMeanSpikeTime = np.mean(allTrialSpikeTimes)
    print 'total mean spike time = %.1fms' % totalMeanSpikeTime
    
    for n in range(len(vmNames)):
        nrSpikeTrials = 0
        nrNoSpikeTrials = 0
        earlyProxSyns = 0
        for i in range(1, len(data[n])):
#        for i in range(1,10):
            print 'Computing reverse correlation for spikes in trial %d of %d\r' % (i, len(data[n])-1),
            sys.stdout.flush()
            trialNr = i-1
            synTrialStr = 'simulation_run%04d_synapses.csv' % trialNr
            synTrialFile = ''
            tmpVmName = vmNames[n]
            for name in synNames:
                if synTrialStr in name and os.path.split(tmpVmName)[0] == os.path.split(name)[0]:
                    synTrialFile = name
                    break
            if synTrialFile == '':
                errstr = 'Could not find synapse activation file for trial nr. %d' % trialNr
                raise RuntimeError(errstr)
#            activeSyns = read_complete_synapse_activation_file(synTrialFile)
            activeSyns = synData[synTrialFile]
            
            synapseTimes = {}
            synapseTimesProx = {}
            synapseTimesDistal = {}
            for excType in excTypes:
                synapseTimes[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
                synapseTimesProx[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
                synapseTimesDistal[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimes['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimesProx['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimesDistal['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimes['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
            synapseTimesProx['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
            synapseTimesDistal['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
            
            for synType in activeSyns.keys():
                preCellType = synType.split('_')[0]
                for excType in excTypes:
                    if excType == preCellType:
                        for syn in activeSyns[synType]:
                            somaDist = syn[1]
                            structure = syn[4]
                            synTimes = syn[5]
                            synapseTimes[excType][structure].extend(synTimes)
                            synapseTimes[excType]['Total'].extend(synTimes)
                            synapseTimes['EXC'][structure].extend(synTimes)
                            synapseTimes['EXC']['Total'].extend(synTimes)
#                            synapseTimesProx[excType][structure].extend(synTimes)
#                            synapseTimesProx[excType]['Total'].extend(synTimes)
#                            synapseTimesProx['EXC'][structure].extend(synTimes)
#                            synapseTimesProx['EXC']['Total'].extend(synTimes)
                            if somaDist < 500.0:
                                synapseTimesProx[excType][structure].extend(synTimes)
                                synapseTimesProx[excType]['Total'].extend(synTimes)
                                synapseTimesProx['EXC'][structure].extend(synTimes)
                                synapseTimesProx['EXC']['Total'].extend(synTimes)
                            else:
                                synapseTimesDistal[excType][structure].extend(synTimes)
                                synapseTimesDistal[excType]['Total'].extend(synTimes)
                                synapseTimesDistal['EXC'][structure].extend(synTimes)
                                synapseTimesDistal['EXC']['Total'].extend(synTimes)
                for inhType in inhTypes:
                    if inhType == preCellType:
                        for syn in activeSyns[synType]:
                            somaDist = syn[1]
                            structure = syn[4]
                            synTimes = syn[5]
                            synapseTimes['INH'][structure].extend(synTimes)
                            synapseTimes['INH']['Total'].extend(synTimes)
#                            synapseTimesProx['INH'][structure].extend(synTimes)
#                            synapseTimesProx['INH']['Total'].extend(synTimes)
                            if somaDist < 500.0:
                                synapseTimesProx['INH'][structure].extend(synTimes)
                                synapseTimesProx['INH']['Total'].extend(synTimes)
                            else:
                                synapseTimesDistal['INH'][structure].extend(synTimes)
                                synapseTimesDistal['INH']['Total'].extend(synTimes)
            
            if not trialWithSpikes[n][trialNr]:
                nrNoSpikeTrials += 1
                for cellType in plotTypes:
                    noSpikeTrialSyns[cellType].append([])
                    noSpikeTrialSynsProx[cellType].append([])
                    noSpikeTrialSynsDistal[cellType].append([])
                    noSpikeTrialSynsEarly[cellType].append([])
                    noSpikeTrialSynsEarlyProx[cellType].append([])
                    noSpikeTrialSynsEarlyDistal[cellType].append([])
                    noSpikeTrialSynsLate[cellType].append([])
                    noSpikeTrialSynsLateProx[cellType].append([])
                    noSpikeTrialSynsLateDistal[cellType].append([])
                    synTimes = synapseTimes[cellType]['Total']
                    synTimesProx = synapseTimesProx[cellType]['Total']
                    synTimesDistal = synapseTimesDistal[cellType]['Total']
                    tSynTmpList = []
                    for tSyn in synTimesProx:
                        # pre-spike mode
#                        if tStim + totalMeanSpikeTime - correlationWindowBegin <= tSyn < tStim + totalMeanSpikeTime - correlationWindowEnd:
#                            noSpikeTrialSynsEarlyProx[cellType][-1].append(tSyn-tStim)
                        # spike trial absolute timing mode
#                        if tOffset + correlationWindowBegin <= tSyn < tOffset + correlationWindowEnd:
#                            noSpikeTrialSynsEarlyProx[cellType][-1].append(tSyn-tOffset)
                        # spike trial 0-50ms all synapses mode
                        if tStim <= tSyn < tStim + offsetWindow50:
                            noSpikeTrialSynsEarlyProx[cellType][-1].append(tSyn-tStim)
                        if tStim <= tSyn < tStim + offsetWindow50:
                            tSynTmpList.append(tSyn)
                    noSpikeTrialSynsTimingVarEarlyProx[cellType].append(np.std(tSynTmpList))
            
            elif trialWithSpikes[n][trialNr]:
                nrSpikeTrials += 1
                tSpikeReference = np.min(trialSpikeTimes[n][trialNr])
                for cellType in plotTypes:
                    spikeTrialSyns[cellType].append([])
                    spikeTrialSynsProx[cellType].append([])
                    spikeTrialSynsDistal[cellType].append([])
                    spikeTrialSynsEarly[cellType].append([])
                    spikeTrialSynsEarlyProx[cellType].append([])
                    spikeTrialSynsEarlyDistal[cellType].append([])
                    spikeTrialSynsLate[cellType].append([])
                    spikeTrialSynsLateProx[cellType].append([])
                    spikeTrialSynsLateDistal[cellType].append([])
                    synTimes = synapseTimes[cellType]['Total']
                    synTimesProx = synapseTimesProx[cellType]['Total']
                    synTimesDistal = synapseTimesDistal[cellType]['Total']
                    tSynTmpList = []
                    for tSyn in synTimesProx:
                        # pre-spike mode
#                        if tStim + tSpikeReference - correlationWindowBegin <= tSyn < tStim + tSpikeReference - correlationWindowEnd:
#                            spikeTrialSynsEarlyProx[cellType][-1].append(tSyn-tStim)
#                            earlyProxSyns += 1
                        # spike trial absolute timing mode
#                        if tOffset + correlationWindowBegin <= tSyn < tOffset + correlationWindowEnd:
#                            spikeTrialSynsEarlyProx[cellType][-1].append(tSyn-tOffset)
#                            earlyProxSyns += 1
                        # spike trial 0-50ms all synapses mode
                        if tStim <= tSyn < tStim + offsetWindow50:
                            spikeTrialSynsEarlyProx[cellType][-1].append(tSyn-tStim)
                            earlyProxSyns += 1
                        if tStim <= tSyn < tStim + offsetWindow50:
                            tSynTmpList.append(tSyn)
                    spikeTrialSynsTimingVarEarlyProx[cellType].append(np.std(tSynTmpList))
            
        print ''
        print 'Nr of trials with spike: %d' % nrSpikeTrials
        print 'Nr of trials without spike: %d' % nrNoSpikeTrials
        print 'earlyProxSyns = %d' % earlyProxSyns
        print 'mean spike time = %.1fms' % np.mean([tSpike for trace in trialSpikeTimes[n] for tSpike in trace])
    
    spikeTrialsName = ''
    noSpikeTrialsName = ''
    allTrialsName = ''
    if not outName.endswith('.csv'):
        spikeTrialsName = outName + '_spike_trials.csv'
        noSpikeTrialsName = outName + '_no_spike_trials.csv'
        allTrialsName = outName + '_all_trials.csv'
    else:
        spikeTrialsName = outName[:-4] + '_spike_trials.csv'
        noSpikeTrialsName = outName[:-4] + '_no_spike_trials.csv'
        allTrialsName = outName[:-4] + '_all_trials.csv'
    
    spikeTrialSynsNr = []
    spikeTrialEvokedSynsNr = []
    noSpikeTrialSynsNr = []
    noSpikeTrialEvokedSynsNr = []
    
    with open(allTrialsName, 'w') as outFile:
        header = 'Trial nr.'
        header += '\tspike (1/0)\tspike time (ms)'
        for cellType in plotTypes:
            header += '\t'
            header += cellType
            header += ' nr. early proximal\t'
            header += cellType
            header += ' nr. evoked early proximal\t'
            header += cellType
            header += ' STD early proximal\t'
            header += cellType
            header += ' MED time early proximal'
        header += '\n'
        outFile.write(header)
        nrOfTrials = len(spikeTrialSyns[plotTypes[0]])
        for i in range(nrOfTrials):
            line = str(i)
            line += '\t1\t'
            line += str(spikeTimesEarly[i])
            totalDiff = 0.0
            evokedDiff = 0.0
            for cellType in plotTypes:
                line += '\t'
                totalSyn = len(spikeTrialSynsEarlyProx[cellType][i])
                evokedSyn = totalSyn
                #if cellType == 'INH':
                    #evokedSyn -= (correlationWindowEnd-correlationWindowBegin)*averageSynsPerCellType[cellType]
                evokedSyn -= offsetWindow50*averageSynsPerCellType[cellType]
                line += str(totalSyn)
                line += '\t'
                line += str(evokedSyn)
                #if cellType == 'EXC':
                #    totalDiff += totalSyn
                #    evokedDiff += evokedSyn
                #if cellType == 'INH':
                totalDiff -= totalSyn
                evokedDiff -= evokedSyn
                # estimate of evoked timing var
                totalHist, bins = np.histogram(spikeTrialSynsEarlyProx[cellType][i], bins=range(int(offsetWindow50+1)))
                evokedHist = totalHist - averageSynsPerCellType[cellType]
                evokedHist = evokedHist*(evokedHist>=0)
                norm = np.sum(evokedHist)
                if norm:
                    meanTime = np.dot(bins[:-1], evokedHist)/norm
                    diff = 0.0
                    for j in range(len(evokedHist)):
                        diff += evokedHist[j]*(j - meanTime)**2
                    SD = np.sqrt(diff/(len(evokedHist)-1))
                    medBin = 0
                    cumSum = 0.0
                    for j in range(len(evokedHist)):
                        cumSum += evokedHist[j]
                        if cumSum >= 0.5*norm:
                            medBin = j
                            break
                    medTime = medBin + 1.0
                else:
                    SD = np.nan
                    medTime = np.nan
                line += '\t'
                line += str(SD)
                line += '\t'
                line += str(medTime)
                #line += str(spikeTrialSynsTimingVarEarlyProx[cellType][i])
            line += '\n'
            outFile.write(line)
            spikeTrialSynsNr.append(totalDiff)
            spikeTrialEvokedSynsNr.append(evokedDiff)
        nrOfTrials = len(noSpikeTrialSyns[plotTypes[0]])
        for i in range(nrOfTrials):
            line = str(i)
            line += '\t0\t'
            line += str(totalMeanSpikeTime)
            totalDiff = 0.0
            evokedDiff = 0.0
            for cellType in plotTypes:
                line += '\t'
                totalSyn = len(noSpikeTrialSynsEarlyProx[cellType][i])
                evokedSyn = totalSyn
                #if cellType == 'INH':
                    #evokedSyn -= (correlationWindowEnd-correlationWindowBegin)*averageSynsPerCellType[cellType]
                evokedSyn -= offsetWindow50*averageSynsPerCellType[cellType]
                line += str(totalSyn)
                line += '\t'
                line += str(evokedSyn)
                #if cellType == 'EXC':
                #    totalDiff += totalSyn
                #    evokedDiff += evokedSyn
                #if cellType == 'INH':
                totalDiff -= totalSyn
                evokedDiff -= evokedSyn
                # estimate of evoked timing var
                totalHist, bins = np.histogram(noSpikeTrialSynsEarlyProx[cellType][i], bins=range(int(offsetWindow50+1)))
                evokedHist = totalHist - averageSynsPerCellType[cellType]
                evokedHist = evokedHist*(evokedHist>=0)
                norm = np.sum(evokedHist)
                if norm:
                    meanTime = np.dot(bins[:-1], evokedHist)/norm
                    diff = 0.0
                    for j in range(len(evokedHist)):
                        diff += evokedHist[j]*(j - meanTime)**2
                    SD = np.sqrt(diff/(len(evokedHist)-1))
                    medBin = 0
                    cumSum = 0.0
                    for j in range(len(evokedHist)):
                        cumSum += evokedHist[j]
                        if cumSum >= 0.5*norm:
                            medBin = j
                            break
                    medTime = medBin + 1.0
                else:
                    SD = np.nan
                    medTime = np.nan
                line += '\t'
                line += str(SD)
                line += '\t'
                line += str(medTime)
                #line += str(noSpikeTrialSynsTimingVarEarlyProx[cellType][i])
            line += '\n'
            outFile.write(line)
            noSpikeTrialSynsNr.append(totalDiff)
            noSpikeTrialEvokedSynsNr.append(evokedDiff)
    

def synaptic_input_PCA(vmData, synData, averageSynsPerCellType, outName, synchronyWindow, synchronyWindowOffset):
    '''
    calculate PCA of temporal synaptic input patterns for all trials
    with/without spikes to identify best discriminating dimensions
    Parameterization: Per trial, E/I syn. proximal/distal in 1ms bins (0-25ms)
    --> 4*25-dimensional space
    '''
    #excTypes = ('Generic',)
    # all
    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    # no L6cc
#    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6ccinv', 'L6ct', 'VPM')
    # no VPM
#    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6cc', 'L6ccinv', 'L6ct')
    # no L6cc, VPM
#    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6ccinv', 'L6ct')
    # no L34
#    excTypes = ('L2', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    inhTypes = ('L1','L23Trans','L45Sym','L45Peak','L56Trans',\
                'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
    cellTypeColorMap = {'L2': 'dodgerblue', 'L34': 'blue', 'L4py': 'palegreen',\
                    'L4sp': 'green', 'L4ss': 'lime', 'L5st': 'yellow', 'L5tt': 'orange',\
                    'L6cc': 'indigo', 'L6ccinv': 'violet', 'L6ct': 'magenta', 'VPM': 'black',\
                    'INH': 'grey', 'EXC': 'red'}
    # Generic excitatory synapse analysis
    plotTypes = ('EXC', 'INH')
    # all
    #plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                #'L6cc', 'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
    # no L6cc
#    plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
    # no VPM
#    plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6cc', 'L6ccinv', 'L6ct', 'EXC', 'INH')
    # no L6cc, VPM
#    plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6ccinv', 'L6ct', 'EXC', 'INH')
    # no L34
#    plotTypes = ('L2', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
#                'L6cc', 'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
    
#    vmNames = []
##    scan_directory(folder, vmNames, vmSuffix)
#    scan_directory2(folder, vmNames, vmSuffix, contains)
#    synNames = []
##    scan_directory(folder, synNames, synapseSuffix)
#    scan_directory2(folder, synNames, synapseSuffix, contains)
#    
#    #===========================================================================
#    # first, determine avg. ongoing nr. of synapses per cell type
#    # later subtracted from total nr. to give evoked nr. of synapses per type
#    #===========================================================================
#    averageSynsPerCellType = load_average_synapse_file(averageSynapseName)
#    
#    print 'Computing reverse correlation of synaptic input from %d files' % len(synNames)
    
    vmNames = vmData.keys()
    synNames = synData.keys()
    
    # look at spikes in early phase (0-20ms)
    # and late phase (20-50ms) separately
    spikeTrialSyns = {}
    spikeTrialSynsProx = {}
    spikeTrialSynsDistal = {}
    spikeTrialSynsEarly = {}
    spikeTrialSynsEarlyProx = {}
    spikeTrialSynsEarlyDistal = {}
    spikeTimesEarly = []
    spikeTrialSynsLate = {}
    spikeTrialSynsLateProx = {}
    spikeTrialSynsLateDistal = {}
    noSpikeTrialSyns = {}
    noSpikeTrialSynsProx = {}
    noSpikeTrialSynsDistal = {}
    noSpikeTrialSynsEarly = {}
    noSpikeTrialSynsEarlyProx = {}
    noSpikeTrialSynsEarlyDistal = {}
    noSpikeTrialSynsLate = {}
    noSpikeTrialSynsLateProx = {}
    noSpikeTrialSynsLateDistal = {}
    spikeTrialSynsTimingVarEarlyProx = {}
    noSpikeTrialSynsTimingVarEarlyProx = {}
    for cellType in plotTypes:
        spikeTrialSyns[cellType] = []
        spikeTrialSynsProx[cellType] = []
        spikeTrialSynsDistal[cellType] = []
        spikeTrialSynsEarly[cellType] = []
        spikeTrialSynsEarlyProx[cellType] = []
        spikeTrialSynsEarlyDistal[cellType] = []
        spikeTrialSynsLate[cellType] = []
        spikeTrialSynsLateProx[cellType] = []
        spikeTrialSynsLateDistal[cellType] = []
        noSpikeTrialSyns[cellType] = []
        noSpikeTrialSynsProx[cellType] = []
        noSpikeTrialSynsDistal[cellType] = []
        noSpikeTrialSynsEarly[cellType] = []
        noSpikeTrialSynsEarlyProx[cellType] = []
        noSpikeTrialSynsEarlyDistal[cellType] = []
        noSpikeTrialSynsLate[cellType] = []
        noSpikeTrialSynsLateProx[cellType] = []
        noSpikeTrialSynsLateDistal[cellType] = []
        spikeTrialSynsTimingVarEarlyProx[cellType] = []
        noSpikeTrialSynsTimingVarEarlyProx[cellType] = []
    
    tOffset = 245.0
#    tStim = 245.0
    tStim = 253.0 # after VPM activation only
    #tStim = 259.0 # PW late phase
    #tStim = 261.0 # SuW late phase
    tStimWindow = 50.0
    earlySynWindow = 25.0
#    earlyWindow = 25.0
    earlyWindow = 17.0 # after VPM activation only: 25-8
    #earlyWindow = 6.0 # after VPM activation early phase PW (14-8)
    #earlyWindow = 11.0 # after VPM activation late phase PW (25-14)
    #earlyWindow = 8.0 # after VPM activation early phase SuW (16-8)
    #earlyWindow = 9.0 # after VPM activation late phase SuW (25-16)
#    correlationWindowBegin = 7.5
#    correlationWindowEnd = 0.0 # i.e. only synapses in [tSpike-correlationWindowBegin,tSpike-correlationWindowEnd]
#    correlationWindowBegin = synchronyWindow + synchronyWindowOffset
#    correlationWindowEnd = synchronyWindowOffset # i.e. only synapses in [tSpike-correlationWindowBegin,tSpike-correlationWindowEnd]
    correlationWindowBegin = synchronyWindowOffset
    correlationWindowEnd = synchronyWindow + synchronyWindowOffset # i.e. only synapses in [correlationWindowBegin,correlationWindowEnd]
    binWidth = 1.0
    data = []
    trialSpikeTimes = [[] for j in range(len(vmNames))]
    trialWithSpikes = {}
    
    for n in range(len(vmNames)):
        fname = vmNames[n]
#        pathName = fname[:fname.rfind('/')]
#        print 'pathName = %s' % pathName
        print 'Analyzing traces in file %s' % fname
#        data_ = np.loadtxt(fname, skiprows=1, unpack=True)
        data_ = vmData[fname]
        data.append(data_)
        t = data_[0]
        trialWithSpikes[n] = []
        for i in range(1, len(data_)):
            trialSpikeTimes[n].append([])
            v = data_[i]
            trialSpikeTimes_ = simple_spike_detection(t, v)
            trialWithSpikes_ = False
            for tSpike in trialSpikeTimes_:
                if tSpike >= tStim and tSpike < tStim + earlyWindow:
                    trialSpikeTimes[n][-1].append(tSpike-tStim)
                    trialWithSpikes_ = True
            if trialWithSpikes_:
                tSpikeMin = np.min(trialSpikeTimes[n][-1])
                spikeTimesEarly.append(tSpikeMin)
            trialWithSpikes[n].append(trialWithSpikes_)
    
    allTrialSpikeTimes = [tSpike for traceFile in trialSpikeTimes for trace in traceFile for tSpike in trace]
    totalMeanSpikeTime = np.mean(allTrialSpikeTimes)
    print 'total mean spike time = %.1fms' % totalMeanSpikeTime
    
    for n in range(len(vmNames)):
        nrSpikeTrials = 0
        nrNoSpikeTrials = 0
        for i in range(1, len(data[n])):
#        for i in range(1,10):
            print 'Computing reverse correlation for spikes in trial %d of %d\r' % (i, len(data[n])-1),
            sys.stdout.flush()
            trialNr = i-1
            synTrialStr = 'simulation_run%04d_synapses.csv' % trialNr
            synTrialFile = ''
            tmpVmName = vmNames[n]
            for name in synNames:
                if synTrialStr in name and os.path.split(tmpVmName)[0] == os.path.split(name)[0]:
                    synTrialFile = name
                    break
            if synTrialFile == '':
                errstr = 'Could not find synapse activation file for trial nr. %d' % trialNr
                raise RuntimeError(errstr)
#            activeSyns = read_complete_synapse_activation_file(synTrialFile)
            activeSyns = synData[synTrialFile]
            
            synapseTimes = {}
            synapseTimesProx = {}
            synapseTimesDistal = {}
            for excType in excTypes:
                synapseTimes[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
                synapseTimesProx[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
                synapseTimesDistal[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimes['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimesProx['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimesDistal['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimes['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
            synapseTimesProx['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
            synapseTimesDistal['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
            
            for synType in activeSyns.keys():
                preCellType = synType.split('_')[0]
                for excType in excTypes:
                    if excType == preCellType:
                        for syn in activeSyns[synType]:
                            somaDist = syn[1]
                            structure = syn[4]
                            synTimes = syn[5]
                            #synapseTimes[excType][structure].extend(synTimes)
                            #synapseTimes[excType]['Total'].extend(synTimes)
                            #synapseTimes['EXC'][structure].extend(synTimes)
                            #synapseTimes['EXC']['Total'].extend(synTimes)
#                            synapseTimesProx[excType][structure].extend(synTimes)
#                            synapseTimesProx[excType]['Total'].extend(synTimes)
#                            synapseTimesProx['EXC'][structure].extend(synTimes)
#                            synapseTimesProx['EXC']['Total'].extend(synTimes)
                            if somaDist < 500.0:
                                #synapseTimesProx[excType][structure].extend(synTimes)
                                #synapseTimesProx[excType]['Total'].extend(synTimes)
                                #synapseTimesProx['EXC'][structure].extend(synTimes)
                                synapseTimesProx['EXC']['Total'].extend(synTimes)
                            else:
                                #synapseTimesDistal[excType][structure].extend(synTimes)
                                #synapseTimesDistal[excType]['Total'].extend(synTimes)
                                #synapseTimesDistal['EXC'][structure].extend(synTimes)
                                synapseTimesDistal['EXC']['Total'].extend(synTimes)
                for inhType in inhTypes:
                    if inhType == preCellType:
                        for syn in activeSyns[synType]:
                            somaDist = syn[1]
                            structure = syn[4]
                            synTimes = syn[5]
                            #synapseTimes['INH'][structure].extend(synTimes)
                            #synapseTimes['INH']['Total'].extend(synTimes)
#                            synapseTimesProx['INH'][structure].extend(synTimes)
#                            synapseTimesProx['INH']['Total'].extend(synTimes)
                            if somaDist < 500.0:
                                #synapseTimesProx['INH'][structure].extend(synTimes)
                                synapseTimesProx['INH']['Total'].extend(synTimes)
                            else:
                                #synapseTimesDistal['INH'][structure].extend(synTimes)
                                synapseTimesDistal['INH']['Total'].extend(synTimes)
            
            if not trialWithSpikes[n][trialNr]:
                nrNoSpikeTrials += 1
                for cellType in plotTypes:
                    #noSpikeTrialSyns[cellType].append([])
                    #noSpikeTrialSynsProx[cellType].append([])
                    #noSpikeTrialSynsDistal[cellType].append([])
                    #noSpikeTrialSynsEarly[cellType].append([])
                    #noSpikeTrialSynsLate[cellType].append([])
                    #noSpikeTrialSynsLateProx[cellType].append([])
                    #noSpikeTrialSynsLateDistal[cellType].append([])
                    #synTimes = synapseTimes[cellType]['Total']
                    synTimesProx = synapseTimesProx[cellType]['Total']
                    synTimesDistal = synapseTimesDistal[cellType]['Total']
                    tSynProxTmpList = []
                    tSynDistalTmpList = []
                    for tSyn in synTimesProx:
                        if tOffset <= tSyn < tOffset + earlySynWindow:
                            tSynProxTmpList.append(tSyn-tOffset)
                    for tSyn in synTimesDistal:
                        if tOffset <= tSyn < tOffset + earlySynWindow:
                            tSynDistalTmpList.append(tSyn-tOffset)
                    proxHist, tmpBins = np.histogram(tSynProxTmpList, bins=range(26))
                    distalHist, tmpBins = np.histogram(tSynDistalTmpList, bins=range(26))
                    noSpikeTrialSynsEarlyProx[cellType].append(proxHist)
                    noSpikeTrialSynsEarlyDistal[cellType].append(distalHist)
            
            elif trialWithSpikes[n][trialNr]:
                nrSpikeTrials += 1
                tSpikeReference = np.min(trialSpikeTimes[n][trialNr])
                for cellType in plotTypes:
                    #spikeTrialSyns[cellType].append([])
                    #spikeTrialSynsProx[cellType].append([])
                    #spikeTrialSynsDistal[cellType].append([])
                    #spikeTrialSynsEarly[cellType].append([])
                    #spikeTrialSynsLate[cellType].append([])
                    #spikeTrialSynsLateProx[cellType].append([])
                    #spikeTrialSynsLateDistal[cellType].append([])
                    #synTimes = synapseTimes[cellType]['Total']
                    synTimesProx = synapseTimesProx[cellType]['Total']
                    synTimesDistal = synapseTimesDistal[cellType]['Total']
                    tSynProxTmpList = []
                    tSynDistalTmpList = []
                    for tSyn in synTimesProx:
                        if tOffset <= tSyn < tOffset + earlySynWindow:
                            tSynProxTmpList.append(tSyn-tOffset)
                    for tSyn in synTimesDistal:
                        if tOffset <= tSyn < tOffset + earlySynWindow:
                            tSynDistalTmpList.append(tSyn-tOffset)
                    proxHist, tmpBins = np.histogram(tSynProxTmpList, bins=range(26))
                    distalHist, tmpBins = np.histogram(tSynDistalTmpList, bins=range(26))
                    spikeTrialSynsEarlyProx[cellType].append(proxHist)
                    spikeTrialSynsEarlyDistal[cellType].append(distalHist)
            
        print ''
        print 'Nr of trials with spike: %d' % nrSpikeTrials
        print 'Nr of trials without spike: %d' % nrNoSpikeTrials
        print 'mean spike time = %.1fms' % np.mean([tSpike for trace in trialSpikeTimes[n] for tSpike in trace])
    
    trialSpikeList = []
    totalHistList = []
    totalSpikeTrials = len(spikeTrialSynsEarlyProx['EXC'])
    totalNoSpikeTrials = len(noSpikeTrialSynsEarlyProx['EXC'])
    # concatenate E/I/prox/distal histograms in the order:
    # prox. E/prox. I/distal E/distal I
    for i in range(totalSpikeTrials):
        trialSpikeList.append(1)
        trialHist = []
        trialHist.extend(spikeTrialSynsEarlyProx['EXC'][i])
        trialHist.extend(spikeTrialSynsEarlyProx['INH'][i])
        trialHist.extend(spikeTrialSynsEarlyDistal['EXC'][i])
        trialHist.extend(spikeTrialSynsEarlyDistal['INH'][i])
        totalHistList.append(trialHist)
    for i in range(totalNoSpikeTrials):
        trialSpikeList.append(0)
        trialHist = []
        trialHist.extend(noSpikeTrialSynsEarlyProx['EXC'][i])
        trialHist.extend(noSpikeTrialSynsEarlyProx['INH'][i])
        trialHist.extend(noSpikeTrialSynsEarlyDistal['EXC'][i])
        trialHist.extend(noSpikeTrialSynsEarlyDistal['INH'][i])
        totalHistList.append(trialHist)
    
    allData = np.array(totalHistList)
    dataMean = np.mean(allData)
    allData = allData - dataMean
    eigenvectors, eigenvals, V = np.linalg.svd(allData.T, full_matrices=False)
    projectedData = np.dot(allData, eigenvectors).transpose()
    PC1LoadVec = eigenvectors.transpose()[0]
    proxELoad = PC1LoadVec[:25]
    proxILoad = PC1LoadVec[25:50]
    distalELoad = PC1LoadVec[50:75]
    distalILoad = PC1LoadVec[75:]
#    plt.figure(1)
#    plt.plot(range(25), proxELoad, label='Prox E')
#    plt.plot(range(25), proxILoad, label='Prox I')
#    plt.plot(range(25), distalELoad, label='Distal E')
#    plt.plot(range(25), distalILoad, label='Distal I')
#    plt.legend()
#    plt.xlabel('Time post-stimulus (ms)')
#    plt.ylabel('Load')
#    plt.savefig(outName+'_PC1Load.pdf')
    PC2LoadVec = eigenvectors.transpose()[1]
    proxELoad2 = PC2LoadVec[:25]
    proxILoad2 = PC2LoadVec[25:50]
    distalELoad2 = PC2LoadVec[50:75]
    distalILoad2 = PC2LoadVec[75:]
    0/0
    
    #with open(outName+'_PC1_PC2_Load.csv', 'w') as outFile1:
        #header = 'time (ms)\tPC1 E prox load\tPC1 I prox load\tPC1 E distal load\tPC1 I distal load\t'
        #header += 'PC2 E prox load\tPC2 I prox load\tPC2 E distal load\tPC2 I distal load\n'
        #outFile1.write(header)
        #for i in range(25):
            #line = str(i+0.5)
            #line += '\t'
            #line += str(proxELoad[i])
            #line += '\t'
            #line += str(proxILoad[i])
            #line += '\t'
            #line += str(distalELoad[i])
            #line += '\t'
            #line += str(distalILoad[i])
            #line += '\t'
            #line += str(proxELoad2[i])
            #line += '\t'
            #line += str(proxILoad2[i])
            #line += '\t'
            #line += str(distalELoad2[i])
            #line += '\t'
            #line += str(distalILoad2[i])
            #line += '\n'
            #outFile1.write(line)
    #with open(outName+'_PC1_PC2.csv', 'w') as outFile2:
        #header = 'spike trial 1/0\tPC1\tPC2\n'
        #outFile2.write(header)
        #for i in range(len(projectedData[0])):
            #line = str(trialSpikeList[i])
            #line += '\t'
            #line += str(projectedData[0][i])
            #line += '\t'
            #line += str(projectedData[1][i])
            #line += '\n'
            #outFile2.write(line)
    
#    plt.figure(2)
#    plt.plot(range(25), proxELoad2, label='Prox E')
#    plt.plot(range(25), proxILoad2, label='Prox I')
#    plt.plot(range(25), distalELoad2, label='Distal E')
#    plt.plot(range(25), distalILoad2, label='Distal I')
#    plt.legend()
#    plt.xlabel('Time post-stimulus (ms)')
#    plt.ylabel('Load')
#    plt.savefig(outName+'_PC1Load.pdf')
#    plt.figure(3)
#    plt.scatter(projectedData[0], projectedData[1], c=trialSpikeList, cmap='Blues')
#    plt.xlabel('PC 1')
#    plt.ylabel('PC 2')
#    plt.savefig(outName+'_PC1_PC2.pdf')
    #plt.show()

def simple_spike_detection(t, v, tBegin=None, tEnd=None, threshold=0.0, mode='regular'):
    '''
    Simple spike detection method. Identify
    spike times within optional window [tBegin, tEnd]
    by determining threshold crossing times from below
    supported modes:
    regular: absolute threshold crossing
    differential: threshold crossing of dv/dt
    '''
    if len(t) != len(v):
        errstr = 'Dimensions of time vector and membrane potential vector not matching'
        raise RuntimeError(errstr)
    
    tSpike = []
    beginIndex = 1
    endIndex = len(t)
    if tBegin is not None:
        for i in range(1,len(t)):
            if t[i-1] < tBegin and t[i] >= tBegin:
                beginIndex = i
                break
    if tEnd is not None:
        for i in range(1,len(t)):
            if t[i-1] < tEnd and t[i] >= tEnd:
                endIndex = i
                break
    
    if mode == 'regular':
        for i in range(beginIndex,endIndex):
            if v[i-1] < threshold and v[i] >= threshold:
                tSpike.append(t[i])
    
    if mode == 'slope':
        dvdt = np.diff(v)/np.diff(t)
        for i in range(beginIndex,endIndex):
            if dvdt[i-1] < threshold and dvdt[i] >= threshold:
                tSpike.append(t[i])
    
    return tSpike

def read_complete_synapse_activation_file(fname):
    '''
    reads list of all functional synapses and their activation times.
    Input: file of format:
        synapse type\\tsynapse ID\\tsoma distance\\tsection ID\\tsection pt ID\\tdendrite label\\tactivation times
    returns: dictionary with cell types as keys and list of synapse locations and activation times,
    coded as tuples: (synapse ID, soma distance, section ID, point ID, structure label, [t1, t2, ... , tn])
    '''
    synapses = {}
    with open(fname, 'r') as synFile:
        for line in synFile:
            line = line.strip()
            if not line:
                continue
            if line[0] == '#':
                continue
            splitLine = line.split('\t')
            cellType = splitLine[0]
            synID = int(splitLine[1])
            somaDist = float(splitLine[2])
            secID = int(splitLine[3])
            ptID = int(splitLine[4])
            structure = splitLine[5]
            synTimes = []
            synTimesStr = splitLine[6].split(',')
            for tStr in synTimesStr:
                if tStr:
                    synTimes.append(float(tStr))
            if not synapses.has_key(cellType):
                synapses[cellType] = [(synID, somaDist, secID, ptID, structure, synTimes)]
            else:
                synapses[cellType].append((synID, somaDist, secID, ptID, structure, synTimes))
    
    return synapses

def load_average_synapse_file(fname):
    print 'Loading average synaptic input from file %s' % fname
    print 'Computing average synaptic input per cell type per ms between 20-120ms'
    averageSynapses = {}
    with open(fname) as averageSynapseFile:
        lineCnt = 0
        for line in averageSynapseFile:
            stripLine = line.strip()
            if not stripLine:
                continue
            splitLine = line.split('\t')
            if lineCnt:
                cellType = splitLine[0]
                tmpSyn = 0.0
                for i in range(5, 25):
                    tmpSyn += float(splitLine[i])
                tmpSyn = 0.01*tmpSyn
                averageSynapses[cellType] = tmpSyn
            lineCnt += 1
    
    excSyn = 0.0
    for cellType in averageSynapses:
        if cellType != 'INH':
            excSyn += averageSynapses[cellType]
    averageSynapses['EXC'] = excSyn
    return averageSynapses

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
#        elif fname.endswith(suffix) and fname.find(contains) != -1 and fname.find('C2_evoked') == -1 and fname.find('E2_evoked') == -1\
#            and fname.find('L6cc_timing') == -1 and fname.find('dynamic_syns') == -1 and fname.find('L5tt_evoked_inact') == -1:
#        elif fname.endswith(suffix) and fname.find(contains) != -1 and fname.find('L6cc_timing') == -1 and fname.find('dynamic_syns') == -1 and fname.find('L5tt_evoked_inact') == -1:
#        elif fname.endswith(suffix) and fname.find(contains) != -1 and fname.find('C2_evoked') == -1 and fname.find('E2_evoked') == -1:
#        elif fname.endswith(suffix) and fname.find(contains) != -1 and fname.find('E2_evoked') == -1:
        elif fname.endswith(suffix) and fname.find(contains) != -1:
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 7:
        folder = sys.argv[1]
        contains = sys.argv[2]
        vmSuffix = sys.argv[3]
        synapseSuffix = sys.argv[4]
        avgSynName = sys.argv[5]
        outName = sys.argv[6]
#        synchronyWindow = float(sys.argv[7])
#        synchronyWindowOffset = float(sys.argv[8])
#        synaptic_input_window_analysis(folder, contains, vmSuffix, synapseSuffix, avgSynName, outName)
#        automated_analysis(folder, contains, vmSuffix, synapseSuffix, avgSynName, outName, synchronyWindow, synchronyWindowOffset)
        automated_analysis(folder, contains, vmSuffix, synapseSuffix, avgSynName, outName)
    else:
        print 'Error! Number of arguments is %d; should be 6' % (len(sys.argv)-1)
    
    
    
    
