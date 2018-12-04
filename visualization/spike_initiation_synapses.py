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
        synapseData[synTrialFile] = scp.read_complete_synapse_activation_file(synTrialFile)
    
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

#def synaptic_input_window_analysis(folder, contains, vmSuffix, synapseSuffix, averageSynapseName, outName, synchronyWindow, synchronyWindowOffset):
#def synaptic_input_window_analysis(folder, contains, vmSuffix, synapseSuffix, averageSynapseName, outName):
def synaptic_input_window_analysis(vmData, synData, averageSynsPerCellType, outName, synchronyWindow, synchronyWindowOffset):
    '''
    load all traces, compute spike times
    and create raster plots
    '''
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
    plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
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
    
    tOffset = 245.0
#    tStim = 245.0
    tStim = 253.0 # after VPM activation only
    tStimWindow = 50.0
#    earlyWindow = 25.0
    earlyWindow = 17.0 # after VPM activation only: 25-8
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
            trialSpikeTimes_ = sca.simple_spike_detection(t, v)
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
#            activeSyns = scp.read_complete_synapse_activation_file(synTrialFile)
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
                    for tSyn in synTimesProx:
                        # pre-spike mode
#                        if tStim + totalMeanSpikeTime - correlationWindowBegin <= tSyn < tStim + totalMeanSpikeTime - correlationWindowEnd:
#                            noSpikeTrialSynsEarlyProx[cellType][-1].append(tSyn-tStim)
                        # spike trial absolute timing mode
                        if tOffset + correlationWindowBegin <= tSyn < tOffset + correlationWindowEnd:
                            noSpikeTrialSynsEarlyProx[cellType][-1].append(tSyn-tOffset)
            
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
                    for tSyn in synTimesProx:
                        # pre-spike mode
#                        if tStim + tSpikeReference - correlationWindowBegin <= tSyn < tStim + tSpikeReference - correlationWindowEnd:
#                            spikeTrialSynsEarlyProx[cellType][-1].append(tSyn-tStim)
#                            earlyProxSyns += 1
                        # spike trial absolute timing mode
                        if tOffset + correlationWindowBegin <= tSyn < tOffset + correlationWindowEnd:
                            spikeTrialSynsEarlyProx[cellType][-1].append(tSyn-tOffset)
                            earlyProxSyns += 1
            
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
            header += ' nr. evoked early proximal'
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
                evokedSyn = totalSyn - (correlationWindowEnd-correlationWindowBegin)*averageSynsPerCellType[cellType]
                line += str(totalSyn)
                line += '\t'
                line += str(evokedSyn)
                if cellType == 'EXC':
                    totalDiff += totalSyn
                    evokedDiff += evokedSyn
                if cellType == 'INH':
                    totalDiff -= totalSyn
                    evokedDiff -= evokedSyn
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
                evokedSyn = totalSyn - (correlationWindowEnd-correlationWindowBegin)*averageSynsPerCellType[cellType]
                line += str(totalSyn)
                line += '\t'
                line += str(evokedSyn)
                if cellType == 'EXC':
                    totalDiff += totalSyn
                    evokedDiff += evokedSyn
                if cellType == 'INH':
                    totalDiff -= totalSyn
                    evokedDiff -= evokedSyn
            line += '\n'
            outFile.write(line)
            noSpikeTrialSynsNr.append(totalDiff)
            noSpikeTrialEvokedSynsNr.append(evokedDiff)
    
    synBinWidth = 8.0
    synMin = np.min([np.min(spikeTrialSynsNr),np.min(noSpikeTrialSynsNr)])
    synMax = np.max([np.max(spikeTrialSynsNr),np.max(noSpikeTrialSynsNr)])
    evokedSynMin = np.min([np.min(spikeTrialEvokedSynsNr),np.min(noSpikeTrialEvokedSynsNr)])
    evokedSynMax = np.max([np.max(spikeTrialEvokedSynsNr),np.max(noSpikeTrialEvokedSynsNr)])
    synMinBin = int(synMin/synBinWidth)*synBinWidth
    synMaxBin = (int(synMax/synBinWidth) + 1)*synBinWidth
    synBinsInput = np.arange(synMinBin, synMaxBin, synBinWidth)
    evokedSynMinBin = int(evokedSynMin/synBinWidth)*synBinWidth
    evokedSynMaxBin = (int(evokedSynMax/synBinWidth) + 1)*synBinWidth
    evokedSynBinsInput = np.arange(evokedSynMinBin, evokedSynMaxBin, synBinWidth)
    
    synHistSpikeTrials, synBins = np.histogram(spikeTrialSynsNr, synBinsInput)
    synHistNoSpikeTrials, synBins = np.histogram(noSpikeTrialSynsNr, synBinsInput)
    evokedSynHistSpikeTrials, evokedSynBins = np.histogram(spikeTrialEvokedSynsNr, evokedSynBinsInput)
    evokedSynHistNoSpikeTrials, evokedSynBins = np.histogram(noSpikeTrialEvokedSynsNr, evokedSynBinsInput)
    
    synHistNorm = synHistSpikeTrials + synHistNoSpikeTrials
    synHistProb = synHistSpikeTrials/(synHistNorm + 1e-6)
    evokedSynHistNorm = evokedSynHistSpikeTrials + evokedSynHistNoSpikeTrials
    evokedSynHistProb = evokedSynHistSpikeTrials/(evokedSynHistNorm + 1e-6)
    
    synProbName = ''
    evokedSynProbName = ''
    if not outName.endswith('.csv'):
        synProbName = outName + '_syn_spike_prob.csv'
        evokedSynProbName = outName + '_syn_spike_prob_evoked.csv'
    else:
        synProbName = outName[:-4] + '_syn_spike_prob.csv'
        evokedSynProbName = outName[:-4] + '_syn_spike_prob_evoked.csv'
    
    with open(synProbName, 'w') as outFile:
        header = 'syn bin begin\tsyn bin end\tsyn bin\tnr spike trials\tnr no spike trials\tspike prob\n'
        outFile.write(header)
        for i in range(len(synHistProb)):
            line = str(synBins[i])
            line += '\t'
            line += str(synBins[i+1])
            line += '\t'
            line += str(0.5*(synBins[i]+synBins[i+1]))
            line += '\t'
            line += str(synHistSpikeTrials[i])
            line += '\t'
            line += str(synHistNoSpikeTrials[i])
            line += '\t'
            line += str(synHistProb[i])
            line += '\n'
            outFile.write(line)
    
    with open(evokedSynProbName, 'w') as outFile:
        header = 'evoked syn bin begin\tevoked syn bin end\tevoked syn bin\tnr spike trials\tnr no spike trials\tspike prob\n'
        outFile.write(header)
        for i in range(len(evokedSynHistProb)):
            line = str(evokedSynBins[i])
            line += '\t'
            line += str(evokedSynBins[i+1])
            line += '\t'
            line += str(0.5*(evokedSynBins[i]+evokedSynBins[i+1]))
            line += '\t'
            line += str(evokedSynHistSpikeTrials[i])
            line += '\t'
            line += str(evokedSynHistNoSpikeTrials[i])
            line += '\t'
            line += str(evokedSynHistProb[i])
            line += '\n'
            outFile.write(line)

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
    
    
    
    