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

def reverse_correlation_synaptic_input(folder, contains, vmSuffix, synapseSuffix):
    '''
    load all traces, compute spike times
    and create raster plots
    '''
    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    inhTypes = ('L1','L23Trans','L45Sym','L45Peak','L56Trans',\
                'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
    cellTypeColorMap = {'L2': 'dodgerblue', 'L34': 'blue', 'L4py': 'palegreen',\
                    'L4sp': 'green', 'L4ss': 'lime', 'L5st': 'yellow', 'L5tt': 'orange',\
                    'L6cc': 'indigo', 'L6ccinv': 'violet', 'L6ct': 'magenta', 'VPM': 'black',\
                    'INH': 'grey', 'EXC': 'red'}
    plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
    vmNames = []
#    scan_directory(folder, vmNames, vmSuffix)
    scan_directory2(folder, vmNames, vmSuffix, contains)
    synNames = []
#    scan_directory(folder, synNames, synapseSuffix)
    scan_directory2(folder, synNames, synapseSuffix, contains)
    
    print 'Computing reverse correlation of synaptic input from %d files' % len(synNames)
    
    # look at spikes in early phase (0-20ms)
    # and late phase (20-50ms) separately
    reverseTimes = {}
    reverseTimesProx = {}
    reverseTimesDistal = {}
    reverseTimesEarly = {}
    reverseTimesEarlyProx = {}
    reverseTimesEarlyDistal = {}
    reverseTimesLate = {}
    reverseTimesLateProx = {}
    reverseTimesLateDistal = {}
    for cellType in plotTypes:
        reverseTimes[cellType] = []
        reverseTimesProx[cellType] = []
        reverseTimesDistal[cellType] = []
        reverseTimesEarly[cellType] = []
        reverseTimesEarlyProx[cellType] = []
        reverseTimesEarlyDistal[cellType] = []
        reverseTimesLate[cellType] = []
        reverseTimesLateProx[cellType] = []
        reverseTimesLateDistal[cellType] = []
    
    tOffset = 100.0
    tStim = 245.0
    tStimWindow = 50.0
    earlyWindow = 25.0
    correlationWindow = 50.0
    binWidth = 1.0
    for n in range(len(vmNames)):
        fname = vmNames[n]
        pathName = fname[:fname.rfind('/')]
        print 'pathName = %s' % pathName
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        t = data[0]
        
        for i in range(1, len(data)):
#        for i in range(1,10):
            print 'Computing reverse correlation for spikes in trial %d of %d\r' % (i, len(data)-1),
            sys.stdout.flush()
            trialNr = i-1
            synTrialStr = 'simulation_run%04d_synapses.csv' % trialNr
            synTrialFile = ''
            for name in synNames:
                if synTrialStr in name and pathName in name:
                    synTrialFile = name
                    break
            if synTrialFile == '':
                errstr = 'Could not find synapse activation file for trial nr. %d' % trialNr
                raise RuntimeError(errstr)
            activeSyns = scp.read_complete_synapse_activation_file(synTrialFile)
            
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
                            if somaDist < 500.0:
                                synapseTimesProx['INH'][structure].extend(synTimes)
                                synapseTimesProx['INH']['Total'].extend(synTimes)
                            else:
                                synapseTimesDistal['INH'][structure].extend(synTimes)
                                synapseTimesDistal['INH']['Total'].extend(synTimes)
            v = data[i]
            trialSpikeTimes = sca.simple_spike_detection(t, v)
            for tSpike in trialSpikeTimes:
                #evoked
                if tSpike <= tStim or tSpike > tStim + tStimWindow:
                #ongoing
                #if tSpike < tOffset or tSpike > tStim:
                    continue
                for cellType in plotTypes:
                    reverseTimes[cellType].append([])
                    reverseTimesProx[cellType].append([])
                    reverseTimesDistal[cellType].append([])
                    if tStim <= tSpike < tStim + earlyWindow:
                        reverseTimesEarly[cellType].append([])
                        reverseTimesEarlyProx[cellType].append([])
                        reverseTimesEarlyDistal[cellType].append([])
                    else:
                        reverseTimesLate[cellType].append([])
                        reverseTimesLateProx[cellType].append([])
                        reverseTimesLateDistal[cellType].append([])
                    synTimes = synapseTimes[cellType]['Total']
                    synTimesProx = synapseTimesProx[cellType]['Total']
                    synTimesDistal = synapseTimesDistal[cellType]['Total']
                    for tSyn in synTimes:
#                    for tSyn in synTimesDistal:
                        if tSpike - correlationWindow <= tSyn < tSpike:
                            reverseTimes[cellType][-1].append(tSyn-tSpike)
                            if tStim <= tSpike < tStim + earlyWindow:
                                reverseTimesEarly[cellType][-1].append(tSyn-tSpike)
                            else:
                                reverseTimesLate[cellType][-1].append(tSyn-tSpike)
                    for tSyn in synTimesProx:
                        if tSpike - correlationWindow <= tSyn < tSpike:
                            reverseTimesProx[cellType][-1].append(tSyn-tSpike)
                            if tStim <= tSpike < tStim + earlyWindow:
                                reverseTimesEarlyProx[cellType][-1].append(tSyn-tSpike)
                            else:
                                reverseTimesLateProx[cellType][-1].append(tSyn-tSpike)
                    for tSyn in synTimesDistal:
                        if tSpike - correlationWindow <= tSyn < tSpike:
                            reverseTimesDistal[cellType][-1].append(tSyn-tSpike)
                            if tStim <= tSpike < tStim + earlyWindow:
                                reverseTimesEarlyDistal[cellType][-1].append(tSyn-tSpike)
                            else:
                                reverseTimesLateDistal[cellType][-1].append(tSyn-tSpike)
            
            #fig.add_subplot(2,1,1)
            #spikes = [i for time in trialSpikeTimes]
            #ax.append(plt.plot(trialSpikeTimes, spikes, 'k|'))
        print ''
    
    ax = []
    fig = plt.figure(1)
    for cellType in plotTypes:
        hist, bins = sca.PSTH_from_spike_times(reverseTimes[cellType], binWidth, -correlationWindow, 0.0)
        offset = 0.5*(bins[1] - bins[0])
        fig.add_subplot(3,1,1)
        ax.append(plt.plot(bins[:-1]+offset, hist, color=cellTypeColorMap[cellType], label=cellType, linewidth=2))
        hist, bins = sca.PSTH_from_spike_times(reverseTimesEarly[cellType], binWidth, -correlationWindow, 0.0)
        offset = 0.5*(bins[1] - bins[0])
        fig.add_subplot(3,1,2)
        ax.append(plt.plot(bins[:-1]+offset, hist, color=cellTypeColorMap[cellType], label=cellType, linewidth=2))
        hist, bins = sca.PSTH_from_spike_times(reverseTimesLate[cellType], binWidth, -correlationWindow, 0.0)
        offset = 0.5*(bins[1] - bins[0])
        fig.add_subplot(3,1,3)
        ax.append(plt.plot(bins[:-1]+offset, hist, color=cellTypeColorMap[cellType], label=cellType, linewidth=2))
    
    #fig.add_subplot(3,1,1)
    plt.xlabel('t [ms]')
    plt.ylabel('All spikes')
    plt.xlim([-correlationWindow, 0.0])
    #fig.add_subplot(3,1,2)
    #plt.xlabel('t [ms]')
    #plt.ylabel('Early spikes')
    #plt.xlim([-correlationWindow, 0.0])
    #fig.add_subplot(3,1,3)
    #plt.xlabel('t [ms]')
    #plt.ylabel('Late spikes')
    #plt.xlim([-correlationWindow, 0.0])
    
    fig = plt.figure(2)
    for cellType in plotTypes:
        hist, bins = sca.PSTH_from_spike_times(reverseTimesProx[cellType], binWidth, -correlationWindow, 0.0)
        offset = 0.5*(bins[1] - bins[0])
        fig.add_subplot(3,1,1)
        ax.append(plt.plot(bins[:-1]+offset, hist, color=cellTypeColorMap[cellType], label=cellType, linewidth=2))
        hist, bins = sca.PSTH_from_spike_times(reverseTimesEarlyProx[cellType], binWidth, -correlationWindow, 0.0)
        offset = 0.5*(bins[1] - bins[0])
        fig.add_subplot(3,1,2)
        ax.append(plt.plot(bins[:-1]+offset, hist, color=cellTypeColorMap[cellType], label=cellType, linewidth=2))
        hist, bins = sca.PSTH_from_spike_times(reverseTimesLateProx[cellType], binWidth, -correlationWindow, 0.0)
        offset = 0.5*(bins[1] - bins[0])
        fig.add_subplot(3,1,3)
        ax.append(plt.plot(bins[:-1]+offset, hist, color=cellTypeColorMap[cellType], label=cellType, linewidth=2))
    plt.xlabel('t [ms]')
    plt.ylabel('All spikes (proximal syns)')
    plt.xlim([-correlationWindow, 0.0])
    
    fig = plt.figure(3)
    for cellType in plotTypes:
        hist, bins = sca.PSTH_from_spike_times(reverseTimesDistal[cellType], binWidth, -correlationWindow, 0.0)
        offset = 0.5*(bins[1] - bins[0])
        fig.add_subplot(3,1,1)
        ax.append(plt.plot(bins[:-1]+offset, hist, color=cellTypeColorMap[cellType], label=cellType, linewidth=2))
        hist, bins = sca.PSTH_from_spike_times(reverseTimesEarlyDistal[cellType], binWidth, -correlationWindow, 0.0)
        offset = 0.5*(bins[1] - bins[0])
        fig.add_subplot(3,1,2)
        ax.append(plt.plot(bins[:-1]+offset, hist, color=cellTypeColorMap[cellType], label=cellType, linewidth=2))
        hist, bins = sca.PSTH_from_spike_times(reverseTimesLateDistal[cellType], binWidth, -correlationWindow, 0.0)
        offset = 0.5*(bins[1] - bins[0])
        fig.add_subplot(3,1,3)
        ax.append(plt.plot(bins[:-1]+offset, hist, color=cellTypeColorMap[cellType], label=cellType, linewidth=2))
    plt.xlabel('t [ms]')
    plt.ylabel('All spikes (distal syns)')
    plt.xlim([-correlationWindow, 0.0])
    
#        outName = fname[:-4]
#        outName += '_spike_raster_plot.pdf'
#        plt.savefig(outName)
    plt.show() 

def spike_time_correlation_synaptic_input(folder, contains, vmSuffix, synapseSuffix, averageSynapseName, outName):
    '''
    relative contributions of different
    presynaptic cell types as a function
    of spike time after stimulus
    '''
    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    inhTypes = ('L1','L23Trans','L45Sym','L45Peak','L56Trans',\
                'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
    cellTypeColorMap = {'L2': 'dodgerblue', 'L34': 'blue', 'L4py': 'palegreen',\
                    'L4sp': 'green', 'L4ss': 'lime', 'L5st': 'yellow', 'L5tt': 'orange',\
                    'L6cc': 'indigo', 'L6ccinv': 'violet', 'L6ct': 'magenta', 'VPM': 'black',\
                    'INH': 'grey', 'EXC': 'red'}
    plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
    
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
    
    print 'Computing reverse correlation of synaptic input from %d files' % len(synNames)
    
#    vmNames = vmData.keys()
#    synNames = synapseData.keys()
    
    # look at spikes in early phase (0-25ms)
    # and late phase (25-50ms) separately
    tOffset = 245.0
#    tStim = 245.0
    tStim = 253.0 # after VPM activation only
    tStimWindow = 50.0
#    earlyWindow = 25.0
    earlyWindow = 17.0 # after VPM activation only: 25-8
    binWidth = 1.0
    
    spikeTrialSyns = {}
    spikeTrialSynsProx = {}
    spikeTrialSynsDistal = {}
    spikeTrialSynsEarly = {}
    spikeTrialSynsEarlyProx = {}
    spikeTrialSynsEarlyDistal = {}
    # contains bins (i.e. int) of spike times at resolution binWidth (in ms)
    spikeTimesEarly = []
    spikeTrialSynsLate = {}
    spikeTrialSynsLateProx = {}
    spikeTrialSynsLateDistal = {}
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
                spikeTimesEarly.append(int(tSpikeMin/binWidth))
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
            activeSyns = synapseData[synTrialFile]
            
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
                            if somaDist < 500.0:
                                synapseTimesProx['INH'][structure].extend(synTimes)
                                synapseTimesProx['INH']['Total'].extend(synTimes)
                            else:
                                synapseTimesDistal['INH'][structure].extend(synTimes)
                                synapseTimesDistal['INH']['Total'].extend(synTimes)
            
            if not trialWithSpikes[n][trialNr]:
                nrNoSpikeTrials += 1
                continue
            
            elif trialWithSpikes[n][trialNr]:
                nrSpikeTrials += 1
                tSpikeReference = np.min(trialSpikeTimes[n][trialNr])
                for cellType in plotTypes:
                    spikeTrialSyns[cellType].append(0)
                    spikeTrialSynsProx[cellType].append(0)
                    spikeTrialSynsDistal[cellType].append(0)
                    spikeTrialSynsEarly[cellType].append(0)
                    spikeTrialSynsEarlyProx[cellType].append(0)
                    spikeTrialSynsEarlyDistal[cellType].append(0)
                    spikeTrialSynsLate[cellType].append(0)
                    spikeTrialSynsLateProx[cellType].append(0)
                    spikeTrialSynsLateDistal[cellType].append(0)
                    synTimes = synapseTimes[cellType]['Total']
                    synTimesProx = synapseTimesProx[cellType]['Total']
                    synTimesDistal = synapseTimesDistal[cellType]['Total']
                    for tSyn in synTimesProx:
                        # 8.0: VPM offset (tOffset vs. tStim)
                        if tOffset <= tSyn < tOffset + tSpikeReference + 8.0:
                            spikeTrialSynsEarlyProx[cellType][-1] += 1
                            earlyProxSyns += 1
                
            
        print ''
        print 'Nr of trials with spike: %d' % nrSpikeTrials
        print 'Nr of trials without spike: %d' % nrNoSpikeTrials
        print 'earlyProxSyns = %d' % earlyProxSyns
        print 'mean spike time = %.1fms' % np.mean([tSpike for trace in trialSpikeTimes[n] for tSpike in trace])
    
    cellTypeTrialContributions = {}
    for cellType in plotTypes:
        cellTypeTrialContributions[cellType] = [[] for i in range(int(earlyWindow))]
    
    for i in range(len(spikeTimesEarly)):
        evokedSyn = {}
#        totalNr = spikeTrialSynsEarlyProx['EXC'][i] - (spikeTimesEarly[i]+9.0)*averageSynsPerCellType['EXC']
#        if not totalNr:
#            continue
#        for cellType in excTypes:
#            # 9.0 = 1.0 + tStim-tOffset
#            evokedSyn_ = spikeTrialSynsEarlyProx[cellType][i] - (spikeTimesEarly[i]+9.0)*averageSynsPerCellType[cellType]
#            evokedSyn[cellType] = np.max([0.0, evokedSyn_])
#            cellTypeTrialContributions[cellType][spikeTimesEarly[i]].append(evokedSyn[cellType]/totalNr)
        norm = 0.0
        for cellType in excTypes:
            # 9.0 = 1.0 + tStim-tOffset
            evokedSyn_ = spikeTrialSynsEarlyProx[cellType][i] - (spikeTimesEarly[i]+9.0)*averageSynsPerCellType[cellType]
            evokedSyn[cellType] = np.max([0.0, evokedSyn_])
            norm += evokedSyn[cellType]
        for cellType in excTypes:
            cellTypeTrialContributions[cellType][spikeTimesEarly[i]].append(evokedSyn[cellType]/norm)
    
    if not outName.endswith('.csv'):
        outName += '.csv'
    
    with open(outName, 'w') as outFile:
        header = 'spike time\tnr. of spikes'
        for cellType in excTypes:
            header += '\t'
            header += cellType
            header += ' AVG\t'
            header += cellType
            header += ' STD'
        header += '\n'
        outFile.write(header)
        for i in range(int(earlyWindow)):
            line = str(i+0.5)
            line += '\t'
            line += str(len(cellTypeTrialContributions['VPM'][i]))
            for cellType in excTypes:
                avg = np.mean(cellTypeTrialContributions[cellType][i])
                std = np.std(cellTypeTrialContributions[cellType][i])
                line += '\t'
                line += str(avg)
                line += '\t'
                line += str(std)
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
#        elif fname.endswith(suffix) and fname.find(contains) != -1 and fname.find('C2_evoked') == -1 and fname.find('E2_evoked') == -1:
        elif fname.endswith(suffix) and fname.find(contains) != -1:
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 5:
        folder = sys.argv[1]
        contains = sys.argv[2]
        vmSuffix = sys.argv[3]
        synapseSuffix = sys.argv[4]
        reverse_correlation_synaptic_input(folder, contains, vmSuffix, synapseSuffix)
    if len(sys.argv) == 7:
        folder = sys.argv[1]
        contains = sys.argv[2]
        vmSuffix = sys.argv[3]
        synapseSuffix = sys.argv[4]
        averageSynapseName = sys.argv[5]
        outName = sys.argv[6]
        spike_time_correlation_synaptic_input(folder, contains, vmSuffix, synapseSuffix, averageSynapseName, outName)
    else:
        print 'Error! Number of arguments is %d; should be 4' % (len(sys.argv)-1)
    
    
    
    