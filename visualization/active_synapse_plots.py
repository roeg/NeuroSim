'''
Created on Jan 28, 2013

ongoing activity L2 neuron model

@author: robert
'''

import sys
import os, os.path
import glob
#import single_cell_analyzer as sca
#import single_cell_parser as scp
import numpy as np
import matplotlib.pyplot as plt

#def create_active_synapse_histogram(folder, suffix, contains, outName):
    #'''
    #load all traces, compute spike times
    #and create raster plots
    #'''
    #excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                #'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    #inhTypes = ('L1','L23Trans','L45Sym','L45Peak','L56Trans',\
                #'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
    #cellTypeColorMap = {'L2': 'dodgerblue', 'L34': 'blue', 'L4py': 'palegreen',\
                    #'L4sp': 'green', 'L4ss': 'lime', 'L5st': 'yellow', 'L5tt': 'orange',\
                    #'L6cc': 'indigo', 'L6ccinv': 'violet', 'L6ct': 'magenta', 'VPM': 'black',\
                    #'INH': 'grey', 'EXC': 'red'}
    #plotTypes = ('INH', 'L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                #'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    #columns = ('A1','A2','A3','A4','Alpha','B1','B2','B3','B4','Beta',\
               #'C1','C2','C3','C4','D1','D2','D3','D4','Delta','E1','E2','E3','E4','Gamma')
    
    #fnames = []
##    scan_directory(folder, fnames, suffix)
    #scan_directory2(folder, fnames, suffix, contains)
    #nrOfFiles = len(fnames)
    #print 'Creating active synapse plots from %d files' % nrOfFiles
    
    #synapseTimes = {}
    #synapseTimesProximal = {}
    #synapseTimesDistal = {}
    #for col in columns:
        #synapseTimes[col] = {}
        #synapseTimesProximal[col] = {}
        #synapseTimesDistal[col] = {}
        #for excType in excTypes:
            #synapseTimes[col][excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            #synapseTimesProximal[col][excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            #synapseTimesDistal[col][excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        #synapseTimes[col]['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        #synapseTimes[col]['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
        #synapseTimesProximal[col]['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        #synapseTimesProximal[col]['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
        #synapseTimesDistal[col]['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        #synapseTimesDistal[col]['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
    
    #for n in range(len(fnames)):
        #fname = fnames[n]
        #print 'Loading synapse activation times from file %s (file %d of %d)\r' % (fname,n+1,nrOfFiles) ,
        #sys.stdout.flush()
        #activeSyns = scp.read_complete_synapse_activation_file(fname)
        
        #for col in columns:
            #for excType in excTypes:
                #synapseTimes[col][excType]['ApicalDendrite'].append([])
                #synapseTimes[col][excType]['Dendrite'].append([])
                #synapseTimes[col][excType]['Total'].append([])
                #synapseTimesProximal[col][excType]['ApicalDendrite'].append([])
                #synapseTimesProximal[col][excType]['Dendrite'].append([])
                #synapseTimesProximal[col][excType]['Total'].append([])
                #synapseTimesDistal[col][excType]['ApicalDendrite'].append([])
                #synapseTimesDistal[col][excType]['Dendrite'].append([])
                #synapseTimesDistal[col][excType]['Total'].append([])
            #synapseTimes[col]['EXC']['ApicalDendrite'].append([])
            #synapseTimes[col]['EXC']['Dendrite'].append([])
            #synapseTimes[col]['EXC']['Total'].append([])
            #synapseTimesProximal[col]['EXC']['ApicalDendrite'].append([])
            #synapseTimesProximal[col]['EXC']['Dendrite'].append([])
            #synapseTimesProximal[col]['EXC']['Total'].append([])
            #synapseTimesDistal[col]['EXC']['ApicalDendrite'].append([])
            #synapseTimesDistal[col]['EXC']['Dendrite'].append([])
            #synapseTimesDistal[col]['EXC']['Total'].append([])
            #synapseTimes[col]['INH']['ApicalDendrite'].append([])
            #synapseTimes[col]['INH']['Dendrite'].append([])
            #synapseTimes[col]['INH']['Soma'].append([])
            #synapseTimes[col]['INH']['Total'].append([])
            #synapseTimesProximal[col]['INH']['ApicalDendrite'].append([])
            #synapseTimesProximal[col]['INH']['Dendrite'].append([])
            #synapseTimesProximal[col]['INH']['Soma'].append([])
            #synapseTimesProximal[col]['INH']['Total'].append([])
            #synapseTimesDistal[col]['INH']['ApicalDendrite'].append([])
            #synapseTimesDistal[col]['INH']['Dendrite'].append([])
            #synapseTimesDistal[col]['INH']['Soma'].append([])
            #synapseTimesDistal[col]['INH']['Total'].append([])
        
        #for synType in activeSyns.keys():
            #preCellType = synType.split('_')[0]
            #preCol = synType.split('_')[1]
            #for col in columns:
                #if col == preCol:
                    #for excType in excTypes:
                        #if excType == preCellType:
                            #for syn in activeSyns[synType]:
                                #somaDist = syn[1]
                                #structure = syn[4]
                                #synTimes = syn[5]
                                #synapseTimes[col][excType][structure][n].extend(synTimes)
                                #synapseTimes[col][excType]['Total'][n].extend(synTimes)
                                #synapseTimes[col]['EXC'][structure][n].extend(synTimes)
                                #synapseTimes[col]['EXC']['Total'][n].extend(synTimes)
                                #if somaDist < 500.0:
                                    #synapseTimesProximal[col][excType][structure][n].extend(synTimes)
                                    #synapseTimesProximal[col][excType]['Total'][n].extend(synTimes)
                                    #synapseTimesProximal[col]['EXC'][structure][n].extend(synTimes)
                                    #synapseTimesProximal[col]['EXC']['Total'][n].extend(synTimes)
                                #if somaDist >= 500.0:
                                    #synapseTimesDistal[col][excType][structure][n].extend(synTimes)
                                    #synapseTimesDistal[col][excType]['Total'][n].extend(synTimes)
                                    #synapseTimesDistal[col]['EXC'][structure][n].extend(synTimes)
                                    #synapseTimesDistal[col]['EXC']['Total'][n].extend(synTimes)
                    #for inhType in inhTypes:
                        #if inhType == preCellType:
                            #for syn in activeSyns[synType]:
                                #somaDist = syn[1]
                                #structure = syn[4]
                                #synTimes = syn[5]
                                #synapseTimes[col]['INH'][structure][n].extend(synTimes)
                                #synapseTimes[col]['INH']['Total'][n].extend(synTimes)
                                #if somaDist < 500.0:
                                    #synapseTimesProximal[col]['INH'][structure][n].extend(synTimes)
                                    #synapseTimesProximal[col]['INH']['Total'][n].extend(synTimes)
                                #if somaDist >= 500.0:
                                    #synapseTimesDistal[col]['INH'][structure][n].extend(synTimes)
                                    #synapseTimesDistal[col]['INH']['Total'][n].extend(synTimes)
    
    #print ''
    #tOffset = 100.0
    #tOngoing = 120.0
    #tOngoingWindow = 100.0
    #tStim = 245.0
    #tStimWindow = 50.0
##    tStim = 195.0
##    tStimWindow = 100.0
    #binWidth = 5.0
    #maxCount = 0
    #synapseHistogramsEvoked = {}
    #synapseHistogramsOngoing = {}
    #synapseHistogramsEvokedProximal = {}
    #synapseHistogramsOngoingProximal = {}
    #synapseHistogramsEvokedDistal = {}
    #synapseHistogramsOngoingDistal = {}
    #for col in columns:
        #synapseHistogramsEvoked[col] = {}
        #synapseHistogramsOngoing[col] = {}
        #synapseHistogramsEvokedProximal[col] = {}
        #synapseHistogramsOngoingProximal[col] = {}
        #synapseHistogramsEvokedDistal[col] = {}
        #synapseHistogramsOngoingDistal[col] = {}
        #for cellType in synapseTimes[col].keys():
            #synapseHistogramsEvoked[col][cellType] = {}
            #synapseHistogramsOngoing[col][cellType] = {}
            #synapseHistogramsEvokedProximal[col][cellType] = {}
            #synapseHistogramsOngoingProximal[col][cellType] = {}
            #synapseHistogramsEvokedDistal[col][cellType] = {}
            #synapseHistogramsOngoingDistal[col][cellType] = {}
            #for structure in synapseTimes[col][cellType].keys():
                #synTimes1 = synapseTimes[col][cellType][structure]
                #hist1, bins1 = sca.PSTH_from_spike_times(synTimes1, binWidth, tStim, tStim+tStimWindow)
                #synapseHistogramsEvoked[col][cellType][structure] = np.sum(hist1)
                #hist2, bins2 = sca.PSTH_from_spike_times(synTimes1, binWidth, tOngoing, tOngoing+tOngoingWindow)
                #synapseHistogramsOngoing[col][cellType][structure] = np.sum(hist2)
                
                #synTimes2 = synapseTimesProximal[col][cellType][structure]
                #hist3, bins3 = sca.PSTH_from_spike_times(synTimes2, binWidth, tStim, tStim+tStimWindow)
                #synapseHistogramsEvokedProximal[col][cellType][structure] = np.sum(hist3)
                #hist4, bins4 = sca.PSTH_from_spike_times(synTimes2, binWidth, tOngoing, tOngoing+tOngoingWindow)
                #synapseHistogramsOngoingProximal[col][cellType][structure] = np.sum(hist4)
                
                #synTimes3 = synapseTimesDistal[col][cellType][structure]
                #hist5, bins5 = sca.PSTH_from_spike_times(synTimes3, binWidth, tStim, tStim+tStimWindow)
                #synapseHistogramsEvokedDistal[col][cellType][structure] = np.sum(hist5)
                #hist6, bins6 = sca.PSTH_from_spike_times(synTimes3, binWidth, tOngoing, tOngoing+tOngoingWindow)
                #synapseHistogramsOngoingDistal[col][cellType][structure] = np.sum(hist6)
##                maxBin = np.max(hist)
##                if cellType != 'INH' and cellType != 'EXC' and maxBin > maxCount:
##                    maxCount = maxBin
    
    #ongoingOutName = outName + '_ongoing_syn.csv'
    #with open(ongoingOutName, 'w') as outputTable:
        #header = 'Column\tcell type\tnr. of active synapses (100ms)\n'
        #outputTable.write(header)
        #for cellType in plotTypes:
            #for col in columns:
                #line = col
                #line += '\t'
                #line += cellType
                #line += '\t'
                #line += str(synapseHistogramsOngoing[col][cellType]['Total'])
                #line += '\n'
                #outputTable.write(line)
    
    #writeHeader = True
    #evokedOutName = outName + '_total_evoked_syn.csv'
    #with open(evokedOutName, 'w') as outputTable:
        #header = 'Column\tcell type\tnr. of active synapses (50ms post-stimulus)\n'
        #outputTable.write(header)
        #for cellType in plotTypes:
            #for col in columns:
                #line = col
                #line += '\t'
                #line += cellType
                #line += '\t'
                #line += str(synapseHistogramsEvoked[col][cellType]['Total'])
                #line += '\n'
                #outputTable.write(line)
    
    #ongoingProximalOutName = outName + '_ongoing_proximal_syn.csv'
    #with open(ongoingProximalOutName, 'w') as outputTable:
        #header = 'Column\tcell type\tnr. of active synapses (100ms)\n'
        #outputTable.write(header)
        #for cellType in plotTypes:
            #for col in columns:
                #line = col
                #line += '\t'
                #line += cellType
                #line += '\t'
                #line += str(synapseHistogramsOngoingProximal[col][cellType]['Total'])
                #line += '\n'
                #outputTable.write(line)
    
    #writeHeader = True
    #evokedProximalOutName = outName + '_total_evoked_proximal_syn.csv'
    #with open(evokedProximalOutName, 'w') as outputTable:
        #header = 'Column\tcell type\tnr. of active synapses (50ms post-stimulus)\n'
        #outputTable.write(header)
        #for cellType in plotTypes:
            #for col in columns:
                #line = col
                #line += '\t'
                #line += cellType
                #line += '\t'
                #line += str(synapseHistogramsEvokedProximal[col][cellType]['Total'])
                #line += '\n'
                #outputTable.write(line)
    
    #ongoingDistalOutName = outName + '_ongoing_distal_syn.csv'
    #with open(ongoingDistalOutName, 'w') as outputTable:
        #header = 'Column\tcell type\tnr. of active synapses (100ms)\n'
        #outputTable.write(header)
        #for cellType in plotTypes:
            #for col in columns:
                #line = col
                #line += '\t'
                #line += cellType
                #line += '\t'
                #line += str(synapseHistogramsOngoingDistal[col][cellType]['Total'])
                #line += '\n'
                #outputTable.write(line)
    
    #writeHeader = True
    #evokedDistalOutName = outName + '_total_evoked_distal_syn.csv'
    #with open(evokedDistalOutName, 'w') as outputTable:
        #header = 'Column\tcell type\tnr. of active synapses (50ms post-stimulus)\n'
        #outputTable.write(header)
        #for cellType in plotTypes:
            #for col in columns:
                #line = col
                #line += '\t'
                #line += cellType
                #line += '\t'
                #line += str(synapseHistogramsEvokedDistal[col][cellType]['Total'])
                #line += '\n'
                #outputTable.write(line)
    
###    figCount = 0
###    for cellType in synapseHistograms.keys():
###        for structure in synapseHistograms[cellType].keys():
###            hist, bins = synapseHistograms[cellType][structure]
###            offset = 0.5*(bins[1] - bins[0])
###            fig = plt.figure(figCount)
###            plt.bar(bins[:-1], hist, color=cellTypeColorMap[cellType], width=binWidth)
###            if cellType != 'INH' and cellType != 'EXC':
###                plt.ylim([0, maxCount+1])
###            plt.xlabel('t [ms]')
###            yLabelStr = 'active synapses per %.1f ms' % binWidth
###            plt.ylabel(yLabelStr)
###            synName = outName
###            synName += '_' + cellType
###            synName += '_active_syns_' + structure
###            synName += '_%.1fms.pdf' % binWidth
###            plt.savefig(synName)
###            scp.write_PSTH(synName[:-4]+'.csv', hist, bins)
###            figCount += 1
##    
##    figCount = 1
##    lastBars = None
###    plt.figure(figCount)
##    tableOutName = outName + '.csv'
##    with open(tableOutName, 'w') as outputTable:
##        for cellType in plotTypes:
###            for structure in synapseHistograms[cellType].keys():
##            hist, bins = synapseHistograms[cellType]['Total']
##            
##            if lastBars is None:
##                header = 'bin'
##                for i in range(len(bins)-1):
##                    header += '\t'
##                    header += str(bins[i])
##                header += '\n'
##                outputTable.write(header)
##            line = cellType
##            for i in range(len(hist)):
##                line += '\t'
##                line += str(hist[i])
##            line += '\n'
##            outputTable.write(line)
##            
###            offset = 0.5*(bins[1] - bins[0])
###            fig = plt.figure(figCount)
###            plt.bar(bins[:-1], hist, color=cellTypeColorMap[cellType], width=binWidth, bottom=lastBars, label=cellType)
####            if cellType != 'INH' and cellType != 'EXC':
####                plt.ylim([0, maxCount+1])
####            synName = outName
####            synName += '_' + cellType
####            synName += '_active_syns_' + structure
####            synName += '_%.1fms.pdf' % binWidth
####            plt.savefig(synName)
####            scp.write_PSTH(synName[:-4]+'.csv', hist, bins)
####            figCount += 1
##            if lastBars is None:
##                lastBars = np.array(hist)
##            else:
##                lastBars += hist
###    plt.xlabel('t [ms]')
###    plt.xlim([tStim-tStimWindow, tStim+tStimWindow])
###    yLabelStr = 'active synapses per %.1f ms' % binWidth
###    plt.ylabel(yLabelStr)
###    plt.legend(loc='best')
####    plt.show()
###    plt.savefig(outName+'.pdf')

def create_active_synapse_histogram2(folder, suffix, contains, window, outName):
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
    plotTypes = ('INH', 'L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    columns = ('A1','A2','A3','A4','Alpha','B1','B2','B3','B4','Beta',\
               'C1','C2','C3','C4','D1','D2','D3','D4','Delta','E1','E2','E3','E4','Gamma')
    
    fnames = []
#    scan_directory(folder, fnames, suffix)
    scan_directory2(folder, fnames, suffix, contains)
    nrOfFiles = len(fnames)
    print 'Creating active synapse plots from %d files' % nrOfFiles
    
    synapseTimes = {}
    synapseTimesProximal = {}
    synapseTimesDistal = {}
    for col in columns:
        synapseTimes[col] = {}
        synapseTimesProximal[col] = {}
        synapseTimesDistal[col] = {}
        for excType in excTypes:
            synapseTimes[col][excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimesProximal[col][excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
            synapseTimesDistal[col][excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        synapseTimes[col]['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        synapseTimes[col]['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
        synapseTimesProximal[col]['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        synapseTimesProximal[col]['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
        synapseTimesDistal[col]['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        synapseTimesDistal[col]['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
    
    for n in range(len(fnames)):
        fname = fnames[n]
        print 'Loading synapse activation times from file %s (file %d of %d)\r' % (fname,n+1,nrOfFiles) ,
        sys.stdout.flush()
        activeSyns = read_complete_synapse_activation_file(fname)
        
        for col in columns:
            for excType in excTypes:
                synapseTimes[col][excType]['ApicalDendrite'].append([])
                synapseTimes[col][excType]['Dendrite'].append([])
                synapseTimes[col][excType]['Total'].append([])
                synapseTimesProximal[col][excType]['ApicalDendrite'].append([])
                synapseTimesProximal[col][excType]['Dendrite'].append([])
                synapseTimesProximal[col][excType]['Total'].append([])
                synapseTimesDistal[col][excType]['ApicalDendrite'].append([])
                synapseTimesDistal[col][excType]['Dendrite'].append([])
                synapseTimesDistal[col][excType]['Total'].append([])
            synapseTimes[col]['EXC']['ApicalDendrite'].append([])
            synapseTimes[col]['EXC']['Dendrite'].append([])
            synapseTimes[col]['EXC']['Total'].append([])
            synapseTimesProximal[col]['EXC']['ApicalDendrite'].append([])
            synapseTimesProximal[col]['EXC']['Dendrite'].append([])
            synapseTimesProximal[col]['EXC']['Total'].append([])
            synapseTimesDistal[col]['EXC']['ApicalDendrite'].append([])
            synapseTimesDistal[col]['EXC']['Dendrite'].append([])
            synapseTimesDistal[col]['EXC']['Total'].append([])
            synapseTimes[col]['INH']['ApicalDendrite'].append([])
            synapseTimes[col]['INH']['Dendrite'].append([])
            synapseTimes[col]['INH']['Soma'].append([])
            synapseTimes[col]['INH']['Total'].append([])
            synapseTimesProximal[col]['INH']['ApicalDendrite'].append([])
            synapseTimesProximal[col]['INH']['Dendrite'].append([])
            synapseTimesProximal[col]['INH']['Soma'].append([])
            synapseTimesProximal[col]['INH']['Total'].append([])
            synapseTimesDistal[col]['INH']['ApicalDendrite'].append([])
            synapseTimesDistal[col]['INH']['Dendrite'].append([])
            synapseTimesDistal[col]['INH']['Soma'].append([])
            synapseTimesDistal[col]['INH']['Total'].append([])
        
        for synType in activeSyns.keys():
            preCellType = synType.split('_')[0]
            preCol = synType.split('_')[1]
            for col in columns:
                if col == preCol:
                    for excType in excTypes:
                        if excType == preCellType:
                            for syn in activeSyns[synType]:
                                somaDist = syn[1]
                                structure = syn[4]
                                synTimes = syn[5]
                                synapseTimes[col][excType][structure][n].extend(synTimes)
                                synapseTimes[col][excType]['Total'][n].extend(synTimes)
                                synapseTimes[col]['EXC'][structure][n].extend(synTimes)
                                synapseTimes[col]['EXC']['Total'][n].extend(synTimes)
                                if somaDist < 500.0:
                                    synapseTimesProximal[col][excType][structure][n].extend(synTimes)
                                    synapseTimesProximal[col][excType]['Total'][n].extend(synTimes)
                                    synapseTimesProximal[col]['EXC'][structure][n].extend(synTimes)
                                    synapseTimesProximal[col]['EXC']['Total'][n].extend(synTimes)
                                if somaDist >= 500.0:
                                    synapseTimesDistal[col][excType][structure][n].extend(synTimes)
                                    synapseTimesDistal[col][excType]['Total'][n].extend(synTimes)
                                    synapseTimesDistal[col]['EXC'][structure][n].extend(synTimes)
                                    synapseTimesDistal[col]['EXC']['Total'][n].extend(synTimes)
                    for inhType in inhTypes:
                        if inhType == preCellType:
                            for syn in activeSyns[synType]:
                                somaDist = syn[1]
                                structure = syn[4]
                                synTimes = syn[5]
                                synapseTimes[col]['INH'][structure][n].extend(synTimes)
                                synapseTimes[col]['INH']['Total'][n].extend(synTimes)
                                if somaDist < 500.0:
                                    synapseTimesProximal[col]['INH'][structure][n].extend(synTimes)
                                    synapseTimesProximal[col]['INH']['Total'][n].extend(synTimes)
                                if somaDist >= 500.0:
                                    synapseTimesDistal[col]['INH'][structure][n].extend(synTimes)
                                    synapseTimesDistal[col]['INH']['Total'][n].extend(synTimes)
    
    print ''
    tOffset = 100.0
    tOngoing = 120.0
    tOngoingWindow = 100.0
    tStim = 245.0
    tStimWindow = window
#    tStim = 195.0
#    tStimWindow = 100.0
    binWidth = 5.0
    maxCount = 0
    synapseHistogramsEvoked = {}
    synapseHistogramsOngoing = {}
    synapseHistogramsEvokedProximal = {}
    synapseHistogramsOngoingProximal = {}
    synapseHistogramsEvokedDistal = {}
    synapseHistogramsOngoingDistal = {}
    for col in columns:
        synapseHistogramsEvoked[col] = {}
        synapseHistogramsOngoing[col] = {}
        synapseHistogramsEvokedProximal[col] = {}
        synapseHistogramsOngoingProximal[col] = {}
        synapseHistogramsEvokedDistal[col] = {}
        synapseHistogramsOngoingDistal[col] = {}
        for cellType in synapseTimes[col].keys():
            synapseHistogramsEvoked[col][cellType] = {}
            synapseHistogramsOngoing[col][cellType] = {}
            synapseHistogramsEvokedProximal[col][cellType] = {}
            synapseHistogramsOngoingProximal[col][cellType] = {}
            synapseHistogramsEvokedDistal[col][cellType] = {}
            synapseHistogramsOngoingDistal[col][cellType] = {}
            for structure in synapseTimes[col][cellType].keys():
                synTimes1 = synapseTimes[col][cellType][structure]
                synsPerTrialInWindow1 = []
                for times in synTimes1:
                    times_ = np.array(times)
                    synsPerTrialInWindow1.append(np.sum((tStim <= times_)*(times_ <= tStim+tStimWindow)))
                synapseHistogramsEvoked[col][cellType][structure] = np.mean(synsPerTrialInWindow1), np.std(synsPerTrialInWindow1)
                synsPerTrialInWindow2 = []
                for times in synTimes1:
                    times_ = np.array(times)
                    synsPerTrialInWindow2.append(np.sum((tOngoing <= times_)*(times_ <= tOngoing+tOngoingWindow)))
                synapseHistogramsOngoing[col][cellType][structure] = np.mean(synsPerTrialInWindow2), np.std(synsPerTrialInWindow2)
                
                synTimes2 = synapseTimesProximal[col][cellType][structure]
                synsPerTrialInWindow3 = []
                for times in synTimes2:
                    times_ = np.array(times)
                    synsPerTrialInWindow3.append(np.sum((tStim <= times_)*(times_ <= tStim+tStimWindow)))
                synapseHistogramsEvokedProximal[col][cellType][structure] = np.mean(synsPerTrialInWindow3), np.std(synsPerTrialInWindow3)
                synsPerTrialInWindow4 = []
                for times in synTimes2:
                    times_ = np.array(times)
                    synsPerTrialInWindow4.append(np.sum((tOngoing <= times_)*(times_ <= tOngoing+tOngoingWindow)))
                synapseHistogramsOngoingProximal[col][cellType][structure] = np.mean(synsPerTrialInWindow4), np.std(synsPerTrialInWindow4)
                
                synTimes3 = synapseTimesDistal[col][cellType][structure]
                synsPerTrialInWindow5 = []
                for times in synTimes3:
                    times_ = np.array(times)
                    synsPerTrialInWindow5.append(np.sum((tStim <= times_)*(times_ <= tStim+tStimWindow)))
                synapseHistogramsEvokedDistal[col][cellType][structure] = np.mean(synsPerTrialInWindow5), np.std(synsPerTrialInWindow5)
                synsPerTrialInWindow6 = []
                for times in synTimes3:
                    times_ = np.array(times)
                    synsPerTrialInWindow6.append(np.sum((tOngoing <= times_)*(times_ <= tOngoing+tOngoingWindow)))
                synapseHistogramsOngoingDistal[col][cellType][structure] = np.mean(synsPerTrialInWindow6), np.std(synsPerTrialInWindow6)
#                maxBin = np.max(hist)
#                if cellType != 'INH' and cellType != 'EXC' and maxBin > maxCount:
#                    maxCount = maxBin
    
    ongoingOutName = outName + '_ongoing_syn.csv'
    with open(ongoingOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (100ms)\tSTD\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsOngoing[col][cellType]['Total'][0])
                line += '\t'
                line += str(synapseHistogramsOngoing[col][cellType]['Total'][1])
                line += '\n'
                outputTable.write(line)
    
    writeHeader = True
    evokedOutName = outName + '_total_evoked_syn.csv'
    with open(evokedOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (' + str(window) + 'ms post-stimulus)\tSTD\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsEvoked[col][cellType]['Total'][0])
                line += '\t'
                line += str(synapseHistogramsEvoked[col][cellType]['Total'][1])
                line += '\n'
                outputTable.write(line)
    
    ongoingProximalOutName = outName + '_ongoing_proximal_syn.csv'
    with open(ongoingProximalOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (100ms)\tSTD\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsOngoingProximal[col][cellType]['Total'][0])
                line += '\t'
                line += str(synapseHistogramsOngoingProximal[col][cellType]['Total'][1])
                line += '\n'
                outputTable.write(line)
    
    writeHeader = True
    evokedProximalOutName = outName + '_total_evoked_proximal_syn.csv'
    with open(evokedProximalOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (' + str(window) + 'ms post-stimulus)\tSTD\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsEvokedProximal[col][cellType]['Total'][0])
                line += '\t'
                line += str(synapseHistogramsEvokedProximal[col][cellType]['Total'][1])
                line += '\n'
                outputTable.write(line)
    
    ongoingDistalOutName = outName + '_ongoing_distal_syn.csv'
    with open(ongoingDistalOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (100ms)\tSTD\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsOngoingDistal[col][cellType]['Total'][0])
                line += '\t'
                line += str(synapseHistogramsOngoingDistal[col][cellType]['Total'][1])
                line += '\n'
                outputTable.write(line)
    
    writeHeader = True
    evokedDistalOutName = outName + '_total_evoked_distal_syn.csv'
    with open(evokedDistalOutName, 'w') as outputTable:
        header = 'Column\tcell type\tnr. of active synapses (' + str(window) + 'ms post-stimulus)\tSTD\n'
        outputTable.write(header)
        for cellType in plotTypes:
            for col in columns:
                line = col
                line += '\t'
                line += cellType
                line += '\t'
                line += str(synapseHistogramsEvokedDistal[col][cellType]['Total'][0])
                line += '\t'
                line += str(synapseHistogramsEvokedDistal[col][cellType]['Total'][1])
                line += '\n'
                outputTable.write(line)
    
##    figCount = 0
##    for cellType in synapseHistograms.keys():
##        for structure in synapseHistograms[cellType].keys():
##            hist, bins = synapseHistograms[cellType][structure]
##            offset = 0.5*(bins[1] - bins[0])
##            fig = plt.figure(figCount)
##            plt.bar(bins[:-1], hist, color=cellTypeColorMap[cellType], width=binWidth)
##            if cellType != 'INH' and cellType != 'EXC':
##                plt.ylim([0, maxCount+1])
##            plt.xlabel('t [ms]')
##            yLabelStr = 'active synapses per %.1f ms' % binWidth
##            plt.ylabel(yLabelStr)
##            synName = outName
##            synName += '_' + cellType
##            synName += '_active_syns_' + structure
##            synName += '_%.1fms.pdf' % binWidth
##            plt.savefig(synName)
##            scp.write_PSTH(synName[:-4]+'.csv', hist, bins)
##            figCount += 1
#    
#    figCount = 1
#    lastBars = None
##    plt.figure(figCount)
#    tableOutName = outName + '.csv'
#    with open(tableOutName, 'w') as outputTable:
#        for cellType in plotTypes:
##            for structure in synapseHistograms[cellType].keys():
#            hist, bins = synapseHistograms[cellType]['Total']
#            
#            if lastBars is None:
#                header = 'bin'
#                for i in range(len(bins)-1):
#                    header += '\t'
#                    header += str(bins[i])
#                header += '\n'
#                outputTable.write(header)
#            line = cellType
#            for i in range(len(hist)):
#                line += '\t'
#                line += str(hist[i])
#            line += '\n'
#            outputTable.write(line)
#            
##            offset = 0.5*(bins[1] - bins[0])
##            fig = plt.figure(figCount)
##            plt.bar(bins[:-1], hist, color=cellTypeColorMap[cellType], width=binWidth, bottom=lastBars, label=cellType)
###            if cellType != 'INH' and cellType != 'EXC':
###                plt.ylim([0, maxCount+1])
###            synName = outName
###            synName += '_' + cellType
###            synName += '_active_syns_' + structure
###            synName += '_%.1fms.pdf' % binWidth
###            plt.savefig(synName)
###            scp.write_PSTH(synName[:-4]+'.csv', hist, bins)
###            figCount += 1
#            if lastBars is None:
#                lastBars = np.array(hist)
#            else:
#                lastBars += hist
##    plt.xlabel('t [ms]')
##    plt.xlim([tStim-tStimWindow, tStim+tStimWindow])
##    yLabelStr = 'active synapses per %.1f ms' % binWidth
##    plt.ylabel(yLabelStr)
##    plt.legend(loc='best')
###    plt.show()
##    plt.savefig(outName+'.pdf')

def create_active_synapse_histogram_spike_no_spike(folder, suffix, vmSuffix, contains, outName):
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
#    scan_directory(folder, synNames, suffix)
    scan_directory2(folder, synNames, suffix, contains)
    nrOfFiles = len(synNames)
    print 'Creating active synapse plots from %d files' % nrOfFiles
    synData = {}
    for fname in synNames:
        activeSyns = scp.read_complete_synapse_activation_file(fname)
        synData[fname] = activeSyns
    
    synapseTimes = {}
    spikeTrialSyns = {}
    noSpikeTrialSyns = {}
    for excType in excTypes:
        synapseTimes[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        spikeTrialSyns[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
        noSpikeTrialSyns[excType] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
    synapseTimes['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
    synapseTimes['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
    spikeTrialSyns['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
    spikeTrialSyns['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
    noSpikeTrialSyns['EXC'] = {'ApicalDendrite': [], 'Dendrite': [], 'Total': []}
    noSpikeTrialSyns['INH'] = {'ApicalDendrite': [], 'Dendrite': [], 'Soma': [], 'Total': []}
    
    data = []
    trialSpikeTimes = [[] for j in range(len(vmNames))]
    trialWithSpikes = {}
    
    tStim = 253.0
    earlyWindow = 17.0
    
    for n in range(len(vmNames)):
        fname = vmNames[n]
#        pathName = fname[:fname.rfind('/')]
#        print 'pathName = %s' % pathName
        print 'Analyzing traces in file %s' % fname
        data_ = np.loadtxt(fname, skiprows=1, unpack=True)
#        data_ = vmData[fname]
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
#            if trialWithSpikes_:
#                tSpikeMin = np.min(trialSpikeTimes[n][-1])
#                spikeTimesEarly.append(tSpikeMin)
            trialWithSpikes[n].append(trialWithSpikes_)
    
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
            
            for excType in excTypes:
                synapseTimes[excType]['ApicalDendrite'].append([])
                synapseTimes[excType]['Dendrite'].append([])
                synapseTimes[excType]['Total'].append([])
                if trialWithSpikes[n][trialNr]:
                    spikeTrialSyns[excType]['ApicalDendrite'].append([])
                    spikeTrialSyns[excType]['Dendrite'].append([])
                    spikeTrialSyns[excType]['Total'].append([])
                else:
                    noSpikeTrialSyns[excType]['ApicalDendrite'].append([])
                    noSpikeTrialSyns[excType]['Dendrite'].append([])
                    noSpikeTrialSyns[excType]['Total'].append([])
            synapseTimes['EXC']['ApicalDendrite'].append([])
            synapseTimes['EXC']['Dendrite'].append([])
            synapseTimes['EXC']['Total'].append([])
            synapseTimes['INH']['ApicalDendrite'].append([])
            synapseTimes['INH']['Dendrite'].append([])
            synapseTimes['INH']['Soma'].append([])
            synapseTimes['INH']['Total'].append([])
            if trialWithSpikes[n][trialNr]:
                spikeTrialSyns['EXC']['ApicalDendrite'].append([])
                spikeTrialSyns['EXC']['Dendrite'].append([])
                spikeTrialSyns['EXC']['Total'].append([])
                spikeTrialSyns['INH']['ApicalDendrite'].append([])
                spikeTrialSyns['INH']['Dendrite'].append([])
                spikeTrialSyns['INH']['Soma'].append([])
                spikeTrialSyns['INH']['Total'].append([])
            else:
                noSpikeTrialSyns['EXC']['ApicalDendrite'].append([])
                noSpikeTrialSyns['EXC']['Dendrite'].append([])
                noSpikeTrialSyns['EXC']['Total'].append([])
                noSpikeTrialSyns['INH']['ApicalDendrite'].append([])
                noSpikeTrialSyns['INH']['Dendrite'].append([])
                noSpikeTrialSyns['INH']['Soma'].append([])
                noSpikeTrialSyns['INH']['Total'].append([])
            
            for synType in activeSyns.keys():
                preCellType = synType.split('_')[0]
                for excType in excTypes:
                    if excType == preCellType:
                        for syn in activeSyns[synType]:
                            somaDist = syn[1]
                            structure = syn[4]
                            synTimes = syn[5]
                            if somaDist < 500.0:
                                synapseTimes[excType][structure][n].extend(synTimes)
                                synapseTimes[excType]['Total'][n].extend(synTimes)
                                synapseTimes['EXC'][structure][n].extend(synTimes)
                                synapseTimes['EXC']['Total'][n].extend(synTimes)
                                if trialWithSpikes[n][trialNr]:
                                    spikeTrialSyns[excType][structure][n].extend(synTimes)
                                    spikeTrialSyns[excType]['Total'][n].extend(synTimes)
                                    spikeTrialSyns['EXC'][structure][n].extend(synTimes)
                                    spikeTrialSyns['EXC']['Total'][n].extend(synTimes)
                                else:
                                    noSpikeTrialSyns[excType][structure][n].extend(synTimes)
                                    noSpikeTrialSyns[excType]['Total'][n].extend(synTimes)
                                    noSpikeTrialSyns['EXC'][structure][n].extend(synTimes)
                                    noSpikeTrialSyns['EXC']['Total'][n].extend(synTimes)
                for inhType in inhTypes:
                    if inhType == preCellType:
                        for syn in activeSyns[synType]:
                            somaDist = syn[1]
                            structure = syn[4]
                            synTimes = syn[5]
                            if somaDist < 500.0:
                                synapseTimes['INH'][structure][n].extend(synTimes)
                                synapseTimes['INH']['Total'][n].extend(synTimes)
                                if trialWithSpikes[n][trialNr]:
                                    spikeTrialSyns['INH'][structure][n].extend(synTimes)
                                    spikeTrialSyns['INH']['Total'][n].extend(synTimes)
                                else:
                                    noSpikeTrialSyns['INH'][structure][n].extend(synTimes)
                                    noSpikeTrialSyns['INH']['Total'][n].extend(synTimes)
        print ''
    
    tOffset = 100.0
#    tStim = 245.0
#    tStimWindow = 50.0
    tPlotBegin = 195.0
    tPlotBeginWindow = 100.0
#    binWidth = 5.0
    binWidth = 1.0
    maxCount = 0
    synapseHistograms = {}
    spikeTrialHistograms = {}
    noSpikeTrialHistograms = {}
    for cellType in synapseTimes.keys():
        synapseHistograms[cellType] = {}
        spikeTrialHistograms[cellType] = {}
        noSpikeTrialHistograms[cellType] = {}
        for structure in synapseTimes[cellType].keys():
            synTimes = synapseTimes[cellType][structure]
            hist, bins = sca.PSTH_from_spike_times(synTimes, binWidth, tOffset, tPlotBegin+tPlotBeginWindow)
            synapseHistograms[cellType][structure] = hist, bins
            spikeTrialSynTimes = spikeTrialSyns[cellType][structure]
            hist2, bins2 = sca.PSTH_from_spike_times(spikeTrialSynTimes, binWidth, tOffset, tPlotBegin+tPlotBeginWindow)
            spikeTrialHistograms[cellType][structure] = hist2, bins2
            noSpikeTrialSynTimes = noSpikeTrialSyns[cellType][structure]
            hist3, bins3 = sca.PSTH_from_spike_times(noSpikeTrialSynTimes, binWidth, tOffset, tPlotBegin+tPlotBeginWindow)
            noSpikeTrialHistograms[cellType][structure] = hist3, bins3
#            maxBin = np.max(hist)
#            if cellType != 'INH' and cellType != 'EXC' and maxBin > maxCount:
#                maxCount = maxBin
    
#    figCount = 0
#    for cellType in synapseHistograms.keys():
#        for structure in synapseHistograms[cellType].keys():
#            hist, bins = synapseHistograms[cellType][structure]
#            offset = 0.5*(bins[1] - bins[0])
#            fig = plt.figure(figCount)
#            plt.bar(bins[:-1], hist, color=cellTypeColorMap[cellType], width=binWidth)
#            if cellType != 'INH' and cellType != 'EXC':
#                plt.ylim([0, maxCount+1])
#            plt.xlabel('t [ms]')
#            yLabelStr = 'active synapses per %.1f ms' % binWidth
#            plt.ylabel(yLabelStr)
#            synName = outName
#            synName += '_' + cellType
#            synName += '_active_syns_' + structure
#            synName += '_%.1fms.pdf' % binWidth
#            plt.savefig(synName)
#            scp.write_PSTH(synName[:-4]+'.csv', hist, bins)
#            figCount += 1
    
    figCount = 1
#    lastBars = None
#    plt.figure(figCount)
    tableOutName = outName + '_all_trials.csv'
    with open(tableOutName, 'w') as outputTable:
        hist, bins = synapseHistograms['EXC']['Total']
        header = 'Bin start\tbin end\tbin center'
        for cellType in plotTypes:
            header += '\t'
            header += cellType
        header += '\n'
        outputTable.write(header)
        for i in range(len(bins)-1):
            line = str(bins[i])
            line += '\t'
            line += str(bins[i+1])
            line += '\t'
            line += str(0.5*(bins[i]+bins[i+1]))
            for cellType in plotTypes:
#                for structure in synapseHistograms[cellType].keys():
                line += '\t'
                line += str(synapseHistograms[cellType]['Total'][0][i])
            line += '\n'
            outputTable.write(line)
    
    tableOutName2 = outName + '_spike_trials.csv'
    with open(tableOutName2, 'w') as outputTable:
        hist, bins = spikeTrialHistograms['EXC']['Total']
        header = 'Bin start\tbin end\tbin center'
        for cellType in plotTypes:
            header += '\t'
            header += cellType
        header += '\n'
        outputTable.write(header)
        for i in range(len(bins)-1):
            line = str(bins[i])
            line += '\t'
            line += str(bins[i+1])
            line += '\t'
            line += str(0.5*(bins[i]+bins[i+1]))
            for cellType in plotTypes:
#                for structure in spikeTrialHistograms[cellType].keys():
                line += '\t'
                line += str(spikeTrialHistograms[cellType]['Total'][0][i])
            line += '\n'
            outputTable.write(line)
    
    tableOutName3 = outName + '_no_spike_trials.csv'
    with open(tableOutName3, 'w') as outputTable:
        hist, bins = noSpikeTrialHistograms['EXC']['Total']
        header = 'Bin start\tbin end\tbin center'
        for cellType in plotTypes:
            header += '\t'
            header += cellType
        header += '\n'
        outputTable.write(header)
        for i in range(len(bins)-1):
            line = str(bins[i])
            line += '\t'
            line += str(bins[i+1])
            line += '\t'
            line += str(0.5*(bins[i]+bins[i+1]))
            for cellType in plotTypes:
#                for structure in noSpikeTrialHistograms[cellType].keys():
                line += '\t'
                line += str(noSpikeTrialHistograms[cellType]['Total'][0][i])
            line += '\n'
            outputTable.write(line)
            
#            offset = 0.5*(bins[1] - bins[0])
#            fig = plt.figure(figCount)
#            plt.bar(bins[:-1], hist, color=cellTypeColorMap[cellType], width=binWidth, bottom=lastBars, label=cellType)
##            if cellType != 'INH' and cellType != 'EXC':
##                plt.ylim([0, maxCount+1])
##            synName = outName
##            synName += '_' + cellType
##            synName += '_active_syns_' + structure
##            synName += '_%.1fms.pdf' % binWidth
##            plt.savefig(synName)
##            scp.write_PSTH(synName[:-4]+'.csv', hist, bins)
##            figCount += 1
#            if lastBars is None:
#                lastBars = np.array(hist)
#            else:
#                lastBars += hist
#    plt.xlabel('t [ms]')
#    plt.xlim([tStim-tStimWindow, tStim+tStimWindow])
#    yLabelStr = 'active synapses per %.1f ms' % binWidth
#    plt.ylabel(yLabelStr)
#    plt.legend(loc='best')
##    plt.show()
#    plt.savefig(outName+'.pdf')

def create_active_synapse_summary_files(folder, outName):
    '''
    create summary tables and plots for 3x3 active synapese
    for average population and individual locations,
    sorted by total/per cell type
    '''
    # parse directory, load csv files, subtract ongoing activity, summarize in 5x10ms bins
    # add attributes to each PSTH:
    # -postsynaptic location/avg population
    # -subcellular structure
    # -presynaptic cell type
    # -deflected whisker
    
    locations = ('average','B1border','B2border','B3border','C1border',\
                 'C2center','C3border','D1border','D2border','D3border')
    structures = ('ApicalDendrite','Dendrite','Soma','Total')
#    activeStructures = ('ApicalDendrite','Dendrite','Total')
    activeStructures = ('ApicalDendrite','Dendrite','Soma','Total')
    preCellTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                   'L6cc', 'L6ccinv', 'L6ct', 'VPM','EXC','INH')
#    activeTypes = ('EXC', 'L34', 'L4sp', 'L4ss', 'L5tt', 'L6cc', 'VPM')
    activeTypes = ('INH',)
    whiskers = ('B1','B2','B3','C1','C2','C3','D1','D2','D3')
    whiskerLUT = {'B1': (2,2), 'B2': (2,1), 'B3': (2,0),\
                  'C1': (1,2), 'C2': (1,1), 'C3': (1,0),\
                  'D1': (0,2), 'D2': (0,1), 'D3': (0,0)}
    rowLUT = {0: 'D row', 1: 'C row', 2: 'B row'}
    
    class ActiveSynapseContainer(object):
        PSTH = None
        location = None
        structure = None
        preCellType = None
        whisker = None
        
        def __init__(self, PSTH):
            self.PSTH = PSTH
        
        def set_properties(self, filename):
            splitName = filename.split('/')
            csvName = splitName[-1]
            locationName = splitName[-2]
            whiskerName = splitName[-3]
            csvSplitName = csvName.split('_')
            structureName = csvSplitName[-2]
            preCellTypeName = csvSplitName[-5]
            for location in locations:
                if location == locationName:
                    self.location = location
                    break
            for structure in structures:
                if structure == structureName:
                    self.structure = structure
                    break
            for preCellType in preCellTypes:
                if preCellType == preCellTypeName:
                    self.preCellType = preCellType
                    break
            for whisker in whiskers:
                if whisker in whiskerName:
                    self.whisker = whisker
                    break
            if self.location is None:
                errstr = 'Could not determine location from filename %s' % filename
                raise RuntimeError(errstr)
            if self.structure is None:
                errstr = 'Could not determine structure from filename %s' % filename
                raise RuntimeError(errstr)
            if self.preCellType is None:
                errstr = 'Could not determine preCellType from filename %s' % filename
                raise RuntimeError(errstr)
            if self.whisker is None:
                errstr = 'Could not determine whisker from filename %s' % filename
                raise RuntimeError(errstr)
    
    allFileNames = []
    allActiveSynHist = []
    scan_directory(folder, allFileNames, '.csv')
    for fileName in allFileNames:
        tBegin, tEnd, synCount = np.loadtxt(fileName, skiprows=1, unpack=True)
        ongoing = np.mean(synCount[:20])
        evokedSyn = synCount - ongoing
        synHist = []
        for i in range(5):
            synHist.append(evokedSyn[20+2*i]+evokedSyn[20+2*i+1])
        synContainer = ActiveSynapseContainer(synHist)
        synContainer.set_properties(fileName)
        allActiveSynHist.append(synContainer)
    
    # summarize 3x3 plots/tables, i.e. sort by:
    # -location
    # -structure
    # -total/active cell types
    # (for now)
    summaryData = {}
    for location in locations:
        summaryData[location] = {}
        for structure in activeStructures:
            summaryData[location][structure] = {}
            for activeType in activeTypes:
                summaryData[location][structure][activeType] = [np.zeros(shape=(3,3)) for i in range(5)]
    
    for synHist in allActiveSynHist:
        location = synHist.location
        structure = synHist.structure
        cellType = synHist.preCellType
        index = whiskerLUT[synHist.whisker]
        if cellType not in activeTypes or structure not in activeStructures:
            continue
        for i in range(len(synHist.PSTH)):
            summaryData[location][structure][cellType][i][index] += synHist.PSTH[i]
    
    for location in summaryData:
        for structure in summaryData[location]:
            locName = outName
            locName += '_' + location + '_' + structure + '.csv'
            with open(locName, 'w') as outFile:
                for i in range(5):
                    line = '%d-%dms\n' % (i*10, (i+1)*10)
                    outFile.write(line)
                    line = ''
                    for cellType in activeTypes:
                        line += cellType
                        line += '\t3\t2\t1\t\t'
                    line += '\n'
                    outFile.write(line)
                    for j in range(3):
                        line = ''
                        for cellType in activeTypes:
                            line += rowLUT[j]
                            line += '\t'
                            for k in range(3):
                                line += str(summaryData[location][structure][cellType][i][(j,k)])
                                line += '\t'
                            line += '\t'
                        line += '\n'
                        outFile.write(line)
                    line = '\n'
                    outFile.write(line)

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
        #elif fname.endswith(suffix) and fname.find(contains) != -1 and fname.find('C2_evoked') == -1 and fname.find('E2_evoked') == -1:
        elif fname.endswith(suffix) and fname.find(contains) != -1:
            fnames.append(fname)
        else:
            continue

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

if __name__ == '__main__':
    if len(sys.argv) == 5:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        contains = sys.argv[3]
        outName = sys.argv[4]
        #create_active_synapse_histogram(folder, suffix, contains, outName)
    elif len(sys.argv) == 6:
        #folder = sys.argv[1]
        #suffix = sys.argv[2]
        #vmSuffix = sys.argv[3]
        #contains = sys.argv[4]
        #outName = sys.argv[5]
        #create_active_synapse_histogram_spike_no_spike(folder, suffix, vmSuffix, contains, outName)
        folder = sys.argv[1]
        suffix = sys.argv[2]
        contains = sys.argv[3]
        window = float(sys.argv[4])
        outName = sys.argv[5]
        create_active_synapse_histogram2(folder, suffix, contains, window, outName)
    elif len(sys.argv) == 3:
        folder = sys.argv[1]
        outName = sys.argv[2]
        create_active_synapse_summary_files(folder, outName)
    else:
        print 'Error! Number of arguments is %d; should be 2, 4 or 5' % (len(sys.argv)-1)
    
    
    
    
