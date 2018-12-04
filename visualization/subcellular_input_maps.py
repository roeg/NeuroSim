'''
@author: robert
'''

import sys
import os, os.path
import glob
import numpy as np
import matplotlib.pyplot as plt
import single_cell_parser as scp
import single_cell_analyzer as sca

def subcellular_input_map(cellName, networkParamName, synapseFilename, outputFilename, cellType=None, tBegin=None, tEnd=None):
    '''
    visualization of subcellular input distributions:
    -per cell type and time window possible
    -output: per cell type landmark files, list of distances
    '''
    neuronParameters = scp.build_parameters(cellName)
    networkParameters = scp.build_parameters(networkParamName)
    scp.load_NMODL_parameters(neuronParameters)
    scp.load_NMODL_parameters(networkParameters)
    cellParam = neuronParameters.neuron
    
    print 'Loading cell morphology...'
    parser = scp.CellParser(cellParam.filename)
    parser.spatialgraph_to_cell()
    cell = parser.cell
    
    synLocations = []
    synDistances = []
    if tBegin is None and tEnd is None:
        evokedNW = scp.NetworkMapper(cell, networkParameters.network)
        evokedNW._assign_anatomical_synapses()
        synTypes = cell.synapses.keys()
        for synType in synTypes:
            preCellType, preColumn = synType.split('_')
            if cellType is not None and preCellType != cellType:
                continue
            print 'Computing synapse distances for cell type %s' % synType
            synDistances.extend(sca.compute_syn_distances(cell, synType))
            for syn in cell.synapses[synType]:
                secID = syn.secID
                ptID = syn.ptID
                x = syn.x
                synLoc3D = cell.sections[secID].pts[ptID]
                synLocations.append(synLoc3D)
    else:
        synInfo = scp.read_synapse_activation_file(synapseFilename)
        synTypes = synInfo.keys()
        for synType in synTypes:
            preCellType, preColumn = synType.split('_')
            if cellType is not None and preCellType != cellType:
                continue
            print 'Loading synapse distances and activation times for cell type %s' % synType
            for syn in synInfo[synType]:
                synID, secID, ptID, synTimes, somaDist = syn
                activeSyn = False
                for t in synTimes:
                    if tBegin <= t <= tEnd:
                        activeSyn = True
                        break
                if not activeSyn:
                    continue
                synLoc3D = cell.sections[secID].pts[ptID]
                synLocations.append(synLoc3D)
                synDistances.append(somaDist)
    
    scp.write_landmark_file(outputFilename, synLocations)
    write_synapse_distances(outputFilename, synDistances)

def subcellular_distribution_active_synapses(folder, contains, synapseSuffix, outName):
    '''
    load all synapse times from all trials
    and compute average subcellular histograms
    of active synapses at 50 micron resolution
    '''
    synNames = []
#    scan_directory(folder, synNames, synapseSuffix)
    scan_directory2(folder, synNames, synapseSuffix, contains)
    
    print 'Loading %d synapse activation files...' % len(synNames)
    synapseData = {}
    for synTrialFile in synNames:
        synapseData[synTrialFile] = scp.read_complete_synapse_activation_file(synTrialFile)
    
    # all separately
    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    inhTypes = ('L1','L23Trans','L45Sym','L45Peak','L56Trans',\
                'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
    cellTypeColorMap = {'L2': 'dodgerblue', 'L34': 'blue', 'L4py': 'palegreen',\
                    'L4sp': 'green', 'L4ss': 'lime', 'L5st': 'yellow', 'L5tt': 'orange',\
                    'L6cc': 'indigo', 'L6ccinv': 'violet', 'L6ct': 'magenta', 'VPM': 'black',\
                    'INH': 'grey', 'EXC': 'red'}
    # all separately
    #plotTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                #'L6cc', 'L6ccinv', 'L6ct', 'VPM', 'EXC', 'INH')
    # only EXC/INH
    plotTypes = ('EXC', 'INH')
    
    synFileNames = synapseData.keys()
    
    tOffset = 245.0
    tStim = 253.0 # after VPM activation only
    offsetWindow50 = 42.0
    
    bins = [i*50 for i in range(33)] # 0-1600 microns
    synapseHistograms = {}
    for cellType in plotTypes:
        synapseHistograms[cellType] = []
    for n in range(len(synFileNames)):
        print 'Computing subcellular histogram of synapses in trial %d of %d\r' % (n+1, len(synFileNames)),
        synTrialFile = synFileNames[n]
        activeSyns = synapseData[synTrialFile]
        synapseDistances = {}
        for excType in excTypes:
            synapseDistances[excType] = []
        synapseDistances['EXC'] = []
        synapseDistances['INH'] = []
        
        for synType in activeSyns.keys():
            preCellType = synType.split('_')[0]
            for excType in excTypes:
                if excType == preCellType:
                    for syn in activeSyns[synType]:
                        somaDist = syn[1]
                        synTimes = syn[5]
                        for tSyn in synTimes:
                            # spike trial 0-50ms all synapses mode
                            if tStim <= tSyn < tStim + offsetWindow50:
                                synapseDistances[excType].append(somaDist)
                                synapseDistances['EXC'].append(somaDist)
                                break
            for inhType in inhTypes:
                if inhType == preCellType:
                    for syn in activeSyns[synType]:
                        somaDist = syn[1]
                        synTimes = syn[5]
                        for tSyn in synTimes:
                            # spike trial 0-50ms all synapses mode
                            if tStim <= tSyn < tStim + offsetWindow50:
                                synapseDistances['INH'].append(somaDist)
                                break
        
        for cellType in plotTypes:
            hist, _ = np.histogram(synapseDistances[cellType], bins=bins)
            synapseHistograms[cellType].append(hist)
    
    print ''
    
    histAvg, histStd = {}, {}
    for cellType in plotTypes:
        histAvg[cellType] = np.mean(synapseHistograms[cellType], axis=0)
        histStd[cellType] = np.std(synapseHistograms[cellType], axis=0)
    
    histName = ''
    if not outName.endswith('.csv'):
        histName = outName + '_subcellular_synapse_distribution.csv'
    else:
        histName = outName[:-4] + '_subcellular_synapse_distribution.csv'
    
    with open(histName, 'w') as outFile:
        header = '#bin begin (microns)\tbin end'
        for cellType in plotTypes:
            header += '\t'
            header += cellType
            header += ' AVG'
            header += '\t'
            header += cellType
            header += ' STD'
        header += '\n'
        outFile.write(header)
        for i in range(len(bins)-1):
            line = str(bins[i])
            line += '\t'
            line += str(bins[i+1])
            for cellType in plotTypes:
                line += '\t'
                line += str(histAvg[cellType][i])
                line += '\t'
                line += str(histStd[cellType][i])
            line += '\n'
            outFile.write(line)

def write_synapse_distances(fname, distances):
    if not fname.endswith('distances.csv'):
        fname += '_distances.csv'
    with open(fname, 'w') as outputFile:
        header = '# distance to soma [micron]\n\n'
        outputFile.write(header)
        for dist in distances:
            line = '%.3f\n' % dist
            outputFile.write(line)

def scan_directory2(path, fnames, suffix, contains):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory2(fname, fnames, suffix, contains)
#        elif fname.endswith(suffix) and fname.find(contains) != -1 and fname.find('C2_evoked') == -1 and fname.find('E2_evoked') == -1\
#            and fname.find('L6cc_timing') == -1 and fname.find('dynamic_syns') == -1 and fname.find('L5tt_evoked_inact') == -1:
#        elif fname.endswith(suffix) and fname.find(contains) != -1 and fname.find('L6cc_timing') == -1 and fname.find('dynamic_syns') == -1 and fname.find('L5tt_evoked_inact') == -1:
        elif fname.endswith(suffix) and fname.find(contains) != -1 and fname.find('C2_evoked') == -1 and fname.find('E2_evoked') == -1:
#        elif fname.endswith(suffix) and fname.find(contains) != -1 and fname.find('E2_evoked') == -1:
#        elif fname.endswith(suffix) and fname.find(contains) != -1:
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 4:
        folder = sys.argv[1]
        contains = sys.argv[2]
        outName = sys.argv[3]
        synapseSuffix = 'synapses.csv'
        subcellular_distribution_active_synapses(folder, contains, synapseSuffix, outName)
    elif len(sys.argv) >= 5:
        cellName = sys.argv[1]
        networkParamName = sys.argv[2]
        synapseFilename = sys.argv[3]
        outputFilename = sys.argv[4]
        cellType = None
        tBegin = None
        tEnd = None
        if len(sys.argv) == 6:
            cellType = sys.argv[5]
        elif len(sys.argv) == 7:
            tBegin = float(sys.argv[5])
            tEnd = float(sys.argv[6])
        elif len(sys.argv) == 8:
            cellType = sys.argv[5]
            tBegin = float(sys.argv[6])
            tEnd = float(sys.argv[7])
        subcellular_input_map(cellName, networkParamName, synapseFilename, outputFilename, cellType, tBegin, tEnd)
    else:
        print 'Error! Number of arguments is %d; should be at least 4' % (len(sys.argv)-1)
    
    
    
    
