'''
@author: robert
'''

import sys
import numpy as np
import single_cell_parser as scp
    
def visualize_active_soma_locations(somaLocationName, connectionName, synapseActivationName, outputFilename, cellType=None, tBegin=None, tEnd=None, structure=None):
    '''
    write landmark file for all active
    somata within specified time frame and/or
    for specific cell types
    '''
    
    '''
    BUG: cannot look at active presynaptic neurons -
    have to look at active synapses and determine
    presynaptic neurons (prel < 1)!!!
    '''
    
    somaLocations = load_cell_location_realization(somaLocationName)
    connectedNeuronsPerSynapse = load_connected_neuron_per_synapse(connectionName)
    activeSynapses = load_active_synapses(synapseActivationName)
    
    somaLocations3DCC = []
    somaLocations3DVPM = []
    preCellIDs = {}
    preCellTypes = activeSynapses.keys()
    for preCellType in preCellTypes:
        compareCellType = preCellType.split('_')[0]
        if cellType is not None and compareCellType != cellType:
            continue
        print 'Loading locations and activation times for cell type %s' % preCellType
        for syn in activeSynapses[preCellType]:
            synID, somaDist, synStructure, synTimes = syn
            if structure is not None and synStructure != structure:
                continue
            if tBegin is not None and tEnd is not None:
                activeSyn = False
                for t in synTimes:
                    if tBegin <= t <= tEnd:
                        activeSyn = True
                        break
                if not activeSyn:
                    continue
            cellID = connectedNeuronsPerSynapse[preCellType][synID]
            if not preCellIDs.has_key(preCellType):
                preCellIDs[preCellType] = []
            else:
                if cellID in preCellIDs[preCellType]:
                    continue
            preCellIDs[preCellType].append(cellID)
            somaLoc3D = somaLocations[preCellType][cellID]
            if 'VPM' in preCellType:
                somaLocations3DVPM.append(somaLoc3D)
            else:
                somaLocations3DCC.append(somaLoc3D)
    
    CCOutName = outputFilename + '_intracortical'
    VPMOutName = outputFilename + '_VPM'
    scp.write_landmark_file(CCOutName, somaLocations3DCC)
    scp.write_landmark_file(VPMOutName, somaLocations3DVPM)

def load_connected_neuron_per_synapse(fname):
    '''
    returns:
    dict of
    cell type, list of presynpatic IDs (list index = synapse index)
    '''
    connectedNeurons = {}
    with open(fname, 'r') as connectionFile:
        for line in connectionFile:
            line = line.strip()
            if not line:
                continue
            if line[0] == '#':
                continue
            splitLine = line.split()
            cellType = splitLine[0]
            cellID = int(splitLine[1])
            synID = int(splitLine[2])
            if not connectedNeurons.has_key(cellType):
                connectedNeurons[cellType] = [cellID]
            else:
                connectedNeurons[cellType].append(cellID)
    
    return connectedNeurons

def load_active_synapses(fname):
    '''
    reads list of all functional synapses and their activation times.
    Input: file of format:
        synapse type\\tsynapse ID\\tsoma distance\\tsection ID\\tsection pt ID\\tdendrite label\\tactivation times
    returns: dictionary with cell types as keys and list of synapse locations and activation times,
    coded as tuples: (synapse ID, soma distance, structure label, [t1, t2, ... , tn])
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
#            secID = int(splitLine[3])
#            ptID = int(splitLine[4])
            structure = splitLine[5]
            synTimes = []
            synTimesStr = splitLine[6].split(',')
            for tStr in synTimesStr:
                if tStr:
                    synTimes.append(float(tStr))
            if not synapses.has_key(cellType):
                synapses[cellType] = [(synID, somaDist, structure, synTimes)]
            else:
                synapses[cellType].append((synID, somaDist, structure, synTimes))
    
    return synapses

def load_cell_location_realization(fname):
    connectedSomata = {}
    with open(fname, 'r') as connectedSomataFile:
        for line in connectedSomataFile:
            line = line.strip()
            if not line:
                continue
            if line[0] == '#':
                continue
            splitLine = line.split()
            cellType = splitLine[0]
#            cellID = int(splitLine[1])
            x = float(splitLine[2])
            y = float(splitLine[3])
            z = float(splitLine[4])
            somaLoc = x, y, z
            if not connectedSomata.has_key(cellType):
                connectedSomata[cellType] = []
            connectedSomata[cellType].append(somaLoc)
    
    return connectedSomata

if __name__ == '__main__':
    if len(sys.argv) >= 5:
        somaLocationName = sys.argv[1]
        connectionName = sys.argv[2]
        synapseActivationName = sys.argv[3]
        outputFilename = sys.argv[4]
        cellType=None
        tBegin=None
        tEnd=None
        structure=None
        if len(sys.argv) == 6:
            cellType = sys.argv[5]
        elif len(sys.argv) == 7:
            tBegin = float(sys.argv[5])
            tEnd = float(sys.argv[6])
        elif len(sys.argv) == 8:
            try:
                cellType = sys.argv[5]
                tBegin = float(sys.argv[6])
                tEnd = float(sys.argv[7])
            except ValueError:
                tBegin = float(sys.argv[5])
                tEnd = float(sys.argv[6])
                structure = sys.argv[7]
        elif len(sys.argv) == 9:
            cellType = sys.argv[5]
            tBegin = float(sys.argv[6])
            tEnd = float(sys.argv[7])
            structure = sys.argv[8]
        visualize_active_soma_locations(somaLocationName, connectionName, synapseActivationName, outputFilename, cellType, tBegin, tEnd, structure)
    else:
        print 'Error! Number of arguments is %d; should be at least 4' % (len(sys.argv)-1)
    
    
    
    