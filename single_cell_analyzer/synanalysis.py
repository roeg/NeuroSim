'''
Created on Apr 2, 2012

@author: regger
'''

import numpy as np
import single_cell_parser as scp

def compute_synapse_distances_times(fname, cell, t=None, synTypes=None):
    synDistances = {}
    synTimes = {}
    activeSyns = {}
    synCnt = 0
    if synTypes is None:
        synTypes = cell.synapses.keys()
        synTypes.sort
    for synType in synTypes:
        synDistances[synType] = compute_syn_distances(cell, synType)
        synTimes[synType] = []
        activeSyns[synType] = []
        for syn in cell.synapses[synType]:
            if not syn.is_active():
                activeSyns[synType].append(False)
                synTimes[synType].append([])
            else:
                activeSyns[synType].append(True)
                tmpSynTimes = syn.releaseSite.spikeTimes[:]
                synTimes[synType].append(tmpSynTimes)
            synCnt += 1
    
    scp.write_synapse_activation_file(fname, cell, synTypes, synDistances, synTimes, activeSyns)

def synapse_activation_times(tVec, cntVec):
    synTVec = []
    for i in range(1, len(cntVec)):
        if cntVec[i] > cntVec[i-1]:
            synTVec.append(tVec[i])
    return synTVec

def synapse_distances(pname):
    parameters = scp.build_parameters(pname)
    cellParam = parameters.network.post
    preParam = parameters.network.pre
    
    parser = scp.CellParser(cellParam.filename)
    parser.spatialgraph_to_cell()
    cell = parser.cell
    for preType in preParam.keys():
        synapseFName = preParam[preType].synapses.realization
        synDist = scp.read_synapse_realization(synapseFName)
        mapper = scp.SynapseMapper(cell, synDist)
        mapper.map_synapse_realization()
    
    for synType in cell.synapses.keys():
        dist = compute_syn_distances(cell, synType)
        name = parameters.info.outputname
        name += '_'
        name += synType
        name += '_syn_distances.csv'
        with open(name, 'w') as distFile:
            header = 'Distance to soma (micron)\n'
            distFile.write(header)
            for d in dist:
                distFile.write(str(d)+'\n')

def synapse_distances_2D(pname):
    parameters = scp.build_parameters(pname)
    cellParam = parameters.network.post
    preParam = parameters.network.pre
    
    parser = scp.CellParser(cellParam.filename)
    parser.spatialgraph_to_cell()
    cell = parser.cell
    for preType in preParam.keys():
        synapseFName = preParam[preType].synapses.realization
        synDist = scp.read_synapse_realization(synapseFName)
        mapper = scp.SynapseMapper(cell, synDist)
        mapper.map_synapse_realization()
    
    for synType in cell.synapses.keys():
        dist = compute_syn_distances_2Dprojected(cell, synType)
        name = parameters.info.outputname
        name += '_'
        name += synType
        name += '_syn_distances_2Dprojected.csv'
        with open(name, 'w') as distFile:
            header = 'Distance to soma (micron)\n'
            distFile.write(header)
            for d in dist:
                distFile.write(str(d)+'\n')

def compute_syn_distances(cell, synType, label=None):
    '''
    computes distances of all synapses on dendrite w.r.t. soma
    
    cell is cell object with attached synapses
    presynaptic cell type given by synType (string)
    optional: dendrite type given by label (string)
    
    returns 1D numpy array of distances to soma
    '''
    if not cell.synapses.has_key(synType):
        errStr = 'Cell does not have synapses of type %s' % synType
        raise KeyError(errStr)
    
    distances = []
    for syn in cell.synapses[synType]:
        currentSec = cell.sections[syn.secID]
        if label is not None and currentSec.label != label:
            continue
        
        if currentSec.label == 'Soma':
            dist = 0.0
            distances.append(dist)
            continue
        
        parentSec = currentSec.parent
        '''compute distance from synapse location to parent branch first'''
        dist = 0.0
        dist = syn.x*currentSec.L
        parentLabel = parentSec.label
        while parentLabel != 'Soma':
            dist += parentSec.L
            currentSec = parentSec
            parentSec = currentSec.parent
            parentLabel = parentSec.label
        distances.append(dist)
    
    return np.array(distances)

def compute_syn_distances_2Dprojected(cell, synType, label=None):
    '''
    computes distances of all synapses on dendrite w.r.t. soma
    projected on 2D plane as seen during 2-photon spine imaging
    
    cell is cell object with attached synapses
    presynaptic cell type given by synType (string)
    optional: dendrite type given by label (string)
    
    returns 1D numpy array of distances to soma
    '''
    if not cell.synapses.has_key(synType):
        errStr = 'Cell does not have synapses of type %s' % synType
        raise KeyError(errStr)
    
    somaLoc = cell.soma.pts[cell.soma.nrOfPts//2]
    distances = []
    for syn in cell.synapses[synType]:
        currentSec = cell.sections[syn.secID]
        if label is not None and currentSec.label != label:
            continue
        
        synLoc = syn.coordinates
        diff = synLoc - somaLoc
        dist = np.sqrt(np.dot(diff[:2], diff[:2]))
        distances.append(dist)
    
    return np.array(distances)

