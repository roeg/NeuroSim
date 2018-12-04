'''
@author: robert

convert input map of single neuron to
CortexInSilico format (test import to CATMAID)
'''

import sys
import os, os.path
import glob
import numpy as np
import matplotlib.pyplot as plt
import single_cell_parser as scp

def single_neuron_embedding_to_catmaid(cellName, networkParamName, somaDistributionName, connectedNeuronsName, outputName):
    '''
    randomly assign 3D locations of corresponding cell type/column
    Assumes there are >= somata than connected somata!!!
    (i.e. best to use soma distribution from corresponding NeuroNet output)
    '''
    neuronParameters = scp.build_parameters(cellName)
    networkParameters = scp.build_parameters(networkParamName)
    scp.load_NMODL_parameters(neuronParameters)
    scp.load_NMODL_parameters(networkParameters)
    cellParam = neuronParameters.neuron
    
    somaDistribution = load_cell_location_realization(somaDistributionName)
    synapseConnectivity = load_connected_neuron_per_synapse(connectedNeuronsName)
    
    print 'Loading cell morphology...'
    parser = scp.CellParser(cellParam.filename)
    parser.spatialgraph_to_cell()
    cell = parser.cell
    
    evokedNW = scp.NetworkMapper(cell, networkParameters.network)
    evokedNW._assign_anatomical_synapses()
    evokedNW._create_presyn_cells()
    evokedNW._map_complete_anatomical_realization()
    
    # create continuous morphology IDs ('DUMMYIDXXXXX')
    IDCounter = 0
    cellTypeColumns = somaDistribution.keys()
    for cellTypeColumn in cellTypeColumns:
        cellType = cellTypeColumn.split('_')[0]
        column = cellTypeColumn.split('_')[1]
        if len(somaDistribution[cellTypeColumn]) != len(evokedNW.cells[cellTypeColumn]):
            errstr = 'Number of cells in Network not equal to number of somata'
            raise RuntimeError(errstr)
        for i in range(len(somaDistribution[cellTypeColumn])):
            somaLoc = somaDistribution[cellTypeColumn][i]
            pointCell = evokedNW.cells[cellTypeColumn][i]
            morphologyID = 'DUMMYID%05d' % IDCounter
            pointCell.cellTypeStr = cellType
            pointCell.columnStr = column
            pointCell.morphologyID = morphologyID
            pointCell.somaLoc3D = somaLoc
            IDCounter += 1
        for i in range(len(synapseConnectivity[cellTypeColumn])):
            syn = cell.synapses[cellTypeColumn][i]
            preCellID = synapseConnectivity[cellTypeColumn][i]
            preCell = evokedNW.cells[cellTypeColumn][preCellID]
            syn.preCell = preCell
            create_catmaid_connection(cell, syn)
            
    # write cell file in JSON format
    morphologyName = outputName
    morphologyName += '_morphologies.json'
    write_catmaid_morphology_file(morphologyName, evokedNW.cells)
    # write synapse file in CortexInSilico/Catmaid format
    synapseName = outputName
    synapseName += '_synapses.csv'
    write_catmaid_synapse_file(synapseName, cell.synapses)

def create_catmaid_connection(cell, synapse):
    '''
    create Catmaid Connector Node
    close to synapse location
    '''
    diameter = cell.sections[synapse.secID].diamList[synapse.ptID]
    norm = 0.0
    while not norm:
        direction = np.random.rand(3)
        direction = 2*direction - 1.0
        norm = np.sqrt(np.dot(direction, direction))
    direction /= norm
    connectorNode = synapse.coordinates + (0.5*diameter + 1.0)*direction
    preNode = synapse.coordinates + (0.5*diameter + 2.0)*direction
    synapse.connectorNode = connectorNode
    synapse.preNode = preNode
    synapse.coordinates = synapse.coordinates + 0.5*diameter*direction

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

def write_catmaid_morphology_file(fname, cells):
    with open(fname, 'w') as outputFile:
        line = '{\n'
        outputFile.write(line)
        for cellType in cells.keys():
            for cell in cells[cellType]:
                line = '\t'
                line += '\"' + cell.morphologyID + '\": {\n'
                line += '\t\t\"annotations\": {\n'
                line += '\t\t\t\"celltype\": \"' + cell.cellTypeStr + '\",\n'
                line += '\t\t\t\"column\": \"' + cell.columnStr + '\",\n'
                line += '\t\t\t\"somalocation\": [' + str(cell.somaLoc3D[0]) + ',' + str(cell.somaLoc3D[1]) + ',' + str(cell.somaLoc3D[2]) + ']\n'
                line += '\t\t}\n'
                line += '\t},\n'
                outputFile.write(line)
        line = '}\n'
        outputFile.write(line)

def write_catmaid_synapse_file(fname, synapses):
    with open(fname, 'w') as outputFile:
        header = 'synapse ID\tpostsynaptic MORPHOLOGYID\tpresynaptic MORPHOLOGYID\tpresynaptic celltype\t3D location (x,y,z)\tpost section\tpost node\t3D connector node location\t3D pre node location\n'
        morphologyID = '86_L5_CDK20041214_nr3L5B_dend_PC_neuron_transform_registered_C2.hoc'
        outputFile.write(header)
        IDCounter = 0
        for cellType in synapses.keys():
            for syn in synapses[cellType]:
                presynapticID = syn.preCell.morphologyID
                presynapticCellType = syn.preCell.cellTypeStr
                line = str(IDCounter)
                line += '\t'
                line += morphologyID
                line += '\t'
                line += presynapticID
                line += '\t'
                line += presynapticCellType
                line += '\t'
                line += str(syn.coordinates[0]) + ',' + str(syn.coordinates[1]) + ',' + str(syn.coordinates[2])
                line += '\t'
                line += str(syn.secID)
                line += '\t'
                line += str(syn.ptID)
                line += '\t'
                line += str(syn.connectorNode[0]) + ',' + str(syn.connectorNode[1]) + ',' + str(syn.connectorNode[2])
                line += '\t'
                line += str(syn.preNode[0]) + ',' + str(syn.preNode[1]) + ',' + str(syn.preNode[2])
                line += '\n'
                outputFile.write(line)
                IDCounter += 1

if __name__ == '__main__':
    if len(sys.argv) == 6:
        cellName = sys.argv[1]
        networkParamName = sys.argv[2]
        somaDistributionName = sys.argv[3]
        connectedNeuronsName = sys.argv[4]
        outputName = sys.argv[5]
        single_neuron_embedding_to_catmaid(cellName, networkParamName, somaDistributionName, connectedNeuronsName, outputName)
    else:
        print 'Error! Number of arguments is %d; should be 5' % (len(sys.argv)-1)
    
    
    
    