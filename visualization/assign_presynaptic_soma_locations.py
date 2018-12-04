'''
@author: robert
'''

import sys
import os, os.path
import glob
import numpy as np
import matplotlib.pyplot as plt
import single_cell_parser as scp

def assign_3D_soma_locations(cellName, networkParamName, somaDistributionName, VPMFolderName, connectedSomataName):
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
    
    somaDistribution = load_soma_distribution(somaDistributionName)
    somaDistribution['VPM'] = load_soma_distribution_VPM(VPMFolderName)
#    connectedSomata = load_connected_somata(connectedSomataName)
#    
#    for cellType in connectedSomata.keys():
#        somaLocations3D = []
#        for column in connectedSomata[cellType].keys():
#            numberOfConnectedCells = len(connectedSomata[cellType][column])
#            numberOfTotalCells = len(somaDistribution[cellType][column])
#            somaLocationIDs = np.random.permutation(numberOfTotalCells)
#            for i in range(numberOfConnectedCells):
#                ID = somaLocationIDs[i]
#                loc = somaDistribution[cellType][column][ID]
#                somaLocations3D.append(loc)
#        cellTypeName = connectedSomataName[:-4]
#        cellTypeName += '_'
#        cellTypeName += cellType
#        scp.write_landmark_file(cellTypeName, somaLocations3D)
    
    print 'Loading cell morphology...'
    parser = scp.CellParser(cellParam.filename)
    parser.spatialgraph_to_cell()
    cell = parser.cell
    
    evokedNW = scp.NetworkMapper(cell, networkParameters.network)
    evokedNW._assign_anatomical_synapses()
    evokedNW._create_presyn_cells()
    evokedNW._map_complete_anatomical_realization()
    
    somaLocationsPerCellType3D = {}
    somaLocationsPerSynType3D = {}
    for synType in evokedNW.cells.keys():
        cellType, column = synType.split('_')
        numberOfConnectedCells = len(evokedNW.cells[synType])
        numberOfTotalCells = len(somaDistribution[cellType][column])
        somaLocationIDs = np.random.permutation(numberOfTotalCells)
        for i in range(numberOfConnectedCells):
            ID = somaLocationIDs[i]
            loc = somaDistribution[cellType][column][ID]
            try:
                somaLocationsPerCellType3D[cellType].append(loc)
            except KeyError:
                somaLocationsPerCellType3D[cellType] = [loc]
            try:
                somaLocationsPerSynType3D[synType].append(loc)
            except KeyError:
                somaLocationsPerSynType3D[synType] = [loc]
    
    for cellType in somaLocationsPerCellType3D.keys():
        cellTypeName = connectedSomataName[:-4]
        cellTypeName += '_'
        cellTypeName += cellType
        scp.write_landmark_file(cellTypeName, somaLocationsPerCellType3D[cellType])
    write_cell_location_realization(somaLocationsPerSynType3D, connectedSomataName)

def load_soma_distribution(somaDistributionName):
    '''
    input:
    PSI Format V1.1

    column[0] = "x"
    column[1] = "y"
    column[2] = "z"
    column[3] = "Excitatory"
    column[4] = "NearestColumnId"
    column[5] = "Id"
    column[6] = "InsideColumn"
    column[7] = "CellTypeId"
    
    returns:
    list of
    cell type, column, 3D soma location
    '''
    if not somaDistributionName.endswith('.psi'):
        errstr = 'Wrong input file format: has to be PSI file format!'
        raise RuntimeError(errstr)
    
    columnMap = {1: 'Alpha', 2: 'A1', 3: 'A2', 4: 'A3', 5: 'A4', 6: 'Beta', 7: 'B1', 8: 'B2', 9: 'B3', 10: 'B4',\
                 11: 'Gamma', 12: 'C1', 13: 'C2', 14: 'C3', 15: 'C4', 16: 'Delta', 17: 'D1', 18: 'D2', 19: 'D3', 20: 'D4',\
                 21: 'E1', 22: 'E2', 23: 'E3', 24: 'E4'}
    cellTypeMap = {1: 'L2', 2: 'L34', 3: 'L4py', 4: 'L4sp', 5: 'L4ss', 6: 'L5st', 7: 'L5tt', 8: 'L6cc', 9: 'L6ccinv',\
                   10: 'L6ct', 22: 'SymLocal1', 23: 'SymLocal2', 24: 'SymLocal3', 25: 'SymLocal4', 26: 'SymLocal5',\
                   27: 'SymLocal6', 28: 'L1', 29: 'L23Trans', 30: 'L45Sym', 31: 'L45Peak', 32: 'L56Trans'}
    
    somaDistribution = {}
    with open(somaDistributionName, 'r') as somaDistributionFile:
        for line in somaDistributionFile:
            line = line.strip()
            if not line:
                continue
            if line[0] == '#':
                continue
            splitLine = line.split()
            if len(splitLine) == 3:
                continue
            x = float(splitLine[0])
            y = float(splitLine[1])
            z = float(splitLine[2])
#            exc = int(splitLine[3])
            column = columnMap[int(splitLine[4])]
#            cellID = int(splitLine[5])
#            insideColumn = int(splitLine[6])
            cellTypeNumber = int(splitLine[7])
            if not cellTypeNumber in cellTypeMap.keys():
                continue
            cellType = cellTypeMap[cellTypeNumber]
            somaLocation = x, y, z
            try:
                somaDistribution[cellType][column].append(somaLocation)
            except KeyError:
                try:
                    somaDistribution[cellType][column] = [somaLocation]
                except KeyError:
                    somaDistribution[cellType] = {column: [somaLocation]}
    
    return somaDistribution

def load_soma_distribution_VPM(VPMFolderName):
    '''
    input:
    folder name where counted VPM somata are located,
    separately for each barreloid
    
    returns:
    list of
    cell type, barreloid (=column), 3D soma location (in VPM coordinates)
    '''
    barreloidNames = ['Alpha','A1','A2','A3','A4','Beta','B1','B2','B3','B4','Gamma',\
                      'C1','C2','C3','C4','Delta','D1','D2','D3','D4','E1','E2','E3','E4']
    
    somaDistribution = {}
    for barreloidName in barreloidNames:
        suffix = 'somata_%s.landmarkAscii' % barreloidName
        landmarkName = None
        for fname in glob.glob(os.path.join(VPMFolderName, '*')):
            if fname.endswith(suffix):
                landmarkName = fname
        if landmarkName is None:
            errstr = 'No landmark file for somata in barreloid %s found' % barreloidName
            raise RuntimeError(errstr)
        
        VPMsomata = load_landmark_file(landmarkName)
        somaDistribution[barreloidName] = VPMsomata
    
    return somaDistribution

def load_landmark_file(landmarkFilename):
    '''
    returns list of (x,y,z) points
    '''
    if not landmarkFilename.endswith('.landmarkAscii'):
        errstr = 'Wrong input format: has to be landmarkAscii format'
        raise RuntimeError(errstr)
    
    landmarks = []
    with open(landmarkFilename, 'r') as landmarkFile:
        readPoints = False
        for line in landmarkFile:
            line.strip()
            if not line:
                continue
            if line[:2] == '@1':
                readPoints = True
                continue
            if readPoints:
                splitLine = line.split()
                x = float(splitLine[0])
                y = float(splitLine[1])
                z = float(splitLine[2])
                landmarks.append((x,y,z))
    
    return landmarks

def load_connected_somata(connectedSomataName):
    '''
    returns:
    list of
    cell type, column, unique cell ID
    '''
    connectedSomata = {}
    with open(connectedSomataName, 'r') as connectedSomataFile:
        for line in connectedSomataFile:
            line = line.strip()
            if not line:
                continue
            if line[0] == '#':
                continue
            splitLine = line.split()
            cellTypeColumn = splitLine[0]
            cellType = cellTypeColumn.split('_')[0]
            column = cellTypeColumn.split('_')[1]
            cellID = splitLine[1]
#            if cellType == 'VPM':
#                continue
            if not connectedSomata.has_key(cellType):
                connectedSomata[cellType] = {}
            if not connectedSomata[cellType].has_key(column):
                connectedSomata[cellType][column] = []
            if cellID in connectedSomata[cellType][column]:
                continue
            connectedSomata[cellType][column].append(cellID)
    
    return connectedSomata

def write_cell_location_realization(somaLocations, connectedSomataName):
    outputName = connectedSomataName[:-4]
    outputName += '_soma_locations.csv'
    with open(outputName, 'w') as outputFile:
        header = '# 3D presynaptic soma location file; only valid with input map realization:\n'
        header += '# '
        header += connectedSomataName
        header += '\n'
        header += '# Type - cell ID - soma x - soma y - soma z\n\n'
        outputFile.write(header)
        for synType in somaLocations.keys():
            for i in range(len(somaLocations[synType])):
                loc = somaLocations[synType][i]
                line = synType
                line += '\t'
                line += str(i)
                line += '\t'
                line += '%.3f\t%.3f\t%.3f\n' % (loc[0],loc[1],loc[2])
                outputFile.write(line)

if __name__ == '__main__':
    if len(sys.argv) == 6:
        cellName = sys.argv[1]
        networkParamName = sys.argv[2]
        somaDistributionName = sys.argv[3]
        VPMFolderName = sys.argv[4]
        connectedSomataName = sys.argv[5]
        assign_3D_soma_locations(cellName, networkParamName, somaDistributionName, VPMFolderName, connectedSomataName)
    else:
        print 'Error! Number of arguments is %d; should be 5' % (len(sys.argv)-1)
    
    
    
    