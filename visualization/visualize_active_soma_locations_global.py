'''
@author: robert
'''

import sys
import numpy as np
import glob
import os

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
# correspondence between anatomical column
# and whisker PSTH relative to PW whisker
# (e.g, C2 whisker deflection in B1
# looks like D3 whisker deflection in C2)
surroundPSTHLookup = {1: 'D3', 2: 'D2', 3: 'D1', 4: 'C3', 5: 'C2',\
                        6: 'C1', 7: 'B3', 8: 'B2', 9: 'B1'}
surroundPSTHLookup2 = {1: 8, 2: 7, 3: 6, 4: 5, 5: 4, 6: 3, 7: 2, 8: 1, 9: 0}

# 0-25ms post-stimulus
L6ccPSTH = {'B1': 0.042, 'B2': 0.033, 'B3': 0.033, 'C1': 0.358, 'C2': 0.425, 'C3': 0.117,\
            'D1': 0.000, 'D2': 0.017, 'D3': 0.017}
cellTypePSTH = {'L2': [0.0214,0.0086,0.0057,0.0086,0.0229,0.0271,0.0243,0.0143,0.0086],\
    'L34': [0.0143,0.0000,0.0029,0.0000,0.1314,0.0229,0.0286,0.0029,0.0000],\
    'L4py': [0.0000,0.0750,0.0750,0.0750,0.0250,0.0250,0.0750,0.0250,0.0250],\
    'L4sp': [0.0000,0.0088,0.0125,0.1525,0.1338,0.0175,0.0088,0.0000,0.0063],\
    'L4ss': [0.0000,0.0088,0.0125,0.1525,0.1338,0.0175,0.0088,0.0000,0.0063],\
    'L5st': [0.0377,0.0377,0.0215,0.0085,0.0254,0.0269,0.0238,0.0262,0.0046],\
    'L5tt': [0.1556,0.2000,0.1000,0.3389,0.3333,0.2167,0.3556,0.2167,0.0889],\
    'L6cc': [0.0417,0.0333,0.0333,0.3583,0.4250,0.1167,0.0000,0.0167,0.0167],\
    'L6ccinv': [0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0500,0.0000,0.0000],\
    'L6ct': [0.0000,0.0100,0.0000,0.0100,0.0400,0.0000,0.0000,0.0000,0.0600],\
    'VPM': [0.0252,0.0332,0.0000,0.0492,0.5449,0.0540,0.0000,0.0120,0.0000]}
cellTypeSpont25 = {'L2': 0.0115,\
    'L34': 0.0077,\
    'L4py': 0.016,\
    'L4sp': 0.015,\
    'L4ss': 0.0132,\
    'L5st': 0.0275,\
    'L5tt': 0.0882,\
    'L6cc': 0.0238,\
    'L6ccinv': 0.0107,\
    'L6ct': 0.002,\
    'VPM': 0.0053}
                    
def visualize_active_soma_locations(cellType, PWName, somaDistributionName, VPMFolderName, outName):
    '''
    write landmark file for all active
    somata of specific cell type for a given PSTH
    and deflected whisker
    '''
    
    somaLocations = load_soma_distribution(somaDistributionName)
    somaLocations['VPM'] = load_soma_distribution_VPM(VPMFolderName)
    
    activeSomaLocations = []
    for column in somaLocations[cellType]:
        if not surroundColumns[PWName].has_key(column):
            #continue
            spikeProb = cellTypeSpont25[cellType]
        #spikeProb = L6ccPSTH[surroundPSTHLookup[surroundColumns[PWName][column]]]
        else:
            PSTH = cellTypePSTH[cellType]
            spikeProb = PSTH[surroundPSTHLookup2[surroundColumns[PWName][column]]]
        for loc in somaLocations[cellType][column]:
            if np.random.rand() < spikeProb:
                activeSomaLocations.append(loc)
    
    write_landmark_file(outName, activeSomaLocations)

def active_neurons_per_celltype(PWName, somaDistributionName, VPMFolderName, outName):
    '''
    write landmark file for all active
    somata of specific cell type for a given PSTH
    and deflected whisker
    '''
    
    somaLocations = load_soma_distribution(somaDistributionName)
    somaLocations['VPM'] = load_soma_distribution_VPM(VPMFolderName)
    
    activeSomaLocations = []
    cellTypes = cellTypePSTH.keys()
    cellTypes.sort()
    PCCelltypeNumbers = {}
    SuCCelltypeNumbers = {}
    for cellType in cellTypes:
        for column in somaLocations[cellType]:
            if not surroundColumns[PWName].has_key(column):
                #spikeProb = cellTypeSpont25[cellType]
                continue
            else:
                PSTH = cellTypePSTH[cellType]
                spikeProb = PSTH[surroundPSTHLookup2[surroundColumns[PWName][column]]]
            activeSomata = int(spikeProb*len(somaLocations[cellType][column]) + 0.5)
            if PWName == column:
                PCCelltypeNumbers[cellType] = activeSomata
            else:
                if not SuCCelltypeNumbers.has_key(cellType):
                    SuCCelltypeNumbers[cellType] = activeSomata
                else:
                    SuCCelltypeNumbers[cellType] += activeSomata
    
    if not outName.endswith('.csv'):
        outName += '.csv'
    with open(outName, 'w') as outFile:
        header = 'Cell type\tPC\tSuCs\n'
        outFile.write(header)
        for cellType in cellTypes:
            line = cellType
            line += '\t'
            line += str(PCCelltypeNumbers[cellType])
            line += '\t'
            line += str(SuCCelltypeNumbers[cellType])
            line += '\n'
            outFile.write(line)

def celltype_population_centerofmass(somaDistributionName, outName):
    '''
    write table of format celltype_column [x y z]
    '''
    somaLocations = load_soma_distribution(somaDistributionName)
    avgLocations = {}
    for cellType in somaLocations:
        for column in somaLocations[cellType]:
            locations = somaLocations[cellType][column]
            populationName = cellType + '_' + column
            avgLocations[populationName] = np.mean(locations, axis=0)
    
    if not outName.endswith('.csv'):
        outName += '.csv'
    with open(outName, 'w') as outFile:
        header = 'population\tcenter of mass\n'
        outFile.write(header)
        for population in avgLocations:
            loc = avgLocations[population]
            line = population
            line += '\t'
            line += str(loc[0])
            line += '\t'
            line += str(loc[1])
            line += '\t'
            line += str(loc[2])
            line += '\n'
            outFile.write(line)

def whisker_evoked_PSTH(column, deflectedWhisker, cellType):
    columns = surroundColumns[deflectedWhisker].keys()
    evokedTypes = evokedTemplates.keys()
    if column not in columns or cellType not in evokedTypes:
        return None
    evokedTemplate = evokedTemplates[cellType]
    PSTHwhisker = surroundPSTHLookup[surroundColumns[deflectedWhisker][column]]
    PSTHstr = cellType + '_' + PSTHwhisker
    PSTH = evokedTemplate[PSTHstr]
    return PSTH

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

def write_landmark_file(fname=None, landmarkList=None):
    '''
    write Amira landmark file
    landmarkList has to be iterable of tuples,
    each of which holds 3 float coordinates
    '''
    if fname is None:
        err_str = 'No landmark output file name given'
        raise RuntimeError(err_str)
    
    if not landmarkList:
        print 'Landmark list empty!'
        return
    nrCoords = len(landmarkList[0])
    if nrCoords != 3:
        err_str = 'Landmarks have wrong format! Number of coordinates is ' + str(nrCoords) + ', should be 3'
        raise RuntimeError(err_str)
    
    if not fname.endswith('.landmarkAscii'):
        fname += '.landmarkAscii'
    
    with open(fname, 'w') as landmarkFile:
        nrOfLandmarks = len(landmarkList)
        header = '# AmiraMesh 3D ASCII 2.0\n\n'\
                'define Markers ' + str(nrOfLandmarks) + '\n\n'\
                'Parameters {\n'\
                '\tNumSets 1,\n'\
                '\tContentType \"LandmarkSet\"\n'\
                '}\n\n'\
                'Markers { float[3] Coordinates } @1\n\n'\
                '# Data section follows\n'\
                '@1\n'
        landmarkFile.write(header)
        for pt in landmarkList:
            line = '%.6f %.6f %.6f\n' % (pt[0], pt[1], pt[2])
            landmarkFile.write(line)

if __name__ == '__main__':
    if len(sys.argv) == 6:
        cellType = sys.argv[1]
        PWName = sys.argv[2]
        somaDistributionName = sys.argv[3]
        VPMFolderName = sys.argv[4]
        outName = sys.argv[5]
        visualize_active_soma_locations(cellType, PWName, somaDistributionName, VPMFolderName, outName)
    elif len(sys.argv) == 5:
        PWName = sys.argv[1]
        somaDistributionName = sys.argv[2]
        VPMFolderName = sys.argv[3]
        outName = sys.argv[4]
        active_neurons_per_celltype(PWName, somaDistributionName, VPMFolderName, outName)
    elif len(sys.argv) == 3:
        somaDistributionName = sys.argv[1]
        outName = sys.argv[2]
        celltype_population_centerofmass(somaDistributionName, outName)
    else:
        print 'Error! Number of arguments is %d; should be at least 2' % (len(sys.argv)-1)
    
    
    
    