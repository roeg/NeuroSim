import sys
import numpy as np
import singlecell_input_mapper as sm
import os, os.path
import glob

#1. define cell type names
cellTypes = ['L2axon','L34axon','L4pyaxon','L4spaxon','L4ssaxon',\
            'L5staxon','L5ttaxon','L6ccaxon','L6ccinvaxon','L6ctaxon','VPM']
#2. define cell type PSTHs
cellTypePSTHs = {'L2axon': [0.0214,0.0086,0.0057,0.0086,0.0229,0.0271,0.0243,0.0143,0.0086],\
                'L34axon': [0.0143,0.0000,0.0029,0.0000,0.1314,0.0229,0.0286,0.0029,0.0000],\
                'L4pyaxon': [0.0000,0.0750,0.0750,0.0750,0.0250,0.0250,0.0750,0.0250,0.0250],\
                'L4spaxon': [0.0000,0.0088,0.0125,0.1525,0.1338,0.0175,0.0088,0.0000,0.0063],\
                'L4ssaxon': [0.0000,0.0088,0.0125,0.1525,0.1338,0.0175,0.0088,0.0000,0.0063],\
                'L5staxon': [0.0377,0.0377,0.0215,0.0085,0.0254,0.0269,0.0238,0.0262,0.0046],\
                'L5ttaxon': [0.1556,0.2000,0.1000,0.3389,0.3333,0.2167,0.3556,0.2167,0.0889],\
                'L6ccaxon': [0.0417,0.0333,0.0333,0.3583,0.4250,0.1167,0.0000,0.0167,0.0167],\
                'L6ccinvaxon': [0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0500,0.0000,0.0000],\
                'L6ctaxon': [0.0000,0.0100,0.0000,0.0100,0.0400,0.0000,0.0000,0.0000,0.0600],\
                'VPM': [0.0252,0.0332,0,0.0492,0.5449,0.054,0,0.012,0]}
#3. define SuC LUT
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
surroundPSTHLookup = {1: 8, 2: 7, 3: 6, 4: 5, 5: 4, 6: 3, 7: 2, 8: 1, 9: 0}

# approximate release probability of exc. synapses
pGlut = 0.6

def main(whisker, boutonFolder, outName):
    fnames = []
    scan_directory(boutonFolder, fnames, 'Boutons.am')
    #4. go through all cell types
    for cellType in cellTypes:
        cellTypeDensities = []
        cellTypeBounds = [1e6,-1e6,1e6,-1e6,1e6,-1e6] # min, max
        cellTypeSpacing = None
        activeColumns = []
    #5. go through all surrounding columns
        SuCs = surroundColumns[whisker].keys()
        for column in SuCs:
            spikeProb = cellTypePSTHs[cellType][surroundPSTHLookup[surroundColumns[whisker][column]]]
            if not spikeProb:
                continue
            activeColumns.append(column)
    #6. load bouton density
            fileName = None
            for fname in fnames:
                if cellType in fname and column in fname:
                    fileName = fname
                    break
            if fileName is None:
                errstr = 'Bouton density of cell type %s in column %s not found' % (cellType, column)
                raise IOError(errstr)
            tmpDensity = sm.read_scalar_field(fileName)
    #7. scale bouton density by PSTH
            tmpDensity.mesh *= spikeProb # without pGlut for now; does not really change anything
            cellTypeDensities.append(tmpDensity)
    #8. add to cell type bouton density
        for i in range(len(cellTypeDensities)):
            if not i:
                cellTypeSpacing = cellTypeDensities[i].spacing
            tmpBounds = cellTypeDensities[i].boundingBox
            for j in range(3):
                if tmpBounds[2*j] < cellTypeBounds[2*j]:
                    cellTypeBounds[2*j] = tmpBounds[2*j]
                if tmpBounds[2*j+1] > cellTypeBounds[2*j+1]:
                    cellTypeBounds[2*j+1] = tmpBounds[2*j+1]
        dims, extent, origin = [], [], []
        for i in range(3):
            dims.append(int((cellTypeBounds[2*i+1]-cellTypeBounds[2*i])/cellTypeSpacing[i]+1.001))
            extent.append(0)
            extent.append(dims[i]-1)
            origin.append(cellTypeBounds[2*i])
        cellTypeMesh = np.zeros(shape=dims)
        for i in range(len(cellTypeDensities)):
            tmpBounds = cellTypeDensities[i].boundingBox
            offset = []
            for j in range(3):
                minOffset = int((tmpBounds[2*j]-cellTypeBounds[2*j])/cellTypeSpacing[j]+0.001)
                #maxOffset = extent[2*j+1] - int((cellTypeBounds[2*j+1]-tmpBounds[2*j+1])/cellTypeSpacing[j]+0.001)
                maxOffset = minOffset + cellTypeDensities[i].extent[2*j+1]
                diff = (maxOffset - minOffset) - (cellTypeDensities[i].extent[2*j+1])
                if diff:
                    errstr = 'Offset difference does not match size of cell type density: maxOffset = %d - minOffset = %d - dim = %d' % (maxOffset, minOffset, cellTypeDensities[i].extent[2*j+1])
                    raise RuntimeError(errstr)
                offset.append(minOffset)
                offset.append(maxOffset)
            cellTypeMesh[offset[0]:offset[1]+1,offset[2]:offset[3]+1,offset[4]:offset[5]+1] += cellTypeDensities[i].mesh
    #9. write cell type-specific bouton density
        cellTypeScalarField = sm.ScalarField(cellTypeMesh, origin, extent, cellTypeSpacing, cellTypeBounds)
        cellTypeOutName = outName
        cellTypeOutName += '_'
        cellTypeOutName += whisker
        cellTypeOutName += '_'
        cellTypeOutName += cellType
        cellTypeOutName += '_active_synapses'
        sm.write_scalar_field(cellTypeOutName, cellTypeScalarField)

def scan_directory(path, fnames, suffix):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory(fname, fnames, suffix)
        elif fname.endswith(suffix):
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 4:
        whisker = sys.argv[1]
        boutonFolder = sys.argv[2]
        outName = sys.argv[3]
        main(whisker, boutonFolder, outName)
    else:
        print 'Wrong number of arguments: [whisker] [bouton folder] [output filename]'