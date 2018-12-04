'''
Created on Mar 8, 2012

@author: regger
'''

labels2int = {\
    "Neuron":                 2,\
    "Dendrite":               3,\
    "ApicalDendrite":         4,\
    "BasalDendrite":          5,\
    "Axon":                   6,\
    "AIS":                    6,\
    "Myelin":                 6,\
    "Node":                   6,\
    "Soma":                   7,\
    "Landmark":               8,\
    "Pia":                    9,\
    "WhiteMatter":           48,\
    "Vessel":                10,\
    "Barrel":                11,\
    "ZAxis":                 50,\
    "aRow":                  12,\
    "A1":                    13,\
    "A2":                    14,\
    "A3":                    15,\
    "A4":                    16,\
    "bRow":                  17,\
    "B1":                    18,\
    "B2":                    19,\
    "B3":                    20,\
    "B4":                    21,\
    "cRow":                  22,\
    "C1":                    23,\
    "C2":                    24,\
    "C3":                    25,\
    "C4":                    26,\
    "C5":                    27,\
    "C6":                    28,\
    "dRow":                  29,\
    "D1":                    30,\
    "D2":                    31,\
    "D3":                    32,\
    "D4":                    33,\
    "D5":                    34,\
    "D6":                    35,\
    "eRow":                  36,\
    "E1":                    37,\
    "E2":                    38,\
    "E3":                    39,\
    "E4":                    40,\
    "E5":                    41,\
    "E6":                    42,\
    "greekRow":              43,\
    "Alpha":                 44,\
    "Beta":                  45,\
    "Gamma":                 46,\
    "Delta":                 47,\
    "Septum":                 0,\
              }

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

def write_sim_results(fname, t, v):
    with open(fname, 'w') as outputFile:
        header = '# t\tvsoma'
        header += '\n\n'
        outputFile.write(header)
        for i in range(len(t)):
            line = str(t[i])
            line += '\t'
            line += str(v[i])
            line += '\n'
            outputFile.write(line)

def write_all_traces(fname, t, vTraces):
    with open(fname, 'w') as outputFile:
        header = 't'
        for i in range(len(vTraces)):
            header += '\tVm run %02d' % i
        header += '\n'
        outputFile.write(header)
        for i in range(len(t)):
            line = str(t[i])
            for j in range(len(vTraces)):
                line += '\t'
                line += str(vTraces[j][i])
            line += '\n'
            outputFile.write(line)
    
def write_cell_synapse_locations(fname=None, synapses=None, cellID=None):
    '''
    writes list of all synapses with the locations
    coded by section ID and section x of cell with ID 'cellID'
    '''
    if fname is None or synapses is None or cellID is None:
        err_str = 'Incomplete data! Cannot write synapse location file'
        raise RuntimeError(err_str)
    
    with open(fname, 'w') as outputFile:
        header = '# Synapse distribution file\n'
        header += '# corresponding to cell: '
        header += cellID
        header += '\n'
        header += '# Type - section - section.x\n\n'
        outputFile.write(header)
        for synType in synapses.keys():
            for syn in synapses[synType]:
                line = syn.preCellType
                line += '\t'
                line += str(syn.secID)
                line += '\t'
                if syn.x > 1.0:
                    syn.x = 1.0
                if syn.x < 0.0:
                    syn.x = 0.0
                line += str(syn.x)
                line += '\n'
                outputFile.write(line)
    
def write_pruned_synapse_locations(fname=None, synapses=None, cellID=None):
    '''
    writes list of all synapses with the locations
    coded by section ID and section x of cell with ID 'cellID'
    and a pruned flag (1 or 0)
    '''
    if fname is None or synapses is None or cellID is None:
        err_str = 'Incomplete data! Cannot write synapse location file'
        raise RuntimeError(err_str)
    
    with open(fname, 'w') as outputFile:
        header = '# Synapse distribution file\n'
        header += '# corresponding to cell: '
        header += cellID
        header += '\n'
        header += '# Type - section - section.x - pruned\n\n'
        outputFile.write(header)
        for synType in synapses.keys():
            for syn in synapses[synType]:
                line = syn.preCellType
                line += '\t'
                line += str(syn.secID)
                line += '\t'
                if syn.x > 1.0:
                    syn.x = 1.0
                if syn.x < 0.0:
                    syn.x = 0.0
                line += str(syn.x)
                line += '\t'
                line += str(syn.pruned)
                line += '\n'
                outputFile.write(line)

def write_functional_realization_map(fname=None, functionalMap=None, anatomicalID=None):
    '''
    writes list of all functional connections
    coded by tuples (cell type, presynaptic cell index, synapse index).
    Only valid for anatomical synapse realization given by anatomicalID
    '''
    if fname is None or functionalMap is None or anatomicalID is None:
        err_str = 'Incomplete data! Cannot write functional realization file'
        raise RuntimeError(err_str)
    
    if not fname.endswith('.con') and not fname.endswith('.CON'):
        fname += '.con'
    
    with open(fname, 'w') as outputFile:
        header = '# Functional realization file; only valid with synapse realization:\n'
        header += '# ' + anatomicalID
        header += '\n'
        header += '# Type - cell ID - synapse ID\n\n'
        outputFile.write(header)
        for con in functionalMap:
            line = con[0]
            line += '\t'
            line += str(con[1])
            line += '\t'
            line += str(con[2])
            line += '\n'
            outputFile.write(line)

def write_synapse_activation_file(fname=None, cell=None, synTypes=None, synDistances=None, synTimes=None, activeSyns=None):
    if fname is None or cell is None or synTypes is None or synDistances is None or synTimes is None or activeSyns is None:
        err_str = 'Incomplete data! Cannot write functional realization file'
        raise RuntimeError(err_str)
    
    with open(fname, 'w') as outputFile:
        header = '# synapse type\t'
        header += 'synapse ID\t'
        header += 'soma distance\t'
        header += 'section ID\t'
        header += 'section pt ID\t'
        header += 'dendrite label\t'
        header += 'activation times\n'
        outputFile.write(header)
        for synType in synTypes:
            for i in range(len(cell.synapses[synType])):
                if not activeSyns[synType][i]:
                    continue
                secID = cell.synapses[synType][i].secID
                ptID = cell.synapses[synType][i].ptID
                dendLabel = cell.sections[secID].label
                line = synType
                line += '\t'
                line += str(i)
                line += '\t'
                line += str(synDistances[synType][i])
                line += '\t'
                line += str(secID)
                line += '\t'
                line += str(ptID)
                line += '\t'
                line += str(dendLabel)
                line += '\t'
                for t in synTimes[synType][i]:
                    line += str(t)
                    line += ','
                line += '\n'
                outputFile.write(line)

def write_synapse_weight_file(fname=None, cell=None):
    if fname is None or cell is None:
        err_str = 'Incomplete data! Cannot write functional realization file'
        raise RuntimeError(err_str)
    
    with open(fname, 'w') as outputFile:
        header = '# synapse type\t'
        header += 'synapse ID\t'
        header += 'section ID\t'
        header += 'section pt ID\t'
        header += 'receptor type\t'
        header += 'synapse weights\n'
        outputFile.write(header)
        for synType in cell.synapses.keys():
            for i in range(len(cell.synapses[synType])):
                for recepStr in cell.synapses[synType][i].weight.keys():
                    secID = cell.synapses[synType][i].secID
                    ptID = cell.synapses[synType][i].ptID
                    line = synType
                    line += '\t'
                    line += str(i)
                    line += '\t'
                    line += str(secID)
                    line += '\t'
                    line += str(ptID)
                    line += '\t'
                    line += recepStr
                    line += '\t'
                    for g in cell.synapses[synType][i].weight[recepStr]:
                        line += str(g)
                        line += ','
                    line += '\n'
                    outputFile.write(line)

def write_PSTH(fname=None, PSTH=None, bins=None):
    '''
    Write PSTH and time bins of PSTH, where
    bins contain left and right end of each bin,
    i.e. len(bins) = len(PSTH) + 1
    '''
    if fname is None or PSTH is None or bins is None:
        err_str = 'Incomplete data! Cannot write PSTH'
        raise RuntimeError(err_str)
    
    with open(fname, 'w') as outputFile:
        header = '# bin begin\t'
        header += 'bin end\t'
        header += 'APs/trial/bin\n'
        outputFile.write(header)
        for i in range(len(PSTH)):
            line = str(bins[i])
            line += '\t'
            line += str(bins[i+1])
            line += '\t'
            line += str(PSTH[i])
            line += '\n'
            outputFile.write(line)

def write_spike_times_file(fname=None, spikeTimes=None):
    '''
    Write file containing trial numbers and all spike
    times in each trial (may be empty).
    spikeTimes should be dictionary with trial numbers as keys (integers),
    and tuples of spike times in each trial as values.
    '''
    if fname is None or spikeTimes is None:
        err_str = 'Incomplete data! Cannot write spike times file'
        raise RuntimeError(err_str)
    
    with open(fname, 'w') as outFile:
        header = '# trial\tspike times\n'
        outFile.write(header)
        trials = spikeTimes.keys()
        trials.sort()
        for trial in trials:
            line = str(trial)
            line += '\t'
            for tSpike in spikeTimes[trial]:
                line += str(tSpike)
                line += ','
            line += '\n'
            outFile.write(line)

def write_presynaptic_spike_times(fname=None, cells=None):
    '''
    Write cell type, PC and spike times of all connected
    presynaptic point cells
    '''
    if fname is None or cells is None:
        err_str = 'Incomplete data! Cannot write presynaptic spike times'
        raise RuntimeError(err_str)
    
    with open(fname, 'w') as outputFile:
        header = '# presynaptic cell type\tcell ID\tspike times\n'
        outputFile.write(header)
        preTypes = cells.keys()
        preTypes.sort()
        for preType in preTypes:
            for i in range(len(cells[preType])):
                cell = cells[preType][i]
                spikeTimes = cell.spikeTimes
                if not len(cell.spikeTimes):
                    continue
                line = preType
                line += '\t'
                line += str(i)
                line += '\t'
                spikeTimes.sort()
                for t in spikeTimes:
                    line += str(t)
                    line += ','
                line += '\n'
                outputFile.write(line)

def write_cell_simulation(fname=None, cell=None, traces=None, tVec=None, allPoints=False):
    '''
    write Amira SpatialGraph files corresponding to time steps
    of entire simulation run. Recorded quantities are passed
    in tuple traces with strings and recorded in Vectors
    attached to Sections of cell
    TODO: Performs interpolation if nseg != nrOfPts for a Section
    '''
    if fname is None or cell is None or traces is None or tVec is None:
        err_str = 'Incomplete data! Cannot write SpatialGraph simulation results'
        raise RuntimeError(err_str)
    
    axonLabels = ['Axon', 'AIS', 'Myelin']
    
    totalNrPts = 0
    nrOfEdges = 0
    for sec in cell.tree:
        if sec.label in axonLabels:
            continue
        nrOfEdges += 1
        if not allPoints:
            totalNrPts += sec.nseg
        else:
            totalNrPts += sec.nrOfPts
    
    header = "# AmiraMesh 3D ASCII 2.0" + "\n"
    header += "# This SpatialGraph file was created by the Neuron Registration Tool NeuroMap " + "\n"
    header += "# NeuroMap was programmed by Robert Egger," + "\n"
    header += "# Max-Planck-Florida Institute, Jupiter, Florida " + "\n"
    header += "" + "\n"
    header += "define VERTEX " + str(nrOfEdges*2) + "\n"
    header += "define EDGE " + str(nrOfEdges) + "\n"
    header += "define POINT " + str(totalNrPts) + "\n"
    header += "" + "\n"
    header += "Parameters {GraphLabels {" + "\n"
    header += "        Neuron { " + "\n"
    header += "            Dendrite {" + "\n"
    header += "                ApicalDendrite {" + "\n"
    header += "                    Color 1 0.5 0.5," + "\n"
    header += "                    Id 4 }" + "\n"
    header += "                BasalDendrite {" + "\n"
    header += "                    Color 0.8 0.4 0.4," + "\n"
    header += "                    Id 5 }" + "\n"
    header += "                Color 1 0 0," + "\n"
    header += "                Id 3 }" + "\n"
    header += "            Axon {" + "\n"
    header += "                Color 0 0 1," + "\n"
    header += "                Id 6 }" + "\n"
    header += "            Soma {" + "\n"
    header += "                Color 1 0 0," + "\n"
    header += "                Id 7 }" + "\n"
    header += "            Color 1 0 0," + "\n"
    header += "            Id 2 }" + "\n"
    header += "        Landmark {" + "\n"
    header += "            Pia {" + "\n"
    header += "                Color 0 1 0.5," + "\n"
    header += "                Id 9 }" + "\n"
    header += "            Vessel {" + "\n"
    header += "                Color 1 0.5 0," + "\n"
    header += "                Id 10 }" + "\n"
    header += "            Barrel {" + "\n"
    header += "                aRow {" + "\n"
    header += "                    A1 {" + "\n"
    header += "                        Color 1 0.2 0.2," + "\n"
    header += "                        Id 13 }" + "\n"
    header += "                    A2 {" + "\n"
    header += "                        Color 1 0.2 0.2," + "\n"
    header += "                        Id 14 }" + "\n"
    header += "                    A3 {" + "\n"
    header += "                        Color 1 0.2 0.2," + "\n"
    header += "                        Id 15 }" + "\n"
    header += "                    A4 {" + "\n"
    header += "                        Color 1 0.2 0.2," + "\n"
    header += "                        Id 16 }" + "\n"
    header += "                    Color 1 0.2 0.2," + "\n"
    header += "                    Id 12 }" + "\n"
    header += "                bRow {" + "\n"
    header += "                    B1 {" + "\n"
    header += "                        Color 1 0.25 0.25," + "\n"
    header += "                        Id 18 }" + "\n"
    header += "                    B2 {" + "\n"
    header += "                        Color 1 0.25 0.25," + "\n"
    header += "                        Id 19 }" + "\n"
    header += "                    B3 {" + "\n"
    header += "                        Color 1 0.25 0.25," + "\n"
    header += "                        Id 20 }" + "\n"
    header += "                    B4 {" + "\n"
    header += "                        Color 1 0.25 0.25," + "\n"
    header += "                        Id 21 }" + "\n"
    header += "                    Color 1 0.25 0.25," + "\n"
    header += "                    Id 17 }" + "\n"
    header += "                cRow {" + "\n"
    header += "                    C1 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 23 }" + "\n"
    header += "                    C2 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 24 }" + "\n"
    header += "                    C3 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 25 }" + "\n"
    header += "                    C4 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 26 }" + "\n"
    header += "                    C5 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 27 }" + "\n"
    header += "                    C6 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 28 }" + "\n"
    header += "                    Color 1 0.3 0.3," + "\n"
    header += "                    Id 22 }" + "\n"
    header += "                dRow {" + "\n"
    header += "                    D1 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 30 }" + "\n"
    header += "                    D2 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 31 }" + "\n"
    header += "                    D3 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 32 }" + "\n"
    header += "                    D4 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 33 }" + "\n"
    header += "                    D5 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 34 }" + "\n"
    header += "                    D6 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 35 }" + "\n"
    header += "                    Color 1 0.35 0.35," + "\n"
    header += "                    Id 29 }" + "\n"
    header += "                eRow {" + "\n"
    header += "                    E1 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 37 }" + "\n"
    header += "                    E2 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 38 }" + "\n"
    header += "                    E3 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 39 }" + "\n"
    header += "                    E4 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 40 }" + "\n"
    header += "                    E5 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 41 }" + "\n"
    header += "                    E6 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 42 }" + "\n"
    header += "                    Color 1 0.4 0.4," + "\n"
    header += "                    Id 36 }" + "\n"
    header += "                greekRow {" + "\n"
    header += "                    Alpha {" + "\n"
    header += "                        Color 1 0.1 0.1," + "\n"
    header += "                        Id 44 }" + "\n"
    header += "                    Beta {" + "\n"
    header += "                        Color 1 0.1 0.1," + "\n"
    header += "                        Id 45 }" + "\n"
    header += "                    Gamma {" + "\n"
    header += "                        Color 1 0.1 0.1," + "\n"
    header += "                        Id 46 }" + "\n"
    header += "                    Delta {" + "\n"
    header += "                        Color 1 0.1 0.1," + "\n"
    header += "                        Id 47 }" + "\n"
    header += "                    Color 1 0.1 0.1," + "\n"
    header += "                    Id 43 }" + "\n"
    header += "                Color 0 1 0," + "\n"
    header += "                Id 11 }" + "\n"
    header += "            WhiteMatter {" + "\n"
    header += "                Color 0.5 1 0.75," + "\n"
    header += "                Id 48 }" + "\n"
    header += "            OtherBarrels {" + "\n"
    header += "                Color 1 0 1," + "\n"
    header += "                Id 49 }" + "\n"
    header += "            ZAxis {" + "\n"
    header += "                Color 0 0 0," + "\n"
    header += "                Id 50 }" + "\n"
    header += "            Color 0 1 1," + "\n"
    header += "            Id 8 }" + "\n"
    header += "        Id 0," + "\n"
    header += "        Color 0 0 0 }" + "\n"
    header += "ContentType \"HxSpatialGraph\" }" + "\n"
    header += "" + "\n"
    header += "VERTEX { float[3] VertexCoordinates } @1 " + "\n"
    header += "VERTEX {int GraphLabels } @2 " + "\n"
    header += "" + "\n"
    header += "EDGE { int[2] EdgeConnectivity } @3 " + "\n"
    header += "EDGE { int NumEdgePoints } @4 " + "\n"
    header += "EDGE { int GraphLabels } @5 " + "\n"
    header += "" + "\n"
    header += "POINT { float[3] EdgePointCoordinates } @6 " + "\n"
    header += "POINT { float Diameter } @7 " + "\n"
    
    dataIndex = []
    for i in range(len(traces)):
#        dataIndex.append(str(i + 7))
#        header += "POINT { float " + traces[i] + " } @" + str(i + 7) + "\n"
        dataIndex.append(str(i + 8))
        header += "POINT { float " + traces[i] + " } @" + str(i + 8) + "\n"
    
    for i in range(len(tVec)):
#        only write every 10th time step for visualization
#        (step size for vis. will then be 0.25ms)
        if i%10:
            continue
        stepFName = fname
        stepFName += '_'
        stepFName += '%07.3f' % tVec[i]
        stepFName += '.am'
        with open(stepFName, 'w') as outFile:
            outFile.write(header)
            
            outFile.write('\n@1 # Vertices xyz coordinates\n')
            for sec in cell.tree:
                if sec.label in axonLabels:
                    continue
                if not allPoints:
                    v1 = sec.segPts[0]
                    v2 = sec.segPts[-1]
                else:
                    v1 = sec.pts[0]
                    v2 = sec.pts[-1]
                line1 = '%.6f %.6f %.6f\n' % (v1[0], v1[1], v1[2])
                line2 = '%.6f %.6f %.6f\n' % (v2[0], v2[1], v2[2])
                outFile.write(line1)
                outFile.write(line2)
            
            outFile.write('\n@2 # Vertex Graph Label\n')
            for sec in cell.tree:
                if sec.label in axonLabels:
                    continue
                line = '%d\n' % labels2int[sec.label]
                outFile.write(line)
                outFile.write(line)
            
            outFile.write('\n@3 # Edge Identifiers\n')
            j = 0
            for sec in cell.tree:
                if sec.label in axonLabels:
                    continue
                line = '%d %d\n' % (2*j, 2*j+1)
                outFile.write(line)
                j += 1
            
            outFile.write('\n@4 # Number of Points per Edge\n')
            for sec in cell.tree:
                if sec.label in axonLabels:
                    continue
                if not allPoints:
                    line = '%d\n' % sec.nseg
                else:
                    line = '%d\n' % sec.nrOfPts
                outFile.write(line)
            
            outFile.write('\n@5 # Edge Graph Labels\n')
            for sec in cell.tree:
                if sec.label in axonLabels:
                    continue
                line = '%d\n' % labels2int[sec.label]
                outFile.write(line)
            
            outFile.write('\n@6 # Point xyz coordinates\n')
            for sec in cell.tree:
                if sec.label in axonLabels:
                    continue
                if not allPoints:
                    for pt in sec.segPts:
                        line = '%.6f %.6f %.6f\n' % (pt[0], pt[1], pt[2])
                        outFile.write(line)
                else:
                    for pt in sec.pts:
                        line = '%.6f %.6f %.6f\n' % (pt[0], pt[1], pt[2])
                        outFile.write(line)
            
            outFile.write('\n@7 # Diameter at Point\n')
            for sec in cell.tree:
                if sec.label in axonLabels:
                    continue
                for diam in sec.segDiams:
                    line = '%.6f\n' % diam
                    outFile.write(line)
            
            outFile.write('\n@8 # Vm at Point\n')
            for sec in cell.tree:
                if sec.label in axonLabels:
                    continue
                for vec in sec.recVList:
                    line = '%.6f\n' % vec[i]
                    outFile.write(line)
            
            if len(traces) > 1:
                for j in range(len(traces))[1:]:
                    var = traces[j]
                    outFile.write('\n@%d # %s at Point\n' % (j+8, var))
                    for sec in cell.tree:
                        if sec.label in axonLabels:
                            continue
                        for vec in sec.recordVars[var]:
                            line = '%.6f\n' % vec[i]
                            outFile.write(line)
            
#            outFile.write('\n@7 # Vm at Point\n')
#            for sec in cell.tree:
#                for vec in sec.recVList:
#                    line = '%.6f\n' % vec[i]
#                    outFile.write(line)
            
    
    
    
    
def write_functional_map(fname, functionalMap):
    totalNrPts = 0
    for key in functionalMap.keys():
        totalNrPts += len(functionalMap[key])
    
    header = "# AmiraMesh 3D ASCII 2.0" + "\n"
    header += "# This SpatialGraph file was created by the Neuron Registration Tool NeuroMap " + "\n"
    header += "# NeuroMap was programmed by Robert Egger," + "\n"
    header += "# Max-Planck-Florida Institute, Jupiter, Florida " + "\n"
    header += "" + "\n"
    header += "define VERTEX " + str(totalNrPts*2) + "\n"
    header += "define EDGE " + str(totalNrPts) + "\n"
    header += "define POINT " + str(totalNrPts*2) + "\n"
    header += "" + "\n"
    header += "Parameters {GraphLabels {" + "\n"
    header += "        Neuron { " + "\n"
    header += "            Dendrite {" + "\n"
    header += "                ApicalDendrite {" + "\n"
    header += "                    Color 1 0.5 0.5," + "\n"
    header += "                    Id 4 }" + "\n"
    header += "                BasalDendrite {" + "\n"
    header += "                    Color 0.8 0.4 0.4," + "\n"
    header += "                    Id 5 }" + "\n"
    header += "                Color 1 0 0," + "\n"
    header += "                Id 3 }" + "\n"
    header += "            Axon {" + "\n"
    header += "                Color 0 0 1," + "\n"
    header += "                Id 6 }" + "\n"
    header += "            Soma {" + "\n"
    header += "                Color 1 0 0," + "\n"
    header += "                Id 7 }" + "\n"
    header += "            Color 1 0 0," + "\n"
    header += "            Id 2 }" + "\n"
    header += "        Landmark {" + "\n"
    header += "            Pia {" + "\n"
    header += "                Color 0 1 0.5," + "\n"
    header += "                Id 9 }" + "\n"
    header += "            Vessel {" + "\n"
    header += "                Color 1 0.5 0," + "\n"
    header += "                Id 10 }" + "\n"
    header += "            Barrel {" + "\n"
    header += "                aRow {" + "\n"
    header += "                    A1 {" + "\n"
    header += "                        Color 1 0.2 0.2," + "\n"
    header += "                        Id 13 }" + "\n"
    header += "                    A2 {" + "\n"
    header += "                        Color 1 0.2 0.2," + "\n"
    header += "                        Id 14 }" + "\n"
    header += "                    A3 {" + "\n"
    header += "                        Color 1 0.2 0.2," + "\n"
    header += "                        Id 15 }" + "\n"
    header += "                    A4 {" + "\n"
    header += "                        Color 1 0.2 0.2," + "\n"
    header += "                        Id 16 }" + "\n"
    header += "                Color 1 0.2 0.2," + "\n"
    header += "                Id 12 }" + "\n"
    header += "                bRow {" + "\n"
    header += "                    B1 {" + "\n"
    header += "                        Color 1 0.25 0.25," + "\n"
    header += "                        Id 18 }" + "\n"
    header += "                    B2 {" + "\n"
    header += "                        Color 1 0.25 0.25," + "\n"
    header += "                        Id 19 }" + "\n"
    header += "                    B3 {" + "\n"
    header += "                        Color 1 0.25 0.25," + "\n"
    header += "                        Id 20 }" + "\n"
    header += "                    B4 {" + "\n"
    header += "                        Color 1 0.25 0.25," + "\n"
    header += "                        Id 21 }" + "\n"
    header += "                    Color 1 0.25 0.25," + "\n"
    header += "                    Id 17 }" + "\n"
    header += "                cRow {" + "\n"
    header += "                    C1 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 23 }" + "\n"
    header += "                    C2 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 24 }" + "\n"
    header += "                    C3 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 25 }" + "\n"
    header += "                    C4 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 26 }" + "\n"
    header += "                    C5 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 27 }" + "\n"
    header += "                    C6 {" + "\n"
    header += "                        Color 1 0.3 0.3," + "\n"
    header += "                        Id 28 }" + "\n"
    header += "                    Color 1 0.3 0.3," + "\n"
    header += "                    Id 22 }" + "\n"
    header += "                dRow {" + "\n"
    header += "                    D1 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 30 }" + "\n"
    header += "                    D2 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 31 }" + "\n"
    header += "                    D3 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 32 }" + "\n"
    header += "                    D4 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 33 }" + "\n"
    header += "                    D5 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 34 }" + "\n"
    header += "                    D6 {" + "\n"
    header += "                        Color 1 0.35 0.35," + "\n"
    header += "                        Id 35 }" + "\n"
    header += "                    Color 1 0.35 0.35," + "\n"
    header += "                    Id 29 }" + "\n"
    header += "                eRow {" + "\n"
    header += "                    E1 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 37 }" + "\n"
    header += "                    E2 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 38 }" + "\n"
    header += "                    E3 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 39 }" + "\n"
    header += "                    E4 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 40 }" + "\n"
    header += "                    E5 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 41 }" + "\n"
    header += "                    E6 {" + "\n"
    header += "                        Color 1 0.4 0.4," + "\n"
    header += "                        Id 42 }" + "\n"
    header += "                    Color 1 0.4 0.4," + "\n"
    header += "                    Id 36 }" + "\n"
    header += "                greekRow {" + "\n"
    header += "                    Alpha {" + "\n"
    header += "                        Color 1 0.1 0.1," + "\n"
    header += "                        Id 44 }" + "\n"
    header += "                    Beta {" + "\n"
    header += "                        Color 1 0.1 0.1," + "\n"
    header += "                        Id 45 }" + "\n"
    header += "                    Gamma {" + "\n"
    header += "                        Color 1 0.1 0.1," + "\n"
    header += "                        Id 46 }" + "\n"
    header += "                    Delta {" + "\n"
    header += "                        Color 1 0.1 0.1," + "\n"
    header += "                        Id 47 }" + "\n"
    header += "                    Color 1 0.1 0.1," + "\n"
    header += "                    Id 43 }" + "\n"
    header += "                Color 0 1 0," + "\n"
    header += "                Id 11 }" + "\n"
    header += "            WhiteMatter {" + "\n"
    header += "                Color 0.5 1 0.75," + "\n"
    header += "                Id 48 }" + "\n"
    header += "            OtherBarrels {" + "\n"
    header += "                Color 1 0 1," + "\n"
    header += "                Id 49 }" + "\n"
    header += "            ZAxis {" + "\n"
    header += "                Color 0 0 0," + "\n"
    header += "                Id 50 }" + "\n"
    header += "            Color 0 1 1," + "\n"
    header += "            Id 8 }" + "\n"
    header += "        Id 0," + "\n"
    header += "        Color 0 0 0 }" + "\n"
    header += "ContentType \"HxSpatialGraph\" }" + "\n"
    header += "" + "\n"
    header += "VERTEX { float[3] VertexCoordinates } @1 " + "\n"
    header += "VERTEX {int GraphLabels } @2 " + "\n"
    header += "" + "\n"
    header += "EDGE { int[2] EdgeConnectivity } @3 " + "\n"
    header += "EDGE { int NumEdgePoints } @4 " + "\n"
    header += "EDGE { int GraphLabels } @5 " + "\n"
    header += "" + "\n"
    header += "POINT { float[3] EdgePointCoordinates } @6 " + "\n"
    
    with open(fname, 'w') as outFile:
        outFile.write(header)
        
        outFile.write('\n@1 # Vertices xyz coordinates\n')
        for key in functionalMap.keys():
            for pair in functionalMap[key]:
                v1, v2 = pair
                line1 = '%.6f %.6f %.6f\n' % (v1[0], v1[1], v1[2])
                line2 = '%.6f %.6f %.6f\n' % (v2[0], v2[1], v2[2])
                outFile.write(line1)
                outFile.write(line2)
        
        outFile.write('\n@2 # Vertex Graph Label\n')
        for key in functionalMap.keys():
            for pair in functionalMap[key]:
                line = '3\n'
                outFile.write(line)
                outFile.write(line)
        
        outFile.write('\n@3 # Edge Identifiers\n')
        j = 0
        for key in functionalMap.keys():
            for pair in functionalMap[key]:
                line = '%d %d\n' % (2*j, 2*j+1)
                outFile.write(line)
                j += 1
        
        outFile.write('\n@4 # Number of Points per Edge\n')
        for key in functionalMap.keys():
            for pair in functionalMap[key]:
                line = '2\n'
                outFile.write(line)
        
        outFile.write('\n@5 # Edge Graph Labels\n')
        for key in functionalMap.keys():
            for pair in functionalMap[key]:
                line = '3\n'
                outFile.write(line)
        
        outFile.write('\n@6 # Point xyz coordinates\n')
        for key in functionalMap.keys():
            for pair in functionalMap[key]:
                for pt in pair:
                    line = '%.6f %.6f %.6f\n' % (pt[0], pt[1], pt[2])
                    outFile.write(line)
    
    
    
