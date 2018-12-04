'''
@author: robert
'''

import sys
import numpy as np
import os.path
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
    
def plot_trace_frames(traceName):
    outFolder = traceName[:-4] + '_individual_frames'
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
    
    t, v = np.loadtxt(traceName, skiprows=1, unpack=True)
    stepSize = 10
    maxFrames = len(t)//stepSize
    xMin = np.min(t)
    xMax = np.max(t)
    vMin = np.min(v-5.0)
    vMax = np.max(v+5.0)
    tWhisker = [200.0, 200.0]
    
#    for i in range(801,810):
    for i in range(maxFrames):
        plotT = t[:i*stepSize]
        plotV = v[:i*stepSize]
        verts = [(200, vMin),(250, vMin),(250, vMax),(200, vMax)]
        poly = Polygon(verts, facecolor='0.9', edgecolor='0.9')
        fig, ax = plt.subplots()
        ax.add_patch(poly)
        plt.plot(plotT, plotV, 'k')
        plt.xlim([xMin, xMax])
        plt.ylim([vMin, vMax])
        frameName = 'Vm_frame%04d.png' % i
        outName = '/'.join([outFolder, frameName])
        plt.savefig(outName)
#    plt.show()
    
def plot_multiple_trace_frames(trace1Name, trace2Name, trace3Name):
    outFolder = trace1Name[:-4] + '_individual_frames'
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
    
    t1, v1 = np.loadtxt(trace1Name, skiprows=1, unpack=True)
    t2, v2 = np.loadtxt(trace2Name, skiprows=1, unpack=True)
    t3, v3 = np.loadtxt(trace3Name, skiprows=1, unpack=True)
    stepSize = 10
    maxFrames = len(t1)//stepSize
    xMin = np.min(t1)
    xMax = np.max(t1)
    vMin = np.min((np.min(v1), np.min(v2), np.min(v3))) - 5.0
    vMax = np.max((np.max(v1), np.max(v2), np.max(v3))) + 5.0
    
#    for i in range(801,810):
    for i in range(maxFrames):
        plotT = t1[:i*stepSize]
        plotV = v1[:i*stepSize]
        verts = [(200, vMin),(250, vMin),(250, vMax),(200, vMax)]
        poly = Polygon(verts, facecolor='0.9', edgecolor='0.9')
        fig, ax = plt.subplots()
        ax.add_patch(poly)
        plt.plot(plotT, plotV, 'k')
        plt.xlim([xMin, xMax])
        plt.ylim([vMin, vMax])
        frameName = 'Vm_trace1_frame%04d.png' % i
        outName = '/'.join([outFolder, frameName])
        plt.savefig(outName)
    
    for i in range(maxFrames):
        plotT = t2[:i*stepSize]
        plotV = v2[:i*stepSize]
        verts = [(200, vMin),(250, vMin),(250, vMax),(200, vMax)]
        poly = Polygon(verts, facecolor='0.9', edgecolor='0.9')
        fig, ax = plt.subplots()
        ax.add_patch(poly)
        plt.plot(t1, v1, 'lightgrey')
        plt.plot(plotT, plotV, 'k')
        plt.xlim([xMin, xMax])
        plt.ylim([vMin, vMax])
        frameName = 'Vm_trace2_frame%04d.png' % i
        outName = '/'.join([outFolder, frameName])
        plt.savefig(outName)
    
    for i in range(maxFrames):
        plotT = t3[:i*stepSize]
        plotV = v3[:i*stepSize]
        verts = [(200, vMin),(250, vMin),(250, vMax),(200, vMax)]
        poly = Polygon(verts, facecolor='0.9', edgecolor='0.9')
        fig, ax = plt.subplots()
        ax.add_patch(poly)
        plt.plot(t1, v1, 'lightgrey')
        plt.plot(t2, v2, 'lightgrey')
        plt.plot(plotT, plotV, 'k')
        plt.xlim([xMin, xMax])
        plt.ylim([vMin, vMax])
        frameName = 'Vm_trace3_frame%04d.png' % i
        outName = '/'.join([outFolder, frameName])
        plt.savefig(outName)
#    plt.show()

def plot_active_synapse_frames(activeSynapseName, traceName):
#    outFolder = traceName[:-4] + '_individual_frames'
#    if not os.path.exists(outFolder):
#        os.makedirs(outFolder)
    
    excTypes = ('L2', 'L34', 'L4py', 'L4sp', 'L4ss', 'L5st', 'L5tt',\
                'L6cc', 'L6ccinv', 'L6ct', 'VPM')
    inhTypes = ('L1','L23Trans','L45Sym','L45Peak','L56Trans',\
                'SymLocal1','SymLocal2','SymLocal3','SymLocal4','SymLocal5','SymLocal6')
    cellTypeColorMap = {'L2': 'dodgerblue', 'L34': 'blue', 'L4py': 'palegreen',\
                    'L4sp': 'green', 'L4ss': 'lime', 'L5st': 'yellow', 'L5tt': 'orange',\
                    'L6cc': 'indigo', 'L6ccinv': 'violet', 'L6ct': 'magenta', 'VPM': 'black',\
                    'INH': 'grey'}
    plotTypes = ('L34', 'L5tt', 'L6cc', 'VPM')
#    plotTypes = ('L6cc', 'VPM')
#    maxTypes = ('L6cc', 'VPM')
    maxTypes = ('L34', 'L5tt', 'L6cc', 'VPM')
#    plotTypes = ('VPM',)
    
    activeSyns = load_active_synapses(activeSynapseName)
    data = np.loadtxt(traceName, skiprows=1, unpack=True)
    t = data[0]
    v = data[1]
    stepSize = 10
    maxFrames = len(t)//stepSize
    xMin = np.min(t)
    xMax = np.max(t)
    vMin = np.min(v-5.0)
    vMax = np.max(v+5.0)
    tWhisker = [200.0, 200.0]
    
    synapseTimes = {}
    for excType in excTypes:
        synapseTimes[excType] = []
    synapseTimes['INH'] = []
    
    for synType in activeSyns.keys():
        for excType in excTypes:
            if excType in synType:
                for syn in activeSyns[synType]:
#                    synapseTimes[excType].extend(syn[3])
                    if syn[2] == 'Dendrite':
#                    if syn[2] == 'ApicalDendrite':
                        synapseTimes[excType].extend(syn[3])
        for inhType in inhTypes:
            if inhType in synType:
                continue
#                for syn in activeSyns[synType]:
#                    synapseTimes['INH'].extend(syn[3])
    
    synStepSize = 4*stepSize # 1ms
    synFrames = maxFrames//4
    tSyns = [t[i*synStepSize] for i in range(synFrames)]
    nrOfBins = len(tSyns)
    
    nMax = 0
    synapseHist = {}
#    for cellType in synapseTimes.keys():
    for cellType in maxTypes:
        hist, bins = np.histogram(synapseTimes[cellType], bins=nrOfBins, range=(xMin, xMax))
        synapseHist[cellType] = hist
        cellTypeMax = np.max(hist)
        if cellTypeMax > nMax:
            nMax = cellTypeMax
    nMax += 1
    
    for cellType in plotTypes:
        plt.plot(tSyns, synapseHist[cellType], color=cellTypeColorMap[cellType])
    plt.xlim([xMin, xMax])
    plt.ylim([0, nMax])
    plt.show()
    
##    for i in range(maxFrames-1, maxFrames):
#    for i in range(maxFrames):
#        plotT = t[:i*stepSize]
#        plotV = v[:i*stepSize]
#        synIndex = i*stepSize//synStepSize + 1
#        verts = [(200, 0),(250, 0),(250, nMax),(200, nMax)]
#        poly = Polygon(verts, facecolor='0.9', edgecolor='0.9')
#        fig, ax = plt.subplots()
#        ax.add_patch(poly)
#        for cellType in plotTypes:
#            plt.plot(tSyns[:synIndex], synapseHist[cellType][:synIndex], color=cellTypeColorMap[cellType])
#        plt.xlim([xMin, xMax])
#        plt.ylim([0, nMax])
##        frameName = 'active_syns_frame%04d.png' % i
##        frameName = 'active_VPM_L6cc_basal_syns_frame%04d.png' % i
#        frameName = 'active_VPM_basal_syns_frame%04d.png' % i
##        frameName = 'active_apical_syns_frame%04d.png' % i
#        outName = '/'.join([outFolder, frameName])
#        plt.savefig(outName)
##    plt.show()

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

if __name__ == '__main__':
    if len(sys.argv) == 2:
        traceName = sys.argv[1]
        plot_trace_frames(traceName)
    elif len(sys.argv) == 3:
        traceName = sys.argv[1]
        synapseName = sys.argv[2]
        plot_active_synapse_frames(synapseName, traceName)
    elif len(sys.argv) == 4:
        trace1Name = sys.argv[1]
        trace2Name = sys.argv[2]
        trace3Name = sys.argv[3]
        plot_multiple_trace_frames(trace1Name, trace2Name, trace3Name)
    else:
        print 'Error! Number of arguments is %d; should be 1, 2 or 3' % (len(sys.argv)-1)
    
    
    
    