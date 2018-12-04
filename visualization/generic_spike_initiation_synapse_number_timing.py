import sys
import os, os.path
import glob
import numpy as np
import matplotlib.pyplot as plt
import single_cell_analyzer as sca
import single_cell_parser as scp


def PSTH_timing_number_analysis(folder, suffix, tBegin, window, outName):
    fnames = []
    scan_directory(folder, fnames, suffix)
    
    synNrTimingFilenames = {}
    for fname in fnames:
        basePath = fname.split('/')[1]
        nSyn = int(basePath.split('_')[1])
        synTiming = float(basePath.split('_')[5])
        if not synNrTimingFilenames.has_key(nSyn):
            synNrTimingFilenames[nSyn] = {}
        synNrTimingFilenames[nSyn][synTiming] = fname
    
    #synNrTimingProbs = {}
    #binWidth = 5.0
    #offset = 145.0 # 245 - 100
    #beginBin = int((offset+tBegin)/binWidth)
    #endBin = beginBin + int(window/binWidth) + 1
    #for nSyn in synNrTimingFilenames.keys():
        #synNrTimingProbs[nSyn] = {}
        #for synTiming in synNrTimingFilenames[nSyn].keys():
            #fname = synNrTimingFilenames[nSyn][synTiming]
            #data = np.loadtxt(fname, skiprows=1, unpack=True)
            #PSTH = data[2]
            #spikeProb = np.sum(PSTH[beginBin:endBin])
            #synNrTimingProbs[nSyn][synTiming] = spikeProb
    
    #synNrTimingProbs = {}
    #offset = 245.0
    #for nSyn in synNrTimingFilenames.keys():
    #    synNrTimingProbs[nSyn] = {}
    #    for synTiming in synNrTimingFilenames[nSyn].keys():
    #        fname = synNrTimingFilenames[nSyn][synTiming]
    #        print 'Analyzing traces in file %s' % fname
    #        vmData = np.loadtxt(fname, skiprows=1, unpack=True)
    #        t = vmData[0]
    #        nrOfTrials = len(vmData) - 1
    #        spikeTrials = 0.0
    #        for i in range(nrOfTrials):
    #            v = vmData[i+1]
    #            trialSpikeTimes = sca.simple_spike_detection(t, v)
    #            for tSpike in trialSpikeTimes:
    #                if tBegin <= tSpike - offset < tBegin + window:
    #                    spikeTrials += 1
    #                    break
    #        spikeProb = spikeTrials/nrOfTrials
    #        synNrTimingProbs[nSyn][synTiming] = spikeProb
    #        print 'Window: %.1f - %.1f ms: spike prob = %.2f' % (tBegin, tBegin+window, spikeProb)
    
    synNrTimingProbs = {}
    offset = 245.0
    for nSyn in synNrTimingFilenames.keys():
        synNrTimingProbs[nSyn] = {}
        for synTiming in synNrTimingFilenames[nSyn].keys():
            fname = synNrTimingFilenames[nSyn][synTiming]
            print 'Analyzing spike times in file %s' % fname
            trialSpikeTimes = scp.read_spike_times_file(fname)
            nrOfTrials = len(trialSpikeTimes.keys())
            spikeTrials = 0.0
            for trial in trialSpikeTimes.keys():
                for tSpike in trialSpikeTimes[trial]:
                    if tBegin <= tSpike - offset < tBegin + window:
                        spikeTrials += 1
                        break
            spikeProb = spikeTrials/nrOfTrials
            synNrTimingProbs[nSyn][synTiming] = spikeProb
            print 'Window: %.1f - %.1f ms: spike prob = %.2f' % (tBegin, tBegin+window, spikeProb)
            
    numbers = np.array(range(50,400,25), dtype=np.float64)
    timings = np.array(range(2,11,1), dtype=np.float64)
    numberTimingMesh = np.meshgrid(numbers,timings)
    spikeProbMesh = np.zeros_like(numberTimingMesh[0])
    for i in range(len(spikeProbMesh)-1):
        for j in range(len(spikeProbMesh[i])-1):
            number = int(numberTimingMesh[0][i][j])
            timing = numberTimingMesh[1][i][j]
            spikeProbMesh[i][j] = synNrTimingProbs[number][timing]
    
    with open(outName+'_grid_probs.csv', 'w') as outFile:
        synNumbers = synNrTimingProbs.keys()
        synNumbers.sort()
        synTimings = synNrTimingProbs[synNumbers[0]].keys()
        synTimings.sort()
        header = '#syn nr.\\timing'
        for synTiming in synTimings:
            header += '\t'
            header += '%.1f' % synTiming
        header += '\n'
        outFile.write(header)
        for synNr in synNumbers:
            line = str(synNr)
            for synTiming in synTimings:
                line += '\t'
                line += str(synNrTimingProbs[synNr][synTiming])
            line += '\n'
            outFile.write(line)
    
    plt.figure(1)
    plt.pcolormesh(numberTimingMesh[0], numberTimingMesh[1], spikeProbMesh, cmap='hot')
    plt.xlabel('Synapse number')
    plt.ylabel('Synapse timing (median; ms)')
    plt.xlim(50,350)
    plt.ylim(2.0,10.0)
    plt.xticks(np.arange(62.5,362.5,25.0), [str(0.76*i) for i in range(50,350,25)]) # PW: 76% syn <500 microns
    plt.yticks(np.arange(2.5,10.5,1.0), [str(i) for i in range(2,10,1)])
    cbar = plt.colorbar()
    cbar.set_label('Spike prob')
    plt.savefig(outName+'_heatmap.pdf')
    #plt.show()

def scan_directory(path, fnames, suffix):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory(fname, fnames, suffix)
        elif fname.endswith(suffix):
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 6:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        tBegin = float(sys.argv[3])
        window = float(sys.argv[4])
        outName = sys.argv[5]
        PSTH_timing_number_analysis(folder, suffix, tBegin, window, outName)
    else:
        print 'Error! Number of arguments is %d; should be 5' % (len(sys.argv)-1)