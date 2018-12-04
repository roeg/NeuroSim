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
    #    synNrTimingProbs[nSyn] = {}
    #    for synTiming in synNrTimingFilenames[nSyn].keys():
    #        fname = synNrTimingFilenames[nSyn][synTiming]
    #        data = np.loadtxt(fname, skiprows=1, unpack=True)
    #        PSTH = data[2]
    #        spikeProb = np.sum(PSTH[beginBin:endBin])
    #        synNrTimingProbs[nSyn][synTiming] = spikeProb
     
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
    
    numbers = np.array(range(25,350,25), dtype=np.float64)
    timings = np.array(range(2,25,1), dtype=np.float64)
    #timings = np.array(range(2,11,1), dtype=np.float64)
    numberTimingMesh = np.meshgrid(numbers,timings)
    spikeProbMesh = np.zeros_like(numberTimingMesh[0])
    for i in range(len(spikeProbMesh)-1):
        for j in range(len(spikeProbMesh[i])-1):
            number = int(numberTimingMesh[0][i][j])
            timing = numberTimingMesh[1][i][j]
            try:
                spikeProbMesh[i][j] = synNrTimingProbs[number][timing]
            except KeyError:
                spikeProbMesh[i][j] = 0.0
            ## for SuW use 0.5 as max:
            #if spikeProbMesh[i][j] > 0.5:
                #spikeProbMesh[i][j] = 0.5
    
    #with open(outName+'_grid_probs.csv', 'w') as outFile:
        #synNumbers = synNrTimingProbs.keys()
        #synNumbers.sort()
        #synTimings = synNrTimingProbs[synNumbers[0]].keys()
        #synTimings.sort()
        #header = '#syn nr.\\timing'
        #for synTiming in synTimings:
            #header += '\t'
            #header += '%.1f' % synTiming
        #header += '\n'
        #outFile.write(header)
        #for synNr in synNumbers:
            ##line = str(synNr)
            ##line = str(0.76*synNr) # PW
            #line = str(0.69*synNr) # SuW
            #for synTiming in synTimings:
                #line += '\t'
                #try:
                    #line += str(synNrTimingProbs[synNr][synTiming])
                #except KeyError:
                    #line += '-1.0'
            #line += '\n'
            #outFile.write(line)
    
    # PW
    # L3, L4ss/sp, L5tt, L6cc, VPM, L6cc+VPM
    #numberMeans = 1/0.76*np.array([74.4, 89.9, 62.2, 54.8, 54.6, 109.4])
    #numberSTD = 1/0.76*np.array([14.3, 24.8, 18.2, 12.5, 14.1, 17.7])
    #timingMeans = np.array([16.0, 20.3, 14.9, 7.0, 2.8, 4.1])
    #timingSTD = np.array([2.8, 2.3, 3.3, 3.6, 0.7, 0.6])
    
    # SuW
    # L3, L4ss/sp, L5tt, L6cc, VPM
    #numberMeans = 1/0.69*np.array([24.4, 23.0, 40.0, 37.7, 6.1])
    #numberSTD = 1/0.69*np.array([15.1, 22.2, 19.9, 12.2, 6.8])
    #timingMeans = np.array([20.5, 21.0, 18.2, 10.1, 8.1])
    #timingSTD = np.array([5.2, 6.0, 4.2, 5.1, 7.2])
    
    # E2
    # L3/4, L5tt, L6cc, VPM
    numberMeans = 1/0.68*np.array([2.2, 9.0, 19.1, 0.1])
    numberSTD = 1/0.68*np.array([9.6, 11.8, 10.2, 1.4])
    timingMeans = np.array([21.6, 20.8, 14.1, 18.3])
    timingSTD = np.array([5.3, 5.1, 7.0, 11.7])
    
    plt.figure(1)
    
    plt.errorbar(numberMeans, timingMeans, xerr=numberSTD, yerr=timingSTD, fmt='ro')
    
    #plt.pcolormesh(numberTimingMesh[0], numberTimingMesh[1], spikeProbMesh, cmap='hot')
    # PW
    #isocontours = [0.02, 0.04, 0.08, 0.15, 0.18, 0.3, 0.6, 0.9]
    # SuW
    isocontours = [0.01, 0.02, 0.04, 0.08, 0.15, 0.3, 0.6, 0.9]
    CS = plt.contour(numberTimingMesh[0], numberTimingMesh[1], spikeProbMesh, isocontours)
    plt.clabel(CS, inline=1, fontsize=10, inline_spacing=-15)
    plt.title('Simplest default with labels')
    plt.xlabel('Synapse number')
    plt.ylabel('Synapse timing (median; ms)')
    plt.xlim(0,325)
    plt.ylim(2.0,24.0)
    #plt.xticks(np.arange(25,325,25.0), [str(0.76*i) for i in range(25,325,25)]) # PW: 76% syn <500 microns
    #plt.xticks(np.arange(25,325,25.0), [str(0.69*i) for i in range(25,325,25)]) # SuW: 69% syn <500 microns
    plt.xticks(np.arange(25,325,25.0), [str(0.68*i) for i in range(25,325,25)]) # SuW: 68% syn <500 microns
    #plt.yticks(np.arange(2,23,1.0), [str(i) for i in range(2,23,1)])
    plt.yticks(np.arange(2,24,1.0), [str(i) for i in range(2,24,1)])
    cbar = plt.colorbar()
    cbar.set_label('Spike prob')
    plt.savefig(outName+'_isocontours.pdf')
    #plt.show()
    
    #numbers = np.array(range(25,325,25), dtype=np.float64)
    #timings = np.array(range(2,24,1), dtype=np.float64) # SuW
    ##timings = np.array(range(2,23,1), dtype=np.float64) # PW
    ##timings = np.array(range(2,11,1), dtype=np.float64)
    #numberTimingMesh = np.meshgrid(numbers,timings)
    #spikeProbMesh = np.zeros_like(numberTimingMesh[0])
    #for i in range(len(spikeProbMesh)):
    #    for j in range(len(spikeProbMesh[i])):
    #        number = int(numberTimingMesh[0][i][j])
    #        timing = numberTimingMesh[1][i][j]
    #        try:
    #            # invert colormap
    #            #spikeProbMesh[i][j] = 1.0 - synNrTimingProbs[number][timing]
    #            spikeProbMesh[i][j] = synNrTimingProbs[number][timing]
    #        except KeyError:
    #            spikeProbMesh[i][j] = 0.0
    #        # for SuW use 0.5 as max:
    #        ## invert colormap
    #        #if spikeProbMesh[i][j] < 0.5:
    #        # invert colormap
    #        if spikeProbMesh[i][j] > 0.5:
    #            spikeProbMesh[i][j] = 0.5
    #
    #plt.figure(1)
    ##plt.pcolormesh(numberTimingMesh[0], numberTimingMesh[1], spikeProbMesh, cmap='hot')
    #plt.imshow(spikeProbMesh, cmap='Greys', interpolation='bilinear')
    #plt.xlabel('Synapse number')
    #plt.ylabel('Synapse timing (median; ms)')
    ##plt.xlim(25,325)
    ##plt.ylim(2.0,24.0)
    ##plt.xticks([i for i in range(len(numbers))], [str(0.76*i) for i in numbers]) # PW: 76% syn <500 microns
    #plt.xticks([i for i in range(len(numbers))], [str(0.69*i) for i in numbers]) # SuW: 69% syn <500 microns
    ##plt.yticks(np.arange(2.5,23.5,1.0), [str(i) for i in range(2,23,1)])
    #plt.yticks([i for i in range(len(timings))], [str(i) for i in timings])
    #cbar = plt.colorbar()
    #cbar.set_label('Spike prob')
    ##plt.savefig(outName+'_heatmap.pdf')
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
