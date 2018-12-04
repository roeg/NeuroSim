'''
Created on Jan 28, 2013

ongoing activity L2 neuron model

@author: robert
'''

import sys
import os, os.path
import glob
import single_cell_analyzer as sca
import single_cell_parser as scp
import numpy as np
import matplotlib.pyplot as plt

def create_spike_times_files(folder):
    '''
    automatically creates spike time files
    from all somatic recording files found below folder
    '''
    suffix = 'vm_all_traces.csv'
    fnames = []
    scan_directory(folder, fnames, suffix)
    
    for fname in fnames:
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        t = data[0]
        allSpikeTimes = {}
        for i in range(1, len(data)):
            v = data[i]
            spikeTimes = sca.simple_spike_detection(t, v)
            allSpikeTimes[i-1] = spikeTimes
        
        outName = fname[:-4] + '_spike_times.csv'
        scp.write_spike_times_file(outName, allSpikeTimes)

def create_spike_raster_plots(folder, suffix):
    '''
    load all traces, compute spike times
    and create raster plots
    '''
    fnames = []
    scan_directory(folder, fnames, suffix)
    
    print 'Creating spike raster plots from %d files' % len(fnames)
    
    tOffset = 100.0
    tStop = 345.0
    for n in range(len(fnames)):
        fname = fnames[n]
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        t = data[0]
        allSpikeTimes = []
        ax = []
        fig = plt.figure(2*n+1)
        for i in range(1, len(data)):
            fig.add_subplot(2,1,1)
            ax.append(plt.plot(t, data[i]))
            fig.add_subplot(2,1,2)
            v = data[i]
            spikeTimes = sca.simple_spike_detection(t, v)
            allSpikeTimes.append(spikeTimes)
            spikes = [i for time in spikeTimes]
            ax.append(plt.plot(spikeTimes, spikes, 'k|'))
#            plt.plot(spikeTimes, spikes, 'k|')
        fig.add_subplot(2,1,1)
        plt.xlabel('t [ms]')
        plt.ylabel('Vm [mV]')
        plt.xlim([tOffset,tStop])
        fig.add_subplot(2,1,2)
        plt.xlabel('t [ms]')
        plt.ylabel('trial nr.')
        plt.xlim([tOffset,tStop])
        
        outName = fname[:-4]
        outName += '_spike_raster_plot.pdf'
        plt.savefig(outName)
        
#        hist, bins = sca.PSTH_from_spike_times(allSpikeTimes, 1.0, 200.0, 250.0)
#        hist, bins = sca.PSTH_from_spike_times(allSpikeTimes, 1.0, 100.0, 250.0)
#        hist, bins = sca.PSTH_from_spike_times(allSpikeTimes, 5.0, 100.0, 350.0)
        hist, bins = sca.PSTH_from_spike_times(allSpikeTimes, 5.0, tOffset, tStop)
#        hist, bins = sca.PSTH_from_spike_times(allSpikeTimes, 10.0, 100.0, 250.0)
        offset = 0.5*(bins[1] - bins[0])
        fig = plt.figure(2*n+2)
#        plt.bar(bins[:-1]+offset, hist, color='b')
        plt.bar(bins[:-1], hist, color='b', width=5)
#        plt.bar(bins[:-1], hist, color='b', width=10)
        plt.xlabel('t [ms]')
        plt.ylabel('AP/stim')
        plt.xlim([tOffset,tStop])
        outName = fname[:-4]
#        outName += '_PSTH_total_10ms.pdf'
        outName += '_PSTH_total_5ms.pdf'
        plt.savefig(outName)
        scp.write_PSTH(outName[:-4]+'.csv', hist, bins)

def create_PSTH(folder, suffix, contains, outName):
    '''
    load all traces, compute spike times
    and create raster plots
    '''
    fnames = []
    scan_directory2(folder, fnames, suffix, contains)
    
    print 'Creating PSTH from %d files' % len(fnames)
    
    tOffset = 100.0
    tStop = 345.0
    allSpikeTimes = []
    for n in range(len(fnames)):
        fname = fnames[n]
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        t = data[0]
        for i in range(1, len(data)):
            v = data[i]
            spikeTimes = sca.simple_spike_detection(t, v)
            allSpikeTimes.append(spikeTimes)
    
    binWidth = 1.0
    hist, bins = sca.PSTH_from_spike_times(allSpikeTimes, binWidth, tOffset, tStop)
    offset = 0.5*(bins[1] - bins[0])
    fig = plt.figure(2*n+2)
    plt.bar(bins[:-1], hist, color='b', width=binWidth)
    plt.xlabel('t [ms]')
    plt.ylabel('AP/stim')
    plt.xlim([tOffset,tStop])
    outName += '_PSTH_total_%.1fms' % binWidth
    plt.savefig(outName+'.pdf')
    scp.write_PSTH(outName+'.csv', hist, bins)

def compute_latencies(folder, suffix, contains, outName):
    '''
    load all traces, compute spike times
    and create raster plots
    '''
    fnames = []
    scan_directory2(folder, fnames, suffix, contains)
    
    print 'Computing latencies from %d files' % len(fnames)
    
    tOffset = 100.0
    tStim = 245.0
    tWindow = 25.0
    tStop = 345.0
    allLatencies = []
    nTrials = 0
    for n in range(len(fnames)):
        fname = fnames[n]
        print 'Loading traces from file %s' % fname
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        t = data[0]
        for i in range(1, len(data)):
            nTrials += 1
            v = data[i]
            spikeTimes = sca.simple_spike_detection(t, v)
            for tSpike in spikeTimes:
                if tStim < tSpike <= tStim + tWindow:
                    allLatencies.append(tSpike-tStim)
                    break
    
    if not outName.endswith('.csv'):
        outName += '.csv'
    with open(outName, 'w') as outFile:
        line = '# Latencies in %d trials\n' % nTrials
        for latency in allLatencies:
            line += str(latency)
            line += '\n'
        outFile.write(line)
    
#    binWidth = 1.0
#    hist, bins = sca.PSTH_from_spike_times(allLatencies, binWidth, tOffset, tStop)
#    offset = 0.5*(bins[1] - bins[0])
##    fig = plt.figure(2*n+2)
##    plt.bar(bins[:-1], hist, color='b', width=binWidth)
##    plt.xlabel('t [ms]')
##    plt.ylabel('AP/stim')
##    plt.xlim([tOffset,tStop])
#    outName += '_latency_hist_%.1fms' % binWidth
##    plt.savefig(outName+'.pdf')
#    scp.write_PSTH(outName+'.csv', hist, bins)

def scan_directory(path, fnames, suffix):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory(fname, fnames, suffix)
        elif fname.endswith(suffix):
            fnames.append(fname)
        else:
            continue

def scan_directory2(path, fnames, suffix, contains):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory2(fname, fnames, suffix, contains)
        elif fname.endswith(suffix) and fname.find(contains) != -1:
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 2:
        folder = sys.argv[1]
        create_spike_times_files(folder)
    elif len(sys.argv) == 3:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        create_spike_raster_plots(folder, suffix)
    elif len(sys.argv) == 5:
        folder = sys.argv[1]
        suffix = sys.argv[2]
        contains = sys.argv[3]
        outName = sys.argv[4]
        create_PSTH(folder, suffix, contains, outName)
        #compute_latencies(folder, suffix, contains, outName)
    else:
        print 'Error! Number of arguments is %d; should be 1, 2 or 4' % (len(sys.argv)-1)
    
    
    
    