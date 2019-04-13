import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

def get_lfp(fname, baselineTime=60000, airpuffDuration=700,\
    ISI=2500, nTrials=30, lfpTime=(-100.0, 500), lowpassFrequency=500.0):
    
    # load and filter signal
    samplingFrequency = 32000.0
    dt = 1000.0/samplingFrequency # in ms
    signal = read_labview_dat_file(fname)
    filterOrder = 3 # Daniel used 6th order with forward only filter
    b, a = butter(filterOrder, lowpassFrequency/(samplingFrequency/2), 'lowpass')
    print 'Filtering input signal with %d order Butterworth LP filter @ %.1f Hz' % (filterOrder, lowpassFrequency)
    filteredSignal = filtfilt(b, a, signal)
    
    # check length of signal
    recTime = dt*len(signal)
    expTime = baselineTime + nTrials*(airpuffDuration + ISI)
    if abs(recTime - expTime) > 0:
        diff = recTime - expTime
        print 'Something is wrong with signal. Recording time %.1f does not match stimulus setup time %.1f' % (recTime/1000, expTime/1000)
        print 'Difference is %.5f ms' % diff
    
    # compute trial-averaged LFP
    nSamplesLFP = int((lfpTime[1] - lfpTime[0])/dt) + 1
    avgSignal = np.zeros(nSamplesLFP)
    avgSignalFiltered = np.zeros(nSamplesLFP)
    individualTraces = []
    individualTracesFiltered = []
    lfpTimeAxis = np.linspace(lfpTime[0], lfpTime[1], nSamplesLFP)
    for i in range(nTrials):
        tmpT = baselineTime + i*(airpuffDuration + ISI)
        tmpT += lfpTime[0]
        stimStartIndex = int(tmpT/dt)
        stimEndIndex = stimStartIndex + nSamplesLFP
        
        s = signal[stimStartIndex:stimEndIndex]
        sFiltered = filteredSignal[stimStartIndex:stimEndIndex]
        
        avgSignal += s
        avgSignalFiltered += sFiltered
        individualTraces.append(s)
        individualTracesFiltered.append(sFiltered)
    
    avgSignal /= nTrials
    avgSignalFiltered /= nTrials
    # latency: minimum within 100ms after stimulus time
    latency = dt*np.argmin(avgSignalFiltered[(lfpTimeAxis > 0)*(lfpTimeAxis <= 100)])
    print 'Latency = %.1f ms' % latency
    
    # plot LFP individual traces and avg
    plt.figure(1)
    for i in range(nTrials):
        plt.plot(lfpTimeAxis, individualTraces[i], 'grey')
    plt.plot(lfpTimeAxis, avgSignal, 'r')
    plt.title('Raw traces')
    plt.xlabel('Time (ms)')
    
    plt.figure(2)
    for i in range(nTrials):
        plt.plot(lfpTimeAxis, individualTracesFiltered[i], 'grey')
    plt.plot(lfpTimeAxis, avgSignalFiltered, 'r')
    plt.title('Filtered traces')
    plt.xlabel('Time (ms)')
    
    plt.show()

def read_labview_dat_file(fname):
    print 'Reading LabView dat file %s' % fname
    scalingFactor = 100.0
    with open(fname, 'rb') as f:
        inHeader = True
        offset = 0
        while inHeader:
            line = f.readline()
            if 'ENDHEADER' in line:
                inHeader = False
                offset = f.tell()
        # skip one number before recording
        # number size: float32 = 4 bytes
        offset += 4
        f.seek(offset)
        # files saved with LabView in Windows in big endian format
        signal = np.fromfile(f, np.dtype('>f4'))
        return scalingFactor*signal

if __name__ == '__main__':
    fname = sys.argv[1]
    get_lfp(fname)