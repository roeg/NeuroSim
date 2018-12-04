#!/usr/bin/python

import sys
import numpy as np
import os, os.path
import matplotlib.pyplot as plt

def plot_cumulative_latency(fnames):
    for fname in fnames:
        PWLatencies, SWLatencies = load_latency_file(fname)
        PWLatencies.sort()
        PWLatencyHist = np.array(range(1, len(PWLatencies) + 1))/float(len(PWLatencies))
        plt.figure(1)
        plt.plot(PWLatencies, PWLatencyHist)
        SWLatencies.sort()
        SWLatencyHist = np.array(range(1, len(SWLatencies) + 1))/float(len(SWLatencies))
        plt.figure(2)
        plt.plot(SWLatencies, SWLatencyHist)
    
    plt.figure(1)
    plt.xlim([0, 100])
    plt.ylim([0, 1])
    plt.figure(2)
    plt.xlim([0, 100])
    plt.ylim([0, 1])
    plt.show()

def load_latency_file(fname):
    PWLatencies = []
    SWLatencies = []
    with open(fname, 'r') as latencyFile:
        lineCnt = 0
        for line in latencyFile:
            if not line:
                continue
            lineCnt += 1
            if lineCnt == 1:
                continue
            splitLine = line.strip().split('\t')
            if splitLine[1]:
                PWLatencies.append(float(splitLine[1]))
            if splitLine[2]:
                SWLatencies.append(float(splitLine[2]))
    
    return PWLatencies, SWLatencies

if __name__ == '__main__':
    fnames = []
    for i in range(1, len(sys.argv)):
        fnames.append(sys.argv[i])
    plot_cumulative_latency(fnames)
    