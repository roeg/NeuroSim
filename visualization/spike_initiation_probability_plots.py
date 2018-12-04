'''

batch fitting of sigmoid curve to spiking probability
curves. Caution: no way to check for p=0 artifacts due
to empty bins in trial histograms
=> manual check of fit results necessary!!!

@author: robert
'''

import sys
import os, os.path
import glob
import single_cell_analyzer as sca
import single_cell_parser as scp
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def automated_analysis(folder, contains, suffix, outName):
    '''
    batch fitting of sigmoid curve to spiking probability
    curves. Caution: no way to check for p=0 artifacts due
    to empty bins in trial histograms
    => manual check of fit results necessary!!!
    '''
    probHistNames = []
    scan_directory2(folder, probHistNames, suffix, contains)
        
    print 'Batch fitting %d spike probability files...' % len(probHistNames)
    probHistFiles = {}
    trialHistFiles = {}
    sigmoidFitParam = {}
    gaussFitParam = {}
    for fname in probHistNames:
        data = np.loadtxt(fname, skiprows=1, unpack=True)
        synchWindow = fname.split('_')[-6]
        windowOffset = fname.split('_')[-4]
        if not probHistFiles.has_key(synchWindow):
            probHistFiles[synchWindow] = {}
            trialHistFiles[synchWindow] = {}
            sigmoidFitParam[synchWindow] = {}
            gaussFitParam[synchWindow] = {}
        if probHistFiles[synchWindow].has_key(windowOffset):
            errstr = 'Duplicated window/offset file %s' % fname
            raise RuntimeError(errstr)
        probHistFiles[synchWindow][windowOffset] = (data[2], data[5])
        trialHistFiles[synchWindow][windowOffset] = (data[2], data[3], data[4])
    
    figCount = 1
    for window in probHistFiles:
        fig = plt.figure(figCount)
        for offset in probHistFiles[window]:
            x, y = probHistFiles[window][offset]
            try:
                p0 = [0.0, 10.0]
                pOpt, pCov = curve_fit(sigmoid, x, y, p0)
            except:
                print 'Error fitting sigmoid window %s / offset %s' % (window, offset)
                print 'Setting default (bad) fit values'
                pOpt = (0.0,1000.0)
            sigmoidFitParam[window][offset] = pOpt
            labelStr = 'w: %s; o: %s' % (window, offset)
            plt.plot(x, y, 'o', label=labelStr)
        plt.gca().set_color_cycle(None)
        for offset in probHistFiles[window]:
            x, y = probHistFiles[window][offset]
            pOpt = sigmoidFitParam[window][offset]
            xFit = np.linspace(np.min(x), np.max(x), 200)
            yFit = sigmoid(xFit, pOpt[0], pOpt[1])
            plt.plot(xFit, yFit, '-')
        plt.legend()
        figName = outName + '_prob_sigmoid_fit_synchWindow_' + window + '.pdf'
        plt.savefig(figName)
        figCount += 1
    
    ax = plt.gca()
    for window in trialHistFiles:
        fig = plt.figure(figCount)
        for offset in trialHistFiles[window]:
            x, y1, y2 = trialHistFiles[window][offset]
#            p0 = [50.0, 300.0, 500.0]
            p0 = [0.0, 300.0, 80.0]
            try:
                pOpt1, pCov1 = curve_fit(gauss, x, y1, p0)
            except:
                print 'Error fitting gauss1 window %s / offset %s' % (window, offset)
                print 'Setting default (bad) fit values'
                pOpt1 = (0.0, 1.0e6, 80.0)
            try:
                pOpt2, pCov2 = curve_fit(gauss, x, y2, p0)
            except:
                print 'Error fitting gauss2 window %s / offset %s' % (window, offset)
                print 'Setting default (bad) fit values'
                pOpt2 = (0.0, 1.0e6, 80.0)
            gaussFitParam[window][offset] = pOpt1, pOpt2
            labelStr = 'w: %s; o: %s' % (window, offset)
            color = next(ax._get_lines.color_cycle)
            plt.plot(x, y1, 'o', label=labelStr, color=color)
            plt.plot(x, y2, 'o', color=color)
            xFit = np.linspace(np.min(x), np.max(x), 200)
            yFit1 = gauss(xFit, pOpt1[0], pOpt1[1], pOpt1[2])
            yFit2 = gauss(xFit, pOpt2[0], pOpt2[1], pOpt2[2])
#            color = next(ax._get_lines.color_cycle)
            plt.plot(xFit, yFit1, '-', color=color)
            plt.plot(xFit, yFit2, '-', color=color)
#        plt.gca().set_color_cycle(None)
#        for offset in probHistFiles[window]:
#            x, y1, y2 = trialHistFiles[window][offset]
#            pOpt1, pOpt2 = gaussFitParam[window][offset]
#            xFit = np.linspace(np.min(x), np.max(x), 200)
#            yFit1 = gauss(xFit, pOpt1[0], pOpt1[1], pOpt1[2])
#            yFit2 = gauss(xFit, pOpt2[0], pOpt2[1], pOpt2[2])
#            color = next(ax._get_lines.color_cycle)
#            plt.plot(xFit, yFit1, '-', color=color)
#            plt.plot(xFit, yFit2, '-', color=color)
        plt.legend()
        figName = outName + '_trials_gauss_fit_synchWindow_' + window + '.pdf'
        plt.savefig(figName)
        figCount += 1
    
    outParamName = outName + '_all_fit_parameters.csv'
    with open(outParamName, 'w') as outFile:
        header = 'synchrony window\twindow offset\tsigmoid center\tsigmoid width\tgauss center1\tgauss width1\tgauss center2\tgauss width2\tseparation\n'
        outFile.write(header)
        for window in sigmoidFitParam:
            for offset in sigmoidFitParam[window]:
                line = window
                line += '\t'
                line += offset
                line += '\t'
                line += str(sigmoidFitParam[window][offset][0])
                line += '\t'
                line += str(sigmoidFitParam[window][offset][1])
                line += '\t'
                center1 = gaussFitParam[window][offset][0][0]
                widthsquared1 = gaussFitParam[window][offset][0][1]
                center2 = gaussFitParam[window][offset][1][0]
                widthsquared2 = gaussFitParam[window][offset][1][1]
                separation = np.abs(center1-center2)/np.sqrt(widthsquared1 + widthsquared2)
                line += str(center1)
                line += '\t'
                line += str(np.sqrt(widthsquared1))
                line += '\t'
                line += str(center2)
                line += '\t'
                line += str(np.sqrt(widthsquared2))
                line += '\t'
                line += str(separation)
                line += '\n'
                outFile.write(line)

def probability_slope_diagram(fname):
    '''
    2D diagram of probability slope
    as a function of (window, offset)
    '''
    data = np.loadtxt(fname, skiprows=1, unpack=True)
    windows = np.array(range(1,27,1), dtype=np.float64)
    offsets = np.array(range(0,26,1), dtype=np.float64)
    windowOffsetMesh = np.meshgrid(windows,offsets)
    slopeMesh = np.zeros_like(windowOffsetMesh[0])
    for i in range(len(slopeMesh)):
        for j in range(len(slopeMesh[i])):
            tmpWindow = windowOffsetMesh[0][i][j]
            tmpOffset = windowOffsetMesh[1][i][j]
            for n in range(len(data[0])):
                if np.abs(tmpWindow-data[0][n])<1e-6 and np.abs(tmpOffset-data[1][n])<1e-6:
                    slopeMesh[i][j] = 1.0/data[3][n]
    
    plt.figure(1)
    plt.pcolormesh(windowOffsetMesh[0], windowOffsetMesh[1], slopeMesh, cmap='hot')
    plt.xlabel('Synchrony window (ms)')
    plt.ylabel('Window offset (ms)')
    plt.xlim(1,26)
    plt.ylim(0,25)
    plt.xticks(np.arange(1.5,26.5,1.0), [str(i) for i in range(1,26,1)])
    plt.yticks(np.arange(0.5,25.5,1.0), [str(i) for i in range(0,25,1)])
    cbar = plt.colorbar()
    cbar.set_label('Slope (1/syn)')
    plt.show()

def sigmoid(x, center, width):
    return 1.0/(1.0 + np.exp(-1.0/width*(x - center)))

def gauss(x, center, widthsquared, amp):
    return amp*np.exp(-0.5/widthsquared*(x - center)**2)

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
    if len(sys.argv) == 5:
        folder = sys.argv[1]
        contains = sys.argv[2]
        suffix = sys.argv[3]
        outName = sys.argv[4]
        automated_analysis(folder, contains, suffix, outName)
    elif len(sys.argv) == 2:
        fname = sys.argv[1]
        probability_slope_diagram(fname)
    else:
        print 'Error! Number of arguments is %d; should be 1 or 4' % (len(sys.argv)-1)
    
    
    
    