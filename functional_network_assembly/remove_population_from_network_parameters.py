#!/usr/bin/python

import sys
import single_cell_parser as scp
import os, os.path
import glob

def remove_population(folderName, contains, suffix, population):
    
    fnames = []
    get_files_in_directory(folderName, fnames, contains, suffix)
    for tmpName in fnames:
        print '*************'
        print 'Removing population %s from file %s' % (population, tmpName)
        print '*************'
        nwParamCluster = scp.build_parameters(tmpName)
        nwParamCluster.network.pop(population)
        
        insertIndex = tmpName.find('_UpState')
        insertIndex += 8
        clusterOutFileName = tmpName[:insertIndex]
        clusterOutFileName += '_'
        clusterOutFileName += population
        clusterOutFileName += '_inact'
        clusterOutFileName += tmpName[insertIndex:]
        nwParamCluster.save(clusterOutFileName)

def get_files_in_directory(path, fnames, contains, suffix):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            continue
        elif fname.endswith(suffix) and fname.find(contains) > -1:
            fnames.append(fname)
        else:
            continue

if __name__ == '__main__':
    if len(sys.argv) == 5:
        folderName = sys.argv[1]
        contains = sys.argv[2]
        suffix = sys.argv[3]
        population = sys.argv[4]
        remove_population(folderName, contains, suffix, population)
    else:
        print 'parameters: [folderName] [contains] [suffix] [population]'
    