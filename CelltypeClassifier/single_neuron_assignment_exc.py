#!/usr/bin/python

import sys
import numpy as np
from scipy.stats import chi2

#===============================================================================
# define global variables (input filenames, cell types etc...)
#===============================================================================
clusterPath = 'ground_truth_exc/'
L2Name = clusterPath + 'L2_unambiguous_cluster.csv'
L3Name = clusterPath + 'L3_unambiguous_cluster.csv'
L3_4Name = clusterPath + 'L3_4_unambiguous_cluster.csv'
L4_3Name = clusterPath + 'L4_3_unambiguous_cluster.csv'
L4ssName = clusterPath + 'L4ss_unambiguous_cluster.csv'
L4spaName = clusterPath + 'L4sp-py_unambiguous_cluster.csv'
L4spbName = clusterPath + 'L4sp2_unambiguous_cluster.csv'
L4pyName = clusterPath + 'L4py_unambiguous_cluster.csv'
L5stName = clusterPath + 'L5st_unambiguous_cluster.csv'
L5ttName = clusterPath + 'L5tt_unambiguous_cluster.csv'
L6ccName = clusterPath + 'L6cc_unambiguous_cluster.csv'
L6ctName = clusterPath + 'L6ct_unambiguous_cluster.csv'
supraTypes = ['L2','L3','L3_4','L4_3']
L3_4TypeNames = ['L3','L3_4','L4_3']
granTypes = ['L4ss','L4sp-py','L4sp-ss','L4py']
L4spTypeNames = ['L4sp-py','L4sp-ss']
infraTypes = ['L5st','L5tt','L6cc','L6ct']
L6Types = ['L6cc','L6ct']
L6Border = 1400.0
#pia is at the end of DB output files...
parametersSupra = [-1,13,14]
parametersGran = [-1,1,14]
parametersInfra = [-1,0,1,2,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21]
DOF = {}
DOF['L2'] = 3
DOF['L3'] = 3
DOF['L3_4'] = 3
DOF['L4_3'] = 3
DOF['L4ss'] = 3
DOF['L4sp-py'] = 3
DOF['L4sp-ss'] = 3
DOF['L4py'] = 3
DOF['L5st'] = 21
DOF['L5tt'] = 21
DOF['L6cc'] = 21
DOF['L6ct'] = 21
#------------------------------------------------------------------------------ 

def main(fname):
    clusters = load_cluster_covariances()
    neuronValues = load_single_neuron_values_DB(fname)
    
    clusterDistances, clusterProbs = compute_cluster_probabilities(neuronValues, clusters)
    assignment = assign_celltype(neuronValues, clusterDistances, clusterProbs)
    
    assignedType = assignment[0]
    assignedDistance = assignment[1]
    assignedPercentile = assignment[2]
    assignedProb = assignment[3]
    
    ofname = fname[:-4]
    ofname += '_assigned_celltype.txt'
    with open(ofname, 'w') as outFile:
        line = '%s\t%f\t%f\n' % (assignedType, assignedProb, assignedPercentile)
        outFile.write(line)

def compute_cluster_probabilities(neuronValues, clusters):
    clusterIDs = clusters.keys()
    clusterIDs.sort()
    clusterDistances = {}
    clusterProbs = {}
    for ID in clusterIDs:
        layer = 0
        if ID in supraTypes:
            layer = 1
        elif ID in granTypes:
            layer = 2
        elif ID in infraTypes:
            layer = 3
        dist = compute_cluster_distance(neuronValues, layer, clusters[ID])
        dof = DOF[ID]
        prob = 1.0 - chi2.cdf(dist**2, dof)
        clusterDistances[ID] = dist
        clusterProbs[ID] = prob
    return clusterDistances, clusterProbs

def assign_celltype(neuronValues, clusterDistances, clusterProbs):
    clusterIDs = clusterDistances.keys()
    clusterIDs.sort()
    normClusterProbs = {}
    norm = 0.0
    for ID in clusterIDs:
        prob = clusterProbs[ID]
        norm += prob
    for ID in clusterIDs:
        prob = clusterProbs[ID]
        tmp = prob/norm
        normClusterProbs[ID] = tmp
    assignedType = 'none'
    assignedDistance = np.inf
    assignedPercentile = -1.0
    assignedProb = -1.0
    for ID in clusterIDs:
        tmpDist = clusterDistances[ID]
        tmpPercentile = clusterProbs[ID]
        tmpProb = normClusterProbs[ID]
        if tmpProb == np.nan:
            continue
        if tmpProb > assignedProb:
            assignedType = ID
            assignedDistance = tmpDist
            assignedPercentile = tmpPercentile
            assignedProb = tmpProb
    if assignedType == 'none':
        somaDepth = neuronValues[parametersSupra[0]]
        if somaDepth >= L6Border:
            clusterIDs = L6Types
        for ID in clusterIDs:
            tmpDist = clusterDistances[ID]
            if tmpDist < assignedDistance:
                assignedType = ID
                assignedDistance = tmpDist
                assignedPercentile = np.nan
                assignedProb = np.nan
    if assignedType in L3_4TypeNames:
        assignedType = 'L3-4'
    if assignedType in L4spTypeNames:
        assignedType = 'L4sp'
    return assignedType, assignedDistance, assignedPercentile, assignedProb

def compute_cluster_distance(neuronValues, layer, cluster):
    dist = 0.0
    mean = cluster[0]
    covInv = cluster[1]
    if layer == 1:
        x = []
        for i in range(len(parametersSupra)):
            ID = parametersSupra[i]
            x.append(neuronValues[ID])
        x = np.array(x)
        diff = x-mean
        dist = np.dot(diff, np.dot(covInv, diff))
    elif layer == 2:
        x = []
        for i in range(len(parametersGran)):
            ID = parametersGran[i]
            x.append(neuronValues[ID])
        x = np.array(x)
        diff = x-mean
        dist = np.dot(diff, np.dot(covInv, diff))
    elif layer == 3:
        x = []
        for i in range(len(parametersInfra)):
            ID = parametersInfra[i]
            x.append(neuronValues[ID])
        x = np.array(x)
        diff = x-mean
        dist = np.dot(diff, np.dot(covInv, diff))
    else:
        errstr = 'Invalid layer: %d' % layer
        raise RuntimeError(errstr)
    
    return np.sqrt(dist)

def load_single_neuron_values(fname):
    return np.loadtxt(fname, delimiter=',')

def load_single_neuron_values_DB(fname):
    neuronFeatures = []
    with open(fname, 'r') as neuronFile:
        lineCnt = 0
        for line in neuronFile:
            if lineCnt == 1:
                splitLine = line.split(',')
                for i in range(1, len(splitLine)-3):
                    neuronFeatures.append(float(splitLine[i].strip('\"')))
                neuronFeatures.append(float(splitLine[-1].strip('\"\n\r')))
                break
            lineCnt += 1
    return np.array(neuronFeatures)

def load_cluster_covariances():
    #===========================================================================
    # load ground truth sample for each cell type
    # and compute mean and (pseudo-)inverse covariance matrix
    #===========================================================================
    L2Cluster = np.loadtxt(L2Name, skiprows=1, delimiter=',')
    L3Cluster = np.loadtxt(L3Name, skiprows=1, delimiter=',')
    L3_4Cluster = np.loadtxt(L3_4Name, skiprows=1, delimiter=',')
    L4_3Cluster = np.loadtxt(L4_3Name, skiprows=1, delimiter=',')
    L4ssCluster = np.loadtxt(L4ssName, skiprows=1, delimiter=',')
    L4spaCluster = np.loadtxt(L4spaName, skiprows=1, delimiter=',')
    L4spbCluster = np.loadtxt(L4spbName, skiprows=1, delimiter=',')
    L4pyCluster = np.loadtxt(L4pyName, skiprows=1, delimiter=',')
    L5stCluster = np.loadtxt(L5stName, skiprows=1, delimiter=',')
    L5ttCluster = np.loadtxt(L5ttName, skiprows=1, delimiter=',')
    L6ccCluster = np.loadtxt(L6ccName, skiprows=1, delimiter=',')
    L6ctCluster = np.loadtxt(L6ctName, skiprows=1, delimiter=',')
    
    L2ClusterMean = np.mean(L2Cluster, axis=0)
    diff = L2Cluster[0]-L2ClusterMean
    L2ClusterCov = np.outer(diff, diff)
    for i in range(1, len(L2Cluster)):
        diff = L2Cluster[i]-L2ClusterMean
        L2ClusterCov += np.outer(diff, diff)
    L2ClusterCov /= len(L2Cluster) - 1
    L2ClusterCovInv = np.linalg.pinv(L2ClusterCov)
    
    L3ClusterMean = np.mean(L3Cluster, axis=0)
    diff = L3Cluster[0]-L3ClusterMean
    L3ClusterCov = np.outer(diff, diff)
    for i in range(1, len(L3Cluster)):
        diff = L3Cluster[i]-L3ClusterMean
        L3ClusterCov += np.outer(diff, diff)
    L3ClusterCov /= len(L3Cluster) - 1
    L3ClusterCovInv = np.linalg.pinv(L3ClusterCov)
    
    L3_4ClusterMean = np.mean(L3_4Cluster, axis=0)
    diff = L3_4Cluster[0]-L3_4ClusterMean
    L3_4ClusterCov = np.outer(diff, diff)
    for i in range(1, len(L3_4Cluster)):
        diff = L3_4Cluster[i]-L3_4ClusterMean
        L3_4ClusterCov += np.outer(diff, diff)
    L3_4ClusterCov /= len(L3_4Cluster) - 1
    L3_4ClusterCovInv = np.linalg.pinv(L3_4ClusterCov)
    
    L4_3ClusterMean = np.mean(L4_3Cluster, axis=0)
    diff = L4_3Cluster[0]-L4_3ClusterMean
    L4_3ClusterCov = np.outer(diff, diff)
    for i in range(1, len(L4_3Cluster)):
        diff = L4_3Cluster[i]-L4_3ClusterMean
        L4_3ClusterCov += np.outer(diff, diff)
    L4_3ClusterCov /= len(L4_3Cluster) - 1
    L4_3ClusterCovInv = np.linalg.pinv(L4_3ClusterCov)

    L4ssClusterMean = np.mean(L4ssCluster, axis=0)
    diff = L4ssCluster[0]-L4ssClusterMean
    L4ssClusterCov = np.outer(diff, diff)
    for i in range(1, len(L4ssCluster)):
        diff = L4ssCluster[i]-L4ssClusterMean
        L4ssClusterCov += np.outer(diff, diff)
    L4ssClusterCov /= len(L4ssCluster) - 1
    L4ssClusterCovInv = np.linalg.pinv(L4ssClusterCov)
    
    L4spaClusterMean = np.mean(L4spaCluster, axis=0)
    diff = L4spaCluster[0]-L4spaClusterMean
    L4spaClusterCov = np.outer(diff, diff)
    for i in range(1, len(L4spaCluster)):
        diff = L4spaCluster[i]-L4spaClusterMean
        L4spaClusterCov += np.outer(diff, diff)
    L4spaClusterCov /= len(L4spaCluster) - 1
    L4spaClusterCovInv = np.linalg.pinv(L4spaClusterCov)
    
    L4spbClusterMean = np.mean(L4spbCluster, axis=0)
    diff = L4spbCluster[0]-L4spbClusterMean
    L4spbClusterCov = np.outer(diff, diff)
    for i in range(1, len(L4spbCluster)):
        diff = L4spbCluster[i]-L4spbClusterMean
        L4spbClusterCov += np.outer(diff, diff)
    L4spbClusterCov /= len(L4spbCluster) - 1
    L4spbClusterCovInv = np.linalg.pinv(L4spbClusterCov)

    L4pyClusterMean = np.mean(L4pyCluster, axis=0)
    diff = L4pyCluster[0]-L4pyClusterMean
    L4pyClusterCov = np.outer(diff, diff)
    for i in range(1, len(L4pyCluster)):
        diff = L4pyCluster[i]-L4pyClusterMean
        L4pyClusterCov += np.outer(diff, diff)
    L4pyClusterCov /= len(L4pyCluster) - 1
    L4pyClusterCovInv = np.linalg.pinv(L4pyClusterCov)
    
    L5stClusterMean = np.mean(L5stCluster, axis=0)
    diff = L5stCluster[0]-L5stClusterMean
    L5stClusterCov = np.outer(diff, diff)
    for i in range(1, len(L5stCluster)):
        diff = L5stCluster[i]-L5stClusterMean
        L5stClusterCov += np.outer(diff, diff)
    L5stClusterCov /= len(L5stCluster) - 1
    L5stClusterCovInv = np.linalg.pinv(L5stClusterCov)
    
    L5ttClusterMean = np.mean(L5ttCluster, axis=0)
    diff = L5ttCluster[0]-L5ttClusterMean
    L5ttClusterCov = np.outer(diff, diff)
    for i in range(1, len(L5ttCluster)):
        diff = L5ttCluster[i]-L5ttClusterMean
        L5ttClusterCov += np.outer(diff, diff)
    L5ttClusterCov /= len(L5ttCluster) - 1
    L5ttClusterCovInv = np.linalg.pinv(L5ttClusterCov)
    
    L6ccClusterMean = np.mean(L6ccCluster, axis=0)
    diff = L6ccCluster[0]-L6ccClusterMean
    L6ccClusterCov = np.outer(diff, diff)
    for i in range(1, len(L6ccCluster)):
        diff = L6ccCluster[i]-L6ccClusterMean
        L6ccClusterCov += np.outer(diff, diff)
    L6ccClusterCov /= len(L6ccCluster) - 1
    L6ccClusterCovInv = np.linalg.pinv(L6ccClusterCov)
   
    L6ctClusterMean = np.mean(L6ctCluster, axis=0)
    diff = L6ctCluster[0]-L6ctClusterMean
    L6ctClusterCov = np.outer(diff, diff)
    for i in range(1, len(L6ctCluster)):
        diff = L6ctCluster[i]-L6ctClusterMean
        L6ctClusterCov += np.outer(diff, diff)
    L6ctClusterCov /= len(L6ctCluster) - 1
    L6ctClusterCovInv = np.linalg.pinv(L6ctClusterCov)
    
    #===========================================================================
    # estimate rank of covariance matrix:
    # rank = n-1 where n = number of samples in a cluster
    #===========================================================================
    DOF['L2'] = min(3, len(L2Cluster)-1)
    DOF['L3'] = min(3, len(L3Cluster)-1)
    DOF['L3_4'] = min(3, len(L3_4Cluster)-1)
    DOF['L4_3'] = min(3, len(L4_3Cluster)-1)
    DOF['L4ss'] = min(3, len(L4ssCluster)-1)
    DOF['L4sp-py'] = min(3, len(L4spaCluster)-1)
    DOF['L4sp-ss'] = min(3, len(L4spbCluster)-1)
    DOF['L4py'] = min(3, len(L4pyCluster)-1)
    DOF['L5st'] = min(21, len(L5stCluster)-1)
    DOF['L5tt'] = min(21, len(L5ttCluster)-1)
    DOF['L6cc'] = min(21, len(L6ccCluster)-1)
    DOF['L6ct'] = min(21, len(L6ctCluster)-1)
    
    clusters = {}
    clusters['L2'] = L2ClusterMean, L2ClusterCovInv
    clusters['L3'] = L3ClusterMean, L3ClusterCovInv
    clusters['L3_4'] = L3_4ClusterMean, L3_4ClusterCovInv
    clusters['L4_3'] = L4_3ClusterMean, L4_3ClusterCovInv
    clusters['L4ss'] = L4ssClusterMean, L4ssClusterCovInv
    clusters['L4sp-py'] = L4spaClusterMean, L4spaClusterCovInv
    clusters['L4sp-ss'] = L4spbClusterMean, L4spbClusterCovInv
    clusters['L4py'] = L4pyClusterMean, L4pyClusterCovInv
    clusters['L5st'] = L5stClusterMean, L5stClusterCovInv
    clusters['L5tt'] = L5ttClusterMean, L5ttClusterCovInv
    clusters['L6cc'] = L6ccClusterMean, L6ccClusterCovInv
    clusters['L6ct'] = L6ctClusterMean, L6ctClusterCovInv
    return clusters

if __name__ == "__main__":
    inputFName = sys.argv[1]
    main(inputFName)
