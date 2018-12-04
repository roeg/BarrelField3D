#!/usr/bin/python

import sys
import numpy as np

def crop_landmarks(landmarkName, bbox):
    print '*************'
    print 'cropping landmark file %s' % landmarkName
    print '*************'
    
    eps = 1E-06
    insideLandmarks = []
    outsideLandmarks = 0
    readData = 0
    with open(landmarkName, 'r') as inputFile:
        for line in inputFile:
            if line:
                if not readData and line[:2] == '@1':
                    readData = 1
                    continue
                if readData:
                    splitLine = line.strip().split(' ')
                    if len(splitLine) == 3:
                        pt = []
                        for i in range(3):
                            pt.append(float(splitLine[i]))
                        inside = 1
                        for i in range(3):
                            if bbox[2*i] + eps > pt[i] or bbox[2*i+1] - eps < pt[i]:
                                inside = 0
                                break
                        if not inside:
                            outsideLandmarks += 1
                        else:
                            insideLandmarks.append(pt)
    
    print '%d landmarks inside bounding box' % len(insideLandmarks)
    print '%d landmarks outside bounding box' % outsideLandmarks
    
    outName = landmarkName[:-14] + '_cropped.landmarkAscii'
    print '*************'
    print 'writing cropped landmark file %s' % outName
    print '*************'
    with open(outName, 'w') as outputFile:
        header = '# AmiraMesh 3D ASCII 2.0\n\n\n'
        header += 'define Markers %d\n\n' % (len(insideLandmarks)+2)
        header += 'Parameters {\n'
        header += '    NumSets 1,\n'
        header += '    ContentType \"LandmarkSet\"\n'
        header += '}\n\n'
        header += 'Markers { float[3] Coordinates } @1\n\n'
        header += '# Landmark set cropped to bounding box %.3f %.3f %.3f %.3f %.3f %.3f\n' % (bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5])
        header += '# Two landmarks corresponding to bounding box added\n'
        header += '# Data section follows\n'
        header += '@1\n'
        outputFile.write(header)
        bboxPt1 = str(bbox[0]) + ' ' + str(bbox[2]) + ' ' + str(bbox[4]) + '\n'
        bboxPt2 = str(bbox[1]) + ' ' + str(bbox[3]) + ' ' + str(bbox[5]) + '\n'
        outputFile.write(bboxPt1)
        outputFile.write(bboxPt2)
        for pt in insideLandmarks:
            line = str(pt[0]) + ' ' + str(pt[1]) + ' ' + str(pt[2]) + '\n'
            outputFile.write(line)

if __name__ == '__main__':
    landmarkName = sys.argv[1]
    bbox = []
    for i in range(6):
        bbox.append(float(sys.argv[2+i]))
    else:
        crop_landmarks(landmarkName, bbox)