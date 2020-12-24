/****************************************************************************/
/*                                                                          */
/* Program:   3DCellCountAnalysis                                           */
/*                                                                          */
/* File:      pipeline.cpp                                                  */
/*                                                                          */
/* Purpose:   pipeline for analysis of cell counts with respect to          */
/*            barreloids in VPM reconstructed in 3D                         */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*            5353 Parkside Drive                                           */
/*            Jupiter, Florida 33458                                        */
/*            USA                                                           */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   22.12.2010                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "../../common/typedefs.h"
#include "../../common/amiraReader.h"
#include "../../common/inputcheckpoint.h"
#include "../../common/inputparameters.h"
#include "../../common/profile.h"
#include <gsl/gsl_rng.h>
// #include <string>

#define DEBUG

double getMeanPresynapticDistance(AmiraSpatialGraph * presynapticSG, int edgeID, std::vector< int > edgePoints);
void getMeanPoint3D(double meanPoint[3], AmiraSpatialGraph * presynapticSG, int edgeID, std::vector< int > edgePoints);
double getRASynapseProbability(double somaLocation[3], double synapseLocation[3]);

int main( int argc , char * argv[])
{
	// set up RNG
	const gsl_rng_type * T;
	gsl_rng * rng;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc(T);
	
	if(argc == 4)
	{
		const char * sgListFilename = argv[1];
		const char * avgHVCFilename = argv[2];
		const char * outputFilename = argv[3];
		
		std::vector< char * > individualSpatialGraphNames;
		std::vector< AmiraSpatialGraph * > individualSpatialGraphs;
		std::vector< double * > somaLocations;
		std::ifstream inputStream(sgListFilename);
		if(!inputStream.fail())
		{
			std::cout << "Loading SpatialGraph files from " << sgListFilename << std::endl;
			std::string currentLine;
			while(!std::getline(inputStream, currentLine).eof())
				if(currentLine.size())
				{
					char * tmpChar = new char[128];
					sscanf(currentLine.c_str(), " %s ", tmpChar);
					individualSpatialGraphNames.push_back(tmpChar);
					
					std::string sgStr(tmpChar);
					Reader * sgReader = new Reader(tmpChar, outputFilename);
					if(sgStr.find(".hoc") != std::string::npos)
					{
						std::cout << "Loading SpatialGraph file " << tmpChar << std::endl;
						sgReader->readHocFile();
						AmiraSpatialGraph * thisSG = sgReader->getSpatialGraph();
						individualSpatialGraphs.push_back(thisSG);
						PolyDataPointerType thisSoma = PolyDataPointerType::New();
						if(thisSG->extractLandmark(Soma, thisSoma))
						{
							int subId;
							double * somaLoc = new double[3];
							double pCenter[3], * weights = new double[thisSoma->GetNumberOfPoints()];
							thisSoma->GetCell(0)->GetParametricCenter(pCenter);
							thisSoma->GetCell(0)->EvaluateLocation(subId, pCenter, somaLoc, weights);
							somaLocations.push_back(somaLoc);
							std::cout << "\tSoma location @ [" << somaLoc[0] << "," << somaLoc[1] << "," << somaLoc[2] << "]" << std::endl;
							delete [] weights;
						}
						else
						{
							std::cout << "Error! SpatialGraph file " << tmpChar << " does not contain a soma!" << std::endl;
							delete sgReader;
							return 0;
						}
					}
					else
					{
						std::cout << "Error! SpatialGraph file " << tmpChar << " has to be Amira '.hoc' ascii file!" << std::endl;
						delete sgReader;
						return 0;
					}
					
					delete sgReader;
				}
		}
		
		Reader * surfReader = new Reader(avgHVCFilename, avgHVCFilename);
		PolyDataPointerType surface;
		std::string surfName(avgHVCFilename);
		if(surfName.find(".surf") != std::string::npos)
		{
			surface = surfReader->readAmiraSurfaceFile();
		}
		else
		{
			std::cout << "Error! Surface file has to be Amira '.surf' file!" << std::endl;
			delete surfReader;
			return 0;
		}
		ClosedSurface * HVCSurface = new ClosedSurface(surface);
		
		// determine bounding box of registered somata
		double sampleBoundingBox[6];
		double minX = 1e6, maxX = -1e6, minY = 1e6, maxY = -1e6, minZ = 1e6, maxZ = -1e6;
		for(int i = 0; i < somaLocations.size(); ++i)
		{
			if(somaLocations[i][0] < minX)
			{
				minX = somaLocations[i][0];
			}
			if(somaLocations[i][0] > maxX)
			{
				maxX = somaLocations[i][0];
			}
			if(somaLocations[i][1] < minY)
			{
				minY = somaLocations[i][1];
			}
			if(somaLocations[i][1] > maxY)
			{
				maxY = somaLocations[i][1];
			}
			if(somaLocations[i][2] < minZ)
			{
				minZ = somaLocations[i][2];
			}
			if(somaLocations[i][2] > maxZ)
			{
				maxZ = somaLocations[i][2];
			}
		}
		sampleBoundingBox[0] = minX;
		sampleBoundingBox[1] = maxX;
		sampleBoundingBox[2] = minY;
		sampleBoundingBox[3] = maxY;
		sampleBoundingBox[4] = minZ;
		sampleBoundingBox[5] = maxZ;
// 		// sample bounding box all registered:
// 		double sampleBoundingBox[] = {-523.999,598.7,-368.841,244.961,-43.3674,274.065};
#ifdef DEBUG
		std::cout << "sampleBoundingBox = [" << sampleBoundingBox[0] << "," << sampleBoundingBox[1] << "," << sampleBoundingBox[2] << "," << sampleBoundingBox[3] << "," << sampleBoundingBox[4] << "," << sampleBoundingBox[5] << "]" << std::endl;
#endif
		// place grid for virtual spheres into bounding box
		// however, in the ventral direction the grid should be bounded
		// by the ventral avg HVC extent (sample is limited to dorsal half of HVC)
		double HVCBounds[6];
		HVCSurface->ptr()->GetBounds(HVCBounds);
// 		HVCSurface->ptr()->Print(std::cout);
		const double sphereGridSpacing = 20;
		const double sphereRadius = 100;
		int sphereGrid[3];
// 		sphereGrid[0] = int((sampleBoundingBox[1] - sampleBoundingBox[0])/sphereGridSpacing + 1);
// 		sphereGrid[1] = int((sampleBoundingBox[3] - sampleBoundingBox[2])/sphereGridSpacing + 1);
// // 		sphereGrid[2] = int((sampleBoundingBox[5] - sampleBoundingBox[4])/sphereGridSpacing + 1);
// 		sphereGrid[2] = int((sampleBoundingBox[5] - HVCBounds[4])/sphereGridSpacing + 1);
		sphereGrid[0] = int((HVCBounds[1] - HVCBounds[0])/sphereGridSpacing + 1);
		sphereGrid[1] = int((HVCBounds[3] - HVCBounds[2])/sphereGridSpacing + 1);
		sphereGrid[2] = int((HVCBounds[5] - HVCBounds[4])/sphereGridSpacing + 1);
#ifdef DEBUG
		std::cout << "HVCBounds = [" << HVCBounds[0] << "," << HVCBounds[1] << "," << HVCBounds[2] << "," << HVCBounds[3] << "," << HVCBounds[4] << "," << HVCBounds[5] << "]" << std::endl;
		std::cout << "sphereGrid = [" << sphereGrid[0] << "," << sphereGrid[1] << "," << sphereGrid[2] << "]" << std::endl;
#endif
		// iterate over cells and spheres:
		// for each cell, store identity of innervated spheres
		// for each sphere, store soma locations of innervating cells
		std::vector< std::vector< int > > individualSphereInnervations;
		std::vector< double > presynapticInnervationDistances;
		std::vector< double > presynapticSomaDistances;
		
// #ifdef DEBUG
// 		std::cout << "Processing morphology " << individualSpatialGraphNames[n] << std::endl;
// #endif
		int totalSpheresInsideHVC = 0;
		ImageDataPointerType innervatedSpheresGrid = ImageDataPointerType::New();
		double gridOrigin[3];
		gridOrigin[0] = sampleBoundingBox[0];
		gridOrigin[1] = sampleBoundingBox[2];
		gridOrigin[2] = HVCBounds[4];
		innervatedSpheresGrid->SetOrigin(gridOrigin);
		innervatedSpheresGrid->SetSpacing(sphereGridSpacing, sphereGridSpacing, sphereGridSpacing);
		innervatedSpheresGrid->SetDimensions(sphereGrid);
#ifdef DEBUG
		int axonInsideCount = 0;
#endif
		for(int i = 0; i < sphereGrid[0]; ++i)
			for(int j = 0; j < sphereGrid[1]; ++j)
				for(int k = 0; k < sphereGrid[2]; ++k)
				{
					int index = k + (sphereGrid[2]*(j + sphereGrid[1]*i));
					double sphereCenter[3];
					sphereCenter[0] = HVCBounds[0] + i*sphereGridSpacing;
					sphereCenter[1] = HVCBounds[2] + j*sphereGridSpacing;
					sphereCenter[2] = HVCBounds[4] + k*sphereGridSpacing;
					int ijk[3];
					ijk[0] = i;
					ijk[1] = j;
					ijk[2] = k;
					if(!HVCSurface->isPointInsideSurface(sphereCenter))
					{
						double * densVal = static_cast< double * >(innervatedSpheresGrid->GetScalarPointer(ijk));
						*densVal = -1;
						continue;
					}
					++totalSpheresInsideHVC;
					std::vector< int > innervatingNeuronIDs;
					for(int n = 0; n < individualSpatialGraphs.size(); ++n)
					{
						std::map< int, std::vector< int > > pointsInSphere;
						AmiraSpatialGraph * neuronMorphology = individualSpatialGraphs[n];
						bool axonInsideSphere = false;
						//main loop
						std::vector< Edge * >::const_iterator edgeIt;
#ifdef DEBUG
						int axonCnt = 1;
#endif
						int edgeCnt = 0;
						for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt, ++edgeCnt)
						{
							if((*edgeIt)->label != Axon)
							{
								continue;
							}
							int edgePtCnt = 0;
							std::list< double * >::const_iterator edgePtIt;
							for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++ edgePtIt, ++edgePtCnt)
							{
								double * currentPt = *edgePtIt;
								double currentDist = sqrt(vtkMath::Distance2BetweenPoints(currentPt, sphereCenter));
								if(currentDist < sphereRadius)
								{
									axonInsideSphere = true;
#ifdef DEBUG
									++axonInsideCount;
#endif
									if(pointsInSphere.find(edgeCnt) == pointsInSphere.end())
									{
										std::vector< int > tmpVec;
										tmpVec.push_back(edgePtCnt);
										pointsInSphere[edgeCnt] = tmpVec;
									}
									else
									{
										pointsInSphere[edgeCnt].push_back(edgePtCnt);
									}
								}
							}
							if(axonInsideSphere)
							{
								break;
							}
#ifdef DEBUG
							++axonCnt;
#endif
						}
						
						if(axonInsideSphere)
						{
							innervatingNeuronIDs.push_back(n);
							std::map< int, std::vector< int > >::const_iterator pointsInSphereIt;
							for(pointsInSphereIt = pointsInSphere.begin(); pointsInSphereIt != pointsInSphere.end(); ++pointsInSphereIt)
							{
								int edgeID = pointsInSphereIt->first;
								std::vector< int > edgePoints = pointsInSphereIt->second;
								double meanDist = getMeanPresynapticDistance(neuronMorphology, edgeID, edgePoints);
								double meanPoint[3];
								getMeanPoint3D(meanPoint, neuronMorphology, edgeID, edgePoints);
								double prob = getRASynapseProbability(somaLocations[n], meanPoint);
								if(gsl_rng_uniform(rng) < prob)
								{
									presynapticInnervationDistances.push_back(meanDist);
									presynapticSomaDistances.push_back(sqrt(vtkMath::Distance2BetweenPoints(somaLocations[n], meanPoint)));
								}
							}
						}
						else
						{
							// nothing to be done
						}
					}
					if(innervatingNeuronIDs.size())
					{
						individualSphereInnervations.push_back(innervatingNeuronIDs);
					}
				}
#ifdef DEBUG
		std::cout << "axonInsideCount = " << axonInsideCount << std::endl;
#endif
		std::cout << "Total number of spheres inside avg. HVC contour = " << totalSpheresInsideHVC << std::endl;
		
		std::string outName(outputFilename);
		outName += "_data_presynaptic_distances_completeHVC.csv";
		std::cout << "outName = " << outName.c_str() << std::endl;

		std::ofstream outFile(outName.c_str());
		outFile << "Path length distance\tEuclidean distance" << std::endl;
		for(int i = 0; i < presynapticInnervationDistances.size(); ++i)
		{
			outFile << presynapticInnervationDistances[i] << "\t" << presynapticSomaDistances[i] << std::endl;
		}
		outFile.close();
	}
	
	// find sphere close to specified 3D location (for visualization)
	if(argc == 7)
	{
		const char * sgListFilename = argv[1];
		const char * avgHVCFilename = argv[2];
		const char * outputFilename = argv[3];
		double sphereLocation[3];
		sphereLocation[0] = atof(argv[4]);
		sphereLocation[1] = atof(argv[5]);
		sphereLocation[2] = atof(argv[6]);
		
		std::vector< char *  > individualSpatialGraphNames;
		std::vector< AmiraSpatialGraph *   > individualSpatialGraphs;
		std::vector< double * > somaLocations;
		std::ifstream inputStream(sgListFilename);
		if(!inputStream.fail())
		{
			std::cout << "Loading SpatialGraph files from " << sgListFilename << std::endl;
			std::string currentLine;
			while(!std::getline(inputStream, currentLine).eof())
				if(currentLine.size())
				{
					char * tmpChar = new char[128];
					sscanf(currentLine.c_str(), " %s ", tmpChar);
					individualSpatialGraphNames.push_back(tmpChar);
					
					std::string sgStr(tmpChar);
					Reader * sgReader = new Reader(tmpChar, tmpChar);
					if(sgStr.find(".am") != std::string::npos)
					{
						std::cout << "Loading SpatialGraph file " << tmpChar << std::endl;
						sgReader->readSpatialGraphFile(1);
						AmiraSpatialGraph * thisSG = sgReader->getSpatialGraph();
						individualSpatialGraphs.push_back(thisSG);
						PolyDataPointerType thisSoma = PolyDataPointerType::New();
						if(thisSG->extractLandmark(Soma, thisSoma))
						{
							int subId;
							double * somaLoc = new double[3];
							double pCenter[3], * weights = new double[thisSoma->GetNumberOfPoints()];
							thisSoma->GetCell(0)->GetParametricCenter(pCenter);
							thisSoma->GetCell(0)->EvaluateLocation(subId, pCenter, somaLoc, weights);
							somaLocations.push_back(somaLoc);
							std::cout << "\tSoma location @ [" << somaLoc[0] << "," << somaLoc[1] << "," << somaLoc[2] << "]" << std::endl;
							delete [] weights;
						}
						else
						{
							std::cout << "Error! SpatialGraph file " << tmpChar << " does not contain a soma!" << std::endl;
							delete sgReader;
							return 0;
						}
					}
					else
					{
						std::cout << "Error! SpatialGraph file " << tmpChar << " has to be Amira '.am' ascii file!" << std::endl;
						delete sgReader;
						return 0;
					}
					
					delete sgReader;
				}
		}
		
		Reader * surfReader = new Reader(avgHVCFilename, avgHVCFilename);
		PolyDataPointerType surface;
		std::string surfName(avgHVCFilename);
		if(surfName.find(".surf") != std::string::npos)
		{
			surface = surfReader->readAmiraSurfaceFile();
		}
		else
		{
			std::cout << "Error! Surface file has to be Amira '.surf' file!" << std::endl;
			delete surfReader;
			return 0;
		}
		ClosedSurface * HVCSurface = new ClosedSurface(surface);
		
		// determine bounding box of registered somata
		double sampleBoundingBox[6];
		double minX = 1e6, maxX = -1e6, minY = 1e6, maxY = -1e6, minZ = 1e6, maxZ = -1e6;
		for(int i = 0; i < somaLocations.size(); ++i)
		{
			if(somaLocations[i][0] < minX)
			{
				minX = somaLocations[i][0];
			}
			if(somaLocations[i][0] > maxX)
			{
				maxX = somaLocations[i][0];
			}
			if(somaLocations[i][1] < minY)
			{
				minY = somaLocations[i][1];
			}
			if(somaLocations[i][1] > maxY)
			{
				maxY = somaLocations[i][1];
			}
			if(somaLocations[i][2] < minZ)
			{
				minZ = somaLocations[i][2];
			}
			if(somaLocations[i][2] > maxZ)
			{
				maxZ = somaLocations[i][2];
			}
		}
		sampleBoundingBox[0] = minX;
		sampleBoundingBox[1] = maxX;
		sampleBoundingBox[2] = minY;
		sampleBoundingBox[3] = maxY;
		sampleBoundingBox[4] = minZ;
		sampleBoundingBox[5] = maxZ;
// 		// sample bounding box all registered:
// 		double sampleBoundingBox[] = {-523.999,598.7,-368.841,244.961,-43.3674,274.065};
#ifdef DEBUG
		std::cout << "sampleBoundingBox = [" << sampleBoundingBox[0] << "," << sampleBoundingBox[1] << "," << sampleBoundingBox[2] << "," << sampleBoundingBox[3] << "," << sampleBoundingBox[4] << "," << sampleBoundingBox[5] << "]" << std::endl;
#endif
		// place grid for virtual spheres into bounding box
		// however, in the ventral direction the grid should be bounded
		// by the ventral avg HVC extent (sample is limited to dorsal half of HVC)
		double HVCBounds[6];
		HVCSurface->ptr()->GetBounds(HVCBounds);
// 		HVCSurface->ptr()->Print(std::cout);
		const double sphereGridSpacing = 20;
		const double sphereRadius = 100;
		int sphereGrid[3];
// 		sphereGrid[0] = int((sampleBoundingBox[1] - sampleBoundingBox[0])/sphereGridSpacing + 1);
// 		sphereGrid[1] = int((sampleBoundingBox[3] - sampleBoundingBox[2])/sphereGridSpacing + 1);
// // 		sphereGrid[2] = int((sampleBoundingBox[5] - sampleBoundingBox[4])/sphereGridSpacing + 1);
// 		sphereGrid[2] = int((sampleBoundingBox[5] - HVCBounds[4])/sphereGridSpacing + 1);
		sphereGrid[0] = int((HVCBounds[1] - HVCBounds[0])/sphereGridSpacing + 1);
		sphereGrid[1] = int((HVCBounds[3] - HVCBounds[2])/sphereGridSpacing + 1);
		sphereGrid[2] = int((HVCBounds[5] - HVCBounds[4])/sphereGridSpacing + 1);
#ifdef DEBUG
		std::cout << "HVCBounds = [" << HVCBounds[0] << "," << HVCBounds[1] << "," << HVCBounds[2] << "," << HVCBounds[3] << "," << HVCBounds[4] << "," << HVCBounds[5] << "]" << std::endl;
		std::cout << "sphereGrid = [" << sphereGrid[0] << "," << sphereGrid[1] << "," << sphereGrid[2] << "]" << std::endl;
#endif
		// iterate over cells and spheres:
		// for each cell, store identity of innervated spheres
		// for each sphere, store soma locations of innervating cells
		std::vector< std::vector< int > > individualSphereInnervations;
		std::vector< std::vector< double > > individualSphereLocations;
// #ifdef DEBUG
// 		std::cout << "Processing morphology " << individualSpatialGraphNames[n] << std::endl;
// #endif
		int totalSpheresInsideHVC = 0;
#ifdef DEBUG
		int axonInsideCount = 0;
#endif
		for(int i = 0; i < sphereGrid[0]; ++i)
			for(int j = 0; j < sphereGrid[1]; ++j)
				for(int k = 0; k < sphereGrid[2]; ++k)
				{
					int index = k + (sphereGrid[2]*(j + sphereGrid[1]*i));
					double sphereCenter[3];
					sphereCenter[0] = HVCBounds[0] + i*sphereGridSpacing;
					sphereCenter[1] = HVCBounds[2] + j*sphereGridSpacing;
					sphereCenter[2] = HVCBounds[4] + k*sphereGridSpacing;
					int ijk[3];
					ijk[0] = i;
					ijk[1] = j;
					ijk[2] = k;
					if(!HVCSurface->isPointInsideSurface(sphereCenter))
					{
						continue;
					}
					++totalSpheresInsideHVC;
					std::vector< int > innervatingNeuronIDs;
					for(int n = 0; n < individualSpatialGraphs.size(); ++n)
					{
						AmiraSpatialGraph * neuronMorphology = individualSpatialGraphs[n];
						bool axonInsideSphere = false;
						//main loop
						std::vector< Edge * >::const_iterator edgeIt;
#ifdef DEBUG
						int axonCnt = 1;
#endif
						for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
						{
							if((*edgeIt)->label != Axon)
							{
								continue;
							}
							std::list< double * >::const_iterator edgePtIt;
							for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
							{
								double * currentPt = *edgePtIt;
								double currentDist = sqrt(vtkMath::Distance2BetweenPoints(currentPt, sphereCenter));
								if(currentDist < sphereRadius)
								{
									axonInsideSphere = true;
#ifdef DEBUG
									++axonInsideCount;
#endif
									break;
								}
							}
#ifdef DEBUG
							++axonCnt;
#endif
						}
						
						if(axonInsideSphere)
						{
							innervatingNeuronIDs.push_back(n);
							break;
						}
						else
						{
							// nothing to be done
						}
					}
					if(innervatingNeuronIDs.size())
					{
						individualSphereInnervations.push_back(innervatingNeuronIDs);
						std::vector< double > tmpSphereLocation;
						tmpSphereLocation.push_back(sphereCenter[0]);
						tmpSphereLocation.push_back(sphereCenter[1]);
						tmpSphereLocation.push_back(sphereCenter[2]);
						individualSphereLocations.push_back(tmpSphereLocation);
					}
				}
#ifdef DEBUG
		std::cout << "axonInsideCount = " << axonInsideCount << std::endl;
#endif
		std::cout << "Total number of spheres inside avg. HVC contour = " << totalSpheresInsideHVC << std::endl;
		int usedInnervatedSpheres = 0;
		bool foundTargetSphere = 0;
		double targetThreshold = sphereGridSpacing;
		double minDist = 1.0e6;
		for(int i = 0; i < individualSphereInnervations.size(); ++i)
		{
			int clusterSize = individualSphereInnervations[i].size();
			double tmpRatio = double(clusterSize)/double(individualSpatialGraphs.size());
			if(clusterSize < 2)
			{
				continue;
			}
			++usedInnervatedSpheres;
			
			double tmpSphereLocation[3];
			tmpSphereLocation[0] = individualSphereLocations[i][0];
			tmpSphereLocation[1] = individualSphereLocations[i][1];
			tmpSphereLocation[2] = individualSphereLocations[i][2];
			double sphereTargetDist = sqrt(vtkMath::Distance2BetweenPoints(tmpSphereLocation, sphereLocation));
			if(sphereTargetDist < minDist)
			{
				minDist = sphereTargetDist;
			}
			if(sphereTargetDist > targetThreshold)
			{
				continue;
			}
			
			foundTargetSphere = 1;
			double totalDist = 0;
			PointsPointerType clusterNeurons = PointsPointerType::New();
			clusterNeurons->SetDataTypeToFloat();
			clusterNeurons->Allocate(1);
			for(int j = 0; j < clusterSize; ++j)
			{
				int neuronID1 = individualSphereInnervations[i][j];
				double * somaLoc1 = somaLocations[neuronID1];
				clusterNeurons->InsertNextPoint(somaLoc1);
				for(int k = j + 1; k < clusterSize; ++k)
				{
					int neuronID2 = individualSphereInnervations[i][k];
					double * somaLoc2 = somaLocations[neuronID2];
					double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(somaLoc1, somaLoc2));
					totalDist += tmpDist;
				}
			}
			int pairwiseComparisons = clusterSize*(clusterSize - 1)/2;
			double avgDist = totalDist/pairwiseComparisons;
			
			std::cout << "*********************" << std::endl;
			std::cout << "Found close location!" << std::endl;
			std::cout << "distance to target = " << sphereTargetDist << std::endl;
			std::cout << "Target location @ [" << sphereLocation[0] << "," << sphereLocation[1] << "," << sphereLocation[2] << "]" << std::endl;
			std::cout << "Sphere location @ [" << tmpSphereLocation[0] << "," << tmpSphereLocation[1] << "," << tmpSphereLocation[2] << "]" << std::endl;
			std::cout << "cluster size = " << clusterSize << std::endl;
			std::cout << "average pairwise dist = " << avgDist << std::endl;
			
			std::string outStr(outputFilename);
			char suffix[64];
			sprintf(suffix, "_sphere_%d_connected_cluster", usedInnervatedSpheres);
			outStr += suffix;
			Reader * landmarkWriter = new Reader(outStr.c_str(), outStr.c_str());
			landmarkWriter->writeLandmarkFile(clusterNeurons);
			delete landmarkWriter;
		}
		
		if(!foundTargetSphere)
		{
			std::cout << "Could not find target sphere, consider increasing the distance threshold (currently " << targetThreshold << " microns)" << std::endl;
			std::cout << "minDist = " << minDist << std::endl;
		}
		
		std::cout << "Total number of spheres used (i.e. >= 2 innervating neurons) = " << usedInnervatedSpheres << std::endl;
		
	}
	
	// find sphere location close to target avg. pairwise distance
	if(argc == 6)
	{
		const char * sgListFilename = argv[1];
		const char * avgHVCFilename = argv[2];
		const char * outputFilename = argv[3];
		double minMLCoordinate = atof(argv[4]);
		double clusterTargetDistance = atof(argv[5]);
		
		std::vector< char *  > individualSpatialGraphNames;
		std::vector< AmiraSpatialGraph *   > individualSpatialGraphs;
		std::vector< double * > somaLocations;
		std::ifstream inputStream(sgListFilename);
		if(!inputStream.fail())
		{
			std::cout << "Loading SpatialGraph files from " << sgListFilename << std::endl;
			std::string currentLine;
			while(!std::getline(inputStream, currentLine).eof())
				if(currentLine.size())
				{
					char * tmpChar = new char[128];
					sscanf(currentLine.c_str(), " %s ", tmpChar);
					individualSpatialGraphNames.push_back(tmpChar);
					
					std::string sgStr(tmpChar);
					Reader * sgReader = new Reader(tmpChar, tmpChar);
					if(sgStr.find(".am") != std::string::npos)
					{
						std::cout << "Loading SpatialGraph file " << tmpChar << std::endl;
						sgReader->readSpatialGraphFile(1);
						AmiraSpatialGraph * thisSG = sgReader->getSpatialGraph();
						individualSpatialGraphs.push_back(thisSG);
						PolyDataPointerType thisSoma = PolyDataPointerType::New();
						if(thisSG->extractLandmark(Soma, thisSoma))
						{
							int subId;
							double * somaLoc = new double[3];
							double pCenter[3], * weights = new double[thisSoma->GetNumberOfPoints()];
							thisSoma->GetCell(0)->GetParametricCenter(pCenter);
							thisSoma->GetCell(0)->EvaluateLocation(subId, pCenter, somaLoc, weights);
							somaLocations.push_back(somaLoc);
							std::cout << "\tSoma location @ [" << somaLoc[0] << "," << somaLoc[1] << "," << somaLoc[2] << "]" << std::endl;
							delete [] weights;
						}
						else
						{
							std::cout << "Error! SpatialGraph file " << tmpChar << " does not contain a soma!" << std::endl;
							delete sgReader;
							return 0;
						}
					}
					else
					{
						std::cout << "Error! SpatialGraph file " << tmpChar << " has to be Amira '.am' ascii file!" << std::endl;
						delete sgReader;
						return 0;
					}
					
					delete sgReader;
				}
		}
		
		Reader * surfReader = new Reader(avgHVCFilename, avgHVCFilename);
		PolyDataPointerType surface;
		std::string surfName(avgHVCFilename);
		if(surfName.find(".surf") != std::string::npos)
		{
			surface = surfReader->readAmiraSurfaceFile();
		}
		else
		{
			std::cout << "Error! Surface file has to be Amira '.surf' file!" << std::endl;
			delete surfReader;
			return 0;
		}
		ClosedSurface * HVCSurface = new ClosedSurface(surface);
		
		// determine bounding box of registered somata
		double sampleBoundingBox[6];
		double minX = 1e6, maxX = -1e6, minY = 1e6, maxY = -1e6, minZ = 1e6, maxZ = -1e6;
		for(int i = 0; i < somaLocations.size(); ++i)
		{
			if(somaLocations[i][0] < minX)
			{
				minX = somaLocations[i][0];
			}
			if(somaLocations[i][0] > maxX)
			{
				maxX = somaLocations[i][0];
			}
			if(somaLocations[i][1] < minY)
			{
				minY = somaLocations[i][1];
			}
			if(somaLocations[i][1] > maxY)
			{
				maxY = somaLocations[i][1];
			}
			if(somaLocations[i][2] < minZ)
			{
				minZ = somaLocations[i][2];
			}
			if(somaLocations[i][2] > maxZ)
			{
				maxZ = somaLocations[i][2];
			}
		}
		sampleBoundingBox[0] = minX;
		sampleBoundingBox[1] = maxX;
		sampleBoundingBox[2] = minY;
		sampleBoundingBox[3] = maxY;
		sampleBoundingBox[4] = minZ;
		sampleBoundingBox[5] = maxZ;
// 		// sample bounding box all registered:
// 		double sampleBoundingBox[] = {-523.999,598.7,-368.841,244.961,-43.3674,274.065};
#ifdef DEBUG
		std::cout << "sampleBoundingBox = [" << sampleBoundingBox[0] << "," << sampleBoundingBox[1] << "," << sampleBoundingBox[2] << "," << sampleBoundingBox[3] << "," << sampleBoundingBox[4] << "," << sampleBoundingBox[5] << "]" << std::endl;
#endif
		// place grid for virtual spheres into bounding box
		// however, in the ventral direction the grid should be bounded
		// by the ventral avg HVC extent (sample is limited to dorsal half of HVC)
		double HVCBounds[6];
		HVCSurface->ptr()->GetBounds(HVCBounds);
// 		HVCSurface->ptr()->Print(std::cout);
		const double sphereGridSpacing = 20;
		const double sphereRadius = 100;
		int sphereGrid[3];
// 		sphereGrid[0] = int((sampleBoundingBox[1] - sampleBoundingBox[0])/sphereGridSpacing + 1);
// 		sphereGrid[1] = int((sampleBoundingBox[3] - sampleBoundingBox[2])/sphereGridSpacing + 1);
// // 		sphereGrid[2] = int((sampleBoundingBox[5] - sampleBoundingBox[4])/sphereGridSpacing + 1);
// 		sphereGrid[2] = int((sampleBoundingBox[5] - HVCBounds[4])/sphereGridSpacing + 1);
		sphereGrid[0] = int((HVCBounds[1] - HVCBounds[0])/sphereGridSpacing + 1);
		sphereGrid[1] = int((HVCBounds[3] - HVCBounds[2])/sphereGridSpacing + 1);
		sphereGrid[2] = int((HVCBounds[5] - HVCBounds[4])/sphereGridSpacing + 1);
#ifdef DEBUG
		std::cout << "HVCBounds = [" << HVCBounds[0] << "," << HVCBounds[1] << "," << HVCBounds[2] << "," << HVCBounds[3] << "," << HVCBounds[4] << "," << HVCBounds[5] << "]" << std::endl;
		std::cout << "sphereGrid = [" << sphereGrid[0] << "," << sphereGrid[1] << "," << sphereGrid[2] << "]" << std::endl;
#endif
		// iterate over cells and spheres:
		// for each cell, store identity of innervated spheres
		// for each sphere, store soma locations of innervating cells
		std::vector< std::vector< int > > individualSphereInnervations;
		std::vector< std::vector< double > > individualSphereLocations;
		
// #ifdef DEBUG
// 		std::cout << "Processing morphology " << individualSpatialGraphNames[n] << std::endl;
// #endif
		int totalSpheresInsideHVC = 0;
		ImageDataPointerType innervatedSpheresGrid = ImageDataPointerType::New();
		double gridOrigin[3];
		gridOrigin[0] = sampleBoundingBox[0];
		gridOrigin[1] = sampleBoundingBox[2];
		gridOrigin[2] = HVCBounds[4];
		innervatedSpheresGrid->SetOrigin(gridOrigin);
		innervatedSpheresGrid->SetSpacing(sphereGridSpacing, sphereGridSpacing, sphereGridSpacing);
		innervatedSpheresGrid->SetDimensions(sphereGrid);
#ifdef DEBUG
		int axonInsideCount = 0;
#endif
		for(int i = 0; i < sphereGrid[0]; ++i)
			for(int j = 0; j < sphereGrid[1]; ++j)
				for(int k = 0; k < sphereGrid[2]; ++k)
				{
					int index = k + (sphereGrid[2]*(j + sphereGrid[1]*i));
					double sphereCenter[3];
					sphereCenter[0] = HVCBounds[0] + i*sphereGridSpacing;
					sphereCenter[1] = HVCBounds[2] + j*sphereGridSpacing;
					sphereCenter[2] = HVCBounds[4] + k*sphereGridSpacing;
					int ijk[3];
					ijk[0] = i;
					ijk[1] = j;
					ijk[2] = k;
					if(!HVCSurface->isPointInsideSurface(sphereCenter))
					{
						double * densVal = static_cast< double * >(innervatedSpheresGrid->GetScalarPointer(ijk));
						*densVal = -1;
						continue;
					}
					++totalSpheresInsideHVC;
					std::vector< int > innervatingNeuronIDs;
					for(int n = 0; n < individualSpatialGraphs.size(); ++n)
					{
						AmiraSpatialGraph * neuronMorphology = individualSpatialGraphs[n];
						bool axonInsideSphere = false;
						//main loop
						std::vector< Edge * >::const_iterator edgeIt;
#ifdef DEBUG
						int axonCnt = 1;
#endif
						for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
						{
							if((*edgeIt)->label != Axon)
							{
								continue;
							}
							std::list< double * >::const_iterator edgePtIt;
							for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++ edgePtIt)
							{
								double * currentPt = *edgePtIt;
								double currentDist = sqrt(vtkMath::Distance2BetweenPoints(currentPt, sphereCenter));
								if(currentDist < sphereRadius)
								{
									axonInsideSphere = true;
#ifdef DEBUG
									++axonInsideCount;
#endif
									break;
								}
							}
							if(axonInsideSphere)
							{
								break;
							}
#ifdef DEBUG
							++axonCnt;
#endif
						}
						
						if(axonInsideSphere)
						{
							innervatingNeuronIDs.push_back(n);
						}
						else
						{
							// nothing to be done
						}
					}
					if(innervatingNeuronIDs.size())
					{
						individualSphereInnervations.push_back(innervatingNeuronIDs);
						std::vector< double > tmpSphereLocation;
						tmpSphereLocation.push_back(sphereCenter[0]);
						tmpSphereLocation.push_back(sphereCenter[1]);
						tmpSphereLocation.push_back(sphereCenter[2]);
						individualSphereLocations.push_back(tmpSphereLocation);
					}
				}
#ifdef DEBUG
		std::cout << "axonInsideCount = " << axonInsideCount << std::endl;
#endif
		std::cout << "Total number of spheres inside avg. HVC contour = " << totalSpheresInsideHVC << std::endl;
		int usedInnervatedSpheres = 0;
// 		bool foundTargetSphere = 0;
// 		double targetThreshold = sphereGridSpacing;
		double minDist = 1.0e6;
		int closestLocationID = -1;
		double closestAvgPairwiseDist = -1;
		for(int i = 0; i < individualSphereInnervations.size(); ++i)
		{
			int clusterSize = individualSphereInnervations[i].size();
			if(clusterSize < 2)
			{
				continue;
			}
			
			double tmpSphereLocation[3];
			tmpSphereLocation[0] = individualSphereLocations[i][0];
			tmpSphereLocation[1] = individualSphereLocations[i][1];
			tmpSphereLocation[2] = individualSphereLocations[i][2];
			
			double totalDist = 0;
			for(int j = 0; j < clusterSize; ++j)
			{
				int neuronID1 = individualSphereInnervations[i][j];
				double * somaLoc1 = somaLocations[neuronID1];
				for(int k = j + 1; k < clusterSize; ++k)
				{
					int neuronID2 = individualSphereInnervations[i][k];
					double * somaLoc2 = somaLocations[neuronID2];
					double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(somaLoc1, somaLoc2));
					totalDist += tmpDist;
				}
			}
			int pairwiseComparisons = clusterSize*(clusterSize - 1)/2;
			double avgDist = totalDist/pairwiseComparisons;
			double diff = abs(clusterTargetDistance - avgDist);
			if(diff < minDist && tmpSphereLocation[0] > minMLCoordinate)
			{
				minDist = diff;
				closestLocationID = i;
				closestAvgPairwiseDist = avgDist;
			}
		}
		
		if(closestLocationID == -1)
		{
			std::cout << "Could not find cluster where sphere is located at x > " << minMLCoordinate << std::endl;
		}
		else
		{
			int clusterSize = individualSphereInnervations[closestLocationID].size();
			double tmpSphereLocation[3];
			tmpSphereLocation[0] = individualSphereLocations[closestLocationID][0];
			tmpSphereLocation[1] = individualSphereLocations[closestLocationID][1];
			tmpSphereLocation[2] = individualSphereLocations[closestLocationID][2];
			double totalDist = 0;
			PointsPointerType clusterNeurons = PointsPointerType::New();
			clusterNeurons->SetDataTypeToFloat();
			clusterNeurons->Allocate(1);
			for(int j = 0; j < clusterSize; ++j)
			{
				int neuronID1 = individualSphereInnervations[closestLocationID][j];
				double * somaLoc1 = somaLocations[neuronID1];
				clusterNeurons->InsertNextPoint(somaLoc1);
				for(int k = j + 1; k < clusterSize; ++k)
				{
					int neuronID2 = individualSphereInnervations[closestLocationID][k];
					double * somaLoc2 = somaLocations[neuronID2];
					double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(somaLoc1, somaLoc2));
					totalDist += tmpDist;
				}
			}
			int pairwiseComparisons = clusterSize*(clusterSize - 1)/2;
			double avgDist = totalDist/pairwiseComparisons;
			
			std::cout << "*********************" << std::endl;
			std::cout << "Found close location!" << std::endl;
			std::cout << "distance to target = " << minDist << std::endl;
			std::cout << "Target avg dist = " << clusterTargetDistance << std::endl;
			std::cout << "Sphere location @ [" << tmpSphereLocation[0] << "," << tmpSphereLocation[1] << "," << tmpSphereLocation[2] << "]" << std::endl;
			std::cout << "cluster size = " << clusterSize << std::endl;
			std::cout << "average pairwise dist = " << avgDist << std::endl;
			
			std::string outStr(outputFilename);
			char suffix[64];
			sprintf(suffix, "_sphere_%d_connected_cluster_avgDist_%.1f", closestLocationID, avgDist);
			outStr += suffix;
			Reader * landmarkWriter = new Reader(outStr.c_str(), outStr.c_str());
			landmarkWriter->writeLandmarkFile(clusterNeurons);
			delete landmarkWriter;
		}
		
// 		if(!foundTargetSphere)
// 		{
// 			std::cout << "Could not find target sphere, consider increasing the distance threshold (currently " << targetThreshold << " microns)" << std::endl;
// 			std::cout << "minDist = " << minDist << std::endl;
// 		}
		
		std::cout << "Total number of spheres used (i.e. >= 2 innervating neurons) = " << usedInnervatedSpheres << std::endl;
		
	}
	
	return 0;
}


double getMeanPresynapticDistance(AmiraSpatialGraph* presynapticSG, int edgeID, std::vector< int > edgePoints)
{
	Edge * tmpEdge = presynapticSG->edgesPointer()->at(edgeID);
	double meanDist = 0.0;
	double nrOfPoints = edgePoints.size();
	for(int i = 0; i < nrOfPoints; ++i)
	{
		int targetEdgePt = edgePoints[i];
		int edgePtCnt = 0;
		double tmpDist = 0.0;
		std::list< double * >::const_iterator edgePtIt1, edgePtIt2;
		edgePtIt1 = tmpEdge->edgePointCoordinates.begin();
		edgePtIt2 = tmpEdge->edgePointCoordinates.begin();
		++edgePtIt2;
		for( ; edgePtIt2 != tmpEdge->edgePointCoordinates.end() && edgePtCnt < targetEdgePt; ++edgePtIt1, ++edgePtIt2, ++edgePtCnt)
		{
			double * pt1 = *edgePtIt1;
			double * pt2 = *edgePtIt2;
			tmpDist += sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
		}
		tmpDist += presynapticSG->totalSegmentLength(tmpEdge->fatherID);
		meanDist += tmpDist;
	}
	meanDist /= nrOfPoints;
	
	return meanDist;
}

void getMeanPoint3D(double meanPoint[3], AmiraSpatialGraph* presynapticSG, int edgeID, std::vector< int > edgePoints)
{
	Edge * tmpEdge = presynapticSG->edgesPointer()->at(edgeID);
	meanPoint[0] = 0.0;
	meanPoint[1] = 0.0;
	meanPoint[2] = 0.0;
	double nrOfPoints = edgePoints.size();
	for(int i = 0; i < nrOfPoints; ++i)
	{
		int targetEdgePt = edgePoints[i];
		std::list< double * >::const_iterator edgePtIt = tmpEdge->edgePointCoordinates.begin();
		for(int j = 0; edgePtIt != tmpEdge->edgePointCoordinates.end() && j < targetEdgePt; ++edgePtIt, ++j)
		{
			// just incrementing...
		}
		double * pt = *edgePtIt;
		for(int j = 0; j < 3; ++j)
		{
			meanPoint[j] += pt[j];
		}
	}
	for(int j = 0; j < 3; ++j)
	{
		meanPoint[j] /= nrOfPoints;
	}
}

double getRASynapseProbability(double somaLocation[3], double synapseLocation[3])
{
	double distance = sqrt(vtkMath::Distance2BetweenPoints(somaLocation, synapseLocation));
	return 1.0/(1 + 30.72*exp(-0.04577*distance));
}



