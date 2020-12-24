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

int main( int argc , char * argv[])
{
	if(argc == 4)
	{
		const char * sgListFilename = argv[1];
		const char * avgHVCFilename = argv[2];
		const char * outputFilename = argv[3];
		
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
					Reader * sgReader = new Reader(tmpChar, outputFilename);
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
		HVCSurface->ptr()->Print(std::cout);
		const double sphereGridSpacing = 20;
		const double sphereRadius = 100;
		int sphereGrid[3];
		sphereGrid[0] = int((sampleBoundingBox[1] - sampleBoundingBox[0])/sphereGridSpacing + 1);
		sphereGrid[1] = int((sampleBoundingBox[3] - sampleBoundingBox[2])/sphereGridSpacing + 1);
// 		sphereGrid[2] = int((sampleBoundingBox[5] - sampleBoundingBox[4])/sphereGridSpacing + 1);
		sphereGrid[2] = int((sampleBoundingBox[5] - HVCBounds[4])/sphereGridSpacing + 1);
#ifdef DEBUG
		std::cout << "HVCBounds = [" << HVCBounds[0] << "," << HVCBounds[1] << "," << HVCBounds[2] << "," << HVCBounds[3] << "," << HVCBounds[4] << "," << HVCBounds[5] << "]" << std::endl;
		std::cout << "sphereGrid = [" << sphereGrid[0] << "," << sphereGrid[1] << "," << sphereGrid[2] << "]" << std::endl;
#endif
		// iterate over cells and spheres:
		// for each cell, store identity of innervated spheres
		// for each sphere, store soma locations of innervating cells
		std::vector< ImageDataPointerType > individualNeuronInnervations;
		std::vector< std::vector< double * > > individualSphereInnervations;
		std::vector< std::vector< double * > > individualSphereNonInnervations;
		std::vector< std::vector< double > > individualSphereInnervationDistances;
		std::vector< std::vector< double > > individualSphereNonInnervationDistances;
		for(int i = 0; i < sphereGrid[0]; ++i)
			for(int j = 0; j < sphereGrid[1]; ++j)
				for(int k = 0; k < sphereGrid[2]; ++k)
				{
					std::vector< double * > tmpVec;
					individualSphereInnervations.push_back(tmpVec);
					individualSphereNonInnervations.push_back(tmpVec);
					std::vector< double > tmpVec2;
					individualSphereInnervationDistances.push_back(tmpVec2);
					individualSphereNonInnervationDistances.push_back(tmpVec2);
				}
		int totalSpheres = individualSphereInnervations.size();
#ifdef DEBUG
		std::cout << "Total number of spheres (on grid) = " << totalSpheres << std::endl;
#endif
		
		for(int n = 0; n < individualSpatialGraphs.size(); ++n)
		{
#ifdef DEBUG
			std::cout << "Processing morphology " << individualSpatialGraphNames[n] << std::endl;
#endif
			int totalSpheresInsideHVC = 0;
			AmiraSpatialGraph * neuronMorphology = individualSpatialGraphs[n];
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
						sphereCenter[0] = gridOrigin[0] + i*sphereGridSpacing;
						sphereCenter[1] = gridOrigin[1] + j*sphereGridSpacing;
						sphereCenter[2] = gridOrigin[2] + k*sphereGridSpacing;
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
							double * densVal = static_cast< double * >(innervatedSpheresGrid->GetScalarPointer(ijk));
							*densVal = 1;
							double * somaDirection = new double[3];
							vtkMath::Subtract(somaLocations[n], sphereCenter, somaDirection);
							double somaDirection2D[2];
							somaDirection2D[0] = somaDirection[0];
							somaDirection2D[1] = somaDirection[1];
							double somaDist = vtkMath::Normalize(somaDirection);
							double somaDist2D = vtkMath::Normalize(somaDirection2D);
							if(index >= totalSpheres)
							{
								std::cout << "Error: index = " << index << " - total nr. of spheres = " << totalSpheres << std::endl;
							}
							else
							{
								individualSphereInnervations[index].push_back(somaDirection);
								individualSphereInnervationDistances[index].push_back(somaDist2D);
							}
						}
						else
						{
							double * densVal = static_cast< double * >(innervatedSpheresGrid->GetScalarPointer(ijk));
							*densVal = 0;
							double * somaDirection = new double[3];
							vtkMath::Subtract(somaLocations[n], sphereCenter, somaDirection);
							double somaDirection2D[2];
							somaDirection2D[0] = somaDirection[0];
							somaDirection2D[1] = somaDirection[1];
							double somaDist = vtkMath::Normalize(somaDirection);
							double somaDist2D = vtkMath::Normalize(somaDirection2D);
							if(index >= totalSpheres)
							{
								std::cout << "Error: index = " << index << " - total nr. of spheres = " << totalSpheres << std::endl;
							}
							else
							{
								individualSphereNonInnervations[index].push_back(somaDirection);
								individualSphereNonInnervationDistances[index].push_back(somaDist2D);
							}
						}
					}
			individualNeuronInnervations.push_back(innervatedSpheresGrid);
#ifdef DEBUG
			std::cout << "axonInsideCount = " << axonInsideCount << std::endl;
#endif
			std::cout << "Total number of spheres inside avg. HVC contour = " << totalSpheresInsideHVC << std::endl;
		}
		
		// output: per-cell 3D pattern of innervated spheres as AmiraMesh
// 		for(int n = 0; n < individualSpatialGraphs.size(); ++n)
// 		{
// 			std::string cellOutputName(outputFilename);
// 			cellOutputName += "_";
// 			cellOutputName += individualSpatialGraphNames[n];
// 			Reader * innervationMeshWriter = new Reader(cellOutputName.c_str(), cellOutputName.c_str());
// 			innervationMeshWriter->writeScalarField(individualNeuronInnervations[n]);
// 			delete innervationMeshWriter;
// 		}
		// output: per sphere soma locations and average as AmiraMesh
// 		double totalAvgAngle = 0;
// 		int totalNrSpheres = 0;
// 		std::vector< double > spheresAngleVec;
// 		ImageDataPointerType spheresAverageGrid = ImageDataPointerType::New();
// 		double gridOrigin[3];
// 		gridOrigin[0] = sampleBoundingBox[0];
// 		gridOrigin[1] = sampleBoundingBox[2];
// 		gridOrigin[2] = HVCBounds[4];
// 		spheresAverageGrid->SetOrigin(gridOrigin);
// 		spheresAverageGrid->SetSpacing(sphereGridSpacing, sphereGridSpacing, sphereGridSpacing);
// 		spheresAverageGrid->SetDimensions(sphereGrid);
// 		for(int index = 0; index < individualSphereInnervations.size(); ++index)
// 		{
// 			int x, y, z;
// 			x = int(index/(sphereGrid[2]*sphereGrid[1]));
// 			y = int(index/sphereGrid[2]) - sphereGrid[1]*x;
// 			z = index - (sphereGrid[2]*(y + sphereGrid[1]*x));
// 			double avgAngle = -1;
// 			if(individualSphereInnervations[index].size())
// 			{
// 				avgAngle = 0;
// 				for(int n = 0; n < individualSphereInnervations[index].size(); ++n)
// 				{
// 					double * tmpVec = individualSphereInnervations[index][n];
// 					double tmpAngle = atan2(tmpVec[1],tmpVec[0]);
// 					avgAngle += tmpAngle;
// 				}
// 				avgAngle /= individualSphereInnervations[index].size();
// 			}
// 			double * angleVal = static_cast< double * >(spheresAverageGrid->GetScalarPointer(x, y, z));
// 			*angleVal = avgAngle;
// 			totalAvgAngle += avgAngle;
// 			++totalNrSpheres;
// 			spheresAngleVec.push_back(avgAngle);
// 		}
// 		std::string spheresAvgOutputname(outputFilename);
// 		spheresAvgOutputname += "_spheres_avg";
// 		Reader * spheresAvgMeshWriter = new Reader(spheresAvgOutputname.c_str(), spheresAvgOutputname.c_str());
// 		spheresAvgMeshWriter->writeScalarField(spheresAverageGrid);
// 		delete spheresAvgMeshWriter;
// 		// output: average across spheres
// 		totalAvgAngle /= totalNrSpheres;
// 		double totalStdAngle = 0;
// 		for(int i = 0; i < spheresAngleVec.size(); ++i)
// 		{
// 			totalStdAngle += (spheresAngleVec[i] - totalAvgAngle)*(spheresAngleVec[i] - totalAvgAngle);
// 		}
// 		totalStdAngle = sqrt(totalStdAngle/totalNrSpheres);
// 		std::cout << "Avg angle = " << totalAvgAngle*180/PI << " +- " << totalStdAngle*180/PI << std::endl;
// 		
// 		std::string spheresAllDirectionsOutputname(outputFilename);
// 		spheresAllDirectionsOutputname += "_spheres_all_directions.csv";
// 		std::ofstream allDirectionsFile(spheresAllDirectionsOutputname.c_str());
// 		allDirectionsFile << "All angles across all spheres" << std::endl;
// 		for(int i = 0; i < individualSphereInnervations.size(); ++i)
// 			for(int j = 0; j < individualSphereInnervations[i].size(); ++j)
// 			{
// 				double * tmpVec = individualSphereInnervations[i][j];
// 				double tmpAngle = atan2(tmpVec[1],tmpVec[0]);
// 				allDirectionsFile << tmpAngle << std::endl;
// 			}
// 		allDirectionsFile.close();
// 		
// 		std::string spheresNonInnervationAllDirectionsOutputname(outputFilename);
// 		spheresNonInnervationAllDirectionsOutputname += "_spheres_all_directions_noInnervation.csv";
// 		std::ofstream noInnervationAllDirectionsFile(spheresNonInnervationAllDirectionsOutputname.c_str());
// 		noInnervationAllDirectionsFile << "All angles across all spheres" << std::endl;
// 		for(int i = 0; i < individualSphereNonInnervations.size(); ++i)
// 			for(int j = 0; j < individualSphereNonInnervations[i].size(); ++j)
// 			{
// 				double * tmpVec = individualSphereNonInnervations[i][j];
// 				double tmpAngle = atan2(tmpVec[1],tmpVec[0]);
// 				noInnervationAllDirectionsFile << tmpAngle << std::endl;
// 			}
// 		noInnervationAllDirectionsFile.close();
		
		double distanceBinSize = 20;
		Profile * innervatedProfile = new Profile(distanceBinSize);
		for(int i = 0; i < individualSphereInnervationDistances.size(); ++i)
			for(int j = 0; j < individualSphereInnervationDistances[i].size(); ++j)
			{
				double tmpDist = individualSphereInnervationDistances[i][j];
				int tmpBin = int(tmpDist/distanceBinSize);
				innervatedProfile->incrementBin(tmpBin);
			}
		std::string spheresInnervationDistancesOutputname(outputFilename);
		spheresInnervationDistancesOutputname += "_spheres_all_distances_innervation.csv";
		std::ofstream spheresInnervationDistancesFile(spheresInnervationDistancesOutputname.c_str());
		spheresInnervationDistancesFile << "Bin begin\tbin end\tcount\n";
		for(int i = 0; i < innervatedProfile->getProfile()->size(); ++i)
		{
			spheresInnervationDistancesFile << i*distanceBinSize << "\t" << (i+1)*distanceBinSize << "\t" << innervatedProfile->getProfile()->at(i) << std::endl;
		}
		spheresInnervationDistancesFile.close();
		
		Profile * nonInnervatedProfile = new Profile(distanceBinSize);
		for(int i = 0; i < individualSphereNonInnervationDistances.size(); ++i)
			for(int j = 0; j < individualSphereNonInnervationDistances[i].size(); ++j)
			{
				double tmpDist = individualSphereNonInnervationDistances[i][j];
				int tmpBin = int(tmpDist/distanceBinSize);
				nonInnervatedProfile->incrementBin(tmpBin);
			}
		std::string spheresNonInnervationDistancesOutputname(outputFilename);
		spheresNonInnervationDistancesOutputname += "_spheres_all_distances_noInnervation.csv";
		std::ofstream spheresNonInnervationDistancesFile(spheresNonInnervationDistancesOutputname.c_str());
		spheresNonInnervationDistancesFile << "Bin begin\tbin end\tcount\n";
		for(int i = 0; i < nonInnervatedProfile->getProfile()->size(); ++i)
		{
			spheresNonInnervationDistancesFile << i*distanceBinSize << "\t" << (i+1)*distanceBinSize << "\t" << nonInnervatedProfile->getProfile()->at(i) << std::endl;
		}
		spheresNonInnervationDistancesFile.close();
	}
	
	if(argc == 5)
	{
		const char * sgListFilename = argv[1];
		const char * avgHVCFilename = argv[2];
		const char * outputFilename = argv[3];
		int mode = atoi(argv[4]);
// defines for different modes
#define DATA 0
#define CLUSTERS 1 // synthetic
#define UNIFORM 2  // synthetic
#define CLOSEST 3  // synthetic
		
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
					Reader * sgReader = new Reader(tmpChar, outputFilename);
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
					}
				}
#ifdef DEBUG
		std::cout << "axonInsideCount = " << axonInsideCount << std::endl;
#endif
		std::cout << "Total number of spheres inside avg. HVC contour = " << totalSpheresInsideHVC << std::endl;
		const gsl_rng_type * T;
		gsl_rng * r;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		int usedInnervatedSpheres = 0;
		std::vector< double > avgPairwiseDistances;
		std::vector< double > innervationRatio;
		std::vector< int > clusterSizes;
		int bootStrapSamples = 10;
		std::cout << "Using mode " << mode << " (0 - DATA; 1 - CLUSTERS; 2 - UNIFORM; 3 - CLOSEST)" << std::endl;
		std::cout << "bootStrapSamples = " << bootStrapSamples << std::endl;
		for(int i = 0; i < individualSphereInnervations.size(); ++i)
		{
			int clusterSize = individualSphereInnervations[i].size();
			double tmpRatio = double(clusterSize)/double(individualSpatialGraphs.size());
			innervationRatio.push_back(tmpRatio);
			if(clusterSize < 2)
			{
				continue;
			}
			++usedInnervatedSpheres;
			if(mode == DATA)
			{
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
				avgPairwiseDistances.push_back(avgDist);
				clusterSizes.push_back(clusterSize);
			}
			else if(mode == CLUSTERS)
			{
				for(int n = 0; n < bootStrapSamples; ++n)
				{
					int neuronID1 = gsl_rng_uniform_int(r, clusterSize);
					double * somaLoc1 = somaLocations[neuronID1];
					std::map< double, int > distanceNeuronIDs;
					for(int j = 0; j < somaLocations.size(); ++j)
					{
						if(j == neuronID1)
						{
							continue;
						}
						double * somaLoc2 = somaLocations[j];
						double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(somaLoc1, somaLoc2));
						distanceNeuronIDs[tmpDist] = j;
					}
					std::vector< int > closestNeurons;
					closestNeurons.push_back(neuronID1);
					for(std::map< double, int >::const_iterator distanceNeuronIDsIt = distanceNeuronIDs.begin(); 
						distanceNeuronIDsIt != distanceNeuronIDs.end() && closestNeurons.size() < clusterSize;
						++distanceNeuronIDsIt)
					{
						closestNeurons.push_back(distanceNeuronIDsIt->second);
					}
					double totalDist = 0;
					for(int j = 0; j < clusterSize; ++j)
					{
						neuronID1 = closestNeurons[j];
						somaLoc1 = somaLocations[neuronID1];
						for(int k = j + 1; k < clusterSize; ++k)
						{
							int neuronID2 = closestNeurons[k];
							double * somaLoc2 = somaLocations[neuronID2];
							double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(somaLoc1, somaLoc2));
							totalDist += tmpDist;
						}
					}
					int pairwiseComparisons = clusterSize*(clusterSize - 1)/2;
					double avgDist = totalDist/pairwiseComparisons;
					avgPairwiseDistances.push_back(avgDist);
					clusterSizes.push_back(clusterSize);
				}
			}
			else if(mode == UNIFORM)
			{
				for(int n = 0; n < bootStrapSamples; ++n)
				{
					int * tmpNeuronIDs = new int[somaLocations.size()];
					for(int j = 0; j < somaLocations.size(); ++j)
					{
						tmpNeuronIDs[j] = j;
					}
					int * randomNeuronIDs = new int[clusterSize];
					gsl_ran_choose(r, randomNeuronIDs, clusterSize, tmpNeuronIDs, somaLocations.size(), sizeof(int));
					
					double totalDist = 0;
					for(int j = 0; j < clusterSize; ++j)
					{
						int neuronID1 = randomNeuronIDs[j];
						double * somaLoc1 = somaLocations[neuronID1];
						for(int k = j + 1; k < clusterSize; ++k)
						{
							int neuronID2 = randomNeuronIDs[k];
							double * somaLoc2 = somaLocations[neuronID2];
							double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(somaLoc1, somaLoc2));
							totalDist += tmpDist;
						}
					}
					int pairwiseComparisons = clusterSize*(clusterSize - 1)/2;
					double avgDist = totalDist/pairwiseComparisons;
					avgPairwiseDistances.push_back(avgDist);
					clusterSizes.push_back(clusterSize);
				}
			}
			else if(mode == CLOSEST)
			{
				
			}
		}
		
		std::cout << "Total number of spheres used (i.e. >= 2 innervating neurons) = " << usedInnervatedSpheres << std::endl;
		std::string outName(outputFilename);
		std::string outName2(outputFilename);
		std::string outName3(outputFilename);
		if(mode == DATA)
		{
			outName += "_data_pairwise_distances_completeHVC.csv";
			outName2 += "_data_innervationRatio_completeHVC.csv";
			outName3 += "_data_clusterSizes_completeHVC.csv";
			std::cout << "outName = " << outName.c_str() << std::endl;
			std::cout << "outName2 = " << outName2.c_str() << std::endl;
			std::cout << "outName3 = " << outName3.c_str() << std::endl;
		}
		else if(mode == CLUSTERS)
		{
			char suffix[64];
			sprintf(suffix, "_clusters_pairwise_distances_bootstrap_%d_completeHVC.csv", bootStrapSamples);
			outName += suffix;
// 			outName += "_clusters_pairwise_distances_bootstrap_";
// 			outName += std::to_string(bootStrapSamples);
// 			outName += ".csv";
			std::cout << "outName = " << outName.c_str() << std::endl;
			outName3 += "_clusters_clusterSizes_completeHVC.csv";
		}
		else if(mode == UNIFORM)
		{
			char suffix[64];
			sprintf(suffix, "_uniform_pairwise_distances_bootstrap_%d_completeHVC.csv", bootStrapSamples);
			outName += suffix;
// 			outName += "_uniform_pairwise_distances_bootstrap_";
// 			outName += std::to_string(bootStrapSamples);
// 			outName += ".csv";
			std::cout << "outName = " << outName.c_str() << std::endl;
			outName3 += "_uniform_clusterSizes_completeHVC.csv";
		}
		else if(mode == CLOSEST)
		{
			char suffix[64];
			sprintf(suffix, "_closest_pairwise_distances_bootstrap_%d_completeHVC.csv", bootStrapSamples);
			outName += suffix;
// 			outName += "_closest_pairwise_distances_bootstrap_";
// 			outName += std::to_string(bootStrapSamples);
// 			outName += ".csv";
			std::cout << "outName = " << outName.c_str() << std::endl;
		}
		std::ofstream outFile(outName.c_str());
		for(int i = 0; i < avgPairwiseDistances.size(); ++i)
		{
			outFile << avgPairwiseDistances[i] << std::endl;
		}
		outFile.close();
			
		std::ofstream outFile3(outName3.c_str());
		for(int i = 0; i < clusterSizes.size(); ++i)
		{
			outFile3 << clusterSizes[i] << std::endl;
		}
		outFile3.close();
		
		if(mode == DATA)
		{
			std::ofstream outFile2(outName2.c_str());
			for(int i = 0; i < innervationRatio.size(); ++i)
			{
				outFile2 << innervationRatio[i] << std::endl;
			}
			outFile2.close();
		}
	}
	
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


