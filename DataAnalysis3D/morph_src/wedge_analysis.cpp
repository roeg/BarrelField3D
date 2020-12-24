/****************************************************************************/
/*                                                                          */
/* Program:   MorphAnalyzer                                                 */
/*                                                                          */
/* File:      main.cpp                                                      */
/*                                                                          */
/* Purpose:   Program for analysis of registered neuron morphologies with   */
/*            respect to columns, septa and layers in standardized barrel   */
/*            cortex. E.g., computes distribution of axon in different      */
/*            columns/layers, and computes z-profiles of axons taking local */
/*            orientation into account                                      */
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
// #include "morph_analyzer.h"
#include "../../common/inputcheckpoint.h"
#include "../../common/amiraReader.h"
#include "../../common/basics.h"
#include "../../common/profile.h"

#define HORIZONTAL 0
#define CORONAL 1
#define SAGITTAL 2

// #define DEBUG

void computeAxonWedgeLengths(AmiraSpatialGraph * neuronMorphology, std::vector< ClosedSurface * > wedgeSurfaces, int orientation, const char * outputFilename);
void computeHVCWedgeVolumes(AmiraSpatialGraph * neuronMorphology, PolyDataPointerType HVCSurface, std::vector< ClosedSurface * > wedgeSurfaces, double voxelSize, const char * outputFilename);
void computeLowpassFilteredAngularDensity(AmiraSpatialGraph* neuronMorphology, PolyDataPointerType HVCSurface, std::vector< ClosedSurface* > wedgeSurfaces, int orientation, double voxelSize, const char* outputFilename);
void computeLowpassFilteredAngularDensity2(AmiraSpatialGraph* neuronMorphology, PolyDataPointerType HVCSurface, std::vector< ClosedSurface* > wedgeSurfaces, int orientation, double voxelSize, const char* outputFilename);
void computeAxonInsideOutsideSphere(AmiraSpatialGraph * neuronMorphology, double proximalBorder, const char * outputFilename);

int main( int argc , char * argv[])
{
	if(argc == 3)
	{
		const char * inputFilename = argv[1];
		double proximalBorder = atof(argv[2]);
		
		std::string ifName(inputFilename);
		size_t suffix;
		Reader * hocReader = new Reader(inputFilename, inputFilename);
		if(ifName.find(".am") != std::string::npos)
		{
			hocReader->readSpatialGraphFile(0);
			suffix = ifName.find(".am");
		}
		else
		{
			std::cout << "Error! Can only analyze .hoc or .am files!" << std::endl;
			delete hocReader;
			return 0;
		}
		
#ifdef DEBUG
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		if(checkPoint->getParameters().axonFlag)
		{
			computeAxonInsideOutsideSphere(hocReader->getSpatialGraph(), proximalBorder, inputFilename);
		}
#endif
#ifndef DEBUG
		computeAxonInsideOutsideSphere(hocReader->getSpatialGraph(), proximalBorder, inputFilename);
#endif
	}
	
	else if(argc == 5)
	{
		const char * inputFilename = argv[1];
		const char * wedgeSurfacesFilename = argv[2];
		int orientation = atoi(argv[3]);
		const char * outputFilename = argv[4];
		double binSize = 10;
		
		std::string ifName(inputFilename);
		size_t suffix;
		
		Reader * hocReader = new Reader(inputFilename, inputFilename);
		if(ifName.find(".hoc") != std::string::npos)
		{
			hocReader->readHocFile();
			suffix = ifName.find(".hoc");
		}
		else if(ifName.find(".am") != std::string::npos)
		{
			hocReader->readSpatialGraphFile(0);
			suffix = ifName.find(".am");
		}
		else
		{
			std::cout << "Error! Can only analyze .hoc or .am files!" << std::endl;
			return 0;
		}
		
		std::vector< ClosedSurface * > individualSurfaces;
		std::ifstream inputStream(wedgeSurfacesFilename);
		if(!inputStream.fail())
		{
			std::cout << "Loading surface files from " << wedgeSurfacesFilename << std::endl;
			std::string currentLine;
			while(!std::getline(inputStream, currentLine).eof())
				if(currentLine.size())
				{
					char * tmpChar = new char[128];
					sscanf(currentLine.c_str(), " %s ", tmpChar);
					
					std::string surfaceStr(tmpChar);
					Reader * wedgeSurfaceReader = new Reader(tmpChar, outputFilename);
					if(surfaceStr.find(".surf") != std::string::npos)
					{
						std::cout << "Loading surface file " << tmpChar << std::endl;
						PolyDataPointerType wedgeSurface_ = wedgeSurfaceReader->readAmiraSurfaceFile();
						ClosedSurface * wedgeSurface = new ClosedSurface(wedgeSurface_);
						individualSurfaces.push_back(wedgeSurface);
					}
					else
					{
						std::cout << "Error! Surface file" << tmpChar << " has to be Amira '.surf' ascii file!" << std::endl;
						return 0;
					}
					
					delete wedgeSurfaceReader;
					delete [] tmpChar;
				}
		}
		
		std::string ofName(ifName, 0, suffix);
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		
		if(checkPoint->getParameters().axonFlag)
		{
			computeAxonWedgeLengths(hocReader->getSpatialGraph(), individualSurfaces, orientation, outputFilename);
			
// 			std::vector< double > zProfile = morphAnalyzer->compute1DProfileLocally("Axon", 50);
// 			ofName += "_axon_1D_profile.csv";
// 			if(zProfile.size())
// 			{
// 				std::ofstream ProfileWriter;
// 				ProfileWriter.open(ofName.c_str());
// 				ProfileWriter << "# 1D Profile of axon" << std::endl;
// 				ProfileWriter << "Depth [um]\tlength per bin [um]" << std::endl;
// 				for(int ii = 0; ii < zProfile.size(); ++ii)
// 					ProfileWriter << ii*50+25 << "\t" << zProfile[ii] << std::endl;
// 				ProfileWriter.close();
// 			}
// 			else
// 				std::cout << "Error! Z Profile is empty!" << std::endl;
		}
		
		delete checkPoint;
		delete hocReader;
	}
	
	else if(argc == 6)
	{
		const char * inputFilename = argv[1];
		const char * HVCSurfaceFilename = argv[2];
		const char * wedgeSurfacesFilename = argv[3];
		double voxelSize = atof(argv[4]);
		const char * outputFilename = argv[5];
		
		std::string ifName(inputFilename);
		size_t suffix;
		
		Reader * hocReader = new Reader(inputFilename, inputFilename);
		if(ifName.find(".hoc") != std::string::npos)
		{
			hocReader->readHocFile();
			suffix = ifName.find(".hoc");
		}
		else if(ifName.find(".am") != std::string::npos)
		{
			hocReader->readSpatialGraphFile(1);
			suffix = ifName.find(".am");
		}
		else
		{
			std::cout << "Error! Can only analyze .hoc or .am files!" << std::endl;
			return 0;
		}
		
		std::string surfaceStr(HVCSurfaceFilename);
		Reader * hvcSurfaceReader = new Reader(HVCSurfaceFilename, outputFilename);
		PolyDataPointerType HVCSurface = PolyDataPointerType::New();
		if(surfaceStr.find(".surf") != std::string::npos)
		{
			std::cout << "Loading surface file " << HVCSurfaceFilename << std::endl;
			HVCSurface = hvcSurfaceReader->readAmiraSurfaceFile();
		}
		else
		{
			std::cout << "Error! Surface file" << HVCSurfaceFilename << " has to be Amira '.surf' ascii file!" << std::endl;
			return 0;
		}
		
		std::vector< ClosedSurface * > individualSurfaces;
		std::ifstream inputStream(wedgeSurfacesFilename);
		if(!inputStream.fail())
		{
			std::cout << "Loading surface files from " << wedgeSurfacesFilename << std::endl;
			std::string currentLine;
			while(!std::getline(inputStream, currentLine).eof())
				if(currentLine.size())
				{
					char * tmpChar = new char[128];
					sscanf(currentLine.c_str(), " %s ", tmpChar);
					
					std::string surfaceStr(tmpChar);
					Reader * wedgeSurfaceReader = new Reader(tmpChar, outputFilename);
					if(surfaceStr.find(".surf") != std::string::npos)
					{
						std::cout << "Loading surface file " << tmpChar << std::endl;
						PolyDataPointerType wedgeSurface_ = wedgeSurfaceReader->readAmiraSurfaceFile();
						ClosedSurface * wedgeSurface = new ClosedSurface(wedgeSurface_);
						individualSurfaces.push_back(wedgeSurface);
					}
					else
					{
						std::cout << "Error! Surface file" << tmpChar << " has to be Amira '.surf' ascii file!" << std::endl;
						return 0;
					}
					
					delete wedgeSurfaceReader;
					delete [] tmpChar;
				}
		}
		
		std::string ofName(ifName, 0, suffix);
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		
		if(checkPoint->getParameters().somaFlag)
		{
			computeHVCWedgeVolumes(hocReader->getSpatialGraph(), HVCSurface, individualSurfaces, voxelSize, outputFilename);
		}
		else
		{
			std::cout << "Error: Neuron morphology does not contain a soma!" << std::endl;
		}
		
		delete checkPoint;
		delete hocReader, delete hvcSurfaceReader;
	}
	
	else if(argc == 7)
	{
		const char * inputFilename = argv[1];
		const char * HVCSurfaceFilename = argv[2];
		const char * wedgeSurfacesFilename = argv[3];
		int orientation = atoi(argv[4]);
		double voxelSize = atof(argv[5]);
		const char * outputFilename = argv[6];
		
		std::string ifName(inputFilename);
		size_t suffix;
		
		Reader * hocReader = new Reader(inputFilename, inputFilename);
		if(ifName.find(".hoc") != std::string::npos)
		{
			hocReader->readHocFile();
			suffix = ifName.find(".hoc");
		}
		else if(ifName.find(".am") != std::string::npos)
		{
			hocReader->readSpatialGraphFile(1);
			suffix = ifName.find(".am");
		}
		else
		{
			std::cout << "Error! Can only analyze .hoc or .am files!" << std::endl;
			return 0;
		}
		
		std::string surfaceStr(HVCSurfaceFilename);
		Reader * hvcSurfaceReader = new Reader(HVCSurfaceFilename, outputFilename);
		PolyDataPointerType HVCSurface = PolyDataPointerType::New();
		if(surfaceStr.find(".surf") != std::string::npos)
		{
			std::cout << "Loading surface file " << HVCSurfaceFilename << std::endl;
			HVCSurface = hvcSurfaceReader->readAmiraSurfaceFile();
		}
		else
		{
			std::cout << "Error! Surface file" << HVCSurfaceFilename << " has to be Amira '.surf' ascii file!" << std::endl;
			return 0;
		}
		
		std::vector< ClosedSurface * > individualSurfaces;
		std::ifstream inputStream(wedgeSurfacesFilename);
		if(!inputStream.fail())
		{
			std::cout << "Loading surface files from " << wedgeSurfacesFilename << std::endl;
			std::string currentLine;
			while(!std::getline(inputStream, currentLine).eof())
				if(currentLine.size())
				{
					char * tmpChar = new char[128];
					sscanf(currentLine.c_str(), " %s ", tmpChar);
					
					std::string surfaceStr(tmpChar);
					Reader * wedgeSurfaceReader = new Reader(tmpChar, outputFilename);
					if(surfaceStr.find(".surf") != std::string::npos)
					{
						std::cout << "Loading surface file " << tmpChar << std::endl;
						PolyDataPointerType wedgeSurface_ = wedgeSurfaceReader->readAmiraSurfaceFile();
						ClosedSurface * wedgeSurface = new ClosedSurface(wedgeSurface_);
						individualSurfaces.push_back(wedgeSurface);
					}
					else
					{
						std::cout << "Error! Surface file" << tmpChar << " has to be Amira '.surf' ascii file!" << std::endl;
						return 0;
					}
					
					delete wedgeSurfaceReader;
					delete [] tmpChar;
					// only load first wedge surface file
					break;
				}
		}
		
		std::string ofName(ifName, 0, suffix);
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		
		if(checkPoint->getParameters().somaFlag)
		{
			computeLowpassFilteredAngularDensity(hocReader->getSpatialGraph(), HVCSurface, individualSurfaces, orientation, voxelSize, outputFilename);
// 			computeLowpassFilteredAngularDensity2(hocReader->getSpatialGraph(), HVCSurface, individualSurfaces, orientation, voxelSize, outputFilename);
		}
		else
		{
			std::cout << "Error: Neuron morphology does not contain a soma!" << std::endl;
		}
		
		delete checkPoint;
		delete hocReader, delete hvcSurfaceReader;
	}
	
	else
	{
		std::cout << "Usage: MorphAnalysis [Input filename] [Output filename] [Label] [Binsize in micron] [Mode: Global - 1 ; Local - 2]" << std::endl;
	}
	return 0;
}

void computeAxonWedgeLengths(AmiraSpatialGraph* neuronMorphology, std::vector< ClosedSurface * > wedgeSurfaces, int orientation, const char* outputFilename)
{
	double proximalBorder = 200; // microns
	std::vector< double > wedgeLengths;
	std::vector< double > wedgeLengthsProximal;
	std::vector< Profile * > wedgeDistanceLengths;
	for(int i = 0; i < wedgeSurfaces.size(); ++i)
	{
		wedgeLengths.push_back(0);
		wedgeLengthsProximal.push_back(0);
		Profile * tmpProfile = new Profile(proximalBorder);
		wedgeDistanceLengths.push_back(tmpProfile);
// 		if(i == 0)
// 		{
// 			wedgeSurfaces[i]->ptr()->Print(std::cout);
// 		}
	}
	
	double somaLoc[3];
	PolyDataPointerType thisSoma = PolyDataPointerType::New();
	if(neuronMorphology->extractLandmark(Soma, thisSoma))
	{
		int subId;
		double pCenter[3], * weights = new double[thisSoma->GetNumberOfPoints()];
		thisSoma->GetCell(0)->GetParametricCenter(pCenter);
		thisSoma->GetCell(0)->EvaluateLocation(subId, pCenter, somaLoc, weights);
		delete [] weights;
	}
	else
	{
		std::cout << "Error! SpatialGraph file does not contain a soma!" << std::endl;
		return;
	}
	
	double somaShift[3];
	somaShift[0] = -somaLoc[0];
	somaShift[1] = -somaLoc[1];
	somaShift[2] = -somaLoc[2];
	
	TransformPointerType shiftTransform = TransformPointerType::New();
	shiftTransform->Translate(somaShift);
	neuronMorphology->setTransformation(shiftTransform);
	neuronMorphology->applyTransformation();
	
	// main loop
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
	{
// 		if((*edgeIt)->label != Axon)
		if((*edgeIt)->label != Dendrite && (*edgeIt)->label != BasalDendrite && (*edgeIt)->label != ApicalDendrite)
			continue;
		
		std::list< double * >::const_iterator ptIt1, ptIt2;
		std::list< double * > tmpEdgePointList((*edgeIt)->edgePointCoordinates);
		ptIt1 = tmpEdgePointList.begin();
		ptIt2 = tmpEdgePointList.begin();
		++ptIt1;
		while(ptIt1 != tmpEdgePointList.end())
		{
#ifdef DEBUG
			std::flush(std::cout << "Processing point on edge " << *edgeIt << " @ [" << (*ptIt1)[0] << "," << (*ptIt1)[1] << "," << (*ptIt1)[2] << "]\r");
#endif
			double * pt1, * pt2;
			int pt1WedgeIndex = -1;
			int pt2WedgeIndex = -1;
			
			pt1 = *ptIt1;
			pt2 = *ptIt2;
			
			double dist1, dist2;
			dist1 = sqrt(vtkMath::Distance2BetweenPoints(pt1, somaLoc));
			dist2 = sqrt(vtkMath::Distance2BetweenPoints(pt2, somaLoc));
// 			if(orientation == HORIZONTAL)
// 			{
// 				dist1 = sqrt((pt1[0]-somaLoc[0])*(pt1[0]-somaLoc[0]) + (pt1[1]-somaLoc[1])*(pt1[1]-somaLoc[1]));
// 				dist2 = sqrt((pt2[0]-somaLoc[0])*(pt2[0]-somaLoc[0]) + (pt2[1]-somaLoc[1])*(pt2[1]-somaLoc[1]));
// 			}
// 			if(orientation == CORONAL)
// 			{
// 				dist1 = sqrt((pt1[0]-somaLoc[0])*(pt1[0]-somaLoc[0]) + (pt1[2]-somaLoc[2])*(pt1[2]-somaLoc[2]));
// 				dist2 = sqrt((pt2[0]-somaLoc[0])*(pt2[0]-somaLoc[0]) + (pt2[2]-somaLoc[2])*(pt2[2]-somaLoc[2]));
// 			}
// 			if(orientation == SAGITTAL)
// 			{
// 				dist1 = sqrt((pt1[1]-somaLoc[1])*(pt1[1]-somaLoc[1]) + (pt1[2]-somaLoc[2])*(pt1[2]-somaLoc[2]));
// 				dist2 = sqrt((pt2[1]-somaLoc[1])*(pt2[1]-somaLoc[1]) + (pt2[2]-somaLoc[2])*(pt2[2]-somaLoc[2]));
// 			}
			
			// if points on different sides of proximal border:
			// virtually insert point 
			if((dist1 <= proximalBorder && dist2 > proximalBorder) || (dist1 > proximalBorder && dist2 <= proximalBorder))
			{
				// TBD - for now just lump to proximal point
			}
			
			for(int i = 0; i < wedgeSurfaces.size(); ++i)
			{
				ClosedSurface *  tmpSurface = wedgeSurfaces[i];
				if(tmpSurface->isPointInsideSurface(pt1))
				{
					pt1WedgeIndex = i;
// 					std::cout << "pt1WedgeIndex = " << pt1WedgeIndex << std::endl;
				}
				if(tmpSurface->isPointInsideSurface(pt2))
				{
					pt2WedgeIndex = i;
// 					std::cout << "pt2WedgeIndex = " << pt2WedgeIndex << std::endl;
				}
			}
			
// 			#ifdef DEBUG
// 			assignDistanceToLayers(pt1, pt2, pt1Layer, pt2Layer, onlyLayerProfiles, onlyLayerProfiles);
// 			#endif
			
			// sanity check
			if(pt1WedgeIndex == -1 || pt2WedgeIndex == -1)
			{
				std::cout << "Error: Point outside of wedges: " << std::endl;
				std::cout << "pt1WedgeIndex = " << pt1WedgeIndex << std::endl;
				std::cout << "pt2WedgeIndex = " << pt2WedgeIndex << std::endl;
				std::cout << "pt1  @ [" << pt1[0] << "," << pt1[1] << "," << pt1[2] << "]" << std::endl;
				std::cout << "pt2  @ [" << pt2[0] << "," << pt2[1] << "," << pt2[2] << "]" << std::endl;
				
				++ptIt1, ++ptIt2;
				continue;
			}
			
			// case 1 : both in same wedge
			if(pt1WedgeIndex == pt2WedgeIndex)
			{
				double length = sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
				if(dist1 <= proximalBorder)
				{
					wedgeLengthsProximal[pt1WedgeIndex] += length;
				}
				wedgeLengths[pt1WedgeIndex] += length;
				int bin = int(dist1/proximalBorder);
				wedgeDistanceLengths[pt1WedgeIndex]->addSegment(length, bin);
			}
			// case 2: points in neighboring wedges (assume that axon sampling distance is sufficiently small)
			else
			{
				double line12[3];
				ClosedSurface * pt1WedgeSurface = wedgeSurfaces[pt1WedgeIndex];
				vtkMath::Subtract(pt2, pt1, line12);
// 				pt1WedgeSurface->intersectLineInDirection(line12, pt1);
				pt1WedgeSurface->intersectLineSegment(pt1, pt2);
				double intersectPt[3];
				if(pt1WedgeSurface->isIntersectionFound())
				{
					pt1WedgeSurface->getLastIntersectPoint(intersectPt);
				}
				else
				{
					std::cout << "Error: Could not find intersection point with wedge " << pt1WedgeIndex << std::endl;
					++ptIt1, ++ptIt2;
					continue;
				}
				
				double distance1 = sqrt(vtkMath::Distance2BetweenPoints(pt1, intersectPt));
				double distance2 = sqrt(vtkMath::Distance2BetweenPoints(pt2, intersectPt));
				if(dist1 <= proximalBorder)
				{
					wedgeLengthsProximal[pt1WedgeIndex] += distance1;
				}
				wedgeLengths[pt1WedgeIndex] += distance1;
				int bin1 = int(dist1/proximalBorder);
				wedgeDistanceLengths[pt1WedgeIndex]->addSegment(distance1, bin1);
				if(dist2 <= proximalBorder)
				{
					wedgeLengthsProximal[pt2WedgeIndex] += distance2;
				}
				wedgeLengths[pt2WedgeIndex] += distance2;
				int bin2 = int(dist2/proximalBorder);
				wedgeDistanceLengths[pt2WedgeIndex]->addSegment(distance2, bin2);
			}
			
			++ptIt1, ++ptIt2;
		}
	}
#ifdef DEBUG
		std::cout << "\nDone processing axon length!" << std::endl;
#endif
	
	#ifndef DEBUG
	// write output files
	
	std::string axonOutName(outputFilename);
// 	axonOutName += "_axon_wedge_length.csv";
	axonOutName += "_dendrite_wedge_length.csv";
	std::ofstream AxonFile;
	AxonFile.open(axonOutName.c_str());
	AxonFile << "Wedge\tLength (microns)" << std::endl;
	for(int i = 0; i < wedgeLengths.size(); ++i)
	{
		AxonFile << i << "\t" << wedgeLengths[i] << std::endl;
	}
	AxonFile.close();
	
// 	std::string axonOutName2(outputFilename);
// 	axonOutName2 += "_axon_wedge_length_proximal.csv";
// 	std::ofstream AxonFile2;
// 	AxonFile2.open(axonOutName2.c_str());
// 	AxonFile2 << "Wedge\tLength (microns)" << std::endl;
// 	for(int i = 0; i < wedgeLengthsProximal.size(); ++i)
// 	{
// 		AxonFile2 << i << "\t" << wedgeLengthsProximal[i] << std::endl;
// 	}
// 	AxonFile2.close();
// 	
// 	std::string axonOutName3(outputFilename);
// 	axonOutName3 += "_axon_wedge_length_distal.csv";
// 	std::ofstream AxonFile3;
// 	AxonFile3.open(axonOutName3.c_str());
// 	AxonFile3 << "Wedge\tLength (microns)" << std::endl;
// 	for(int i = 0; i < wedgeLengths.size(); ++i)
// 	{
// 		AxonFile3 << i << "\t" << wedgeLengths[i] - wedgeLengthsProximal[i] << std::endl;
// 	}
// 	AxonFile3.close();
// 	
// 	int maxBins = 0;
// 	for(int i = 0; i < wedgeSurfaces.size(); ++i)
// 	{
// 		std::vector< double > * profilePtr = wedgeDistanceLengths[i]->getProfile();
// 		if(profilePtr->size() > maxBins)
// 		{
// 			maxBins = profilePtr->size();
// 		}
// 	}
// 	std::string axonOutName4(outputFilename);
// 	axonOutName4 += "_axon_wedge_length_distance_bins.csv";
// 	std::ofstream AxonFile4;
// 	AxonFile4.open(axonOutName4.c_str());
// 	AxonFile4 << "Wedge\t";
// 	for(int i = 0; i < maxBins; ++i)
// 	{
// 		AxonFile4 << i*proximalBorder << "\t";
// 	}
// 	AxonFile4 << std::endl;
// 	for(int i = 0; i < wedgeDistanceLengths.size(); ++i)
// 	{
// 		AxonFile4 << i << "\t";
// 		std::vector< double > * profilePtr = wedgeDistanceLengths[i]->getProfile();
// 		int binCount = 0;
// 		for(int j = 0; j < profilePtr->size(); ++j)
// 		{
// 			AxonFile4 << profilePtr->at(j) << "\t";
// 			++binCount;
// 		}
// 		for(int j = 0; j < maxBins - binCount; ++j)
// 		{
// 			AxonFile4 << "0\t";
// 		}
// 		AxonFile4 << std::endl;
// 	}
// 	AxonFile4.close();
	#endif
}

void computeHVCWedgeVolumes(AmiraSpatialGraph* neuronMorphology, PolyDataPointerType HVCSurface, std::vector< ClosedSurface* > wedgeSurfaces, double voxelSize, const char* outputFilename)
{
	double proximalBorder = 200; // microns
	std::vector< double > wedgeVolumes;
	std::vector< double > wedgeVolumesProximal;
	std::vector< Profile * > wedgeDistanceVolumes;
	for(int i = 0; i < wedgeSurfaces.size(); ++i)
	{
		wedgeVolumes.push_back(0);
		wedgeVolumesProximal.push_back(0);
		Profile * tmpProfile = new Profile(proximalBorder);
		wedgeDistanceVolumes.push_back(tmpProfile);
// 		if(i == 0)
// 		{
// 			wedgeSurfaces[i]->ptr()->Print(std::cout);
// 		}
	}
	
	double somaLoc[3];
	PolyDataPointerType thisSoma = PolyDataPointerType::New();
	if(neuronMorphology->extractLandmark(Soma, thisSoma))
	{
		int subId;
		double pCenter[3], * weights = new double[thisSoma->GetNumberOfPoints()];
		thisSoma->GetCell(0)->GetParametricCenter(pCenter);
		thisSoma->GetCell(0)->EvaluateLocation(subId, pCenter, somaLoc, weights);
		delete [] weights;
	}
	else
	{
		std::cout << "Error! SpatialGraph file does not contain a soma!" << std::endl;
		return;
	}
	
	double surfaceShift[3];
	surfaceShift[0] = -somaLoc[0];
	surfaceShift[1] = -somaLoc[1];
	surfaceShift[2] = -somaLoc[2];
	
	TransformPointerType shiftTransform = TransformPointerType::New();
	shiftTransform->Translate(surfaceShift);
	TransformFilterType surfaceShiftFilter = TransformFilterType::New();
	surfaceShiftFilter->SetTransform(shiftTransform);
	surfaceShiftFilter->SetInput(HVCSurface);
	surfaceShiftFilter->Update();
	
	ClosedSurface * closedHVCSurface = new ClosedSurface(surfaceShiftFilter->GetOutput());
	double voxelVolume = voxelSize*voxelSize*voxelSize;
	double HVCBounds[6];
	closedHVCSurface->ptr()->GetBounds(HVCBounds);
// 	closedHVCSurface->ptr()->Print(std::cout);
	
	PointsPointerType wedge0Points = PointsPointerType::New();
	wedge0Points->SetDataTypeToFloat();
	wedge0Points->Allocate(1);
	
	for(int i = 0; i < wedgeSurfaces.size(); ++i)
	{
		std::cout << "Checking wedge " << i << std::endl;
		double wedgeBounds[6];
		wedgeSurfaces[i]->ptr()->GetBounds(wedgeBounds);
// 		if(i == 0)
// 		{
// 			wedgeSurfaces[i]->ptr()->Print(std::cout);
// 		}
		for(double x = HVCBounds[0]; x <= HVCBounds[1]; x += voxelSize)
			for(double y = HVCBounds[2]; y <= HVCBounds[3]; y += voxelSize)
				for(double z = HVCBounds[4]; z <= HVCBounds[5]; z += voxelSize)
				{
					std::cout << "Checking point X = " << x << " Y = " << y << " Z = " << z << "\r";
					
// 					double testPt[] = {32,-166,-227};
					double xyz[3];
					xyz[0] = x;
					xyz[1] = y;
					xyz[2] = z;
// 					bool testPtNow = false;
// 					if(i == 0)
// 					{
// 						double testDist = sqrt(vtkMath::Distance2BetweenPoints(testPt, xyz));
// 						if(testDist < 10)
// 						{
// 							std::cout << "Checking test point" << std::endl;
// 							testPtNow = true;
// 						}
// 					}
					
					if(x < wedgeBounds[0] || x > wedgeBounds[1] || y < wedgeBounds[2] || y > wedgeBounds[3] || z < wedgeBounds[4] || z > wedgeBounds[5])
					{
// 						if(testPtNow)
// 						{
// 							std::cout << "Test point outside of wedge bounds!?" << std::endl;
// 						}
						continue;
					}
					if(closedHVCSurface->isPointInsideSurface(xyz) && wedgeSurfaces[i]->isPointInsideSurface(xyz))
					{
						if(i == 0)
						{
							wedge0Points->InsertNextPoint(xyz);
						}
						wedgeVolumes[i] += voxelVolume;
						double dist = sqrt(x*x + y*y + z*z);
						if(dist < proximalBorder)
						{
							wedgeVolumesProximal[i] += voxelVolume;
						}
						int bin = int(dist/proximalBorder);
						wedgeDistanceVolumes[i]->addSegment(voxelVolume, bin);
					}
				}
		std::cout << std::endl;
	}
	
	// write output files
	
	std::string volumeOutName(outputFilename);
	volumeOutName += "_HVC_wedge_overlap_volumes.csv";
	std::ofstream VolumeFile;
	VolumeFile.open(volumeOutName.c_str());
	VolumeFile << "Wedge\tTotal volume (cubic microns)\tProximal volume (cubic microns)\tDistal volume (cubic microns)" << std::endl;
	for(int i = 0; i < wedgeVolumes.size(); ++i)
	{
		VolumeFile << i << "\t" << wedgeVolumes[i] << "\t" << wedgeVolumesProximal[i] << "\t" << wedgeVolumes[i] - wedgeVolumesProximal[i] << std::endl;
	}
	VolumeFile.close();
	
	int maxBins = 0;
	for(int i = 0; i < wedgeSurfaces.size(); ++i)
	{
		std::vector< double > * profilePtr = wedgeDistanceVolumes[i]->getProfile();
		if(profilePtr->size() > maxBins)
		{
			maxBins = profilePtr->size();
		}
	}
	std::string volumeOutName2(outputFilename);
	volumeOutName2 += "_HVC_wedge_overlap_distance_volumes.csv";
	std::ofstream VolumeFile2;
	VolumeFile2.open(volumeOutName2.c_str());
	VolumeFile2 << "Wedge\t";
	for(int i = 0; i < maxBins; ++i)
	{
		VolumeFile2 << i*proximalBorder << "\t";
	}
	VolumeFile2 << std::endl;
	for(int i = 0; i < wedgeDistanceVolumes.size(); ++i)
	{
		VolumeFile2 << i << "\t";
		std::vector< double > * profilePtr = wedgeDistanceVolumes[i]->getProfile();
		int binCount = 0;
		for(int j = 0; j < profilePtr->size(); ++j)
		{
			VolumeFile2 << profilePtr->at(j) << "\t";
			++binCount;
		}
		for(int j = 0; j < maxBins - binCount; ++j)
		{
			VolumeFile2 << "0\t";
		}
		VolumeFile2 << std::endl;
	}
	VolumeFile2.close();
	
// 	std::string wedgePointsName(outputFilename);
// 	wedgePointsName += "_wedge0Points.landmarkAscii";
// 	Reader * pointWriter = new Reader(wedgePointsName.c_str(), wedgePointsName.c_str());
// 	pointWriter->writeLandmarkFile(wedge0Points);
// 	delete pointWriter;
	
	delete closedHVCSurface;
}

void computeLowpassFilteredAngularDensity(AmiraSpatialGraph* neuronMorphology, PolyDataPointerType HVCSurface, std::vector< ClosedSurface* > wedgeSurfaces, int orientation, double voxelSize, const char* outputFilename)
{
	const double proximalBorder = 200; // microns
	double degreeSteps = 360;
// 	int degreeStepSize = 10;
	double degreeStepSize = 22.5;
	double degreeOffset = -28.5;
	std::vector< double > wedgeVolumes;
	std::vector< double > wedgeVolumesDistal;
	std::vector< double > wedgeLengths;
	std::vector< double > wedgeLengthsDistal;
	std::vector< double > wedgeDensities;
	std::vector< double > wedgeDensitiesDistal;
// 	for(int i = 0; i < degreeSteps; i += degreeStepSize)
// 	{
// 		wedgeDensities.push_back(0);
// 		wedgeDensitiesDistal.push_back(0);
// // 		if(i == 0)
// // 		{
// // 			wedgeSurfaces[i]->ptr()->Print(std::cout);
// // 		}
// 	}
	
	double somaLoc[3];
	PolyDataPointerType thisSoma = PolyDataPointerType::New();
	if(neuronMorphology->extractLandmark(Soma, thisSoma))
	{
		int subId;
		double pCenter[3], * weights = new double[thisSoma->GetNumberOfPoints()];
		thisSoma->GetCell(0)->GetParametricCenter(pCenter);
		thisSoma->GetCell(0)->EvaluateLocation(subId, pCenter, somaLoc, weights);
		delete [] weights;
	}
	else
	{
		std::cout << "Error! SpatialGraph file does not contain a soma!" << std::endl;
		return;
	}
	
	double surfaceShift[3];
	surfaceShift[0] = -somaLoc[0];
	surfaceShift[1] = -somaLoc[1];
	surfaceShift[2] = -somaLoc[2];
	
	TransformPointerType shiftTransform = TransformPointerType::New();
	shiftTransform->Translate(surfaceShift);
	TransformFilterType surfaceShiftFilter = TransformFilterType::New();
	surfaceShiftFilter->SetTransform(shiftTransform);
	surfaceShiftFilter->SetInput(HVCSurface);
	surfaceShiftFilter->Update();
	neuronMorphology->setTransformation(shiftTransform);
	neuronMorphology->applyTransformation();
	
	double maxDist = 0;
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label != Axon)
			continue;
		
		std::list< double * >::const_iterator ptIt;
		std::list< double * > tmpEdgePointList((*edgeIt)->edgePointCoordinates);
		ptIt = tmpEdgePointList.begin();
		while(ptIt != tmpEdgePointList.end())
		{
// #ifdef DEBUG
// 			std::flush(std::cout << "Processing point on edge " << *edgeIt << " @ [" << (*ptIt1)[0] << "," << (*ptIt1)[1] << "," << (*ptIt1)[2] << "]\r");
// #endif
			double * pt1 = *ptIt;
			double tmpDist;
			if(orientation == HORIZONTAL)
			{
				tmpDist = sqrt(pt1[0]*pt1[0] + pt1[1]*pt1[1]);
			}
			else if(orientation == CORONAL)
			{
				tmpDist = sqrt(pt1[0]*pt1[0] + pt1[2]*pt1[2]);
			}
			else if(orientation == SAGITTAL)
			{
				tmpDist = sqrt(pt1[1]*pt1[1] + pt1[2]*pt1[2]);
			}
			if(tmpDist > maxDist)
			{
				maxDist = tmpDist;
			}
			++ptIt;
		}
	}
	std::cout << "Maximum 2D-projected radial extent of axon morphology: " << maxDist << std::endl;
	
	ClosedSurface * closedHVCSurface = new ClosedSurface(surfaceShiftFilter->GetOutput());
	double voxelVolume = voxelSize*voxelSize*voxelSize;
	double HVCBounds[6];
	closedHVCSurface->ptr()->GetBounds(HVCBounds);
// 	closedHVCSurface->ptr()->Print(std::cout);
	
// 	PointsPointerType debugPoints = PointsPointerType::New();
// 	debugPoints->SetDataTypeToFloat();
// 	debugPoints->Allocate(1);
	
// 	AmiraSpatialGraph * debugSG = new AmiraSpatialGraph;
	
// 	for(int i = 0; i < 1; ++i)
// 	for(int i = 187; i < 188; ++i)
	for(double i = 0; i < degreeSteps; i += degreeStepSize)
	{
		double tmpAngle = i + degreeOffset;
		std::cout << "Checking wedge " << tmpAngle << std::endl;
		
		PolyDataPointerType tmpWedge = wedgeSurfaces[0]->ptr();
		TransformPointerType tmpRotation = TransformPointerType::New();
		if(orientation == HORIZONTAL)
		{
			tmpRotation->RotateZ(tmpAngle);
		}
		else if(orientation == CORONAL)
		{
			tmpRotation->RotateY(tmpAngle);
		}
		else if(orientation == SAGITTAL)
		{
			tmpRotation->RotateX(tmpAngle);
		}
		TransformFilterType tmpRotateFilter = TransformFilterType::New();
		tmpRotateFilter->SetInput(tmpWedge);
		tmpRotateFilter->SetTransform(tmpRotation);
		tmpRotateFilter->Update();
		ClosedSurface * tmpRotatingWedge = new ClosedSurface(tmpRotateFilter->GetOutput());
		double wedgeBounds[6];
		tmpRotatingWedge->ptr()->GetBounds(wedgeBounds);
// 		tmpRotatingWedge->ptr()->Print(std::cout);
// 		if(i == 0)
// 		{
// 			wedgeSurfaces[i]->ptr()->Print(std::cout);
// 		}
		
		double tmpVolume = 0;
		double tmpDistalVolume = 0;
		double tmpLength = 0;
		double tmpDistalLength = 0;
		
// 		for(double x = HVCBounds[0]; x <= HVCBounds[1]; x += voxelSize)
// 			for(double y = HVCBounds[2]; y <= HVCBounds[3]; y += voxelSize)
// 				for(double z = HVCBounds[4]; z <= HVCBounds[5]; z += voxelSize)
// 				{
// 					std::cout << "Checking point X = " << x << " Y = " << y << " Z = " << z << "\r";
// 					
// // 					double testPt[] = {32,-166,-227};
// 					double xyz[3];
// 					xyz[0] = x;
// 					xyz[1] = y;
// 					xyz[2] = z;
// // 					bool testPtNow = false;
// // 					if(i == 0)
// // 					{
// // 						double testDist = sqrt(vtkMath::Distance2BetweenPoints(testPt, xyz));
// // 						if(testDist < 10)
// // 						{
// // 							std::cout << "Checking test point" << std::endl;
// // 							testPtNow = true;
// // 						}
// // 					}
// 					
// 					if(x < wedgeBounds[0] || x > wedgeBounds[1] || y < wedgeBounds[2] || y > wedgeBounds[3] || z < wedgeBounds[4] || z > wedgeBounds[5])
// 					{
// // 						if(testPtNow)
// // 						{
// // 							std::cout << "Test point outside of wedge bounds!?" << std::endl;
// // 						}
// 						continue;
// 					}
// 					if(closedHVCSurface->isPointInsideSurface(xyz) && tmpRotatingWedge->isPointInsideSurface(xyz))
// 					{
// // 						if(i == 0)
// // 						{
// // 							wedge0Points->InsertNextPoint(xyz);
// // 						}
// 						double dist2D;
// 						if(orientation == HORIZONTAL)
// 						{
// 							dist2D = sqrt(x*x + y*y);
// 						}
// 						else if(orientation == CORONAL)
// 						{
// 							dist2D = sqrt(x*x + z*z);
// 						}
// 						else if(orientation == SAGITTAL)
// 						{
// 							dist2D = sqrt(y*y + z*z);
// 						}
// 						if(dist2D < maxDist)
// 						{
// 							tmpVolume += voxelVolume;
// // 							double dist = sqrt(x*x + y*y + z*z);
// 							if(dist2D > proximalBorder)
// 							{
// 								tmpDistalVolume += voxelVolume;
// 							}
// 						}
// 					}
// 				}
// 		std::cout << std::endl;
		wedgeVolumes.push_back(tmpVolume);
		wedgeVolumesDistal.push_back(tmpDistalVolume);
		
		// main loop
		for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label != Axon)
				continue;
			
			std::list< double * >::const_iterator ptIt1, ptIt2;
			std::list< double * > tmpEdgePointList((*edgeIt)->edgePointCoordinates);
			ptIt1 = tmpEdgePointList.begin();
			ptIt2 = tmpEdgePointList.begin();
			++ptIt1;
			while(ptIt1 != tmpEdgePointList.end())
			{
// 	#ifdef DEBUG
// 				std::flush(std::cout << "Processing point on edge " << *edgeIt << " @ [" << (*ptIt1)[0] << "," << (*ptIt1)[1] << "," << (*ptIt1)[2] << "]\r");
// 	#endif
				double * pt1, * pt2;
				int pt1WedgeIndex = -1;
				int pt2WedgeIndex = -1;
				
				pt1 = *ptIt1;
				pt2 = *ptIt2;
				
				double dist1, dist2;
				if(orientation == HORIZONTAL)
				{
					dist1 = sqrt(pt1[0]*pt1[0] + pt1[1]*pt1[1]);
					dist2 = sqrt(pt2[0]*pt2[0] + pt2[1]*pt2[1]);
				}
				if(orientation == CORONAL)
				{
					dist1 = sqrt(pt1[0]*pt1[0] + pt1[2]*pt1[2]);
					dist2 = sqrt(pt2[0]*pt2[0] + pt2[2]*pt2[2]);
				}
				if(orientation == SAGITTAL)
				{
					dist1 = sqrt(pt1[1]*pt1[1] + pt1[2]*pt1[2]);
					dist2 = sqrt(pt2[1]*pt2[1] + pt2[2]*pt2[2]);
				}
				
				// if points on different sides of proximal border:
				// virtually insert point 
				if((dist1 <= proximalBorder && dist2 > proximalBorder) || (dist1 > proximalBorder && dist2 <= proximalBorder))
				{
					// TBD - for now just lump to proximal point
				}
				
				if(tmpRotatingWedge->isPointInsideSurface(pt1))
				{
					pt1WedgeIndex = 1;
// 					std::cout << "pt1WedgeIndex = " << pt1WedgeIndex << std::endl;
				}
				if(tmpRotatingWedge->isPointInsideSurface(pt2))
				{
					pt2WedgeIndex = 1;
// 					std::cout << "pt2WedgeIndex = " << pt2WedgeIndex << std::endl;
				}
				
// 				#ifdef DEBUG
// 				assignDistanceToLayers(pt1, pt2, pt1Layer, pt2Layer, onlyLayerProfiles, onlyLayerProfiles);
// 				#endif
				
// 				// sanity check
// 				if(pt1WedgeIndex == -1 || pt2WedgeIndex == -1)
// 				{
// 					std::cout << "Error: Point outside of wedges: " << std::endl;
// 					std::cout << "pt1WedgeIndex = " << pt1WedgeIndex << std::endl;
// 					std::cout << "pt2WedgeIndex = " << pt2WedgeIndex << std::endl;
// 					std::cout << "pt1  @ [" << pt1[0] << "," << pt1[1] << "," << pt1[2] << "]" << std::endl;
// 					std::cout << "pt2  @ [" << pt2[0] << "," << pt2[1] << "," << pt2[2] << "]" << std::endl;
// 					
// 					++ptIt1, ++ptIt2;
// 					continue;
// 				}
				
				// case 0: both not in wedge:
				if(pt1WedgeIndex == -1 && pt2WedgeIndex == -1)
				{
					// do nothing
				}
				// case 1 : both in same wedge
				else if(pt1WedgeIndex == 1 && pt2WedgeIndex == 1)
				{
					double length = sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
					if(dist1 > proximalBorder)
					{
						tmpDistalLength += length;
					}
					tmpLength += length;
// 					if(i==187)
// 					{
// // 						debugPoints->InsertNextPoint(pt1);
// // 						debugPoints->InsertNextPoint(pt2);
// 					}
				}
				// case 2: points in neighboring wedges (assume that axon sampling distance is sufficiently small)
				else
				{
					double line12[3];
					tmpRotatingWedge->intersectLineSegment(pt1, pt2);
					double intersectPt[3];
					if(tmpRotatingWedge->isIntersectionFound())
					{
						tmpRotatingWedge->getLastIntersectPoint(intersectPt);
					}
					else
					{
						std::cout << "Error: Could not find intersection point with wedge " << pt1WedgeIndex << std::endl;
						++ptIt1, ++ptIt2;
						continue;
					}
					
					double distance1 = sqrt(vtkMath::Distance2BetweenPoints(pt1, intersectPt));
					double distance2 = sqrt(vtkMath::Distance2BetweenPoints(pt2, intersectPt));
					if(pt1WedgeIndex == 1)
					{
						if(dist1 > proximalBorder)
						{
							tmpDistalLength += distance1;
						}
						tmpLength += distance1;
// 						if(i==187)
// 						{
// // 							debugPoints->InsertNextPoint(pt1);
// 							debugPoints->InsertNextPoint(intersectPt);
// 							double * v1Pt = new double[3];
// 							double * v2Pt = new double[3];
// 							double * intPt = new double[3];
// 							for(int k = 0; k < 3; ++k)
// 							{
// 								v1Pt[k] = pt1[k];
// 								v2Pt[k] = pt2[k];
// 								intPt[k] = intersectPt[k];
// 							}
// 							Vertex * newVertex1 = new Vertex(v1Pt, Axon);
// 							Vertex * newVertex2 = new Vertex(v2Pt, Axon);
// 							std::list< double * > edgePtCoords;
// 							edgePtCoords.push_back(v1Pt);
// 							edgePtCoords.push_back(intPt);
// 							edgePtCoords.push_back(v2Pt);
// 							int * edgeConnectivity = new int[2];
// 							edgeConnectivity[0] = debugSG->getNumberOfVertices();
// 							edgeConnectivity[1] = debugSG->getNumberOfVertices() + 1;
// 							Edge * newEdge = new Edge(edgeConnectivity, 3, Axon, edgePtCoords, 1);
// 							debugSG->addVertex(newVertex1);
// 							debugSG->addVertex(newVertex2);
// 							debugSG->addEdge(newEdge);
// 						}
					}
					if(pt2WedgeIndex == 1)
					{
						if(dist2 > proximalBorder)
						{
							tmpDistalLength += distance2;
						}
						tmpLength += distance2;
// 						if(i==187)
// 						{
// // 							debugPoints->InsertNextPoint(pt2);
// 							debugPoints->InsertNextPoint(intersectPt);
// 							double * v1Pt = new double[3];
// 							double * v2Pt = new double[3];
// 							double * intPt = new double[3];
// 							for(int k = 0; k < 3; ++k)
// 							{
// 								v1Pt[k] = pt1[k];
// 								v2Pt[k] = pt2[k];
// 								intPt[k] = intersectPt[k];
// 							}
// 							Vertex * newVertex1 = new Vertex(v1Pt, Axon);
// 							Vertex * newVertex2 = new Vertex(v2Pt, Axon);
// 							std::list< double * > edgePtCoords;
// 							edgePtCoords.push_back(v1Pt);
// 							edgePtCoords.push_back(intPt);
// 							edgePtCoords.push_back(v2Pt);
// 							int * edgeConnectivity = new int[2];
// 							edgeConnectivity[0] = debugSG->getNumberOfVertices();
// 							edgeConnectivity[1] = debugSG->getNumberOfVertices() + 1;
// 							Edge * newEdge = new Edge(edgeConnectivity, 3, Axon, edgePtCoords, 1);
// 							debugSG->addVertex(newVertex1);
// 							debugSG->addVertex(newVertex2);
// 							debugSG->addEdge(newEdge);
// 						}
					}
				}
				
				++ptIt1, ++ptIt2;
			}
		}
		wedgeLengths.push_back(tmpLength);
		wedgeLengthsDistal.push_back(tmpDistalLength);
		
		double tmpDensity = -1;
		double tmpDistalDensity = -1;
		if(tmpVolume > 0)
		{
			tmpDensity = tmpLength/tmpVolume;
			std::cout << "Non-zero volume, calculating density: " << tmpDensity << std::endl;
		}
		if(tmpDistalVolume > 0)
		{
			tmpDistalDensity = tmpDistalLength/tmpDistalVolume;
		}
		wedgeDensities.push_back(tmpDensity);
		wedgeDensitiesDistal.push_back(tmpDistalDensity);
		
		delete tmpRotatingWedge;
	}
	
	// write output files
	
	std::string volumeOutName(outputFilename);
	volumeOutName += "_HVC_angular_density.csv";
	std::ofstream VolumeFile;
	VolumeFile.open(volumeOutName.c_str());
	VolumeFile << "Degree\tLength\tDistal length\tVolume\tDensity [microns/cubic microns]" << std::endl;
	for(int i = 0; i < wedgeLengths.size(); ++i)
	{
		VolumeFile << i << "\t" << wedgeLengths[i] << "\t" << wedgeLengthsDistal[i] << "\t" << wedgeVolumes[i] << "\t" << wedgeDensities[i] << std::endl;
	}
	VolumeFile.close();
	
// 	std::string volumeOutName(outputFilename);
// 	volumeOutName += "_HVC_lowpass_filtered_angular_density.csv";
// 	std::ofstream VolumeFile;
// 	VolumeFile.open(volumeOutName.c_str());
// 	VolumeFile << "Degree\tLength\tVolume\tDensity [microns/cubic microns]" << std::endl;
// 	for(int i = 0; i < wedgeLengths.size(); ++i)
// 	{
// 		VolumeFile << i << "\t" << wedgeLengths[i] << "\t" << wedgeVolumes[i] << "\t" << wedgeDensities[i] << std::endl;
// 	}
// 	VolumeFile.close();
	
// 	std::string volumeOutName2(outputFilename);
// 	volumeOutName2 += "_HVC_lowpass_filtered_angular_density_distal.csv";
// 	std::ofstream VolumeFile2;
// 	VolumeFile2.open(volumeOutName2.c_str());
// 	VolumeFile2 << "Degree\tLength\tVolume\tDensity [microns/cubic microns]" << std::endl;
// 	for(int i = 0; i < wedgeDensitiesDistal.size(); ++i)
// 	{
// 		VolumeFile2 << i << "\t" << wedgeLengthsDistal[i] << "\t" << wedgeVolumesDistal[i] << "\t" << wedgeDensitiesDistal[i] << std::endl;
// 	}
// 	VolumeFile2.close();
	
// 	std::string wedgePointsName(outputFilename);
// 	wedgePointsName += "_wedge187Points.landmarkAscii";
// 	Reader * pointWriter = new Reader(wedgePointsName.c_str(), wedgePointsName.c_str());
// 	pointWriter->writeLandmarkFile(debugPoints);
// 	delete pointWriter;
// 	
// 	std::string debugSGName(outputFilename);
// 	debugSGName += "_debugSG.am";
// 	Reader * debugSGWriter = new Reader(debugSGName.c_str(), debugSGName.c_str());
// 	debugSGWriter->setSpatialGraph(debugSG);
// 	pointWriter->writeSpatialGraphFile();
// 	delete pointWriter;
	
	delete closedHVCSurface;
}

void computeLowpassFilteredAngularDensity2(AmiraSpatialGraph* neuronMorphology, PolyDataPointerType HVCSurface, std::vector< ClosedSurface* > wedgeSurfaces, int orientation, double voxelSize, const char* outputFilename)
{
	const double proximalBorder = 200; // microns
	int degreeSteps = 360;
	int degreeStepSize = 10;
	double degreeOffset = -28.5;
	std::vector< double > wedgeVolumes;
	std::vector< double > wedgeVolumesDistal;
	std::vector< double > wedgeLengths;
	std::vector< double > wedgeLengthsDistal;
	std::vector< double > wedgeDensities;
	std::vector< double > wedgeDensitiesDistal;
// 	for(int i = 0; i < degreeSteps; i += degreeStepSize)
// 	{
// 		wedgeDensities.push_back(0);
// 		wedgeDensitiesDistal.push_back(0);
// // 		if(i == 0)
// // 		{
// // 			wedgeSurfaces[i]->ptr()->Print(std::cout);
// // 		}
// 	}
	
	double somaLoc[3];
	PolyDataPointerType thisSoma = PolyDataPointerType::New();
	if(neuronMorphology->extractLandmark(Soma, thisSoma))
	{
		int subId;
		double pCenter[3], * weights = new double[thisSoma->GetNumberOfPoints()];
		thisSoma->GetCell(0)->GetParametricCenter(pCenter);
		thisSoma->GetCell(0)->EvaluateLocation(subId, pCenter, somaLoc, weights);
		delete [] weights;
	}
	else
	{
		std::cout << "Error! SpatialGraph file does not contain a soma!" << std::endl;
		return;
	}
	
	double surfaceShift[3];
	surfaceShift[0] = -somaLoc[0];
	surfaceShift[1] = -somaLoc[1];
	surfaceShift[2] = -somaLoc[2];
	
	TransformPointerType shiftTransform = TransformPointerType::New();
	shiftTransform->Translate(surfaceShift);
	TransformFilterType surfaceShiftFilter = TransformFilterType::New();
	surfaceShiftFilter->SetTransform(shiftTransform);
	surfaceShiftFilter->SetInput(HVCSurface);
	surfaceShiftFilter->Update();
	neuronMorphology->setTransformation(shiftTransform);
	neuronMorphology->applyTransformation();
	
	double maxDist = 0;
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label != Axon)
// 		if((*edgeIt)->label != Dendrite &&(*edgeIt)->label != BasalDendrite && (*edgeIt)->label != ApicalDendrite)
			continue;
		
		std::list< double * >::const_iterator ptIt;
		std::list< double * > tmpEdgePointList((*edgeIt)->edgePointCoordinates);
		ptIt = tmpEdgePointList.begin();
		while(ptIt != tmpEdgePointList.end())
		{
// #ifdef DEBUG
// 			std::flush(std::cout << "Processing point on edge " << *edgeIt << " @ [" << (*ptIt1)[0] << "," << (*ptIt1)[1] << "," << (*ptIt1)[2] << "]\r");
// #endif
			double * pt1 = *ptIt;
			double tmpDist;
			if(orientation == HORIZONTAL)
			{
				tmpDist = sqrt(pt1[0]*pt1[0] + pt1[1]*pt1[1]);
			}
			else if(orientation == CORONAL)
			{
				tmpDist = sqrt(pt1[0]*pt1[0] + pt1[2]*pt1[2]);
			}
			else if(orientation == SAGITTAL)
			{
				tmpDist = sqrt(pt1[1]*pt1[1] + pt1[2]*pt1[2]);
			}
			if(tmpDist > maxDist)
			{
				maxDist = tmpDist;
			}
			++ptIt;
		}
	}
	std::cout << "Maximum 2D-projected radial extent of axon morphology: " << maxDist << std::endl;
	
	ClosedSurface * closedHVCSurface = new ClosedSurface(surfaceShiftFilter->GetOutput());
	double voxelVolume = voxelSize*voxelSize*voxelSize;
	double HVCBounds[6];
	closedHVCSurface->ptr()->GetBounds(HVCBounds);
	for(int i = 0; i < degreeSteps; i += degreeStepSize)
	{
		double rotationAngle = i + degreeOffset;
		std::cout << "Checking wedge " << i << " at angle " << rotationAngle << std::endl;
		PolyDataPointerType tmpWedge = wedgeSurfaces[0]->ptr();
		TransformPointerType tmpRotation = TransformPointerType::New();
		if(orientation == HORIZONTAL)
		{
			tmpRotation->RotateZ(rotationAngle);
		}
		else if(orientation == CORONAL)
		{
			tmpRotation->RotateY(rotationAngle);
		}
		else if(orientation == SAGITTAL)
		{
			tmpRotation->RotateX(rotationAngle);
		}
		TransformFilterType tmpRotateFilter = TransformFilterType::New();
		tmpRotateFilter->SetInput(tmpWedge);
		tmpRotateFilter->SetTransform(tmpRotation);
		tmpRotateFilter->Update();
		ClosedSurface * tmpRotatingWedge = new ClosedSurface(tmpRotateFilter->GetOutput());
		double wedgeBounds[6];
		tmpRotatingWedge->ptr()->GetBounds(wedgeBounds);
// 		tmpRotatingWedge->ptr()->Print(std::cout);
// 		if(i == 0)
// 		{
// 			wedgeSurfaces[i]->ptr()->Print(std::cout);
// 		}
		
		double tmpVolume = 0;
		double tmpDistalVolume = 0;
		double tmpLength = 0;
		double tmpDistalLength = 0;
		
		for(double x = HVCBounds[0]; x <= HVCBounds[1]; x += voxelSize)
			for(double y = HVCBounds[2]; y <= HVCBounds[3]; y += voxelSize)
				for(double z = HVCBounds[4]; z <= HVCBounds[5]; z += voxelSize)
				{
					std::cout << "Checking point X = " << x << " Y = " << y << " Z = " << z << "\r";
					
// 					double testPt[] = {32,-166,-227};
					double xyz[3];
					xyz[0] = x + 0.5*voxelSize;
					xyz[1] = y + 0.5*voxelSize;
					xyz[2] = z + 0.5*voxelSize;
// 					bool testPtNow = false;
// 					if(i == 0)
// 					{
// 						double testDist = sqrt(vtkMath::Distance2BetweenPoints(testPt, xyz));
// 						if(testDist < 10)
// 						{
// 							std::cout << "Checking test point" << std::endl;
// 							testPtNow = true;
// 						}
// 					}
					
					if(x < wedgeBounds[0] || x > wedgeBounds[1] || y < wedgeBounds[2] || y > wedgeBounds[3] || z < wedgeBounds[4] || z > wedgeBounds[5])
					{
// 						if(testPtNow)
// 						{
// 							std::cout << "Test point outside of wedge bounds!?" << std::endl;
// 						}
						continue;
					}
					if(closedHVCSurface->isPointInsideSurface(xyz) && tmpRotatingWedge->isPointInsideSurface(xyz))
					{
// 						if(i == 0)
// 						{
// 							wedge0Points->InsertNextPoint(xyz);
// 						}
// 						double dist2D;
// 						if(orientation == HORIZONTAL)
// 						{
// 							dist2D = sqrt(x*x + y*y);
// 						}
// 						else if(orientation == CORONAL)
// 						{
// 							dist2D = sqrt(x*x + z*z);
// 						}
// 						else if(orientation == SAGITTAL)
// 						{
// 							dist2D = sqrt(y*y + z*z);
// 						}
						double dist = sqrt(x*x + y*y + z*z);
// 						if(dist2D < maxDist)
						if(dist < maxDist)
						{
							tmpVolume += voxelVolume;
// 							double dist = sqrt(x*x + y*y + z*z);
// 							if(dist2D > proximalBorder)
							if(dist > proximalBorder)
							{
								tmpDistalVolume += voxelVolume;
							}
							double clippingBox[6];
							clippingBox[0] = x;
							clippingBox[1] = x + voxelSize;
							clippingBox[2] = y;
							clippingBox[3] = y + voxelSize;
							clippingBox[4] = z;
							clippingBox[5] = z + voxelSize;
							AmiraSpatialGraph * clippedGraph = neuronMorphology->clipSpatialGraph(clippingBox);
							if(clippedGraph)
							{
								for(edgeIt = clippedGraph->edgesBegin(); edgeIt != clippedGraph->edgesEnd(); ++edgeIt)
								{
									if((*edgeIt)->label != Axon)
// 									if((*edgeIt)->label != Dendrite &&(*edgeIt)->label != BasalDendrite && (*edgeIt)->label != ApicalDendrite)
										continue;
									tmpLength += (*edgeIt)->segmentLength();
// 									if(dist2D > proximalBorder)
									if(dist > proximalBorder)
									{
										tmpDistalLength += (*edgeIt)->segmentLength();
									}
								}
							}
						}
					}
				}
		std::cout << std::endl;
		wedgeVolumes.push_back(tmpVolume);
		wedgeVolumesDistal.push_back(tmpDistalVolume);
		wedgeLengths.push_back(tmpLength);
		wedgeLengthsDistal.push_back(tmpDistalLength);
		double tmpDensity = -1, tmpDistalDensity = -1;
		if(tmpVolume > 0)
		{
			tmpDensity = tmpLength/tmpVolume;
		}
		if(tmpDistalVolume > 0)
		{
			tmpDistalDensity = tmpDistalLength/tmpDistalVolume;
		}
		wedgeDensities.push_back(tmpDensity);
		wedgeDensitiesDistal.push_back(tmpDistalDensity);
	}
	
	// write output files
	
	std::string volumeOutName(outputFilename);
	volumeOutName += "_HVC_lowpass_filtered_voxel_wedge_density.csv";
	std::ofstream VolumeFile;
	VolumeFile.open(volumeOutName.c_str());
	VolumeFile << "Degree\tLength\tVolume\tDensity [microns/cubic microns]" << std::endl;
	for(int i = 0; i < wedgeLengths.size(); ++i)
	{
		VolumeFile << i << "\t" << wedgeLengths[i] << "\t" << wedgeVolumes[i] << "\t" << wedgeDensities[i] << std::endl;
	}
	VolumeFile.close();
	
	std::string volumeOutName2(outputFilename);
	volumeOutName2 += "_HVC_lowpass_filtered_voxel_wedge_density_distal.csv";
	std::ofstream VolumeFile2;
	VolumeFile2.open(volumeOutName2.c_str());
	VolumeFile2 << "Degree\tLength\tVolume\tDensity [microns/cubic microns]" << std::endl;
	for(int i = 0; i < wedgeDensitiesDistal.size(); ++i)
	{
		VolumeFile2 << i << "\t" << wedgeLengthsDistal[i] << "\t" << wedgeVolumesDistal[i] << "\t" << wedgeDensitiesDistal[i] << std::endl;
	}
	VolumeFile2.close();
}

void computeAxonInsideOutsideSphere(AmiraSpatialGraph* neuronMorphology, double proximalBorder, const char* outputFilename)
{
	double somaLoc[3];
	PolyDataPointerType thisSoma = PolyDataPointerType::New();
	if(neuronMorphology->extractLandmark(Soma, thisSoma))
	{
		int subId;
		double pCenter[3], * weights = new double[thisSoma->GetNumberOfPoints()];
		thisSoma->GetCell(0)->GetParametricCenter(pCenter);
		thisSoma->GetCell(0)->EvaluateLocation(subId, pCenter, somaLoc, weights);
		delete [] weights;
	}
	else
	{
		std::cout << "Error! SpatialGraph file does not contain a soma!" << std::endl;
		return;
	}
	
	double proximalLength = 0;
	double distalLength = 0;
	// main loop
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label != Axon && (*edgeIt)->label != ProjectionAxon)
			continue;
		
		std::list< double * >::const_iterator ptIt1, ptIt2;
		std::list< double * > tmpEdgePointList((*edgeIt)->edgePointCoordinates);
		ptIt1 = tmpEdgePointList.begin();
		ptIt2 = tmpEdgePointList.begin();
		++ptIt1;
		while(ptIt1 != tmpEdgePointList.end())
		{
// #ifdef DEBUG
// 			std::flush(std::cout << "Processing point on edge " << *edgeIt << " @ [" << (*ptIt1)[0] << "," << (*ptIt1)[1] << "," << (*ptIt1)[2] << "]\r");
// #endif
			double * pt1, * pt2;
			pt1 = *ptIt1;
			pt2 = *ptIt2;
			
			double dist1, dist2;
			dist1 = sqrt((pt1[0]-somaLoc[0])*(pt1[0]-somaLoc[0]) + (pt1[1]-somaLoc[1])*(pt1[1]-somaLoc[1]) + (pt1[2]-somaLoc[2])*(pt1[2]-somaLoc[2]));
			dist2 = sqrt((pt2[0]-somaLoc[0])*(pt2[0]-somaLoc[0]) + (pt2[1]-somaLoc[1])*(pt2[1]-somaLoc[1]) + (pt2[2]-somaLoc[2])*(pt2[2]-somaLoc[2]));
			
			// if points on different sides of proximal border:
			// virtually insert point on sphere
			if(dist1 < proximalBorder && dist2 >= proximalBorder)
			{
				double p2p1Diff[3], p1SomaDiff[3];
				double d12Square, alpha, r2, t, p1Square;
				vtkMath::Subtract(pt2, pt1, p2p1Diff);
				vtkMath::Subtract(pt1, somaLoc, p1SomaDiff);
				d12Square = vtkMath::Distance2BetweenPoints(pt2, pt1);
				r2 = proximalBorder*proximalBorder;
				alpha = vtkMath::Dot(p1SomaDiff, p2p1Diff);
				p1Square = vtkMath::Distance2BetweenPoints(pt1, somaLoc);
				t = -alpha/d12Square + sqrt(alpha*alpha/(d12Square*d12Square) + (r2 - p1Square)/d12Square);
#ifdef DEBUG
				std::flush(std::cout << "d12Square = " << d12Square << std::endl);
				std::flush(std::cout << "r2 = " << r2 << std::endl);
				std::flush(std::cout << "alpha = " << alpha << std::endl);
				std::flush(std::cout << "p1Square = " << p1Square << std::endl);
				std::flush(std::cout << "pt1 inside, pt2 outside - t = " << t << std::endl);
#endif
				
				double intersectPt[3];
				for(int i = 0; i < 3; ++i)
				{
					intersectPt[i] = pt1[i] + t*p2p1Diff[i];
				}
				
				proximalLength += sqrt(vtkMath::Distance2BetweenPoints(pt1, intersectPt));
				distalLength += sqrt(vtkMath::Distance2BetweenPoints(intersectPt, pt2));
			}
			else if(dist1 >= proximalBorder && dist2 < proximalBorder)
			{
				double p1p2Diff[3], p2SomaDiff[3];
				double d12Square, alpha, r2, t, p2Square;
				vtkMath::Subtract(pt1, pt2, p1p2Diff);
				vtkMath::Subtract(pt2, somaLoc, p2SomaDiff);
				d12Square = vtkMath::Distance2BetweenPoints(pt2, pt1);
				r2 = proximalBorder*proximalBorder;
				alpha = vtkMath::Dot(p2SomaDiff, p1p2Diff);
				p2Square = vtkMath::Distance2BetweenPoints(pt2, somaLoc);
				t = -alpha/d12Square + sqrt(alpha*alpha/(d12Square*d12Square) + (r2 - p2Square)/d12Square);
#ifdef DEBUG
				std::flush(std::cout << "d12Square = " << d12Square << std::endl);
				std::flush(std::cout << "r2 = " << r2 << std::endl);
				std::flush(std::cout << "alpha = " << alpha << std::endl);
				std::flush(std::cout << "p2Square = " << p2Square << std::endl);
				std::flush(std::cout << "pt2 inside, pt1 outside - t = " << t << std::endl);
#endif
				
				double intersectPt[3];
				for(int i = 0; i < 3; ++i)
				{
					intersectPt[i] = pt2[i] + t*p1p2Diff[i];
				}
				
				proximalLength += sqrt(vtkMath::Distance2BetweenPoints(pt2, intersectPt));
				distalLength += sqrt(vtkMath::Distance2BetweenPoints(intersectPt, pt1));
			}
			else if(dist1 >= proximalBorder && dist2 >= proximalBorder)
			{
				distalLength += sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
			}
			else if(dist1 < proximalBorder && dist2 < proximalBorder)
			{
				proximalLength += sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
			}
			
			++ptIt1, ++ptIt2;
		}
	}
	// write output files
	
	char * outName = new char[128];
	sprintf(outName, "%s_proximal_distal_%.1f.csv", outputFilename, proximalBorder);
	std::ofstream ProxDistFile;
	ProxDistFile.open(outName);
	ProxDistFile << "Proximal length (microns)\tDistal length (microns)" << std::endl;
	ProxDistFile << proximalLength << "\t" << distalLength << std::endl;
	ProxDistFile.close();
	
	std::cout << outputFilename << "\t" << proximalLength << "\t" << distalLength << std::endl;
}





