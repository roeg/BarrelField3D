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
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/amiraReader.h"
#include "../../common/inputcheckpoint.h"
#include <gsl/gsl_rng.h>

double RASynapseProbability(double somaDistance);

int main( int argc , char * argv[])
{
	if(argc == 2)
	{
		const char * inputFilename = argv[1];
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
		
		std::string ofName(ifName, 0, suffix);
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		
		AmiraSpatialGraph * neuronMorphology = hocReader->getSpatialGraph();
        //std::cout << "Number of vertices: " << neuronMorphology->getNumberOfVertices() << std::endl;
        //std::cout << "Number of edges: " << neuronMorphology->getNumberOfEdges() << std::endl;
        //std::cout << "Number of points: " << neuronMorphology->getNumberOfPoints() << std::endl;
		
		if(checkPoint->getParameters().axonFlag)
		{
            
            //Sam's cells with only one vertex at soma location:
            double somaCenter[3];
        	std::vector< Vertex * >::const_iterator vertexIt;
        	for(vertexIt = neuronMorphology->verticesBegin(); vertexIt != neuronMorphology->verticesEnd(); ++vertexIt)
        	{
                if((*vertexIt)->label == Soma)
                {
                    somaCenter[0] = (*vertexIt)->coordinates[0];
                    somaCenter[1] = (*vertexIt)->coordinates[1];
                    somaCenter[2] = (*vertexIt)->coordinates[2];
                    break;
                }
        	}
            //PolyDataPointerType structure = PolyDataPointerType::New();
        	//if(!neuronMorphology->extractLandmark(Soma, structure))
        	//{
        	//	std::cout << "Error! Could not find structure with ID Soma in SpatialGraph!" << std::endl;
        	//	return 0;
        	//}
        	//int subID;
        	//double pCoords[3], * weights;
        	//weights = new double[structure->GetCell(0)->GetNumberOfPoints()];
        	//structure->GetCell(0)->GetParametricCenter(pCoords);
        	//structure->GetCell(0)->EvaluateLocation(subID, pCoords, somaCenter, weights);
        	//delete [] weights;
        	
            std::vector< double > segmentLengths;
            std::vector< double > segmentMidpointDistances;
			std::vector< double > segmentHorizontalAngles;
        	// start iteration through all edges
        	std::vector< Edge * >::const_iterator edgeIt;
        	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
        	{
        		if((*edgeIt)->label != Axon)
        			continue;
        		
                double tmpLength = 0.0;
        		std::list< double * >::const_iterator ptIt1, ptIt2;
        		ptIt1 = (*edgeIt)->edgePointCoordinates.begin();
        		ptIt2 = ptIt1;
        		++ptIt1;
        		while(ptIt1 != (*edgeIt)->edgePointCoordinates.end())
        		{
                    double * pt1 = *ptIt1;
                    double * pt2 = *ptIt2;
                    tmpLength += sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
        			++ptIt1, ++ptIt2;
        		}
                segmentLengths.push_back(tmpLength);
        		
                double centerPt[3];
                double * beginNode = (*edgeIt)->edgePointCoordinates.front();
                double * endNode = (*edgeIt)->edgePointCoordinates.back();
                vtkMath::Add(beginNode, endNode, centerPt);
                vtkMath::MultiplyScalar(centerPt, 0.5);
                double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(centerPt, somaCenter));
                segmentMidpointDistances.push_back(tmpDist);
				
				double diffVev[3];
				vtkMath::Subtract(endNode, beginNode, diffVev);
				double segmentAngle = 180/PI*atan2(diffVev[1], diffVev[0]) - 90;
				if(segmentAngle < 0)
				{
					segmentAngle += 360;
				}
				if(segmentAngle > 180)
				{
					segmentAngle -= 180;
				}
				segmentHorizontalAngles.push_back(segmentAngle);
    		}
    		
//             std::string lengthOutName(ofName);
//             lengthOutName += "_length_distance_list.csv";
//         	std::ofstream Table;
//         	Table.open(lengthOutName.c_str());
//         	Table << "Segment length (microns)\tSegment midpoint-soma distance (microns)\n";
//             for(int i = 0; i < segmentLengths.size(); ++i)
//             {
//                 Table << segmentLengths[i] << "\t" << segmentMidpointDistances[i] << std::endl;
//             }
//             Table.close();
    		
            std::string angleOutName(ofName);
            angleOutName += "_segment_angle_horizontal_list.csv";
        	std::ofstream Table2;
        	Table2.open(angleOutName.c_str());
        	Table2 << "Segment angle horizontal (degrees)\tsegment length (microns)\n";
            for(int i = 0; i < segmentHorizontalAngles.size(); ++i)
            {
                Table2 << segmentHorizontalAngles[i] << "\t" << segmentLengths[i] << std::endl;
            }
            Table2.close();
		}
		else
		{
            std::cout << "Error: Morphology does not have any structure with label Axon or Soma" << std::endl;
		}
		
		delete checkPoint;
		delete hocReader;
	}
	
	else if(argc == 3)
	{
		const char * inputFilename = argv[1];
		double meanSynapseDistance = atof(argv[2]);
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
		
		std::string ofName(ifName, 0, suffix);
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		
		AmiraSpatialGraph * neuronMorphology = hocReader->getSpatialGraph();
        //std::cout << "Number of vertices: " << neuronMorphology->getNumberOfVertices() << std::endl;
        //std::cout << "Number of edges: " << neuronMorphology->getNumberOfEdges() << std::endl;
        //std::cout << "Number of points: " << neuronMorphology->getNumberOfPoints() << std::endl;
		std::vector< std::vector< double > > allPoints;
		
		if(checkPoint->getParameters().axonFlag)
		{
            
            //Sam's cells with only one vertex at soma location:
            double somaCenter[3];
        	std::vector< Vertex * >::const_iterator vertexIt;
        	for(vertexIt = neuronMorphology->verticesBegin(); vertexIt != neuronMorphology->verticesEnd(); ++vertexIt)
        	{
                if((*vertexIt)->label == Soma)
                {
                    somaCenter[0] = (*vertexIt)->coordinates[0];
                    somaCenter[1] = (*vertexIt)->coordinates[1];
                    somaCenter[2] = (*vertexIt)->coordinates[2];
                    break;
                }
        	}
            //PolyDataPointerType structure = PolyDataPointerType::New();
        	//if(!neuronMorphology->extractLandmark(Soma, structure))
        	//{
        	//	std::cout << "Error! Could not find structure with ID Soma in SpatialGraph!" << std::endl;
        	//	return 0;
        	//}
        	//int subID;
        	//double pCoords[3], * weights;
        	//weights = new double[structure->GetCell(0)->GetNumberOfPoints()];
        	//structure->GetCell(0)->GetParametricCenter(pCoords);
        	//structure->GetCell(0)->EvaluateLocation(subID, pCoords, somaCenter, weights);
        	//delete [] weights;
        	
        	AmiraSpatialGraph * pathLengthDistanceSG = new AmiraSpatialGraph;
			AmiraSpatialGraph * delayRatioSG = new AmiraSpatialGraph;
			pathLengthDistanceSG->mergeSpatialGraph(neuronMorphology);
			delayRatioSG->mergeSpatialGraph(neuronMorphology);
        	
        	double meanPointDistance = 0.0;
			unsigned int nrOfPoints = 0;
			unsigned int nrOfEdges = 1;
            std::vector< double > tipPathLengths;
            std::vector< double > pointPathLengths;
            std::vector< double > pointEuclideanDistances;
			PointsPointerType endNodes = PointsPointerType::New();
			endNodes->SetDataTypeToFloat();
			endNodes->Allocate(1);
        	// start iteration through all edges
        	std::vector< Edge * >::const_iterator edgeIt;
        	std::vector< Edge * >::const_iterator pathLengthEdgeIt;
        	std::vector< Edge * >::const_iterator delayEdgeIt;
			
			std::vector< int > intermediateEdges;
        	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
			{
				Edge * currentEdge = *edgeIt;
				if(currentEdge->label != Axon)
				{
					continue;
				}
				int fatherID = currentEdge->fatherID;
				if(fatherID != -1)
				{
					intermediateEdges.push_back(fatherID);
				}
			}
        	for(int i = 0; i < neuronMorphology->edgesPointer()->size(); ++i)
			{
				Edge * currentEdge = neuronMorphology->edgesPointer()->at(i);
				if(currentEdge->label != Axon)
				{
					continue;
				}
				if(std::find(intermediateEdges.begin(), intermediateEdges.end(), i) == intermediateEdges.end())
				{
					double tipDistance = currentEdge->segmentLength();
					int fatherID = currentEdge->fatherID;
					Edge * tmpEdge;
					while(fatherID != -1)
					{
						tmpEdge = neuronMorphology->edgesPointer()->at(fatherID);
						tipDistance += tmpEdge->segmentLength();
						fatherID = tmpEdge->fatherID;
					}
					tipPathLengths.push_back(tipDistance);
					double * tmpPt = currentEdge->edgePointCoordinates.back();
					endNodes->InsertNextPoint(tmpPt);
				}
			}
			
        	for(edgeIt = neuronMorphology->edgesBegin(), pathLengthEdgeIt = pathLengthDistanceSG->edgesBegin(), delayEdgeIt = delayRatioSG->edgesBegin();
				edgeIt != neuronMorphology->edgesEnd() && pathLengthEdgeIt != pathLengthDistanceSG->edgesEnd() && delayEdgeIt != delayRatioSG->edgesEnd();
				++edgeIt, ++pathLengthEdgeIt, ++delayEdgeIt)
        	{
        		if((*edgeIt)->label != Axon)
				{
					std::list< double >::iterator pathLengthEdgePtIt;
					for(pathLengthEdgePtIt = (*pathLengthEdgeIt)->radiusList.begin(); pathLengthEdgePtIt != (*pathLengthEdgeIt)->radiusList.end(); ++pathLengthEdgePtIt)
					{
						*pathLengthEdgePtIt = 0;
					}
					std::list< double >::iterator delayRatioEdgePtIt;
					for(delayRatioEdgePtIt = (*delayEdgeIt)->radiusList.begin(); delayRatioEdgePtIt != (*delayEdgeIt)->radiusList.end(); ++delayRatioEdgePtIt)
					{
						*delayRatioEdgePtIt = 0;
					}
        			continue;
				}
        		
				std::flush(std::cout << "Processing axon edge " << nrOfEdges << std::endl);
				std::flush(std::cout << "Number of points = " << (*edgeIt)->numEdgePoints << std::endl);
				
				meanPointDistance += (*edgeIt)->segmentLength();
				nrOfPoints += (*edgeIt)->numEdgePoints;
				
				std::list< double > edgePathLengthDistances;
				std::list< double > edgeDelayRatios;
                double tmpLength = 0.0;
        		std::list< double * >::const_iterator ptIt1, ptIt2;
        		ptIt1 = (*edgeIt)->edgePointCoordinates.begin();
        		ptIt2 = ptIt1;
				++ptIt2;
				unsigned int pointCnt = 1;
        		while(ptIt2 != (*edgeIt)->edgePointCoordinates.end())
        		{
					std::flush(std::cout << "\tPoint nr. " << pointCnt << "\r");
                    double * pt1 = *ptIt1;
                    double * pt2 = *ptIt2;
                    tmpLength += sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
					double tmpLength2 = tmpLength;
					int fatherID = (*edgeIt)->fatherID;
// 					std::flush(std::cout << "\tFather ID = " << fatherID << "\r");
					int currID;
					Edge * currentEdge;
					if(fatherID == -1)
					{
						currentEdge = *edgeIt;
						double * tmpPt = currentEdge->edgePointCoordinates.front();
						double pt1SomaDist = sqrt(vtkMath::Distance2BetweenPoints(tmpPt, somaCenter));
						tmpLength2 += pt1SomaDist;
// 						std::cout << "Pt1 Soma dist = " << pt1SomaDist << "\r";
					}
					while(fatherID != -1)
					{
						currID = fatherID;
						currentEdge = neuronMorphology->edgesPointer()->at(currID);
						tmpLength2 += currentEdge->segmentLength();
						fatherID = currentEdge->fatherID;
						if(fatherID == -1)
						{
							double * tmpPt = currentEdge->edgePointCoordinates.front();
							double pt1SomaDist = sqrt(vtkMath::Distance2BetweenPoints(tmpPt, somaCenter));
							tmpLength2 += pt1SomaDist;
// 							std::cout << "Pt1 Soma dist = " << pt1SomaDist << "\r";
						}
					}
        		
					pointPathLengths.push_back(tmpLength2);
					double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(pt2, somaCenter));
					pointEuclideanDistances.push_back(tmpDist);
					edgePathLengthDistances.push_back(tmpLength2);
					edgeDelayRatios.push_back(tmpLength2/tmpDist);
					std::vector< double > tmpPt;
					for(int n = 0; n < 3; ++n)
					{
						tmpPt.push_back(pt2[n]);
					}
					if(edgePathLengthDistances.size() == 1)
					{
						edgePathLengthDistances.push_back(tmpLength2);
						edgeDelayRatios.push_back(tmpLength2/tmpDist);
						for(int n = 0; n < 3; ++n)
						{
							tmpPt.push_back(pt2[n]);
						}
					}
					allPoints.push_back(tmpPt);
					
        			++ptIt1, ++ptIt2;
					++pointCnt;
					
					if(tmpLength2 < tmpDist)
					{
						std::cout << std::endl;
						std::cout << "Path length soma distance smaller than euclidean soma distance!" << std::endl;
						std::cout << "tmpLength2 = " << tmpLength2 << std::endl;
						std::cout << "tmpDist = " << tmpDist << std::endl;
						std::cout << "pt1 @ [" << pt1[0] << "," << pt1[1] << "," << pt1[2] << "]" << std::endl;
						std::cout << "pt2 @ [" << pt2[0] << "," << pt2[1] << "," << pt2[2] << "]" << std::endl;
					}
        		}
        		std::cout << std::endl;
				++nrOfEdges;
				
				(*pathLengthEdgeIt)->radiusList = edgePathLengthDistances;
				(*delayEdgeIt)->radiusList = edgeDelayRatios;
    		}
    		
    		meanPointDistance /= nrOfPoints;
			int downSampleFactor = int(meanSynapseDistance/meanPointDistance + 0.5);
			std::cout << "meanPointDistance = " << meanPointDistance << std::endl;
			std::cout << "downSampleFactor = " << downSampleFactor << std::endl;
    		
			std::vector< std::vector< double > > ptAttributeList;
			
			std::string TableName = ofName + "_euclidean_distance_path_length_list_%.1f_mu.csv";
			char * TableNameChar = new char[256];
			sprintf(TableNameChar, TableName.c_str(), meanSynapseDistance);
        	std::ofstream Table;
        	Table.open(TableNameChar);
        	Table << "Euclidean distance (microns)\tPath length distance (microns)\n";
            for(int i = 0; i < pointPathLengths.size(); ++i)
            {
				if(!(i%downSampleFactor))
				{
					Table << pointEuclideanDistances[i] << "\t" << pointPathLengths[i] << std::endl;
					std::vector< double > tmpVec;
					for(int j = 0; j < 3; ++j)
					{
						tmpVec.push_back(allPoints[i][j]);
					}
					tmpVec.push_back(pointPathLengths[i]);
					ptAttributeList.push_back(tmpVec);
				}
            }
            Table.close();
			
			std::string sampledSpheresName = ofName + "_sampledSpheres_%.1f_mu";
			char * sampledSpheresNameChar = new char[256];
			sprintf(sampledSpheresNameChar, sampledSpheresName.c_str(), meanSynapseDistance);
			Reader * sampledSpheresWriter = new Reader(sampledSpheresNameChar, sampledSpheresNameChar);
			sampledSpheresWriter->writeLandmarkAttributeFile(ptAttributeList);
			delete sampledSpheresWriter;
			
			std::string endNodesName = ofName + "_endNodes_%.1f_mu";
			char * endNodesNameChar = new char[256];
			sprintf(endNodesNameChar, endNodesName.c_str(), meanSynapseDistance);
			Reader * endNodesWriter = new Reader(endNodesNameChar, endNodesNameChar);
			endNodesWriter->writeLandmarkFile(endNodes);
			delete endNodesWriter;
			
			std::string pathLengthName = ofName + "_path_length_distances_%.1f_mu";
			char * pathLengthNameChar = new char[256];
			sprintf(pathLengthNameChar, pathLengthName.c_str(), meanSynapseDistance);
			Reader * pathLengthWriter = new Reader(pathLengthNameChar, pathLengthNameChar);
			pathLengthWriter->setSpatialGraph(pathLengthDistanceSG);
			pathLengthWriter->writeSpatialGraphFile();
			
			std::string delayRatioName = ofName + "_delay_ratios_%.1f_mu";
			char * delayRatioNameChar = new char[256];
			sprintf(delayRatioNameChar, delayRatioName.c_str(), meanSynapseDistance);
			Reader * delayRatioWriter = new Reader(delayRatioNameChar, delayRatioNameChar);
			delayRatioWriter->setSpatialGraph(delayRatioSG);
			delayRatioWriter->writeSpatialGraphFile();
			
			std::string Table2Name = ofName + "_endnode_path_length_list_%.1f_mu.csv";
			char * Table2NameChar = new char[256];
			sprintf(Table2NameChar, Table2Name.c_str(), meanSynapseDistance);
        	std::ofstream Table2;
        	Table2.open(Table2NameChar);
        	Table2 << "End node path length distance (microns)\n";
            for(int i = 0; i < tipPathLengths.size(); ++i)
            {
				Table2 << tipPathLengths[i] << std::endl;
            }
            Table2.close();
		}
		else
		{
            std::cout << "Error: Morphology does not have any structure with label Axon or Soma" << std::endl;
		}
		
		delete checkPoint;
		delete hocReader;
	}
	
	else if(argc == 4)
	{
		const char * inputFilename = argv[1];
		double meanSynapseDistance = atof(argv[2]);
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
		
		std::string ofName(ifName, 0, suffix);
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		
		AmiraSpatialGraph * neuronMorphology = hocReader->getSpatialGraph();
		std::vector< std::vector< double > > allPoints;
		
		if(checkPoint->getParameters().axonFlag)
		{
            
            //Sam's cells with only one vertex at soma location:
            double somaCenter[3];
        	std::vector< Vertex * >::const_iterator vertexIt;
        	for(vertexIt = neuronMorphology->verticesBegin(); vertexIt != neuronMorphology->verticesEnd(); ++vertexIt)
        	{
                if((*vertexIt)->label == Soma)
                {
                    somaCenter[0] = (*vertexIt)->coordinates[0];
                    somaCenter[1] = (*vertexIt)->coordinates[1];
                    somaCenter[2] = (*vertexIt)->coordinates[2];
                    break;
                }
        	}
        	
        	double meanPointDistance = 0.0;
			unsigned int nrOfPoints = 0;
			unsigned int nrOfEdges = 1;
            std::vector< double > pointPathLengths;
            std::vector< double > pointEuclideanDistances;
			PointsPointerType sampledSynapseLocations = PointsPointerType::New();
			sampledSynapseLocations->SetDataTypeToFloat();
			sampledSynapseLocations->Allocate(1);
        	// start iteration through all edges
        	std::vector< Edge * >::const_iterator edgeIt;
        	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
        	{
        		if((*edgeIt)->label != Axon)
				{
        			continue;
				}
        		
				std::flush(std::cout << "Processing axon edge " << nrOfEdges << std::endl);
				std::flush(std::cout << "Number of points = " << (*edgeIt)->numEdgePoints << std::endl);
				
				meanPointDistance += (*edgeIt)->segmentLength();
				nrOfPoints += (*edgeIt)->numEdgePoints;
				
				std::list< double > edgePathLengthDistances;
				std::list< double > edgeDelayRatios;
                double tmpLength = 0.0;
        		std::list< double * >::const_iterator ptIt1, ptIt2;
        		ptIt1 = (*edgeIt)->edgePointCoordinates.begin();
        		ptIt2 = ptIt1;
				++ptIt2;
				unsigned int pointCnt = 1;
        		while(ptIt2 != (*edgeIt)->edgePointCoordinates.end())
        		{
					std::flush(std::cout << "\tPoint nr. " << pointCnt << "\r");
                    double * pt1 = *ptIt1;
                    double * pt2 = *ptIt2;
                    tmpLength += sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
					double tmpLength2 = tmpLength;
					int fatherID = (*edgeIt)->fatherID;
// 					std::flush(std::cout << "\tFather ID = " << fatherID << "\r");
					int currID;
					Edge * currentEdge;
					if(fatherID == -1)
					{
						currentEdge = *edgeIt;
						double * tmpPt = currentEdge->edgePointCoordinates.front();
						double pt1SomaDist = sqrt(vtkMath::Distance2BetweenPoints(tmpPt, somaCenter));
						tmpLength2 += pt1SomaDist;
// 						std::cout << "Pt1 Soma dist = " << pt1SomaDist << "\r";
					}
					while(fatherID != -1)
					{
						currID = fatherID;
						currentEdge = neuronMorphology->edgesPointer()->at(currID);
						tmpLength2 += currentEdge->segmentLength();
						fatherID = currentEdge->fatherID;
						if(fatherID == -1)
						{
							double * tmpPt = currentEdge->edgePointCoordinates.front();
							double pt1SomaDist = sqrt(vtkMath::Distance2BetweenPoints(tmpPt, somaCenter));
							tmpLength2 += pt1SomaDist;
// 							std::cout << "Pt1 Soma dist = " << pt1SomaDist << "\r";
						}
					}
        		
					pointPathLengths.push_back(tmpLength2);
					double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(pt2, somaCenter));
					pointEuclideanDistances.push_back(tmpDist);
					edgePathLengthDistances.push_back(tmpLength2);
					edgeDelayRatios.push_back(tmpLength2/tmpDist);
					std::vector< double > tmpPt;
					for(int n = 0; n < 3; ++n)
					{
						tmpPt.push_back(pt2[n]);
					}
					if(edgePathLengthDistances.size() == 1)
					{
						edgePathLengthDistances.push_back(tmpLength2);
						edgeDelayRatios.push_back(tmpLength2/tmpDist);
						for(int n = 0; n < 3; ++n)
						{
							tmpPt.push_back(pt2[n]);
						}
					}
					allPoints.push_back(tmpPt);
					
        			++ptIt1, ++ptIt2;
					++pointCnt;
					
					if(tmpLength2 < tmpDist)
					{
						std::cout << std::endl;
						std::cout << "Path length soma distance smaller than euclidean soma distance!" << std::endl;
						std::cout << "tmpLength2 = " << tmpLength2 << std::endl;
						std::cout << "tmpDist = " << tmpDist << std::endl;
						std::cout << "pt1 @ [" << pt1[0] << "," << pt1[1] << "," << pt1[2] << "]" << std::endl;
						std::cout << "pt2 @ [" << pt2[0] << "," << pt2[1] << "," << pt2[2] << "]" << std::endl;
					}
        		}
        		std::cout << std::endl;
				++nrOfEdges;
    		}
    		
    		meanPointDistance /= nrOfPoints;
			int downSampleFactor = int(meanSynapseDistance/meanPointDistance + 0.5);
			std::cout << "meanPointDistance = " << meanPointDistance << std::endl;
			std::cout << "downSampleFactor = " << downSampleFactor << std::endl;
    		
			const gsl_rng_type * T;
			gsl_rng * r;
			gsl_rng_env_setup();
			T = gsl_rng_default;
			r = gsl_rng_alloc(T);
			
			std::vector< std::vector< double > > ptAttributeList;
			std::string TableName = ofName + "_euclidean_distance_path_length_list_RAsynapses.csv";
        	std::ofstream Table;
        	Table.open(TableName.c_str());
        	Table << "Euclidean distance (microns)\tPath length distance (microns)\n";
            for(int i = 0; i < pointPathLengths.size(); ++i)
            {
				if(!(i%downSampleFactor) && gsl_rng_uniform(r) < RASynapseProbability(pointEuclideanDistances[i]))
				{
					Table << pointEuclideanDistances[i] << "\t" << pointPathLengths[i] << std::endl;
					std::vector< double > tmpVec;
					for(int j = 0; j < 3; ++j)
					{
						tmpVec.push_back(allPoints[i][j]);
					}
					tmpVec.push_back(pointPathLengths[i]);
					ptAttributeList.push_back(tmpVec);
				}
            }
            Table.close();
			
			std::string sampledSpheresName = ofName + "_sampledSpheres_RAsynapses";
			char * sampledSpheresNameChar = new char[256];
			sprintf(sampledSpheresNameChar, sampledSpheresName.c_str(), meanSynapseDistance);
			Reader * sampledSpheresWriter = new Reader(sampledSpheresNameChar, sampledSpheresNameChar);
			sampledSpheresWriter->writeLandmarkAttributeFile(ptAttributeList);
			delete sampledSpheresWriter;
		}
		else
		{
            std::cout << "Error: Morphology does not have any structure with label Axon or Soma" << std::endl;
		}
		
		delete checkPoint;
		delete hocReader;
	}
	
	if(argc == 5)
	{
		const char * inputFilename = argv[1];
		double sigmoidA = atof(argv[2]);
		double sigmoidB = atof(argv[3]);
		double sigmoidC = atof(argv[4]);
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
		
		std::string ofName(ifName, 0, suffix);
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		
		AmiraSpatialGraph * neuronMorphology = hocReader->getSpatialGraph();
        //std::cout << "Number of vertices: " << neuronMorphology->getNumberOfVertices() << std::endl;
        //std::cout << "Number of edges: " << neuronMorphology->getNumberOfEdges() << std::endl;
        //std::cout << "Number of points: " << neuronMorphology->getNumberOfPoints() << std::endl;
		std::vector< std::vector< double > > allPoints;
		
		if(checkPoint->getParameters().axonFlag)
		{
            
            //Sam's cells with only one vertex at soma location:
            double somaCenter[3];
        	std::vector< Vertex * >::const_iterator vertexIt;
        	for(vertexIt = neuronMorphology->verticesBegin(); vertexIt != neuronMorphology->verticesEnd(); ++vertexIt)
        	{
                if((*vertexIt)->label == Soma)
                {
                    somaCenter[0] = (*vertexIt)->coordinates[0];
                    somaCenter[1] = (*vertexIt)->coordinates[1];
                    somaCenter[2] = (*vertexIt)->coordinates[2];
                    break;
                }
        	}

        	AmiraSpatialGraph * somaDistanceSG = new AmiraSpatialGraph;
			somaDistanceSG->mergeSpatialGraph(neuronMorphology);

        	// start iteration through all edges
			int nrOfEdges = 0;
        	std::vector< Edge * >::const_iterator edgeIt;
        	std::vector< Edge * >::const_iterator pathLengthEdgeIt;
        	for(edgeIt = neuronMorphology->edgesBegin(), pathLengthEdgeIt = somaDistanceSG->edgesBegin();
				edgeIt != neuronMorphology->edgesEnd() && pathLengthEdgeIt != somaDistanceSG->edgesEnd();
				++edgeIt, ++pathLengthEdgeIt)
        	{
        		if((*edgeIt)->label != Axon)
				{
					std::list< double >::iterator pathLengthEdgePtIt;
					for(pathLengthEdgePtIt = (*pathLengthEdgeIt)->radiusList.begin(); pathLengthEdgePtIt != (*pathLengthEdgeIt)->radiusList.end(); ++pathLengthEdgePtIt)
					{
						*pathLengthEdgePtIt = 0;
					}
				}
        		nrOfEdges += 1;
				std::flush(std::cout << "Processing axon edge " << nrOfEdges << std::endl);
				std::flush(std::cout << "Number of points = " << (*edgeIt)->numEdgePoints << std::endl);
				
				std::list< double > edgeSomaDistances;
        		std::list< double * >::const_iterator ptIt = (*edgeIt)->edgePointCoordinates.begin();
        		while(ptIt != (*edgeIt)->edgePointCoordinates.end())
        		{
                    double * pt = *ptIt;
                    double somaDist = sqrt(vtkMath::Distance2BetweenPoints(pt, somaCenter));
					double proportion = sigmoidA/(1.0 + sigmoidB*exp(sigmoidC*somaDist));
					edgeSomaDistances.push_back(proportion);
					++ptIt;
        		}
        		std::cout << std::endl;
				
				(*pathLengthEdgeIt)->radiusList = edgeSomaDistances;
    		}
			
			std::string somaDistanceName = ofName + "_target_proportion";
			Reader * somaDistanceWriter = new Reader(somaDistanceName.c_str(), somaDistanceName.c_str());
			somaDistanceWriter->setSpatialGraph(somaDistanceSG);
			somaDistanceWriter->writeSpatialGraphFile();
			delete somaDistanceWriter;
		}
		else
		{
            std::cout << "Error: Morphology does not have any structure with label Axon or Soma" << std::endl;
		}
		
		delete checkPoint;
		delete hocReader;
	}
	
 	else
 	{
 		std::cout << "Usage: AxonBranchLengthDistance [Input filename] or AxonBranchLengthDistance [Input filename] [avg. inter-synapse distance]" << std::endl;
 	}
	return 0;
}

double RASynapseProbability(double somaDistance)
{
	double A = 0.157;
	double B = 30.72;
	double C = -0.04577;
	return A/(1 + B*exp(C*somaDistance));
}


