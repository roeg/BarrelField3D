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
#include "../../common/profile.h"
#include <gsl/gsl_rng.h>

double RASynapseProbability(double somaDistance);
void rotateVector2D(double vector[2], double rotatedVector[2], double angle);
double getMaximumPathlength(AmiraSpatialGraph * sg, int label);
Profile * getPathLengthHistogram(int nrOfSynapses, double mean, double std, double binSize, double maxPathlength);
double getRootPathlengthDistance(AmiraSpatialGraph * sg, int edgeID);

int main( int argc , char * argv[])
{
	if(argc == 3)
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
			double HVCAngle = 28.5;
			std::vector< std::vector< double > > pointLocationsSomaCenter;
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
					double pt2SomaCentered[2], pt2HVCAligned[2];
					pt2SomaCentered[0] = pt2[0] - somaCenter[0];
					pt2SomaCentered[1] = pt2[1] - somaCenter[1];
					rotateVector2D(pt2SomaCentered, pt2HVCAligned, HVCAngle);
					std::vector< double > vec2SomaCentered;
					vec2SomaCentered.push_back(pt2HVCAligned[0]);
					vec2SomaCentered.push_back(pt2HVCAligned[1]);
					pointLocationsSomaCenter.push_back(vec2SomaCentered);
					
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
			std::string TableName = ofName + "_euclidean_distance_path_length_list_RAsynapses_HVCRotated.csv";
        	std::ofstream Table;
        	Table.open(TableName.c_str());
        	Table << "Euclidean distance (microns)\tPath length distance (microns)\tLocation (x)\tLocation (y)\n";
            for(int i = 0; i < pointPathLengths.size(); ++i)
            {
				if(!(i%downSampleFactor) && gsl_rng_uniform(r) < RASynapseProbability(pointEuclideanDistances[i]))
				{
					Table << pointEuclideanDistances[i] << "\t" << pointPathLengths[i] << "\t" << pointLocationsSomaCenter[i][0] << "\t"  << pointLocationsSomaCenter[i][1] << std::endl;
					std::vector< double > tmpVec;
					// HVC reference frame; 3D location on axon
// 					for(int j = 0; j < 3; ++j)
// 					{
// 						tmpVec.push_back(allPoints[i][j]);
// 					}
					// HVC rotated reference frame; projected 2D location
					tmpVec.push_back(pointLocationsSomaCenter[i][0]);
					tmpVec.push_back(pointLocationsSomaCenter[i][1]);
					tmpVec.push_back(0);
					tmpVec.push_back(pointPathLengths[i]);
					ptAttributeList.push_back(tmpVec);
				}
            }
            Table.close();
			
			std::string sampledSpheresName = ofName + "_sampledSpheres_RAsynapses_HVCRotated";
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
	
	if(argc == 6)
	{
		// map arbitrary mean/STD log-normal delay map onto axonal arbors
		// collapse into radial distribution and compare with 2p imaging data
		
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		const int nrOfSynapses = atoi(argv[3]);
		const double mean_ = atof(argv[4]);
		const double std_ = atof(argv[5]);
		// convert ms delay to microns path length
		const double mean = mean_*187.0;
		const double std = std_*187.0;
		
		// load neuron (hoc file!)
		std::string ifName(inputFilename);
		size_t suffix;
		
		Reader * hocReader = new Reader(inputFilename, inputFilename);
		if(ifName.find(".hoc") != std::string::npos)
		{
			hocReader->readHocFile();
			suffix = ifName.find(".hoc");
		}
		else
		{
			std::cout << "Error! Can only analyze .hoc files!" << std::endl;
			return 0;
		}
		
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		AmiraSpatialGraph * neuronMorphology = hocReader->getSpatialGraph();
		
		// determine maximum path length
		const double maxPathlength = getMaximumPathlength(neuronMorphology, Axon);
		
		// sample N path length distances from log-normal distribution; 
		// replace if > max. path length
		// compute path length histogram in 50 or 100 micron bins (imaging radial extent ~400 microns)
		const double binSize = 50.0;
		Profile * pathLengthHistogram = getPathLengthHistogram(nrOfSynapses, mean, std, binSize, maxPathlength);
		
		// data structure: map: int (path length histogram bin) -> std::vector of points
		// along axon (std::vector< double>) in each bin
		std::map< int, std::vector< std::vector< double > > > binnedPoints;
		for(int i = 0; i < pathLengthHistogram->getProfile()->size(); ++i)
		{
			std::vector< std::vector< double > > tmpVec;
			binnedPoints[i] = tmpVec;
		}
		
		// iterate over axon edges
		for(int i = 0; i < neuronMorphology->edgesPointer()->size(); ++i)
		{
			Edge * currentEdge = neuronMorphology->edgesPointer()->at(i);
			if(currentEdge->label != Axon)
			{
				continue;
			}
			
			const double rootPathlength = getRootPathlengthDistance(neuronMorphology, i);
			std::list< double * >::const_iterator edgePtListIt1 = currentEdge->edgePointCoordinates.begin();
			std::list< double * >::const_iterator edgePtListIt2 = currentEdge->edgePointCoordinates.begin();
			++edgePtListIt2;
			double tmpPathlength = rootPathlength;
			for(; edgePtListIt2 != currentEdge->edgePointCoordinates.end(); ++edgePtListIt1, ++edgePtListIt2)
			{
				// if successive points are > 1 micron apart insert points along
				// connecting line at 0.5 micron intervals
				// store each point location in path length histogram (same resolution as above)
				double * pt1 = *edgePtListIt1;
				double * pt2 = *edgePtListIt2;
				const double tmpDist = sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
				if(tmpDist <= 1.0)
				{
					std::vector< double > vec1;
					std::vector< double > vec2;
					for(int j = 0; j < 3; ++j)
					{
						vec1.push_back(pt1[j]);
						vec2.push_back(pt2[j]);
					}
					int bin1 = int(tmpPathlength/binSize);
					int bin2 = int((tmpPathlength + tmpDist)/binSize);
					binnedPoints[bin1].push_back(vec1);
					binnedPoints[bin2].push_back(vec2);
				}
				else
				{
					// add first point
					std::vector< double > vec1;
					for(int j = 0; j < 3; ++j)
					{
						vec1.push_back(pt1[j]);
					}
					int bin1 = int(tmpPathlength/binSize);
					binnedPoints[bin1].push_back(vec1);
					// do the interpolation here
					double diff[3];
					vtkMath::Subtract(pt2, pt1, diff);
					vtkMath::Normalize(diff);
					const double stepSize = 0.5;
					double tmpStepLength = stepSize;
					while(tmpStepLength < tmpDist)
					{
						double tmpVec[3], tmpPt[3];
						for(int j = 0; j < 3; ++j)
						{
							tmpVec[j] = tmpStepLength * diff[j];
						}
						vtkMath::Add(pt1, tmpVec, tmpPt);
						
						std::vector< double > vecN;
						for(int j = 0; j < 3; ++j)
						{
							vecN.push_back(tmpPt[j]);
						}
						int binN = int((tmpPathlength + tmpStepLength)/binSize);
						binnedPoints[binN].push_back(vecN);
						
						tmpStepLength += stepSize;
					}
					// add last point
					std::vector< double > vec2;
					for(int j = 0; j < 3; ++j)
					{
						vec2.push_back(pt2[j]);
					}
					int bin2 = int((tmpPathlength + tmpDist)/binSize);
					binnedPoints[bin2].push_back(vec2);
				}
				tmpPathlength += tmpDist;
			}
		}
		
		// randomly sample number of these points in each bin according to above histogram
		const gsl_rng_type * T;
		gsl_rng * rng;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		rng = gsl_rng_alloc(T);
		
		PointsPointerType sampledPoints = PointsPointerType::New();
		sampledPoints->SetDataTypeToFloat();
		sampledPoints->Allocate(1);
		int missingPoints = 0;
		for(int i = 0; i < pathLengthHistogram->getProfile()->size(); ++i)
		{
			int nrSamples = int(pathLengthHistogram->getProfile()->at(i));
			int nrPoints = binnedPoints[i].size();
			std::cout << "Bin " << i << " - nr samples: " << nrSamples << " - nr points: " << nrPoints << std::endl;
			int visualizedPoints;
			if(nrPoints < nrSamples)
			{
				missingPoints += nrSamples - nrPoints;
				visualizedPoints = nrPoints;
			}
			else
			{
				visualizedPoints = nrSamples;
			}
			
			if(visualizedPoints)
			{
				int * shuffledIndices = new int[nrPoints];
				for(int j = 0; j < nrPoints; ++j)
				{
					shuffledIndices[j] = j;
				}
				gsl_ran_shuffle(rng, shuffledIndices, nrPoints, sizeof(int));
				for(int j = 0; j < visualizedPoints; ++j)
				{
					int pointIndex = shuffledIndices[j];
					std::vector< double > tmpPt = binnedPoints[i].at(pointIndex);
					double newPt[3];
					newPt[0] = tmpPt[0];
					newPt[1] = tmpPt[1];
					newPt[2] = tmpPt[2];
					sampledPoints->InsertNextPoint(newPt);
				}
				delete [] shuffledIndices;
			}
		}
		std::cout << "Missing points: " << missingPoints << std::endl;
		
		// write out point locations relative to axon morphology for visualization as well as 
		// euclidean soma distance (2D projected for comparison with 2p imaging data)
		std::string ofName(outputFilename);
		std::string pointOutname = ofName + "_sampled_points";
		Reader * landmarkWriter = new Reader(pointOutname.c_str(), pointOutname.c_str());
		landmarkWriter->writeLandmarkFile(sampledPoints);
		
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
		std::cout << "\tSoma location @ [" << somaLoc[0] << "," << somaLoc[1] << "," << somaLoc[2] << "]" << std::endl;
		
		std::string distanceOutname = ofName + "_sampled_points_distances.csv";
		std::ofstream distanceOutfile(distanceOutname.c_str());
		distanceOutfile << "Soma distance (microns)" << std::endl;
		for(int i = 0; i < sampledPoints->GetNumberOfPoints(); ++i)
		{
			double tmpPt[3];
			sampledPoints->GetPoint(i, tmpPt);
// 			std::cout << "Point " << i << " @ [" << tmpPt[0] << "," << tmpPt[1] << "," << tmpPt[2] << "]" << std::endl;
			double tmpDist = sqrt((tmpPt[0] - somaLoc[0]) * (tmpPt[0] - somaLoc[0]) + (tmpPt[1] - somaLoc[1]) * (tmpPt[1] - somaLoc[1]));
// 			std::cout << "tmpDist = " << tmpDist << std::endl;
			distanceOutfile << tmpDist << std::endl;
		}
		distanceOutfile.close();
		
		std::string locationOutname = ofName + "_sampled_points_2Dlocations.csv";
		std::ofstream locationOutfile(locationOutname.c_str());
		locationOutfile << "Point x\tPoint y" << std::endl;
		for(int i = 0; i < sampledPoints->GetNumberOfPoints(); ++i)
		{
			double tmpPt[3];
			sampledPoints->GetPoint(i, tmpPt);
			locationOutfile << tmpPt[0] - somaLoc[0] << "\t" << tmpPt[1] - somaLoc[1] << std::endl;
		}
		locationOutfile.close();
		
		delete landmarkWriter;
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

void rotateVector2D(double vector[2], double rotatedVector[2], double angle)
{
	double tmpVec[3];
	tmpVec[0] = vector[0];
	tmpVec[1] = vector[1];
	tmpVec[2] = 0;
	TransformPointerType rotation = TransformPointerType::New();
	rotation->RotateZ(angle);
	rotation->TransformPoint(tmpVec, tmpVec);
// 	rotation->Print(std::cout);
	rotatedVector[0] = tmpVec[0];
	rotatedVector[1] = tmpVec[1];
}

double getMaximumPathlength(AmiraSpatialGraph * sg, int label)
{
	double maxPathlength = 0.0;
	for(int i = 0; i < sg->edgesPointer()->size(); ++i)
	{
		Edge * currentEdge = sg->edgesPointer()->at(i);
		if(currentEdge->label != label)
		{
			continue;
		}
		
		double tmpPathlength = currentEdge->segmentLength() + getRootPathlengthDistance(sg, i);
		if(tmpPathlength > maxPathlength)
		{
			maxPathlength = tmpPathlength;
		}
	}
	
	return maxPathlength;
}

Profile * getPathLengthHistogram(int nrOfSynapses, double mean, double std, double binSize, double maxPathlength)
{
		const gsl_rng_type * T;
		gsl_rng * rng;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		rng = gsl_rng_alloc(T);
		
		double mu = log(mean / sqrt(1 + std * std / mean / mean));
		double sigma = sqrt(log(1 + std * std / mean / mean));
		
		Profile * histogram = new Profile(binSize);
		
		for(int i = 0; i < nrOfSynapses; ++i)
		{
			double sample = gsl_ran_lognormal(rng, mu, sigma);
			while(sample > maxPathlength)
			{
				sample = gsl_ran_lognormal(rng, mu, sigma);
			}
			
			int bin = int(sample/binSize);
			histogram->incrementBin(bin);
		}
		
		return histogram;
}

double getRootPathlengthDistance(AmiraSpatialGraph * sg, int edgeID)
{
	double rootPathlength = 0.0;
	Edge * currentEdge = sg->edgesPointer()->at(edgeID);
	while(currentEdge->fatherID != -1)
	{
		Edge * fatherEdge = sg->edgesPointer()->at(currentEdge->fatherID);
		if(fatherEdge->label != Soma)
		{
			rootPathlength += fatherEdge->segmentLength();
		}
		currentEdge = fatherEdge;
	}
	
	return rootPathlength;
}




