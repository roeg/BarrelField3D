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

// #define DEBUG

std::map< unsigned int, std::list< unsigned int > > getConnectedComponents(AmiraSpatialGraph * sg);
void getCenterOfMass(AmiraSpatialGraph * sg, std::list< unsigned int > connectedComponent, double COM[3]);
void getCenterOfMass(Edge * edge, double COM[3]);

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
		else
		{
			std::cout << "Error! Can only analyze .hoc files!" << std::endl;
			return 0;
		}
		
		std::string ofName(ifName, 0, suffix);
		InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
		checkPoint->checkNeuronMorphology();
		
		AmiraSpatialGraph * neuronMorphology = hocReader->getSpatialGraph();
        //std::cout << "Number of vertices: " << neuronMorphology->getNumberOfVertices() << std::endl;
        //std::cout << "Number of edges: " << neuronMorphology->getNumberOfEdges() << std::endl;
        //std::cout << "Number of points: " << neuronMorphology->getNumberOfPoints() << std::endl;
		
		if(checkPoint->getParameters().axonFlag && checkPoint->getParameters().somaFlag)
		{
			// per-branch analysis
			std::vector< double > BPSomaDistancesEuclidean;
			std::vector< double > BPSomaDistancesPathLength;
			std::vector< double > EndNodeSomaDistancesEuclidean;
			std::vector< double > EndNodeSomaDistancesPathLength;
			std::vector< double > correspondingBranchLengthBP;
			std::vector< double > correspondingBranchLengthEndNode;
			std::vector< double > branchSomaDistanceEuclidean;
			std::vector< double > branchSomaDistancePathLength;
			std::vector< double > branchLength;
			std::vector< double > branchNodeEuclideanDistances;
			
            double somaCenter[3];
            PolyDataPointerType structure = PolyDataPointerType::New();
        	if(!neuronMorphology->extractLandmark(Soma, structure))
        	{
        		std::cout << "Error! Could not find structure with ID Soma in SpatialGraph!" << std::endl;
        		return 0;
        	}
        	int subID;
        	double pCoords[3], * weights;
        	weights = new double[structure->GetCell(0)->GetNumberOfPoints()];
        	structure->GetCell(0)->GetParametricCenter(pCoords);
        	structure->GetCell(0)->EvaluateLocation(subID, pCoords, somaCenter, weights);
        	delete [] weights;
			
#ifdef DEBUG
			PointsPointerType AxonBranchPoints = PointsPointerType::New();
			AxonBranchPoints->SetDataTypeToFloat();
#endif
			
// 			std::vector< double * > branchPoints;
			std::vector< int > visitedEdgeIDs;
        	std::vector< Edge * >::const_iterator edgeIt;
        	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
        	{
        		if((*edgeIt)->label != Axon)
        			continue;
				
				int fatherID = (*edgeIt)->fatherID;
				if(fatherID > -1 && std::find(visitedEdgeIDs.begin(), visitedEdgeIDs.end(), fatherID) == visitedEdgeIDs.end())
				{
// 					double * tmpBPCoordinates = neuronMorphology->edgesPointer()->at(fatherID)->edgePointCoordinates.back();
// 					branchPoints.push_back(tmpBPCoordinates);
					visitedEdgeIDs.push_back(fatherID);
#ifdef DEBUG
					AxonBranchPoints->InsertNextPoint(tmpBPCoordinates);
#endif
				}
			}
			// calculate path length and euclidean distance of BPs/end nodes
        	for(int branchID = 0; branchID < neuronMorphology->edgesPointer()->size(); ++branchID)
        	{
				// branch points
				if(std::find(visitedEdgeIDs.begin(), visitedEdgeIDs.end(), branchID) != visitedEdgeIDs.end())
				{
					Edge * branch = neuronMorphology->edgesPointer()->at(branchID);
					double * tmpBP = branch->edgePointCoordinates.back();
					double * tmpBP2 = branch->edgePointCoordinates.front();
					double euclideanDistance = sqrt(vtkMath::Distance2BetweenPoints(tmpBP, somaCenter));
					double branchNodeDistance = sqrt(vtkMath::Distance2BetweenPoints(tmpBP, tmpBP2));
					double pathLengthDistance = branch->segmentLength();
					correspondingBranchLengthBP.push_back(pathLengthDistance);
					branchLength.push_back(pathLengthDistance);
					while(branch->fatherID > -1)
					{
						int tmpFatherID = branch->fatherID;
						branch = neuronMorphology->edgesPointer()->at(tmpFatherID);
						pathLengthDistance += branch->segmentLength();
					}
					BPSomaDistancesEuclidean.push_back(euclideanDistance);
					BPSomaDistancesPathLength.push_back(pathLengthDistance);
					branchSomaDistanceEuclidean.push_back(euclideanDistance);
					branchSomaDistancePathLength.push_back(pathLengthDistance);
					branchNodeEuclideanDistances.push_back(branchNodeDistance);
				}
				// end nodes
				else
				{
					Edge * branch = neuronMorphology->edgesPointer()->at(branchID);
					double * tmpBP = branch->edgePointCoordinates.back();
					double * tmpBP2 = branch->edgePointCoordinates.front();
					double euclideanDistance = sqrt(vtkMath::Distance2BetweenPoints(tmpBP, somaCenter));
					double branchNodeDistance = sqrt(vtkMath::Distance2BetweenPoints(tmpBP, tmpBP2));
					double pathLengthDistance = branch->segmentLength();
					correspondingBranchLengthEndNode.push_back(pathLengthDistance);
					branchLength.push_back(pathLengthDistance);
					while(branch->fatherID > -1)
					{
						int tmpFatherID = branch->fatherID;
						branch = neuronMorphology->edgesPointer()->at(tmpFatherID);
						pathLengthDistance += branch->segmentLength();
					}
					EndNodeSomaDistancesEuclidean.push_back(euclideanDistance);
					EndNodeSomaDistancesPathLength.push_back(pathLengthDistance);
					branchSomaDistanceEuclidean.push_back(euclideanDistance);
					branchSomaDistancePathLength.push_back(pathLengthDistance);
					branchNodeEuclideanDistances.push_back(branchNodeDistance);
				}
        	}
			
// 			std::string BPOutName(ofName);
// 			BPOutName += "_branch_nodes.csv";
// 			std::ofstream TableBPs;
// 			TableBPs.open(BPOutName.c_str());
// 			TableBPs << "Branch length (microns)\tBP-soma euclidean distance\tBP-soma path length distance\n";
// 			for(int i = 0; i < correspondingBranchLengthBP.size(); ++i)
// 			{
// 				TableBPs << correspondingBranchLengthBP[i] << "\t" << BPSomaDistancesEuclidean[i] << "\t" << BPSomaDistancesPathLength[i] << std::endl;
// 			}
// 			TableBPs.close();
// 			
// 			std::string EndNodeOutName(ofName);
// 			EndNodeOutName += "_end_nodes.csv";
// 			std::ofstream TableEndNodes;
// 			TableEndNodes.open(EndNodeOutName.c_str());
// 			TableEndNodes << "Branch length (microns)\tEnd node-soma euclidean distance\tEnd node-soma path length distance\n";
// 			for(int i = 0; i < correspondingBranchLengthEndNode.size(); ++i)
// 			{
// 				TableEndNodes << correspondingBranchLengthEndNode[i] << "\t" << EndNodeSomaDistancesEuclidean[i] << "\t" << EndNodeSomaDistancesPathLength[i] << std::endl;
// 			}
// 			TableEndNodes.close();
			
			std::string TotalOutName(ofName);
			TotalOutName += "_branch_nodes.csv";
			std::ofstream Table;
			Table.open(TotalOutName.c_str());
			Table << "Branch length (microns)\tbranch node euclidean distance\tBP-soma euclidean distance\tBP-soma path length distance\tTotal Tortoisity\tSegment Tortoisity\n";
			for(int i = 0; i < branchLength.size(); ++i)
			{
				double totalTortoisity = branchSomaDistancePathLength[i]/branchSomaDistanceEuclidean[i];
				double segmentTortoisity = branchLength[i]/branchNodeEuclideanDistances[i];
				Table << branchLength[i] << "\t" << branchNodeEuclideanDistances[i] << "\t" << branchSomaDistanceEuclidean[i] << "\t" << branchSomaDistancePathLength[i] << "\t" << totalTortoisity << "\t" << segmentTortoisity << std::endl;
			}
			Table.close();
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
 		std::cout << "Usage: AxonBoxAnalysis [Input filename] [box size X] [box size Y] [box size Z]" << std::endl;
 	}
	return 0;
}

std::map< unsigned int, std::list< unsigned int > > getConnectedComponents(AmiraSpatialGraph* sg)
{
	// construct dual graph from SpatialGraph 
	// in dual graph, SpatialGraph edges are vertices,
	// and edges describe which SpatialGraph edges are connected
	MutableUndirectedGraphPointerType dualGraph = MutableUndirectedGraphPointerType::New();
	std::vector< vtkIdType > dualVertexIDs;
	for(int i = 0; i < sg->edgesPointer()->size(); ++i)
	{
		dualVertexIDs.push_back(dualGraph->AddVertex());
	}
	for(int i = 0; i < sg->edgesPointer()->size(); ++i)
	{
		vtkIdType ID1 = dualVertexIDs[i];
		double * pt1 = sg->edgesPointer()->at(i)->edgePointCoordinates.front();
		double * pt2 = sg->edgesPointer()->at(i)->edgePointCoordinates.back();
		for(int j = i+1; j < sg->edgesPointer()->size(); ++j)
		{
			vtkIdType ID2 = dualVertexIDs[j];
			double * pt3 = sg->edgesPointer()->at(j)->edgePointCoordinates.front();
			double * pt4 = sg->edgesPointer()->at(j)->edgePointCoordinates.back();
			double dist13 = sqrt(vtkMath::Distance2BetweenPoints(pt1, pt3));
			double dist14 = sqrt(vtkMath::Distance2BetweenPoints(pt1, pt4));
			double dist23 = sqrt(vtkMath::Distance2BetweenPoints(pt2, pt3));
			double dist24 = sqrt(vtkMath::Distance2BetweenPoints(pt2, pt4));
			if(dist13 < 1e-6 || dist14 < 1e-6 || dist23 < 1e-6 || dist24 < 1e-6)
			{
				dualGraph->AddEdge(ID1, ID2);
			}
		}
	}
	
	vtkBoostConnectedComponents * connectedComponentFilter = vtkBoostConnectedComponents::New();
	connectedComponentFilter->SetInput(dualGraph);
	connectedComponentFilter->Update();
	GraphPointerType outGraph = connectedComponentFilter->GetOutput();
	IntArrayPointerType components = vtkIntArray::SafeDownCast(outGraph->GetVertexData()->GetArray("component"));
	
	std::map< unsigned int, std::list< unsigned int > > edgesPerComponent;
	for(int i = 0; i < sg->edgesPointer()->size(); ++i)
	{
		int component = components->GetValue(dualVertexIDs[i]);
		if(edgesPerComponent.find(component) == edgesPerComponent.end())
		{
			std::list< unsigned int > componentEdges;
			componentEdges.push_back(i);
			edgesPerComponent[component] = componentEdges;
		}
		else
		{
			edgesPerComponent[component].push_back(i);
		}
	}
	
	return edgesPerComponent;
};

void getCenterOfMass(AmiraSpatialGraph* sg, std::list< unsigned int > connectedComponent, double COM[3])
{
	COM[0] = 0;
	COM[1] = 0;
	COM[2] = 0;
	
	unsigned int totalNrEdges = 0;
	std::list< unsigned int >::const_iterator connectedComponentIt;
	for(connectedComponentIt = connectedComponent.begin(); connectedComponentIt != connectedComponent.end(); ++connectedComponentIt)
	{
		Edge * edge = sg->edgesPointer()->at(*connectedComponentIt);
		std::list< double * >::const_iterator edgePtIt;
		for(edgePtIt = edge->edgePointCoordinates.begin(); edgePtIt != edge->edgePointCoordinates.end(); ++edgePtIt)
		{
			double * pt = *edgePtIt;
			COM[0] += pt[0];
			COM[1] += pt[1];
			COM[2] += pt[2];
		}
		totalNrEdges += edge->numEdgePoints;
	}
	
	COM[0] /= totalNrEdges;
	COM[1] /= totalNrEdges;
	COM[2] /= totalNrEdges;
};

void getCenterOfMass(Edge* edge, double COM[3])
{
	COM[0] = 0;
	COM[1] = 0;
	COM[2] = 0;
	
	std::list< double * >::const_iterator edgePtIt;
	for(edgePtIt = edge->edgePointCoordinates.begin(); edgePtIt != edge->edgePointCoordinates.end(); ++edgePtIt)
	{
		double * pt = *edgePtIt;
		COM[0] += pt[0];
		COM[1] += pt[1];
		COM[2] += pt[2];
	}
	
	COM[0] /= edge->numEdgePoints;
	COM[1] /= edge->numEdgePoints;
	COM[2] /= edge->numEdgePoints;
};














