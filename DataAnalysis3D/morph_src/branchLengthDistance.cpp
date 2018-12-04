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
    		}
    		
            ofName += "_length_distance_list.csv";
        	std::ofstream Table;
        	Table.open(ofName.c_str());
        	Table << "Segment length (microns)\tSegment midpoint-soma distance (microns)\n";
            for(int i = 0; i < segmentLengths.size(); ++i)
            {
                Table << segmentLengths[i] << "\t" << segmentMidpointDistances[i] << std::endl;
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
 		std::cout << "Usage: AxonBranchLengthDistance [Input filename]" << std::endl;
 	}
	return 0;
}
