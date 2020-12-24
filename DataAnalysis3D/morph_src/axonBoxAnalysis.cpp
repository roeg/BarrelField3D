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
	if(argc == 5 || argc == 8)
	{
		const char * inputFilename = argv[1];
		unsigned int boxSize[3];
		boxSize[0] = atoi(argv[2]);
		boxSize[1] = atoi(argv[3]);
		boxSize[2] = atoi(argv[4]);
		double shift[] = {0,0,0};
		if(argc == 8)
		{
			shift[X_COORD] = atof(argv[5]);
			shift[Y_COORD] = atof(argv[6]);
			shift[Z_COORD] = atof(argv[7]);
		}
		
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
		
		if(checkPoint->getParameters().axonFlag && checkPoint->getParameters().somaFlag)
		{
			// per-box analysis
            std::vector< double > boxSomaDistances;
            std::vector< unsigned int > boxBPCounts;
            std::vector< std::vector< double > > boxBranchLengths;
            std::vector< std::vector< unsigned int > > boxBranchBPCounts;
			// per-piece analysis
            std::vector< double > pieceSomaDistances;
            std::vector< double > pieceLengths;
			
			// per-branch analysis
            std::vector< double > fragmentSomaDistances;
            std::vector< double > fragmentLengths;
			std::vector< double > fragmentEuclideanLengths;
            std::vector< unsigned int > fragmentBPs;
			std::vector< double > longTerminalSomaDistances;
			std::vector< double > longTerminalTruncatedLength;
			std::vector< double > longTerminalNonTruncatedLength;
			std::vector< double > isolatedTerminalSomaDistances;
			std::vector< double > isolatedTerminalTruncatedLength;
			std::vector< double > isolatedTerminalNonTruncatedLength;
			
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
			// get 3D location of all branch/end nodes for later comparison
			// with truncated end node branches
			std::vector< double * > branchPoints;
			std::vector< double * > endNodes;
			std::vector< int > terminalBranchIDs;
			std::vector< int > visitedEdgeIDs;
        	std::vector< Edge * >::const_iterator edgeIt;
        	for(edgeIt = neuronMorphology->edgesBegin(); edgeIt != neuronMorphology->edgesEnd(); ++edgeIt)
        	{
        		if((*edgeIt)->label != Axon)
        			continue;
				
				int fatherID = (*edgeIt)->fatherID;
				if(fatherID > -1 && std::find(visitedEdgeIDs.begin(), visitedEdgeIDs.end(), fatherID) == visitedEdgeIDs.end())
				{
					visitedEdgeIDs.push_back(fatherID);
#ifdef DEBUG
					AxonBranchPoints->InsertNextPoint(tmpBPCoordinates);
#endif
				}
			}
			// calculate path length and euclidean distance of BPs/end nodes
        	for(int branchID = 0; branchID < neuronMorphology->edgesPointer()->size(); ++branchID)
        	{
				Edge * branch = neuronMorphology->edgesPointer()->at(branchID);
				double * tmpBP = branch->edgePointCoordinates.back();
				// branch points
				if(std::find(visitedEdgeIDs.begin(), visitedEdgeIDs.end(), branchID) != visitedEdgeIDs.end())
				{
					branchPoints.push_back(tmpBP);
				}
				// end nodes
				else
				{
					endNodes.push_back(tmpBP);
					terminalBranchIDs.push_back(branchID);
				}
        	}
			
#ifdef DEBUG
			std::cout << "Number of BPs: " << branchPoints.size() << std::endl;
			std::string BPOutName = ofName + "_axon_BPs";
			Reader * BPWriter = new Reader(BPOutName.c_str(), BPOutName.c_str());
			BPWriter->writeLandmarkFile(AxonBranchPoints);
			delete BPWriter;
			
			unsigned int SGCount = 0;
#endif
			double totalAxonLength = 0.0;
			unsigned int totalBPCount = 0;
			
			double neuronBounds[6];
			neuronMorphology->getBoundingBox(neuronBounds);
			for(double xMin = neuronBounds[0] - shift[X_COORD]; xMin < neuronBounds[1]; xMin += boxSize[0])
				for(double yMin = neuronBounds[2] - shift[Y_COORD]; yMin < neuronBounds[3]; yMin += boxSize[1])
					for(double zMin = neuronBounds[4] - shift[Z_COORD]; zMin < neuronBounds[5]; zMin += boxSize[2])
					{
						double tmpBox[6];
						tmpBox[0] = xMin;
						tmpBox[1] = xMin + boxSize[0];
						tmpBox[2] = yMin;
						tmpBox[3] = yMin + boxSize[1];
						tmpBox[4] = zMin;
						tmpBox[5] = zMin + boxSize[2];
						
						double boxCenter[3];
						for(int i = 0; i < 3; ++i)
						{
							boxCenter[i] = 0.5*(tmpBox[2*i] + tmpBox[2*i+1]);
						}
						double boxSomaDistance = sqrt(vtkMath::Distance2BetweenPoints(somaCenter, boxCenter));
						
						// 	determine nr. of BPs in box
						unsigned int tmpBPCount = 0;
						double delta[] = {1e-6, 1e-6, 1e-6};
						for(int i = 0; i < branchPoints.size(); ++i)
						{
							if(vtkMath::PointIsWithinBounds(branchPoints[i], tmpBox, delta))
							{
								++tmpBPCount;
							}
						}
						
						// 	truncate SpatialGraph with box dimensions
						// 	if nothing remains: continue with next box
						// 	else: determine length of remaining segments
						AmiraSpatialGraph * tmpSG = neuronMorphology->clipSpatialGraph(tmpBox);
						
// #ifdef DEBUG
// 						std::cout << "clipping box: (" << tmpBox[0] << "," << tmpBox[1] << ") , (" << tmpBox[2] << "," << tmpBox[3] << ") , (" << tmpBox[4] << "," << tmpBox[5] << ")" << std::endl;
// #endif
						if(!tmpSG)
						{
// #ifdef DEBUG
// 							std::cout << "Nothing in box..." << std::endl;
// #endif
							continue;
						}
						
#ifdef DEBUG
						std::cout << "Analyzing remaining segments in box!" << std::endl;
						std::cout << "Remaining segments: " << tmpSG->edgesPointer()->size() << std::endl;
						std::cout << "Remaining vertices: " << tmpSG->verticesPointer()->size() << std::endl;
						char * tmpSuffix = new char[64];
						sprintf(tmpSuffix, "_SG_%d", SGCount);
						++SGCount;
						std::string tmpSGOutName = ofName + std::string(tmpSuffix);
						Reader * tmpSGWriter = new Reader(tmpSGOutName.c_str(), tmpSGOutName.c_str());
						tmpSGWriter->setSpatialGraph(tmpSG);
						tmpSGWriter->writeSpatialGraphFile();
						delete tmpSGWriter;
#endif
						// identify connected components
						std::map< unsigned int, std::list< unsigned int > > connectedComponents = getConnectedComponents(tmpSG);
#ifdef DEBUG
						std::cout << "Nr. connected components: " << connectedComponents.size() << std::endl;
						unsigned int individualComponents = 0;
						unsigned int edgeIDSum = 0;
						double connectedComponentsLength = 0;
						double individualEdgesLength = 0;
						for(edgeIt = tmpSG->edgesBegin(); edgeIt != tmpSG->edgesEnd(); ++edgeIt)
						{
							if((*edgeIt)->label != Axon)
								continue;
							
							individualEdgesLength += (*edgeIt)->segmentLength();
						}
#endif
						std::vector< double > tmpFragmentLengths;
						std::vector< unsigned int > tmpFragmentBPs;
						std::map< unsigned int, std::list< unsigned int > >::const_iterator connectedComponentsIt;
						for(connectedComponentsIt = connectedComponents.begin(); connectedComponentsIt != connectedComponents.end(); ++connectedComponentsIt)
						{
#ifdef DEBUG
							std::cout << "Connected component: " << connectedComponentsIt->first << std::endl;
#endif
							unsigned int somaCorrection = 0;
							double tmpLength = 0;
							double tmpEuclideanDistance = 0;
							std::list< unsigned int >::const_iterator componentEdgesIt;
							for(componentEdgesIt = connectedComponentsIt->second.begin(); componentEdgesIt != connectedComponentsIt->second.end(); ++componentEdgesIt)
							{
								unsigned int edgeID = *componentEdgesIt;
								if(tmpSG->edgesPointer()->at(edgeID)->label == Soma)
								{
									somaCorrection = 1;
									continue;
								}
								tmpLength += tmpSG->edgesPointer()->at(edgeID)->segmentLength();
								totalAxonLength += tmpSG->edgesPointer()->at(edgeID)->segmentLength();
								tmpEuclideanDistance += sqrt(vtkMath::Distance2BetweenPoints(tmpSG->edgesPointer()->at(edgeID)->edgePointCoordinates.front(), tmpSG->edgesPointer()->at(edgeID)->edgePointCoordinates.back()));
								
								double edgeCOM[3];
								getCenterOfMass(tmpSG->edgesPointer()->at(edgeID), edgeCOM);
								double edgeSomaDist = sqrt(vtkMath::Distance2BetweenPoints(edgeCOM, somaCenter));
								pieceSomaDistances.push_back(edgeSomaDist);
								pieceLengths.push_back(tmpSG->edgesPointer()->at(edgeID)->segmentLength());
#ifdef DEBUG
								connectedComponentsLength += tmpSG->edgesPointer()->at(edgeID)->segmentLength();
								++individualComponents;
								edgeIDSum += edgeID;
								std::cout << "\tedge ID: " << edgeID << std::endl;
#endif
							}
							tmpFragmentLengths.push_back(tmpLength);
							unsigned int tmpBPs = (connectedComponentsIt->second.size()-somaCorrection-1)/2;
							tmpFragmentBPs.push_back(tmpBPs);
							totalBPCount += tmpBPs;
							
							double branchCOM[3];
							getCenterOfMass(tmpSG, connectedComponentsIt->second, branchCOM);
							double branchSomaDist = sqrt(vtkMath::Distance2BetweenPoints(branchCOM, somaCenter));
							fragmentSomaDistances.push_back(branchSomaDist);
							fragmentLengths.push_back(tmpLength);
							fragmentEuclideanLengths.push_back(tmpEuclideanDistance);
							fragmentBPs.push_back(tmpBPs);
							
							// terminal nodes
							for(componentEdgesIt = connectedComponentsIt->second.begin(); componentEdgesIt != connectedComponentsIt->second.end(); ++componentEdgesIt)
							{
								unsigned int edgeID = *componentEdgesIt;
								if(tmpSG->edgesPointer()->at(edgeID)->label == Soma)
								{
									continue;
								}
								
// 								double * pt1 = tmpSG->edgesPointer()->at(edgeID)->edgePointCoordinates.front();
// 								double * pt2 = tmpSG->edgesPointer()->at(edgeID)->edgePointCoordinates.back();
								double * tmpPt = tmpSG->edgesPointer()->at(edgeID)->edgePointCoordinates.back();
								bool isTerminalNode = true;
								std::list< unsigned int >::const_iterator otherComponentEdgesIt;
								for(otherComponentEdgesIt = connectedComponentsIt->second.begin(); otherComponentEdgesIt != connectedComponentsIt->second.end(); ++otherComponentEdgesIt)
								{
									if(componentEdgesIt == otherComponentEdgesIt)
									{
										continue;
									}
									unsigned int otherEdgeID = *otherComponentEdgesIt;
									double * otherNode1 = tmpSG->edgesPointer()->at(otherEdgeID)->edgePointCoordinates.front();
									double * otherNode2 = tmpSG->edgesPointer()->at(otherEdgeID)->edgePointCoordinates.back();
									if(sqrt(vtkMath::Distance2BetweenPoints(tmpPt, otherNode1)) < 1e-6 || sqrt(vtkMath::Distance2BetweenPoints(tmpPt, otherNode2)) < 1e-6)
									{
										isTerminalNode = false;
										break;
									}
								}
								
								if(isTerminalNode)
								{
									int terminalBranch = -1;
									for(int i = 0; i < endNodes.size(); ++i)
									{
	// 									if(sqrt(vtkMath::Distance2BetweenPoints(endNodes[i], pt1)) < 1e-6 || sqrt(vtkMath::Distance2BetweenPoints(endNodes[i], pt2)) < 1e-6)
										if(sqrt(vtkMath::Distance2BetweenPoints(endNodes[i], tmpPt)) < 1e-6)
										{
											terminalBranch = terminalBranchIDs[i];
											break;
										}
									}
									if(terminalBranch > -1/* && tmpSG->edgesPointer()->at(edgeID)->segmentLength() > 200.0*/)
									{
										double terminalCOM[3];
										getCenterOfMass(tmpSG->edgesPointer()->at(edgeID), terminalCOM);
										longTerminalSomaDistances.push_back(sqrt(vtkMath::Distance2BetweenPoints(terminalCOM, somaCenter)));
										longTerminalTruncatedLength.push_back(tmpSG->edgesPointer()->at(edgeID)->segmentLength());
										longTerminalNonTruncatedLength.push_back(neuronMorphology->edgesPointer()->at(terminalBranch)->segmentLength());
										// isolated terminal branches
										if(connectedComponentsIt->second.size() == 1)
										{
											isolatedTerminalSomaDistances.push_back(sqrt(vtkMath::Distance2BetweenPoints(terminalCOM, somaCenter)));
											isolatedTerminalTruncatedLength.push_back(tmpSG->edgesPointer()->at(edgeID)->segmentLength());
											isolatedTerminalNonTruncatedLength.push_back(neuronMorphology->edgesPointer()->at(terminalBranch)->segmentLength());
										}
									}
								}
							}
						}
#ifdef DEBUG
						if(fabs(connectedComponentsLength-individualEdgesLength) > 1e-6)
						{
							std::cout << "*********************************************" << std::endl;
							std::cout << "Connected components length = " << connectedComponentsLength << std::endl;
							std::cout << "Individual edges length = " << individualEdgesLength << std::endl;
							std::cout << "*********************************************" << std::endl;
						}
						if(individualComponents != tmpSG->edgesPointer()->size())
						{
							std::cout << "*********************************************" << std::endl;
							std::cout << "WARNING: individualComponents (" << individualComponents << ") != nr. edges (" << tmpSG->edgesPointer()->size() << ")" << std::endl;
							std::cout << "*********************************************" << std::endl;
						}
						unsigned int tmpN = tmpSG->edgesPointer()->size();
						if(edgeIDSum != tmpN*(tmpN-1)/2)
						{
							std::cout << "*********************************************" << std::endl;
							std::cout << "WARNING: not all edges used/duplicate edges used" << std::endl;
							std::cout << "*********************************************" << std::endl;
						}
#endif
						boxSomaDistances.push_back(boxSomaDistance);
						boxBranchBPCounts.push_back(tmpFragmentBPs);
						boxBranchLengths.push_back(tmpFragmentLengths);
						
// 						// version ignoring connected components
// 						std::vector< double > tmpBranchLengths;
// 						for(edgeIt = tmpSG->edgesBegin(); edgeIt != tmpSG->edgesEnd(); ++edgeIt)
// 						{
// 							if((*edgeIt)->label != Axon)
// 								continue;
// 							
// 							tmpBranchLengths.push_back((*edgeIt)->segmentLength());
// #ifdef DEBUG
// 							totalAxonLength += (*edgeIt)->segmentLength();
// #endif
// 						}
// 						boxSomaDistances.push_back(boxSomaDistance);
// 						boxBPCounts.push_back(tmpBPCount);
// 						boxBranchLengths.push_back(tmpBranchLengths);
						
						delete tmpSG;
					}
			
			// 	store # BPs and length with distance box center-soma
			
			std::cout << "Total axon length from boxes: " << totalAxonLength << std::endl;
			std::cout << "Total BPs from boxes: " << totalBPCount << std::endl;
			
// 			// per-piece output
// 			std::string pieceOFName(ofName);
// 			pieceOFName += "_box_analysis_pieces.csv";
// 			std::ofstream TablePerPiece;
// 			TablePerPiece.open(pieceOFName.c_str());
// 			TablePerPiece << "Piece COM-soma distance (microns)\tPiece length (microns)\n";
// 			for(int i = 0; i < pieceSomaDistances.size(); ++i)
// 			{
// 				TablePerPiece << pieceSomaDistances[i] << "\t" << pieceLengths[i] << std::endl;
// 			}
// 			TablePerPiece.close();
			
			// per-branch output
			std::string branchOFName(ofName);
			char * tmpBranchName = new char[64];
			sprintf(tmpBranchName, "_box_analysis_branches_dx_%.0f_dy_%.0f_dz_%.0f.csv", shift[X_COORD], shift[Y_COORD], shift[Z_COORD]);
			branchOFName += std::string(tmpBranchName);
			std::ofstream TablePerBranch;
			TablePerBranch.open(branchOFName.c_str());
			TablePerBranch << "Branch COM-soma distance (microns)\tBranch length (microns)\tEuclidean distance between end/branch nodes\tBPs\n";
			for(int i = 0; i < fragmentSomaDistances.size(); ++i)
			{
				TablePerBranch << fragmentSomaDistances[i] << "\t" << fragmentLengths[i] << "\t" << fragmentEuclideanLengths[i] << "\t" << fragmentBPs[i] << std::endl;
			}
			TablePerBranch.close();
			
			// long terminal branches output
			std::string terminalOFName(ofName);
			char * tmpTermName = new char[64];
			sprintf(tmpTermName, "_box_analysis_all_terminals_dx_%.0f_dy_%.0f_dz_%.0f.csv", shift[X_COORD], shift[Y_COORD], shift[Z_COORD]);
			terminalOFName += std::string(tmpTermName);
			std::ofstream TablePerTerminal;
			TablePerTerminal.open(terminalOFName.c_str());
			TablePerTerminal << "Terminal COM-soma distance (microns)\tTerminal truncated length (microns)\tTerminal non-truncated length (microns)\n";
			for(int i = 0; i < longTerminalSomaDistances.size(); ++i)
			{
				TablePerTerminal << longTerminalSomaDistances[i] << "\t" << longTerminalTruncatedLength[i] << "\t" << longTerminalNonTruncatedLength[i] << std::endl;
			}
			TablePerTerminal.close();
			
			// long terminal branches output
			std::string isolatedTerminalOFName(ofName);
			char * tmpTermName2 = new char[64];
			sprintf(tmpTermName2, "_box_analysis_isolated_terminals_dx_%.0f_dy_%.0f_dz_%.0f.csv", shift[X_COORD], shift[Y_COORD], shift[Z_COORD]);
			isolatedTerminalOFName += std::string(tmpTermName2);
			std::ofstream TablePerIsolatedTerminal;
			TablePerIsolatedTerminal.open(isolatedTerminalOFName.c_str());
			TablePerIsolatedTerminal << "Terminal COM-soma distance (microns)\tTerminal truncated length (microns)\tTerminal non-truncated length (microns)\n";
			for(int i = 0; i < isolatedTerminalSomaDistances.size(); ++i)
			{
				TablePerIsolatedTerminal << isolatedTerminalSomaDistances[i] << "\t" << isolatedTerminalTruncatedLength[i] << "\t" << isolatedTerminalNonTruncatedLength[i] << std::endl;
			}
			TablePerIsolatedTerminal.close();
			
			// per-box output
// 			ofName += "_box_analysis.csv";
// 			std::ofstream Table;
// 			Table.open(ofName.c_str());
// 			Table << "Box center-soma distance (microns)\tNr. of pieces in box\tPiece lengths\tPiece BPs\n";
// 			for(int i = 0; i < boxSomaDistances.size(); ++i)
// 			{
// 				Table << boxSomaDistances[i] << "\t" << boxBranchLengths[i].size() << "\t";
// 				for(int j = 0; j < boxBranchLengths[i].size(); ++j)
// 				{
// 					if(j > 0)
// 					{
// 						Table << ",";
// 					}
// 					Table << boxBranchLengths[i][j];
// 				}
// 				Table << "\t";
// 				for(int j = 0; j < boxBranchBPCounts[i].size(); ++j)
// 				{
// 					if(j > 0)
// 					{
// 						Table << ",";
// 					}
// 					Table << boxBranchBPCounts[i][j];
// 				}
// 				Table << std::endl;
// 			}
// 			Table.close();
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














