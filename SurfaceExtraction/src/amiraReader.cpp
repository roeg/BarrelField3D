/****************************************************************************/
/*                                                                          */
/* File:      amiraReader.cpp                                                */
/*                                                                          */
/* Purpose:   this program seeks to extract the exact number ond position   */
/*            of neuronal cellbodies labelled with a fluorescent dye        */
/*            from confocal or 2-Photon stacks                              */
/*                                                                          */
/* Author:    Marcel Oberlaender                                            */
/*            Max-Planck-Institute of Neurobiology                          */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/* Co-Author: Stefan Reissl                                                 */
/*            Max-Planck-Institute of Neurobiology                          */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/* Co-Author: Robert Egger                                                  */
/*            Max-Planck-Institute for Medical Research                     */
/*            Jahnstrasse 19                                                */
/*            D-69120 Heidelberg                                            */
/*                                                                          */
/* EMail:     oberlaender@neuro.mpg.de                                      */
/*                                                                          */
/* History:   03.01.2010                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "amiraReader.h"


Reader::Reader(const char * filename)
{
	this->inputFilename = filename;
	letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	numbers = "0123456789";
	signs = "+-";
	otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
	whitespace = "\t ";
};

Reader::Reader(const char * filename, const char * outputFilename)
{
	this->inputFilename = filename;
	this->outputFilename = outputFilename;
	letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	numbers = "0123456789";
	signs = "+-";
	otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
	whitespace = "\t ";
};

Reader::~Reader()
{
	//tbd
};

void Reader::readSpatialGraphFile()
{
	std::ifstream inputStream(inputFilename);
	
	if(!inputStream.fail())
	{
		std::string currentLine;
		unsigned int line = 0;
		
		bool parameters = 0;
		bool transform = 0;
		unsigned int brackets = 0;
		unsigned int vertexTransformIndex = 1000000, edgeTransformIndex = 1000000;
		unsigned int vertexCoordIndex = 1000000, vertexLabelIndex = 1000000, edgeConnectivityIndex = 1000000, edgePointIndex = 1000000, edgeLabelIndex = 1000000, edgePointCoordIndex = 1000000, edgeRadiusIndex = 1000000;
		unsigned int currentIndex = 0;
		unsigned int vertex = 0, edge = 0, point = 0;
		
		std::list< Vertex* > inputVertices;
		std::list< Edge* > inputEdges;
		std::list< double * > tmpVertices;
		std::list< int > tmpVertexLabels;
		std::list< int * > tmpEdgeConnections;
		std::list< int > tmpNoEdgePoints;
		std::list< int > tmpEdgeLabels;
		std::list< double * > edgePoints;
		
		while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
		{
// 			if(!parameters)
// 			{
// 				std::cout << currentLine << std::endl;
// 				++line;
// 			}
			
			if(currentLine.size())
			{
				if(currentLine.find("@", 0) == 0)
				{
					char * tmp = new char[currentLine.size() - 1];
					currentLine.copy(tmp, currentLine.size() - 1, 1);
					currentIndex = atoi(tmp);
					std::cout << "Reading data section " << currentIndex << std::endl;
					delete [] tmp;
					continue;
				}
				
				if(currentIndex == 0)
				{
					std::string::size_type loc = currentLine.find("define", 0);
					if(loc != std::string::npos)
					{
						if(currentLine.find("VERTEX", 7) != std::string::npos)
						{
							char * tmp = new char[currentLine.size() - 14];
							currentLine.copy(tmp, currentLine.size() - 14, 14);
							vertex = atoi(tmp);
// 							std::cout << "vertex = " << vertex << std::endl;
// 							inputVertices.resize(vertex);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("EDGE", 7) != std::string::npos)
						{
							char * tmp = new char[currentLine.size() - 12];
							currentLine.copy(tmp, currentLine.size() - 12, 12);
							edge = atoi(tmp);
// 							std::cout << "edges = " << edge << std::endl;
// 							inputEdges.resize(edge);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("POINT", 7) != std::string::npos)
						{
							char * tmp = new char[currentLine.size() - 13];
							currentLine.copy(tmp, currentLine.size() - 13, 13);
							point = atoi(tmp);
// 							std::cout << "points = " << point << std::endl;
							delete [] tmp;
							continue;
						}
					}
					
					loc = currentLine.find("Parameters", 0);
					if(loc == 0)
					{
						parameters = 1;
						brackets = 1;
						continue;
					}
					if(parameters && currentLine.find("{", 0) == std::string::npos && currentLine.find("}", 0) == std::string::npos)
						continue;
					if(parameters && currentLine.find("{", 0) != std::string::npos)
					{
						++brackets;
						continue;
					}
					if(parameters && currentLine.find("}", 0) != std::string::npos)
					{
						--brackets;
						if(!brackets)
							parameters = 0;
						continue;
					}
					
					loc = currentLine.find("VERTEX", 0);
					if(loc == 0)
					{
						if(currentLine.find("VertexCoordinates", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							vertexCoordIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("GraphLabels", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							vertexLabelIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("TransformInfo", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							vertexTransformIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
					}
					
					loc = currentLine.find("EDGE", 0);
					if(loc == 0)
					{
						if(currentLine.find("EdgeConnectivity", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeConnectivityIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("NumEdgePoints", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgePointIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("GraphLabels", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeLabelIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("TransformInfo", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeTransformIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
					}
					
					loc = currentLine.find("POINT", 0);
					if(loc == 0)
					{
						if(currentLine.find("EdgePointCoordinates", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgePointCoordIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("Radius", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeRadiusIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
					}
				}
				
				if(currentIndex == vertexTransformIndex || currentIndex == edgeTransformIndex || currentIndex == edgeRadiusIndex)
					continue;
				
				if(currentIndex == vertexCoordIndex)
				{
					std::string::size_type loc1 = currentLine.find(" ", 0);
					std::string::size_type loc2 = currentLine.find(" ", loc1 + 1);
					char * tmp1 = new char[loc1];
					char * tmp2 = new char[loc2 - loc1];
					char * tmp3 = new char[currentLine.size() - loc2 - 1];
					currentLine.copy(tmp1, loc1, 0);
					currentLine.copy(tmp2, loc2 - loc1, loc1 + 1);
					currentLine.copy(tmp3, currentLine.size() - loc2 - 1, loc2 + 1);
					double * tmpCoords = new double[3];
					tmpCoords[0] = atof(tmp1);
					tmpCoords[1] = atof(tmp2);
					tmpCoords[2] = atof(tmp3);
					
					tmpVertices.push_back(tmpCoords);
				}
				
				if(currentIndex == vertexLabelIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					int tmplabel = atoi(tmp);
					
					tmpVertexLabels.push_back(tmplabel);
				}
				
				if(currentIndex == edgeConnectivityIndex)
				{
					std::string::size_type loc = currentLine.find(" ", 0);
					char * tmp1 = new char[loc];
					char * tmp2 = new char[currentLine.size() - loc - 1];
					currentLine.copy(tmp1, loc, 0);
					currentLine.copy(tmp2, currentLine.size() - loc - 1, loc + 1);
					int * tmpCoords = new int[2];
					tmpCoords[0] = atoi(tmp1);
					tmpCoords[1] = atoi(tmp2);
					
					tmpEdgeConnections.push_back(tmpCoords);
				}
				
				if(currentIndex == edgePointIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					int tmplabel = atoi(tmp);
					
					tmpNoEdgePoints.push_back(tmplabel);
				}
				
				if(currentIndex == edgeLabelIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					int tmplabel = atoi(tmp);
					
					tmpEdgeLabels.push_back(tmplabel);
				}
				
				if(currentIndex == edgePointCoordIndex)
				{
					std::string::size_type loc1 = currentLine.find(" ", 0);
					std::string::size_type loc2 = currentLine.find(" ", loc1 + 1);
					char * tmp1 = new char[loc1];
					char * tmp2 = new char[loc2 - loc1];
					char * tmp3 = new char[currentLine.size() - loc2 - 1];
					currentLine.copy(tmp1, loc1, 0);
					currentLine.copy(tmp2, loc2 - loc1, loc1 + 1);
					currentLine.copy(tmp3, currentLine.size() - loc2 - 1, loc2 + 1);
					double * tmpCoords = new double[3];
					tmpCoords[0] = atof(tmp1);
					tmpCoords[1] = atof(tmp2);
					tmpCoords[2] = atof(tmp3);
					
					edgePoints.push_back(tmpCoords);
				}
			}
		}
		
		std::list< double * >::iterator tmpvertexit;
		std::list< int >::iterator tmpvertexlabelit;
		std::list< int * >::iterator tmpedgeconnectivityit;
		std::list< int >::iterator tmpnumberedgepointsit;
		std::list< int >::iterator tmpedgelabelit;
		std::list< double * >::iterator edgepointit;
		
		for(tmpvertexit = tmpVertices.begin(), tmpvertexlabelit = tmpVertexLabels.begin(); tmpvertexit != tmpVertices.end() && tmpvertexlabelit != tmpVertexLabels.end(); ++tmpvertexit, ++tmpvertexlabelit)
		{
			Vertex * tmpVertex = new Vertex(*tmpvertexit, *tmpvertexlabelit);
			inputVertices.push_back(tmpVertex);
		}
		
		edgepointit = edgePoints.begin();
		for(tmpedgeconnectivityit = tmpEdgeConnections.begin(), tmpnumberedgepointsit = tmpNoEdgePoints.begin(), tmpedgelabelit = tmpEdgeLabels.begin();
		tmpedgeconnectivityit != tmpEdgeConnections.end() && tmpnumberedgepointsit != tmpNoEdgePoints.end() && tmpedgelabelit != tmpEdgeLabels.end();
		++tmpedgeconnectivityit, ++tmpnumberedgepointsit, ++tmpedgelabelit)
		{
			std::list< double * >::iterator tmpit = edgepointit;
			for(int ii = 0; ii < *tmpnumberedgepointsit; ++ii)
				++tmpit;
			std::list< double * > tmpPoints(edgepointit, tmpit);
			Edge * tmpEdge = new Edge(*tmpedgeconnectivityit, *tmpnumberedgepointsit, *tmpedgelabelit, tmpPoints);
			inputEdges.push_back(tmpEdge);
			edgepointit = tmpit;
		}
		
		inputSpatialGraph = new AmiraSpatialGraph;
		
		std::list< Vertex* >::iterator vertexIter;
		std::list< Edge* >::iterator edgeIter;
		
		for(vertexIter = inputVertices.begin(); vertexIter != inputVertices.end(); ++ vertexIter)
			inputSpatialGraph->addVertex(*vertexIter);
		
		for(edgeIter = inputEdges.begin(); edgeIter != inputEdges.end(); ++edgeIter)
			inputSpatialGraph->addEdge(*edgeIter);
		
// 		std::cout << "Vertex number = " << inputSpatialGraph->getNumberOfVertices() << std::endl;
// 		std::cout << "Edge number = " << inputSpatialGraph->getNumberOfEdges() << std::endl;
// 		std::cout << "Point number = " << inputSpatialGraph->getNumberOfPoints() << std::endl;
		
// 		std::cout << "VertexCoordinates @" << vertexCoordIndex << std::endl;
// 		std::cout << "Vertex GraphLabels @" << vertexLabelIndex << std::endl;
// 		std::cout << "Vertex TransformInfo @" << vertexTransformIndex << std::endl;
// 		std::cout << "EdgeConnectivity @" << edgeConnectivityIndex << std::endl;
// 		std::cout << "NumEdgePoints @" << edgePointIndex << std::endl;
// 		std::cout << "Edge GraphLabels @" << edgeLabelIndex << std::endl;
// 		std::cout << "Edge TransformInfo @" << edgeTransformIndex << std::endl;
// 		std::cout << "EdgePointCoordinates @" << edgePointCoordIndex << std::endl;
// 		std::cout << "Edge Radius @" << edgeRadiusIndex << std::endl;
	}
	
	inputStream.close();
};

void Reader::writeSpatialGraphFile()
{
// 	std::list<std::list<Compartment * > >::iterator edge_list_it;
// 	std::list<Compartment * >::iterator edge_it;
	std::list< Vertex* >::iterator vertexIt;
	std::list< Edge* >::iterator edgeIt;
	
	int number_of_edge_points = inputSpatialGraph->getNumberOfPoints();
	
// 	for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
// 	{
// 		number_of_edge_points += (*edge_list_it).size();
// 	}
// 	for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
// 	{
// 		number_of_edge_points += (*edge_list_contour_it).size();
// 	}
// 	for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
// 	{
// 		number_of_edge_points += (*edge_list_contour_it).size();
// 	}
	
	std::string format = outputFilename;
	format += ".am";
	
	#ifdef DEBUG
	std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
	//std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
	#endif
	
	std::ofstream NeuroMorphData( format.c_str() );
	
	NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	NeuroMorphData << "# This SpatilaGraph file was created by the Neuron Reconstruction Tool NeuroMorph " << std::endl;
	NeuroMorphData << "# NeuroMorph was programmed by Marcel Oberlaender and Philip J. Broser," << std::endl;
	NeuroMorphData << "# Max-Planck-Institute for Medical Research Heidelberg, Germany " << std::endl;
	
	NeuroMorphData << "define VERTEX " << /*amira_spatial_graph->vertice_list.size() +*/ /*amira_contour_graph->vertice_list.size() +*/ inputSpatialGraph->getNumberOfVertices() << std::endl;
	NeuroMorphData << "define EDGE " << /*amira_spatial_graph->edge_list.size() +*/ /*amira_contour_graph->edge_list.size() +*/ inputSpatialGraph->getNumberOfEdges()  << std::endl;
// 	NeuroMorphData << "define GRAPH " << /*amira_spatial_graph->vertice_list.size() +*/ /*amira_contour_graph->vertice_list.size() +*/ inputSpatialGraph->getNumberOfVertices() + /*amira_spatial_graph->edge_list.size() +*/ /*amira_contour_graph->edge_list.size() +*/ inputSpatialGraph->getNumberOfEdges() << std::endl;
	NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
	
	NeuroMorphData << "Parameters {VertexLabels"<< std::endl;
	NeuroMorphData << "           {UndefinedVertexLabel{Color 1 1 1}"  	<<std::endl;
	NeuroMorphData << "            LowEndingVertexLabel{Color 0 0 1}"   	<<std::endl;
	NeuroMorphData << "            HighEndingVertexLabel{Color 0 1 0}"  	<<std::endl;
	NeuroMorphData << "            IntersecVertexLabel{Color 1 1 0}"    	<<std::endl; 
	NeuroMorphData << "            NormalEndingVertexLabel{Color 1 0 0}}"	<<std::endl; 
	NeuroMorphData << "            EdgeLabels"                          	<<std::endl;
	NeuroMorphData << "           {UndefinedEdgeLabel{Color 1 1 1}"     	<<std::endl;
	NeuroMorphData << "            BottomEdgeLabel{Color 0 0 1}"        	<<std::endl;
	NeuroMorphData << "            TopEdgeLabel{Color 0 1 0}"           	<<std::endl;
	NeuroMorphData << "            IntermediateEdgeLabel{Color 1 1 0}"  	<<std::endl;
	NeuroMorphData << "            TopToBottomEdgeLabel{Color 1 0 0}}"  	<<std::endl;
	NeuroMorphData << "    GraphLabels {"                                	<<std::endl;
	NeuroMorphData << "        Neuron { "                                	<<std::endl;
	NeuroMorphData << "            Dendrite {"                           	<<std::endl;
	NeuroMorphData << "                ApicalDendrite {"                 	<<std::endl;
	NeuroMorphData << "                    Color 1 0.5 0.5,"          	<<std::endl;
	NeuroMorphData << "                    Id 4 }"                     	<<std::endl;
	NeuroMorphData << "                BasalDendrite {"         		<<std::endl;
	NeuroMorphData << "                    Color 0.8 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                    Id 5 }"				<<std::endl;
	NeuroMorphData << "                Color 1 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 3 }"				<<std::endl;
	NeuroMorphData << "            Axon {"					<<std::endl;
	NeuroMorphData << "                Color 0 0 1,"			<<std::endl;
	NeuroMorphData << "                Id 6 }"				<<std::endl;
	NeuroMorphData << "            Soma {"					<<std::endl;
	NeuroMorphData << "                Color 1 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 7 }"				<<std::endl;
	NeuroMorphData << "            Color 1 0 0,"				<<std::endl;
	NeuroMorphData << "            Id 2 }"					<<std::endl;
	NeuroMorphData << "        Landmark {"					<<std::endl;
	NeuroMorphData << "            Pia {"					<<std::endl;
	NeuroMorphData << "                Color 0 1 0.5,"			<<std::endl;
	NeuroMorphData << "                Id 9 }"				<<std::endl;
	NeuroMorphData << "            Vessel {"				<<std::endl;
	NeuroMorphData << "                Color 1 0.5 0,"			<<std::endl;
	NeuroMorphData << "                Id 10 }"				<<std::endl;
	NeuroMorphData << "            Barrel {"				<<std::endl;
	NeuroMorphData << "                aRow {"				<<std::endl;
	NeuroMorphData << "                    A1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 13 }"			<<std::endl;
	NeuroMorphData << "                    A2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 14 }"			<<std::endl;
	NeuroMorphData << "                    A3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 15 }"			<<std::endl;
	NeuroMorphData << "                    A4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 16 }"			<<std::endl;
	NeuroMorphData << "                Color 1 0.2 0.2,"			<<std::endl;
	NeuroMorphData << "                Id 12 }"				<<std::endl;
	NeuroMorphData << "                bRow {"				<<std::endl;
	NeuroMorphData << "                    B1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 18 }"			<<std::endl;
	NeuroMorphData << "                    B2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 19 }"			<<std::endl;
	NeuroMorphData << "                    B3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 20 }"			<<std::endl;
	NeuroMorphData << "                    B4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 21 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                    Id 17 }"				<<std::endl;
	NeuroMorphData << "                cRow {"				<<std::endl;
	NeuroMorphData << "                    C1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 23 }"			<<std::endl;
	NeuroMorphData << "                    C2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 24 }"			<<std::endl;
	NeuroMorphData << "                    C3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 25 }"			<<std::endl;
	NeuroMorphData << "                    C4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 26 }"			<<std::endl;
	NeuroMorphData << "                    C5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 27 }"			<<std::endl;
	NeuroMorphData << "                    C6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 28 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                    Id 22 }"				<<std::endl;
	NeuroMorphData << "                dRow {"				<<std::endl;
	NeuroMorphData << "                    D1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 30 }"			<<std::endl;
	NeuroMorphData << "                    D2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 31 }"			<<std::endl;
	NeuroMorphData << "                    D3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 32 }"			<<std::endl;
	NeuroMorphData << "                    D4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 33 }"			<<std::endl;
	NeuroMorphData << "                    D5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 34 }"			<<std::endl;
	NeuroMorphData << "                    D6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 35 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                    Id 29 }"				<<std::endl;
	NeuroMorphData << "                eRow {"				<<std::endl;
	NeuroMorphData << "                    E1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 37 }"			<<std::endl;
	NeuroMorphData << "                    E2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 38 }"			<<std::endl;
	NeuroMorphData << "                    E3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 39 }"			<<std::endl;
	NeuroMorphData << "                    E4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 40 }"			<<std::endl;
	NeuroMorphData << "                    E5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 41 }"			<<std::endl;
	NeuroMorphData << "                    E6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 42 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                    Id 36 }"				<<std::endl;
	NeuroMorphData << "                greekRow {"				<<std::endl;
	NeuroMorphData << "                    Alpha {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 44 }"			<<std::endl;
	NeuroMorphData << "                    Beta {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 45 }"			<<std::endl;
	NeuroMorphData << "                    Gamma {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 46 }"			<<std::endl;
	NeuroMorphData << "                    Delta {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 47 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                    Id 43 }"				<<std::endl;
	NeuroMorphData << "                Color 0 1 0,"			<<std::endl;
	NeuroMorphData << "                Id 11 }"				<<std::endl;
	NeuroMorphData << "            WhiteMatter {"				<<std::endl;
	NeuroMorphData << "                Color 0.5 1 0.75,"			<<std::endl;
	NeuroMorphData << "                Id 48 }"				<<std::endl;
	NeuroMorphData << "            OtherBarrels {"				<<std::endl;
	NeuroMorphData << "                Color 1 0 1,"			<<std::endl;
	NeuroMorphData << "                Id 49 }"				<<std::endl;
	NeuroMorphData << "            Color 0 1 1,"				<<std::endl;
	NeuroMorphData << "            Id 8 }"					<<std::endl;
	NeuroMorphData << "        Id 0,"					<<std::endl;
	NeuroMorphData << "        Color 0 0 0 }"				<<std::endl;
	NeuroMorphData << "ContentType \"HxSpatialGraph\" }"                 	<<std::endl;
	
	NeuroMorphData << "VERTEX { float[3] VertexCoordinates } = @1 " 	<< std::endl;
	NeuroMorphData << "VERTEX { int VertexLabels } @2 " 			<< std::endl;
	NeuroMorphData << "VERTEX {int GraphLabels } = @3 " 			<< std::endl;
	
	NeuroMorphData << "EDGE { int[2] EdgeConnectivity } = @4 " 		<< std::endl;
	NeuroMorphData << "EDGE { int NumEdgePoints } = @5 " 			<< std::endl;
	NeuroMorphData << "EDGE { int EdgeLabels } @6 " 			<< std::endl;
	NeuroMorphData << "EDGE { int GraphLabels } = @7 " 			<< std::endl;
	
	NeuroMorphData << "POINT { float[3] EdgePointCoordinates } = @8 " 	<< std::endl;
	NeuroMorphData << "POINT { float Radius } = @9 " 			<< std::endl;
	
// 	NeuroMorphData << "GRAPH { int GraphLabel} = @8 " 			<< std::endl;
	
	NeuroMorphData << "\n@1 # Vertices xyz coordinates" 			<< std::endl;
// 	for(int i = 0; i < amira_spatial_graph->vertice_list.size(); i++)
// 		NeuroMorphData << XYSAMPLING * amira_spatial_graph->vertice_list[i]->GetMeasurementVector()[X_COORD] << " " << XYSAMPLING * amira_spatial_graph->vertice_list[i]->GetMeasurementVector()[Y_COORD]  << " " << ZSAMPLING * amira_spatial_graph->vertice_list[i]->GetMeasurementVector()[Z_COORD]  << std::endl;
// 	for(contour_it=amira_contour_graph->vertice_list.begin(); contour_it!=amira_contour_graph->vertice_list.end(); ++contour_it)
// 		NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD]  << " " << (*contour_it)[Z_COORD]  << std::endl;
	for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
		NeuroMorphData << (*vertexIt)->coordinates[X_COORD] << " " << (*vertexIt)->coordinates[Y_COORD]  << " " << (*vertexIt)->coordinates[Z_COORD]  << std::endl;
	
	NeuroMorphData << "\n@2 # Vertex Label" << std::endl;
// 	for(int i = 0; i < amira_spatial_graph->vertice_list.size(); i++)
// 	{
// 		if(amira_spatial_graph->vertice_list[i]->FLAG(LOWENDING))
// 			NeuroMorphData << 1 << std::endl;
// 		else
// 		{
// 			if(amira_spatial_graph->vertice_list[i]->FLAG(HIGHENDING))
// 				NeuroMorphData << 2 << std::endl;
// 			else
// 				if(amira_spatial_graph->vertice_list[i]->FLAG(INTERSECTION))
// 					NeuroMorphData << 3 << std::endl;
// 				else
// 					if(amira_spatial_graph->vertice_list[i]->FLAG(IRREDUCEABLEPOINTS))
// 						NeuroMorphData << 4 << std::endl;
// 		}
// 	}
// 	for(int i = 0; i < amira_contour_graph->vertice_list.size(); i++)
// 	{
// 		NeuroMorphData << 1 << std::endl;
// 	}
	for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
	{
		NeuroMorphData << 0 << std::endl;
	}
	
	NeuroMorphData << "\n@3 # Vertex Graph Label" << std::endl;
// 	for(int i = 0; i < amira_spatial_graph->vertice_list.size(); i++)
// 	{
// 		if(amira_spatial_graph->vertice_list[i]->FLAG(LOWENDING))
// 			NeuroMorphData << 1 << std::endl;
// 		else
// 		{
// 			if(amira_spatial_graph->vertice_list[i]->FLAG(HIGHENDING))
// 				NeuroMorphData << 2 << std::endl;
// 			else
// 				if(amira_spatial_graph->vertice_list[i]->FLAG(INTERSECTION))
// 					NeuroMorphData << 3 << std::endl;
// 				else
// 					if(amira_spatial_graph->vertice_list[i]->FLAG(IRREDUCEABLEPOINTS))
// 						NeuroMorphData << 4 << std::endl;
// 		}
// 	}
// 	for(int i = 0; i < amira_contour_graph->vertice_list.size(); i++)
// 	{
// 		NeuroMorphData << 1 << std::endl;
// 	}
	for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
	{
		NeuroMorphData << (*vertexIt)->label << std::endl;
	}
	
	NeuroMorphData << "\n@4 # Edge Identifiers" << std::endl;
	int last_index = 0;
// 	for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
// 	{	  
// 		int ID1 = 0;
// 		int ID2 = 0;
// 		int index1 = 0; 
// 		int index2 = 0;
// 		
// 		edge_it = (*edge_list_it).begin();
// 		ID1 = (int)((*edge_it)->GetMeasurementVector()[IDENTIFIER]);
// 		
// 		edge_it = (*edge_list_it).end();
// 		edge_it--;
// 		
// 		ID2 = (int)((*edge_it)->GetMeasurementVector()[IDENTIFIER]);
// 		//std::cout<< "Edge Indices: " << ID1 << " " << ID2 << std::endl;
// 		
// 		
// 		for(int i = 0; i < amira_spatial_graph->vertice_list.size(); i++)
// 		{
// 			if((int)(amira_spatial_graph->vertice_list[i]->GetMeasurementVector()[IDENTIFIER])==ID1)
// 				index1 = i;
// 		}
// 		for(int j = 0; j < amira_spatial_graph->vertice_list.size(); j++)
// 		{
// 			if((int)(amira_spatial_graph->vertice_list[j]->GetMeasurementVector()[IDENTIFIER])==ID2)
// 				index2 = j;
// 		}
// 		NeuroMorphData <<  index1 << " " << index2 <<std::endl;
// 		last_index=index2;
// 	}
	
// 	int last_index = amira_contour_graph->vertice_list.size();
// 	for(int i=0; i < amira_contour_graph->vertice_list.size(); i++)
// 	{
// 		NeuroMorphData << /*last_index+i+1*/ i << " " << /*last_index+i+1*/ i << std::endl;
// 	}
	for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
	{
		NeuroMorphData << (*edgeIt)->edgeConnectivity[0] << " " << (*edgeIt)->edgeConnectivity[1] << std::endl;
	}
	
	NeuroMorphData << "\n@5 # Number of Points per Edge" << std::endl;
// 	for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
// 	{
// 		int number_of_points = 0;
// 		
// 		number_of_points = (*edge_list_it).size();
// 		
// 		NeuroMorphData <<  number_of_points <<std::endl;
// 	}
// 	for(edge_list_contour_it=amira_contour_graph->edge_list.begin(); edge_list_contour_it!=amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
// 	{
// 		int number_of_contur_points = 0;
// 		
// 		number_of_contur_points = (*edge_list_contour_it).size();
// 		
// 		NeuroMorphData <<  number_of_contur_points <<std::endl;
// 	}
	for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
	{
// 		int number_of_contur_points = 0;
// 		
// 		number_of_contur_points = (*edge_list_contour_it).size();
		
		NeuroMorphData <<  (*edgeIt)->numEdgePoints <<std::endl;
	}
	
	NeuroMorphData << "\n@6 # EdgeLabels" << std::endl;
// 	for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
// 	{	 
// 		edge_it = (*edge_list_it).begin();
// 		
// 		switch((*edge_it)->GetColorTag())
// 		{
// 			case blue:   NeuroMorphData <<  1 <<std::endl;break; //low
// 			case green:  NeuroMorphData <<  2 <<std::endl;break; //high
// 			case yellow: NeuroMorphData <<  3 <<std::endl;break; //nothing
// 			case red:    NeuroMorphData <<  4 <<std::endl;break; //top to bottom
// 			default:     NeuroMorphData <<  0 <<std::endl;break; //nothing
// 		}
// 	}
// 	for(edge_list_contour_it=amira_contour_graph->edge_list.begin(); edge_list_contour_it!=amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
// 	{
// 		NeuroMorphData << 1 << std::endl;
// 	}
	for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
	{
		NeuroMorphData << 0 << std::endl;
	}
	
	NeuroMorphData << "\n@7 # Edge Graph Labels" << std::endl;
// 	for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
// 	{	 
// 		edge_it = (*edge_list_it).begin();
// 		
// 		switch((*edge_it)->GetColorTag())
// 		{
// 			case blue:   NeuroMorphData <<  1 <<std::endl;break; //low
// 			case green:  NeuroMorphData <<  2 <<std::endl;break; //high
// 			case yellow: NeuroMorphData <<  3 <<std::endl;break; //nothing
// 			case red:    NeuroMorphData <<  4 <<std::endl;break; //top to bottom
// 			default:     NeuroMorphData <<  0 <<std::endl;break; //nothing
// 		}
// 	}
// 	for(edge_list_contour_it=amira_contour_graph->edge_list.begin(); edge_list_contour_it!=amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
// 	{
// 		NeuroMorphData << 1 << std::endl;
// 	}
	for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
	{
		NeuroMorphData << (*edgeIt)->label << std::endl;
	}
	
	NeuroMorphData << "\n@8 # Point xyz coordinates" << std::endl;
// 	for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
// 	{	  
// 		for(edge_it = (*edge_list_it).begin(); edge_it != (*edge_list_it).end(); ++edge_it)
// 		{
// 			NeuroMorphData << XYSAMPLING * (*edge_it)->GetMeasurementVector()[X_COORD] << " " << XYSAMPLING * (*edge_it)->GetMeasurementVector()[Y_COORD] << " " << ZSAMPLING * (*edge_it)->GetMeasurementVector()[Z_COORD] << std::endl;
// 		}
// 		
// 	}
// 	for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
// 	{	  
// 		for(contour_it = (*edge_list_contour_it).begin(); contour_it != (*edge_list_contour_it).end(); ++contour_it)
// 		{
// 			NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
// 		}
// 		
// 	}
	for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
	{
		std::list< double * >::iterator contour_it;
		for(contour_it = (*edgeIt)->edgePointCoordinates.begin(); contour_it != (*edgeIt)->edgePointCoordinates.end(); ++contour_it)
		{
			NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
		}
		
	}
	
	NeuroMorphData << "\n@9 # Radius at Point" << std::endl;
// 	for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
// 	{	  
// 		for(edge_it = (*edge_list_it).begin(); edge_it != (*edge_list_it).end(); ++edge_it)
// 		{
// 			NeuroMorphData << (*edge_it)->GetMeasurementVector()[SURFACE] << std::endl;
// 		}
// 		
// 	}
// 	for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
// 	{	  
// 		for(contour_it = (*edge_list_contour_it).begin(); contour_it != (*edge_list_contour_it).end(); ++contour_it)
// 		{
// 			NeuroMorphData << 1 << std::endl;
// 		}
// 		
// 	}
	for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
	{	  
		for(int ii = 0; ii < (*edgeIt)->edgePointCoordinates.size(); ++ii)
		{
			NeuroMorphData << (*edgeIt)->radius << std::endl;
		}
		
	}
	
// 	NeuroMorphData << "\n@8 # Graph Property" << std::endl;
// // 	for(int i = 0; i < amira_spatial_graph->vertice_list.size(); i++)
// // 		NeuroMorphData << 6 <<std::endl;
// // 	for(contour_it=amira_contour_graph->vertice_list.begin(); contour_it!=amira_contour_graph->vertice_list.end(); ++contour_it)
// // 		NeuroMorphData << 10 << std::endl;
// 	for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
// 		NeuroMorphData << (*vertexIt)->label << std::endl;
// // 	for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
// // 		NeuroMorphData << 6 << std::endl;
// // 	for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
// // 		NeuroMorphData << 10 << std::endl;
// 	for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
// 		NeuroMorphData << (*edgeIt)->label << std::endl;
	
// 	return 0;
	NeuroMorphData.close();
};

void getTransformationsFromSpatialGraphFile()
{
	//tbd
};

void Reader::writeAmiraSurfaceFile(PolyDataPointerType triangleData)
{
	std::string format = outputFilename;
	format += ".surf";
	
	std::ofstream SurfaceData( format.c_str() );
	
	SurfaceData << "# HyperSurface 0.1 ASCII" << std::endl;
	SurfaceData << "" << std::endl;
	SurfaceData << "Parameters {" << std::endl;
	SurfaceData << "\tMaterials {" << std::endl;
	SurfaceData << "\t\tExterior {" << std::endl;
	SurfaceData << "\t\t\tid 0," << std::endl;
	SurfaceData << "\t\t\tColor 1 1 1" << std::endl;
	SurfaceData << "\t\t}" << std::endl;
	SurfaceData << "\t\tUnsortedContours0 {" << std::endl;
	SurfaceData << "\t\t\tid 1," << std::endl;
	SurfaceData << "\t\t\tColor 1 0 0" << std::endl;
	SurfaceData << "\t\t}" << std::endl;
	SurfaceData << "\t}" << std::endl;
	SurfaceData << "\tBoundaryIds {" << std::endl;
	SurfaceData << "\t\tname \"BoundaryConditions\"" << std::endl;
	SurfaceData << "\t}" << std::endl;
	SurfaceData << "}" << std::endl;
	SurfaceData << "" << std::endl;
	
	SurfaceData << "Vertices " << triangleData->GetNumberOfPoints() << std::endl;
	for(int ii = 0; ii < triangleData->GetNumberOfPoints(); ++ii)
	{
		double * point = triangleData->GetPoint(ii);
		SurfaceData << std::fixed << "\t" << point[0] << " " << point[1] << " " << point[2] << std::endl;
	}
	
	SurfaceData << "NBranchingPoints " << 0 << std::endl;
	SurfaceData << "NVerticesOnCurves " << 0 << std::endl;
	SurfaceData << "BoundaryCurves " << 0 << std::endl;
	SurfaceData << "Patches " << 1 << std::endl;
	SurfaceData << "{" << std::endl;
	SurfaceData << "InnerRegion UnsortedContours" << 0 << std::endl;
	SurfaceData << "OuterRegion Exterior" << std::endl;
	SurfaceData << "BoundaryID " << 0 << std::endl;
	SurfaceData << "BranchingPoints " << 0 << std::endl;
	SurfaceData << "" << std::endl;
	SurfaceData << "Triangles " << triangleData->GetNumberOfPolys() << std::endl;
	for(unsigned int ii = 0; ii < triangleData->GetNumberOfPolys(); ++ii)
	{
		vtkIdList * pointIDs = triangleData->GetCell(ii)->GetPointIds();
		SurfaceData << "\t" << pointIDs->GetId(0) + 1 << " " << pointIDs->GetId(1) + 1 << " " << pointIDs->GetId(2) + 1 << std::endl;
	}
	SurfaceData << "}" << std::endl;
	
	SurfaceData.close();
};

PolyDataPointerType Reader::readAmiraSurfaceFile()
{
	std::ifstream inputStream(inputFilename);
	PolyDataPointerType surface = PolyDataPointerType::New();
	PointsPointerType points = PointsPointerType::New();
	
	if(!inputStream.fail())
	{
		std::string currentLine;
		
		unsigned int currentIndex = 0;
		const unsigned int point = 1, cell = 2;
		unsigned int pointID = 0, cellID = 0;
		
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				std::string::size_type loc1, loc2, loc3;
				
				if(!currentIndex)
				{
					loc1 = currentLine.find("Vertices", 0);
					if(loc1 == 0)
					{
						
						currentIndex = point;
						char * tmp = new char[currentLine.size() - 9];
						currentLine.copy(tmp, currentLine.size() - 9, 9);
						int noOfPoints = atoi(tmp);
						points->SetDataTypeToFloat();
						points->SetNumberOfPoints(noOfPoints);
						delete [] tmp;
					}
					
					loc1 = currentLine.find("Triangles", 0);
					if(loc1 == 0)
					{
						currentIndex = cell;
						char * tmp = new char[currentLine.size() - 10];
						currentLine.copy(tmp, currentLine.size() - 10, 10);
						int noOfCells = atoi(tmp);
						surface->Allocate(noOfCells);
						delete [] tmp;
					}
				}
				
				else if(currentIndex == point || currentIndex == cell)
				{
					if(currentLine.find_first_of(letters, 0) != std::string::npos || currentLine.find_first_of(otherChars, 0) != std::string::npos)
					{
						currentIndex = 0;
						continue;
					}
					
					if(currentIndex == point)
					{
						loc1 = currentLine.find_first_of(numbers, 0);
						loc2 = currentLine.find_first_of(signs, 0);
						if(loc2 != std::string::npos)
							if(loc2 < loc1)
								loc1 = loc2;
						loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
						loc3 = currentLine.find_first_of(whitespace, loc2 + 1);
						char * tmp1 = new char[loc2 - loc1];
						char * tmp2 = new char[loc3 - loc2 - 1];
						char * tmp3 = new char[currentLine.size() - loc3 - 1];
						currentLine.copy(tmp1, loc2 - loc1, loc1);
						currentLine.copy(tmp2, loc3 - loc2 - 1, loc2 + 1);
						currentLine.copy(tmp3, currentLine.size() - loc3 - 1, loc3 + 1);
						double pointCoords[] = {atof(tmp1), atof(tmp2), atof(tmp3)};
						points->SetPoint(pointID, pointCoords);
						++pointID;
						delete [] tmp1;
						delete [] tmp2;
						delete [] tmp3;
					}
					if(currentIndex == cell)
					{
						loc1 = currentLine.find_first_of(numbers, 0);
						loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
						loc3 = currentLine.find_first_of(whitespace, loc2 + 1);
						char * tmp1 = new char[loc2 - loc1];
						char * tmp2 = new char[loc3 - loc2 - 1];
						char * tmp3 = new char[currentLine.size() - loc3 - 1];
						currentLine.copy(tmp1, loc2 - loc1, loc1);
						currentLine.copy(tmp2, loc3 - loc2 - 1, loc2 + 1);
						currentLine.copy(tmp3, currentLine.size() - loc3 - 1, loc3 + 1);
						int cellPoints[] = {atoi(tmp1) - 1, atoi(tmp2) - 1, atoi(tmp3) - 1};
						PolygonPointerType poly = PolygonPointerType::New();
						poly->GetPointIds()->SetNumberOfIds(3);
						for(int ii = 0; ii < 3; ++ii)
							poly->GetPointIds()->SetId(ii, cellPoints[ii]);
						surface->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
						delete [] tmp1;
						delete [] tmp2;
						delete [] tmp3;
					}
				}
			}
		}
		surface->SetPoints(points);
		surface->Update();
	}
	inputStream.close();
	return surface;
};


void Reader::getTransformationsFromSpatialGraphFile()
{
	
};

























