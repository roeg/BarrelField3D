/****************************************************************************/
/*                                                                          */
/* File:      basics.h										                */
/*                                                                          */
/* Purpose:   class for basic functions and	common variables				*/
/*			  especially writing and reading functions						*/
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

#pragma once
#include "typedefs.h"

#ifndef SPATIALGRAPH
#define SPATIALGRAPH


#define WATERSHED_COLOR_STEP	10
#define KMEANS_COLOR_STEP		13
#define CRITICAL_DENSITY		0.00008	// [N/�m^3]

class Vertex
{
	public:
		Vertex(float * _coordinates, int _label);
		Vertex(double * _coordinates, int _label);
		~Vertex();
// 	private:
		double coordinates[3];
		int label;
};

class Edge
{
	public:
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates, float _radius);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates, float _radius);
		~Edge();
// 	private:
		int edgeConnectivity[2];
		int numEdgePoints;
		int label;
		std::list< double * > edgePointCoordinates;
		float radius;
};

// class Basics
// {
// 	public:
// 		Basics();
// 		~Basics();
// 		
// 	protected:
// 		const char * inputFilename;						// filename of the input images
// 		const char * outputFilename;					// filename of CellClusterList.csv, image files (*.sr), landmark file and histo stuff 
// };

class AmiraSpatialGraph
{
	public:
		AmiraSpatialGraph();
		~AmiraSpatialGraph();
		
		void addVertex( Vertex * newVertex );
		void addEdge( Edge * newEdge );
		
		std::list< Vertex * >::iterator verticesBegin();
		std::list< Edge * >::iterator edgesBegin();
		std::list< Vertex * >::iterator verticesEnd();
		std::list< Edge * >::iterator edgesEnd();
		
		//extract all points of landmark 'label' planewise; return their plane indices in zIndexList; return false if empty
		bool extractLandmark(int label, std::list< std::list< double * > >& planewisePointList, std::list< int >& zIndexList);
		bool extractLandmark(int label, PolyDataPointerType polyData, std::list< int >& zIndexList);
		
		unsigned int getNumberOfVertices() { return vertices.size(); }
		unsigned int getNumberOfEdges() { return edges.size(); }
		unsigned int getNumberOfPoints();
		
		void vesselsToPoints();
	
	private:
// 		unsigned int numberOfVertices, numberOfEdges, numberOfPoints;
		
		std::list< Vertex * > vertices;
		std::list< Edge * > edges;
		
		void removeDuplicatePoints();
};

#endif