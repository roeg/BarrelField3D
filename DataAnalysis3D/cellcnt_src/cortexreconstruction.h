/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      cortexsreconstruction.h                                       */
/*                                                                          */
/* Purpose:   Class for 3D reconstruction of cortical landmarks from 2D     */
/*            contour stacks. Same methods as used in CortexCoordinates     */
/*            for reconstruction of the automatic landmarks                 */
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

#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/barrel_field.h"
#include "../../common/amiraReader.h"
#include "../../common/inputparameters.h"
#include "cortexsurfacereconstruct.h"

#ifndef CORTEXRECONSTRUCTION_H
#define CORTEXRECONSTRUCTION_H

class CortexReconstruction
{
	public:
		CortexReconstruction(AmiraSpatialGraph * inputSG, InputParameters parameters);
		CortexReconstruction(AmiraSpatialGraph * inputSG);
		~CortexReconstruction();
		
		void startReconstruction(const char * outputFilename);
		void convertSGToBarrelField();
		
		void vesselReconstruction(const char * outputFilename);
		
		std::map< int, Column * > getBarrelField();
		Surface * getPiaSurface();
		Surface * getWMSurface();
		
	private:
		InputParameters parameters;
		AmiraSpatialGraph * spatialGraph;
		PolyDataPointerType piaSurface;
		PolyDataPointerType wmSurface;
		std::map< int, Column * > finalBarrels;
		BarrelField * SBF;
		
		void barrelReconstruction();
		void writeBarrelParameters(const char * ofName);
		PolyDataPointerType computeAverageBarrelContour(PolyDataPointerType completeBarrel, double barrelCentroid[3], double barrelAxis[3]);
		PolyDataPointerType sampleBarrelContour(PolyDataPointerType completeBarrel);
		void computeMaxBarrelContours(std::map< int, PolyDataPointerType > barrels, std::map< int, double * > barrelAxes, std::map< int, double * > barrelCenters, std::map< int, PolyDataPointerType >& avgBarrels, std::map< int, std::vector< double * > >& endPointMap);
		PolyDataPointerType smoothBarrelAlongAxis(double newAxis[3], double barrelCentroid[3], PolyDataPointerType barrel, std::vector< double * > endPoints);
		void enforceOverlapConstraint(std::map< int, double * > barrelAxes, std::map< int, PolyDataPointerType >& avgBarrels, std::map< int, std::vector< double * > >& endPointMap);
		PolyDataPointerType createTopBottomContours(PolyDataPointerType avgContour, std::vector< double * > endPoints);
		void newBarrelAxis(PolyDataPointerType barrel, PolyDataPointerType piaSurface, std::list< unsigned int > vessels, double alpha, double axis[3]);
		std::multimap< double, double * > barrelAxisScores(PolyDataPointerType piaSurface, double * barrelCentroid, double radius, double alpha);
		PolyDataPointerType selectSurfaceRoi(PolyDataPointerType surface, double * center, double radius);
		void enforceAxisDivergence(std::map< int, double * > barrelAxes, std::map< int, double * > barrelCenters);
		std::map< int, std::list< int > > createBarrelGrid(std::map< int, double * > barrelAxes);
		void calculateBarrelCentroid(PolyDataPointerType barrel, double centroid[3]);
		void closeBarrelAlongNewAxis(double newAxis[3], double barrelCentroid[3], PolyDataPointerType barrel, std::vector< double * >& endPoints);
		void getLandmarkMinMaxIDs(PolyDataPointerType landmark, int& minID, int& maxID);
		int intersectConvexCellsInPlane(CellPointerType cell1, CellPointerType cell2, double tol, double p0[3], double p1[3]);
		
		void computeBloodVessels(bool constrained, double minRadius);
		std::multimap< double, Edge* > vesselDistances3D(double origin[3], std::pair< std::multimap< double, Edge* >::iterator, std::multimap< double, Edge* >::iterator > vessels);
		std::multimap< double, Edge* > vesselDistances2D(double origin[3], std::pair< std::multimap< double, Edge* >::iterator, std::multimap< double, Edge* >::iterator > vessels);
		std::list< unsigned int > computeConstrainingVessels(PolyDataPointerType piaSurface);
		unsigned int vesselsAroundBarrel(double barrelCentroid[3], std::list< unsigned int > vessels);
		
// 		PolyDataPointerType surfaceReconstruction(int label);
// 		ImageDataPointerType addTop2(int label, int additionalSections, bool zReversed = 0, double zSpacing = 0);
// 		ImageDataPointerType piaVolume(int label, int additionalSections, bool zReversed = 0, double zSpacing = 0);
// 		ImageDataPointerType createImageVolumeFromPolyData(PolyDataPointerType poly, int label, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax, double zSpacing = 0);
// 		ImageDataPointerType distanceTransform(ImageDataPointerType volume);
// 		int * calculateExtent(int minX, int maxX, int minY, int maxY, int minZ, int maxZ, int label);
// 		int * calculateExtent(int minX, int maxX, int minY, int maxY, int minZ, int maxZ, int label, double spacing[3]);
// 		NeighborhoodOffsetVectorType CreateLookUpTable();
// 		PolyDataPointerType smoothSurface(PolyDataPointerType surface);
};

#endif // CORTEXRECONSTRUCTION_H
