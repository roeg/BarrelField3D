/****************************************************************************/
/*                                                                          */
/* Program:   CortexCoordinates                                             */
/*                                                                          */
/* File:      geometry.h                                                    */
/*                                                                          */
/* Purpose:   class providing all methods for computing geometric           */
/*            parameters and structures: 3D vessels, surfaces, barrel       */
/*            parameters                                                    */
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

#pragma once
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/amiraReader.h"

#ifndef GEOMETRY
#define GEOMETRY

class Geometry
{
	public:
		Geometry(AmiraSpatialGraph * inputSpatialGraph);
		Geometry();
		~Geometry();
		
		void computeVesselToNormalAngleHistogram(const char * graphFName, const char * surfFName, const char * oFName, double binSize);
		void computeSurfaces(const char * filename, int interpolParam);
		PolyDataPointerType smoothSurface(PolyDataPointerType surface);
		
		void computeBarrelSurfaces(const char * filename, double alpha);
		void computeTotalVolumes(const char * filename);
		
		void reconstructionError(const char * filename);
		
		AmiraSpatialGraph * getSpatialGraph() { return spatialGraph; }
		AmiraSpatialGraph * getSpatialGraph2() { return spatialGraph2; }
		
	private:
		std::list< int > barrelLabels;
		std::map< int, const char * > int2Labels;
		AmiraSpatialGraph * spatialGraph;
		AmiraSpatialGraph * spatialGraph2;
		
		void computeBloodVessels(bool constrained);
		std::multimap< double, Edge* > vesselDistances3D(double * origin, std::pair< std::multimap< double, Edge* >::iterator, std::multimap< double, Edge* >::iterator > vessels);
		std::multimap< double, Edge* > vesselDistances2D(double * origin, std::pair< std::multimap< double, Edge* >::iterator, std::multimap< double, Edge* >::iterator > vessels);
		
		ImageDataPointerType piaVolume(int label, const char * filename);
		ImageDataPointerType piaVolume(int label, int additionalSections);
		ImageDataPointerType distanceTransform(ImageDataPointerType volume);
		int * calculateExtent(int minX, int maxX, int minY, int maxY, int minZ, int maxZ, int label);
		NeighborhoodOffsetVectorType CreateLookUpTable();
		
		PolyDataPointerType smoothBarrelInZ();
		PolyDataPointerType smoothBarrelInZ2(int label);
		PolyDataPointerType correctOutlierContours(PolyDataPointerType barrel);
// 		std::vector< double > newBarrelAxis(int label, std::list< unsigned int > vessels, PolyDataPointerType piaSurface, const char * outputFilename, double alpha);
		double * newBarrelAxis(PolyDataPointerType barrel, std::list< unsigned int > vessels, PolyDataPointerType piaSurface, const char * outputFilename, double alpha);
		double * simpleBarrelAxis(PolyDataPointerType barrel, PolyDataPointerType piaSurface, double alpha);
		std::list< unsigned int > computeConstrainingVessels(PolyDataPointerType piaSurface);
		std::multimap< double, double * > barrelAxisScores(PolyDataPointerType piaSurface, double * barrelCentroid, double radius, double alpha);
		PolyDataPointerType selectSurfaceRoi(PolyDataPointerType surface, double * center, double radius);
		unsigned int vesselsAroundBarrel(double * barrelCentroid, std::list< unsigned int > vessels);
		void enforceAxisDivergence(std::map< int, double * > barrelAxes, std::map< int, double * > barrelCenters);
		void calculateAvgContours(std::map< int, PolyDataPointerType > barrels, std::map< int, double * > barrelAxes, std::map< int, double * > barrelCenters, std::map< int, PolyDataPointerType >& avgBarrels, std::map< int, std::vector< double * > >& endPoints);
		double * calculateBarrelCentroid(PolyDataPointerType barrel);
		double * calculateBarrelRadius(PolyDataPointerType barrel, double * barrelCenter);
		void transformBarrelZAxis(double * newAxis, double * barrelCentroid, PolyDataPointerType barrel);
		PolyDataPointerType smoothBarrelAlongNewAxis(double * newAxis, double * barrelCentroid, PolyDataPointerType barrel, std::vector< double * > endPoints);
		void computeAverageHomeBarrel(PolyDataPointerType completeBarrel, double barrelCentroid[3], double barrelAxis[3]);
		void closeBarrelAlongNewAxis(double * newAxis, double * barrelCentroid, PolyDataPointerType barrel, std::vector< double * >& endPoints);
		std::vector< double > computeBarrelParameters(PolyDataPointerType barrel, PolyDataPointerType pia, PolyDataPointerType wm, double * newAxis, double * barrelCenter, std::vector< double * > endPoints, int label, std::map< int, Column * >& barrelColumns, std::map< int, PolyDataPointerType >& avgBarrels);
		std::vector< double > computeManualBarrelParameters(PolyDataPointerType barrel, PolyDataPointerType pia, PolyDataPointerType wm, double * newAxis, double * barrelCenter, std::vector< double * > endPoints);
		PolyDataPointerType maxBarrelContour(PolyDataPointerType barrel, double * newAxis, double * barrelCenter, std::vector< double * > endPoints);
		double * axisSurfaceIntersection(PolyDataPointerType surface, double * axis, double * center);
		void writeBestBarrelAxes(std::multimap< double, double * > axes, double * barrelCentroid);
		void normalize(double * vec);
		
		int intersectConvexCellsInPlane(CellPointerType cell1, CellPointerType cell2, double tol, double p0[3], double p1[3]);
		
		Column * createBarrelColumn(PolyDataPointerType avgBarrel, double * top, double * bottom);
		std::map< int, std::vector< double > > calculateColumnOverlap(std::map< int, Column * > barrelColumns, int label);
		std::map< int, std::vector< double > > calculateSeptalDistances(std::map< int, Column * > barrelColumns, int label);
		ImageDataPointerType createColumnVoxels(Column * column, std::vector< double * >& contour, double& minDist, double& maxDist);
		std::map< int, std::list< int > > computeBarrelVesselCorrelations(std::map< int, PolyDataPointerType > avgBarrels, std::map< int, double* > barrelCenters, std::vector< Edge * > vesselVec);
		
		void septalDistances(std::map< int, PolyDataPointerType > avgBarrels, std::map< int, double* > barrelCenters);
		
		std::map< int, std::list< int > > createBarrelGrid(std::map< int, double * > barrelAxes);
		std::map< int, std::list< int > > createBarrelGrid(std::map< int, Column * > barrelColumns);
		std::map< int, std::list< int > > createNearestNeighborBarrelGrid(std::map< int, double * > barrelCenters);
		
		ImageDataPointerType addTop(int label, int interpolationSlices);
		ImageDataPointerType addTop2(int label, int additionalSections);
		double * samplingRayDerivatives(PointsPointerType);
		PolyDataPointerType calculateTop(std::vector< double * > derivatives, std::vector< PointsPointerType > existingPoints);
		
		ImageDataPointerType createImageVolumeFromPolyData(PolyDataPointerType poly, int label, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax);
		ImageDataPointerType createImageVolume(int label, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax);
		PolyDataPointerType createPolyDataFromPlanePointList(std::list< std::list< double * > > planewisePointList);
		ImageDataPointerType mergeStructures(ImageDataPointerType structure1, ImageDataPointerType structure2);
};



#endif
