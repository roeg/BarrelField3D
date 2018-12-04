/****************************************************************************/
/*                                                                          */
/* Program:   NeuroRegistration                                             */
/*                                                                          */
/* File:      utilities.h                                                   */
/*                                                                          */
/* Purpose:   class providing methods for common computations and data      */
/*            structure handling                                            */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*            5353 Parkside Drive                                           */
/*            Jupiter, Florida 33458                                        */
/*            USA                                                           */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   10.02.2011                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#pragma once
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/barrel_field.h"
#include "../../common/amiraReader.h"

#ifndef UTILITIES
#define UTILITIES

// #define REG_ACCURACY

#ifdef REG_ACCURACY
//variables for parameter evaluation
extern double var_gamma;
extern double var_alpha;
#endif

class Utilities
{
	public:
		Utilities(AmiraSpatialGraph * inputSpatialGraph);
		Utilities();
		~Utilities();
		
		AmiraSpatialGraph * getSpatialGraph() { return spatialGraph; }
		
		// set flag explicitly; default 0
		void setL1flag(bool flag){L1flag = flag;}
		
	protected:
		AmiraSpatialGraph * spatialGraph;
		BarrelField * SBF;
		std::list< int > inputLabels;
		std::map< int, PolyDataPointerType > avgBarrelContours;
		// cell structure flags
		bool somaFlag, dendriteFlag, apicalFlag, basalFlag, axonFlag;
		// landmark flags
		bool piaFlag, wmFlag;
		// L1 neuron flag
		bool L1flag;
		// thickness of sections; should be detected automatically
		double piaSpacing;
		double wmSpacing;
		// z-direction of contours; automatic detection if possible
		bool zReversed;
		
// 		void initializeConstants();
// 		void readStandardBarrelField();
		void inputConsistencyCheck();
		
		void normalize(double * vec);
		double L2Distance3D(double x[3], double y[3]);
		double L2Distance2D(double x[2], double y[2]);
		
		void getPCenterOfStructure(AmiraSpatialGraph * sg, int ID, double centerPt[3]);
		
		std::map< int, std::list< int > > createBarrelGrid(std::map< int, double * > barrelAxes);
		HomogeneousMatrixPointerType getLocalBarrelCoordinates(double * newAxis);
		HomogeneousMatrixPointerType transformToBarrelCoordinates(double * newAxis);
		HomogeneousMatrixPointerType transformToBarrelCoordinates(double * oldAxis, double * newAxis);
		void calculateBarrelFieldCentroid(std::map< int, Column * > barrels, double centroid[3]);
		
		void detectSectionThickness();
		void detectInputZDirection();
		void alignSpatialGraphGlobalZ(std::list< double > zIndexList);
		void getLandmarkMinMaxIDs(PolyDataPointerType landmark, int& minID, int& maxID);
		
		PolyDataPointerType surfaceReconstruction(int label);
		double * newBarrelAxis(PolyDataPointerType barrel, PolyDataPointerType piaSurface, double alpha);
		std::multimap< double, double * > barrelAxisScores(PolyDataPointerType piaSurface, double * barrelCentroid, double radius, double alpha);
		PolyDataPointerType selectSurfaceRoi(PolyDataPointerType surface, double * center, double radius);
		double * calculateBarrelCentroid(PolyDataPointerType barrel);
		void closeBarrelAlongNewAxis(double * newAxis, double * barrelCentroid, PolyDataPointerType barrel, std::vector< double * >& endPoints, bool HBRecon);
		
		ImageDataPointerType piaVolume(int label, int additionalSections, bool zReversed = 0, double zSpacing = 0);
		ImageDataPointerType distanceTransform(ImageDataPointerType volume);
		int * calculateExtent(int minX, int maxX, int minY, int maxY, int minZ, int maxZ, int label);
		int * calculateExtent(int minX, int maxX, int minY, int maxY, int minZ, int maxZ, int label, double spacing[3]);
		NeighborhoodOffsetVectorType CreateLookUpTable();
		ImageDataPointerType addTop2(int label, int additionalSections, bool zReversed = 0, double zSpacing = 0);
		PolyDataPointerType smoothSurface(PolyDataPointerType surface);
		
		ImageDataPointerType createImageVolumeFromPolyData(PolyDataPointerType poly, int label, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax, double zSpacing = 0);
		ImageDataPointerType createImageVolume(int label, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax);
		PolyDataPointerType createPolyDataFromPlanePointList(std::list< std::list< double * > > planewisePointList);
		
		void alignBarrelFieldCentroids(std::map< int, Column * > * barrels, int nrOfBarrelFields);
		gsl_matrix * computeOptimalRotation(std::map< int, Column * > refBF, std::map< int, Column * > matchBF);
		gsl_matrix * computeOptimalRotation(std::map< int, double * > refBF, std::map< int, double * > matchBF);
		gsl_matrix * computeOptimalRotation2D(std::map< int, Column * > refBF, std::map< int, Column * > matchBF);
		HomogeneousMatrixPointerType gsl2VtkMatrix(gsl_matrix * mIn);
};

#endif