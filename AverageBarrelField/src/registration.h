/****************************************************************************/
/*                                                                          */
/* Program:   AverageBarrelField                                            */
/*                                                                          */
/* File:      registration.h                                                */
/*                                                                          */
/* Purpose:   class providing all methods for computing average             */
/*            barrel field                                                  */
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
#include "../../common/amiraReader.h"

#ifndef REGISTRATION
#define REGISTRATION

class Registration
{
	public:
		Registration(AmiraSpatialGraph * inputSpatialGraph);
		Registration();
		~Registration();
		
		std::list< AmiraSpatialGraph * > averageBarrelField(int nrOfBarrelFields, std::list< AmiraSpatialGraph * > allBarrelFields, std::vector< const char * > fileNames);
		
		AmiraSpatialGraph * getSpatialGraph() { return spatialGraph; }
		AmiraSpatialGraph * getSpatialGraph2() { return spatialGraph2; }
		
	private:
		bool zReversed;	// in the sense that the top usually has a lower z-value than the bottom (actually, z <- -z b/c it is really the depth below pia)
		int nrOfPointSets;
		std::list< int > barrelLabels;
		std::list< int > borderBarrels;
		std::map< int, double > avgTopDist;
		std::map< int, double > avgPiaWMDist;
		std::map< int, double > avgBarrelHeight;
		std::map< int, double > avgBarrelArea;
		std::list< double > cellTypeRatioDepths;
		std::list< double > cellTypeRatioDepthsSupra;
		std::list< double > cellTypeRatioDepthsGran;
		std::list< double > cellTypeRatioDepthsInfra;
		std::map< int, const char * > int2Labels;
		AmiraSpatialGraph * spatialGraph;
		AmiraSpatialGraph * spatialGraph2;
		std::map< int, Column * > avgColumns;
		std::map< int, Column * > avgBarrels;
		std::map< int, double * > avgAxes;
		std::map< int, double * > avgCenters;
		std::vector< TransformPointerType > transformVec;
		
		void initializeConstants();
		//helper methods
		void normalize(double * vec);
		double L2Distance3D(double x[3], double y[3]);
		double L2Distance2D(double x[2], double y[2]);
		std::map< int, std::list< int > > createBarrelGrid(std::map< int, double * > barrelAxes);
		HomogeneousMatrixPointerType getLocalBarrelCoordinates(double * newAxis);
		HomogeneousMatrixPointerType getLocalCoordinatesForVariation(double * newAxis);
		HomogeneousMatrixPointerType transformToBarrelCoordinates(double * newAxis);
		HomogeneousMatrixPointerType getNewXAxis(double angle);
		HomogeneousMatrixPointerType getNewZAxis(double angle);
		
		// removes centroids of barrel fields (optimal translation)
		void alignBarrelFieldCentroids(std::map< int, Column * > * barrels, int nrOfBarrelFields);
		// computes optimal rotation matrix matching two point sets
		// with known correspondences
		gsl_matrix * computeOptimalRotation(std::map< int, Column * > refBF, std::map< int, Column * > matchBF);
		// computes mean BT/BB points after one iteration of
		// finding optimal translations/rotations
		void computeMeanShape(std::map< int, Column * > * meanShape, std::map< int, Column * > * barrels);
		// computes residuals after one iteration of
		// finding optimal translations/rotations
		double getResiduals(std::map< int, Column * > * barrels);
		
		// creates standard barrel field from mean BT/BB points
		// using average values for laminar thickness, barrel area etc.
		void computeAverageBarrelField(std::map< int, Column * > * matchedBarrels, std::vector< const char * > fileNames);
		// transformation of standard barrel field to defined
		// coordinates (usually D2 center/axis aligned)
		void alignToReferenceFrame(std::map< int, Column * > * matchedBarrels);
		// make sure the columns are not converging towards the pia
		void enforceAxisDivergence(std::map< int, double * > barrelAxes, std::map< int, double * > barrelCenters);
		// quantification of the residual variability between corresponding
		// BT/BB points by computing eigenvalues of the covariance matrix
		// for each set of points
		void alignmentQuality(std::map< int, Column * > * barrels, std::vector< const char * > fileNames);
		void writeAvgBarrelData(std::map< int, Column * > * matchedBarrels, std::vector< const char * > fileNames);
		
		// create standardized surfaces by triangulation and interpolation
		// of the pia/WM/L4 intersection points along each column
		void computeAverageSurfaces(const char * surfaceFilename);
		// create standardized axis field by interpolation of the
		// standardized column axes. (refine to 50 micron voxels separately)
		void computeAverageAxesField();
		// helper methods for computation of the average surfaces/axes
		double * findExtremePoint(PointsPointerType pts, double center[3], double com[3]);
		double * findExtremeAxis(int thisID, double extrapolatedPt[3]);
};



#endif
