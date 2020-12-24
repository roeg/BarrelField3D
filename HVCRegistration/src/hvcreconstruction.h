/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      barreloidsurfacereconstruct.h                                 */
/*                                                                          */
/* Purpose:   Class for 3D reconstruction and measurement of VPM barreloids */
/*            from aligned 2D contour stacks                                */
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

#ifndef HVCRECONSTRUCTION_H
#define HVCRECONSTRUCTION_H

#include "../../common/typedefs.h"
#include "../../common/surface3dreconstruct.h"

class HVCSurfaceReconstruct
{
	public:
		HVCSurfaceReconstruct ( InputParameters parameters, PolyDataPointerType HVCSurface );
		HVCSurfaceReconstruct ( InputParameters parameters, AmiraSpatialGraph * barreloidContours );
		// implementation of 3D surface reconstruction from 2D contours
		PolyDataPointerType surfaceReconstruction ( int label );
		
		// fit ellipsoid to reconstructed surface and return center
		std::vector< double > getEllipsoidCenter();
		// computes principal axis vectors of ellipsoid approximation of barreloid
		std::vector< std::vector< double > > getPrincipalAxes();
		// measure length of barreloid from fitted center along fitted axes
		std::vector< double > getHalfAxisLengths();
		
	private:
		bool ellipsoidDone;
		
		PolyDataPointerType HVCContours;
		PolyDataPointerType HVCSurface;
		UnstructuredGridPointerType HVCVolume;
		
		std::vector< double > ellipsoidCenter;
		std::vector< std::vector< double > > ellipsoidAxes;
		std::vector< double > ellipsoidAxisLengths;
		
		// implementation of ellipsoid fitting to surface points
		void computeEllipsoidFromContours();
		// use center of mass and principal axis decomposition
		void computePrincipalAxes();
		//helper methods
		void initializeData ( gsl_matrix * mData );
		void setEllipsoidMatrix ( gsl_matrix * mA, gsl_vector * vX );
		void initializeCovarianceMatrix(gsl_matrix* mCov, double mean[3]);
		
		// helper methods for surface triangulation
		std::vector< int > getOrderedContourIndices();
		UnstructuredGridPointerType triangulateAdjacentContours ( int index1, int index2 );
		UnstructuredGridPointerType triangulateEndContours ( int index, bool lastContour );

};

#endif // HVCRECONSTRUCTION_H
