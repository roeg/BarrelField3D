/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      landmarkregistration.h                                        */
/*                                                                          */
/* Purpose:   Class for registration of two barrel fields. Implements       */
/*            method to obtain optimal transformation as described in the   */
/*            paper, but only for two barrel fields                         */
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

#ifndef LANDMARKREGISTRATION_H
#define LANDMARKREGISTRATION_H

class LandmarkRegistration
{
	public:
		LandmarkRegistration();
		LandmarkRegistration(std::map< int, Column * > inputBarrelField);
		LandmarkRegistration(AmiraSpatialGraph * inputSG);
		~LandmarkRegistration();
		
		// main method for computation of registration
		// using corresponding BT/BB points
		void startRegistration();
		// main method for computation of registration
		// using corresponding Column top/bottom points
		void startColumnRegistration();
		// method for registration including
		// anisotropic scaling
		void startRegAnisotropicScale();
		
		// returns optimal transformation
		TransformPointerType getTransform();
		void writeTransform ( const char * outputFilename );
		
	private:
		BarrelField * SBF;
		std::map< int, Column * > inputBarrelField;
		TransformPointerType regTransform;
		
		// methods implementing calculation of the optimal transformations
		void alignBarrelFieldCentroids(std::map< int, Column * > inputBarrels, std::map< int, Column * > refBarrels, double shift[3]);
		gsl_matrix * computeOptimalRotation(std::map< int, Column * > refBF, std::map< int, Column * > matchBF);
		void computeOptimalRotationScale(gsl_matrix * mX, gsl_matrix * mY, gsl_matrix * mU, gsl_matrix * mLambda);
		void computeOptimalScale(gsl_matrix * mX, gsl_matrix * mY, gsl_matrix * mU, gsl_matrix * mLambda, bool constrained=0);
		gsl_matrix * createPointPositionMatrix(std::map< int, Column * > barrels);
		// computes residuals after one iteration of
		// finding optimal translations/rotations
		double getResiduals(std::map< int, Column* > barrels1, std::map< int, Column* > barrels2);
		
		HomogeneousMatrixPointerType gsl2VtkMatrix(gsl_matrix * mIn);
};

#endif // LANDMARKREGISTRATION_H
