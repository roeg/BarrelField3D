/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      vpmreconstruction.h                                           */
/*                                                                          */
/* Purpose:   Class providing methods for 3D reconstruction of VPM          */
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
#include "barreloidsurfacereconstruct.h"

#ifndef VPMRECONSTRUCTION_H
#define VPMRECONSTRUCTION_H

class VPMReconstruction
{
	public:
		VPMReconstruction ( AmiraSpatialGraph * inputSG, InputParameters parameters );
		~VPMReconstruction();
		
		// main function starting the 3D reconstruction
		void startReconstruction ( const char * outputFilename );
		
		std::map< int, ClosedSurface * > getBarreloids();
		
	private:
		InputParameters parameters;
		AmiraSpatialGraph * spatialGraph;
		BarrelField * SBF;
		
		std::map< int, ClosedSurface * > barreloidField;
		// order: z, row, arc
		std::map< int, std::vector< double > > barreloidDimensions;
		std::map< int, std::vector< double > > barreloidCenters;
		std::map< int, std::vector< std::vector< double > > > barreloidEVecs;
		
		void writeBarreloidParameters ( const char * outputFilename );
		// computation of barreloid volume from closed surface
		// in the same fashion as Amira (add signed volumes of all
		// Tetraeders made up of surface triangles and origin)
		double barreloidVolume ( PolyDataPointerType barreloid );
		// computes dimension of the barreloids from ellipsoid
		// fitting results and barreloid surface
		void computeBarreloidDimensions();
};

#endif // VPMRECONSTRUCTION_H
