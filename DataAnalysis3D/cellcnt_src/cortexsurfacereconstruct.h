/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      cortexsurfacereconstruct.h                                    */
/*                                                                          */
/* Purpose:   Class for 3D reconstruction of cortical surfaces from 2D      */
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

#ifndef CORTEXSURFACERECONSTRUCT_H
#define CORTEXSURFACERECONSTRUCT_H

#include "../../common/surface3dreconstruct.h"


class CortexSurfaceReconstruct : public Surface3DReconstruct
{
	public:
		CortexSurfaceReconstruct ( InputParameters parameters, PolyDataPointerType structure ) :
		Surface3DReconstruct(parameters) {this->structure = structure;}
		
		// implementation of the surface reconstruction
		virtual PolyDataPointerType surfaceReconstruction ( int label );

	protected:
		virtual PolyDataPointerType smoothSurface ( PolyDataPointerType surface );
		virtual ImageDataPointerType createImageVolumeFromPolyData ( PolyDataPointerType poly, double bounds[6], int label, double zSpacing = 0 );
		
	private:
		PolyDataPointerType structure;
		
		ImageDataPointerType addTop2(int label, int additionalSections, bool zReversed = 0, double zSpacing = 0);
		ImageDataPointerType piaVolume(int label, int additionalSections, bool zReversed = 0, double zSpacing = 0);
};

#endif // CORTEXSURFACERECONSTRUCT_H
