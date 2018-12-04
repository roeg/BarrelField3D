#include "typedefs.h"
#include "basics.h"
#include "barrel_field.h"
#include "inputparameters.h"

#ifndef SURFACE3DRECONSTRUCT_H
#define SURFACE3DRECONSTRUCT_H

// virtual base class for all surface reconstruction methods
// defines standard interface
class Surface3DReconstruct
{
	public:
		Surface3DReconstruct(InputParameters parameters);
		virtual ~Surface3DReconstruct();
		
		virtual PolyDataPointerType surfaceReconstruction(int label) = 0;
	
	protected:
		InputParameters parameters;
		BarrelField * SBF;
		
		// methods used instandard implementation
		// for isosurface creation from 2D contour stacks
		NeighborhoodOffsetVectorType CreateLookUpTable();
		virtual ImageDataPointerType createImageVolumeFromPolyData(PolyDataPointerType poly, double bounds[6], int label, double zSpacing = 0);
		virtual ImageDataPointerType distanceTransform(ImageDataPointerType volume);
		virtual void calculateExtent(double bounds[6], int extent[6], double spacing[3]);
		virtual PolyDataPointerType smoothSurface(PolyDataPointerType surface) = 0;
};

#endif // SURFACE3DRECONSTRUCT_H
