#include "cortexsurfacereconstruct.h"

PolyDataPointerType CortexSurfaceReconstruct::surfaceReconstruction ( int label )
{
	
	double spacing = 0;
	if(label == Pia)
		spacing = parameters.piaSpacing;
	else if(label == WhiteMatter)
		spacing = parameters.wmSpacing;
	ImageDataPointerType completeSurface = addTop2(label, 20, parameters.zReversed, spacing);
	MarchingCubesPointerType mcSurfaceFilter1 = MarchingCubesPointerType::New();
	mcSurfaceFilter1->SetInput(completeSurface);
	mcSurfaceFilter1->SetValue(0, 0);
	mcSurfaceFilter1->ComputeScalarsOff();
	mcSurfaceFilter1->ComputeGradientsOff();
	mcSurfaceFilter1->ComputeNormalsOff();
	mcSurfaceFilter1->Update();
	PolyDataPointerType smoothedSurface = smoothSurface(mcSurfaceFilter1->GetOutput());
	return smoothedSurface;
}

PolyDataPointerType CortexSurfaceReconstruct::smoothSurface ( PolyDataPointerType surface )
{
// 	AveragePolyDataFilterType smoothingFilter = AveragePolyDataFilterType::New();
// 	smoothingFilter->BoundarySmoothingOff();
// 	smoothingFilter->FeatureEdgeSmoothingOff();
// // 	smoothingFilter->FeatureEdgeSmoothingOn();
// // 	smoothingFilter->SetFeatureAngle(90);
// 	smoothingFilter->SetConvergence(0.001);
// 	smoothingFilter->SetNumberOfIterations(200);
// 	smoothingFilter->SetInput(surface);
// 	smoothingFilter->Update();
// 	return smoothingFilter->GetOutput();
	
	LowpassPolyDataFilterType smoothingFilter = LowpassPolyDataFilterType::New();
	smoothingFilter->BoundarySmoothingOff();
	smoothingFilter->FeatureEdgeSmoothingOff();
	smoothingFilter->NormalizeCoordinatesOn();
	smoothingFilter->SetNumberOfIterations(20);
	smoothingFilter->SetPassBand(0.01);
	smoothingFilter->SetInput(surface);
	smoothingFilter->Update();
	return smoothingFilter->GetOutput();
}

ImageDataPointerType CortexSurfaceReconstruct::createImageVolumeFromPolyData ( PolyDataPointerType poly, double bounds[6], int label, double zSpacing )
{
	if(poly->GetNumberOfCells())
	{
		ImageDataPointerType volume = ImageDataPointerType::New();
		double spacing[3];
		if(!zSpacing)
			switch(label)
			{
				case Pia:
					spacing[0] = spacing[1] = spacing[2] = 50;
					break;
					
				case WhiteMatter:
					spacing[0] = spacing[1] = spacing[2] = 50;
					break;
					
				default:
					spacing[0] = spacing[1] = spacing[2] = 1;
					break;
			}
		else
			spacing[0] = spacing[1] = spacing[2] = zSpacing;
		volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
		
		int dims[6];
		calculateExtent(bounds, dims, spacing);
// 		std::flush(std::cout << "max extent of input: [" << bounds[0] << "," << bounds[1] << "], [" << bounds[2] << "," << bounds[3] << "],[" << bounds[4] << "," << bounds[5] << "]" << std::endl);
// 		std::flush(std::cout << "Allocating memory for image  with dimensions [" << dims[0] << "," << dims[1] << "], [" << dims[2] << "," << dims[3] << "],[" << dims[4] << "," << dims[5] << "]" << std::endl);
		volume->SetExtent(dims);
		volume->SetNumberOfScalarComponents(1);
		volume->SetScalarTypeToUnsignedChar();
		
// 		poly->Print(std::cout);
// 		volume->Print(std::cout);
		
		volume->AllocateScalars();
		for(int z = dims[4]; z <= dims[5]; ++z)
			for(int y = dims[2]; y <= dims[3]; ++y)
				for(int x = dims[0]; x <= dims[1]; ++x)
				{
					unsigned char * px = static_cast< unsigned char * >(volume->GetScalarPointer(x, y, z));
					*px = 0;
				}
		volume->Update();
// 		std::flush(std::cout << "Calculating pixels inside polygons!" << std::endl);
		unsigned long insidePoints = 0;
		unsigned long outsidePoints = 0;
		for(int currPlane = 0; currPlane < poly->GetNumberOfCells(); ++currPlane)
		{
			double closestPoint[3];
			int subId;
			double pCoords[3];
			double * weights = new double[poly->GetCell(currPlane)->GetNumberOfPoints()];
			double dist2;
			double tmp[3];
			poly->GetCell(currPlane)->GetPoints()->GetPoint(0, tmp);
			double * cellBounds = poly->GetCell(currPlane)->GetBounds();
// 			std::cout << "point tmp @ [" << tmp[0] << "," << tmp[1] << "," << tmp[2] << "]" << std::endl;
			int z = lround(tmp[2]/spacing[2]);
// 			std::flush(std::cout << "Determining inside/outside polygon for all points plane z = " << z << std::endl);
			
			// Force all points of polygon to lie in the same plane
			// sometimes necessary b/c of some Amira imprecision artifacts...
			PolyDataPointerType thisPlanePoly = PolyDataPointerType::New();
			PolygonPointerType thisPolygon = PolygonPointerType::New();
			PointsPointerType thisPolyPoints = PointsPointerType::New();
			thisPlanePoly->Allocate(1);
			thisPolygon->GetPointIds()->SetNumberOfIds(poly->GetCell(currPlane)->GetNumberOfPoints());
			thisPolyPoints->SetDataTypeToFloat();
			thisPolyPoints->SetNumberOfPoints(poly->GetCell(currPlane)->GetNumberOfPoints());
			for(int ii = 0; ii < poly->GetCell(currPlane)->GetNumberOfPoints(); ++ii)
			{
				double tmp2[3];
				poly->GetCell(currPlane)->GetPoints()->GetPoint(ii, tmp2);
				tmp2[2] = 0;	// all in plane z = 0
				thisPolygon->GetPointIds()->SetId(ii, ii);
				thisPolyPoints->InsertPoint(ii, tmp2);
			}
			thisPlanePoly->InsertNextCell(thisPolygon->GetCellType(), thisPolygon->GetPointIds());
			thisPlanePoly->SetPoints(thisPolyPoints);
			thisPlanePoly->Update();
			
// 			std::flush(std::cout << "physical z = " << tmp[2] << std::endl);
// 			std::flush(std::cout << "z spacing = " << spacing[2] << std::endl);
// 			std::flush(std::cout << "cellBounds = [" << cellBounds[0] << "," << cellBounds[1] << "], [" << cellBounds[2] << "," << cellBounds[3] << "],[" << cellBounds[4] << "," << cellBounds[5] << "]" << std::endl);
// 			#pragma omp parallel for
			for(int y = dims[2]; y <= dims[3]; ++y)
			{
				for(int x = dims[0]; x <= dims[1]; ++x)
				{
					unsigned char * px = static_cast< unsigned char * >(volume->GetScalarPointer(x, y, z));
// 					double tmpCoord[] = {x*spacing[0], y*spacing[1], z*spacing[2]};
					// brute force method in case pia/WM have different spacing/offset
// 					double tmpCoord[] = {x*spacing[0], y*spacing[1], tmp[2]};
					double tmpCoord[] = {x*spacing[0], y*spacing[1], 0};	// check only in x-y direction!!!
					if(tmpCoord[0] < cellBounds[0] || tmpCoord[0] > cellBounds[1] 
						|| tmpCoord[1] < cellBounds[2] || tmpCoord[1] > cellBounds[3] 
						/*|| tmpCoord[2] < cellBounds[4] || tmpCoord[2] > cellBounds[5]*/)
					{
// 						std::flush(std::cout << "out of bounds." << std::endl);
						++outsidePoints;
						*px = 0;
						continue;
					}
// 					int insidePolygon = poly->GetCell(currPlane)->EvaluatePosition(tmpCoord, closestPoint, subId, pCoords, dist2, weights);
					int insidePolygon = thisPlanePoly->GetCell(0)->EvaluatePosition(tmpCoord, closestPoint, subId, pCoords, dist2, weights);
					if(insidePolygon == 1)
					{
// 						std::flush(std::cout << "hit!" << std::endl);
						*px = 255;
						++insidePoints;
					}
					else
					{
						++outsidePoints;
						*px = 0;
					}
				}
			}
// 			std::flush(std::cout << "outsidePoints = " << outsidePoints << " --- insidePoints = " << insidePoints << " --- volume: " << (dims[1]-dims[0]+1)*(dims[3]-dims[2]+1)*(dims[5]-dims[4]+1) << " points" << std::endl);
		}
// 		std::flush(std::cout << "outsidePoints = " << outsidePoints << " --- insidePoints = " << insidePoints << " --- volume: " << (dims[1]-dims[0]+1)*(dims[3]-dims[2]+1)*(dims[5]-dims[4]+1) << " points" << std::endl);
		volume->Update();
		return volume;
	}
	else
	{
		std::cout << "Error! Empty PolyData! Cannot convert to vtkImageData!" << std::endl;
		return NULL;
	}
}

ImageDataPointerType CortexSurfaceReconstruct::addTop2 ( int label, int additionalSections, bool zReversed, double zSpacing )
{
	ImageDataPointerType pia = piaVolume(label, additionalSections, zReversed, zSpacing);
	int extent[6] = {0, 0, 0, 0, 0, 0};	//minX, maxX, minY, maxY, minZ, maxZ
	pia->GetExtent(extent);
	int startZ;	//first section w/o data
	int stopZ = zReversed ? extent[4] : extent[5];
// 	std::cout << "stopZ = " << stopZ << std::endl;
	if(zReversed)
	{
		startZ = extent[4] + additionalSections;
// 		std::cout << "startZ = " << startZ << std::endl;
		if(startZ > extent[5] - 2)
		{
			std::cout << "Error! Pia volume incorrect: vertical extent smaller than data extent" << std::endl;
			return pia;
		}
	}
	else
	{
		startZ = extent[5] - additionalSections;
// 		std::cout << "startZ = " << startZ << std::endl;
		if(startZ < extent[4] + 2)
		{
			std::cout << "Error! Pia volume incorrect: vertical extent smaller than data extent" << std::endl;
			return pia;
		}
	}
	
	// for NeuroCount: zCount starts at 1
	int zCount = 1;
	int sign = zReversed ? -1 : 1;
	double curvature = zSpacing ? zSpacing/50.0 : 1;
	
	#ifdef REG_ACCURACY
	curvature *= var_gamma;
	#endif
	
// 	if(L1flag)
// 	{
// 		// flat pia top
// 		curvature *= 10000;
// 		zCount += 1;
// 		std::cout << ">>> Computing Pia top for L1 neurons..." << std::endl;
// 	}
	
// 	std::flush(std::cout << "sign = " << sign << std::endl);
// 	std::flush(std::cout << "curvature = " << curvature << std::endl);
// 	std::flush(std::cout << "startZ = " << startZ << std::endl);
// 	std::flush(std::cout << "stopZ = " << stopZ << std::endl);
	for(int z = startZ; z != stopZ; z += sign)
	{
		for(int y = extent[2]; y <= extent[3]; ++y)
			for(int x = extent[0]; x <= extent[1]; ++x)
			{
				float * dist0 = static_cast< float * >(pia->GetScalarPointer(x, y, z));
				float * dist1 = static_cast< float * >(pia->GetScalarPointer(x, y, z-sign*1));
				float * dist2 = static_cast< float * >(pia->GetScalarPointer(x, y, z-sign*2));
				*dist0 = *dist1 + (*dist1 - *dist2) + (curvature*zCount)*(curvature*zCount);	// effectively calculates gradient of the distance image plus an artificial curvature
			}
		++zCount;
	}
	
	return pia;
}

ImageDataPointerType CortexSurfaceReconstruct::piaVolume ( int label, int additionalSections, bool zReversed, double zSpacing )
{
// 		std::flush(std::cout << "Calculating isosurface..." << std::endl);
		double * bounds = structure->GetBounds();
		for(int ii = 0; ii < 3; ++ii)
		{
			bounds[2*ii] = round(bounds[2*ii]);
			bounds[2*ii+1] = round(bounds[2*ii+1]);
		}
		double zOffset = zSpacing ? zSpacing : 50;
		if(zReversed)
		{
			bounds[5] -= zOffset;	//open @ bottom for pia & WM
			bounds[4] -= zOffset*additionalSections;
		}
		else
		{
			bounds[5] += zOffset*additionalSections;
			bounds[4] += zOffset;	//open @ bottom for pia & WM
		}
		ImageDataPointerType volume = createImageVolumeFromPolyData(structure, bounds, label, zSpacing);
		ImageDataPointerType distVolume = distanceTransform(volume);
		distVolume->SetSpacing(volume->GetSpacing());
		return distVolume;
}
