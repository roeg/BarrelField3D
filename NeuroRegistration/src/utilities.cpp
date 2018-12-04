/****************************************************************************/
/*                                                                          */
/* Program:   NeuroRegistration                                             */
/*                                                                          */
/* File:      utilities.cpp                                                 */
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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "utilities.h"

// #define DEBUG

Utilities::Utilities(AmiraSpatialGraph * inputSpatialGraph)
{
	spatialGraph = inputSpatialGraph;
	SBF = new BarrelField;
	inputConsistencyCheck();
	L1flag = 0;
};

Utilities::Utilities()
{
	spatialGraph = NULL;
	SBF = new BarrelField;
	L1flag = 0;
};

Utilities::~Utilities()
{
	if(spatialGraph) delete spatialGraph;
	if(SBF) delete SBF;
};

ImageDataPointerType Utilities::piaVolume(int label, int additionalSections, bool zReversed, double zSpacing)
{
	std::list< std::list< double * > > planeEdgePointList;
	std::list< int > zIndexList;
	PolyDataPointerType polyData = PolyDataPointerType::New();
// 	if(spatialGraph->extractLandmark(label, planeEdgePointList, zIndexList))
	if(spatialGraph->extractLandmark(label, polyData))
	{
// 		std::flush(std::cout << "Calculating isosurface..." << std::endl);
		double * bounds = polyData->GetBounds();
		int xMin = round(bounds[0]), xMax = round(bounds[1]);
		int yMin = round(bounds[2]), yMax = round(bounds[3]);
		int zMin = round(bounds[4]), zMax = round(bounds[5]);
		double zOffset = zSpacing ? zSpacing : 50;
		if(zReversed)
		{
			zMax -= zOffset;	//open @ bottom for pia & WM
			zMin -= zOffset*additionalSections;
		}
		else
		{
			zMax += zOffset*additionalSections;
			zMin += zOffset;	//open @ bottom for pia & WM
		}
// 		PolyDataPointerType polyData = createPolyDataFromPlanePointList(planeEdgePointList);
		ImageDataPointerType volume = createImageVolumeFromPolyData(polyData, label, xMin, xMax, yMin, yMax, zMin, zMax, zSpacing);
		ImageDataPointerType distVolume = distanceTransform(volume);
		distVolume->SetSpacing(volume->GetSpacing());
		return distVolume;
	}
	
	else
	{
		std::cout << "Error! Empty SpatialGraph!" << std::endl;
		return NULL;
	}
};

int * Utilities::calculateExtent(int minX, int maxX, int minY, int maxY, int minZ, int maxZ, int label)
{
	int * extent = new int[6];
	double spacing[3];
	switch(label)
	{
		case Pia:
			spacing[0] = spacing[1] = spacing[2] = 50;
			break;
			
		case WhiteMatter:
			spacing[0] = spacing[1] = spacing[2] = 50;
			break;
			
		case Barrel:
			spacing[0] = spacing[1] = spacing[2] = 10;
			break;

		default:
			spacing[0] = spacing[1] = spacing[2] = 1;
			break;
	}
	
	//make sure that maxCoordinates are inside an integer number of cells defined by spacing
	extent[0] = ((double)minX - spacing[0])/spacing[0]/* - 0.5*/;
	extent[1] = ((double)maxX + spacing[0])/spacing[0]/* + 0.5*/;
	extent[2] = ((double)minY - spacing[1])/spacing[1]/* - 0.5*/;
	extent[3] = ((double)maxY + spacing[1])/spacing[1]/* + 0.5*/;
	extent[4] = ((double)minZ - spacing[2])/spacing[2]/* - 0.5*/;
	extent[5] = ((double)maxZ + spacing[2])/spacing[2]/* + 0.5*/;
	
	return extent;
};

int * Utilities::calculateExtent(int minX, int maxX, int minY, int maxY, int minZ, int maxZ, int label, double spacing[3])
{
	int * extent = new int[6];
	//make sure that maxCoordinates are inside an integer number of cells defined by spacing
	extent[0] = ((double)minX - spacing[0])/spacing[0]/* - 0.5*/;
	extent[1] = ((double)maxX + spacing[0])/spacing[0]/* + 0.5*/;
	extent[2] = ((double)minY - spacing[1])/spacing[1]/* - 0.5*/;
	extent[3] = ((double)maxY + spacing[1])/spacing[1]/* + 0.5*/;
	extent[4] = ((double)minZ - spacing[2])/spacing[2]/* - 0.5*/;
	extent[5] = ((double)maxZ + spacing[2])/spacing[2]/* + 0.5*/;
	// ugly fix in case minZ < 0... not sure if it's 100% bullet proof...
	if(minZ < 0 && abs(minZ)%lround(spacing[2]))
		extent[4] -= 1;
	
	return extent;
};

/******************************************************************************/
/*determine contour pixels for each plane, then compute plane-wise            */
/*distance-transform of image. Implemented in ITK.                            */
/******************************************************************************/
ImageDataPointerType Utilities::distanceTransform(ImageDataPointerType volume)
{
	VTK2ITKImageExportPointerType v2iExport = VTK2ITKImageExportPointerType::New();
	ITK2VTKImageImportPointerType i2vImport = ITK2VTKImageImportPointerType::New();
	VTK2ITKImageImportType::Pointer v2iImport = VTK2ITKImageImportType::New();
	ITK2VTKCalcImageExportType::Pointer i2vExport = ITK2VTKCalcImageExportType::New();
	
	ImageType::Pointer volumeImage = ImageType::New();
	ImageType::Pointer contourImage = ImageType::New();
	CalcImageType::Pointer distanceMap = CalcImageType::New();
	DistanceMapImageFilterType::Pointer distanceMapFilter = DistanceMapImageFilterType::New();
	
	v2iExport->SetInput(volume);
	v2iImport->SetCallbackUserData(v2iExport->GetCallbackUserData());
	v2iImport->SetDataExtentCallback(v2iExport->GetDataExtentCallback());
	v2iImport->SetNumberOfComponentsCallback(v2iExport->GetNumberOfComponentsCallback());
	v2iImport->SetOriginCallback(v2iExport->GetOriginCallback());
	v2iImport->SetPipelineModifiedCallback(v2iExport->GetPipelineModifiedCallback());
	v2iImport->SetPropagateUpdateExtentCallback(v2iExport->GetPropagateUpdateExtentCallback());
	v2iImport->SetScalarTypeCallback(v2iExport->GetScalarTypeCallback());
	v2iImport->SetSpacingCallback(v2iExport->GetSpacingCallback());
	v2iImport->SetUpdateDataCallback(v2iExport->GetUpdateDataCallback());
	v2iImport->SetUpdateInformationCallback(v2iExport->GetUpdateInformationCallback());
	v2iImport->SetWholeExtentCallback(v2iExport->GetWholeExtentCallback());
	v2iImport->SetBufferPointerCallback(v2iExport->GetBufferPointerCallback());
	
	volumeImage = v2iImport->GetOutput();
	volumeImage->Update();
	contourImage->SetRegions(volumeImage->GetLargestPossibleRegion());
	contourImage->Allocate();
	contourImage->FillBuffer(0);
	distanceMap->SetRegions(volumeImage->GetLargestPossibleRegion());
	distanceMap->Allocate();
	ImageType::SizeType radius;
	radius.Fill(1);
	SegNeighborhoodIteratorType volumeIter(radius, volumeImage, volumeImage->GetLargestPossibleRegion());
	SegNeighborhoodIteratorType contourIter(radius, contourImage, contourImage->GetLargestPossibleRegion());
	NeighborhoodOffsetVectorType offset = CreateLookUpTable();
	int neighborhood4[] = {10, 12, 13, 15};
	for(volumeIter.GoToBegin(), contourIter.GoToBegin(); !volumeIter.IsAtEnd() && !contourIter.IsAtEnd(); ++volumeIter, ++contourIter)
		if(volumeIter.GetCenterPixel())
		{
			contourIter.SetCenterPixel(volumeIter.GetCenterPixel());
			bool border[] = {0, 0, 0, 0};
			for(int ii = 0; ii < 4; ++ii)
				if(!volumeIter.GetPixel(offset[neighborhood4[ii]]))
					contourIter.SetPixel(offset[neighborhood4[ii]], volumeIter.GetCenterPixel());
		}
	contourImage->Update();
	
	int startZ = contourImage->GetLargestPossibleRegion().GetIndex()[2];
	int deltaZ = contourImage->GetLargestPossibleRegion().GetSize()[2];
	for(int z = startZ; z < (startZ + deltaZ); ++z)
	{
		ImageType::Pointer planeImage = ImageType::New();
		CalcImageType::Pointer planeDistanceMap = CalcImageType::New();
		ImageType::RegionType planeRegion;
		ImageType::RegionType planeImageRegion;
		ImageType::IndexType planeIndex;
		ImageType::SizeType planeSize;
		
		planeIndex = contourImage->GetLargestPossibleRegion().GetIndex();
		planeIndex[2] = z;
		planeSize = contourImage->GetLargestPossibleRegion().GetSize();
		planeSize[2] = 1;
		planeRegion.SetIndex(planeIndex);
		planeRegion.SetSize(planeSize);
		planeImageRegion.SetIndex(contourImage->GetLargestPossibleRegion().GetIndex());
		planeImageRegion.SetSize(planeSize);
		planeImage->SetRegions(planeImageRegion);
		planeImage->Allocate();
		
		ConstIteratorType copyIter(contourImage, planeRegion);
		IteratorType2 pasteIter(planeImage, planeImageRegion);
		for(copyIter.GoToBegin(), pasteIter.GoToBegin(); !copyIter.IsAtEnd() && !pasteIter.IsAtEnd(); ++copyIter, ++pasteIter)
			pasteIter.Set(copyIter.Get());
		planeImage->Update();
		
		distanceMapFilter->UseImageSpacingOff();
		distanceMapFilter->InsideIsPositiveOff();
		distanceMapFilter->SetInput(planeImage);
		planeDistanceMap = distanceMapFilter->GetDistanceMap();
		planeDistanceMap->Update();
		
		ConstCalcIteratorType copyIter2(planeDistanceMap, planeImageRegion);
		CalcIteratorType pasteIter2(distanceMap, planeRegion);
		for(copyIter2.GoToBegin(), pasteIter2.GoToBegin(); !copyIter2.IsAtEnd() && !pasteIter2.IsAtEnd(); ++copyIter2, ++pasteIter2)
			pasteIter2.Set(copyIter2.Get());
		distanceMap->Update();
	}
	
	i2vExport->SetInput(distanceMap);
	i2vImport->SetCallbackUserData(i2vExport->GetCallbackUserData());
	i2vImport->SetDataExtentCallback(i2vExport->GetDataExtentCallback());
	i2vImport->SetNumberOfComponentsCallback(i2vExport->GetNumberOfComponentsCallback());
	i2vImport->SetOriginCallback(i2vExport->GetOriginCallback());
	i2vImport->SetPipelineModifiedCallback(i2vExport->GetPipelineModifiedCallback());
	i2vImport->SetPropagateUpdateExtentCallback(i2vExport->GetPropagateUpdateExtentCallback());
	i2vImport->SetScalarTypeCallback(i2vExport->GetScalarTypeCallback());
	i2vImport->SetSpacingCallback(i2vExport->GetSpacingCallback());
	i2vImport->SetUpdateDataCallback(i2vExport->GetUpdateDataCallback());
	i2vImport->SetUpdateInformationCallback(i2vExport->GetUpdateInformationCallback());
	i2vImport->SetWholeExtentCallback(i2vExport->GetWholeExtentCallback());
	i2vImport->SetBufferPointerCallback(i2vExport->GetBufferPointerCallback());
	i2vImport->Update();
	
	ImageDataPointerType distanceVolume = ImageDataPointerType::New();
	distanceVolume->DeepCopy(i2vImport->GetOutput());
	distanceVolume->Update();
	return distanceVolume;
};

/******************************************************************************/
/*interpolate barrel/pia/WM top by following the gradient of the distance     */
/*transform of the existing landmark                                          */
/******************************************************************************/
ImageDataPointerType Utilities::addTop2(int label, int additionalSections, bool zReversed, double zSpacing)
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
	int zCount = 0;
	int sign = zReversed ? -1 : 1;
	double curvature = zSpacing ? zSpacing/50.0 : 1;
	
	#ifdef REG_ACCURACY
	curvature *= var_gamma;
	#endif
	
	if(L1flag)
	{
		// flat pia top
		curvature *= 10000;
		zCount += 1;
		std::cout << ">>> Computing Pia top for L1 neurons..." << std::endl;
	}
	
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
};

PolyDataPointerType Utilities::smoothSurface(PolyDataPointerType surface)
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
};

ImageDataPointerType Utilities::createImageVolumeFromPolyData(PolyDataPointerType poly, int label, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax, double zSpacing)
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
		int * dims = calculateExtent(xMin, xMax, yMin, yMax, zMin, zMax, label, spacing);
// 		std::flush(std::cout << "max extent of input: [" << xMin << "," << xMax << "], [" << yMin << "," << yMax << "],[" << zMin << "," << zMax << "]" << std::endl);
// 		std::flush(std::cout << "Allocating memory for image  with dimensions [" << dims[0] << "," << dims[1] << "], [" << dims[2] << "," << dims[3] << "],[" << dims[4] << "," << dims[5] << "]" << std::endl);
		volume->SetExtent(dims);
		volume->SetNumberOfScalarComponents(1);
		volume->SetScalarTypeToUnsignedChar();
		
		// try this instead of global shift:
// 		double zOffset = zMin - spacing[2];
		
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
		
		// try this instead of global shift:
// 		double origin[3];
// 		volume->GetOrigin(origin);
// 		std::cout << "origin @ [" << origin[0] << "," << origin[1] << "," << origin[2] << "]" << std::endl;
// 		origin[2] = zOffset;
// 		volume->SetOrigin(origin);
// 		volume->GetOrigin(origin);
// 		std::cout << "now origin @ [" << origin[0] << "," << origin[1] << "," << origin[2] << "]" << std::endl;
// 		volume->Update();
		
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
			double * bounds = poly->GetCell(currPlane)->GetBounds();
// 			std::cout << "point tmp @ [" << tmp[0] << "," << tmp[1] << "," << tmp[2] << "]" << std::endl;
			// try this instead of global shift:
// 			int z = lround((tmp[2] - zOffset)/spacing[2]);
			int z = lround(tmp[2]/spacing[2]);
			
// 			std::flush(std::cout << "Determining inside/outside polygon for all points plane z = " << z << std::endl);
// 			std::flush(std::cout << "physical z = " << tmp[2] << std::endl);
// 			std::flush(std::cout << "z spacing = " << spacing[2] << std::endl);
// 			std::flush(std::cout << "bounds = [" << bounds[0] << "," << bounds[1] << "], [" << bounds[2] << "," << bounds[3] << "],[" << bounds[4] << "," << bounds[5] << "]" << std::endl);
// 			#pragma omp parallel for
			for(int y = dims[2]; y <= dims[3]; ++y)
			{
				for(int x = dims[0]; x <= dims[1]; ++x)
				{
					unsigned char * px = static_cast< unsigned char * >(volume->GetScalarPointer(x, y, z));
// 					double tmpCoord[] = {x*spacing[0], y*spacing[1], z*spacing[2]};
					// brute force method in case pia/WM have different spacing/offset
					double tmpCoord[] = {x*spacing[0], y*spacing[1], tmp[2]};
					if(tmpCoord[0] < bounds[0] || tmpCoord[0] > bounds[1] 
						|| tmpCoord[1] < bounds[2] || tmpCoord[1] > bounds[3] 
						|| tmpCoord[2] < bounds[4] || tmpCoord[2] > bounds[5])
					{
// 						std::flush(std::cout << "out of bounds." << std::endl);
						++outsidePoints;
						*px = 0;
						continue;
					}
					int insidePolygon = poly->GetCell(currPlane)->EvaluatePosition(tmpCoord, closestPoint, subId, pCoords, dist2, weights);
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
};

ImageDataPointerType Utilities::createImageVolume(int label, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax)
{
	
	ImageDataPointerType volume = ImageDataPointerType::New();
	double spacing[3];
	if(Barrel <= label && label <= E6)
		label = Barrel;
	switch(label)
	{
		case Pia:
			spacing[0] = spacing[1] = spacing[2] = 50;
			break;
			
		case WhiteMatter:
			spacing[0] = spacing[1] = spacing[2] = 50;
			break;
			
		case Barrel:
			spacing[0] = spacing[1] = spacing[2] = 10;
			break;
			
		default:
			spacing[0] = spacing[1] = spacing[2] = 1;
			break;
	}
	
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	int * dims = calculateExtent(xMin, xMax, yMin, yMax, zMin, zMax, label);
	// 		std::flush(std::cout << "max extent of input: [" << xMin << "," << xMax << "], [" << yMin << "," << yMax << "],[" << zMin << "," << zMax << "]" << std::endl);
	// 		std::flush(std::cout << "Allocating memory for image  with dimensions [" << dims[0] << "," << dims[1] << "], [" << dims[2] << "," << dims[3] << "],[" << dims[4] << "," << dims[5] << "]" << std::endl);
	volume->SetExtent(dims);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToUnsignedChar();
	// 	volume->Print(std::cout);
	volume->AllocateScalars();
	for(int z = dims[4]; z <= dims[5]; ++z)
		for(int y = dims[2]; y <= dims[3]; ++y)
			for(int x = dims[0]; x <= dims[1]; ++x)
			{
				unsigned char * px = static_cast< unsigned char * >(volume->GetScalarPointer(x, y, z));
				*px = 0;
			}
	
	volume->Update();
	return volume;
};


/****************************************************************************/
/*moves all centroids to the origin                                         */
/****************************************************************************/
void Utilities::alignBarrelFieldCentroids(std::map< int, Column * > * barrels, int nrOfBarrelFields)
{
	for(int ii = 0; ii < nrOfBarrelFields; ++ii)
	{
		double shift[] = {0, 0, 0};
		std::map< int, Column * >::const_iterator bfIt;
		for(bfIt = barrels[ii].begin(); bfIt != barrels[ii].end(); ++bfIt)
		{
			for(int jj = 0; jj < 3; ++jj)
			{
				shift[jj] += bfIt->second->top[jj];
				shift[jj] += bfIt->second->bottom[jj];
			}
		}
		for(int jj = 0; jj < 3; ++jj)
		{
			if(barrels[ii].size())
				shift[jj] = -shift[jj]/(double)(2*barrels[ii].size());
		}
		for(bfIt = barrels[ii].begin(); bfIt != barrels[ii].end(); ++bfIt)
			bfIt->second->translateColumn(shift);
	}
};

/****************************************************************************/
/*compute optimal rotation matrix in the least squares sense                */
/*using the Kabsch algorithm:                                               */
/*minimize sum_i (U x_i - y_i)^2 where U = optimal rotation matrix,         */
/*x_i = ith landmark of matchBF, y_i = ith landmark of refBF                */
/****************************************************************************/
gsl_matrix * Utilities::computeOptimalRotation(std::map< int, Column * > refBF, std::map< int, Column * > matchBF)
{
	std::list< int > commonLandmarks;
	std::map< int, Column * >::const_iterator refBFIt, matchBFIt;
	for(refBFIt = refBF.begin(); refBFIt != refBF.end(); ++refBFIt)
		for(matchBFIt = matchBF.begin(); matchBFIt != matchBF.end(); ++matchBFIt)
			if(refBFIt->first == matchBFIt->first)
			{
				commonLandmarks.push_back(refBFIt->first);
				break;
			}
	
	gsl_matrix * mX = gsl_matrix_alloc(3, 2*commonLandmarks.size());
	gsl_matrix * mY = gsl_matrix_alloc(3, 2*commonLandmarks.size());
	gsl_matrix * mCov = gsl_matrix_alloc(3, 3);
	gsl_matrix * mLU = gsl_matrix_alloc(3, 3);
	gsl_matrix * mU = gsl_matrix_alloc(3, 3);
	gsl_matrix * mV = gsl_matrix_alloc(3, 3);
	gsl_matrix * mId = gsl_matrix_alloc(3, 3);
	gsl_vector * vS = gsl_vector_alloc(3);
	gsl_permutation * permLU = gsl_permutation_alloc(3);
	if(!mX || !mY || !mCov || !mLU || !mU || !mV || !mId || !vS || !permLU)
	{
		std::cout << "Error! Could not allocate enough memory for optimal rotation. Aborting..." << std::endl;
		return NULL;
	}
	
	std::list< int >::const_iterator commonLandmarksIt;
	int ii = 0;
	for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt, ++ii)
	{
		int ID = *commonLandmarksIt;
		double * refTop, * refBottom, * matchTop, * matchBottom;
		refTop = refBF[ID]->top, refBottom = refBF[ID]->bottom;
		matchTop = matchBF[ID]->top, matchBottom = matchBF[ID]->bottom;
		for(int jj = 0; jj < 3; ++jj)
		{
			double * mXTopPtr = gsl_matrix_ptr(mX, jj, 2*ii);
			double * mXBottomPtr = gsl_matrix_ptr(mX, jj, 2*ii+1);
			double * mYTopPtr = gsl_matrix_ptr(mY, jj, 2*ii);
			double * mYBottomPtr = gsl_matrix_ptr(mY, jj, 2*ii+1);
			if(!mXTopPtr || !mXBottomPtr || !mYTopPtr || !mYBottomPtr)
			{
				std::cout << "Error! Invalid memory access during landmark matrix allocation. Aborting..." << std::endl;
				return NULL;
			}
			*mXTopPtr = matchTop[jj], *mXBottomPtr = matchBottom[jj];
			*mYTopPtr = refTop[jj], *mYBottomPtr = refBottom[jj];
		}
	}
	int status;
	// compute X Y^t
	status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mX, mY, 0.0, mCov);
	
// 	std::cout << "mCov:" << status << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
// 	}
	// compute sign of determinant of mCov (right-/left-handed)
	// does not work for 3 points! rank of mCov = 2 in that case
	// -> det(mCov) = 0 -> not useful
// 	gsl_matrix_memcpy(mLU, mCov);
// 	int sign, detSign;
// 	status = gsl_linalg_LU_decomp(mLU, permLU, &sign);
// 	detSign = gsl_linalg_LU_sgndet(mLU, sign);
// 	
// 	std::cout << "mLU: " << status << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mLU, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mLU, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "detSign = " << detSign << std::endl;
	
	// compute SVD of mCov
	gsl_vector * vTmp = gsl_vector_alloc(3);
	status = gsl_linalg_SV_decomp(mCov, mV, vS, vTmp);
	
// 	std::cout << "mCov after SVD: " << status << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "mV after SVD: " << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mV, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mV, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "vS after SVD: " << std::endl;
// 	std::cout << "[" << gsl_vector_get(vS, 0) << "," << gsl_vector_get(vS, 1) << "," << gsl_vector_get(vS, 2) << "]" << std::endl;
	
	gsl_matrix * mUTmp, * mVTmp;
	gsl_permutation * permU, * permV;
	mUTmp = gsl_matrix_alloc(3, 3);
	mVTmp = gsl_matrix_alloc(3, 3);
	permU = gsl_permutation_alloc(3);
	permV = gsl_permutation_alloc(3);
	gsl_matrix_memcpy(mUTmp, mCov);
	gsl_matrix_memcpy(mVTmp, mV);
	int signumU, signumV, signU, signV, detSign;
	status = gsl_linalg_LU_decomp(mUTmp, permU, &signumU);
	signU = gsl_linalg_LU_sgndet(mUTmp, signumU);
	status = gsl_linalg_LU_decomp(mVTmp, permV, &signumV);
	signV = gsl_linalg_LU_sgndet(mVTmp, signumV);
	detSign = signU*signV;
	
// 	std::cout << "mUTmp: " << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mUTmp, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mUTmp, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "mVTmp: " << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mVTmp, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mVTmp, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "detSign = " << detSign << std::endl;
	
	// compute mU = mV * diag(1,1,detSign) * mCov^t
	// note: gsl already computes mV so that mCov = U * S * V^t
	gsl_matrix_set_identity(mId);
	gsl_matrix_set(mId, 2, 2, detSign);
	gsl_matrix * mTmp = gsl_matrix_alloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mId, mCov, 0.0, mTmp);
	status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mV, mTmp, 0.0, mU);
	
// 	std::cout << "mU: " << status << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mU, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mU, ii, 2) << "]" << std::endl;
// 	}
	
	gsl_matrix_free(mUTmp), gsl_matrix_free(mVTmp);
	gsl_permutation_free(permU), gsl_permutation_free(permV);
	gsl_matrix_free(mTmp), gsl_matrix_free(mCov), gsl_matrix_free(mId), gsl_matrix_free(mLU), gsl_matrix_free(mV), gsl_matrix_free(mX), gsl_matrix_free(mY);
	gsl_permutation_free(permLU);
	gsl_vector_free(vS), gsl_vector_free(vTmp);
	
	return mU;
};

/****************************************************************************/
/*compute optimal rotation matrix in the least squares sense                */
/*using the Kabsch algorithm:                                               */
/*minimize sum_i (U x_i - y_i)^2 where U = optimal rotation matrix,         */
/*x_i = ith landmark of matchBF, y_i = ith landmark of refBF                */
/*landmarks are only the barrel centers in this version!                    */
/****************************************************************************/
gsl_matrix * Utilities::computeOptimalRotation(std::map< int, double * > refBF, std::map< int, double * > matchBF)
{
	std::list< int > commonLandmarks;
	std::map< int, double * >::const_iterator refBFIt, matchBFIt;
	for(refBFIt = refBF.begin(); refBFIt != refBF.end(); ++refBFIt)
		for(matchBFIt = matchBF.begin(); matchBFIt != matchBF.end(); ++matchBFIt)
			if(refBFIt->first == matchBFIt->first)
			{
				commonLandmarks.push_back(refBFIt->first);
				break;
			}
	
	gsl_matrix * mX = gsl_matrix_alloc(3, commonLandmarks.size());
	gsl_matrix * mY = gsl_matrix_alloc(3, commonLandmarks.size());
	gsl_matrix * mCov = gsl_matrix_alloc(3, 3);
	gsl_matrix * mLU = gsl_matrix_alloc(3, 3);
	gsl_matrix * mU = gsl_matrix_alloc(3, 3);
	gsl_matrix * mV = gsl_matrix_alloc(3, 3);
	gsl_matrix * mId = gsl_matrix_alloc(3, 3);
	gsl_vector * vS = gsl_vector_alloc(3);
	gsl_permutation * permLU = gsl_permutation_alloc(3);
	if(!mX || !mY || !mCov || !mLU || !mU || !mV || !mId || !vS || !permLU)
	{
		std::cout << "Error! Could not allocate enough memory for optimal rotation. Aborting..." << std::endl;
		return NULL;
	}
	
	std::list< int >::const_iterator commonLandmarksIt;
	int ii = 0;
	for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt, ++ii)
	{
		int ID = *commonLandmarksIt;
		double * refCenter, * matchCenter;
		refCenter = refBF[ID];
		matchCenter = matchBF[ID];
		for(int jj = 0; jj < 3; ++jj)
		{
			double * mXPtr = gsl_matrix_ptr(mX, jj, ii);
			double * mYPtr = gsl_matrix_ptr(mY, jj, ii);
			if(!mXPtr || !mYPtr)
			{
				std::cout << "Error! Invalid memory access during landmark matrix allocation. Aborting..." << std::endl;
				return NULL;
			}
			*mXPtr = matchCenter[jj];
			*mYPtr = refCenter[jj];
		}
	}
	int status;
	// compute X Y^t
	status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mX, mY, 0.0, mCov);
	
// 	std::cout << "mCov:" << status << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
// 	}
	
	// compute sign of determinant of mCov (right-/left-handed)
	// does not work for 3 points! rank of mCov = 2 in that case
	// -> det(mCov) = 0 -> not useful
// 	gsl_matrix_memcpy(mLU, mCov);
// 	int sign, detSign;
// 	status = gsl_linalg_LU_decomp(mLU, permLU, &sign);
// 	detSign = gsl_linalg_LU_sgndet(mLU, sign);
// 	
// 	std::cout << "mLU: " << status << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mLU, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mLU, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "detSign = " << detSign << std::endl;
	
	// compute SVD of mCov
	gsl_vector * vTmp = gsl_vector_alloc(3);
	status = gsl_linalg_SV_decomp(mCov, mV, vS, vTmp);
	
// 	std::cout << "mCov after SVD: " << status << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "mV after SVD: " << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mV, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mV, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "vS after SVD: " << std::endl;
// 	std::cout << "[" << gsl_vector_get(vS, 0) << "," << gsl_vector_get(vS, 1) << "," << gsl_vector_get(vS, 2) << "]" << std::endl;
	
	gsl_matrix * mUTmp, * mVTmp;
	gsl_permutation * permU, * permV;
	mUTmp = gsl_matrix_alloc(3, 3);
	mVTmp = gsl_matrix_alloc(3, 3);
	permU = gsl_permutation_alloc(3);
	permV = gsl_permutation_alloc(3);
	gsl_matrix_memcpy(mUTmp, mCov);
	gsl_matrix_memcpy(mVTmp, mV);
	int signumU, signumV, signU, signV, detSign;
	status = gsl_linalg_LU_decomp(mUTmp, permU, &signumU);
	signU = gsl_linalg_LU_sgndet(mUTmp, signumU);
	status = gsl_linalg_LU_decomp(mVTmp, permV, &signumV);
	signV = gsl_linalg_LU_sgndet(mVTmp, signumV);
	detSign = signU*signV;
	
// 	std::cout << "mUTmp: " << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mUTmp, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mUTmp, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "mVTmp: " << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mVTmp, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mVTmp, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "detSign = " << detSign << std::endl;
	
	// compute mU = mV * diag(1,1,detSign) * mCov^t
	// note: gsl already computes mV so that mCov = U * S * V^t
	gsl_matrix_set_identity(mId);
	gsl_matrix_set(mId, 2, 2, detSign);
	gsl_matrix * mTmp = gsl_matrix_alloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mId, mCov, 0.0, mTmp);
	status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mV, mTmp, 0.0, mU);
	
// 	std::cout << "mU: " << status << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mU, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mU, ii, 2) << "]" << std::endl;
// 	}
	
	gsl_matrix_free(mUTmp), gsl_matrix_free(mVTmp);
	gsl_permutation_free(permU), gsl_permutation_free(permV);
	gsl_matrix_free(mTmp), gsl_matrix_free(mCov), gsl_matrix_free(mId), gsl_matrix_free(mLU), gsl_matrix_free(mV), gsl_matrix_free(mX), gsl_matrix_free(mY);
	gsl_permutation_free(permLU);
	gsl_vector_free(vS), gsl_vector_free(vTmp);
	
	return mU;
};

/****************************************************************************/
/*compute optimal rotation matrix in the least squares sense                */
/*using the Kabsch algorithm:                                               */
/*minimize sum_i (U x_i - y_i)^2 where U = optimal rotation matrix,         */
/*x_i = ith landmark of matchBF, y_i = ith landmark of refBF                */
/****************************************************************************/
gsl_matrix * Utilities::computeOptimalRotation2D(std::map< int, Column * > refBF, std::map< int, Column * > matchBF)
{
	std::list< int > commonLandmarks;
	std::map< int, Column * >::const_iterator refBFIt, matchBFIt;
	for(refBFIt = refBF.begin(); refBFIt != refBF.end(); ++refBFIt)
		for(matchBFIt = matchBF.begin(); matchBFIt != matchBF.end(); ++matchBFIt)
			if(refBFIt->first == matchBFIt->first)
			{
				commonLandmarks.push_back(refBFIt->first);
				break;
			}
	
	gsl_matrix * mX = gsl_matrix_alloc(2, 2*commonLandmarks.size());
	gsl_matrix * mY = gsl_matrix_alloc(2, 2*commonLandmarks.size());
	gsl_matrix * mCov = gsl_matrix_alloc(2, 2);
	gsl_matrix * mLU = gsl_matrix_alloc(2, 2);
	gsl_matrix * mU = gsl_matrix_alloc(2, 2);
	gsl_matrix * mV = gsl_matrix_alloc(2, 2);
	gsl_matrix * mId = gsl_matrix_alloc(2, 2);
	gsl_vector * vS = gsl_vector_alloc(2);
	gsl_permutation * permLU = gsl_permutation_alloc(2);
	if(!mX || !mY || !mCov || !mLU || !mU || !mV || !mId || !vS || !permLU)
	{
		std::cout << "Error! Could not allocate enough memory for optimal rotation. Aborting..." << std::endl;
		return NULL;
	}
	
	std::list< int >::const_iterator commonLandmarksIt;
	int ii = 0;
	for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt, ++ii)
	{
		int ID = *commonLandmarksIt;
		double * refTop, * refBottom, * matchTop, * matchBottom;
		refTop = refBF[ID]->top, refBottom = refBF[ID]->bottom;
		matchTop = matchBF[ID]->top, matchBottom = matchBF[ID]->bottom;
		for(int jj = 0; jj < 2; ++jj)
		{
			double * mXTopPtr = gsl_matrix_ptr(mX, jj, 2*ii);
			double * mXBottomPtr = gsl_matrix_ptr(mX, jj, 2*ii+1);
			double * mYTopPtr = gsl_matrix_ptr(mY, jj, 2*ii);
			double * mYBottomPtr = gsl_matrix_ptr(mY, jj, 2*ii+1);
			if(!mXTopPtr || !mXBottomPtr || !mYTopPtr || !mYBottomPtr)
			{
				std::cout << "Error! Invalid memory access during landmark matrix allocation. Aborting..." << std::endl;
				return NULL;
			}
			*mXTopPtr = matchTop[jj], *mXBottomPtr = matchBottom[jj];
			*mYTopPtr = refTop[jj], *mYBottomPtr = refBottom[jj];
		}
	}
	int status;
	// compute X Y^t
	status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mX, mY, 0.0, mCov);
// 	std::cout << "mCov:" << status << std::endl;
// 	for(int ii = 0; ii < 2; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 1; ++jj)
// 			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mCov, ii, 1) << "]" << std::endl;
// 	}
	// compute sign of determinant of mCov (right-/left-handed)
	gsl_matrix_memcpy(mLU, mCov);
	int sign, detSign;
	status = gsl_linalg_LU_decomp(mLU, permLU, &sign);
// 	std::cout << "mLU: " << status << std::endl;
// 	for(int ii = 0; ii < 2; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 1; ++jj)
// 			std::cout << gsl_matrix_get(mLU, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mLU, ii, 1) << "]" << std::endl;
// 	}
	detSign = gsl_linalg_LU_sgndet(mLU, sign);
// 	std::cout << "detSign = " << detSign << std::endl;
	// compute SVD of mCov
	gsl_vector * vTmp = gsl_vector_alloc(2);
	status = gsl_linalg_SV_decomp(mCov, mV, vS, vTmp);
	// compute mU = mV * diag(1,1,detSign) * mCov^t
	// note: gsl already computes mV so that mCov = U * S * V^t
// 	std::cout << "mCov after SVD: " << status << std::endl;
// 	for(int ii = 0; ii < 2; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 1; ++jj)
// 			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mCov, ii, 1) << "]" << std::endl;
// 	}
// 	std::cout << "mV after SVD: " << std::endl;
// 	for(int ii = 0; ii < 2; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 1; ++jj)
// 			std::cout << gsl_matrix_get(mV, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mV, ii, 1) << "]" << std::endl;
// 	}
	gsl_matrix_set_identity(mId);
	gsl_matrix_set(mId, 1, 1, detSign);
	gsl_matrix * mTmp = gsl_matrix_alloc(2, 2);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mId, mCov, 0.0, mTmp);
	status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mV, mTmp, 0.0, mU);
// 	std::cout << "mU: " << status << std::endl;
// 	for(int ii = 0; ii < 2; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 1; ++jj)
// 			std::cout << gsl_matrix_get(mU, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mU, ii, 1) << "]" << std::endl;
// 	}
	
	gsl_matrix_free(mTmp), gsl_matrix_free(mCov), gsl_matrix_free(mId), gsl_matrix_free(mLU), gsl_matrix_free(mV), gsl_matrix_free(mX), gsl_matrix_free(mY);
	gsl_permutation_free(permLU);
	gsl_vector_free(vS), gsl_vector_free(vTmp);
	
	return mU;
};

//deprecated; use AmiraSpatialGraph::extractLandmark() directly instead
PolyDataPointerType Utilities::createPolyDataFromPlanePointList(std::list< std::list< double * > > planewisePointList)
{
	if(planewisePointList.size())
	{
		PolyDataPointerType polyData = PolyDataPointerType::New();
		polyData->Allocate(1);
		PointsPointerType points = PointsPointerType::New();
		points->SetDataTypeToFloat();
		int lastID = 0;
		std::list< std::list< double * > >::iterator planeEdgePointListIt;
		for(planeEdgePointListIt = planewisePointList.begin(); planeEdgePointListIt != planewisePointList.end(); ++planeEdgePointListIt)
		{
			// 			std::flush(std::cout << "Allocating memory for " << planeEdgePointListIt->size() - 1 << " points in plane " << planeEdgePointListIt->back()[2] << "..." << std::endl);
			int end = planeEdgePointListIt->size() - 1;	// vtkPolygon does NOT use the same point twice on a contour as a SpatialGraph does
			PolygonPointerType poly = PolygonPointerType::New();
			poly->GetPointIds()->SetNumberOfIds(end);
			
			std::list< double * >::iterator planeVertexIt = planeEdgePointListIt->begin();
			for(int ii = 0; ii < end; ++planeVertexIt, ++ii)
			{
				double * tmp = *planeVertexIt;
				points->InsertNextPoint(tmp);
				poly->GetPointIds()->SetId(ii, ii + lastID);
			}
			lastID += end;
			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
		}
		polyData->SetPoints(points);
		polyData->Update();
		return polyData;
	}
	else
	{
		std::cout << "Error! Empty PointList! Cannot create vtkPolyData!" << std::endl;
		return NULL;
	}
};

NeighborhoodOffsetVectorType Utilities::CreateLookUpTable()
{
	SegNeighborhoodIteratorType::OffsetType offset;
	NeighborhoodOffsetVectorType look_up_table(0);
	
	for(int z = -1; z <= 1 ; z++)
	{
		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				if( (x == 0) && (y == 0) && (z == 0) )
				{}
				else
				{
					offset[0]=x;
					offset[1]=y;
					offset[2]=z;
					look_up_table.push_back(offset);
				}
			}
		}
	}
	
	return look_up_table;
};

/******************************************************************************/
/*returns matrix that gives coordinates of vectors from data coordinate       */
/*system in local barrel coordinate system determined by newAxis              */
/******************************************************************************/
HomogeneousMatrixPointerType Utilities::getLocalBarrelCoordinates(double * newAxis)
{
	TransformPointerType rotate = TransformPointerType::New();
	//old axis == z axis
	//rotation axis = newAxis x oldAxis
	//then, angle of rotation = acos(oldAxis*newAxis)
	double rotationAxis[3];
	rotationAxis[0] = -newAxis[1];
	rotationAxis[1] = newAxis[0];
	rotationAxis[2] = 0;
	double angle = acos(-newAxis[2])*180/PI;
	
	rotate->RotateWXYZ(angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);
// 	rotate->Print(std::cout);
	
	return rotate->GetMatrix();
};

/******************************************************************************/
/*returns matrix that rotates vector into local barrel coordinate system      */
/*determined by newAxis                                                       */
/******************************************************************************/
HomogeneousMatrixPointerType Utilities::transformToBarrelCoordinates(double * newAxis)
{
	TransformPointerType rotate = TransformPointerType::New();
	//old axis == z axis
	//rotation axis = newAxis x oldAxis
	//then, angle of rotation = acos(oldAxis*newAxis)
	double rotationAxis[3];
	rotationAxis[0] = newAxis[1];
	rotationAxis[1] = -newAxis[0];
	rotationAxis[2] = 0;
	double angle = acos(newAxis[2])*180/PI;
	
	rotate->RotateWXYZ(-angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);
// 	rotate->Print(std::cout);
	
	return rotate->GetMatrix();
};

/******************************************************************************/
/*returns matrix that rotates vector into local barrel coordinate system      */
/*determined by newAxis                                                       */
/******************************************************************************/
HomogeneousMatrixPointerType Utilities::transformToBarrelCoordinates(double * oldAxis, double * newAxis)
{
	TransformPointerType rotate = TransformPointerType::New();
	//old axis == z axis
	//rotation axis = oldAxis x newAxis
	//angle of rotation = acos(oldAxis*newAxis)
	double rotationAxis[3];
	rotationAxis[0] = oldAxis[1]*newAxis[2] - oldAxis[2]*newAxis[1];
	rotationAxis[1] = oldAxis[2]*newAxis[0] - oldAxis[0]*newAxis[2];
	rotationAxis[2] = oldAxis[0]*newAxis[1] - oldAxis[1]*newAxis[0];
	double angle = 0;
	for(int ii = 0; ii < 3; ++ii)
		angle += oldAxis[ii]*newAxis[ii];
	angle = acos(angle)*180/PI;
	
	rotate->RotateWXYZ(angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);
// 	rotate->Print(std::cout);
	
	return rotate->GetMatrix();
};

void Utilities::calculateBarrelFieldCentroid(std::map< int, Column* > barrels, double centroid[3])
{
	std::map< int, Column * >::const_iterator barrelIt;
	centroid[0] = centroid[1] = centroid[2] = 0;
	for(barrelIt = barrels.begin(); barrelIt != barrels.end(); ++barrelIt)
	{
		for(int ii = 0; ii < 3; ++ii)
		{
			centroid[ii] += barrelIt->second->top[ii];
			centroid[ii] += barrelIt->second->bottom[ii];
		}
	}
	for(int ii = 0; ii < 3; ++ii)
		centroid[ii] /= double(2*barrels.size());
};

/******************************************************************************/
/*assumes either barrel or pia contours are present and all have the same     */
/*spacing. also assumes that no landmark has only a single contour            */
/******************************************************************************/
void Utilities::detectSectionThickness()
{
	if(spatialGraph)
	{
		PolyDataPointerType pia = PolyDataPointerType::New();
		if(spatialGraph->extractLandmark(Pia, pia))
		{
			std::list< double > spacingList;
			std::list< double > zIndexList;
			for(int ii = 0; ii < pia->GetNumberOfCells(); ++ii)
			{
				double currZIndex = pia->GetCell(ii)->GetBounds()[4];
// 				std::cout << "currZIndex = " << currZIndex << std::endl;
				zIndexList.push_back(round(currZIndex));
			}
			zIndexList.sort();
			std::list< double >::const_iterator zIndexListIt1 = zIndexList.begin();
			std::list< double >::const_iterator zIndexListIt2;
			++zIndexListIt1;
			for(zIndexListIt2 = zIndexList.begin() ; zIndexListIt1 != zIndexList.end(); ++zIndexListIt1, ++zIndexListIt2)
				spacingList.push_back(round(*zIndexListIt1 - *zIndexListIt2));
			spacingList.sort();
			spacingList.unique();
			if(spacingList.size() > 1)
			{
				std::cout << "Warning! Input SpatialGraph has non-uniform pia section thickness." << std::endl;
				std::cout << "Run 'check_hoc_file.py' on input .hoc file first." << std::endl;
				std::cout << "Output invalid!" << std::endl;
			}
			else
			{
				piaSpacing = spacingList.front();
// 				std::cout << "piaSpacing = " << piaSpacing << std::endl;
				alignSpatialGraphGlobalZ(zIndexList);
			}
		}
		
		PolyDataPointerType wm = PolyDataPointerType::New();
		if(spatialGraph->extractLandmark(WhiteMatter, wm))
		{
			std::list< double > spacingList;
			std::list< double > zIndexList;
			for(int ii = 0; ii < wm->GetNumberOfCells(); ++ii)
			{
				double currZIndex = wm->GetCell(ii)->GetBounds()[4];
// 				std::cout << "currZIndex = " << currZIndex << std::endl;
				zIndexList.push_back(round(currZIndex));
			}
			zIndexList.sort();
			std::list< double >::const_iterator zIndexListIt1 = zIndexList.begin();
			std::list< double >::const_iterator zIndexListIt2;
			++zIndexListIt1;
			for(zIndexListIt2 = zIndexList.begin() ; zIndexListIt1 != zIndexList.end(); ++zIndexListIt1, ++zIndexListIt2)
				spacingList.push_back(round(*zIndexListIt1 - *zIndexListIt2));
			spacingList.sort();
			spacingList.unique();
			if(spacingList.size() > 1)
			{
				std::cout << "Warning! Input SpatialGraph has non-uniform WM section thickness." << std::endl;
				std::cout << "Run 'check_hoc_file.py' on input .hoc file first." << std::endl;
				std::cout << "Output invalid!" << std::endl;
			}
			else
			{
				wmSpacing = spacingList.front();
// 				std::cout << "wmSpacing = " << wmSpacing << std::endl;
// 				alignSpatialGraphGlobalZ(zIndexList);
			}
		}
	}
	else
		std::cout << "Error! Input SpatialGraph is empty" << std::endl;
};

/******************************************************************************/
/*determines z direction from pia contours:                                   */
/* z(smallest contour) > z(largest contour) ? zReversed = 0 : zReversed = 1   */
/******************************************************************************/
void Utilities::detectInputZDirection()
{
	if(spatialGraph)
	{
		if(!spatialGraph->isLabelInSpatialGraph(Pia))
		{
			std::cout << "Pia not found in Input SpatialGraph. Could not determine z-direction." << std::endl;
			std::cout << "Setting zReversed = 0" << std::endl;
			zReversed = 0;
			return;
		}
		
		PolyDataPointerType landmark = PolyDataPointerType::New();
		if(spatialGraph->extractLandmark(Pia, landmark))
		{
			double minArea, maxArea = landmark->GetCell(0)->GetLength2();
			minArea = maxArea;
			int minID = 0, maxID = 0;
			for(int ii = 1; ii < landmark->GetNumberOfCells(); ++ii)
			{
// 				landmark->GetCell(ii)->Print(std::cout);
				double tmpArea = landmark->GetCell(ii)->GetLength2();
// 				std::cout << "tmpArea = " << tmpArea << std::endl;
				if(tmpArea > maxArea)
				{
					maxArea = tmpArea;
					maxID = ii;
				}
				if(tmpArea < minArea)
				{
					minArea = tmpArea;
					minID = ii;
				}
			}
			double zMinID = landmark->GetCell(minID)->GetBounds()[4], zMaxID = landmark->GetCell(maxID)->GetBounds()[4];
			zReversed = zMinID > zMaxID ? 0 : 1;
		}
	}
	else
		std::cout << "Error! Input SpatialGraph is empty" << std::endl;
};

/******************************************************************************/
/*align z values of SpatialGraph such that they are aligned with z spacing    */
/*to give integer dimensions for vtkImageData (i.e., shift z_min of Pia to 0) */
/******************************************************************************/
void Utilities::alignSpatialGraphGlobalZ(std::list< double > zIndexList)
{
	if(zIndexList.size() && spatialGraph)
	{
		zIndexList.sort();
		double minZ = zIndexList.front();
		double ** transform = new double *[4];
		for(int ii = 0; ii < 4; ++ii)
		{
			transform[ii] = new double [4];
			for(int jj = 0; jj < 4; ++jj)
			{
				if(ii != jj)
					transform[ii][jj] = 0;
				else
					transform[ii][jj] = 1;
			}
		}
		transform[2][3] = -minZ;
		spatialGraph->setTransformation(transform);
		spatialGraph->applyTransformation();
		for(int ii = 0; ii < 4; ++ii)
			delete [] transform[ii];
		delete [] transform;
	}
	else
		std::cout << "Error! Input SpatialGraph is empty" << std::endl;
};

/******************************************************************************/
/*align z values of SpatialGraph such that they are aligned with z spacing    */
/*to give integer dimensions for vtkImageData (i.e., shift z_min of Pia to 0) */
/******************************************************************************/
void Utilities::getLandmarkMinMaxIDs(PolyDataPointerType landmark, int& minID, int& maxID)
{
	double minZ, maxZ = landmark->GetCell(0)->GetBounds()[4];
	minZ = maxZ;
	minID = 0, maxID = 0;
	for(int ii = 1; ii < landmark->GetNumberOfCells(); ++ii)
	{
		double tmpZ = landmark->GetCell(ii)->GetBounds()[4];
		if(tmpZ > maxZ)
		{
			maxZ = tmpZ;
			maxID = ii;
		}
		if(tmpZ < minZ)
		{
			minZ = tmpZ;
			minID = ii;
		}
	}
#ifdef DEBUG
	std::cout << "minZ = " << minZ << " -> minID = " << minID << std::endl;
	std::cout << "maxZ = " << maxZ << " -> maxID = " << maxID << std::endl;
#endif
	if(zReversed)
	{
		int tmp = minID;
		minID = maxID;
		maxID = tmp;
#ifdef DEBUG
		std::cout << "zReversed = 1:" << std::endl;
		std::cout << "new minID = " << minID << std::endl;
		std::cout << "new maxID = " << maxID << std::endl;
#endif
	}
};

PolyDataPointerType Utilities::surfaceReconstruction(int label)
{
	double spacing = 0;
	if(label == Pia)
		spacing = piaSpacing;
	else if(label == WhiteMatter)
		spacing = wmSpacing;
	ImageDataPointerType completeSurface = addTop2(label, 20, zReversed, spacing);
	MarchingCubesPointerType mcSurfaceFilter1 = MarchingCubesPointerType::New();
	mcSurfaceFilter1->SetInput(completeSurface);
	mcSurfaceFilter1->SetValue(0, 0);
	mcSurfaceFilter1->ComputeScalarsOff();
	mcSurfaceFilter1->ComputeGradientsOff();
	mcSurfaceFilter1->ComputeNormalsOff();
	mcSurfaceFilter1->Update();
	PolyDataPointerType smoothedSurface = smoothSurface(mcSurfaceFilter1->GetOutput());
	return smoothedSurface;
};

/****************************************************************************/
/*computes new barrel axis from possible candidate axes in a x-y radius of  */
/*2000micron. rates candidate axes based on distance to Pia, angle to Pia   */
/*normal at intersection point                                              */
/****************************************************************************/
double * Utilities::newBarrelAxis(PolyDataPointerType barrel, PolyDataPointerType piaSurface, double alpha)
{
	double * barrelCOM = calculateBarrelCentroid(barrel);
	std::multimap< double, double * > scores = barrelAxisScores(piaSurface, barrelCOM, 2000, alpha);
	
	/*** BEGIN MAKE NICE PLOTS OF AXIS BEFORE/AFTER ***/
// 		std::cout << "Writing 5 highest scoring axes out of " << scores.size() << std::endl;
// 		std::multimap< double, double * >::const_reverse_iterator axesIt;
// 		int count = 0;
// 		int nrOfAxisPoints = 2;
// 		for(axesIt = scores.rbegin(); axesIt != scores.rend() && count < 5; ++axesIt, ++count)
// 		{
// 			double * endPoint = new double[3];
// 			double * bottomPoint = new double[3];
// 			for(int ii = 0; ii < 3; ++ii)
// 			{
// 				endPoint[ii] = barrelCOM[ii] + axesIt->second[ii];
// 				bottomPoint[ii] = barrelCOM[ii] - 3*axesIt->second[ii];
// 			}
// 			int label;
// 			if(count)
// 				label = WhiteMatter;
// 			else
// 				label = ZAxis;
// 			Vertex * newVert1 = new Vertex(endPoint, label);
// 			Vertex * newVert2 = new Vertex(barrelCOM, label);
// 			workingSG->addVertex(newVert1);
// 			workingSG->addVertex(newVert2);
// 			int connectionIndex[2];
// 			if(!workingSG->getNumberOfVertices())
// 			{
// 				connectionIndex[0] = 0;
// 				connectionIndex[1] = 1;
// 			}
// 			else
// 			{
// 				connectionIndex[0] = workingSG->getNumberOfVertices() - 2;
// 				connectionIndex[1] = workingSG->getNumberOfVertices() - 1;
// 			}
// // 			int noOfAxisPoints = 2;
// 			std::list< double * > axisCoords;
// 			axisCoords.push_back(endPoint);
// 			axisCoords.push_back(barrelCOM);
// 			Edge * newAxis = new Edge(connectionIndex, nrOfAxisPoints, label, axisCoords);
// 			workingSG->addEdge(newAxis);
// 		}
	/*** END MAKE NICE PLOTS OF AXIS BEFORE/AFTER ***/
	
	double * newAxis = new double[3];
	newAxis[0] = scores.rbegin()->second[0], newAxis[1] = scores.rbegin()->second[1], newAxis[2] = scores.rbegin()->second[2];
	return newAxis;
};

/******************************************************************************/
/*returns pairs of (double score, double * vec) where score is the score of   */
/*the potential barrel axis in direction of vec from barrelCentroid           */
/*to the Pia                                                                  */
/*radius defines the radius of the cylinder defining the ROI                  */
/*of Pia surface patches                                                      */
/******************************************************************************/
std::multimap< double, double * > Utilities::barrelAxisScores(PolyDataPointerType piaSurface, double * barrelCentroid, double radius, double alpha)
{
	if(barrelCentroid != NULL)
	{
		PolyDataPointerType piaROI = selectSurfaceRoi(piaSurface, barrelCentroid, radius);
		if(piaROI->GetNumberOfCells())
		{
			std::vector< double > dotP;
			std::vector< double > dist;
			std::vector< double * > axis;
			std::multimap< double, double * > score;
			for(int ii = 0; ii < piaROI->GetNumberOfCells(); ++ii)
			{
				CellPointerType currentCell = piaROI->GetCell(ii);
				if(currentCell->GetNumberOfPoints() != 3)
				{
					std::cout << "Warning! vtkExtractPolyDataGeometry ouput is not triangles!" << std::endl;
				}
				PointsPointerType tmpPoints = currentCell->GetPoints();
				double centerPoint[] = {0, 0, 0};
				double axisNorm = 0;
				double * tmpAxis = new double[3];
				int nrOfPoints = tmpPoints->GetNumberOfPoints();
				for(int jj = 0; jj < nrOfPoints; ++jj)
				{
					double * tmpCoords = new double[3];
					tmpPoints->GetPoint(jj, tmpCoords);
					for(int kk = 0; kk < 3; ++kk)
						centerPoint[kk] += tmpCoords[kk];
					delete [] tmpCoords;
				}
				for(int jj = 0; jj < 3; ++jj)
					centerPoint[jj] = centerPoint[jj]/(double)nrOfPoints;
				
				for(int jj = 0; jj < 3; ++jj)
				{
					tmpAxis[jj] = centerPoint[jj] - barrelCentroid[jj];
					axisNorm += tmpAxis[jj]*tmpAxis[jj];
				}
				axisNorm = sqrt(axisNorm);
				
				PointsPointerType cellPoints = currentCell->GetPoints();
				double p1[3], p2[3], p3[3], normal[3];
				cellPoints->GetPoint(0, p1), cellPoints->GetPoint(1, p2), cellPoints->GetPoint(2, p3);
				for(int jj = 0; jj < 3; ++jj)
				{
					p3[jj] -= p1[jj];
					p2[jj] -= p1[jj];
				}
				normal[0] = p2[1]*p3[2] - p2[2]*p3[1];
				normal[1] = p2[2]*p3[0] - p2[0]*p3[2];
				normal[2] = p2[0]*p3[1] - p2[1]*p3[0];
				
				double tmpDotP = 0, normN = 0;
				for(int jj = 0; jj < 3; ++jj)
				{
					tmpDotP += tmpAxis[jj]*normal[jj];
					normN += normal[jj]*normal[jj];
				}
				normN = sqrt(normN);
				tmpDotP /= normN;
				tmpDotP /= axisNorm;
				tmpDotP = std::abs(tmpDotP);
				dotP.push_back(tmpDotP);
				dist.push_back(axisNorm);
				axis.push_back(tmpAxis);
			}
			
			if(dotP.size() == dist.size() && dotP.size() == axis.size())
			{
				if(score.size())
					score.clear();
				double maxDist = 0;
				for(int ii = 0; ii < dist.size(); ++ii)
					if(dist[ii] > maxDist)
						maxDist = dist[ii];
				for(int ii = 0; ii < dist.size(); ++ii)
				{
					double tmpScore = alpha*dotP[ii] + (1 - alpha)*(1 - dist[ii]/maxDist);
					score.insert(std::pair< double, double * >(tmpScore, axis[ii]));
				}
				return score;
			}
			else
			{
				std::cout << "Error! dotP.size() != dist.size() . Cannot calculate best axis." << std::endl;
				std::multimap< double, double * > emptyScore;
				return emptyScore;
			}
		}
	}
	std::multimap< double, double * > emptyScore;
	return emptyScore;
};

PolyDataPointerType Utilities::selectSurfaceRoi(PolyDataPointerType surface, double * center, double radius)
{
	CylinderPointerType cylinder = CylinderPointerType::New();
	TransformPointerType cylinderTransform = TransformPointerType::New();
	TransformPointerType cylinderTranslate = TransformPointerType::New();
	ExtractPolyDataGeometryPointerType roiClipper = ExtractPolyDataGeometryPointerType::New();
	CellLocatorPointerType locator = CellLocatorPointerType::New();
	//Cylinder is centered at Center and axes of rotation is along the y-axis
	//The transformation transforms a point into the space of the implicit function (i.e., the model space).
	//=> rotation about x-axis by 90 degrees! has to be in homogeneous coordinates
	double offset[] = {0, 0, 0};
	offset[0] = -center[0];
	offset[1] = -center[1];
	offset[2] = -center[2];
	cylinderTranslate->Translate(offset[0], offset[1], offset[2]);
	cylinderTransform->RotateX(90);
	cylinderTransform->Concatenate(cylinderTranslate);
	
	cylinder->SetRadius(radius);
	cylinder->SetCenter(0, 0, 0);
	cylinder->SetTransform(cylinderTransform);
	
	roiClipper->SetImplicitFunction(cylinder);
	roiClipper->ExtractInsideOn();
	roiClipper->ExtractBoundaryCellsOn();
	roiClipper->SetInput(surface);
	roiClipper->Update();
	return roiClipper->GetOutput();
};

/******************************************************************************/
/*adds polygons on top and bottom to make barrel extend up to extreme points  */
/*of existing contours measured along the new z axis.                         */
/*Stores the extreme points in vector endpoints in the order top, bottom      */
/******************************************************************************/
void Utilities::closeBarrelAlongNewAxis(double* newAxis, double* barrelCentroid, PolyDataPointerType barrel, std::vector< double* >& endPoints, bool HBRecon)
{
	if(barrel->GetNumberOfCells())
	{
		PolyDataPointerType barrelTop = PolyDataPointerType::New();
		PolyDataPointerType barrelBottom = PolyDataPointerType::New();
		IdListPointerType barrelTopCellIds = IdListPointerType::New();
		IdListPointerType barrelBottomCellIds = IdListPointerType::New();
		barrelTop->Allocate();
		barrelBottom->Allocate();
		
		int minID, maxID;
		getLandmarkMinMaxIDs(barrel, minID, maxID);
		barrelTopCellIds->InsertId(0, maxID);
		barrelBottomCellIds->InsertId(0, minID);
		barrelTop->CopyCells(barrel, barrelTopCellIds);
		barrelBottom->CopyCells(barrel, barrelBottomCellIds);
		
		// case first reconstruction:
		// find extreme bottom/top point by transforming the coordinate system
		// for every point in lowest/highest polygon and looking for extreme z values
		double normA = 0;
		double alpha = 0;
		double newZAxis[3];
		newZAxis[0] = newAxis[0], newZAxis[1] = newAxis[1], newZAxis[2] = newAxis[2];
		normalize(newZAxis);
		alpha = acos(std::abs(newZAxis[2]));
		TransformPointerType translate = TransformPointerType::New();
		TransformPointerType rotate = TransformPointerType::New();
		TransformPointerType inverseTranslate = TransformPointerType::New();
		
		//old axis == z axis
		//rotation axis = newAxus x oldAxis
		//then, angle of rotation = acos(oldAxis*newAxis)
		double rotationAxis[3];
		rotationAxis[0] = -newZAxis[1];
		rotationAxis[1] = newZAxis[0];
		rotationAxis[2] = 0;
		double sign = zReversed ? -1 : 1;
// 		double rotationAxis[3];
// 		rotationAxis[0] = -sign*newZAxis[1];
// 		rotationAxis[1] = sign*newZAxis[0];
// 		rotationAxis[2] = 0;
		
// 		double angle = sign*alpha*180/PI;	// -1 b/c we want to rotate coordinate system, not points
		double angle = alpha*180/PI;	// -1 b/c we want to rotate coordinate system, not points
		
		translate->Translate(-barrelCentroid[0], -barrelCentroid[1], -barrelCentroid[2]);
		rotate->RotateWXYZ(-angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);
		inverseTranslate->Translate(barrelCentroid[0], barrelCentroid[1], barrelCentroid[2]);
		
		rotate->Concatenate(translate);
		inverseTranslate->Concatenate(rotate);
		inverseTranslate->Update();
		
		// case reconstruction of registered home barrel along standard axis:
		// simply project onto new axis and look at extreme z-values
		// (all handled by bool flag HBRecon)
		
		int extremeBottomID = 0, extremeTopID = 0;
		double bottomMin = sign*1E06, topMax = -sign*1E06;
		PointsPointerType topPts = barrelTop->GetCell(0)->GetPoints();
		PointsPointerType bottomPts = barrelBottom->GetCell(0)->GetPoints();
		for(int ii = 0; ii < topPts->GetNumberOfPoints(); ++ii)
		{
			double thisZ;
			double * tmpPt = new double[3];
			topPts->GetPoint(ii, tmpPt);
			if(HBRecon)
			{
				double linePt1[3], linePt2[3], projPt[3], t;
				for(int jj = 0; jj < 3; ++jj)
				{
					linePt1[jj] = barrelCentroid[jj] + 1000*newZAxis[jj];
					linePt2[jj] = barrelCentroid[jj] - 1000*newZAxis[jj];
				}
				vtkLine::DistanceToLine(tmpPt, linePt1, linePt2, t, projPt);
				thisZ = projPt[2];
			}
			else
			{
				double homPt[4], transPt[4];
				homPt[0] = tmpPt[0], homPt[1] = tmpPt[1], homPt[2] = tmpPt[2], homPt[3] = 1;
				inverseTranslate->MultiplyPoint(homPt, transPt);
				thisZ = transPt[2];
			}
			if(zReversed ? thisZ < topMax : thisZ > topMax)
			{
				topMax = thisZ;
				extremeTopID = ii;
			}
			delete [] tmpPt;
		}
		for(int ii = 0; ii < bottomPts->GetNumberOfPoints(); ++ii)
		{
			double thisZ;
			double * tmpPt = new double[3];
			bottomPts->GetPoint(ii, tmpPt);
			if(HBRecon)
			{
				double linePt1[3], linePt2[3], projPt[3], t;
				for(int jj = 0; jj < 3; ++jj)
				{
					linePt1[jj] = barrelCentroid[jj] + 1000*newZAxis[jj];
					linePt2[jj] = barrelCentroid[jj] - 1000*newZAxis[jj];
				}
				vtkLine::DistanceToLine(tmpPt, linePt1, linePt2, t, projPt);
				thisZ = projPt[2];
			}
			else
			{
				double homPt[4], transPt[4];
				homPt[0] = tmpPt[0], homPt[1] = tmpPt[1], homPt[2] = tmpPt[2], homPt[3] = 1;
				inverseTranslate->MultiplyPoint(homPt, transPt);
				thisZ = transPt[2];
			}
			if(zReversed ? thisZ > bottomMin : thisZ < bottomMin)
			{
				bottomMin = thisZ;
				extremeBottomID = ii;
			}
			delete [] tmpPt;
		}
		
		LinePointerType newZAxisLine = LinePointerType::New();
		double * extremeTopPt = new double[3];
		double * extremeBottomPt = new double[3];
		double * finalTopPt = new double[3];
		double * finalBottomPt = new double[3];
		topPts->GetPoint(extremeTopID, extremeTopPt);
		bottomPts->GetPoint(extremeBottomID, extremeBottomPt);
		double linePt1[3], linePt2[3];
		linePt1[0] = barrelCentroid[0] + 1000*newZAxis[0], linePt1[1] = barrelCentroid[1] + 1000*newZAxis[1], linePt1[2] = barrelCentroid[2] + 1000*newZAxis[2];
		linePt2[0] = barrelCentroid[0] - 1000*newZAxis[0], linePt2[1] = barrelCentroid[1] - 1000*newZAxis[1], linePt2[2] = barrelCentroid[2] - 1000*newZAxis[2];
		
		//check whether new axis actually passes through top contour!!!
		double * topBounds = barrelTop->GetCell(0)->GetBounds();
		double * bottomBounds = barrelBottom->GetCell(0)->GetBounds();
		double * axisTopCoord = new double[3];
		{axisTopCoord[0] = barrelCentroid[0] + std::abs(barrelCentroid[2]-topBounds[4])*newZAxis[0]/cos(alpha), axisTopCoord[1] = barrelCentroid[1] + std::abs(barrelCentroid[2]-topBounds[4])*newZAxis[1]/cos(alpha), axisTopCoord[2] = topBounds[4];}
		double * axisBottomCoord = new double[3];
		{axisBottomCoord[0] = barrelCentroid[0] - std::abs(bottomBounds[5] - barrelCentroid[2])*newZAxis[0]/cos(alpha), axisBottomCoord[1] = barrelCentroid[1] - std::abs(bottomBounds[5] - barrelCentroid[2])*newZAxis[1]/cos(alpha), axisBottomCoord[2] = bottomBounds[5];}
		int insideTop = 0, insideBottom = 0;
		double closestPoint1[3], closestPoint2[3];
		int subId1, subId2;
		double pCoords1[3], pCoords2[3], * weights1, * weights2, dist1, dist2, tmp1[3], tmp2[3];
		weights1 = new double[barrelTop->GetCell(0)->GetNumberOfPoints()];
		weights2 = new double[barrelBottom->GetCell(0)->GetNumberOfPoints()];
		if(axisTopCoord[0] >= topBounds[0] && axisTopCoord[0] <= topBounds[1] 
			&& axisTopCoord[1] >= topBounds[2] && axisTopCoord[1] <= topBounds[3] 
			&& axisTopCoord[2] >= topBounds[4] && axisTopCoord[2] <= topBounds[5])
		{
			insideTop = barrelTop->GetCell(0)->EvaluatePosition(axisTopCoord, closestPoint1, subId1, pCoords1, dist1, weights1);
		}
		if(axisBottomCoord[0] >= bottomBounds[0] && axisBottomCoord[0] <= bottomBounds[1] 
			&& axisBottomCoord[1] >= bottomBounds[2] && axisBottomCoord[1] <= bottomBounds[3] 
			&& axisBottomCoord[2] >= bottomBounds[4] && axisBottomCoord[2] <= bottomBounds[5])
		{
			insideBottom = barrelBottom->GetCell(0)->EvaluatePosition(axisBottomCoord, closestPoint2, subId2, pCoords2, dist2, weights2);
		}
		if(insideTop)
		{
			double t = 0;
			double topPtDist = newZAxisLine->DistanceToLine(extremeTopPt, linePt1, linePt2, t, finalTopPt);
			endPoints.push_back(finalTopPt);
		}
		else if(!insideTop)
		{
			double t = 0;
			double topPtDist = newZAxisLine->DistanceToLine(extremeTopPt, linePt1, linePt2, t, finalTopPt);
			double * finalTopPt2 = new double[3];
			finalTopPt2[0] = finalTopPt[0], finalTopPt2[1] = finalTopPt[1], finalTopPt2[2] = finalTopPt[2];
			endPoints.push_back(finalTopPt2);
			double * centerPoint = new double[3];
			int subID;
			double pCoords[3], * weights;
			weights = new double[barrelTop->GetCell(0)->GetNumberOfPoints()];
			barrelTop->GetCell(0)->GetParametricCenter(pCoords);
			barrelTop->GetCell(0)->EvaluateLocation(subID, pCoords, centerPoint, weights);
			PlanePointerType topPlane = PlanePointerType::New();
			topPlane->SetOrigin(finalTopPt);
			topPlane->SetNormal(newZAxis);
			double planeDist = std::abs(topPlane->DistanceToPlane(centerPoint));
			finalTopPt[0] = centerPoint[0], finalTopPt[1] = centerPoint[1], finalTopPt[2] = centerPoint[2];
			if(cos(alpha) && !HBRecon)
				finalTopPt[2] -= planeDist/cos(alpha);
		}
		if(insideBottom)
		{
			double t = 0;
			double bottomPtDist = newZAxisLine->DistanceToLine(extremeBottomPt, linePt1, linePt2, t, finalBottomPt);
			endPoints.push_back(finalBottomPt);
		}
		else if(!insideBottom)
		{
			double t = 0;
			double bottomPtDist = newZAxisLine->DistanceToLine(extremeBottomPt, linePt1, linePt2, t, finalBottomPt);
			double * finalBottomPt2 = new double[3];
			finalBottomPt2[0] = finalBottomPt[0], finalBottomPt2[1] = finalBottomPt[1], finalBottomPt2[2] = finalBottomPt[2];
			endPoints.push_back(finalBottomPt2);
			double * centerPoint = new double[3];
			int subID;
			double pCoords[3], * weights;
			weights = new double[barrelBottom->GetCell(0)->GetNumberOfPoints()];
			barrelBottom->GetCell(0)->GetParametricCenter(pCoords);
			barrelBottom->GetCell(0)->EvaluateLocation(subID, pCoords, centerPoint, weights);
			PlanePointerType bottomPlane = PlanePointerType::New();
			bottomPlane->SetOrigin(finalBottomPt);
			bottomPlane->SetNormal(newZAxis);
			double planeDist = std::abs(bottomPlane->DistanceToPlane(centerPoint));
			finalBottomPt[0] = centerPoint[0], finalBottomPt[1] = centerPoint[1], finalBottomPt[2] = centerPoint[2];
			if(cos(alpha) && !HBRecon)
				finalBottomPt[2] += planeDist/cos(alpha);
		}
		// start linear interpolation of contours up to final points
		int nrOfTopPolys = 0, nrOfBottomPolys = 0;
		if(HBRecon)
		{
			nrOfTopPolys = int(L2Distance3D(finalTopPt, extremeTopPt));
			nrOfBottomPolys = int(L2Distance3D(finalBottomPt, extremeBottomPt));
		}
		else
		{
			nrOfTopPolys = int(std::abs(finalTopPt[2] - extremeTopPt[2]));
			nrOfBottomPolys = int(std::abs(finalBottomPt[2] - extremeBottomPt[2]));
		}
		std::vector< double * > topContourSlopes;
		std::vector< double * > bottomContourSlopes;
		for(int ii = 0; ii < topPts->GetNumberOfPoints(); ++ii)
		{
			double * tmpPt = new double[3];
			topPts->GetPoint(ii, tmpPt);
			for(int jj = 0; jj < 3; ++jj)
				tmpPt[jj] = (finalTopPt[jj] - tmpPt[jj])/(double)(nrOfTopPolys);
			topContourSlopes.push_back(tmpPt);
		}
		for(int ii = 0; ii < bottomPts->GetNumberOfPoints(); ++ii)
		{
			double * tmpPt = new double[3];
			bottomPts->GetPoint(ii, tmpPt);
			for(int jj = 0; jj < 3; ++jj)
				tmpPt[jj] = (finalBottomPt[jj] - tmpPt[jj])/(double)(nrOfBottomPolys);
			bottomContourSlopes.push_back(tmpPt);
		}
		
		// first top
		PolyDataPointerType newBarrelTop = PolyDataPointerType::New();
		PointsPointerType newBarrelTopPoints = PointsPointerType::New();
		newBarrelTop->Allocate();
		newBarrelTopPoints->SetDataTypeToFloat();
		int lastID = 0;
		double startZ = 0;
		for(int ii = 1; ii < nrOfTopPolys; ++ii)	// leave out last polygon on purpose to make sure we do not have some self-intersecting mess
		{
			PolygonPointerType newPoly = PolygonPointerType::New();
			newPoly->GetPointIds()->SetNumberOfIds(topPts->GetNumberOfPoints());
			for(int jj = 0; jj < topPts->GetNumberOfPoints(); ++jj)
			{
				double * tmpPt = new double[3];
				topPts->GetPoint(jj, tmpPt);
				for(int kk = 0; kk < 3; ++kk)
					tmpPt[kk] += ii*topContourSlopes[jj][kk];
				if(!HBRecon)
				{
					if(ii == 1)
					{
						tmpPt[2] = round(tmpPt[2]);
						startZ = tmpPt[2];
					}
					else
						tmpPt[2] = startZ + sign*(ii - 1);
				}
				newBarrelTopPoints->InsertNextPoint(tmpPt);
				newPoly->GetPointIds()->SetId(jj, jj + lastID);
			}
			lastID += topPts->GetNumberOfPoints();
			newBarrelTop->InsertNextCell(newPoly->GetCellType(), newPoly->GetPointIds());
		}
		newBarrelTop->SetPoints(newBarrelTopPoints);
		newBarrelTop->Update();
		
		//now bottom
		PolyDataPointerType newBarrelBottom = PolyDataPointerType::New();
		PointsPointerType newBarrelBottomPoints = PointsPointerType::New();
		newBarrelBottom->Allocate();
		newBarrelBottomPoints->SetDataTypeToFloat();
		lastID = 0;
		for(int ii = 1; ii < nrOfBottomPolys; ++ii)	// leave out last polygon on purpose to make sure we do not have some self-intersecting mess
		{
			PolygonPointerType newPoly = PolygonPointerType::New();
			newPoly->GetPointIds()->SetNumberOfIds(bottomPts->GetNumberOfPoints());
			for(int jj = 0; jj < bottomPts->GetNumberOfPoints(); ++jj)
			{
				double * tmpPt = new double[3];
				bottomPts->GetPoint(jj, tmpPt);
				for(int kk = 0; kk < 3; ++kk)
					tmpPt[kk] += ii*bottomContourSlopes[jj][kk];
				if(!HBRecon)
				{
					if(ii == 1)
					{
						tmpPt[2] = round(tmpPt[2]);
						startZ = tmpPt[2];
					}
					else
						tmpPt[2] = startZ - sign*(ii - 1);
				}
				newBarrelBottomPoints->InsertNextPoint(tmpPt);
				newPoly->GetPointIds()->SetId(jj, jj + lastID);
			}
			lastID += bottomPts->GetNumberOfPoints();
			newBarrelBottom->InsertNextCell(newPoly->GetCellType(), newPoly->GetPointIds());
		}
		newBarrelBottom->SetPoints(newBarrelBottomPoints);
		newBarrelBottom->Update();
		AppendPolyDataPointerType appendToBarrel = AppendPolyDataPointerType::New();
		appendToBarrel->AddInput(barrel);
		appendToBarrel->AddInput(newBarrelTop);
		appendToBarrel->AddInput(newBarrelBottom);
		appendToBarrel->Update();
		barrel->DeepCopy(appendToBarrel->GetOutput());
		barrel->Update();
	}
	else
	{
		std::cout << "Error! PolyData barrel is empty! Could not calculate barrel caps." << std::endl;
		return;
	}
};

/******************************************************************************/
/*calculate barrel centroid as average of cell centers                        */
/*do not use contour points; they might be spaced very irregularly            */
/******************************************************************************/
double * Utilities::calculateBarrelCentroid(PolyDataPointerType barrel)
{
	double * centroid = new double[3];
	for(int ii = 0; ii < 3; ++ii)
		centroid[ii] = 0;
	int nrOfCells = barrel->GetNumberOfCells();
	for(int ii = 0; ii < nrOfCells; ++ii)
	{
		double centerPoint[3], paramCenter[3];
		//parametric center
		int subID;
		double pCoords[3], * weights;
		weights = new double[barrel->GetCell(ii)->GetNumberOfPoints()];
		barrel->GetCell(ii)->GetParametricCenter(pCoords);
		barrel->GetCell(ii)->EvaluateLocation(subID, pCoords, centerPoint, weights);
		for(int jj = 0; jj < 3; ++jj)
			centroid[jj] += centerPoint[jj];
		delete [] weights;
	}
	for(int jj = 0; jj < 3; ++jj)
		centroid[jj] = centroid[jj]/(double)nrOfCells;
	return centroid;
};

/******************************************************************************/
/*create grid of barrels actually present in data: each barrel label has a    */
/*LUT of neighbors associated with it (i.e., create a graph)                  */
/******************************************************************************/
std::map< int, std::list< int > > Utilities::createBarrelGrid(std::map< int, double * > barrelAxes)
{
	std::map< int, std::list< int > > grid;
	std::list< int > barrelIDs;
	std::map< int, double * >::const_iterator barrelAxesIt;
	for(barrelAxesIt = barrelAxes.begin(); barrelAxesIt != barrelAxes.end(); ++barrelAxesIt)
		barrelIDs.push_back(barrelAxesIt->first);
	barrelIDs.sort();
	
	std::list< int > alpha;
	alpha.push_back(A1), alpha.push_back(B1), alpha.push_back(Beta);
	std::list< int > beta;
	beta.push_back(Alpha), beta.push_back(B1), beta.push_back(C1), beta.push_back(Gamma);
	std::list< int > gamma;
	gamma.push_back(Beta), gamma.push_back(C1), gamma.push_back(D1), gamma.push_back(Delta);
	std::list< int > delta;
	delta.push_back(Gamma), delta.push_back(D1), delta.push_back(E1);
	std::list< int > a1;
	a1.push_back(Alpha), a1.push_back(B1), a1.push_back(B2), a1.push_back(A2);
	std::list< int > a2;
	a2.push_back(A1), a2.push_back(B1), a2.push_back(B2), a2.push_back(B3), a2.push_back(A3);
	std::list< int > a3;
	a3.push_back(A2), a3.push_back(B2), a3.push_back(B3), a3.push_back(B4), a3.push_back(A4);
	std::list< int > a4;
	a4.push_back(A3), a4.push_back(B3), a4.push_back(B4);
	std::list< int > b1;
	b1.push_back(Beta), b1.push_back(C1), b1.push_back(C2), b1.push_back(B2), b1.push_back(A2), b1.push_back(A1), b1.push_back(Alpha);
	std::list< int > b2;
	b2.push_back(B1), b2.push_back(C1), b2.push_back(C2), b2.push_back(C3), b2.push_back(B3), b2.push_back(A3), b2.push_back(A2), b2.push_back(A1);
	std::list< int > b3;
	b3.push_back(B2), b3.push_back(C2), b3.push_back(C3), b3.push_back(C4), b3.push_back(B4), b3.push_back(A4), b3.push_back(A3), b3.push_back(A2);
	std::list< int > b4;
	b4.push_back(B3), b4.push_back(C3), b4.push_back(C4), b4.push_back(C5), b4.push_back(A4), b4.push_back(A3);
	std::list< int > c1;
	c1.push_back(Gamma), c1.push_back(D1), c1.push_back(D2), c1.push_back(C2), c1.push_back(B2), c1.push_back(B1), c1.push_back(Beta);
	std::list< int > c2;
	c2.push_back(C1), c2.push_back(D1), c2.push_back(D2), c2.push_back(D3), c2.push_back(C3), c2.push_back(B3), c2.push_back(B2), c2.push_back(B1);
	std::list< int > c3;
	c3.push_back(C2), c3.push_back(D2), c3.push_back(D3), c3.push_back(D4), c3.push_back(C4), c3.push_back(B4), c3.push_back(B3), c3.push_back(B2);
	std::list< int > c4;
	c4.push_back(C3), c4.push_back(D3), c4.push_back(D4), c4.push_back(D5), c4.push_back(C5), c4.push_back(B4), c4.push_back(B3);
	std::list< int > c5;
	c5.push_back(C4), c5.push_back(D4), c5.push_back(D5), c5.push_back(D6), c5.push_back(C6), c5.push_back(B4);
	std::list< int > c6;
	c6.push_back(C5), c6.push_back(D5), c6.push_back(D6);
	std::list< int > d1;
	d1.push_back(Delta), d1.push_back(E1), d1.push_back(E2), d1.push_back(D2), d1.push_back(C2), d1.push_back(C1), d1.push_back(Gamma);
	std::list< int > d2;
	d2.push_back(D1), d2.push_back(E1), d2.push_back(E2), d2.push_back(E3), d2.push_back(D3), d2.push_back(C3), d2.push_back(C2), d2.push_back(C1);
	std::list< int > d3;
	d3.push_back(D2), d3.push_back(E2), d3.push_back(E3), d3.push_back(E4), d3.push_back(D4), d3.push_back(C4), d3.push_back(C3), d3.push_back(C2);
	std::list< int > d4;
	d4.push_back(D3), d4.push_back(E3), d4.push_back(E4), d4.push_back(E5), d4.push_back(D5), d4.push_back(C5), d4.push_back(C4), d4.push_back(C3);
	std::list< int > d5;
	d5.push_back(D4), d5.push_back(E4), d5.push_back(E5), d5.push_back(E6), d5.push_back(D6), d5.push_back(C6), d5.push_back(C5), d5.push_back(C4);
	std::list< int > d6;
	d6.push_back(D5), d6.push_back(E5), d6.push_back(E6), d6.push_back(C6), d6.push_back(C5);
	std::list< int > e1;
	e1.push_back(Delta), e1.push_back(E2), e1.push_back(D2), e1.push_back(D1);
	std::list< int > e2;
	e2.push_back(E1), e2.push_back(E3), e2.push_back(D3), e2.push_back(D2), e2.push_back(D1);
	std::list< int > e3;
	e3.push_back(E2), e3.push_back(E4), e3.push_back(D4), e3.push_back(D3), e3.push_back(D2);
	std::list< int > e4;
	e4.push_back(E3), e4.push_back(E5), e4.push_back(D5), e4.push_back(D4), e4.push_back(D3);
	std::list< int > e5;
	e5.push_back(E4), e5.push_back(E6), e5.push_back(D6), e5.push_back(D5), e5.push_back(D4);
	std::list< int > e6;
	e6.push_back(E5), e6.push_back(D6), e6.push_back(D5);
	
	grid.insert(std::pair< int, std::list< int > >(Alpha, alpha));
	grid.insert(std::pair< int, std::list< int > >(Beta, beta));
	grid.insert(std::pair< int, std::list< int > >(Gamma, gamma));
	grid.insert(std::pair< int, std::list< int > >(Delta, delta));
	grid.insert(std::pair< int, std::list< int > >(A1, a1));
	grid.insert(std::pair< int, std::list< int > >(A2, a2));
	grid.insert(std::pair< int, std::list< int > >(A3, a3));
	grid.insert(std::pair< int, std::list< int > >(A4, a4));
	grid.insert(std::pair< int, std::list< int > >(B1, b1));
	grid.insert(std::pair< int, std::list< int > >(B2, b2));
	grid.insert(std::pair< int, std::list< int > >(B3, b3));
	grid.insert(std::pair< int, std::list< int > >(B4, b4));
	grid.insert(std::pair< int, std::list< int > >(C1, c1));
	grid.insert(std::pair< int, std::list< int > >(C2, c2));
	grid.insert(std::pair< int, std::list< int > >(C3, c3));
	grid.insert(std::pair< int, std::list< int > >(C4, c4));
	grid.insert(std::pair< int, std::list< int > >(C5, c5));
	grid.insert(std::pair< int, std::list< int > >(C6, c6));
	grid.insert(std::pair< int, std::list< int > >(D1, d1));
	grid.insert(std::pair< int, std::list< int > >(D2, d2));
	grid.insert(std::pair< int, std::list< int > >(D3, d3));
	grid.insert(std::pair< int, std::list< int > >(D4, d4));
	grid.insert(std::pair< int, std::list< int > >(D5, d5));
	grid.insert(std::pair< int, std::list< int > >(D6, d6));
	grid.insert(std::pair< int, std::list< int > >(E1, e1));
	grid.insert(std::pair< int, std::list< int > >(E2, e2));
	grid.insert(std::pair< int, std::list< int > >(E3, e3));
	grid.insert(std::pair< int, std::list< int > >(E4, e4));
	grid.insert(std::pair< int, std::list< int > >(E5, e5));
	grid.insert(std::pair< int, std::list< int > >(E6, e6));
	
	std::map< int, std::list< int > >::iterator gridIt;
	std::list< int >::iterator neighborIt;
	for(gridIt = grid.begin(); gridIt != grid.end(); )
	{
		if(std::find(barrelIDs.begin(), barrelIDs.end(), gridIt->first) != barrelIDs.end())
		{
			for(neighborIt = gridIt->second.begin(); neighborIt != gridIt->second.end(); )
			{
				if(std::find(barrelIDs.begin(), barrelIDs.end(), *neighborIt) != barrelIDs.end())
					++neighborIt;
				else
					neighborIt = gridIt->second.erase(neighborIt);
			}
			++gridIt;
		}
		else
		{
			std::map< int, std::list< int > >::iterator eraseIt = gridIt;
			++gridIt;
			grid.erase(eraseIt);
		}
	}
	
	return grid;
};

// does what you think it does
void Utilities::normalize(double * vec)
{
	double norm = 0;
	for(int ii = 0; ii < 3; ++ii)
		norm += vec[ii]*vec[ii];
	norm = sqrt(norm);
	if(norm)
		for(int ii = 0; ii < 3; ++ii)
			vec[ii] = vec[ii]/norm;
};

double Utilities::L2Distance3D(double x[3], double y[3])
{
	return sqrt((x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]) + (x[2] - y[2])*(x[2] - y[2]));
};

double Utilities::L2Distance2D(double x[2], double y[2])
{
	return sqrt((x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]));
};

void Utilities::getPCenterOfStructure(AmiraSpatialGraph * sg, int ID, double centerPt[3])
{
	PolyDataPointerType structure = PolyDataPointerType::New();
	if(!sg->extractLandmark(ID, structure))
	{
		std::cout << "Error! Could not find structure with ID " << ID <<" in SpatialGraph!" << std::endl;
		return;
	}
	int subID;
	double pCoords[3], * weights;
	weights = new double[structure->GetCell(0)->GetNumberOfPoints()];
	structure->GetCell(0)->GetParametricCenter(pCoords);
	structure->GetCell(0)->EvaluateLocation(subID, pCoords, centerPt, weights);
	delete [] weights;
};

// turns 3x3 rotation matrix into 4x4 homogeneous VTK transformation matrix
// Warning: performs no checking whatsoever
HomogeneousMatrixPointerType Utilities::gsl2VtkMatrix(gsl_matrix * mIn)
{
	HomogeneousMatrixPointerType mOut = HomogeneousMatrixPointerType::New();
	mOut->Identity();
	for(int ii = 0; ii < 3; ++ii)
		for(int jj = 0; jj < 3; ++jj)
			mOut->SetElement(ii, jj, gsl_matrix_get(mIn, ii, jj));
	return mOut;
};

void Utilities::inputConsistencyCheck()
{
	somaFlag = dendriteFlag = apicalFlag = basalFlag = axonFlag = 0;
	piaFlag = wmFlag = 0;
	if(spatialGraph)
	{
		for(int ii = Neuron; ii <= ZAxis; ++ii)
			if(spatialGraph->isLabelInSpatialGraph(ii))
				inputLabels.push_back(ii);
		
		if(std::find(inputLabels.begin(), inputLabels.end(), Soma) != inputLabels.end())
		{
			std::cout << "Soma found in input hoc file!" << std::endl;
			somaFlag = 1;
		}
		if(std::find(inputLabels.begin(), inputLabels.end(), Dendrite) != inputLabels.end())
		{
			std::cout << "Dendrite found in input hoc file!" << std::endl;
			dendriteFlag = 1;
		}
		if(std::find(inputLabels.begin(), inputLabels.end(), ApicalDendrite) != inputLabels.end())
		{
			std::cout << "Apical dendrite found in input hoc file!" << std::endl;
			apicalFlag = 1;
		}
		if(std::find(inputLabels.begin(), inputLabels.end(), BasalDendrite) != inputLabels.end())
		{
			std::cout << "Basal dendrite found in input hoc file!" << std::endl;
			basalFlag = 1;
		}
		if(std::find(inputLabels.begin(), inputLabels.end(), Axon) != inputLabels.end())
		{
			std::cout << "Axon found in input hoc file!" << std::endl;
			axonFlag = 1;
		}
		if(std::find(inputLabels.begin(), inputLabels.end(), Pia) != inputLabels.end())
		{
			std::cout << "Pia found in input hoc file!" << std::endl;
			piaFlag = 1;
		}
		if(std::find(inputLabels.begin(), inputLabels.end(), WhiteMatter) != inputLabels.end())
		{
			std::cout << "WM found in input hoc file!" << std::endl;
			wmFlag = 1;
		}
		
		//round Pia and WM z values to integers!
		// (sometimes they are systematically off by some decimal places...)
		//makes everything A LOT easier and more robust
		
// 		std::vector< Vertex * >::iterator vertexIt;
// 		for(vertexIt = spatialGraph->verticesBegin(); vertexIt != spatialGraph->verticesEnd(); ++vertexIt)
// 			if((*vertexIt)->label > Landmark)
// 			{
// 				double tmpZ = (*vertexIt)->coordinates[2];
// 				tmpZ = round(tmpZ);
// 				if(tmpZ-(*vertexIt)->coordinates[2])
// 					(*vertexIt)->coordinates[2] = tmpZ;
// 				
// 			}
// 		std::vector< Edge * >::iterator edgeIt;
// 		for(edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt)
// 		{
// 			if((*edgeIt)->label > Landmark)
// 			{
// 				std::list< double * >::iterator edgePtIt;
// 				for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
// 				{
// 					double tmpZ = (*edgePtIt)[2];
// 					tmpZ = round(tmpZ);
// 					if(tmpZ-(*edgePtIt)[2])
// 						(*edgePtIt)[2] = tmpZ;
// 				}
// 				
// 			}
// 		}
	}
};
