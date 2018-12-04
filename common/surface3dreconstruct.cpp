#include "surface3dreconstruct.h"

Surface3DReconstruct::Surface3DReconstruct ( InputParameters parameters )
{
	this->parameters = parameters;
	SBF = new BarrelField(false);
}

Surface3DReconstruct::~Surface3DReconstruct()
{
	if(SBF) delete SBF;
}

/****************************************************************************/
/* private methods                                                          */
/****************************************************************************/

ImageDataPointerType Surface3DReconstruct::createImageVolumeFromPolyData ( PolyDataPointerType poly, double bounds[6], int label, double zSpacing )
{
	if(poly->GetNumberOfCells())
	{
		ImageDataPointerType volume = ImageDataPointerType::New();
		double spacing[3];
		if(!zSpacing)
			spacing[0] = spacing[1] = spacing[2] = 1;
		else
			spacing[0] = spacing[1] = spacing[2] = zSpacing;
		volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
		
		int dims[6];
		calculateExtent(bounds, dims, spacing);
// 		std::flush(std::cout << "max extent of input: [" << xMin << "," << xMax << "], [" << yMin << "," << yMax << "],[" << zMin << "," << zMax << "]" << std::endl);
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
						/*|| tmpCoord[2] < bounds[4] || tmpCoord[2] > bounds[5]*/)
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

ImageDataPointerType Surface3DReconstruct::distanceTransform ( ImageDataPointerType volume )
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
}

void Surface3DReconstruct::calculateExtent ( double bounds[6], int extent[6], double spacing[3] )
{
	//make sure that maxCoordinates are inside an integer number of cells defined by spacing
	extent[0] = (bounds[0] - spacing[0])/spacing[0]/* - 0.5*/;
	extent[1] = (bounds[1] + spacing[0])/spacing[0]/* + 0.5*/;
	extent[2] = (bounds[2] - spacing[1])/spacing[1]/* - 0.5*/;
	extent[3] = (bounds[3] + spacing[1])/spacing[1]/* + 0.5*/;
	extent[4] = (bounds[4] - spacing[2])/spacing[2]/* - 0.5*/;
	extent[5] = (bounds[5] + spacing[2])/spacing[2]/* + 0.5*/;
	// ugly fix in case minZ < 0... not sure if it's 100% bullet proof...
	if(bounds[4] < 0 && abs(bounds[4])%lround(spacing[2]))
		extent[4] -= 1;
}

NeighborhoodOffsetVectorType Surface3DReconstruct::CreateLookUpTable()
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
}

