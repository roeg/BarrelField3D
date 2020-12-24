#include "hvcreconstruction.h"

// #define DEBUG

HVCSurfaceReconstruct::HVCSurfaceReconstruct ( InputParameters parameters, PolyDataPointerType HVCSurface )
{
	this->HVCSurface = HVCSurface;
	ellipsoidDone = 0;
}

HVCSurfaceReconstruct::HVCSurfaceReconstruct(InputParameters parameters, AmiraSpatialGraph * barreloidContours)
{
	this->HVCContours = PolyDataPointerType::New();
	if(!barreloidContours->extractLandmark(D2, this->HVCContours))
	{
		std::flush(std::cout << "Error: Could not find contour landmarks (label D2) in SpatialGraph" << std::endl);
	}
// #ifdef DEBUG
// 	this->HVCContours->Print(std::cout);
// #endif
	ellipsoidDone = 0;
}

std::vector< double > HVCSurfaceReconstruct::getEllipsoidCenter()
{
	if(!ellipsoidDone)
	{
// 		computeEllipsoidFromContours();
		computePrincipalAxes();
	}
	
	return ellipsoidCenter;
}

std::vector< std::vector< double > > HVCSurfaceReconstruct::getPrincipalAxes()
{
	if(!ellipsoidDone)
	{
// 		computeEllipsoidFromContours();
		computePrincipalAxes();
	}
	
	return ellipsoidAxes;
}

std::vector< double > HVCSurfaceReconstruct::getHalfAxisLengths()
{
	if(!ellipsoidDone)
	{
// 		computeEllipsoidFromContours();
		computePrincipalAxes();
	}
	
	return ellipsoidAxisLengths;
}

// use 3D Delaunay and find
// outer surface of tetrahedrons
PolyDataPointerType HVCSurfaceReconstruct::surfaceReconstruction ( int label )
{
	AppendFilterPointerType mergeTetrasFilter = AppendFilterPointerType::New();
	mergeTetrasFilter->MergePointsOn();
	
	std::vector< int > contourIndices = getOrderedContourIndices();
	mergeTetrasFilter->AddInput(triangulateEndContours(contourIndices[0], 0));
	for(int ii = 0; ii < contourIndices.size()-1; ++ii)
// 	for(int ii = 0; ii < 1; ++ii)
	{
		mergeTetrasFilter->AddInput(triangulateAdjacentContours(contourIndices[ii], contourIndices[ii+1]));
	}
	mergeTetrasFilter->AddInput(triangulateEndContours(contourIndices.back(), 1));
	mergeTetrasFilter->Update();
	
	HVCVolume = mergeTetrasFilter->GetOutput();
#ifdef DEBUG
	HVCVolume->Print(std::cout);
#endif
	
// 	Delaunay3DFilterPointerType del3DTriangulation = Delaunay3DFilterPointerType::New();
// 	del3DTriangulation->SetInput(HVCContours);
// 	del3DTriangulation->Update();
// 	HVCVolume = del3DTriangulation->GetOutput();
	
// 	return smoothSurface(mergeSurfacesFilter->GetOutput());
// 	DataSetSurfaceFilterPointerType surfaceExtractor = DataSetSurfaceFilterPointerType::New();
// 	surfaceExtractor->SetInput(mergeTetrasFilter->GetOutput());
// 	surfaceExtractor->Update();
// 	#ifdef DEBUG
// 	surfaceExtractor->GetOutput()->Print(std::cout);
// 	#endif
	GeometryFilterPointerType surfaceExtractor = GeometryFilterPointerType::New();
	surfaceExtractor->SetInput(mergeTetrasFilter->GetOutput());
// 	surfaceExtractor->SetInput(del3DTriangulation->GetOutput());
	surfaceExtractor->MergingOn();
	surfaceExtractor->Update();
	#ifdef DEBUG
	surfaceExtractor->GetOutput()->Print(std::cout);
	#endif
	
	PolyDataPointerType hvcSurface_ = PolyDataPointerType::New();
	hvcSurface_ = surfaceExtractor->GetOutput();
	double contourCorrectionShift[] = {0, 0, -50};
	TransformPointerType contourCorrectionTransform = TransformPointerType::New();
	contourCorrectionTransform->Translate(contourCorrectionShift);
	TransformFilterType contourCorrectionFilter = TransformFilterType::New();
	contourCorrectionFilter->SetTransform(contourCorrectionTransform);
	contourCorrectionFilter->SetInput(hvcSurface_);
	contourCorrectionFilter->Update();
	HVCSurface = contourCorrectionFilter->GetOutput();
	return HVCSurface;
// 	return smoothSurface(surfaceExtractor->GetOutput());
}

void HVCSurfaceReconstruct::computePrincipalAxes()
{
	if(HVCSurface)
	{
		#ifdef DEBUG
		std::cout << "Fitting principal axes to " << HVCSurface->GetNumberOfCells() << " points" << std::endl;
		#endif
		int status;
		gsl_matrix * mCov = gsl_matrix_alloc(3, 3);
		double center[3];
		initializeCovarianceMatrix(mCov, center);
		
		// diagonalization
		gsl_vector * eigenVals = gsl_vector_alloc(3);
		gsl_matrix * eigenVecs = gsl_matrix_alloc(3, 3);
		gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(3);
		status = gsl_eigen_symmv(mCov, eigenVals, eigenVecs, w);
		
		// eigenvalues are unordered, but eigenvectors
		// are in (column) order corresponding
		// to their respective eigenvalues.
		if(ellipsoidCenter.size()) ellipsoidCenter.clear();
		if(ellipsoidAxes.size()) ellipsoidAxes.clear();
		if(ellipsoidAxisLengths.size()) ellipsoidAxisLengths.clear();
		for(int ii = 0; ii < 3; ++ii)
		{
			std::vector< double > eigenVec;
			for(int jj = 0; jj < 3; ++jj)
				eigenVec.push_back(gsl_matrix_get(eigenVecs, jj, ii));
			ellipsoidAxes.push_back(eigenVec);
			ellipsoidAxisLengths.push_back(sqrt(gsl_vector_get(eigenVals, ii)));
			ellipsoidCenter.push_back(center[ii]);
			
			#ifdef DEBUG
			std::cout << "Eigenvector " << ii << " = ";
			std::cout << "[" << eigenVec[0] << "," << eigenVec[1] << "," << eigenVec[2] << "]" << std::endl;
			std::cout << "Diameter " << ii << " = " << 2*ellipsoidAxisLengths[ii] << std::endl;
			#endif
		}
		
		ellipsoidDone = 1;
		
		// clean up
		gsl_matrix_free(mCov);
		gsl_matrix_free(eigenVecs);
		gsl_vector_free(eigenVals);
		gsl_eigen_symmv_free(w);
	}
}

// fit ellipsoid in the form
// Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz = 1
// (nothing fancy, may therefore be unstable; based in parts on
// http://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit)
void HVCSurfaceReconstruct::computeEllipsoidFromContours()
{
	if(HVCSurface)
	{
		#ifdef DEBUG
		std::cout << "Fitting ellipsoid to " << HVCSurface->GetNumberOfCells() << " points" << std::endl;
		#endif
		int status;
		gsl_matrix * mData = gsl_matrix_alloc(HVCSurface->GetNumberOfCells(), 9);
		gsl_matrix * mU = gsl_matrix_alloc(HVCSurface->GetNumberOfCells(), 9);
		gsl_matrix * mV = gsl_matrix_alloc(9, 9);
		gsl_vector * vS = gsl_vector_alloc(9);
		gsl_vector * vWork = gsl_vector_alloc(9);
		gsl_vector * vB1 = gsl_vector_alloc(HVCSurface->GetNumberOfCells());
		gsl_vector * vX = gsl_vector_alloc(9);
		
		initializeData(mData);
		// initialize vector vB1 = (1,1,...,1)^T
		gsl_vector_set_all(vB1, 1.0);
		
		// solve linear system Ax = vB1
		// where A = mData
		status = gsl_linalg_SV_decomp(mData, mV, vS, vWork);
		status = gsl_linalg_SV_solve(mData, mV, vS, vB1, vX);
		
		// recover ellipsoid parameters
		// from algebraic parameters
		// (see lab book for explanation)
		gsl_matrix * mA = gsl_matrix_alloc(3,3);
		gsl_matrix * mE = gsl_matrix_alloc(3,3);
		gsl_vector * vB2 = gsl_vector_alloc(3);
		gsl_vector * vC = gsl_vector_alloc(3);
		gsl_permutation * perm = gsl_permutation_alloc(3);
		double c = -1.0, k;
		int signum;
		
		setEllipsoidMatrix(mA, vX);
		setEllipsoidMatrix(mE, vX);
		gsl_vector_set(vB2, 0, 2*gsl_vector_get(vX, 6));
		gsl_vector_set(vB2, 1, 2*gsl_vector_get(vX, 7));
		gsl_vector_set(vB2, 2, 2*gsl_vector_get(vX, 8));
		
		#ifdef DEBUG
		std::cout << std::endl;
		std::cout << "Algebraic ellipsoid fitting results" << std::endl;
		std::cout << "mE:" << std::endl;
		for(int ii = 0; ii < 3; ++ii)
		{
			std::cout << "[";
			for(int jj = 0; jj < 2; ++jj)
				std::cout << gsl_matrix_get(mE, ii, jj) << ", ";
			std::cout << gsl_matrix_get(mE, ii, 2) << "]" << std::endl;
		}
		std::cout << "vB2:" << std::endl;
		std::cout << "[";
		for(int ii = 0; ii < 2; ++ii)
			std::cout << gsl_vector_get(vB2, ii) << ", ";
		std::cout << gsl_vector_get(vB2, 2) << "]" << std::endl;
		std::cout << std::endl;
		#endif
		
		// solve for center
		status = gsl_linalg_LU_decomp(mA, perm, &signum);
		status = gsl_linalg_LU_solve(mA, perm, vB2, vC);
		gsl_vector_scale(vC, -0.5);
		
		// compute scaling factor of ellipsoid matrix
		gsl_blas_ddot(vB2, vC, &k);
		k = 1/(-0.5*k - c);
		gsl_matrix_scale(mE, k);
		
		#ifdef DEBUG
		std::cout << "Recovering ellipsoid parameters" << std::endl;
		std::cout << "k:" << std::endl;
		std::cout << k << std::endl;
		std::cout << "mE:" << std::endl;
		for(int ii = 0; ii < 3; ++ii)
		{
			std::cout << "[";
			for(int jj = 0; jj < 2; ++jj)
				std::cout << gsl_matrix_get(mE, ii, jj) << ", ";
			std::cout << gsl_matrix_get(mE, ii, 2) << "]" << std::endl;
		}
		std::cout << "vB2:" << std::endl;
		std::cout << "[";
		for(int ii = 0; ii < 2; ++ii)
			std::cout << gsl_vector_get(vB2, ii) << ", ";
		std::cout << gsl_vector_get(vB2, 2) << "]" << std::endl;
		std::cout << "vC:" << std::endl;
		std::cout << "[";
		for(int ii = 0; ii < 2; ++ii)
			std::cout << gsl_vector_get(vC, ii) << ", ";
		std::cout << gsl_vector_get(vC, 2) << "]" << std::endl;
		std::cout << std::endl;
		#endif
		
		// diagonalization
		gsl_vector * eigenVals = gsl_vector_alloc(3);
		gsl_matrix * eigenVecs = gsl_matrix_alloc(3, 3);
		gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(3);
		status = gsl_eigen_symmv(mE, eigenVals, eigenVecs, w);
		
		// eigenvalues are unordered, but eigenvectors
		// are in (column) order corresponding
		// to their respective eigenvalues.
		if(ellipsoidAxes.size()) ellipsoidAxes.clear();
		if(ellipsoidAxisLengths.size()) ellipsoidAxisLengths.clear();
		for(int ii = 0; ii < 3; ++ii)
		{
			std::vector< double > eigenVec;
			for(int jj = 0; jj < 3; ++jj)
				eigenVec.push_back(gsl_matrix_get(eigenVecs, jj, ii));
			ellipsoidAxes.push_back(eigenVec);
			ellipsoidAxisLengths.push_back(sqrt(1/gsl_vector_get(eigenVals, ii)));
			ellipsoidCenter.push_back(gsl_vector_get(vC, ii));
			
			#ifdef DEBUG
			std::cout << "Eigenvector " << ii << " = ";
			std::cout << "[" << eigenVec[0] << "," << eigenVec[1] << "," << eigenVec[2] << "]" << std::endl;
			std::cout << "Diameter " << ii << " = " << 2*ellipsoidAxisLengths[ii] << std::endl;
			#endif
		}
		
		ellipsoidDone = 1;
		
		// clean up
		gsl_matrix_free(mData), gsl_matrix_free(mU), gsl_matrix_free(mV);
		gsl_vector_free(vS), gsl_vector_free(vWork), gsl_vector_free(vB1), gsl_vector_free(vX);
		
		gsl_matrix_free(mA), gsl_matrix_free(mE);
		gsl_permutation_free(perm);
		gsl_vector_free(vB2), gsl_vector_free(vC);
		
		gsl_matrix_free(eigenVecs);
		gsl_vector_free(eigenVals);
		gsl_eigen_symmv_free(w);
	}
}

std::vector< int > HVCSurfaceReconstruct::getOrderedContourIndices()
{
	std::map< double, int > contourIndexMap;
	for(int ii = 0; ii < HVCContours->GetNumberOfCells(); ++ii)
	{
		double pt[3];
		HVCContours->GetCell(ii)->GetPoints()->GetPoint(0, pt);
		contourIndexMap.insert(std::pair< double, int >(pt[2], ii));
	}
	
	#ifdef DEBUG
	std::cout << " Barreloid contour order:" << std::endl;
	#endif
	std::vector< int > orderedIndices;
	std::map< double, int >::const_iterator contourIndexMapIt;
	for(contourIndexMapIt = contourIndexMap.begin(); contourIndexMapIt != contourIndexMap.end(); ++contourIndexMapIt)
	{
		orderedIndices.push_back(contourIndexMapIt->second);
		#ifdef DEBUG
		std::cout << " z = " << contourIndexMapIt->first << "\tindex = " << contourIndexMapIt->second << std::endl;
		#endif
	}
	
	return orderedIndices;
}

UnstructuredGridPointerType HVCSurfaceReconstruct::triangulateAdjacentContours ( int index1, int index2 )
{
	Delaunay3DFilterPointerType del3DTriangulation = Delaunay3DFilterPointerType::New();
	PolyDataPointerType contourSubset = PolyDataPointerType::New();
	IdListPointerType subsetIDs = IdListPointerType::New();
	
	subsetIDs->InsertId(0, index1);
	subsetIDs->InsertId(1, index2);
	contourSubset->Allocate(1);
	contourSubset->CopyCells(HVCContours, subsetIDs);
	contourSubset->Update();
	
	del3DTriangulation->SetInput(contourSubset);
	del3DTriangulation->Update();
	
	return del3DTriangulation->GetOutput();
}

// triangulate first/last polygon including center point
// to ensure a closed surface
UnstructuredGridPointerType HVCSurfaceReconstruct::triangulateEndContours ( int index, bool lastContour )
{
	PolyDataPointerType endPolyData = PolyDataPointerType::New();
	PointsPointerType endPolyPts = PointsPointerType::New();
	PolygonPointerType endPoly = PolygonPointerType::New();
	endPolyData->Allocate(1);
	endPolyPts->SetDataTypeToFloat();
	endPoly->GetPointIds()->SetNumberOfIds(HVCContours->GetCell(index)->GetNumberOfPoints()+1);
	
	PointsPointerType cellPts = HVCContours->GetCell(index)->GetPoints();
	unsigned int cellPtNr = cellPts->GetNumberOfPoints();
	for(int ii = 0; ii < cellPtNr; ++ii)
	{
		double pt[3];
		cellPts->GetPoint(ii, pt);
		endPolyPts->InsertPoint(ii, pt);
		endPoly->GetPointIds()->InsertId(ii, ii);
	}
	
	int subId;
	double center[3], pCenter[3], * weights = new double[cellPtNr];
	HVCContours->GetCell(index)->GetParametricCenter(pCenter);
	HVCContours->GetCell(index)->EvaluateLocation(subId, pCenter, center, weights);
	delete [] weights;
	
	// put center pt at end of section
	// to close surface in 3D
	if(!lastContour)
		center[2] -= 50;
	else
		center[2] += 50;
	#ifdef DEBUG
	if(lastContour)
	{
		std::cout << " Last contour! Center @ [" << center[0] << "," << center[1] << "," << center[2] << "]" << std::endl;
	}
	else if(!lastContour)
	{
		std::cout << " First contour! Center @ [" << center[0] << "," << center[1] << "," << center[2] << "]" << std::endl;
	}
	#endif
	
	endPolyPts->InsertPoint(cellPtNr, center);
	endPoly->GetPointIds()->InsertId(cellPtNr, cellPtNr);
	endPolyData->InsertNextCell(endPoly->GetCellType(), endPoly->GetPointIds());
	endPolyData->SetPoints(endPolyPts);
	
	Delaunay3DFilterPointerType del3DTriangulation = Delaunay3DFilterPointerType::New();
	del3DTriangulation->SetInput(endPolyData);
	del3DTriangulation->Update();
	
	return del3DTriangulation->GetOutput();
}

void HVCSurfaceReconstruct::initializeData ( gsl_matrix* mData )
{
	// initialize data matrix
	for(int ii = 0; ii < HVCSurface->GetNumberOfCells(); ++ii)
	{
		PointsPointerType cellPts = HVCSurface->GetCell(ii)->GetPoints();
		double * pt = new double[9], x = 0, y = 0, z = 0;
		cellPts->GetPoint(0, pt);
		cellPts->GetPoint(1, pt+3);
		cellPts->GetPoint(2, pt+6);
		for(int jj = 0; jj < 3; ++jj)
		{
			x += pt[jj*3]/3;
			y += pt[jj*3+1]/3;
			z += pt[jj*3+2]/3;
		}
		delete [] pt;
		
		gsl_matrix_set(mData, ii, 0, x*x);
		gsl_matrix_set(mData, ii, 1, y*y);
		gsl_matrix_set(mData, ii, 2, z*z);
		gsl_matrix_set(mData, ii, 3, 2*x*y);
		gsl_matrix_set(mData, ii, 4, 2*x*z);
		gsl_matrix_set(mData, ii, 5, 2*y*z);
		gsl_matrix_set(mData, ii, 6, 2*x);
		gsl_matrix_set(mData, ii, 7, 2*y);
		gsl_matrix_set(mData, ii, 8, 2*z);
		
		#ifdef DEBUG
		// print row
		#endif
	}
}

void HVCSurfaceReconstruct::setEllipsoidMatrix ( gsl_matrix* mA, gsl_vector* vX )
{
	gsl_matrix_set(mA, 0, 0, gsl_vector_get(vX, 0));
	gsl_matrix_set(mA, 1, 1, gsl_vector_get(vX, 1));
	gsl_matrix_set(mA, 2, 2, gsl_vector_get(vX, 2));
	gsl_matrix_set(mA, 0, 1, gsl_vector_get(vX, 3));
	gsl_matrix_set(mA, 1, 0, gsl_vector_get(vX, 3));
	gsl_matrix_set(mA, 0, 2, gsl_vector_get(vX, 4));
	gsl_matrix_set(mA, 2, 0, gsl_vector_get(vX, 4));
	gsl_matrix_set(mA, 1, 2, gsl_vector_get(vX, 5));
	gsl_matrix_set(mA, 2, 1, gsl_vector_get(vX, 5));
}

void HVCSurfaceReconstruct::initializeCovarianceMatrix(gsl_matrix* mCov, double mean[3])
{
	// initialize covariance matrix
	std::vector < double * > dataPts;
	mean[0] = 0;
	mean[1] = 0;
	mean[2] = 0;
	for(int ii = 0; ii < HVCSurface->GetNumberOfCells(); ++ii)
	{
		PointsPointerType cellPts = HVCSurface->GetCell(ii)->GetPoints();
		double * pt = new double[9], x = 0, y = 0, z = 0;
		cellPts->GetPoint(0, pt);
		cellPts->GetPoint(1, pt+3);
		cellPts->GetPoint(2, pt+6);
		for(int jj = 0; jj < 3; ++jj)
		{
			x += pt[jj*3]/3;
			y += pt[jj*3+1]/3;
			z += pt[jj*3+2]/3;
		}
		delete [] pt;
		double * tmpPt = new double[3];
		tmpPt[0] = x;
		tmpPt[1] = y;
		tmpPt[2] = z;
		mean[0] += x;
		mean[1] += y;
		mean[2] += z;
		dataPts.push_back(tmpPt);
		
		#ifdef DEBUG
		// print row
		#endif
	}
	
	mean[0] /= dataPts.size();
	mean[1] /= dataPts.size();
	mean[2] /= dataPts.size();
	
	gsl_matrix_set_zero(mCov);
	for(int i = 0; i < dataPts.size(); ++i)
	{
		for(int j = 0; j < 3; ++j)
			for(int k = 0; k < 3; ++k)
			{
				double * element = gsl_matrix_ptr(mCov, j, k);
				*element += (dataPts[i][j] - mean[j])*(dataPts[i][k] - mean[k]);
			}
	}
	for(int j = 0; j < 3; ++j)
		for(int k = 0; k < 3; ++k)
		{
			double * element = gsl_matrix_ptr(mCov, j, k);
			*element /= (dataPts.size()-1);
		}
}
