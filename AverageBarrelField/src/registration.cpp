/****************************************************************************/
/*                                                                          */
/* Program:   AverageBarrelField                                            */
/*                                                                          */
/* File:      registration.cpp                                              */
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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "registration.h"

#define ORIGINAL_RECONS
// #define REG_ACCURACY
#define GLOBAL_CENTER_D2
// #define GLOBAL_CENTER_C2
// #define CONVEX_HULL_SURFS
// Macro defs for quantifying other registration procedures
// #define REG_OPTIMIZATION
// #define CENTER_REG
// #define PIA_WM_REG

Registration::Registration(AmiraSpatialGraph * inputSpatialGraph)
{
	spatialGraph = inputSpatialGraph;
	initializeConstants();
};

Registration::Registration()
{
	spatialGraph = NULL;
	spatialGraph2 = NULL;
	initializeConstants();
};

Registration::~Registration()
{
	barrelLabels.clear();
	int2Labels.clear();
	avgTopDist.clear();
	avgPiaWMDist.clear();
	avgBarrelHeight.clear();
	avgBarrelArea.clear();
};

std::list< AmiraSpatialGraph * > Registration::averageBarrelField(int nrOfBarrelFields, std::list< AmiraSpatialGraph * > allBarrelFields, std::vector< const char * > fileNames)
{
	nrOfPointSets = nrOfBarrelFields;
	bool checkedZDirection = 0;
	std::map< int, Column * > * barrels = new std::map< int, Column * >[nrOfBarrelFields];
	#ifdef PIA_WM_REG
	std::map< int, Column * > * trueBarrels = new std::map< int, Column * >[nrOfBarrelFields];
	#endif
	std::list< int >::const_iterator labelIt;
	std::list< AmiraSpatialGraph * >::const_iterator allBarrelFieldsIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ii = 0;
		for(allBarrelFieldsIt = allBarrelFields.begin(); allBarrelFieldsIt != allBarrelFields.end(); ++allBarrelFieldsIt, ++ii)
		{
			PolyDataPointerType contour = PolyDataPointerType::New();
			if((*allBarrelFieldsIt)->extractLandmark(*labelIt, contour))
			{
				#ifndef PIA_WM_REG
				double * top = new double[3], * bottom = new double[3];
				int subId, subId2;
				double pCoords[3], pCoords2[3], * weights = new double[contour->GetCell(0)->GetNumberOfPoints()], * weights2 = new double[contour->GetCell(1)->GetNumberOfPoints()];
				contour->GetCell(0)->GetParametricCenter(pCoords);
				contour->GetCell(0)->EvaluateLocation(subId, pCoords, top, weights);
				contour->GetCell(1)->GetParametricCenter(pCoords2);
				contour->GetCell(1)->EvaluateLocation(subId2, pCoords2, bottom, weights2);
				Column * thisCol;
				if(top[2] < bottom[2])
					thisCol = new Column(contour, top, bottom);
				else
					thisCol = new Column(contour, bottom, top);
				barrels[ii].insert(std::pair< int, Column * >(*labelIt, thisCol));
				if(!checkedZDirection)
				{
#ifndef ORIGINAL_RECONS
					if(thisCol->top[2] >= thisCol->bottom[2])
						zReversed = 1;
					else
						zReversed = 0;
#elif defined ORIGINAL_RECONS
					zReversed = 1;
#endif
					checkedZDirection = 1;
					std::cout << "zReversed = " << zReversed << std::endl;
				}
				#elif defined PIA_WM_REG
// 				std::cout << "Loading barrel and column " << int2Labels[*labelIt] << std::endl;
// 				contour->Print(std::cout);
				double * center1 = new double[3], * center2 = new double[3], * center3 = new double[3], * center4 = new double[3];
				int subId1, subId2, subId3, subId4;
				double pCoords1[3], pCoords2[3], pCoords3[3], pCoords4[3];
				double * weights1 = new double[contour->GetCell(0)->GetNumberOfPoints()];
				double * weights2 = new double[contour->GetCell(1)->GetNumberOfPoints()];
				double * weights3 = new double[contour->GetCell(2)->GetNumberOfPoints()];
				double * weights4 = new double[contour->GetCell(3)->GetNumberOfPoints()];
				contour->GetCell(0)->GetParametricCenter(pCoords1);
				contour->GetCell(0)->EvaluateLocation(subId1, pCoords1, center1, weights1);
				contour->GetCell(1)->GetParametricCenter(pCoords2);
				contour->GetCell(1)->EvaluateLocation(subId2, pCoords2, center2, weights2);
				contour->GetCell(2)->GetParametricCenter(pCoords3);
				contour->GetCell(2)->EvaluateLocation(subId3, pCoords3, center3, weights3);
				contour->GetCell(3)->GetParametricCenter(pCoords4);
				contour->GetCell(3)->EvaluateLocation(subId4, pCoords4, center4, weights4);
				Column * thisBarrel, * thisColumn;
				thisBarrel = new Column(contour, center1, center2);
				thisColumn = new Column(contour, center3, center4);
// 				std::cout << "Column top @[" << thisColumn->top[0] << "," << thisColumn->top[1] << "," << thisColumn->top[2] << "]" << std::endl;
// 				std::cout << "Column bottom @[" << thisColumn->bottom[0] << "," << thisColumn->bottom[1] << "," << thisColumn->bottom[2] << "]" << std::endl;
// 				std::cout << "Barrel top @[" << thisBarrel->top[0] << "," << thisBarrel->top[1] << "," << thisBarrel->top[2] << "]" << std::endl;
// 				std::cout << "Barrel bottom @[" << thisBarrel->bottom[0] << "," << thisBarrel->bottom[1] << "," << thisBarrel->bottom[2] << "]" << std::endl;
				barrels[ii].insert(std::pair< int, Column * >(*labelIt, thisColumn));
				trueBarrels[ii].insert(std::pair< int, Column * >(*labelIt, thisBarrel));
				zReversed = 0;
				checkedZDirection = 1;
				#endif
			}
		}
	}
	for(int ii = 0; ii < nrOfBarrelFields; ++ii)
	{
		TransformPointerType bfTrans = TransformPointerType::New();
		bfTrans->Identity();
		transformVec.push_back(bfTrans);
	}
	
	// move all centroids to the origin
	// then compute optimal rotation matrix
	// for each barrel field wrt 1st barrel field
	// using the Kabsch algorithm
	alignBarrelFieldCentroids(barrels, nrOfBarrelFields);
	
	std::list< AmiraSpatialGraph * > regBarrelFields;
	allBarrelFieldsIt = allBarrelFields.begin();
	bool success = 0;
	// one-step version
	// not necessarily (even locally) optimal yet
	for(int ii = 0; ii < nrOfBarrelFields; ++ii, ++allBarrelFieldsIt)
	{
		#ifdef PIA_WM_REG
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(trueBarrels[ii].find(ID) != trueBarrels[ii].end())
				trueBarrels[ii][ID]->translateColumn(transformVec[ii]->GetPosition());
		}
		#endif
		if(ii)
		{
			gsl_matrix * mU = computeOptimalRotation(barrels[0], barrels[ii]);
			if(mU)
			{
				success = 1;
				double * transformation[4];
				for(int jj = 0; jj < 4; ++jj)
					transformation[jj] = new double[4];
				for(int jj = 0; jj < 3; ++jj)
					for(int kk = 0; kk < 3; ++kk)
						transformation[jj][kk] = gsl_matrix_get(mU, jj, kk);
				for(int jj = 0; jj < 3; ++jj)
				{
					transformation[jj][3] = 0;
					transformation[3][jj] = 0;
				}
				transformation[3][3] = 1;
				
				TransformPointerType rot = TransformPointerType::New();
				HomogeneousMatrixPointerType mRot = HomogeneousMatrixPointerType::New();
				mRot->Identity();
				for(int jj = 0; jj < 3; ++jj)
					for(int kk = 0; kk < 3; ++kk)
						mRot->SetElement(jj, kk, gsl_matrix_get(mU, jj, kk));
				rot->SetMatrix(mRot);
				rot->Concatenate(transformVec[ii]);
				rot->Update();
				transformVec[ii] = rot;
				
				AmiraSpatialGraph * regSpatialGraph = new AmiraSpatialGraph();
				regSpatialGraph->setTransformation(transformation);
				for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
					if((*allBarrelFieldsIt)->isLabelInSpatialGraph(*labelIt))
					{
						regSpatialGraph->addPolyDataObject(barrels[ii][*labelIt]->contours, *labelIt);
						barrels[ii][*labelIt]->rotateColumn(mU); // order important b/c at the moment transformation is applied to SpatialGraph separately!
						#ifdef PIA_WM_REG
						trueBarrels[ii][*labelIt]->rotateColumn(mU);
						#endif
					}
				regSpatialGraph->applyTransformation();
				regBarrelFields.push_back(regSpatialGraph);
				gsl_matrix_free(mU);
			}
		}
		else
		{
			AmiraSpatialGraph * regSpatialGraph = new AmiraSpatialGraph();
			for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
				if((*allBarrelFieldsIt)->isLabelInSpatialGraph(*labelIt))
					regSpatialGraph->addPolyDataObject(barrels[ii][*labelIt]->contours, *labelIt);
			regBarrelFields.push_back(regSpatialGraph);
		}
	}
	
	// Orthogonal Procrustes alignment
	// during first iteration, arbitrarily choose 
	// first BF as reference; after that, use mean shape
// 	std::map< int, Column * > meanShape;
// 	std::vector< AmiraSpatialGraph * > regSpatialGraphVec;
// 	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
// 	{
// 		int ID = *labelIt;
// 		if(barrels[0].find(ID) != barrels[0].end())
// 		{
// 			Column * newCol = new Column();
// 			for(int ii = 0; ii < 3; ++ii)
// 			{
// 				newCol->bottom[ii] = barrels[0][ID]->bottom[ii];
// 				newCol->top[ii] = barrels[0][ID]->top[ii];
// 			}
// 			meanShape.insert(std::pair< int, Column * >(ID, newCol));
// 		}
// 	}
// 	
// 	double residuals = getResiduals(barrels);
// 	double tol = 1E-03, delta = 0;
// 	unsigned long iterCnt = 1;
// 	do
// 	{
// 		std::cout << std::endl;
// 		std::cout << "Iteration " << iterCnt << std::endl;
// 		bool firstPass = 1;
// 		if(regSpatialGraphVec.size())
// 			firstPass = 0;
// 		allBarrelFieldsIt = allBarrelFields.begin();
// 		for(int ii = 0; ii < nrOfBarrelFields; ++ii, ++allBarrelFieldsIt)
// 		{
// 			gsl_matrix * mU = computeOptimalRotation(meanShape, barrels[ii]);
// 			if(mU)
// 			{
// 				success = 1;
// 				double * transformation[4];
// 				for(int jj = 0; jj < 4; ++jj)
// 					transformation[jj] = new double[4];
// 				for(int jj = 0; jj < 3; ++jj)
// 					for(int kk = 0; kk < 3; ++kk)
// 						transformation[jj][kk] = gsl_matrix_get(mU, jj, kk);
// 				for(int jj = 0; jj < 3; ++jj)
// 				{
// 					transformation[jj][3] = 0;
// 					transformation[3][jj] = 0;
// 				}
// 				transformation[3][3] = 1;
// 				
// 				TransformPointerType rot = TransformPointerType::New();
// 				HomogeneousMatrixPointerType mRot = HomogeneousMatrixPointerType::New();
// 				mRot->Identity();
// 				for(int jj = 0; jj < 3; ++jj)
// 					for(int kk = 0; kk < 3; ++kk)
// 						mRot->SetElement(jj, kk, gsl_matrix_get(mU, jj, kk));
// 				
// 				std::cout << "\tBarrel field ID " << ii << std::endl;
// 				std::cout << "\tRotation matrix:" << std::endl;
// 				mRot->Print(std::cout);
// 				std::flush(std::cout);
// 				
// 				rot->SetMatrix(mRot);
// 				rot->Concatenate(transformVec[ii]);
// 				rot->Update();
// 				transformVec[ii] = rot;
// 				
// 				AmiraSpatialGraph * regSpatialGraph;
// 				if(firstPass)
// 					regSpatialGraph = new AmiraSpatialGraph();
// 				else
// 					regSpatialGraph = regSpatialGraphVec[ii];
// 				regSpatialGraph->setTransformation(transformation);
// 				for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
// 					if((*allBarrelFieldsIt)->isLabelInSpatialGraph(*labelIt))
// 					{
// 						if(firstPass)
// 							regSpatialGraph->addPolyDataObject(barrels[ii][*labelIt]->contours, *labelIt);
// 						barrels[ii][*labelIt]->rotateColumn(mU); // order important b/c at the moment transformation is applied to SpatialGraph separately!
// 					}
// 				regSpatialGraph->applyTransformation();
// 				if(firstPass)
// 					regSpatialGraphVec.push_back(regSpatialGraph);
// 				gsl_matrix_free(mU);
// 			}
// 		}
// 		
// 		std::cout << "Computing new mean shape..." << std::endl;
// 		computeMeanShape(&meanShape, barrels);
// 		
// 		double newResiduals = getResiduals(barrels);
// 		delta = fabs(residuals-newResiduals)/residuals;
// 		
// 		std::cout << "old residuals =\t" << residuals << std::endl;
// 		std::cout << "new residuals =\t" << newResiduals << std::endl;
// 		std::cout << "delta =\t" << delta << std::endl;
// 		std::cout << "tol =\t" << tol << std::endl;
// 		std::cout << "*********************************" << std::endl;
// 		
// 		residuals = newResiduals;
// 		++iterCnt;
// 	} while(delta > tol && iterCnt <= 10);
// 	for(int ii = 0; ii < regSpatialGraphVec.size(); ++ii)
// 		regBarrelFields.push_back(regSpatialGraphVec[ii]);
	
	if(success)
	{
		computeAverageBarrelField(barrels, fileNames);
		computeAverageSurfaces(fileNames.back());
		computeAverageAxesField();
		#ifdef REG_ACCURACY
		#ifndef PIA_WM_REG
		alignmentQuality(barrels, fileNames);
		#elif defined PIA_WM_REG
		alignmentQuality(trueBarrels, fileNames);
		#endif
		#endif
	}
	return regBarrelFields;
};

/****************************************************************************/
/*moves all centroids to the origin                                         */
/****************************************************************************/
void Registration::alignBarrelFieldCentroids(std::map< int, Column * > * barrels, int nrOfBarrelFields)
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
		
		TransformPointerType shiftTrans = TransformPointerType::New();
		shiftTrans->Translate(shift);
		transformVec[ii] = shiftTrans;
	}
};

/****************************************************************************/
/*compute optimal rotation matrix in the least squares sense                */
/*using the Kabsch algorithm:                                               */
/*minimize sum_i (U x_i - y_i)^2 where U = optimal rotation matrix,         */
/*x_i = ith landmark of matchBF, y_i = ith landmark of refBF                */
/****************************************************************************/
gsl_matrix * Registration::computeOptimalRotation(std::map< int, Column * > refBF, std::map< int, Column * > matchBF)
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
	#ifndef CENTER_REG
	gsl_matrix * mX = gsl_matrix_alloc(3, 2*commonLandmarks.size());
	gsl_matrix * mY = gsl_matrix_alloc(3, 2*commonLandmarks.size());
	#endif
	#ifdef CENTER_REG
	gsl_matrix * mX = gsl_matrix_alloc(3, commonLandmarks.size());
	gsl_matrix * mY = gsl_matrix_alloc(3, commonLandmarks.size());
	#endif
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
		#ifndef CENTER_REG
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
		#endif
		#ifdef CENTER_REG
		double refCenter[3], matchCenter[3];
		vtkMath::Add(refBF[ID]->top, refBF[ID]->bottom, refCenter);
		vtkMath::MultiplyScalar(refCenter, 0.5);
		vtkMath::Add(matchBF[ID]->top, matchBF[ID]->bottom, matchCenter);
		vtkMath::MultiplyScalar(matchCenter, 0.5);
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
		#endif
	}
	
	// compute X Y^t
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mX, mY, 0.0, mCov);
// 	std::cout << "mCov:" << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
// 	}
	// compute sign of determinant of mCov (right-/left-handed)
	gsl_matrix_memcpy(mLU, mCov);
	int sign, detSign;
	gsl_linalg_LU_decomp(mLU, permLU, &sign);
// 	std::cout << "mLU:" << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mLU, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mLU, ii, 2) << "]" << std::endl;
// 	}
	detSign = gsl_linalg_LU_sgndet(mLU, sign);
// 	std::cout << "detSign = " << detSign << std::endl;
	// compute SVD of mCov
	gsl_vector * vTmp = gsl_vector_alloc(3);
	gsl_linalg_SV_decomp(mCov, mV, vS, vTmp);
	// compute mU = mV * diag(1,1,detSign) * mCov^t
	// note: gsl already computes mV so that mCov = U * S * V^t
// 	std::cout << "mCov after SVD:" << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
// 	}
// 	std::cout << "mV after SVD:" << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mV, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mV, ii, 2) << "]" << std::endl;
// 	}
	gsl_matrix_set_identity(mId);
	gsl_matrix_set(mId, 2, 2, detSign);
	gsl_matrix * mTmp = gsl_matrix_alloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mId, mCov, 0.0, mTmp);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mV, mTmp, 0.0, mU);
// 	std::cout << "mU:" << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mU, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mU, ii, 2) << "]" << std::endl;
// 	}
	
	gsl_matrix_free(mTmp), gsl_matrix_free(mCov), gsl_matrix_free(mId), gsl_matrix_free(mLU), gsl_matrix_free(mV), gsl_matrix_free(mX), gsl_matrix_free(mY);
	gsl_permutation_free(permLU);
	gsl_vector_free(vS), gsl_vector_free(vTmp);
	
	return mU;
};

void Registration::computeMeanShape(std::map< int, Column * > * meanShape, std::map< int, Column * > * barrels)
{
	std::map< int, double[3] > meanTop, meanBottom;
	std::map< int, int > barrelFrequency;
	std::list< int >::const_iterator labelIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrels[0].find(ID) != barrels[0].end() && meanShape->find(ID) != meanShape->end())
		{
			int cnt = 0;
			double avgTop[3], avgBottom[3];
			avgTop[0] = avgTop[1] = avgTop[2] = 0;
			avgBottom[0] = avgBottom[1] = avgBottom[2] = 0;
			for(int ii = 0; ii < nrOfPointSets; ++ii)
				if(barrels[ii].find(ID) != barrels[ii].end())
				{
					++cnt;
					for(int jj = 0; jj < 3; ++jj)
					{
						avgTop[jj] += barrels[ii][ID]->top[jj];
						avgBottom[jj] += barrels[ii][ID]->bottom[jj];
					}
				}
			vtkMath::MultiplyScalar(avgTop, 1/(double)cnt);
			vtkMath::MultiplyScalar(avgBottom, 1/(double)cnt);
			for(int ii = 0; ii < 3; ++ii)
			{
				meanShape[0][ID]->top[ii] = avgTop[ii];
				meanShape[0][ID]->bottom[ii] = avgBottom[ii];
			}
		}
	}
};

double Registration::getResiduals(std::map< int, Column * > * barrels)
{
	double res = 0;
	std::list< int >::const_iterator labelIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		for(int ii = 0; ii < nrOfPointSets; ++ii)
			for(int jj = 0; jj < nrOfPointSets; ++jj)
				if(barrels[ii].find(ID) != barrels[ii].end() && barrels[jj].find(ID) != barrels[jj].end())
				{
					res += vtkMath::Distance2BetweenPoints(barrels[ii][ID]->top, barrels[jj][ID]->top);
					res += vtkMath::Distance2BetweenPoints(barrels[ii][ID]->bottom, barrels[jj][ID]->bottom);
				}
	}
	
	return res;
};

void Registration::computeAverageBarrelField(std::map< int, Column * > * matchedBarrels, std::vector< const char * > fileNames)
{
	std::flush(std::cout << "computing avg barrel field" << std::endl);
	std::list< int >::const_iterator labelIt;
	spatialGraph = new AmiraSpatialGraph;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		#ifndef REG_OPTIMIZATION
		double * avgTop = new double[3], * avgBottom = new double[3];
		avgTop[0] = avgTop[1] = avgTop[2] = 0;
		avgBottom[0] = avgBottom[1] = avgBottom[2] = 0;
		int avgCnt = 0;
		for(int ii = 0; ii < nrOfPointSets; ++ii)
			if(matchedBarrels[ii].find(ID) != matchedBarrels[ii].end())
			{
				++avgCnt;
				for(int jj = 0; jj < 3; ++jj)
				{
					avgTop[jj] += matchedBarrels[ii][ID]->top[jj];
					avgBottom[jj] += matchedBarrels[ii][ID]->bottom[jj];
				}
			}
		
		if(avgCnt)
		{
			double * axis = new double[3], * center = new double[3];
			for(int ii = 0; ii < 3; ++ii)
			{
				avgTop[ii] = avgTop[ii]/avgCnt;
				avgBottom[ii] = avgBottom[ii]/avgCnt;
				axis[ii] = avgTop[ii] - avgBottom[ii];
				center[ii] = (avgTop[ii] + avgBottom[ii])*0.5;
			}
			Column * newAvgBarrel = new Column;
			newAvgBarrel->top = avgTop;
			newAvgBarrel->bottom = avgBottom;
			avgBarrels.insert(std::pair< int, Column * >(ID, newAvgBarrel));
			avgAxes.insert(std::pair< int, double * >(ID, axis));
			avgCenters.insert(std::pair< int, double * >(ID, center));
		}
		else
		{
			delete [] avgTop, delete [] avgBottom;
		}
		#endif
		
		#ifdef REG_OPTIMIZATION
		double * avgTop = new double[3], * avgBottom = new double[3];
		double * axis = new double[3], * center = new double[3];
		center[0] = center[1] = center[2] = 0;
		axis[0] = axis[1] = axis[2] = 0;
		avgTop[0] = avgTop[1] = avgTop[2] = 0;
		avgBottom[0] = avgBottom[1] = avgBottom[2] = 0;
		int avgCnt = 0;
		for(int ii = 0; ii < nrOfPointSets; ++ii)
			if(matchedBarrels[ii].find(ID) != matchedBarrels[ii].end())
			{
				++avgCnt;
				#ifdef CENTER_REG
				double tmpAxis[3];
				for(int jj = 0; jj < 3; ++jj)
				{
					center[jj] += 0.5*(matchedBarrels[ii][ID]->top[jj] + matchedBarrels[ii][ID]->bottom[jj]);
					tmpAxis[jj] = matchedBarrels[ii][ID]->top[jj] - matchedBarrels[ii][ID]->bottom[jj];
				}
				normalize(tmpAxis);
				for(int jj = 0; jj < 3; ++jj)
					axis[jj] += tmpAxis[jj];
				#endif
				#ifdef PIA_WM_REG
				for(int jj = 0; jj < 3; ++jj)
				{
					avgTop[jj] += matchedBarrels[ii][ID]->top[jj];
					avgBottom[jj] += matchedBarrels[ii][ID]->bottom[jj];
				}
				#endif
			}
		
		if(avgCnt)
		{
			#ifdef CENTER_REG
			normalize(axis);
			for(int ii = 0; ii < 3; ++ii)
			{
				center[ii] /= avgCnt;
				avgTop[ii] = center[ii] + 0.5*avgBarrelHeight[ID]*axis[ii];
				avgBottom[ii] = center[ii] - 0.5*avgBarrelHeight[ID]*axis[ii];
			}
			#endif
			#ifdef PIA_WM_REG
			for(int ii = 0; ii < 3; ++ii)
			{
				avgTop[ii] = avgTop[ii]/avgCnt;
				avgBottom[ii] = avgBottom[ii]/avgCnt;
				axis[ii] = avgTop[ii] - avgBottom[ii];
				center[ii] = (avgTop[ii] + avgBottom[ii])*0.5;
			}
			#endif
			
			#ifdef CENTER_REG
			Column * newAvgBarrel = new Column;
			newAvgBarrel->top = avgTop;
			newAvgBarrel->bottom = avgBottom;
			avgBarrels.insert(std::pair< int, Column * >(ID, newAvgBarrel));
			#endif
			#ifdef PIA_WM_REG
			Column * newAvgCol = new Column;
			newAvgCol->top = avgTop;
			newAvgCol->bottom = avgBottom;
			avgColumns.insert(std::pair< int, Column * >(ID, newAvgCol));
			#endif
			avgAxes.insert(std::pair< int, double * >(ID, axis));
			avgCenters.insert(std::pair< int, double * >(ID, center));
		}
		#endif
	}
	// make sure avg axes are divergent towards pia
	enforceAxisDivergence(avgAxes, avgCenters);
	
	#ifndef PIA_WM_REG
	#ifndef REG_ACCURACY
	alignToReferenceFrame(matchedBarrels);
	#endif
	
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(avgAxes.find(ID) != avgAxes.end() && avgCenters.find(ID) != avgCenters.end() && avgBarrels.find(ID) != avgBarrels.end())
		{
			double * colTop = new double[3], * colBottom = new double[3];
			double height = avgBarrels[ID]->getHeight();
			normalize(avgAxes[ID]);
			for(int ii = 0; ii < 3; ++ii)
			{
				avgBarrels[ID]->top[ii] = avgCenters[ID][ii] + 0.5*height*avgAxes[ID][ii];
				avgBarrels[ID]->bottom[ii] = avgCenters[ID][ii] - 0.5*height*avgAxes[ID][ii];
				colTop[ii] = avgBarrels[ID]->top[ii] + avgTopDist[ID]*avgAxes[ID][ii];
				colBottom[ii] = colTop[ii] - avgPiaWMDist[ID]*avgAxes[ID][ii];
			}
			Column * newAvgCol = new Column;
			newAvgCol->top = colTop;
			newAvgCol->bottom = colBottom;
			avgColumns.insert(std::pair< int, Column * >(ID, newAvgCol));
		}
	}
	#elif defined PIA_WM_REG
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(avgAxes.find(ID) != avgAxes.end() && avgCenters.find(ID) != avgCenters.end() && avgColumns.find(ID) != avgColumns.end())
		{
			normalize(avgAxes[ID]);
			for(int ii = 0; ii < 3; ++ii)
			{
				avgCenters[ID][ii] = avgColumns[ID]->top[ii] - (avgTopDist[ID] + 0.5*avgBarrelHeight[ID])*avgAxes[ID][ii];
			}
		}
	}
	#ifndef REG_ACCURACY
	alignToReferenceFrame(matchedBarrels);
	#endif
	//need to loop once more b/c only now everything is aligned to the barrel center, not the column center
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(avgAxes.find(ID) != avgAxes.end() && avgCenters.find(ID) != avgCenters.end() && avgColumns.find(ID) != avgColumns.end())
		{
			double * barrelTop = new double[3], * barrelBottom = new double[3];
			normalize(avgAxes[ID]);
			for(int ii = 0; ii < 3; ++ii)
			{
				barrelTop[ii] = avgCenters[ID][ii] + 0.5*avgBarrelHeight[ID]*avgAxes[ID][ii];
				barrelBottom[ii] = avgCenters[ID][ii] - 0.5*avgBarrelHeight[ID]*avgAxes[ID][ii];
				avgColumns[ID]->top[ii] = barrelTop[ii] + avgTopDist[ID]*avgAxes[ID][ii];
				avgColumns[ID]->bottom[ii] = avgColumns[ID]->top[ii] - avgPiaWMDist[ID]*avgAxes[ID][ii];
			}
			Column * newAvgBarrel = new Column;
			newAvgBarrel->top = barrelTop;
			newAvgBarrel->bottom = barrelBottom;
			avgBarrels.insert(std::pair< int, Column * >(ID, newAvgBarrel));
		}
	}
	#endif
	
	writeAvgBarrelData(matchedBarrels, fileNames);
};

/****************************************************************************/
/*align coordinate system to D2 column:                                     */
/*D2 barrel center @ (0,0,0)                                                */
/*D2 axis = z-axis                                                          */
/*D2-D3 axis = x-axis                                                       */
/*Warning: this version only works correctly when used after calculating    */
/*the average axes and enforcing the divergence constraint                  */
/****************************************************************************/
void Registration::alignToReferenceFrame(std::map< int, Column * > * matchedBarrels)
{
	#ifndef PIA_WM_REG
	if(avgCenters.find(D2) != avgCenters.end() && avgCenters.find(D3) != avgCenters.end() && avgAxes.find(D2) != avgAxes.end() && avgBarrels.find(D2) != avgBarrels.end())
	#endif
	#ifdef PIA_WM_REG
	if(avgCenters.find(D2) != avgCenters.end() && avgCenters.find(D3) != avgCenters.end() && avgAxes.find(D2) != avgAxes.end() && avgColumns.find(D2) != avgColumns.end())
	#endif
	{
		#ifdef GLOBAL_CENTER_D2
		double shiftD2[3];
		shiftD2[0] = -avgCenters[D2][0], shiftD2[1] = -avgCenters[D2][1], shiftD2[2] = -avgCenters[D2][2];
		// rotation so D2 axis is new z axis
		HomogeneousMatrixPointerType rMat = getLocalBarrelCoordinates(avgAxes[D2]);
		// rotation about new x axis by 180 degrees makes z = (0,0,1) (from (0,0,-1) before)
		HomogeneousMatrixPointerType zAxisMat;
		if(zReversed)
			zAxisMat = getNewZAxis(180);
		else
			zAxisMat = getNewZAxis(0);
		// rotation so D2-D3 axis is new x axis
		double currXAxis[3];
		currXAxis[0] = avgCenters[D3][0] - avgCenters[D2][0], currXAxis[1] = avgCenters[D3][1] - avgCenters[D2][1], currXAxis[2] = 0;
		normalize(currXAxis);
		HomogeneousMatrixPointerType xAxisMat = getNewXAxis(acos(currXAxis[0])*180/PI);
		
		#elif defined GLOBAL_CENTER_C2
 		// C2 @ center
 		double shiftD2[3];
 		shiftD2[0] = -avgCenters[C2][0], shiftD2[1] = -avgCenters[C2][1], shiftD2[2] = -avgCenters[C2][2];
 		// rotation so C2 axis is new z axis
 		HomogeneousMatrixPointerType rMat = getLocalBarrelCoordinates(avgAxes[C2]);
 		// rotation about new x axis by 180 degrees makes z = (0,0,1) (from (0,0,-1) before)
		HomogeneousMatrixPointerType zAxisMat;
		if(zReversed)
			zAxisMat = getNewZAxis(180);
		else
			zAxisMat = getNewZAxis(0);
 		// rotation so D2-D3 axis is new x axis
 		double currXAxis[3];
 		currXAxis[0] = avgCenters[C3][0] - avgCenters[C2][0], currXAxis[1] = avgCenters[C3][1] - avgCenters[C2][1], currXAxis[2] = 0;
 		normalize(currXAxis);
 		HomogeneousMatrixPointerType xAxisMat = getNewXAxis(acos(currXAxis[0])*180/PI);
 		std::vector< double * > C2Center, C3Center, C2Axis;
 		for(int ii = 0; ii < nrOfPointSets; ++ii)
 		{
 			double * newC2Center = new double[3], * newC3Center = new double[3], * newC2Axis = new double[3];
 			for(int jj = 0; jj <3; ++jj)
 			{
 				newC2Center[jj] = 0.5*(matchedBarrels[ii][C2]->top[jj] + matchedBarrels[ii][C2]->bottom[jj]);
 				newC3Center[jj] = 0.5*(matchedBarrels[ii][C3]->top[jj] + matchedBarrels[ii][C3]->bottom[jj]);
 				newC2Axis[jj] = matchedBarrels[ii][C2]->top[jj] - matchedBarrels[ii][C2]->bottom[jj];
 			}
 			normalize(newC2Axis);
 			C2Center.push_back(newC2Center);
 			C3Center.push_back(newC3Center);
 			C2Axis.push_back(newC2Axis);
 		}
 		#endif
		
		std::list< int >::const_iterator labelIt;
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(avgCenters.find(ID) != avgCenters.end() && avgAxes.find(ID) != avgAxes.end())
			{
				double hCenter[4], hAxis[4];
				for(int ii = 0; ii < 3; ++ii)
				{
					avgCenters[ID][ii] += shiftD2[ii];
					hCenter[ii] = avgCenters[ID][ii];
					hAxis[ii] = hCenter[ii] + avgAxes[ID][ii];
				}
				hCenter[3] = 1, hAxis[3] = 1;
				rMat->MultiplyPoint(hCenter, hCenter);
				rMat->MultiplyPoint(hAxis, hAxis);
				zAxisMat->MultiplyPoint(hCenter, hCenter);
				zAxisMat->MultiplyPoint(hAxis, hAxis);
				xAxisMat->MultiplyPoint(hCenter, hCenter);
				xAxisMat->MultiplyPoint(hAxis, hAxis);
				for(int ii = 0; ii < 3; ++ii)
				{
					avgCenters[ID][ii] = hCenter[ii];
					avgAxes[ID][ii] = hAxis[ii] - hCenter[ii];
				}
			}
			
			#ifdef GLOBAL_CENTER_D2
			for(int ii = 0; ii < nrOfPointSets; ++ii)
				if(matchedBarrels[ii].find(ID) != matchedBarrels[ii].end())
				{
// 					std::cout << std::endl;
// 					std::cout << "matched barrel " << int2Labels[ID] << "; BF nr. " << ii << std::endl;
// 					std::cout << "registered:" << std::endl;
// 					std::cout << "top\tbottom" << std::endl;
					double hTop[4], hBottom[4];
					for(int jj = 0; jj < 3; ++jj)
					{
// 						std::flush(std::cout << matchedBarrels[ii][ID]->top[jj] << "\t" << matchedBarrels[ii][ID]->bottom[jj] << std::endl);
						matchedBarrels[ii][ID]->top[jj] += shiftD2[jj];
						matchedBarrels[ii][ID]->bottom[jj] += shiftD2[jj];
						hTop[jj] = matchedBarrels[ii][ID]->top[jj];
						hBottom[jj] = matchedBarrels[ii][ID]->bottom[jj];
					}
					hTop[3] = 1, hBottom[3] = 1;
					rMat->MultiplyPoint(hTop, hTop);
					rMat->MultiplyPoint(hBottom, hBottom);
					zAxisMat->MultiplyPoint(hTop, hTop);
					zAxisMat->MultiplyPoint(hBottom, hBottom);
					xAxisMat->MultiplyPoint(hTop, hTop);
					xAxisMat->MultiplyPoint(hBottom, hBottom);
// 					std::cout << "registered to D2:" << std::endl;
// 					std::cout << "top\tbottom" << std::endl;
					for(int jj = 0; jj < 3; ++jj)
					{
						matchedBarrels[ii][ID]->top[jj] = hTop[jj];
						matchedBarrels[ii][ID]->bottom[jj] = hBottom[jj];
// 						std::flush(std::cout << matchedBarrels[ii][ID]->top[jj] << "\t" << matchedBarrels[ii][ID]->bottom[jj] << std::endl);
					}
				}
				
			#elif defined GLOBAL_CENTER_C2
 			for(int ii = 0; ii < nrOfPointSets; ++ii)
 				if(matchedBarrels[ii].find(ID) != matchedBarrels[ii].end())
 				{
 					double shiftInputBFC2[3];
 					shiftInputBFC2[0] = -C2Center[ii][0], shiftInputBFC2[1] = -C2Center[ii][1], shiftInputBFC2[2] = -C2Center[ii][2];
 					// rotation so D2 axis is new z axis
 					HomogeneousMatrixPointerType rMatInputBF = getLocalBarrelCoordinates(C2Axis[ii]);
 					// rotation about new x axis by 180 degrees makes z = (0,0,1) (from (0,0,-1) before)
 					HomogeneousMatrixPointerType zAxisMatInputBF;
					if(zReversed)
						zAxisMatInputBF = getNewZAxis(180);
					else
						zAxisMatInputBF = getNewZAxis(0);
 					// rotation so D2-D3 axis is new x axis
 					double xAxisInputBF[3];
 					xAxisInputBF[0] = C3Center[ii][0] - C2Center[ii][0], xAxisInputBF[1] = C3Center[ii][1] - C2Center[ii][1], xAxisInputBF[2] = 0;
 					normalize(xAxisInputBF);
 					HomogeneousMatrixPointerType xAxisMatInputBF = getNewXAxis(acos(xAxisInputBF[0])*180/PI);
 					
 // 					std::cout << std::endl;
 // 					std::cout << "matched barrel " << int2Labels[ID] << "; BF nr. " << ii << std::endl;
 // 					std::cout << "registered:" << std::endl;
 // 					std::cout << "top\tbottom" << std::endl;
 					double hTop[4], hBottom[4];
 					for(int jj = 0; jj < 3; ++jj)
 					{
 // 						std::flush(std::cout << matchedBarrels[ii][ID]->top[jj] << "\t" << matchedBarrels[ii][ID]->bottom[jj] << std::endl);
 						matchedBarrels[ii][ID]->top[jj] += shiftInputBFC2[jj];
 						matchedBarrels[ii][ID]->bottom[jj] += shiftInputBFC2[jj];
 						hTop[jj] = matchedBarrels[ii][ID]->top[jj];
 						hBottom[jj] = matchedBarrels[ii][ID]->bottom[jj];
 					}
 					hTop[3] = 1, hBottom[3] = 1;
 					rMatInputBF->MultiplyPoint(hTop, hTop);
 					rMatInputBF->MultiplyPoint(hBottom, hBottom);
 					zAxisMatInputBF->MultiplyPoint(hTop, hTop);
 					zAxisMatInputBF->MultiplyPoint(hBottom, hBottom);
 					xAxisMatInputBF->MultiplyPoint(hTop, hTop);
 					xAxisMatInputBF->MultiplyPoint(hBottom, hBottom);
 // 					std::cout << "registered to D2:" << std::endl;
 // 					std::cout << "top\tbottom" << std::endl;
 					for(int jj = 0; jj < 3; ++jj)
 					{
 						matchedBarrels[ii][ID]->top[jj] = hTop[jj];
 						matchedBarrels[ii][ID]->bottom[jj] = hBottom[jj];
 // 						std::flush(std::cout << matchedBarrels[ii][ID]->top[jj] << "\t" << matchedBarrels[ii][ID]->bottom[jj] << std::endl);
 					}
 				}
			#endif
		}
		
		for(int ii = 0; ii < nrOfPointSets; ++ii)
		{
			TransformPointerType alignCenter = TransformPointerType::New();
			TransformPointerType rot1 = TransformPointerType::New();
			TransformPointerType rot2 = TransformPointerType::New();
			TransformPointerType rot3 = TransformPointerType::New();
			alignCenter->Translate(shiftD2);
			rot1->SetMatrix(rMat);
			rot2->SetMatrix(zAxisMat);
			rot3->SetMatrix(xAxisMat);
			
			alignCenter->Concatenate(transformVec[ii]);
			rot1->Concatenate(alignCenter);
			rot2->Concatenate(rot1);
			rot3->Concatenate(rot2);
			rot3->Update();
			transformVec[ii] = rot3;
		}
	}
};

/****************************************************************************/
/*calculate avg deviation of individual barrel fields from average barrel   */
/*field (in new local coordinates though)                                   */
/****************************************************************************/
void Registration::alignmentQuality(std::map< int, Column * > * barrels, std::vector< const char * > fileNames)
{
	if(nrOfPointSets > 1 && avgBarrels.size())
	{
		// create avg row direction vectors to analyze row/arc x-y std dev
		std::map< int, double * > rowDirections;
		double * aRowVec = new double[3];
		double * bRowVec = new double[3];
		double * cRowVec = new double[3];
		double * dRowVec = new double[3];
		double * eRowVec = new double[3];
		double * greekRowVec = new double[3];
		for(int ii = 0; ii < 3; ++ii)
		{
			aRowVec[ii] = avgBarrels[A4]->top[ii] - avgBarrels[A1]->top[ii];
			bRowVec[ii] = avgBarrels[B4]->top[ii] - avgBarrels[B1]->top[ii];
			cRowVec[ii] = avgBarrels[C4]->top[ii] - avgBarrels[C1]->top[ii];
			dRowVec[ii] = avgBarrels[D4]->top[ii] - avgBarrels[D1]->top[ii];
			eRowVec[ii] = avgBarrels[E4]->top[ii] - avgBarrels[E1]->top[ii];
			greekRowVec[ii] = avgBarrels[Delta]->top[ii] - avgBarrels[Alpha]->top[ii];
		}                                     
		normalize(aRowVec);                   
		normalize(bRowVec);                   
		normalize(cRowVec);                   
		normalize(dRowVec);
		normalize(eRowVec);
		normalize(greekRowVec);
		rowDirections.insert(std::pair< int, double * >(aRow, aRowVec));
		rowDirections.insert(std::pair< int, double * >(bRow, bRowVec));
		rowDirections.insert(std::pair< int, double * >(cRow, cRowVec));
		rowDirections.insert(std::pair< int, double * >(dRow, dRowVec));
		rowDirections.insert(std::pair< int, double * >(eRow, eRowVec));
		rowDirections.insert(std::pair< int, double * >(greekRow, greekRowVec));
		
		// analyze std dev in x-y plane and z -  all w.r.t. local barrel axis
		std::string AlignmentName(fileNames[0]);
		AlignmentName += "_avg_deviation.csv";
		std::ofstream AlignmentData(AlignmentName.c_str());
		AlignmentData << "Barrel field ID\tFile name" << std::endl;
		for(int ii = 0; ii < fileNames.size(); ++ii)
			AlignmentData << ii << "\t" << fileNames[ii] << std::endl;
		AlignmentData << std::endl;
		AlignmentData << "deviation from avg barrel field" << std::endl;
// 		AlignmentData << "Barrel\tavg height\ttop x-y variability\ttop z variability\ttop row variability\ttop arc variability\ttop total variability";
// 		AlignmentData << "\tbottom x-y variability\tbottom z variability\tbottom total variability" << std::endl;
		AlignmentData << "Barrel\tavg height\ttop z variability\ttop row variability\ttop arc variability\ttop total variability";
		AlignmentData << "\tbottom z variability\tbottom row variability\tbottom arc variability\tbottom total variability" << std::endl;
		
		std::list< int >::const_iterator labelIt;
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(avgBarrels.find(ID) != avgBarrels.end() && avgAxes.find(ID) != avgAxes.end())
			{
				gsl_matrix * mCovTop = gsl_matrix_calloc(3,3);
				gsl_matrix * mCovBottom = gsl_matrix_calloc(3,3);
				
				HomogeneousMatrixPointerType localCoords = getLocalCoordinatesForVariation(avgAxes[ID]);
				double hRowAxis[4], rowAxis[2];
				int row;
				if(aRow < ID && ID < bRow)
					row = aRow;
				else if(bRow < ID && ID < cRow)
					row = bRow;
				else if(cRow < ID && ID < dRow)
					row = cRow;
				else if(dRow < ID && ID < eRow)
					row = dRow;
				else if(eRow < ID && ID < greekRow)
					row = eRow;
				else if(greekRow < ID)
					row = greekRow;
				for(int ii = 0; ii < 3; ++ii)
					hRowAxis[ii] = rowDirections[row][ii];
				hRowAxis[3] = 1;
				localCoords->MultiplyPoint(hRowAxis, hRowAxis);
				rowAxis[0] = hRowAxis[0], rowAxis[1] = hRowAxis[1];
				double tmpRowAxNorm = sqrt(rowAxis[0]*rowAxis[0] + rowAxis[1]*rowAxis[1]);
				rowAxis[0] /= tmpRowAxNorm, rowAxis[1] /= tmpRowAxNorm;
				
// 				double hAxisDebug[4];
// 				for(int ii = 0; ii < 3; ++ii)
// 					hAxisDebug[ii] = avgAxes[ID][ii];
// 				hAxisDebug[3] = 1;
// 				localCoords->MultiplyPoint(hAxisDebug, hAxisDebug);
// 				std::cout << "barrel " << int2Labels[ID] << " local z axis: [" << hAxisDebug[0] << "," << hAxisDebug[1] << "," << hAxisDebug[2] << "]" << std::endl;
				
				double topPlaneDist = 0, bottomPlaneDist = 0, topZDist = 0, bottomZDist = 0;
				double topRowDist = 0, topArcDist = 0;
				double topTotalDist = 0, bottomTotalDist = 0;
				double avgCnt = 0;
				for(int ii = 0; ii < nrOfPointSets; ++ii)
					if(barrels[ii].find(ID) != barrels[ii].end())
					{
						++avgCnt;
						
						// slightly awkward version
						double hTopDisplacement[4], hBottomDisplacement[4];
						double tmpTopDist = 0, tmpBottomDist = 0;
						double tmpTopRowDist = 0, tmpTopArcDist = 0;
						for(int jj = 0; jj < 3; ++jj)
						{
							hTopDisplacement[jj] = barrels[ii][ID]->top[jj] - avgBarrels[ID]->top[jj];
							hBottomDisplacement[jj] = barrels[ii][ID]->bottom[jj] - avgBarrels[ID]->bottom[jj];
							tmpTopDist += hTopDisplacement[jj]*hTopDisplacement[jj];
							tmpBottomDist += hBottomDisplacement[jj]*hBottomDisplacement[jj];
						}
						hTopDisplacement[3] = hBottomDisplacement[3] = 1;
						localCoords->MultiplyPoint(hTopDisplacement, hTopDisplacement), localCoords->MultiplyPoint(hBottomDisplacement, hBottomDisplacement);
						
						topTotalDist += tmpTopDist;
						bottomTotalDist += tmpBottomDist;
// 						topPlaneDist += hTopDisplacement[0]*hTopDisplacement[0] + hTopDisplacement[1]*hTopDisplacement[1];
// 						bottomPlaneDist += hBottomDisplacement[0]*hBottomDisplacement[0] + hBottomDisplacement[1]*hBottomDisplacement[1];
// 						topZDist += hTopDisplacement[2]*hTopDisplacement[2];
// 						bottomZDist += hBottomDisplacement[2]*hBottomDisplacement[2];
// 						
// 						double tmpTopPlaneDist2 = hTopDisplacement[0]*hTopDisplacement[0] + hTopDisplacement[1]*hTopDisplacement[1];
// 						if(row == greekRow)
// 						{
// 							for(int jj = 0; jj < 2; ++jj)
// 								tmpTopArcDist += hTopDisplacement[jj]*rowAxis[jj];
// 							tmpTopRowDist = sqrt(tmpTopPlaneDist2 - tmpTopArcDist*tmpTopArcDist);
// 						}
// 						else
// 						{
// 							for(int jj = 0; jj < 2; ++jj)
// 								tmpTopRowDist += hTopDisplacement[jj]*rowAxis[jj];
// 							tmpTopArcDist = sqrt(tmpTopPlaneDist2 - tmpTopRowDist*tmpTopRowDist);
// 						}
// 						topRowDist += tmpTopRowDist*tmpTopRowDist;
// 						topArcDist += tmpTopArcDist*tmpTopArcDist;
						
						// covariance matrix version
						for(int jj = 0; jj < 3; ++jj)
							for(int kk = 0; kk < 3; ++kk)
							{
								double topDev = (barrels[ii][ID]->top[jj] - avgBarrels[ID]->top[jj])*(barrels[ii][ID]->top[kk] - avgBarrels[ID]->top[kk]);
								double bottomDev = (barrels[ii][ID]->bottom[jj] - avgBarrels[ID]->bottom[jj])*(barrels[ii][ID]->bottom[kk] - avgBarrels[ID]->bottom[kk]);
								double oldTopDev = gsl_matrix_get(mCovTop, jj, kk);
								double oldBottomDev = gsl_matrix_get(mCovBottom, jj, kk);
								gsl_matrix_set(mCovTop, jj, kk, oldTopDev+topDev);
								gsl_matrix_set(mCovBottom, jj, kk, oldBottomDev+bottomDev);
							}
					}
				if(avgCnt)
				{
					// slightly awkward version
					topTotalDist = sqrt(topTotalDist/avgCnt), bottomTotalDist = sqrt(bottomTotalDist/avgCnt);
// 					topZDist = sqrt(topZDist/avgCnt), bottomZDist = sqrt(bottomZDist/avgCnt);
// 					topPlaneDist = sqrt(topPlaneDist/avgCnt), bottomPlaneDist = sqrt(bottomPlaneDist/avgCnt);
// 					topRowDist = sqrt(topRowDist/avgCnt), topArcDist = sqrt(topArcDist/avgCnt);
/*					
					double avgHeight = L2Distance3D(avgBarrels[ID]->top, avgBarrels[ID]->bottom);
					AlignmentData << int2Labels[ID] << "\t" << avgHeight << "\t" << topPlaneDist << "\t" << topZDist << "\t" << topRowDist << "\t" << topArcDist << "\t" << topTotalDist;
					AlignmentData << "\t" << bottomPlaneDist << "\t" << bottomZDist << "\t" << bottomTotalDist << std::endl;*/
					
					// covariance matrix version
					gsl_matrix_scale(mCovTop, 1/double(avgCnt-1));
					gsl_matrix_scale(mCovBottom, 1/double(avgCnt-1));
					// diagonalization
					gsl_vector * eigenValsTop = gsl_vector_alloc(3);
					gsl_matrix * eigenVecsTop = gsl_matrix_alloc(3, 3);
					gsl_eigen_symmv_workspace * wTop = gsl_eigen_symmv_alloc(3);
					gsl_vector * eigenValsBottom = gsl_vector_alloc(3);
					gsl_matrix * eigenVecsBottom = gsl_matrix_alloc(3, 3);
					gsl_eigen_symmv_workspace * wBottom = gsl_eigen_symmv_alloc(3);
					int status = gsl_eigen_symmv(mCovTop, eigenValsTop, eigenVecsTop, wTop);
					status = gsl_eigen_symmv(mCovBottom, eigenValsBottom, eigenVecsBottom, wBottom);
					// sort eigenvalues:
					// compute row/arc/zaxis closest
					// to corresponding eigenvector
					int topRow = -1, topArc = -1, topZ = -1;
					int bottomRow = -1, bottomArc = -1, bottomZ = -1;
					double topEVec[9];
					double bottomEVec[9];
					for(int jj = 0; jj < 3; ++jj)
					{
						topEVec[jj] = gsl_matrix_get(eigenVecsTop, jj, 0);
						topEVec[jj+3] = gsl_matrix_get(eigenVecsTop, jj, 1);
						topEVec[jj+6] = gsl_matrix_get(eigenVecsTop, jj, 2);
						bottomEVec[jj] = gsl_matrix_get(eigenVecsBottom, jj, 0);
						bottomEVec[jj+3] = gsl_matrix_get(eigenVecsBottom, jj, 1);
						bottomEVec[jj+6] = gsl_matrix_get(eigenVecsBottom, jj, 2);
					}
					double topRowScore = 0, topZScore = 0;
					double bottomRowScore = 0, bottomZScore = 0;
					for(int jj = 0; jj < 3; ++jj)
					{
						double tmpTopRow, tmpTopZ;
						double tmpBottomRow, tmpBottomZ;
						tmpTopRow = fabs(vtkMath::Dot(rowDirections[row], topEVec+3*jj));
						tmpTopZ = fabs(vtkMath::Dot(avgAxes[ID], topEVec+3*jj));
						tmpBottomRow = fabs(vtkMath::Dot(rowDirections[row], bottomEVec+3*jj));
						tmpBottomZ = fabs(vtkMath::Dot(avgAxes[ID], bottomEVec+3*jj));
						if(tmpTopRow > topRowScore)
						{
							topRow = jj;
							topRowScore = tmpTopRow;
						}
						if(tmpTopZ > topZScore)
						{
							topZ = jj;
							topZScore = tmpTopZ;
						}
						if(tmpBottomRow > bottomRowScore)
						{
							bottomRow = jj;
							bottomRowScore = tmpBottomRow;
						}
						if(tmpBottomZ > bottomZScore)
						{
							bottomZ = jj;
							bottomZScore = tmpBottomZ;
						}
					}
					if(topRow == topZ)
					{
						if(bottomRow == bottomZ)
						{
							std::cout << "Error! Could not determine row/arc/z closest to evec!" << std::endl;
							std::cout << "topRow = " << topRow << " - topZ = " << topZ << std::endl;
							std::cout << "bottomRow = " << bottomRow << " - bottomZ = " << bottomZ << std::endl;
							std::cout << "Cannot determine variability of avg barrel " << int2Labels[ID] << std::endl;
							continue;
						}
						else
						{
							topRow = bottomRow;
							topZ = bottomZ;
						}
					}
					if(bottomRow == bottomZ)
					{
						if(topRow == topZ)
						{
							std::cout << "Error! Could not determine row/arc/z closest to evec!" << std::endl;
							std::cout << "topRow = " << topRow << " - topZ = " << topZ << std::endl;
							std::cout << "bottomRow = " << bottomRow << " - bottomZ = " << bottomZ << std::endl;
							std::cout << "Cannot determine variability of avg barrel " << int2Labels[ID] << std::endl;
							continue;
						}
						else
						{
							bottomRow = topRow;
							bottomZ = topZ;
						}
					}
					if(row == greekRow)
					{
						// switch row/arc
						topArc = 3 - topRow - topZ;
						bottomArc = 3 - bottomRow - bottomZ;
						int tmp = topRow;
						topRow = topArc;
						topArc = tmp;
						tmp = bottomRow;
						bottomRow = bottomArc;
						bottomArc = tmp;
						
					}
					else
					{
						topArc = 3 - topRow - topZ;
						bottomArc = 3 - bottomRow - bottomZ;
					}
					
					double avgHeight = L2Distance3D(avgBarrels[ID]->top, avgBarrels[ID]->bottom);
					AlignmentData << int2Labels[ID] << "\t" << avgHeight << "\t";
					AlignmentData << sqrt(gsl_vector_get(eigenValsTop, topZ)) << "\t" << sqrt(gsl_vector_get(eigenValsTop, topRow)) << "\t";
					AlignmentData << sqrt(gsl_vector_get(eigenValsTop, topArc)) << "\t" << topTotalDist << "\t";
					AlignmentData << sqrt(gsl_vector_get(eigenValsBottom, bottomZ)) << "\t" << sqrt(gsl_vector_get(eigenValsBottom, bottomRow)) << "\t";
					AlignmentData << sqrt(gsl_vector_get(eigenValsBottom, bottomArc)) << "\t" << bottomTotalDist << std::endl;
					
					gsl_matrix_free(eigenVecsTop);
					gsl_vector_free(eigenValsTop);
					gsl_eigen_symmv_free(wTop);
					gsl_matrix_free(eigenVecsBottom);
					gsl_vector_free(eigenValsBottom);
					gsl_eigen_symmv_free(wBottom);
				}
				
				gsl_matrix_free(mCovTop);
				gsl_matrix_free(mCovBottom);
			}
		}
		
		AlignmentData.close();
	}
};

void Registration::writeAvgBarrelData(std::map< int, Column * > * matchedBarrels, std::vector< const char * > fileNames)
{
	std::string AvgBarrelParameters(fileNames[0]);
	AvgBarrelParameters += "_avg_barrel_field_parameters.csv";
	std::ofstream AvgBarrelData(AvgBarrelParameters.c_str());
	AvgBarrelData << "Barrel\tcenter\tz axis unit vec\tradius" << std::endl;
	std::list< int >::const_iterator labelIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		#ifndef PIA_WM_REG
		if(avgAxes.find(ID) != avgAxes.end() && avgCenters.find(ID) != avgCenters.end() && avgBarrels.find(ID) != avgBarrels.end())
		#elif defined PIA_WM_REG
		if(avgAxes.find(ID) != avgAxes.end() && avgCenters.find(ID) != avgCenters.end() && avgColumns.find(ID) != avgColumns.end())
		#endif
		{
			// create avg contour: use circle as 1st approximation
			// better contour possible w/ voronoi region growth???
			PolyDataPointerType barrelContours = PolyDataPointerType::New();
			PolyDataPointerType colContours = PolyDataPointerType::New();
			PointsPointerType barrelPoints = PointsPointerType::New();
			PointsPointerType colPoints = PointsPointerType::New();
			PolygonPointerType barrelTopPoly = PolygonPointerType::New();
			PolygonPointerType barrelBottomPoly = PolygonPointerType::New();
			PolygonPointerType colTopPoly = PolygonPointerType::New();
			PolygonPointerType colBottomPoly = PolygonPointerType::New();
			barrelContours->Allocate(1);
			colContours->Allocate(1);
			barrelPoints->SetDataTypeToFloat();
			colPoints->SetDataTypeToFloat();
			barrelTopPoly->GetPointIds()->SetNumberOfIds(36), barrelBottomPoly->GetPointIds()->SetNumberOfIds(36);
			colTopPoly->GetPointIds()->SetNumberOfIds(36), colBottomPoly->GetPointIds()->SetNumberOfIds(36);
			
			//for illustration
			PolyDataPointerType centerContours = PolyDataPointerType::New();
			PointsPointerType centerPoints = PointsPointerType::New();
			PolygonPointerType centerTopPoly = PolygonPointerType::New();
			PolygonPointerType centerBottomPoly = PolygonPointerType::New();
			centerContours->Allocate(1);
			centerPoints->SetDataTypeToFloat();
			centerTopPoly->GetPointIds()->SetNumberOfIds(36), centerBottomPoly->GetPointIds()->SetNumberOfIds(36);
			std::vector< double * > centerPtVec;
			
			double radius = sqrt(avgBarrelArea[ID]/PI);
			HomogeneousMatrixPointerType mat = transformToBarrelCoordinates(avgAxes[ID]);
			std::vector< double * > contourPtVec;
			for(int ii = 0; ii < 36; ++ii)
			{
				double * pt = new double[3];
				double hPt[4];
				pt[0] = radius*cos(ii*PI/18), pt[1] = radius*sin(ii*PI/18), pt[2] = 0;
				hPt[0] = pt[0], hPt[1] = pt[1], hPt[2] = pt[2], hPt[3] = 1;
				mat->MultiplyPoint(hPt, hPt);
				pt[0] = hPt[0], pt[1] = hPt[1], pt[2] = hPt[2];
				
				double * pt2 = new double[3];
				double hPt2[4];
				pt2[0] = 5*cos(ii*PI/18), pt2[1] = 5*sin(ii*PI/18), pt2[2] = 0;
				hPt2[0] = pt2[0], hPt2[1] = pt2[1], hPt2[2] = pt2[2], hPt2[3] = 1;
				mat->MultiplyPoint(hPt2, hPt2);
				pt2[0] = hPt2[0], pt2[1] = hPt2[1], pt2[2] = hPt2[2];
				
				double * barrelTopPt = new double[3], * colTopPt = new double[3];
				double * centerTopPt = new double[3];
				for(int jj = 0; jj < 3; ++jj)
				{
					barrelTopPt[jj] = avgBarrels[ID]->top[jj] + pt[jj];
					colTopPt[jj] = avgColumns[ID]->top[jj] + pt[jj];
					centerTopPt[jj] = avgBarrels[ID]->top[jj] + pt2[jj];
				}
				contourPtVec.push_back(pt);
				barrelPoints->InsertNextPoint(barrelTopPt);
				colPoints->InsertNextPoint(colTopPt);
				barrelTopPoly->GetPointIds()->SetId(ii, ii);
				colTopPoly->GetPointIds()->SetId(ii, ii);
				
				centerPtVec.push_back(pt2);
				centerPoints->InsertNextPoint(centerTopPt);
				centerTopPoly->GetPointIds()->SetId(ii, ii);
			}
			for(int ii = 0; ii < 36; ++ii)
			{
				double * barrelBottomPt = new double[3], * colBottomPt = new double[3];
				double * centerBottomPt = new double[3];
				for(int jj = 0; jj < 3; ++jj)
				{
					barrelBottomPt[jj] = avgBarrels[ID]->bottom[jj] + contourPtVec[ii][jj];
					colBottomPt[jj] = avgColumns[ID]->bottom[jj] + contourPtVec[ii][jj];
					centerBottomPt[jj] = avgBarrels[ID]->bottom[jj] + centerPtVec[ii][jj];
				}
				barrelPoints->InsertNextPoint(barrelBottomPt);
				colPoints->InsertNextPoint(colBottomPt);
				barrelBottomPoly->GetPointIds()->SetId(ii, ii+36);
				colBottomPoly->GetPointIds()->SetId(ii, ii+36);
				
				centerPoints->InsertNextPoint(centerBottomPt);
				centerBottomPoly->GetPointIds()->SetId(ii, ii+36);
				
				delete contourPtVec[ii];
				delete centerPtVec[ii];
			}
			contourPtVec.clear();
			barrelContours->InsertNextCell(barrelTopPoly->GetCellType(), barrelTopPoly->GetPointIds());
			barrelContours->InsertNextCell(barrelBottomPoly->GetCellType(), barrelBottomPoly->GetPointIds());
			barrelContours->SetPoints(barrelPoints);
			colContours->InsertNextCell(colTopPoly->GetCellType(), colTopPoly->GetPointIds());
			colContours->InsertNextCell(colBottomPoly->GetCellType(), colBottomPoly->GetPointIds());
			colContours->SetPoints(colPoints);
			avgBarrels[ID]->contours = barrelContours;
			avgColumns[ID]->contours = colContours;
			
			spatialGraph->addPolyDataObject(avgBarrels[ID]->contours, ID);
			spatialGraph->addPolyDataObject(avgColumns[ID]->contours, ID);
			
			centerContours->InsertNextCell(centerTopPoly->GetCellType(), centerTopPoly->GetPointIds());
			centerContours->InsertNextCell(centerBottomPoly->GetCellType(), centerBottomPoly->GetPointIds());
			centerContours->SetPoints(centerPoints);
			spatialGraph->addPolyDataObject(centerContours, ID);
			centerPtVec.clear();
			
			AvgBarrelData << int2Labels[ID] << "\t" << avgCenters[ID][0] << "," << avgCenters[ID][1] << "," << avgCenters[ID][2] << "\t";
			AvgBarrelData << avgAxes[ID][0] << "," << avgAxes[ID][1] << "," << avgAxes[ID][2] << "\t" << radius << std::endl;
			
// 			Vertex * newVert1 = new Vertex(avgColumns[ID]->top, ID);
// 			Vertex * newVert2 = new Vertex(avgColumns[ID]->bottom, ID);
// 			spatialGraph->addVertex(newVert1), spatialGraph->addVertex(newVert2);
// 			int connectivity[2];
// 			connectivity[0] = spatialGraph->getNumberOfVertices()-2, connectivity[1] = spatialGraph->getNumberOfVertices()-1;
// 			std::list< double * > edgePtList;
// 			edgePtList.push_back(avgColumns[ID]->top), edgePtList.push_back(avgColumns[ID]->bottom);
// 			Edge * newEdge = new Edge(connectivity, 2, ID, edgePtList);
// 			spatialGraph->addEdge(newEdge);
		}
	}
	AvgBarrelData.close();
	
	for(int n = 0; n < nrOfPointSets; ++n)
	{
		std::string TransformParameters(fileNames[n]);
		TransformParameters += "_transform.log";
		std::ofstream TransformData(TransformParameters.c_str());
		TransformData << "# Transformation matrix" << std::endl;
		TransformData << "# for use with Amira" << std::endl;
		TransformData << "# Format: column 1 column 2 column 3 column 4" << std::endl;
		for(int jj = 0; jj < 4; ++jj)
			for(int ii = 0; ii < 4; ++ii)
				TransformData << transformVec[n]->GetMatrix()->GetElement(ii, jj) << " ";
			TransformData << std::endl;
		TransformData.close();
		
		std::string BarrelParameters(fileNames[n]);
		BarrelParameters += "_reg_barrel_centers.log";
		std::ofstream BarrelParameterData(BarrelParameters.c_str());
		BarrelParameterData << "# Registered barrel center coordinates" << std::endl;
		BarrelParameterData << "Barrel\tcenter\tz axis unit vec\tangle with std axis" << std::endl;
//  		BarrelParameterData << "Barrel\tcenter\tz axis unit vec" << std::endl;
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(matchedBarrels[n].find(ID) != matchedBarrels[n].end())
			{
				double center[3], axis[3], angle = 0;
				for(int ii = 0; ii < 3; ++ii)
				{
					center[ii] = 0.5*(matchedBarrels[n][ID]->top[ii] + matchedBarrels[n][ID]->bottom[ii]);
					axis[ii] = matchedBarrels[n][ID]->top[ii] - matchedBarrels[n][ID]->bottom[ii];
				}
				normalize(axis);
				for(int ii = 0; ii < 3; ++ii)
					angle += axis[ii]*avgAxes[ID][ii];
				angle = std::acos(angle)*180/PI;
				BarrelParameterData << int2Labels[ID] << "\t" << center[0] << "\t" << center[1] << "\t" << center[2] << "\t";
				BarrelParameterData << axis[0] << "\t" << axis[1] << "\t" << axis[2] << "\t" << angle << std::endl;
			}
		}
		BarrelParameterData.close();
	}
};

/******************************************************************************/
/*create average surfaces by triangulating isosurface points along individual */
/*barrel axes                                                                 */
/******************************************************************************/
void Registration::computeAverageSurfaces(const char * surfaceFilename)
{
	if(avgBarrels.size() && avgColumns.size())
	{
		std::cout << "Computing initial triangulation of average barrel field surfaces..." << std::endl;
		Delaunay2DFilterPointerType delTriangulation1 = Delaunay2DFilterPointerType::New();
		Delaunay2DFilterPointerType delTriangulation2 = Delaunay2DFilterPointerType::New();
		Delaunay2DFilterPointerType delTriangulation3 = Delaunay2DFilterPointerType::New();
		Delaunay2DFilterPointerType delTriangulation4 = Delaunay2DFilterPointerType::New();
		PointsPointerType piaPoints = PointsPointerType::New();
		PointsPointerType L4UpperPoints = PointsPointerType::New();
		PointsPointerType L4LowerPoints = PointsPointerType::New();
		PointsPointerType wmPoints = PointsPointerType::New();
		PolyDataPointerType piaSampling = PolyDataPointerType::New(); // Delaunay filter needs DataOpbject as input...
		PolyDataPointerType L4UpperSampling = PolyDataPointerType::New();
		PolyDataPointerType L4LowerSampling = PolyDataPointerType::New();
		PolyDataPointerType wmSampling = PolyDataPointerType::New();
		piaPoints->SetDataTypeToFloat(), L4UpperPoints->SetDataTypeToFloat(), L4LowerPoints->SetDataTypeToFloat(), wmPoints->SetDataTypeToFloat();
		piaSampling->Allocate(1), L4UpperSampling->Allocate(1), L4LowerSampling->Allocate(1), wmSampling->Allocate(1);
		std::list< int >::const_iterator labelIt;
		
		#ifndef CONVEX_HULL_SURFS
		double barrelFieldCOM[] = {0,0,0};
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(avgBarrels.find(ID) != avgBarrels.end())
			{
				for(int ii = 0; ii < 3; ++ii)
				{
					barrelFieldCOM[ii] += avgBarrels[ID]->top[ii];
					barrelFieldCOM[ii] += avgBarrels[ID]->bottom[ii];
				}
			}
		}
		for(int ii = 0; ii < 3; ++ii)
			barrelFieldCOM[ii] /= double(2*avgBarrels.size());
		#endif
		
		#ifdef CONVEX_HULL_SURFS
		// set up convex hull pts
		ConvexHull2DFilterPointerType piaProjectedHull = ConvexHull2DFilterPointerType::New();
		PointsPointerType piaBorderBarrelPts = PointsPointerType::New();
		piaProjectedHull->SetDataTypeToFloat();
		piaBorderBarrelPts->SetDataTypeToFloat();
		ConvexHull2DFilterPointerType L4UProjectedHull = ConvexHull2DFilterPointerType::New();
		PointsPointerType L4UBorderBarrelPts = PointsPointerType::New();
		L4UProjectedHull->SetDataTypeToFloat();
		L4UBorderBarrelPts->SetDataTypeToFloat();
		ConvexHull2DFilterPointerType L4LProjectedHull = ConvexHull2DFilterPointerType::New();
		PointsPointerType L4LBorderBarrelPts = PointsPointerType::New();
		L4LProjectedHull->SetDataTypeToFloat();
		L4LBorderBarrelPts->SetDataTypeToFloat();
		ConvexHull2DFilterPointerType wmProjectedHull = ConvexHull2DFilterPointerType::New();
		PointsPointerType wmBorderBarrelPts = PointsPointerType::New();
		wmProjectedHull->SetDataTypeToFloat();
		wmBorderBarrelPts->SetDataTypeToFloat();
		#endif
		
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(avgColumns.find(ID) != avgColumns.end())
			{
				piaPoints->InsertNextPoint(avgColumns[ID]->top);
				L4UpperPoints->InsertNextPoint(avgBarrels[ID]->top);
				L4LowerPoints->InsertNextPoint(avgBarrels[ID]->bottom);
				wmPoints->InsertNextPoint(avgColumns[ID]->bottom);
				
				#ifdef CONVEX_HULL_SURFS
				if(std::find(borderBarrels.begin(), borderBarrels.end(), ID) != borderBarrels.end())
				{
					for(int ii = 0; ii < avgColumns[ID]->contours->GetCell(0)->GetNumberOfPoints(); ++ii)
					{
						double piaPt[3];
						avgColumns[ID]->contours->GetCell(0)->GetPoints()->GetPoint(ii, piaPt);
						piaProjectedHull->InsertNextPoint(piaPt);
						piaBorderBarrelPts->InsertNextPoint(piaPt);
					}
					for(int ii = 0; ii < avgBarrels[ID]->contours->GetCell(0)->GetNumberOfPoints(); ++ii)
					{
						double L4UPt[3];
						avgBarrels[ID]->contours->GetCell(0)->GetPoints()->GetPoint(ii, L4UPt);
						L4UProjectedHull->InsertNextPoint(L4UPt);
						L4UBorderBarrelPts->InsertNextPoint(L4UPt);
					}
					for(int ii = 0; ii < avgBarrels[ID]->contours->GetCell(1)->GetNumberOfPoints(); ++ii)
					{
						double L4LPt[3];
						avgBarrels[ID]->contours->GetCell(1)->GetPoints()->GetPoint(ii, L4LPt);
						L4LProjectedHull->InsertNextPoint(L4LPt);
						L4LBorderBarrelPts->InsertNextPoint(L4LPt);
					}
					for(int ii = 0; ii < avgColumns[ID]->contours->GetCell(1)->GetNumberOfPoints(); ++ii)
					{
						double wmPt[3];
						avgColumns[ID]->contours->GetCell(1)->GetPoints()->GetPoint(ii, wmPt);
						wmProjectedHull->InsertNextPoint(wmPt);
						wmBorderBarrelPts->InsertNextPoint(wmPt);
					}
				}
				#endif
				
				#ifndef CONVEX_HULL_SURFS
				if(std::find(borderBarrels.begin(), borderBarrels.end(), ID) != borderBarrels.end())
					piaPoints->InsertNextPoint(findExtremePoint(avgColumns[ID]->contours->GetCell(0)->GetPoints(), avgColumns[ID]->top, barrelFieldCOM));
				if(std::find(borderBarrels.begin(), borderBarrels.end(), ID) != borderBarrels.end())
					L4UpperPoints->InsertNextPoint(findExtremePoint(avgBarrels[ID]->contours->GetCell(0)->GetPoints(), avgBarrels[ID]->top, barrelFieldCOM));
				if(std::find(borderBarrels.begin(), borderBarrels.end(), ID) != borderBarrels.end())
					L4LowerPoints->InsertNextPoint(findExtremePoint(avgBarrels[ID]->contours->GetCell(1)->GetPoints(), avgBarrels[ID]->bottom, barrelFieldCOM));
				if(std::find(borderBarrels.begin(), borderBarrels.end(), ID) != borderBarrels.end())
					wmPoints->InsertNextPoint(findExtremePoint(avgColumns[ID]->contours->GetCell(1)->GetPoints(), avgColumns[ID]->bottom, barrelFieldCOM));
				#endif
			}
		}
		
		#ifdef CONVEX_HULL_SURFS
		// Pia
		piaProjectedHull->Update();
		int nrPiaHullPts = piaProjectedHull->GetSizeCCWHullZ();
		double * piaHullPts = new double[2*nrPiaHullPts];
		piaProjectedHull->GetCCWHullZ(piaHullPts, nrPiaHullPts);
		for(int ii = 0; ii < nrPiaHullPts; ++ii)
		{
			double minDist = 1E09;
			int minID = -1;
			for(int jj = 0; jj < piaBorderBarrelPts->GetNumberOfPoints(); ++jj)
			{
				double pt2[3], dist2D;
				piaBorderBarrelPts->GetPoint(jj, pt2);
				dist2D = L2Distance2D(pt2, piaHullPts+2*ii);
				if(dist2D < minDist)
				{
					minDist = dist2D;
					minID = jj;
				}
			}
			
			if(minID >= 0)
			{
				double thisHullPt[3];
				piaBorderBarrelPts->GetPoint(minID, thisHullPt);
				piaPoints->InsertNextPoint(thisHullPt);
			}
			else
				std::cout << "Oops! Something is wrong calculating the convex hull of Pia..." << std::endl;
		}
		// L4U
		L4UProjectedHull->Update();
		int nrL4UHullPts = L4UProjectedHull->GetSizeCCWHullZ();
		double * L4UHullPts = new double[2*nrL4UHullPts];
		L4UProjectedHull->GetCCWHullZ(L4UHullPts, nrL4UHullPts);
		for(int ii = 0; ii < nrL4UHullPts; ++ii)
		{
			double minDist = 1E09;
			int minID = -1;
			for(int jj = 0; jj < L4UBorderBarrelPts->GetNumberOfPoints(); ++jj)
			{
				double pt2[3], dist2D;
				L4UBorderBarrelPts->GetPoint(jj, pt2);
				dist2D = L2Distance2D(pt2, L4UHullPts+2*ii);
				if(dist2D < minDist)
				{
					minDist = dist2D;
					minID = jj;
				}
			}
			
			if(minID >= 0)
			{
				double thisHullPt[3];
				L4UBorderBarrelPts->GetPoint(minID, thisHullPt);
				L4UpperPoints->InsertNextPoint(thisHullPt);
			}
			else
				std::cout << "Oops! Something is wrong calculating the convex hull of L4U..." << std::endl;
		}
		// L4L
		L4LProjectedHull->Update();
		int nrL4LHullPts = L4LProjectedHull->GetSizeCCWHullZ();
		double * L4LHullPts = new double[2*nrL4LHullPts];
		L4LProjectedHull->GetCCWHullZ(L4LHullPts, nrL4LHullPts);
		for(int ii = 0; ii < nrL4LHullPts; ++ii)
		{
			double minDist = 1E09;
			int minID = -1;
			for(int jj = 0; jj < L4LBorderBarrelPts->GetNumberOfPoints(); ++jj)
			{
				double pt2[3], dist2D;
				L4LBorderBarrelPts->GetPoint(jj, pt2);
				dist2D = L2Distance2D(pt2, L4LHullPts+2*ii);
				if(dist2D < minDist)
				{
					minDist = dist2D;
					minID = jj;
				}
			}
			
			if(minID >= 0)
			{
				double thisHullPt[3];
				L4LBorderBarrelPts->GetPoint(minID, thisHullPt);
				L4LowerPoints->InsertNextPoint(thisHullPt);
			}
			else
				std::cout << "Oops! Something is wrong calculating the convex hull of L4L..." << std::endl;
		}
		// WM
		wmProjectedHull->Update();
		int nrwmHullPts = wmProjectedHull->GetSizeCCWHullZ();
		double * wmHullPts = new double[2*nrwmHullPts];
		wmProjectedHull->GetCCWHullZ(wmHullPts, nrwmHullPts);
		for(int ii = 0; ii < nrwmHullPts; ++ii)
		{
			double minDist = 1E09;
			int minID = -1;
			for(int jj = 0; jj < wmBorderBarrelPts->GetNumberOfPoints(); ++jj)
			{
				double pt2[3], dist2D;
				wmBorderBarrelPts->GetPoint(jj, pt2);
				dist2D = L2Distance2D(pt2, wmHullPts+2*ii);
				if(dist2D < minDist)
				{
					minDist = dist2D;
					minID = jj;
				}
			}
			
			if(minID >= 0)
			{
				double thisHullPt[3];
				wmBorderBarrelPts->GetPoint(minID, thisHullPt);
				wmPoints->InsertNextPoint(thisHullPt);
			}
			else
				std::cout << "Oops! Something is wrong calculating the convex hull of WM..." << std::endl;
		}
		#endif
		
		piaSampling->SetPoints(piaPoints);
		delTriangulation1->SetInput(piaSampling);
		delTriangulation1->Update();
		PolyDataPointerType piaTriangulation = delTriangulation1->GetOutput();
		L4UpperSampling->SetPoints(L4UpperPoints);
		delTriangulation2->SetInput(L4UpperSampling);
		delTriangulation2->Update();
		PolyDataPointerType L4UpperTriangulation = delTriangulation2->GetOutput();
		L4LowerSampling->SetPoints(L4LowerPoints);
		delTriangulation3->SetInput(L4LowerSampling);
		delTriangulation3->Update();
		PolyDataPointerType L4LowerTriangulation = delTriangulation3->GetOutput();
		wmSampling->SetPoints(wmPoints);
		delTriangulation4->SetInput(wmSampling);
		delTriangulation4->Update();
		PolyDataPointerType wmTriangulation = delTriangulation4->GetOutput();
		
		std::flush(std::cout << "Refining mesh of average barrel field surfaces..." << std::endl);
		MeshRefinementFilterPointerType piaMeshRefinement = MeshRefinementFilterPointerType::New();
		piaMeshRefinement->SetNumberOfSubdivisions(1);
		piaMeshRefinement->SetInput(piaTriangulation);
		piaMeshRefinement->Update();
		PolyDataPointerType piaMeshFine = piaMeshRefinement->GetOutput();
		MeshRefinementFilterPointerType L4UpperMeshRefinement = MeshRefinementFilterPointerType::New();
		L4UpperMeshRefinement->SetNumberOfSubdivisions(1);
		L4UpperMeshRefinement->SetInput(L4UpperTriangulation);
		L4UpperMeshRefinement->Update();
		PolyDataPointerType L4UpperMeshFine = L4UpperMeshRefinement->GetOutput();
		MeshRefinementFilterPointerType L4LowerMeshRefinement = MeshRefinementFilterPointerType::New();
		L4LowerMeshRefinement->SetNumberOfSubdivisions(1);
		L4LowerMeshRefinement->SetInput(L4LowerTriangulation);
		L4LowerMeshRefinement->Update();
		PolyDataPointerType L4LowerMeshFine = L4LowerMeshRefinement->GetOutput();
		MeshRefinementFilterPointerType WMMeshRefinement = MeshRefinementFilterPointerType::New();
		WMMeshRefinement->SetNumberOfSubdivisions(1);
		WMMeshRefinement->SetInput(wmTriangulation);
		WMMeshRefinement->Update();
		PolyDataPointerType wmMeshFine = WMMeshRefinement->GetOutput();
		
		// comment following lines for variability calculation
		// otherwise, produces too much output...
		#ifndef REG_ACCURACY
		std::string piaFilename(surfaceFilename);
		piaFilename += "_pia";
		Reader * piaWriter = new Reader(piaFilename.c_str(), piaFilename.c_str());
		piaWriter->writeAmiraSurfaceFile(piaMeshFine);
		delete piaWriter;
		std::string L4UpperFilename(surfaceFilename);
		L4UpperFilename += "_L4Upper";
		Reader * L4UpperWriter = new Reader(L4UpperFilename.c_str(), L4UpperFilename.c_str());
		L4UpperWriter->writeAmiraSurfaceFile(L4UpperMeshFine);
		delete L4UpperWriter;
		std::string L4LowerFilename(surfaceFilename);
		L4LowerFilename += "_L4Lower";
		Reader * L4LowerWriter = new Reader(L4LowerFilename.c_str(), L4LowerFilename.c_str());
		L4LowerWriter->writeAmiraSurfaceFile(L4LowerMeshFine);
		delete L4LowerWriter;
		std::string WMFilename(surfaceFilename);
		WMFilename += "_WM";
		Reader * wmWriter = new Reader(WMFilename.c_str(), WMFilename.c_str());
		wmWriter->writeAmiraSurfaceFile(wmMeshFine);
		delete wmWriter;
		#endif
		
		// calculate and write all cell type boundary surfaces
		std::list< double >::const_iterator ratioDepthIt;
		int surfCnt = 0;
		for(ratioDepthIt = cellTypeRatioDepthsSupra.begin(); ratioDepthIt != cellTypeRatioDepthsSupra.end(); ++ratioDepthIt, ++surfCnt)
		{
			Delaunay2DFilterPointerType ratioTriangulation1 = Delaunay2DFilterPointerType::New();
			PointsPointerType ratioPoints = PointsPointerType::New();
			PolyDataPointerType ratioSampling = PolyDataPointerType::New(); // Delaunay filter needs DataOpbject as input...
			ratioPoints->SetDataTypeToFloat();
			ratioSampling->Allocate(1);
			
			#ifndef CONVEX_HULL_SURFS
			for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
			{
				int ID = *labelIt;
				if(avgColumns.find(ID) != avgColumns.end())
				{
					double * newPt = new double[3];
					double barrelTop[3];
					for(int ii = 0; ii < 3; ++ii)
						barrelTop[ii] = avgBarrels[ID]->top[ii];
					double dist = avgTopDist[ID];
					for(int ii = 0; ii < 3; ++ii)
						newPt[ii] = barrelTop[ii] + *ratioDepthIt*dist*avgAxes[ID][ii];
					ratioPoints->InsertNextPoint(newPt);
					if(std::find(borderBarrels.begin(), borderBarrels.end(), ID) != borderBarrels.end())
					{
						PolyDataPointerType tmpContour = PolyDataPointerType::New();
						tmpContour->Allocate(1);
						tmpContour->DeepCopy(avgColumns[ID]->contours);
						Column * tmpCol = new Column;
						tmpCol->contours = tmpContour;
						double shiftVec[3];
						for(int ii = 0; ii < 3; ++ii)
							shiftVec[ii] = newPt[ii] - avgColumns[ID]->top[ii];
						tmpCol->translateColumn(shiftVec);
						ratioPoints->InsertNextPoint(findExtremePoint(tmpCol->contours->GetCell(0)->GetPoints(), newPt, barrelFieldCOM));
						delete tmpCol;
					}
				}
			}
			#endif
			
			#ifdef CONVEX_HULL_SURFS
			for(int ii = 0; ii < L4UpperPoints->GetNumberOfPoints(); ++ii)
			{
				double pt[3];
				L4UpperPoints->GetPoint(ii, pt);
				// brute force way finding the corresponding columns
				// of all sampled L4L pts:
				double minDist = 1E09;
				int minID = -1;
				for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
				{
					int ID = *labelIt;
					if(avgBarrels.find(ID) != avgBarrels.end())
					{
						double centerDist = L2Distance3D(pt, avgBarrels[ID]->bottom);
						if(centerDist < minDist)
						{
							minDist = centerDist;
							minID = ID;
						}
						for(int n = 0; n < avgBarrels[ID]->contours->GetCell(1)->GetNumberOfPoints(); ++n)
						{
							double pt2[3];
							avgBarrels[ID]->contours->GetCell(1)->GetPoints()->GetPoint(n, pt2);
							double edgeDist =  L2Distance3D(pt, pt2);
							if(edgeDist < minDist)
							{
								minDist = edgeDist;
								minID = ID;
							}
						}
					}
				}
				if(minID >= 0)
				{
					int ID = minID;
					double barrelTop[3], newPt[3];
					for(int jj = 0; jj < 3; ++jj)
						barrelTop[jj] = pt[jj];
					double dist = avgTopDist[ID];
					for(int jj = 0; jj < 3; ++jj)
						newPt[jj] = barrelTop[jj] + *ratioDepthIt*dist*avgAxes[ID][jj];
					ratioPoints->InsertNextPoint(newPt);
				}
				else
					std::cout << "Oops! Something is wrong calculating the convex hull of the ratio surfs..." << std::endl;
			}
			#endif
			
			ratioSampling->SetPoints(ratioPoints);
			ratioTriangulation1->SetInput(ratioSampling);
			ratioTriangulation1->Update();
			PolyDataPointerType ratioTriangulation = ratioTriangulation1->GetOutput();
			MeshRefinementFilterPointerType ratioMeshRefinement = MeshRefinementFilterPointerType::New();
			ratioMeshRefinement->SetNumberOfSubdivisions(1);
			ratioMeshRefinement->SetInput(ratioTriangulation);
			ratioMeshRefinement->Update();
			PolyDataPointerType ratioMeshFine = ratioMeshRefinement->GetOutput();
			std::string ratioFilename(surfaceFilename);
			ratioFilename += "_celltype_border_%02d";
			// comment remaining lines for variability calculation
			// otherwise, produces too much output...
			#ifndef REG_ACCURACY
			char * borderFilename = new char[256];
			sprintf(borderFilename, ratioFilename.c_str(), surfCnt);
			Reader * ratioWriter = new Reader(borderFilename, borderFilename);
			ratioWriter->writeAmiraSurfaceFile(ratioMeshFine);
			delete ratioWriter;
			#endif
		}
		for(ratioDepthIt = cellTypeRatioDepthsGran.begin(); ratioDepthIt != cellTypeRatioDepthsGran.end(); ++ratioDepthIt, ++surfCnt)
		{
			Delaunay2DFilterPointerType ratioTriangulation1 = Delaunay2DFilterPointerType::New();
			PointsPointerType ratioPoints = PointsPointerType::New();
			PolyDataPointerType ratioSampling = PolyDataPointerType::New(); // Delaunay filter needs DataOpbject as input...
			ratioPoints->SetDataTypeToFloat();
			ratioSampling->Allocate(1);
			
			#ifndef CONVEX_HULL_SURFS
			for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
			{
				int ID = *labelIt;
				if(avgColumns.find(ID) != avgColumns.end())
				{
					double * newPt = new double[3];
					double barrelBottom[3];
					for(int ii = 0; ii < 3; ++ii)
						barrelBottom[ii] = avgBarrels[ID]->bottom[ii];
					double dist = avgBarrelHeight[ID];
					for(int ii = 0; ii < 3; ++ii)
						newPt[ii] = barrelBottom[ii] + *ratioDepthIt*dist*avgAxes[ID][ii];
					ratioPoints->InsertNextPoint(newPt);
					if(std::find(borderBarrels.begin(), borderBarrels.end(), ID) != borderBarrels.end())
					{
						PolyDataPointerType tmpContour = PolyDataPointerType::New();
						tmpContour->Allocate(1);
						tmpContour->DeepCopy(avgColumns[ID]->contours);
						Column * tmpCol = new Column;
						tmpCol->contours = tmpContour;
						double shiftVec[3];
						for(int ii = 0; ii < 3; ++ii)
							shiftVec[ii] = newPt[ii] - avgColumns[ID]->top[ii];
						tmpCol->translateColumn(shiftVec);
						ratioPoints->InsertNextPoint(findExtremePoint(tmpCol->contours->GetCell(0)->GetPoints(), newPt, barrelFieldCOM));
						delete tmpCol;
					}
				}
			}
			#endif
			
			#ifdef CONVEX_HULL_SURFS
			for(int ii = 0; ii < L4LowerPoints->GetNumberOfPoints(); ++ii)
			{
				double pt[3];
				L4LowerPoints->GetPoint(ii, pt);
				// brute force way finding the corresponding columns
				// of all sampled L4L pts:
				double minDist = 1E09;
				int minID = -1;
				for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
				{
					int ID = *labelIt;
					if(avgBarrels.find(ID) != avgBarrels.end())
					{
						double centerDist = L2Distance3D(pt, avgBarrels[ID]->bottom);
						if(centerDist < minDist)
						{
							minDist = centerDist;
							minID = ID;
						}
						for(int n = 0; n < avgBarrels[ID]->contours->GetCell(1)->GetNumberOfPoints(); ++n)
						{
							double pt2[3];
							avgBarrels[ID]->contours->GetCell(1)->GetPoints()->GetPoint(n, pt2);
							double edgeDist =  L2Distance3D(pt, pt2);
							if(edgeDist < minDist)
							{
								minDist = edgeDist;
								minID = ID;
							}
						}
					}
				}
				if(minID >= 0)
				{
					int ID = minID;
					double barrelBottom[3], newPt[3];
					for(int jj = 0; jj < 3; ++jj)
						barrelBottom[jj] = pt[jj];
					double dist = avgBarrelHeight[ID];
					for(int jj = 0; jj < 3; ++jj)
						newPt[jj] = barrelBottom[jj] + *ratioDepthIt*dist*avgAxes[ID][jj];
					ratioPoints->InsertNextPoint(newPt);
				}
				else
					std::cout << "Oops! Something is wrong calculating the convex hull of the ratio surfs..." << std::endl;
			}
			#endif
			
			ratioSampling->SetPoints(ratioPoints);
			ratioTriangulation1->SetInput(ratioSampling);
			ratioTriangulation1->Update();
			PolyDataPointerType ratioTriangulation = ratioTriangulation1->GetOutput();
			MeshRefinementFilterPointerType ratioMeshRefinement = MeshRefinementFilterPointerType::New();
			ratioMeshRefinement->SetNumberOfSubdivisions(1);
			ratioMeshRefinement->SetInput(ratioTriangulation);
			ratioMeshRefinement->Update();
			PolyDataPointerType ratioMeshFine = ratioMeshRefinement->GetOutput();
			std::string ratioFilename(surfaceFilename);
			ratioFilename += "_celltype_border_%02d";
			// comment remaining lines for variability calculation
			// otherwise, produces too much output...
			#ifndef REG_ACCURACY
			char * borderFilename = new char[256];
			sprintf(borderFilename, ratioFilename.c_str(), surfCnt);
			Reader * ratioWriter = new Reader(borderFilename, borderFilename);
			ratioWriter->writeAmiraSurfaceFile(ratioMeshFine);
			delete ratioWriter;
			#endif
		}
		for(ratioDepthIt = cellTypeRatioDepthsInfra.begin(); ratioDepthIt != cellTypeRatioDepthsInfra.end(); ++ratioDepthIt, ++surfCnt)
		{
			Delaunay2DFilterPointerType ratioTriangulation1 = Delaunay2DFilterPointerType::New();
			PointsPointerType ratioPoints = PointsPointerType::New();
			PolyDataPointerType ratioSampling = PolyDataPointerType::New(); // Delaunay filter needs DataOpbject as input...
			ratioPoints->SetDataTypeToFloat();
			ratioSampling->Allocate(1);
			
			#ifndef CONVEX_HULL_SURFS
			for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
			{
				int ID = *labelIt;
				if(avgColumns.find(ID) != avgColumns.end())
				{
					double * newPt = new double[3];
					double barrelBottom[3];
					for(int ii = 0; ii < 3; ++ii)
						barrelBottom[ii] = avgBarrels[ID]->bottom[ii];
					double dist = avgPiaWMDist[ID] - avgTopDist[ID] - avgBarrelHeight[ID];
					for(int ii = 0; ii < 3; ++ii)
						newPt[ii] = barrelBottom[ii] + *ratioDepthIt*dist*avgAxes[ID][ii];
					ratioPoints->InsertNextPoint(newPt);
					if(std::find(borderBarrels.begin(), borderBarrels.end(), ID) != borderBarrels.end())
					{
						PolyDataPointerType tmpContour = PolyDataPointerType::New();
						tmpContour->Allocate(1);
						tmpContour->DeepCopy(avgColumns[ID]->contours);
						Column * tmpCol = new Column;
						tmpCol->contours = tmpContour;
						double shiftVec[3];
						for(int ii = 0; ii < 3; ++ii)
							shiftVec[ii] = newPt[ii] - avgColumns[ID]->top[ii];
						tmpCol->translateColumn(shiftVec);
						ratioPoints->InsertNextPoint(findExtremePoint(tmpCol->contours->GetCell(0)->GetPoints(), newPt, barrelFieldCOM));
						delete tmpCol;
					}
				}
			}
			#endif
			
			#ifdef CONVEX_HULL_SURFS
			for(int ii = 0; ii < L4LowerPoints->GetNumberOfPoints(); ++ii)
			{
				double pt[3];
				L4LowerPoints->GetPoint(ii, pt);
				// brute force way finding the corresponding columns
				// of all sampled L4L pts:
				double minDist = 1E09;
				int minID = -1;
				for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
				{
					int ID = *labelIt;
					if(avgBarrels.find(ID) != avgBarrels.end())
					{
						double centerDist = L2Distance3D(pt, avgBarrels[ID]->bottom);
						if(centerDist < minDist)
						{
							minDist = centerDist;
							minID = ID;
						}
						for(int n = 0; n < avgBarrels[ID]->contours->GetCell(1)->GetNumberOfPoints(); ++n)
						{
							double pt2[3];
							avgBarrels[ID]->contours->GetCell(1)->GetPoints()->GetPoint(n, pt2);
							double edgeDist =  L2Distance3D(pt, pt2);
							if(edgeDist < minDist)
							{
								minDist = edgeDist;
								minID = ID;
							}
						}
					}
				}
				if(minID >= 0)
				{
					int ID = minID;
					double barrelBottom[3], newPt[3];
					for(int jj = 0; jj < 3; ++jj)
						barrelBottom[jj] = pt[jj];
					double dist = avgPiaWMDist[ID] - avgTopDist[ID] - avgBarrelHeight[ID];
					for(int jj = 0; jj < 3; ++jj)
						newPt[jj] = barrelBottom[jj] + *ratioDepthIt*dist*avgAxes[ID][jj];
					ratioPoints->InsertNextPoint(newPt);
				}
				else
					std::cout << "Oops! Something is wrong calculating the convex hull of the ratio surfs..." << std::endl;
			}
			#endif
			
			ratioSampling->SetPoints(ratioPoints);
			ratioTriangulation1->SetInput(ratioSampling);
			ratioTriangulation1->Update();
			PolyDataPointerType ratioTriangulation = ratioTriangulation1->GetOutput();
			MeshRefinementFilterPointerType ratioMeshRefinement = MeshRefinementFilterPointerType::New();
			ratioMeshRefinement->SetNumberOfSubdivisions(1);
			ratioMeshRefinement->SetInput(ratioTriangulation);
			ratioMeshRefinement->Update();
			PolyDataPointerType ratioMeshFine = ratioMeshRefinement->GetOutput();
			std::string ratioFilename(surfaceFilename);
			ratioFilename += "_celltype_border_%02d";
			// comment remaining lines for variability calculation
			// otherwise, produces too much output...
			#ifndef REG_ACCURACY
			char * borderFilename = new char[256];
			sprintf(borderFilename, ratioFilename.c_str(), surfCnt);
			Reader * ratioWriter = new Reader(borderFilename, borderFilename);
			ratioWriter->writeAmiraSurfaceFile(ratioMeshFine);
			delete ratioWriter;
			#endif
		}
		// old version: linear scaling Barrel center-Pia and Barrel center-WM
// 		std::list< double >::const_iterator ratioDepthIt;
// 		int surfCnt = 0;
// 		for(ratioDepthIt = cellTypeRatioDepths.begin(); ratioDepthIt != cellTypeRatioDepths.end(); ++ratioDepthIt, ++surfCnt)
// 		{
// 			Delaunay2DFilterPointerType ratioTriangulation1 = Delaunay2DFilterPointerType::New();
// 			PointsPointerType ratioPoints = PointsPointerType::New();
// 			PolyDataPointerType ratioSampling = PolyDataPointerType::New(); // Delaunay filter needs DataOpbject as input...
// 			ratioPoints->SetDataTypeToFloat();
// 			ratioSampling->Allocate(1);
// 			
// 			#ifndef CONVEX_HULL_SURFS
// 			for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
// 			{
// 				int ID = *labelIt;
// 				if(avgColumns.find(ID) != avgColumns.end())
// 				{
// 					double * newPt = new double[3];
// 					double barrelCenter[3];
// 					for(int ii = 0; ii < 3; ++ii)
// 						barrelCenter[ii] = avgBarrels[ID]->bottom[ii] + 0.5*avgBarrels[ID]->getHeight()*avgAxes[ID][ii];
// 					double dist;
// 					if(*ratioDepthIt >= 0)
// 						dist = avgTopDist[ID] + 0.5*avgBarrels[ID]->getHeight();
// 					else
// 						dist = avgPiaWMDist[ID] - avgTopDist[ID] - 0.5*avgBarrels[ID]->getHeight();
// 					for(int ii = 0; ii < 3; ++ii)
// 						newPt[ii] = barrelCenter[ii] + *ratioDepthIt*dist*avgAxes[ID][ii];
// 					ratioPoints->InsertNextPoint(newPt);
// 					if(std::find(borderBarrels.begin(), borderBarrels.end(), ID) != borderBarrels.end())
// 					{
// 						PolyDataPointerType tmpContour = PolyDataPointerType::New();
// 						tmpContour->Allocate(1);
// 						tmpContour->DeepCopy(avgColumns[ID]->contours);
// 						Column * tmpCol = new Column;
// 						tmpCol->contours = tmpContour;
// 						double shiftVec[3];
// 						for(int ii = 0; ii < 3; ++ii)
// 							shiftVec[ii] = newPt[ii] - avgColumns[ID]->top[ii];
// 						tmpCol->translateColumn(shiftVec);
// 						ratioPoints->InsertNextPoint(findExtremePoint(tmpCol->contours->GetCell(0)->GetPoints(), newPt, barrelFieldCOM));
// 						delete tmpCol;
// 					}
// 				}
// 			}
// 			#endif
// 			
// 			#ifdef CONVEX_HULL_SURFS
// 			for(int ii = 0; ii < L4LowerPoints->GetNumberOfPoints(); ++ii)
// 			{
// 				double pt[3];
// 				L4LowerPoints->GetPoint(ii, pt);
// 				// brute force way finding the corresponding columns
// 				// of all sampled L4L pts:
// 				double minDist = 1E09;
// 				int minID = -1;
// 				for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
// 				{
// 					int ID = *labelIt;
// 					if(avgBarrels.find(ID) != avgBarrels.end())
// 					{
// 						double centerDist = L2Distance3D(pt, avgBarrels[ID]->bottom);
// 						if(centerDist < minDist)
// 						{
// 							minDist = centerDist;
// 							minID = ID;
// 						}
// 						for(int n = 0; n < avgBarrels[ID]->contours->GetCell(1)->GetNumberOfPoints(); ++n)
// 						{
// 							double pt2[3];
// 							avgBarrels[ID]->contours->GetCell(1)->GetPoints()->GetPoint(n, pt2);
// 							double edgeDist =  L2Distance3D(pt, pt2);
// 							if(edgeDist < minDist)
// 							{
// 								minDist = edgeDist;
// 								minID = ID;
// 							}
// 						}
// 					}
// 				}
// 				if(minID >= 0)
// 				{
// 					int ID = minID;
// 					double barrelCenter[3], newPt[3];
// // 					BUG?!? don't use ii here! it is already declared in the main for-loop
// 					for(int jj = 0; jj < 3; ++jj)
// 						barrelCenter[jj] = pt[jj] + 0.5*avgBarrels[ID]->getHeight()*avgAxes[ID][jj];
// 					double dist;
// 					if(*ratioDepthIt >= 0)
// 						dist = avgTopDist[ID] + 0.5*avgBarrels[ID]->getHeight();
// 					else
// 						dist = avgPiaWMDist[ID] - avgTopDist[ID] - 0.5*avgBarrels[ID]->getHeight();
// // 					BUG?!? don't use ii here! it is already declared in the main for-loop
// 					for(int jj = 0; jj < 3; ++jj)
// 						newPt[jj] = barrelCenter[jj] + *ratioDepthIt*dist*avgAxes[ID][jj];
// 					ratioPoints->InsertNextPoint(newPt);
// 				}
// 				else
// 					std::cout << "Oops! Something is wrong calculating the convex hull of the ratio surfs..." << std::endl;
// 			}
// 			#endif
// 			
// 			ratioSampling->SetPoints(ratioPoints);
// 			ratioTriangulation1->SetInput(ratioSampling);
// 			ratioTriangulation1->Update();
// 			PolyDataPointerType ratioTriangulation = ratioTriangulation1->GetOutput();
// 			MeshRefinementFilterPointerType ratioMeshRefinement = MeshRefinementFilterPointerType::New();
// 			ratioMeshRefinement->SetNumberOfSubdivisions(1);
// 			ratioMeshRefinement->SetInput(ratioTriangulation);
// 			ratioMeshRefinement->Update();
// 			PolyDataPointerType ratioMeshFine = ratioMeshRefinement->GetOutput();
// 			std::string ratioFilename(surfaceFilename);
// 			ratioFilename += "_celltype_border_%02d";
// 			// comment remaining lines for variability calculation
// 			// otherwise, produces too much output...
// 			#ifndef REG_ACCURACY
// 			char * borderFilename = new char[256];
// 			sprintf(borderFilename, ratioFilename.c_str(), surfCnt);
// 			Reader * ratioWriter = new Reader(borderFilename, borderFilename);
// 			ratioWriter->writeAmiraSurfaceFile(ratioMeshFine);
// 			delete ratioWriter;
// 			#endif
// 		}
	}
	else
		std::cout << "Error! Cannot calculate average surfaces before average Barrel Field and average Columns are calculated!" << std::endl;
};

/******************************************************************************/
/*create average axis field by triangulating axis directions measured at      */
/*lower L4 surface b/c uncertainty is minimal there                           */
/******************************************************************************/
void Registration::computeAverageAxesField()
{
	if(avgBarrels.size() && avgAxes.size())
	{
		std::cout << "Computing initial triangulation of average axis field..." << std::endl;
		PointsPointerType axisPts = PointsPointerType::New();
		PolyDataPointerType axisSampling = PolyDataPointerType::New();
		DoubleArrayPointerType axisVals = DoubleArrayPointerType::New();
		axisSampling->Allocate(1);
		axisPts->SetDataTypeToDouble();
		axisVals->SetNumberOfComponents(3);
		std::list< int >::const_iterator labelIt;
		double barrelFieldCOM[] = {0,0,0};
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(avgBarrels.find(ID) != avgBarrels.end())
			{
				for(int ii = 0; ii < 3; ++ii)
				{
					barrelFieldCOM[ii] += avgBarrels[ID]->top[ii];
					barrelFieldCOM[ii] += avgBarrels[ID]->bottom[ii];
				}
			}
		}
		for(int ii = 0; ii < 3; ++ii)
			barrelFieldCOM[ii] /= double(2*avgBarrels.size());
		
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(avgBarrels.find(ID) != avgBarrels.end() && avgAxes.find(ID) != avgAxes.end())
			{
				axisVals->InsertNextTupleValue(avgAxes[ID]);
				axisPts->InsertNextPoint(avgBarrels[ID]->bottom);
				if(std::find(borderBarrels.begin(), borderBarrels.end(), ID) != borderBarrels.end())
				{
					double * extrapolAxisPt = findExtremePoint(avgBarrels[ID]->contours->GetCell(1)->GetPoints(), avgBarrels[ID]->bottom, barrelFieldCOM);
					axisPts->InsertNextPoint(extrapolAxisPt);
					axisVals->InsertNextTupleValue(findExtremeAxis(ID, extrapolAxisPt));
				}
			}
		}
		axisSampling->SetPoints(axisPts);
		axisSampling->GetPointData()->SetVectors(axisVals);
		Delaunay2DFilterPointerType axisTriangulation = Delaunay2DFilterPointerType::New();
		axisTriangulation->SetInput(axisSampling);
		axisTriangulation->Update();
		PolyDataPointerType axisMesh = axisTriangulation->GetOutput();
		
		std::flush(std::cout << "Refining mesh for axis field..." << std::endl);
		// re-normalize axes when subdividing recursively??? b/c norm prob not conserved during interpolation
		MeshRefinementFilterPointerType axisMeshRefinement = MeshRefinementFilterPointerType::New();
		axisMeshRefinement->SetNumberOfSubdivisions(1);
		axisMeshRefinement->SetInput(axisMesh);
		axisMeshRefinement->Update();
		PolyDataPointerType axisMeshFine = axisMeshRefinement->GetOutput();
		
// 		PointsPointerType meshPts = axisMesh->GetPoints();
// 		DataArrayPointerType axisMeshData = axisMesh->GetPointData()->GetVectors();
		PointsPointerType meshPts = axisMeshFine->GetPoints();
		DataArrayPointerType axisMeshData = axisMeshFine->GetPointData()->GetVectors();
		for(int ii = 0; ii < meshPts->GetNumberOfPoints() && ii < axisMeshData->GetNumberOfTuples(); ++ii)
		{
			double * top = new double[3], * bottom = new double[3], axis[3], phi, theta;
			meshPts->GetPoint(ii, bottom);
			axis[0] = axisMeshData->GetComponent(ii, 0), axis[1] = axisMeshData->GetComponent(ii, 1), axis[2] = axisMeshData->GetComponent(ii, 2);
			normalize(axis);
			for(int jj = 0; jj < 3; ++jj)
			{
				top[jj] = bottom[jj] + 600*axis[jj];
				bottom[jj] -= 500*axis[jj];
			}
			
			if(!spatialGraph)
				spatialGraph = new AmiraSpatialGraph;
			Vertex * newVert1 = new Vertex(top, ZAxis);
			Vertex * newVert2 = new Vertex(bottom, ZAxis);
			spatialGraph->addVertex(newVert1), spatialGraph->addVertex(newVert2);
			int connectivity[2];
			connectivity[0] = spatialGraph->getNumberOfVertices()-2, connectivity[1] = spatialGraph->getNumberOfVertices()-1;
			std::list< double * > edgePtList;
			edgePtList.push_back(top), edgePtList.push_back(bottom);
			Edge * newEdge = new Edge(connectivity, 2, ZAxis, edgePtList);
			spatialGraph->addEdge(newEdge);
		}
	}
	else
		std::cout << "Error! Cannot calculate average axis field before average barrel field and average axes are calculated!" << std::endl;
};

void Registration::enforceAxisDivergence(std::map< int, double * > barrelAxes, std::map< int, double * > barrelCenters)
{
	std::map< int, double * >::iterator barrelAxesIt;
	std::map< int, double * >::reverse_iterator barrelAxesRevIt;
// 	std::map< int, std::list< int > >::iterator barrelGridIt;
// 	std::map< int, std::list< int > >::reverse_iterator barrelGridRevIt;
	for(barrelAxesIt = barrelAxes.begin(); barrelAxesIt != barrelAxes.end(); ++barrelAxesIt)
		normalize(barrelAxesIt->second);
	std::map< int, std::list< int > > barrelGrid = createBarrelGrid(barrelAxes);
	
	int maxLoops = 0;
	bool change = 0;
	do
	{
// 		std::cout << "loop nr " << maxLoops << std::endl;
		change = 0;
		for(barrelAxesIt = barrelAxes.begin(); barrelAxesIt != barrelAxes.end(); ++barrelAxesIt)
		{
			int thisBarrel = barrelAxesIt->first;
			double thisAxis[3], cumulatedChange[3];
			thisAxis[0] = barrelAxesIt->second[0], thisAxis[1] = barrelAxesIt->second[1], thisAxis[2] = barrelAxesIt->second[2];
			cumulatedChange[0] = cumulatedChange[1] = cumulatedChange[2] = 0;
			std::list< int >::iterator neighborIt;
			for(neighborIt = barrelGrid[barrelAxesIt->first].begin(); neighborIt != barrelGrid[barrelAxesIt->first].end(); ++neighborIt)
			{
				double neighborAxis[3];
				neighborAxis[0] = barrelAxes[*neighborIt][0], neighborAxis[1] = barrelAxes[*neighborIt][1], neighborAxis[2] = barrelAxes[*neighborIt][2];
				double direction[3];
				for(int ii = 0; ii < 3; ++ii)
					direction[ii] = barrelCenters[*neighborIt][ii] - barrelCenters[thisBarrel][ii];
				normalize(direction);
				//eliminate part of neighboraxis that is parallel to this axis
				double prod = 0;
				for(int ii = 0; ii < 3; ++ii)
					prod += thisAxis[ii]*neighborAxis[ii];
				for(int ii = 0; ii < 3; ++ii)
					neighborAxis[ii] -= prod*thisAxis[ii];
				//check if other axis is pointing out of our local coordinate system
				double localDirection = 0;
				for(int ii = 0; ii < 3; ++ii)
					localDirection += direction[ii]*neighborAxis[ii];
				// < 0 means other axis is pointing towards us
				if(localDirection < 0)
				{
					change = 1;
					for(int ii = 0; ii < 3; ++ii)
					{
						cumulatedChange[ii] += 0.5*localDirection*direction[ii];
						barrelAxes[*neighborIt][ii] -= 0.5*localDirection*direction[ii];
					}
					normalize(barrelAxes[*neighborIt]);
				}
			}
			for(int ii = 0; ii < 3; ++ii)
				barrelAxesIt->second[ii] += cumulatedChange[ii];
			normalize(barrelAxesIt->second);
		}
		if(change)
		{
			for(barrelAxesRevIt = barrelAxes.rbegin(); barrelAxesRevIt != barrelAxes.rend(); ++barrelAxesRevIt)
			{
				int thisBarrel = barrelAxesRevIt->first;
				double thisAxis[3], cumulatedChange[3];
				thisAxis[0] = barrelAxesRevIt->second[0], thisAxis[1] = barrelAxesRevIt->second[1], thisAxis[2] = barrelAxesRevIt->second[2];
				cumulatedChange[0] = cumulatedChange[1] = cumulatedChange[2] = 0;
				std::list< int >::iterator neighborIt;
				for(neighborIt = barrelGrid[barrelAxesRevIt->first].begin(); neighborIt != barrelGrid[barrelAxesRevIt->first].end(); ++neighborIt)
				{
					double neighborAxis[3];
					neighborAxis[0] = barrelAxes[*neighborIt][0], neighborAxis[1] = barrelAxes[*neighborIt][1], neighborAxis[2] = barrelAxes[*neighborIt][2];
					double direction[3];
					for(int ii = 0; ii < 3; ++ii)
						direction[ii] = barrelCenters[*neighborIt][ii] - barrelCenters[thisBarrel][ii];
					normalize(direction);
					//eliminate part of neighboraxis that is parallel to this axis
					double prod = 0;
					for(int ii = 0; ii < 3; ++ii)
						prod += thisAxis[ii]*neighborAxis[ii];
					for(int ii = 0; ii < 3; ++ii)
						neighborAxis[ii] -= prod*thisAxis[ii];
					//check if other axis is pointing out of our local coordinate system
					double localDirection = 0;
					for(int ii = 0; ii < 3; ++ii)
						localDirection += direction[ii]*neighborAxis[ii];
					// < 0 means other axis is pointing towards us
					if(localDirection < 0)
					{
						change = 1;
						for(int ii = 0; ii < 3; ++ii)
						{
							cumulatedChange[ii] += 0.5*localDirection*direction[ii];
							barrelAxes[*neighborIt][ii] -= 0.5*localDirection*direction[ii];
						}
						normalize(barrelAxes[*neighborIt]);
					}
				}
				for(int ii = 0; ii < 3; ++ii)
					barrelAxesRevIt->second[ii] += cumulatedChange[ii];
				normalize(barrelAxesRevIt->second);
			}
		}
		++maxLoops;
	} while(change && maxLoops < 10);
	
};

/******************************************************************************/
/*returns matrix that gives coordinates of vectors from data coordinate       */
/*system in local barrel coordinate system determined by newAxis              */
/******************************************************************************/
HomogeneousMatrixPointerType Registration::getLocalBarrelCoordinates(double * newAxis)
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
/*returns matrix that gives coordinates of vectors from data coordinate       */
/*system in local barrel coordinate system determined by newAxis              */
/******************************************************************************/
HomogeneousMatrixPointerType Registration::getLocalCoordinatesForVariation(double * newAxis)
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
HomogeneousMatrixPointerType Registration::transformToBarrelCoordinates(double * newAxis)
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
/*returns rotation matrix in xy-plane to get new x axis                       */
/******************************************************************************/
HomogeneousMatrixPointerType Registration::getNewXAxis(double angle)
{
	TransformPointerType rotate = TransformPointerType::New();
	double rotationAxis[3];
	rotationAxis[0] = 0;
	rotationAxis[1] = 0;
	rotationAxis[2] = -1;
	
	rotate->RotateWXYZ(angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);
	
	return rotate->GetMatrix();
};

/******************************************************************************/
/*returns rotation matrix in yz-plane to get new z axis                       */
/******************************************************************************/
HomogeneousMatrixPointerType Registration::getNewZAxis(double angle)
{
	TransformPointerType rotate = TransformPointerType::New();
	double rotationAxis[3];
	rotationAxis[0] = 1;
	rotationAxis[1] = 0;
	rotationAxis[2] = 0;
	
	rotate->RotateWXYZ(angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);
	
	return rotate->GetMatrix();
};

/******************************************************************************/
/*returns point that has greatest distance from origin measured in xy-plane   */
/*and "dilates" that point by moving it radially away from center             */
/******************************************************************************/
double * Registration::findExtremePoint(PointsPointerType pts, double center[3], double com[3])
{
	double * maxPt = new double[3];
	double maxDist = 0;
	for(int ii = 0; ii < pts->GetNumberOfPoints(); ++ii)
	{
		double dist = 0, pt[3];
		pts->GetPoint(ii, pt);
		double pt2D[] = {pt[0], pt[1]}, com2D[] = {com[0], com[1]};
		dist = L2Distance2D(pt2D, com2D);
		if(dist > maxDist)
		{
			maxDist = dist;
			maxPt[0] = pt[0], maxPt[1] = pt[1], maxPt[2] = pt[2];
		}
	}
// 	double radial[3];
// 	for(int ii = 0; ii < 3; ++ii)
// 		radial[ii] = maxPt[ii] - center[ii];
// 	normalize(radial);
// 	for(int ii = 0; ii < 3; ++ii)
// 		maxPt[ii] += 400*radial[ii];
	
	double dilation = 400.0;
	double radial[3];
	for(int ii = 0; ii < 3; ++ii)
		radial[ii] = maxPt[ii] - center[ii];
	normalize(radial);
	for(int ii = 0; ii < 3; ++ii)
		maxPt[ii] += dilation*radial[ii];
	
	std::list< int >::const_iterator labelIt;
	for(labelIt = borderBarrels.begin(); labelIt != borderBarrels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(avgColumns.find(ID) != avgColumns.end())
		{
			double radius, distance, diff;
			double t, closestPoint[3];
			radius = sqrt(avgBarrelArea[ID]/PI);
			distance = vtkLine::DistanceToLine(maxPt, avgColumns[ID]->top, avgColumns[ID]->bottom, t, closestPoint);
			distance = sqrt(distance) - radius;
			diff = distance - dilation;
			if(diff < 0)
			{
				for(int ii = 0; ii < 3; ++ii)
					radial[ii] = maxPt[ii] - closestPoint[ii];
				normalize(radial);
				for(int ii = 0; ii < 3; ++ii)
					maxPt[ii] += -diff*radial[ii];
			}
		}
	}
	return maxPt;
};

/******************************************************************************/
/*returns axis vector that is linearly extrapolated from Barrel thisID by     */
/*approximating the axis normal gradient from the neighboring barrel axes     */
/*but: gradient is approximated as grad(axis_i) = unitvec_i \partial_i axis_i */
/*only way to get reasonable results b/c barrels are not on a grid!!!         */
/******************************************************************************/
double * Registration::findExtremeAxis(int thisID, double extrapolatedPt[3])
{
	double * newAxis = new double[3], ** axisGradient = new double* [3], direction[3];
	if(avgAxes.size() && avgBarrels.size())
	{
		std::map< int, std::list< int > > barrelGrid = createBarrelGrid(avgAxes);
		std::list< int >::const_iterator neighborIt;
		if(barrelGrid.find(thisID) != barrelGrid.end())
		{
			axisGradient[0] = new double[3], axisGradient[1] = new double[3], axisGradient[2] = new double[3];
			double * thisPt = avgBarrels[thisID]->bottom, * thisAxis = avgAxes[thisID];
			for(int ii = 0; ii < 3; ++ii)
			{
				direction[ii] = extrapolatedPt[ii] - thisPt[ii];
				axisGradient[0][ii] = axisGradient[1][ii] = axisGradient[2][ii] = 0;
			}
			normalize(direction);
// 			if(thisID == Alpha)
// 				std::cout << "direction = [" << direction[0] << "," << direction[1] << "," << direction[2] << "]" << std::endl;
			int neighborCnt = barrelGrid[thisID].size();
			double * weights = new double[neighborCnt], weightNorm = 0;
			int cnt = 0;
			for(neighborIt = barrelGrid[thisID].begin(); neighborIt != barrelGrid[thisID].end(); ++neighborIt, ++cnt)
			{
				weights[cnt] = L2Distance3D(thisPt, avgBarrels[*neighborIt]->bottom);
				if(weights[cnt])
				{
					weights[cnt] = 1/weights[cnt];
					weightNorm += weights[cnt];
				}
			}
			if(weightNorm)
				for(int ii = 0; ii < cnt; ++ii)
				{
					weights[ii] /= weightNorm;
// 					if(thisID == Alpha)
// 						std::cout << "weights[" << ii << "] = " << weights[ii] << std::endl;
				}
			cnt = 0;
			for(neighborIt = barrelGrid[thisID].begin(); neighborIt != barrelGrid[thisID].end(); ++neighborIt, ++cnt)
			{
				int nextID = *neighborIt;
				double * nextPt = avgBarrels[nextID]->bottom, * nextAxis = avgAxes[nextID];
				double tmpSlope[3];
				for(int ii = 0; ii < 3; ++ii)
				{
					tmpSlope[ii] = thisAxis[ii] - nextAxis[ii];
					for(int jj = 0; jj < 3; ++jj)
					{
						double dist = thisPt[jj]-nextPt[jj];
						dist += std::abs(dist) > 150 ? 0 : dist > 0 ? 150 : -150;	// regularization
						if(dist)
							axisGradient[ii][jj] += weights[cnt]*tmpSlope[ii]/dist;
// 						if(thisID == Alpha)
// 							std::cout << "axisGradient[" << ii << "][" << jj << "] = " << axisGradient[ii][jj] << std::endl;
					}
				}
			}
			delete [] weights;
			// actual slope in direction: direction * gradient
			double slope[] = {0, 0, 0};
			for(int ii = 0; ii < 3; ++ii)
				for(int jj = 0; jj < 3; ++jj)
				{
					if(ii != jj)
						axisGradient[ii][jj] = 0;
					slope[ii] += direction[jj]*axisGradient[ii][jj];
// 					if(thisID == Alpha && ii == jj)
// 						std::cout << "slope[" << ii << "] = " << slope[ii] << std::endl;
				}
			double dist = L2Distance3D(extrapolatedPt, thisPt);
// 			if(thisID == Alpha)
// 				std::cout << "dist = " << dist << std::endl;
			for(int ii = 0; ii < 3; ++ii)
				newAxis[ii] = thisAxis[ii] + dist*slope[ii];
			normalize(newAxis);
		}
	}
	else
		std::cout << "Error! Cannot extrapolate axis before average barrel field and average axes are calculated!" << std::endl;
	delete [] axisGradient;
	return newAxis;
};

/******************************************************************************/
/*create grid of barrels actually present in data: each barrel label has a    */
/*LUT of neighbors associated with it (i.e., create a graph)                  */
/******************************************************************************/
std::map< int, std::list< int > > Registration::createBarrelGrid(std::map< int, double * > barrelAxes)
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
void Registration::normalize(double * vec)
{
	double norm = 0;
	for(int ii = 0; ii < 3; ++ii)
		norm += vec[ii]*vec[ii];
	norm = sqrt(norm);
	if(norm)
		for(int ii = 0; ii < 3; ++ii)
			vec[ii] = vec[ii]/norm;
};

double Registration::L2Distance3D(double x[3], double y[3])
{
	return sqrt((x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]) + (x[2] - y[2])*(x[2] - y[2]));
};

double Registration::L2Distance2D(double x[2], double y[2])
{
	return sqrt((x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]));
};

void Registration::initializeConstants()
{
	if(barrelLabels.size())
		barrelLabels.clear();
	barrelLabels.push_back(Alpha);
	barrelLabels.push_back(A1);
	barrelLabels.push_back(A2);
	barrelLabels.push_back(A3);
	barrelLabels.push_back(A4);
	barrelLabels.push_back(Beta);
	barrelLabels.push_back(B1);
	barrelLabels.push_back(B2);
	barrelLabels.push_back(B3);
	barrelLabels.push_back(B4);
	barrelLabels.push_back(Gamma);
	barrelLabels.push_back(C1);
	barrelLabels.push_back(C2);
	barrelLabels.push_back(C3);
	barrelLabels.push_back(C4);
	barrelLabels.push_back(C5);
	barrelLabels.push_back(C6);
	barrelLabels.push_back(Delta);
	barrelLabels.push_back(D1);
	barrelLabels.push_back(D2);
	barrelLabels.push_back(D3);
	barrelLabels.push_back(D4);
	barrelLabels.push_back(D5);
	barrelLabels.push_back(D6);
	barrelLabels.push_back(E1);
	barrelLabels.push_back(E2);
	barrelLabels.push_back(E3);
	barrelLabels.push_back(E4);
	barrelLabels.push_back(E5);
	barrelLabels.push_back(E6);
	if(borderBarrels.size())
		borderBarrels.clear();
	borderBarrels.push_back(Alpha);
	borderBarrels.push_back(A1);
	borderBarrels.push_back(A2);
	borderBarrels.push_back(A3);
	borderBarrels.push_back(A4);
	borderBarrels.push_back(Beta);
	borderBarrels.push_back(B4);
	borderBarrels.push_back(Gamma);
	borderBarrels.push_back(C4);
	borderBarrels.push_back(Delta);
	borderBarrels.push_back(D4);
	borderBarrels.push_back(E1);
	borderBarrels.push_back(E2);
	borderBarrels.push_back(E3);
	borderBarrels.push_back(E4);
	if(int2Labels.size())
		int2Labels.clear();
	int2Labels.insert(std::pair< int, const char * >(Alpha, "Alpha"));
	int2Labels.insert(std::pair< int, const char * >(A1, "A1"));
	int2Labels.insert(std::pair< int, const char * >(A2, "A2"));
	int2Labels.insert(std::pair< int, const char * >(A3, "A3"));
	int2Labels.insert(std::pair< int, const char * >(A4, "A4"));
	int2Labels.insert(std::pair< int, const char * >(Beta, "Beta"));
	int2Labels.insert(std::pair< int, const char * >(B1, "B1"));
	int2Labels.insert(std::pair< int, const char * >(B2, "B2"));
	int2Labels.insert(std::pair< int, const char * >(B3, "B3"));
	int2Labels.insert(std::pair< int, const char * >(B4, "B4"));
	int2Labels.insert(std::pair< int, const char * >(Gamma, "Gamma"));
	int2Labels.insert(std::pair< int, const char * >(C1, "C1"));
	int2Labels.insert(std::pair< int, const char * >(C2, "C2"));
	int2Labels.insert(std::pair< int, const char * >(C3, "C3"));
	int2Labels.insert(std::pair< int, const char * >(C4, "C4"));
	int2Labels.insert(std::pair< int, const char * >(C5, "C5"));
	int2Labels.insert(std::pair< int, const char * >(C6, "C6"));
	int2Labels.insert(std::pair< int, const char * >(Delta, "Delta"));
	int2Labels.insert(std::pair< int, const char * >(D1, "D1"));
	int2Labels.insert(std::pair< int, const char * >(D2, "D2"));
	int2Labels.insert(std::pair< int, const char * >(D3, "D3"));
	int2Labels.insert(std::pair< int, const char * >(D4, "D4"));
	int2Labels.insert(std::pair< int, const char * >(D5, "D5"));
	int2Labels.insert(std::pair< int, const char * >(D6, "D6"));
	int2Labels.insert(std::pair< int, const char * >(E1, "E1"));
	int2Labels.insert(std::pair< int, const char * >(E2, "E2"));
	int2Labels.insert(std::pair< int, const char * >(E3, "E3"));
	int2Labels.insert(std::pair< int, const char * >(E4, "E4"));
	int2Labels.insert(std::pair< int, const char * >(E5, "E5"));
	int2Labels.insert(std::pair< int, const char * >(E6, "E6"));
	
	// 12 BF:
	if(avgBarrelHeight.size())
		avgBarrelHeight.clear();
	avgBarrelHeight.insert(std::pair< int, double >(Alpha, 308));
	avgBarrelHeight.insert(std::pair< int, double >(A1, 322));
	avgBarrelHeight.insert(std::pair< int, double >(A2, 337));
	avgBarrelHeight.insert(std::pair< int, double >(A3, 315));
	avgBarrelHeight.insert(std::pair< int, double >(A4, 316));
	avgBarrelHeight.insert(std::pair< int, double >(Beta, 338));
	avgBarrelHeight.insert(std::pair< int, double >(B1, 346));
	avgBarrelHeight.insert(std::pair< int, double >(B2, 344));
	avgBarrelHeight.insert(std::pair< int, double >(B3, 354));
	avgBarrelHeight.insert(std::pair< int, double >(B4, 352));
	avgBarrelHeight.insert(std::pair< int, double >(Gamma, 355));
	avgBarrelHeight.insert(std::pair< int, double >(C1, 358));
	avgBarrelHeight.insert(std::pair< int, double >(C2, 360));
	avgBarrelHeight.insert(std::pair< int, double >(C3, 354));
	avgBarrelHeight.insert(std::pair< int, double >(C4, 363));
	avgBarrelHeight.insert(std::pair< int, double >(C5, 343));
	avgBarrelHeight.insert(std::pair< int, double >(C6, 352));
	avgBarrelHeight.insert(std::pair< int, double >(Delta, 354));
	avgBarrelHeight.insert(std::pair< int, double >(D1, 371));
	avgBarrelHeight.insert(std::pair< int, double >(D2, 362));
	avgBarrelHeight.insert(std::pair< int, double >(D3, 368));
	avgBarrelHeight.insert(std::pair< int, double >(D4, 363));
	avgBarrelHeight.insert(std::pair< int, double >(D5, 349));
	avgBarrelHeight.insert(std::pair< int, double >(D6, 352));
	avgBarrelHeight.insert(std::pair< int, double >(E1, 343));
	avgBarrelHeight.insert(std::pair< int, double >(E2, 355));
	avgBarrelHeight.insert(std::pair< int, double >(E3, 364));
	avgBarrelHeight.insert(std::pair< int, double >(E4, 344));
	avgBarrelHeight.insert(std::pair< int, double >(E5, 339));
	avgBarrelHeight.insert(std::pair< int, double >(E6, 365));
	if(avgTopDist.size())
		avgTopDist.clear();
	avgTopDist.insert(std::pair< int, double >(Alpha, 479));
	avgTopDist.insert(std::pair< int, double >(A1, 455));
	avgTopDist.insert(std::pair< int, double >(A2, 467));
	avgTopDist.insert(std::pair< int, double >(A3, 485));
	avgTopDist.insert(std::pair< int, double >(A4, 489));
	avgTopDist.insert(std::pair< int, double >(Beta, 472));
	avgTopDist.insert(std::pair< int, double >(B1, 481));
	avgTopDist.insert(std::pair< int, double >(B2, 490));
	avgTopDist.insert(std::pair< int, double >(B3, 490));
	avgTopDist.insert(std::pair< int, double >(B4, 501));
	avgTopDist.insert(std::pair< int, double >(Gamma, 478));
	avgTopDist.insert(std::pair< int, double >(C1, 478));
	avgTopDist.insert(std::pair< int, double >(C2, 496));
	avgTopDist.insert(std::pair< int, double >(C3, 534));
	avgTopDist.insert(std::pair< int, double >(C4, 557));
	avgTopDist.insert(std::pair< int, double >(C5, 549));
	avgTopDist.insert(std::pair< int, double >(C6, 526));
	avgTopDist.insert(std::pair< int, double >(Delta, 506));
	avgTopDist.insert(std::pair< int, double >(D1, 491));
	avgTopDist.insert(std::pair< int, double >(D2, 526));
	avgTopDist.insert(std::pair< int, double >(D3, 552));
	avgTopDist.insert(std::pair< int, double >(D4, 556));
	avgTopDist.insert(std::pair< int, double >(D5, 567));
	avgTopDist.insert(std::pair< int, double >(D6, 549));
	avgTopDist.insert(std::pair< int, double >(E1, 546));
	avgTopDist.insert(std::pair< int, double >(E2, 557));
	avgTopDist.insert(std::pair< int, double >(E3, 549));
	avgTopDist.insert(std::pair< int, double >(E4, 580));
	avgTopDist.insert(std::pair< int, double >(E5, 566));
	avgTopDist.insert(std::pair< int, double >(E6, 545));
	if(avgPiaWMDist.size())
		avgPiaWMDist.clear();
	avgPiaWMDist.insert(std::pair< int, double >(Alpha, 1600));
	avgPiaWMDist.insert(std::pair< int, double >(A1, 1651));
	avgPiaWMDist.insert(std::pair< int, double >(A2, 1759));
	avgPiaWMDist.insert(std::pair< int, double >(A3, 1825));
	avgPiaWMDist.insert(std::pair< int, double >(A4, 1916));
	avgPiaWMDist.insert(std::pair< int, double >(Beta, 1623));
	avgPiaWMDist.insert(std::pair< int, double >(B1, 1736));
	avgPiaWMDist.insert(std::pair< int, double >(B2, 1815));
	avgPiaWMDist.insert(std::pair< int, double >(B3, 1899));
	avgPiaWMDist.insert(std::pair< int, double >(B4, 1961));
	avgPiaWMDist.insert(std::pair< int, double >(Gamma, 1713));
	avgPiaWMDist.insert(std::pair< int, double >(C1, 1800));
	avgPiaWMDist.insert(std::pair< int, double >(C2, 1892));
	avgPiaWMDist.insert(std::pair< int, double >(C3, 1985));
	avgPiaWMDist.insert(std::pair< int, double >(C4, 2038));
	avgPiaWMDist.insert(std::pair< int, double >(C5, 2036));
	avgPiaWMDist.insert(std::pair< int, double >(C6, 2064));
	avgPiaWMDist.insert(std::pair< int, double >(Delta, 1845));
	avgPiaWMDist.insert(std::pair< int, double >(D1, 1865));
	avgPiaWMDist.insert(std::pair< int, double >(D2, 1957));
	avgPiaWMDist.insert(std::pair< int, double >(D3, 2046));
	avgPiaWMDist.insert(std::pair< int, double >(D4, 2081));
	avgPiaWMDist.insert(std::pair< int, double >(D5, 2071));
	avgPiaWMDist.insert(std::pair< int, double >(D6, 2087));
	avgPiaWMDist.insert(std::pair< int, double >(E1, 1977));
	avgPiaWMDist.insert(std::pair< int, double >(E2, 2096));
	avgPiaWMDist.insert(std::pair< int, double >(E3, 2117));
	avgPiaWMDist.insert(std::pair< int, double >(E4, 2111));
	avgPiaWMDist.insert(std::pair< int, double >(E5, 2087));
	avgPiaWMDist.insert(std::pair< int, double >(E6, 2119));
	if(avgBarrelArea.size())
		avgBarrelArea.clear();
	avgBarrelArea.insert(std::pair< int, double >(Alpha, 8.79E4));
	avgBarrelArea.insert(std::pair< int, double >(A1, 8.55E4));
	avgBarrelArea.insert(std::pair< int, double >(A2, 8.30E4));
	avgBarrelArea.insert(std::pair< int, double >(A3, 6.48E4));
	avgBarrelArea.insert(std::pair< int, double >(A4, 6.50E4));
	avgBarrelArea.insert(std::pair< int, double >(Beta, 10.2E4));
	avgBarrelArea.insert(std::pair< int, double >(B1, 8.74E4));
	avgBarrelArea.insert(std::pair< int, double >(B2, 9.01E4));
	avgBarrelArea.insert(std::pair< int, double >(B3, 7.66E4));
	avgBarrelArea.insert(std::pair< int, double >(B4, 7.89E4));
	avgBarrelArea.insert(std::pair< int, double >(Gamma, 12.7E4));
	avgBarrelArea.insert(std::pair< int, double >(C1, 10.1E4));
	avgBarrelArea.insert(std::pair< int, double >(C2, 10.3E4));
	avgBarrelArea.insert(std::pair< int, double >(C3, 10.4E4));
	avgBarrelArea.insert(std::pair< int, double >(C4, 9.03E4));
	avgBarrelArea.insert(std::pair< int, double >(C5, 7.06E4));
	avgBarrelArea.insert(std::pair< int, double >(C6, 5.63E4));
	avgBarrelArea.insert(std::pair< int, double >(Delta, 14.3E4));
	avgBarrelArea.insert(std::pair< int, double >(D1, 11.2E4));
	avgBarrelArea.insert(std::pair< int, double >(D2, 12.4E4));
	avgBarrelArea.insert(std::pair< int, double >(D3, 11.6E4));
	avgBarrelArea.insert(std::pair< int, double >(D4, 10.3E4));
	avgBarrelArea.insert(std::pair< int, double >(D5, 8.54E4));
	avgBarrelArea.insert(std::pair< int, double >(D6, 6.40E4));
	avgBarrelArea.insert(std::pair< int, double >(E1, 14.6E4));
	avgBarrelArea.insert(std::pair< int, double >(E2, 15.9E4));
	avgBarrelArea.insert(std::pair< int, double >(E3, 15.8E4));
	avgBarrelArea.insert(std::pair< int, double >(E4, 12.3E4));
	avgBarrelArea.insert(std::pair< int, double >(E5, 8.89E4));
	avgBarrelArea.insert(std::pair< int, double >(E6, 7.31E4));
	
	// 6 BF:
// 	if(avgTopDist.size())
// 		avgTopDist.clear();
// 	avgTopDist.insert(std::pair< int, double >(Alpha, 496));
// 	avgTopDist.insert(std::pair< int, double >(A1, 459));
// 	avgTopDist.insert(std::pair< int, double >(A2, 480));
// 	avgTopDist.insert(std::pair< int, double >(A3, 507));
// 	avgTopDist.insert(std::pair< int, double >(A4, 501));
// 	avgTopDist.insert(std::pair< int, double >(Beta, 492));
// 	avgTopDist.insert(std::pair< int, double >(B1, 489));
// 	avgTopDist.insert(std::pair< int, double >(B2, 491));
// 	avgTopDist.insert(std::pair< int, double >(B3, 469));
// 	avgTopDist.insert(std::pair< int, double >(B4, 491));
// 	avgTopDist.insert(std::pair< int, double >(Gamma, 487));
// 	avgTopDist.insert(std::pair< int, double >(C1, 485));
// 	avgTopDist.insert(std::pair< int, double >(C2, 487));
// 	avgTopDist.insert(std::pair< int, double >(C3, 508));
// 	avgTopDist.insert(std::pair< int, double >(C4, 531));
// 	avgTopDist.insert(std::pair< int, double >(C5, 544));
// 	avgTopDist.insert(std::pair< int, double >(C6, 506));
// 	avgTopDist.insert(std::pair< int, double >(Delta, 525));
// 	avgTopDist.insert(std::pair< int, double >(D1, 491));
// 	avgTopDist.insert(std::pair< int, double >(D2, 544));
// 	avgTopDist.insert(std::pair< int, double >(D3, 547));
// 	avgTopDist.insert(std::pair< int, double >(D4, 549));
// 	avgTopDist.insert(std::pair< int, double >(D5, 541));
// 	avgTopDist.insert(std::pair< int, double >(D6, 559));
// 	avgTopDist.insert(std::pair< int, double >(E1, 552));
// 	avgTopDist.insert(std::pair< int, double >(E2, 568));
// 	avgTopDist.insert(std::pair< int, double >(E3, 530));
// 	avgTopDist.insert(std::pair< int, double >(E4, 571));
// 	avgTopDist.insert(std::pair< int, double >(E5, 533));
// 	avgTopDist.insert(std::pair< int, double >(E6, 538));
// 	if(avgPiaWMDist.size())
// 		avgPiaWMDist.clear();
// 	avgPiaWMDist.insert(std::pair< int, double >(Alpha, 1573));
// 	avgPiaWMDist.insert(std::pair< int, double >(A1, 1630));
// 	avgPiaWMDist.insert(std::pair< int, double >(A2, 1754));
// 	avgPiaWMDist.insert(std::pair< int, double >(A3, 1805));
// 	avgPiaWMDist.insert(std::pair< int, double >(A4, 1856));
// 	avgPiaWMDist.insert(std::pair< int, double >(Beta, 1578));
// 	avgPiaWMDist.insert(std::pair< int, double >(B1, 1726));
// 	avgPiaWMDist.insert(std::pair< int, double >(B2, 1792));
// 	avgPiaWMDist.insert(std::pair< int, double >(B3, 1828));
// 	avgPiaWMDist.insert(std::pair< int, double >(B4, 1911));
// 	avgPiaWMDist.insert(std::pair< int, double >(Gamma, 1681));
// 	avgPiaWMDist.insert(std::pair< int, double >(C1, 1786));
// 	avgPiaWMDist.insert(std::pair< int, double >(C2, 1851));
// 	avgPiaWMDist.insert(std::pair< int, double >(C3, 1929));
// 	avgPiaWMDist.insert(std::pair< int, double >(C4, 1966));
// 	avgPiaWMDist.insert(std::pair< int, double >(C5, 1997));
// 	avgPiaWMDist.insert(std::pair< int, double >(C6, 2017));
// 	avgPiaWMDist.insert(std::pair< int, double >(Delta, 1853));
// 	avgPiaWMDist.insert(std::pair< int, double >(D1, 1825));
// 	avgPiaWMDist.insert(std::pair< int, double >(D2, 1938));
// 	avgPiaWMDist.insert(std::pair< int, double >(D3, 2006));
// 	avgPiaWMDist.insert(std::pair< int, double >(D4, 2054));
// 	avgPiaWMDist.insert(std::pair< int, double >(D5, 2017));
// 	avgPiaWMDist.insert(std::pair< int, double >(D6, 2076));
// 	avgPiaWMDist.insert(std::pair< int, double >(E1, 1957));
// 	avgPiaWMDist.insert(std::pair< int, double >(E2, 2079));
// 	avgPiaWMDist.insert(std::pair< int, double >(E3, 2077));
// 	avgPiaWMDist.insert(std::pair< int, double >(E4, 2086));
// 	avgPiaWMDist.insert(std::pair< int, double >(E5, 2031));
// 	avgPiaWMDist.insert(std::pair< int, double >(E6, 2071));
// 	if(avgBarrelArea.size())
// 		avgBarrelArea.clear();
// 	avgBarrelArea.insert(std::pair< int, double >(Alpha, 7.99E4));
// 	avgBarrelArea.insert(std::pair< int, double >(A1, 7.49E4));
// 	avgBarrelArea.insert(std::pair< int, double >(A2, 8.13E4));
// 	avgBarrelArea.insert(std::pair< int, double >(A3, 6.15E4));
// 	avgBarrelArea.insert(std::pair< int, double >(A4, 6.14E4));
// 	avgBarrelArea.insert(std::pair< int, double >(Beta, 9.52E4));
// 	avgBarrelArea.insert(std::pair< int, double >(B1, 7.96E4));
// 	avgBarrelArea.insert(std::pair< int, double >(B2, 8.16E4));
// 	avgBarrelArea.insert(std::pair< int, double >(B3, 7.60E4));
// 	avgBarrelArea.insert(std::pair< int, double >(B4, 7.84E4));
// 	avgBarrelArea.insert(std::pair< int, double >(Gamma, 11.0E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C1, 9.73E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C2, 10.2E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C3, 9.77E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C4, 8.61E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C5, 7.33E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C6, 5.82E4));
// 	avgBarrelArea.insert(std::pair< int, double >(Delta, 12.0E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D1, 10.8E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D2, 11.7E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D3, 11.0E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D4, 10.1E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D5, 8.10E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D6, 6.26E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E1, 13.3E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E2, 15.1E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E3, 15.0E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E4, 11.9E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E5, 8.76E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E6, 6.70E4));
	
// 	21Jan2013: old, WRONG version
// 	(assumed Standard Column BC == SBF D2 BC, which is off by 100mu)
// 	if(cellTypeRatioDepths.size())
// 		cellTypeRatioDepths.clear();
// 	cellTypeRatioDepths.push_back(0.929);
// 	cellTypeRatioDepths.push_back(0.714);
// 	cellTypeRatioDepths.push_back(0.571);
// 	cellTypeRatioDepths.push_back(0.429);
// 	cellTypeRatioDepths.push_back(0.357);
// 	cellTypeRatioDepths.push_back(0.289);
// 	cellTypeRatioDepths.push_back(0.214);
// 	cellTypeRatioDepths.push_back(0.143);
// 	cellTypeRatioDepths.push_back(0.071);
// 	cellTypeRatioDepths.push_back(0.0);
// 	cellTypeRatioDepths.push_back(-0.04);
// 	cellTypeRatioDepths.push_back(-0.08);
// 	cellTypeRatioDepths.push_back(-0.16);
// 	cellTypeRatioDepths.push_back(-0.24);
// 	cellTypeRatioDepths.push_back(-0.4);
// 	cellTypeRatioDepths.push_back(-0.48);
// 	cellTypeRatioDepths.push_back(-0.6);
// 	cellTypeRatioDepths.push_back(-0.64);
// 	cellTypeRatioDepths.push_back(-0.8);
	
// // 	Excitatory boundary surfaces:
// // 	13Jan2014 AxonColumn final:
// // 	cell type surfaces scaled 
// // 	separately supra-/granular/infra-;
// // 	as during neuron registration
// 	if(cellTypeRatioDepthsSupra.size())
// 		cellTypeRatioDepthsSupra.clear();
// // 	L1 0-125mu (see also L1-L2 project)
// 	cellTypeRatioDepthsSupra.push_back(0.773);
// // 	L1 0-150mu
// // 	cellTypeRatioDepthsSupra.push_back(0.727);
// 	cellTypeRatioDepthsSupra.push_back(0.455);
// 	cellTypeRatioDepthsSupra.push_back(0.364);
// 	cellTypeRatioDepthsSupra.push_back(0.273);
// 	cellTypeRatioDepthsSupra.push_back(0.182);
// 	cellTypeRatioDepthsSupra.push_back(0.091);
// 	if(cellTypeRatioDepthsGran.size())
// 		cellTypeRatioDepthsGran.clear();
// 	cellTypeRatioDepthsGran.push_back(1.0);
// 	cellTypeRatioDepthsGran.push_back(0.857);
// 	cellTypeRatioDepthsGran.push_back(0.714);
// 	cellTypeRatioDepthsGran.push_back(0.571);
// 	cellTypeRatioDepthsGran.push_back(0.429);
// 	cellTypeRatioDepthsGran.push_back(0.286);
// 	cellTypeRatioDepthsGran.push_back(0.143);
// 	cellTypeRatioDepthsGran.push_back(0.0);
// 	if(cellTypeRatioDepthsInfra.size())
// 		cellTypeRatioDepthsInfra.clear();
// 	cellTypeRatioDepthsInfra.push_back(-0.048);
// 	cellTypeRatioDepthsInfra.push_back(-0.095);
// 	cellTypeRatioDepthsInfra.push_back(-0.143);
// 	cellTypeRatioDepthsInfra.push_back(-0.190);
// 	cellTypeRatioDepthsInfra.push_back(-0.238);
// 	cellTypeRatioDepthsInfra.push_back(-0.286);
// 	cellTypeRatioDepthsInfra.push_back(-0.333);
// 	cellTypeRatioDepthsInfra.push_back(-0.381);
// 	cellTypeRatioDepthsInfra.push_back(-0.429);
// 	cellTypeRatioDepthsInfra.push_back(-0.476);
// 	cellTypeRatioDepthsInfra.push_back(-0.524);
// 	cellTypeRatioDepthsInfra.push_back(-0.571);
// 	cellTypeRatioDepthsInfra.push_back(-0.619);
// 	cellTypeRatioDepthsInfra.push_back(-0.667);
	
// // 	Inhibitory boundary surfaces:
// // 	06May2014 INColumn V4 (Daniel Master Thesis):
// // 	cell type surfaces scaled 
// // 	separately supra-/granular/infra-;
// // 	as during neuron registration
// 	if(cellTypeRatioDepthsSupra.size())
// 		cellTypeRatioDepthsSupra.clear();
// // 	L1 0-125mu (see also L1-L2 project)
// 	cellTypeRatioDepthsSupra.push_back(0.773);
// // 	L1 0-150mu
// // 	cellTypeRatioDepthsSupra.push_back(0.727);
// 	cellTypeRatioDepthsSupra.push_back(0.455);
// 	cellTypeRatioDepthsSupra.push_back(0.273);
// 	cellTypeRatioDepthsSupra.push_back(0.091);
// 	if(cellTypeRatioDepthsGran.size())
// 		cellTypeRatioDepthsGran.clear();
// 	cellTypeRatioDepthsGran.push_back(1.0);
// 	cellTypeRatioDepthsGran.push_back(0.857);
// 	cellTypeRatioDepthsGran.push_back(0.714);
// 	cellTypeRatioDepthsGran.push_back(0.571);
// 	cellTypeRatioDepthsGran.push_back(0.143);
// 	cellTypeRatioDepthsGran.push_back(0.0);
// 	if(cellTypeRatioDepthsInfra.size())
// 		cellTypeRatioDepthsInfra.clear();
// 	cellTypeRatioDepthsInfra.push_back(-0.048);
// 	cellTypeRatioDepthsInfra.push_back(-0.095);
// 	cellTypeRatioDepthsInfra.push_back(-0.143);
// 	cellTypeRatioDepthsInfra.push_back(-0.190);
// 	cellTypeRatioDepthsInfra.push_back(-0.238);
// 	cellTypeRatioDepthsInfra.push_back(-0.286);
// 	cellTypeRatioDepthsInfra.push_back(-0.381);
// 	cellTypeRatioDepthsInfra.push_back(-0.429);
// 	cellTypeRatioDepthsInfra.push_back(-0.476);
// 	cellTypeRatioDepthsInfra.push_back(-0.524);
// 	cellTypeRatioDepthsInfra.push_back(-0.571);
// 	cellTypeRatioDepthsInfra.push_back(-0.619);
// 	cellTypeRatioDepthsInfra.push_back(-0.667);
// 	cellTypeRatioDepthsInfra.push_back(-0.714);
// 	cellTypeRatioDepthsInfra.push_back(-0.762);
// 	cellTypeRatioDepthsInfra.push_back(-0.810);
	
// 	Inhibitory boundary surfaces:
// 	18Dec2014 INColumn V6 (IN paper):
// 	cell type surfaces scaled 
// 	separately supra-/granular/infra-;
// 	as during neuron registration
	if(cellTypeRatioDepthsSupra.size())
		cellTypeRatioDepthsSupra.clear();
// 	L1 0-125mu (see also L1-L2 project)
	cellTypeRatioDepthsSupra.push_back(0.773);
// 	L1 0-150mu
// 	cellTypeRatioDepthsSupra.push_back(0.727);
	cellTypeRatioDepthsSupra.push_back(0.455);
	cellTypeRatioDepthsSupra.push_back(0.273);
	cellTypeRatioDepthsSupra.push_back(0.091);
	if(cellTypeRatioDepthsGran.size())
		cellTypeRatioDepthsGran.clear();
	cellTypeRatioDepthsGran.push_back(1.0);
	cellTypeRatioDepthsGran.push_back(0.857);
	cellTypeRatioDepthsGran.push_back(0.714);
	cellTypeRatioDepthsGran.push_back(0.571);
	cellTypeRatioDepthsGran.push_back(0.143);
	cellTypeRatioDepthsGran.push_back(0.0);
	if(cellTypeRatioDepthsInfra.size())
		cellTypeRatioDepthsInfra.clear();
	cellTypeRatioDepthsInfra.push_back(-0.048);
	cellTypeRatioDepthsInfra.push_back(-0.095);
	cellTypeRatioDepthsInfra.push_back(-0.190);
	cellTypeRatioDepthsInfra.push_back(-0.238);
	cellTypeRatioDepthsInfra.push_back(-0.286);
	cellTypeRatioDepthsInfra.push_back(-0.333);
	cellTypeRatioDepthsInfra.push_back(-0.381);
	cellTypeRatioDepthsInfra.push_back(-0.429);
	cellTypeRatioDepthsInfra.push_back(-0.476);
	cellTypeRatioDepthsInfra.push_back(-0.524);
	cellTypeRatioDepthsInfra.push_back(-0.571);
	cellTypeRatioDepthsInfra.push_back(-0.619);
	cellTypeRatioDepthsInfra.push_back(-0.667);
	cellTypeRatioDepthsInfra.push_back(-0.714);
	cellTypeRatioDepthsInfra.push_back(-0.762);
	cellTypeRatioDepthsInfra.push_back(-0.810);
	
// // 	Inhibitory boundary surfaces:
// // 	17Mar2014 InColumn V3 (Daniel):
// // 	cell type surfaces scaled 
// // 	separately supra-/granular/infra-;
// // 	as during neuron registration
// 	if(cellTypeRatioDepthsSupra.size())
// 		cellTypeRatioDepthsSupra.clear();
// // 	L1 0-125mu (see also L1-L2 project)
// 	cellTypeRatioDepthsSupra.push_back(0.773);
// // 	L1 0-150mu
// // 	cellTypeRatioDepthsSupra.push_back(0.727);
// 	cellTypeRatioDepthsSupra.push_back(0.091);
// 	if(cellTypeRatioDepthsGran.size())
// 		cellTypeRatioDepthsGran.clear();
// 	cellTypeRatioDepthsGran.push_back(1.0);
// 	cellTypeRatioDepthsGran.push_back(0.857);
// 	cellTypeRatioDepthsGran.push_back(0.714);
// 	cellTypeRatioDepthsGran.push_back(0.571);
// 	cellTypeRatioDepthsGran.push_back(0.0);
// 	if(cellTypeRatioDepthsInfra.size())
// 		cellTypeRatioDepthsInfra.clear();
// 	cellTypeRatioDepthsInfra.push_back(-0.048);
// 	cellTypeRatioDepthsInfra.push_back(-0.095);
// 	cellTypeRatioDepthsInfra.push_back(-0.143);
// 	cellTypeRatioDepthsInfra.push_back(-0.190);
// 	cellTypeRatioDepthsInfra.push_back(-0.238);
// 	cellTypeRatioDepthsInfra.push_back(-0.286);
// 	cellTypeRatioDepthsInfra.push_back(-0.429);
// 	cellTypeRatioDepthsInfra.push_back(-0.476);
// 	cellTypeRatioDepthsInfra.push_back(-0.524);
// 	cellTypeRatioDepthsInfra.push_back(-0.571);
// 	cellTypeRatioDepthsInfra.push_back(-0.619);
// 	cellTypeRatioDepthsInfra.push_back(-0.667);
// 	cellTypeRatioDepthsInfra.push_back(-0.714);
// 	cellTypeRatioDepthsInfra.push_back(-0.762);
// 	cellTypeRatioDepthsInfra.push_back(-0.810);
	
// // 	27Nov2013 AxonColumn V1:
// // 	cell type surfaces scaled 
// // 	separately supra-/granular/infra-;
// // 	as during neuron registration
// 	if(cellTypeRatioDepthsSupra.size())
// 		cellTypeRatioDepthsSupra.clear();
// // 	L1 0-125mu (see also L1-L2 project)
// 	cellTypeRatioDepthsSupra.push_back(0.773);
// // 	L1 0-150mu
// // 	cellTypeRatioDepthsSupra.push_back(0.727);
// 	cellTypeRatioDepthsSupra.push_back(0.455);
// 	cellTypeRatioDepthsSupra.push_back(0.364);
// 	cellTypeRatioDepthsSupra.push_back(0.273);
// 	cellTypeRatioDepthsSupra.push_back(0.182);
// 	cellTypeRatioDepthsSupra.push_back(0.091);
// 	if(cellTypeRatioDepthsGran.size())
// 		cellTypeRatioDepthsGran.clear();
// 	cellTypeRatioDepthsGran.push_back(1.0);
// 	cellTypeRatioDepthsGran.push_back(0.857);
// 	cellTypeRatioDepthsGran.push_back(0.714);
// 	cellTypeRatioDepthsGran.push_back(0.571);
// 	cellTypeRatioDepthsGran.push_back(0.429);
// 	cellTypeRatioDepthsGran.push_back(0.286);
// 	cellTypeRatioDepthsGran.push_back(0.0);
// 	if(cellTypeRatioDepthsInfra.size())
// 		cellTypeRatioDepthsInfra.clear();
// 	cellTypeRatioDepthsInfra.push_back(-0.048);
// 	cellTypeRatioDepthsInfra.push_back(-0.095);
// 	cellTypeRatioDepthsInfra.push_back(-0.190);
// 	cellTypeRatioDepthsInfra.push_back(-0.238);
// 	cellTypeRatioDepthsInfra.push_back(-0.286);
// 	cellTypeRatioDepthsInfra.push_back(-0.333);
// 	cellTypeRatioDepthsInfra.push_back(-0.381);
// 	cellTypeRatioDepthsInfra.push_back(-0.429);
// 	cellTypeRatioDepthsInfra.push_back(-0.476);
// 	cellTypeRatioDepthsInfra.push_back(-0.524);
// 	cellTypeRatioDepthsInfra.push_back(-0.571);
// 	cellTypeRatioDepthsInfra.push_back(-0.619);
// 	cellTypeRatioDepthsInfra.push_back(-0.667);
// 	cellTypeRatioDepthsInfra.push_back(-0.714);
	
// 	27Nov2013: AxonColumn V1 (linear interpolation; deprecated)
	if(cellTypeRatioDepths.size())
		cellTypeRatioDepths.clear();
// 	L1 for Arno 0-125mu
	cellTypeRatioDepths.push_back(0.821);
// 	cellTypeRatioDepths.push_back(0.786);
	cellTypeRatioDepths.push_back(0.571);
	cellTypeRatioDepths.push_back(0.500);
	cellTypeRatioDepths.push_back(0.429);
	cellTypeRatioDepths.push_back(0.357);
	cellTypeRatioDepths.push_back(0.286);
	cellTypeRatioDepths.push_back(0.214);
	cellTypeRatioDepths.push_back(0.143);
	cellTypeRatioDepths.push_back(0.071);
	cellTypeRatioDepths.push_back(0.0);
	cellTypeRatioDepths.push_back(-0.04);
	cellTypeRatioDepths.push_back(-0.08);
	cellTypeRatioDepths.push_back(-0.16);
	cellTypeRatioDepths.push_back(-0.20);
	cellTypeRatioDepths.push_back(-0.24);
	cellTypeRatioDepths.push_back(-0.32);
	cellTypeRatioDepths.push_back(-0.36);
	cellTypeRatioDepths.push_back(-0.40);
	cellTypeRatioDepths.push_back(-0.44);
	cellTypeRatioDepths.push_back(-0.48);
	cellTypeRatioDepths.push_back(-0.52);
	cellTypeRatioDepths.push_back(-0.56);
	cellTypeRatioDepths.push_back(-0.60);
	cellTypeRatioDepths.push_back(-0.64);
	cellTypeRatioDepths.push_back(-0.68);
	cellTypeRatioDepths.push_back(-0.72);
	cellTypeRatioDepths.push_back(-0.76);
	
// 	21Jan2013: new, correct(er) version
// 	if(cellTypeRatioDepths.size())
// 		cellTypeRatioDepths.clear();
// // 	L1 for Arno 0-125mu
// 	cellTypeRatioDepths.push_back(0.821);
// // 	cellTypeRatioDepths.push_back(0.786);
// 	cellTypeRatioDepths.push_back(0.571);
// 	cellTypeRatioDepths.push_back(0.429);
// 	cellTypeRatioDepths.push_back(0.286);
// 	cellTypeRatioDepths.push_back(0.214);
// 	cellTypeRatioDepths.push_back(0.143);
// 	cellTypeRatioDepths.push_back(0.071);
// 	cellTypeRatioDepths.push_back(0.0);
// 	cellTypeRatioDepths.push_back(-0.04);
// 	cellTypeRatioDepths.push_back(-0.08);
// 	cellTypeRatioDepths.push_back(-0.12);
// 	cellTypeRatioDepths.push_back(-0.16);
// 	cellTypeRatioDepths.push_back(-0.24);
// 	cellTypeRatioDepths.push_back(-0.32);
// 	cellTypeRatioDepths.push_back(-0.48);
// 	cellTypeRatioDepths.push_back(-0.56);
// 	cellTypeRatioDepths.push_back(-0.68);
// 	cellTypeRatioDepths.push_back(-0.72);
// 	cellTypeRatioDepths.push_back(-0.88);
	
// 	18Jun2013: cell type surfaces scaled 
// 	separately supra-/granular/infra-;
// 	as during neuron registration
// 	if(cellTypeRatioDepthsSupra.size())
// 		cellTypeRatioDepthsSupra.clear();
// 	cellTypeRatioDepthsSupra.push_back(0.727);
// 	cellTypeRatioDepthsSupra.push_back(0.455);
// 	cellTypeRatioDepthsSupra.push_back(0.273);
// 	cellTypeRatioDepthsSupra.push_back(0.091);
// 	if(cellTypeRatioDepthsGran.size())
// 		cellTypeRatioDepthsGran.clear();
// 	cellTypeRatioDepthsGran.push_back(1.0);
// 	cellTypeRatioDepthsGran.push_back(0.857);
// 	cellTypeRatioDepthsGran.push_back(0.714);
// 	cellTypeRatioDepthsGran.push_back(0.571);
// 	cellTypeRatioDepthsGran.push_back(0.429);
// 	cellTypeRatioDepthsGran.push_back(0.286);
// 	cellTypeRatioDepthsGran.push_back(0.143);
// 	cellTypeRatioDepthsGran.push_back(0.0);
// 	if(cellTypeRatioDepthsInfra.size())
// 		cellTypeRatioDepthsInfra.clear();
// 	cellTypeRatioDepthsInfra.push_back(-0.095);
// 	cellTypeRatioDepthsInfra.push_back(-0.190);
// 	cellTypeRatioDepthsInfra.push_back(-0.381);
// 	cellTypeRatioDepthsInfra.push_back(-0.476);
// 	cellTypeRatioDepthsInfra.push_back(-0.619);
// 	cellTypeRatioDepthsInfra.push_back(-0.667);
}








