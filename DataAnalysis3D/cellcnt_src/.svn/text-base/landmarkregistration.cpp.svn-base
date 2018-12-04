#include "landmarkregistration.h"
// #define DEBUG

LandmarkRegistration::LandmarkRegistration()
{

}

LandmarkRegistration::LandmarkRegistration ( std::map< int, Column* > inputBarrelField )
{
	SBF = new BarrelField(false);
	this->inputBarrelField = inputBarrelField;
}

LandmarkRegistration::LandmarkRegistration ( AmiraSpatialGraph* inputSG )
{
	SBF = new BarrelField(false);
	if(inputSG)
	{
		if(inputBarrelField.size()) inputBarrelField.clear();
		
		std::list< int >::const_iterator labelIt;
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			PolyDataPointerType barrel = PolyDataPointerType::New();
			if(inputSG->extractLandmark(ID, barrel))
			{
				double top[3] = {0,0,0}, bottom[3] = {0,0,0};
				
				#ifdef DEBUG
				PointsPointerType barrelPts = barrel->GetPoints();
				for(int ii = 0; ii < barrel->GetCell(0)->GetNumberOfPoints(); ++ii)
				{
					double pt[3];
					barrelPts->GetPoint(barrel->GetCell(0)->GetPointId(ii), pt);
					top[0] += pt[0];
					top[1] += pt[1];
					top[2] += pt[2];
				}
				for(int ii = 0; ii < barrel->GetCell(1)->GetNumberOfPoints(); ++ii)
				{
					double pt[3];
					barrelPts->GetPoint(barrel->GetCell(1)->GetPointId(ii), pt);
					bottom[0] += pt[0];
					bottom[1] += pt[1];
					bottom[2] += pt[2];
				}
				vtkMath::MultiplyScalar(top, 1/double(barrel->GetCell(0)->GetNumberOfPoints()));
				vtkMath::MultiplyScalar(bottom, 1/double(barrel->GetCell(1)->GetNumberOfPoints()));
				#endif
				
				// parametric center:
				double pCenter1[3], pCenter2[3], * weights1, * weights2;
				int subID1, subID2;
				weights1 = new double[barrel->GetCell(0)->GetNumberOfPoints()];
				weights2 = new double[barrel->GetCell(1)->GetNumberOfPoints()];
				barrel->GetCell(0)->GetParametricCenter(pCenter1);
				barrel->GetCell(0)->EvaluateLocation(subID1, pCenter1, top, weights1);
				barrel->GetCell(1)->GetParametricCenter(pCenter2);
				barrel->GetCell(1)->EvaluateLocation(subID2, pCenter2, bottom, weights2);
				delete [] weights1, delete [] weights2;
				
				Column * thisBarrel = new Column(barrel, top, bottom);
				inputBarrelField.insert(std::pair< int, Column * >(ID, thisBarrel));
				
				#ifdef DEBUG
				std::cout << "Barrel " << SBF->int2Labels[ID] << ":" << std::endl;
				barrelPts->Print(std::cout);
				std::cout << "top @ [" << thisBarrel->top[0] << "," << thisBarrel->top[1] << "," << thisBarrel->top[2] << "]" << std::endl;
				std::cout << "bottom @ [" << thisBarrel->bottom[0] << "," << thisBarrel->bottom[1] << "," << thisBarrel->bottom[2] << "]" << std::endl;
				#endif
			}
		}
	}
}

LandmarkRegistration::~LandmarkRegistration()
{
	inputBarrelField.clear();
	if(SBF) delete SBF;
}

void LandmarkRegistration::startRegistration()
{
	std::list< int > commonLandmarks;
	std::map< int, Column * >::const_iterator refBFIt, matchBFIt;
	for(refBFIt = SBF->avgBarrels.begin(); refBFIt != SBF->avgBarrels.end(); ++refBFIt)
		for(matchBFIt = inputBarrelField.begin(); matchBFIt != inputBarrelField.end(); ++matchBFIt)
			if(refBFIt->first == matchBFIt->first)
			{
				commonLandmarks.push_back(refBFIt->first);
				#ifdef DEBUG
				std::cout << ">>>        Found common landmark: " << SBF->int2Labels[refBFIt->first] << std::endl;
				#endif
				break;
			}
	
	std::map< int, Column * > commonRefBarrels;
	std::map< int, Column * > commonInputBarrels;
	std::list< int >::const_iterator commonLandmarksIt;
	for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt)
	{
		int ID = *commonLandmarksIt;
		Column * newRefBarrel = new Column(SBF->avgBarrels[ID]);
		Column * newInputBarrel = new Column(inputBarrelField[ID]);
		commonRefBarrels.insert(std::pair< int, Column * >(ID, newRefBarrel));
		commonInputBarrels.insert(std::pair< int, Column * >(ID, newInputBarrel));
	}
	
	double shift[3], invShift[3];
	alignBarrelFieldCentroids(commonInputBarrels, commonRefBarrels, shift);
	gsl_matrix * mU = computeOptimalRotation(commonRefBarrels, commonInputBarrels);
	if(mU)
	{
		TransformPointerType shiftTrans = TransformPointerType::New();
		TransformPointerType invShiftTrans = TransformPointerType::New();
		shiftTrans->Translate(shift);
		commonRefBarrels[D2]->getCenter(invShift);
		invShift[0] = -invShift[0];
		invShift[1] = -invShift[1];
		invShift[2] = -invShift[2];
		invShiftTrans->Translate(invShift);
		HomogeneousMatrixPointerType regMatrix = gsl2VtkMatrix(mU);
		regTransform = invShiftTrans;
		regTransform->Concatenate(regMatrix);
		regTransform->Concatenate(shiftTrans);
// 		regTransform = shiftTrans;
// 		regTransform->Concatenate(regMatrix);
// 		regTransform->Concatenate(invShiftTrans);
// 		regTransform->Update();
	}
	else
	{
		std::cout << "Error! Registration to Standard Barrel Field failed!" << std::endl;
	}
	
	for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt)
	{
		int ID = *commonLandmarksIt;
		delete commonRefBarrels[ID];
		delete commonInputBarrels[ID];
	}
}

void LandmarkRegistration::startColumnRegistration()
{
	std::list< int > commonLandmarks;
	std::map< int, Column * >::const_iterator refBFIt, matchBFIt;
	for(refBFIt = SBF->avgColumns.begin(); refBFIt != SBF->avgColumns.end(); ++refBFIt)
		for(matchBFIt = inputBarrelField.begin(); matchBFIt != inputBarrelField.end(); ++matchBFIt)
			if(refBFIt->first == matchBFIt->first)
			{
				commonLandmarks.push_back(refBFIt->first);
				#ifdef DEBUG
				std::cout << ">>>        Found common landmark: " << SBF->int2Labels[refBFIt->first] << std::endl;
				#endif
				break;
			}
	
	std::map< int, Column * > commonRefBarrels;
	std::map< int, Column * > commonInputBarrels;
	std::list< int >::const_iterator commonLandmarksIt;
	for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt)
	{
		int ID = *commonLandmarksIt;
		Column * newRefBarrel = new Column(SBF->avgColumns[ID]);
		Column * newInputBarrel = new Column(inputBarrelField[ID]);
		commonRefBarrels.insert(std::pair< int, Column * >(ID, newRefBarrel));
		commonInputBarrels.insert(std::pair< int, Column * >(ID, newInputBarrel));
	}
	
	double shift[3], invShift[3], oldD2ColCenter[3];
	commonRefBarrels[D2]->getCenter(oldD2ColCenter);
	alignBarrelFieldCentroids(commonInputBarrels, commonRefBarrels, shift);
	gsl_matrix * mU = computeOptimalRotation(commonRefBarrels, commonInputBarrels);
	if(mU)
	{
		TransformPointerType shiftTrans = TransformPointerType::New();
		TransformPointerType invShiftTrans = TransformPointerType::New();
		shiftTrans->Translate(shift);
		commonRefBarrels[D2]->getCenter(invShift);
		invShift[0] = oldD2ColCenter[0] - invShift[0];
		invShift[1] = oldD2ColCenter[1] - invShift[1];
		invShift[2] = oldD2ColCenter[2] - invShift[2];
		invShiftTrans->Translate(invShift);
		HomogeneousMatrixPointerType regMatrix = gsl2VtkMatrix(mU);
		regTransform = invShiftTrans;
		regTransform->Concatenate(regMatrix);
		regTransform->Concatenate(shiftTrans);
// 		regTransform = shiftTrans;
// 		regTransform->Concatenate(regMatrix);
// 		regTransform->Concatenate(invShiftTrans);
// 		regTransform->Update();
	}
	else
	{
		std::cout << "Error! Registration to Standard Barrel Field failed!" << std::endl;
	}
	
	for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt)
	{
		int ID = *commonLandmarksIt;
		delete commonRefBarrels[ID];
		delete commonInputBarrels[ID];
	}
}

void LandmarkRegistration::startRegAnisotropicScale()
{
	std::list< int > commonLandmarks;
	std::map< int, Column * >::const_iterator refBFIt, matchBFIt;
	for(refBFIt = SBF->avgBarrels.begin(); refBFIt != SBF->avgBarrels.end(); ++refBFIt)
		for(matchBFIt = inputBarrelField.begin(); matchBFIt != inputBarrelField.end(); ++matchBFIt)
			if(refBFIt->first == matchBFIt->first)
			{
				commonLandmarks.push_back(refBFIt->first);
				#ifdef DEBUG
				std::cout << ">>>        Found common landmark: " << SBF->int2Labels[refBFIt->first] << std::endl;
				#endif
				break;
			}
	
	std::map< int, Column * > commonRefBarrels;
	std::map< int, Column * > commonInputBarrels;
	std::list< int >::const_iterator commonLandmarksIt;
	for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt)
	{
		int ID = *commonLandmarksIt;
		Column * newRefBarrel = new Column(SBF->avgBarrels[ID]);
		Column * newInputBarrel = new Column(inputBarrelField[ID]);
		commonRefBarrels.insert(std::pair< int, Column * >(ID, newRefBarrel));
		commonInputBarrels.insert(std::pair< int, Column * >(ID, newInputBarrel));
	}
	
	// optimal translation
	double shift[3], invShift[3];
	alignBarrelFieldCentroids(commonInputBarrels, commonRefBarrels, shift);
	
	// compute optimal rotation and anisotropic scaling
	TransformPointerType totalTransform = TransformPointerType::New();
	gsl_matrix * mU = gsl_matrix_alloc(3, 3);
	gsl_matrix * mLambda = gsl_matrix_alloc(3, 3);
	gsl_matrix_set_identity(mLambda);
	gsl_matrix * mOpt = gsl_matrix_alloc(3, 3);
	double residuals = getResiduals(commonInputBarrels, commonRefBarrels);
	double tol = 1E-03, delta = 0;
	unsigned long iterCnt = 1;
	std::cout << "*********************************" << std::endl;
	do
	{
		std::cout << std::endl;
		std::cout << "Iteration " << iterCnt << std::endl;
		std::cout << std::endl;
		
		gsl_matrix * mX = createPointPositionMatrix(commonInputBarrels);
		gsl_matrix * mY = createPointPositionMatrix(commonRefBarrels);
		
		computeOptimalRotationScale(mX, mY, mU, mLambda);
		bool zScaleConstraint = 1;
		computeOptimalScale(mX, mY, mU, mLambda, zScaleConstraint);
		
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mU, mLambda, 0.0, mOpt);
		
// 		double * transformation[4];
// 		for(int jj = 0; jj < 4; ++jj)
// 			transformation[jj] = new double[4];
// 		for(int jj = 0; jj < 3; ++jj)
// 			for(int kk = 0; kk < 3; ++kk)
// 				transformation[jj][kk] = gsl_matrix_get(mOpt, jj, kk);
// 		for(int jj = 0; jj < 3; ++jj)
// 		{
// 			transformation[jj][3] = 0;
// 			transformation[3][jj] = 0;
// 		}
// 		transformation[3][3] = 1;
		
		HomogeneousMatrixPointerType mRot = HomogeneousMatrixPointerType::New();
		mRot->Identity();
		for(int jj = 0; jj < 3; ++jj)
			for(int kk = 0; kk < 3; ++kk)
				mRot->SetElement(jj, kk, gsl_matrix_get(mOpt, jj, kk));
		
		if(iterCnt == 1)
		{
			totalTransform->SetMatrix(mRot);
		}
		else
		{
			totalTransform->Concatenate(mRot);
		}
		std::cout << "Current optimal scaling and rotation matrix:" << std::endl;
		totalTransform->GetMatrix()->Print(std::cout);
		std::flush(std::cout);
		
		std::map< int, Column * >::const_iterator matchBarrelIt;
		for(matchBarrelIt = commonInputBarrels.begin(); matchBarrelIt != commonInputBarrels.end(); ++matchBarrelIt)
		{
			matchBarrelIt->second->rotateColumn(mOpt);
		}
		
		double newResiduals = getResiduals(commonInputBarrels, commonRefBarrels);
		delta = fabs(residuals-newResiduals)/residuals;
		
		std::cout << std::endl;
		std::cout << "old residuals =\t" << residuals << std::endl;
		std::cout << "new residuals =\t" << newResiduals << std::endl;
		std::cout << "delta =\t" << delta << std::endl;
		std::cout << "tol =\t" << tol << std::endl;
		std::cout << "*********************************" << std::endl;
		
		residuals = newResiduals;
		++iterCnt;
	} while(delta > tol && iterCnt <= 10);
	
	// inverse translation to D2 origin
	TransformPointerType shiftTrans = TransformPointerType::New();
	TransformPointerType invShiftTrans = TransformPointerType::New();
	shiftTrans->Translate(shift);
	commonRefBarrels[D2]->getCenter(invShift);
	invShift[0] = -invShift[0];
	invShift[1] = -invShift[1];
	invShift[2] = -invShift[2];
	invShiftTrans->Translate(invShift);
	HomogeneousMatrixPointerType regMatrix = totalTransform->GetMatrix();
	regTransform = invShiftTrans;
	regTransform->Concatenate(regMatrix);
	regTransform->Concatenate(shiftTrans);
	
	gsl_matrix_free(mU), gsl_matrix_free(mLambda), gsl_matrix_free(mOpt);
	for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt)
	{
		int ID = *commonLandmarksIt;
		delete commonRefBarrels[ID];
		delete commonInputBarrels[ID];
	}
}

TransformPointerType LandmarkRegistration::getTransform()
{
	if(regTransform) return regTransform;
	else
	{
		TransformPointerType id = TransformPointerType::New();
		id->Identity();
		return id;
	}
}

void LandmarkRegistration::writeTransform ( const char* outputFilename )
{
	if(regTransform)
	{
		std::string ofName(outputFilename);
		ofName += "_registered_transform.log";
		std::ofstream TransformStream;
		TransformStream.open(ofName.c_str());
		
		TransformStream << "# Transformation matrix 4x4" << std::endl;
		for(int ii = 0; ii < 4; ++ii)
		{
			TransformStream << "[";
			for(int jj = 0; jj < 4; ++jj)
			{
				if(jj < 3)
					TransformStream << regTransform->GetMatrix()->GetElement(ii, jj) << " ";
				else
					TransformStream << regTransform->GetMatrix()->GetElement(ii, jj);
			}
			TransformStream << "]" << std::endl;
		}
		TransformStream << std::endl;
		
		TransformStream << "# Transformation matrix Amira format" << std::endl;
		for(int ii = 0; ii < 4; ++ii)
			for(int jj = 0; jj < 4; ++jj)
				TransformStream << regTransform->GetMatrix()->GetElement(jj, ii) << " ";
		TransformStream << std::endl;
		
		TransformStream.close();
	}
}

void LandmarkRegistration::alignBarrelFieldCentroids ( std::map< int, Column* > inputBarrels, std::map< int, Column* > refBarrels, double shift[3] )
{
	double inputShift[] = {0,0,0}, refShift[] = {0,0,0};
	std::map< int, Column * >::const_iterator bfIt;
	for(bfIt = inputBarrels.begin(); bfIt != inputBarrels.end(); ++bfIt)
	{
		for(int jj = 0; jj < 3; ++jj)
		{
			inputShift[jj] += bfIt->second->top[jj];
			inputShift[jj] += bfIt->second->bottom[jj];
		}
	}
	for(bfIt = refBarrels.begin(); bfIt != refBarrels.end(); ++bfIt)
	{
		for(int jj = 0; jj < 3; ++jj)
		{
			refShift[jj] += bfIt->second->top[jj];
			refShift[jj] += bfIt->second->bottom[jj];
		}
	}
	for(int jj = 0; jj < 3; ++jj)
	{
		if(inputBarrels.size())
		{
			inputShift[jj] = -inputShift[jj]/(double)(2*inputBarrels.size());
			shift[jj] = inputShift[jj];
		}
		if(refBarrels.size())
		{
			refShift[jj] = -refShift[jj]/(double)(2*refBarrels.size());
		}
	}
	for(bfIt = inputBarrels.begin(); bfIt != inputBarrels.end(); ++bfIt)
		bfIt->second->translateColumn(inputShift);
	for(bfIt = refBarrels.begin(); bfIt != refBarrels.end(); ++bfIt)
		bfIt->second->translateColumn(refShift);
}

gsl_matrix* LandmarkRegistration::computeOptimalRotation ( std::map< int, Column* > refBF, std::map< int, Column* > matchBF )
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
	
	#ifdef DEBUG
	std::cout << "mCov:" << status << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
	}
	#endif
	
	// compute SVD of mCov
	gsl_vector * vTmp = gsl_vector_alloc(3);
	status = gsl_linalg_SV_decomp(mCov, mV, vS, vTmp);
	
	#ifdef DEBUG
	std::cout << "mCov after SVD: " << status << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
	}
	std::cout << "mV after SVD: " << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mV, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mV, ii, 2) << "]" << std::endl;
	}
	std::cout << "vS after SVD: " << std::endl;
	std::cout << "[" << gsl_vector_get(vS, 0) << "," << gsl_vector_get(vS, 1) << "," << gsl_vector_get(vS, 2) << "]" << std::endl;
	#endif
	
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
	
	#ifdef DEBUG
	std::cout << "mUTmp: " << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mUTmp, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mUTmp, ii, 2) << "]" << std::endl;
	}
	std::cout << "mVTmp: " << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mVTmp, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mVTmp, ii, 2) << "]" << std::endl;
	}
	std::cout << "detSign = " << detSign << std::endl;
	#endif
	
	// compute mU = mV * diag(1,1,detSign) * mCov^t
	// note: gsl already computes mV so that mCov = U * S * V^t
	gsl_matrix_set_identity(mId);
	gsl_matrix_set(mId, 2, 2, detSign);
	gsl_matrix * mTmp = gsl_matrix_alloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mId, mCov, 0.0, mTmp);
	status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mV, mTmp, 0.0, mU);
	
	#ifdef DEBUG
	std::cout << "mU: " << status << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mU, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mU, ii, 2) << "]" << std::endl;
	}
	#endif
	
	gsl_matrix_free(mUTmp), gsl_matrix_free(mVTmp);
	gsl_permutation_free(permU), gsl_permutation_free(permV);
	gsl_matrix_free(mTmp), gsl_matrix_free(mCov), gsl_matrix_free(mId), gsl_matrix_free(mLU), gsl_matrix_free(mV), gsl_matrix_free(mX), gsl_matrix_free(mY);
	gsl_permutation_free(permLU);
	gsl_vector_free(vS), gsl_vector_free(vTmp);
	
	return mU;
}

void LandmarkRegistration::computeOptimalRotationScale(gsl_matrix* mX, gsl_matrix* mY, gsl_matrix* mU, gsl_matrix* mLambda)
{
	gsl_matrix * mCov = gsl_matrix_alloc(3, 3);
	gsl_matrix * mTmp = gsl_matrix_alloc(3, 3);
	gsl_matrix * mLU = gsl_matrix_alloc(3, 3);
	gsl_matrix * mV = gsl_matrix_alloc(3, 3);
	gsl_matrix * mId = gsl_matrix_alloc(3, 3);
	gsl_vector * vS = gsl_vector_alloc(3);
	gsl_permutation * permLU = gsl_permutation_alloc(3);
	if(!mX || !mY || !mCov || !mLU || !mU || !mV || !mId || !vS || !permLU)
	{
		std::cout << "Error! Could not allocate enough memory for optimal rotation. Aborting..." << std::endl;
		return;
	}
	
	int status;
	// compute Lambda X Y^t
	status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mX, mY, 0.0, mTmp);
	status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mLambda, mTmp, 0.0, mCov);
	
	#ifdef DEBUG
	std::cout << "mCov:" << status << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
	}
	#endif
	
	// compute SVD of mCov
	gsl_vector * vTmp = gsl_vector_alloc(3);
	status = gsl_linalg_SV_decomp(mCov, mV, vS, vTmp);
	
	#ifdef DEBUG
	std::cout << "mCov after SVD: " << status << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
	}
	std::cout << "mV after SVD: " << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mV, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mV, ii, 2) << "]" << std::endl;
	}
	std::cout << "vS after SVD: " << std::endl;
	std::cout << "[" << gsl_vector_get(vS, 0) << "," << gsl_vector_get(vS, 1) << "," << gsl_vector_get(vS, 2) << "]" << std::endl;
	#endif
	
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
	
	#ifdef DEBUG
	std::cout << "mUTmp: " << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mUTmp, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mUTmp, ii, 2) << "]" << std::endl;
	}
	std::cout << "mVTmp: " << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mVTmp, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mVTmp, ii, 2) << "]" << std::endl;
	}
	std::cout << "detSign = " << detSign << std::endl;
	#endif
	
	// compute mU = mV * diag(1,1,detSign) * mCov^t
	// note: gsl already computes mV so that mCov = U * S * V^t
	gsl_matrix_set_identity(mId);
	gsl_matrix_set(mId, 2, 2, detSign);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mId, mCov, 0.0, mTmp);
	status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mV, mTmp, 0.0, mU);
	
	#ifdef DEBUG
	std::cout << "mU: " << status << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mU, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mU, ii, 2) << "]" << std::endl;
	}
	#endif
	
	gsl_matrix_free(mUTmp), gsl_matrix_free(mVTmp);
	gsl_permutation_free(permU), gsl_permutation_free(permV);
	gsl_matrix_free(mTmp), gsl_matrix_free(mCov), gsl_matrix_free(mId), gsl_matrix_free(mLU), gsl_matrix_free(mV);
	gsl_permutation_free(permLU);
	gsl_vector_free(vS), gsl_vector_free(vTmp);
}

void LandmarkRegistration::computeOptimalScale(gsl_matrix* mX, gsl_matrix* mY, gsl_matrix* mU, gsl_matrix* mLambda, bool constrained)
{
	gsl_matrix * mCov = gsl_matrix_alloc(3, 3);
	gsl_matrix * mCovX = gsl_matrix_alloc(3, 3);
	gsl_matrix * mTmp = gsl_matrix_alloc(3, 3);
	if(!mCov || !mCovX || !mTmp)
	{
		std::cout << "Error! Could not allocate enough memory for optimal rotation. Aborting..." << std::endl;
		return;
	}
	
	int status;
	// compute X Y^t U
	status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mX, mY, 0.0, mTmp);
	status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mTmp, mU, 0.0, mCov);
	
	#ifdef DEBUG
	std::cout << "mCov:" << status << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mCov, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mCov, ii, 2) << "]" << std::endl;
	}
	#endif
	
	// compute X X^t
	status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mX, mX, 0.0, mCovX);
	
	#ifdef DEBUG
	std::cout << "mCovX:" << status << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mCovX, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mCovX, ii, 2) << "]" << std::endl;
	}
	#endif
	
	// compute lambda_ii
	double num, denom, lambda;
	for(int ii = 0; ii < 3; ++ii)
	{
		num = gsl_matrix_get(mCov, ii, ii);
		denom = gsl_matrix_get(mCovX, ii, ii);
		lambda = num/denom;
		gsl_matrix_set(mLambda, ii, ii, lambda);
	}
	// constraint for shrinkage correction:
	// not necessary in z-direction
	if(constrained)
	{
		gsl_matrix_set(mLambda, 2, 2, 1.0);
	}
	
	#ifdef DEBUG
	std::cout << "mLambda:" << std::endl;
	for(int ii = 0; ii < 3; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 2; ++jj)
			std::cout << gsl_matrix_get(mLambda, ii, jj) << ", ";
		std::cout << gsl_matrix_get(mLambda, ii, 2) << "]" << std::endl;
	}
	#endif
	
	gsl_matrix_free(mTmp), gsl_matrix_free(mCov), gsl_matrix_free(mCovX);
}

gsl_matrix* LandmarkRegistration::createPointPositionMatrix(std::map< int, Column* > barrels)
{
	gsl_matrix * mX = gsl_matrix_alloc(3, 2*barrels.size());
	std::map< int, Column* >::const_iterator barrelIt;
	int ii = 0;
	for(barrelIt = barrels.begin(); barrelIt != barrels.end(); ++barrelIt, ++ii)
	{
		double * top, * bottom;
		top = barrelIt->second->top, bottom = barrelIt->second->bottom;
		for(int jj = 0; jj < 3; ++jj)
		{
			double * mXTopPtr = gsl_matrix_ptr(mX, jj, 2*ii);
			double * mXBottomPtr = gsl_matrix_ptr(mX, jj, 2*ii+1);
			if(!mXTopPtr || !mXBottomPtr)
			{
				std::cout << "Error! Invalid memory access during landmark matrix allocation. Aborting..." << std::endl;
				return NULL;
			}
			*mXTopPtr = top[jj];
			*mXBottomPtr = bottom[jj];
		}
	}
	
	return mX;
}

double LandmarkRegistration::getResiduals(std::map< int, Column* > barrels1, std::map< int, Column* > barrels2)
{
	double res = 0;
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrels1.find(ID) != barrels1.end() && barrels2.find(ID) != barrels2.end())
		{
			res += vtkMath::Distance2BetweenPoints(barrels1[ID]->top, barrels2[ID]->top);
			res += vtkMath::Distance2BetweenPoints(barrels1[ID]->bottom, barrels2[ID]->bottom);
		}
	}
	
	return res;
};

HomogeneousMatrixPointerType LandmarkRegistration::gsl2VtkMatrix ( gsl_matrix* mIn )
{
	HomogeneousMatrixPointerType mOut = HomogeneousMatrixPointerType::New();
	mOut->Identity();
	for(int ii = 0; ii < 3; ++ii)
		for(int jj = 0; jj < 3; ++jj)
			mOut->SetElement(ii, jj, gsl_matrix_get(mIn, ii, jj));
	return mOut;
}
