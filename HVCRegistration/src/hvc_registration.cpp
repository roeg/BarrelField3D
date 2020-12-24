/****************************************************************************/
/*                                                                          */
/* Program:   3DCellCountAnalysis                                           */
/*                                                                          */
/* File:      pipeline.cpp                                                  */
/*                                                                          */
/* Purpose:   pipeline for analysis of cell counts with respect to          */
/*            barreloids in VPM reconstructed in 3D                         */
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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "../../common/typedefs.h"
#include "../../common/amiraReader.h"
#include "../../common/inputcheckpoint.h"
#include "../../common/inputparameters.h"
#include "hvcreconstruction.h"

// measure length along "cardinal" anatomical axes (in aligned coordinates)
std::vector< double > getCardinalAxesLengths(PolyDataPointerType alignedSurface, double measuringPoint[3]);

int main( int argc , char * argv[])
{
	if(argc == 3)
	{
// 		const char * surfaceFilename = argv[1];
// 		const char * outputFilename = argv[2];
		
// 		PolyDataPointerType hvcSurface;
// 		
// 		std::string surfaceStr(surfaceFilename);
// 		Reader * hvcSurfaceReader = new Reader(surfaceFilename, outputFilename);
// 		if(surfaceStr.find(".surf") != std::string::npos)
// 			hvcSurface = hvcSurfaceReader->readAmiraSurfaceFile();
// 		else
// 		{
// 			std::cout << "Error! Surface file has to be Amira '.surf' ascii file!" << std::endl;
// 			return 0;
// 		}
		
// 		// test barreloid surface reconstruction:
// 		PolyDataPointerType C1barreloid = PolyDataPointerType::New();
// 		if(landmarkSG->extractLandmark(E3, C1barreloid))
// 		{
// 			BarreloidSurfaceReconstruct * barreloidReconstruct = new BarreloidSurfaceReconstruct(checkpoint->getParameters(), C1barreloid);
// 			PolyDataPointerType C1barreloidSurface = barreloidReconstruct->surfaceReconstruction(E3);
// 			Reader * surfWriter = new Reader(outputFilename, outputFilename);
// 			surfWriter->writeAmiraSurfaceFile(C1barreloidSurface);
// 			delete surfWriter;
// 			delete barreloidReconstruct;
// 		}
		
// 		// smooth surface
// 		LowpassPolyDataFilterType smoothingFilter = LowpassPolyDataFilterType::New();
// 		smoothingFilter->BoundarySmoothingOff();
// 		smoothingFilter->FeatureEdgeSmoothingOff();
// 		smoothingFilter->NormalizeCoordinatesOn();
// 		smoothingFilter->SetNumberOfIterations(10);
// 		smoothingFilter->SetPassBand(0.01);
// 		smoothingFilter->SetInput(hvcSurface);
// 		smoothingFilter->Update();
// 		PolyDataPointerType hvcSurfaceSmooth = smoothingFilter->GetOutput();
		
		const char * contourFilename = argv[1];
		const char * outputFilename = argv[2];
		AmiraSpatialGraph * contourSG;
		std::string contourStr(contourFilename);
		Reader * contourFileReader = new Reader(contourFilename, outputFilename);
		if(contourStr.find(".am") != std::string::npos)
		{
			contourFileReader->readSpatialGraphFile(0);
			contourSG = contourFileReader->getSpatialGraph();
		}
		else
		{
			std::cout << "Error! Landmark file has to be Amira '.am' file!" << std::endl;
			return 0;
		}
		
		// check if input is ok and set flags
		InputCheckpoint * checkpoint = new InputCheckpoint(contourSG);
		checkpoint->run();
// 		if(!checkpoint->getInputOK())
// 		{
// 			std::cout << "Error! Landmark file corrupt. Aborting..." << std::endl;
// 			delete checkpoint;
// 			delete contourFileReader;
// 			delete contourSG;
// 			return 0;
// 		}
		
		std::flush(std::cout << "Starting HVC surface reconstruction..." << std::endl);
		// test ellipsoid approximation:
// 		HVCSurfaceReconstruct * hvcReconstruct = new HVCSurfaceReconstruct(hvcSurfaceSmooth);
		HVCSurfaceReconstruct * hvcReconstruct = new HVCSurfaceReconstruct(checkpoint->getParameters(), contourSG);
		PolyDataPointerType hvcSurface = PolyDataPointerType::New();
		hvcSurface = hvcReconstruct->surfaceReconstruction(D2);
		std::vector< double > HVCCenterVec = hvcReconstruct->getEllipsoidCenter();
		std::vector< double > HVCAxesLengthsVec = hvcReconstruct->getHalfAxisLengths();
		std::vector< std::vector< double > > HVCAxesVec = hvcReconstruct->getPrincipalAxes();
		
		double HVCCenterShift[3];
		double HVCAxesLengths[3];
		double HVCAxis1[3], HVCAxis2[3], HVCAxis3[3];
		double tmpAxis1[3], tmpAxis2[3], tmpAxis3[3];
		int MLIndex = -1, APIndex = -1, DVIndex = -1;
		double MLTmpLength = 0, APTmpLength = 0, DVTmpLength = 0;
		for(int i = 0; i < 3; ++i)
		{
			if(fabs(HVCAxesVec[i][2]) > MLTmpLength)
			{
				MLTmpLength = fabs(HVCAxesVec[i][2]);
				MLIndex = i;
			}
			if(fabs(HVCAxesVec[i][0]) > APTmpLength)
			{
				APTmpLength = fabs(HVCAxesVec[i][0]);
				APIndex = i;
			}
			if(fabs(HVCAxesVec[i][1]) > DVTmpLength)
			{
				DVTmpLength = fabs(HVCAxesVec[i][1]);
				DVIndex = i;
			}
		}
		std::cout << "HVCAxesVec[0][0]: " << HVCAxesVec[0][0] << std::endl;
		std::cout << "HVCAxesVec[1][0]: " << HVCAxesVec[1][0] << std::endl;
		std::cout << "HVCAxesVec[2][0]: " << HVCAxesVec[2][0] << std::endl;
		std::cout << "HVCAxesVec[0][1]: " << HVCAxesVec[0][1] << std::endl;
		std::cout << "HVCAxesVec[1][1]: " << HVCAxesVec[1][1] << std::endl;
		std::cout << "HVCAxesVec[2][1]: " << HVCAxesVec[2][1] << std::endl;
		std::cout << "HVCAxesVec[0][2]: " << HVCAxesVec[0][2] << std::endl;
		std::cout << "HVCAxesVec[1][2]: " << HVCAxesVec[1][2] << std::endl;
		std::cout << "HVCAxesVec[2][2]: " << HVCAxesVec[2][2] << std::endl;
		std::cout << "MLIndex: " << MLIndex << std::endl;
		std::cout << "MLTmpLength: " << MLTmpLength << std::endl;
		std::cout << "APIndex: " << APIndex << std::endl;
		std::cout << "APTmpLength: " << APTmpLength << std::endl;
		std::cout << "DVIndex: " << DVIndex << std::endl;
		std::cout << "DVTmpLength: " << DVTmpLength << std::endl;
		for(int i = 0; i < 3; ++i)
		{
			HVCCenterShift[i] = -HVCCenterVec[i];
			HVCAxesLengths[i] = HVCAxesLengthsVec[i];
			tmpAxis1[i] = HVCAxesVec[MLIndex][i];
			tmpAxis2[i] = HVCAxesVec[APIndex][i];
			tmpAxis3[i] = HVCAxesVec[DVIndex][i];
		}
		// Pick direction for eigenvectors
		// aligned with Sam's choice of global coordinate system
		// M-L in reconstruction: negative z -> lateral; final orientaion: positive x: lateral
		double axis1Sign = tmpAxis1[2] > 0 ? -1 : 1;
		// A-P in reconstruction: positive x -> anterior; final orientaion: positive y: anterior
		double axis2Sign = tmpAxis2[0] > 0 ? 1 : -1;
		// D-V in reconstruction: negative y -> dorsal; final orientaion: positive z: dorsal
		double axis3Sign = tmpAxis3[1] > 0 ? -1 : 1;
// 		double axis1Sign = tmpAxis1[0] > 0 ? 1 : -1;
// 		double axis2Sign = tmpAxis2[1] > 0 ? 1 : -1;
// 		double axis3Sign = tmpAxis3[2] > 0 ? 1 : -1;
		std::cout << "Axis 1 sign: " << axis1Sign << std::endl;
		std::cout << "Axis 2 sign: " << axis2Sign << std::endl;
		std::cout << "Axis 3 sign: " << axis3Sign << std::endl;
		for(int i = 0; i < 3; ++i)
		{
			HVCAxis1[i] = axis1Sign*tmpAxis1[i];
			HVCAxis2[i] = axis2Sign*tmpAxis2[i];
			HVCAxis3[i] = axis3Sign*tmpAxis3[i];
		}
		
		// Transformation:
		// Step 1: shift center to origin
		// Step 2: align x' with x (calculate <x,x'> and x x x' for angle and rotation axis)
		// Step 3: align z'' with z, where z'' is z' transformed into the x', x - aligned coordinate system
		double xAxis[] = {1,0,0};
		double yAxis[] = {0,1,0};
		double zAxis[] = {0,0,1};
		// hack for 106 (dorsal-ventral extent wider than rostral-caudal)
// 		if(surfaceStr.find("106") != std::string::npos)
// 		{
// 			yAxis[0] = 0;
// 			yAxis[1] = 0;
// 			yAxis[2] = 1;
// 			zAxis[0] = 0;
// 			zAxis[1] = 1;
// 			zAxis[2] = 0;
// 		}
		
		// Step 1
		TransformPointerType HVCShift = TransformPointerType::New();
		HVCShift->Translate(HVCCenterShift);
		// Step 2
		double angle1, axis1[3];
		angle1 = -acos(HVCAxis1[0])*180/PI;
		vtkMath::Cross(xAxis, HVCAxis1, axis1);
		TransformPointerType alignXAxis = TransformPointerType::New();
		alignXAxis->RotateWXYZ(angle1, axis1);
		// Step 3
		double zDoublePrime[3];
		alignXAxis->TransformPoint(HVCAxis3, zDoublePrime);
		double angle2, axis2[3];
		angle2 = -acos(zDoublePrime[2])*180/PI;
		vtkMath::Cross(zAxis, zDoublePrime, axis2);
// 		vtkMath::Cross(zDoublePrime, zAxis, axis2);
		TransformPointerType alignZDoublePrimeAxis = TransformPointerType::New();
		alignZDoublePrimeAxis->RotateWXYZ(angle2, axis2);
		// Concatenate
		
		TransformPointerType hvcTransform = TransformPointerType::New();
		hvcTransform = alignZDoublePrimeAxis;
		hvcTransform->Concatenate(alignXAxis);
		hvcTransform->Concatenate(HVCShift);
		hvcTransform->Update();
		
		// transform surface and measure extent 
		// along x/y/z axes
		TransformFilterType alignSurface = TransformFilterType::New();
		alignSurface->SetTransform(hvcTransform);
// 		alignSurface->SetInput(hvcSurfaceSmooth);
		alignSurface->SetInput(hvcSurface);
		alignSurface->Update();
		PolyDataPointerType transformedHVCSurfaceData = alignSurface->GetOutput();
		
		double xAxisNeg[] = {-1,0,0};
		double yAxisNeg[] = {0,-1,0};
		double zAxisNeg[] = {0,0,-1};
		double origin[] = {0,0,0};
		double xPoint1[3];
		double xPoint2[3];
		double yPoint1[3];
		double yPoint2[3];
		double zPoint1[3];
		double zPoint2[3];
		
		Surface * transformedHVCSurface = new Surface(transformedHVCSurfaceData);
		transformedHVCSurface->intersectLineInDirection(xAxis, origin);
		transformedHVCSurface->getLastIntersectPoint(xPoint1);
		transformedHVCSurface->intersectLineInDirection(xAxisNeg, origin);
		transformedHVCSurface->getLastIntersectPoint(xPoint2);
		transformedHVCSurface->intersectLineInDirection(yAxis, origin);
		transformedHVCSurface->getLastIntersectPoint(yPoint1);
		transformedHVCSurface->intersectLineInDirection(yAxisNeg, origin);
		transformedHVCSurface->getLastIntersectPoint(yPoint2);
		transformedHVCSurface->intersectLineInDirection(zAxis, origin);
		transformedHVCSurface->getLastIntersectPoint(zPoint1);
		transformedHVCSurface->intersectLineInDirection(zAxisNeg, origin);
		transformedHVCSurface->getLastIntersectPoint(zPoint2);
		double MLDistance = sqrt(vtkMath::Distance2BetweenPoints(xPoint1, xPoint2));
		double APDistance = sqrt(vtkMath::Distance2BetweenPoints(yPoint1, yPoint2));
		double DVDistance = sqrt(vtkMath::Distance2BetweenPoints(zPoint1, zPoint2));
		
		std::vector< double > cardinalExtent = getCardinalAxesLengths(transformedHVCSurfaceData, origin);
		
		AmiraSpatialGraph * regSG = new AmiraSpatialGraph;
		regSG->mergeSpatialGraph(contourSG);
		regSG->removeLabel(D2);
		regSG->setTransformation(hvcTransform);
		regSG->applyTransformation();
		std::string sgOutName(outputFilename);
		sgOutName += "_registered_SpatialGraph";
		Reader * regSGWriter = new Reader(sgOutName.c_str(), sgOutName.c_str());
		regSGWriter->setSpatialGraph(regSG);
		regSGWriter->writeSpatialGraphFile();
		delete regSGWriter;
		
		// Output
		std::string ofName(outputFilename);
		ofName += "_registered_transform.log";
		std::ofstream TransformStream;
		TransformStream.open(ofName.c_str());
		TransformStream << "# Transformation matrix Amira format" << std::endl;
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j)
				TransformStream << hvcTransform->GetMatrix()->GetElement(j, i) << " ";
		TransformStream << std::endl;
		
		TransformStream << "# Principal axes lengths" << std::endl;
		TransformStream << HVCAxesLengths[0] << std::endl;
		TransformStream << HVCAxesLengths[1] << std::endl;
		TransformStream << HVCAxesLengths[2] << std::endl;
		TransformStream << std::endl;
		
		TransformStream << "# Extent along principal axes through origin" << std::endl;
		TransformStream << "M-L\t" << MLDistance << std::endl;
		TransformStream << "A-P\t" << APDistance << std::endl;
		TransformStream << "D-V\t" << DVDistance << std::endl;
		TransformStream << std::endl;
		
// 		distances.push_back(APHorizontalDistance);
// 		distances.push_back(APSagittalDistance);
// 		distances.push_back(MLHorizontalDistance);
// 		distances.push_back(MLCoronalDistance);
// 		distances.push_back(DVSagittalDistance);
// 		distances.push_back(DVCoronalDistance);
		TransformStream << "# Maximum extent along canonical axes" << std::endl;
		TransformStream << "M-L Horizontal\t" << cardinalExtent[2] << std::endl;
		TransformStream << "M-L Coronal\t" << cardinalExtent[3] << std::endl;
		TransformStream << "A-P Horizontal\t" << cardinalExtent[0] << std::endl;
		TransformStream << "A-P Sagittal\t" << cardinalExtent[1] << std::endl;
		TransformStream << "D-V Sagittal\t" << cardinalExtent[4] << std::endl;
		TransformStream << "D-V Coronal\t" << cardinalExtent[5] << std::endl;
		
		std::string ofName2(outputFilename);
		ofName2 += "_registered_surface";
		Reader * surfaceWriter = new Reader(ofName2.c_str(), ofName2.c_str());
		surfaceWriter->writeAmiraSurfaceFile(transformedHVCSurfaceData);
		
// 		delete hvcReconstruct, delete hvcSurfaceReader, delete surfaceWriter;
		delete hvcReconstruct, delete contourFileReader, delete contourSG, delete surfaceWriter;
	}
	
	// enforce axon/dendrites within canonical HVC surface constraint
	if(argc == 4)
	{
		const char * spatialGraphFilename = argv[1];
		const char * averageHVCFilename = argv[2];
		const char * outputFilename = argv[3];
	}
	
	return 0;
}

std::vector< double > getCardinalAxesLengths(PolyDataPointerType alignedSurface, double measuringPoint[3])
{
	// M-L: x-axis; measure in both coronal and horizontal plane
	// A-P: y-axis; measure in both horizontal and sagittal plane
	// D-V: z-axis; measure in both sagittal and coronal plane
	
	double MLAxis[] = {1,0,0};
	double APAxis[] = {0,1,0};
	double DVAxis[] = {0,0,1};
	
	PlanePointerType coronalPlane = PlanePointerType::New();
	coronalPlane->SetOrigin(measuringPoint);
	coronalPlane->SetNormal(APAxis);
	
	PlanePointerType sagittalPlane = PlanePointerType::New();
	sagittalPlane->SetOrigin(measuringPoint);
	sagittalPlane->SetNormal(MLAxis);
	
	PlanePointerType horizontalPlane = PlanePointerType::New();
	horizontalPlane->SetOrigin(measuringPoint);
	horizontalPlane->SetNormal(DVAxis);
	
	CutterPointerType cutSurfaceCoronalFilter = CutterPointerType::New();
	cutSurfaceCoronalFilter->SetCutFunction(coronalPlane);
	cutSurfaceCoronalFilter->SetInput(alignedSurface);
	cutSurfaceCoronalFilter->Update();
	PolyDataPointerType coronalSurfaceCut = cutSurfaceCoronalFilter->GetOutput();
// 	coronalSurfaceCut->Print(std::cout);
	
	CutterPointerType cutSurfaceSagittalFilter = CutterPointerType::New();
	cutSurfaceSagittalFilter->SetCutFunction(sagittalPlane);
	cutSurfaceSagittalFilter->SetInput(alignedSurface);
	cutSurfaceSagittalFilter->Update();
	PolyDataPointerType sagittalSurfaceCut = cutSurfaceSagittalFilter->GetOutput();
	
	CutterPointerType cutSurfaceHorizontalFilter = CutterPointerType::New();
	cutSurfaceHorizontalFilter->SetCutFunction(horizontalPlane);
	cutSurfaceHorizontalFilter->SetInput(alignedSurface);
	cutSurfaceHorizontalFilter->Update();
	PolyDataPointerType horizontalSurfaceCut = cutSurfaceHorizontalFilter->GetOutput();
	
	double coronalBounds[6];
	coronalSurfaceCut->GetBounds(coronalBounds);
	double MLCoronalDistance = coronalBounds[1] - coronalBounds[0];
	double DVCoronalDistance = coronalBounds[5] - coronalBounds[4];
	
	double sagittalBounds[6];
	sagittalSurfaceCut->GetBounds(sagittalBounds);
	double APSagittalDistance = sagittalBounds[3] - sagittalBounds[2];
	double DVSagittalDistance = sagittalBounds[5] - sagittalBounds[4];
	
	double horizontalBounds[6];
	horizontalSurfaceCut->GetBounds(horizontalBounds);
	double MLHorizontalDistance = horizontalBounds[1] - horizontalBounds[0];
	double APHorizontalDistance = horizontalBounds[3] - horizontalBounds[2];
	
	std::vector< double > distances;
	distances.push_back(APHorizontalDistance);
	distances.push_back(APSagittalDistance);
	distances.push_back(MLHorizontalDistance);
	distances.push_back(MLCoronalDistance);
	distances.push_back(DVSagittalDistance);
	distances.push_back(DVCoronalDistance);
	
	return distances;
}

