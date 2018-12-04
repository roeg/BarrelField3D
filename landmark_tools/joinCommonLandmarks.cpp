/****************************************************************************/
/*                                                                          */
/* Program:   JoinCommonLandmarks                                           */
/*                                                                          */
/* File:      joinCommonLandmarks.cpp                                       */
/*                                                                          */
/* Purpose:   calculates correspondence between two landmark sets           */
/*            and keeps all landmarks present within tolerance radius       */
/*                                                                          */
/* EMail:     Robert.Egger@mpfi.org                                         */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "../common/typedefs.h"
#include "../common/amiraReader.h"
#include "../common/basics.h"

PointsPointerType assignCorrespondingPoints(PointsPointerType refPoints, PointsPointerType matchPoints, double tolRadius);
PointsPointerType mergePoints(PointsPointerType refPoints, PointsPointerType matchPoints, double tolRadius);

int main( int argc , char * argv[])
{
	if(argc == 5)
	{
		const char * inputFilename1 = argv[1];
		const char * inputFilename2 = argv[2];
		const char * outputFilename = argv[3];
		double tolRadius = atof(argv[4]);
		std::string outName(outputFilename);
		outName += "_correspondence_radius_";
		outName += argv[4];
		std::string outName1(outName);
		outName1 += "_set1-2";
		Reader * landmarkReader1 = new Reader(inputFilename1, outName1.c_str());
		PointsPointerType matchPoints1 = landmarkReader1->readLandmarkFile();
		std::string outName2(outName);
		outName2 += "_set1-3";
		Reader * landmarkReader2 = new Reader(inputFilename2, outName2.c_str());
		PointsPointerType matchPoints2 = landmarkReader2->readLandmarkFile();
		
		std::cout << "Set 1 size" << "\t" << matchPoints1->GetNumberOfPoints() << std::endl;
		std::cout << "Set 2 size" << "\t" << matchPoints2->GetNumberOfPoints() << std::endl;
		
		PointsPointerType correspondenceSet1Set2 = assignCorrespondingPoints(matchPoints1, matchPoints2, tolRadius);
		
		std::string summaryName(outName);
		summaryName += "_summary.csv";
		std::ofstream summaryOutStream(summaryName.c_str());
		summaryOutStream << "Correspondence radius (microns)" << "\t" << tolRadius << std::endl;
		summaryOutStream << "Set 1 size" << "\t" << matchPoints1->GetNumberOfPoints() << std::endl;
		summaryOutStream << "Set 2 size" << "\t" << matchPoints2->GetNumberOfPoints() << std::endl;
		summaryOutStream << "Correspondence set 1/2" << "\t" << correspondenceSet1Set2->GetNumberOfPoints() << std::endl;
		summaryOutStream.close();
		
		landmarkReader1->writeLandmarkFile(correspondenceSet1Set2);
		
		delete landmarkReader1, delete landmarkReader2;
	}
	else
		std::cout << "Usage: ./LandmarkCorrespondence [landmark set 1 filename] [landmark set 2 filename] [output filename] [tolerance radius]" << std::endl;
	return 0;
	
// 	if(argc >= 4)
// 	{
// 		std::vector< const char * > inputFilenameVec;
// 		for(int i = 1; i < argc - 2; ++i)
// 		{
// 			inputFilenameVec.push_back(argv[i]);
// 		}
// 		const char * outputFilename = argv[argc-2];
// 		double tolRadius = atof(argv[argc-1]);
// 		std::string outName(outputFilename);
// 		outName += "_correspondence_radius_";
// 		outName += argv[argc-1];
// 		std::vector< PointsPointerType > pointSetVec;
// 		for(int i = 1; i < argc - 2; ++i)
// 		{
// 			Reader * landmarkReader = new Reader(inputFilenameVec[i-1], inputFilenameVec[i-1]);
// 			PointsPointerType pointSet = landmarkReader->readLandmarkFile();
// 			pointSetVec.push_back(pointSet);
// 			std::cout << "Set " << i << " size" << "\t" << pointSet->GetNumberOfPoints() << std::endl;
// 			delete landmarkReader;
// 		}
// 		
// 		std::vector< PointsPointerType > correspondingSets;
// 		for(int i = 0; i < pointSetVec.size() - 1; ++i)
// 			for(int j = i + 1; j < pointSetVec.size(); ++j)
// 			{
// 				PointsPointerType correspondence = assignCorrespondingPoints(pointSetVec[i], pointSetVec[j], tolRadius);
// 				std::cout << "Correspondence set " << i << "/" << j << "\t" << correspondence->GetNumberOfPoints() << std::endl;
// 				correspondingSets.push_back(correspondence);
// 			}
// 		
// 		PointsPointerType mergedPointSet = mergePoints(correspondingSets[0], correspondingSets[1], tolRadius);
// 		for(int i = 2; i < correspondingSets.size(); ++i)
// 		{
// 			mergedPointSet = mergePoints(mergedPointSet, correspondingSets[i], tolRadius);
// 		}
// 		
// 		std::string summaryName(outName);
// 		summaryName += "_summary.csv";
// 		std::ofstream summaryOutStream(summaryName.c_str());
// 		summaryOutStream << "Correspondence radius (microns)" << "\t" << tolRadius << std::endl;
// 		summaryOutStream << "Correspondence set" << "\t" << mergedPointSet->GetNumberOfPoints() << std::endl;
// 		summaryOutStream.close();
// 		
// 		Reader * landmarkWriter = new Reader(outName.c_str(), outName.c_str());
// 		landmarkWriter->writeLandmarkFile(mergedPointSet);
// 		delete landmarkWriter;
// 	}
// 	else
// 		std::cout << "Usage: ./LandmarkCorrespondence [landmark set 1 filename] [landmark set 2 filename] [...] [landmark set n filename] [tolerance radius]" << std::endl;
// 	return 0;
}

PointsPointerType assignCorrespondingPoints(PointsPointerType refPoints, PointsPointerType matchPoints, double tolRadius)
{
	std::list< int > assignedRefPts;
	std::list< int > assignedMatchPts;
	bool change = 0;
	do
	{
// 			std::cout << "next iteration..." << std::endl;
		int oldNrPts = assignedMatchPts.size();
		std::map< int, int > candidateCorrespondences;
		std::list< int > candidateMatchPts;
		for(int ii = 0; ii < matchPoints->GetNumberOfPoints(); ++ii)
		{
			if(std::find(assignedMatchPts.begin(), assignedMatchPts.end(), ii) != assignedMatchPts.end())
				continue;
			
			double pt[3];
			matchPoints->GetPoint(ii, pt);
			bool hasCorrespondence = 0;
			double minDist = 1E10;
			int minID = -1;
			for(int jj = 0; jj < refPoints->GetNumberOfPoints(); ++jj)
			{
				if(std::find(assignedRefPts.begin(), assignedRefPts.end(), jj) != assignedRefPts.end())
					continue;
				
				double refPt[3];
				refPoints->GetPoint(jj, refPt);
// 				double dist = vtkMath::Distance2BetweenPoints(pt, refPt);
				double dist = (pt[0] - refPt[0])*(pt[0] - refPt[0]) + (pt[1] - refPt[1])*(pt[1] - refPt[1]);
				if(dist <= tolRadius*tolRadius && dist < minDist)
				{
					hasCorrespondence = 1;
					minDist = dist;
					minID = jj;
				}
			}
			if(hasCorrespondence)
			{
				candidateCorrespondences.insert(std::pair< int, int >(ii, minID));
				candidateMatchPts.push_back(minID);
			}
		}
		for(int ii = 0; ii < matchPoints->GetNumberOfPoints(); ++ii)
		{
			if(candidateCorrespondences.find(ii) != candidateCorrespondences.end())
			{
				int candID = candidateCorrespondences[ii];
				int multiplicity = 0;
				std::list< int >::const_iterator candidateMatchPtsIt;
				for(candidateMatchPtsIt = candidateMatchPts.begin(); candidateMatchPtsIt != candidateMatchPts.end(); ++candidateMatchPtsIt)
					if(*candidateMatchPtsIt == candID)
						++multiplicity;
				if(multiplicity > 1)
				{
// 						std::cout << "Pt " << candID << " has multiplicity " << multiplicity << std::endl;
					//pick closest of all ambiguous pts
					double pt[3];
					refPoints->GetPoint(candID, pt);
					double minDist = 1E10;
					int minID = -1;
					for(int jj = 0; jj < matchPoints->GetNumberOfPoints(); ++jj)
					{
						if(std::find(assignedMatchPts.begin(), assignedMatchPts.end(), jj) != assignedMatchPts.end())
							continue;
						
						double matchPt[3];
						matchPoints->GetPoint(jj, matchPt);
// 						double dist = vtkMath::Distance2BetweenPoints(pt, matchPt);
						double dist = (pt[0] - matchPt[0])*(pt[0] - matchPt[0]) + (pt[1] - matchPt[1])*(pt[1] - matchPt[1]);
						if(dist <= tolRadius*tolRadius && dist < minDist)
						{
							minDist = dist;
							minID = jj;
						}
					}
// 						std::cout << "\tminDist = " << sqrt(minDist) << std::endl;
// 						std::cout << "\tminID = " << minID << std::endl;
					assignedMatchPts.push_back(minID);
					assignedRefPts.push_back(candID);
				}
				else if(multiplicity == 1)
				{
					assignedMatchPts.push_back(ii);
					assignedRefPts.push_back(candID);
				}
				else
				{
					std::cout << "Warning: candidate corresponding point without multiplicity!" << std::endl;
				}
			}
		}
		
		candidateCorrespondences.clear();
		candidateMatchPts.clear();
		int newNrPts = assignedMatchPts.size();
		if(newNrPts != oldNrPts) change = 1;
		else change = 0;
	} while(change);
	
	std::cout << "Found " << assignedMatchPts.size() << " matching points!" << std::endl;
	
	PointsPointerType correspondingPts = PointsPointerType::New();
	correspondingPts->SetDataTypeToFloat();
	int fnPtCount = 0;
	for(int ii = 0; ii < matchPoints->GetNumberOfPoints(); ++ii)
	{
		if(std::find(assignedMatchPts.begin(), assignedMatchPts.end(), ii) != assignedMatchPts.end())
		{
			double * refPt = new double[3];
			matchPoints->GetPoint(ii, refPt);
			correspondingPts->InsertPoint(fnPtCount, refPt);
			++fnPtCount;
		}
	}
	
	return correspondingPts;
}

PointsPointerType mergePoints(PointsPointerType refPoints, PointsPointerType matchPoints, double tolRadius)
{
	std::list< int > assignedRefPts;
	std::list< int > assignedMatchPts;
	bool change = 0;
	do
	{
// 			std::cout << "next iteration..." << std::endl;
		int oldNrPts = assignedMatchPts.size();
		std::map< int, int > candidateCorrespondences;
		std::list< int > candidateMatchPts;
		for(int ii = 0; ii < matchPoints->GetNumberOfPoints(); ++ii)
		{
			if(std::find(assignedMatchPts.begin(), assignedMatchPts.end(), ii) != assignedMatchPts.end())
				continue;
			
			double pt[3];
			matchPoints->GetPoint(ii, pt);
			bool hasCorrespondence = 0;
			double minDist = 1E10;
			int minID = -1;
			for(int jj = 0; jj < refPoints->GetNumberOfPoints(); ++jj)
			{
				if(std::find(assignedRefPts.begin(), assignedRefPts.end(), jj) != assignedRefPts.end())
					continue;
				
				double refPt[3];
				refPoints->GetPoint(jj, refPt);
// 				double dist = vtkMath::Distance2BetweenPoints(pt, refPt);
				double dist = (pt[0] - refPt[0])*(pt[0] - refPt[0]) + (pt[1] - refPt[1])*(pt[1] - refPt[1]);
				if(dist <= tolRadius*tolRadius && dist < minDist)
				{
					hasCorrespondence = 1;
					minDist = dist;
					minID = jj;
				}
			}
			if(hasCorrespondence)
			{
				candidateCorrespondences.insert(std::pair< int, int >(ii, minID));
				candidateMatchPts.push_back(minID);
			}
		}
		for(int ii = 0; ii < matchPoints->GetNumberOfPoints(); ++ii)
		{
			if(candidateCorrespondences.find(ii) != candidateCorrespondences.end())
			{
				int candID = candidateCorrespondences[ii];
				int multiplicity = 0;
				std::list< int >::const_iterator candidateMatchPtsIt;
				for(candidateMatchPtsIt = candidateMatchPts.begin(); candidateMatchPtsIt != candidateMatchPts.end(); ++candidateMatchPtsIt)
					if(*candidateMatchPtsIt == candID)
						++multiplicity;
				if(multiplicity > 1)
				{
// 						std::cout << "Pt " << candID << " has multiplicity " << multiplicity << std::endl;
					//pick closest of all ambiguous pts
					double pt[3];
					refPoints->GetPoint(candID, pt);
					double minDist = 1E10;
					int minID = -1;
					for(int jj = 0; jj < matchPoints->GetNumberOfPoints(); ++jj)
					{
						if(std::find(assignedMatchPts.begin(), assignedMatchPts.end(), jj) != assignedMatchPts.end())
							continue;
						
						double matchPt[3];
						matchPoints->GetPoint(jj, matchPt);
// 						double dist = vtkMath::Distance2BetweenPoints(pt, matchPt);
						double dist = (pt[0] - matchPt[0])*(pt[0] - matchPt[0]) + (pt[1] - matchPt[1])*(pt[1] - matchPt[1]);
						if(dist <= tolRadius*tolRadius && dist < minDist)
						{
							minDist = dist;
							minID = jj;
						}
					}
// 						std::cout << "\tminDist = " << sqrt(minDist) << std::endl;
// 						std::cout << "\tminID = " << minID << std::endl;
					assignedMatchPts.push_back(minID);
					assignedRefPts.push_back(candID);
				}
				else if(multiplicity == 1)
				{
					assignedMatchPts.push_back(ii);
					assignedRefPts.push_back(candID);
				}
				else
				{
					std::cout << "Warning: candidate corresponding point without multiplicity!" << std::endl;
				}
			}
		}
		
		candidateCorrespondences.clear();
		candidateMatchPts.clear();
		int newNrPts = assignedMatchPts.size();
		if(newNrPts != oldNrPts) change = 1;
		else change = 0;
	} while(change);
	
	std::cout << "Found " << assignedMatchPts.size() << " matching points!" << std::endl;
	
	PointsPointerType mergePts = PointsPointerType::New();
	mergePts->SetDataTypeToFloat();
	int fnPtCount = 0;
	for(int ii = 0; ii < refPoints->GetNumberOfPoints(); ++ii)
	{
		double * refPt = new double[3];
		refPoints->GetPoint(ii, refPt);
		mergePts->InsertPoint(fnPtCount, refPt);
		++fnPtCount;
	}
	for(int ii = 0; ii < matchPoints->GetNumberOfPoints(); ++ii)
	{
		// only add points that do NOT have correspondence (i.e., we create set 1 OR set 2)
		if(std::find(assignedMatchPts.begin(), assignedMatchPts.end(), ii) == assignedMatchPts.end())
		{
			double * matchPt = new double[3];
			matchPoints->GetPoint(ii, matchPt);
			mergePts->InsertPoint(fnPtCount, matchPt);
			++fnPtCount;
		}
	}
	
	return mergePts;
};




