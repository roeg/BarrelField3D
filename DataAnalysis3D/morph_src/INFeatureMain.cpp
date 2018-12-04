/****************************************************************************/
/*                                                                          */
/* Program:   MorphAnalyzer                                                 */
/*                                                                          */
/* File:      main.cpp                                                      */
/*                                                                          */
/* Purpose:   Program for analysis of registered neuron morphologies with   */
/*            respect to columns, septa and layers in standardized barrel   */
/*            cortex. E.g., computes distribution of axon in different      */
/*            columns/layers, and computes z-profiles of axons taking local */
/*            orientation into account                                      */
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
#include "morph_analyzer.h"
#include "../../common/inputcheckpoint.h"

/* ALL */
//int main( int argc , char * argv[])
//{
////	const char * inputFilename = argv[1];
////	const char * outputFilename = argv[2];
//
//	//std::ifstream inputStream("/home/dudvary/testCell.txt");
//	//std::ifstream inputStream("/home/dudvary/Interneurons/Lists/registratedCells.txt");
//	//std::string rootStr = "/home/dudvary/Interneurons/hoc/registrated/";
//	//std::string outrootStr = "/home/dudvary/Interneurons/ClassificationTab/";
//
//	std::ifstream inputStream("/home/dudvary/Interneurons/L1/L1axon_D2/L1.txt");
//	std::string rootStr = "/home/dudvary/Interneurons/L1/L1axon_D2/";
//	std::string outrootStr = "/home/dudvary/Interneurons/ClassificationTab/L1/";
//
//	std::cout << "ID";
//	for(int ii = 0; ii < 9; ++ii)
//		std::cout << ",Feature_" << ii;
//	std::cout << "," << std::endl;
//	std::vector< std::string > printStr;
//	printStr.push_back("ID,extent,soma,CMtrans,CM-soma,CMtrans-soma,ratio,Z4-Z5,Z1-Z4,Z1-Z2");
//
//	if(!inputStream.fail())
//	{
//		std::string currentName;
//		while(!inputStream.eof())
//		{
//			getline(inputStream,currentName);
//			if (currentName.empty())
//				break;
//
//			std::string tmpStr = rootStr + currentName;
//			const char * inputFilename = tmpStr.c_str();
//			const char * outputFilename = tmpStr.c_str();
//
//			const char * label = "Axon";
//			int binSize = 50;
//			int binSize2 = 25;
//			double epsilon = 1E-6;
//			int extentDel = 6; // range for deleting local peak
//			int maxlen = 39; // 39*50 = 1950um (WM)
//
//			std::string ofName(outputFilename);
//			std::string ifName(inputFilename);
//
//			Reader * hocReader = new Reader(inputFilename, outputFilename);
//			if(ifName.find(".hoc") != std::string::npos)
//				hocReader->readHocFile();
//			else if(ifName.find(".am") != std::string::npos)
//				hocReader->readSpatialGraphFile(0);
//			else
//			{
//				std::cout << "Error! Can only analyze .hoc or .am files!" << std::endl;
//				return 0;
//			}
//
//			InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
//			checkPoint->checkNeuronMorphology();
//
//			Analyzer * morphAnalyzer = new Analyzer(hocReader->getSpatialGraph(), checkPoint->getParameters());
//			std::vector< double >* zProfile;
//			std::vector< double >* zProfile_25;
//
//			zProfile = morphAnalyzer->compute1DProfileInsideS1(label, binSize);
//			zProfile_25 = morphAnalyzer->compute1DProfileInsideS1(label, binSize2);
//
//			/* For CHECKING */
////			std::cout << "25um bin z profile" << std::endl;
////			for(int ii = 0; ii < zProfile_25->size(); ++ii)
////			{
////				std::cout << zProfile_25->at(ii) << std::endl;
////			}
////			std::cout << "50um bin z profile" << std::endl;
////			for(int ii = 0; ii < zProfile->size(); ++ii)
////			{
////				std::cout << zProfile->at(ii) << std::endl;
////			}
//			/* ------------ */
//			if(zProfile->size()==1)
//			{
//				std::cout << "Warning! Width is below 50um!" << std::endl;
//			}
//			if(!zProfile->size())
//			{
//				std::cout << "Error! Z Profile is empty! (bin size = " << binSize << ")" << std::endl;
//				return 0;
//			}
//			if(!zProfile_25->size())
//			{
//				std::cout << "Error! Z Profile is empty! (bin size = " << binSize2 << ")" << std::endl;
//				return 0;
//			}
//
//			// Normalize Length to 1 and restrict length to 39 from pia to WM (0 to 1950 = 39*50)
//			double totalLength_S1 = 0;
//			for(int ii = 0; ii < zProfile->size(); ++ii)
//			{
//				if (ii>=maxlen)
//					break;
//				totalLength_S1 += zProfile->at(ii);
//			}
//			for(int ii = 0; ii < maxlen; ++ii)
//			{
//				if (ii>=zProfile->size())
//					(*zProfile).push_back(0.0);
//				else
//					zProfile->at(ii) = zProfile->at(ii)/totalLength_S1;
//			}
//			zProfile->erase(zProfile->begin()+maxlen,zProfile->end());
//
//			/* Initialize Features */
//			std::vector< double > features;
//			features.assign(9,-1);
//
//			/* Feature 2: Soma */
//			bool hasSoma = checkPoint->getParameters().somaFlag;
//			double somaPt[] = {0,0,0};
//
//			if (hasSoma)
//			{
//				AmiraSpatialGraph * neuronMorphology = morphAnalyzer->getNeuronMorphology();
//				PolyDataPointerType structure = PolyDataPointerType::New();
//				if(!neuronMorphology->extractLandmark(Soma, structure))
//				{
//					std::cout << "Error! Could not find structure with ID " << Soma <<" in SpatialGraph!" << std::endl;
//					return 0;
//				}
//
//				PointsPointerType pts = structure->GetCell(0)->GetPoints();
//				for (int ii=0; ii < pts->GetNumberOfPoints(); ii++)
//				{
//					double tmp[3];
//					pts->GetPoint(ii,tmp);
//					somaPt[0] += tmp[0];
//					somaPt[1] += tmp[1];
//					somaPt[2] += tmp[2];
//				}
//
//				somaPt[0] /= pts->GetNumberOfPoints();
//				somaPt[1] /= pts->GetNumberOfPoints();
//				somaPt[2] /= pts->GetNumberOfPoints();
//			}
//			else
//			{
//				std::cout << "Error! No soma label found!" << std::endl;
//				return 0;
//			}
//
//			somaPt[2] = 706-somaPt[2];
//			features.at(1) = somaPt[2];
//
//			/* Feature 1: 90% Columnar Extent */
//			double mintmp = 1E6;
//			int idxSoma[2] = { -1, -1 };
//			double ratio = 0.9;
//
//			// Find Index of Soma & Compute total length of morphology
//			idxSoma[0] = int((somaPt[2]-epsilon)/binSize2);
//			idxSoma[1] = idxSoma[0];
//
//			// Sanity Check: Check total length of 25um and 50um binsz profiles
//			double totalLength25 = 0.0;
//			double totalLength25_S1 = 0.0;
//			for(int ii = 0; ii < zProfile_25->size(); ++ii)
//			{
//				if (ii>=maxlen*2)
//					totalLength25 = totalLength25 + zProfile_25->at(ii);
//				else
//					totalLength25_S1 = totalLength25_S1 + zProfile_25->at(ii);
//			}
//			totalLength25 = totalLength25 + totalLength25_S1;
//			if (fabs(totalLength25_S1-totalLength_S1) > epsilon)
//			{
//				std::cout << "Warning! Length of profiles is not the same! Length (bin size=" << binSize << ") = " << totalLength_S1;
//				std::cout << " Length2 (bin size=" << binSize2 << ") = " << totalLength25_S1 << " Diff = " << fabs(totalLength25_S1-totalLength_S1) << std::endl;
//
//				if (fabs(totalLength25_S1-totalLength_S1) > 2)
//				{
//					std::cout << "Error! Length of profiles is not the same! Length (bin size=" << binSize << ") = " << totalLength_S1;
//					std::cout << " Length2 (bin size=" << binSize2 << ") = " << totalLength25_S1 << " Diff = " << fabs(totalLength25_S1-totalLength_S1) << std::endl;
//					return 0;
//				}
//			}
//
//			if (idxSoma[0]<0 || idxSoma[0] != idxSoma[1])
//			{
//				std::cout << "Error! Soma position is outside of Z Profile!" << std::endl;
//				return 0;
//			}
//			if (totalLength25_S1<=0)
//			{
//				std::cout << "Error! Total axonal length (len=" << totalLength25_S1 << ") is equal or below zero! " << std::endl;
//				return 0;
//			}
//
//			double tmpLength = zProfile_25->at(idxSoma[0]);
//
//			while (tmpLength<ratio*totalLength25)
//			{
//				// update range
//				idxSoma[0] = idxSoma[0]-1;
//				idxSoma[1] = idxSoma[1]+1;
//
//				// check whether range is out of bound
//				if (idxSoma[0]<0)
//					idxSoma[0]=0;
//				if (idxSoma[1]>=zProfile_25->size())
//					idxSoma[1]=zProfile_25->size()-1;
//
//				// check whether values above (0 to idxSoma[0]) are all zero
//				tmpLength = 0;
//				for(int ii = 0; ii < idxSoma[0]+1; ++ii)
//				{
//					tmpLength += zProfile_25->at(ii);
//				}
//				if (tmpLength==0)
//					idxSoma[0]=idxSoma[0]+1;
//
//				// check whether values below (idxSoma[1] to end) are all zero
//				tmpLength = 0;
//				for(int ii = idxSoma[1]; ii < zProfile_25->size(); ++ii)
//				{
//					tmpLength += zProfile_25->at(ii);
//				}
//				if (tmpLength==0)
//					idxSoma[1]=idxSoma[1]-1;
//
//				// compute length in current range (idxSoma[0] to idxSoma[1])
//				tmpLength = 0;
//				for(int ii = idxSoma[0]; ii < idxSoma[1]+1; ++ii)
//				{
//					tmpLength += zProfile_25->at(ii);
//				}
//				/* CHECK each loop */
//				//std::cout << idxSoma[0] << "," << idxSoma[1] << "," << tmpLength << "/" << ratio*totalLength25 << std::endl;
//				/* ----- */
//			}
//
//			// 90% columnar extent
//			features.at(0) = (idxSoma[1]-idxSoma[0])*binSize2;
//
//			if (features.at(0)<binSize2)
//			{
//				std::cout << "Error! 90% z Extent " << features.at(0) << " is smaller than " << binSize2 << std::endl;
//				return 0;
//			}
//
//			/* Feature 3: Center of non-local Mass */
//			// Find peaks in z profile
//			std::vector< int > loc;
//			std::vector< double > locVal;
//			for(int ii = 1; ii < zProfile->size()-1; ++ii)
//			{
//				if (zProfile->at(ii)>zProfile->at(ii-1) && zProfile->at(ii)>zProfile->at(ii+1))
//				{
//					loc.push_back(ii);
//					locVal.push_back(zProfile->at(ii));
//				}
//			}
//
//			// Position of maximum closest to soma (soma peak)
//			int idxMax = -1;
//			// if no local maximum found take maximum value
//			if (loc.size()==0)
//			{
//				std::cout << "Warning! No Local Maximum found in z profile, take Maximum instead" << std::endl;
//
//				// Get maximum value/position
//				double maxtmp = 0;
//				for(int ii = 0; ii < zProfile->size(); ++ii)
//				{
//					 if (zProfile->at(ii) > maxtmp)
//					{
//						maxtmp = zProfile->at(ii);
//						idxMax = ii;
//					}
//				}
//			}
//			else // Find local peak within 230um arround the soma
//			{
//				std::vector< int > ix;
//				for(int ii = 0; ii < loc.size(); ++ii)
//				{
//					double tmpDist = (loc.at(ii) * binSize + binSize/2) - somaPt[2];
//					if (fabs(tmpDist)<230)
//						ix.push_back(ii);
//				}
//
//				if (ix.size()==0)
//				{
//					std::cout << "Warning! No Local Maximum in z profile within 230um around the soma" << std::endl;
//					double mintmp = 1E6;
//					for(int ii = 0; ii < loc.size(); ++ii)
//					{
//						double tmpDist = (loc.at(ii) * binSize + binSize/2) - somaPt[2];
//						if (fabs(tmpDist)<mintmp)
//						{
//							mintmp = fabs(tmpDist);
//							ix.clear();
//							ix.push_back(ii);
//						}
//					}
//				}
//				// Get Index with Maximum Peak within 230um
//				int ix2 = -1;
//				double maxtmp = 0;
//				for(int ii = 0; ii < ix.size(); ++ii)
//				{
//					if (locVal.at(ix.at(ii)) > maxtmp)
//					{
//						maxtmp = locVal.at(ix.at(ii));
//						ix2 = ii;
//					}
//				}
//
//				if (ix2<0)
//				{
//					std::cout << "Warning! No local maximum within 230um around soma" << std::endl;
//					return 0;
//				}
//				idxMax = loc.at(ix.at(ix2));
//			}
//
//			// Get non local Length Profile
//			int cutPts[2] = { idxMax-extentDel, idxMax+extentDel};
//			if (cutPts[0]<0)
//				cutPts[0]=0;
//			if (cutPts[1]>maxlen-1)
//				cutPts[1]=maxlen-1;
//
//			std::vector< double > zProfileClean;
//			for(int ii = 0; ii < maxlen; ++ii)
//			{
//				if (ii<cutPts[0] || ii>cutPts[1])
//					zProfileClean.push_back(zProfile->at(ii));
//				else
//					zProfileClean.push_back(0);
//				/* CHECK non local profile */
//				std::cout << zProfileClean.at(ii) << "," << zProfile->at(ii) << std::endl;
//				/* ---- */
//			}
//			// Get non local length ratio
//			double nonlocalratio = 0.0;
//			for(int ii = 0; ii < zProfileClean.size(); ++ii)
//			{
//				nonlocalratio = nonlocalratio+zProfileClean.at(ii);
//			}
//			/* CHECK trans profile */
//			//std::cout << cutPts[0] << "," << cutPts[1] << "," << ix2 << "," << ix.at(ix2) << "," << loc.at(ix.at(ix2)) << std::endl;
//			/* --- */
//
//			/* Feature 4: Center of Mass - Soma */
//			double tmpCM = 0;
//			for(int ii = 0; ii < zProfile->size(); ++ii)
//			{
//				tmpCM = tmpCM + zProfile->at(ii) * (ii * binSize + binSize/2);
//			}
//			features.at(3) = somaPt[2]-tmpCM;
//
//			/* Feature 3: Center of Mass Non Local */
//			tmpCM = 0;
//			if (nonlocalratio>0)
//			{
//				for(int ii = 0; ii < zProfileClean.size(); ++ii)
//				{
//					tmpCM = tmpCM + zProfileClean.at(ii) * (ii * binSize + binSize/2);
//
//					//std::cout << tmpCM << std::endl;
//				}
//				//std::cout << ">> " << tmpCM << " " << nonlocalratio << std::endl;
//				features.at(2) = 706-tmpCM/nonlocalratio;
//			}
//			else // use Center of Mass if nonlocal profile is all zero
//			{
//				features.at(2) = -somaPt[2] + features.at(3) + 706;
//			}
//
//			/* Feature 5: Center of non-local Mass - Soma */
//			features.at(4) = (features.at(2) + features.at(1)) - 706;
//
//			/* Feature 6: Non-Local Length Ratio */
//			features.at(5) = nonlocalratio;
//
//			/* Feature 7: L45HotZone-WM */
//			double f6 = 0.0;
//			if (nonlocalratio>0)
//			{
//				for(int ii = 13; ii < zProfileClean.size(); ++ii)
//				{
//					f6 = f6+zProfileClean.at(ii);
//				}
//				f6 = f6/nonlocalratio;
//			}
//			features.at(6) = f6;
//
//			/* Feature 8: Pia-L45Hotzone */
//			double f7 = 0.0;
//			if (nonlocalratio>0)
//			{
//				for(int ii = 0; ii < 21; ++ii)
//				{
//					f7 = f7+zProfileClean.at(ii);
//				}
//				f7 = f7/nonlocalratio;
//			}
//			features.at(7) = f7;
//
//			/* Feature 9: L2HotZone */
//			double f8 = 0.0;
//			for(int ii = 0; ii < 8; ++ii)
//			{
//				f8 = f8+zProfileClean.at(ii);
//			}
//			features.at(8) = f8;
//
//			// store in string
//			size_t idx = currentName.find(".hoc");
//			std::ostringstream strtmp;
//			strtmp << currentName.substr(0,idx);
//			for(int ii = 0; ii < features.size(); ++ii)
//				strtmp << std::setprecision(15) << "," << features.at(ii);
//			printStr.push_back(strtmp.str());
//
//			// Store in .csv file
//			std::string fname = outrootStr + currentName.substr(0,idx) + ".csv";
//			std::ofstream Writer;
//			Writer.open(fname.c_str());
//			Writer << "\"extent\",\"soma\",\"CMtrans\",\"CM-soma\",\"CMtrans-soma\",\"ratio\",\"Z4-Z5\",\"Z1-Z4\",\"Z1-Z2\"" << std::endl;
//			for(int ii = 0; ii < features.size()-1; ++ii)
//			{
//				Writer << std::setprecision(15) << "\"" << features.at(ii) << "\",";
//			}
//			Writer << std::setprecision(15) << "\"" << features.at(features.size()-1) << "\"\n";
//			Writer.close();
//
//			delete morphAnalyzer, delete hocReader;
//		}
//
//		// Store all features of all cells in .csv file
//		std::string fname_all = outrootStr + "ground_truth_inh/features_all.csv";
//		std::ofstream WriterAll;
//		WriterAll.open(fname_all.c_str());
//		for(int ii = 0; ii < printStr.size(); ++ii)
//		{
//			WriterAll << printStr.at(ii) << std::endl;
//		}
//		WriterAll.close();
//	}
//	return 0;
//}

int main( int argc , char * argv[])
{
	const char * inputFilename = argv[1];
	const char * outputFilename = argv[2];

	//std::string ofName(outputFilename);
	std::string ifName(inputFilename);

	const char * label = "Axon";
	int binSize = 50;
	int binSize2 = 25;
	double epsilon = 1E-6;
	int extentDel = 6; // range for deleting local peak
	int maxlen = 39; // 39*50 = 1950um (WM)

	Reader * hocReader = new Reader(inputFilename, outputFilename);
	if(ifName.find(".hoc") != std::string::npos)
		hocReader->readHocFile();
	else if(ifName.find(".am") != std::string::npos)
		hocReader->readSpatialGraphFile(0);
	else
	{
		std::cout << "Error! Can only analyze .hoc or .am files!" << std::endl;
		return 0;
	}

	InputCheckpoint * checkPoint = new InputCheckpoint(hocReader->getSpatialGraph());
	checkPoint->checkNeuronMorphology();

	Analyzer * morphAnalyzer = new Analyzer(hocReader->getSpatialGraph(), checkPoint->getParameters());
	std::vector< double >* zProfile;
	std::vector< double >* zProfile_25;

	zProfile = morphAnalyzer->compute1DProfileInsideS1(label, binSize);
	zProfile_25 = morphAnalyzer->compute1DProfileInsideS1(label, binSize2);

	/* For CHECKING */
//			std::cout << "25um bin z profile" << std::endl;
//			for(int ii = 0; ii < zProfile_25->size(); ++ii)
//			{
//				std::cout << zProfile_25->at(ii) << std::endl;
//			}
//			std::cout << "50um bin z profile" << std::endl;
//			for(int ii = 0; ii < zProfile->size(); ++ii)
//			{
//				std::cout << zProfile->at(ii) << std::endl;
//			}
	/* ------------ */
	if(zProfile->size()==1)
	{
		std::cout << "Warning! Width is below 50um!" << std::endl;
	}
	if(!zProfile->size())
	{
		std::cout << "Error! Z Profile is empty! (bin size = " << binSize << ")" << std::endl;
		return 0;
	}
	if(!zProfile_25->size())
	{
		std::cout << "Error! Z Profile is empty! (bin size = " << binSize2 << ")" << std::endl;
		return 0;
	}

	// Normalize Length to 1 and restrict length to 39 from pia to WM (0 to 1950 = 39*50)
	double totalLength_S1 = 0;
	for(int ii = 0; ii < zProfile->size(); ++ii)
	{
		if (ii>=maxlen)
			break;
		totalLength_S1 += zProfile->at(ii);
	}
	for(int ii = 0; ii < maxlen; ++ii)
	{
		if (ii>=zProfile->size())
			(*zProfile).push_back(0.0);
		else
			zProfile->at(ii) = zProfile->at(ii)/totalLength_S1;
	}
	zProfile->erase(zProfile->begin()+maxlen,zProfile->end());

	/* Initialize Features */
	std::vector< double > features;
	features.assign(9,-1);

	/* Feature 2: Soma */
	bool hasSoma = checkPoint->getParameters().somaFlag;
	double somaPt[] = {0,0,0};

	if (hasSoma)
	{
		AmiraSpatialGraph * neuronMorphology = morphAnalyzer->getNeuronMorphology();
		PolyDataPointerType structure = PolyDataPointerType::New();
		if(!neuronMorphology->extractLandmark(Soma, structure))
		{
			std::cout << "Error! Could not find structure with ID " << Soma <<" in SpatialGraph!" << std::endl;
			return 0;
		}

		PointsPointerType pts = structure->GetCell(0)->GetPoints();
		for (int ii=0; ii < pts->GetNumberOfPoints(); ii++)
		{
			double tmp[3];
			pts->GetPoint(ii,tmp);
			somaPt[0] += tmp[0];
			somaPt[1] += tmp[1];
			somaPt[2] += tmp[2];
		}

		somaPt[0] /= pts->GetNumberOfPoints();
		somaPt[1] /= pts->GetNumberOfPoints();
		somaPt[2] /= pts->GetNumberOfPoints();
	}
	else
	{
		std::cout << "Error! No soma label found!" << std::endl;
		return 0;
	}

	somaPt[2] = 706-somaPt[2];
	features.at(1) = somaPt[2];

	/* Feature 1: 90% Columnar Extent */
	double mintmp = 1E6;
	int idxSoma[2] = { -1, -1 };
	double ratio = 0.9;

	// Find Index of Soma & Compute total length of morphology
	idxSoma[0] = int((somaPt[2]-epsilon)/binSize2);
	idxSoma[1] = idxSoma[0];

	// Sanity Check: Check total length of 25um and 50um binsz profiles
	double totalLength25 = 0.0;
	double totalLength25_S1 = 0.0;
	for(int ii = 0; ii < zProfile_25->size(); ++ii)
	{
		if (ii>=maxlen*2)
			totalLength25 = totalLength25 + zProfile_25->at(ii);
		else
			totalLength25_S1 = totalLength25_S1 + zProfile_25->at(ii);
	}
	totalLength25 = totalLength25 + totalLength25_S1;
	if (fabs(totalLength25_S1-totalLength_S1) > epsilon)
	{
		std::cout << "Warning! Length of profiles is not the same! Length (bin size=" << binSize << ") = " << totalLength_S1;
		std::cout << " Length2 (bin size=" << binSize2 << ") = " << totalLength25_S1 << " Diff = " << fabs(totalLength25_S1-totalLength_S1) << std::endl;

		if (fabs(totalLength25_S1-totalLength_S1) > 2)
		{
			std::cout << "Error! Length of profiles is not the same! Length (bin size=" << binSize << ") = " << totalLength_S1;
			std::cout << " Length2 (bin size=" << binSize2 << ") = " << totalLength25_S1 << " Diff = " << fabs(totalLength25_S1-totalLength_S1) << std::endl;
			return 0;
		}
	}

	if (idxSoma[0]<0 || idxSoma[0] != idxSoma[1])
	{
		std::cout << "Error! Soma position is outside of Z Profile!" << std::endl;
		return 0;
	}
	if (totalLength25_S1<=0)
	{
		std::cout << "Error! Total axonal length (len=" << totalLength25_S1 << ") is equal or below zero! " << std::endl;
		return 0;
	}

	double tmpLength = zProfile_25->at(idxSoma[0]);

	while (tmpLength<ratio*totalLength25)
	{
		// update range
		idxSoma[0] = idxSoma[0]-1;
		idxSoma[1] = idxSoma[1]+1;

		// check whether range is out of bound
		if (idxSoma[0]<0)
			idxSoma[0]=0;
		if (idxSoma[1]>=zProfile_25->size())
			idxSoma[1]=zProfile_25->size()-1;

		// check whether values above (0 to idxSoma[0]) are all zero
		tmpLength = 0;
		for(int ii = 0; ii < idxSoma[0]+1; ++ii)
		{
			tmpLength += zProfile_25->at(ii);
		}
		if (tmpLength==0)
			idxSoma[0]=idxSoma[0]+1;

		// check whether values below (idxSoma[1] to end) are all zero
		tmpLength = 0;
		for(int ii = idxSoma[1]; ii < zProfile_25->size(); ++ii)
		{
			tmpLength += zProfile_25->at(ii);
		}
		if (tmpLength==0)
			idxSoma[1]=idxSoma[1]-1;

		// compute length in current range (idxSoma[0] to idxSoma[1])
		tmpLength = 0;
		for(int ii = idxSoma[0]; ii < idxSoma[1]+1; ++ii)
		{
			tmpLength += zProfile_25->at(ii);
		}
		/* CHECK each loop */
		//std::cout << idxSoma[0] << "," << idxSoma[1] << "," << tmpLength << "/" << ratio*totalLength25 << std::endl;
		/* ----- */
	}

	// 90% columnar extent
	features.at(0) = (idxSoma[1]-idxSoma[0])*binSize2;

	if (features.at(0)<binSize2)
	{
		std::cout << "Error! 90% z Extent " << features.at(0) << " is smaller than " << binSize2 << std::endl;
		return 0;
	}

	/* Feature 3: Center of non-local Mass */
	// Find peaks in z profile
	std::vector< int > loc;
	std::vector< double > locVal;
	for(int ii = 1; ii < zProfile->size()-1; ++ii)
	{
		if (zProfile->at(ii)>zProfile->at(ii-1) && zProfile->at(ii)>zProfile->at(ii+1))
		{
			loc.push_back(ii);
			locVal.push_back(zProfile->at(ii));
		}
	}

	// Position of maximum closest to soma (soma peak)
	int idxMax = -1;
	// if no local maximum found take maximum value
	if (loc.size()==0)
	{
		std::cout << "Warning! No Local Maximum found in z profile, take Maximum instead" << std::endl;

		// Get maximum value/position
		double maxtmp = 0;
		for(int ii = 0; ii < zProfile->size(); ++ii)
		{
			 if (zProfile->at(ii) > maxtmp)
			{
				maxtmp = zProfile->at(ii);
				idxMax = ii;
			}
		}
	}
	else // Find local peak within 230um arround the soma
	{
		std::vector< int > ix;
		for(int ii = 0; ii < loc.size(); ++ii)
		{
			double tmpDist = (loc.at(ii) * binSize + binSize/2) - somaPt[2];
			if (fabs(tmpDist)<230)
				ix.push_back(ii);
		}

		if (ix.size()==0)
		{
			std::cout << "Warning! No Local Maximum in z profile within 230um around the soma" << std::endl;
			double mintmp = 1E6;
			for(int ii = 0; ii < loc.size(); ++ii)
			{
				double tmpDist = (loc.at(ii) * binSize + binSize/2) - somaPt[2];
				if (fabs(tmpDist)<mintmp)
				{
					mintmp = fabs(tmpDist);
					ix.clear();
					ix.push_back(ii);
				}
			}
		}
		// Get Index with Maximum Peak within 230um
		int ix2 = -1;
		double maxtmp = 0;
		for(int ii = 0; ii < ix.size(); ++ii)
		{
			if (locVal.at(ix.at(ii)) > maxtmp)
			{
				maxtmp = locVal.at(ix.at(ii));
				ix2 = ii;
			}
		}

		if (ix2<0)
		{
			std::cout << "Warning! No local maximum within 230um around soma" << std::endl;
			return 0;
		}
		idxMax = loc.at(ix.at(ix2));
	}

	// Get non local Length Profile
	int cutPts[2] = { idxMax-extentDel, idxMax+extentDel};
	if (cutPts[0]<0)
		cutPts[0]=0;
	if (cutPts[1]>maxlen-1)
		cutPts[1]=maxlen-1;

	std::vector< double > zProfileClean;
	for(int ii = 0; ii < maxlen; ++ii)
	{
		if (ii<cutPts[0] || ii>cutPts[1])
			zProfileClean.push_back(zProfile->at(ii));
		else
			zProfileClean.push_back(0);
		/* CHECK non local profile */
		//std::cout << zProfileClean.at(ii) << "," << zProfile->at(ii) << std::endl;
		/* ---- */
	}
	// Get non local length ratio
	double nonlocalratio = 0.0;
	for(int ii = 0; ii < zProfileClean.size(); ++ii)
	{
		nonlocalratio = nonlocalratio+zProfileClean.at(ii);
	}
	/* CHECK trans profile */
	//std::cout << cutPts[0] << "," << cutPts[1] << "," << ix2 << "," << ix.at(ix2) << "," << loc.at(ix.at(ix2)) << std::endl;
	/* --- */

	/* Feature 4: Center of Mass - Soma */
	double tmpCM = 0;
	for(int ii = 0; ii < zProfile->size(); ++ii)
	{
		tmpCM = tmpCM + zProfile->at(ii) * (ii * binSize + binSize/2);
	}
	features.at(3) = somaPt[2]-tmpCM;

	/* Feature 3: Center of Mass Non Local */
	tmpCM = 0;
	if (nonlocalratio>0)
	{
		for(int ii = 0; ii < zProfileClean.size(); ++ii)
		{
			tmpCM = tmpCM + zProfileClean.at(ii) * (ii * binSize + binSize/2);

			//std::cout << tmpCM << std::endl;
		}
		//std::cout << ">> " << tmpCM << " " << nonlocalratio << std::endl;
		features.at(2) = 706-tmpCM/nonlocalratio;
	}
	else // use Center of Mass if nonlocal profile is all zero
	{
		features.at(2) = -somaPt[2] + features.at(3) + 706;
	}

	/* Feature 5: Center of non-local Mass - Soma */
	features.at(4) = (features.at(2) + features.at(1)) - 706;

	/* Feature 6: Non-Local Length Ratio */
	features.at(5) = nonlocalratio;

	/* Feature 7: L45HotZone-WM */
	double f6 = 0.0;
	if (nonlocalratio>0)
	{
		for(int ii = 13; ii < zProfileClean.size(); ++ii)
		{
			f6 = f6+zProfileClean.at(ii);
		}
		f6 = f6/nonlocalratio;
	}
	features.at(6) = f6;

	/* Feature 8: Pia-L45Hotzone */
	double f7 = 0.0;
	if (nonlocalratio>0)
	{
		for(int ii = 0; ii < 21; ++ii)
		{
			f7 = f7+zProfileClean.at(ii);
		}
		f7 = f7/nonlocalratio;
	}
	features.at(7) = f7;

	/* Feature 9: L2HotZone */
	double f8 = 0.0;
	for(int ii = 0; ii < 8; ++ii)
	{
		f8 = f8+zProfileClean.at(ii);
	}
	features.at(8) = f8;

	// Store in .csv file
	std::ofstream Writer;
	Writer.open(outputFilename);
	Writer << "\"extent\",\"soma\",\"CMtrans\",\"CM-soma\",\"CMtrans-soma\",\"ratio\",\"Z4-Z5\",\"Z1-Z4\",\"Z1-Z2\"" << std::endl;
	for(int ii = 0; ii < features.size()-1; ++ii)
	{
		Writer << std::setprecision(15) << "\"" << features.at(ii) << "\",";
	}
	Writer << std::setprecision(15) << "\"" << features.at(features.size()-1) << "\"\n";
	Writer.close();

	delete morphAnalyzer, delete hocReader;
	return 0;
}

