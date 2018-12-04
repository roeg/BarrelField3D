/****************************************************************************/
/*                                                                          */
/* Program:   AverageBarrelField                                            */
/*                                                                          */
/* File:      main.cpp                                                      */
/*                                                                          */
/* Purpose:  Program for registering barrel field reconstructions to common */
/*           coordinate system using orthogonal transformations and         */
/*           creating standard barrel field from average position of BT/BB  */
/*           and measured average values for laminar thicknesses. Also      */
/*           creates Pia/WM/L4 border surfaces and 3D z-axis field.         */
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
#include "registration.h"

int main( int argc , char * argv[])
{
	if(argc == 3)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		Reader * hocReader = new Reader(inputFilename, outputFilename);
		hocReader->readHocFile();
		hocReader->writeSpatialGraphFile();
		delete hocReader;
	}
	if(argc > 3)
	{
		std::vector< const char * > inputNameVec;
		std::vector< std::string > outputNameVec;
		std::vector< const char * > outFileNames;
		std::vector< Reader * > amiraReaderVec;
		std::list< AmiraSpatialGraph * > barrelFields;
		
		const char * outputFilename = argv[argc-1];
		std::string append(outputFilename);
		
		for(int ii = 1; ii < argc-1; ++ii)
		{
			inputNameVec.push_back(argv[ii]);
			std::string newStr(argv[ii]);
			newStr += "_" + append;
			outputNameVec.push_back(newStr);
			outFileNames.push_back(newStr.c_str());
			Reader * newReader = new Reader(inputNameVec[ii-1], outFileNames[ii-1]);
			newReader->readSpatialGraphFile(0);
			barrelFields.push_back(newReader->getSpatialGraph());
			amiraReaderVec.push_back(newReader);
		}
		std::string finalOutputName(inputNameVec[0]);
		finalOutputName += "_" + append + "_avg_barrel_field";
		outFileNames.push_back(finalOutputName.c_str());
		
		Registration * spatialGraphProcessor = new Registration();
		std::list< AmiraSpatialGraph * > regBarrelFields = spatialGraphProcessor->averageBarrelField(barrelFields.size(), barrelFields, outFileNames);
		Reader * avgWriter = new Reader(finalOutputName.c_str(), finalOutputName.c_str());
		avgWriter->setSpatialGraph(spatialGraphProcessor->getSpatialGraph());
		avgWriter->writeSpatialGraphFile();
		
		std::list< AmiraSpatialGraph * >::const_iterator bfIt = regBarrelFields.begin();
		for(int ii = 0; ii < barrelFields.size(); ++ii, ++bfIt)
		{
			amiraReaderVec[ii]->setSpatialGraph(*bfIt);
			amiraReaderVec[ii]->writeSpatialGraphFile();
			delete amiraReaderVec[ii];
		}
		
		for(bfIt = regBarrelFields.begin(); bfIt != regBarrelFields.end(); ++bfIt)
			delete *bfIt;
		for(bfIt = barrelFields.begin(); bfIt != barrelFields.end(); ++bfIt)
			delete *bfIt;
		regBarrelFields.clear(), barrelFields.clear(), outFileNames.clear();
		delete spatialGraphProcessor;
	}
// 	if(argc == 5)
// 	{
// 		const char * inputFilename1 = argv[1];
// 		const char * inputFilename2 = argv[2];
// 		const char * inputFilename3 = argv[3];
// 		const char * outputFilename = argv[4];
// 		std::string ofName1(inputFilename1);
// 		std::string ofName2(inputFilename2);
// 		std::string ofName3(inputFilename3);
// 		std::string ofName4(inputFilename1);
// 		std::string append(outputFilename);
// 		ofName1 += "_" + append;
// 		ofName2 += "_" + append;
// 		ofName3 += "_" + append;
// 		ofName4 += "_" + append + "_avg_barrel_field";
// 		
// 		std::vector< const char * > fileNames;
// 		fileNames.push_back(ofName1.c_str());
// 		fileNames.push_back(ofName2.c_str());
// 		fileNames.push_back(ofName3.c_str());
// 		fileNames.push_back(ofName4.c_str());
// 		std::list< AmiraSpatialGraph * > barrelFields;
// 		
// 		Reader * amiraReader1 = new Reader(inputFilename1, ofName1.c_str());
// 		amiraReader1->readSpatialGraphFile(0);
// 		barrelFields.push_back(amiraReader1->getSpatialGraph());
// 		
// 		Reader * amiraReader2 = new Reader(inputFilename2, ofName2.c_str());
// 		amiraReader2->readSpatialGraphFile(0);
// 		barrelFields.push_back(amiraReader2->getSpatialGraph());
// 		
// 		Reader * amiraReader3 = new Reader(inputFilename3, ofName3.c_str());
// 		amiraReader3->readSpatialGraphFile(0);
// 		barrelFields.push_back(amiraReader3->getSpatialGraph());
// 		
// 		Registration * spatialGraphProcessor = new Registration();
// 		std::list< AmiraSpatialGraph * > regBarrelFields = spatialGraphProcessor->averageBarrelField(barrelFields.size(), barrelFields, fileNames);
// 		Reader * avgWriter = new Reader(ofName4.c_str(), ofName4.c_str());
// 		avgWriter->setSpatialGraph(spatialGraphProcessor->getSpatialGraph());
// 		avgWriter->writeSpatialGraphFile();
// 		
// 		std::list< AmiraSpatialGraph * >::const_iterator bfIt = regBarrelFields.begin();
// 		amiraReader1->setSpatialGraph(*bfIt);
// 		amiraReader1->writeSpatialGraphFile();
// 		++bfIt;
// 		amiraReader2->setSpatialGraph(*bfIt);
// 		amiraReader2->writeSpatialGraphFile();
// 		++bfIt;
// 		amiraReader3->setSpatialGraph(*bfIt);
// 		amiraReader3->writeSpatialGraphFile();
// 		
// 		for(bfIt = regBarrelFields.begin(); bfIt != regBarrelFields.end(); ++bfIt)
// 			delete *bfIt;
// 		for(bfIt = barrelFields.begin(); bfIt != barrelFields.end(); ++bfIt)
// 			delete *bfIt;
// 		regBarrelFields.clear(), barrelFields.clear(), fileNames.clear();
// 		delete spatialGraphProcessor;
// 		delete amiraReader1;
// 		delete amiraReader2;
// 		delete amiraReader3;
// 	}
	
	return 0;
}