/****************************************************************************/
/*                                                                          */
/* Program:   SingleNeuronSynapseMapper                                     */
/*                                                                          */
/* File:      singleneuronmain.cpp                                          */
/*                                                                          */
/* Purpose:   single neuron synapse mapping                                 */
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

std::list< unsigned int > getCellTypeIndices(unsigned int cellTypeID, std::vector< unsigned int > cellTypeIDs);
void writeTransformedSpatialGraphSet(const char * outputFilename, std::vector< unsigned int > originalGraphIndices, std::vector< unsigned int > cellTypeIDs,
									 std::vector< double * > spatialGraphTransforms, std::vector< std::string > originalGraphFiles, std::map< unsigned int,
									 std::string > cellTypeIDLabels, double boundingBox[6]);
void writePrePostIDMap(const char * outputFilename, std::map< unsigned int, unsigned int > prePostIDMap);
void writeMorphologyShiftData(const char * outputFilename, std::map< unsigned int, std::vector< double > > originalGraphShiftValues, std::vector< std::string > originalGraphFiles);
TransformPointerType amiraToVTKTransform(double * amiraTransform);
void vtkToAmiraTransform(TransformPointerType vtkTransform, double * amiraTransform);
TransformPointerType dendriteToAxonTransform(TransformPointerType dendriteTransform, AmiraSpatialGraph* dendriteSG, AmiraSpatialGraph* axonSG, double* distBefore, double* distAfter);
void getPCenterOfStructure(AmiraSpatialGraph * sg, int ID, double centerPt[3]);
void getCOMOfStructure(AmiraSpatialGraph * sg, int ID, double centerPt[3]);

int main( int argc , char * argv[])
{
	if(argc == 3)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		
		const gsl_rng_type * T;
		gsl_rng * r;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		
		std::map< unsigned int, unsigned int > prePostIDMap;
		std::vector< unsigned int > originalGraphIndices;
		std::vector< unsigned int > cellTypeIDs;
		std::vector< double * > spatialGraphTransforms;
		std::vector< std::string > originalGraphFiles;
		std::map< unsigned int, std::string > cellTypeIDLabels;
		std::vector< AmiraSpatialGraph * > readSpatialGraphs;
		Reader::readSpatialGraphSetFile(inputFilename, originalGraphIndices, cellTypeIDs, spatialGraphTransforms, originalGraphFiles, cellTypeIDLabels);
		
		for(int i = 0; i < originalGraphFiles.size(); ++i)
		{
			std::flush(std::cout << "Loading SpatialGraph " << i << " of " << originalGraphFiles.size() << "\r");
			std::string originalName = originalGraphFiles[i];
			std::string loadName = originalName.substr(1, originalName.size()-2);
			Reader * spatialGraphReader = new Reader(loadName.c_str(), loadName.c_str());
			spatialGraphReader->readSpatialGraphFile(0);
			AmiraSpatialGraph * loadSG = spatialGraphReader->getSpatialGraph();
			readSpatialGraphs.push_back(loadSG);
			delete spatialGraphReader;
		}
		std::cout << std::endl;
		
		std::list< std::string > dendriteTypes;
		dendriteTypes.push_back(std::string("L2"));
		dendriteTypes.push_back(std::string("L34"));
		dendriteTypes.push_back(std::string("L4py"));
		dendriteTypes.push_back(std::string("L4sp"));
		dendriteTypes.push_back(std::string("L4ss"));
		dendriteTypes.push_back(std::string("L5st"));
		dendriteTypes.push_back(std::string("L5tt"));
		dendriteTypes.push_back(std::string("L6cc"));
		dendriteTypes.push_back(std::string("L6ccinv"));
		dendriteTypes.push_back(std::string("L6ct"));
		dendriteTypes.push_back(std::string("SymLocal1"));
		dendriteTypes.push_back(std::string("SymLocal2"));
		dendriteTypes.push_back(std::string("SymLocal3"));
		dendriteTypes.push_back(std::string("SymLocal4"));
		dendriteTypes.push_back(std::string("SymLocal5"));
		dendriteTypes.push_back(std::string("SymLocal6"));
		dendriteTypes.push_back(std::string("L1"));
		dendriteTypes.push_back(std::string("L23Trans"));
		dendriteTypes.push_back(std::string("L45Sym"));
		dendriteTypes.push_back(std::string("L45Peak"));
		dendriteTypes.push_back(std::string("L56Trans"));
		std::list< std::string > moveAxonTypes;
		moveAxonTypes.push_back(std::string("L2axon"));
		moveAxonTypes.push_back(std::string("L34axon"));
		moveAxonTypes.push_back(std::string("L4pyaxon"));
		moveAxonTypes.push_back(std::string("L4spaxon"));
		moveAxonTypes.push_back(std::string("L4ssaxon"));
		moveAxonTypes.push_back(std::string("L5staxon"));
		moveAxonTypes.push_back(std::string("L5ttaxon"));
		moveAxonTypes.push_back(std::string("L6ccaxon"));
		moveAxonTypes.push_back(std::string("L6ccinvaxon"));
		moveAxonTypes.push_back(std::string("L6ctaxon"));
		moveAxonTypes.push_back(std::string("SymLocal1axon"));
		moveAxonTypes.push_back(std::string("SymLocal2axon"));
		moveAxonTypes.push_back(std::string("SymLocal3axon"));
		moveAxonTypes.push_back(std::string("SymLocal4axon"));
		moveAxonTypes.push_back(std::string("SymLocal5axon"));
		moveAxonTypes.push_back(std::string("SymLocal6axon"));
		moveAxonTypes.push_back(std::string("L1axon"));
		moveAxonTypes.push_back(std::string("L23Transaxon"));
		moveAxonTypes.push_back(std::string("L45Symaxon"));
		moveAxonTypes.push_back(std::string("L45Peakaxon"));
		moveAxonTypes.push_back(std::string("L56Transaxon"));
		
// 		algorithm: determine all dendrites within D2 column,
// 		randomly pick an axon of the same type registered to the
// 		same column (->check filename!), and copy the dendrite
// 		transformation to the axon.
		std::cout << "Cell type (D2 registered)\tcount" << std::endl;
		std::map< unsigned int, std::string >::const_iterator cellTypeIDLabelsIt;
		std::map< unsigned int, std::vector< double > > originalGraphShiftValues;
		
		std::list< std::string >::const_iterator dendriteTypesIt;
		for(dendriteTypesIt = dendriteTypes.begin(); dendriteTypesIt != dendriteTypes.end(); ++dendriteTypesIt)
		{
			unsigned int dendriteType = -1;
			for(cellTypeIDLabelsIt = cellTypeIDLabels.begin(); cellTypeIDLabelsIt != cellTypeIDLabels.end(); ++cellTypeIDLabelsIt)
			{
				if(cellTypeIDLabelsIt->second == *dendriteTypesIt)
				{
					dendriteType = cellTypeIDLabelsIt->first;
					continue;
				}
			}
			std::list< unsigned int > dendriteTypeIDs = getCellTypeIndices(dendriteType, cellTypeIDs);
			unsigned int dendriteTypeCount = 0;
			std::list< unsigned int >::const_iterator dendriteTypeIDsIt;
			for(dendriteTypeIDsIt = dendriteTypeIDs.begin(); dendriteTypeIDsIt != dendriteTypeIDs.end(); ++dendriteTypeIDsIt)
			{
				if(originalGraphFiles[originalGraphIndices[*dendriteTypeIDsIt]].find("registered_D2") != std::string::npos)
				{
					++dendriteTypeCount;
				}
			}
			std::cout << *dendriteTypesIt << "\t" << dendriteTypeCount << std::endl;
		}
		
		std::list< std::string >::const_iterator axonTypesIt;
		for(axonTypesIt = moveAxonTypes.begin(); axonTypesIt != moveAxonTypes.end(); ++axonTypesIt)
		{
			unsigned int axonType = -1;
			for(cellTypeIDLabelsIt = cellTypeIDLabels.begin(); cellTypeIDLabelsIt != cellTypeIDLabels.end(); ++cellTypeIDLabelsIt)
			{
				if(cellTypeIDLabelsIt->second == *axonTypesIt)
				{
					axonType = cellTypeIDLabelsIt->first;
					continue;
				}
			}
			std::list< unsigned int > axonTypeIDs = getCellTypeIndices(axonType, cellTypeIDs);
			unsigned int axonTypeCount = 0;
			std::list< unsigned int >::const_iterator axonTypeIDsIt;
			for(axonTypeIDsIt = axonTypeIDs.begin(); axonTypeIDsIt != axonTypeIDs.end(); ++axonTypeIDsIt)
			{
				if(originalGraphFiles[originalGraphIndices[*axonTypeIDsIt]].find("registered_D2") != std::string::npos)
				{
					++axonTypeCount;
				}
			}
			std::cout << *axonTypesIt << "\t" << axonTypeCount << std::endl;
		}
		
		// HACK for Itamar's barrel: L4ss axon morphology ID239 is projecting 
		// only to the septum when shifted around, and we therefore replace it
		// by a randomly chosen different L4ss axon morphology
		// here: store L4ss dend/axon correspondence
		std::vector< unsigned int > dendriteIDsPermutationL4ss;
		
		for(dendriteTypesIt = dendriteTypes.begin(), axonTypesIt = moveAxonTypes.begin();
			dendriteTypesIt != dendriteTypes.end() && axonTypesIt != moveAxonTypes.end(); ++dendriteTypesIt, ++axonTypesIt)
		{
			unsigned int dendriteType = -1;
			for(cellTypeIDLabelsIt = cellTypeIDLabels.begin(); cellTypeIDLabelsIt != cellTypeIDLabels.end(); ++cellTypeIDLabelsIt)
			{
				if(cellTypeIDLabelsIt->second == *dendriteTypesIt)
				{
					dendriteType = cellTypeIDLabelsIt->first;
					continue;
				}
			}
			unsigned int axonType = -1;
			for(cellTypeIDLabelsIt = cellTypeIDLabels.begin(); cellTypeIDLabelsIt != cellTypeIDLabels.end(); ++cellTypeIDLabelsIt)
			{
				if(cellTypeIDLabelsIt->second == *axonTypesIt)
				{
					axonType = cellTypeIDLabelsIt->first;
					continue;
				}
			}
			
			std::list< unsigned int > dendriteTypeD2IDs;
			for(unsigned int i = 0; i < cellTypeIDs.size(); ++i)
			{
				if(cellTypeIDs[i] == dendriteType
					&& originalGraphFiles[originalGraphIndices[i]].find("registered_D2") != std::string::npos)
				{
					dendriteTypeD2IDs.push_back(i);
				}
			}
			std::list< unsigned int > axonTypeD2IDs;
			for(unsigned int i = 0; i < cellTypeIDs.size(); ++i)
			{
				if(cellTypeIDs[i] == axonType
					&& originalGraphFiles[originalGraphIndices[i]].find("registered_D2") != std::string::npos)
				{
					axonTypeD2IDs.push_back(i);
				}
			}
			unsigned int * dendriteIDsPermutation = new unsigned int[dendriteTypeD2IDs.size()];
			std::list< unsigned int >::const_iterator dendriteTypeIDsIt = dendriteTypeD2IDs.begin();
			for(int i = 0; i < dendriteTypeD2IDs.size(); ++i, ++dendriteTypeIDsIt)
			{
				dendriteIDsPermutation[i] = *dendriteTypeIDsIt;
			}
			std::cout << "Randomly assigning dendrite (" << *dendriteTypesIt << ") transformations to axon (" << *axonTypesIt << ") morphologies..." << std::endl;
			gsl_ran_shuffle(r, dendriteIDsPermutation, dendriteTypeD2IDs.size(), sizeof(unsigned int));
			
			std::list< unsigned int >::const_iterator axonTypeIDsIt = axonTypeD2IDs.begin();
			for(int i = 0; i < axonTypeD2IDs.size(); ++i, ++axonTypeIDsIt)
			{
// 				std::flush(std::cout << "Processing axon/dendrite pair " << i+1 << " of " << axonTypeD2IDs.size() << "\r");
				
				prePostIDMap.insert(std::pair< unsigned int, unsigned int >(*axonTypeIDsIt, dendriteIDsPermutation[i]));
				double * dendTransform = spatialGraphTransforms[dendriteIDsPermutation[i]];
				// HACK for Itamar's barrel
				if(*axonTypesIt == std::string("L4ssaxon"))
				{
					dendriteIDsPermutationL4ss.push_back(dendriteIDsPermutation[i]);
				}
				
				// only load SpatialGraphs here
				unsigned int dendIndex = dendriteIDsPermutation[i];
				unsigned int axonIndex = *axonTypeIDsIt;
				AmiraSpatialGraph * dendriteSG = new AmiraSpatialGraph;
				AmiraSpatialGraph * axonSG = new AmiraSpatialGraph;
				dendriteSG->mergeSpatialGraph(readSpatialGraphs[originalGraphIndices[dendIndex]]);
				axonSG->mergeSpatialGraph(readSpatialGraphs[originalGraphIndices[axonIndex]]);
				
				double distBefore, distAfter;
				TransformPointerType axonTransform = dendriteToAxonTransform(amiraToVTKTransform(dendTransform), dendriteSG, axonSG, &distBefore, &distAfter);
				vtkToAmiraTransform(axonTransform, spatialGraphTransforms[*axonTypeIDsIt]);
				if(distAfter > 0.1)
				{
					std::cout << "Problem: distAfter = " << distAfter << std::endl;
					std::cout << "dendIndex = " << dendIndex << std::endl;
					std::cout << "axonIndex = " << axonIndex << std::endl;
					std::cout << "originalGraphIndices[dendIndex] = " << originalGraphIndices[dendIndex] << std::endl;
					std::cout << "originalGraphIndices[axonIndex] = " << originalGraphIndices[axonIndex] << std::endl;
					std::cout << "cellTypeIDs[dendIndex] = " << cellTypeIDLabels[cellTypeIDs[dendIndex]] << std::endl;
					std::cout << "cellTypeIDs[axonIndex] = " << cellTypeIDLabels[cellTypeIDs[axonIndex]] << std::endl;
				}
				if(originalGraphShiftValues.find(originalGraphIndices[axonIndex]) != originalGraphShiftValues.end())
				{
					originalGraphShiftValues[originalGraphIndices[axonIndex]].push_back(distBefore);
				}
				else
				{
					std::vector< double > distBeforeVec;
					distBeforeVec.push_back(distBefore);
					originalGraphShiftValues[originalGraphIndices[axonIndex]] = distBeforeVec;
				}
				
				delete dendriteSG, delete axonSG;
			}
			std::cout << std::endl;
			delete [] dendriteIDsPermutation;
		}
		
// 		VPM: Move randomly within +-50 microns in each dimension
		unsigned int vpmID;
		for(cellTypeIDLabelsIt = cellTypeIDLabels.begin(); cellTypeIDLabelsIt != cellTypeIDLabels.end(); ++cellTypeIDLabelsIt)
		{
			if(cellTypeIDLabelsIt->second == std::string("VPM"))
			{
				vpmID = cellTypeIDLabelsIt->first;
				break;
			}
		}
		for(int i = 0, vpmCount = 0; i < cellTypeIDs.size(); ++i)
		{
			if(cellTypeIDs[i] == vpmID && originalGraphFiles[originalGraphIndices[i]].find("registered_D2") != std::string::npos)
			{
				++vpmCount;
				float dx, dy, dz;
				dx = 100*gsl_rng_uniform(r) - 50;
				dy = 100*gsl_rng_uniform(r) - 50;
				dz = 100*gsl_rng_uniform(r) - 50;
				double * axonTrans = spatialGraphTransforms[i];
				axonTrans[12] = dx;
				axonTrans[13] = dy;
				axonTrans[14] = dz;
				std::flush(std::cout << "Processing VPM axon " << vpmCount << ": sqrt(dx^2) = " << sqrt(dx*dx + dy*dy + dz*dz) << "\r");
			}
		}
		std::cout << std::endl;
		
		// HACK for Itamar's barrel: L4ss axon morphology ID239 is projecting 
		// only to the septum when shifted around, and we therefore replace it
		// by a randomly chosen different L4ss axon morphology
		
		unsigned int replaceCount = 0;
		for(dendriteTypesIt = dendriteTypes.begin(), axonTypesIt = moveAxonTypes.begin();
			dendriteTypesIt != dendriteTypes.end() && axonTypesIt != moveAxonTypes.end(); ++dendriteTypesIt, ++axonTypesIt)
		{
			unsigned int dendriteType = -1;
			for(cellTypeIDLabelsIt = cellTypeIDLabels.begin(); cellTypeIDLabelsIt != cellTypeIDLabels.end(); ++cellTypeIDLabelsIt)
			{
				if(cellTypeIDLabelsIt->second == *dendriteTypesIt)
				{
					dendriteType = cellTypeIDLabelsIt->first;
					continue;
				}
			}
			unsigned int axonType = -1;
			for(cellTypeIDLabelsIt = cellTypeIDLabels.begin(); cellTypeIDLabelsIt != cellTypeIDLabels.end(); ++cellTypeIDLabelsIt)
			{
				if(cellTypeIDLabelsIt->second == *axonTypesIt)
				{
					axonType = cellTypeIDLabelsIt->first;
					continue;
				}
			}
			
			std::list< unsigned int > dendriteTypeD2IDs;
			for(unsigned int i = 0; i < cellTypeIDs.size(); ++i)
			{
				if(cellTypeIDs[i] == dendriteType
					&& originalGraphFiles[originalGraphIndices[i]].find("registered_D2") != std::string::npos)
				{
					dendriteTypeD2IDs.push_back(i);
				}
			}
			std::list< unsigned int > axonTypeD2IDs;
			for(unsigned int i = 0; i < cellTypeIDs.size(); ++i)
			{
				if(cellTypeIDs[i] == axonType
					&& originalGraphFiles[originalGraphIndices[i]].find("registered_D2") != std::string::npos)
				{
					axonTypeD2IDs.push_back(i);
				}
			}
			
			std::cout << "Re-assigning axon morphologies of L4ss ID 239..." << std::endl;
			std::cout << "Current cell type: " << *axonTypesIt << std::endl;
			unsigned int goodL4ssIDs[] = {7403, 7404, 7405, 7406, 7408, 7409, 7410, 7411};
			std::list< unsigned int >::const_iterator axonTypeIDsIt = axonTypeD2IDs.begin();
			unsigned int L4ssCount = 0;
			for(int i = 0; i < axonTypeD2IDs.size(); ++i, ++axonTypeIDsIt)
			{
// 				std::flush(std::cout << "Processing axon/dendrite pair " << i+1 << " of " << axonTypeD2IDs.size() << "\r");
				if(*axonTypesIt == std::string("L4ssaxon"))
				{
					unsigned int dendritePermutationIndex = dendriteIDsPermutationL4ss[L4ssCount];
// 					prePostIDMap.insert(std::pair< unsigned int, unsigned int >(*axonTypeIDsIt, dendriteIDsPermutation[i]));
					double * dendTransform = spatialGraphTransforms[dendritePermutationIndex];
					
					// only load SpatialGraphs here
					unsigned int dendIndex = dendritePermutationIndex;
					unsigned int axonIndex = *axonTypeIDsIt;
					if(originalGraphIndices[axonIndex] == 7407)
					{
						std::cout << "Replacing axon morphology corresponding to dendrite ID " << dendIndex << std::endl;
						++replaceCount;
						unsigned int newAxonID = goodL4ssIDs[gsl_rng_uniform_int (r, 8)];
						originalGraphIndices[axonIndex] = newAxonID;
						AmiraSpatialGraph * dendriteSG = new AmiraSpatialGraph;
						AmiraSpatialGraph * axonSG = new AmiraSpatialGraph;
						dendriteSG->mergeSpatialGraph(readSpatialGraphs[originalGraphIndices[dendIndex]]);
						axonSG->mergeSpatialGraph(readSpatialGraphs[newAxonID]);
						
						double distBefore, distAfter;
						TransformPointerType axonTransform = dendriteToAxonTransform(amiraToVTKTransform(dendTransform), dendriteSG, axonSG, &distBefore, &distAfter);
						vtkToAmiraTransform(axonTransform, spatialGraphTransforms[*axonTypeIDsIt]);
						if(distAfter > 0.1)
						{
							std::cout << "Problem: distAfter = " << distAfter << std::endl;
							std::cout << "dendIndex = " << dendIndex << std::endl;
							std::cout << "axonIndex = " << axonIndex << std::endl;
							std::cout << "originalGraphIndices[dendIndex] = " << originalGraphIndices[dendIndex] << std::endl;
							std::cout << "originalGraphIndices[axonIndex] = " << originalGraphIndices[axonIndex] << std::endl;
							std::cout << "cellTypeIDs[dendIndex] = " << cellTypeIDLabels[cellTypeIDs[dendIndex]] << std::endl;
							std::cout << "cellTypeIDs[axonIndex] = " << cellTypeIDLabels[cellTypeIDs[axonIndex]] << std::endl;
						}
						if(originalGraphShiftValues.find(originalGraphIndices[axonIndex]) != originalGraphShiftValues.end())
						{
							originalGraphShiftValues[originalGraphIndices[axonIndex]].push_back(distBefore);
						}
						else
						{
							std::vector< double > distBeforeVec;
							distBeforeVec.push_back(distBefore);
							originalGraphShiftValues[originalGraphIndices[axonIndex]] = distBeforeVec;
						}
						delete dendriteSG, delete axonSG;
					}
					
					++L4ssCount;
				}
			}
			std::cout << std::endl;
		}
		std::cout << "Re-assigned " << replaceCount << " axon morphologies of L4ss ID 239 (should be 258)" << std::endl;
		
		// END HACK for Itamar's barrel
		
// 		// check new bounding box - this can be done more efficiently in Amira!
		double globalBounds[6];
		globalBounds[0] = -5000;
		globalBounds[1] = 5000;
		globalBounds[2] = -5000;
		globalBounds[3] = 5000;
		globalBounds[4] = -5000;
		globalBounds[5] = 5000;
		writeTransformedSpatialGraphSet(outputFilename, originalGraphIndices, cellTypeIDs, spatialGraphTransforms, originalGraphFiles, cellTypeIDLabels, globalBounds);
		
		std::string prePostName(outputFilename);
		prePostName += "_pre_post_IDs.csv";
		writePrePostIDMap(prePostName.c_str(), prePostIDMap);
		
		std::string morphologyShiftName(outputFilename);
		morphologyShiftName += "_morphology_shift_values.csv";
		writeMorphologyShiftData(morphologyShiftName.c_str(), originalGraphShiftValues, originalGraphFiles);
	}
	
	return 0;
}

std::list< unsigned int > getCellTypeIndices(unsigned int cellTypeID, std::vector< unsigned int > cellTypeIDs)
{
	std::list< unsigned int > cellTypeIndexList;
	for(unsigned int i = 0; i < cellTypeIDs.size(); ++i)
	{
		if(cellTypeIDs[i] == cellTypeID)
		{
			cellTypeIndexList.push_back(i);
		}
	}
	return cellTypeIndexList;
};

TransformPointerType amiraToVTKTransform(double* amiraTransform)
{
	HomogeneousMatrixPointerType mat = HomogeneousMatrixPointerType::New();
	for(int i = 0; i < 4; ++i)
		for(int j = 0; j < 4; ++j)
		{
			mat->SetElement(j, i, amiraTransform[i*4+j]);
		}
	
	TransformPointerType vtkTransform = TransformPointerType::New();
	vtkTransform->SetMatrix(mat);
	return vtkTransform;
};

void vtkToAmiraTransform(TransformPointerType vtkTransform, double* amiraTransform)
{
	HomogeneousMatrixPointerType mat = vtkTransform->GetMatrix();
	for(int i = 0; i < 4; ++i)
		for(int j = 0; j < 4; ++j)
		{
			amiraTransform[i*4+j] = mat->GetElement(j, i);
		}
};

TransformPointerType dendriteToAxonTransform(TransformPointerType dendriteTransform, AmiraSpatialGraph* dendriteSG, AmiraSpatialGraph* axonSG, double * distBefore, double * distAfter)
{
	double dendriteSoma[3], axonSoma[3], somaDiff[3];
// 	getPCenterOfStructure(dendriteSG, Soma, dendriteSoma);
// 	getPCenterOfStructure(axonSG, Soma, axonSoma);
	getCOMOfStructure(dendriteSG, Soma, dendriteSoma);
	getCOMOfStructure(axonSG, Soma, axonSoma);
	vtkMath::Subtract(dendriteSoma, axonSoma, somaDiff);
// 	somaDiff[0] = -axonSoma[0];
// 	somaDiff[1] = -axonSoma[1];
// 	somaDiff[2] = -axonSoma[2];
	
	TransformPointerType axonSomaShift = TransformPointerType::New();
	axonSomaShift->Translate(somaDiff);
	// Amira: T1 * T2 if T1 is first transform, and T2 second...
// 	dendriteTransform->PostMultiply();
// 	dendriteTransform->Concatenate(axonSomaShift);
// 	dendriteTransform->Print(std::cout);
// 	axonSomaShift->Concatenate(dendriteTransform);
// 	axonSomaShift->Print(std::cout);
	
	HomogeneousMatrixPointerType completeTransMat = HomogeneousMatrixPointerType::New();
	HomogeneousMatrixPointerType step1Mat = HomogeneousMatrixPointerType::New();
	HomogeneousMatrixPointerType dendRotMat = HomogeneousMatrixPointerType::New();
	HomogeneousMatrixPointerType dendTransMat = HomogeneousMatrixPointerType::New();
	HomogeneousMatrixPointerType dendMat = dendriteTransform->GetMatrix();
	dendRotMat->Identity();
	dendTransMat->Identity();
	for(int i = 0; i < 4; ++i)
		for(int j = 0; j < 4; ++j)
		{
			if(i < 3 && j < 3)
			{
				dendRotMat->SetElement(i, j, dendMat->GetElement(i, j));
			}
			if(j == 3)
			{
				dendTransMat->SetElement(i, j, dendMat->GetElement(i, j));
			}
		}
// 	dendRotMat->Print(std::cout);
// 	dendTransMat->Print(std::cout);
	vtkMatrix4x4::Multiply4x4(dendRotMat, axonSomaShift->GetMatrix(), step1Mat);
	vtkMatrix4x4::Multiply4x4(dendTransMat, step1Mat, completeTransMat);
	
	double axonSomaHom[4], axonSomaHomTrans[4], axonSomaTrans[3];
	double dendSomaHom[4], dendSomaHomTrans[4], dendSomaTrans[3];
	for(int i = 0; i < 3; ++i)
	{
		axonSomaHom[i] = axonSoma[i];
		dendSomaHom[i] = dendriteSoma[i];
	}
	axonSomaHom[3] = 1;
	dendSomaHom[3] = 1;
// 	HomogeneousMatrixPointerType mat = axonSomaShift->GetMatrix();
// 	HomogeneousMatrixPointerType mat = dendriteTransform->GetMatrix();
// 	mat->Print(std::cout);
// 	for(int i = 0; i < 4; ++i)
// 	{
// 		axonSomaHomTrans[i] = 0;
// 		for(int j = 0; j < 4; ++j)
// 		{
// // 			axonSomaHomTrans[i] += axonSomaHom[j]*mat->GetElement(i, j);
// 			axonSomaHomTrans[i] += axonSomaHom[j]*completeTransMat->GetElement(i, j);
// 		}
// 	}
	completeTransMat->MultiplyPoint(axonSomaHom, axonSomaHomTrans);
	dendMat->MultiplyPoint(dendSomaHom, dendSomaHomTrans);
	for(int i = 0; i < 3; ++i)
	{
		axonSomaTrans[i] = axonSomaHomTrans[i];
		dendSomaTrans[i] = dendSomaHomTrans[i];
	}
// 	double distBefore = sqrt(vtkMath::Distance2BetweenPoints(dendSomaTrans, axonSoma));
// 	double distAfter = sqrt(vtkMath::Distance2BetweenPoints(dendSomaTrans, axonSomaTrans));
	*distBefore = sqrt(vtkMath::Distance2BetweenPoints(dendSomaTrans, axonSoma));
	*distAfter = sqrt(vtkMath::Distance2BetweenPoints(dendSomaTrans, axonSomaTrans));
// 	std::flush(std::cout << "Distance before = " << distBefore << std::endl);
// 	std::flush(std::cout << "Distance after = " << distAfter << std::endl);
	
	TransformPointerType completeTrans = TransformPointerType::New();
	completeTrans->SetMatrix(completeTransMat);
// 	completeTrans->Print(std::cout);
	return completeTrans;
// 	return axonSomaShift;
// 	return dendriteTransform;
};

void getPCenterOfStructure(AmiraSpatialGraph * sg, int ID, double centerPt[3])
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

void getCOMOfStructure(AmiraSpatialGraph* sg, int ID, double centerPt[3])
{
	PolyDataPointerType structure = PolyDataPointerType::New();
	if(!sg->extractLandmark(ID, structure))
	{
		std::cout << "Error! Could not find structure with ID " << ID <<" in SpatialGraph!" << std::endl;
		return;
	}
	centerPt[0] = 0;
	centerPt[1] = 0;
	centerPt[2] = 0;
	unsigned int nrOfPoints = structure->GetCell(0)->GetNumberOfPoints();
	for(int i = 0; i < nrOfPoints; ++i)
	{
		double pt[3];
		structure->GetCell(0)->GetPoints()->GetPoint(i, pt);
		for(int j = 0; j < 3; ++j)
		{
			centerPt[j] += pt[j];
		}
	}
	centerPt[0] = centerPt[0]/(double)nrOfPoints;
	centerPt[1] = centerPt[1]/(double)nrOfPoints;
	centerPt[2] = centerPt[2]/(double)nrOfPoints;
};

void writeTransformedSpatialGraphSet(const char* outputFilename, std::vector< unsigned int > originalGraphIndices, std::vector< unsigned int > cellTypeIDs,
									 std::vector< double* > spatialGraphTransforms, std::vector< std::string > originalGraphFiles,
									 std::map< unsigned int, std::string > cellTypeIDLabels, double boundingBox[6])
{
	std::string format = outputFilename;
	format += ".am";
	std::cout << "Writing SpatialGraphSet file " << format << std::endl;
	std::ofstream SGSet(format.c_str());
	
	SGSet << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	SGSet << "# Axons shifted" << std::endl;
	SGSet << std::endl;
	SGSet << std::endl;
	SGSet << "define GRAPH " << originalGraphIndices.size() << std::endl;
	SGSet << std::endl;
	SGSet << "Parameters {" << std::endl;
	SGSet << "    SpatialGraphSetLabels {" << std::endl;
	std::map< unsigned int, std::string >::const_iterator cellTypeIDLabelsIt;
	for(cellTypeIDLabelsIt = cellTypeIDLabels.begin(); cellTypeIDLabelsIt != cellTypeIDLabels.end(); ++cellTypeIDLabelsIt)
	{
		SGSet << "        " << cellTypeIDLabelsIt->second << " {" << std::endl;
		SGSet << "            Color 0.890625 0.101562 0.109375," << std::endl;
		SGSet << "            Id " << cellTypeIDLabelsIt->first << std::endl;
		SGSet << "        }" << std::endl;
	}
	SGSet << "        Id 0," << std::endl;
	SGSet << "        Color 0 0 0" << std::endl;
	SGSet << "    }" << std::endl;
	SGSet << "    Files {" << std::endl;
	for(int i = 0; i < originalGraphFiles.size(); ++i)
	{
		char * fileChar = new char[64];
		sprintf(fileChar, "        File%05d {", i);
		SGSet << fileChar << std::endl;
		SGSet << "            Path " << originalGraphFiles[i] << std::endl;
		SGSet << "        }" << std::endl;
	}
	SGSet << "    }" << std::endl;
	SGSet << "    ContentType \"HxSpatialGraphSet\"," << std::endl;
	SGSet << "    BoundingBox " << boundingBox[0] << " " << boundingBox[1] << " " << boundingBox[2] << " " << boundingBox[3] << " " << boundingBox[4] << " " << boundingBox[5] << std::endl;
	SGSet << "}" << std::endl;
	SGSet << std::endl;
	SGSet << "GRAPH { int OriginalGraphIndices } @1" << std::endl;
	SGSet << "GRAPH { int SpatialGraphSetLabels } @2" << std::endl;
	SGSet << "GRAPH { float[16] AdditionalTransforms } @3" << std::endl;
	SGSet << std::endl;
	SGSet << "# Data section follows" << std::endl;
	SGSet << "@1" << std::endl;
	for(int i = 0; i < originalGraphIndices.size(); ++i)
	{
		SGSet << originalGraphIndices[i] << " " << std::endl;
	}
	SGSet << std::endl;
	SGSet << "@2" << std::endl;
	for(int i = 0; i < cellTypeIDs.size(); ++i)
	{
		SGSet << cellTypeIDs[i] << " " << std::endl;
	}
	SGSet << std::endl;
	SGSet << "@3" << std::endl;
	for(int i = 0; i < spatialGraphTransforms.size(); ++i)
	{
		SGSet << std::scientific << std::setprecision(15);
		for(int j = 0; j < 16; ++j)
		{
			SGSet << spatialGraphTransforms[i][j] << " ";
		}
		SGSet << std::endl;
	}
	SGSet << std::endl;
	
	SGSet.close();
};

void writePrePostIDMap(const char* outputFilename, std::map< unsigned int, unsigned int > prePostIDMap)
{
	std::ofstream PrePostMap(outputFilename);
	PrePostMap << "#Presynaptic neuron ID\tPostsynaptic neuron ID\n";
	std::map< unsigned int, unsigned int >::const_iterator prePostIDMapIt;
	for(prePostIDMapIt = prePostIDMap.begin(); prePostIDMapIt != prePostIDMap.end(); ++prePostIDMapIt)
	{
		PrePostMap << prePostIDMapIt->first << "\t" << prePostIDMapIt->second << std::endl;
	}
	
	PrePostMap.close();
};

void writeMorphologyShiftData(const char* outputFilename, std::map< unsigned int, std::vector< double > > originalGraphShiftValues, std::vector< std::string > originalGraphFiles)
{
	std::ofstream MorphologyShiftFile(outputFilename);
	MorphologyShiftFile << "Morphology filename\tMean shift (microns)\tSTD shift (microns)" << std::endl;
	std::map< unsigned int, std::vector< double > >::const_iterator originalGraphShiftValuesIt;
	for(originalGraphShiftValuesIt = originalGraphShiftValues.begin(); originalGraphShiftValuesIt != originalGraphShiftValues.end(); ++originalGraphShiftValuesIt)
	{
		std::string morphologyFilename = originalGraphFiles[originalGraphShiftValuesIt->first];
		std::vector< double > morphologyShiftValues = originalGraphShiftValuesIt->second;
		
		double meanShift = 0, STDShift = 0;
		for(int i = 0; i < morphologyShiftValues.size(); ++i)
		{
			double dr = morphologyShiftValues[i];
			meanShift += dr;
		}
		meanShift = meanShift/morphologyShiftValues.size();
		for(int i = 0; i < morphologyShiftValues.size(); ++i)
		{
			double dr = morphologyShiftValues[i];
			STDShift += (dr - meanShift)*(dr - meanShift);
		}
		STDShift = sqrt(STDShift/(morphologyShiftValues.size()-1));
		
		MorphologyShiftFile << morphologyFilename << "\t" << meanShift << "\t" << STDShift << std::endl;
	}
	
	MorphologyShiftFile.close();
};


