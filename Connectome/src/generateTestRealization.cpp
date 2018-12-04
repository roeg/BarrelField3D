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
#include "../../common/barrel_field.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// #define DEBUG

// Vincent Definitions 11/2015

struct SynapseInfo {
	unsigned int preID;
	unsigned int postID;
// 	unsigned int preCellType;
// 	unsigned int postCellType;
	std::vector< float > synapseLocation;
	float preSomaDistance;
	float postSomaDistance;
	std::string postsynapticStructure;
};

struct CellInfo {
	unsigned int cellID;
	unsigned int cellType;
	std::vector< float > somaLocation;
// 	bool insideS1;
// 	unsigned int column;
// 	bool insideColumn;
	unsigned int region;
	unsigned int nearestColumn;
	std::string laminarLocation;
// 	unsigned int principalWhisker;
};

std::map< unsigned int, unsigned int > cortexRegionLUT;
std::map< unsigned int, unsigned int > VPMRegionLUT;

std::list< unsigned int > getCellTypeIndices(unsigned int cellTypeID, std::vector< unsigned int > cellTypeIDs);
void getPCenterOfStructure(AmiraSpatialGraph * sg, int ID, double centerPt[3]);
TransformPointerType amiraToVTKTransform(double * amiraTransform);
std::pair< double *, unsigned int > placeSynapse(AmiraSpatialGraph * morphology, unsigned int preType, unsigned int postType, gsl_rng * rng);

std::map< unsigned int, unsigned int > readPrePostIDMap(const char * filename);
void writeSynapseInfoFile(const char * outputFilename, std::list< SynapseInfo * > synInfoList);
void writeCellInfoFile(const char * outputFilename, std::list< CellInfo * > cellInfoList);
void initializeLUTs();

int main( int argc , char * argv[])
{
	if(argc == 5)
	{
		initializeLUTs();
		const char * spatialGraphSetFilename = argv[1];
		const char * prePostIDMapFilename = argv[2];
		const char * innervationMatrixFilename = argv[3];
		const char * outputFilename = argv[4];
		
		const gsl_rng_type * T;
		gsl_rng * randomNumberGenerator;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		randomNumberGenerator = gsl_rng_alloc(T);
		
		BarrelField * SBF = new BarrelField;
		
		std::cout << "Reading pre-/postsynaptic ID pairs" << std::endl;
		std::map< unsigned int, unsigned int > prePostIDMap = readPrePostIDMap(prePostIDMapFilename);
		
		std::vector< unsigned int > originalGraphIndices;
		std::vector< unsigned int > cellTypeIDs;
		std::vector< double * > spatialGraphTransforms;
		std::vector< std::string > originalGraphFiles;
		std::map< unsigned int, std::string > cellTypeIDLabels;
		std::vector< AmiraSpatialGraph * > readSpatialGraphs;
		Reader::readSpatialGraphSetFile(spatialGraphSetFilename, originalGraphIndices, cellTypeIDs, spatialGraphTransforms, originalGraphFiles, cellTypeIDLabels);
		
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
		
		ConnectionMatrix * connectome = new ConnectionMatrix;
		Reader * matrixReader = new Reader(innervationMatrixFilename);
		matrixReader->readConnectionMatrix(connectome);
		
		// go through all presynaptic neurons in connection matrix
		// check if corresponding postsynaptic ID is also in connection matrix
		// (unless preType = VPM)
		// store all thus determined valid IDs ("complete" cells)
		// load all valid postsynapic dendrite morphologies
		// collect cell information
		std::vector< unsigned int > validPreIDs, validPostIDs;
		std::map< unsigned int, unsigned int >::const_iterator prePostIDMapIt;
		for(prePostIDMapIt = prePostIDMap.begin(); prePostIDMapIt != prePostIDMap.end(); ++prePostIDMapIt)
		{
			unsigned int preID = prePostIDMapIt->first;
			unsigned int postID = prePostIDMapIt->second;
			if((std::find(connectome->preColumnIDs[D2].begin(), connectome->preColumnIDs[D2].end(), preID) != connectome->preColumnIDs[D2].end())
				&& (std::find(connectome->postColumnIDs[D2].begin(), connectome->postColumnIDs[D2].end(), postID) != connectome->postColumnIDs[D2].end()))
			{
				validPreIDs.push_back(preID);
				validPostIDs.push_back(postID);
			}
		}
		unsigned int completeCells = validPreIDs.size();
		std::cout << "Found " << completeCells << " valid (complete) cells in connectome!" << std::endl;
		std::cout << "Adding VPM neurons to allowed presynaptic cells" << std::endl;
		std::vector< unsigned int > validPreIDsVPM = connectome->getPreColumnCelltypeSelection(D2, VPM);
		
		std::map< unsigned int, AmiraSpatialGraph * > postMorphologies;
		for(int i = 0; i < validPostIDs.size(); ++i)
		{
			unsigned int postID = validPostIDs[i];
			AmiraSpatialGraph * newMorphology = new AmiraSpatialGraph;
			newMorphology->mergeSpatialGraph(readSpatialGraphs[originalGraphIndices[postID]]);
			newMorphology->setTransformation(amiraToVTKTransform(spatialGraphTransforms[postID]));
			newMorphology->applyTransformation();
			postMorphologies.insert(std::pair< unsigned int, AmiraSpatialGraph * >(postID, newMorphology));
		}
		
		std::list< CellInfo * > cellInfoList;
		for(int i = 0; i < validPostIDs.size(); ++i)
		{
			CellInfo * newCell = new CellInfo;
			unsigned int cellID = validPostIDs[i];
			unsigned int cellType = cellTypeIDs[cellID];
			AmiraSpatialGraph * cellMorphology = postMorphologies[cellID];
			double somaPt[3];
			getPCenterOfStructure(cellMorphology, Soma, somaPt);
			std::vector< float > somaLocation;
			somaLocation.push_back(somaPt[0]);
			somaLocation.push_back(somaPt[1]);
			somaLocation.push_back(somaPt[2]);
			
			bool insideS1 = SBF->isInsideS1(somaPt);
			unsigned int column = -1;
			bool insideColumn = 0;
			std::string laminarLocation("N/A");
// 			if(insideS1)
// 			{
				column = SBF->closestBarrel(somaPt);
				insideColumn = SBF->insideColumn(somaPt);
				unsigned int layer = SBF->laminarPosition(somaPt);
				if(layer == SUPRA)
				{
					laminarLocation = std::string("S");
				}
				if(layer == GRAN)
				{
					laminarLocation = std::string("G");
				}
				if(layer == INFRA)
				{
					laminarLocation = std::string("I");
				}
				newCell->cellID = cellID;
				newCell->cellType = cellType;
				newCell->somaLocation = somaLocation;
// 				newCell->insideS1 = insideS1;
// 				newCell->column = column;
// 				newCell->insideColumn = insideColumn;
				if(insideColumn)
				{
					newCell->region = cortexRegionLUT[column];
				}
				else
				{
					newCell->region = cortexRegionLUT[Septum];
				}
				newCell->nearestColumn = cortexRegionLUT[column];
				newCell->laminarLocation = laminarLocation;
				cellInfoList.push_back(newCell);
// 			}
// 			else
// 			{
// 				column = SBF->closestBarrel(somaPt);
// 				unsigned int layer = SBF->laminarPosition(somaPt);
// 				std::cout << "Oops! Found a non-VPM cell outside of S1..." << std::endl;
// 				std::cout << "CellID:" << cellID << std::endl;
// 				std::cout << "Cell type:" << cellType << std::endl;
// 				std::cout << "Soma loc:" << somaLocation[0] << "," << somaLocation[1] << "," << somaLocation[2] << std::endl;
// 				std::cout << "Nearest column:" << SBF->int2Labels[column] << std::endl;
// 				std::cout << "Layer:" << layer << std::endl;
// 			}
		}
		for(int i = 0; i < validPreIDsVPM.size(); ++i)
		{
			CellInfo * newCell = new CellInfo;
			unsigned int cellID = validPreIDsVPM[i];
			unsigned int cellType = cellTypeIDs[cellID];
			AmiraSpatialGraph * cellMorphology = new AmiraSpatialGraph;
			cellMorphology->mergeSpatialGraph(readSpatialGraphs[originalGraphIndices[cellID]]);
			cellMorphology->setTransformation(amiraToVTKTransform(spatialGraphTransforms[cellID]));
			cellMorphology->applyTransformation();
			double somaPt[3];
			getPCenterOfStructure(cellMorphology, Soma, somaPt);
			std::vector< float > somaLocation;
			somaLocation.push_back(somaPt[0]);
			somaLocation.push_back(somaPt[1]);
			somaLocation.push_back(somaPt[2]);
			
			bool insideS1 = 0;
			unsigned int column = -1;
			bool insideColumn = 0;
			std::string laminarLocation("N/A");
			newCell->cellID = cellID;
			newCell->cellType = cellType;
			newCell->somaLocation = somaLocation;
// 			newCell->insideS1 = insideS1;
// 			newCell->column = column;
// 			newCell->insideColumn = insideColumn;
			newCell->region = VPMRegionLUT[D2];
			newCell->nearestColumn = -1;
			newCell->laminarLocation = laminarLocation;
			cellInfoList.push_back(newCell);
			
			delete cellMorphology;
		}
		
		// iterate through all complete cells presynaptic:
		// iterate through all postsynaptic IDs:
		// draw number of synapses from Poisson(I_ij)
		// randomly choose 3D location of synapses on postsynaptic neuron
		// for now: use euclidian distance synapse-pre/post soma
		// insert synapses into synapse data structure
		for(int i = 0; i < validPreIDsVPM.size(); ++i)
		{
			validPreIDs.push_back(validPreIDsVPM[i]);
		}
// 		std::list< SynapseInfo * > synInfoList;
		std::string synInfoName(outputFilename);
		synInfoName += "_Synapses.csv";
		std::cout << "Writing synapse info file " << synInfoName.c_str() << " ..." << std::endl;
		
		std::ofstream SynInfoFile(synInfoName.c_str());
// 		SynInfoFile << "Pre Cell ID\tPost Cell ID\tPre Cell Type\tPost Cell type\tSynapse location x\tSynapse location y\tSynapse location z\tPath length distance presynaptic soma\tPath length distance postsynaptic soma\tSubstructure" << std::endl;
		SynInfoFile << "PreNeuronID\tPostNeuronID\tPosX\tPosY\tPosZ\tDistanceToSomaPre\tDistanceToSomaPost\tSubcellularLocation" << std::endl;
	
		unsigned int totalNrSyn = 0;
		for(int i = 0; i < validPreIDs.size(); ++i)
		{
			bool isVPM = false;
			unsigned int preID = validPreIDs[i];
			if(cellTypeIDs[preID] == VPM)
			{
				isVPM = true;
			}
			unsigned int correspondingPostID = -1;
			if(!isVPM)
			{
				correspondingPostID = prePostIDMap[preID];
			}
			std::cout << "Mapping synapses from cell " << preID << " (postID = " << correspondingPostID << ")" << std::endl;
			if(isVPM)
			{
				std::cout << "\tPresynaptic cell is VPM (i.e. axon only)" << std::endl;
			}
			unsigned int preCellType;
			if(!isVPM)
			{
				preCellType = cellTypeIDs[correspondingPostID];
			}
			else
			{
				preCellType = cellTypeIDs[preID];
			}
			AmiraSpatialGraph * preMorphology;
			if(!isVPM)
			{
				preMorphology = postMorphologies[correspondingPostID];
			}
			unsigned int cellNrSyn = 0;
			unsigned int connectedCells = 0;
			for(int j = 0; j < validPostIDs.size(); ++j)
			{
				unsigned int postID = validPostIDs[j];
				if(postID == correspondingPostID)
				{
					continue;
				}
				MatrixIndexType connectionIndex(preID, postID);
				if(connectome->matrix.find(connectionIndex) != connectome->matrix.end())
				{
					float innervation = connectome->matrix[connectionIndex];
					unsigned int nrSyn = gsl_ran_poisson(randomNumberGenerator, innervation);
					if(nrSyn)
					{
#ifdef DEBUG
						std::cout << "\tPlacing " << nrSyn << " synapses on postsynaptic cell " << postID << std::endl;
#endif
						//unsigned int preID = correspondingPostID;
						//unsigned int postID;
						unsigned int postCellType = cellTypeIDs[postID];
						AmiraSpatialGraph * postMorphology = postMorphologies[postID];
						for(unsigned int k = 0; k < nrSyn; ++k)
						{
							std::pair< double *, unsigned int > mappedSynapse = placeSynapse(postMorphology, preCellType, postCellType, randomNumberGenerator);
							
							double * synapseLocation = mappedSynapse.first;
							std::vector< float > synLocVector;
							synLocVector.push_back(synapseLocation[0]);
							synLocVector.push_back(synapseLocation[1]);
							synLocVector.push_back(synapseLocation[2]);
							
							unsigned int postsynapticStructureLabel = mappedSynapse.second;
							std::string postsynapticStructure("N/A");
							if(postsynapticStructureLabel == Soma)
							{
								postsynapticStructure = std::string("Soma");
							}
							if(postsynapticStructureLabel == ApicalDendrite)
							{
								postsynapticStructure = std::string("ApicalDendrite");
							}
							if(postsynapticStructureLabel == Dendrite)
							{
								postsynapticStructure = std::string("Dendrite");
							}
							if(postsynapticStructureLabel == BasalDendrite)
							{
								postsynapticStructure = std::string("Dendrite");
							}
							
							double postSomaPt[3];
							getPCenterOfStructure(postMorphology, Soma, postSomaPt);
							float postSomaDistance = sqrt(vtkMath::Distance2BetweenPoints(synapseLocation, postSomaPt));
							float preSomaDistance = -1;
							if(!isVPM)
							{
								double preSomaPt[3];
								getPCenterOfStructure(preMorphology, Soma, preSomaPt);
								preSomaDistance = sqrt(vtkMath::Distance2BetweenPoints(synapseLocation, preSomaPt));
							}
							
							SynapseInfo * newSyn = new SynapseInfo;
							if(!isVPM)
							{
								newSyn->preID = correspondingPostID;
							}
							else
							{
								newSyn->preID = preID;
							}
							newSyn->postID = postID;
// 							newSyn->preCellType = preCellType;
// 							newSyn->postCellType = postCellType;
							newSyn->synapseLocation = synLocVector;
							newSyn->preSomaDistance = preSomaDistance;
							newSyn->postSomaDistance = postSomaDistance;
							newSyn->postsynapticStructure = postsynapticStructure;
							
							SynInfoFile << newSyn->preID << "\t";
							SynInfoFile << newSyn->postID << "\t";
// 							SynInfoFile << newSyn->preCellType << "\t";
// 							SynInfoFile << newSyn->postCellType << "\t";
							SynInfoFile << newSyn->synapseLocation[0] << "\t";
							SynInfoFile << newSyn->synapseLocation[1] << "\t";
							SynInfoFile << newSyn->synapseLocation[2] << "\t";
							SynInfoFile << newSyn->preSomaDistance << "\t";
							SynInfoFile << newSyn->postSomaDistance << "\t";
							SynInfoFile << newSyn->postsynapticStructure.c_str() << "\n";
// 							SynInfoFile.flush();
							delete newSyn;
// 							synInfoList.push_back(newSyn);
						}
						totalNrSyn += nrSyn;
						cellNrSyn += nrSyn;
						++connectedCells;
					}
				}
			}
			std::cout << "\tnumber of synapses: " << cellNrSyn << std::endl;
			std::cout << "\tconnected postsynaptic cells: " << connectedCells << std::endl;
		}
		std::cout << "Total number of synapses: " << totalNrSyn << std::endl;
		SynInfoFile.close();
		 
		// write soma location file
		writeCellInfoFile(outputFilename, cellInfoList);
		// write synapse file
// 		writeSynapseInfoFile(outputFilename, synInfoList);
		
		
// 		std::list< SynapseInfo * >::const_iterator synInfoListIt;
// 		for(synInfoListIt = synInfoList.begin(); synInfoListIt != synInfoList.end(); ++synInfoListIt)
// 		{
// 			delete *synInfoListIt;
// 		}
		std::list< CellInfo * >::const_iterator cellInfoListIt;
		for(cellInfoListIt = cellInfoList.begin(); cellInfoListIt != cellInfoList.end(); ++cellInfoListIt)
		{
			delete *cellInfoListIt;
		}
		
		gsl_rng_free(randomNumberGenerator);
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

std::pair< double*, unsigned int > placeSynapse(AmiraSpatialGraph* morphology, unsigned int preType, unsigned int postType, gsl_rng* rng)
{
	std::list< unsigned int > excTypes;
	excTypes.push_back(L2);
	excTypes.push_back(L34);
	excTypes.push_back(L4py);
	excTypes.push_back(L4sp);
	excTypes.push_back(L4ss);
	excTypes.push_back(L5st);
	excTypes.push_back(L5tt);
	excTypes.push_back(L6cc);
	excTypes.push_back(L6ccinv);
	excTypes.push_back(L6ct);
	excTypes.push_back(VPM);
	std::list< unsigned int > inhTypes;
	inhTypes.push_back(SymLocal1);
	inhTypes.push_back(SymLocal2);
	inhTypes.push_back(SymLocal3);
	inhTypes.push_back(SymLocal4);
	inhTypes.push_back(SymLocal5);
	inhTypes.push_back(SymLocal6);
	inhTypes.push_back(L1);
	inhTypes.push_back(L23Trans);
	inhTypes.push_back(L45Sym);
	inhTypes.push_back(L45Peak);
	inhTypes.push_back(L56Trans);
	
#define EXC 1
#define INH 0
	unsigned int preFunctionalType = -1, postFunctionalType = -1;
	if(std::find(excTypes.begin(), excTypes.end(), preType) != excTypes.end())
	{
		preFunctionalType = EXC;
	}
	else if(std::find(inhTypes.begin(), inhTypes.end(), preType) != inhTypes.end())
	{
		preFunctionalType = INH;
	}
	if(std::find(excTypes.begin(), excTypes.end(), postType) != excTypes.end())
	{
		postFunctionalType = EXC;
	}
	else if(std::find(inhTypes.begin(), inhTypes.end(), postType) != inhTypes.end())
	{
		postFunctionalType = INH;
	}
	
	std::vector< unsigned int > targetEdgeIndices;
	if(postFunctionalType == EXC)
	{
		if(preFunctionalType == EXC)
		{
			// place syns on apical/basal
			for(int i = 0; i < morphology->edgesPointer()->size(); ++i)
			{
				if(morphology->edgesPointer()->at(i)->label == ApicalDendrite || morphology->edgesPointer()->at(i)->label == Dendrite
					|| morphology->edgesPointer()->at(i)->label == BasalDendrite)
				{
					targetEdgeIndices.push_back(i);
				}
			}
			
		}
		if(preFunctionalType == INH)
		{
			// place syns on soma/apical/basal
			for(int i = 0; i < morphology->edgesPointer()->size(); ++i)
			{
				if(morphology->edgesPointer()->at(i)->label == ApicalDendrite || morphology->edgesPointer()->at(i)->label == Dendrite
					|| morphology->edgesPointer()->at(i)->label == BasalDendrite || morphology->edgesPointer()->at(i)->label == Soma)
				{
					targetEdgeIndices.push_back(i);
				}
			}
		}
	}
	if(postFunctionalType == INH)
	{
		// always target all structures
		for(int i = 0; i < morphology->edgesPointer()->size(); ++i)
		{
			if(morphology->edgesPointer()->at(i)->label == ApicalDendrite || morphology->edgesPointer()->at(i)->label == Dendrite
				|| morphology->edgesPointer()->at(i)->label == BasalDendrite || morphology->edgesPointer()->at(i)->label == Soma)
			{
				targetEdgeIndices.push_back(i);
			}
		}
	}
	
#ifdef DEBUG
	std::cout << "preType = " << preType << std::endl;
	std::cout << "preFunctionalType = " << preFunctionalType << std::endl;
	std::cout << "postType = " << postType << std::endl;
	std::cout << "postFunctionalType = " << postFunctionalType << std::endl;
	std::cout << "\tFound " << targetEdgeIndices.size() << " potential target edges" << std::endl;
#endif
	
	unsigned int * targetEdgeIndicesTmp = new unsigned int[targetEdgeIndices.size()];
	unsigned int pickedEgdeIndex;
	for(int i = 0; i < targetEdgeIndices.size(); ++i)
	{
		targetEdgeIndicesTmp[i] = i;
	}
	gsl_ran_choose(rng, &pickedEgdeIndex, 1, targetEdgeIndicesTmp, targetEdgeIndices.size(), sizeof(unsigned int));
	unsigned int nrOfEdgePoints = morphology->edgesPointer()->at(pickedEgdeIndex)->numEdgePoints;
	unsigned int * targetPointIndicesTmp = new unsigned int[nrOfEdgePoints];
	unsigned int pickedPointIndex;
	for(int i = 0; i < nrOfEdgePoints; ++i)
	{
		targetPointIndicesTmp[i] = i;
	}
	gsl_ran_choose(rng, &pickedPointIndex, 1, targetPointIndicesTmp, nrOfEdgePoints, sizeof(unsigned int));
	
	unsigned int structureLabel = morphology->edgesPointer()->at(pickedEgdeIndex)->label;
	std::list< double * >::const_iterator edgePtIt = morphology->edgesPointer()->at(pickedEgdeIndex)->edgePointCoordinates.begin();
	for(int i = 0; i < pickedPointIndex; ++i)
	{
		++edgePtIt;
	}
	double * synLocation = *edgePtIt;
	
#undef EXC
#undef INH
	
	return std::pair< double * , unsigned int>(synLocation, structureLabel);
};

std::map< unsigned int, unsigned int > readPrePostIDMap(const char* filename)
{
	std::map< unsigned int, unsigned int > prePostIDMap;
	std::ifstream PrePostMapStream(filename);
	if(!PrePostMapStream.fail())
	{
		std::string currentLine;
		while(!std::getline(PrePostMapStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				if(currentLine.find("#") != std::string::npos)
				{
					continue;
				}
				unsigned int preID, postID;
				std::sscanf(currentLine.c_str(), "%d %d", &preID, &postID);
				prePostIDMap.insert(std::pair< unsigned int, unsigned int >(preID, postID));
#ifdef DEBUG
				std::cout << "preID: " << preID << " - postID: " << postID << std::endl;
#endif
			}
		}
	}
	else
	{
		std::cout << "Error reading pre/post file " << filename << std::endl;
	}
	
	return prePostIDMap;
};

void writeSynapseInfoFile(const char* outputFilename, std::list< SynapseInfo* > synInfoList)
{
	std::string synInfoName(outputFilename);
	synInfoName += "_Synapses.csv";
	std::cout << "Writing synapse info file " << synInfoName.c_str() << " ..." << std::endl;
	
	std::ofstream SynInfoFile(synInfoName.c_str());
// 	SynInfoFile << "Pre Cell ID\tPost Cell ID\tPre Cell Type\tPost Cell type\tSynapse location x\tSynapse location y\tSynapse location z\tPath length distance presynaptic soma\tPath length distance postsynaptic soma\tSubstructure" << std::endl;
	SynInfoFile << "PreNeuronID\tPostNeuronID\tPosX\tPosY\tPosZ\tDistanceToSomaPre\tDistanceToSomaPost\tSubcellularLocation" << std::endl;
	
	std::list< SynapseInfo * >::const_iterator synInfoListIt;
	for(synInfoListIt = synInfoList.begin(); synInfoListIt != synInfoList.end(); ++synInfoListIt)
	{
		SynapseInfo * synInformation = *synInfoListIt;
		SynInfoFile << synInformation->preID << "\t";
		SynInfoFile << synInformation->postID << "\t";
// 		SynInfoFile << synInformation->preCellType << "\t";
// 		SynInfoFile << synInformation->postCellType << "\t";
		SynInfoFile << synInformation->synapseLocation[0] << "\t";
		SynInfoFile << synInformation->synapseLocation[1] << "\t";
		SynInfoFile << synInformation->synapseLocation[2] << "\t";
		SynInfoFile << synInformation->preSomaDistance << "\t";
		SynInfoFile << synInformation->postSomaDistance << "\t";
		SynInfoFile << synInformation->postsynapticStructure.c_str() << "\n";
	}
	
	SynInfoFile.close();
};

void writeCellInfoFile(const char* outputFilename, std::list< CellInfo* > cellInfoList)
{
	std::string cellInfoName(outputFilename);
	cellInfoName += "_SomaDistribution.csv";
	std::cout << "Writing cell info file " << cellInfoName.c_str() << " ..." << std::endl;
	
	std::ofstream CellInfoFile(cellInfoName.c_str());
// 	CellInfoFile << "Cell ID\tCell type\tSoma location x\tSoma location y\tSoma location z\tInsideS1\tColumn\tInside Column\tLaminarLocation" << std::endl;
	CellInfoFile << "ID\tCellTypeID\tSomaX\tSomaY\tSomaZ\tRegionID\tNearestColumnID\tLaminarLocation" << std::endl;
	
	std::list< CellInfo * >::const_iterator cellInfoListIt;
	for(cellInfoListIt = cellInfoList.begin(); cellInfoListIt != cellInfoList.end(); ++cellInfoListIt)
	{
		CellInfo * cellInformation = *cellInfoListIt;
		CellInfoFile << cellInformation->cellID << "\t";
		CellInfoFile << cellInformation->cellType << "\t";
		CellInfoFile << cellInformation->somaLocation[0] << "\t";
		CellInfoFile << cellInformation->somaLocation[1] << "\t";
		CellInfoFile << cellInformation->somaLocation[2] << "\t";
// 		CellInfoFile << cellInformation->insideS1 << "\t";
// 		CellInfoFile << cellInformation->column << "\t";
// 		CellInfoFile << cellInformation->insideColumn << "\t";
		CellInfoFile << cellInformation->region << "\t";
		if(cellInformation->nearestColumn == -1)
		{
			CellInfoFile << "N/A" << "\t";
		}
		else
		{
			CellInfoFile << cellInformation->nearestColumn << "\t";
		}
		CellInfoFile << cellInformation->laminarLocation.c_str() << "\n";
	}
	
	CellInfoFile.close();
};

void initializeLUTs()
{
	cortexRegionLUT[Septum] = 5;
	cortexRegionLUT[Alpha] = 6;
	cortexRegionLUT[A1] = 7;
	cortexRegionLUT[A2] = 8;
	cortexRegionLUT[A3] = 9;
	cortexRegionLUT[A4] = 10;
	cortexRegionLUT[Beta] = 11;
	cortexRegionLUT[B1] = 12;
	cortexRegionLUT[B2] = 13;
	cortexRegionLUT[B3] = 14;
	cortexRegionLUT[B4] = 15;
	cortexRegionLUT[Gamma] = 16;
	cortexRegionLUT[C1] = 17;
	cortexRegionLUT[C2] = 18;
	cortexRegionLUT[C3] = 19;
	cortexRegionLUT[C4] = 20;
	cortexRegionLUT[Delta] = 21;
	cortexRegionLUT[D1] = 22;
	cortexRegionLUT[D2] = 23;
	cortexRegionLUT[D3] = 24;
	cortexRegionLUT[D4] = 25;
	cortexRegionLUT[E1] = 26;
	cortexRegionLUT[E2] = 27;
	cortexRegionLUT[E3] = 28;
	cortexRegionLUT[E4] = 29;
	
	VPMRegionLUT[Alpha] = 30;
	VPMRegionLUT[A1] = 31;
	VPMRegionLUT[A2] = 32;
	VPMRegionLUT[A3] = 33;
	VPMRegionLUT[A4] = 34;
	VPMRegionLUT[Beta] = 35;
	VPMRegionLUT[B1] = 36;
	VPMRegionLUT[B2] = 37;
	VPMRegionLUT[B3] = 38;
	VPMRegionLUT[B4] = 39;
	VPMRegionLUT[Gamma] = 40;
	VPMRegionLUT[C1] = 41;
	VPMRegionLUT[C2] = 42;
	VPMRegionLUT[C3] = 43;
	VPMRegionLUT[C4] = 44;
	VPMRegionLUT[Delta] = 50;
	VPMRegionLUT[D1] = 46;
	VPMRegionLUT[D2] = 47;
	VPMRegionLUT[D3] = 48;
	VPMRegionLUT[D4] = 49;
	VPMRegionLUT[E1] = 50;
	VPMRegionLUT[E2] = 51;
	VPMRegionLUT[E3] = 52;
	VPMRegionLUT[E4] = 53;
};



