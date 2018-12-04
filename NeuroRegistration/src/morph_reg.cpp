/****************************************************************************/
/*                                                                          */
/* Program:   NeuroRegistration                                             */
/*                                                                          */
/* File:      morph_reg.cpp                                                 */
/*                                                                          */
/* Purpose:   class providing the methods for registration of 3D neuron     */
/*            morphologies to the standard barrel field                     */
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
#include "morph_reg.h"

// #define PIPELINE_DOC
// #define REG_ACCURACY
// #define DEBUG

Registration::Registration(AmiraSpatialGraph * inputSpatialGraph, bool registerLandmarks, int landmarkMode) : Utilities(inputSpatialGraph)
{
	if(registerLandmarks)
	{
		detectSectionThickness();
		detectInputZDirection();
#ifdef DEBUG
		std::cout << "Input SpatialGraph z-reversed detected: " << zReversed << std::endl;
#endif
		setUpContourCorrection(landmarkMode);
	}
	else	// use this when registering already registered morphology to new home column
	{
		piaSpacing = 0;
		zReversed = 0;
	}
	workingSG = NULL;
};

Registration::Registration()
{
	detectSectionThickness();
	detectInputZDirection();
	workingSG = NULL;
};

Registration::~Registration()
{
	if(workingSG) delete workingSG;
};

/******************************************************************************/
/*start of the morphology registration pipeline                               */
/******************************************************************************/
int Registration::startNeuronRegistration(const char* outputFilename, int axisSelect, int cellType, int goodApical,
										  bool piaCorrection, bool onlyAxon, std::string manualHB)
{
	if(spatialGraph)
	{
		#ifdef PIPELINE_DOC
		pipelineDocName = outputFilename;
		#endif
		#ifdef REG_ACCURACY
		logFileName = outputFilename;
		#endif
		
		// set flags for registration mode
// 		ignoreBarrels = noBarrels;
		ignoreBarrels = 0;	// for historical reasons, this implements a dV registration
		correctPia = piaCorrection;
		if(onlyAxon)
		{
			if(!axonFlag)
			{
				std::cout << "Error! Trying to force axon-only mode but there is no axon in input hoc file!" << std::endl;
				std::cout << "Aborting registration..." << std::endl;
				return 0;
			}
			somaFlag = 0;
			dendriteFlag = 0;
			apicalFlag = 0;
			basalFlag = 0;
			axisSelect = 0;
			cellType = 0;
			goodApical = 0;
			
			std::cout << std::endl;
			if(manualHB != std::string("0"))
			{
				int manualHomeBarrel = 0;
				std::map< std::string, int >::const_iterator compareIt;
				for(compareIt = SBF->labels2Int.begin(); compareIt != SBF->labels2Int.end(); ++compareIt)
				{
					if(manualHB.compare(compareIt->first) == 0)
					{
						manualHomeBarrel = compareIt->second;
						break;
					}
				}
				if(compareIt == SBF->labels2Int.end() || SBF->avgBarrels.find(manualHomeBarrel) == SBF->avgBarrels.end())
				{
					std::cout << "Manual home barrel not recognized. Has to be Greek Arc or A-E row, Arc 1-4" << std::endl;
					return 0;
				}
				spatialGraph->setHomeBarrel(manualHomeBarrel);
				std::cout << "Registration optimized for manual home barrel: " << manualHB << std::endl;
			}
			
			std::cout << "Using axon registration mode; ignoring any other cell structures!" << std::endl;
		}
		
		std::string transformName(outputFilename);
		transformName += "_transform.log";
		TransformLog.open(transformName.c_str());
		TransformLog << "NeuroMap v3.2.1" << std::endl;
		TransformLog << "Updated for BuildingBrains-3D.org" << std::endl;
		TransformLog << "Build " << __DATE__ << " " << __TIME__ << std::endl;
		TransformLog << std::endl;
		
		double * neuronAxis = NULL;
		if(!workingSG)
		{
			workingSG = new AmiraSpatialGraph;
		}
		workingSG->mergeSpatialGraph(spatialGraph);
		workingSG->setHomeBarrel(spatialGraph->getHomeBarrel());
		if(!workingSG->getHomeBarrel())
		{
			workingSG->setHomeBarrel(closestMorphBarrel());
			// if there's still no home barrel, do registration to global barrel field optimum
			if(!workingSG->getHomeBarrel())
			{
// 				std::cout << "Error! No home barrel in hoc file and not able to determine home barrel automatically!" << std::endl;
// 				std::cout << "Aborting registration..." << std::endl;
				std::cout << "WARNING: No home barrel in hoc file and not able to determine home barrel automatically!" << std::endl;
			}
		}
		
		std::cout << std::endl;
		std::cout << "**********************************************************" << std::endl;
		std::cout << "**********************************************************" << std::endl;
		std::cout << ">>> Starting registration!" << std::endl;
		if(workingSG->getHomeBarrel())
		{
			std::cout << ">>> Manual home barrel: " << SBF->int2Labels[workingSG->getHomeBarrel()] << std::endl;
		}
		else
		{
			std::cout << ">>> WARNING: No home barrel determined, registration without home barrel optimization!" << std::endl;
		}
		std::cout << std::endl;
		TransformLog << "Starting registration!" << std::endl;
// 		TransformLog << "Registration mode: " << ignoreBarrels << "\t (0 - use barrels for orientation; 1 - daVinci mode)" << std::endl;
		if(L1flag)
		{
			TransformLog << "L1 neuron registration mode" << std::endl;
		}
		if(workingSG->getHomeBarrel())
		{
			TransformLog << "Home barrel: " << SBF->int2Labels[workingSG->getHomeBarrel()] << std::endl;
		}
		else
		{
			TransformLog << "WARNING: No home barrel determined, registration without home barrel optimization!" << std::endl;
		}
		
		PolyDataPointerType pia = NULL;
		PolyDataPointerType wm = NULL;
		std::cout << ">>> Computing landmark reconstruction..." << std::endl;
		std::cout << std::endl;
		if(piaFlag)
		{
			pia = surfaceReconstruction(Pia);
#ifdef DEBUG
			std::string piaDebugName(outputFilename);
			piaDebugName += "_pia_reconstruction";
			Reader * piaDebugWriter = new Reader(piaDebugName.c_str(), piaDebugName.c_str());
			piaDebugWriter->writeAmiraSurfaceFile(pia);
#endif
		}
		if(wmFlag)
		{
			wm = surfaceReconstruction(WhiteMatter);
#ifdef DEBUG
			std::string wmDebugName(outputFilename);
			wmDebugName += "_WM_reconstruction";
			Reader * wmDebugWriter = new Reader(wmDebugName.c_str(), wmDebugName.c_str());
			wmDebugWriter->writeAmiraSurfaceFile(wm);
#endif
		}
		constantLandmarkOffset(pia, wm, -39);
		
		#ifdef PIPELINE_DOC
		if(piaFlag)
		{
// 			constantPiaOffset(pia, -100);
			globalPia = PolyDataPointerType::New();
			globalPia->DeepCopy(pia);
			registeredSomaPosition("slicing direction", globalPia);
			double minDendDist = unregisteredPiaDendriteDistance(globalPia);
			TransformLog << "Unregistered minimum Pia-dendrite distance (micron):\t" << minDendDist << std::endl;
			TransformLog << std::endl;
		}
		#endif
		
		std::map< int, double * > barrelAxes;
		std::map< int, Column * > * barrels = new std::map< int, Column * >[2]; // input barrels = barrels[0], reference barrels = barrels[1]
		barrelReconstruction(pia, neuronAxis, &barrelAxes, barrels);
#ifdef DEBUG
		if(workingSG->getHomeBarrel())
		{
			std::cout << "Home barrel height: " << barrels[0][workingSG->getHomeBarrel()]->getHeight() << std::endl;
		}
#endif
		
		std::cout << ">>> Computing landmark registration..." << std::endl;
		TransformPointerType landmarkTransform = landmarkRegistration(barrels);
		
		int regHomeBarrel;
		if(onlyAxon)
		{
			regHomeBarrel = workingSG->getHomeBarrel();
		}
		else
		{
			regHomeBarrel = registeredHomeBarrel();
		}
		
		// in case we have a home barrel, check whether the registration went ok
		if(regHomeBarrel)
		{
			// case: registered to wrong home barrel
			if(regHomeBarrel != workingSG->getHomeBarrel() && workingSG->isLabelInSpatialGraph(regHomeBarrel))
			{
				#ifdef PIPELINE_DOC
				globalPia->DeepCopy(pia);
				#endif
				workingSG->clear();
				workingSG->mergeSpatialGraph(spatialGraph);
				workingSG->setHomeBarrel(regHomeBarrel);
				barrelAxes.clear(), barrels[0].clear(), barrels[1].clear();
				std::map< int, double * > newBarrelAxes;
				std::map< int, Column * > * newBarrels = new std::map< int, Column * >[2]; // input barrels = barrels[0], reference barrels = barrels[1]
				std::cout << ">>> Detected new home barrel: " << SBF->int2Labels[workingSG->getHomeBarrel()] << std::endl;
				std::cout << ">>> Repeating landmark registration with " << SBF->int2Labels[workingSG->getHomeBarrel()] << " as reference..." << std::endl;
				TransformLog << std::endl;
				TransformLog << "WARNING!!! Detected new home barrel: " << SBF->int2Labels[workingSG->getHomeBarrel()] << std::endl;
				TransformLog << "Repeating landmark registration with " << SBF->int2Labels[workingSG->getHomeBarrel()] << " as reference..." << std::endl;
				barrelReconstruction(pia, neuronAxis, &newBarrelAxes, newBarrels);
				landmarkTransform = landmarkRegistration(newBarrels);
				barrelAxes = newBarrelAxes, barrels = newBarrels;
			}
			// case: registered to wrong home barrel, but true home barrel is not available
			else if(regHomeBarrel != workingSG->getHomeBarrel() && !workingSG->isLabelInSpatialGraph(regHomeBarrel))
			{
				int nextHomeBarrel = closestAvgBarrel();
				if(!nextHomeBarrel)
				{
					std::cout << "Error! No home barrel in hoc file and not able to determine home barrel automatically!" << std::endl;
					std::cout << "Aborting registration..." << std::endl;
					return 0;
				}
				#ifdef PIPELINE_DOC
				globalPia->DeepCopy(pia);
				#endif
				workingSG->clear();
				workingSG->mergeSpatialGraph(spatialGraph);
				workingSG->setHomeBarrel(nextHomeBarrel);
				barrelAxes.clear(), barrels[0].clear(), barrels[1].clear();
				std::map< int, double * > newBarrelAxes;
				std::map< int, Column * > * newBarrels = new std::map< int, Column * >[2]; // input barrels = barrels[0], reference barrels = barrels[1]
				std::cout << ">>> Detected new home barrel: " << SBF->int2Labels[regHomeBarrel] << std::endl;
				std::cout << ">>> Warning! Home barrel contours not available!" << std::endl;
				std::cout << ">>> Repeating landmark registration with " << SBF->int2Labels[nextHomeBarrel] << " as reference..." << std::endl;
				TransformLog << std::endl;
				TransformLog << "WARNING!!! Detected new home barrel: " << SBF->int2Labels[regHomeBarrel] << std::endl;
				TransformLog << "WARNING! Home barrel contours not available!" << std::endl;
				TransformLog << "Repeating landmark registration with " << SBF->int2Labels[nextHomeBarrel] << " as reference..." << std::endl;
				barrelReconstruction(pia, neuronAxis, &newBarrelAxes, newBarrels);
				landmarkTransform = landmarkRegistration(newBarrels);
				barrelAxes = newBarrelAxes, barrels = newBarrels;
			}
			
			// after landmark registration aligned with z-axis
			// -> set zReversed to 0 for other method calls
			if(zReversed)
				zReversed = 0;
			alignHomeBarrelAxis(barrels, workingSG->getHomeBarrel());
#ifdef DEBUG
			if(workingSG->getHomeBarrel())
			{
				std::cout << "Home barrel height after refinement: " << barrels[0][workingSG->getHomeBarrel()]->getHeight() << std::endl;
			}
#endif
		}
			
		// after landmark registration aligned with z-axis
		// -> set zReversed to 0 for other method calls
		if(zReversed)
			zReversed = 0;
		
		if(piaFlag)
		{
			TransformFilterType piaTransform = TransformFilterType::New();
			piaTransform->SetTransform(landmarkTransform);
			piaTransform->SetInput(pia);
			piaTransform->Update();
			pia = piaTransform->GetOutput();
		}
		if(wmFlag)
		{
			TransformFilterType wmTransform = TransformFilterType::New();
			wmTransform->SetTransform(landmarkTransform);
			wmTransform->SetInput(wm);
			wmTransform->Update();
			wm = wmTransform->GetOutput();
		}
		
		TransformLog << std::endl;
		TransformLog << "Global landmark registration:" << std::endl;
// 		for(int ii = 0; ii < 4; ++ii)
// 		{
// 			TransformLog << "[";
// 			for(int jj = 0; jj < 4; ++jj)
// 			{
// 				if(jj < 3)
// 					TransformLog << landmarkTransform->GetMatrix()->GetElement(ii, jj) << " ";
// 				else
// 					TransformLog << landmarkTransform->GetMatrix()->GetElement(ii, jj);
// 			}
// 			TransformLog << "]" << std::endl;
// 		}
// 		TransformLog << std::endl;
		
		TransformLog << "# Transformation matrix Amira format" << std::endl;
		for(int ii = 0; ii < 4; ++ii)
			for(int jj = 0; jj < 4; ++jj)
				TransformLog << landmarkTransform->GetMatrix()->GetElement(jj, ii) << " ";
		TransformLog << std::endl;
		TransformLog << "# Transformation matrix python/numpy format:" << std::endl;
		TransformLog << "[";
		for(int ii = 0; ii < 4; ++ii)
		{
			TransformLog << "[";
			for(int jj = 0; jj < 4; ++jj)
			{
				TransformLog << landmarkTransform->GetMatrix()->GetElement(ii, jj);
				if(jj < 3) TransformLog << ",";
			}
			TransformLog << "]";
			if(ii < 3) TransformLog << ",";
		}
		TransformLog << "]" << std::endl;
		TransformLog << std::endl;
		TransformLog << std::endl;
		
		TransformLog << "Parameters used for neuron z-axis determination:" << std::endl;
		TransformLog << "neuronAxis = " << axisSelect << std::endl;
		TransformLog << "cellType = " << cellType << std::endl;
		TransformLog << "goodApicalDend = " << goodApical << std::endl;
		switch(axisSelect)
		{
			case 1:
				if(cellType == SUPRA || cellType == GRAN)
				{
					if(goodApical)
						neuronAxis = principalAxes(ApicalDendrite, workingSG, goodApical);
					else
						neuronAxis = detectMorphZAxis(ApicalDendrite, cellType, workingSG);
				}
				else
					neuronAxis = principalAxes(ApicalDendrite, workingSG, goodApical);
				break;
				
			case 2:
				neuronAxis = detectMorphZAxis(Axon, cellType, workingSG);
				break;
		}
		if(neuronAxis)
		{
			TransformLog << "Neuron z-axis = [" << neuronAxis[0] << "," << neuronAxis[1] << "," << neuronAxis[2] <<"]" << std::endl;
			TransformLog << std::endl;
			std::cout << std::endl;
			std::cout << ">>> Neuron axis unit vector: [" << neuronAxis[0] << "," << neuronAxis[1] << "," << neuronAxis[2] <<"]" << std::endl;
			std::cout << std::endl;
		}
		else
		{
			TransformLog << "WARNING!!! Neuron z-axis invalid or not used." << std::endl;
			TransformLog << std::endl;
			std::cout << std::endl;
			std::cout << "Neuron axis invalid!" << std::endl;
			std::cout << std::endl;
		}
		
		std::cout << ">>> Computing neuron morphology transformations..." << std::endl;
		morphologyRegistration(neuronAxis, pia, wm, barrels);
		
		if(somaFlag)
		{
			int finalHomeBarrel = registeredHomeBarrel();
			if(finalHomeBarrel)
				workingSG->setHomeBarrel(finalHomeBarrel);
		}
		
		if(piaFlag && somaFlag)
		{
			#ifndef PIPELINE_DOC
			double piaDistance = registeredPiaDistance();
			double radialDist = somaDistanceToColumnAxis();
			double minDendDist = registeredPiaDendriteDistance();
			int colSepFlag = getColumnSeptumFlag();
			int regLayer = getLaminarPosition();
			double somaPos[3], localOrientation[3] = {0,0,1};
			getPCenterOfStructure(workingSG, Soma, somaPos);
			SBF->localZAxis(somaPos, localOrientation);
			TransformLog << "Registered Pia-soma distance (micron):\t" << piaDistance << std::endl;
			TransformLog << std::endl;
			TransformLog << "Registered Soma-column axis distance (micron):\t" << radialDist << std::endl;
			TransformLog << std::endl;
			TransformLog << "Registered minimum Pia-dendrite distance (micron):\t" << minDendDist << std::endl;
			TransformLog << std::endl;
			TransformLog << "Registered column (1)/septum (0) cell:\t" << colSepFlag << std::endl;
			TransformLog << std::endl;
			TransformLog << "Registered layer ( (s)upra, (g)ran, (i)nfra ):\t";
			if(regLayer == SUPRA) TransformLog << "s" << std::endl;
			if(regLayer == GRAN) TransformLog << "g" << std::endl;
			if(regLayer == INFRA) TransformLog << "i" << std::endl;
			if(regLayer == 0) TransformLog << "N/A" << std::endl;
			TransformLog << std::endl;
			TransformLog << "Local orientation:\t[" << localOrientation[0] << "," << localOrientation[1] << "," << localOrientation[2] << "]" << std::endl;
			TransformLog << std::endl;
			
			#else
			registeredSomaPosition("z scaling");
			double minDendDist = registeredPiaDendriteDistance();
			int colSepFlag = getColumnSeptumFlag();
			TransformLog << "Registered minimum Pia-dendrite distance (micron):\t" << minDendDist << std::endl;
			TransformLog << "Registered column (1)/septum (0) cell:\t" << colSepFlag << std::endl;
			TransformLog << std::endl;
			#endif
		}
		if(!somaFlag)
		{
			double localOrientation[3] = {0,0,1};
			SBF->localZAxis(workingSG->getHomeBarrel(), localOrientation);
			TransformLog << "Local orientation:\t[" << localOrientation[0] << "," << localOrientation[1] << "," << localOrientation[2] << "]" << std::endl;
			TransformLog << std::endl;
		}
		
		#ifdef PIPELINE_DOC
		if(piaFlag)
		{
			std::string piaOutputName(outputFilename);
			piaOutputName += "_pia";
			Reader * piaWriter = new Reader(piaOutputName.c_str(), piaOutputName.c_str());
			piaWriter->writeAmiraSurfaceFile(pia);
			delete piaWriter;
		}
		if(wmFlag)
		{
			std::string wmOutputName(outputFilename);
			wmOutputName += "_WM";
			Reader * wmWriter = new Reader(wmOutputName.c_str(), wmOutputName.c_str());
			wmWriter->writeAmiraSurfaceFile(wm);
			delete wmWriter;
		}
		#endif
		
		#ifdef REG_ACCURACY
		regVariability();
		#endif
		
		std::cout << std::endl;
		if(workingSG->getHomeBarrel())
		{
			std::cout << ">>> Final home barrel: " << SBF->int2Labels[workingSG->getHomeBarrel()] << std::endl;
			TransformLog << "Final home barrel: " << SBF->int2Labels[workingSG->getHomeBarrel()] << std::endl;
		}
		std::cout << "**********************************************************" << std::endl;
		std::cout << "**********************************************************" << std::endl;
		TransformLog.close();
		spatialGraph->clear();
		spatialGraph->mergeSpatialGraph(workingSG);
		spatialGraph->setHomeBarrel(workingSG->getHomeBarrel());
		
		if(neuronAxis) delete [] neuronAxis;
// 		delete [] barrels;
		return 1;
	}
	else
	{
		std::cout << "Error! Input SpatialGraph is empty!" << std::endl;
		return 0;
	}
};

void Registration::neuronMorphologyZAxis()
{
	double * neuronAxis;
	int axisSelect = 0, cellType = 0;
	std::cout << std::endl;
	std::cout << "Use neuron z-Axis:" << std::endl;
	std::cout << "No axis         -> 0" << std::endl;
	std::cout << "Apical dendrite -> 1" << std::endl;
	std::cout << "Axon            -> 2" << std::endl;
	std::cin >> axisSelect;
	while(cellType != SUPRA && cellType != GRAN && cellType != INFRA)
	{
		std::cout << "celltype (supra = 1, gran = 2, infra = 3): ";
		std::cin >> cellType;
	}
	switch(axisSelect)
	{
		case 1:
			std::cout << "Non-stellate cell (0/1)? ";
			bool goodApical;
			std::cin >> goodApical;
			if(cellType == SUPRA || cellType == GRAN)
			{
				if(goodApical)
					neuronAxis = principalAxes(ApicalDendrite, spatialGraph, goodApical);
				else
					neuronAxis = detectMorphZAxis(ApicalDendrite, cellType, spatialGraph);
			}
			else
				neuronAxis = principalAxes(ApicalDendrite, spatialGraph, goodApical);
			break;
			
		case 2:
			neuronAxis = detectMorphZAxis(Axon, cellType, spatialGraph);
			break;
	}
	if(neuronAxis)
		std::flush(std::cout << "neuron z axis = [" << neuronAxis[0] << "," << neuronAxis[1] << "," << neuronAxis[2] <<"]" << std::endl);
	else
		std::flush(std::cout << "neuron z axis invalid!" << std::endl);
};

int Registration::registerToDifferentColumn(const char * outputFilename, const char * newHomeBarrelLabel, bool onlyAxon)
{
	int newHomeBarrel;
	std::string tmpStr(newHomeBarrelLabel);
	std::map< std::string, int >::const_iterator compareIt;
	for(compareIt = SBF->labels2Int.begin(); compareIt != SBF->labels2Int.end(); ++compareIt)
	{
		if(tmpStr.compare(compareIt->first) == 0)
		{
			newHomeBarrel = compareIt->second;
			break;
		}
	}
	if(compareIt == SBF->labels2Int.end())
	{
		std::cout << "New home barrel not recognized. Has to be Greek Arc or A-E row, Arc 1-4" << std::endl;
		return 0;
	}
	if(spatialGraph)
	{
		if(onlyAxon)
		{
			if(!axonFlag)
			{
				std::cout << "Error! Trying to force axon-only mode but there is no axon in input hoc file!" << std::endl;
				std::cout << "Aborting registration..." << std::endl;
				return 0;
			}
			somaFlag = 0;
			dendriteFlag = 0;
			apicalFlag = 0;
			basalFlag = 0;
			std::cout << std::endl;
			std::cout << "Using axon registration mode; ignoring any other cell structures!" << std::endl;
		}
		
		if(!workingSG)
			workingSG = new AmiraSpatialGraph;
		workingSG->mergeSpatialGraph(spatialGraph);
		
		std::cout << "**********************************************************" << std::endl;
		std::cout << "**********************************************************" << std::endl;
		std::cout << ">>> Starting registration to column " << SBF->int2Labels[newHomeBarrel] << std::endl;
		
		std::string transformName(outputFilename);
		transformName += "_newHB_";
		transformName += SBF->int2Labels[newHomeBarrel];
		transformName += "_transform.log";
		TransformLog.open(transformName.c_str());
		TransformLog << "NeuroMap v3.2.1" << std::endl;
		TransformLog << "Updated for BuildingBrains-3D.org" << std::endl;
		TransformLog << "Build " << __DATE__ << " " << __TIME__ << std::endl;
		TransformLog << std::endl;
		
		int oldHomeBarrel;
		if(somaFlag)
			oldHomeBarrel = registeredHomeBarrel();
		else
		{
			oldHomeBarrel = spatialGraph->getHomeBarrel();
			if(!oldHomeBarrel)
			{
				std::cout << "Error! Cannot determine old home barrel. Aborting registration..." << std::endl;
				return 0;
			}
		}
		if(oldHomeBarrel == newHomeBarrel)
		{
			std::cout << "Already registered to " << SBF->int2Labels[newHomeBarrel] << "!" << std::endl;
			double piaDistance = registeredPiaDistance();
			TransformLog << "Home column is identical to " << SBF->int2Labels[newHomeBarrel] << "!" << std::endl;
			TransformLog << std::endl;
			TransformLog << "Registered Pia-soma distance (micron):\t" << piaDistance << std::endl;
			TransformLog << std::endl;
			if(somaFlag)
			{
				double somaPos[3], localOrientation[3] = {0,0,1};
				getPCenterOfStructure(workingSG, Soma, somaPos);
				SBF->localZAxis(somaPos, localOrientation);
				TransformLog << "Local orientation:\t[" << localOrientation[0] << "," << localOrientation[1] << "," << localOrientation[2] << "]" << std::endl;
				TransformLog << std::endl;
			}
			else
			{
				double localOrientation[3] = {0,0,1};
				SBF->localZAxis(newHomeBarrel, localOrientation);
				TransformLog << "Local orientation:\t[" << localOrientation[0] << "," << localOrientation[1] << "," << localOrientation[2] << "]" << std::endl;
				TransformLog << std::endl;
			}
			TransformLog.close();
			// in case there's still a wrong home barrel somewhere in the hoc file...
			spatialGraph->setHomeBarrel(newHomeBarrel);
			return 1;
		}
		workingSG->setHomeBarrel(oldHomeBarrel);
		
		// assume that neuron is already registered, i.e. z-axis is parallel to local z-axis!!!
		double somaPt[3];
		if(somaFlag)
			getPCenterOfStructure(workingSG, Soma, somaPt);
		else
			for(int ii = 0; ii < 3; ++ii)
				somaPt[ii] = SBF->avgCenters[oldHomeBarrel][ii];
// 		double * neuronAxis = localZAxis(somaPt);
		double neuronAxis[3];
		SBF->localZAxis(somaPt, neuronAxis);
		
		// translation and rotation to new location
		TransformPointerType trans = transformToNewColumn(neuronAxis, newHomeBarrel);
		workingSG->setHomeBarrel(newHomeBarrel);
		
		if(somaFlag)
		{
			#ifndef PIPELINE_DOC
			double piaDistance = registeredPiaDistance();
			double somaPos[3], localOrientation[3] = {0,0,1};
			getPCenterOfStructure(workingSG, Soma, somaPos);
			SBF->localZAxis(somaPos, localOrientation);
			TransformLog << "Registered Pia-soma distance (micron):\t" << piaDistance << std::endl;
			TransformLog << std::endl;
			TransformLog << "Local orientation:\t[" << localOrientation[0] << "," << localOrientation[1] << "," << localOrientation[2] << "]" << std::endl;
			TransformLog << std::endl;
			#else
			registeredSomaPosition("z scaling");
			double minDendDist = registeredPiaDendriteDistance();
			TransformLog << "Registered minimum Pia-dendrite distance (micron):\t" << minDendDist << std::endl;
			TransformLog << std::endl;
			#endif
		}
		else
		{
			double localOrientation[3] = {0,0,1};
			SBF->localZAxis(newHomeBarrel, localOrientation);
			TransformLog << "Local orientation:\t[" << localOrientation[0] << "," << localOrientation[1] << "," << localOrientation[2] << "]" << std::endl;
			TransformLog << std::endl;
		}
		
		std::cout << ">>> Final home barrel: " << SBF->int2Labels[newHomeBarrel] << std::endl;
		std::cout << "**********************************************************" << std::endl;
		std::cout << "**********************************************************" << std::endl;
		TransformLog << "Final home barrel: " << SBF->int2Labels[newHomeBarrel] << std::endl;
		TransformLog.close();
		spatialGraph->clear();
		spatialGraph->mergeSpatialGraph(workingSG);
		spatialGraph->setHomeBarrel(newHomeBarrel);
		
		return 1;
	}
	else
		return 0;
};

void Registration::barrelFieldRegistration(const char * outputFilename, int mode, const char * refBarrel)
{
	if(spatialGraph)
	{
		if(!workingSG)
			workingSG = new AmiraSpatialGraph;
		workingSG->mergeSpatialGraph(spatialGraph);
		workingSG->setHomeBarrel(spatialGraph->getHomeBarrel());
		if(!workingSG->getHomeBarrel())
		{
			workingSG->setHomeBarrel(closestMorphBarrel());
			// if there's still no home barrel, abort registration
			if(!workingSG->getHomeBarrel())
			{
				std::cout << "Error! No home barrel in hoc file and not able to determine home barrel automatically!" << std::endl;
				std::cout << "Aborting registration..." << std::endl;
				return;
			}
		}
		
		int refBarrelID = 0;
		if(/*mode == 2 &&*/ refBarrel)
		{
			std::string tmpStr(refBarrel);
			std::map< std::string, int >::const_iterator compareIt;
			for(compareIt = SBF->labels2Int.begin(); compareIt != SBF->labels2Int.end(); ++compareIt)
			{
				if(tmpStr.compare(compareIt->first) == 0)
				{
					refBarrelID = compareIt->second;
					break;
				}
			}
			if(compareIt == SBF->labels2Int.end())
			{
				std::cout << "New home barrel not recognized. Has to be Greek Arc or A-E row, Arc 1-4" << std::endl;
				return;
			}
		}
		else if(mode == 2 && !refBarrel)
		{
			std::cout << "Error! Reference barrel has to be specified for local registration!" << std::endl;
			return;
		}
		PolyDataPointerType pia = NULL;
		PolyDataPointerType wm = NULL;
		std::cout << ">>> Computing landmark reconstruction..." << std::endl;
		std::cout << std::endl;
		if(piaFlag)
			pia = surfaceReconstruction(Pia);
		if(wmFlag)
			wm = surfaceReconstruction(WhiteMatter);
// 		constantLandmarkOffset(pia, wm, -39);
		
		std::map< int, double * > barrelAxes;
		std::map< int, Column * > * barrels = new std::map< int, Column * >[2]; // input barrels = barrels[0], reference barrels = barrels[1]
		double neuronAxis[3];
		barrelReconstruction(pia, neuronAxis, &barrelAxes, barrels);
		//now, add the barrel contours to SpatialGraph for nice visualization
		std::list< int >::const_iterator labelIt;
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(barrels[0].find(ID) != barrels[0].end())
				workingSG->addPolyDataObject(barrels[0][ID]->contours, ID);
		}
		
		std::cout << ">>> Computing landmark registration..." << std::endl;
		TransformPointerType landmarkTransform = landmarkRegistration(barrels, mode, refBarrelID);
		HomogeneousMatrixPointerType mTrans = landmarkTransform->GetMatrix();
		
		// test whether this reduces the RMSE
		localLandmarkOffset(barrels);
		
		std::string ofName(outputFilename);
		std::string ofName2(outputFilename);
		ofName += ".csv";
		ofName2 += "_HB.csv";
		
		std::ofstream HBErrorFile;
		HBErrorFile.open(ofName2.c_str());
		HBErrorFile << "# HB squared 3D error:" << std::endl;
		
		std::ofstream TransformFile;
		TransformFile.open(ofName.c_str());
		TransformFile << "# Landmark transformation matrix" << std::endl;
		TransformFile << "# Mode: " << mode << " (0: global; 1: local)" << std::endl;
		TransformFile << "# Amira format:" << std::endl;
		for(int ii = 0; ii < 4; ++ii)
			for(int jj = 0; jj < 4; ++jj)
				TransformFile << mTrans->GetElement(jj, ii) << " ";
		TransformFile << std::endl;
		TransformFile << "# python/numpy format:" << std::endl;
		TransformFile << "[";
		for(int ii = 0; ii < 4; ++ii)
		{
			TransformFile << "[";
			for(int jj = 0; jj < 4; ++jj)
			{
				TransformFile << mTrans->GetElement(ii, jj);
				if(jj < 3) TransformFile << ",";
			}
			TransformFile << "]";
			if(ii < 3) TransformFile << ",";
		}
		TransformFile << "]" << std::endl;
		TransformFile << std::endl;
		TransformFile << "Barrel\tTop squared distance\tBottom squared distance" << std::endl;
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			TransformFile << SBF->int2Labels[ID] << "\t";
			if(mode == 2)
			{
				// only measure RMSE of neighboring barrels
				if(barrels[0].find(ID) != barrels[0].end() && SBF->avgBarrels.find(ID) != SBF->avgBarrels.end()
					&& std::find(SBF->barrelGrid[refBarrelID].begin(), SBF->barrelGrid[refBarrelID].end(), ID)
					!= SBF->barrelGrid[refBarrelID].end())
				{
					double topDist2 = vtkMath::Distance2BetweenPoints(barrels[0][ID]->top, SBF->avgBarrels[ID]->top);
					double bottomDist2 = vtkMath::Distance2BetweenPoints(barrels[0][ID]->bottom, SBF->avgBarrels[ID]->bottom);
					TransformFile << topDist2 << "\t" << bottomDist2;
				}
			}
			else
			{
// 				if(barrels[0].find(ID) != barrels[0].end() && SBF->avgBarrels.find(ID) != SBF->avgBarrels.end())
// 				{
// 					double topDist2 = vtkMath::Distance2BetweenPoints(barrels[0][ID]->top, SBF->avgBarrels[ID]->top);
// 					double bottomDist2 = vtkMath::Distance2BetweenPoints(barrels[0][ID]->bottom, SBF->avgBarrels[ID]->bottom);
// 					TransformFile << topDist2 << "\t" << bottomDist2;
// 				}
				// only MSE of neighboring barrels
				if(barrels[0].find(ID) != barrels[0].end() && SBF->avgBarrels.find(ID) != SBF->avgBarrels.end()
					&& std::find(SBF->barrelGrid[refBarrelID].begin(), SBF->barrelGrid[refBarrelID].end(), ID)
					!= SBF->barrelGrid[refBarrelID].end())
				{
					double topDist2 = vtkMath::Distance2BetweenPoints(barrels[0][ID]->top, SBF->avgBarrels[ID]->top);
					double bottomDist2 = vtkMath::Distance2BetweenPoints(barrels[0][ID]->bottom, SBF->avgBarrels[ID]->bottom);
					TransformFile << topDist2 << "\t" << bottomDist2;
				}
				// only MSE of home barrel
				if(barrels[0].find(ID) != barrels[0].end() && SBF->avgBarrels.find(ID) != SBF->avgBarrels.end()
					&& ID == refBarrelID)
				{
					double topDist2 = vtkMath::Distance2BetweenPoints(barrels[0][ID]->top, SBF->avgBarrels[ID]->top);
					double bottomDist2 = vtkMath::Distance2BetweenPoints(barrels[0][ID]->bottom, SBF->avgBarrels[ID]->bottom);
					HBErrorFile << topDist2 << "\t" << bottomDist2;
				}
			}
			TransformFile << std::endl;
		}
		
		HBErrorFile.close();
		TransformFile.close();
		
		spatialGraph->clear();
		spatialGraph->mergeSpatialGraph(workingSG);
	}
};

//calculate deviation of automatic barrel detection from manual barrel outlines
//absolute deviation of barrel centers/centroids
void Registration::reconstructionError(const char * filename)
{
	if(!workingSG) workingSG = new AmiraSpatialGraph;
	workingSG->mergeSpatialGraph(spatialGraph);
	
	PolyDataPointerType pia = NULL;
	PolyDataPointerType wm = NULL;
	if(piaFlag)
		pia = surfaceReconstruction(Pia);
	if(wmFlag)
		wm = surfaceReconstruction(WhiteMatter);
	
	std::map< int, PolyDataPointerType > barrels;
	std::map< int, PolyDataPointerType > avgBarrels;
	std::map< int, double * > barrelAxes;
	std::map< int, double * > barrelCenters;
	std::map< int, std::vector< double * > > endPointMap;
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		PolyDataPointerType barrel = PolyDataPointerType::New();
		if(workingSG->extractLandmark(*labelIt, barrel))
		{
			std::cout << "Computing parameters of Barrel " << SBF->int2Labels[*labelIt] << std::endl;
			double alpha = 0.5;
			double * newAxis = newBarrelAxis(barrel, pia, alpha);
			double * barrelCenter = calculateBarrelCentroid(barrel);
			barrels.insert(std::pair< int, PolyDataPointerType >(*labelIt, barrel));
			barrelCenters.insert(std::pair< int, double * >(*labelIt, barrelCenter));
			barrelAxes.insert(std::pair< int, double * >(*labelIt, newAxis));
			
			std::vector< double * > endPoints;
			computeAverageHomeBarrel(barrels[*labelIt], barrelCenters[*labelIt], barrelAxes[*labelIt], *labelIt);
			closeBarrelAlongNewAxis(barrelAxes[*labelIt], barrelCenters[*labelIt], barrels[*labelIt], endPoints, 0);
			endPointMap.insert(std::pair< int, std::vector< double * > >(*labelIt, endPoints));
		}
		else
			std::cout << "Barrel " << SBF->int2Labels[*labelIt] << " not found in SpatialGraph!" << std::endl;
	}
	// apply divergence constraint on barrel axis vector field
// 	enforceAxisDivergence(barrelAxes, barrelCenters);
	
// 	// calculate avg barrel contours and ensure mutual non-overlap
// 	// of those in barrel field (they can overlap in deeper layers...)
// 	std::map< int, std::vector< double * > > endPointMap;
// 	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
// 		if(workingSG->isLabelInSpatialGraph(*labelIt))
// 		{
// 			std::vector< double * > endPoints;
// 			closeBarrelAlongNewAxis(barrelAxes[*labelIt], barrelCenters[*labelIt], barrels[*labelIt], endPoints, 0);
// 			endPointMap.insert(std::pair< int, std::vector< double * > >(*labelIt, endPoints));
// 			computeAverageHomeBarrel(barrels[*labelIt], barrelCenters[*labelIt], barrelAxes[*labelIt]);
// // 			avgBarrels.insert(std::pair< int, PolyDataPointerType >(*labelIt, smoothBarrelAlongNewAxis(barrelAxes[*labelIt], barrelCenters[*labelIt], barrels[*labelIt], endPointMap[*labelIt])));
// 		}
	
	std::string outFile(filename);
	outFile += "_parameter.csv";
	std::ofstream ParameterFile(outFile.c_str());
	ParameterFile << "Barrel\tbottom\ttop\theight\tWM\tarea" << std::endl;
	std::map< int, std::vector< double > > barrelParameters;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		ParameterFile << SBF->int2Labels[*labelIt];
		if(workingSG->isLabelInSpatialGraph(*labelIt))
		{
			PolyDataPointerType barrel = barrels[*labelIt];
			double * finalAxis = barrelAxes[*labelIt];
			double * barrelCenter = barrelCenters[*labelIt];
			std::vector< double * > endPoints = endPointMap[*labelIt];
// 			std::vector< double * > endPoints;
// 			closeBarrelAlongNewAxis(finalAxis, barrelCenter, barrel, endPoints);
			std::vector< double > parameters = computeManualBarrelParameters(barrel, pia, wm, finalAxis, barrelCenter, endPoints, *labelIt, avgBarrelContours);
			for(int ii = 0; ii < parameters.size(); ++ii)
				ParameterFile << "\t" << parameters[ii];
			
			barrelParameters.insert(std::pair< int, std::vector< double > >(*labelIt, parameters));
			workingSG->addPolyDataObject(barrel, *labelIt);
		}
		ParameterFile << std::endl;
	}
	ParameterFile.close();
	
	if(piaFlag)
	{
		std::string piaSurfFile(filename);
		piaSurfFile += "_pia";
		Reader * piaSurfaceWriter = new Reader(piaSurfFile.c_str(), piaSurfFile.c_str());
		piaSurfaceWriter->writeAmiraSurfaceFile(pia);
		delete piaSurfaceWriter;
	}
	
	if(wmFlag)
	{
		std::string wmSurfFile(filename);
		wmSurfFile += "_wm";
		Reader * wmSurfaceWriter = new Reader(wmSurfFile.c_str(), wmSurfFile.c_str());
		wmSurfaceWriter->writeAmiraSurfaceFile(wm);
		delete wmSurfaceWriter;
	}
	
	spatialGraph->clear();
	spatialGraph->mergeSpatialGraph(workingSG);
};

/******************************************************************************/
/*optional: detect z-axis of neuron reconstruction using either the apical    */
/*dendrite or the axon                                                        */
/*useage: apical: 1; axon: 0                                                  */
/******************************************************************************/
double * Registration::detectMorphZAxis(int label, int cellType, AmiraSpatialGraph* morphology)
{
	std::map< double, int > segmentScores;
	double startPt[3];
	PolyDataPointerType soma = PolyDataPointerType::New();
	if(!morphology->extractLandmark(Soma, soma))
	{
		std::cout << "Error! Could not find soma in input SpatialGraph! Cannot determine neuron z-axis." << std::endl;
		return NULL;
	}
	int subID;
	double pCoords[3], * weights;
	weights = new double[soma->GetCell(0)->GetNumberOfPoints()];
	soma->GetCell(0)->GetParametricCenter(pCoords);
	soma->GetCell(0)->EvaluateLocation(subID, pCoords, startPt, weights);
	
	// case: use apical
	if(label == ApicalDendrite)
	{
		if(!morphology->isLabelInSpatialGraph(ApicalDendrite))
		{
			std::cout << "Error! Could not find apical dendrite in input SpatialGraph! Cannot determine neuron z-axis." << std::endl;
			return NULL;
		}
		double alpha = 0.1;
// 		std::cout << "neuron axis alpha = ";
// 		std::cin >> alpha;
		while(cellType != SUPRA && cellType != GRAN && cellType != INFRA)
		{
			std::cout << "celltype (supra = 1, gran = 2, infra = 3): ";
			std::cin >> cellType;
		}
		
		double maxDist = 0;
		for(int ii = 0; ii < (*(morphology->edgesPointer())).size(); ++ii)
			if((*(morphology->edgesPointer()))[ii]->label == ApicalDendrite)
			{
				double tmpDist = L2Distance3D(startPt, (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.back());
				if(tmpDist > maxDist)
					maxDist = tmpDist;
			}
		
		for(int ii = 0; ii < (*(morphology->edgesPointer())).size(); ++ii)
			if((*(morphology->edgesPointer()))[ii]->label == ApicalDendrite)
			{
				double segmentLength = L2Distance3D((*(morphology->edgesPointer()))[ii]->edgePointCoordinates.front(), (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.back());
				double directLength = L2Distance3D(startPt, (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.back());
				double length = morphology->totalSegmentLength(ii);
				double angle = morphology->cumulatedSegmentAngle(ii);
				if(angle && length > 50)
				{
					double score;
					if(cellType == SUPRA)
						score = length/directLength + alpha*angle;
					else if(cellType == GRAN)
						score = length/directLength + maxDist/length + alpha*angle;
					else if(cellType == INFRA)
						score = segmentLength;
// 					std::cout << "ID " << ii << std::endl;
// 					std::cout << "directLength = " << directLength << std::endl;
// 					std::cout << "maxDist = " << maxDist << std::endl;
// 					std::cout << "length = " << length << std::endl;
// 					std::cout << "cumulative angle to soma = " << angle << std::endl;
// 					std::cout << "score = " << score << std::endl;
					segmentScores.insert(std::pair< double, int >(score, ii));
				}
			}
		
		int axisID1 = -1, axisID2 = -1;
		if(cellType == SUPRA || cellType == GRAN)
		{
			std::map< double, int >::const_iterator scoreIt = segmentScores.begin();
			axisID1 = scoreIt->second;
			++scoreIt;
			if(scoreIt != segmentScores.begin())
				axisID2 = scoreIt->second;
		}
		if(cellType == INFRA)
		{
			std::map< double, int >::const_reverse_iterator scoreRIt = segmentScores.rbegin();
			axisID1 = scoreRIt->second;
			++scoreRIt;
			if(scoreRIt != segmentScores.rend())
				axisID2 = scoreRIt->second;
		}
		double * zAxis = new double[3], axis1[3], axis2[3];
		for(int ii = 0; ii < 3; ++ii)
		{
			axis1[ii] = (*(morphology->edgesPointer()))[axisID1]->edgePointCoordinates.back()[ii] - startPt[ii];
			if(axisID2 != -1)
				axis2[ii] = (*(morphology->edgesPointer()))[axisID2]->edgePointCoordinates.back()[ii] - startPt[ii];
		}
		bool averageAxes = 0;
		if(axisID2 != -1)
		{
			double tmpAngle = 0, tmpAxis1[3], tmpAxis2[3];
			for(int ii = 0; ii < 3; ++ii)
			{
				tmpAxis1[ii] = (*(morphology->edgesPointer()))[axisID1]->edgePointCoordinates.front()[ii] - (*(morphology->edgesPointer()))[axisID1]->edgePointCoordinates.back()[ii];
				tmpAxis2[ii] = (*(morphology->edgesPointer()))[axisID2]->edgePointCoordinates.front()[ii] - (*(morphology->edgesPointer()))[axisID2]->edgePointCoordinates.back()[ii];
			}
			normalize(tmpAxis1), normalize(tmpAxis2);
			for(int ii = 0; ii < 3; ++ii)
				tmpAngle += tmpAxis1[ii]*tmpAxis2[ii];
			tmpAngle = acos(tmpAngle)*180/PI;
			if(tmpAngle < 10)
				averageAxes = 1;
		}
		
		morphology->addLine(startPt, (*(morphology->edgesPointer()))[axisID1]->edgePointCoordinates.back());
		
		for(int ii = 0; ii < 3; ++ii)
		{
			zAxis[ii] = axis1[ii];
			if(averageAxes)
				zAxis[ii] += axis2[ii];
		}
		normalize(zAxis);
		
		if(averageAxes)
		{
			morphology->addLine(startPt, (*(morphology->edgesPointer()))[axisID2]->edgePointCoordinates.back());
			double axisEndPt[3];
			for(int ii = 0; ii < 3; ++ii)
				axisEndPt[ii] = startPt[ii] + 1000*zAxis[ii];
			morphology->addLine(startPt, axisEndPt);
		}
		
		return zAxis;
	}
	// case: use axon
	else if(label == Axon)
	{
		if(!morphology->isLabelInSpatialGraph(Axon))
		{
			std::cout << "Error! Could not find axon in input SpatialGraph! Cannot determine neuron z-axis." << std::endl;
			return NULL;
		}
		
		double maxDist = 0, alpha = 0.1, beta = 0.1;
// 		std::cout << "neuron axis alpha = ";
// 		std::cin >> alpha;
// 		std::cout << "neuron axis beta = ";
// 		std::cin >> beta;
// 		std::cout << "maxdist = ";
// 		std::cin >> maxDist;
		for(int ii = 0; ii < (*(morphology->edgesPointer())).size(); ++ii)
			if((*(morphology->edgesPointer()))[ii]->label == Axon)
			{
				double tmpDist = L2Distance3D(startPt, (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.back());
				if(tmpDist > maxDist)
					maxDist = tmpDist;
			}
			
		for(int ii = 0; ii < (*(morphology->edgesPointer())).size(); ++ii)
			if((*(morphology->edgesPointer()))[ii]->label == Axon)
			{
				double segmentLength = L2Distance3D((*(morphology->edgesPointer()))[ii]->edgePointCoordinates.front(), (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.back());
				double directLength = L2Distance3D(startPt, (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.back());
				double length = morphology->totalSegmentLength(ii);
				double angle = morphology->cumulatedSegmentAngle(ii);
				if(angle && length > 200 && directLength < maxDist)
				{
					double score;
// 					score = length/directLength + alpha*angle;
					score = length/directLength + beta*maxDist/length /*+ alpha*angle*/;
// 					std::cout << "ID " << ii << std::endl;
// 					std::cout << "directLength = " << directLength << std::endl;
// 					std::cout << "maxDist = " << maxDist << std::endl;
// 					std::cout << "length = " << length << std::endl;
// 					std::cout << "cumulative angle to soma = " << angle << std::endl;
// 					std::cout << "score = " << score << std::endl;
					segmentScores.insert(std::pair< double, int >(score, ii));
				}
			}
		
		int axisID = segmentScores.begin()->second;
		double * zAxis = new double[3], axisEndPt[3];
		for(int ii = 0; ii < 3; ++ii)
			zAxis[ii] = (*(morphology->edgesPointer()))[axisID]->edgePointCoordinates.back()[ii] - startPt[ii];
		normalize(zAxis);
		// this is called after landmark registration, so the z-axis has to point upwards:
		if(zAxis[2] < 0)
			zAxis[0] = -zAxis[0], zAxis[1] = -zAxis[1], zAxis[2] = -zAxis[2];
// 		if(!zReversed && zAxis[2] < 0)
// 			zAxis[0] = -zAxis[0], zAxis[1] = -zAxis[1], zAxis[2] = -zAxis[2];
// 		else if(zReversed && zAxis[2] > 0)
// 			zAxis[0] = -zAxis[0], zAxis[1] = -zAxis[1], zAxis[2] = -zAxis[2];
		morphology->addLine(startPt, (*(morphology->edgesPointer()))[axisID]->edgePointCoordinates.back());
		
		return zAxis;
	}
	else
	{
		std::cout << "Error! Use axon or apical dendrite to determine neuron z-axis." << std::endl;
		return NULL;
	}
};

double * Registration::principalAxes(int label, AmiraSpatialGraph * morphology, int biTufted)
{
	double startPt[3];
	PolyDataPointerType soma = PolyDataPointerType::New();
	if(!morphology->extractLandmark(Soma, soma))
	{
		std::cout << "Error! Could not find soma in input SpatialGraph! Cannot determine neuron z-axis." << std::endl;
		return NULL;
	}
	int subID;
	double pCoords[3], * weights;
	weights = new double[soma->GetCell(0)->GetNumberOfPoints()];
	soma->GetCell(0)->GetParametricCenter(pCoords);
	soma->GetCell(0)->EvaluateLocation(subID, pCoords, startPt, weights);
	
	std::list< double * > ptList;
	for(int ii = 0; ii < (*(morphology->edgesPointer())).size(); ++ii)
	{
		if(label == ApicalDendrite)
		{
			if(biTufted != 2 && ((*(morphology->edgesPointer()))[ii]->label == label || (*(morphology->edgesPointer()))[ii]->label == label-1))
			{
				std::list< double * >::const_iterator edgePtIt;
				for(edgePtIt = (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.begin(); edgePtIt != (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.end(); ++edgePtIt)
				{
					double * newPt = new double[3];
					newPt[0] = (*edgePtIt)[0], newPt[1] = (*edgePtIt)[1], newPt[2] = (*edgePtIt)[2];
					ptList.push_back(newPt);
				}
			}
			else if(biTufted == 2)
			{
				if((*(morphology->edgesPointer()))[ii]->label == label)
				{
					std::list< double * >::const_iterator edgePtIt;
					for(edgePtIt = (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.begin(); edgePtIt != (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.end(); ++edgePtIt)
					{
						double * newPt = new double[3];
						newPt[0] = (*edgePtIt)[0], newPt[1] = (*edgePtIt)[1], newPt[2] = (*edgePtIt)[2];
						ptList.push_back(newPt);
					}
				}
			}
		}
		else
			if((*(morphology->edgesPointer()))[ii]->label == label)
			{
				std::list< double * >::const_iterator edgePtIt;
				for(edgePtIt = (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.begin(); edgePtIt != (*(morphology->edgesPointer()))[ii]->edgePointCoordinates.end(); ++edgePtIt)
				{
					double * newPt = new double[3];
					newPt[0] = (*edgePtIt)[0], newPt[1] = (*edgePtIt)[1], newPt[2] = (*edgePtIt)[2];
					ptList.push_back(newPt);
				}
			}
	}
	double com[] = {0,0,0};
	std::list< double * >::iterator ptListIt;
	for(ptListIt = ptList.begin(); ptListIt != ptList.end(); ++ptListIt)
		for(int ii = 0; ii < 3; ++ii)
			com[ii] += (*ptListIt)[ii];
	if(ptList.size())
	{
		for(int ii = 0; ii < 3; ++ii)
			com[ii] /= double(ptList.size());
	}
	for(ptListIt = ptList.begin(); ptListIt != ptList.end(); ++ptListIt)
		for(int ii = 0; ii < 3; ++ii)
			(*ptListIt)[ii] -= com[ii];
	
	gsl_matrix * mI = gsl_matrix_calloc(3, 3);
	for(ptListIt = ptList.begin(); ptListIt != ptList.end(); ++ptListIt)
		for(int ii = 0; ii < 3; ++ii)
			for(int jj = ii; jj < 3; ++jj)
			{
				double * pt = *ptListIt;
				double tmp = -pt[ii]*pt[jj];
				if(ii == jj)
					tmp += pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2];
				double oldVal = gsl_matrix_get(mI, ii, jj);
				tmp += oldVal;
				gsl_matrix_set(mI, ii, jj, tmp);
			}
	for(int ii = 0; ii < 3; ++ii)
		for(int jj = ii+1; jj < 3; ++jj)
		{
			double tmp = gsl_matrix_get(mI, ii, jj);
			gsl_matrix_set(mI, jj, ii, tmp);
		}
		
// 	std::cout << "mI:" << std::endl;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 2; ++jj)
// 			std::cout << gsl_matrix_get(mI, ii, jj) << ", ";
// 		std::cout << gsl_matrix_get(mI, ii, 2) << "]" << std::endl;
// 	}
	
	// diagonalization
	gsl_vector * eigenVals = gsl_vector_alloc(3);
	gsl_matrix * eigenVecs = gsl_matrix_alloc(3, 3);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(3);
	int status = gsl_eigen_symmv(mI, eigenVals, eigenVecs, w);
	// eigenvalues are unordered, but eigenvectors are in (column) order corresponding to their respective eigenvalues
	double minEVal = gsl_vector_get(eigenVals, 0);
	int minEValID = 0;
	for(int ii = 1; ii < 3; ++ii)
	{
		double tmp = gsl_vector_get(eigenVals, ii);
		if(tmp < minEVal)
		{
			minEVal = tmp;
			minEValID = ii;
		}
	}
	
	for(int ii = 0; ii < 3; ++ii)
	{
		double tmpAxis[3], tmpAxisEndPt[3], tmpLength;
		tmpLength = gsl_vector_get(eigenVals, ii);
		for(int jj = 0; jj < 3; ++jj)
		{
			tmpAxis[jj] = gsl_matrix_get(eigenVecs, jj, ii);
			tmpAxisEndPt[jj] = com[jj] + 500*tmpAxis[jj];
		}
	}
	
	double * pAxis = new double[3];
	for(int ii = 0; ii < 3; ++ii)
		pAxis[ii] = gsl_matrix_get(eigenVecs, ii, minEValID);
	
	gsl_matrix_free(eigenVecs), gsl_matrix_free(mI);
	gsl_vector_free(eigenVals);
	gsl_eigen_symmv_free(w);
	
	// this is called after landmark registration, so the z-axis has to point upwards:
	if(pAxis[2] < 0)
		pAxis[0] = -pAxis[0], pAxis[1] = -pAxis[1], pAxis[2] = -pAxis[2];
// 	if(!zReversed && pAxis[2] < 0)
// 		pAxis[0] = -pAxis[0], pAxis[1] = -pAxis[1], pAxis[2] = -pAxis[2];
// 	else if(zReversed && pAxis[2] > 0)
// 		pAxis[0] = -pAxis[0], pAxis[1] = -pAxis[1], pAxis[2] = -pAxis[2];
	
	double axisEndPt[3];
	for(int ii = 0; ii < 3; ++ii)
		axisEndPt[ii] = startPt[ii] + 1000*pAxis[ii];
	morphology->addLine(startPt, axisEndPt);
	
	return pAxis;
};

/****************************************************************************/
/*pipeline for reconstruction of pia surface, barrels, axes                 */
/*required for registration to standard barrel field                        */
/****************************************************************************/
void Registration::barrelReconstruction(PolyDataPointerType pia, double * neuronAxis, std::map< int, double * > * barrelAxes, std::map< int, Column * > * barrels)
{
	double alpha = 0.5;
	#ifdef REG_ACCURACY
	alpha = var_alpha;
	#endif
// 	std::cout << "alpha = ";
// 	std::cin >> alpha;
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		PolyDataPointerType barrel = PolyDataPointerType::New();
		if(workingSG->extractLandmark(ID, barrel))
		{
			double * barrelCentroid = calculateBarrelCentroid(barrel);
			double * barrelAxis = NULL;
			if(piaFlag)
				barrelAxis = newBarrelAxis(barrel, pia, alpha);
			else
			{
				barrelAxis = new double[3];
				barrelAxis[0] = 0, barrelAxis[1] = 0, barrelAxis[2] = 1;
			}
#ifdef DEBUG
			std::cout << "Barrel " << SBF->int2Labels[ID] << " axis: [" << barrelAxis[0] << "," << barrelAxis[1] << "," << barrelAxis[2] << "]" << std::endl;
#endif
			barrelAxes->insert(std::pair< int, double * >(ID, barrelAxis));
			std::vector< double * > endPts;
			bool HBRecon = 0;
			closeBarrelAlongNewAxis(barrelAxis, barrelCentroid, barrel, endPts, HBRecon);
			Column * newBarrel = new Column;
			newBarrel->top = endPts[0], newBarrel->bottom = endPts[1];
			newBarrel->contours = barrel;
			barrels->insert(std::pair< int, Column * >(ID, newBarrel));
		}
	}
	// following lines are for one-time barrel reconstruction
// 	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
// 	{
// 		int ID = *labelIt;
// 		if(barrels->find(ID) != barrels->end())
// 			workingSG->addPolyDataObject(barrels[0][ID]->contours, ID);
// 	}
};

/****************************************************************************/
/*pipeline for registration to standard barrel field                        */
/****************************************************************************/
TransformPointerType Registration::landmarkRegistration(std::map< int, Column * > * barrels)
{
	if(barrels->size() && SBF->avgBarrels.size())
	{
		std::list< int > commonLandmarks;
		std::map< int, Column * >::const_iterator refBFIt, matchBFIt;
		for(refBFIt = SBF->avgBarrels.begin(); refBFIt != SBF->avgBarrels.end(); ++refBFIt)
			for(matchBFIt = barrels->begin(); matchBFIt != barrels->end(); ++matchBFIt)
				if(refBFIt->first == matchBFIt->first)
				{
					commonLandmarks.push_back(refBFIt->first);
					std::cout << ">>>        Found common landmark: " << SBF->int2Labels[refBFIt->first] << std::endl;
					break;
				}
		
		std::map< int, Column * > commonBarrels;
		std::list< int >::const_iterator commonLandmarksIt;
		for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt)
		{
			barrels[1].insert(std::pair< int, Column * >(*commonLandmarksIt, SBF->avgBarrels[*commonLandmarksIt]));
			commonBarrels.insert(std::pair< int, Column * >(*commonLandmarksIt, barrels[0][*commonLandmarksIt]));
		}
		
		/*** global registration first, then alignment of HB, then global orientation ***/
		if(commonLandmarks.size() > 2 && !ignoreBarrels)
		{
// 			// should really be COM of common barrels
// 			// fix at some point; not big problem though b/c
// 			// home barrel is aligned in next step
// 			// FIXED 08/24/11 (see commonBarrels above)
			
			// for evaluation of registration:
			// align HB first, then compute optimal
			// rotation of remaining barrels
			// w.r.t. HB!!!
			double refCOM[3], matchCOM[3];
			int HBID = workingSG->getHomeBarrel();
			if(HBID)
			{
				for(int ii  = 0; ii < 3; ++ii)
				{
					refCOM[ii] = SBF->avgCenters[HBID][ii];
					matchCOM[ii] = 0.5*(barrels[0][HBID]->top[ii] + barrels[0][HBID]->bottom[ii]);
				}
			}
			else
			{
				calculateBarrelFieldCentroid(barrels[1], refCOM);
				calculateBarrelFieldCentroid(commonBarrels, matchCOM);
			}
			
// 			std::map< int, double * > refBF, matchBF;
			//using barrel top/bottom as reference (preferred for consistency)
			std::map< int, Column * > refBF, matchBF;
			for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt)
			{
				int ID = *commonLandmarksIt;
				Column * refBarrel = new Column;
				Column * matchBarrel = new Column;
				double * refTop = new double[3], * refBottom = new double[3];
				double * matchTop = new double[3], * matchBottom = new double[3];
				// using barrel centers as reference
// 				double * refCenter = new double[3], * matchCenter = new double[3];
				for(int ii = 0; ii < 3; ++ii)
				{
// 					refCenter[ii] = 0.5*(barrels[1][ID]->top[ii] + barrels[1][ID]->bottom[ii]) - refCOM[ii];
// 					matchCenter[ii] = 0.5*(barrels[0][ID]->top[ii] + barrels[0][ID]->bottom[ii]) - matchCOM[ii];
					refTop[ii] = barrels[1][ID]->top[ii] - refCOM[ii];
					refBottom[ii] = barrels[1][ID]->bottom[ii] - refCOM[ii];
					matchTop[ii] = barrels[0][ID]->top[ii] - matchCOM[ii];
					matchBottom[ii] = barrels[0][ID]->bottom[ii] - matchCOM[ii];
				}
				refBarrel->top = refTop, refBarrel->bottom = refBottom;
				matchBarrel->top = matchTop, matchBarrel->bottom = matchBottom;
				refBF.insert(std::pair< int, Column * >(ID, refBarrel));
				matchBF.insert(std::pair< int, Column * >(ID, matchBarrel));
// 				refBF.insert(std::pair< int, double * >(ID, refCenter));
// 				matchBF.insert(std::pair< int, double * >(ID, matchCenter));
			}
			
			gsl_matrix * mU = computeOptimalRotation(refBF, matchBF);
			
			double invShift[3];
			invShift[0] = -matchCOM[0], invShift[1] = -matchCOM[1], invShift[2] = -matchCOM[2];
			std::list< int >::const_iterator labelIt;
			for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
			{
				int ID = *labelIt;
				if(barrels[0].find(ID) != barrels[0].end())
				{
					barrels[0][ID]->translateColumn(invShift);
					barrels[0][ID]->rotateColumn(mU);
					barrels[0][ID]->translateColumn(refCOM);
				}
			}
			
			#ifdef PIPELINE_DOC
			TransformPointerType outTrans1 = TransformPointerType::New();
			TransformPointerType outTrans2 = TransformPointerType::New();
			TransformPointerType outTrans3 = TransformPointerType::New();
			outTrans1->Translate(invShift);
			outTrans2->SetMatrix(gsl2VtkMatrix(mU));
			outTrans3->Translate(refCOM);
			outTrans2->Concatenate(outTrans1);
			outTrans3->Concatenate(outTrans2);
			writeTransformSG(outTrans3, "_global_alignment");
// // 			writeRegStepParameters("_global_alignment");
// 			registeredSomaPosition("global alignment");
// 			double minDendDist = registeredPiaDendriteDistance();
// 			int colSepFlag = getColumnSeptumFlag();
// 			TransformLog << "Registered minimum Pia-dendrite distance (micron):\t" << minDendDist << std::endl;
// 			TransformLog << "Registered column (1)/septum (0) cell:\t" << colSepFlag << std::endl;
// 			TransformLog << std::endl;
			#endif
			
			// for evaluation of registration:
			// align HB first, then compute optimal
			// rotation of remaining barrels
			// w.r.t. HB!!! -> comment next two lines
// 			HomogeneousMatrixPointerType homeBarrelMat = alignHomeBarrel2(barrels);
// 			HomogeneousMatrixPointerType planeRotMat = alignRemainingBarrels(barrels);
			
			TransformPointerType step1 = TransformPointerType::New();
			TransformPointerType step2 = TransformPointerType::New();
			TransformPointerType step3 = TransformPointerType::New();
// 			TransformPointerType step4 = TransformPointerType::New();
// 			TransformPointerType step5 = TransformPointerType::New();
			step1->Translate(invShift);
			step2->SetMatrix(gsl2VtkMatrix(mU));
			step3->Translate(refCOM);
// 			step4->SetMatrix(homeBarrelMat);
// 			step5->SetMatrix(planeRotMat);
			step2->Concatenate(step1);
			step3->Concatenate(step2);
			step3->Update();
// 			step4->Concatenate(step3);
// 			step4->Update();
// 			step5->Concatenate(step4);
// 			step5->Update();
			
// 			workingSG->setTransformation(step5);
// 			workingSG->setTransformation(step4);
			workingSG->setTransformation(step3);
			workingSG->applyTransformation();
			
// 			#ifdef PIPELINE_DOC
// // 			TransformPointerType outTrans4 = TransformPointerType::New();
// // 			TransformPointerType outTrans5 = TransformPointerType::New();
// // 			outTrans4->SetMatrix(homeBarrelMat);
// // 			outTrans5->SetMatrix(planeRotMat);
// // 			outTrans4->Concatenate(outTrans3);
// // 			outTrans5->Concatenate(outTrans4);
// // 			writeTransformSG(outTrans5, "_local_alignment");
// // 			writeRegStepParameters("_local_alignment");
// 			TransformFilterType globalPiaTrans = TransformFilterType::New();
// 			globalPiaTrans->SetTransform(step3);
// 			globalPiaTrans->SetInput(globalPia);
// 			globalPiaTrans->Update();
// 			globalPia->DeepCopy(globalPiaTrans->GetOutput());
// 			registeredSomaPosition("home barrel alignment", globalPia);
// 			double minDendDist = registeredPiaDendriteDistance(globalPia);
// 			int colSepFlag = getColumnSeptumFlag();
// 			TransformLog << "Registered minimum Pia-dendrite distance (micron):\t" << minDendDist << std::endl;
// 			TransformLog << "Registered column (1)/septum (0) cell:\t" << colSepFlag << std::endl;
// 			TransformLog << std::endl;
// 			#endif
			
// 			workingSG->setTransformation(step5);
// 			workingSG->applyTransformation();
// 			
// 			#ifdef PIPELINE_DOC
// 			registeredSomaPosition("home barrel orientation");
// 			minDendDist = registeredPiaDendriteDistance();
// 			colSepFlag = getColumnSeptumFlag();
// 			TransformLog << "Registered minimum Pia-dendrite distance (micron):\t" << minDendDist << std::endl;
// 			TransformLog << "Registered column (1)/septum (0) cell:\t" << colSepFlag << std::endl;
// 			TransformLog << std::endl;
// 			#endif
// 			step5->Concatenate(step4);
// 			step5->Update();
			
// 			delete [] refCOM, delete [] matchCOM;
			gsl_matrix_free(mU);
// 			return step5;
			return step3;
		}
		else
		{
			HomogeneousMatrixPointerType homeBarrelMat = alignHomeBarrel(barrels);
			HomogeneousMatrixPointerType planeRotMat = alignRemainingBarrels(barrels);
			
			TransformPointerType step1 = TransformPointerType::New();
			TransformPointerType step2 = TransformPointerType::New();
			TransformPointerType step3 = TransformPointerType::New();
			step1->SetMatrix(homeBarrelMat);
			step2->SetMatrix(planeRotMat);
			step2->Concatenate(step1);
			step2->Update();
			
			workingSG->setTransformation(step2);
			workingSG->applyTransformation();
			return step2;
		}
	}
	else
	{
		std::cout << "Error! Barrels cannot be empty!" << std::endl;
		return NULL;
	}
};

TransformPointerType Registration::landmarkRegistration(std::map< int, Column * > * barrels, int mode, int refBarrelID)
{
	if(barrels->size() && SBF->avgBarrels.size())
	{
		std::list< int > commonLandmarks;
		std::map< int, Column * >::const_iterator refBFIt, matchBFIt;
		for(refBFIt = SBF->avgBarrels.begin(); refBFIt != SBF->avgBarrels.end(); ++refBFIt)
			for(matchBFIt = barrels->begin(); matchBFIt != barrels->end(); ++matchBFIt)
				if(refBFIt->first == matchBFIt->first)
				{
					commonLandmarks.push_back(refBFIt->first);
					std::cout << ">>>        Found common landmark: " << SBF->int2Labels[refBFIt->first] << std::endl;
					break;
				}
		
		std::map< int, Column * > commonBarrels;
		std::list< int >::const_iterator commonLandmarksIt;
		for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt)
		{
			barrels[1].insert(std::pair< int, Column * >(*commonLandmarksIt, SBF->avgBarrels[*commonLandmarksIt]));
			commonBarrels.insert(std::pair< int, Column * >(*commonLandmarksIt, barrels[0][*commonLandmarksIt]));
		}
		
		double refCOM[3], matchCOM[3];
		if(mode == 1)
		{
// 			double * tmpRef, * tmpMatch;
// 			tmpMatch = calculateBarrelFieldCentroid(commonBarrels);
// 			tmpRef = calculateBarrelFieldCentroid(barrels[1]);
			double tmpRef[3], tmpMatch[3];
			calculateBarrelFieldCentroid(commonBarrels, tmpMatch);
			 calculateBarrelFieldCentroid(barrels[1], tmpRef);
			for(int ii  = 0; ii < 3; ++ii)
			{
				refCOM[ii] = tmpRef[ii];
				matchCOM[ii] = tmpMatch[ii];
			}
// 			delete [] tmpRef, delete [] tmpMatch;
		}
		else if(mode == 2)
		{
			// int HBID = workingSG->getHomeBarrel();
			for(int ii  = 0; ii < 3; ++ii)
			{
				refCOM[ii] = SBF->avgCenters[refBarrelID][ii];
				matchCOM[ii] = 0.5*(barrels[0][refBarrelID]->top[ii] + barrels[0][refBarrelID]->bottom[ii]);
			}
		}
		else
		{
			std::cout << "Error! Only modes 1 (global optimization) or 2 (local optimization) valid!" << std::endl;
			TransformPointerType id = TransformPointerType::New();
			id->Identity();
			return id;
		}
		
		std::map< int, Column * > refBF, matchBF;
		for(commonLandmarksIt = commonLandmarks.begin(); commonLandmarksIt != commonLandmarks.end(); ++commonLandmarksIt)
		{
			int ID = *commonLandmarksIt;
			Column * refBarrel = new Column;
			Column * matchBarrel = new Column;
			double * refTop = new double[3], * refBottom = new double[3];
			double * matchTop = new double[3], * matchBottom = new double[3];
			for(int ii = 0; ii < 3; ++ii)
			{
				refTop[ii] = barrels[1][ID]->top[ii] - refCOM[ii];
				refBottom[ii] = barrels[1][ID]->bottom[ii] - refCOM[ii];
				matchTop[ii] = barrels[0][ID]->top[ii] - matchCOM[ii];
				matchBottom[ii] = barrels[0][ID]->bottom[ii] - matchCOM[ii];
			}
			refBarrel->top = refTop, refBarrel->bottom = refBottom;
			matchBarrel->top = matchTop, matchBarrel->bottom = matchBottom;
			refBF.insert(std::pair< int, Column * >(ID, refBarrel));
			matchBF.insert(std::pair< int, Column * >(ID, matchBarrel));
		}
		
		gsl_matrix * mU = computeOptimalRotation(refBF, matchBF);
		
		double invShift[3];
		invShift[0] = -matchCOM[0], invShift[1] = -matchCOM[1], invShift[2] = -matchCOM[2];
		std::list< int >::const_iterator labelIt;
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(barrels[0].find(ID) != barrels[0].end())
			{
				barrels[0][ID]->translateColumn(invShift);
				barrels[0][ID]->rotateColumn(mU);
				barrels[0][ID]->translateColumn(refCOM);
			}
		}
				
		TransformPointerType step1 = TransformPointerType::New();
		TransformPointerType step2 = TransformPointerType::New();
		TransformPointerType step3 = TransformPointerType::New();
		step1->Translate(invShift);
		step2->SetMatrix(gsl2VtkMatrix(mU));
		step3->Translate(refCOM);
		step2->Concatenate(step1);
		step3->Concatenate(step2);
		step3->Update();
		
		workingSG->setTransformation(step3);
		workingSG->applyTransformation();
		
		gsl_matrix_free(mU);
		return step3;
	}
	else
	{
		std::cout << "Error! Barrels cannot be empty!" << std::endl;
		return NULL;
	}
};

void Registration::morphologyRegistration(double * neuronAxis, PolyDataPointerType pia, PolyDataPointerType wm, std::map< int, Column * > * barrels)
{
	if(workingSG)
	{
		int HBID = workingSG->getHomeBarrel();
#ifdef DEBUG
		std::flush(std::cout << "HBID = " << HBID << std::endl);
#endif
		// rotation if neuron z axis available
		if(neuronAxis)
		{
			double somaPt[3];
			if(somaFlag)
				getPCenterOfStructure(workingSG, Soma, somaPt);
			
			// version 1: rotation around soma
			// possibly systematic error in soma xy-registration?
			double invSomaPt[3];
			invSomaPt[0] = -somaPt[0], invSomaPt[1] = -somaPt[1], invSomaPt[2] = -somaPt[2];
// 			double * newZAxis = localZAxis(somaPt);
			double newZAxis[3];
			SBF->localZAxis(somaPt, newZAxis);
			HomogeneousMatrixPointerType mRot = transformToBarrelCoordinates(neuronAxis, newZAxis);
			TransformPointerType step1 = TransformPointerType::New();
			TransformPointerType step2 = TransformPointerType::New();
			TransformPointerType step3 = TransformPointerType::New();
			step1->Translate(invSomaPt);
			step2->SetMatrix(mRot);
			step3->Translate(somaPt);
			step2->Concatenate(step1);
			step3->Concatenate(step2);
			
			#ifdef PIPELINE_DOC
			TransformPointerType outTrans1 = TransformPointerType::New();
			TransformPointerType outTrans2 = TransformPointerType::New();
			TransformPointerType outTrans3 = TransformPointerType::New();
			outTrans1->Translate(invSomaPt);
			outTrans2->SetMatrix(mRot);
			outTrans3->Translate(somaPt);
			outTrans2->Concatenate(outTrans1);
			outTrans3->Concatenate(outTrans2);
			writeTransformSG(outTrans3, "_neuron_axis_alignment");
// 			writeRegStepParameters("_neuron_axis_alignment");
			#endif
			
			for(int label = Neuron; label <= Soma; ++label)
				if(workingSG->isLabelInSpatialGraph(label))
				{
					workingSG->setTransformation(step3);
					workingSG->applyTransformation(label);
				}
			if(workingSG->isLabelInSpatialGraph(ZAxis))
			{
				workingSG->setTransformation(step3);
				workingSG->applyTransformation(ZAxis);
			}
			
			TransformLog << "Neuron z-axis correction:" << std::endl;
// 			for(int ii = 0; ii < 4; ++ii)
// 			{
// 				TransformLog << "[";
// 				for(int jj = 0; jj < 4; ++jj)
// 				{
// 					if(jj < 3)
// 						TransformLog << step3->GetMatrix()->GetElement(ii, jj) << " ";
// 					else
// 						TransformLog << step3->GetMatrix()->GetElement(ii, jj);
// 				}
// 				TransformLog << "]" << std::endl;
// 			}
// 			TransformLog << std::endl;
			
			TransformLog << "# Transformation matrix Amira format" << std::endl;
			for(int ii = 0; ii < 4; ++ii)
				for(int jj = 0; jj < 4; ++jj)
					TransformLog << step3->GetMatrix()->GetElement(jj, ii) << " ";
				TransformLog << std::endl;
			TransformLog << std::endl;
			TransformLog << "# Transformation matrix python/numpy format:" << std::endl;
			TransformLog << "[";
			for(int ii = 0; ii < 4; ++ii)
			{
				TransformLog << "[";
				for(int jj = 0; jj < 4; ++jj)
				{
					TransformLog << step3->GetMatrix()->GetElement(ii, jj);
					if(jj < 3) TransformLog << ",";
				}
				TransformLog << "]";
				if(ii < 3) TransformLog << ",";
			}
			TransformLog << "]" << std::endl;
			TransformLog << std::endl;
			
// 			#ifdef PIPELINE_DOC
// 			registeredSomaPosition("neuron axis alignment", globalPia);
// 			double minDendDist = registeredPiaDendriteDistance(globalPia);
// 			int colSepFlag = getColumnSeptumFlag();
// 			TransformLog << "Registered minimum Pia-dendrite distance (micron):\t" << minDendDist << std::endl;
// 			TransformLog << "Registered column (1)/septum (0) cell:\t" << colSepFlag << std::endl;
// 			TransformLog << std::endl;
// 			#endif
			
// 			delete [] newZAxis;
		}
		#ifdef PIPELINE_DOC
		if(piaFlag)
		{
			std::string piaOutputName(pipelineDocName);
			piaOutputName += "_pia_global_reg";
			Reader * piaWriter = new Reader(piaOutputName.c_str(), piaOutputName.c_str());
			piaWriter->writeAmiraSurfaceFile(pia);
			delete piaWriter;
		}
		if(wmFlag)
		{
			std::string wmOutputName(pipelineDocName);
			wmOutputName += "_wm_global_reg";
			Reader * wmWriter = new Reader(wmOutputName.c_str(), wmOutputName.c_str());
			wmWriter->writeAmiraSurfaceFile(wm);
			delete wmWriter;
		}
		#endif
		if(correctPia && piaFlag)
		{
			TransformLog << "Artificial Pia correction [micron]:\t" << correctPiaDist << std::endl;
			TransformLog << std::endl;
			correctArtificialPia(pia);
		}
		
		if(contourCorrection.find(HBID) != contourCorrection.end())
			if(contourCorrection[HBID].size() == 3)
			{
				TransformLog << "Manual contour correction distances [micron]:" << std::endl;
// 				TransformLog << "Manual contour correction scale:" << std::endl;
				TransformLog << "Pia-L4 upper\tL4 upper-L4 lower\tL4 lower-WM" << std::endl;
				TransformLog << contourCorrection[HBID][0] << "\t" << contourCorrection[HBID][1] << "\t" << contourCorrection[HBID][2] << std::endl;
				TransformLog << std::endl;
				correctManualContours(pia, wm, barrels);
				
				#ifdef PIPELINE_DOC
				std::string quickName(pipelineDocName);
				quickName += "_sys_contour_correction";
				Reader * quickWriter = new Reader(quickName.c_str(), quickName.c_str());
				quickWriter->setSpatialGraph(workingSG);
				quickWriter->writeSpatialGraphFile();
				delete quickWriter;
// 				registeredSomaPosition("systematic contour correction", globalPia);
// 				double minDendDist = registeredPiaDendriteDistance(globalPia);
// 				int colSepFlag = getColumnSeptumFlag();
// 				TransformLog << "Registered minimum Pia-dendrite distance (micron):\t" << minDendDist << std::endl;
// 				TransformLog << "Registered column (1)/septum (0) cell:\t" << colSepFlag << std::endl;
// 				TransformLog << std::endl;
				#endif
			}
		
		// inhomogeneous z scaling; whole barrel field version; 4 pt registration
		// ref pts = pia, L4upper, L4lower, WM (z coords)
		// scale -> pia-L4u, L4u-L4l, L4l-WM
		double refTopDist, refBarrelHeight, refBottomDist;
		double inputTopDist = 0, inputBarrelHeight = 0, inputBottomDist = 0;
		double zScale[] = {1,1,1};
		double totalScale = 0;
		
		// compute scaling factors based on home barrel
		if(HBID)
		{
			double tmpZAxis[3];
			for(int ii = 0; ii < 3; ++ii)
				tmpZAxis[ii] = barrels[0][HBID]->top[ii] - barrels[0][HBID]->bottom[ii];
			normalize(tmpZAxis);
			double * inputPiaPt = NULL;
			double * inputWMPt = NULL;
			if(piaFlag)
				inputPiaPt = axisSurfaceIntersection(pia, tmpZAxis, barrels[0][HBID]->top);
			if(wmFlag)
				inputWMPt = axisSurfaceIntersection(wm, tmpZAxis, barrels[0][HBID]->bottom);
			
			inputBarrelHeight = barrels[0][HBID]->getHeight();
			if(inputPiaPt)
				inputTopDist = L2Distance3D(inputPiaPt, barrels[0][HBID]->top);
			if(inputWMPt)
				inputBottomDist = L2Distance3D(barrels[0][HBID]->bottom, inputWMPt);
			
			refTopDist = SBF->avgTopDist[HBID];
			refBarrelHeight = SBF->avgBarrels[HBID]->getHeight();
			refBottomDist = SBF->avgPiaWMDist[HBID] - refTopDist - refBarrelHeight;
			
			if(inputTopDist)
			{
				zScale[0] = refTopDist/inputTopDist;
#ifdef DEBUG
				std::cout << "Determine HB scale factors:" << std::endl;
				std::cout << "refTopDist = " << refTopDist << std::endl;
				std::cout << "inputTopDist = " << inputTopDist << std::endl;
				std::cout << "zScale[0] = " << zScale[0] << std::endl;
#endif
			}
			if(inputBarrelHeight)
			{
				zScale[1] = refBarrelHeight/inputBarrelHeight;
			}
			if(inputBottomDist)
			{
				zScale[2] = refBottomDist/inputBottomDist;
			}
		
			if(inputPiaPt) delete [] inputPiaPt;
			if(inputWMPt) delete [] inputWMPt;
		}
		// compute average scaling factors
		else
		{
#ifdef DEBUG
			std::flush(std::cout << "Computing average z-scale parameters..." << std::endl);
#endif
			double supraScaleAvg = 0, granScaleAvg = 0, infraScaleAvg = 0;
			double supraScaleCnt = 0, granScaleCnt = 0, infraScaleCnt = 0;
			std::map< int, Column * >::const_iterator barrelIt;
			for(barrelIt = barrels->begin(); barrelIt != barrels->end(); ++barrelIt)
			{
				int barrelID = barrelIt->first;
				if(SBF->avgBarrels.find(barrelID) != SBF->avgBarrels.end())
				{
#ifdef DEBUG
					std::flush(std::cout << "Checking barrel column " << SBF->int2Labels[barrelIt->first] << std::endl);
#endif
					double tmpRefTopDist, tmpRefBarrelHeight, tmpRefBottomDist;
					double tmpInputTopDist = 0, tmpInputBarrelHeight = 0, tmpInputBottomDist = 0;
					double tmpZAxis[3];
					for(int ii = 0; ii < 3; ++ii)
						tmpZAxis[ii] = barrels[0][barrelID]->top[ii] - barrels[0][barrelID]->bottom[ii];
					normalize(tmpZAxis);
					double * inputPiaPt = NULL;
					double * inputWMPt = NULL;
					if(piaFlag)
						inputPiaPt = axisSurfaceIntersection(pia, tmpZAxis, barrels[0][barrelID]->top);
					if(wmFlag)
						inputWMPt = axisSurfaceIntersection(wm, tmpZAxis, barrels[0][barrelID]->bottom);
					
					tmpInputBarrelHeight = barrels[0][barrelID]->getHeight();
					if(inputPiaPt)
						tmpInputTopDist = L2Distance3D(inputPiaPt, barrels[0][barrelID]->top);
					if(inputWMPt)
						tmpInputBottomDist = L2Distance3D(barrels[0][barrelID]->bottom, inputWMPt);
					
					tmpRefTopDist = SBF->avgTopDist[barrelID];
					tmpRefBarrelHeight = SBF->avgBarrels[barrelID]->getHeight();
					tmpRefBottomDist = SBF->avgPiaWMDist[barrelID] - tmpRefTopDist - tmpRefBarrelHeight;
					
					if(tmpInputTopDist)
					{
						supraScaleAvg += tmpRefTopDist/tmpInputTopDist;
						++supraScaleCnt;
					}
					if(tmpInputBarrelHeight)
					{
						granScaleAvg += tmpRefBarrelHeight/tmpInputBarrelHeight;
						++granScaleCnt;
					}
					if(tmpInputBottomDist)
					{
						infraScaleAvg += tmpRefBottomDist/tmpInputBottomDist;
						++infraScaleCnt;
					}
				
					if(inputPiaPt) delete [] inputPiaPt;
					if(inputWMPt) delete [] inputWMPt;
				}
			}
			if(supraScaleCnt)
			{
				zScale[0] = supraScaleAvg/supraScaleCnt;
			}
			if(granScaleCnt)
			{
				zScale[1] = granScaleAvg/granScaleCnt;
			}
			if(infraScaleCnt)
			{
				zScale[2] = infraScaleAvg/infraScaleCnt;
			}
		}
		
		// somaRefPt is used as "direction" reference for small corrections
		// computed in buildLocalZScaleCorrection() and inhomogeneousZScaling().
		// does not necessarily need to be soma location;
		// central point within SBF is sufficient
		double somaRefPt[3];
		if(somaFlag)
		{
			getPCenterOfStructure(workingSG, Soma, somaRefPt);
		}
		else
		{
			if(HBID)
			{
				for(int ii = 0; ii < 3; ++ii)
					somaRefPt[ii] = SBF->avgCenters[HBID][ii];
			}
			// in case there is no home barrel, 
			else
			{
				for(int ii = 0; ii < 3; ++ii)
					somaRefPt[ii] = SBF->avgCenters[C2][ii];
			}
		}
		
		std::flush(std::cout << std::endl);
		std::flush(std::cout << ">>>        Checking morphology for z-scale correction" << std::endl);
		bool zScaleModified = 0;
		double oldZScale = zScale[0];
		std::map< vtkIdType, std::vector< double > > scaleCorrection = buildLocalZScaleCorrection(zScale, somaRefPt);
		// for now: apply z-scale correction only to dendrite 
		std::map< vtkIdType, std::vector< double > >::const_iterator scaleCorrectionIt;
		for(scaleCorrectionIt = scaleCorrection.begin(); scaleCorrectionIt != scaleCorrection.end(); ++scaleCorrectionIt)
			if(scaleCorrectionIt->second[0] < zScale[0])
			{
				zScale[0] = scaleCorrectionIt->second[0];
				zScaleModified = 1;
			}
		
		std::cout << ">>>        Z-scale factors: [" << zScale[0] << "," << zScale[1] << "," << zScale[2] << "]" << std::endl;
		if(HBID)
		{
			// refTopDist = inputTopDist*zScale[0] in case z scale was globally adjusted
			totalScale = ((inputTopDist ? inputTopDist*zScale[0] : refTopDist) + refBarrelHeight + refBottomDist);
			totalScale /= ((inputTopDist ? inputTopDist : refTopDist) + inputBarrelHeight + (inputBottomDist ? inputBottomDist : refBottomDist));
			
			std::cout << ">>>        Total column scale: " << totalScale << std::endl;
			TransformLog << "Unscaled distances:" << std::endl;
			TransformLog << "Pia-L4 upper\tL4 upper-L4 lower\tL4 lower-WM" << std::endl;
			TransformLog << inputTopDist << "\t" << inputBarrelHeight << "\t" << inputBottomDist << std::endl;
			TransformLog << "Target distances:" << std::endl;
			TransformLog << "Pia-L4 upper\tL4 upper-L4 lower\tL4 lower-WM" << std::endl;
			TransformLog << refTopDist << "\t" << refBarrelHeight << "\t" << refBottomDist << std::endl;
			TransformLog << "Morphology scale factors:" << std::endl;
			TransformLog << "Pia-L4 upper\tL4 upper-L4 lower\tL4 lower-WM\tTotal column scale" << std::endl;
			TransformLog << zScale[0] << "\t" << zScale[1] << "\t" << zScale[2] << "\t" << totalScale << std::endl;
		}
		else
		{
			TransformLog << "Using average scale factors!" << std::endl;
			TransformLog << "Morphology scale factors:" << std::endl;
			TransformLog << "Pia-L4 upper\tL4 upper-L4 lower\tL4 lower-WM" << std::endl;
			TransformLog << zScale[0] << "\t" << zScale[1] << "\t" << zScale[2] << std::endl;
		}
		if(zScaleModified)
		{
			TransformLog << "WARNING: Pia-L4 upper scale was modified to prevent parts of the cell from sticking out of pia." << std::endl;
			TransformLog << "Pia-L4 upper scale before was:\t" << oldZScale << std::endl;
		}
		TransformLog << std::endl;
		
		inhomogeneousZScaling(zScale, somaRefPt);
// 		inhomogeneousZScalingPRL(zScale, somaRefPt);
		
		// inhomogeneous z scaling; single barrel version
		/*** for visualization purposes only: BEGIN scale landmarks globally ***/
		if(piaFlag && HBID)
		{
			tmpCoordinateSystem(pia, wm, barrels);
			double globalZAxis[] = {0,0,1};
			double * landmarkPiaPt = axisSurfaceIntersection(pia, globalZAxis, barrels[0][HBID]->top);
			double * landmarkWMPt = NULL;
			if(wm)
				landmarkWMPt = axisSurfaceIntersection(wm, globalZAxis, barrels[0][HBID]->bottom);
			
			// ref pts = pia, L4upper, L4lower, WM (z coords)
			// scale -> pia-L4u, L4u-L4l, L4l-WM
			double inRefPts[4] = {0,0,0,0};
			double refPoints[4] = {0,0,0,0};
			
			inRefPts[0] = landmarkPiaPt[2];
			inRefPts[1] = barrels[0][HBID]->top[2];
			inRefPts[2] = barrels[0][HBID]->bottom[2];
			if(landmarkWMPt)
				inRefPts[3] = landmarkWMPt[2];
			else
				inRefPts[3] = inRefPts[2];
			
			refPoints[0] = 0.5*refBarrelHeight + refTopDist;
			refPoints[1] = 0.5*refBarrelHeight;
			refPoints[2] = -0.5*refBarrelHeight;
			if(landmarkWMPt)
				refPoints[3] = refPoints[2] - refBottomDist;
			else
				refPoints[3] = refPoints[2];
			
			inhomogeneousZScaling(zScale, refPoints, inRefPts, pia, wm);
			
			tmpCoordinateSystemInv(pia, wm, barrels);
			
			delete [] landmarkPiaPt;
			if(landmarkWMPt) delete [] landmarkWMPt;
		}
		/*** END scale landmarks globally ***/
	}
};

TransformPointerType Registration::transformToNewColumn(double * neuronAxis, int newHomeBarrel)
{
	if(workingSG)
	{
		int oldHomeBarrel = workingSG->getHomeBarrel();
		
		// vertically, preserve soma position in relative position
		// in supra-/granular/infragranular region
		// because the relative thicknesses of those layers
		// are not preserved across columns
// 		double newPiaL4Dist = L2Distance3D(SBF->avgColumns[newHomeBarrel]->top, SBF->avgBarrels[newHomeBarrel]->top);
// 		double oldPiaL4Dist = L2Distance3D(SBF->avgColumns[oldHomeBarrel]->top, SBF->avgBarrels[oldHomeBarrel]->top);
// 		double newL4Height = SBF->avgBarrels[newHomeBarrel]->getHeight();
// 		double oldL4Height = SBF->avgBarrels[oldHomeBarrel]->getHeight();
// 		double newL4WMDist = L2Distance3D(SBF->avgColumns[newHomeBarrel]->bottom, SBF->avgBarrels[newHomeBarrel]->bottom);
// 		double oldL4WMDist = L2Distance3D(SBF->avgColumns[oldHomeBarrel]->bottom, SBF->avgBarrels[oldHomeBarrel]->bottom);
// 		double newPiaWMDist = SBF->avgColumns[newHomeBarrel]->getHeight();
// 		double oldPiaWMDist = SBF->avgColumns[oldHomeBarrel]->getHeight();
		double newPiaL4Dist, oldPiaL4Dist;
		double newL4Height, oldL4Height;
		double newL4WMDist, oldL4WMDist;
		double newPiaWMDist, oldPiaWMDist;
		
// 		// inhomogeneous z-scale
// 		// ref pts = pia, L4upper, L4lower, WM (z coords)
// 		// scale -> pia-L4u, L4u-L4l, L4l-WM
// 		double zScale[] = {1,1,1};
// 		double totalScale = 0;
// 		zScale[0] = newPiaL4Dist/oldPiaL4Dist;
// 		zScale[1] = newL4Height/oldL4Height;
// 		zScale[2] = newL4WMDist/oldL4WMDist;
// 		totalScale = SBF->avgColumns[newHomeBarrel]->getHeight()/SBF->avgColumns[oldHomeBarrel]->getHeight();
		
		double somaCenter[3], oldZAxis[3];
		double newSomaPt[3], invSomaPt[3];
		
		if(somaFlag)
		{
			double piaIntersectPt[3], wmIntersectPt[3], L4UIntersectPt[3], L4LIntersectPt[3];
			getPCenterOfStructure(workingSG, Soma, somaCenter);
			
			SBF->localZAxis(somaCenter, oldZAxis);
			SBF->avgPiaSurface->intersectLine(oldZAxis, somaCenter);
			SBF->avgPiaSurface->getLastIntersectPoint(piaIntersectPt);
			SBF->avgL4UpperSurface->intersectLine(oldZAxis, somaCenter);
			SBF->avgL4UpperSurface->getLastIntersectPoint(L4UIntersectPt);
			SBF->avgL4LowerSurface->intersectLine(oldZAxis, somaCenter);
			SBF->avgL4LowerSurface->getLastIntersectPoint(L4LIntersectPt);
			SBF->avgWMSurface->intersectLine(oldZAxis, somaCenter);
			SBF->avgWMSurface->getLastIntersectPoint(wmIntersectPt);
			oldPiaL4Dist = sqrt(vtkMath::Distance2BetweenPoints(piaIntersectPt, L4UIntersectPt));
			oldL4Height = sqrt(vtkMath::Distance2BetweenPoints(L4UIntersectPt, L4LIntersectPt));
			oldL4WMDist = sqrt(vtkMath::Distance2BetweenPoints(L4LIntersectPt, wmIntersectPt));
			oldPiaWMDist = sqrt(vtkMath::Distance2BetweenPoints(piaIntersectPt, wmIntersectPt));
			
			double projSomaPt[3], radialSomaPt[3], newProjSomaPt[3], t;
			double radialDist = vtkLine::DistanceToLine(somaCenter, SBF->avgColumns[oldHomeBarrel]->top, SBF->avgColumns[oldHomeBarrel]->bottom, t, projSomaPt);
			radialDist = sqrt(radialDist);
			double radialRatio = radialDist/sqrt(SBF->avgBarrelArea[oldHomeBarrel]/PI);
			
			double verticalDist, newVerticalDist;
			if(projSomaPt[2] > SBF->avgBarrels[oldHomeBarrel]->top[2])
			{
				verticalDist = L2Distance3D(projSomaPt, SBF->avgBarrels[oldHomeBarrel]->top);
				// + necessary correction b/c the barrel centers are made to overlap virtually
				newVerticalDist = verticalDist + 0.5*oldL4Height;
			}
			else if(projSomaPt[2] > SBF->avgBarrels[oldHomeBarrel]->bottom[2])
			{
				verticalDist = L2Distance3D(projSomaPt, SBF->avgCenters[oldHomeBarrel]);
				// + necessary correction b/c the barrel centers are made to overlap virtually
				newVerticalDist = verticalDist;
// 				// addtl correction enforcing that soma stays in correct layer
// 				if(newL4Height < oldL4Height)
// 					newVerticalDist -= 0.5*(oldL4Height - newL4Height);
				if(projSomaPt[2] < SBF->avgCenters[oldHomeBarrel][2])
					newVerticalDist = -newVerticalDist;
			}
			else
			{
				verticalDist = L2Distance3D(projSomaPt, SBF->avgBarrels[oldHomeBarrel]->bottom);
				// + necessary correction b/c the barrel centers are made to overlap virtually
				newVerticalDist = -(verticalDist + 0.5*oldL4Height);
			}
			
			for(int ii = 0; ii < 3; ++ii)
			{
				radialSomaPt[ii] = somaCenter[ii] - projSomaPt[ii];
				newProjSomaPt[ii] = SBF->avgCenters[newHomeBarrel][ii] + newVerticalDist*SBF->avgAxes[newHomeBarrel][ii];
			}
			normalize(radialSomaPt);
			
			// move soma radially only so far that the position
			// relative to (circular) column boundary is preserved
			// otherwise, could be a disaster to map e.g. an E2 cell to A1
			// This is still a little inaccurate for cells far from the
			// column axis depending on whether they are moved along
			// directions of low or high curvature
			// Better solution = ???
			double newRadius = radialRatio*sqrt(SBF->avgBarrelArea[newHomeBarrel]/PI);
			for(int ii = 0; ii < 3; ++ii)
			{
				invSomaPt[ii] = -somaCenter[ii];
				newSomaPt[ii] = newProjSomaPt[ii] + newRadius*radialSomaPt[ii];
			}
		}
		else
		{
			for(int ii = 0; ii < 3; ++ii)
			{
				somaCenter[ii] = SBF->avgCenters[oldHomeBarrel][ii];
				invSomaPt[ii] = -somaCenter[ii];
				newSomaPt[ii] = SBF->avgCenters[newHomeBarrel][ii];
			}
			oldPiaL4Dist = L2Distance3D(SBF->avgColumns[oldHomeBarrel]->top, SBF->avgBarrels[oldHomeBarrel]->top);
			oldL4Height = SBF->avgBarrels[oldHomeBarrel]->getHeight();
			oldL4WMDist = L2Distance3D(SBF->avgColumns[oldHomeBarrel]->bottom, SBF->avgBarrels[oldHomeBarrel]->bottom);
			oldPiaWMDist = SBF->avgColumns[oldHomeBarrel]->getHeight();
		}
		
// 		double * newZAxis = localZAxis(newSomaPt);
		double newZAxis[3];
		SBF->localZAxis(newSomaPt, newZAxis);
		
		double newPiaIntersectPt[3], newWMIntersectPt[3], newL4UIntersectPt[3], newL4LIntersectPt[3];
		SBF->avgPiaSurface->intersectLine(newZAxis, newSomaPt);
		SBF->avgPiaSurface->getLastIntersectPoint(newPiaIntersectPt);
		SBF->avgL4UpperSurface->intersectLine(newZAxis, newSomaPt);
		SBF->avgL4UpperSurface->getLastIntersectPoint(newL4UIntersectPt);
		SBF->avgL4LowerSurface->intersectLine(newZAxis, newSomaPt);
		SBF->avgL4LowerSurface->getLastIntersectPoint(newL4LIntersectPt);
		SBF->avgWMSurface->intersectLine(newZAxis, newSomaPt);
		SBF->avgWMSurface->getLastIntersectPoint(newWMIntersectPt);
		newPiaL4Dist = sqrt(vtkMath::Distance2BetweenPoints(newPiaIntersectPt, newL4UIntersectPt));
		newL4Height = sqrt(vtkMath::Distance2BetweenPoints(newL4UIntersectPt, newL4LIntersectPt));
		newL4WMDist = sqrt(vtkMath::Distance2BetweenPoints(newL4LIntersectPt, newWMIntersectPt));
		newPiaWMDist = sqrt(vtkMath::Distance2BetweenPoints(newPiaIntersectPt, newWMIntersectPt));
		
		// inhomogeneous z-scale
		// ref pts = pia, L4upper, L4lower, WM (z coords)
		// scale -> pia-L4u, L4u-L4l, L4l-WM
		double zScale[] = {1,1,1};
		double totalScale = 0;
		zScale[0] = newPiaL4Dist/oldPiaL4Dist;
		zScale[1] = newL4Height/oldL4Height;
		zScale[2] = newL4WMDist/oldL4WMDist;
		totalScale = newPiaWMDist/oldPiaWMDist;
		
		HomogeneousMatrixPointerType mRot = transformToBarrelCoordinates(neuronAxis, newZAxis);
		TransformPointerType step1 = TransformPointerType::New();
		TransformPointerType step2 = TransformPointerType::New();
		TransformPointerType step3 = TransformPointerType::New();
		step1->Translate(invSomaPt);
		step2->SetMatrix(mRot);
		step3->Translate(newSomaPt);
		step2->Concatenate(step1);
		step3->Concatenate(step2);
		step3->Update();
		
		for(int label = Neuron; label < Landmark; ++label)
		{
			workingSG->setTransformation(step3);
			workingSG->applyTransformation(label);
		}
		
		TransformLog << "Neuron morphology transformation:" << std::endl;
// 		for(int ii = 0; ii < 4; ++ii)
// 		{
// 			TransformLog << "[";
// 			for(int jj = 0; jj < 4; ++jj)
// 			{
// 				if(jj < 3)
// 					TransformLog << step3->GetMatrix()->GetElement(ii, jj) << " ";
// 				else
// 					TransformLog << step3->GetMatrix()->GetElement(ii, jj);
// 			}
// 			TransformLog << "]" << std::endl;
// 		}
		
		TransformLog << "# Transformation matrix Amira format" << std::endl;
		for(int ii = 0; ii < 4; ++ii)
			for(int jj = 0; jj < 4; ++jj)
				TransformLog << step3->GetMatrix()->GetElement(jj, ii) << " ";
			TransformLog << std::endl;
		TransformLog << std::endl;
		TransformLog << "# Transformation matrix python/numpy format:" << std::endl;
		TransformLog << "[";
		for(int ii = 0; ii < 4; ++ii)
		{
			TransformLog << "[";
			for(int jj = 0; jj < 4; ++jj)
			{
				TransformLog << step3->GetMatrix()->GetElement(ii, jj);
				if(jj < 3) TransformLog << ",";
			}
			TransformLog << "]";
			if(ii < 3) TransformLog << ",";
		}
		TransformLog << "]" << std::endl;
		TransformLog << std::endl;
		
		#ifdef PIPELINE_DOC
		registeredSomaPosition("neuron translation and rotation");
		double minDendDist = registeredPiaDendriteDistance();
		TransformLog << "Registered minimum Pia-dendrite distance (micron):\t" << minDendDist << std::endl;
		TransformLog << std::endl;
		#endif
		
		std::flush(std::cout << std::endl);
		std::flush(std::cout << ">>>        Checking morphology for z-scale correction" << std::endl);
		bool zScaleModified = 0;
		double oldZScale;
		std::map< vtkIdType, std::vector< double > > scaleCorrection = buildLocalZScaleCorrection(zScale, newSomaPt);
		// for now: apply z-scale correction globally to get rid of local artifacts,
		// e.g. steps in morphology
		std::map< vtkIdType, std::vector< double > >::const_iterator scaleCorrectionIt;
		for(scaleCorrectionIt = scaleCorrection.begin(); scaleCorrectionIt != scaleCorrection.end(); ++scaleCorrectionIt)
			if(scaleCorrectionIt->second[0] < zScale[0])
			{
				oldZScale = zScale[0];
				zScale[0] = scaleCorrectionIt->second[0];
				zScaleModified = 1;
			}
		
		std::cout << std::endl;
		std::cout << ">>>        Z scale factors: [" << zScale[0] << "," << zScale[1] << "," << zScale[2] << "]" << std::endl;
		std::cout << ">>>        Total column scale: " << totalScale << std::endl;
		std::cout << std::endl;
		
		inhomogeneousZScaling(zScale, newSomaPt);
		
		TransformLog << std::endl;
		TransformLog << "Old home column parameters:" << std::endl;
		TransformLog << "Pia-L4 upper\tL4 upper-L4 lower\tL4 lower-WM" << std::endl;
		TransformLog << oldPiaL4Dist << "\t" << oldL4Height << "\t" << oldL4WMDist << std::endl;
		TransformLog << "New home column parameters:" << std::endl;
		TransformLog << "Pia-L4 upper\tL4 upper-L4 lower\tL4 lower-WM" << std::endl;
		TransformLog << newPiaL4Dist << "\t" << newL4Height << "\t" << newL4WMDist << std::endl;
		TransformLog << "Morphology scale factors:" << std::endl;
		TransformLog << "Pia-L4 upper\tL4 upper-L4 lower\tL4 lower-WM\tTotal column scale" << std::endl;
		TransformLog << zScale[0] << "\t" << zScale[1] << "\t" << zScale[2] << "\t" << totalScale << std::endl;
		if(zScaleModified)
		{
			TransformLog << "WARNING: Pia-L4 upper scale was modified to prevent parts of the cell from sticking out of pia." << std::endl;
			TransformLog << "Pia-L4 upper scale befor was:\t" << oldZScale << std::endl;
		}
		TransformLog << std::endl;
		
		return step3;
	}
	else
	{
		std::cout << "Error! workingSG NULL!" << std::endl;
		TransformPointerType emptyTrans = TransformPointerType::New();
		return emptyTrans;
	}
};

/****************************************************************************/
/*returns soma-pia distance in micron                                       */
/*use whenever you want to, but result only meaningful after registration   */
/*and z scaling are finished                                                */
/****************************************************************************/
double Registration::registeredPiaDistance()
{
	if(workingSG && SBF->avgPiaSurface)
	{
		PolyDataPointerType soma = PolyDataPointerType::New();
		if(workingSG->extractLandmark(Soma, soma))
		{
			int subID;
			double somaCenter[3], pCoords[3], * weights;
			weights = new double[soma->GetCell(0)->GetNumberOfPoints()];
			soma->GetCell(0)->GetParametricCenter(pCoords);
			soma->GetCell(0)->EvaluateLocation(subID, pCoords, somaCenter, weights);
			delete [] weights;
			
			double dist = -1;
// 			double * zAxis = localZAxis(somaCenter);
			double zAxis[3];
			SBF->localZAxis(somaCenter, zAxis);
			SBF->avgPiaSurface->intersectLine(zAxis, somaCenter);
			double * piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
			if(piaIntersectPt)
			{
				dist = L2Distance3D(somaCenter, piaIntersectPt);
				delete [] piaIntersectPt;
			}
// 			delete [] zAxis;
			
// 			#ifdef PIPELINE_DOC
// 			
// 			//compute coordinates in cylindrical coordinates
// 			//measured in home column
// 			//phi measured relative to neighboring barrel along row
// 			double rPos, phiPos, zPos = dist;
// 			double somaPt[3], radialSomaPt[3];
// 			getPCenterOfStructure(workingSG, Soma, somaPt);
// 			int HBID = workingSG->getHomeBarrel();
// 			int NBID = neighborBarrel[HBID];
// 			
// 			double t, closestPt[3], colAxis[3], xAxis[3], yAxis[3];
// 			rPos = vtkLine::DistanceToLine(somaPt, avgColumns[HBID]->top, avgColumns[HBID]->bottom, t, closestPt);
// 			rPos = sqrt(rPos);
// 			
// 			for(int ii = 0; ii < 3; ++ii)
// 			{
// 				radialSomaPt[ii] = somaPt[ii] - closestPt[ii];
// 				colAxis[ii] = avgColumns[HBID]->top[ii] - avgColumns[HBID]->bottom[ii];
// 				xAxis[ii] = avgCenters[NBID][ii] - avgCenters[HBID][ii];
// 			}
// 			normalize(colAxis);
// 			double correction = vtkMath::Dot(colAxis, xAxis);
// 			for(int ii = 0; ii < 3; ++ii)
// 				xAxis[ii] -= correction*colAxis[ii];
// 			normalize(xAxis);
// 			vtkMath::Cross(colAxis, xAxis, yAxis);
// 			double xNew = vtkMath::Dot(radialSomaPt, xAxis);
// 			double yNew = vtkMath::Dot(radialSomaPt, yAxis);
// 			phiPos = std::atan2(yNew, xNew) + PI;
// 			
// 			std::string somaVarFilename(pipelineDocName);
// 			somaVarFilename += "_soma_pos.csv";
// 			std::ofstream somaVarFile;
// 			somaVarFile.open(somaVarFilename.c_str());
// 			somaVarFile << "# Soma location in global cartesian and local cylindrical coordinates" << std::endl;
// 			somaVarFile << "x y z\t" << somaPt[0] << "\t" << somaPt[1] << "\t" << somaPt[2] << std::endl;
// 			somaVarFile << "r phi z\t" << rPos << "\t" << phiPos << "\t" << zPos << std::endl;
// 			somaVarFile.close();
// 			#endif
			
			return dist;
		}
		return -1;
	}
	return -1;
};

void Registration::registeredSomaPosition(const char * label)
{
	if(workingSG && SBF->avgPiaSurface)
	{
		double somaCenter[3];
		if(somaFlag)
			getPCenterOfStructure(workingSG, Soma, somaCenter);
		
		double dist = -1;
// 		double * zAxis = localZAxis(somaCenter);
		double zAxis[3];
		SBF->localZAxis(somaCenter, zAxis);
		SBF->avgPiaSurface->intersectLine(zAxis, somaCenter);
		double * piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
		if(piaIntersectPt)
		{
			dist = L2Distance3D(somaCenter, piaIntersectPt);
			delete [] piaIntersectPt;
		}
// 		delete [] zAxis;
		
		#ifdef PIPELINE_DOC
		
		//compute coordinates in cylindrical coordinates
		//measured in home column
		//phi measured relative to neighboring barrel along row
		double rPos, phiPos, zPos = dist;
		double somaPt[3], radialSomaPt[3];
		getPCenterOfStructure(workingSG, Soma, somaPt);
		int HBID = workingSG->getHomeBarrel();
		int NBID = SBF->neighborBarrel[HBID];
		
		double t, closestPt[3], colAxis[3], xAxis[3], yAxis[3];
		rPos = vtkLine::DistanceToLine(somaPt, SBF->avgColumns[HBID]->top, SBF->avgColumns[HBID]->bottom, t, closestPt);
		rPos = sqrt(rPos);
		
		for(int ii = 0; ii < 3; ++ii)
		{
			radialSomaPt[ii] = somaPt[ii] - closestPt[ii];
			colAxis[ii] = SBF->avgColumns[HBID]->top[ii] - SBF->avgColumns[HBID]->bottom[ii];
			xAxis[ii] = SBF->avgCenters[NBID][ii] - SBF->avgCenters[HBID][ii];
		}
		normalize(colAxis);
		double correction = vtkMath::Dot(colAxis, xAxis);
		for(int ii = 0; ii < 3; ++ii)
			xAxis[ii] -= correction*colAxis[ii];
		normalize(xAxis);
		vtkMath::Cross(colAxis, xAxis, yAxis);
		double xNew = vtkMath::Dot(radialSomaPt, xAxis);
		double yNew = vtkMath::Dot(radialSomaPt, yAxis);
		phiPos = std::atan2(yNew, xNew) + PI;
		
// 		std::string somaVarFilename(pipelineDocName);
// 		somaVarFilename += label;
// 		somaVarFilename += "_soma_pos.csv";
// 		std::ofstream somaVarFile;
// 		somaVarFile.open(somaVarFilename.c_str());
// 		somaVarFile << "# Soma location in global cartesian and local cylindrical coordinates" << std::endl;
// 		somaVarFile << "x y z\t" << somaPt[0] << "\t" << somaPt[1] << "\t" << somaPt[2] << std::endl;
// 		somaVarFile << "r phi z\t" << rPos << "\t" << phiPos << "\t" << zPos << std::endl;
// 		somaVarFile.close();
		TransformLog << "Soma location step: " << label << std::endl;
		TransformLog << "x\ty\tz" << std::endl;
		TransformLog << somaPt[0] << "\t" << somaPt[1] << "\t" << somaPt[2] << std::endl;
		TransformLog << "r\tphi\tz" << std::endl;
		TransformLog << rPos << "\t" << phiPos << "\t" << zPos << std::endl;
		TransformLog << std::endl;
		#endif
	}
};

void Registration::registeredSomaPosition(const char * label, PolyDataPointerType pia)
{	if(workingSG && piaFlag)
	{
		double somaCenter[3];
		if(somaFlag)
			getPCenterOfStructure(workingSG, Soma, somaCenter);
		
		double dist = -1;
// 		double * zAxis = localZAxis(somaCenter);
		double zAxis[3];
		SBF->localZAxis(somaCenter, zAxis);
		double * piaIntersectPt = axisSurfaceIntersection(pia, zAxis, somaCenter);
		if(piaIntersectPt)
		{
			dist = L2Distance3D(somaCenter, piaIntersectPt);
			delete [] piaIntersectPt;
		}
		
		#ifdef PIPELINE_DOC
		
		//compute coordinates in cylindrical coordinates
		//measured in home column
		//phi measured relative to neighboring barrel along row
		double rPos, phiPos, zPos = dist;
		double somaPt[3], radialSomaPt[3];
		getPCenterOfStructure(workingSG, Soma, somaPt);
		int HBID = workingSG->getHomeBarrel();
		int NBID = SBF->neighborBarrel[HBID];
		
		double t, closestPt[3], colAxis[3], xAxis[3], yAxis[3];
		rPos = vtkLine::DistanceToLine(somaPt, SBF->avgColumns[HBID]->top, SBF->avgColumns[HBID]->bottom, t, closestPt);
		rPos = sqrt(rPos);
		
		for(int ii = 0; ii < 3; ++ii)
		{
			radialSomaPt[ii] = somaPt[ii] - closestPt[ii];
			colAxis[ii] = SBF->avgColumns[HBID]->top[ii] - SBF->avgColumns[HBID]->bottom[ii];
			xAxis[ii] = SBF->avgCenters[NBID][ii] - SBF->avgCenters[HBID][ii];
		}
		normalize(colAxis);
		double correction = vtkMath::Dot(colAxis, xAxis);
		for(int ii = 0; ii < 3; ++ii)
			xAxis[ii] -= correction*colAxis[ii];
		normalize(xAxis);
		vtkMath::Cross(colAxis, xAxis, yAxis);
		double xNew = vtkMath::Dot(radialSomaPt, xAxis);
		double yNew = vtkMath::Dot(radialSomaPt, yAxis);
		phiPos = std::atan2(yNew, xNew) + PI;
		
// 		std::string somaVarFilename(pipelineDocName);
// 		somaVarFilename += label;
// 		somaVarFilename += "_soma_pos.csv";
// 		std::ofstream somaVarFile;
// 		somaVarFile.open(somaVarFilename.c_str());
// 		somaVarFile << "# Soma location in global cartesian and local cylindrical coordinates" << std::endl;
// 		somaVarFile << "x y z\t" << somaPt[0] << "\t" << somaPt[1] << "\t" << somaPt[2] << std::endl;
// 		somaVarFile << "r phi z\t" << rPos << "\t" << phiPos << "\t" << zPos << std::endl;
// 		somaVarFile.close();
		TransformLog << "Soma location step: " << label << std::endl;
		TransformLog << "x\ty\tz" << std::endl;
		TransformLog << somaPt[0] << "\t" << somaPt[1] << "\t" << somaPt[2] << std::endl;
		TransformLog << "r\tphi\tz" << std::endl;
		TransformLog << rPos << "\t" << phiPos << "\t" << zPos << std::endl;
		TransformLog << std::endl;
		#endif
	}
};

/****************************************************************************/
/*returns soma-column axis distance in micron                               */
/*use whenever you want to, but result only meaningful after registration   */
/*and z scaling are finished                                                */
/****************************************************************************/
double Registration::somaDistanceToColumnAxis()
{
	if(workingSG)
	{
		int HBID = workingSG->getHomeBarrel();
		double projSomaPt[3], somaPt[3], t;
		getPCenterOfStructure(workingSG, Soma, somaPt);
		double radialDist = vtkLine::DistanceToLine(somaPt, SBF->avgColumns[HBID]->top, SBF->avgColumns[HBID]->bottom, t, projSomaPt);
		return sqrt(radialDist);
	}
	return -1;
};

/****************************************************************************/
/*returns closest dendrite-pia distance in micron                           */
/****************************************************************************/
double Registration::registeredPiaDendriteDistance()
{
	if(workingSG && SBF->avgPiaSurface)
	{
		double minDist = 1E8;
		std::vector< Edge * >::const_iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label != Dendrite && (*edgeIt)->label != ApicalDendrite && (*edgeIt)->label != BasalDendrite)
				continue;
			
			std::list< double * >::const_iterator edgePtIt;
			for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
			{
				double * dendPt = *edgePtIt;
// 				double * zAxis = localZAxis(dendPt);
				double zAxis[3];
				SBF->localZAxis(dendPt, zAxis);
				SBF->avgPiaSurface->intersectLine(zAxis, dendPt);
				double * piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
				if(piaIntersectPt)
				{
					double tmpDist = L2Distance3D(dendPt, piaIntersectPt);
					if(tmpDist < minDist)
						minDist = tmpDist;
					delete [] piaIntersectPt;
				}
// 				delete [] zAxis;
			}
		}
		// in case Dendrite really sticks out of Pia
		if(minDist < 5)
		{
			double maxDendPt[3], zAxis[3];
			maxDendPoint(maxDendPt);
			SBF->localZAxis(maxDendPt, zAxis);
			SBF->avgPiaSurface->intersectLine(zAxis, maxDendPt);
			double * piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
			if(piaIntersectPt)
			{
				double tmpDist = L2Distance3D(maxDendPt, piaIntersectPt);
				if(maxDendPt[2] > piaIntersectPt[2])
					minDist = -tmpDist;
				else
					minDist = tmpDist;
				delete [] piaIntersectPt;
			}
		}
		return minDist;
	}
	return -1;
};

double Registration::registeredPiaDendriteDistance(PolyDataPointerType pia)
{
	if(workingSG && piaFlag)
	{
		CellLocatorPointerType locator = CellLocatorPointerType::New();
		locator->AutomaticOn();
		locator->SetDataSet(pia);
		locator->BuildLocator();
		double minDist = 1E8;
		std::vector< Edge * >::const_iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label != Dendrite && (*edgeIt)->label != ApicalDendrite && (*edgeIt)->label != BasalDendrite)
				continue;
			
			std::list< double * >::const_iterator edgePtIt;
			for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
			{
				double * dendPt = *edgePtIt;
// 				double * zAxis = localZAxis(dendPt);
				double zAxis[3];
				SBF->localZAxis(dendPt, zAxis);
				double * piaIntersectPt = axisSurfaceIntersection(locator, zAxis, dendPt);
				if(piaIntersectPt)
				{
					double tmpDist = L2Distance3D(dendPt, piaIntersectPt);
					if(tmpDist < minDist)
						minDist = tmpDist;
					delete [] piaIntersectPt;
				}
// 				delete [] zAxis;
			}
		}
		// in case Dendrite really sticks out of Pia
		if(minDist < 5)
		{
			double maxDendPt[3], zAxis[3];
			maxDendPoint(maxDendPt);
			SBF->localZAxis(maxDendPt, zAxis);
			double * piaIntersectPt = axisSurfaceIntersection(locator, zAxis, maxDendPt);
			if(piaIntersectPt)
			{
				double tmpDist = L2Distance3D(maxDendPt, piaIntersectPt);
				if(maxDendPt[2] > piaIntersectPt[2])
					minDist = -tmpDist;
				else
					minDist = tmpDist;
				delete [] piaIntersectPt;
			}
		}
		return minDist;
	}
	return -1;
};

double Registration::unregisteredPiaDendriteDistance(PolyDataPointerType pia)
{
	if(workingSG && piaFlag)
	{
		CellLocatorPointerType locator = CellLocatorPointerType::New();
		locator->AutomaticOn();
		locator->SetDataSet(pia);
		locator->BuildLocator();
		double minDist = 1E8;
		std::vector< Edge * >::const_iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label != Dendrite && (*edgeIt)->label != ApicalDendrite && (*edgeIt)->label != BasalDendrite)
				continue;
			
			std::list< double * >::const_iterator edgePtIt;
			for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
			{
				double * dendPt = *edgePtIt;
// 				double * zAxis = localZAxis(dendPt);
				double zAxis[] = {0,0,zReversed ? -1 : 1};
				double * piaIntersectPt = axisSurfaceIntersection(locator, zAxis, dendPt);
				if(piaIntersectPt)
				{
					double tmpDist = L2Distance3D(dendPt, piaIntersectPt);
					if(tmpDist < minDist)
						minDist = tmpDist;
					delete [] piaIntersectPt;
				}
// 				delete [] zAxis;
			}
		}
		return minDist;
	}
	return -1;
};

double Registration::registeredPiaSegmentDistance(int label, int edgeNr)
{
	if(workingSG && SBF->avgPiaSurface)
	{
		double minDist = 1E09;
		unsigned int labelCnt = 1;
		std::vector< Edge * >::const_iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label != label)
				continue;
			
			if(labelCnt != edgeNr)
			{
				++labelCnt;
				continue;
			}
			++labelCnt;
			
			// debug
			(*edgeIt)->label = Landmark;
			
			std::list< double * >::const_iterator edgePtIt;
			for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
			{
// 				double * pt = *edgePtIt, * tmpZAxis = localZAxis(pt), intersectPt[3];
				double * pt = *edgePtIt, tmpZAxis[3], intersectPt[3];
				SBF->localZAxis(pt, tmpZAxis);
				SBF->avgPiaSurface->intersectLine(tmpZAxis, pt);
				if(SBF->avgPiaSurface->isValid())
				{
					SBF->avgPiaSurface->getLastIntersectPoint(intersectPt);
					double tmpDist = L2Distance3D(pt, intersectPt);
					if(tmpDist < minDist)
						minDist = tmpDist;
				}
// 				delete [] tmpZAxis;
			}
		}
		return minDist;
	}
	return -1;
};

double Registration::globalPiaSegmentDistance(int label, int edgeNr, PolyDataPointerType pia)
{
	if(workingSG && piaFlag)
	{
		double minDist = 1E09;
		unsigned int labelCnt = 1;
		std::vector< Edge * >::const_iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label != label)
				continue;
			
			if(labelCnt != edgeNr)
			{
				++labelCnt;
				continue;
			}
			++labelCnt;
			
			std::list< double * >::const_iterator edgePtIt;
			for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
			{
				double * pt = *edgePtIt, * intersectPt, tmpZAxis[] = {0,0,zReversed ? -1 : 1};
				intersectPt = axisSurfaceIntersection(pia, tmpZAxis, pt);
				if(intersectPt)
				{
					double tmpDist = L2Distance3D(pt, intersectPt);
					if(tmpDist < minDist)
						minDist = tmpDist;
					delete [] intersectPt;
				}
			}
		}
		return minDist;
	}
	return -1;
};

/****************************************************************************/
/*determine column/septum based on registered soma position                 */
/*WARNING: assumes that avg home barrel contour has been calculated!!!      */
/****************************************************************************/
int Registration::getColumnSeptumFlag()
{
	if(workingSG)
	{
		int HBID = workingSG->getHomeBarrel();
		if(avgBarrelContours.find(HBID) != avgBarrelContours.end())
		{
			double somaPt[3], zAxis[3];
			getPCenterOfStructure(workingSG, Soma, somaPt);
			SBF->localZAxis(HBID, zAxis);
			
			double tol = 1E-3, t, intersectPt[3], pCoords[3], pt1[3], pt2[3];
			int subID;
			for(int ii = 0; ii < 3; ++ii)
			{
				pt1[ii] = somaPt[ii] + 4000*zAxis[ii];
				pt2[ii] = somaPt[ii] - 4000*zAxis[ii];
			}
			if(avgBarrelContours[HBID]->GetCell(0)->IntersectWithLine(pt1, pt2, tol, t, intersectPt, pCoords, subID))
			{
				return 1;
			}
			return 0;
		}
	}
	return -1;
};

/******************************************************************************/
/*get laminar position (spura/gran/infra) of registered neuron                */
/******************************************************************************/
int Registration::getLaminarPosition()
{
	if(somaFlag)
	{
		int HBID = workingSG->getHomeBarrel();
		double somaPt[3];
		getPCenterOfStructure(workingSG, Soma, somaPt);
		
		return SBF->laminarPosition(somaPt);
//		if(somaPt[2] > SBF->avgBarrels[HBID]->top[2])
//			return 0;
//		else if(somaPt[2] > SBF->avgBarrels[HBID]->bottom[2])
//			return 1;
//		else
//			return 2;
	}
	return 0;
//	return -1;
}

/******************************************************************************/
/*simple avg contour by sampling barrel contours at regular angular intervals */
/*and then averaging in z                                                     */
/******************************************************************************/
void Registration::computeAverageHomeBarrel(PolyDataPointerType completeBarrel, double barrelCentroid[3], double barrelAxis[3], int ID)
{
	int nrAngles = 36;
	std::vector< std::vector< double * > > thisPlaneSamplingPoints;	// one vector for each sampling angle
	for(int ii = 0; ii < nrAngles; ++ii)
	{
		std::vector< double * > angleSamplingVec;
		thisPlaneSamplingPoints.push_back(angleSamplingVec);
	}
// 	std::cout << std::endl;
// 	std::cout << "nr. of home barrel contours: " << completeBarrel->GetNumberOfCells() << std::endl;
	for(int ii = 0; ii < completeBarrel->GetNumberOfCells(); ++ii)
	{
		double * centerPoint = new double[3];
		double paramCenter[3];
		
		//parametric center
		int subID;
		double pCoords[3], * weights;
		weights = new double[completeBarrel->GetCell(ii)->GetNumberOfPoints()];
		completeBarrel->GetCell(ii)->GetParametricCenter(pCoords);
		completeBarrel->GetCell(ii)->EvaluateLocation(subID, pCoords, centerPoint, weights);
		
		// identify polygon normals so each contour can be traversed in same direction
		// this method is only called AFTER landmark registration
		// -> already aligned to D2 coordinate system -> ignore zReversed
		double normal[3], rotAxis[3], zUnitVec[3] = {0,0,1};
// 		double normal[3], rotAxis[3], zUnitVec[3] = {0,0,zReversed ? -1 : 1};
		double rotAngle = 0;
		vtkPolygon::ComputeNormal(completeBarrel->GetCell(ii)->GetPoints(), normal);
		bool reverseDirection;
		if(normal[2] < 0)
		{
			reverseDirection = 1;
			zUnitVec[2] *= -1;
		}
		else
			reverseDirection = 0;
		
		// rotate sampling plane so it is perpendicular to 
		// the plane defined by the contour
		vtkMath::Cross(normal, zUnitVec, rotAxis);
		normalize(rotAxis);
		rotAngle = std::acos(vtkMath::Dot(normal, zUnitVec))*180/PI;
		TransformPointerType rot = TransformPointerType::New();
		rot->RotateWXYZ(rotAngle, rotAxis);
		
// 		std::cout << std::endl;
// 		std::cout << "home barrel cell " << ii << std::endl;
// 		std::cout << "normal = [" << normal[0] << "," << normal[1] << "," << normal[2] << "]" << std::endl;
// 		std::cout << "rotAxis = [" << rotAxis[0] << "," << rotAxis[1] << "," << rotAxis[2] << "]" << std::endl;
// 		std::cout << "zUnitVec = [" << zUnitVec[0] << "," << zUnitVec[1] << "," << zUnitVec[2] << "]" << std::endl;
// 		std::cout << "rotAngle = " << rotAngle << std::endl;
		
		PlanePointerType thisPlane = PlanePointerType::New();
		PointsPointerType thisCellPoints = completeBarrel->GetCell(ii)->GetPoints();
		for(int jj = 0; jj < nrAngles; ++jj)	// has to be unique b/c order of the points is not clear
		{
// 			std::cout << std::endl;
// 			std::cout << "sampling angle " << jj << std::endl;
			
			double angle = jj*2*PI/nrAngles;
			double planeNormal[3], hPlaneNormal[4];
			if(reverseDirection)
			{
				hPlaneNormal[0] = -sin(angle);
				hPlaneNormal[1] = -cos(angle);
			}
			else
			{
				hPlaneNormal[0] = sin(angle);
				hPlaneNormal[1] = cos(angle);
			}
			hPlaneNormal[2] = 0, hPlaneNormal[3] = 1;
			rot->MultiplyPoint(hPlaneNormal, hPlaneNormal);
			planeNormal[0] = hPlaneNormal[0];
			planeNormal[1] = hPlaneNormal[1];
			planeNormal[2] = hPlaneNormal[2];
			
			thisPlane->SetOrigin(centerPoint);
			thisPlane->SetNormal(planeNormal);
			
			double * firstPt = new double[3];
			thisCellPoints->GetPoint(0, firstPt);
			double lastVal = thisPlane->EvaluateFunction(firstPt);
			delete [] firstPt;
			for(int kk = 1; kk <= thisCellPoints->GetNumberOfPoints(); ++kk)
			{
// 				std::cout << "cell " << ii << " pt " << kk << std::endl;
				
				double * thisPt = new double[3];
				thisCellPoints->GetPoint(kk%thisCellPoints->GetNumberOfPoints(), thisPt);
				double ptVal = thisPlane->EvaluateFunction(thisPt);
				if((ptVal < 0 && lastVal > 0))	// has to be unique b/c order of the points is not clear
				{
					// regular case
					double * samplePt =  new double[3];
					double * lastPt = new double[3];
					thisCellPoints->GetPoint(kk-1, lastPt);
					double direction[3];
					double normD = 0;
					for(int ll = 0; ll < 3; ++ll)
					{
						direction[ll] = lastPt[ll] - thisPt[ll];
						normD += direction[ll]*direction[ll];
					}
					normD = sqrt(normD);
					if(normD)
					{
						double distAlongLine = std::abs(ptVal);
						double correctionAngle = std::abs(direction[0]*planeNormal[0]/normD + direction[1]*planeNormal[1]/normD);
						if(correctionAngle)
							distAlongLine /= correctionAngle;
						for(int ll = 0; ll < 3; ++ll)
							samplePt[ll] = thisPt[ll] + direction[ll]*distAlongLine/normD;
					}
					else
						for(int ll = 0; ll < 3; ++ll)
							samplePt[ll] = thisPt[ll];
						thisPlaneSamplingPoints[jj].push_back(samplePt);
					delete [] lastPt;
// 					std::cout << "regular order" << std::endl;
// 					std::cout << "pt " << kk << " @  [" << samplePt[0] << "," << samplePt[1] << "," << samplePt[2] << "]" << std::endl;
				}
				else if(ptVal == 0 && lastVal > 0)
				{
					// case: this point is directly on plane
					double * samplePt =  new double[3];
					for(int ll = 0; ll < 3; ++ll)
						samplePt[ll] = thisPt[ll];
					thisPlaneSamplingPoints[jj].push_back(samplePt);
// 					std::cout << "pt on plane" << std::endl;
// 					std::cout << "pt " << kk << " @  [" << samplePt[0] << "," << samplePt[1] << "," << samplePt[2] << "]" << std::endl;
				}
// 				else if(ptVal != 0 && lastVal == 0)
// 				{
// 					;// do nothing
// 				}
// 				else if(ptVal == 0 && lastVal == 0)
// 				{
// 					;// either an error or both points lie directly on the plane (extremely unlikely)
// 				}
				lastVal = ptVal;
				delete [] thisPt;
			} // all cell points
		} // all angles
	} // all barrel cells
	
	PolyDataPointerType avgContour = PolyDataPointerType::New();
	PointsPointerType contourPts = PointsPointerType::New();
	PolygonPointerType contourPoly = PolygonPointerType::New();
	avgContour->Allocate(1);
	contourPts->SetDataTypeToFloat();
	contourPts->SetNumberOfPoints(nrAngles);
	contourPoly->GetPointIds()->SetNumberOfIds(nrAngles);
	
	for(int ii = 0; ii < nrAngles; ++ii)
	{
		if(thisPlaneSamplingPoints[ii].size())
		{
// 			std::cout << std::endl;
// 			std::cout << "sampling angle " << ii << std::endl;
			
			double * avgPt = new double[3];
			avgPt[0] = 0, avgPt[1] = 0, avgPt[2] = 0;
			for(int jj = 0; jj < thisPlaneSamplingPoints[ii].size(); ++jj)
			{
				double * tmpPt = thisPlaneSamplingPoints[ii][jj];
				for(int kk = 0; kk < 3; ++kk)
					avgPt[kk] += tmpPt[kk];
				
// 				std::cout << "contour " << jj << std::endl;
// 				std::cout << "pt = [" << tmpPt[0] << "," << tmpPt[1] << "," << tmpPt[2] << "]" << std::endl;
			}
			for(int jj = 0; jj < 3; ++jj)
				avgPt[jj] = avgPt[jj]/double(thisPlaneSamplingPoints[ii].size());
			
			contourPts->InsertPoint(ii, avgPt);
			contourPoly->GetPointIds()->SetId(ii, ii);
		}
	}
	avgContour->InsertNextCell(contourPoly->GetCellType(), contourPoly->GetPointIds());
	avgContour->SetPoints(contourPts);
	avgContour->Update();
	
	avgBarrelContours.insert(std::pair< int, PolyDataPointerType >(ID, avgContour));
	
	// debug only:
	if(workingSG)
		workingSG->addPolyDataObject(avgContour, Barrel);
};

std::vector< double > Registration::computeManualBarrelParameters(PolyDataPointerType barrel, PolyDataPointerType pia, PolyDataPointerType wm, double* newAxis, double* barrelCenter, std::vector< double* > endPoints, int label, std::map< int, PolyDataPointerType >& avgBarrels)
{
	std::vector< double > params;
	if(barrel->GetNumberOfPoints())
	{
		double height, topPiaDist, bottomPiaDist, piaWMDist = 0, area;
		double * topPt = new double[3];
		double * bottomPt = new double[3];
		topPt[0] = endPoints[0][0], topPt[1] = endPoints[0][1], topPt[2] = endPoints[0][2];
		bottomPt[0] = endPoints[1][0], bottomPt[1] = endPoints[1][1], bottomPt[2] = endPoints[1][2];
		height = sqrt((topPt[0] - bottomPt[0])*(topPt[0] - bottomPt[0]) + (topPt[1] - bottomPt[1])*(topPt[1] - bottomPt[1]) + (topPt[2] - bottomPt[2])*(topPt[2] - bottomPt[2]));
		
		double * piaIntersection = axisSurfaceIntersection(pia, newAxis, barrelCenter);
		bottomPiaDist = sqrt((piaIntersection[0] - bottomPt[0])*(piaIntersection[0] - bottomPt[0]) + (piaIntersection[1] - bottomPt[1])*(piaIntersection[1] - bottomPt[1]) + (piaIntersection[2] - bottomPt[2])*(piaIntersection[2] - bottomPt[2]));
		topPiaDist = sqrt((piaIntersection[0] - topPt[0])*(piaIntersection[0] - topPt[0]) + (piaIntersection[1] - topPt[1])*(piaIntersection[1] - topPt[1]) + (piaIntersection[2] - topPt[2])*(piaIntersection[2] - topPt[2]));
		
		double * wmIntersection;
		if(wmFlag)
		{
			wmIntersection = axisSurfaceIntersection(wm, newAxis, barrelCenter);
			if(wmIntersection)
				piaWMDist = sqrt((piaIntersection[0] - wmIntersection[0])*(piaIntersection[0] - wmIntersection[0]) + (piaIntersection[1] - wmIntersection[1])*(piaIntersection[1] - wmIntersection[1]) + (piaIntersection[2] - wmIntersection[2])*(piaIntersection[2] - wmIntersection[2]));
		}
		PolyDataPointerType avgContours = avgBarrels[label];
		double normal[3];
		// don't use GetCell(0)->GetPoints b/c points in PointIDList are not sorted in strictly ascending order for overlap-corrected contours!
		area = vtkPolygon::ComputeArea(avgContours->GetPoints(), avgContours->GetCell(0)->GetNumberOfPoints(), avgContours->GetCell(0)->GetPointIds()->GetPointer(0), normal);
		
		params.push_back(bottomPiaDist);
		params.push_back(topPiaDist);
		params.push_back(height);
		params.push_back(piaWMDist);
		params.push_back(area);
		
		delete [] topPt, delete [] bottomPt, delete [] piaIntersection;
		if(wm) delete [] wmIntersection;
		return params;
	}
	else
	{
		std::cout << "Error! PolyData barrel is empty! Could not calculate barrel parameters." << std::endl;
		return params;
	}
};

/******************************************************************************/
/*align home barrel endpoints by translation and rotation                     */
/*barrels[0] : input; barrels[1] : standard                                   */
/******************************************************************************/
HomogeneousMatrixPointerType Registration::alignHomeBarrel(std::map< int, Column * > * barrels)
{
	int HBID = workingSG->getHomeBarrel();
	double shift[3], invShift[3], inputAxis[3], standardAxis[3];
	for(int ii = 0; ii < 3; ++ii)
	{
		shift[ii] = -0.5*(barrels[0][HBID]->top[ii] + barrels[0][HBID]->bottom[ii]);
		invShift[ii] = SBF->avgCenters[workingSG->getHomeBarrel()][ii];
		inputAxis[ii] = barrels[0][HBID]->top[ii] - barrels[0][HBID]->bottom[ii];
		standardAxis[ii] = SBF->avgAxes[HBID][ii];
	}
	normalize(inputAxis), normalize(standardAxis);
	
	TransformPointerType trans = TransformPointerType::New();
	TransformPointerType rot = TransformPointerType::New();
	TransformPointerType invTrans = TransformPointerType::New();
	trans->Translate(shift);
	rot->SetMatrix(transformToBarrelCoordinates(inputAxis, standardAxis));
	invTrans->Translate(invShift);
	
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrels[0].find(ID) != barrels[0].end())
		{
			barrels[0][ID]->translateColumn(shift);
			barrels[0][ID]->rotateColumn(rot->GetMatrix());
			barrels[0][ID]->translateColumn(invShift);
		}
	}
	rot->Concatenate(trans);
	invTrans->Concatenate(rot);
	
	return invTrans->GetMatrix();
};

/******************************************************************************/
/*align home barrel endpoints only by translation                             */
/*barrels[0] : input; barrels[1] : standard                                   */
/******************************************************************************/
HomogeneousMatrixPointerType Registration::alignHomeBarrel2(std::map< int, Column * > * barrels)
{
	int HBID = workingSG->getHomeBarrel();
	double shift[3];
	for(int ii = 0; ii < 3; ++ii)
	{
		shift[ii] = SBF->avgCenters[HBID][ii] - 0.5*(barrels[0][HBID]->top[ii] + barrels[0][HBID]->bottom[ii]);
	}
	
	TransformPointerType trans = TransformPointerType::New();
	trans->Translate(shift);
	
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrels[0].find(ID) != barrels[0].end())
			barrels[0][ID]->translateColumn(shift);
	}
	
	return trans->GetMatrix();
};

/******************************************************************************/
/*align remaining landmarks by rotation in plane perpendicular to the aligned */
/*home barrel axis                                                            */
/*rotation angle is average signed offset angle between top/bottom of input   */
/*and reference barrels (except for home barrel, which is already aligned)    */
/*barrels[0] : input; barrels[1] : standard                                   */
/******************************************************************************/
HomogeneousMatrixPointerType Registration::alignRemainingBarrels(std::map< int, Column * > * barrels)
{
	int HBID = workingSG->getHomeBarrel();
	PlanePointerType topPlane = PlanePointerType::New();
	PlanePointerType bottomPlane = PlanePointerType::New();
	topPlane->SetNormal(SBF->avgAxes[HBID]);
	topPlane->SetOrigin(SBF->avgBarrels[HBID]->top);
	bottomPlane->SetNormal(SBF->avgAxes[HBID]);
	bottomPlane->SetOrigin(SBF->avgBarrels[HBID]->bottom);
	double avgOffsetAngle = 0;
	double avgAngleCnt = 0;
	
	#ifdef REG_ACCURACY
	std::vector< double > angleOffsetVec;
	#endif
	
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrels[0].find(ID) != barrels[0].end() && barrels[1].find(ID) != barrels[1].end() && ID != HBID)
		{
			double topProj1[4], bottomProj1[4];
			double topProj2[4], bottomProj2[4];
			double topDist1 = topPlane->EvaluateFunction(barrels[0][ID]->top), bottomDist1 = bottomPlane->EvaluateFunction(barrels[0][ID]->bottom);
			double topDist2 = topPlane->EvaluateFunction(barrels[1][ID]->top), bottomDist2 = bottomPlane->EvaluateFunction(barrels[1][ID]->bottom);
			for(int ii = 0; ii < 3; ++ii)
			{
				topProj1[ii] = barrels[0][ID]->top[ii] - topDist1*SBF->avgAxes[HBID][ii] - SBF->avgBarrels[HBID]->top[ii];
				bottomProj1[ii] = barrels[0][ID]->bottom[ii] - bottomDist1*SBF->avgAxes[HBID][ii] - SBF->avgBarrels[HBID]->bottom[ii];
				topProj2[ii] = barrels[1][ID]->top[ii] - topDist2*SBF->avgAxes[HBID][ii] - SBF->avgBarrels[HBID]->top[ii];
				bottomProj2[ii] = barrels[1][ID]->bottom[ii] - bottomDist2*SBF->avgAxes[HBID][ii] - SBF->avgBarrels[HBID]->bottom[ii];
			}
			topProj1[3] = bottomProj1[3] = 1;
			topProj2[3] = bottomProj2[3] = 1;
			double tmpTopRotAngle = atan2(topProj2[1], topProj2[0])*180/PI;
			double tmpBottomRotAngle = atan2(bottomProj2[1], bottomProj2[0])*180/PI;
			TransformPointerType tmpTopRot = TransformPointerType::New();
			TransformPointerType tmpBottomRot = TransformPointerType::New();
			tmpTopRot->RotateWXYZ(-tmpTopRotAngle, 0, 0, 1);
			tmpBottomRot->RotateWXYZ(-tmpBottomRotAngle, 0, 0, 1);
			tmpTopRot->MultiplyPoint(topProj1, topProj1);
			tmpBottomRot->MultiplyPoint(bottomProj1, bottomProj1);
			double tmpTopDiff = atan2(topProj1[1], topProj1[0]);
			double tmpBottomDiff = atan2(bottomProj1[1], bottomProj1[0]);
			avgOffsetAngle += tmpTopDiff;
			avgOffsetAngle += tmpBottomDiff;
			avgAngleCnt += 2;
			#ifdef REG_ACCURACY
			angleOffsetVec.push_back(tmpTopDiff);
			angleOffsetVec.push_back(tmpBottomDiff);
			#endif
		}
	}
	if(avgAngleCnt)
		avgOffsetAngle = avgOffsetAngle*180/(PI*avgAngleCnt);
	
	#ifdef REG_ACCURACY
	double angleStd = 0;
	for(int ii = 0; ii < angleOffsetVec.size(); ++ii)
		angleStd += (angleOffsetVec[ii] - avgOffsetAngle*PI/180)*(angleOffsetVec[ii] - avgOffsetAngle*PI/180);
	if(angleOffsetVec.size() > 1)
		angleStd = sqrt(angleStd/(angleOffsetVec.size()-1))*180/PI;
	else
		angleStd = avgOffsetAngle;
	angleUncertainty = angleStd;
	#endif
	
	double shift[3], inverseShift[3];
	for(int ii = 0; ii < 3; ++ii)
	{
		shift[ii] = -SBF->avgCenters[HBID][ii];
		inverseShift[ii] = SBF->avgCenters[HBID][ii];
	}
	
	TransformPointerType step1 = TransformPointerType::New();
	TransformPointerType step2 = TransformPointerType::New();
	TransformPointerType step3 = TransformPointerType::New();
	step1->Translate(shift);
	step2->RotateWXYZ(-avgOffsetAngle, SBF->avgAxes[HBID]);
	step3->Translate(inverseShift);
	
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrels[0].find(ID) != barrels[0].end())
		{
			barrels[0][ID]->translateColumn(shift);
			barrels[0][ID]->rotateColumn(step2->GetMatrix());
			barrels[0][ID]->translateColumn(inverseShift);
		}
	}
	step2->Concatenate(step1);
	step3->Concatenate(step2);
	
	return step3->GetMatrix();
};

/******************************************************************************/
/*re-calculate home barrel reconstruction w.r.t. local z axis                 */
/*important to avoid systematic error when  calculating                       */
/*correct z-scaling factors!!!                                                */
/******************************************************************************/
void Registration::alignHomeBarrelAxis(std::map< int, Column * > * barrels, int HBID)
{
	// one-time barrel reconstruction
// 	double * HBAxis = localZAxis(HBID);
// 	double oldAxis[3], HBCenter[3], oldHeight = barrels[0][HBID]->getHeight();
// 	for(int ii = 0; ii < 3; ++ii)
// 		oldAxis[ii] = barrels[0][HBID]->top[ii] - barrels[0][HBID]->bottom[ii];
// 	normalize(oldAxis);
// 	for(int ii = 0; ii < 3; ++ii)
// 		HBCenter[ii] = barrels[0][HBID]->bottom[ii] + 0.5*oldHeight*oldAxis[ii];
// 	
// 	double newTop[3], newBottom[3], topLinePt[3], bottomLinePt[3], t;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		topLinePt[ii] = HBCenter[ii] + 1000*HBAxis[ii];
// 		bottomLinePt[ii] = HBCenter[ii] - 1000*HBAxis[ii];
// 	}
// 	delete [] HBAxis;
// 	
// 	vtkLine::DistanceToLine(barrels[0][HBID]->top, topLinePt, bottomLinePt, t, newTop);
// 	vtkLine::DistanceToLine(barrels[0][HBID]->bottom, topLinePt, bottomLinePt, t, newBottom);
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		barrels[0][HBID]->top[ii] = newTop[ii];
// 		barrels[0][HBID]->bottom[ii] = newBottom[ii];
// 	}
	
		
	//now, add the barrel contours to SpatialGraph for nice visualization
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrels[0].find(ID) != barrels[0].end())
		{
			double HBAxis[3];
			double barrelCentroid[3];
			barrels[0][ID]->getCenter(barrelCentroid);
			// reconstruction of home barrel along standard axis
			if(ID == HBID)
				SBF->localZAxis(HBID, HBAxis);
			else
			{
				for(int ii = 0; ii < 3; ++ii)
					HBAxis[ii] = barrels[0][ID]->top[ii] - barrels[0][ID]->bottom[ii];
				vtkMath::Normalize(HBAxis);
			}
			
			PolyDataPointerType currentBarrel = PolyDataPointerType::New();
			if(workingSG->extractLandmark(ID, currentBarrel))
			{
				// compute avg barrel contour to determine
				// column/septum cell later
				computeAverageHomeBarrel(currentBarrel, barrelCentroid, HBAxis, ID);
				
				// recompute HB top and bottom points
				// along standard axis
				if(ID == HBID)
				{
					std::vector< double * > endPts;
					bool HBRecon = 1;
					closeBarrelAlongNewAxis(HBAxis, barrelCentroid, currentBarrel, endPts, HBRecon);
					barrels[0][HBID]->top = endPts[0], barrels[0][HBID]->bottom = endPts[1];
					barrels[0][HBID]->contours = currentBarrel;
				}
			}
			
			// now, add contours to SpatialGraph for visualization
			workingSG->addPolyDataObject(barrels[0][ID]->contours, ID);
		}
	}
};

/******************************************************************************/
/*optional correction for systematic offset of all barrels along their local  */
/*registration axis                                                           */
/******************************************************************************/
void Registration::localLandmarkOffset(std::map< int, Column * > * barrels)
{
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrels->find(ID) != barrels->end() && contourCorrection.find(ID) != contourCorrection.end())
		{
			
			double zAxis[3], barrelCentroid[3];
			for(int ii = 0; ii < 3; ++ii)
				barrelCentroid[ii] = 0.5*(barrels[0][ID]->top[ii] + barrels[0][ID]->bottom[ii]);
			SBF->localZAxis(barrelCentroid, zAxis);
			
			// this part for absolute mode
			double dSupra = contourCorrection[ID][0];
			double dL4 = contourCorrection[ID][1];
			double dInfra = contourCorrection[ID][2];
			
			// shift barrel top/bottom
			for(int ii = 0; ii < 3; ++ii)
			{
				barrels[0][ID]->top[ii] += 0.5*dL4*zAxis[ii];
				barrels[0][ID]->bottom[ii] -= 0.5*dL4*zAxis[ii];
			}
		}
	}
};

/******************************************************************************/
/*optional correction for systematic offset in artificial Pia                 */
/******************************************************************************/
void Registration::correctArtificialPia(PolyDataPointerType pia)
{
	if(workingSG)
	{
		int HBID = workingSG->getHomeBarrel();
		if(!HBID)
		{
			std::cout << "Error! No home barrel defined. Cannot correct manual contours." << std::endl;
			return;
		}
// 		double somaPt[3], * zAxis, piaShift[3];
		double somaPt[3], zAxis[3], piaShift[3];
		getPCenterOfStructure(workingSG, Soma, somaPt);
// 		zAxis = localZAxis(somaPt);
		SBF->localZAxis(somaPt, zAxis);
		for(int ii = 0; ii < 3; ++ii)
			piaShift[ii] = correctPiaDist*zAxis[ii];
// 		delete [] zAxis;
		
		TransformPointerType piaTrans = TransformPointerType::New();
		piaTrans->Translate(piaShift);
		if(piaFlag)
		{
			workingSG->setTransformation(piaTrans);
			workingSG->applyTransformation(Pia);
			TransformFilterType piaTransFilter = TransformFilterType::New();
			piaTransFilter->SetTransform(piaTrans);
			piaTransFilter->SetInput(pia);
			piaTransFilter->Update();
			pia->DeepCopy(piaTransFilter->GetOutput());
			#ifdef PIPELINE_DOC
			globalPia->DeepCopy(pia);
			#endif
		}
	}
};

/******************************************************************************/
/*optional correction for systematic offset in artificial Pia                 */
/******************************************************************************/
void Registration::constantPiaOffset(PolyDataPointerType pia, double offset)
{
	if(workingSG)
	{
		double zAxis[] = {0,0,zReversed ? -1 : 1}, piaShift[3];
		for(int ii = 0; ii < 3; ++ii)
			piaShift[ii] = offset*zAxis[ii];
		
		TransformPointerType piaTrans = TransformPointerType::New();
		piaTrans->Translate(piaShift);
		if(piaFlag)
		{
			workingSG->setTransformation(piaTrans);
			workingSG->applyTransformation(Pia);
			TransformFilterType piaTransFilter = TransformFilterType::New();
			piaTransFilter->SetTransform(piaTrans);
			piaTransFilter->SetInput(pia);
			piaTransFilter->Update();
			pia->DeepCopy(piaTransFilter->GetOutput());
		}
	}
};

void Registration::constantLandmarkOffset(PolyDataPointerType pia, PolyDataPointerType wm, double offset)
{
	if(workingSG)
	{
		double zAxis[] = {0,0,zReversed ? -1 : 1}, shiftVec[3];
		for(int ii = 0; ii < 3; ++ii)
			shiftVec[ii] = offset*zAxis[ii];
		
		TransformPointerType translation = TransformPointerType::New();
		translation->Translate(shiftVec);
		
		std::list< int >::const_iterator labelIt;
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(ID >= Landmark)
			{
				workingSG->setTransformation(translation);
				workingSG->applyTransformation(ID);
			}
		}
		if(piaFlag)
		{
			TransformFilterType piaTransFilter = TransformFilterType::New();
			piaTransFilter->SetTransform(translation);
			piaTransFilter->SetInput(pia);
			piaTransFilter->Update();
			pia->DeepCopy(piaTransFilter->GetOutput());
		}
		if(wmFlag)
		{
			TransformFilterType wmTransFilter = TransformFilterType::New();
			wmTransFilter->SetTransform(translation);
			wmTransFilter->SetInput(wm);
			wmTransFilter->Update();
			wm->DeepCopy(wmTransFilter->GetOutput());
		}
	}
};

/******************************************************************************/
/*optional correction for systematic offsets in manual barrel contours        */
/*corrects position of landmarks such that z-scaling can be performed without */
/*any changes in implementation (with corrected scaling factors)              */
/******************************************************************************/
void Registration::correctManualContours(PolyDataPointerType pia, PolyDataPointerType wm, std::map< int, Column * > * barrels)
{
	if(workingSG)
	{
		int HBID = workingSG->getHomeBarrel();
		if(!HBID)
		{
			std::cout << "Error! No home barrel defined. Cannot correct manual contours." << std::endl;
			return;
		}
		double zAxis[3];
		for(int ii = 0; ii < 3; ++ii)
			zAxis[ii] = barrels[0][HBID]->top[ii] - barrels[0][HBID]->bottom[ii];
		normalize(zAxis);
		
		// this part for absolute mode
		double dSupra = contourCorrection[HBID][0];
		double dL4 = contourCorrection[HBID][1];
		double dInfra = contourCorrection[HBID][2];
		
		// this part for relative mode
// 		double * inputPiaPt = NULL;
// 		double * inputWMPt = NULL;
// 		if(piaFlag)
// 			inputPiaPt = axisSurfaceIntersection(pia, zAxis, barrels[0][HBID]->top);
// 		if(wmFlag)
// 			inputWMPt = axisSurfaceIntersection(wm, zAxis, barrels[0][HBID]->bottom);
// 		double inputTopDist = 0, inputBarrelHeight = 0, inputBottomDist = 0;
// 		inputBarrelHeight = barrels[0][HBID]->getHeight();
// 		if(inputPiaPt)
// 			inputTopDist = L2Distance3D(inputPiaPt, barrels[0][HBID]->top);
// 		if(inputWMPt)
// 			inputBottomDist = L2Distance3D(barrels[0][HBID]->bottom, inputWMPt);
// 		
// 		double dSupra = (contourCorrection[HBID][0] - 1)*inputTopDist;
// 		double dL4 = (contourCorrection[HBID][1] - 1)*inputBarrelHeight;
// 		double dInfra = (contourCorrection[HBID][2] - 1)*inputBottomDist;
		
		// step 1: shift home barrel top/bottom
		for(int ii = 0; ii < 3; ++ii)
		{
			barrels[0][HBID]->top[ii] += 0.5*dL4*zAxis[ii];
			barrels[0][HBID]->bottom[ii] -= 0.5*dL4*zAxis[ii];
		}
		
		// step 2: shift Pia & neuron
		double supraShift[3];
		for(int ii  = 0; ii < 3; ++ii)
			supraShift[ii] = (dSupra + 0.5*dL4)*zAxis[ii];
		TransformPointerType supraTrans = TransformPointerType::New();
		supraTrans->Translate(supraShift);
		
		for(int label = Neuron; label <= Soma; ++label)
			if(workingSG->isLabelInSpatialGraph(label))
			{
				workingSG->setTransformation(supraTrans);
				workingSG->applyTransformation(label);
			}
		if(piaFlag)
		{
			workingSG->setTransformation(supraTrans);
			workingSG->applyTransformation(Pia);
			if(pia)
			{
				TransformFilterType piaTransFilter = TransformFilterType::New();
				piaTransFilter->SetTransform(supraTrans);
				piaTransFilter->SetInput(pia);
				piaTransFilter->Update();
				pia->DeepCopy(piaTransFilter->GetOutput());
				#ifdef PIPELINE_DOC
				globalPia->DeepCopy(pia);
				#endif
			}
		}
		
		// step 3: shift WM
		if(wmFlag)
		{
			double infraShift[3];
			for(int ii  = 0; ii < 3; ++ii)
				infraShift[ii] = -1*(dInfra + 0.5*dL4)*zAxis[ii];
			TransformPointerType infraTrans = TransformPointerType::New();
			infraTrans->Translate(infraShift);
			workingSG->setTransformation(infraTrans);
			workingSG->applyTransformation(WhiteMatter);
			if(wm)
			{
				TransformFilterType WMTransFilter = TransformFilterType::New();
				WMTransFilter->SetTransform(infraTrans);
				WMTransFilter->SetInput(wm);
				WMTransFilter->Update();
				wm->DeepCopy(WMTransFilter->GetOutput());
			}
		}
	}
};

/******************************************************************************/
/*computes new z coordinates for all structures; strictly valid only for      */
/*home barrel. surfaces are only shifted, since they are only written         */
/*for visualization purposes anyways                                          */
/******************************************************************************/
void Registration::inhomogeneousZScaling(double zScaling[3], double refPts[4], double inputRefPts[4], PolyDataPointerType pia, PolyDataPointerType wm)
{
	if(workingSG)
	{
		std::vector< Vertex * >::iterator vertexIt;
		for(vertexIt = workingSG->verticesBegin(); vertexIt != workingSG->verticesEnd(); ++vertexIt)
		{
			if((*vertexIt)->label < Landmark || (*vertexIt)->label > WhiteMatter)
				continue;
			double oldZ = (*vertexIt)->coordinates[2];
			double refZ, inZ, scale;
			for(int ii = 1; ii < 4; ++ii)
			{
				if(ii == 3)
				{
					refZ = refPts[ii];
					inZ = inputRefPts[ii];
					scale = zScaling[ii-1];
					break;
				}
				if(oldZ >= inputRefPts[ii])
				{
					refZ = refPts[ii];
					inZ = inputRefPts[ii];
					scale = zScaling[ii-1];
					break;
				}
			}
			double newZ = refZ + scale*(oldZ - inZ);
			(*vertexIt)->coordinates[2] = newZ;
		}
		std::vector< Edge * >::iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label < Landmark || (*edgeIt)->label > WhiteMatter)
				continue;
			std::list< double * >::iterator edgePtIt;
			for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
			{
				double oldZ = (*edgePtIt)[2];
				double refZ, inZ, scale;
				for(int ii = 1; ii < 4; ++ii)
				{
					if(ii == 3)
					{
						refZ = refPts[ii];
						inZ = inputRefPts[ii];
						scale = zScaling[ii-1];
						break;
					}
					if(oldZ >= inputRefPts[ii])
					{
						refZ = refPts[ii];
						inZ = inputRefPts[ii];
						scale = zScaling[ii-1];
						break;
					}
				}
				double newZ = refZ + scale*(oldZ - inZ);
				(*edgePtIt)[2] = newZ;
			}
		}
	}
// 	double piaShift[] = {0,0,refPts[0]-inputRefPts[0]};
// 	TransformPointerType piaTrans = TransformPointerType::New();
// 	piaTrans->Translate(piaShift);
// 	TransformFilterType piaTransFilter = TransformFilterType::New();
// 	piaTransFilter->SetTransform(piaTrans);
// 	piaTransFilter->SetInput(pia);
// 	piaTransFilter->Update();
// 	pia->DeepCopy(piaTransFilter->GetOutput());
// 	
// 	if(wm)
// 	{
// 		double wmShift[] = {0,0,refPts[3]-inputRefPts[3]};
// 		TransformPointerType wmTrans = TransformPointerType::New();
// 		wmTrans->Translate(wmShift);
// 		TransformFilterType wmTransFilter = TransformFilterType::New();
// 		wmTransFilter->SetTransform(wmTrans);
// 		wmTransFilter->SetInput(wm);
// 		wmTransFilter->Update();
// 		wm->DeepCopy(wmTransFilter->GetOutput());
// 	}
};

/******************************************************************************/
/*computes new z coordinates for all structures; scaling factors are          */
/*calculated for the home barrel and assumed constant in the barrel field.    */
/*scaling is however computed with those scaling factors relative to local    */
/*z axis and local pia, wm and L4 surfaces                                    */
/******************************************************************************/
void Registration::inhomogeneousZScaling(double zScaling[3], double somaPt[3])
{
	if(SBF->avgPiaSurface && SBF->avgWMSurface && SBF->avgL4UpperSurface && SBF->avgL4LowerSurface)
	{
		std::flush(std::cout << ">>>        Applying inhomogeneous z-scaling" << std::endl);
		std::cout << std::endl;
		
		std::vector< Vertex * >::iterator vertexIt;
		for(vertexIt = workingSG->verticesBegin(); vertexIt != workingSG->verticesEnd(); ++vertexIt)
		{
			if((*vertexIt)->label < Neuron || (*vertexIt)->label > Soma)
				continue;
			
			computeNewPosition((*vertexIt)->coordinates, somaPt, zScaling);
		}
		
		std::vector< Edge * >::iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label < Neuron || (*edgeIt)->label > Soma)
				continue;
			
			std::list< double * >::iterator edgePtIt;
			for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
			{
				computeNewPosition(*edgePtIt, somaPt, zScaling);
			}
		}
	}
	else
		std::cout << "Error! No z-scaling possible because standard surfaces are incomplete!" << std::endl;
};

/******************************************************************************/
/*computes new z coordinates for all structures; scaling factors are          */
/*calculated for the home barrel and assumed constant in the barrel field.    */
/*scaling is however computed with those scaling factors relative to local    */
/*z axis and local pia, wm and L4 surfaces. Parallel version.                 */
/******************************************************************************/
// void Registration::inhomogeneousZScalingPRL(double zScaling[3], double somaPt[3])
// {
// 	if(avgPia->GetNumberOfCells() && avgWM->GetNumberOfCells() && avgL4Upper->GetNumberOfCells() && avgL4Lower->GetNumberOfCells())
// 	{
// 		std::flush(std::cout << ">>>        Applying inhomogeneous z-scaling" << std::endl);
// 		std::cout << std::endl;
// 		
// 		unsigned int stepsize = 8;
// // 		std::vector< Vertex * >::iterator vertexIt;
// 		for(unsigned int ii = 0; ii < workingSG->verticesPointer()->size(); ii += stepsize)
// 		{
// 			#pragma omp parallel for
// 			for(unsigned int jj = 0; jj < stepsize; ++jj)
// 			{
// 				unsigned int n = ii + jj;
// 				if(n >= workingSG->verticesPointer()->size())
// 					continue;
// 				if((*(workingSG->verticesPointer()))[n]->label < Neuron || (*(workingSG->verticesPointer()))[n]->label > Soma)
// 					continue;
// 				
// 				computeNewPosition((*(workingSG->verticesPointer()))[n]->coordinates, somaPt, zScaling);
// 			}
// 		}
// 		
// // 		std::vector< Edge * >::iterator edgeIt;
// 		for(unsigned int ii = 0; ii < workingSG->edgesPointer()->size(); ii += stepsize)
// 		{
// 			#pragma omp parallel for
// 			for(unsigned int jj = 0; jj < stepsize; ++jj)
// 			{
// 				unsigned int n = ii + jj;
// 				if(n >= workingSG->edgesPointer()->size())
// 					continue;
// 				if((*(workingSG->edgesPointer()))[n]->label < Neuron || (*(workingSG->edgesPointer()))[n]->label > Soma)
// 					continue;
// 				
// 				std::list< double * >::iterator edgePtIt;
// 				for(edgePtIt = (*(workingSG->edgesPointer()))[n]->edgePointCoordinates.begin();
// 				edgePtIt != (*(workingSG->edgesPointer()))[n]->edgePointCoordinates.end();
// 				++edgePtIt)
// 				{
// 					computeNewPosition(*edgePtIt, somaPt, zScaling);
// 				}
// 			}
// 		}
// 	}
// 	else
// 		std::cout << "Error! No z-scaling possible because standard surfaces are incomplete!" << std::endl;
// };

/******************************************************************************/
/*computes new z coordinates for all structures; scaling factors are          */
/*calculated for the home barrel and assumed constant in the barrel field.    */
/*scaling is however computed with those scaling factors relative to local    */
/*z axis and local pia, wm and L4 surfaces; newSomaPt is fixed pt!!!          */
/******************************************************************************/
void Registration::inhomogeneousZScaling(double zScaling[3], double relSomaPt[3], double oldSomaPt[3])
{
	if(SBF->avgPiaSurface && SBF->avgWMSurface && SBF->avgL4UpperSurface && SBF->avgL4LowerSurface)
	{
		std::flush(std::cout << ">>>        Applying inhomogeneous z scaling" << std::endl);
		std::cout << std::endl;
		
		std::vector< Vertex * >::iterator vertexIt;
		for(vertexIt = workingSG->verticesBegin(); vertexIt != workingSG->verticesEnd(); ++vertexIt)
		{
			if((*vertexIt)->label < Neuron || (*vertexIt)->label > Soma)
				continue;
			
			computeNewPosition((*vertexIt)->coordinates, relSomaPt, oldSomaPt, zScaling);
		}
		
		std::vector< Edge * >::iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label < Neuron || (*edgeIt)->label > Soma)
				continue;
			
			std::list< double * >::iterator edgePtIt;
			for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
			{
				computeNewPosition(*edgePtIt, relSomaPt, oldSomaPt, zScaling);
			}
		}
	}
	else
		std::cout << "Error! No z-scaling possible because standard surfaces are incomplete!" << std::endl;
};

/******************************************************************************/
/*computes new coordinates for all structures; scaling factors are            */
/*calculated for the home barrel and assumed constant in the barrel field.    */
/*scaling is however computed with those scaling factors relative to local    */
/*z axis and local pia, wm and L4 surfaces                                    */
/******************************************************************************/
void Registration::computeNewPosition(double * oldPt, double somaPt[3], double zScaling[3])
{
	double closestPt[3], correction[3];	//in case there is structure outside of standard surfaces -> use closest pt
// 	double * zAxis = localZAxis(oldPt);
	double zAxis[3];
	SBF->localZAxis(oldPt, zAxis);
	//use closestPt for all calculations and simply shift everything back with correction vector at the end
	closestPt[0] = oldPt[0], closestPt[1] = oldPt[1], closestPt[2] = oldPt[2];
	correction[0] = oldPt[0] - closestPt[0], correction[1] = oldPt[1] - closestPt[1], correction[2] = oldPt[2] - closestPt[2];
	
	SBF->avgPiaSurface->intersectLine(zAxis, closestPt);
	SBF->avgL4UpperSurface->intersectLine(zAxis, closestPt);
	SBF->avgL4LowerSurface->intersectLine(zAxis, closestPt);
	double * piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
	double * L4UIntersectPt = SBF->avgL4UpperSurface->getLastIntersectPoint();
	double * L4LIntersectPt = SBF->avgL4LowerSurface->getLastIntersectPoint();
	double * WMIntersectPt = NULL;
	//heuristic safeguard in case axis at edge of axis field is just outside of interpolated surfaces
	//maybe fix better by interpolating surfaces further than axes
	int safeCnt = 0;
	while((!piaIntersectPt || !L4UIntersectPt || !L4LIntersectPt) && safeCnt < 10)
	{
		++safeCnt;
		double delta[3], tmp[3];
		SBF->localZAxis(oldPt, closestPt, tmp);
		if(safeCnt == 1)
			correction[0] = oldPt[0] - closestPt[0], correction[1] = oldPt[1] - closestPt[1], correction[2] = oldPt[2] - closestPt[2];
		delta[0] = somaPt[0] - closestPt[0], delta[1] = somaPt[1] - closestPt[1], delta[2] = somaPt[2] - closestPt[2];
		normalize(delta);
		delta[0] = 10*delta[0], delta[1] = 10*delta[1], delta[2] = 10*delta[2];
		correction[0] -= delta[0], correction[1] -= delta[1], correction[2] -= delta[2];
		closestPt[0] += delta[0], closestPt[1] += delta[1], closestPt[2] += delta[2];
		SBF->avgPiaSurface->intersectLine(zAxis, closestPt);
		SBF->avgL4UpperSurface->intersectLine(zAxis, closestPt);
		SBF->avgL4LowerSurface->intersectLine(zAxis, closestPt);
		piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
		L4UIntersectPt = SBF->avgL4UpperSurface->getLastIntersectPoint();
		L4LIntersectPt = SBF->avgL4LowerSurface->getLastIntersectPoint();
		if(zScaling[2] != 1)
		{
			SBF->avgWMSurface->intersectLine(zAxis, closestPt);
			WMIntersectPt = SBF->avgWMSurface->getLastIntersectPoint();
		}
	}
	
	// TEMPORARY!!!! ignoring the problem does not solve it...
	if(piaIntersectPt && L4UIntersectPt && L4LIntersectPt)
	{
		double L4CenterPt[3], localZScale[3];
		for(int ii = 0; ii < 3; ++ii)
		{
			L4CenterPt[ii] = 0.5*(L4UIntersectPt[ii] + L4LIntersectPt[ii]);
			localZScale[ii] = zScaling[ii];
		}
		
		double refTopDist = L2Distance3D(piaIntersectPt, L4UIntersectPt);
		double refBarrelHeight = L2Distance3D(L4UIntersectPt, L4LIntersectPt);
		double refBottomDist = 0;
		if(localZScale[2] != 1 && WMIntersectPt)
			refBottomDist = L2Distance3D(L4LIntersectPt, WMIntersectPt);
		
		std::vector< double * > refPts, inputRefPts;
		refPts.push_back(piaIntersectPt), refPts.push_back(L4UIntersectPt), refPts.push_back(L4LIntersectPt);
		if(localZScale[2] != 1 && WMIntersectPt)
			refPts.push_back(WMIntersectPt);
		else
			refPts.push_back(L4LIntersectPt);
		
		//now, reconstruct virtual input reference points
		double inTopDist, inBarrelHeight, inBottomDist;
		double inPiaPt[3], inL4UPt[3], inL4LPt[3], inWMPt[3];
		inTopDist = refTopDist/localZScale[0];
		inBarrelHeight = refBarrelHeight/localZScale[1];
		inBottomDist = refBottomDist/localZScale[2];
		for(int ii = 0; ii < 3; ++ii)
		{
			inL4UPt[ii] = L4CenterPt[ii] + 0.5*inBarrelHeight*zAxis[ii];
			inL4LPt[ii] = L4CenterPt[ii] - 0.5*inBarrelHeight*zAxis[ii];
			inPiaPt[ii] = inL4UPt[ii] + inTopDist*zAxis[ii];
			inWMPt[ii] = inL4LPt[ii] - inBottomDist*zAxis[ii];
		}
		inputRefPts.push_back(inPiaPt);
		inputRefPts.push_back(inL4UPt);
		inputRefPts.push_back(inL4LPt);
		inputRefPts.push_back(inWMPt);
		
		double * refPt, * inPt, scale;
		for(int ii = 1; ii < 4; ++ii)
		{
			if(ii == 3)
			{
				refPt = refPts[ii];
				inPt = inputRefPts[ii];
				scale = localZScale[ii-1];
				break;
			}
			if(closestPt[2] >= inputRefPts[ii][2])
			{
				refPt = refPts[ii];
				inPt = inputRefPts[ii];
				scale = localZScale[ii-1];
				break;
			}
		}
		for(int ii = 0; ii < 3; ++ii)
			oldPt[ii] = refPt[ii] + scale*(closestPt[ii] - inPt[ii]) + correction[ii];
		
		refPts.clear(), inputRefPts.clear();
		delete [] piaIntersectPt, delete [] L4UIntersectPt, delete [] L4LIntersectPt;
	}
	if(WMIntersectPt) delete [] WMIntersectPt;
// 	delete [] zAxis;
};

/******************************************************************************/
/*computes new coordinates for all structures; scaling factors are            */
/*calculated for the home barrel and assumed constant in the barrel field.    */
/*scaling is however computed with those scaling factors relative to local    */
/*z axis and local pia, wm and L4 surfaces                                    */
/******************************************************************************/
void Registration::computeNewPosition(double * oldPt, double somaPt[3], double zScaling[3], bool verbose)
{
	double closestPt[3], correction[3];	//in case there is structure outside of standard surfaces -> use closest pt
// 	double * zAxis = localZAxis(oldPt);
	double zAxis[3];
	SBF->localZAxis(oldPt, zAxis);
	//use closestPt for all calculations and simply shift everything back with correction vector at the end
	closestPt[0] = oldPt[0], closestPt[1] = oldPt[1], closestPt[2] = oldPt[2];
	correction[0] = oldPt[0] - closestPt[0], correction[1] = oldPt[1] - closestPt[1], correction[2] = oldPt[2] - closestPt[2];
	
	SBF->avgPiaSurface->intersectLine(zAxis, closestPt);
	SBF->avgL4UpperSurface->intersectLine(zAxis, closestPt);
	SBF->avgL4LowerSurface->intersectLine(zAxis, closestPt);
	double * piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
	double * L4UIntersectPt = SBF->avgL4UpperSurface->getLastIntersectPoint();
	double * L4LIntersectPt = SBF->avgL4LowerSurface->getLastIntersectPoint();
	double * WMIntersectPt = NULL;
	//heuristic safeguard in case axis at edge of axis field is just outside of interpolated surfaces
	//maybe fix better by interpolating surfaces further than axes
	int safeCnt = 0;
	while((!piaIntersectPt || !L4UIntersectPt || !L4LIntersectPt) && safeCnt < 10)
	{
		++safeCnt;
		double delta[3], tmp[3];
		SBF->localZAxis(oldPt, closestPt, tmp);
		if(safeCnt == 1)
			correction[0] = oldPt[0] - closestPt[0], correction[1] = oldPt[1] - closestPt[1], correction[2] = oldPt[2] - closestPt[2];
		delta[0] = somaPt[0] - closestPt[0], delta[1] = somaPt[1] - closestPt[1], delta[2] = somaPt[2] - closestPt[2];
		normalize(delta);
		delta[0] = 10*delta[0], delta[1] = 10*delta[1], delta[2] = 10*delta[2];
		correction[0] -= delta[0], correction[1] -= delta[1], correction[2] -= delta[2];
		closestPt[0] += delta[0], closestPt[1] += delta[1], closestPt[2] += delta[2];
		SBF->avgPiaSurface->intersectLine(zAxis, closestPt);
		SBF->avgL4UpperSurface->intersectLine(zAxis, closestPt);
		SBF->avgL4LowerSurface->intersectLine(zAxis, closestPt);
		piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
		L4UIntersectPt = SBF->avgL4UpperSurface->getLastIntersectPoint();
		L4LIntersectPt = SBF->avgL4LowerSurface->getLastIntersectPoint();
		if(zScaling[2] != 1)
		{
			SBF->avgWMSurface->intersectLine(zAxis, closestPt);
			WMIntersectPt = SBF->avgWMSurface->getLastIntersectPoint();
		}
	}
	
	// TEMPORARY!!!! ignoring the problem does not solve it...
	if(piaIntersectPt && L4UIntersectPt && L4LIntersectPt)
	{
		double L4CenterPt[3], localZScale[3];
		for(int ii = 0; ii < 3; ++ii)
		{
			L4CenterPt[ii] = 0.5*(L4UIntersectPt[ii] + L4LIntersectPt[ii]);
			localZScale[ii] = zScaling[ii];
		}
		
		double refTopDist = L2Distance3D(piaIntersectPt, L4UIntersectPt);
		double refBarrelHeight = L2Distance3D(L4UIntersectPt, L4LIntersectPt);
		double refBottomDist = 0;
		if(localZScale[2] != 1 && WMIntersectPt)
			refBottomDist = L2Distance3D(L4LIntersectPt, WMIntersectPt);
		
		std::vector< double * > refPts, inputRefPts;
		refPts.push_back(piaIntersectPt), refPts.push_back(L4UIntersectPt), refPts.push_back(L4LIntersectPt);
		if(localZScale[2] != 1 && WMIntersectPt)
			refPts.push_back(WMIntersectPt);
		else
			refPts.push_back(L4LIntersectPt);
		
		//now, reconstruct virtual input reference points
		double inTopDist, inBarrelHeight, inBottomDist;
		double inPiaPt[3], inL4UPt[3], inL4LPt[3], inWMPt[3];
		inTopDist = refTopDist/localZScale[0];
		inBarrelHeight = refBarrelHeight/localZScale[1];
		inBottomDist = refBottomDist/localZScale[2];
		for(int ii = 0; ii < 3; ++ii)
		{
			inL4UPt[ii] = L4CenterPt[ii] + 0.5*inBarrelHeight*zAxis[ii];
			inL4LPt[ii] = L4CenterPt[ii] - 0.5*inBarrelHeight*zAxis[ii];
			inPiaPt[ii] = inL4UPt[ii] + inTopDist*zAxis[ii];
			inWMPt[ii] = inL4LPt[ii] - inBottomDist*zAxis[ii];
		}
		inputRefPts.push_back(inPiaPt);
		inputRefPts.push_back(inL4UPt);
		inputRefPts.push_back(inL4LPt);
		inputRefPts.push_back(inWMPt);
		
		// verbose output for debugging
		std::cout << "oldPt @ [" << oldPt[0] << "," << oldPt[1] << "," << oldPt[2] << "]" << std::endl;
		std::cout << "closestPt @ [" << closestPt[0] << "," << closestPt[1] << "," << closestPt[2] << "]" << std::endl;
		std::cout << "z axis = [" << zAxis[0] << "," << zAxis[1] << "," << zAxis[2] << "]" << std::endl;
		if(piaIntersectPt)
			std::cout << "piaIntersectPt @ [" << piaIntersectPt[0] << "," << piaIntersectPt[1] << "," << piaIntersectPt[2] << "]" << std::endl;
		if(L4UIntersectPt)
			std::cout << "L4UIntersectPt @ [" << L4UIntersectPt[0] << "," << L4UIntersectPt[1] << "," << L4UIntersectPt[2] << "]" << std::endl;
		if(L4LIntersectPt)
			std::cout << "L4LIntersectPt @ [" << L4LIntersectPt[0] << "," << L4LIntersectPt[1] << "," << L4LIntersectPt[2] << "]" << std::endl;
		if(WMIntersectPt)
			std::cout << "WMIntersectPt @ [" << WMIntersectPt[0] << "," << WMIntersectPt[1] << "," << WMIntersectPt[2] << "]" << std::endl;
		std::cout << std::endl;
		
		double * refPt, * inPt, scale;
		for(int ii = 1; ii < 4; ++ii)
		{
			if(ii == 3)
			{
				refPt = refPts[ii];
				inPt = inputRefPts[ii];
				scale = localZScale[ii-1];
				break;
			}
			if(closestPt[2] >= inputRefPts[ii][2])
			{
				refPt = refPts[ii];
				inPt = inputRefPts[ii];
				scale = localZScale[ii-1];
				break;
			}
		}
		for(int ii = 0; ii < 3; ++ii)
			oldPt[ii] = refPt[ii] + scale*(closestPt[ii] - inPt[ii]) + correction[ii];
		
		refPts.clear(), inputRefPts.clear();
		delete [] piaIntersectPt, delete [] L4UIntersectPt, delete [] L4LIntersectPt;
	}
	else
	{
		if(!piaIntersectPt)
			std::cout << "Warning! piaIntersectPt not found!" << std::endl;
		if(!L4UIntersectPt)
			std::cout << "Warning! L4UIntersectPt not found!" << std::endl;
		if(!L4LIntersectPt)
			std::cout << "Warning! L4LIntersectPt not found!" << std::endl;
	}
	if(WMIntersectPt) delete [] WMIntersectPt;
// 	delete [] zAxis;
};

/******************************************************************************/
/*computes new coordinates for all structures; scaling factors are            */
/*calculated for the home barrel and assumed constant in the barrel field.    */
/*scaling is however computed with those scaling factors relative to local    */
/*z axis and local pia, wm and L4 surfaces. one reference pt is replace with  */
/*the soma point -> for use in re-registering cell to different column        */
/*(preserves relative soma position along vertical axis; that's why the new)  */
/*soma pt is a double[3]: 0 <= newSomaPt[i] <= 1 -> position along column     */
/*for each layer separately                                                   */
/******************************************************************************/
void Registration::computeNewPosition(double * oldPt, double relSomaPt[3], double oldSomaPt[3], double zScaling[3])
{
	double closestPt[3], correction[3];	//in case there is structure outside of standard surfaces -> use closest pt
// 	double * zAxis = localZAxis(oldPt);
	double zAxis[3];
	SBF->localZAxis(oldPt, zAxis);
	//use closestPt for all calculations and simply shift everything back with correction vector at the end
	closestPt[0] = oldPt[0], closestPt[1] = oldPt[1], closestPt[2] = oldPt[2];
	correction[0] = oldPt[0] - closestPt[0], correction[1] = oldPt[1] - closestPt[1], correction[2] = oldPt[2] - closestPt[2];
	
	SBF->avgPiaSurface->intersectLine(zAxis, closestPt);
	SBF->avgL4UpperSurface->intersectLine(zAxis, closestPt);
	SBF->avgL4LowerSurface->intersectLine(zAxis, closestPt);
	double * piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
	double * L4UIntersectPt = SBF->avgL4UpperSurface->getLastIntersectPoint();
	double * L4LIntersectPt = SBF->avgL4LowerSurface->getLastIntersectPoint();
	double * WMIntersectPt = NULL;
	//heuristic safeguard in case axis at edge of axis field is just outside of interpolated surfaces
	//maybe fix better by interpolating surfaces further than axes
	int safeCnt = 0;
	while((!piaIntersectPt || !L4UIntersectPt || !L4LIntersectPt) && safeCnt < 10)
	{
		++safeCnt;
		double delta[3], tmp[3];
		SBF->localZAxis(oldPt, closestPt, tmp);
		if(safeCnt == 1)
			correction[0] = oldPt[0] - closestPt[0], correction[1] = oldPt[1] - closestPt[1], correction[2] = oldPt[2] - closestPt[2];
		delta[0] = oldSomaPt[0] - closestPt[0], delta[1] = oldSomaPt[1] - closestPt[1], delta[2] = oldSomaPt[2] - closestPt[2];
		normalize(delta);
		delta[0] = 10*delta[0], delta[1] = 10*delta[1], delta[2] = 10*delta[2];
		correction[0] -= delta[0], correction[1] -= delta[1], correction[2] -= delta[2];
		closestPt[0] += delta[0], closestPt[1] += delta[1], closestPt[2] += delta[2];
		SBF->avgPiaSurface->intersectLine(zAxis, closestPt);
		SBF->avgL4UpperSurface->intersectLine(zAxis, closestPt);
		SBF->avgL4LowerSurface->intersectLine(zAxis, closestPt);
		piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
		L4UIntersectPt = SBF->avgL4UpperSurface->getLastIntersectPoint();
		L4LIntersectPt = SBF->avgL4LowerSurface->getLastIntersectPoint();
		if(zScaling[2] != 1)
		{
			SBF->avgWMSurface->intersectLine(zAxis, closestPt);
			WMIntersectPt = SBF->avgWMSurface->getLastIntersectPoint();
		}
	}
	
	// TEMPORARY!!!! ignoring the problem does not solve it...
	if(piaIntersectPt && L4UIntersectPt && L4LIntersectPt)
	{
		double L4CenterPt[3], localZScale[3];
		for(int ii = 0; ii < 3; ++ii)
		{
			L4CenterPt[ii] = 0.5*(L4UIntersectPt[ii] + L4LIntersectPt[ii]);
			localZScale[ii] = zScaling[ii];
		}
		
		double refTopDist = L2Distance3D(piaIntersectPt, L4UIntersectPt);
		double refBarrelHeight = L2Distance3D(L4UIntersectPt, L4LIntersectPt);
		double refBottomDist = 0;
		if(localZScale[2] != 1 && WMIntersectPt)
			refBottomDist = L2Distance3D(L4LIntersectPt, WMIntersectPt);
		
		std::vector< double * > refPts, inputRefPts;
		refPts.push_back(piaIntersectPt), refPts.push_back(L4UIntersectPt), refPts.push_back(L4LIntersectPt);
		if(localZScale[2] != 1 && WMIntersectPt)
			refPts.push_back(WMIntersectPt);
		else
			refPts.push_back(L4LIntersectPt);
		
		//now, reconstruct virtual input reference points
		double inTopDist, inBarrelHeight, inBottomDist;
		double inPiaPt[3], inL4UPt[3], inL4LPt[3], inWMPt[3];
		inTopDist = refTopDist/localZScale[0];
		inBarrelHeight = refBarrelHeight/localZScale[1];
		inBottomDist = refBottomDist/localZScale[2];
		for(int ii = 0; ii < 3; ++ii)
		{
			inL4UPt[ii] = L4CenterPt[ii] + 0.5*inBarrelHeight*zAxis[ii];
			inL4LPt[ii] = L4CenterPt[ii] - 0.5*inBarrelHeight*zAxis[ii];
			inPiaPt[ii] = inL4UPt[ii] + inTopDist*zAxis[ii];
			inWMPt[ii] = inL4LPt[ii] - inBottomDist*zAxis[ii];
		}
		inputRefPts.push_back(inPiaPt);
		inputRefPts.push_back(inL4UPt);
		inputRefPts.push_back(inL4LPt);
		inputRefPts.push_back(inWMPt);
		
		double * refPt, * inPt, scale;
		int layer;
		for(int ii = 1; ii < 4; ++ii)
		{
			if(ii == 3)
			{
				refPt = refPts[ii];
				inPt = inputRefPts[ii];
				scale = localZScale[ii-1];
				layer = ii;
				break;
			}
			if(closestPt[2] >= inputRefPts[ii][2])
			{
				refPt = refPts[ii];
				inPt = inputRefPts[ii];
				scale = localZScale[ii-1];
				layer = ii;
				break;
			}
		}
		// check if reference pts need to be adjusted
		// to soma position (only if in same layer)
		if(relSomaPt[0] < 1 && layer == SUPRA)
			for(int ii = 0; ii < 3; ++ii)
			{
				refPt[ii] = piaIntersectPt[ii] - relSomaPt[0]*refTopDist*zAxis[ii];
				inPt[ii] = inPiaPt[ii] - relSomaPt[0]*inTopDist*zAxis[ii];
			}
		else if(relSomaPt[1] < 1 && layer == GRAN)
			for(int ii = 0; ii < 3; ++ii)
			{
				refPt[ii] = piaIntersectPt[ii] - (refTopDist + relSomaPt[1]*refBarrelHeight)*zAxis[ii];
				inPt[ii] = inPiaPt[ii] - (inTopDist + relSomaPt[1]*inBarrelHeight)*zAxis[ii];
			}
		else if(relSomaPt[2] <= 1 && layer == INFRA && WMIntersectPt)
			for(int ii = 0; ii < 3; ++ii)
			{
				refPt[ii] = piaIntersectPt[ii] - (refTopDist + refBarrelHeight + relSomaPt[2]*refBottomDist)*zAxis[ii];
				inPt[ii] = inPiaPt[ii] - (inTopDist + inBarrelHeight + relSomaPt[2]*inBottomDist)*zAxis[ii];
			}
		
		for(int ii = 0; ii < 3; ++ii)
			oldPt[ii] = refPt[ii] + scale*(closestPt[ii] - inPt[ii]) + correction[ii];
		
		refPts.clear(), inputRefPts.clear();
		delete [] piaIntersectPt, delete [] L4UIntersectPt, delete [] L4LIntersectPt;
	}
	if(WMIntersectPt) delete [] WMIntersectPt;
// 	delete [] zAxis;
};

/******************************************************************************/
/*check the whole cell morphology to make sure that no structure is going to  */
/*"stick out" of the pia.                                                     */
/*returns a look-up table of z-scale factors where the key corresponds to the */
/*cell ID of local cell in avgL4Lower -> to be used in conjunction with the   */
/*cell locator of this surface                                                */
/******************************************************************************/
std::map< vtkIdType, std::vector< double > > Registration::buildLocalZScaleCorrection(double zScaling[3], double somaPt[3])
{
	if(SBF->avgPiaSurface && SBF->avgWMSurface && SBF->avgL4UpperSurface && SBF->avgL4LowerSurface)
	{
		std::map< vtkIdType, double > maxZTable;
		std::map< vtkIdType, double > piaTable;
		std::vector< Edge * >::const_iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			// for now: apply z-scale correction only to dendrite 
			if((*edgeIt)->label < Neuron || (*edgeIt)->label > Soma /*|| (*edgeIt)->label == Axon*/)
				continue;
			
			std::list< double * >::const_iterator edgePtIt;
			for(edgePtIt = (*edgeIt)->edgePointCoordinates.begin(); edgePtIt != (*edgeIt)->edgePointCoordinates.end(); ++edgePtIt)
			{
				double * oldPt = *edgePtIt;
				double closestPt[3], correction[3];	//in case there is structure outside of standard surfaces -> use closest pt
// 				double * zAxis = localZAxis(oldPt);
				double zAxis[3];
				SBF->localZAxis(oldPt, zAxis);
				//use closestPt for all calculations and simply shift everything back with correction vector at the end
				closestPt[0] = oldPt[0], closestPt[1] = oldPt[1], closestPt[2] = oldPt[2];
				correction[0] = oldPt[0] - closestPt[0], correction[1] = oldPt[1] - closestPt[1], correction[2] = oldPt[2] - closestPt[2];
				
				SBF->avgPiaSurface->intersectLine(zAxis, closestPt);
				SBF->avgL4UpperSurface->intersectLine(zAxis, closestPt);
				SBF->avgL4LowerSurface->intersectLine(zAxis, closestPt);
				double * piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
				double * L4UIntersectPt = SBF->avgL4UpperSurface->getLastIntersectPoint();
				double * L4LIntersectPt = SBF->avgL4LowerSurface->getLastIntersectPoint();
				vtkIdType intersectID = SBF->avgL4UpperSurface->getLastIntersectCellID();
				double * WMIntersectPt = NULL;
				//heuristic safeguard in case axis at edge of axis field is just outside of interpolated surfaces
				//maybe fix better by interpolating surfaces further than axes
				int safeCnt = 0;
				while((!piaIntersectPt || !L4UIntersectPt || !L4LIntersectPt) && safeCnt < 10)
				{
					++safeCnt;
					double delta[3];
					delta[0] = somaPt[0] - closestPt[0], delta[1] = somaPt[1] - closestPt[1], delta[2] = somaPt[2] - closestPt[2];
					normalize(delta);
					delta[0] = 10*delta[0], delta[1] = 10*delta[1], delta[2] = 10*delta[2];
					correction[0] -= delta[0], correction[1] -= delta[1], correction[2] -= delta[2];
					closestPt[0] += delta[0], closestPt[1] += delta[1], closestPt[2] += delta[2];
					SBF->avgPiaSurface->intersectLine(zAxis, closestPt);
					SBF->avgL4UpperSurface->intersectLine(zAxis, closestPt);
					SBF->avgL4LowerSurface->intersectLine(zAxis, closestPt);
					SBF->avgWMSurface->intersectLine(zAxis, closestPt);
					piaIntersectPt = SBF->avgPiaSurface->getLastIntersectPoint();
					L4UIntersectPt = SBF->avgL4UpperSurface->getLastIntersectPoint();
					L4LIntersectPt = SBF->avgL4LowerSurface->getLastIntersectPoint();
					intersectID = SBF->avgL4UpperSurface->getLastIntersectCellID();
					if(zScaling[2] != 1)
					{
						SBF->avgWMSurface->intersectLine(zAxis, closestPt);
						WMIntersectPt = SBF->avgWMSurface->getLastIntersectPoint();
					}
				}
				
				// TEMPORARY!!!! ignoring the problem does not solve it...
				if(piaIntersectPt && L4UIntersectPt && L4LIntersectPt)
				{
					double L4CenterPt[3];
					for(int ii = 0; ii < 3; ++ii)
						L4CenterPt[ii] = 0.5*(L4UIntersectPt[ii] + L4LIntersectPt[ii]);
					
					double refTopDist = L2Distance3D(piaIntersectPt, L4UIntersectPt);
					double refBarrelHeight = L2Distance3D(L4UIntersectPt, L4LIntersectPt);
					double refBottomDist = 0;
					if(zScaling[2] != 1 && WMIntersectPt)
						refBottomDist = L2Distance3D(L4LIntersectPt, WMIntersectPt);
					
					std::vector< double * > refPts, inputRefPts;
					refPts.push_back(piaIntersectPt), refPts.push_back(L4UIntersectPt), refPts.push_back(L4LIntersectPt);
					if(zScaling[2] != 1 && WMIntersectPt)
						refPts.push_back(WMIntersectPt);
					else
						refPts.push_back(L4LIntersectPt);
					
					//now, reconstruct virtual input reference points
					double inTopDist, inBarrelHeight, inBottomDist;
					double inPiaPt[3], inL4UPt[3], inL4LPt[3], inWMPt[3];
					inTopDist = refTopDist/zScaling[0];
					inBarrelHeight = refBarrelHeight/zScaling[1];
					inBottomDist = refBottomDist/zScaling[2];
					for(int ii = 0; ii < 3; ++ii)
					{
						inL4UPt[ii] = L4CenterPt[ii] + 0.5*inBarrelHeight*zAxis[ii];
						inL4LPt[ii] = L4CenterPt[ii] - 0.5*inBarrelHeight*zAxis[ii];
						inPiaPt[ii] = inL4UPt[ii] + inTopDist*zAxis[ii];
						inWMPt[ii] = inL4LPt[ii] - refBottomDist*zAxis[ii];
					}
					inputRefPts.push_back(inPiaPt);
					inputRefPts.push_back(inL4UPt);
					inputRefPts.push_back(inL4LPt);
					inputRefPts.push_back(inWMPt);
					
					double * refPt, * inPt, scale;
					for(int ii = 1; ii < 4; ++ii)
					{
						if(ii == 3)
						{
							refPt = refPts[ii];
							inPt = inputRefPts[ii];
							scale = zScaling[ii-1];
							break;
						}
						if(closestPt[2] >= inputRefPts[ii][2])
						{
							refPt = refPts[ii];
							inPt = inputRefPts[ii];
							scale = zScaling[ii-1];
							break;
						}
					}
					double tmpRegPt[3];
					for(int ii = 0; ii < 3; ++ii)
						tmpRegPt[ii] = refPt[ii] + scale*(closestPt[ii] - inPt[ii]);
					if(tmpRegPt[2] > piaIntersectPt[2])
					{
						if(maxZTable.find(intersectID) == maxZTable.end())
						{
							maxZTable.insert(std::pair< vtkIdType, double >(intersectID, L2Distance3D(tmpRegPt, L4UIntersectPt)));
							piaTable.insert(std::pair< vtkIdType, double >(intersectID, L2Distance3D(piaIntersectPt, L4UIntersectPt)));
						}
						else if(L2Distance3D(tmpRegPt, L4UIntersectPt) > maxZTable[intersectID])
							maxZTable[intersectID] = L2Distance3D(tmpRegPt, L4UIntersectPt);
					}
					
					refPts.clear(), inputRefPts.clear();
					/*delete [] zAxis,*/ delete [] piaIntersectPt, delete [] L4UIntersectPt, delete [] L4LIntersectPt;
				}
				if(WMIntersectPt) delete [] WMIntersectPt;
			}
		}
		
		std::map< vtkIdType, std::vector< double > > correctionTable;
		for(int ii = 0; ii < SBF->avgL4LowerSurface->ptr()->GetNumberOfCells(); ++ii)
			if(maxZTable.find(ii) != maxZTable.end() && piaTable.find(ii) != piaTable.end())
			{
				std::vector< double > correctZScale;
				double piaCorrection = piaTable[ii]/maxZTable[ii]*zScaling[0];
				correctZScale.push_back(piaCorrection);
				correctZScale.push_back(zScaling[1]);
				correctZScale.push_back(zScaling[2]);
				correctionTable.insert(std::pair< vtkIdType, std::vector< double > >(ii, correctZScale));
			}
		return correctionTable;
	}
	else
	{
		std::cout << "Error! No z-scaling possible because standard surfaces are incomplete!" << std::endl;
		std::map< vtkIdType, std::vector< double > > emptyMap;
		return emptyMap;
	}
};

/******************************************************************************/
/*check whether new cell position still corresponds to manually marked home   */
/*barrel. if not, repeat registration w.r.t. new home barrel                  */
/******************************************************************************/
int Registration::registeredHomeBarrel()
{
	if(workingSG && somaFlag)
	{
		int HBID = 0;
		double somaCenter[3];
		getPCenterOfStructure(workingSG, Soma, somaCenter);
		HBID = SBF->closestBarrel(somaCenter);
		
		return HBID;
	}
	else
	{
		std::cout << "WARNING! Cannot determine home barrel (soma = False)" << std::endl;
		return 0;
	}
};

/******************************************************************************/
/*closest barrel in case home barrel is not available in input as reference   */
/******************************************************************************/
int Registration::closestAvgBarrel()
{
	if(workingSG && somaFlag)
	{
		PolyDataPointerType soma = PolyDataPointerType::New();
		if(workingSG->extractLandmark(Soma, soma))
		{
			int HBID = 0;
			double somaCenter[3], minDist = 1E08;
			//parametric center
			int subID;
			double pCoords[3], * weights;
			weights = new double[soma->GetCell(0)->GetNumberOfPoints()];
			soma->GetCell(0)->GetParametricCenter(pCoords);
			soma->GetCell(0)->EvaluateLocation(subID, pCoords, somaCenter, weights);
			delete [] weights;
			std::list< int >::const_iterator labelIt;
			for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
			{
				int ID = *labelIt;
				if(SBF->avgCenters.find(ID) != SBF->avgCenters.end() && SBF->avgAxes.find(ID) != SBF->avgAxes.end() && workingSG->isLabelInSpatialGraph(ID))
				{
					double linePt1[3], linePt2[3], closestPt[3], t;
					for(int ii = 0; ii < 3; ++ii)
					{
						linePt1[ii] = SBF->avgCenters[ID][ii] + 1500*SBF->avgAxes[ID][ii];
						linePt2[ii] = SBF->avgCenters[ID][ii] - 1500*SBF->avgAxes[ID][ii];
					}
					double dist = vtkLine::DistanceToLine(somaCenter, linePt1, linePt2, t, closestPt);
					dist = sqrt(dist);
					if(dist < minDist)
					{
						minDist = dist;
						HBID = ID;
					}
				}
			}
			return HBID;
		}
		else
		{
			std::cout << "Error! Soma not in SpatialGraph! Could not determine registered home barrel" << std::endl;
			return 0;
		}
	}
	else
	{
		std::cout << "WARNING! Cannot determine home barrel (soma = False)" << std::endl;
		return 0;
	}
};

/******************************************************************************/
/*closest barrel in case home barrel is not available in input as reference   */
/******************************************************************************/
int Registration::closestMorphBarrel()
{
	if(workingSG && somaFlag)
	{
		PolyDataPointerType soma = PolyDataPointerType::New();
		if(workingSG->extractLandmark(Soma, soma))
		{
			int HBID = 0;
			double somaCenter[3], minDist = 1E08;
			//parametric center
			int subID;
			double pCoords[3], * weights;
			weights = new double[soma->GetCell(0)->GetNumberOfPoints()];
			soma->GetCell(0)->GetParametricCenter(pCoords);
			soma->GetCell(0)->EvaluateLocation(subID, pCoords, somaCenter, weights);
			delete [] weights;
			std::list< int >::const_iterator labelIt;
			for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
			{
				int ID = *labelIt;
				PolyDataPointerType barrel = PolyDataPointerType::New();
				if(workingSG->extractLandmark(ID, barrel))
				{
					double * barrelCentroid = calculateBarrelCentroid(barrel);
					double dist = L2Distance3D(barrelCentroid, somaCenter);
					if(dist < minDist)
					{
						minDist = dist;
						HBID = ID;
					}
				}
			}
			return HBID;
		}
		else
		{
			std::cout << "Error! Soma not in SpatialGraph! Could not determine home barrel" << std::endl;
			return 0;
		}
	}
	else
	{
		std::cout << "WARNING! Cannot determine home barrel (soma = False)" << std::endl;
		return 0;
	}
};

/******************************************************************************/
/*transform everything into tmp coordinate system with home barrel @ origin   */
/*and home barrel axis vertical; makes z-scaling easier to implement          */
/******************************************************************************/
void Registration::tmpCoordinateSystem(PolyDataPointerType pia, PolyDataPointerType wm, std::map< int, Column * > * barrels)
{
	int HBID = workingSG->getHomeBarrel();
	double tmpAxis[] = {0,0,1}, shift[3];
	for(int ii = 0; ii < 3; ++ii)
		shift[ii] = -SBF->avgCenters[HBID][ii];
	HomogeneousMatrixPointerType tmpRot = transformToBarrelCoordinates(SBF->avgAxes[HBID], tmpAxis);
	TransformPointerType tmpShift = TransformPointerType::New();
	TransformPointerType tmpRotTrans = TransformPointerType::New();
	tmpShift->Translate(shift);
	tmpRotTrans->SetMatrix(tmpRot);
	
	barrels[0][HBID]->translateColumn(shift);
	barrels[0][HBID]->rotateColumn(tmpRot);
	
	tmpRotTrans->Concatenate(tmpShift);
	tmpRotTrans->Update();
	
	workingSG->setTransformation(tmpRotTrans);
	workingSG->applyTransformation();
	
	TransformFilterType piaTmpTrans = TransformFilterType::New();
	piaTmpTrans->SetTransform(tmpRotTrans);
	piaTmpTrans->SetInput(pia);
	piaTmpTrans->Update();
	pia->DeepCopy(piaTmpTrans->GetOutput());
	if(wm)
	{
		TransformFilterType wmTmpTrans = TransformFilterType::New();
		wmTmpTrans->SetTransform(tmpRotTrans);
		wmTmpTrans->SetInput(wm);
		wmTmpTrans->Update();
		wm->DeepCopy(wmTmpTrans->GetOutput());
	}
};

/******************************************************************************/
/*transform everything back from tmp coordinate system into reference system  */
/******************************************************************************/
void Registration::tmpCoordinateSystemInv(PolyDataPointerType pia, PolyDataPointerType wm, std::map< int, Column * > * barrels)
{
	int HBID = workingSG->getHomeBarrel();
	double tmpAxis[] = {0,0,1}, invShift[3];
	for(int ii = 0; ii < 3; ++ii)
		invShift[ii] = SBF->avgCenters[HBID][ii];
	HomogeneousMatrixPointerType tmpRot = transformToBarrelCoordinates(tmpAxis, SBF->avgAxes[HBID]);
	TransformPointerType tmpRotTrans = TransformPointerType::New();
	TransformPointerType tmpInvShift = TransformPointerType::New();
	tmpRotTrans->SetMatrix(tmpRot);
	tmpInvShift->Translate(invShift);
	
	barrels[0][HBID]->rotateColumn(tmpRot);
	barrels[0][HBID]->translateColumn(invShift);
	
	tmpInvShift->Concatenate(tmpRotTrans);
	tmpInvShift->Update();
	
	workingSG->setTransformation(tmpInvShift);
	workingSG->applyTransformation();
	
	TransformFilterType piaTmpTrans = TransformFilterType::New();
	piaTmpTrans->SetTransform(tmpInvShift);
	piaTmpTrans->SetInput(pia);
	piaTmpTrans->Update();
	pia->DeepCopy(piaTmpTrans->GetOutput());
	if(wm)
	{
		TransformFilterType wmTmpTrans = TransformFilterType::New();
		wmTmpTrans->SetTransform(tmpInvShift);
		wmTmpTrans->SetInput(wm);
		wmTmpTrans->Update();
		wm->DeepCopy(wmTmpTrans->GetOutput());
	}
};

void Registration::setUpContourCorrection(int mode)
{
	// mode = 1: brightfield + DAB/Cytochrome staining
	if(mode == 1)
	{
		std::vector< double > A1corr, A2corr, A3corr, A4corr;
		std::vector< double > B1corr, B2corr, B3corr, B4corr;
		std::vector< double > C1corr, C2corr, C3corr, C4corr;
		std::vector< double > D1corr, D2corr, D3corr, D4corr;
		std::vector< double > E1corr, E2corr, E3corr, E4corr;
		// absolute values
		A1corr.push_back(0), A1corr.push_back(0), A1corr.push_back(0);
		A2corr.push_back(-215), A2corr.push_back(115), A2corr.push_back(0);
		A3corr.push_back(0), A3corr.push_back(0), A3corr.push_back(0);
		A4corr.push_back(0), A4corr.push_back(0), A4corr.push_back(0);
		B1corr.push_back(-151), B1corr.push_back(45), B1corr.push_back(0);
		B2corr.push_back(-148), B2corr.push_back(68), B2corr.push_back(0);
		B3corr.push_back(-147), B3corr.push_back(113), B3corr.push_back(0);
		B4corr.push_back(-165), B4corr.push_back(146), B4corr.push_back(0);
		C1corr.push_back(-75), C1corr.push_back(93), C1corr.push_back(87);
		C2corr.push_back(-102), C2corr.push_back(74), C2corr.push_back(80);
		C3corr.push_back(-4), C3corr.push_back(19), C3corr.push_back(94);
		C4corr.push_back(-20), C4corr.push_back(143), C4corr.push_back(82);
		D1corr.push_back(-75), D1corr.push_back(54), D1corr.push_back(42);
		D2corr.push_back(-68), D2corr.push_back(9), D2corr.push_back(64);
		D3corr.push_back(-80), D3corr.push_back(36), D3corr.push_back(109);
		D4corr.push_back(-112), D4corr.push_back(55), D4corr.push_back(138);
		E1corr.push_back(37), E1corr.push_back(-12), E1corr.push_back(40);
		E2corr.push_back(38), E2corr.push_back(48), E2corr.push_back(101);
		E3corr.push_back(-77), E3corr.push_back(116), E3corr.push_back(148);
		E4corr.push_back(-90), E4corr.push_back(76), E4corr.push_back(134);
		// relative values
	// 	A1corr.push_back(1), A1corr.push_back(1), A1corr.push_back(1);
	// 	A2corr.push_back(0.685), A2corr.push_back(1.519), A2corr.push_back(1);
	// 	A3corr.push_back(1), A3corr.push_back(1), A3corr.push_back(1);
	// 	A4corr.push_back(1), A4corr.push_back(1), A4corr.push_back(1);
	// 	B1corr.push_back(0.761), B1corr.push_back(1.148), B1corr.push_back(1);
	// 	B2corr.push_back(0.768), B2corr.push_back(1.248), B2corr.push_back(1);
	// 	B3corr.push_back(0.770), B3corr.push_back(1.471), B3corr.push_back(1);
	// 	B4corr.push_back(0.752), B4corr.push_back(1.710), B4corr.push_back(1);
	// 	C1corr.push_back(0.864), C1corr.push_back(1.352), C1corr.push_back(1.099);
	// 	C2corr.push_back(0.829), C2corr.push_back(1.260), C2corr.push_back(1.084);
	// 	C3corr.push_back(0.993), C3corr.push_back(1.057), C3corr.push_back(1.093);
	// 	C4corr.push_back(0.965), C4corr.push_back(1.653), C4corr.push_back(1.079);
	// 	D1corr.push_back(0.846), D1corr.push_back(1.170), D1corr.push_back(1.043);
	// 	D2corr.push_back(0.885), D2corr.push_back(1.025), D2corr.push_back(1.064);
	// 	D3corr.push_back(0.873), D3corr.push_back(1.107), D3corr.push_back(1.107);
	// 	D4corr.push_back(0.832), D4corr.push_back(1.180), D4corr.push_back(1.135);
	// 	E1corr.push_back(1.072), E1corr.push_back(0.966), E1corr.push_back(1.038);
	// 	E2corr.push_back(1.073), E2corr.push_back(1.158), E2corr.push_back(1.093);
	// 	E3corr.push_back(0.844), E3corr.push_back(1.469), E3corr.push_back(1.140);
	// 	E4corr.push_back(0.866), E4corr.push_back(1.283), E4corr.push_back(1.127);
		contourCorrection.insert(std::pair< int, std::vector< double > >(A1, A1corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(A2, A2corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(A3, A3corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(A4, A4corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(B1, B1corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(B2, B2corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(B3, B3corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(B4, B4corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(C1, C1corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(C2, C2corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(C3, C3corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(C4, C4corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(D1, D1corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(D2, D2corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(D3, D3corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(D4, D4corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(E1, E1corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(E2, E2corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(E3, E3corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(E4, E4corr));
	}
	
	// mode = 2: confocal + epifluorescence
	if(mode == 2)
	{
		std::vector< double > A1corr, A2corr, A3corr, A4corr;
		std::vector< double > B1corr, B2corr, B3corr, B4corr;
		std::vector< double > C1corr, C2corr, C3corr, C4corr;
		std::vector< double > D1corr, D2corr, D3corr, D4corr;
		std::vector< double > E1corr, E2corr, E3corr, E4corr;
		std::vector< double > Alphacorr, Betacorr, Gammacorr, Deltacorr;
		// absolute values
		A1corr.push_back(-162), A1corr.push_back(56), A1corr.push_back(63);
		A2corr.push_back(-91), A2corr.push_back(57), A2corr.push_back(88);
		A3corr.push_back(-126), A3corr.push_back(65), A3corr.push_back(99);
		A4corr.push_back(-91), A4corr.push_back(72), A4corr.push_back(147);
		B1corr.push_back(-104), B1corr.push_back(86), B1corr.push_back(25);
		B2corr.push_back(-91), B2corr.push_back(78), B2corr.push_back(87);
		B3corr.push_back(-104), B3corr.push_back(87), B3corr.push_back(88);
		B4corr.push_back(-84), B4corr.push_back(84), B4corr.push_back(74);
		C1corr.push_back(-130), C1corr.push_back(89), C1corr.push_back(64);
		C2corr.push_back(-104), C2corr.push_back(85), C2corr.push_back(109);
		C3corr.push_back(-65), C3corr.push_back(52), C3corr.push_back(150);
		C4corr.push_back(-34), C4corr.push_back(43), C4corr.push_back(181);
		D1corr.push_back(-120), D1corr.push_back(93), D1corr.push_back(32);
		D2corr.push_back(-101), D2corr.push_back(69), D2corr.push_back(79);
		D3corr.push_back(-49), D3corr.push_back(64), D3corr.push_back(88);
		D4corr.push_back(-18), D4corr.push_back(31), D4corr.push_back(139);
		E1corr.push_back(-82), E1corr.push_back(94), E1corr.push_back(30);
		E2corr.push_back(-72), E2corr.push_back(109), E2corr.push_back(105);
		E3corr.push_back(-91), E3corr.push_back(119), E3corr.push_back(126);
		E4corr.push_back(-36), E4corr.push_back(85), E4corr.push_back(157);
		Alphacorr.push_back(-100), Alphacorr.push_back(93), Alphacorr.push_back(-16);
		Betacorr.push_back(-147), Betacorr.push_back(133), Betacorr.push_back(-33);
		Gammacorr.push_back(-145), Gammacorr.push_back(121), Gammacorr.push_back(-17);
		Deltacorr.push_back(-112), Deltacorr.push_back(98), Deltacorr.push_back(-30);
		contourCorrection.insert(std::pair< int, std::vector< double > >(A1, A1corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(A2, A2corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(A3, A3corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(A4, A4corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(B1, B1corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(B2, B2corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(B3, B3corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(B4, B4corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(C1, C1corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(C2, C2corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(C3, C3corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(C4, C4corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(D1, D1corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(D2, D2corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(D3, D3corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(D4, D4corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(E1, E1corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(E2, E2corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(E3, E3corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(E4, E4corr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(Alpha, Alphacorr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(Beta, Betacorr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(Gamma, Gammacorr));
		contourCorrection.insert(std::pair< int, std::vector< double > >(Delta, Deltacorr));
	}
	
	// values for Pia correction [micron]
	correctPiaDist = -25;
};

double * Registration::axisSurfaceIntersection(PolyDataPointerType surface, double * axis, double * center)
{
	PolyDataPointerType surfRoi = selectSurfaceRoi(surface, center, 2000);
	if(surfRoi->GetNumberOfPoints())
	{
		CellLocatorPointerType locator = CellLocatorPointerType::New();
		locator->AutomaticOn();
		locator->SetDataSet(surfRoi);
		locator->BuildLocator();
		double * intersectPt = new double[3];
		double a0[3], a1[3], tol = 0.1, t, pcoords[3];
		int subId;
		vtkIdType cellID;
		GenericCellPointerType intersectCell = GenericCellPointerType::New();
		for(int jj = 0; jj < 3; ++jj)
		{
			a0[jj] = center[jj];
			a1[jj] = center[jj];
		}
		double * direction = new double[3];
		for(int jj = 0; jj < 3; ++jj)
			direction[jj] = axis[jj];
		normalize(direction);
		for(int jj = 0; jj < 3; ++jj)
		{
			a0[jj] += direction[jj]*2000;
			a1[jj] -= direction[jj]*2000;
		}
		delete [] direction;
		int intersection = locator->IntersectWithLine(a0, a1, tol, t, intersectPt, pcoords, subId, cellID, intersectCell);
		if(intersection)
		{
			return intersectPt;
		}
		else
		{
			std::cout << "Warning! Could not find intersection point of axis with surface!" << std::endl;
			return NULL;
		}
	}
	else
	{
		std::flush(std::cout << "Error! Could not select surface ROI for calculation of intersection point." << std::endl);
		return NULL;
	}
};

double * Registration::axisSurfaceIntersection(CellLocatorPointerType surfaceLocator, double * axis, double * center)
{
	double * intersectPt = new double[3];
	double a0[3], a1[3], tol = 1E-03, t, pcoords[3];
	int subId;
	vtkIdType cellID;
	GenericCellPointerType intersectCell = GenericCellPointerType::New();
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] = center[jj];
		a1[jj] = center[jj];
	}
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] += axis[jj]*2000;
		a1[jj] -= axis[jj]*2000;
	}
	int intersection = surfaceLocator->IntersectWithLine(a0, a1, tol, t, intersectPt, pcoords, subId, cellID, intersectCell);
	if(intersection)
	{
		return intersectPt;
	}
	else
	{
		return NULL;
	}
};

double * Registration::axisSurfaceIntersection(CellLocatorPointerType surfaceLocator, double * axis, double * center, vtkIdType& intersectCellID)
{
	double * intersectPt = new double[3];
	double a0[3], a1[3], tol = 1E-03, t, pcoords[3];
	int subId;
	vtkIdType cellID;
	GenericCellPointerType intersectCell = GenericCellPointerType::New();
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] = center[jj];
		a1[jj] = center[jj];
	}
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] += axis[jj]*2000;
		a1[jj] -= axis[jj]*2000;
	}
	int intersection = surfaceLocator->IntersectWithLine(a0, a1, tol, t, intersectPt, pcoords, subId, cellID, intersectCell);
	if(intersection)
	{
		intersectCellID = cellID;
		return intersectPt;
	}
	else
	{
		return NULL;
	}
};

/******************************************************************************/
/*return max z coordinate of reconstruction to ensure that nothing will be    */
/*"sticking out" of the pia                                                   */
/******************************************************************************/
double Registration::morphologyZExtent()
{
	double maxZ = 0;
	if(workingSG)
	{
		maxZ = (*(workingSG->edgesBegin()))->edgePointCoordinates.front()[2];
		std::vector< Edge * >::iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label >= Neuron && (*edgeIt)->label <= Soma)
			{
				std::list< double * >::iterator edgeListIt;
				for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
				{
					double tmpZ = (*edgeListIt)[2];
					if(!zReversed && tmpZ > maxZ)
						maxZ = tmpZ;
					else if(zReversed && tmpZ < maxZ)
						maxZ = tmpZ;
				}
			}
		}
	}
	return maxZ;
};

/******************************************************************************/
/*return max z coordinate of reconstruction to ensure that nothing will be    */
/*"sticking out" of the pia                                                   */
/******************************************************************************/
void Registration::maxDendPoint(double maxPoint[3])
{
	double maxZ = 0, tmpPt[] = {0,0,0};
	if(workingSG)
	{
		maxZ = (*(workingSG->edgesBegin()))->edgePointCoordinates.front()[2];
		std::vector< Edge * >::iterator edgeIt;
		for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label != Dendrite && (*edgeIt)->label != ApicalDendrite && (*edgeIt)->label != BasalDendrite)
				continue;
			
			std::list< double * >::iterator edgeListIt;
			for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{
				double tmpZ = (*edgeListIt)[2];
				if(!zReversed && tmpZ > maxZ)
				{
					maxZ = tmpZ;
					tmpPt[0] = (*edgeListIt)[0];
					tmpPt[1] = (*edgeListIt)[1];
					tmpPt[2] = (*edgeListIt)[2];
				}
				else if(zReversed && tmpZ < maxZ)
				{
					maxZ = tmpZ;
					tmpPt[0] = (*edgeListIt)[0];
					tmpPt[1] = (*edgeListIt)[1];
					tmpPt[2] = (*edgeListIt)[2];
				}
			}
		}
	}
	maxPoint[0] = tmpPt[0];
	maxPoint[1] = tmpPt[1];
	maxPoint[2] = tmpPt[2];
};

#ifdef PIPELINE_DOC
void Registration::writeTransformSG(TransformPointerType transf, const char * label)
{
// 	std::flush(std::cout << "allocating SpatialGraph" << std::endl);
	AmiraSpatialGraph * outSG = new AmiraSpatialGraph;
// 	std::flush(std::cout << "setting SpatialGraph" << std::endl);
	outSG->mergeSpatialGraph(workingSG);
	outSG->setTransformation(transf);
// 	std::flush(std::cout << "applying transformation to SpatialGraph" << std::endl);
	outSG->applyTransformation();
	std::string ofName(pipelineDocName);
	ofName += label;
// 	std::flush(std::cout << "allocating SpatialGraph writer" << std::endl);
	Reader * transWriter = new Reader(ofName.c_str(), ofName.c_str());
// 	std::flush(std::cout << "setting SpatialGraph writer data" << std::endl);
	transWriter->setSpatialGraph(outSG);
// 	std::flush(std::cout << "writing transformed SpatialGraph" << std::endl);
	transWriter->writeSpatialGraphFile();
	delete transWriter, delete outSG;
};

void Registration::writeRegStepParameters(const char * label)
{
	std::string ofName(pipelineDocName);
	ofName += label;
	ofName += ".csv";
	std::ofstream parameterFile(ofName.c_str());
	if(workingSG && SBF->avgPiaSurface)
	{
		double piaDist = registeredPiaDistance();
		double somaAxisDist = somaDistanceToColumnAxis();
		double dendDist = registeredPiaDendriteDistance();
		int colSepFlag = getColumnSeptumFlag();
		parameterFile << "Pia-soma distance:\t" << piaDist << std::endl;
		parameterFile << "Soma distance to column axis:\t" << somaAxisDist << std::endl;
		parameterFile << "Minimum Pia-dendrite distance:\t" << dendDist << std::endl;
		parameterFile << "Column (1)/ Septum (0):\t" << colSepFlag << std::endl;
		parameterFile.close();
	}
	else
	{
		parameterFile << "Error! Internal data incomplete/inaccessible!" << std::endl;
		parameterFile.close();
	}
};
#endif

#ifdef REG_ACCURACY
void Registration::regVariability()
{
	//compute coordinates in cylindrical coordinates
	//measured in home column
	//phi measured relative to neighboring barrel along row
	double rPos, phiPos, zPos;
	double somaPt[3], radialSomaPt[3];
	getPCenterOfStructure(workingSG, Soma, somaPt);
	int HBID = workingSG->getHomeBarrel();
	int NBID = neighborBarrel[HBID];
	
	zPos = registeredPiaDistance();
	
	double t, closestPt[3], zAxis[3], xAxis[3], yAxis[3];
	rPos = vtkLine::DistanceToLine(somaPt, avgColumns[HBID]->top, avgColumns[HBID]->bottom, t, closestPt);
	rPos = sqrt(rPos);
	
	for(int ii = 0; ii < 3; ++ii)
	{
		radialSomaPt[ii] = somaPt[ii] - closestPt[ii];
		zAxis[ii] = avgColumns[HBID]->top[ii] - avgColumns[HBID]->bottom[ii];
		xAxis[ii] = avgCenters[NBID][ii] - avgCenters[HBID][ii];
	}
	normalize(zAxis);
	double correction = vtkMath::Dot(zAxis, xAxis);
	for(int ii = 0; ii < 3; ++ii)
		xAxis[ii] -= correction*zAxis[ii];
	normalize(xAxis);
	vtkMath::Cross(zAxis, xAxis, yAxis);
	double xNew = vtkMath::Dot(radialSomaPt, xAxis);
	double yNew = vtkMath::Dot(radialSomaPt, yAxis);
	phiPos = std::atan2(yNew, xNew) + PI;
	
	std::string somaVarFilename(logFileName);
	somaVarFilename += "_soma_var.csv";
	std::ofstream somaVarFile;
	somaVarFile.open(somaVarFilename.c_str());
	somaVarFile << "# Soma location in global cartesian and local cylindrical coordinates" << std::endl;
	somaVarFile << "# var_alpha = " << var_alpha << std::endl;
	somaVarFile << "# var_gamma = " << var_gamma << std::endl;
	somaVarFile << "# dphi = " << angleUncertainty << std::endl;
	somaVarFile << "x y z\t" << somaPt[0] << "\t" << somaPt[1] << "\t" << somaPt[2] << std::endl;
	somaVarFile << "r phi z\t" << rPos << "\t" << phiPos << "\t" << zPos << std::endl;
	somaVarFile.close();
	
	//compute total axon and dendrite length
	double axonTotal = 0, dendriteTotal = 0, apicalTotal = 0, basalTotal = 0;
	std::vector< Edge * >::const_iterator edgeIt;
	for(edgeIt = workingSG->edgesBegin(); edgeIt != workingSG->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label == Axon)
			axonTotal += (*edgeIt)->segmentLength();
		else if((*edgeIt)->label == Dendrite)
			dendriteTotal += (*edgeIt)->segmentLength();
		else if((*edgeIt)->label == ApicalDendrite)
			apicalTotal += (*edgeIt)->segmentLength();
		else if((*edgeIt)->label == BasalDendrite)
			basalTotal += (*edgeIt)->segmentLength();
	}
	
	std::string lengthVarFilename(logFileName);
	lengthVarFilename += "_length_var.csv";
	std::ofstream lengthVarFile;
	lengthVarFile.open(lengthVarFilename.c_str());
	lengthVarFile << "# Axon/Dendrite/Apical/Basal length" << std::endl;
	lengthVarFile << "# var_alpha = " << var_alpha << std::endl;
	lengthVarFile << "# var_gamma = " << var_gamma << std::endl;
	lengthVarFile << axonTotal << "\t" << dendriteTotal << "\t" << apicalTotal << "\t" << basalTotal << std::endl;
	lengthVarFile.close();
};
#endif



