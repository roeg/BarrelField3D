#include "inputcheckpoint.h"
// #define DEBUG

InputCheckpoint::InputCheckpoint()
{
	sg = NULL;
	inputOK = 0;
}

InputCheckpoint::InputCheckpoint ( AmiraSpatialGraph* inputSG )
{
	sg = new AmiraSpatialGraph;
	sg->mergeSpatialGraph(inputSG);
	inputOK = 1;
	piaFlag = 0, wmFlag = 0, zReversed = 0;
	piaSpacing = 0, wmSpacing = 0;
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
}

InputCheckpoint::~InputCheckpoint()
{
	if(sg) delete sg;
}

void InputCheckpoint::run()
{
	#ifdef DEBUG
	std::cout << "Checking input SpatialGraph!" << std::endl;
	#endif
	detectSectionThickness();
	detectInputZDirection();
}

void InputCheckpoint::checkBarrelField()
{
	if(sg)
	{
		#ifdef DEBUG
		std::cout << "Checking input SpatialGraph!" << std::endl;
		#endif
		std::list< int >::const_iterator labelIt;
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			PolyDataPointerType barrel = PolyDataPointerType::New();
			if(sg->extractLandmark(ID, barrel))
				if(barrel->GetNumberOfCells() > 2)
				{
					std::cout << "Error! Barrel/barrel column " << ID << " has wrong format.";
					std::cout << "At most two contours allowed for each barrel." << std::endl;
					inputOK = 0;
					break;
				}
		}
		#ifdef DEBUG
		std::cout << "inputOK = " << inputOK << std::endl;
		#endif
	}
}

void InputCheckpoint::checkNeuronMorphology()
{
	if(sg)
	{
		std::list< int > inputLabels;
		for(int ii = Neuron; ii <= ZAxis; ++ii)
			if(sg->isLabelInSpatialGraph(ii))
				inputLabels.push_back(ii);
		if(std::find(inputLabels.begin(), inputLabels.end(), Soma) != inputLabels.end())
		{
			std::cout << "Soma found in input hoc file!" << std::endl;
			somaFlag = 1;
		}
		if(std::find(inputLabels.begin(), inputLabels.end(), Dendrite) != inputLabels.end())
		{
			std::cout << "Dendrite found in input hoc file!" << std::endl;
			dendriteFlag = 1;
		}
		if(std::find(inputLabels.begin(), inputLabels.end(), ApicalDendrite) != inputLabels.end())
		{
			std::cout << "Apical dendrite found in input hoc file!" << std::endl;
			apicalFlag = 1;
		}
		if(std::find(inputLabels.begin(), inputLabels.end(), BasalDendrite) != inputLabels.end())
		{
			std::cout << "Basal dendrite found in input hoc file!" << std::endl;
			basalFlag = 1;
		}
		if(std::find(inputLabels.begin(), inputLabels.end(), Axon) != inputLabels.end())
		{
			std::cout << "Axon found in input hoc file!" << std::endl;
			axonFlag = 1;
		}
	}
}

void InputCheckpoint::detectSectionThickness()
{
	if(sg)
	{
		PolyDataPointerType pia = PolyDataPointerType::New();
		if(sg->extractLandmark(Pia, pia))
		{
			std::list< double > spacingList;
			std::list< double > zIndexList;
			for(int ii = 0; ii < pia->GetNumberOfCells(); ++ii)
			{
				double currZIndex = pia->GetCell(ii)->GetBounds()[4];
// 				std::cout << "currZIndex = " << currZIndex << std::endl;
				zIndexList.push_back(round(currZIndex));
			}
			zIndexList.sort();
			std::list< double >::const_iterator zIndexListIt1 = zIndexList.begin();
			std::list< double >::const_iterator zIndexListIt2;
			++zIndexListIt1;
			for(zIndexListIt2 = zIndexList.begin() ; zIndexListIt1 != zIndexList.end(); ++zIndexListIt1, ++zIndexListIt2)
				spacingList.push_back(round(*zIndexListIt1 - *zIndexListIt2));
			spacingList.sort();
			spacingList.unique();
			if(spacingList.size() > 1)
			{
				inputOK = 0;
				std::cout << "Warning! Input SpatialGraph has non-uniform pia section thickness." << std::endl;
// 				std::cout << "Run 'check_hoc_file.py' on input .hoc file first." << std::endl;
// 				std::cout << "Output invalid!" << std::endl;
			}
			else
			{
				piaFlag = 1;
				piaSpacing = spacingList.front();
				#ifdef DEBUG
				std::cout << "Pia z spacing = " << piaSpacing << "um" << std::endl;
				#endif
// 				alignSpatialGraphGlobalZ(zIndexList);
			}
		}
		
		PolyDataPointerType wm = PolyDataPointerType::New();
		if(sg->extractLandmark(WhiteMatter, wm))
		{
			std::list< double > spacingList;
			std::list< double > zIndexList;
			for(int ii = 0; ii < wm->GetNumberOfCells(); ++ii)
			{
				double currZIndex = wm->GetCell(ii)->GetBounds()[4];
// 				std::cout << "currZIndex = " << currZIndex << std::endl;
				zIndexList.push_back(round(currZIndex));
			}
			zIndexList.sort();
			std::list< double >::const_iterator zIndexListIt1 = zIndexList.begin();
			std::list< double >::const_iterator zIndexListIt2;
			++zIndexListIt1;
			for(zIndexListIt2 = zIndexList.begin() ; zIndexListIt1 != zIndexList.end(); ++zIndexListIt1, ++zIndexListIt2)
				spacingList.push_back(round(*zIndexListIt1 - *zIndexListIt2));
			spacingList.sort();
			spacingList.unique();
			if(spacingList.size() > 1)
			{
				inputOK = 0;
				std::cout << "Warning! Input SpatialGraph has non-uniform WM section thickness." << std::endl;
// 				std::cout << "Run 'check_hoc_file.py' on input .hoc file first." << std::endl;
// 				std::cout << "Output invalid!" << std::endl;
			}
			else
			{
				wmFlag = 1;
				wmSpacing = spacingList.front();
				#ifdef DEBUG
				std::cout << "WM z spacing = " << wmSpacing << "um" << std::endl;
				#endif
// 				alignSpatialGraphGlobalZ(zIndexList);
			}
		}
	}
	else
		std::cout << "Error! Input SpatialGraph is empty" << std::endl;
}

void InputCheckpoint::detectInputZDirection()
{
	if(sg)
	{
		if(!sg->isLabelInSpatialGraph(Pia))
		{
			std::cout << "Pia not found in Input SpatialGraph. Could not determine z-direction." << std::endl;
			std::cout << "Setting zReversed = 0" << std::endl;
			zReversed = 0;
			return;
		}
		
		PolyDataPointerType landmark = PolyDataPointerType::New();
		if(sg->extractLandmark(Pia, landmark))
		{
			double minArea, maxArea = landmark->GetCell(0)->GetLength2();
			minArea = maxArea;
			int minID = 0, maxID = 0;
			for(int ii = 1; ii < landmark->GetNumberOfCells(); ++ii)
			{
// 				landmark->GetCell(ii)->Print(std::cout);
				double tmpArea = landmark->GetCell(ii)->GetLength2();
// 				std::cout << "tmpArea = " << tmpArea << std::endl;
				if(tmpArea > maxArea)
				{
					maxArea = tmpArea;
					maxID = ii;
				}
				if(tmpArea < minArea)
				{
					minArea = tmpArea;
					minID = ii;
				}
			}
			double zMinID = landmark->GetCell(minID)->GetBounds()[4], zMaxID = landmark->GetCell(maxID)->GetBounds()[4];
			zReversed = zMinID > zMaxID ? 0 : 1;
			#ifdef DEBUG
			std::cout << "zReversed = " << zReversed << std::endl;
			#endif
		}
	}
	else
		std::cout << "Error! Input SpatialGraph is empty" << std::endl;
}

// void InputCheckpoint::getLandmarkMinMaxIDs ( PolyDataPointerType landmark, int& minID, int& maxID )
// {
// 
// }

InputParameters InputCheckpoint::getParameters()
{
	InputParameters parameters;
	if(inputOK)
	{
		parameters.piaFlag = this->piaFlag;
		parameters.wmFlag = this->wmFlag;
		parameters.zReversed = this->zReversed;
		parameters.piaSpacing = this->piaSpacing;
		parameters.wmSpacing = this->wmSpacing;
		
		parameters.somaFlag = this->somaFlag;
		parameters.dendriteFlag = this->dendriteFlag;
		parameters.apicalFlag = this->apicalFlag;
		parameters.basalFlag = this->basalFlag;
		parameters.axonFlag = this->axonFlag;
	}
	return parameters;
}
