#include "landmarkanalyzer.h"
// #define DEBUG

#ifdef L1
#undef L1
#endif
#define L1 1
#define L23 2
#define L4 3
#define L5 4
#define L6 5

LandmarkAnalyzer::LandmarkAnalyzer()
{
	piaFlag = 0;
	wmFlag = 0;
	barrelFlag = 0;
	barreloidFlag = 0;
	cellFlag = 0;
	barrelLayerFlag = 0;
	supraFlag = 0;
	granFlag = 0;
	infraFlag = 0;
	SBF = new BarrelField(false);
	SPACING = 10;	// voxel side length in um, CHANGE ONLY HERE!!!
}

LandmarkAnalyzer::~LandmarkAnalyzer()
{
	if(SBF) delete SBF;
	if(columnSomaProfiles.size())
	{
		std::map< int, Profile * >::iterator colProfileIt;
		for(colProfileIt = columnSomaProfiles.begin(); colProfileIt != columnSomaProfiles.end(); ++colProfileIt)
			delete colProfileIt->second;
	}
}

void LandmarkAnalyzer::setPiaSurface ( Surface* piaSurface )
{
	pia = piaSurface;
	if(pia) piaFlag = 1;
}

void LandmarkAnalyzer::setWMSurface ( Surface* wmSurface )
{
	WM = wmSurface;
	if(WM) wmFlag = 1;
}

void LandmarkAnalyzer::setBarrelField ( std::map< int, Column* > barrelField )
{
	this->barrelField = barrelField;
	if(this->barrelField.size()) barrelFlag = 1;
}

void LandmarkAnalyzer::setBarrelFieldLayer ( std::map< int, Column* > barrelFieldLayer )
{
	this->barrelFieldLayer = barrelFieldLayer;
	if(this->barrelFieldLayer.size()) barrelLayerFlag = 1;
}

void LandmarkAnalyzer::setSupragranularLayer(std::map< int, Column* > barrelField)
{
	this->supragranularLayer = barrelField;
	if(this->supragranularLayer.size()) supraFlag = 1;
}

void LandmarkAnalyzer::setGranularLayer(std::map< int, Column* > barrelField)
{
	this->granularLayer = barrelField;
	if(this->granularLayer.size()) granFlag = 1;
}

void LandmarkAnalyzer::setInfragranularLayer(std::map< int, Column* > barrelField)
{
	this->infragranularLayer = barrelField;
	if(this->infragranularLayer.size()) infraFlag = 1;
}

void LandmarkAnalyzer::setLandmarkSet ( PointsPointerType cellLandmarks )
{
	cellSomata = cellLandmarks;
	if(cellSomata->GetNumberOfPoints()) cellFlag = 1;
}

void LandmarkAnalyzer::setBarreloidField ( std::map< int, ClosedSurface* > barreloidField )
{
	this->barreloidField = barreloidField;
	if(this->barreloidField.size()) barreloidFlag = 1;
}

void LandmarkAnalyzer::standardBFAnalysis()
{
	this->barrelField = SBF->avgColumns;
	barrelFlag = 1;
	this->pia = SBF->avgPiaSurface;
	piaFlag = 1;
	this->WM = SBF->avgWMSurface;
	wmFlag = 1;
}

// compute profiles for each column; in ambiguous cases, assign soma
// to nearest column (takes overlap into account)
void LandmarkAnalyzer::computeColumnProfiles(const char* outputFilename)
{
	if(/*!piaFlag || !wmFlag ||*/ !barrelFlag || !cellFlag)
	{
		std::cout << "Error! Cannot compute column profiles. Incomplete input data!" << std::endl;
		return;
	}
	outFilenameGlobal = outputFilename;
	
	std::cout << "Computing unambiguous profiles..." << std::endl;
	if(columnSomaProfiles.size()) columnSomaProfiles.clear();
	
	// create columns from barrel contours
	// and intersection points with Pia/WM
	// allocate data structures for storing
	// somata inside columns
	std::map< int, PointsPointerType > somaColumns;
	std::map< int, Column * > barrelColumns;
	std::map< int, std::vector< std::vector< double > > > radialContours;
	std::map< int, double > minDistances;
	std::map< int, double > maxDistances;
	std::map< int, double > voxelNumbers;
	
	createCountDataStructures(somaColumns, barrelColumns, radialContours, minDistances, maxDistances);
	
	assignSomataToColumns(somaColumns, barrelColumns, radialContours, minDistances, maxDistances);
		
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrelField.find(ID) != barrelField.end())
		{
			Profile * zProfile = computeColumnSomaProfile(somaColumns[ID], barrelColumns[ID]->top, barrelColumns[ID]->bottom);
			columnSomaProfiles.insert(std::pair< int, Profile * >(ID, zProfile));
			Profile * voxelProfile = computeColumnVolumeProfile(barrelColumns, minDistances, maxDistances, radialContours, ID);
			voxelNumbers.insert(std::pair< int, double >(ID, voxelProfile->getIntegral()));
			Profile * densityProfile = new Profile(50);
			for(unsigned int ii = 0; ii < zProfile->getProfile()->size() && ii < voxelProfile->getProfile()->size(); ++ii)
			{
				double voxelVol = SPACING*SPACING*SPACING*1E-9;
				if(voxelProfile->getProfile()->at(ii))
				{
					double binDens = zProfile->getProfile()->at(ii)/(voxelProfile->getProfile()->at(ii)*voxelVol);
					densityProfile->addSegment(binDens, ii);
				}
				else
					densityProfile->addSegment(0, ii);
			}
			
			std::string colVoxelName(outputFilename);
			colVoxelName += "_%03d";
			colVoxelName += "micron";
			colVoxelName += "_voxels_";
			colVoxelName += SBF->int2Labels[ID];
			char * colVoxelStr = new char[256];
			std::sprintf(colVoxelStr, colVoxelName.c_str(), (int)SPACING);
			writeZProfile(voxelProfile, colVoxelStr);
			delete voxelProfile;
			
			std::string colDensName(outputFilename);
			colDensName += "_%03d";
			colDensName += "micron";
			colDensName += "_density_";
			colDensName += SBF->int2Labels[ID];
			char * colDensStr = new char[256];
			std::sprintf(colDensStr, colDensName.c_str(), (int)SPACING);
			writeZProfile(densityProfile, colDensStr);
			delete densityProfile;
		}
	}
	
	// write separate landmark files for each column
	#ifdef DEBUG
	std::cout << "Writing landmark files for all soma columns..." << std::endl;
	#endif
	std::map< int, PointsPointerType >::const_iterator somaColIt;
	for(somaColIt = somaColumns.begin(); somaColIt != somaColumns.end(); ++somaColIt)
	{
		int ID = somaColIt->first;
		std::string ofName(outputFilename);
		ofName += "_somata_";
		ofName += SBF->int2Labels[ID];
		Reader * somaColWriter = new Reader(ofName.c_str(), ofName.c_str());
		somaColWriter->writeLandmarkFile(somaColIt->second);
		delete somaColWriter;
		delete barrelColumns[ID];
	}
	
	writeColumnSomaProfiles(outputFilename);
	
	std::string ofName(outputFilename);
	ofName += "_all_somata.csv";
	std::ofstream SomaList;
	SomaList.open(ofName.c_str());
	SomaList << "# Column\tNr of somata" << std::endl;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		SomaList << SBF->int2Labels[ID];
		if(somaColumns.find(ID) != somaColumns.end())
		{
			SomaList << "\t" << somaColumns[ID]->GetNumberOfPoints();
		}
		SomaList << std::endl;
	}
	SomaList.close();
	
	std::string ofName2(outputFilename);
	ofName2 += "_all_volumes.csv";
	std::ofstream VoxelList;
	VoxelList.open(ofName2.c_str());
	VoxelList << "# Column\tVolume [mm^3]" << std::endl;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		VoxelList << SBF->int2Labels[ID];
		if(voxelNumbers.find(ID) != voxelNumbers.end())
		{
			double vol = voxelNumbers[ID]*SPACING*SPACING*SPACING*1E-9;
			VoxelList << "\t" << vol;
		}
		VoxelList << std::endl;
	}
	VoxelList.close();
}

void LandmarkAnalyzer::computeCorrectedColumnProfiles(const char* outputFilename)
{
	if(/*!piaFlag || !wmFlag ||*/ !barrelFlag || !cellFlag)
	{
		std::cout << "Error! Cannot compute column profiles. Incomplete input data!" << std::endl;
		return;
	}
	outFilenameGlobal = outputFilename;
	
	std::cout << "Computing unambiguous profiles..." << std::endl;
	if(columnSomaProfiles.size()) columnSomaProfiles.clear();
	
	// create columns from barrel contours
	// and intersection points with Pia/WM
	// allocate data structures for storing
	// somata inside columns
	std::map< int, PointsPointerType > somaColumns;
	std::map< int, Column * > barrelColumns;
	std::map< int, std::vector< std::vector< double > > > radialContours;
	std::map< int, double > minDistances;
	std::map< int, double > maxDistances;
	
	createCountDataStructures(somaColumns, barrelColumns, radialContours, minDistances, maxDistances);
	
	assignSomataToColumns(somaColumns, barrelColumns, radialContours, minDistances, maxDistances);
	
	// raw column soma profiles
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrelField.find(ID) != barrelField.end())
		{
#ifdef DEBUG
std::cout << "computing soma z profile for column " << ID << std::endl;
#endif
			Profile * zProfile = computeColumnSomaProfile(somaColumns[ID], barrelColumns[ID]->top, barrelColumns[ID]->bottom);
			columnSomaProfiles.insert(std::pair< int, Profile * >(ID, zProfile));
		}
	}
	writeColumnSomaProfiles(outputFilename);
		
	correctSomaProfiles(columnSomaProfiles);
	std::string correctedOutName(outputFilename);
	correctedOutName += "_corrected";
	writeColumnSomaProfiles(correctedOutName.c_str());
	
	// write separate landmark files for each column
	#ifdef DEBUG
	std::cout << "Writing landmark files for all soma columns..." << std::endl;
	#endif
	
	std::map< int, PointsPointerType >::const_iterator somaColIt;
	for(somaColIt = somaColumns.begin(); somaColIt != somaColumns.end(); ++somaColIt)
	{
		int ID = somaColIt->first;
		std::string ofName(outputFilename);
		ofName += "_somata_";
		ofName += SBF->int2Labels[ID];
		Reader * somaColWriter = new Reader(ofName.c_str(), ofName.c_str());
		somaColWriter->writeLandmarkFile(somaColIt->second);
		delete somaColWriter;
		delete barrelColumns[ID];
	}
	
	std::string ofName(outputFilename);
	ofName += "_all_somata.csv";
	std::ofstream SomaList;
	SomaList.open(ofName.c_str());
	SomaList << "# Column\tNr of somata" << std::endl;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		SomaList << SBF->int2Labels[ID];
		if(columnSomaProfiles.find(ID) != columnSomaProfiles.end())
		{
			SomaList << "\t" << columnSomaProfiles[ID]->getIntegral();
		}
		SomaList << std::endl;
	}
	SomaList.close();
}

void LandmarkAnalyzer::computeCorrectedColumnProfilesInLayer(const char* outputFilename)
{
	if(!barrelFlag || !barrelLayerFlag || !cellFlag)
	{
		std::cout << "Error! Cannot compute column profiles. Incomplete input data!" << std::endl;
		return;
	}
	outFilenameGlobal = outputFilename;
	
	std::cout << "Computing unambiguous profiles..." << std::endl;
	if(columnSomaProfiles.size()) columnSomaProfiles.clear();
	
	// create columns from barrel contours
	// and intersection points with Pia/WM
	// allocate data structures for storing
	// somata inside columns
	std::map< int, PointsPointerType > somaColumns;
	std::map< int, Column * > layerColumns;
	std::map< int, std::vector< std::vector< double > > > radialContours;
	std::map< int, double > minDistances;
	std::map< int, double > maxDistances;
	
	createLayerCountDataStructures(somaColumns, layerColumns, radialContours, minDistances, maxDistances);
	
	assignSomataToColumns(somaColumns, layerColumns, radialContours, minDistances, maxDistances);
	
	// raw column soma profiles
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrelField.find(ID) != barrelField.end())
		{
#ifdef DEBUG
std::cout << "computing soma z profile for column " << ID << std::endl;
#endif
			Profile * zProfile = computeColumnSomaProfile(somaColumns[ID], barrelField[ID]->top, barrelField[ID]->bottom);
			columnSomaProfiles.insert(std::pair< int, Profile * >(ID, zProfile));
		}
	}
	writeColumnSomaProfiles(outputFilename);
		
	correctSomaProfiles(columnSomaProfiles);
	std::string correctedOutName(outputFilename);
	correctedOutName += "_corrected";
	writeColumnSomaProfiles(correctedOutName.c_str());
	
	// write separate landmark files for each column
	#ifdef DEBUG
	std::cout << "Writing landmark files for all soma columns..." << std::endl;
	#endif
	
	std::map< int, PointsPointerType >::const_iterator somaColIt;
	for(somaColIt = somaColumns.begin(); somaColIt != somaColumns.end(); ++somaColIt)
	{
		int ID = somaColIt->first;
		std::string ofName(outputFilename);
		ofName += "_somata_";
		ofName += SBF->int2Labels[ID];
		Reader * somaColWriter = new Reader(ofName.c_str(), ofName.c_str());
		somaColWriter->writeLandmarkFile(somaColIt->second);
		delete somaColWriter;
		delete layerColumns[ID];
	}
	
	std::string ofName(outputFilename);
	ofName += "_all_somata.csv";
	std::ofstream SomaList;
	SomaList.open(ofName.c_str());
	SomaList << "# Column\tNr of somata" << std::endl;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		SomaList << SBF->int2Labels[ID];
		if(columnSomaProfiles.find(ID) != columnSomaProfiles.end())
		{
			SomaList << "\t" << columnSomaProfiles[ID]->getIntegral();
		}
		SomaList << std::endl;
	}
	SomaList.close();
}

void LandmarkAnalyzer::countLaminarNeuronNumbers(const char* outputFilename, bool IN)
{
	if(/*!piaFlag || !wmFlag ||*/ !barrelFlag || !cellFlag)
	{
		std::cout << "Error! Cannot compute column profiles. Incomplete input data!" << std::endl;
		return;
	}
	outFilenameGlobal = outputFilename;
	
	std::cout << "Computing unambiguous profiles..." << std::endl;
	if(columnSomaProfiles.size()) columnSomaProfiles.clear();
	
	std::list< int >::const_iterator labelIt;
	// create columns from barrel contours
	// and intersection points with Pia/WM
	// allocate data structures for storing
	// somata inside columns
	std::map< int, PointsPointerType > somaColumns;
	std::map< int, Column * > barrelColumns;
	std::map< int, std::vector< std::vector< double > > > radialContours;
	std::map< int, double > minDistances;
	std::map< int, double > maxDistances;
	std::map< int, double > voxelNumbers;
	
	createCountDataStructures(somaColumns, barrelColumns, radialContours, minDistances, maxDistances);
	
	// set up data structures for all layers in all columns
	// layer depths Oberlaender dendrite column
	std::map< int, std::string > layerNamesMO2011;
	std::map< int, std::map< int, PointsPointerType > > layerNumbersMO2011;
	std::map< int, std::map< int, Profile* > > layerProfilesCorrectedMO2011;
	layerNamesMO2011.insert(std::pair< int, std::string >(1, std::string("L1")));
	layerNamesMO2011.insert(std::pair< int, std::string >(2, std::string("L2_3")));
	layerNamesMO2011.insert(std::pair< int, std::string >(3, std::string("L4")));
	layerNamesMO2011.insert(std::pair< int, std::string >(4, std::string("L5")));
	layerNamesMO2011.insert(std::pair< int, std::string >(5, std::string("L6")));
	
	std::map< int, std::map< int, PointsPointerType > >::iterator layerIt;
	std::map< int, std::string >::const_iterator layerNameIt;
	for(layerNameIt = layerNamesMO2011.begin(); layerNameIt != layerNamesMO2011.end(); ++layerNameIt)
	{
		int layer = layerNameIt->first;
		std::map< int, PointsPointerType > layerSomata;
		layerNumbersMO2011.insert(std::pair< int, std::map< int, PointsPointerType > >(layer, layerSomata));
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			PointsPointerType somata = PointsPointerType::New();
			somata->Allocate(1);
			somata->SetDataTypeToFloat();
			layerNumbersMO2011[layer].insert(std::pair< int, PointsPointerType >(ID, somata));
		}
		PointsPointerType septumSomata = PointsPointerType::New();
		septumSomata->Allocate(1);
		septumSomata->SetDataTypeToFloat();
		layerNumbersMO2011[layer].insert(std::pair< int, PointsPointerType >(Septum, septumSomata));
	}
	
	// now process every cell soma
	// and assign landmarks to nearest
	// column if in doubt
	#ifdef DEBUG
	std::cout << "Checking " << cellSomata->GetNumberOfPoints() << " somata!" << std::endl;
	#endif
	for(unsigned int ii = 0; ii < cellSomata->GetNumberOfPoints(); ++ii)
	{
		double somaLocation[3];
		cellSomata->GetPoint(ii, somaLocation);
		
		// check distance to all columns
		double depth = 1E09;
		double minDistance = 1E09;
		double outsideMinDistance = 1E09;
		int closestColumn = 0;
		int outsideClosestColumn = 0;
		std::map< int, Column * >::const_iterator columnIt;
		for(columnIt = barrelColumns.begin(); columnIt != barrelColumns.end(); ++columnIt)
		{
			int ID = columnIt->first;
			double t, projectedPt[3];
			double dist = vtkLine::DistanceToLine(somaLocation, barrelColumns[ID]->top, barrelColumns[ID]->bottom, t, projectedPt);
			dist = sqrt(dist);
			if(/*dist > maxDistances[ID] ||*/ t < 0 || t > 1)
				continue;
			if(dist < minDistances[ID])
			{
				if(dist < minDistance)
				{
					minDistance = dist;
					closestColumn = ID;
					depth = t;
				}
				continue;
			}
			if(dist > maxDistances[ID])
			{
				if(dist < outsideMinDistance)
				{
					outsideMinDistance = dist;
					outsideClosestColumn = ID;
					depth = t;
				}
				continue;
			}
			
			PolyDataPointerType polyData = PolyDataPointerType::New();
			polyData->Allocate(1);
			PointsPointerType points = PointsPointerType::New();
			points->SetDataTypeToFloat();
			PolygonPointerType poly = PolygonPointerType::New();
			poly->GetPointIds()->SetNumberOfIds(radialContours[ID].size());
			for(int jj = 0; jj < radialContours[ID].size(); ++jj)
			{
				double tmp[3];
				for(int kk = 0; kk < 3; ++kk)
					tmp[kk] = projectedPt[kk] + radialContours[ID][jj][kk];
				points->InsertNextPoint(tmp);
				poly->GetPointIds()->SetId(jj, jj);
			}
			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
			polyData->SetPoints(points);
			polyData->Update();
			double closestPoint[3], pCoords[3], dist2;
			int subId;
			double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
			int insidePolygon = polyData->GetCell(0)->EvaluatePosition(somaLocation, closestPoint, subId, pCoords, dist2, weights);
			if(insidePolygon == 1)
			{
				if(dist < minDistance)
				{
					minDistance = dist;
					closestColumn = ID;
					depth = t;
				}
			}
			delete [] weights;
		}
		// assign landmark to closest soma column
		if(closestColumn)
		{
#ifdef DEBUG
			std::cout << "found closest column: " << closestColumn << std::endl;
#endif
			int layer = computeLaminarPosition(somaLocation, closestColumn);
#ifdef DEBUG
			std::cout << "found closest layer MO: " << layer << std::endl;
#endif
			if(layer == -1 || layer == 0)
			{
				std::cout << "Warning! Could not find laminar position for point [" << somaLocation[0] << "," << somaLocation[1] << "," << somaLocation[2] << "] ";
				std::cout << "in column " << SBF->int2Labels[closestColumn] << std::endl;
				continue;
			}
			int oldNrOfPoints = layerNumbersMO2011[layer][closestColumn]->GetNumberOfPoints();
			layerNumbersMO2011[layer][closestColumn]->InsertPoint(oldNrOfPoints, somaLocation);
		}
		else if(outsideClosestColumn)
		{
#ifdef DEBUG
			std::cout << "found closest column: " << closestColumn << std::endl;
#endif
			int layer = computeLaminarPosition(somaLocation, outsideClosestColumn);
#ifdef DEBUG
			std::cout << "found closest layer MO: " << layer << std::endl;
#endif
			if(layer == -1 || layer == 0)
			{
				std::cout << "Warning! Could not find laminar position for point [" << somaLocation[0] << "," << somaLocation[1] << "," << somaLocation[2] << "] ";
				std::cout << "outside column " << SBF->int2Labels[outsideClosestColumn] << std::endl;
				continue;
			}
			int oldNrOfPoints = layerNumbersMO2011[layer][Septum]->GetNumberOfPoints();
			layerNumbersMO2011[layer][Septum]->InsertPoint(oldNrOfPoints, somaLocation);
		}
	}
	
	// correction if interneuron counts
	if(IN)
	{
		std::map< int, PointsPointerType >::const_iterator somaColIt;
		for(layerNameIt = layerNamesMO2011.begin(); layerNameIt != layerNamesMO2011.end(); ++layerNameIt)
		{
			int layer = layerNameIt->first;
			for(somaColIt = layerNumbersMO2011[layer].begin(); somaColIt != layerNumbersMO2011[layer].end(); ++somaColIt)
			{
				int ID = somaColIt->first;
				if(barrelField.find(ID) != barrelField.end() && SBF->avgBarrels.find(ID) != SBF->avgBarrels.end())
				{
#ifdef DEBUG
std::cout << "computing soma z profile for column " << ID << std::endl;
#endif
					Profile * zProfile = computeColumnSomaProfile(somaColIt->second, barrelColumns[ID]->top, barrelColumns[ID]->bottom);
					layerProfilesCorrectedMO2011[layer].insert(std::pair< int, Profile * >(ID, zProfile));
				}
			}
			correctSomaProfiles(layerProfilesCorrectedMO2011[layer]);
		}
	}
	
	
	// write separate landmark files for each column
	#ifdef DEBUG
	std::cout << "Writing landmark files for all soma columns..." << std::endl;
	#endif
	std::map< int, PointsPointerType >::const_iterator somaColIt;
	for(layerNameIt = layerNamesMO2011.begin(); layerNameIt != layerNamesMO2011.end(); ++layerNameIt)
	{
		int layer = layerNameIt->first;
		for(somaColIt = layerNumbersMO2011[layer].begin(); somaColIt != layerNumbersMO2011[layer].end(); ++somaColIt)
		{
			int ID = somaColIt->first;
			std::string ofName(outputFilename);
			ofName += "_somata_MO_";
			ofName += layerNameIt->second;
			ofName += "_";
			if(ID)
			{
				ofName += SBF->int2Labels[ID];
			}
			else
			{
				ofName += "Septum";
			}
			Reader * somaColWriter = new Reader(ofName.c_str(), ofName.c_str());
			somaColWriter->writeLandmarkFile(somaColIt->second);
			delete somaColWriter;
			if(layer == 1) delete barrelColumns[ID];
		}
	}
	
	std::string ofNameAll(outputFilename);
	ofNameAll += "_all_somata.csv";
	std::ofstream SomaList;
	SomaList.open(ofNameAll.c_str());
	SomaList << "# Layer\tColumn\tNr of somata" << std::endl;
	SomaList << "# MO 2011 layer borders" << std::endl;
	for(layerNameIt = layerNamesMO2011.begin(); layerNameIt != layerNamesMO2011.end(); ++layerNameIt)
	{
		int layer = layerNameIt->first;
		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			SomaList << layerNameIt->second.c_str() << "\t" << SBF->int2Labels[ID];
			if(layerNumbersMO2011[layer].find(ID) != layerNumbersMO2011[layer].end())
			{
				SomaList << "\t" << layerNumbersMO2011[layer][ID]->GetNumberOfPoints();
			}
			SomaList << std::endl;
		}
		SomaList << layerNameIt->second.c_str() << "\t" << "Septum";
		if(layerNumbersMO2011[layer].find(Septum) != layerNumbersMO2011[layer].end())
		{
			SomaList << "\t" << layerNumbersMO2011[layer][Septum]->GetNumberOfPoints();
		}
		SomaList << std::endl;
	}
	SomaList.close();
	
	if(IN)
	{
		std::string ofNameAllIN(outputFilename);
		ofNameAllIN += "_all_somata_corrected.csv";
		std::ofstream SomaList2;
		SomaList2.open(ofNameAllIN.c_str());
		SomaList2 << "# Layer\tColumn\tNr of somata (corrected)" << std::endl;
		SomaList2 << "# MO 2011 layer borders" << std::endl;
		for(layerNameIt = layerNamesMO2011.begin(); layerNameIt != layerNamesMO2011.end(); ++layerNameIt)
		{
			int layer = layerNameIt->first;
			for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
			{
				int ID = *labelIt;
				SomaList2 << layerNameIt->second.c_str() << "\t" << SBF->int2Labels[ID];
				if(layerProfilesCorrectedMO2011[layer].find(ID) != layerProfilesCorrectedMO2011[layer].end())
				{
					SomaList2 << "\t" << layerProfilesCorrectedMO2011[layer][ID]->getIntegral();
				}
				SomaList2 << std::endl;
			}
		}
		SomaList2.close();
	}
}

// void LandmarkAnalyzer::countLaminarNeuronNumbers(const char* outputFilename, bool IN)
// {
// 	if(/*!piaFlag || !wmFlag ||*/ !barrelFlag || !cellFlag)
// 	{
// 		std::cout << "Error! Cannot compute column profiles. Incomplete input data!" << std::endl;
// 		return;
// 	}
// 	outFilenameGlobal = outputFilename;
// 	
// 	std::cout << "Computing unambiguous profiles..." << std::endl;
// 	if(columnSomaProfiles.size()) columnSomaProfiles.clear();
// 	
// 	std::list< int >::const_iterator labelIt;
// 	// create columns from barrel contours
// 	// and intersection points with Pia/WM
// 	// allocate data structures for storing
// 	// somata inside columns
// 	std::map< int, PointsPointerType > somaColumns;
// 	std::map< int, Column * > barrelColumns;
// 	std::map< int, std::vector< std::vector< double > > > radialContours;
// 	std::map< int, double > minDistances;
// 	std::map< int, double > maxDistances;
// 	std::map< int, double > voxelNumbers;
// 	
// 	createCountDataStructures(somaColumns, barrelColumns, radialContours, minDistances, maxDistances);
// 	
// 	// set up data structures for all layers in all columns
// 	// layer depths Oberlaender dendrite column
// 	std::map< int, std::string > layerNamesMO2011;
// 	std::map< double, int > layerDepthsMO2011;
// 	std::map< int, std::map< int, PointsPointerType > > layerNumbersMO2011;
// 	std::map< int, std::map< int, Profile* > > layerProfilesCorrectedMO2011;
// 	layerNamesMO2011.insert(std::pair< int, std::string >(1, std::string("L1")));
// 	layerNamesMO2011.insert(std::pair< int, std::string >(2, std::string("L2")));
// 	layerNamesMO2011.insert(std::pair< int, std::string >(3, std::string("L3")));
// 	layerNamesMO2011.insert(std::pair< int, std::string >(4, std::string("L4")));
// 	layerNamesMO2011.insert(std::pair< int, std::string >(5, std::string("L5")));
// 	layerNamesMO2011.insert(std::pair< int, std::string >(6, std::string("L6")));
// 	layerDepthsMO2011.insert(std::pair< double, int >(0.08, 1));
// 	layerDepthsMO2011.insert(std::pair< double, int >(0.18, 2));
// 	layerDepthsMO2011.insert(std::pair< double, int >(0.27, 3));
// 	layerDepthsMO2011.insert(std::pair< double, int >(0.51, 4));
// 	layerDepthsMO2011.insert(std::pair< double, int >(0.71, 5));
// 	layerDepthsMO2011.insert(std::pair< double, int >(1.0, 6));
// 	// layer depths Meyer standard column
// 	std::map< int, std::string > layerNamesHSM2010;
// 	std::map< double, int > layerDepthsHSM2010;
// 	std::map< int, std::map< int, PointsPointerType > > layerNumbersHSM2010;
// 	std::map< int, std::map< int, Profile* > > layerProfilesCorrectedHSM2010;
// 	layerNamesHSM2010.insert(std::pair< int, std::string >(1, std::string("L1")));
// 	layerNamesHSM2010.insert(std::pair< int, std::string >(2, std::string("L2")));
// 	layerNamesHSM2010.insert(std::pair< int, std::string >(3, std::string("L3")));
// 	layerNamesHSM2010.insert(std::pair< int, std::string >(4, std::string("L4")));
// 	layerNamesHSM2010.insert(std::pair< int, std::string >(5, std::string("L5A")));
// 	layerNamesHSM2010.insert(std::pair< int, std::string >(6, std::string("L5B")));
// 	layerNamesHSM2010.insert(std::pair< int, std::string >(7, std::string("L6A")));
// 	layerNamesHSM2010.insert(std::pair< int, std::string >(8, std::string("L6B")));
// 	layerDepthsHSM2010.insert(std::pair< double, int >(0.04, 1));
// 	layerDepthsHSM2010.insert(std::pair< double, int >(0.13, 2));
// 	layerDepthsHSM2010.insert(std::pair< double, int >(0.26, 3));
// 	layerDepthsHSM2010.insert(std::pair< double, int >(0.40, 4));
// 	layerDepthsHSM2010.insert(std::pair< double, int >(0.50, 5));
// 	layerDepthsHSM2010.insert(std::pair< double, int >(0.66, 6));
// 	layerDepthsHSM2010.insert(std::pair< double, int >(0.81, 7));
// 	layerDepthsHSM2010.insert(std::pair< double, int >(1.0, 8));
// 	
// 	std::map< int, std::map< int, PointsPointerType > >::iterator layerIt;
// 	std::map< int, std::string >::const_iterator layerNameIt;
// 	for(layerNameIt = layerNamesMO2011.begin(); layerNameIt != layerNamesMO2011.end(); ++layerNameIt)
// 	{
// 		int layer = layerNameIt->first;
// 		std::map< int, PointsPointerType > layerSomata;
// 		layerNumbersMO2011.insert(std::pair< int, std::map< int, PointsPointerType > >(layer, layerSomata));
// 		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
// 		{
// 			int ID = *labelIt;
// 			PointsPointerType somata = PointsPointerType::New();
// 			somata->Allocate(1);
// 			somata->SetDataTypeToFloat();
// 			layerNumbersMO2011[layer].insert(std::pair< int, PointsPointerType >(ID, somata));
// 		}
// 	}
// 	for(layerNameIt = layerNamesHSM2010.begin(); layerNameIt != layerNamesHSM2010.end(); ++layerNameIt)
// 	{
// 		int layer = layerNameIt->first;
// 		std::map< int, PointsPointerType > layerSomata;
// 		layerNumbersHSM2010.insert(std::pair< int, std::map< int, PointsPointerType > >(layer, layerSomata));
// 		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
// 		{
// 			int ID = *labelIt;
// 			PointsPointerType somata = PointsPointerType::New();
// 			somata->Allocate(1);
// 			somata->SetDataTypeToFloat();
// 			layerNumbersHSM2010[layer].insert(std::pair< int, PointsPointerType >(ID, somata));
// 		}
// 	}
// 	
// 	// now process every cell soma
// 	// and assign landmarks to nearest
// 	// column if in doubt
// 	#ifdef DEBUG
// 	std::cout << "Checking " << cellSomata->GetNumberOfPoints() << " somata!" << std::endl;
// 	#endif
// 	for(unsigned int ii = 0; ii < cellSomata->GetNumberOfPoints(); ++ii)
// 	{
// 		double somaLocation[3];
// 		cellSomata->GetPoint(ii, somaLocation);
// 		
// 		// check distance to all columns
// 		double depth = 1E09;
// 		double minDistance = 1E09;
// 		int closestColumn = 0;
// 		std::map< int, Column * >::const_iterator columnIt;
// 		for(columnIt = barrelColumns.begin(); columnIt != barrelColumns.end(); ++columnIt)
// 		{
// 			int ID = columnIt->first;
// 			double t, projectedPt[3];
// 			double dist = vtkLine::DistanceToLine(somaLocation, barrelColumns[ID]->top, barrelColumns[ID]->bottom, t, projectedPt);
// 			dist = sqrt(dist);
// 			if(dist > maxDistances[ID] || t < 0 || t > 1)
// 				continue;
// 			if(dist < minDistances[ID])
// 			{
// 				if(dist < minDistance)
// 				{
// 					minDistance = dist;
// 					closestColumn = ID;
// 					depth = t;
// 				}
// 				continue;
// 			}
// 			
// 			PolyDataPointerType polyData = PolyDataPointerType::New();
// 			polyData->Allocate(1);
// 			PointsPointerType points = PointsPointerType::New();
// 			points->SetDataTypeToFloat();
// 			PolygonPointerType poly = PolygonPointerType::New();
// 			poly->GetPointIds()->SetNumberOfIds(radialContours[ID].size());
// 			for(int jj = 0; jj < radialContours[ID].size(); ++jj)
// 			{
// 				double tmp[3];
// 				for(int kk = 0; kk < 3; ++kk)
// 					tmp[kk] = projectedPt[kk] + radialContours[ID][jj][kk];
// 				points->InsertNextPoint(tmp);
// 				poly->GetPointIds()->SetId(jj, jj);
// 			}
// 			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
// 			polyData->SetPoints(points);
// 			polyData->Update();
// 			double closestPoint[3], pCoords[3], dist2;
// 			int subId;
// 			double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
// 			int insidePolygon = polyData->GetCell(0)->EvaluatePosition(somaLocation, closestPoint, subId, pCoords, dist2, weights);
// 			if(insidePolygon == 1)
// 			{
// 				if(dist < minDistance)
// 				{
// 					minDistance = dist;
// 					closestColumn = ID;
// 					depth = t;
// 				}
// 			}
// 			delete [] weights;
// 		}
// 		// assign landmark to closest soma column
// 		if(closestColumn)
// 		{
// #ifdef DEBUG
// 			std::cout << "found closest column: " << closestColumn << std::endl;
// 			std::cout << "depth: " << depth << std::endl;
// #endif
// 			std::map< double, int >::const_iterator depthIt;
// 			for(depthIt = layerDepthsMO2011.begin(); depthIt != layerDepthsMO2011.end(); ++depthIt)
// 			{
// 				if(depth <= depthIt->first)
// 				{
// #ifdef DEBUG
// 					std::cout << "found closest layer MO: " << depthIt->second << std::endl;
// #endif
// 					int layer = depthIt->second;
// 					int oldNrOfPoints = layerNumbersMO2011[layer][closestColumn]->GetNumberOfPoints();
// 					layerNumbersMO2011[layer][closestColumn]->InsertPoint(oldNrOfPoints, somaLocation);
// 					break;
// 				}
// 			}
// 			for(depthIt = layerDepthsHSM2010.begin(); depthIt != layerDepthsHSM2010.end(); ++depthIt)
// 			{
// 				if(depth <= depthIt->first)
// 				{
// #ifdef DEBUG
// 					std::cout << "found closest layer HSM: " << depthIt->second << std::endl;
// #endif
// 					int layer = depthIt->second;
// 					int oldNrOfPoints = layerNumbersHSM2010[layer][closestColumn]->GetNumberOfPoints();
// 					layerNumbersHSM2010[layer][closestColumn]->InsertPoint(oldNrOfPoints, somaLocation);
// 					break;
// 				}
// 			}
// 		}
// 	}
// 	
// 	// correction if interneuron counts
// 	if(IN)
// 	{
// 		std::map< int, PointsPointerType >::const_iterator somaColIt;
// 		for(layerNameIt = layerNamesMO2011.begin(); layerNameIt != layerNamesMO2011.end(); ++layerNameIt)
// 		{
// 			int layer = layerNameIt->first;
// 			for(somaColIt = layerNumbersMO2011[layer].begin(); somaColIt != layerNumbersMO2011[layer].end(); ++somaColIt)
// 			{
// 				int ID = somaColIt->first;
// 				if(barrelField.find(ID) != barrelField.end() && SBF->avgBarrels.find(ID) != SBF->avgBarrels.end())
// 				{
// #ifdef DEBUG
// std::cout << "computing soma z profile for column " << ID << std::endl;
// #endif
// 					Profile * zProfile = computeColumnSomaProfile(somaColIt->second, barrelColumns[ID]->top, barrelColumns[ID]->bottom);
// 					layerProfilesCorrectedMO2011[layer].insert(std::pair< int, Profile * >(ID, zProfile));
// 				}
// 			}
// 			correctSomaProfiles(layerProfilesCorrectedMO2011[layer]);
// 		}
// 		for(layerNameIt = layerNamesHSM2010.begin(); layerNameIt != layerNamesHSM2010.end(); ++layerNameIt)
// 		{
// 			int layer = layerNameIt->first;
// 			for(somaColIt = layerNumbersHSM2010[layer].begin(); somaColIt != layerNumbersHSM2010[layer].end(); ++somaColIt)
// 			{
// 				int ID = somaColIt->first;
// 				if(barrelField.find(ID) != barrelField.end() && SBF->avgBarrels.find(ID) != SBF->avgBarrels.end())
// 				{
// #ifdef DEBUG
// std::cout << "computing soma z profile for column " << ID << std::endl;
// #endif
// 					Profile * zProfile = computeColumnSomaProfile(somaColIt->second, barrelColumns[ID]->top, barrelColumns[ID]->bottom);
// 					layerProfilesCorrectedHSM2010[layer].insert(std::pair< int, Profile * >(ID, zProfile));
// 				}
// 			}
// 			correctSomaProfiles(layerProfilesCorrectedHSM2010[layer]);
// 		}
// 	}
// 	
// 	
// 	// write separate landmark files for each column
// 	#ifdef DEBUG
// 	std::cout << "Writing landmark files for all soma columns..." << std::endl;
// 	#endif
// 	std::map< int, PointsPointerType >::const_iterator somaColIt;
// 	for(layerNameIt = layerNamesMO2011.begin(); layerNameIt != layerNamesMO2011.end(); ++layerNameIt)
// 	{
// 		int layer = layerNameIt->first;
// 		for(somaColIt = layerNumbersMO2011[layer].begin(); somaColIt != layerNumbersMO2011[layer].end(); ++somaColIt)
// 		{
// 			int ID = somaColIt->first;
// 			std::string ofName(outputFilename);
// 			ofName += "_somata_MO_";
// 			ofName += layerNameIt->second;
// 			ofName += "_";
// 			ofName += SBF->int2Labels[ID];
// 			Reader * somaColWriter = new Reader(ofName.c_str(), ofName.c_str());
// 			somaColWriter->writeLandmarkFile(somaColIt->second);
// 			delete somaColWriter;
// 			if(layer == 1) delete barrelColumns[ID];
// 		}
// 	}
// 	for(layerNameIt = layerNamesHSM2010.begin(); layerNameIt != layerNamesHSM2010.end(); ++layerNameIt)
// 	{
// 		int layer = layerNameIt->first;
// 		for(somaColIt = layerNumbersHSM2010[layer].begin(); somaColIt != layerNumbersHSM2010[layer].end(); ++somaColIt)
// 		{
// 			int ID = somaColIt->first;
// 			std::string ofName(outputFilename);
// 			ofName += "_somata_HSM_";
// 			ofName += layerNameIt->second;
// 			ofName += "_";
// 			ofName += SBF->int2Labels[ID];
// 			Reader * somaColWriter = new Reader(ofName.c_str(), ofName.c_str());
// 			somaColWriter->writeLandmarkFile(somaColIt->second);
// 			delete somaColWriter;
// 		}
// 	}
// 	
// 	std::string ofNameAll(outputFilename);
// 	ofNameAll += "_all_somata.csv";
// 	std::ofstream SomaList;
// 	SomaList.open(ofNameAll.c_str());
// 	SomaList << "# Layer\tColumn\tNr of somata" << std::endl;
// 	SomaList << "# MO 2011 layer borders" << std::endl;
// 	for(layerNameIt = layerNamesMO2011.begin(); layerNameIt != layerNamesMO2011.end(); ++layerNameIt)
// 	{
// 		int layer = layerNameIt->first;
// 		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
// 		{
// 			int ID = *labelIt;
// 			SomaList << layerNameIt->second.c_str() << "\t" << SBF->int2Labels[ID];
// 			if(layerNumbersMO2011[layer].find(ID) != layerNumbersMO2011[layer].end())
// 			{
// 				SomaList << "\t" << layerNumbersMO2011[layer][ID]->GetNumberOfPoints();
// 			}
// 			SomaList << std::endl;
// 		}
// 	}
// 	SomaList << std::endl;
// 	SomaList << "# HSM 2010 layer borders" << std::endl;
// 	for(layerNameIt = layerNamesHSM2010.begin(); layerNameIt != layerNamesHSM2010.end(); ++layerNameIt)
// 	{
// 		int layer = layerNameIt->first;
// 		for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
// 		{
// 			int ID = *labelIt;
// 			SomaList << layerNameIt->second.c_str() << "\t" << SBF->int2Labels[ID];
// 			if(layerNumbersHSM2010[layer].find(ID) != layerNumbersHSM2010[layer].end())
// 			{
// 				SomaList << "\t" << layerNumbersHSM2010[layer][ID]->GetNumberOfPoints();
// 			}
// 			SomaList << std::endl;
// 		}
// 	}
// 	SomaList.close();
// 	
// 	if(IN)
// 	{
// 		std::string ofNameAllIN(outputFilename);
// 		ofNameAllIN += "_all_somata_corrected.csv";
// 		std::ofstream SomaList2;
// 		SomaList2.open(ofNameAllIN.c_str());
// 		SomaList2 << "# Layer\tColumn\tNr of somata (corrected)" << std::endl;
// 		SomaList2 << "# MO 2011 layer borders" << std::endl;
// 		for(layerNameIt = layerNamesMO2011.begin(); layerNameIt != layerNamesMO2011.end(); ++layerNameIt)
// 		{
// 			int layer = layerNameIt->first;
// 			for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
// 			{
// 				int ID = *labelIt;
// 				SomaList2 << layerNameIt->second.c_str() << "\t" << SBF->int2Labels[ID];
// 				if(layerProfilesCorrectedMO2011[layer].find(ID) != layerProfilesCorrectedMO2011[layer].end())
// 				{
// 					SomaList2 << "\t" << layerProfilesCorrectedMO2011[layer][ID]->getIntegral();
// 				}
// 				SomaList2 << std::endl;
// 			}
// 		}
// 		SomaList2 << std::endl;
// 		SomaList2 << "# HSM 2010 layer borders" << std::endl;
// 		for(layerNameIt = layerNamesHSM2010.begin(); layerNameIt != layerNamesHSM2010.end(); ++layerNameIt)
// 		{
// 			int layer = layerNameIt->first;
// 			for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
// 			{
// 				int ID = *labelIt;
// 				SomaList2 << layerNameIt->second.c_str() << "\t" << SBF->int2Labels[ID];
// 				if(layerProfilesCorrectedHSM2010[layer].find(ID) != layerProfilesCorrectedHSM2010[layer].end())
// 				{
// 					SomaList2 << "\t" << layerProfilesCorrectedHSM2010[layer][ID]->getIntegral();
// 				}
// 				SomaList2 << std::endl;
// 			}
// 		}
// 		SomaList2.close();
// 	}
// }

//compute profiles separately for each column; i.e. do not take overlap into account
void LandmarkAnalyzer::computeSeparateColumnProfiles ( const char* outputFilename )
{
	if(!piaFlag || !wmFlag || !barrelFlag || !cellFlag)
	{
		std::cout << "Error! Cannot compute column profiles. Incomplete input data!" << std::endl;
		return;
	}
	outFilenameGlobal = outputFilename;
	
	std::cout << "Computing separate column profiles..." << std::endl;
	if(columnSomaProfiles.size()) columnSomaProfiles.clear();
	
	// create columns from barrel contours
	// and intersection points with Pia/WM
	// allocate data structures for storing
	// somata inside columns
	std::map< int, PointsPointerType > somaColumns;
	std::map< int, Column * > barrelColumns;
	std::map< int, std::vector< std::vector< double > > > radialContours;
	std::map< int, double > minDistances;
	std::map< int, double > maxDistances;
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrelField.find(ID) != barrelField.end())
		{
			PointsPointerType somaColumn = PointsPointerType::New();
			somaColumn->Allocate(1);
			somaColumn->SetDataTypeToFloat();
			somaColumns.insert(std::pair< int, PointsPointerType >(ID, somaColumn));
			
			double columnAxis[3], columnCenter[3], columnTop[3], columnBottom[3];
			for(int ii = 0; ii < 3; ++ii)
			{
				columnAxis[ii] = barrelField[ID]->top[ii] - barrelField[ID]->bottom[ii];
				columnCenter[ii] = 0.5*(barrelField[ID]->top[ii] + barrelField[ID]->bottom[ii]);
			}
			vtkMath::Normalize(columnAxis);
			
			pia->intersectLine(columnAxis, columnCenter);
			if(!pia->isValid())
			{
				std::cout << "Error! Could not compute intersection of column axis with Pia!" << std::endl;
				continue;
			}
			WM->intersectLine(columnAxis, columnCenter);
			if(!WM->isValid())
			{
				std::cout << "Error! Could not compute intersection of column axis with WM!" << std::endl;
				continue;
			}
			pia->getLastIntersectPoint(columnTop);
			WM->getLastIntersectPoint(columnBottom);
			Column * thisCol = new Column(barrelField[ID]->contours, columnTop, columnBottom);
// 			Column * thisCol = new Column;
// 			thisCol->top = columnTop;
// 			thisCol->bottom = columnBottom;
// 			thisCol->contours = PolyDataPointerType::New();
// 			thisCol->contours->DeepCopy(barrelField[ID]->contours);
			barrelColumns.insert(std::pair< int, Column * >(ID, thisCol));
			std::vector< std::vector< double > > contourPts;
			double minDist = 1E09, maxDist = 0;
			for(int ii = 0; ii < thisCol->contours->GetCell(0)->GetNumberOfPoints(); ++ii)
			{
				double pt[3], dist = 0, t, closestPt[3];
				std::vector< double > radialPt;
				thisCol->contours->GetCell(0)->GetPoints()->GetPoint(ii, pt);
				dist = vtkLine::DistanceToLine(pt, thisCol->top, thisCol->bottom, t, closestPt);
				dist = sqrt(dist);
				if(dist < minDist)
					minDist = dist;
				if(dist > maxDist)
					maxDist = dist;
				for(int jj = 0; jj < 3; ++jj)
					radialPt.push_back(pt[jj] - closestPt[jj]);
				contourPts.push_back(radialPt);
			}
			radialContours.insert(std::pair< int, std::vector< std::vector< double > > >(ID, contourPts));
			minDistances.insert(std::pair< int, double >(ID, minDist));
			maxDistances.insert(std::pair< int, double >(ID, maxDist));
// 			#ifdef DEBUG
// 			double boundingBox[6];
// 			barrelField[ID]->contours->GetBounds(boundingBox);
// 			std::cout << "*************************************" << std::endl;
// 			std::cout << "Barrel " << SBF->int2Labels[ID] << ":" << std::endl;
// 			std::cout << "Top @ [" << barrelField[ID]->top[0] << "," << barrelField[ID]->top[1] << "," << barrelField[ID]->top[2] << "]" << std::endl;
// 			std::cout << "Bottom @ [" << barrelField[ID]->bottom[0] << "," << barrelField[ID]->bottom[1] << "," << barrelField[ID]->bottom[2] << "]" << std::endl;
// 			std::cout << "Center @ [" << columnCenter[0] << "," << columnCenter[1] << "," << columnCenter[2] << "]" << std::endl;
// 			std::cout << "Axis = [" << columnAxis[0] << "," << columnAxis[1] << "," << columnAxis[2] << "]" << std::endl;
// 			std::cout << "Bounding box:" << std::endl;
// 			std::cout << "[" << boundingBox[0] << "," << boundingBox[1] << "]" << std::endl;
// 			std::cout << "[" << boundingBox[2] << "," << boundingBox[3] << "]" << std::endl;
// 			std::cout << "[" << boundingBox[4] << "," << boundingBox[5] << "]" << std::endl;
// 			std::cout << "Column " << SBF->int2Labels[ID] << ":" << std::endl;
// 			std::cout << "Top @ [" << barrelColumns[ID]->top[0] << "," << barrelColumns[ID]->top[1] << "," << barrelColumns[ID]->top[2] << "]" << std::endl;
// 			std::cout << "Pia pt @ [" << columnTop[0] << "," << columnTop[1] << "," << columnTop[2] << "]" << std::endl;
// 			std::cout << "Bottom @ [" << barrelColumns[ID]->bottom[0] << "," << barrelColumns[ID]->bottom[1] << "," << barrelColumns[ID]->bottom[2] << "]" << std::endl;
// 			std::cout << "WM pt @ [" << columnBottom[0] << "," << columnBottom[1] << "," << columnBottom[2] << "]" << std::endl;
// 			std::cout << "Min radius = " << minDist << "um" << std::endl;
// 			std::cout << "Max radius = " << maxDist << "um" << std::endl;
// 			#endif
		}
	}
	
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrelField.find(ID) != barrelField.end())
		{
			std::flush(std::cout << "Checking for somata in column " << SBF->int2Labels[ID] << std::endl);
			unsigned int nrOfSomata = cellSomata->GetNumberOfPoints();
			for(unsigned int ii = 0; ii < nrOfSomata; ++ii)
			{
				double somaLocation[3];
				cellSomata->GetPoint(ii, somaLocation);
				double t, projectedPt[3];
				double dist = vtkLine::DistanceToLine(somaLocation, barrelColumns[ID]->top, barrelColumns[ID]->bottom, t, projectedPt);
				dist = sqrt(dist);
				if(dist > maxDistances[ID] || t < 0 || t > 1)
					continue;
				if(dist < minDistances[ID])
				{
					int oldNrOfPoints = somaColumns[ID]->GetNumberOfPoints();
					somaColumns[ID]->InsertPoint(oldNrOfPoints, somaLocation);
					continue;
				}
				
				PolyDataPointerType polyData = PolyDataPointerType::New();
				polyData->Allocate(1);
				PointsPointerType points = PointsPointerType::New();
				points->SetDataTypeToFloat();
				PolygonPointerType poly = PolygonPointerType::New();
				poly->GetPointIds()->SetNumberOfIds(radialContours[ID].size());
				for(int jj = 0; jj < radialContours[ID].size(); ++jj)
				{
					double tmp[3];
					for(int kk = 0; kk < 3; ++kk)
						tmp[kk] = projectedPt[kk] + radialContours[ID][jj][kk];
					points->InsertNextPoint(tmp);
					poly->GetPointIds()->SetId(jj, jj);
				}
				polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
				polyData->SetPoints(points);
				polyData->Update();
				double closestPoint[3], pCoords[3], dist2;
				int subId;
				double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
				int insidePolygon = polyData->GetCell(0)->EvaluatePosition(somaLocation, closestPoint, subId, pCoords, dist2, weights);
				if(insidePolygon == 1)
				{
					int oldNrOfPoints = somaColumns[ID]->GetNumberOfPoints();
					somaColumns[ID]->InsertPoint(oldNrOfPoints, somaLocation);
				}
				delete [] weights;
			}
			
			Profile * zProfile = computeColumnSomaProfile(somaColumns[ID], barrelColumns[ID]->top, barrelColumns[ID]->bottom);
			columnSomaProfiles.insert(std::pair< int, Profile * >(ID, zProfile));
		}
	}
	
	// write separate landmark files for each column
	#ifdef DEBUG
	std::cout << "Writing landmark files for all soma columns..." << std::endl;
	#endif
	std::map< int, PointsPointerType >::const_iterator somaColIt;
	for(somaColIt = somaColumns.begin(); somaColIt != somaColumns.end(); ++somaColIt)
	{
		int ID = somaColIt->first;
		std::string ofName(outputFilename);
		ofName += "_separate_somata_";
		ofName += SBF->int2Labels[ID];
		Reader * somaColWriter = new Reader(ofName.c_str(), ofName.c_str());
		somaColWriter->writeLandmarkFile(somaColIt->second);
		delete somaColWriter;
		delete barrelColumns[ID];
	}
	
	std::string profileName(outputFilename);
	profileName += "_separate";
	writeColumnSomaProfiles(profileName.c_str());
	
	std::string ofName(outputFilename);
	ofName += "_all_somata_separate.csv";
	std::ofstream SomaList;
	SomaList.open(ofName.c_str());
	SomaList << "# Column\tNr of somata" << std::endl;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		SomaList << SBF->int2Labels[ID];
		if(somaColumns.find(ID) != somaColumns.end())
		{
			SomaList << "\t" << somaColumns[ID]->GetNumberOfPoints();
		}
		SomaList << std::endl;
	}
	SomaList.close();
}

// count all somata inside volume defined by barrel columns,
// Pia and WM and assign them to column/septum
void LandmarkAnalyzer::countCellsInSeptum ( Surface* S1Surface, const char* outputFilename )
{
	if(!piaFlag || !wmFlag || !barrelFlag || !cellFlag)
	{
		std::cout << "Error! Cannot compute column profiles. Incomplete input data!" << std::endl;
		std::cout << "Pia:\t" << piaFlag << std::endl;
		std::cout << "WM:\t" << wmFlag << std::endl;
		std::cout << "Barrels:\t" << barrelFlag << std::endl;
		std::cout << "Somata:\t" << cellFlag << std::endl;
		return;
	}
	outFilenameGlobal = outputFilename;
	
	std::cout << "Counting somata in septum..." << std::endl;
	
	// create columns from barrel contours
	// and intersection points with Pia/WM
	// merge all column contours for computation of convex hull
	std::map< int, Column * > barrelColumns;
	std::map< int, std::vector< std::vector< double > > > radialContours;
	std::map< int, double > minDistances;
	std::map< int, double > maxDistances;
	std::list< int >::const_iterator labelIt;
	unsigned int totalPtCnt = 0;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrelField.find(ID) != barrelField.end())
		{
			double columnAxis[3], columnCenter[3], columnTop[3], columnBottom[3];
			for(int ii = 0; ii < 3; ++ii)
			{
				columnAxis[ii] = barrelField[ID]->top[ii] - barrelField[ID]->bottom[ii];
				columnCenter[ii] = 0.5*(barrelField[ID]->top[ii] + barrelField[ID]->bottom[ii]);
			}
			vtkMath::Normalize(columnAxis);
			
			pia->intersectLine(columnAxis, columnCenter);
			if(!pia->isValid())
			{
				std::cout << "Error! Could not compute intersection of column axis with Pia!" << std::endl;
				continue;
			}
			WM->intersectLine(columnAxis, columnCenter);
			if(!WM->isValid())
			{
				std::cout << "Error! Could not compute intersection of column axis with WM!" << std::endl;
				continue;
			}
			pia->getLastIntersectPoint(columnTop);
			WM->getLastIntersectPoint(columnBottom);
			Column * thisCol = new Column(barrelField[ID]->contours, columnTop, columnBottom);
			barrelColumns.insert(std::pair< int, Column * >(ID, thisCol));
			
			PointsPointerType contourPts = barrelField[ID]->contours->GetCell(0)->GetPoints();
			std::vector< std::vector< double > > radialContourPts;
			double minDist = 1E09, maxDist = 0;
			for(int ii = 0; ii < contourPts->GetNumberOfPoints(); ++ii)
			{
				double tmp[3], projectedPt[3], dist = 0, t;
				double * contourTopPt = new double[3], * contourBottomPt = new double[3];
				
				contourPts->GetPoint(ii, tmp);
				dist = vtkLine::DistanceToLine(tmp, columnTop, columnBottom, t, projectedPt);
				dist = sqrt(dist);
				if(dist < minDist)
					minDist = dist;
				if(dist > maxDist)
					maxDist = dist;
				
				vtkMath::Subtract(tmp, projectedPt, tmp);
				vtkMath::Add(columnTop, tmp, contourTopPt);
				vtkMath::Add(columnBottom, tmp, contourBottomPt);
				totalPtCnt += 2;
				
				std::vector< double > radialPt;
				for(int jj = 0; jj < 3; ++jj)
					radialPt.push_back(tmp[jj]);
				radialContourPts.push_back(radialPt);
			}
			radialContours.insert(std::pair< int, std::vector< std::vector< double > > >(ID, radialContourPts));
			minDistances.insert(std::pair< int, double >(ID, minDist));
			maxDistances.insert(std::pair< int, double >(ID, maxDist));
		}
	}
	
	PolyDataPointerType barrelFieldHull = S1Surface->ptr();
	double S1bounds[6], center[3];
	barrelFieldHull->GetBounds(S1bounds);
	for(int ii = 0; ii < 3; ++ii)
		center[ii] = 0.5*(S1bounds[2*ii] + S1bounds[2*ii+1]);
	
	#ifdef DEBUG
	AmiraSpatialGraph * tmpSG = new AmiraSpatialGraph;
	tmpSG->addPolyDataObject(barrelFieldHull, Barrel);
	#endif
	
	std::vector< PlanePointerType > hullPlanes;
	for(int ii = 0; ii < barrelFieldHull->GetNumberOfCells(); ++ii)
	{
		double normal[3], origin[3], pCoords[3], * weights = new double[barrelFieldHull->GetCell(ii)->GetNumberOfPoints()];
		int subId;
		barrelFieldHull->GetCell(ii)->GetParametricCenter(pCoords);
		barrelFieldHull->GetCell(ii)->EvaluateLocation(subId, pCoords, origin, weights);
		delete [] weights;
		vtkPolygon::ComputeNormal(barrelFieldHull->GetCell(ii)->GetPoints(), normal);
		if(vtkMath::Norm(normal))
		{
			PlanePointerType newPlane = PlanePointerType::New();
			newPlane->SetNormal(normal);
			newPlane->SetOrigin(origin);
			double f = newPlane->EvaluateFunction(center);
			// dirty trick to circumvent numerical problems
			// in cases where points nearly lie on a line
			// If center is not on "inside", do not accept
			// this plane
			if(f < 0)
				hullPlanes.push_back(newPlane);
		}
	}
	
	// for performance reasons, implement everything
	// explicitly within the loops here
	// otherwise prohibitively many function calls...
	PointsPointerType septumSomata = PointsPointerType::New();
	PointsPointerType columnSomata = PointsPointerType::New();
	septumSomata->SetDataTypeToFloat();
	columnSomata->SetDataTypeToFloat();
	unsigned int nrOfSomata = cellSomata->GetNumberOfPoints();
	for(int ii = 0; ii < nrOfSomata; ++ii)
	{
		double somaLocation[3];
		cellSomata->GetPoint(ii, somaLocation);
		// no need to waste computing time...
		double tol[3] = {1E-6,1E-6,1E-6};
		if(!vtkMath::PointIsWithinBounds(somaLocation, S1bounds, tol))
			continue;
		
		// candidate points:
		// first check whether point is
		// really inside surface defined
		// by convex hull
		bool inside = 1;
		for(int jj = 0; jj < hullPlanes.size(); ++jj)
		{
			if(hullPlanes[jj]->EvaluateFunction(somaLocation) > 0)
			{
				inside = 0;
				break;
			}
		}
		if(!inside)
			continue;
		
		// Pia-WM: checking along global
		// z axis sufficient
		double zAxis[] = {0,0,1};
		pia->intersectLine(zAxis, somaLocation);
		WM->intersectLine(zAxis, somaLocation);
		// careful: if no intersection with is found
		// with WM, this may well be because the WM
		// is not present below all of S1; thus accept
		// point in these cases
		if(pia->isIntersectionFound() && WM->isIntersectionFound())
		{
			double piaIntersectPt[3], WMIntersectPt[3], projectedPt[3], t;
			pia->getLastIntersectPoint(piaIntersectPt);
			WM->getLastIntersectPoint(WMIntersectPt);
			vtkLine::DistanceToLine(somaLocation, piaIntersectPt, WMIntersectPt, t, projectedPt);
			if(t < 0 || t > 1)
				continue;
		}
		
		// check distance to all columns
		bool isInSeptum = 1;
		std::map< int, Column * >::const_iterator columnIt;
		for(columnIt = barrelColumns.begin(); columnIt != barrelColumns.end(); ++columnIt)
		{
			int ID = columnIt->first;
			double t, projectedPt[3];
			double dist = vtkLine::DistanceToLine(somaLocation, barrelColumns[ID]->top, barrelColumns[ID]->bottom, t, projectedPt);
			dist = sqrt(dist);
			if(dist > maxDistances[ID] || t < 0 || t > 1)
				continue;
			if(dist < minDistances[ID])
			{
				isInSeptum = 0;
				break;
			}
			
			PolyDataPointerType polyData = PolyDataPointerType::New();
			polyData->Allocate(1);
			PointsPointerType points = PointsPointerType::New();
			points->SetDataTypeToFloat();
			PolygonPointerType poly = PolygonPointerType::New();
			poly->GetPointIds()->SetNumberOfIds(radialContours[ID].size());
			for(int jj = 0; jj < radialContours[ID].size(); ++jj)
			{
				double tmp[3];
				for(int kk = 0; kk < 3; ++kk)
					tmp[kk] = projectedPt[kk] + radialContours[ID][jj][kk];
				points->InsertNextPoint(tmp);
				poly->GetPointIds()->SetId(jj, jj);
			}
			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
			polyData->SetPoints(points);
			polyData->Update();
			double closestPoint[3], pCoords[3], dist2;
			int subId;
			double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
			int insidePolygon = polyData->GetCell(0)->EvaluatePosition(somaLocation, closestPoint, subId, pCoords, dist2, weights);
			if(insidePolygon == 1)
			{
				isInSeptum = 0;
				delete [] weights;
				break;
			}
			delete [] weights;
		}
		// assign landmark to septum
		if(isInSeptum)
		{
			unsigned int oldNrOfPoints = septumSomata->GetNumberOfPoints();
			septumSomata->InsertPoint(oldNrOfPoints, somaLocation);
		}
		else
		{
			unsigned int oldNrOfPoints = columnSomata->GetNumberOfPoints();
			columnSomata->InsertPoint(oldNrOfPoints, somaLocation);
		}
	}
	
	std::cout << "Computing total septum volume..." << std::endl;
	double spacing = SPACING;
	PointsPointerType septumVoxels = PointsPointerType::New();
	PointsPointerType columnVoxels = PointsPointerType::New();
	septumVoxels->SetDataTypeToFloat();
	columnVoxels->SetDataTypeToFloat();
	ImageDataPointerType S1Volume = createImageVolume(S1bounds);
	for(int x = S1Volume->GetExtent()[0]; x <= S1Volume->GetExtent()[1]; ++x)
		for(int y = S1Volume->GetExtent()[2]; y <= S1Volume->GetExtent()[3]; ++y)
			for(int z = S1Volume->GetExtent()[4]; z <= S1Volume->GetExtent()[5]; ++z)
			{
				double pt[3];
				pt[0] = x*spacing, pt[1] = y*spacing, pt[2] = z*spacing;
				// no need to waste computing time...
				double tol[3] = {1E-6,1E-6,1E-6};
				if(!vtkMath::PointIsWithinBounds(pt, S1bounds, tol))
					continue;
				
				// candidate points:
				// first check whether point is
				// really inside surface defined
				// by convex hull
				bool inside = 1;
				for(int jj = 0; jj < hullPlanes.size(); ++jj)
				{
					if(hullPlanes[jj]->EvaluateFunction(pt) > 0)
					{
						inside = 0;
						break;
					}
				}
				if(!inside)
					continue;
				
				// Pia-WM: checking along global
				// z axis sufficient
				double zAxis[] = {0,0,1};
				pia->intersectLine(zAxis, pt);
				WM->intersectLine(zAxis, pt);
				// careful: if no intersection with is found
				// with WM, this may well be because the WM
				// is not present below all of S1; thus accept
				// point in these cases
				if(pia->isIntersectionFound() && WM->isIntersectionFound())
				{
					double piaIntersectPt[3], WMIntersectPt[3], projectedPt[3], t;
					pia->getLastIntersectPoint(piaIntersectPt);
					WM->getLastIntersectPoint(WMIntersectPt);
					vtkLine::DistanceToLine(pt, piaIntersectPt, WMIntersectPt, t, projectedPt);
					if(t < 0 || t > 1)
						continue;
				}
				
				// check distance to all columns
				bool isInSeptum = 1;
				std::map< int, Column * >::const_iterator columnIt;
				for(columnIt = barrelColumns.begin(); columnIt != barrelColumns.end(); ++columnIt)
				{
					int ID = columnIt->first;
					double t, projectedPt[3];
					double dist = vtkLine::DistanceToLine(pt, barrelColumns[ID]->top, barrelColumns[ID]->bottom, t, projectedPt);
					dist = sqrt(dist);
					if(dist > maxDistances[ID] || t < 0 || t > 1)
						continue;
					if(dist < minDistances[ID])
					{
						isInSeptum = 0;
						break;
					}
					
					PolyDataPointerType polyData = PolyDataPointerType::New();
					polyData->Allocate(1);
					PointsPointerType points = PointsPointerType::New();
					points->SetDataTypeToFloat();
					PolygonPointerType poly = PolygonPointerType::New();
					poly->GetPointIds()->SetNumberOfIds(radialContours[ID].size());
					for(int jj = 0; jj < radialContours[ID].size(); ++jj)
					{
						double tmp[3];
						for(int kk = 0; kk < 3; ++kk)
							tmp[kk] = projectedPt[kk] + radialContours[ID][jj][kk];
						points->InsertNextPoint(tmp);
						poly->GetPointIds()->SetId(jj, jj);
					}
					polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
					polyData->SetPoints(points);
					polyData->Update();
					double closestPoint[3], pCoords[3], dist2;
					int subId;
					double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
					int insidePolygon = polyData->GetCell(0)->EvaluatePosition(pt, closestPoint, subId, pCoords, dist2, weights);
					if(insidePolygon == 1)
					{
						isInSeptum = 0;
						delete [] weights;
						break;
					}
					delete [] weights;
				}
				// assign landmark to septum
				if(isInSeptum)
				{
					unsigned int oldNrOfPoints = septumVoxels->GetNumberOfPoints();
					septumVoxels->InsertPoint(oldNrOfPoints, pt);
				}
				else
				{
					unsigned int oldNrOfPoints = columnVoxels->GetNumberOfPoints();
					columnVoxels->InsertPoint(oldNrOfPoints, pt);
				}
			}
	
	Profile * septumSomaProfile = computeSeptumProfile(septumSomata, barrelColumns);
	Profile * columnSomaProfile = computeSeptumProfile(columnSomata, barrelColumns);
	Profile * septumVolumeProfile = computeSeptumProfile(septumVoxels, barrelColumns);
	Profile * columnVolumeProfile = computeSeptumProfile(columnVoxels, barrelColumns);
	
	std::string septumSomaName(outputFilename);
	septumSomaName += "_septum_somata";
	writeZProfile(septumSomaProfile, septumSomaName.c_str());
	
	std::string ofName(outputFilename);
	ofName += "_S1_septum_somata";
	Reader * septumSomaWriter = new Reader(ofName.c_str(), ofName.c_str());
	septumSomaWriter->writeLandmarkFile(septumSomata);
	delete septumSomaWriter;
	
	std::string columnSomaName(outputFilename);
	columnSomaName += "_column_somata";
	writeZProfile(columnSomaProfile, columnSomaName.c_str());
	
	std::string ofName2(outputFilename);
	ofName2 += "_S1_column_somata";
	Reader * colSomaWriter = new Reader(ofName2.c_str(), ofName2.c_str());
	colSomaWriter->writeLandmarkFile(columnSomata);
	delete colSomaWriter;
	
	std::string septumVoxelName(outputFilename);
	septumVoxelName += "_septum_voxels_";
	septumVoxelName += "%03d";
	septumVoxelName += "micron";
	char * sepVoxelStr = new char[256];
	std::sprintf(sepVoxelStr, septumVoxelName.c_str(), (int)SPACING);
	writeZProfile(septumVolumeProfile, sepVoxelStr);
	
	std::string columnVoxelName(outputFilename);
	columnVoxelName += "_column_voxels_";
	columnVoxelName += "%03d";
	columnVoxelName += "micron";
	char * colVoxelStr = new char[256];
	std::sprintf(colVoxelStr, columnVoxelName.c_str(), (int)SPACING);
	writeZProfile(columnVolumeProfile, colVoxelStr);
	
	#ifdef DEBUG
	std::string tmpSGName(outputFilename);
	tmpSGName += "_convex_hull";
	Reader * tmpSGWriter = new Reader(tmpSGName.c_str(), tmpSGName.c_str());
	tmpSGWriter->setSpatialGraph(tmpSG);
	tmpSGWriter->writeSpatialGraphFile();
	delete tmpSGWriter, delete tmpSG;
	#endif
	
	delete septumSomaProfile, delete columnSomaProfile;
	delete [] sepVoxelStr, delete [] colVoxelStr;
	delete septumVolumeProfile, delete columnVolumeProfile;
}

void LandmarkAnalyzer::countInhCellsInSeptumInLayer(Surface* S1Surface, const char* outputFilename)
{
	if(!piaFlag || !wmFlag || !barrelFlag || !barrelLayerFlag || !cellFlag)
	{
		std::cout << "Error! Cannot compute column profiles. Incomplete input data!" << std::endl;
		return;
	}
	outFilenameGlobal = outputFilename;
	
	std::cout << "Counting somata in septum in layer, correcting INs..." << std::endl;
	
	// create columns from barrel contours
	// and intersection points with Pia/WM
	// merge all column contours for computation of convex hull
	std::map< int, Column * > barrelColumns;
	std::map< int, std::vector< std::vector< double > > > radialContours;
	std::map< int, double > minDistances;
	std::map< int, double > maxDistances;
	std::list< int >::const_iterator labelIt;
	unsigned int totalPtCnt = 0;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrelFieldLayer.find(ID) != barrelFieldLayer.end())
		{
			double columnAxis[3], columnCenter[3], columnTop[3], columnBottom[3];
			for(int ii = 0; ii < 3; ++ii)
			{
				columnAxis[ii] = barrelFieldLayer[ID]->top[ii] - barrelFieldLayer[ID]->bottom[ii];
				columnCenter[ii] = 0.5*(barrelFieldLayer[ID]->top[ii] + barrelFieldLayer[ID]->bottom[ii]);
			}
			vtkMath::Normalize(columnAxis);
			
			pia->intersectLine(columnAxis, columnCenter);
			if(!pia->isValid())
			{
				std::cout << "Error! Could not compute intersection of column axis with Pia!" << std::endl;
				continue;
			}
			WM->intersectLine(columnAxis, columnCenter);
			if(!WM->isValid())
			{
				std::cout << "Error! Could not compute intersection of column axis with WM!" << std::endl;
				continue;
			}
			pia->getLastIntersectPoint(columnTop);
			WM->getLastIntersectPoint(columnBottom);
			Column * thisCol = new Column(barrelFieldLayer[ID]->contours, columnTop, columnBottom);
			barrelColumns.insert(std::pair< int, Column * >(ID, thisCol));
			
			PointsPointerType contourPts = barrelFieldLayer[ID]->contours->GetCell(0)->GetPoints();
			std::vector< std::vector< double > > radialContourPts;
			double minDist = 1E09, maxDist = 0;
			for(int ii = 0; ii < contourPts->GetNumberOfPoints(); ++ii)
			{
				double tmp[3], projectedPt[3], dist = 0, t;
				double * contourTopPt = new double[3], * contourBottomPt = new double[3];
				
				contourPts->GetPoint(ii, tmp);
				dist = vtkLine::DistanceToLine(tmp, columnTop, columnBottom, t, projectedPt);
				dist = sqrt(dist);
				if(dist < minDist)
					minDist = dist;
				if(dist > maxDist)
					maxDist = dist;
				
				vtkMath::Subtract(tmp, projectedPt, tmp);
				vtkMath::Add(columnTop, tmp, contourTopPt);
				vtkMath::Add(columnBottom, tmp, contourBottomPt);
				totalPtCnt += 2;
				
				std::vector< double > radialPt;
				for(int jj = 0; jj < 3; ++jj)
					radialPt.push_back(tmp[jj]);
				radialContourPts.push_back(radialPt);
			}
			radialContours.insert(std::pair< int, std::vector< std::vector< double > > >(ID, radialContourPts));
			minDistances.insert(std::pair< int, double >(ID, minDist));
			maxDistances.insert(std::pair< int, double >(ID, maxDist));
		}
	}
	
	PolyDataPointerType barrelFieldHull = S1Surface->ptr();
	double S1bounds[6], center[3];
	barrelFieldHull->GetBounds(S1bounds);
	for(int ii = 0; ii < 3; ++ii)
		center[ii] = 0.5*(S1bounds[2*ii] + S1bounds[2*ii+1]);
	
	#ifdef DEBUG
	AmiraSpatialGraph * tmpSG = new AmiraSpatialGraph;
	tmpSG->addPolyDataObject(barrelFieldHull, Barrel);
	#endif
	
	std::vector< PlanePointerType > hullPlanes;
	for(int ii = 0; ii < barrelFieldHull->GetNumberOfCells(); ++ii)
	{
		double normal[3], origin[3], pCoords[3], * weights = new double[barrelFieldHull->GetCell(ii)->GetNumberOfPoints()];
		int subId;
		barrelFieldHull->GetCell(ii)->GetParametricCenter(pCoords);
		barrelFieldHull->GetCell(ii)->EvaluateLocation(subId, pCoords, origin, weights);
		delete [] weights;
		vtkPolygon::ComputeNormal(barrelFieldHull->GetCell(ii)->GetPoints(), normal);
		if(vtkMath::Norm(normal))
		{
			PlanePointerType newPlane = PlanePointerType::New();
			newPlane->SetNormal(normal);
			newPlane->SetOrigin(origin);
			double f = newPlane->EvaluateFunction(center);
			// dirty trick to circumvent numerical problems
			// in cases where points nearly lie on a line
			// If center is not on "inside", do not accept
			// this plane
			if(f < 0)
				hullPlanes.push_back(newPlane);
		}
	}
	
	// for performance reasons, implement everything
	// explicitly within the loops here
	// otherwise prohibitively many function calls...
	PointsPointerType septumSomata = PointsPointerType::New();
	PointsPointerType columnSomata = PointsPointerType::New();
	septumSomata->SetDataTypeToFloat();
	columnSomata->SetDataTypeToFloat();
	unsigned int nrOfSomata = cellSomata->GetNumberOfPoints();
	for(int ii = 0; ii < nrOfSomata; ++ii)
	{
		double somaLocation[3];
		cellSomata->GetPoint(ii, somaLocation);
		// no need to waste computing time...
		double tol[3] = {1E-6,1E-6,1E-6};
		if(!vtkMath::PointIsWithinBounds(somaLocation, S1bounds, tol))
			continue;
		
		// candidate points:
		// first check whether point is
		// really inside surface defined
		// by convex hull
		bool inside = 1;
		for(int jj = 0; jj < hullPlanes.size(); ++jj)
		{
			if(hullPlanes[jj]->EvaluateFunction(somaLocation) > 0)
			{
				inside = 0;
				break;
			}
		}
		if(!inside)
			continue;
		
		// Pia-WM: checking along global
		// z axis sufficient
		double zAxis[] = {0,0,1};
		pia->intersectLine(zAxis, somaLocation);
		WM->intersectLine(zAxis, somaLocation);
		// careful: if no intersection with is found
		// with WM, this may well be because the WM
		// is not present below all of S1; thus accept
		// point in these cases
		if(pia->isIntersectionFound() && WM->isIntersectionFound())
		{
			double piaIntersectPt[3], WMIntersectPt[3], projectedPt[3], t;
			pia->getLastIntersectPoint(piaIntersectPt);
			WM->getLastIntersectPoint(WMIntersectPt);
			vtkLine::DistanceToLine(somaLocation, piaIntersectPt, WMIntersectPt, t, projectedPt);
			if(t < 0 || t > 1)
				continue;
		}
		
		// check distance to all columns
		bool isInSeptum = 1;
		std::map< int, Column * >::const_iterator columnIt;
		for(columnIt = barrelColumns.begin(); columnIt != barrelColumns.end(); ++columnIt)
		{
			int ID = columnIt->first;
			double t, projectedPt[3];
			double dist = vtkLine::DistanceToLine(somaLocation, barrelColumns[ID]->top, barrelColumns[ID]->bottom, t, projectedPt);
			dist = sqrt(dist);
			if(dist > maxDistances[ID] || t < 0 || t > 1)
				continue;
			if(dist < minDistances[ID])
			{
				isInSeptum = 0;
				break;
			}
			
			PolyDataPointerType polyData = PolyDataPointerType::New();
			polyData->Allocate(1);
			PointsPointerType points = PointsPointerType::New();
			points->SetDataTypeToFloat();
			PolygonPointerType poly = PolygonPointerType::New();
			poly->GetPointIds()->SetNumberOfIds(radialContours[ID].size());
			for(int jj = 0; jj < radialContours[ID].size(); ++jj)
			{
				double tmp[3];
				for(int kk = 0; kk < 3; ++kk)
					tmp[kk] = projectedPt[kk] + radialContours[ID][jj][kk];
				points->InsertNextPoint(tmp);
				poly->GetPointIds()->SetId(jj, jj);
			}
			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
			polyData->SetPoints(points);
			polyData->Update();
			double closestPoint[3], pCoords[3], dist2;
			int subId;
			double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
			int insidePolygon = polyData->GetCell(0)->EvaluatePosition(somaLocation, closestPoint, subId, pCoords, dist2, weights);
			if(insidePolygon == 1)
			{
				isInSeptum = 0;
				delete [] weights;
				break;
			}
			delete [] weights;
		}
		// assign landmark to septum
		if(isInSeptum)
		{
			unsigned int oldNrOfPoints = septumSomata->GetNumberOfPoints();
			septumSomata->InsertPoint(oldNrOfPoints, somaLocation);
		}
		else
		{
			unsigned int oldNrOfPoints = columnSomata->GetNumberOfPoints();
			columnSomata->InsertPoint(oldNrOfPoints, somaLocation);
		}
	}
	
	Profile * septumSomaProfile, * columnSomaProfile;
	septumSomaProfile = computeCorrectedSeptumProfile(septumSomata, barrelField, barrelFieldLayer);
	columnSomaProfile = computeCorrectedSeptumProfile(columnSomata, barrelField, barrelFieldLayer);
	
	std::string septumSomaName(outputFilename);
	septumSomaName += "_septum_somata";
	writeZProfile(septumSomaProfile, septumSomaName.c_str());
	
	std::string ofName(outputFilename);
	ofName += "_S1_septum_somata";
	Reader * septumSomaWriter = new Reader(ofName.c_str(), ofName.c_str());
	septumSomaWriter->writeLandmarkFile(septumSomata);
	delete septumSomaWriter;
	
	std::string columnSomaName(outputFilename);
	columnSomaName += "_column_somata";
	writeZProfile(columnSomaProfile, columnSomaName.c_str());
	
	std::string ofName2(outputFilename);
	ofName2 += "_S1_column_somata";
	Reader * colSomaWriter = new Reader(ofName2.c_str(), ofName2.c_str());
	colSomaWriter->writeLandmarkFile(columnSomata);
	delete colSomaWriter;
	
	#ifdef DEBUG
	std::string tmpSGName(outputFilename);
	tmpSGName += "_convex_hull";
	Reader * tmpSGWriter = new Reader(tmpSGName.c_str(), tmpSGName.c_str());
	tmpSGWriter->setSpatialGraph(tmpSG);
	tmpSGWriter->writeSpatialGraphFile();
	delete tmpSGWriter, delete tmpSG;
	#endif
	
	delete septumSomaProfile, delete columnSomaProfile;
}

void LandmarkAnalyzer::countCellsInBarreloids ( const char* outputFilename )
{
	if(!barreloidFlag || !cellFlag)
	{
		std::cout << "Error! Cannot count cells in barreloids. Incomplete input data!" << std::endl;
		return;
	}
	std::cout << "Counting cells in individual barreloids..." << std::endl;
	#ifdef DEBUG
	std::cout << "Checking " << cellSomata->GetNumberOfPoints() << " somata!" << std::endl;
	#endif
	
	std::map< int, PointsPointerType > barreloidFieldSomata;
	double delta[] = {1E-3, 1E-3, 1E-3};
	
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barreloidField.find(ID) != barreloidField.end())
		{
			PointsPointerType barreloidSomata = PointsPointerType::New();
			barreloidSomata->SetDataTypeToFloat();
			
			double bounds[6];
			barreloidField[ID]->ptr()->GetBounds(bounds);
			
			#ifdef DEBUG
			std::cout << "Checking for somata in barreloid " << SBF->int2Labels[ID] << std::endl;
			std::cout << "Bounds: [" << bounds[0] << "," << bounds[1] << "],[";
			std::cout << bounds[2] << "," << bounds[3] << "],[";
			std::cout << bounds[4] << "," << bounds[5] <<"]" << std::endl;
			#endif
			
			unsigned int nrOfSomata = cellSomata->GetNumberOfPoints();
			for(unsigned int ii = 0; ii < nrOfSomata; ++ii)
			{
				double pt[3];
				cellSomata->GetPoint(ii, pt);
				
// 				#ifdef DEBUG
// 				std::cout << "Checking point" << std::endl;
// 				std::cout << "[" << pt[0] << "," << pt[1] << "," << pt[2] << "]" << std::endl;
// 				#endif

				for(int jj = 0; jj < 3; ++jj)
				{
					if(pt[jj]+delta[jj] < bounds[2*jj] || pt[jj]-delta[jj] > bounds[2*jj+1])
						continue;
				}
				
				#ifdef DEBUG
				std::cout << "Found candidate point" << std::endl;
				std::cout << "[" << pt[0] << "," << pt[1] << "," << pt[2] << "]" << std::endl;
				#endif
				
				if(barreloidField[ID]->isPointInsideSurface(pt))
				{
					unsigned int oldNrOfPts = barreloidSomata->GetNumberOfPoints();
					barreloidSomata->InsertPoint(oldNrOfPts, pt);
					#ifdef DEBUG
					std::cout << "Point is inside" << std::endl;
					#endif
				}
			}
			barreloidFieldSomata.insert(std::pair< int, PointsPointerType >(ID, barreloidSomata));
		}
	}
	
	std::string ofName(outputFilename);
	ofName += "_somata.csv";
	std::ofstream SomaList;
	SomaList.open(ofName.c_str());
	SomaList << "# Barreloid\tNr of somata" << std::endl;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barreloidFieldSomata.find(ID) != barreloidFieldSomata.end())
		{
			SomaList << SBF->int2Labels[ID] << "\t" << barreloidFieldSomata[ID]->GetNumberOfPoints() << std::endl;
			std::string somaName(outputFilename);
			somaName += "_somata_";
			somaName += SBF->int2Labels[ID];
			Reader * somaWriter = new Reader(somaName.c_str(), somaName.c_str());
			somaWriter->writeLandmarkFile(barreloidFieldSomata[ID]);
			delete somaWriter;
		}
	}
	SomaList.close();
}

void LandmarkAnalyzer::correctINDensity(ImageDataPointerType density, const char* outputFilename)
{
	if(!density)
	{
		std::cout << "Error! Interneuron density field invalid!" << std::endl;
		return;
	}
	
	Profile * refProfile = createCorrectionProfile();
	double binOffset = 25;
	double refSupra = 525, refGran = 365, refInfra = 1075;
	
	double origin[3], voxelSpacing[3];
	int dimensions[3];
	density->GetOrigin(origin);
	density->GetSpacing(voxelSpacing);
	density->GetDimensions(dimensions);
	for(int kk = 0; kk < dimensions[0]; ++kk)
		for(int ll = 0; ll < dimensions[1]; ++ll)
			for(int mm = 0; mm < dimensions[2]; ++mm)
			{
				int klm[3];
				double xyz[3];
				klm[0] = kk, klm[1] = ll, klm[2] = mm;
				xyz[0] = origin[0] + kk*voxelSpacing[0];
				xyz[1] = origin[1] + ll*voxelSpacing[1];
				xyz[2] = origin[2] + mm*voxelSpacing[2];
				
				int ID = SBF->closestBarrel(xyz);
				double piaDist = SBF->piaDistance(xyz);
				double * densVal = static_cast< double * >(density->GetScalarPointer(klm));
				if(!ID || piaDist == -1 || !(*densVal)/* || !SBF->isInsideS1(xyz)*/)
				{
					*densVal = 0;
					continue;
				}
				
#ifdef DEBUG
std::cout << "klm:\t\t[" << kk << "," << ll << "," << mm << "]" << std::endl;
std::cout << "xyz:\t\t[" << xyz[0] << "," << xyz[1] << "," << xyz[2] << "]" << std::endl;
std::cout << "closest ID:\t" << SBF->int2Labels[ID] << std::endl;
std::cout << "pia distance:\t" << piaDist << std::endl;
std::cout << "densVal:\t" << *densVal << std::endl;
#endif
				
				double colSupra, colGran, colInfra;
				double supraScale, granScale, infraScale;
				colSupra = SBF->avgTopDist[ID];
				colGran = SBF->avgBarrels[ID]->getHeight();
				colInfra = SBF->avgPiaWMDist[ID] - SBF->avgBarrels[ID]->getHeight() - SBF->avgTopDist[ID];
				supraScale = refSupra/colSupra;
				granScale = refGran/colGran;
				infraScale = refInfra/colInfra;
				double refDepth;
				// scale the supragranular part
				if(piaDist <= colSupra)
				{
					refDepth = piaDist*supraScale;
#ifdef DEBUG
std::cout << "supragranular scaling!" << std::endl;
#endif
				}
				// scale the granular part
				else if(piaDist <= colSupra + colGran)
				{
					refDepth = refSupra + (piaDist - colSupra)*granScale;
#ifdef DEBUG
std::cout << "granular scaling!" << std::endl;
#endif
				}
				// scale the infragranular part
				else
				{
					refDepth = refSupra + refGran + (piaDist - colGran - colSupra)*infraScale;
#ifdef DEBUG
std::cout << "infragranular scaling!" << std::endl;
#endif
				}
				// find the scaling bin and interpolate linearly
				// between the two closest bins
				int bin;
				double binDist, binSize;
				double scale;
				binSize = refProfile->getBinSize();
				for(int jj = 0; jj < refProfile->getProfile()->size(); ++jj)
				{
					// assign here in case column is larger than refProfile;
					// then, we automatically scale by 1.0
					bin = jj;
					if((jj+1)*binSize > refDepth)
						break;
				}
				binDist = refDepth - bin*binSize - binOffset;
				if(binDist > 0)
				{
					if(bin < refProfile->getProfile()->size() - 1)
					{
						scale = (1 - binDist/binSize)*refProfile->getProfile()->at(bin) + binDist/binSize*refProfile->getProfile()->at(bin+1);
					}
					else
					{
						scale = refProfile->getProfile()->at(bin);
					}
				}
				if(binDist <= 0)
				{
					if(bin > 0)
					{
						scale = (1 + binDist/binSize)*refProfile->getProfile()->at(bin) - binDist/binSize*refProfile->getProfile()->at(bin-1);
					}
					else
					{
						scale = refProfile->getProfile()->at(bin);
					}
				}
				double val = *densVal;
#ifdef DEBUG
std::cout << "column depth:\t" << piaDist << std::endl;
std::cout << "reference depth:\t" << refDepth << std::endl;
std::cout << "reference bin\t" << bin << std::endl;
std::cout << "bin distance\t" << binDist << std::endl;
std::cout << "old value\t" << val << std::endl;
#endif
				val *= scale;
				*densVal = val;
#ifdef DEBUG
std::cout << "final scale\t" << scale << std::endl;
std::cout << "new value\t" << val << std::endl;
std::cout << "-----------------------------------------------" << std::endl;
#endif
			}
	
	Reader * scaledFieldWriter = new Reader(outputFilename, outputFilename);
	scaledFieldWriter->writeScalarField(density);
	delete scaledFieldWriter;
}

void LandmarkAnalyzer::cropINDensity(ImageDataPointerType density, const char* outputFilename)
{
	ImageDataPointerType CRowDensitySupra = ImageDataPointerType::New();
	CRowDensitySupra->DeepCopy(density);
	ImageDataPointerType CRowDensityGran = ImageDataPointerType::New();
	CRowDensityGran->DeepCopy(density);
	ImageDataPointerType CRowDensityInfra = ImageDataPointerType::New();
	CRowDensityInfra->DeepCopy(density);
	ImageDataPointerType Arc2DensitySupra = ImageDataPointerType::New();
	Arc2DensitySupra->DeepCopy(density);
	ImageDataPointerType Arc2DensityGran = ImageDataPointerType::New();
	Arc2DensityGran->DeepCopy(density);
	ImageDataPointerType Arc2DensityInfra = ImageDataPointerType::New();
	Arc2DensityInfra->DeepCopy(density);
	
	double origin[3], voxelSpacing[3];
	int dimensions[3];
	density->GetOrigin(origin);
	density->GetSpacing(voxelSpacing);
	density->GetDimensions(dimensions);
	
	const double voxelSize = voxelSpacing[0];
	const double voxelVol = voxelSize*voxelSize*voxelSize*1E-9;
	Profile * CRowSupraRawProfile = new Profile(voxelSize);
	Profile * CRowSupraRawColumnProfile = new Profile(voxelSize);
	Profile * CRowSupraRawSeptumProfile = new Profile(voxelSize);
	Profile * CRowSupraVolumeProfile = new Profile(voxelSize);
	Profile * CRowSupraVolumeColumnProfile = new Profile(voxelSize);
	Profile * CRowSupraVolumeSeptumProfile = new Profile(voxelSize);
	Profile * CRowSupraDensityProfile = new Profile(voxelSize);
	Profile * CRowGranRawProfile = new Profile(voxelSize);
	Profile * CRowGranRawColumnProfile = new Profile(voxelSize);
	Profile * CRowGranRawSeptumProfile = new Profile(voxelSize);
	Profile * CRowGranVolumeProfile = new Profile(voxelSize);
	Profile * CRowGranVolumeColumnProfile = new Profile(voxelSize);
	Profile * CRowGranVolumeSeptumProfile = new Profile(voxelSize);
	Profile * CRowGranDensityProfile = new Profile(voxelSize);
	Profile * CRowInfraRawProfile = new Profile(voxelSize);
	Profile * CRowInfraRawColumnProfile = new Profile(voxelSize);
	Profile * CRowInfraRawSeptumProfile = new Profile(voxelSize);
	Profile * CRowInfraVolumeProfile = new Profile(voxelSize);
	Profile * CRowInfraVolumeColumnProfile = new Profile(voxelSize);
	Profile * CRowInfraVolumeSeptumProfile = new Profile(voxelSize);
	Profile * CRowInfraDensityProfile = new Profile(voxelSize);
	
	Profile * Arc2SupraRawProfile = new Profile(voxelSize);
	Profile * Arc2SupraRawColumnProfile = new Profile(voxelSize);
	Profile * Arc2SupraRawSeptumProfile = new Profile(voxelSize);
	Profile * Arc2SupraVolumeProfile = new Profile(voxelSize);
	Profile * Arc2SupraVolumeColumnProfile = new Profile(voxelSize);
	Profile * Arc2SupraVolumeSeptumProfile = new Profile(voxelSize);
	Profile * Arc2SupraDensityProfile = new Profile(voxelSize);
	Profile * Arc2GranRawProfile = new Profile(voxelSize);
	Profile * Arc2GranRawColumnProfile = new Profile(voxelSize);
	Profile * Arc2GranRawSeptumProfile = new Profile(voxelSize);
	Profile * Arc2GranVolumeProfile = new Profile(voxelSize);
	Profile * Arc2GranVolumeColumnProfile = new Profile(voxelSize);
	Profile * Arc2GranVolumeSeptumProfile = new Profile(voxelSize);
	Profile * Arc2GranDensityProfile = new Profile(voxelSize);
	Profile * Arc2InfraRawProfile = new Profile(voxelSize);
	Profile * Arc2InfraRawColumnProfile = new Profile(voxelSize);
	Profile * Arc2InfraRawSeptumProfile = new Profile(voxelSize);
	Profile * Arc2InfraVolumeProfile = new Profile(voxelSize);
	Profile * Arc2InfraVolumeColumnProfile = new Profile(voxelSize);
	Profile * Arc2InfraVolumeSeptumProfile = new Profile(voxelSize);
	Profile * Arc2InfraDensityProfile = new Profile(voxelSize);
	
	std::list< int > CRowLabels, Arc2Labels;
	CRowLabels.push_back(C1);
	CRowLabels.push_back(C2);
	CRowLabels.push_back(C3);
	CRowLabels.push_back(C4);
	Arc2Labels.push_back(A2);
	Arc2Labels.push_back(B2);
	Arc2Labels.push_back(C2);
	Arc2Labels.push_back(D2);
	Arc2Labels.push_back(E2);
	
	for(int kk = 0; kk < dimensions[0]; ++kk)
		for(int ll = 0; ll < dimensions[1]; ++ll)
			for(int mm = 0; mm < dimensions[2]; ++mm)
			{
				int klm[3];
				double xyz[3];
				klm[0] = kk, klm[1] = ll, klm[2] = mm;
				xyz[0] = origin[0] + kk*voxelSpacing[0];
				xyz[1] = origin[1] + ll*voxelSpacing[1];
				xyz[2] = origin[2] + mm*voxelSpacing[2];
				
				int ID = SBF->closestBarrel(xyz);
				int layer = SBF->laminarPosition(xyz);
				
				int rowBin = kk;
				int arcBin = ll;
				
				std::list< int >::const_iterator labelIt;
				
				bool inCRow = false;
				for(labelIt = CRowLabels.begin(); labelIt != CRowLabels.end(); ++labelIt)
				{
					if(*labelIt == ID)
					{
						inCRow = true;
						break;
					}
				}
				
				bool inArc2 = false;
				for(labelIt = Arc2Labels.begin(); labelIt != Arc2Labels.end(); ++labelIt)
				{
					if(*labelIt == ID)
					{
						inArc2 = true;
						break;
					}
				}
				
				if(inCRow && SBF->isInsideS1(xyz))
				{
					switch(layer)
					{
						double * densVal;
						case SUPRA:
							densVal = static_cast< double * >(CRowDensityGran->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(CRowDensityInfra->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(CRowDensitySupra->GetScalarPointer(klm));
							CRowSupraRawProfile->addSegment(*densVal, rowBin);
							CRowSupraVolumeProfile->addSegment(voxelVol, rowBin);
							if(SBF->insideColumn(xyz))
							{
								CRowSupraRawColumnProfile->addSegment(*densVal, rowBin);
								CRowSupraVolumeColumnProfile->addSegment(voxelVol, rowBin);
							}
							else
							{
								CRowSupraRawSeptumProfile->addSegment(*densVal, rowBin);
								CRowSupraVolumeSeptumProfile->addSegment(voxelVol, rowBin);
							}
							break;
						
						case GRAN:
							densVal = static_cast< double * >(CRowDensitySupra->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(CRowDensityInfra->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(CRowDensityGran->GetScalarPointer(klm));
							CRowGranRawProfile->addSegment(*densVal, rowBin);
							CRowGranVolumeProfile->addSegment(voxelVol, rowBin);
							if(SBF->insideColumn(xyz))
							{
								CRowGranRawColumnProfile->addSegment(*densVal, rowBin);
								CRowGranVolumeColumnProfile->addSegment(voxelVol, rowBin);
							}
							else
							{
								CRowGranRawSeptumProfile->addSegment(*densVal, rowBin);
								CRowGranVolumeSeptumProfile->addSegment(voxelVol, rowBin);
							}
							break;
						
						case INFRA:
							densVal = static_cast< double * >(CRowDensitySupra->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(CRowDensityGran->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(CRowDensityInfra->GetScalarPointer(klm));
							CRowInfraRawProfile->addSegment(*densVal, rowBin);
							CRowInfraVolumeProfile->addSegment(voxelVol, rowBin);
							if(SBF->insideColumn(xyz))
							{
								CRowInfraRawColumnProfile->addSegment(*densVal, rowBin);
								CRowInfraVolumeColumnProfile->addSegment(voxelVol, rowBin);
							}
							else
							{
								CRowInfraRawSeptumProfile->addSegment(*densVal, rowBin);
								CRowInfraVolumeSeptumProfile->addSegment(voxelVol, rowBin);
							}
							break;
						
						default:
							densVal = static_cast< double * >(CRowDensitySupra->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(CRowDensityGran->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(CRowDensityInfra->GetScalarPointer(klm));
							*densVal = -1;
					}
				}
				else
				{
					double * densVal = static_cast< double * >(CRowDensitySupra->GetScalarPointer(klm));
					*densVal = -1;
					densVal = static_cast< double * >(CRowDensityGran->GetScalarPointer(klm));
					*densVal = -1;
					densVal = static_cast< double * >(CRowDensityInfra->GetScalarPointer(klm));
					*densVal = -1;
					CRowSupraRawProfile->addSegment(0, rowBin);
					CRowSupraVolumeProfile->addSegment(0, rowBin);
					CRowGranRawProfile->addSegment(0, rowBin);
					CRowGranVolumeProfile->addSegment(0, rowBin);
					CRowInfraRawProfile->addSegment(0, rowBin);
					CRowInfraVolumeProfile->addSegment(0, rowBin);
					CRowGranRawColumnProfile->addSegment(0, rowBin);
					CRowGranVolumeColumnProfile->addSegment(0, rowBin);
					CRowSupraRawSeptumProfile->addSegment(0, rowBin);
					CRowSupraVolumeSeptumProfile->addSegment(0, rowBin);
					CRowGranRawColumnProfile->addSegment(0, rowBin);
					CRowGranVolumeColumnProfile->addSegment(0, rowBin);
					CRowGranRawSeptumProfile->addSegment(0, rowBin);
					CRowGranVolumeSeptumProfile->addSegment(0, rowBin);
					CRowInfraRawColumnProfile->addSegment(0, rowBin);
					CRowInfraVolumeColumnProfile->addSegment(0, rowBin);
					CRowInfraRawSeptumProfile->addSegment(0, rowBin);
					CRowInfraVolumeSeptumProfile->addSegment(0, rowBin);
				}
// 				if(inCRow && !SBF->isInsideS1(xyz))
// 				{
// 					CRowSupraRawProfile->addSegment(0, rowBin);
// 					CRowSupraVolumeProfile->addSegment(0, rowBin);
// 					CRowGranRawProfile->addSegment(0, rowBin);
// 					CRowGranVolumeProfile->addSegment(0, rowBin);
// 					CRowInfraRawProfile->addSegment(0, rowBin);
// 					CRowInfraVolumeProfile->addSegment(0, rowBin);
// 				}
				
				if(inArc2 && SBF->isInsideS1(xyz))
				{
					switch(layer)
					{
						double * densVal;
						case SUPRA:
							densVal = static_cast< double * >(Arc2DensityGran->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(Arc2DensityInfra->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(Arc2DensitySupra->GetScalarPointer(klm));
							Arc2SupraRawProfile->addSegment(*densVal, arcBin);
							Arc2SupraVolumeProfile->addSegment(voxelVol, arcBin);
							if(SBF->insideColumn(xyz))
							{
								Arc2SupraRawColumnProfile->addSegment(*densVal, rowBin);
								Arc2SupraVolumeColumnProfile->addSegment(voxelVol, rowBin);
							}
							else
							{
								Arc2SupraRawSeptumProfile->addSegment(*densVal, rowBin);
								Arc2SupraVolumeSeptumProfile->addSegment(voxelVol, rowBin);
							}
							break;
						
						case GRAN:
							densVal = static_cast< double * >(Arc2DensitySupra->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(Arc2DensityInfra->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(Arc2DensityGran->GetScalarPointer(klm));
							Arc2GranRawProfile->addSegment(*densVal, arcBin);
							Arc2GranVolumeProfile->addSegment(voxelVol, arcBin);
							if(SBF->insideColumn(xyz))
							{
								Arc2GranRawColumnProfile->addSegment(*densVal, rowBin);
								Arc2GranVolumeColumnProfile->addSegment(voxelVol, rowBin);
							}
							else
							{
								Arc2GranRawSeptumProfile->addSegment(*densVal, rowBin);
								Arc2GranVolumeSeptumProfile->addSegment(voxelVol, rowBin);
							}
							break;
						
						case INFRA:
							densVal = static_cast< double * >(Arc2DensitySupra->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(Arc2DensityGran->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(Arc2DensityInfra->GetScalarPointer(klm));
							Arc2InfraRawProfile->addSegment(*densVal, arcBin);
							Arc2InfraVolumeProfile->addSegment(voxelVol, arcBin);
							if(SBF->insideColumn(xyz))
							{
								Arc2InfraRawColumnProfile->addSegment(*densVal, rowBin);
								Arc2InfraVolumeColumnProfile->addSegment(voxelVol, rowBin);
							}
							else
							{
								Arc2InfraRawSeptumProfile->addSegment(*densVal, rowBin);
								Arc2InfraVolumeSeptumProfile->addSegment(voxelVol, rowBin);
							}
							break;
						
						default:
							densVal = static_cast< double * >(Arc2DensitySupra->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(Arc2DensityGran->GetScalarPointer(klm));
							*densVal = -1;
							densVal = static_cast< double * >(Arc2DensityInfra->GetScalarPointer(klm));
							*densVal = -1;
					}
				}
				else
				{
					double * densVal = static_cast< double * >(Arc2DensitySupra->GetScalarPointer(klm));
					*densVal = -1;
					densVal = static_cast< double * >(Arc2DensityGran->GetScalarPointer(klm));
					*densVal = -1;
					densVal = static_cast< double * >(Arc2DensityInfra->GetScalarPointer(klm));
					*densVal = -1;
					Arc2SupraRawProfile->addSegment(0, arcBin);
					Arc2SupraVolumeProfile->addSegment(0, arcBin);
					Arc2GranRawProfile->addSegment(0, arcBin);
					Arc2GranVolumeProfile->addSegment(0, arcBin);
					Arc2InfraRawProfile->addSegment(0, arcBin);
					Arc2InfraVolumeProfile->addSegment(0, arcBin);
					Arc2GranRawColumnProfile->addSegment(0, rowBin);
					Arc2GranVolumeColumnProfile->addSegment(0, rowBin);
					Arc2SupraRawSeptumProfile->addSegment(0, rowBin);
					Arc2SupraVolumeSeptumProfile->addSegment(0, rowBin);
					Arc2GranRawColumnProfile->addSegment(0, rowBin);
					Arc2GranVolumeColumnProfile->addSegment(0, rowBin);
					Arc2GranRawSeptumProfile->addSegment(0, rowBin);
					Arc2GranVolumeSeptumProfile->addSegment(0, rowBin);
					Arc2InfraRawColumnProfile->addSegment(0, rowBin);
					Arc2InfraVolumeColumnProfile->addSegment(0, rowBin);
					Arc2InfraRawSeptumProfile->addSegment(0, rowBin);
					Arc2InfraVolumeSeptumProfile->addSegment(0, rowBin);
				}
// 				if(inArc2 && !SBF->isInsideS1(xyz))
// 				{
// 					Arc2SupraRawProfile->addSegment(0, arcBin);
// 					Arc2SupraVolumeProfile->addSegment(0, arcBin);
// 					Arc2GranRawProfile->addSegment(0, arcBin);
// 					Arc2GranVolumeProfile->addSegment(0, arcBin);
// 					Arc2InfraRawProfile->addSegment(0, arcBin);
// 					Arc2InfraVolumeProfile->addSegment(0, arcBin);
// 				}
			}
	
	computeDensityProfile(CRowSupraRawProfile, CRowSupraVolumeProfile, CRowSupraDensityProfile);
	computeDensityProfile(CRowGranRawProfile, CRowGranVolumeProfile, CRowGranDensityProfile);
	computeDensityProfile(CRowInfraRawProfile, CRowInfraVolumeProfile, CRowInfraDensityProfile);
	computeDensityProfile(Arc2SupraRawProfile, Arc2SupraVolumeProfile, Arc2SupraDensityProfile);
	computeDensityProfile(Arc2GranRawProfile, Arc2GranVolumeProfile, Arc2GranDensityProfile);
	computeDensityProfile(Arc2InfraRawProfile, Arc2InfraVolumeProfile, Arc2InfraDensityProfile);
	
	double CRowSupraColumnDensity, CRowGranColumnDensity, CRowInfraColumnDensity;
	double CRowSupraSeptumDensity, CRowGranSeptumDensity, CRowInfraSeptumDensity;
	double Arc2SupraColumnDensity, Arc2GranColumnDensity, Arc2InfraColumnDensity;
	double Arc2SupraSeptumDensity, Arc2GranSeptumDensity, Arc2InfraSeptumDensity;
	CRowSupraColumnDensity = CRowSupraRawColumnProfile->getIntegral()/CRowSupraVolumeColumnProfile->getIntegral();
	CRowSupraSeptumDensity = CRowSupraRawSeptumProfile->getIntegral()/CRowSupraVolumeSeptumProfile->getIntegral();
	CRowGranColumnDensity = CRowGranRawColumnProfile->getIntegral()/CRowGranVolumeColumnProfile->getIntegral();
	CRowGranSeptumDensity = CRowGranRawSeptumProfile->getIntegral()/CRowGranVolumeSeptumProfile->getIntegral();
	CRowInfraColumnDensity = CRowInfraRawColumnProfile->getIntegral()/CRowInfraVolumeColumnProfile->getIntegral();
	CRowInfraSeptumDensity = CRowInfraRawSeptumProfile->getIntegral()/CRowInfraVolumeSeptumProfile->getIntegral();
	Arc2SupraColumnDensity = Arc2SupraRawColumnProfile->getIntegral()/Arc2SupraVolumeColumnProfile->getIntegral();
	Arc2SupraSeptumDensity = Arc2SupraRawSeptumProfile->getIntegral()/Arc2SupraVolumeSeptumProfile->getIntegral();
	Arc2GranColumnDensity = Arc2GranRawColumnProfile->getIntegral()/Arc2GranVolumeColumnProfile->getIntegral();
	Arc2GranSeptumDensity = Arc2GranRawSeptumProfile->getIntegral()/Arc2GranVolumeSeptumProfile->getIntegral();
	Arc2InfraColumnDensity = Arc2InfraRawColumnProfile->getIntegral()/Arc2InfraVolumeColumnProfile->getIntegral();
	Arc2InfraSeptumDensity = Arc2InfraRawSeptumProfile->getIntegral()/Arc2InfraVolumeSeptumProfile->getIntegral();
	
	std::string summaryName(outputFilename);
	summaryName += "_density_averages.csv";
	std::ofstream summaryWriter;
	summaryWriter.open(summaryName.c_str());
	summaryWriter << "#bouton density [1/mm^3]" << std::endl;
	summaryWriter << "row/arc\tcol/sep\tsupra\tgran\tinfra" << std::endl;
	summaryWriter << "C row\tcol\t" << CRowSupraColumnDensity << "\t" << CRowGranColumnDensity << "\t"<< CRowInfraColumnDensity << std::endl;
	summaryWriter << "C row\tsep\t" << CRowSupraSeptumDensity << "\t" << CRowGranSeptumDensity << "\t"<< CRowInfraSeptumDensity << std::endl;
	summaryWriter << "arc 2\tcol\t" << Arc2SupraColumnDensity << "\t" << Arc2GranColumnDensity << "\t"<< Arc2InfraColumnDensity << std::endl;
	summaryWriter << "arc 2\tsep\t" << Arc2SupraSeptumDensity << "\t" << Arc2GranSeptumDensity << "\t"<< Arc2InfraSeptumDensity << std::endl;
	summaryWriter.close();
	
	std::string CRowNameSupra(outputFilename);
	CRowNameSupra += "_CRowSupra";
	Reader * CRowWriterSupra = new Reader(CRowNameSupra.c_str(), CRowNameSupra.c_str());
	CRowWriterSupra->writeScalarField(CRowDensitySupra);
	std::string CRowNameSupraDens(CRowNameSupra);
	CRowNameSupraDens += "_density";
	writeZProfile(CRowSupraDensityProfile, CRowNameSupraDens.c_str(), origin[0]);
	std::string CRowNameSupraLength(CRowNameSupra);
	CRowNameSupraLength += "_length";
	writeZProfile(CRowSupraRawProfile, CRowNameSupraLength.c_str(), origin[0]);
	std::string CRowNameGran(outputFilename);
	CRowNameGran += "_CRowGran";
	Reader * CRowWriterGran = new Reader(CRowNameGran.c_str(), CRowNameGran.c_str());
	CRowWriterGran->writeScalarField(CRowDensityGran);
	std::string CRowNameGranDens(CRowNameGran);
	CRowNameGranDens += "_density";
	writeZProfile(CRowGranDensityProfile, CRowNameGranDens.c_str(), origin[0]);
	std::string CRowNameGranLength(CRowNameGran);
	CRowNameGranLength += "_length";
	writeZProfile(CRowGranRawProfile, CRowNameGranLength.c_str(), origin[0]);
	std::string CRowNameInfra(outputFilename);
	CRowNameInfra += "_CRowInfra";
	Reader * CRowWriterInfra = new Reader(CRowNameInfra.c_str(), CRowNameInfra.c_str());
	CRowWriterInfra->writeScalarField(CRowDensityInfra);
	std::string CRowNameInfraDens(CRowNameInfra);
	CRowNameInfraDens += "_density";
	writeZProfile(CRowInfraDensityProfile, CRowNameInfraDens.c_str(), origin[0]);
	std::string CRowNameInfraLength(CRowNameInfra);
	CRowNameInfraLength += "_length";
	writeZProfile(CRowInfraRawProfile, CRowNameInfraLength.c_str(), origin[0]);
	
	std::string Arc2NameSupra(outputFilename);
	Arc2NameSupra += "_Arc2Supra";
	Reader * Arc2WriterSupra = new Reader(Arc2NameSupra.c_str(), Arc2NameSupra.c_str());
	Arc2WriterSupra->writeScalarField(Arc2DensitySupra);
	std::string Arc2NameSupraDens(Arc2NameSupra);
	Arc2NameSupraDens += "_density";
	writeZProfile(Arc2SupraDensityProfile, Arc2NameSupraDens.c_str(), origin[1]);
	std::string Arc2NameSupraLength(Arc2NameSupra);
	Arc2NameSupraLength += "_length";
	writeZProfile(Arc2SupraRawProfile, Arc2NameSupraLength.c_str(), origin[1]);
	std::string Arc2NameGran(outputFilename);
	Arc2NameGran += "_Arc2Gran";
	Reader * Arc2WriterGran = new Reader(Arc2NameGran.c_str(), Arc2NameGran.c_str());
	Arc2WriterGran->writeScalarField(Arc2DensityGran);
	std::string Arc2NameGranDens(Arc2NameGran);
	Arc2NameGranDens += "_density";
	writeZProfile(Arc2GranDensityProfile, Arc2NameGranDens.c_str(), origin[1]);
	std::string Arc2NameGranLength(Arc2NameGran);
	Arc2NameGranLength += "_length";
	writeZProfile(Arc2GranRawProfile, Arc2NameGranLength.c_str(), origin[1]);
	std::string Arc2NameInfra(outputFilename);
	Arc2NameInfra += "_Arc2Infra";
	Reader * Arc2WriterInfra = new Reader(Arc2NameInfra.c_str(), Arc2NameInfra.c_str());
	Arc2WriterInfra->writeScalarField(Arc2DensityInfra);
	std::string Arc2NameInfraDens(Arc2NameInfra);
	Arc2NameInfraDens += "_density";
	writeZProfile(Arc2InfraDensityProfile, Arc2NameInfraDens.c_str(), origin[1]);
	std::string Arc2NameInfraLength(Arc2NameInfra);
	Arc2NameInfraLength += "_length";
	writeZProfile(Arc2InfraRawProfile, Arc2NameInfraLength.c_str(), origin[1]);
	
	delete CRowWriterSupra, delete CRowWriterGran, delete CRowWriterInfra;
	delete Arc2WriterSupra, delete Arc2WriterGran, delete Arc2WriterInfra;
	
	delete CRowSupraRawProfile, delete CRowSupraVolumeProfile, delete CRowSupraDensityProfile;
	delete CRowGranRawProfile, delete CRowGranVolumeProfile, delete CRowGranDensityProfile;
	delete CRowInfraRawProfile, delete CRowInfraVolumeProfile, delete CRowInfraDensityProfile;
	delete Arc2SupraRawProfile, delete Arc2SupraVolumeProfile, delete Arc2SupraDensityProfile;
	delete Arc2GranRawProfile, delete Arc2GranVolumeProfile, delete Arc2GranDensityProfile;
	delete Arc2InfraRawProfile, delete Arc2InfraVolumeProfile, delete Arc2InfraDensityProfile;
	delete CRowSupraRawColumnProfile;
	delete CRowSupraRawSeptumProfile;
	delete CRowSupraVolumeColumnProfile;
	delete CRowSupraVolumeSeptumProfile;
	delete CRowGranRawColumnProfile;
	delete CRowGranRawSeptumProfile;
	delete CRowGranVolumeColumnProfile;
	delete CRowGranVolumeSeptumProfile;
	delete CRowInfraRawColumnProfile;
	delete CRowInfraRawSeptumProfile;
	delete CRowInfraVolumeColumnProfile;
	delete CRowInfraVolumeSeptumProfile;
	delete Arc2SupraRawColumnProfile;
	delete Arc2SupraRawSeptumProfile;
	delete Arc2SupraVolumeColumnProfile;
	delete Arc2SupraVolumeSeptumProfile;
	delete Arc2GranRawColumnProfile;
	delete Arc2GranRawSeptumProfile;
	delete Arc2GranVolumeColumnProfile;
	delete Arc2GranVolumeSeptumProfile;
	delete Arc2InfraRawColumnProfile;
	delete Arc2InfraRawSeptumProfile;
	delete Arc2InfraVolumeColumnProfile;
	delete Arc2InfraVolumeSeptumProfile;
}

void LandmarkAnalyzer::cropDensityLayers(ImageDataPointerType density, const char* outputFilename)
{
	ImageDataPointerType DensitySupra = ImageDataPointerType::New();
	DensitySupra->DeepCopy(density);
	ImageDataPointerType DensityGran = ImageDataPointerType::New();
	DensityGran->DeepCopy(density);
	ImageDataPointerType DensityInfra = ImageDataPointerType::New();
	DensityInfra->DeepCopy(density);
	
	double origin[3], voxelSpacing[3];
	int dimensions[3];
	density->GetOrigin(origin);
	density->GetSpacing(voxelSpacing);
	density->GetDimensions(dimensions);
	
	for(int kk = 0; kk < dimensions[0]; ++kk)
		for(int ll = 0; ll < dimensions[1]; ++ll)
			for(int mm = 0; mm < dimensions[2]; ++mm)
			{
				int klm[3];
				double xyz[3];
				klm[0] = kk, klm[1] = ll, klm[2] = mm;
				xyz[0] = origin[0] + kk*voxelSpacing[0];
				xyz[1] = origin[1] + ll*voxelSpacing[1];
				xyz[2] = origin[2] + mm*voxelSpacing[2];
				
				int layer = SBF->laminarPosition(xyz);
				
				if(SBF->isInsideS1(xyz))
				{
					switch(layer)
					{
						double * densVal;
						case SUPRA:
							densVal = static_cast< double * >(DensityGran->GetScalarPointer(klm));
							*densVal = 0;
							densVal = static_cast< double * >(DensityInfra->GetScalarPointer(klm));
							*densVal = 0;
							break;
						
						case GRAN:
							densVal = static_cast< double * >(DensitySupra->GetScalarPointer(klm));
							*densVal = 0;
							densVal = static_cast< double * >(DensityInfra->GetScalarPointer(klm));
							*densVal = 0;
							break;
						
						case INFRA:
							densVal = static_cast< double * >(DensitySupra->GetScalarPointer(klm));
							*densVal = 0;
							densVal = static_cast< double * >(DensityGran->GetScalarPointer(klm));
							*densVal = 0;
							break;
						
						default:
							densVal = static_cast< double * >(DensitySupra->GetScalarPointer(klm));
							*densVal = 0;
							densVal = static_cast< double * >(DensityGran->GetScalarPointer(klm));
							*densVal = 0;
							densVal = static_cast< double * >(DensityInfra->GetScalarPointer(klm));
							*densVal = 0;
					}
				}
				else
				{
					double * densVal = static_cast< double * >(DensitySupra->GetScalarPointer(klm));
					*densVal = 0;
					densVal = static_cast< double * >(DensityGran->GetScalarPointer(klm));
					*densVal = 0;
					densVal = static_cast< double * >(DensityInfra->GetScalarPointer(klm));
					*densVal = 0;
				}
			}
	
	std::string NameSupra(outputFilename);
	NameSupra += "_Supra";
	Reader * WriterSupra = new Reader(NameSupra.c_str(), NameSupra.c_str());
	WriterSupra->writeScalarField(DensitySupra);
	std::string NameGran(outputFilename);
	NameGran += "_Gran";
	Reader * WriterGran = new Reader(NameGran.c_str(), NameGran.c_str());
	WriterGran->writeScalarField(DensityGran);
	std::string NameInfra(outputFilename);
	NameInfra += "_Infra";
	Reader * WriterInfra = new Reader(NameInfra.c_str(), NameInfra.c_str());
	WriterInfra->writeScalarField(DensityInfra);
	
	delete WriterSupra, delete WriterGran, delete WriterInfra;
}

void LandmarkAnalyzer::cropDensityCRowArc2(ImageDataPointerType density, const char* outputFilename)
{
	ImageDataPointerType CRowDensitySupra = ImageDataPointerType::New();
	CRowDensitySupra->DeepCopy(density);
	ImageDataPointerType Arc2DensitySupra = ImageDataPointerType::New();
	Arc2DensitySupra->DeepCopy(density);
	
	double origin[3], voxelSpacing[3];
	int dimensions[3];
	density->GetOrigin(origin);
	density->GetSpacing(voxelSpacing);
	density->GetDimensions(dimensions);
	
	const double voxelSize = voxelSpacing[0];
	const double voxelVol = voxelSize*voxelSize*voxelSize*1E-9;
	Profile * CRowSupraRawProfile = new Profile(voxelSize);
	Profile * CRowSupraRawColumnProfile = new Profile(voxelSize);
	Profile * CRowSupraRawSeptumProfile = new Profile(voxelSize);
	Profile * CRowSupraVolumeProfile = new Profile(voxelSize);
	Profile * CRowSupraVolumeColumnProfile = new Profile(voxelSize);
	Profile * CRowSupraVolumeSeptumProfile = new Profile(voxelSize);
	Profile * CRowSupraDensityProfile = new Profile(voxelSize);
	
	Profile * Arc2SupraRawProfile = new Profile(voxelSize);
	Profile * Arc2SupraRawColumnProfile = new Profile(voxelSize);
	Profile * Arc2SupraRawSeptumProfile = new Profile(voxelSize);
	Profile * Arc2SupraVolumeProfile = new Profile(voxelSize);
	Profile * Arc2SupraVolumeColumnProfile = new Profile(voxelSize);
	Profile * Arc2SupraVolumeSeptumProfile = new Profile(voxelSize);
	Profile * Arc2SupraDensityProfile = new Profile(voxelSize);
	
	std::list< int > CRowLabels, Arc2Labels;
	CRowLabels.push_back(C1);
	CRowLabels.push_back(C2);
	CRowLabels.push_back(C3);
	CRowLabels.push_back(C4);
	Arc2Labels.push_back(A2);
	Arc2Labels.push_back(B2);
	Arc2Labels.push_back(C2);
	Arc2Labels.push_back(D2);
	Arc2Labels.push_back(E2);
	
	for(int kk = 0; kk < dimensions[0]; ++kk)
		for(int ll = 0; ll < dimensions[1]; ++ll)
			for(int mm = 0; mm < dimensions[2]; ++mm)
			{
				int klm[3];
				double xyz[3];
				klm[0] = kk, klm[1] = ll, klm[2] = mm;
				xyz[0] = origin[0] + kk*voxelSpacing[0];
				xyz[1] = origin[1] + ll*voxelSpacing[1];
				xyz[2] = origin[2] + mm*voxelSpacing[2];
				
				int ID = SBF->closestBarrel(xyz);
				
				int rowBin = kk;
				int arcBin = ll;
				
				std::list< int >::const_iterator labelIt;
				
				bool inCRow = false;
				for(labelIt = CRowLabels.begin(); labelIt != CRowLabels.end(); ++labelIt)
				{
					if(*labelIt == ID)
					{
						inCRow = true;
						break;
					}
				}
				
				bool inArc2 = false;
				for(labelIt = Arc2Labels.begin(); labelIt != Arc2Labels.end(); ++labelIt)
				{
					if(*labelIt == ID)
					{
						inArc2 = true;
						break;
					}
				}
				
				if(inCRow && SBF->isInsideS1(xyz))
				{
					double * densVal;
					densVal = static_cast< double * >(CRowDensitySupra->GetScalarPointer(klm));
					CRowSupraRawProfile->addSegment(*densVal, rowBin);
					CRowSupraVolumeProfile->addSegment(voxelVol, rowBin);
					if(SBF->insideColumn(xyz))
					{
						CRowSupraRawColumnProfile->addSegment(*densVal, rowBin);
						CRowSupraVolumeColumnProfile->addSegment(voxelVol, rowBin);
					}
					else
					{
						CRowSupraRawSeptumProfile->addSegment(*densVal, rowBin);
						CRowSupraVolumeSeptumProfile->addSegment(voxelVol, rowBin);
					}
				}
				else
				{
					double * densVal = static_cast< double * >(CRowDensitySupra->GetScalarPointer(klm));
					*densVal = -1;
				}
// 				if(inCRow && !SBF->isInsideS1(xyz))
// 				{
// 					CRowSupraRawProfile->addSegment(0, rowBin);
// 					CRowSupraVolumeProfile->addSegment(0, rowBin);
// 					CRowGranRawProfile->addSegment(0, rowBin);
// 					CRowGranVolumeProfile->addSegment(0, rowBin);
// 					CRowInfraRawProfile->addSegment(0, rowBin);
// 					CRowInfraVolumeProfile->addSegment(0, rowBin);
// 				}
				
				if(inArc2 && SBF->isInsideS1(xyz))
				{
					double * densVal;
					densVal = static_cast< double * >(Arc2DensitySupra->GetScalarPointer(klm));
					Arc2SupraRawProfile->addSegment(*densVal, arcBin);
					Arc2SupraVolumeProfile->addSegment(voxelVol, arcBin);
					if(SBF->insideColumn(xyz))
					{
						Arc2SupraRawColumnProfile->addSegment(*densVal, rowBin);
						Arc2SupraVolumeColumnProfile->addSegment(voxelVol, rowBin);
					}
					else
					{
						Arc2SupraRawSeptumProfile->addSegment(*densVal, rowBin);
						Arc2SupraVolumeSeptumProfile->addSegment(voxelVol, rowBin);
					}
				}
				else
				{
					double * densVal = static_cast< double * >(Arc2DensitySupra->GetScalarPointer(klm));
					*densVal = -1;
					Arc2SupraRawProfile->addSegment(0, arcBin);
					Arc2SupraVolumeProfile->addSegment(0, arcBin);
					Arc2SupraRawSeptumProfile->addSegment(0, rowBin);
					Arc2SupraVolumeSeptumProfile->addSegment(0, rowBin);
				}
// 				if(inArc2 && !SBF->isInsideS1(xyz))
// 				{
// 					Arc2SupraRawProfile->addSegment(0, arcBin);
// 					Arc2SupraVolumeProfile->addSegment(0, arcBin);
// 					Arc2GranRawProfile->addSegment(0, arcBin);
// 					Arc2GranVolumeProfile->addSegment(0, arcBin);
// 					Arc2InfraRawProfile->addSegment(0, arcBin);
// 					Arc2InfraVolumeProfile->addSegment(0, arcBin);
// 				}
			}
	
	computeDensityProfile(CRowSupraRawProfile, CRowSupraVolumeProfile, CRowSupraDensityProfile);
	computeDensityProfile(Arc2SupraRawProfile, Arc2SupraVolumeProfile, Arc2SupraDensityProfile);
	
	double CRowSupraColumnDensity;
	double CRowSupraSeptumDensity;
	double Arc2SupraColumnDensity;
	double Arc2SupraSeptumDensity;
	CRowSupraColumnDensity = CRowSupraRawColumnProfile->getIntegral()/CRowSupraVolumeColumnProfile->getIntegral();
	CRowSupraSeptumDensity = CRowSupraRawSeptumProfile->getIntegral()/CRowSupraVolumeSeptumProfile->getIntegral();
	Arc2SupraColumnDensity = Arc2SupraRawColumnProfile->getIntegral()/Arc2SupraVolumeColumnProfile->getIntegral();
	Arc2SupraSeptumDensity = Arc2SupraRawSeptumProfile->getIntegral()/Arc2SupraVolumeSeptumProfile->getIntegral();
	
	std::string summaryName(outputFilename);
	summaryName += "_density_averages.csv";
	std::ofstream summaryWriter;
	summaryWriter.open(summaryName.c_str());
	summaryWriter << "#bouton density [1/mm^3]" << std::endl;
	summaryWriter << "row/arc\tcol/sep" << std::endl;
	summaryWriter << "C row\tcol\t" << CRowSupraColumnDensity << std::endl;
	summaryWriter << "C row\tsep\t" << CRowSupraSeptumDensity << std::endl;
	summaryWriter << "arc 2\tcol\t" << Arc2SupraColumnDensity << std::endl;
	summaryWriter << "arc 2\tsep\t" << Arc2SupraSeptumDensity << std::endl;
	summaryWriter.close();
	
	std::string CRowNameSupra(outputFilename);
	CRowNameSupra += "_CRow";
	Reader * CRowWriterSupra = new Reader(CRowNameSupra.c_str(), CRowNameSupra.c_str());
	CRowWriterSupra->writeScalarField(CRowDensitySupra);
	std::string CRowNameSupraDens(CRowNameSupra);
	CRowNameSupraDens += "_density";
	writeZProfile(CRowSupraDensityProfile, CRowNameSupraDens.c_str(), origin[0]);
	std::string CRowNameSupraLength(CRowNameSupra);
	CRowNameSupraLength += "_length";
	writeZProfile(CRowSupraRawProfile, CRowNameSupraLength.c_str(), origin[0]);
	
	std::string Arc2NameSupra(outputFilename);
	Arc2NameSupra += "_Arc2";
	Reader * Arc2WriterSupra = new Reader(Arc2NameSupra.c_str(), Arc2NameSupra.c_str());
	Arc2WriterSupra->writeScalarField(Arc2DensitySupra);
	std::string Arc2NameSupraDens(Arc2NameSupra);
	Arc2NameSupraDens += "_density";
	writeZProfile(Arc2SupraDensityProfile, Arc2NameSupraDens.c_str(), origin[1]);
	std::string Arc2NameSupraLength(Arc2NameSupra);
	Arc2NameSupraLength += "_length";
	writeZProfile(Arc2SupraRawProfile, Arc2NameSupraLength.c_str(), origin[1]);
	
	delete CRowWriterSupra;
	delete Arc2WriterSupra;
	
	delete CRowSupraRawProfile, delete CRowSupraVolumeProfile, delete CRowSupraDensityProfile;
	delete Arc2SupraRawProfile, delete Arc2SupraVolumeProfile, delete Arc2SupraDensityProfile;
	delete CRowSupraRawColumnProfile;
	delete CRowSupraRawSeptumProfile;
	delete CRowSupraVolumeColumnProfile;
	delete CRowSupraVolumeSeptumProfile;
	delete Arc2SupraRawColumnProfile;
	delete Arc2SupraRawSeptumProfile;
	delete Arc2SupraVolumeColumnProfile;
	delete Arc2SupraVolumeSeptumProfile;
}


/*************************************************************************/
/* PRIVATE METHODS                                                       */
/*************************************************************************/

void LandmarkAnalyzer::assignSomataToColumns(std::map< int, PointsPointerType >& somaColumns,
						std::map< int, Column * >& barrelColumns,
						std::map< int, std::vector< std::vector< double > > >& radialContours,
						std::map< int, double >& minDistances,
						std::map< int, double >& maxDistances)
{
	// process every cell soma
	// and assign landmarks to nearest
	// column if in doubt
	#ifdef DEBUG
	std::cout << "Checking " << cellSomata->GetNumberOfPoints() << " somata!" << std::endl;
	#endif
	for(unsigned int ii = 0; ii < cellSomata->GetNumberOfPoints(); ++ii)
	{
		double somaLocation[3];
		cellSomata->GetPoint(ii, somaLocation);
		
		// check distance to all columns
		double minDistance = 1E09;
		int closestColumn = 0;
		std::map< int, Column * >::const_iterator columnIt;
		for(columnIt = barrelColumns.begin(); columnIt != barrelColumns.end(); ++columnIt)
		{
			int ID = columnIt->first;
			double t, projectedPt[3];
			double dist = vtkLine::DistanceToLine(somaLocation, barrelColumns[ID]->top, barrelColumns[ID]->bottom, t, projectedPt);
			dist = sqrt(dist);
			if(dist > maxDistances[ID] || t < 0 || t > 1)
				continue;
			if(dist < minDistances[ID])
			{
				if(dist < minDistance)
				{
					minDistance = dist;
					closestColumn = ID;
				}
				continue;
			}
			
			PolyDataPointerType polyData = PolyDataPointerType::New();
			polyData->Allocate(1);
			PointsPointerType points = PointsPointerType::New();
			points->SetDataTypeToFloat();
			PolygonPointerType poly = PolygonPointerType::New();
			poly->GetPointIds()->SetNumberOfIds(radialContours[ID].size());
			for(int jj = 0; jj < radialContours[ID].size(); ++jj)
			{
				double tmp[3];
				for(int kk = 0; kk < 3; ++kk)
					tmp[kk] = projectedPt[kk] + radialContours[ID][jj][kk];
				points->InsertNextPoint(tmp);
				poly->GetPointIds()->SetId(jj, jj);
			}
			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
			polyData->SetPoints(points);
			polyData->Update();
			double closestPoint[3], pCoords[3], dist2;
			int subId;
			double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
			int insidePolygon = polyData->GetCell(0)->EvaluatePosition(somaLocation, closestPoint, subId, pCoords, dist2, weights);
			if(insidePolygon == 1)
			{
				if(dist < minDistance)
				{
					minDistance = dist;
					closestColumn = ID;
				}
			}
			delete [] weights;
		}
		// assign landmark to closest soma column
		if(closestColumn)
		{
			int oldNrOfPoints = somaColumns[closestColumn]->GetNumberOfPoints();
			somaColumns[closestColumn]->InsertPoint(oldNrOfPoints, somaLocation);
		}
	}
}
// assign somata to 50um bins along the column axis
// defined by top and bottom
Profile * LandmarkAnalyzer::computeColumnSomaProfile ( PointsPointerType somata, double top[3], double bottom[3] )
{
	double binSize = 50;
	Profile * zProfile = new Profile(binSize);
	unsigned int nrOfSomata = somata->GetNumberOfPoints();
	#ifdef DEBUG
	std::flush(std::cout << "Computing soma z profile of " << nrOfSomata << " somata...");
	#endif
	for(int ii = 0; ii < nrOfSomata; ++ii)
	{
		double soma[3], tmpProj[3], t;
		somata->GetPoint(ii, soma);
		vtkLine::DistanceToLine(soma, top, bottom, t, tmpProj);
		double depth1 = sqrt(vtkMath::Distance2BetweenPoints(tmpProj, top));
		if(t >= 0 && t <= 1)
		{
			unsigned int bin = (unsigned int)(depth1/binSize);	// round down
			zProfile->incrementBin(bin);
		}
		else
			std::cout << "Warning! Found somata outside of column during calculation of z profile." << std::endl;
	}
	#ifdef DEBUG
	std::flush(std::cout << " done!" << std::endl);
	#endif
	
	return zProfile;
}

// assign voxels to 50um bins along the column axis
// and choose nearest column in ambiguous cases
Profile* LandmarkAnalyzer::computeColumnVolumeProfile ( std::map< int, Column* > barrelColumns, std::map< int, double > minDistances,
							std::map< int, double > maxDistances, 
							std::map< int, std::vector< std::vector< double > > > radialContours, int ID )
{
	std::cout << "Computing voxel z profile of column " << SBF->int2Labels[ID] << std::endl;
	if(barrelColumns.find(ID) != barrelColumns.end())
	{
		double binSize = 50;
		Profile * zProfile = new Profile(binSize);
		double spacing = SPACING;
		
		// Compute transformation matrix
		// in/out of column-aligned coordinates
		// to prevent discretization artifacts
		// along column axis
		HomogeneousMatrixPointerType mGlobalToColumn = HomogeneousMatrixPointerType::New();
		HomogeneousMatrixPointerType mColumnToGlobal = HomogeneousMatrixPointerType::New();
		double columnCenter[3], columnAxis[3], invColumnCenter[3], tmp[3], t, zAxis[] = {0,0,1};
		vtkMath::Add(barrelColumns[ID]->top, barrelColumns[ID]->bottom, tmp);
		vtkMath::MultiplyScalar(tmp, 0.5);
		vtkLine::DistanceToLine(tmp, barrelColumns[ID]->top, barrelColumns[ID]->bottom, t, columnCenter);
		vtkMath::Subtract(barrelColumns[ID]->top, barrelColumns[ID]->bottom, columnAxis);
		vtkMath::Normalize(columnAxis);
		for(int ii = 0; ii < 3; ++ii)
			invColumnCenter[ii] = -1*columnCenter[ii];
		
		if(columnAxis[2] == 1)
		{
			mGlobalToColumn->Identity();
			mColumnToGlobal->Identity();
		}
		else
		{
			double rotAxis[3], rotAngle;
			vtkMath::Cross(zAxis, columnAxis, rotAxis);
			vtkMath::Normalize(rotAxis);
			rotAngle = acos(vtkMath::Dot(columnAxis, zAxis))*180/PI;
			
			TransformPointerType trans1 = TransformPointerType::New();
			TransformPointerType trans2 = TransformPointerType::New();
			TransformPointerType rot1 = TransformPointerType::New();
			TransformPointerType rot2 = TransformPointerType::New();
			TransformPointerType invTrans1 = TransformPointerType::New();
			TransformPointerType invTrans2 = TransformPointerType::New();
			trans1->Translate(invColumnCenter);
			trans2->Translate(invColumnCenter);
			rot1->RotateWXYZ(rotAngle, rotAxis);
			rot2->RotateWXYZ(-1*rotAngle, rotAxis);
			invTrans1->Translate(columnCenter);
			invTrans2->Translate(columnCenter);
			
			rot1->Concatenate(trans1);
			invTrans1->Concatenate(rot1);
			invTrans1->Update();
			mGlobalToColumn = invTrans1->GetMatrix();
			
			rot2->Concatenate(trans2);
			invTrans2->Concatenate(rot2);
			invTrans2->Update();
			mColumnToGlobal = invTrans2->GetMatrix();
		}
		
		// careful: barrelColumn contours are really
		// still only the barrel contours
		// compute column bounds explicitly
		double colBounds[6];
		for(int ii = 0; ii < 3; ++ii)
		{
			colBounds[2*ii] = 1E9;
			colBounds[2*ii+1] = 1E-9;
		}
		std::vector< std::vector< double > >::const_iterator radialContourIt;
		for(radialContourIt = radialContours[ID].begin(); radialContourIt != radialContours[ID].end(); ++radialContourIt)
		{
			double hTop[4], hBottom[4];
			for(int ii = 0; ii < 3; ++ii)
			{
				hTop[ii] = (*radialContourIt)[ii] + barrelColumns[ID]->top[ii];
				hBottom[ii] = (*radialContourIt)[ii] + barrelColumns[ID]->bottom[ii];
			}
			hTop[3] = 1, hBottom[3] = 1;
			mColumnToGlobal->MultiplyPoint(hTop, hTop);
			mColumnToGlobal->MultiplyPoint(hBottom, hBottom);
			
			for(int ii = 0; ii < 3; ++ii)
			{
				double tmp1 = hTop[ii];
				double tmp2 = hBottom[ii];
				if(tmp1 < colBounds[2*ii])
					colBounds[2*ii] = tmp1;
				if(tmp2 < colBounds[2*ii])
					colBounds[2*ii] = tmp2;
				if(tmp1 > colBounds[2*ii+1])
					colBounds[2*ii+1] = tmp1;
				if(tmp2 > colBounds[2*ii+1])
					colBounds[2*ii+1] = tmp2;
			}
// 			for(int ii = 0; ii < 3; ++ii)
// 			{
// 				double tmp1 = (*radialContourIt)[ii] + barrelColumns[ID]->top[ii];
// 				double tmp2 = (*radialContourIt)[ii] + barrelColumns[ID]->bottom[ii];
// 				if(tmp1 < colBounds[2*ii])
// 					colBounds[2*ii] = tmp1;
// 				if(tmp2 < colBounds[2*ii])
// 					colBounds[2*ii] = tmp2;
// 				if(tmp1 > colBounds[2*ii+1])
// 					colBounds[2*ii+1] = tmp1;
// 				if(tmp2 > colBounds[2*ii+1])
// 					colBounds[2*ii+1] = tmp2;
// 			}
		}
		
		#ifdef DEBUG
		PointsPointerType voxelsInsideColumn = PointsPointerType::New();
		voxelsInsideColumn->SetDataTypeToFloat();
		#endif
		
		ImageDataPointerType volume = createImageVolume(colBounds);
		for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
			for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
				for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
				{
					bool isInColumn = 0;
// 					double pt[3];
// 					pt[0] = x*spacing, pt[1] = y*spacing, pt[2] = z*spacing;
					double somaPt3D[3], hPt[4];
					hPt[0] = x*spacing, hPt[1] = y*spacing, hPt[2] = z*spacing, hPt[3] = 1;
					mGlobalToColumn->MultiplyPoint(hPt, hPt);
					somaPt3D[0] = hPt[0], somaPt3D[1] = hPt[1], somaPt3D[2] = hPt[2];
					
					double t, projectedPt[3];
					double dist = vtkLine::DistanceToLine(somaPt3D, barrelColumns[ID]->top, barrelColumns[ID]->bottom, t, projectedPt);
					dist = sqrt(dist);
					if(dist > maxDistances[ID] || t < 0 || t > 1)
						continue;
					if(dist < minDistances[ID])
						isInColumn = 1;
					// check ambiguous contour cases
					if(!isInColumn)
					{
						PolyDataPointerType polyData = PolyDataPointerType::New();
						polyData->Allocate(1);
						PointsPointerType points = PointsPointerType::New();
						points->SetDataTypeToFloat();
						PolygonPointerType poly = PolygonPointerType::New();
						poly->GetPointIds()->SetNumberOfIds(radialContours[ID].size());
						for(int jj = 0; jj < radialContours[ID].size(); ++jj)
						{
							double tmp[3];
							for(int kk = 0; kk < 3; ++kk)
								tmp[kk] = projectedPt[kk] + radialContours[ID][jj][kk];
							points->InsertNextPoint(tmp);
							poly->GetPointIds()->SetId(jj, jj);
						}
						polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
						polyData->SetPoints(points);
						polyData->Update();
						double closestPoint[3], pCoords[3], dist2;
						int subId;
						double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
						int insidePolygon = polyData->GetCell(0)->EvaluatePosition(somaPt3D, closestPoint, subId, pCoords, dist2, weights);
						if(insidePolygon == 1)
							isInColumn = 1;
						delete [] weights;
					}
					// if still not in column, go to next pt
					if(!isInColumn)
						continue;
					// check whether there is a closer column
					std::map< int, Column * >::const_iterator columnIt;
					for(columnIt = barrelColumns.begin(); columnIt != barrelColumns.end(); ++columnIt)
					{
						int otherID = columnIt->first;
						if(otherID != ID)
						{
							double _t, _projectedPt[3];
							double otherDist = vtkLine::DistanceToLine(somaPt3D, barrelColumns[otherID]->top, barrelColumns[otherID]->bottom, _t, _projectedPt);
							otherDist = sqrt(otherDist);
							if(otherDist < dist)
							{
								PolyDataPointerType polyData = PolyDataPointerType::New();
								polyData->Allocate(1);
								PointsPointerType points = PointsPointerType::New();
								points->SetDataTypeToFloat();
								PolygonPointerType poly = PolygonPointerType::New();
								poly->GetPointIds()->SetNumberOfIds(radialContours[otherID].size());
								for(int jj = 0; jj < radialContours[otherID].size(); ++jj)
								{
									double tmp[3];
									for(int kk = 0; kk < 3; ++kk)
										tmp[kk] = _projectedPt[kk] + radialContours[otherID][jj][kk];
									points->InsertNextPoint(tmp);
									poly->GetPointIds()->SetId(jj, jj);
								}
								polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
								polyData->SetPoints(points);
								polyData->Update();
								double closestPoint[3], pCoords[3], dist2;
								int subId;
								double * weights = new double[polyData->GetCell(0)->GetNumberOfPoints()];
								int insidePolygon = polyData->GetCell(0)->EvaluatePosition(somaPt3D, closestPoint, subId, pCoords, dist2, weights);
								if(insidePolygon == 1)
								{
									delete [] weights;
									isInColumn = 0;
									break;
								}
								delete [] weights;
							}
						}
					}
					// if still inside column, add voxel to profile
					if(isInColumn)
					{
						double depth = sqrt(vtkMath::Distance2BetweenPoints(projectedPt, barrelColumns[ID]->top));
						if(t >= 0 && t <= 1)
						{
							unsigned int bin = (unsigned int)(depth/binSize);	// round down
							zProfile->incrementBin(bin);
						}
						#ifdef DEBUG
						unsigned int oldNrOfPoints = voxelsInsideColumn->GetNumberOfPoints();
						voxelsInsideColumn->InsertPoint(oldNrOfPoints, somaPt3D);
						#endif
					}
				}
		
		#ifdef DEBUG
		std::string colVoxelStr(outFilenameGlobal);
		colVoxelStr += "_voxels_inside_";
		colVoxelStr += SBF->int2Labels[ID];
		Reader * colVoxelWriter = new Reader(colVoxelStr.c_str(), colVoxelStr.c_str());
		colVoxelWriter->writeLandmarkFile(voxelsInsideColumn);
		delete colVoxelWriter;
		#endif
		
		return zProfile;
	}
	else
	{
		std::cout << "Error! Barrel column " << SBF->int2Labels[ID] << " is not present in input!";
		std::cout << "Cannot compute column profile." << std::endl;
		Profile * emptyProfile = new Profile;
		return emptyProfile;
	}
}

// assign points to 50um bins along the column axis
// choose nearest column axis
// works for somata as well as for voxels
Profile* LandmarkAnalyzer::computeSeptumProfile ( PointsPointerType points, std::map< int, Column* > barrelColumns )
{
	double binSize = 50;
	Profile * zProfile = new Profile(binSize);
	unsigned int nrOfPoints = points->GetNumberOfPoints();
	#ifdef DEBUG
	std::flush(std::cout << "Computing soma z profile of " << nrOfPoints << " points..." << std::endl);
	#endif
	for(int ii = 0; ii < nrOfPoints; ++ii)
	{
		int closestID = 0;
		double pt[3], closestPt[3], minDist = 1E9;
		points->GetPoint(ii, pt);
		std::map< int, Column * >::const_iterator barrelColumnIt;
		for(barrelColumnIt = barrelColumns.begin(); barrelColumnIt != barrelColumns.end(); ++barrelColumnIt)
		{
			double projectedPt[3], t, dist;
			dist = vtkLine::DistanceToLine(pt, barrelColumnIt->second->top, barrelColumnIt->second->bottom, t, projectedPt);
			dist = sqrt(dist);
			if(dist < minDist && t >= 0 && t <= 1)
			{
				minDist = dist;
				closestID = barrelColumnIt->first;
				for(int ii = 0; ii < 3; ++ii)
					closestPt[ii] = projectedPt[ii];
			}
		}
		if(closestID)
		{
			double depth = sqrt(vtkMath::Distance2BetweenPoints(closestPt, barrelColumns[closestID]->top));
			unsigned int bin = (unsigned int)(depth/binSize);	// round down
			zProfile->incrementBin(bin);
		}
		#ifdef DEBUG
		else
		{
			std::cout << "Warning! Could not find closest column for point " << ii << ";";
			std::cout << " @ [" << pt[0] << "," << pt[1] << "," << pt[2] << "]" << std::endl;
		}
		#endif
	}
	
	return zProfile;
}

Profile* LandmarkAnalyzer::computeCorrectedSeptumProfile(PointsPointerType points, std::map< int, Column* > barrelColumns, 
							std::map< int, Column* > layerColumns)
{
	double binSize = 50;
	Profile * zProfile = new Profile(binSize);
	unsigned int nrOfPoints = points->GetNumberOfPoints();
	Profile * correctionProfile = createCorrectionProfile();
	#ifdef DEBUG
	std::flush(std::cout << "Computing soma z profile of " << nrOfPoints << " points..." << std::endl);
	#endif
	for(int ii = 0; ii < nrOfPoints; ++ii)
	{
		#ifdef DEBUG
		std::flush(std::cout << "Computing bin for point " << ii+1 << " of " << nrOfPoints << " points..." << std::endl);
		#endif
		int closestID = 0;
		double pt[3], closestPt[3], minDist = 1E9;
		points->GetPoint(ii, pt);
		std::map< int, Column * >::const_iterator layerColumnIt;
		for(layerColumnIt = layerColumns.begin(); layerColumnIt != layerColumns.end(); ++layerColumnIt)
		{
			if(SBF->avgBarrels.find(layerColumnIt->first) != SBF->avgBarrels.end())
			{
				double projectedPt[3], t, dist;
				dist = vtkLine::DistanceToLine(pt, layerColumnIt->second->top, layerColumnIt->second->bottom, t, projectedPt);
				dist = sqrt(dist);
				if(dist < minDist && t >= 0 && t <= 1)
				{
					minDist = dist;
					closestID = layerColumnIt->first;
					for(int jj = 0; jj < 3; ++jj)
						closestPt[jj] = projectedPt[jj];
				}
			}
		}
		if(closestID)
		{
			// compute HSM correction factor scale, then add scale to bin in profile
			double depth = sqrt(vtkMath::Distance2BetweenPoints(closestPt, barrelColumns[closestID]->top));
			unsigned int bin = (unsigned int)(depth/binSize);	// round down
			double binOffset = binSize*0.5;
			double refSupra = 525, refGran = 365, refInfra = 1075;
			double colSupra, colGran, colInfra;
			double supraScale, granScale, infraScale;
			colSupra = SBF->avgTopDist[closestID];
			colGran = SBF->avgBarrels[closestID]->getHeight();
			colInfra = SBF->avgPiaWMDist[closestID] - SBF->avgBarrels[closestID]->getHeight() - SBF->avgTopDist[closestID];
			supraScale = refSupra/colSupra;
			granScale = refGran/colGran;
			infraScale = refInfra/colInfra;
			
			double colDepth = binOffset + bin*binSize;
			double refDepth;
			// scale the supragranular part
			if(colDepth <= colSupra)
			{
				refDepth = colDepth*supraScale;
			}
			// scale the granular part
			else if(colDepth <= colSupra + colGran)
			{
				refDepth = refSupra + (colDepth - colSupra)*granScale;
			}
			// scale the infragranular part
			else
			{
				refDepth = refSupra + refGran + (colDepth - colGran - colSupra)*infraScale;
			}
			// find the scaling bin and interpolate linearly
			// between the two closest bins
			int correctionBin;
			double binDist, binSize;
			double scale;
			binSize = correctionProfile->getBinSize();
			for(int jj = 0; jj < correctionProfile->getProfile()->size(); ++jj)
			{
				// assign here in case column is larger than correctionProfile;
				// then, we automatically scale by 1.0
				correctionBin = jj;
				if((jj+1)*binSize > refDepth)
					break;
			}
			binDist = refDepth - correctionBin*binSize - binOffset;
			if(binDist > 0)
			{
				if(correctionBin < correctionProfile->getProfile()->size() - 1)
				{
					scale = (1 - binDist/binSize)*correctionProfile->getProfile()->at(correctionBin) + binDist/binSize*correctionProfile->getProfile()->at(correctionBin+1);
				}
				else
				{
					scale = correctionProfile->getProfile()->at(correctionBin);
				}
			}
			if(binDist <= 0)
			{
				if(correctionBin > 0)
				{
					scale = (1 + binDist/binSize)*correctionProfile->getProfile()->at(correctionBin) - binDist/binSize*correctionProfile->getProfile()->at(correctionBin-1);
				}
				else
				{
					scale = correctionProfile->getProfile()->at(correctionBin);
				}
			}
			
// #ifdef DEBUG
// std::cout << "column bin\t" << bin << std::endl;
// std::cout << "column depth\t" << colDepth << std::endl;
// std::cout << "reference depth\t" << refDepth << std::endl;
// std::cout << "reference bin\t" << correctionBin << std::endl;
// std::cout << "bin distance\t" << binDist << std::endl;
// std::cout << "final scale\t" << scale << std::endl;
// std::cout << "-----------------------------------------------" << std::endl;
// #endif
			zProfile->addSegment(scale, bin);
		}
		#ifdef DEBUG
		else
		{
			std::cout << "Warning! Could not find closest column for point " << ii << ";";
			std::cout << " @ [" << pt[0] << "," << pt[1] << "," << pt[2] << "]" << std::endl;
		}
		#endif
	}
	
	return zProfile;
}

void LandmarkAnalyzer::correctSomaProfiles(std::map< int, Profile* > columnProfiles)
{
	Profile * refProfile = createCorrectionProfile();
	double binOffset = 25;
	double refSupra = 525, refGran = 365, refInfra = 1075;
	double colSupra, colGran, colInfra;
	double supraScale, granScale, infraScale;
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrelField.find(ID) != barrelField.end() && SBF->avgBarrels.find(ID) != SBF->avgBarrels.end())
		{
#ifdef DEBUG
std::cout << "***********************************************" << std::endl;
std::cout << "computing corrected soma z profile for column " << ID << std::endl;
#endif
			Profile * zProfile = columnProfiles[ID];
			colSupra = SBF->avgTopDist[ID];
			colGran = SBF->avgBarrels[ID]->getHeight();
			colInfra = SBF->avgPiaWMDist[ID] - SBF->avgBarrels[ID]->getHeight() - SBF->avgTopDist[ID];
			supraScale = refSupra/colSupra;
			granScale = refGran/colGran;
			infraScale = refInfra/colInfra;
			for(int ii = 0; ii < zProfile->getProfile()->size(); ++ii)
			{
				double colDepth = binOffset + ii*zProfile->getBinSize();
				double refDepth;
				// scale the supragranular part
				if(colDepth <= colSupra)
				{
					refDepth = colDepth*supraScale;
				}
				// scale the granular part
				else if(colDepth <= colSupra + colGran)
				{
					refDepth = refSupra + (colDepth - colSupra)*granScale;
				}
				// scale the infragranular part
				else
				{
					refDepth = refSupra + refGran + (colDepth - colGran - colSupra)*infraScale;
				}
				// find the scaling bin and interpolate linearly
				// between the two closest bins
				int bin;
				double binDist, binSize;
				double scale;
				binSize = refProfile->getBinSize();
				for(int jj = 0; jj < refProfile->getProfile()->size(); ++jj)
				{
					// assign here in case column is larger than refProfile;
					// then, we automatically scale by 1.0
					bin = jj;
					if((jj+1)*binSize > refDepth)
						break;
				}
				binDist = refDepth - bin*binSize - binOffset;
				if(binDist > 0)
				{
					if(bin < refProfile->getProfile()->size() - 1)
					{
						scale = (1 - binDist/binSize)*refProfile->getProfile()->at(bin) + binDist/binSize*refProfile->getProfile()->at(bin+1);
					}
					else
					{
						scale = refProfile->getProfile()->at(bin);
					}
				}
				if(binDist <= 0)
				{
					if(bin > 0)
					{
						scale = (1 + binDist/binSize)*refProfile->getProfile()->at(bin) - binDist/binSize*refProfile->getProfile()->at(bin-1);
					}
					else
					{
						scale = refProfile->getProfile()->at(bin);
					}
				}
				double val = zProfile->getProfile()->at(ii);
#ifdef DEBUG
std::cout << "column bin\t" << ii << std::endl;
std::cout << "column depth:\t" << colDepth << std::endl;
std::cout << "reference depth:\t" << refDepth << std::endl;
std::cout << "reference bin\t" << bin << std::endl;
std::cout << "bin distance\t" << binDist << std::endl;
std::cout << "old value\t" << val << std::endl;
#endif
				val *= scale;
				zProfile->getProfile()->at(ii) = val;
#ifdef DEBUG
std::cout << "final scale\t" << scale << std::endl;
std::cout << "new value\t" << val << std::endl;
std::cout << "-----------------------------------------------" << std::endl;
#endif
			}
			zProfile->updateIntegral();
		}
	}
	delete refProfile;
}

int LandmarkAnalyzer::computeLaminarPosition(double x[3], int columnID)
{
	if(!supraFlag || !granFlag || !infraFlag)
	{
		std::cout << "Warning! Trying to compute laminar position without barrel field layer data!" << std::endl;
		return -1;
	}
	
	if(supragranularLayer.find(columnID) == supragranularLayer.end()
		|| granularLayer.find(columnID) == granularLayer.end()
		|| infragranularLayer.find(columnID) == infragranularLayer.end()
	)
	{
		std::cout << "Warning! Barrel column " << SBF->int2Labels[columnID] << " not present in barrel field layer data!" << std::endl;
		return -1;
	}
	
	// by trial-and-error testing which contours of supra/infra are "top" in the Amira files...
	double L1ratio = 0.727; // 1-0.273
	double L5ratio = 0.476; // 1-0.524
	int layer = 0;
	
	double dist1 = 0, t1, closestPt1[3];
	dist1 = vtkLine::DistanceToLine(x, supragranularLayer[columnID]->top, supragranularLayer[columnID]->bottom, t1, closestPt1);
	if(t1 >= 0 && t1 <= 1)
	{
		if(t1 >= L1ratio)
			layer =  L1;
		else
			layer = L23;
		return layer;
	}
	
	double dist2 = 0, t2, closestPt2[3];
	dist2 = vtkLine::DistanceToLine(x, granularLayer[columnID]->top, granularLayer[columnID]->bottom, t2, closestPt2);
	if(t2 >= 0 && t2 <= 1)
	{
		layer = L4;
		return layer;
	}
	
	double dist3 = 0, t3, closestPt3[3];
	dist3 = vtkLine::DistanceToLine(x, infragranularLayer[columnID]->top, infragranularLayer[columnID]->bottom, t3, closestPt3);
	if(t3 >= 0 && t3 <= 1)
	{
		if(t3 <= L5ratio)
			layer =  L5;
		else
			layer = L6;
		return layer;
	}
	
	return layer;
}

void LandmarkAnalyzer::computeDensityProfile(Profile* rawProfile, Profile* volumeProfile, Profile* densityProfile)
{
	for(int ii = 0; ii < rawProfile->getProfile()->size() && ii < volumeProfile->getProfile()->size(); ++ii)
	{
		double length = rawProfile->getProfile()->at(ii);
		double vol = volumeProfile->getProfile()->at(ii);
		if(vol)
		{
			double density = length/vol;
			densityProfile->addSegment(density, ii);
		}
	}
}

void LandmarkAnalyzer::createCountDataStructures ( std::map< int, PointsPointerType >& somaColumns,
						   std::map< int, Column* >& barrelColumns,
						   std::map< int, std::vector< std::vector< double > > >& radialContours,
						   std::map< int, double >& minDistances,
						   std::map< int, double >& maxDistances )
{
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrelField.find(ID) != barrelField.end())
		{
			PointsPointerType somaColumn = PointsPointerType::New();
			somaColumn->Allocate(1);
			somaColumn->SetDataTypeToFloat();
			somaColumns.insert(std::pair< int, PointsPointerType >(ID, somaColumn));
			
			double columnTop[3], columnBottom[3];
			
			// create barrel columns by extrapolation of reconstructed barrel contours
			if(piaFlag && wmFlag)
			{
				double columnAxis[3], columnCenter[3];
				for(int ii = 0; ii < 3; ++ii)
				{
					columnAxis[ii] = barrelField[ID]->top[ii] - barrelField[ID]->bottom[ii];
					columnCenter[ii] = 0.5*(barrelField[ID]->top[ii] + barrelField[ID]->bottom[ii]);
				}
				vtkMath::Normalize(columnAxis);
				
				pia->intersectLine(columnAxis, columnCenter);
				if(!pia->isValid())
				{
					std::cout << "Error! Could not compute intersection of column axis with Pia!" << std::endl;
					continue;
				}
				WM->intersectLine(columnAxis, columnCenter);
				if(!WM->isValid())
				{
					std::cout << "Error! Could not compute intersection of column axis with WM!" << std::endl;
					continue;
				}
				pia->getLastIntersectPoint(columnTop);
				WM->getLastIntersectPoint(columnBottom);
			}
			// otherwise use supplied contours (e.g. during further analysis supra/gran/infra)
			else
			{
				for(int ii = 0; ii < 3; ++ii)
				{
					columnTop[ii] = barrelField[ID]->top[ii];
					columnBottom[ii] = barrelField[ID]->bottom[ii];
				}
			}
			
			Column * thisCol = new Column(barrelField[ID]->contours, columnTop, columnBottom);
			barrelColumns.insert(std::pair< int, Column * >(ID, thisCol));
			
			std::vector< std::vector< double > > contourPts;
			double minDist = 1E09, maxDist = 0;
			for(int ii = 0; ii < thisCol->contours->GetCell(0)->GetNumberOfPoints(); ++ii)
			{
				double pt[3], dist = 0, t, closestPt[3];
				std::vector< double > radialPt;
				thisCol->contours->GetCell(0)->GetPoints()->GetPoint(ii, pt);
				dist = vtkLine::DistanceToLine(pt, thisCol->top, thisCol->bottom, t, closestPt);
				dist = sqrt(dist);
				if(dist < minDist)
					minDist = dist;
				if(dist > maxDist)
					maxDist = dist;
				for(int jj = 0; jj < 3; ++jj)
					radialPt.push_back(pt[jj] - closestPt[jj]);
				contourPts.push_back(radialPt);
			}
			radialContours.insert(std::pair< int, std::vector< std::vector< double > > >(ID, contourPts));
			minDistances.insert(std::pair< int, double >(ID, minDist));
			maxDistances.insert(std::pair< int, double >(ID, maxDist));
// 			#ifdef DEBUG
// 			double boundingBox[6], columnAxis2[3];
// 			for(int ii = 0; ii < 3; ++ii)
// 				columnAxis2[ii] = barrelColumns[ID]->top[ii] - barrelColumns[ID]->bottom[ii];
// 			vtkMath::Normalize(columnAxis2);
// 			barrelField[ID]->contours->GetBounds(boundingBox);
// 			std::cout << "*************************************" << std::endl;
// 			std::cout << "Barrel " << SBF->int2Labels[ID] << ":" << std::endl;
// 			std::cout << "Top @ [" << barrelField[ID]->top[0] << "," << barrelField[ID]->top[1] << "," << barrelField[ID]->top[2] << "]" << std::endl;
// 			std::cout << "Bottom @ [" << barrelField[ID]->bottom[0] << "," << barrelField[ID]->bottom[1] << "," << barrelField[ID]->bottom[2] << "]" << std::endl;
// 			std::cout << "Center @ [" << columnCenter[0] << "," << columnCenter[1] << "," << columnCenter[2] << "]" << std::endl;
// 			std::cout << "Axis = [" << columnAxis[0] << "," << columnAxis[1] << "," << columnAxis[2] << "]" << std::endl;
// 			std::cout << "Bounding box:" << std::endl;
// 			std::cout << "[" << boundingBox[0] << "," << boundingBox[1] << "]" << std::endl;
// 			std::cout << "[" << boundingBox[2] << "," << boundingBox[3] << "]" << std::endl;
// 			std::cout << "[" << boundingBox[4] << "," << boundingBox[5] << "]" << std::endl;
// 			std::cout << "Column " << SBF->int2Labels[ID] << ":" << std::endl;
// 			std::cout << "Top @ [" << barrelColumns[ID]->top[0] << "," << barrelColumns[ID]->top[1] << "," << barrelColumns[ID]->top[2] << "]" << std::endl;
// 			std::cout << "Bottom @ [" << barrelColumns[ID]->bottom[0] << "," << barrelColumns[ID]->bottom[1] << "," << barrelColumns[ID]->bottom[2] << "]" << std::endl;
// 			std::cout << "Column axis = [" << columnAxis2[0] << "," << columnAxis2[1] << "," << columnAxis2[2] << "]" << std::endl;
// 			std::cout << "Min radius = " << minDist << "um" << std::endl;
// 			std::cout << "Max radius = " << maxDist << "um" << std::endl;
// 			#endif
		}
	}
}

void LandmarkAnalyzer::createLayerCountDataStructures ( std::map< int, PointsPointerType >& somaColumns,
						   std::map< int, Column* >& barrelColumns,
						   std::map< int, std::vector< std::vector< double > > >& radialContours,
						   std::map< int, double >& minDistances,
						   std::map< int, double >& maxDistances )
{
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(barrelFieldLayer.find(ID) != barrelFieldLayer.end())
		{
			PointsPointerType somaColumn = PointsPointerType::New();
			somaColumn->Allocate(1);
			somaColumn->SetDataTypeToFloat();
			somaColumns.insert(std::pair< int, PointsPointerType >(ID, somaColumn));
			
			double columnTop[3], columnBottom[3];
			for(int ii = 0; ii < 3; ++ii)
			{
				columnTop[ii] = barrelFieldLayer[ID]->top[ii];
				columnBottom[ii] = barrelFieldLayer[ID]->bottom[ii];
			}
			
			Column * thisCol = new Column(barrelFieldLayer[ID]->contours, columnTop, columnBottom);
			barrelColumns.insert(std::pair< int, Column * >(ID, thisCol));
			
			std::vector< std::vector< double > > contourPts;
			double minDist = 1E09, maxDist = 0;
			for(int ii = 0; ii < thisCol->contours->GetCell(0)->GetNumberOfPoints(); ++ii)
			{
				double pt[3], dist = 0, t, closestPt[3];
				std::vector< double > radialPt;
				thisCol->contours->GetCell(0)->GetPoints()->GetPoint(ii, pt);
				dist = vtkLine::DistanceToLine(pt, thisCol->top, thisCol->bottom, t, closestPt);
				dist = sqrt(dist);
				if(dist < minDist)
					minDist = dist;
				if(dist > maxDist)
					maxDist = dist;
				for(int jj = 0; jj < 3; ++jj)
					radialPt.push_back(pt[jj] - closestPt[jj]);
				contourPts.push_back(radialPt);
			}
			radialContours.insert(std::pair< int, std::vector< std::vector< double > > >(ID, contourPts));
			minDistances.insert(std::pair< int, double >(ID, minDist));
			maxDistances.insert(std::pair< int, double >(ID, maxDist));
// 			#ifdef DEBUG
// 			double boundingBox[6], columnAxis2[3];
// 			for(int ii = 0; ii < 3; ++ii)
// 				columnAxis2[ii] = barrelColumns[ID]->top[ii] - barrelColumns[ID]->bottom[ii];
// 			vtkMath::Normalize(columnAxis2);
// 			barrelFieldLayer[ID]->contours->GetBounds(boundingBox);
// 			std::cout << "*************************************" << std::endl;
// 			std::cout << "Barrel " << SBF->int2Labels[ID] << ":" << std::endl;
// 			std::cout << "Top @ [" << barrelFieldLayer[ID]->top[0] << "," << barrelFieldLayer[ID]->top[1] << "," << barrelFieldLayer[ID]->top[2] << "]" << std::endl;
// 			std::cout << "Bottom @ [" << barrelFieldLayer[ID]->bottom[0] << "," << barrelFieldLayer[ID]->bottom[1] << "," << barrelFieldLayer[ID]->bottom[2] << "]" << std::endl;
// 			std::cout << "Center @ [" << columnCenter[0] << "," << columnCenter[1] << "," << columnCenter[2] << "]" << std::endl;
// 			std::cout << "Axis = [" << columnAxis[0] << "," << columnAxis[1] << "," << columnAxis[2] << "]" << std::endl;
// 			std::cout << "Bounding box:" << std::endl;
// 			std::cout << "[" << boundingBox[0] << "," << boundingBox[1] << "]" << std::endl;
// 			std::cout << "[" << boundingBox[2] << "," << boundingBox[3] << "]" << std::endl;
// 			std::cout << "[" << boundingBox[4] << "," << boundingBox[5] << "]" << std::endl;
// 			std::cout << "Column " << SBF->int2Labels[ID] << ":" << std::endl;
// 			std::cout << "Top @ [" << barrelColumns[ID]->top[0] << "," << barrelColumns[ID]->top[1] << "," << barrelColumns[ID]->top[2] << "]" << std::endl;
// 			std::cout << "Bottom @ [" << barrelColumns[ID]->bottom[0] << "," << barrelColumns[ID]->bottom[1] << "," << barrelColumns[ID]->bottom[2] << "]" << std::endl;
// 			std::cout << "Column axis = [" << columnAxis2[0] << "," << columnAxis2[1] << "," << columnAxis2[2] << "]" << std::endl;
// 			std::cout << "Min radius = " << minDist << "um" << std::endl;
// 			std::cout << "Max radius = " << maxDist << "um" << std::endl;
// 			#endif
		}
	}
}

void LandmarkAnalyzer::writeColumnSomaProfiles ( const char* outputFilename )
{
	std::list< int >::const_iterator labelIt;
	for(labelIt = SBF->barrelLabels.begin(); labelIt != SBF->barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(columnSomaProfiles.find(ID) != columnSomaProfiles.end())
		{
			std::string ofName(outputFilename);
			ofName += "_zprofile_";
			ofName += SBF->int2Labels[ID];
			ofName += ".csv";
			double binSize = columnSomaProfiles[ID]->getBinSize();
			double offset = 0.5*binSize;
			std::vector< double > * zProfile = columnSomaProfiles[ID]->getProfile();
			std::ofstream ProfileWriter;
			ProfileWriter.open(ofName.c_str());
			ProfileWriter << "# 1D Profile of " << SBF->int2Labels[ID] << std::endl;
			ProfileWriter << "Depth [um]\tcounts per bin" << std::endl;
			for(int ii = 0; ii < zProfile->size(); ++ii)
				ProfileWriter << ii*binSize + offset << "\t" << (*zProfile)[ii] << std::endl;
			ProfileWriter.close();
		}
	}
}

void LandmarkAnalyzer::writeZProfile ( Profile* zProfile, const char* outputFilename, double globalOffset )
{
	std::string ofName(outputFilename);
	ofName += "_zprofile.csv";
	double binSize = zProfile->getBinSize();
	double offset = 0.5*binSize + globalOffset;
	std::vector< double > * profile = zProfile->getProfile();
	std::ofstream ProfileWriter;
	ProfileWriter.open(ofName.c_str());
	ProfileWriter << "# 1D Profile" << std::endl;
	ProfileWriter << "Depth [um]\tcounts per bin" << std::endl;
	for(int ii = 0; ii < profile->size(); ++ii)
		ProfileWriter << ii*binSize + offset << "\t" << (*profile)[ii] << std::endl;
	ProfileWriter.close();
}

Profile* LandmarkAnalyzer::createCorrectionProfile()
{
	Profile * refProfile = new Profile(50);
	// WR09, WR26, WR27, WR611
	refProfile->getProfile()->push_back(0.88);
	refProfile->getProfile()->push_back(1.22);
	refProfile->getProfile()->push_back(0.68);
	refProfile->getProfile()->push_back(0.87);
	refProfile->getProfile()->push_back(1.21);
	refProfile->getProfile()->push_back(1.17);
	refProfile->getProfile()->push_back(1.63);
	refProfile->getProfile()->push_back(1.66);
	refProfile->getProfile()->push_back(1.06);
	refProfile->getProfile()->push_back(0.96);
	refProfile->getProfile()->push_back(0.66);
	refProfile->getProfile()->push_back(0.94);
	refProfile->getProfile()->push_back(1.18);
	refProfile->getProfile()->push_back(0.89);
	refProfile->getProfile()->push_back(1.63);
	refProfile->getProfile()->push_back(2.35);
	refProfile->getProfile()->push_back(1.99);
	refProfile->getProfile()->push_back(1.97);
	refProfile->getProfile()->push_back(1.87);
	refProfile->getProfile()->push_back(2.01);
	refProfile->getProfile()->push_back(1.93);
	refProfile->getProfile()->push_back(2.12);
	refProfile->getProfile()->push_back(1.35);
	refProfile->getProfile()->push_back(1.18);
	refProfile->getProfile()->push_back(1.55);
	refProfile->getProfile()->push_back(1.59);
	refProfile->getProfile()->push_back(1.66);
	refProfile->getProfile()->push_back(1.99);
	refProfile->getProfile()->push_back(1.29);
	refProfile->getProfile()->push_back(1.59);
	refProfile->getProfile()->push_back(1.56);
	refProfile->getProfile()->push_back(1.35);
	refProfile->getProfile()->push_back(1.29);
	refProfile->getProfile()->push_back(1.47);
	refProfile->getProfile()->push_back(1.37);
	refProfile->getProfile()->push_back(1.25);
	refProfile->getProfile()->push_back(1.60);
	refProfile->getProfile()->push_back(2.04);
	refProfile->getProfile()->push_back(1.17);
	refProfile->getProfile()->push_back(1.01);
	refProfile->getProfile()->push_back(2.85);
	refProfile->getProfile()->push_back(1.0);
	refProfile->getProfile()->push_back(1.0);
	refProfile->getProfile()->push_back(1.0);
	refProfile->getProfile()->push_back(1.0);
	refProfile->getProfile()->push_back(1.0);
	refProfile->getProfile()->push_back(1.0);
	refProfile->getProfile()->push_back(1.0);
	refProfile->getProfile()->push_back(1.0);
	refProfile->getProfile()->push_back(1.0);
	refProfile->getProfile()->push_back(1.0);
	refProfile->getProfile()->push_back(1.0);
	
	// WR09, WR26, WR27
// 	refProfile->getProfile()->push_back(0.83);
// 	refProfile->getProfile()->push_back(1.09);
// 	refProfile->getProfile()->push_back(0.63);
// 	refProfile->getProfile()->push_back(0.86);
// 	refProfile->getProfile()->push_back(1.30);
// 	refProfile->getProfile()->push_back(1.21);
// 	refProfile->getProfile()->push_back(1.60);
// 	refProfile->getProfile()->push_back(1.54);
// 	refProfile->getProfile()->push_back(0.96);
// 	refProfile->getProfile()->push_back(0.87);
// 	refProfile->getProfile()->push_back(0.58);
// 	refProfile->getProfile()->push_back(0.86);
// 	refProfile->getProfile()->push_back(1.17);
// 	refProfile->getProfile()->push_back(0.82);
// 	refProfile->getProfile()->push_back(1.46);
// 	refProfile->getProfile()->push_back(2.10);
// 	refProfile->getProfile()->push_back(1.72);
// 	refProfile->getProfile()->push_back(1.67);
// 	refProfile->getProfile()->push_back(1.60);
// 	refProfile->getProfile()->push_back(1.83);
// 	refProfile->getProfile()->push_back(1.69);
// 	refProfile->getProfile()->push_back(1.85);
// 	refProfile->getProfile()->push_back(1.19);
// 	refProfile->getProfile()->push_back(1.02);
// 	refProfile->getProfile()->push_back(1.41);
// 	refProfile->getProfile()->push_back(1.50);
// 	refProfile->getProfile()->push_back(1.41);
// 	refProfile->getProfile()->push_back(1.77);
// 	refProfile->getProfile()->push_back(1.16);
// 	refProfile->getProfile()->push_back(1.51);
// 	refProfile->getProfile()->push_back(1.47);
// 	refProfile->getProfile()->push_back(1.22);
// 	refProfile->getProfile()->push_back(1.18);
// 	refProfile->getProfile()->push_back(1.27);
// 	refProfile->getProfile()->push_back(1.22);
// 	refProfile->getProfile()->push_back(1.22);
// 	refProfile->getProfile()->push_back(1.55);
// 	refProfile->getProfile()->push_back(2.26);
// 	refProfile->getProfile()->push_back(1.34);
// 	refProfile->getProfile()->push_back(1.13);
// 	refProfile->getProfile()->push_back(2.85);
// 	refProfile->getProfile()->push_back(1.0);
// 	refProfile->getProfile()->push_back(1.0);
// 	refProfile->getProfile()->push_back(1.0);
// 	refProfile->getProfile()->push_back(1.0);
// 	refProfile->getProfile()->push_back(1.0);
// 	refProfile->getProfile()->push_back(1.0);
// 	refProfile->getProfile()->push_back(1.0);
// 	refProfile->getProfile()->push_back(1.0);
// 	refProfile->getProfile()->push_back(1.0);
// 	refProfile->getProfile()->push_back(1.0);
// 	refProfile->getProfile()->push_back(1.0);
	
	return refProfile;
}

ImageDataPointerType LandmarkAnalyzer::createImageVolume ( double bounds[6] )
{
	double spacing[3];
	int dims[6];
	spacing[0] = spacing[1] = spacing[2] = SPACING;
	calculateExtent(bounds, dims);
	
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	volume->SetExtent(dims);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToUnsignedChar();
// 	volume->Print(std::cout);
	volume->AllocateScalars();
// 	std::flush(std::cout << "max extent of input: [" << xMin << "," << xMax << "], [" << yMin << "," << yMax << "],[" << zMin << "," << zMax << "]" << std::endl);
// 	std::flush(std::cout << "Allocating memory for image  with dimensions [" << dims[0] << "," << dims[1] << "], [" << dims[2] << "," << dims[3] << "],[" << dims[4] << "," << dims[5] << "]" << std::endl);

	for(int z = dims[4]; z <= dims[5]; ++z)
	       for(int y = dims[2]; y <= dims[3]; ++y)
			for(int x = dims[0]; x <= dims[1]; ++x)
			{
				unsigned char * px = static_cast< unsigned char * >(volume->GetScalarPointer(x, y, z));
				*px = 0;
			}

	volume->Update();
	return volume;
}

void LandmarkAnalyzer::calculateExtent ( double bounds[6], int extent[6] )
{
	double spacing[3];
	spacing[0] = spacing[1] = spacing[2] = SPACING;
	
	//make sure that maxCoordinates are inside an integer number of cells defined by spacing
	extent[0] = (bounds[0] - spacing[0])/spacing[0]/* - 0.5*/;
	extent[1] = (bounds[1] + spacing[0])/spacing[0]/* + 0.5*/;
	extent[2] = (bounds[2] - spacing[1])/spacing[1]/* - 0.5*/;
	extent[3] = (bounds[3] + spacing[1])/spacing[1]/* + 0.5*/;
	extent[4] = (bounds[4] - spacing[2])/spacing[2]/* - 0.5*/;
	extent[5] = (bounds[5] + spacing[2])/spacing[2]/* + 0.5*/;
}

