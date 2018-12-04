#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/amiraReader.h"
#include "../../common/barrel_field.h"
#include "../../Interneuron/src/helper.h"
#include <set>
#include <utility>

/* Initialization Functions */
void initalizeSomaPositions();
void initializeConstants();
void initializeLists();

/* Structs */
/* Iterative Updating of Sample Size (N), Mean (M) and SD (Q) [SD = sqrt(Q / (N-1))]
 * Calculation for Connection Probability:
 * - CellWeight = #PresynapticCells_in_ColumnCellTypePairSelection
 * - CellNorm = #PresynapticCells_total / (#PresynapticColumns * #PresynapticCellTypes)
 * - C = Convergences from PresynapticColumnCellTypePair onto Postsynaptic Cell
 * - Convergences = C * CellWeigth/CellNorm (e.g., C * 0.1 * #PresynapticColumns * #PresynapticCellTypes)
 * - N = #PostsynapticCells * #PresynapticColumns * #PresynapticCellTypes
 * -> #PresynapticColumns * #PresynapticCellTypes necessary because go through each PresynapticColumnCellTypePair;
 * cancels each other out
 * => Normalized by #PostsynapticCells only */
struct bin
{
	int N;
	float M;
	float Q;
};
// Soma Position
struct somaPosition
{
	somaPosition() : x(0), y(0), z(0) {	}
	somaPosition(double _x, double _y, double _z) : x(_x), y(_y), z(_z) { }
	double x;
	double y;
	double z;
};

/* Variables */
std::map< unsigned int, const char * > int2CelltypeLabels;
std::map< unsigned int, const char * > int2ColumnLabels;
std::map< std::string, unsigned int> celltypeLabels2Int;
std::map< std::string, unsigned int> ColumnLabels2Int;
std::vector <int> PostCelltypeList;
std::vector <int> PreCelltypeList;
std::vector <int> PostColumnList;
std::vector <int> PreColumnList;
std::map< ColumnCellTypePair, somaPosition> celltypeLabels2SomaPos;

/* Functions */
void computeSynapseMeanTable();
void writeConvergenceTableCSV();
void writeSomaPosition();
std::vector< float > computeSynMean(std::vector< CellTableRow * > selection);
void writeSynapseMeanTable(const char * Outputfname, std::vector< float > synapsesMean, std::map< ColumnCellTypePair, unsigned int > header, int sz);
void outputCorrBins(const char * inputfname);
void outputCorrBinsItself(const char * inputfname, bool synapse);
void writeProbabilityMeanTable(std::string compartment, bool LocalSubType);
void writeProbabilityMeanTable(const char * inputfname, const char * outputfname, const char * numCellsfname, bool LocalSubType, int AV);
/* Check This with 1D Profile Plot */
void writeCPMap(std::string spatial);
void writeCPMapEXIN(std::string spatial);
void extractFromTable(const char * outputfname, std::string compartment, std::list< unsigned int > preColumns, std::list< unsigned int > preCelltypes, std::list< unsigned int > postColumns, std::list< unsigned int > postCelltypes);

/* Use these */
void writeCPMap3x3(bool bol3x3, std::string compartment);
void writeCPMap3x3EXIN(bool bol3x3, std::string compartment);
void writeCPMap3x3Single(ColumnCellTypePair PreSynapse, ColumnCellTypePair PostSynapse);
void writeCPMap_C2Output(std::string compartment, std::string spatial);
void writeCPMap_C2OutputEXIN(std::string compartment, std::string spatial);
/* ** */
void writeCPMap3(std::string spatial); // Probably Deprecated!
void writeCP1D(std::string spatial, std::string compartment);
void writeCP1D_Cut(std::string spatial, std::string compartment);
void columnsCutDisplay(std::string spatial, const char * fname);
void cutImgVol(const char * Inputfname, const char * Outputfname, std::string spatial);
void getCorrPerBin(std::vector< bin * > valVector, std::vector< CellTableRow * > row1, std::vector< CellTableRow * > row2, int binsz);
void D2Test();
std::map< ColumnCellTypePair, unsigned int > getCellNumbers(const char* fname);
std::map< ColumnCellTypePair, unsigned int > getCellNumbers();
void convergenceDepth1D();
void convergenceDepth2D();

/* Save/Store Functions */
void storeBinMapAsCSV(const char* fname, std::map< MatrixIndexType, bin *> binMap, double binsz);
void storeBinAsCSV(const char* fname, std::map< int, bin *> bin, double binsz);

/* Helper */
void getConvergence(CellTable * convergence, const char * inputfname, const char * numCellsfname, bool LocalSubTypes);
void getConvergence(CellTable * convergence);
void getConvergence(CellTable * convergence, std::string compartment);
void getConvergence(CellTable * convergence, std::string compartment, bool LocalSubTypes);
void getConnectome(CellTable * connectome);
int getTotalCellNumbers(std::map< ColumnCellTypePair, unsigned int > numCells, std::list< unsigned int > columnsPre, std::list< unsigned int > cellTypesPre);
void getPlane(double n[3], double x[3], std::string spatial);
void connectCellIDMorphologyFile();


int main(int argc, const char** argv)
{
//	writeCPMap3x3EXIN(false,"Apical");
//	writeCPMap3x3EXIN(true,"Apical");

//	writeCPMap3x3(false,"Basal");
//	writeCPMap3x3(true,"Basal");

//	writeCP1D_Cut("arcish", "Total");
//	writeCP1D_Cut("arcish", "Apical");
//	writeCP1D_Cut("arcish", "Basal");
//	writeCPMap_C2OutputEXIN("Total","");
//	writeCPMap_C2OutputEXIN("Total","");

//	writeCPMap_C2Output("Apical");
//	writeCPMap_C2Output("Basal");
//	writeCPMapEXIN("rowish");
//	writeCPMapEXIN("arcish");

//	outputCorrBinsItself(argv[1],false);

//	writeProbabilityMeanTable("Total",true);
//	writeProbabilityMeanTable("Basal",true);
//	writeProbabilityMeanTable("Apical",true);
//
//	writeProbabilityMeanTable("Total",false);
//	writeProbabilityMeanTable("Basal",false);
//	writeProbabilityMeanTable("Apical",false);

//	int AV = 0;
//	const char * numCellsfname = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/nrCells.csv";
//	std::string inputfname = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/convergenceTotal.csv";
//	std::string Outputfname = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/ConnectionProbabilityTotal_";
//
//	bool LocalSubType = true;
//	Outputfname = Outputfname + "LocalSubTypes_";
//	Outputfname = Outputfname + "SD.csv";
//
//	writeProbabilityMeanTable(inputfname.c_str(), Outputfname.c_str(), numCellsfname, LocalSubType, AV);


//	const char * inputfname = "/nas1/Data_regger/NeuroNet/NN_paper/PSD_density/Vincent_2015_V1/Network3x3/cache_uniform/convergenceC2IN.csv";
//	const char * outputfname = "/nas1/Data_regger/NeuroNet/NN_paper/PSD_density/Vincent_2015_V1/Network3x3/cache_uniform/ConnectionProbabilityC2IN_cpp.csv";
//	bool LocalSubType = true;
//	int AV = 1;
//	const char * numCellsfname = "";
//
//	// DOES NOT WORK?!
//	writeProbabilityMeanTable(inputfname, outputfname, numCellsfname, LocalSubType, AV);

	/* Extract L4->L6A */
//	std::list< unsigned int > postCelltypes;
//	postCelltypes.push_back(L6ccinv);
//	postCelltypes.push_back(L6cc);
//	postCelltypes.push_back(L6ct);
//	std::list< unsigned int > postColumns;
//	postColumns.push_back(C2);
//
//	std::list< unsigned int > preCelltypes;
//	preCelltypes.push_back(L4ssaxon);
//	preCelltypes.push_back(L4spaxon);
//	preCelltypes.push_back(L4pyaxon);
//
//	std::list< unsigned int > preColumns;
//	preColumns.push_back(C2);
//
//	std::string compartment = "Apical";
//	const char * outputfname = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/L4-L6_Apical.csv";
//	extractFromTable(outputfname, compartment, preColumns, preCelltypes, postColumns, postCelltypes);
//
//	compartment = "Basal";
//	outputfname = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/L4-L6_Basal.csv";
//	extractFromTable(outputfname, compartment, preColumns, preCelltypes, postColumns, postCelltypes);


	/* Extract C2 connectivity All->C2 */
	std::list< unsigned int > postCelltypes;
	std::list< unsigned int > postColumns;
	postColumns.push_back(C2);

	std::list< unsigned int > preCelltypes;
	std::list< unsigned int > preColumns;
	preColumns.push_back(C2);

	std::string compartment = "Total";
	const char * outputfname = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/convergenceTotal_C2-C2.csv";
	extractFromTable(outputfname, compartment, preColumns, preCelltypes, postColumns, postCelltypes);

	compartment = "Basal";
	outputfname =  "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/convergenceBasal_C2-C2.csv";
	extractFromTable(outputfname, compartment, preColumns, preCelltypes, postColumns, postCelltypes);

//	connectCellIDMorphologyFile();
}

/* Extracts and writes excerpt from Convergences (or Synapse) Table specified by
 * presynaptic Columns and Celltypes, if empty select all Columns/CellTypes
 * postsynaptic Columns and Celltypes, if empty select all Columns/CellTypes
 * compartment "Apical" "Basal" "Total"
 * outputfname path/to/.csv-file
 */
void extractFromTable(const char * outputfname, std::string compartment, std::list< unsigned int > preColumns, std::list< unsigned int > preCelltypes, std::list< unsigned int > postColumns, std::list< unsigned int > postCelltypes)
{
	std::ofstream CSVwriter;
	CSVwriter.open(outputfname);

	if(!CSVwriter.fail())
	{
		CellTable * convergence = new CellTable;
		getConvergence(convergence, compartment);

		std::map< ColumnCellTypePair, unsigned int > header = convergence->header;
		std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(postColumns, postCelltypes);

		CSVwriter << "CELLID,CELLTYPE,COLUMN,SOMA_X,SOMA_Y,SOMA_Z,INSIDE_COLUMN," ;

		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
		{
			bool bolColumn = true;
			if(preColumns.size()>0)
			{
				bolColumn = (std::find(preColumns.begin(), preColumns.end(), headerit->first.first) != preColumns.end());
			}

			bool bolCellType = true;
			if(preCelltypes.size()>0)
			{
				bolCellType = (std::find(preCelltypes.begin(), preCelltypes.end(), headerit->first.second) != preCelltypes.end());
			}

			// Presynaptic Cells
			if (bolColumn && bolCellType)
			{
				CSVwriter << int2ColumnLabels[headerit->first.first] << "_" << int2CelltypeLabels[headerit->first.second] << ",";
			}
		}
		CSVwriter << std::endl;

		for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
		{

			std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

			CSVwriter << (*rowit)->cellID << "," << int2CelltypeLabels[(*rowit)->cellType] << "," << int2ColumnLabels[(*rowit)->column] << ",";
			CSVwriter << (*rowit)->somaLocation[0] << "," << (*rowit)->somaLocation[1] << "," << (*rowit)->somaLocation[2] << ",";
			CSVwriter << (*rowit)->insideColumn << ",";

			for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
			{
				bool bolColumn = true;
				if(preColumns.size()>0)
				{
					bolColumn = (std::find(preColumns.begin(), preColumns.end(), headerit->first.first) != preColumns.end());
				}

				bool bolCellType = true;
				if(preCelltypes.size()>0)
				{
					bolCellType = (std::find(preCelltypes.begin(), preCelltypes.end(), headerit->first.second) != preCelltypes.end());
				}

				// Presynaptic Cells
				if (bolColumn && bolCellType)
				{
					std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
					CSVwriter << (*Convit) << ",";
				}
			}
			CSVwriter << std::endl;
		}
		std::cout << "Wrote " << outputfname << " successfully! " << std::endl;
	}
	else
	{
		std::cout << "ERROR! Writing " << outputfname << " failed! " << std::endl;
	}
	CSVwriter.close();
}

void connectCellIDMorphologyFile()
{
	// Connect Cell IDs in SynapseMatrix with Cell IDs in Morphology file (and filenames)
	const char * fname = "/nas1/Data_daniel/Network/INColumn_v6/MorphologiesFilenameID.csv";
	std::ofstream CSVwriter;
	CSVwriter.open(fname);

	if(!CSVwriter.fail())
	{
		// Read SynapseMatrix
		const char * inputSynapse = "/nas1/Data_daniel/Network/INColumn_v6/cache/synapses3x3All.am";
		Reader * cellTableReader = new Reader(inputSynapse, inputSynapse);
		CellTable * connectome  = new CellTable;
		cellTableReader->readSynapsesPerCellTable(connectome);

		// Extract CellType/Column Combination (postsynaptic)
		std::map< ColumnCellTypePair, unsigned int > header = connectome->header;
		std::list< unsigned int > cellTypes;
		cellTypes.push_back(L2);
		cellTypes.push_back(L34);
		cellTypes.push_back(L4py);
		cellTypes.push_back(L4sp);
		cellTypes.push_back(L4ss);
		cellTypes.push_back(L5st);
		cellTypes.push_back(L5tt);
		cellTypes.push_back(L6cc);
		cellTypes.push_back(L6ccinv);
		cellTypes.push_back(L6ct);
		cellTypes.push_back(L1);
		cellTypes.push_back(SymLocal1);
		cellTypes.push_back(SymLocal2);
		cellTypes.push_back(SymLocal3);
		cellTypes.push_back(SymLocal4);
		cellTypes.push_back(SymLocal5);
		cellTypes.push_back(SymLocal6);
		cellTypes.push_back(L23Trans);
		cellTypes.push_back(L45Sym);
		cellTypes.push_back(L45Peak);
		cellTypes.push_back(L56Trans);
		std::list< unsigned int > columns;
		columns.push_back(C2);
		std::vector< CellTableRow * > rows = connectome->getPostColumnCelltypeRows(columns, cellTypes);

		// Read Morphologies
		const char * inputFilename = "/nas1/Data_daniel/Network/INColumn_v6/morphologies/Morphologies.bb";
		std::map< unsigned int, unsigned int > prePostIDMap;
		std::vector< unsigned int > originalGraphIndices;
		std::vector< unsigned int > cellTypeIDs;
		std::vector< double * > spatialGraphTransforms;
		std::vector< std::string > originalGraphFiles;
		std::map< unsigned int, std::string > cellTypeIDLabels;
		std::vector< AmiraSpatialGraph * > readSpatialGraphs;
		Reader::readSpatialGraphSetFile(inputFilename, originalGraphIndices, cellTypeIDs, spatialGraphTransforms, originalGraphFiles, cellTypeIDLabels);

		CSVwriter << "Filename,ID" << std::endl;
		//std::cout << originalGraphFiles.size() << std::endl;
		//std::cout << originalGraphIndices.size() << std::endl;
		//int mx = 0;

		for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
		{
			if ((*rowit)->insideColumn==0)
				continue;
			CSVwriter << originalGraphFiles[originalGraphIndices[(*rowit)->cellID]] << "," << (*rowit)->cellID << std::endl;
			//if (mx<(*rowit)->cellID)
			//	mx = (*rowit)->cellID;
		}

		//std::cout << mx << std::endl;

//		CSVwriter << "Filename,ID" << std::endl;
//		for(int i = 0; i < originalGraphIndices.size(); ++i)
//		{
//			CSVwriter << originalGraphFiles[originalGraphIndices[i]] << "," << i << std::endl;
//		}

		CSVwriter.close();
		std::cout << "Written " << fname << std::endl;
	}
}

void getConvergence(CellTable * convergence)
{
	getConvergence(convergence,"Total",false);
}

void getConvergence(CellTable * convergence, std::string compartment)
{
	getConvergence(convergence,compartment,false);
}

void getConvergence(CellTable * convergence, std::string compartment, bool LocalSubTypes)
{
	initializeConstants();
	initializeLists();
	initalizeSomaPositions();

	std::string inputFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/convergence" + compartment + ".csv";
	std::string numCellsfname = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/nrCells.csv";

	getConvergence(convergence, inputFilename.c_str(), numCellsfname.c_str(), LocalSubTypes);
}

void getConvergence(CellTable * convergence, const char * inputfname, const char * numCellsfname,  bool LocalSubTypes)
{
	initializeConstants();
	initializeLists();
	initalizeSomaPositions();

	// Read in Cell Table
	Reader * cellTableReader = new Reader(inputfname, inputfname);
	cellTableReader->readConvergenceTableCSV(convergence);

	if (!LocalSubTypes)
	{
		convergence->mergeLocal(numCellsfname);
		std::cout << "Merge Locals! " << std::endl;
	}

	std::cout << "Done Reading Cell Table " << inputfname << std::endl;
	std::cout << " " << std::endl;
}

void getConnectome(CellTable * connectome)
{
	initializeConstants();
	initializeLists();
	initalizeSomaPositions();

	std::cout << "CHECK THIS?!?! getConnectome function" << std::endl;

	const char * inputFilename = ""; //"/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/convergence.csv";

	// Read in Cell Table
	Reader * cellTableReader = new Reader(inputFilename, inputFilename);
	cellTableReader->readSynapsesPerCellTable(connectome);
	const char * inputFilename2 = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/nrCells.csv";
	connectome->mergeLocal(inputFilename2);
}


// Compute convergence for VPM - IN wihtin C2 column along z-axis in 2D
// (vertical column depth z vs xy-horizontal distance to columna axis
void convergenceDepth2D()
{
	CellTable * convergence = new CellTable;
	getConvergence(convergence);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;

	// Extract CellType/Column Combination (postsynaptic)
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;
	std::list< unsigned int > cellTypes;
	cellTypes.push_back(SymLocal);
	cellTypes.push_back(L23Trans);
	cellTypes.push_back(L45Sym);
	cellTypes.push_back(L45Peak);
	cellTypes.push_back(L56Trans);
	std::list< unsigned int > columns;
	columns.push_back(C2);

	std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

	// PreCellType/Column
	int pretype = VPM;
	int precolumn = C2;
	double cellNorm = 0.0;
	std::list< unsigned int > celltypesPre (1,pretype);
	std::list< unsigned int > columnsPre (1,precolumn);
	cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
	cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());

	// Bins
	std::map <MatrixIndexType, bin*> binMap;
	double binsz = 50;

	for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
	{
		//if ((*rowit)->insideColumn==0)
		//	continue;

		double t, projectedPt[3], soma_postsyn[3];
		soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
		soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
		soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z
		double dist1 = 676-soma_postsyn[Z_COORD];
		if (dist1<0)
		{
			std::cout << "dist1 (" << dist1 << ") is negative!" << std::endl;
			return;
		}

		double dist2 = vtkLine::DistanceToLine(soma_postsyn, colRef->top, colRef->bottom, t, projectedPt);
		dist2 = sqrt(dist2);

		int pos1 = (int) dist1/binsz;
		int pos2 = (int) dist2/binsz;
		MatrixIndexType idx = MatrixIndexType(pos1,pos2);

		std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
		{
			// Presynaptic Cells
			if (headerit->first.first==precolumn && headerit->first.second==pretype)
			{
				std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
				double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
				float cval = (*Convit) * cellWeight;

				// Add to histogram if not in histogram yet
				if (binMap.count(idx) == 0)
				{
					bin * b = new bin;
					b->M = cval;
					b->N = 1;
					b->Q = 0;
					binMap[idx] = b;
				}
				else // Update histogram
				{
					// Update Mean, SD, and Number of Samples,
					float tmp1 = cval - binMap[idx]->M;
					// Update Mean of correlation r: (Mean * N + t) / (N+1)
					binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
					// Update N (Number of Samples)
					binMap[idx]->N = binMap[idx]->N+1;
					float tmp2 = (cval - binMap[idx]->M);
					// Update SD (can be calculated from Q)
					binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
				}
			}
		}
	}

	// Output
	std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/2DtestVPM-IN-C2.csv";
	storeBinMapAsCSV(outputFilename.c_str(), binMap, binsz);
}

// Compute convergence for VPM - IN wihtin C2 column along cortical column z-axis in 1D
void convergenceDepth1D()
{
	CellTable * convergence = new CellTable;
	getConvergence(convergence);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	// Extract CellType/Column Combination (postsynaptic)
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;
	std::list< unsigned int > cellTypes;
	cellTypes.push_back(SymLocal);
	cellTypes.push_back(L23Trans);
	cellTypes.push_back(L45Sym);
	cellTypes.push_back(L45Peak);
	cellTypes.push_back(L56Trans);
	std::list< unsigned int > columns;
	columns.push_back(C2);

	std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

	// PreCellType/Column
	int pretype = VPM;
	int precolumn = C2;
	double cellNorm = 0.0;
	std::list< unsigned int > celltypesPre (1,pretype);
	std::list< unsigned int > columnsPre (1,precolumn);
	cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
	cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());

//	std::cout << "cellNorm " << cellNorm << std::endl; // 350

	// Bins
	std::map <int , bin*> binMap;
	double binsz = 50;

	for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
	{
		if ((*rowit)->insideColumn==0)
			continue;

		double soma_z = 676-(*rowit)->somaLocation[Z_COORD];
		int idx = (int) soma_z/binsz;

		if (idx<0)
		{
			std::cout << "Negative Value! " << idx << " " << soma_z << std::endl;
			idx = 0;
		}

		std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
		{
			// Presynaptic Cells
			if (headerit->first.first==precolumn && headerit->first.second==pretype)
			{
				std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
				double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
				float cval = (*Convit) * cellWeight;

//				std::cout << "cellWeight " << cellWeight << std::endl; // 1
//				std::cout << "cval " << cval << std::endl; // 0.028
//				std::cout << "(*Convit) " << (*Convit) << std::endl; // 0.028

				// Add to histogram if not in histogram yet
				if (binMap.count(idx) == 0)
				{
					bin * b = new bin;
					b->M = cval;
					b->N = 1;
					b->Q = 0;
					binMap[idx] = b;
				}
				else // Update histogram
				{
					// Update Mean, SD, and Number of Samples,
					float tmp1 = cval - binMap[idx]->M;
					// Update Mean of correlation r: (Mean * N + t) / (N+1)
					binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
					// Update N (Number of Samples)
					binMap[idx]->N = binMap[idx]->N+1;
					float tmp2 = (cval - binMap[idx]->M);
					// Update SD (can be calculated from Q)
					binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
				}
			}
		}
	}

	// Output
	std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/testVPM-IN-C2.csv";
	std::ofstream CSVwriter;
	CSVwriter.open(outputFilename.c_str());

	if(!CSVwriter.fail())
	{
		CSVwriter << "x,M,SD,N_pre" << std::endl;

		for (std::map< int, bin *>::iterator it = binMap.begin(); it != binMap.end(); ++it)
		{
			bin * b = it->second;
			int N = b->N;
			float M = b->M;
			float Q = b->Q;
			float SD = sqrt(Q/(N-1));

			int idx = (it->first);

			CSVwriter << idx * binsz+binsz/2 << "," << M << "," << SD << "," << N << std::endl;
		}
	}
	else
	{
		std::cout << "Error! Could not write " << outputFilename << std::endl;
	}
	CSVwriter.close();
	std::cout << "Written " << outputFilename << std::endl;
}

/* Connection Probability over
 * - x: position of postsynaptic soma in x (row)
 * - y: position of postsynaptic soma in y (arc)
 * Presynapse:
 * - default: C2 (CellType)
 * - rowish: B2 C2 D2 (CellType)
 * - arcish: C1 C2 C3 (CellType)
 * Postsynapse: 3x3 (EX/IN)
 * Compartment: "Apical", "Basal" or "Total"
 * Spatial: "rowish", "arcish" or something else (=C2)
 * Basal and Apical might make no sense here, use Total
 * Stores CP for each Cell Type Combintion
 * - CPMap for x vs. y */
void writeCPMap_C2OutputEXIN(std::string compartment, std::string spatial)
{
	if ((compartment.compare("Apical")!=0) && (compartment.compare("Basal")!=0) && (compartment.compare("Total")!=0))
	{
		std::cout << "Wrong Input '" << compartment << "' only 'Basal' or 'Apical' or 'Total' allowed!" << std::endl;
		return;
	}

	CellTable * convergence = new CellTable;
	getConvergence(convergence,compartment);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;
	double binsz = 50;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	// Iterate over Postsynaptic Types
	for (int i = 0; i<2; ++i)
	{
		cellTypes.clear();

		if (i==0) // Excitatory
		{
			cellTypes.push_back(L2);
			cellTypes.push_back(L34);
			cellTypes.push_back(L4py);
			cellTypes.push_back(L4sp);
			cellTypes.push_back(L4ss);
			cellTypes.push_back(L5st);
			cellTypes.push_back(L5tt);
			cellTypes.push_back(L6cc);
			cellTypes.push_back(L6ccinv);
			cellTypes.push_back(L6ct);
		}
		else if (i==1) // Inhibitory
		{
			if ((compartment.compare("Apical")==0) || (compartment.compare("Basal")==0))
				continue;

			cellTypes.push_back(SymLocal);
			cellTypes.push_back(L1);
			cellTypes.push_back(L23Trans);
			cellTypes.push_back(L45Sym);
			cellTypes.push_back(L45Peak);
			cellTypes.push_back(L56Trans);
		}

		std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

		// Iterate over Presynaptic Types
		for (std::vector<int>::iterator Preit = PreCelltypeList.begin() ; Preit != PreCelltypeList.end(); ++Preit)
		{
			// Presynaptic Cell Type and its Location
			int pretype = *Preit;

			// Connection Probability Map
			std::map <MatrixIndexType, bin*> binMap;
			std::map <int, bin*> binZ;

			// Only one Celltype and Column C2
			std::list< unsigned int > celltypesPre (1,pretype);
			std::list< unsigned int > columnsPre (1,C2);

			if (spatial.compare("arcish")==0)
			{
				columnsPre.push_back(B2);
				columnsPre.push_back(D2);
			}
			if (spatial.compare("rowish")==0)
			{
				columnsPre.push_back(C1);
				columnsPre.push_back(C3);
			}

			double cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
			cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());

			// Compute Connection Probability, go through all post-synaptic cells
			for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
			{
				double t, projectedPt[3], soma_postsyn[3];
				soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
				soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
				soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

				// Compute Position of Post-Synaptic Cell with respect to anatomical landmarks
				double offset = 2000;
				double distx = soma_postsyn[X_COORD] + offset; // Row
				double disty = soma_postsyn[Y_COORD] + offset; // Arc
				double distz = soma_postsyn[Z_COORD] + offset; // z

				if (distx<0)
				{
					std::cout << "distx (" << distx << ") is negative, increase offset to " << offset-distx << "!" << std::endl;
					return;
				}
				if (disty<0)
				{
					std::cout << "disty (" << disty << ") is negative, increase offset to " << offset-disty << "!" << std::endl;
					return;
				}
				if (distz<0)
				{
					std::cout << "distz (" << distz << ") is negative, increase offset to " << offset-distz << "!" << std::endl;
					return;
				}

				int posx = (int) distx/binsz;
				int posy = (int) disty/binsz;
				int posz = (int) distz/binsz;

				MatrixIndexType idx = MatrixIndexType(posx,posy);
				std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

				for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
				{
					// Presynaptic Cells
					bool bolCol = (std::find(columnsPre.begin(), columnsPre.end(), headerit->first.first) != columnsPre.end());

					if (bolCol && headerit->first.second==pretype)
					{
						double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
						std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
						float cval = (*Convit) * cellWeight;

						// Row Arc Position
						// Add to histogram if not in histogram yet
						if (binMap.count(idx) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binMap[idx] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binMap[idx]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
							// Update N (Number of Samples)
							binMap[idx]->N = binMap[idx]->N+1;
							float tmp2 = (cval - binMap[idx]->M);
							// Update SD (can be calculated from Q)
							binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}

						/* 1D Map: z profile */
						if (binZ.count(posz) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binZ[posz] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binZ[posz]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binZ[posz]->M = (binZ[posz]->M * binZ[posz]->N + cval) / (binZ[posz]->N+1);
							// Update N (Number of Samples)
							binZ[posz]->N = binZ[posz]->N+1;
							float tmp2 = (cval - binZ[posz]->M);
							// Update SD (can be calculated from Q)
							binZ[posz]->Q = binZ[posz]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}
					}
				}
			}

			// Output
			std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/C2Output/" + std::string(int2CelltypeLabels[pretype]);

			if (spatial.compare("arcish")==0)
			{
				outputFilename = outputFilename + "_arc_";
			}
			else if (spatial.compare("rowish")==0)
			{
				outputFilename = outputFilename + "_row_";
			}
			else
			{
				outputFilename = outputFilename + "_C2_";
			}

			if (i==0)
			{
				outputFilename = outputFilename + "EX" + compartment;
			}
			else if (i==1)
			{
				outputFilename = outputFilename + "IN" + compartment;
			}

			std::string outputFilename1 = outputFilename + "_Z.csv";
			outputFilename = outputFilename + ".csv";

			storeBinMapAsCSV(outputFilename.c_str(), binMap, binsz);
			storeBinAsCSV(outputFilename1.c_str(), binZ, binsz); // z profile
		}
	}
}

/* Connection Probability over
 * - x: position of postsynaptic soma in x (row)
 * - y: position of postsynaptic soma in y (arc)
 * Presynapse:
 * - default: C2 (CellType)
 * - rowish: B2 C2 D2 (CellType)
 * - arcish: C1 C2 C3 (CellType)
 * Postsynapse: 3x3 (CellType)
 * Compartment: "Apical", "Basal" or "Total"
 * Spatial: "rowish", "arcish" or something else (=C2)
 * Stores CP for each Cell Type Combintion
 * - CPMap for x vs. y */
void writeCPMap_C2Output(std::string compartment, std::string spatial)
{
	if ((compartment.compare("Apical")!=0) && (compartment.compare("Basal")!=0) && (compartment.compare("Total")!=0))
	{
		std::cout << "Wrong Input '" << compartment << "' only 'Basal' or 'Apical' or 'Total' allowed!" << std::endl;
		return;
	}

	// Read in Cell Table
	CellTable * convergence = new CellTable;
	getConvergence(convergence,compartment);

	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;
	double binsz = 50;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	// Iterate over Postsynaptic Types
	for (std::vector<int>::iterator Postit = PostCelltypeList.begin() ; Postit != PostCelltypeList.end(); ++Postit)
	{
		int posttype = *Postit;

		// Skip if Inhibitory Cell and Apical or Basal is required!
		if (helper::isInhibitory(posttype) && ((compartment.compare("Apical")==0) || (compartment.compare("Basal")==0)))
			continue;

		cellTypes.clear();
		cellTypes.push_back(posttype);
		std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

		// Iterate over Presynaptic Types
		for (std::vector<int>::iterator Preit = PreCelltypeList.begin() ; Preit != PreCelltypeList.end(); ++Preit)
		{
			// Presynaptic Cell Type and its Location
			int pretype = *Preit;

			// Connection Probability Map
			std::map <MatrixIndexType, bin*> binMap;
			std::map <int, bin*> binZ;

			// Only one Celltype and Column C2
			std::list< unsigned int > celltypesPre (1,pretype);
			std::list< unsigned int > columnsPre (1,C2);

			if (spatial.compare("arcish")==0)
			{
				columnsPre.push_back(B2);
				columnsPre.push_back(D2);
			}
			if (spatial.compare("rowish")==0)
			{
				columnsPre.push_back(C1);
				columnsPre.push_back(C3);
			}

			double cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
			cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());

			// Compute Connection Probability, go through all post-synaptic cells
			for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
			{
				double t, projectedPt[3], soma_postsyn[3];
				soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
				soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
				soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

				// Compute Position of Post-Synaptic Cell with respect to anatomical landmarks
				double offset = 2000;
				double distx = soma_postsyn[X_COORD] + offset;
				double disty = soma_postsyn[Y_COORD] + offset;
				double distz = soma_postsyn[Z_COORD] + offset;

				if (distx<0)
				{
					std::cout << "distx (" << distx << ") is negative, increase offset to " << offset-distx << "!" << std::endl;
					return;
				}
				if (disty<0)
				{
					std::cout << "disty (" << disty << ") is negative, increase offset to " << offset-disty << "!" << std::endl;
					return;
				}
				if (distz<0)
				{
					std::cout << "distz (" << distz << ") is negative, increase offset to " << offset-distz << "!" << std::endl;
					return;
				}

				int posx = (int) distx/binsz;
				int posy = (int) disty/binsz;
				int posz = (int) distz/binsz;

				MatrixIndexType idx = MatrixIndexType(posx,posy);
				std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

				for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
				{
					// Presynaptic Cells
					bool bolCol = (std::find(columnsPre.begin(), columnsPre.end(), headerit->first.first) != columnsPre.end());

					if (bolCol && headerit->first.second==pretype)
					{
						double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
						std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
						float cval = (*Convit) * cellWeight;

						// Add to histogram if not in histogram yet
						if (binMap.count(idx) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binMap[idx] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binMap[idx]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
							// Update N (Number of Samples)
							binMap[idx]->N = binMap[idx]->N+1;
							float tmp2 = (cval - binMap[idx]->M);
							// Update SD (can be calculated from Q)
							binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}

						/* 1D Map: z profile */
						if (binZ.count(posz) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binZ[posz] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binZ[posz]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binZ[posz]->M = (binZ[posz]->M * binZ[posz]->N + cval) / (binZ[posz]->N+1);
							// Update N (Number of Samples)
							binZ[posz]->N = binZ[posz]->N+1;
							float tmp2 = (cval - binZ[posz]->M);
							// Update SD (can be calculated from Q)
							binZ[posz]->Q = binZ[posz]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}
					}
				}
			}

			// Output
			std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/C2Output/" + std::string(int2CelltypeLabels[pretype]);

			if (spatial.compare("arcish")==0)
			{
				outputFilename = outputFilename + "_arc_" + std::string(int2CelltypeLabels[posttype]) + compartment;
			}
			else if (spatial.compare("rowish")==0)
			{
				outputFilename = outputFilename + "_row_" + std::string(int2CelltypeLabels[posttype]) + compartment;
			}
			else
			{
				outputFilename = outputFilename + "_C2_" + std::string(int2CelltypeLabels[posttype]) + compartment;
			}

			std::string outputFilename1 = outputFilename + "_Z.csv";
			outputFilename = outputFilename + ".csv";

			storeBinMapAsCSV(outputFilename.c_str(), binMap, binsz);
			storeBinAsCSV(outputFilename1.c_str(), binZ, binsz); // z profile
		}
	}
}

void D2Test()
{
	CellTable * convergence = new CellTable;
	getConvergence(convergence);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;
	double binsz = 50;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	columns.push_back(D2);
	std::list< unsigned int > cellTypes;
	cellTypes.push_back(L4ss);

	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;
	std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

	std::list< unsigned int > celltypesPre (1,L4ssaxon);
	std::list< unsigned int > columnsPre;
//	columnsPre.push_back(C1);
//	columnsPre.push_back(C2);
//	columnsPre.push_back(C3);
	columnsPre.push_back(D1);
	columnsPre.push_back(D2);
	columnsPre.push_back(D3);
	//columnsPre.push_back(B1);
	//columnsPre.push_back(B2);
	//columnsPre.push_back(B3);
	double cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
	cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());

	// Connection Probability Map
	std::map <MatrixIndexType, bin*> binMap;

	// Compute Connection Probability, go through all post-synaptic cells
	for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
	{
		double t, projectedPt[3], soma_postsyn[3];
		soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
		soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
		soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

		// Compute Position of Post-Synaptic Cell with respect to anatomical landmarks
		double offset = 2000;
		double dist1 = soma_postsyn[X_COORD] + offset;
		double dist2 = soma_postsyn[Y_COORD] + offset;

		if (dist1<0)
		{
			std::cout << "dist1 (" << dist1 << ") is negative, increase offset to " << offset-dist1 << "!" << std::endl;
		}
		if (dist2<0)
		{
			std::cout << "dist2 (" << dist2 << ") is negative, increase offset to " << offset-dist2 << "!" << std::endl;
		}

		int pos1 = (int) dist1/binsz;
		int pos2 = (int) dist2/binsz;
		MatrixIndexType idx = MatrixIndexType(pos1,pos2);

		std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
		{
			// Presynaptic Cells
//			bool bolCol = (headerit->first.first==C2) || (headerit->first.first==C1) || (headerit->first.first==C3) ||
//							(headerit->first.first==D2) || (headerit->first.first==D1) || (headerit->first.first==D3) ||
//							(headerit->first.first==B2) || (headerit->first.first==B1) || (headerit->first.first==B3);

//			bool bolCol = (headerit->first.first==D2);
			bool bolCol = (headerit->first.first==D2) || (headerit->first.first==D1) || (headerit->first.first==D3);

			if (headerit->first.second==L4ssaxon && bolCol)
			{
				double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;

				std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);

				float cval = (*Convit) * cellWeight;

				// Add to histogram if not in histogram yet
				if (binMap.count(idx) == 0)
				{
					bin * b = new bin;
					b->M = cval;
					b->N = 1;
					b->Q = 0;
					binMap[idx] = b;
				}
				else // Update histogram
				{
					// Update Mean, SD, and Number of Samples,
					float tmp1 = cval - binMap[idx]->M;
					// Update Mean of correlation r: (Mean * N + t) / (N+1)
					binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
					// Update N (Number of Samples)
					binMap[idx]->N = binMap[idx]->N+1;
					float tmp2 = (cval - binMap[idx]->M);
					// Update SD (can be calculated from Q)
					binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
				}
			}
		}
	}

	// Output
	std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/single/L4ssDRow-L4ssD2.csv";
	std::ofstream CSVwriter;
	CSVwriter.open(outputFilename.c_str());

	if(!CSVwriter.fail())
	{
		CSVwriter << "x,y,M,SD,N" << std::endl;

		for (std::map< MatrixIndexType, bin *>::iterator it = binMap.begin(); it != binMap.end(); ++it)
		{
			bin * b = it->second;
			int N = b->N;
			float M = b->M;
			float Q = b->Q;
			float SD = sqrt(Q/(N-1));

			MatrixIndexType idx = (it->first);

			CSVwriter << idx.first * binsz+binsz/2 << "," << idx.second * binsz+binsz/2 << "," << M << "," << SD << "," << N << std::endl;
		}
	}
	else
	{
		std::cout << "Error! Could not write " << outputFilename << std::endl;
	}
	CSVwriter.close();

	std::cout << "Written " << outputFilename << std::endl;
}

/* Set Values of ImgVol to Zero, that are outside an specified range
 * Inputfname: path/to/ImageVolume.am
 * Outputfname: Output pathname of modified ImgVol
 * spatial: whether along arc or row */
void cutImgVol(const char * Inputfname, const char * Outputfname, std::string spatial)
{
	// Definition of Plane to compute Distance
	double n[3]; // normal vector between centers of C1 and C3, and C2 Column Vector
	double v[3]; // from C1 to C3

	if ((spatial.compare("row")==0) || (spatial.compare("rowish")==0))
	{
		getPlane(n, v, "row");
	}
	else if ((spatial.compare("arc")==0) || (spatial.compare("arcish")==0))
	{
		getPlane(n, v, "arc");
	}
	else
	{
		std::cout << "Wrong Input " << spatial << " not defined. Use 'arc'/'arcish' or 'row'/'rowish'!" << std::endl;
		return;
	}

	Reader * ImgReader = new Reader(Inputfname, Outputfname);
	ImageDataPointerType volume = ImgReader->readScalarField();

	double box[6];
	double VoxelSize = 50;
	volume->GetBounds(box);
	double pt[3];

	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{
		pt[0] = box[0]+x*VoxelSize;

		for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
		{
			pt[1] = box[2]+y*VoxelSize;

			for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
			{
				pt[2] = box[4]+z*VoxelSize;

				double dist = vtkPlane::DistanceToPlane(v,n,pt);

				if (dist>200)
				{
					double *px = static_cast< double * >(volume->GetScalarPointer(x,y,z));
					(*px) = 0;
				}
			}
		}
	}
	ImgReader->writeScalarField(volume);
}


/* Computes Plane along 2-arc or C-row and returns normal vector n and spanning vector x
 * 	n[3] normal vector between centers of C1 and C3, and C2 Column Vector, or B2, C2, D2,
 * 	x[3] from C1 to C3 or B2 to D2
 *  orientation arc (B2,C2,D2) or row (C1,C2,C3) */
void getPlane(double n[3], double x[3], std::string spatial)
{
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;

	if ((spatial.compare("row")==0) || (spatial.compare("rowish")==0))
	{
		Column * C1col = avgColumns.find(C1)->second;
		Column * C2col = avgColumns.find(C2)->second;
		Column * C3col = avgColumns.find(C3)->second;

		double centerC1[3];
		double centerC2[3];
		double centerC3[3];
		C1col->getCenter(centerC1);
		C2col->getCenter(centerC2);
		C3col->getCenter(centerC3);

		double * C2top = C2col->top;
		double * C2bottom = C2col->bottom;
		double C2vec[3];
		C2vec[0] = C2bottom[0]-C2top[0];
		C2vec[1] = C2bottom[1]-C2top[1];
		C2vec[2] = C2bottom[2]-C2top[2];

		double u[3]; // from C1 to C3
		u[0] = centerC3[0]-centerC1[0];
		u[1] = centerC3[1]-centerC1[1];
		u[2] = centerC3[2]-centerC1[2];

		// cross product between points spanning the plane
		// normal vector between centers of C1 and C3, and C2 Column Vector
		n[0] = u[1]*C2vec[2] - u[2]*C2vec[1];
		n[1] = u[2]*C2vec[0] - u[0]*C2vec[2];
		n[2] = u[0]*C2vec[1] - u[1]*C2vec[0];
		double norm = (sqrt(pow(n[0],2.0) + pow(n[1],2.0) + pow(n[2],2.0)));
		n[0] = n[0]/norm;
		n[1] = n[1]/norm;
		n[2] = n[2]/norm;

		// from C1 to C3
		x[0] = (centerC3[0]+centerC1[0]+centerC2[0])/3;
		x[1] = (centerC3[1]+centerC1[1]+centerC2[1])/3;
		x[2] = (centerC3[2]+centerC1[2]+centerC2[2])/3;
	}
	else if ((spatial.compare("arc")==0) || (spatial.compare("arcish")==0)) // 2-arc
	{
		Column * B2col = avgColumns.find(B2)->second;
		Column * C2col = avgColumns.find(C2)->second;
		Column * D2col = avgColumns.find(D2)->second;

		double centerB2[3];
		double centerC2[3];
		double centerD2[3];
		B2col->getCenter(centerB2);
		C2col->getCenter(centerC2);
		D2col->getCenter(centerD2);

		double * C2top = C2col->top;
		double * C2bottom = C2col->bottom;
		double C2vec[3];
		C2vec[0] = C2bottom[0]-C2top[0];
		C2vec[1] = C2bottom[1]-C2top[1];
		C2vec[2] = C2bottom[2]-C2top[2];

		double u[3]; // from D2 to B2
		u[0] = centerD2[0]-centerB2[0];
		u[1] = centerD2[1]-centerB2[1];
		u[2] = centerD2[2]-centerB2[2];

		// cross product between points spanning the plane
		// normal vector between centers of B2 and D2, and C2 Column Vector
		n[0] = u[1]*C2vec[2] - u[2]*C2vec[1];
		n[1] = u[2]*C2vec[0] - u[0]*C2vec[2];
		n[2] = u[0]*C2vec[1] - u[1]*C2vec[0];
		double norm = (sqrt(pow(n[0],2.0) + pow(n[1],2.0) + pow(n[2],2.0)));
		n[0] = n[0]/norm;
		n[1] = n[1]/norm;
		n[2] = n[2]/norm;

		// from C1 to C3
		x[0] = (centerD2[0]+centerB2[0]+centerC2[0])/3;
		x[1] = (centerD2[1]+centerB2[1]+centerC2[1])/3;
		x[2] = (centerD2[2]+centerB2[2]+centerC2[2])/3;
	}
	else
	{
		std::cout << "Wrong Input " << spatial << " not defined. Use 'arc'/'arcish' or 'row'/'rowish'!" << std::endl;
	}
}

/* Write Landmarkfile of all Somata in certain range (to check arcish/rowish)
 * Also Displays ratio of Cells within Distance and outside of Distance for each Column
 * spatial: arc or row
 * fname: path/to/landmark file, where landmark file containing soma positions should be stores
 * 	e.g., "/nas1/Data_daniel/CRow"; "/nas1/Data_daniel/2Arc"; */
void columnsCutDisplay(std::string spatial, const char * fname)
{
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;

	double dist_threshold = 200.0;

	// Definition of Plane to compute Distance
	double n[3]; // normal vector between centers of C1 and C3, and C2 Column Vector
	double x[3]; // from C1 to C3

	if ((spatial.compare("row")==0) || (spatial.compare("rowish")==0))
	{
		getPlane(n, x, "row");
	}
	else if ((spatial.compare("arc")==0) || (spatial.compare("arcish")==0))
	{
		getPlane(n, x, "arc");
	}
	else
	{
		std::cout << "Wrong Input " << spatial << " not defined. Use 'arc'/'arcish' or 'row'/'rowish'!" << std::endl;
		return;
	}

	std::map< unsigned int, unsigned int > countCol;
	std::map< unsigned int, unsigned int > countColTotal;

	/* Test Distance */
	CellTable * convergence = new CellTable;
	getConvergence(convergence);

	std::list< unsigned int > columns;
	std::list< unsigned int > cellTypes;
	std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

	PointsPointerType pt = PointsPointerType::New();
	pt->Allocate(1);
	pt->SetDataTypeToFloat();

	// Compute Connection Probability, go through all post-synaptic cells
	for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
	{
		double t, projectedPt[3], soma_postsyn[3];
		soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
		soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
		soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

		double dist = vtkPlane::DistanceToPlane(x,n,soma_postsyn);

		// Total Cells
		unsigned int currCol = (*rowit)->column;
		if (countColTotal.count(currCol) == 0)
		{
			countColTotal[currCol] = 1;
		}
		else
		{
			countColTotal[currCol] = countColTotal[currCol]+1;
		}

		if (dist<dist_threshold)
		{
			unsigned int oldNrOfPoints = pt->GetNumberOfPoints();
			pt->InsertPoint(oldNrOfPoints, soma_postsyn);

			if (countCol.count(currCol) == 0)
			{
				countCol[currCol] = 1;
			}
			else
			{
				countCol[currCol] = countCol[currCol]+1;
			}
		}
	}

	Reader * Writer = new Reader(fname, fname);
	Writer->writeLandmarkFile(pt);

	std::cout << "Column,Cells in range " << dist_threshold << std::endl;
	for (std::map< unsigned int, unsigned int>::iterator it = countColTotal.begin(); it != countColTotal.end(); ++it)
	{
		double n = 0.0;
		if (countCol.count(it->first) != 0)
		{
			n = countCol[it->first];
			n = n/(it->second);
		}

		std::cout << int2ColumnLabels[it->first] << "," << n << std::endl;
	}
}

/* Write CP Map for whole 3x3 grid
 * Postsynaptic: 3x3 (CellType)
 * Presynaptic: 3x3 or entire barrel cortex (EX/IN) (use wroteCPMap3x3 if no ex or in)
 * bol3x3: normalize presynaptic 3x3 (true) or across entire barrel cortex (false)
 * compartment: Basal, Apical or Total
 * Stores CP for each CellType Combination
 * - CPMap for x vs y */
void writeCPMap3x3EXIN(bool bol3x3, std::string compartment)
{
	if ((compartment.compare("Apical")!=0) && (compartment.compare("Basal")!=0) && (compartment.compare("Total")!=0))
	{
		std::cout << "Wrong Input '" << compartment << "' only 'Basal' or 'Apical' or 'Total' allowed!" << std::endl;
		return;
	}

	CellTable * convergence = new CellTable;
	getConvergence(convergence,compartment);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;
	double binsz = 50;
	double offset = 2000;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	// Iterate over Postsynaptic Types
	for (std::vector<int>::iterator Postit = PostCelltypeList.begin() ; Postit != PostCelltypeList.end(); ++Postit)
	{
		int posttype = *Postit;

		// Skip if Inhibitory Cell and Apical or Basal is required!
		if (helper::isInhibitory(posttype) && ((compartment.compare("Apical")==0) || (compartment.compare("Basal")==0)))
			continue;

		cellTypes.clear();
		cellTypes.push_back(posttype);
		std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

		// Iterate over Exc. Types (i=0) and Inh. Types (i=1)
		for (int i=0; i<2; ++i)
		{
			// Connection Probability Map
			std::map <MatrixIndexType, bin*> binMap;

			std::list< unsigned int > celltypesPre;
			std::list< unsigned int > columnsPre;

			if (bol3x3) // Only 3x3
			{
				columnsPre.push_back(C1);
				columnsPre.push_back(C2);
				columnsPre.push_back(C3);
				columnsPre.push_back(D1);
				columnsPre.push_back(D2);
				columnsPre.push_back(D3);
				columnsPre.push_back(B1);
				columnsPre.push_back(B2);
				columnsPre.push_back(B3);
			}
			else // All Columns
			{
				columnsPre.push_back(A1);
				columnsPre.push_back(A2);
				columnsPre.push_back(A3);
				columnsPre.push_back(A4);
				columnsPre.push_back(B1);
				columnsPre.push_back(B2);
				columnsPre.push_back(B3);
				columnsPre.push_back(B4);
				columnsPre.push_back(C1);
				columnsPre.push_back(C2);
				columnsPre.push_back(C3);
				columnsPre.push_back(C4);
				columnsPre.push_back(D1);
				columnsPre.push_back(D2);
				columnsPre.push_back(D3);
				columnsPre.push_back(D4);
				columnsPre.push_back(E1);
				columnsPre.push_back(E2);
				columnsPre.push_back(E3);
				columnsPre.push_back(E4);
				columnsPre.push_back(Alpha);
				columnsPre.push_back(Beta);
				columnsPre.push_back(Gamma);
				columnsPre.push_back(Delta);
			}

			if (i==0) // Excitatory Presynaptic Cells
			{
				celltypesPre.push_back(L2axon);
				celltypesPre.push_back(L34axon);
				celltypesPre.push_back(L4pyaxon);
				celltypesPre.push_back(L4spaxon);
				celltypesPre.push_back(L4ssaxon);
				celltypesPre.push_back(L5staxon);
				celltypesPre.push_back(L5ttaxon);
				celltypesPre.push_back(L6ccaxon);
				celltypesPre.push_back(L6ccinvaxon);
				celltypesPre.push_back(L6ctaxon);
			}
			else // Inhibitory Presynaptic Cells
			{
				celltypesPre.push_back(SymLocalaxon);
				celltypesPre.push_back(L1axon);
				celltypesPre.push_back(L23Transaxon);
				celltypesPre.push_back(L45Symaxon);
				celltypesPre.push_back(L45Peakaxon);
				celltypesPre.push_back(L56Transaxon);
			}

			double cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
			cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());

			// Compute Connection Probability, go through all post-synaptic cells
			for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
			{
				double t, projectedPt[3], soma_postsyn[3];
				soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
				soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
				soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

				// Compute Position of Post-Synaptic Cell with respect to anatomical landmarks
				double distx = soma_postsyn[X_COORD] + offset;
				double disty = soma_postsyn[Y_COORD] + offset;

				if (distx<0)
				{
					std::cout << "distx (" << distx << ") is negative, increase offset to " << offset-distx << "!" << std::endl;
					return;
				}
				if (disty<0)
				{
					std::cout << "disty (" << disty << ") is negative, increase offset to " << offset-disty << "!" << std::endl;
					return;
				}

				int posx = (int) distx/binsz;
				int posy = (int) disty/binsz;
				MatrixIndexType idx = MatrixIndexType(posx,posy);
				std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

				for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
				{
					// Check presynaptic column
					bool bolCol = true;
					if (bol3x3)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==C1) || (headerit->first.first==C3) ||
									(headerit->first.first==D2) || (headerit->first.first==D1) || (headerit->first.first==D3) ||
									(headerit->first.first==B2) || (headerit->first.first==B1) || (headerit->first.first==B3);
					}

					// Check presynaptic celltype
					bool bolCell = false;
					for (std::list<unsigned int>::iterator it=celltypesPre.begin(); it != celltypesPre.end(); ++it)
					{
						if (headerit->first.second == (*it))
						{
							bolCell = true;
							break;
						}
					}

					// Presynaptic Cells
					if (bolCol && bolCell) //headerit->first.second==pretype)
					{
						double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
						std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
						float cval = (*Convit) * cellWeight;

						// Add to histogram if not in histogram yet
						if (binMap.count(idx) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binMap[idx] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binMap[idx]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
							// Update N (Number of Samples)
							binMap[idx]->N = binMap[idx]->N+1;
							float tmp2 = (cval - binMap[idx]->M);
							// Update SD (can be calculated from Q)
							binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}
					}
				}
			}

			// Output
			std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/3x3/";
			if (i==0) // Excitatory Presynaptic Cells
			{
				outputFilename = outputFilename + "EX_" + std::string(int2CelltypeLabels[posttype]) + compartment;
			}
			else // Inhibitory Presynaptic Cells
			{
				outputFilename = outputFilename + "IN_" + std::string(int2CelltypeLabels[posttype]) + compartment;
			}
			if (bol3x3)
			{
				outputFilename = outputFilename + "_3x3";
			}
			outputFilename = outputFilename + ".csv";
			storeBinMapAsCSV(outputFilename.c_str(), binMap, binsz);
		}
	}
}


/* Write CP Map for whole 3x3 grid
 * Presynaptic: Entire Barrel Cortex or 3x3
 * Postsynaptic: 3x3
 * bol3x3: normalize presynaptic 3x3 (true) or across entire barrel cortex (false)
 * compartment: Basal, Apical or Total */
void writeCPMap3x3(bool bol3x3, std::string compartment)
{
	if ((compartment.compare("Apical")!=0) && (compartment.compare("Basal")!=0) && (compartment.compare("Total")!=0))
	{
		std::cout << "Wrong Input '" << compartment << "' only 'Basal' or 'Apical' or 'Total' allowed!" << std::endl;
		return;
	}

	CellTable * convergence = new CellTable;
	getConvergence(convergence,compartment);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;
	double binsz = 50;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	// Iterate over Postsynaptic Types
	for (std::vector<int>::iterator Postit = PostCelltypeList.begin() ; Postit != PostCelltypeList.end(); ++Postit)
	{
		int posttype = *Postit;

		// Skip if Inhibitory Cell and Apical or Basal is required!
		if (helper::isInhibitory(posttype) && ((compartment.compare("Apical")==0) || (compartment.compare("Basal")==0)))
			continue;

		cellTypes.clear();
		cellTypes.push_back(posttype);
		std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

		// Iterate over Presynaptic Types
		for (std::vector<int>::iterator Preit = PreCelltypeList.begin() ; Preit != PreCelltypeList.end(); ++Preit)
		{
			// Presynaptic Cell Type and its Location
			int pretype = *Preit;

			// Connection Probability Map
			std::map <MatrixIndexType, bin*> binMap;

			double cellNorm = 0.0;
			if (bol3x3)
			{
				std::list< unsigned int > celltypesPre;
				celltypesPre.push_back(pretype);
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(C1);
				columnsPre.push_back(C2);
				columnsPre.push_back(C3);
				columnsPre.push_back(D1);
				columnsPre.push_back(D2);
				columnsPre.push_back(D3);
				columnsPre.push_back(B1);
				columnsPre.push_back(B2);
				columnsPre.push_back(B3);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}
			else
			{
				std::list< unsigned int > celltypesPre;
				celltypesPre.push_back(pretype);
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(A1);
				columnsPre.push_back(A2);
				columnsPre.push_back(A3);
				columnsPre.push_back(A4);
				columnsPre.push_back(B1);
				columnsPre.push_back(B2);
				columnsPre.push_back(B3);
				columnsPre.push_back(B4);
				columnsPre.push_back(C1);
				columnsPre.push_back(C2);
				columnsPre.push_back(C3);
				columnsPre.push_back(C4);
				columnsPre.push_back(D1);
				columnsPre.push_back(D2);
				columnsPre.push_back(D3);
				columnsPre.push_back(D4);
				columnsPre.push_back(E1);
				columnsPre.push_back(E2);
				columnsPre.push_back(E3);
				columnsPre.push_back(E4);
				columnsPre.push_back(Alpha);
				columnsPre.push_back(Beta);
				columnsPre.push_back(Gamma);
				columnsPre.push_back(Delta);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}

			// Compute Connection Probability, go through all post-synaptic cells
			for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
			{
				double t, projectedPt[3], soma_postsyn[3];
				soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
				soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
				soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

				// Compute Position of Post-Synaptic Cell with respect to anatomical landmarks
				double offset = 2000;
				double dist1 = soma_postsyn[X_COORD] + offset;
				double dist2 = soma_postsyn[Y_COORD] + offset;

				if (dist1<0)
				{
					std::cout << "dist1 (" << dist1 << ") is negative, increase offset to " << offset-dist1 << "!" << std::endl;
					return;
				}
				if (dist2<0)
				{
					std::cout << "dist1 (" << dist2 << ") is negative, increase offset to " << offset-dist2 << "!" << std::endl;
					return;
				}

				int pos1 = (int) dist1/binsz;
				int pos2 = (int) dist2/binsz;
				MatrixIndexType idx = MatrixIndexType(pos1,pos2);

				std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

				for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
				{
					bool bolCol = true;
					if (bol3x3)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==C1) || (headerit->first.first==C3) ||
									(headerit->first.first==D2) || (headerit->first.first==D1) || (headerit->first.first==D3) ||
									(headerit->first.first==B2) || (headerit->first.first==B1) || (headerit->first.first==B3);
					}

					// Presynaptic Cells
					if (bolCol && headerit->first.second==pretype)
					{
						double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
						std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
						float cval = (*Convit) * cellWeight;

						// Add to histogram if not in histogram yet
						if (binMap.count(idx) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binMap[idx] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binMap[idx]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
							// Update N (Number of Samples)
							binMap[idx]->N = binMap[idx]->N+1;
							float tmp2 = (cval - binMap[idx]->M);
							// Update SD (can be calculated from Q)
							binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}
					}
				}
			}

			// Output
			std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/3x3/" + std::string(int2CelltypeLabels[pretype])
											+ std::string(int2CelltypeLabels[posttype]) + compartment;
			if (bol3x3)
			{
				outputFilename = outputFilename + "_3x3";
			}
			outputFilename = outputFilename + ".csv";
			storeBinMapAsCSV(outputFilename.c_str(), binMap, binsz);
		}
	}
}

/* Write BinMap as CSV file
 * fname: path/to/filename
 * binMap: Matrix
 * binsz for position within matrix */
void storeBinMapAsCSV(const char* fname, std::map< MatrixIndexType, bin *> binMap, double binsz)
{
	std::ofstream CSVwriter;
	CSVwriter.open(fname);

	if(!CSVwriter.fail())
	{
		CSVwriter << "x,y,M,SD,N" << std::endl;

		for (std::map< MatrixIndexType, bin *>::iterator it = binMap.begin(); it != binMap.end(); ++it)
		{
			bin * b = it->second;
			int N = b->N;
			float M = b->M;
			float Q = b->Q;
			float SD = sqrt(Q/(N-1));

			MatrixIndexType idx = (it->first);

			CSVwriter << idx.first * binsz+binsz/2 << "," << idx.second * binsz+binsz/2 << "," << M << "," << SD << "," << N << std::endl;
		}
	}
	else
	{
		std::cout << "Error! Could not write " << fname << std::endl;
	}
	CSVwriter.close();

	std::cout << "Written " << fname << std::endl;
}

/* Write bin1D as CSV file
 * fname: path/to/filename
 * binMap: 1D bin
 * binsz for position within matrix
 */
void storeBinAsCSV(const char* fname, std::map< int, bin *> bin1D, double binsz)
{
	std::ofstream CSVwriter;
	CSVwriter.open(fname);

	if(!CSVwriter.fail())
	{
		CSVwriter << "x,M,SD,N" << std::endl;

		for (std::map< int, bin *>::iterator it = bin1D.begin(); it != bin1D.end(); ++it)
		{
			bin * b = it->second;
			int N = b->N;
			float M = b->M;
			float Q = b->Q;
			float SD = sqrt(Q/(N-1));

			int idx = (it->first);

			CSVwriter << idx * binsz+binsz/2 << "," << M << "," << SD << "," << N << std::endl;
		}
	}
	else
	{
		std::cout << "Error! Could not write " << fname << std::endl;
	}
	CSVwriter.close();

	std::cout << "Written " << fname << std::endl;
}

/* Write Single Connection Probability Map for specified pre- and postsynaptic ColumnCellType Pair
 * Stores Output over x and y position of postsynaptic ColumnCellType Pair */
void writeCPMap3x3Single(ColumnCellTypePair PreSynapse, ColumnCellTypePair PostSynapse)
{
	CellTable * convergence = new CellTable;
	getConvergence(convergence);

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;
	double binsz = 50;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	columns.push_back(PostSynapse.first);
	std::list< unsigned int > cellTypes;
	cellTypes.push_back(PostSynapse.second);

	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;
	std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

	// Connection Probability Map
	std::map <MatrixIndexType, bin*> binMap;

	// Compute Connection Probability, go through all post-synaptic cells
	for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
	{
		double t, projectedPt[3], soma_postsyn[3];
		soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
		soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
		soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

//		if (!((*rowit)->insideColumn))
//		{
//			continue;
//		}

		// Compute Position of Post-Synaptic Cell with respect to anatomical landmarks
		double offset = 2000;
		double dist1 = soma_postsyn[X_COORD] + offset;
		double dist2 = soma_postsyn[Y_COORD] + offset;

		if (dist1<0)
		{
			std::cout << "dist1 (" << dist1 << ") is negative, increase offset to " << offset-dist1 << "!" << std::endl;
			return;
		}
		if (dist2<0)
		{
			std::cout << "dist1 (" << dist2 << ") is negative, increase offset to " << offset-dist2 << "!" << std::endl;
			return;
		}

		int pos1 = (int) dist1/binsz;
		int pos2 = (int) dist2/binsz;
		MatrixIndexType idx = MatrixIndexType(pos1,pos2);

		std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
		{
			// Presynaptic Cells
			if (headerit->first.second==PreSynapse.second && headerit->first.first==PreSynapse.first)
			{
				std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);

				// Add to histogram if not in histogram yet
				if (binMap.count(idx) == 0)
				{
					bin * b = new bin;
					b->M = (*Convit);
					b->N = 1;
					b->Q = 0;
					binMap[idx] = b;
				}
				else // Update histogram
				{
					// Update Mean, SD, and Number of Samples,
					float tmp1 = (*Convit) - binMap[idx]->M;
					// Update Mean of correlation r: (Mean * N + t) / (N+1)
					binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + (*Convit)) / (binMap[idx]->N+1);
					// Update N (Number of Samples)
					binMap[idx]->N = binMap[idx]->N+1;
					float tmp2 = ((*Convit) - binMap[idx]->M);
					// Update SD (can be calculated from Q)
					binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
				}
			}
		}
	}

	// Output
	std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/single/" +
			std::string(int2CelltypeLabels[PreSynapse.second]) + "_" + std::string(int2ColumnLabels[PreSynapse.first]) + "-" +
			std::string(int2CelltypeLabels[PostSynapse.second]) + "_" + std::string(int2ColumnLabels[PostSynapse.first]) + ".csv";
	storeBinMapAsCSV(outputFilename.c_str(), binMap, binsz);
}

/* spatial: "arcish" or "rowish"
 * compartment: "Apical", "Basal", or empty "Total"
 * */
void writeCP1D_Cut(std::string spatial, std::string compartment)
{
	if ((spatial.compare("arcish")!=0) && (spatial.compare("rowish")!=0))
	{
		std::cout << "Wrong Input '" << spatial << "' only 'arcish' or 'rowish' allowed!" << std::endl;
		return;
	}

	if ((compartment.compare("Apical")!=0) && (compartment.compare("Basal")!=0) && (compartment.compare("Total")!=0))
	{
		std::cout << "Wrong Input '" << compartment << "' only 'Basal' or 'Apical' or 'Total' allowed!" << std::endl;
		return;
	}

	CellTable * convergence = new CellTable;
	getConvergence(convergence,compartment);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	double binsz = 50.0;

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;

	// Get Plane to compute distance from plane
	// Definition of Plane to compute Distance
	double n[3]; // normal vector between centers of C1 and C3, and C2 Column Vector
	double x[3]; // from C1 to C3
	getPlane(n, x, spatial);
	double distThreshold = 200.0;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	// Iterate over Postsynaptic Types
	for (std::vector<int>::iterator Postit = PostCelltypeList.begin() ; Postit != PostCelltypeList.end(); ++Postit)
	{
		int posttype = *Postit;
		cellTypes.clear();
		cellTypes.push_back(posttype);
		std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

		// Iterate over Presynaptic Types
		for (std::vector<int>::iterator Preit = PreCelltypeList.begin() ; Preit != PreCelltypeList.end(); ++Preit)
		{
			// Presynaptic Cell Type
			int pretype = *Preit;
			double cellNorm = 0.0;

			if (spatial.compare("rowish")==0)
			{
				std::list< unsigned int > celltypesPre (1,pretype);
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(C1);
				columnsPre.push_back(C2);
				columnsPre.push_back(C3);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}
			else if (spatial.compare("arcish")==0)
			{
				std::list< unsigned int > celltypesPre (1,pretype);
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(D2);
				columnsPre.push_back(C2);
				columnsPre.push_back(B2);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}

			// Connection Probability Map
			std::map <int , bin*> binMap;
			std::map <int , int> binN;

			// Compute Connection Probability, go through all post-synaptic cells
			for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
			{
				double t, projectedPt[3], soma_postsyn[3];
				soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
				soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
				soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

				double distPlane = vtkPlane::DistanceToPlane(x,n,soma_postsyn);

				if (distPlane>distThreshold)
				{
					continue;
				}

				// Compute Position of Post-Synaptic Cell with respect to anatomical landmarks
				double dist = 0;
				double Offset = 2000;
				if (spatial.compare("rowish")==0) // Position on Row (from D2)
				{
					dist = soma_postsyn[X_COORD] + Offset;
				}
				else if (spatial.compare("arcish")==0) // Position on Arc from (D2)
				{
					dist = soma_postsyn[Y_COORD] + Offset;
				}
				if (dist<0)
				{
					std::cout << dist << " is negative ! Increase artificial offset " << Offset << std::endl;
					return;
				}

				int idx = (int) dist/binsz;

				if (binN.count(idx) == 0)
				{
					binN[idx] = 1;
				}
				else // Update histogram
				{
					binN[idx] = binN[idx]+1;
				}

				std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

				for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
				{
					// Presynaptic Cells
					bool bolCol = false; // (headerit->first.first==C2) || (headerit->first.first==C1) || (headerit->first.first==C3) ||
								//	(headerit->first.first==D2) || (headerit->first.first==D1) || (headerit->first.first==D3) ||
								//	(headerit->first.first==B2) || (headerit->first.first==B1) || (headerit->first.first==B3);

					if (spatial.compare("rowish")==0)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==C1) || (headerit->first.first==C3);
					}
					else if (spatial.compare("arcish")==0)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==B2) || (headerit->first.first==D2);
					}

					// Presynaptic Cells
					if (bolCol && headerit->first.second==pretype)
					{
						std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
						double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
						float cval = (*Convit) * cellWeight;

						// Add to histogram if not in histogram yet
						if (binMap.count(idx) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binMap[idx] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binMap[idx]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
							// Update N (Number of Samples)
							binMap[idx]->N = binMap[idx]->N+1;
							float tmp2 = (cval - binMap[idx]->M);
							// Update SD (can be calculated from Q)
							binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}
					}
				}
			}

			// Output
			std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/1D/" + spatial + "_Cut/" + std::string(int2CelltypeLabels[pretype]) + "_" + std::string(int2CelltypeLabels[posttype]) + compartment + ".csv";
			if (compartment.compare("Total")==0)
				outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/1D/" + spatial + "_Cut/" + std::string(int2CelltypeLabels[pretype]) + "_" + std::string(int2CelltypeLabels[posttype]) + ".csv";

			std::ofstream CSVwriter;
			CSVwriter.open(outputFilename.c_str());

			if(!CSVwriter.fail())
			{
				CSVwriter << "x,M,SD,N_pre,N_post" << std::endl;

				std::map< int, int>::iterator itN = binN.begin();
				for (std::map< int, bin *>::iterator it = binMap.begin(); it != binMap.end(); ++it, ++itN)
				{
					bin * b = it->second;
					int N = b->N;
					float M = b->M;
					float Q = b->Q;
					float SD = sqrt(Q/(N-1));
					int idx = (it->first);

					CSVwriter << idx * binsz+binsz/2 << "," << M << "," << SD << "," << N << "," << itN->second << std::endl;
				}
			}
			else
			{
				std::cout << "Error! Could not write " << outputFilename << std::endl;
			}
			CSVwriter.close();
			std::cout << "Written " << outputFilename << std::endl;
		}
	}
}

/* Not sure whether this works, DEPRECATED
 * spatial: "arcish" or "rowish"
 * compartment: "Apical", "Basal", or empty "Total" */
void writeCP1D(std::string spatial, std::string compartment)
{
	if ((spatial.compare("arcish")!=0) && (spatial.compare("rowish")!=0))
	{
		std::cout << "Wrong Input '" << spatial << "' only 'arcish' or 'rowish' allowed!" << std::endl;
		return;
	}

	if ((compartment.compare("Apical")!=0) && (compartment.compare("Basal")!=0) && (compartment.compare("Total")!=0))
	{
		std::cout << "Wrong Input '" << compartment << "' only 'Basal' or 'Apical' or '' allowed!" << std::endl;
		return;
	}

	CellTable * convergence = new CellTable;
	getConvergence(convergence,compartment);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	double binsz = 50;

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	if (spatial.compare("arcish")==0)
	{
		columns.push_back(B2);
		columns.push_back(C2);
		columns.push_back(D2);
	}
	if (spatial.compare("rowish")==0)
	{
		columns.push_back(C1);
		columns.push_back(C2);
		columns.push_back(C3);
	}

	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	// Iterate over Postsynaptic Types
	for (std::vector<int>::iterator Postit = PostCelltypeList.begin() ; Postit != PostCelltypeList.end(); ++Postit)
	{
		int posttype = *Postit;
		cellTypes.clear();
		cellTypes.push_back(posttype);
		std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

		// Iterate over Presynaptic Types
		for (std::vector<int>::iterator Preit = PreCelltypeList.begin() ; Preit != PreCelltypeList.end(); ++Preit)
		{
			// Presynaptic Cell Type
			int pretype = *Preit;
			double cellNorm = 0.0;

			if (spatial.compare("rowish")==0)
			{
				std::list< unsigned int > celltypesPre (1,pretype);
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(C1);
				columnsPre.push_back(C2);
				columnsPre.push_back(C3);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}
			else if (spatial.compare("arcish")==0)
			{
				std::list< unsigned int > celltypesPre (1,pretype);
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(D2);
				columnsPre.push_back(C2);
				columnsPre.push_back(B2);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}

			// Connection Probability Map
			std::map <int , bin*> binMap;
			std::map <int , int> binN;

			// Compute Connection Probability, go through all post-synaptic cells
			for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
			{
				double t, projectedPt[3], soma_postsyn[3];
				soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
				soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
				soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

				int SGI = BF->laminarPosition(soma_postsyn);

				//if ((SGI !=GRAN) || !((*rowit)->insideColumn))
//				if ((SGI !=SUPRA))
//					continue;

				// Compute Position of Post-Synaptic Cell with respect to anatomical landmarks
				double dist = 0;
				double Offset = 2000;
				if (spatial.compare("rowish")==0) // Position on Row (from D2)
				{
					dist = soma_postsyn[X_COORD] + Offset;
				}
				else if (spatial.compare("arcish")==0) // Position on Arc from (D2)
				{
					dist = soma_postsyn[Y_COORD] + Offset;
				}
				if (dist<0)
				{
					std::cout << dist << " is negative ! Increase artificial offset " << Offset << std::endl;
					return;
				}

				int idx = (int) dist/binsz;

				if (binN.count(idx) == 0)
				{
					binN[idx] = 1;
				}
				else // Update histogram
				{
					binN[idx] = binN[idx]+1;
				}

				std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

				for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
				{
					bool bolCol = false;
					if (spatial.compare("rowish")==0)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==C1) || (headerit->first.first==C3);
					}
					else if (spatial.compare("arcish")==0)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==B2) || (headerit->first.first==D2);
					}

					// Presynaptic Cells
					if (bolCol && headerit->first.second==pretype)
					{
						std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
						double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
						float cval = (*Convit) * cellWeight;

						// Add to histogram if not in histogram yet
						if (binMap.count(idx) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binMap[idx] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binMap[idx]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
							// Update N (Number of Samples)
							binMap[idx]->N = binMap[idx]->N+1;
							float tmp2 = (cval - binMap[idx]->M);
							// Update SD (can be calculated from Q)
							binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}
					}
				}
			}

			// Output
			std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/1D/" + spatial + "/" + std::string(int2CelltypeLabels[pretype]) + "_" + std::string(int2CelltypeLabels[posttype]) + compartment + ".csv";
			if (compartment.compare("Total")==0)
				outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/1D/" + spatial + "/" + std::string(int2CelltypeLabels[pretype]) + "_" + std::string(int2CelltypeLabels[posttype]) + ".csv";

			std::ofstream CSVwriter;
			CSVwriter.open(outputFilename.c_str());

			if(!CSVwriter.fail())
			{
				CSVwriter << "x,M,SD,N_pre,N_post" << std::endl;

				std::map< int, int>::iterator itN = binN.begin();
				for (std::map< int, bin *>::iterator it = binMap.begin(); it != binMap.end(); ++it, ++itN)
				{
					bin * b = it->second;
					int N = b->N;
					float M = b->M;
					float Q = b->Q;
					float SD = sqrt(Q/(N-1));

					int idx = (it->first);

					CSVwriter << idx * binsz+binsz/2 << "," << M << "," << SD << "," << N << "," << itN->second << std::endl;
				}
			}
			else
			{
				std::cout << "Error! Could not write " << outputFilename << std::endl;
			}
			CSVwriter.close();
			std::cout << "Written " << outputFilename << std::endl;
		}
	}
}

/* DEPRECATED */
/* Connection Probability vs. Distance of Postsynaptic Cell to vertical Column axis and Inter-Soma Distance pre/postsynaptic Cells either along row or arc (3x3) */
/* spatial = or "arcish" or "rowish" */
/* Postsynaptic Cell only in certain range around plane (not column, to avoid missing septa */
void writeCPMap3(std::string spatial)
{
	if ((spatial.compare("arcish")!=0) && (spatial.compare("rowish")!=0))
	{
		std::cout << "Wrong Input '" << spatial << "' only 'arcish', and 'rowish' allowed!" << std::endl;
		return;
	}

	CellTable * convergence = new CellTable;
	getConvergence(convergence);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;

	// Get acrish or rowish plane
	// Used to calculate the distance between soma and plane to find all postsynaptic cells within this range
	double n[3]; // normal vector between centers of C1 and C3, and C2 Column Vector
	double x[3]; // from C1 to C3
	getPlane(n, x, spatial);
	double distThreshold = 200.0;
	double binsz = 50;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	// Iterate over Postsynaptic Types
	for (std::vector<int>::iterator Postit = PostCelltypeList.begin() ; Postit != PostCelltypeList.end(); ++Postit)
	{
		int posttype = *Postit;
		cellTypes.clear();
		cellTypes.push_back(posttype);
		std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

		// Iterate over Presynaptic Types
		for (std::vector<int>::iterator Preit = PreCelltypeList.begin() ; Preit != PreCelltypeList.end(); ++Preit)
		{
			// Presynaptic Cell Type and its Location
			int pretype = *Preit;
			double cellNorm = 0.0;

			if (spatial.compare("rowish")==0)
			{
				std::list< unsigned int > celltypesPre (1,pretype);
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(C1);
				columnsPre.push_back(C2);
				columnsPre.push_back(C3);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}
			else if (spatial.compare("arcish")==0)
			{
				std::list< unsigned int > celltypesPre (1,pretype);
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(D2);
				columnsPre.push_back(C2);
				columnsPre.push_back(B2);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}

			// Connection Probability Map
			std::map <MatrixIndexType, bin*> binMap;

			// Compute Connection Probability, go through all post-synaptic cells
			for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
			{
				double t, projectedPt[3], soma_postsyn[3];
				soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
				soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
				soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

				// Checks whether soma is within range to arcish/rowish plane. If not, skip
				double distPlane = vtkPlane::DistanceToPlane(x,n,soma_postsyn);
				if (distPlane>distThreshold)
				{
					continue;
				}

				// Compute Position of Post-Synaptic Cell with respect to anatomical landmarks
				double dist1 = 0;
				double offset = 2000;
				if (spatial.compare("rowish")==0) // Position on Row (from D2)
				{
					dist1 = soma_postsyn[X_COORD] + offset;
				}
				else if (spatial.compare("arcish")==0) // Position on Arc from (D2)
				{
					dist1 = soma_postsyn[Y_COORD] + offset;
				}

				if (dist1<0)
				{
					std::cout << "dist1 (" << dist1 << ") is negative, increase offset to " << offset-dist1 << "!" << std::endl;
				}

				int pos1 = (int) dist1/binsz;
				std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

				for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
				{
					bool bolCol = false;
					if (spatial.compare("rowish")==0)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==C1) || (headerit->first.first==C3);
					}
					else if (spatial.compare("arcish")==0)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==B2) || (headerit->first.first==D2);
					}

					// Presynaptic Cells
					if (bolCol && headerit->first.second==pretype)
					{
						// Compute Distance Post-Synaptic Cell to Average Soma Position of Presynaptic Cell type in respective Column
						somaPosition soma_presyn = celltypeLabels2SomaPos[ColumnCellTypePair(headerit->first.first,pretype)];
						double dist2 = sqrt(pow((soma_presyn.x-soma_postsyn[X_COORD]),2.0)+pow((soma_presyn.y-soma_postsyn[Y_COORD]),2.0));
						int pos2 = (int) dist2/binsz;
						MatrixIndexType idx = MatrixIndexType(pos1,pos2);

						std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
						double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
						float cval = (*Convit) * cellWeight;

						// Add to histogram if not in histogram yet
						if (binMap.count(idx) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binMap[idx] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binMap[idx]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
							// Update N (Number of Samples)
							binMap[idx]->N = binMap[idx]->N+1;
							float tmp2 = (cval - binMap[idx]->M);
							// Update SD (can be calculated from Q)
							binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}
					}
				}
			}

			// Output
			std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/" + spatial + "/" + std::string(int2CelltypeLabels[pretype]) + "_" + std::string(int2CelltypeLabels[posttype]) + "_3D.csv";
			storeBinMapAsCSV(outputFilename.c_str(), binMap, binsz);
		}
	}
}

/* Connection Probability over
 * - x: Distance of Postsynaptic Cell to vertical Column axis
 * - y: Intersomatic Distance pre/postsynaptic Cells (Average, not real!)
 * spatial = "horizontal", "vertical", "arcish" or "rowish"
 * - horizontal/vertical: Presynapse: C2 (CellType); Postsynapse: 3x3 (CellType);
 * - rowish: Presynapse: C1, C2, C3 (CellType); Postsynapse: C1, C2, C3 (CellType);
 * - arcish: Presynapse: B2, C2, D2 (CellType); Postsynapse: B2, C2, D2 (CellType);
 * Stores CP for each Cell Type Combination
 * - CPMap for Distance to C2 vs InterSomaDistance */
void writeCPMap(std::string spatial)
{
	if ((spatial.compare("vertical")!=0) && (spatial.compare("horizontal")!=0) && (spatial.compare("arcish")!=0) && (spatial.compare("rowish")!=0))
	{
		std::cout << "Wrong Input '" << spatial << "' only 'vertical', 'horizontal', 'arcish', and 'rowish' allowed!" << std::endl;
		return;
	}

	CellTable * convergence = new CellTable;
	getConvergence(convergence);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;
	double binsz = 50;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	if (spatial.compare("arcish")==0)
	{
		columns.push_back(B2);
		columns.push_back(C2);
		columns.push_back(D2);
	}
	if (spatial.compare("rowish")==0)
	{
		columns.push_back(C1);
		columns.push_back(C2);
		columns.push_back(C3);
	}

	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	// Iterate over Postsynaptic Types
	for (std::vector<int>::iterator Postit = PostCelltypeList.begin() ; Postit != PostCelltypeList.end(); ++Postit)
	{
		int posttype = *Postit;
		cellTypes.clear();
		cellTypes.push_back(posttype);
		std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

		// Iterate over Presynaptic Types
		for (std::vector<int>::iterator Preit = PreCelltypeList.begin() ; Preit != PreCelltypeList.end(); ++Preit)
		{
			// Presynaptic Cell Type and its Location
			int pretype = *Preit;
			double cellNorm = 0.0;

			if (spatial.compare("rowish")==0)
			{
				std::list< unsigned int > celltypesPre (1,pretype);
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(C1);
				columnsPre.push_back(C2);
				columnsPre.push_back(C3);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}
			else if (spatial.compare("arcish")==0)
			{
				std::list< unsigned int > celltypesPre (1,pretype);
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(D2);
				columnsPre.push_back(C2);
				columnsPre.push_back(B2);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}
			else
			{
				std::list< unsigned int > celltypesPre (1,pretype);
				std::list< unsigned int > columnsPre (1,C2);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}

			somaPosition soma_presyn = celltypeLabels2SomaPos[ColumnCellTypePair(C2,pretype)];

			// Connection Probability Map
			std::map <MatrixIndexType, bin*> binMap;

			// Compute Connection Probability, go through all post-synaptic cells
			for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
			{
				double t, projectedPt[3], soma_postsyn[3];
				soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
				soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
				soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

//				int SGI = BF->laminarPosition(soma_postsyn);
//				if (SGI !=GRAN || !((*rowit)->insideColumn))
//					continue;

				// Compute Position of Post-Synaptic Cell with respect to anatomical landmarks
				double dist1 = 0;
				double Offset = 2000;
				if (spatial.compare("rowish")==0) // Position on Row (from D2)
				{
					dist1 = soma_postsyn[X_COORD] + Offset;
				}
				else if (spatial.compare("arcish")==0) // Position on Arc from (D2)
				{
					dist1 = soma_postsyn[Y_COORD] + Offset;
				}
				else // Compute Distance Post-Synaptic Cell to Vertical Column Axis
				{
					dist1 = vtkLine::DistanceToLine(soma_postsyn, colRef->top, colRef->bottom, t, projectedPt);
					dist1 = sqrt(dist1);
				}
				if (dist1<0)
				{
					std::cout << dist1 << " is negative ! Increase artificial offset " << Offset << std::endl;
					return;
				}
				//dist1 += 25.0;
				int pos1 = (int) dist1/binsz;

				// Compute Distance Post-Synaptic Cell to Average Soma Position of Presynaptic Cell type in C2
				double dist2 = 0.0;
				if (spatial.compare("vertical")==0) // 3D distance in x,y,z
				{
					dist2 = sqrt(pow((soma_presyn.x-soma_postsyn[X_COORD]),2.0)+pow((soma_presyn.y-soma_postsyn[Y_COORD]),2.0)+pow((soma_presyn.z-soma_postsyn[Z_COORD]),2.0));
				}
				else if (spatial.compare("horizontal")==0) // 2D (horizontal) distance in x,y
				{
					dist2 = sqrt(pow((soma_presyn.x-soma_postsyn[X_COORD]),2.0)+pow((soma_presyn.y-soma_postsyn[Y_COORD]),2.0));
				}
				else if ((spatial.compare("arcish")==0) || (spatial.compare("rowish")==0)) // 2D (horizontal) distance in x,y
				{
					dist2 = sqrt(pow((soma_presyn.x-soma_postsyn[X_COORD]),2.0)+pow((soma_presyn.y-soma_postsyn[Y_COORD]),2.0));
				}
				else
				{
					std::cout << "Wrong Input '" << spatial << "' only 'vertical' or 'horizontal' allowed!" << std::endl;
					return;
				}
				int pos2 = (int) dist2/binsz;

				MatrixIndexType idx = MatrixIndexType(pos1,pos2);
				std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

				for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
				{
					bool bolCol = headerit->first.first==C2;
					if (spatial.compare("rowish")==0)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==C1) || (headerit->first.first==C3);
					}
					else if (spatial.compare("arcish")==0)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==B2) || (headerit->first.first==D2);
					}

					// Presynaptic Cells
					if (bolCol && headerit->first.second==pretype)
					{
						std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
						double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
						float cval = (*Convit) * cellWeight;

						// Add to histogram if not in histogram yet
						if (binMap.count(idx) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binMap[idx] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binMap[idx]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
							// Update N (Number of Samples)
							binMap[idx]->N = binMap[idx]->N+1;
							float tmp2 = (cval - binMap[idx]->M);
							// Update SD (can be calculated from Q)
							binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}
					}
				}
			}

			// Output
			std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/" + spatial + "/" + std::string(int2CelltypeLabels[pretype]) + "_" + std::string(int2CelltypeLabels[posttype]) + ".csv";
			storeBinMapAsCSV(outputFilename.c_str(), binMap, binsz);
		}
	}
}

/* Connection Probability over
 * - x: Distance of Postsynaptic Cell to vertical Column axis
 * - y: Intersomatic Distance pre/postsynaptic Cells (Average, not real!)
 * spatial = "horizontal", "vertical", "arcish" or "rowish"
 * - horizontal/vertical: Presynapse: C2 (EX/IN); Postsynapse: 3x3 (CellType);
 * - rowish: Presynapse: C1, C2, C3 (EX/IN); Postsynapse: C1, C2, C3 (CellType);
 * - arcish: Presynapse: B2, C2, D2 (EX/IN); Postsynapse: B2, C2, D2 (CellType);
 * Stores CP for each Cell Type-EX and IN Combination
 * - CPMap for Distance to C2 vs InterSomaDistance
 * - 1D CP plot for Distance to C2 (C2_1D)
 * - 1D CP plot for InterSomaDistance (IS_1D) */
void writeCPMapEXIN(std::string spatial)
{
	if ((spatial.compare("vertical")!=0) && (spatial.compare("horizontal")!=0) && (spatial.compare("arcish")!=0) && (spatial.compare("rowish")!=0))
	{
		std::cout << "Wrong Input '" << spatial << "' only 'vertical', 'horizontal', 'arcish', and 'rowish' allowed!" << std::endl;
		return;
	}

	CellTable * convergence = new CellTable;
	getConvergence(convergence);
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers();

	// Reference Column (C2)
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * colRef = avgColumns.find(C2)->second;

	double binsz = 50;

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	if (spatial.compare("arcish")==0)
	{
		columns.push_back(B2);
		columns.push_back(C2);
		columns.push_back(D2);
	}
	if (spatial.compare("rowish")==0)
	{
		columns.push_back(C1);
		columns.push_back(C2);
		columns.push_back(C3);
	}

	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	// Iterate over Postsynaptic Types
	for (std::vector<int>::iterator Postit = PostCelltypeList.begin() ; Postit != PostCelltypeList.end(); ++Postit)
	{
		int posttype = *Postit;
		cellTypes.clear();
		cellTypes.push_back(posttype);
		std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

		// Iterate over Exc. Types (i=0) and Inh. Types (i=1)
		for (int i=0; i<2; ++i)
		{
			// Connection Probability Map
			std::map <MatrixIndexType, bin*> binMap;
			// Connection Probability 1D
			std::map <int, bin*> bin1; // Relation to BC
			std::map <int, bin*> bin2; // InterSomaDistance

			std::list< unsigned int > celltypesPre;
			double cellNorm = 0.0;

			if (i==0) // Excitatory Presynaptic Cells
			{
				celltypesPre.push_back(L2axon);
				celltypesPre.push_back(L34axon);
				celltypesPre.push_back(L4pyaxon);
				celltypesPre.push_back(L4spaxon);
				celltypesPre.push_back(L4ssaxon);
				celltypesPre.push_back(L5staxon);
				celltypesPre.push_back(L5ttaxon);
				celltypesPre.push_back(L6ccaxon);
				celltypesPre.push_back(L6ccinvaxon);
				celltypesPre.push_back(L6ctaxon);
			}
			else // Inhibitory Presynaptic Cells
			{
				celltypesPre.push_back(SymLocalaxon);
				celltypesPre.push_back(L1axon);
				celltypesPre.push_back(L23Transaxon);
				celltypesPre.push_back(L45Symaxon);
				celltypesPre.push_back(L45Peakaxon);
				celltypesPre.push_back(L56Transaxon);
			}

			if (spatial.compare("rowish")==0)
			{
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(C1);
				columnsPre.push_back(C2);
				columnsPre.push_back(C3);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}
			else if (spatial.compare("arcish")==0)
			{
				std::list< unsigned int > columnsPre;
				columnsPre.push_back(D2);
				columnsPre.push_back(C2);
				columnsPre.push_back(B2);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}
			else // horizontal or vertical
			{
				std::list< unsigned int > columnsPre (1,C2);
				cellNorm = double(getTotalCellNumbers(numCells, columnsPre, celltypesPre));
				cellNorm = cellNorm/double(columnsPre.size()*celltypesPre.size());
			}

			// Compute Connection Probability, go through all post-synaptic cells
			for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
			{
				double t, projectedPt[3], soma_postsyn[3];
				soma_postsyn[X_COORD] = (*rowit)->somaLocation[X_COORD]; // x
				soma_postsyn[Y_COORD] = (*rowit)->somaLocation[Y_COORD]; // y
				soma_postsyn[Z_COORD] = (*rowit)->somaLocation[Z_COORD]; // z

				/* Feature 1:
				 * Compute Position of Post-Synaptic Cell with respect to anatomical landmarks  (C2) */				double dist1 = 0;
				double Offset = 2000;
				if (spatial.compare("rowish")==0) // Position on Row (from D2)
				{
					dist1 = soma_postsyn[X_COORD] + Offset;
				}
				else if (spatial.compare("arcish")==0) // Position on Arc from (D2)
				{
					dist1 = soma_postsyn[Y_COORD] + Offset;
				}
				else // Compute Distance Post-Synaptic Cell to Vertical Column Axis
				{
					dist1 = vtkLine::DistanceToLine(soma_postsyn, colRef->top, colRef->bottom, t, projectedPt);
					dist1 = sqrt(dist1);
				}
				if (dist1<0)
				{
					std::cout << dist1 << " is negative ! Increase artificial offset " << Offset << std::endl;
					return;
				}
				int pos1 = (int) dist1/binsz;
				std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;

				for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
				{
					// Check presynaptic Column
					bool bolCol = headerit->first.first==C2;
					if (spatial.compare("rowish")==0)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==C1) || (headerit->first.first==C3);
					}
					else if (spatial.compare("arcish")==0)
					{
						bolCol = (headerit->first.first==C2) || (headerit->first.first==B2) || (headerit->first.first==D2);
					}

					// Check presynaptic Celltype
					bool bolCell = false;
					for (std::list<unsigned int>::iterator it=celltypesPre.begin(); it != celltypesPre.end(); ++it)
					{
						if (headerit->first.second == (*it))
						{
							bolCell = true;
							break;
						}
					}

					// Presynaptic CellType and Column are correct
					if (bolCol && bolCell)
					{
						/* Feature 2:
						* Compute Distance Post-Synaptic Cell to Average Soma Position of Presynaptic Cell type */
						somaPosition soma_presyn = celltypeLabels2SomaPos[ColumnCellTypePair(C2,headerit->first.second)];
						double dist2 = 0.0;
						if (spatial.compare("vertical")==0) // 3D distance in x,y,z
						{
							dist2 = sqrt(pow((soma_presyn.x-soma_postsyn[X_COORD]),2.0)+pow((soma_presyn.y-soma_postsyn[Y_COORD]),2.0)+pow((soma_presyn.z-soma_postsyn[Z_COORD]),2.0));
						}
						else if (spatial.compare("horizontal")==0) // 2D (horizontal) distance in x,y
						{
							dist2 = sqrt(pow((soma_presyn.x-soma_postsyn[X_COORD]),2.0)+pow((soma_presyn.y-soma_postsyn[Y_COORD]),2.0));
						}
						else if ((spatial.compare("arcish")==0) || (spatial.compare("rowish")==0)) // 2D (horizontal) distance in x,y
						{
							dist2 = sqrt(pow((soma_presyn.x-soma_postsyn[X_COORD]),2.0)+pow((soma_presyn.y-soma_postsyn[Y_COORD]),2.0));
						}
						else
						{
							std::cout << "Wrong Input '" << spatial << "' only 'vertical' or 'horizontal' allowed!" << std::endl;
							return;
						}
						int pos2 = (int) dist2/binsz;
						MatrixIndexType idx = MatrixIndexType(pos1,pos2);

						std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
						double cellWeight = (double(numCells[ColumnCellTypePair(headerit->first.first,headerit->first.second)])) / cellNorm;
						float cval = (*Convit) * cellWeight;

						/* 2D Map */
						// Add to histogram if not in histogram yet
						if (binMap.count(idx) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							binMap[idx] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - binMap[idx]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							binMap[idx]->M = (binMap[idx]->M * binMap[idx]->N + cval) / (binMap[idx]->N+1);
							// Update N (Number of Samples)
							binMap[idx]->N = binMap[idx]->N+1;
							float tmp2 = (cval - binMap[idx]->M);
							// Update SD (can be calculated from Q)
							binMap[idx]->Q = binMap[idx]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}

						/* 1D Map */
						// Feature 1: Relation to BC
						if (bin1.count(pos1) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							bin1[pos1] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - bin1[pos1]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							bin1[pos1]->M = (bin1[pos1]->M * bin1[pos1]->N + cval) / (bin1[pos1]->N+1);
							// Update N (Number of Samples)
							bin1[pos1]->N = bin1[pos1]->N+1;
							float tmp2 = (cval - bin1[pos1]->M);
							// Update SD (can be calculated from Q)
							bin1[pos1]->Q = bin1[pos1]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}

						// Feature 2: InterSomaDistance
						if (bin2.count(pos2) == 0)
						{
							bin * b = new bin;
							b->M = cval;
							b->N = 1;
							b->Q = 0;
							bin2[pos2] = b;
						}
						else // Update histogram
						{
							// Update Mean, SD, and Number of Samples,
							float tmp1 = cval - bin2[pos2]->M;
							// Update Mean of correlation r: (Mean * N + t) / (N+1)
							bin2[pos2]->M = (bin2[pos2]->M * bin2[pos2]->N + cval) / (bin2[pos2]->N+1);
							// Update N (Number of Samples)
							bin2[pos2]->N = bin2[pos2]->N+1;
							float tmp2 = (cval - bin2[pos2]->M);
							// Update SD (can be calculated from Q)
							bin2[pos2]->Q = bin2[pos2]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
						}
					}
				}
			}

			// Output
			std::string outputFilename = "/nas1/Data_daniel/Network/INColumn_v6/ConnectionProbCSV/" + spatial + "/";
			if (i==0) // Excitatory Presynaptic Cells
			{
				outputFilename = outputFilename + "EX_" + std::string(int2CelltypeLabels[posttype]);
			}
			else // Inhibitory Presynaptic Cells
			{
				outputFilename = outputFilename + "IN_" + std::string(int2CelltypeLabels[posttype]);
			}

			std::string outputFilename1 = outputFilename + "_C2_1D.csv";
			std::string outputFilename2 = outputFilename + "_IS_1D.csv";
			outputFilename = outputFilename + ".csv";

			storeBinMapAsCSV(outputFilename.c_str(), binMap, binsz); // 2D
			storeBinAsCSV(outputFilename1.c_str(), bin1, binsz); // C2
			storeBinAsCSV(outputFilename2.c_str(), bin2, binsz); // IS
		}
	}
}

/* Write Mean and SD Soma Position (x,y,z) of each Postsynaptic CellType/Column Combination (3x3)  */
void writeSomaPosition()
{
	initializeConstants();
	initializeLists();

	std::string inputFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/convergenceTotal.csv";

	CellTable * convergence = new CellTable;
	getConvergence(convergence);

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	std::string Outputfname = "/home/dudvary/Documents/SomaPos.csv";
	std::ofstream CSVwriter;
	CSVwriter.open(Outputfname.c_str());

	if(!CSVwriter.fail())
	{
		CSVwriter << "Location of Soma [" << inputFilename << "]" << std::endl;
		CSVwriter << "Column,CellType,M_Soma_x,M_Soma_y,M_Soma_z,SD_Soma_x,SD_Soma_y,SD_Soma_z,N," << std::endl;

		for (std::vector<int>::iterator Colit = PostColumnList.begin() ; Colit != PostColumnList.end(); ++Colit)
		{
			 columns.clear();
			 columns.push_back(*Colit);

			 for (std::vector<int>::iterator Typeit = PostCelltypeList.begin() ; Typeit != PostCelltypeList.end(); ++Typeit)
			 {
				 cellTypes.clear();
				 cellTypes.push_back(*Typeit);
				 std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

				 std::vector< bin *> binList;
				 for (int i = 0; i < 3; i++)
				 {
						bin * b = new bin;
						b->M = 0;
						b->N = 0;
						b->Q = 0;
						binList.push_back(b);
				 }

				// Compute Connection Probability for one combination of celltype/column pair
				for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
				{
					// Extract Soma Location Value of Postsynaptic Col/CellType Pair
					double somaLocation[3];
					somaLocation[0] = (*rowit)->somaLocation[0]; // x
					somaLocation[1] = (*rowit)->somaLocation[1]; // y
					somaLocation[2] = (*rowit)->somaLocation[2]; // z

					int ii = 0;
					for (std::vector< bin *>::iterator Binit = binList.begin(); Binit != binList.end(); ++Binit, ++ii)
					{
						// Update Mean, SD, and Number of Samples,
						float tmp1 = somaLocation[ii] - (*Binit)->M;
						// Update Mean of correlation r: (Mean * N + t) / (N+1)
						(*Binit)->M = ((*Binit)->M * (*Binit)->N + somaLocation[ii]) / ((*Binit)->N+1);
						// Update N (Number of Samples)
						(*Binit)->N = (*Binit)->N+1;
						float tmp2 = (somaLocation[ii] - (*Binit)->M);
						// Update SD (can be calculated from Q)
						(*Binit)->Q = (*Binit)->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
					}
				}

				CSVwriter << int2ColumnLabels[*Colit] << "," << int2CelltypeLabels[*Typeit] << ",";

				if ((binList.at(0)->N != binList.at(1)->N) || (binList.at(1)->N != binList.at(2)->N))
				{
					std::cout << "ERROR! N values are different, something is wrong! " << binList.at(0)->N << " " << binList.at(1)->N << " " << binList.at(2)->N << std::endl;
					//return 0;
				}

				CSVwriter << binList.at(0)->M << "," << binList.at(1)->M << "," << binList.at(2)->M;

				for (int i=0; i != 3; ++i)
				{
					int N = binList.at(i)->N;
					float Q = binList.at(i)->Q;
					float SD = sqrt(Q/(N-1));

					CSVwriter << "," << SD;
				}
				CSVwriter << "," << binList.at(0)->N << "," << std::endl;
			 }
		}
	}
	else
	{
		std::cout << "ERROR! Writing " << Outputfname << " failed! " << std::endl;
	}
	CSVwriter.close();
}

/* Compartment: Apical, Basal or Total ('')
 * Input filename
 * Output filename
 * numCellsfname : path to number of cells.csv (necessary to merge locals)
 * LocalSubtype: 0: Merge Local Subtypes 1: Local Subtypes
 * AV: 0: compute SD, 1: compute Mean/Average
 * */
void writeProbabilityMeanTable(const char * inputfname, const char * outputfname, const char * numCellsfname, bool LocalSubType, int AV)
{
	initializeConstants();
	initializeLists();

	CellTable * convergence = new CellTable;
	getConvergence(convergence,inputfname,numCellsfname,LocalSubType);

	// Extract CellType/Column Combination (postsynaptic)
	std::list< unsigned int > columns;
	std::list< unsigned int > cellTypes;
	std::map< ColumnCellTypePair, unsigned int > header = convergence->header;

	if (LocalSubType)
	{
		PostCelltypeList.clear();
		PostCelltypeList.push_back(L2);
		PostCelltypeList.push_back(L34);
		PostCelltypeList.push_back(L4py);
		PostCelltypeList.push_back(L4sp);
		PostCelltypeList.push_back(L4ss);
		PostCelltypeList.push_back(L5st);
		PostCelltypeList.push_back(L5tt);
		PostCelltypeList.push_back(L6cc);
		PostCelltypeList.push_back(L6ccinv);
		PostCelltypeList.push_back(L6ct);
		PostCelltypeList.push_back(SymLocal1);
		PostCelltypeList.push_back(SymLocal2);
		PostCelltypeList.push_back(SymLocal3);
		PostCelltypeList.push_back(SymLocal4);
		PostCelltypeList.push_back(SymLocal5);
		PostCelltypeList.push_back(SymLocal6);
		PostCelltypeList.push_back(L1);
		PostCelltypeList.push_back(L23Trans);
		PostCelltypeList.push_back(L45Sym);
		PostCelltypeList.push_back(L45Peak);
		PostCelltypeList.push_back(L56Trans);
	}
	else
	{
		PostCelltypeList.clear();
		PostCelltypeList.push_back(L2);
		PostCelltypeList.push_back(L34);
		PostCelltypeList.push_back(L4py);
		PostCelltypeList.push_back(L4sp);
		PostCelltypeList.push_back(L4ss);
		PostCelltypeList.push_back(L5st);
		PostCelltypeList.push_back(L5tt);
		PostCelltypeList.push_back(L6cc);
		PostCelltypeList.push_back(L6ccinv);
		PostCelltypeList.push_back(L6ct);
		PostCelltypeList.push_back(SymLocal);
		PostCelltypeList.push_back(L1);
		PostCelltypeList.push_back(L23Trans);
		PostCelltypeList.push_back(L45Sym);
		PostCelltypeList.push_back(L45Peak);
		PostCelltypeList.push_back(L56Trans);
	}

	std::ofstream CSVwriter;
	CSVwriter.open(outputfname);

	if(!CSVwriter.fail())
	{
		CSVwriter << "CONNECTION PROBABILITY AVERAGE [" << inputfname << "]" << std::endl;
		// Write in csv file
		CSVwriter <<  "COLUMN-CELLTYPE,";
		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
		{
			CSVwriter << int2ColumnLabels[headerit->first.first] << "-" << int2CelltypeLabels[headerit->first.second] << ",";
		}
		CSVwriter << std::endl;

		for (std::vector<int>::iterator Colit = PostColumnList.begin() ; Colit != PostColumnList.end(); ++Colit)
		{
			 columns.clear();
			 columns.push_back(*Colit);

			 for (std::vector<int>::iterator Typeit = PostCelltypeList.begin() ; Typeit != PostCelltypeList.end(); ++Typeit)
			 {
				 cellTypes.clear();
				 cellTypes.push_back(*Typeit);
				 std::vector< CellTableRow * > rows = convergence->getPostColumnCelltypeRows(columns, cellTypes);

				 std::vector< bin *> binList;
				 for (int i = 0; i < header.size(); i++)
				 {
						bin * b = new bin;
						b->M = 0;
						b->N = 0;
						b->Q = 0;
						binList.push_back(b);
				 }

				// Compute Connection Probability for one combination of celltype/column pair
				for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
				{
					// Extract Convergence Value of Presynaptic Col/CellType Pair
					std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;
					std::vector< bin *>::iterator Binit = binList.begin();
					for (std::vector<float>::iterator Convit = convergenceVal.begin(); Convit != convergenceVal.end(); ++Convit, ++Binit)
					{
						// Update Mean, SD, and Number of Samples,
						float tmp1 = (*Convit) - (*Binit)->M;
						// Update Mean of correlation r: (Mean * N + t) / (N+1)
						(*Binit)->M = ((*Binit)->M * (*Binit)->N + (*Convit)) / ((*Binit)->N+1);
						// Update N (Number of Samples)
						(*Binit)->N = (*Binit)->N+1;
						float tmp2 = ((*Convit) - (*Binit)->M);
						// Update SD (can be calculated from Q)
						(*Binit)->Q = (*Binit)->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
					}
				}

				CSVwriter << int2ColumnLabels[*Colit] << "-" << int2CelltypeLabels[*Typeit];
				for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
				{
					if (AV==0)
					{	// SD
						int N = binList.at(headerit->second)->N;
						float Q = binList.at(headerit->second)->Q;
						CSVwriter << "," << sqrt(Q/(N-1));
					}
					else
					{	// AVERAGE
						CSVwriter << "," << binList.at(headerit->second)->M;
					}

				}
				CSVwriter << "," << std::endl;
			 }
		}
	}
	else
	{
		std::cout << "ERROR! Writing " << outputfname << " failed! " << std::endl;
	}
	CSVwriter.close();

	// Reset PostCelltypeList
	if (LocalSubType)
	{
		PostCelltypeList.clear();
		PostCelltypeList.push_back(L2);
		PostCelltypeList.push_back(L34);
		PostCelltypeList.push_back(L4py);
		PostCelltypeList.push_back(L4sp);
		PostCelltypeList.push_back(L4ss);
		PostCelltypeList.push_back(L5st);
		PostCelltypeList.push_back(L5tt);
		PostCelltypeList.push_back(L6cc);
		PostCelltypeList.push_back(L6ccinv);
		PostCelltypeList.push_back(L6ct);
		PostCelltypeList.push_back(SymLocal);
		PostCelltypeList.push_back(L1);
		PostCelltypeList.push_back(L23Trans);
		PostCelltypeList.push_back(L45Sym);
		PostCelltypeList.push_back(L45Peak);
		PostCelltypeList.push_back(L56Trans);
	}
}

void writeProbabilityMeanTable(std::string compartment, bool LocalSubType)
{
	int AV = 1;
	const char * numCellsfname = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/nrCells.csv";
	std::string inputfname = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/convergence" + compartment + ".csv";
	std::string Outputfname = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/ConnectionProbability" + compartment + "_";

	if (LocalSubType)
		Outputfname = Outputfname + "LocalSubTypes_";
	if (AV==0)
		Outputfname = Outputfname + "SD.csv";
	else
		Outputfname = Outputfname + "AV.csv";

	writeProbabilityMeanTable(inputfname.c_str(), Outputfname.c_str(), numCellsfname, LocalSubType, AV);
}

/* Compute Correlation between Cell Types */
void outputCorrBins(const char * inputfname)
{
	if (celltypeLabels2Int.count(std::string(inputfname)) > 0)
	{
		std::cout << "ERROR! Wrong Cell Type! " << inputfname << " does not exist!" << std::endl;
		//return NULL;
	}

	CellTable * connectome = new CellTable;
	getConnectome(connectome);

	// Create Map with key (corresponding to bin position, with a pair: Mean,Number of samples
	int binsz = 50;

	std::list< unsigned int > columns;
	//columns.push_back(C2);
	std::list< unsigned int > cellTypes;
	cellTypes.push_back(celltypeLabels2Int[std::string(inputfname)]);

	std::vector< CellTableRow * > row1 = connectome->getPostColumnCelltypeRows(columns, cellTypes);

	// Go through all Cell types
	#pragma omp parallel for
	for (int i = 0; i < PostCelltypeList.size(); ++i)
	{
		// Initialize with 0s
		std::vector< bin * > valVector;
		for (int c=0; c!= 100; ++c)
		{
			bin * b = new bin;
			b->M = 0;
			b->N = 0;
			b->Q = 0;
			valVector.push_back(b);
		}

		// Get two rows
		std::list< unsigned int > cellTypes2;
		cellTypes2.push_back(PostCelltypeList.at(i));
		std::vector< CellTableRow * > row2 = connectome->getPostColumnCelltypeRows(columns, cellTypes2);
		getCorrPerBin(valVector, row1, row2, binsz);

		// Output
		std::string s1 = int2CelltypeLabels[cellTypes.front()];
		std::string s2 = int2CelltypeLabels[PostCelltypeList.at(i)];

		std::string Outputfname = "/home/dudvary/Documents/Corr/" + s1 + "/MeanCorrValues_" + s1 + "-" + s2 + ".csv";
		std::ofstream CSVwriter;
		CSVwriter.open(Outputfname.c_str());

		if(!CSVwriter.fail())
		{
			CSVwriter << "bin,MeanR,N,SD" << std::endl;
			// Write in csv file
			int b = 0;
			for (std::vector< bin * >::iterator it = valVector.begin(); it!=valVector.end(); ++it, ++b)
			{
				float Q = (*it)->Q;
				int N = (*it)->N;
				CSVwriter << b*binsz+binsz/2 << "," << (*it)->M << "," << N << "," << sqrt(Q/(N-1)) << std::endl;
			}
		}
		else
		{
			std::cout << "ERROR! Writing " << Outputfname << " failed! " << std::endl;
		}
		CSVwriter.close();

		std::cout << s1 << "-" << s2 << std::endl;
	}

	// Output file
//	const char * Outputfname = "/home/dudvary/Documents/CorrValues.csv";
//	std::ofstream CSVwriter;
//	CSVwriter.open(Outputfname);
//
//	if(!CSVwriter.fail())
//	{
//		CSVwriter <<  "r,d,celltype,column" << std::endl;
//
//		for (std::vector< CellTableRow * >::iterator row1it = row1.begin(); row1it != row1.end(); ++row1it)
//		{
//			std::vector< float > somaLocation1 = (*row1it)->somaLocation;
//			std::list<double> synapsesPerPreTypeColumn1((*row1it)->synapsesPerPreTypeColumn.begin(),(*row1it)->synapsesPerPreTypeColumn.end());
//			cellIDlist.insert((*row1it)->cellID);
//
//			for (std::vector< CellTableRow * >::iterator row2it = row2.begin(); row2it != row2.end(); ++row2it)
//			{
//				if (cellIDlist.find((*row2it)->cellID) != cellIDlist.end())
//				{
//					continue;
//				}
//
//				std::vector< float > somaLocation2 = (*row2it)->somaLocation;
//				std::list<double> synapsesPerPreTypeColumn2((*row2it)->synapsesPerPreTypeColumn.begin(),(*row2it)->synapsesPerPreTypeColumn.end());
//
//				// Compute Correlation Coefficient r between Input Vector of Postsynaptic Cells
//				//r.push_back(helper::computeCorr(synapsesPerPreTypeColumn1, synapsesPerPreTypeColumn2));
//				// Compute Distance between Postsynaptic Cells
//				//dist.push_back(sqrt(pow((somaLocation2[X_COORD]-somaLocation2[X_COORD]),2.0)+pow((somaLocation2[Y_COORD]-somaLocation2[Y_COORD]),2.0)+pow((somaLocation2[Z_COORD]-somaLocation2[Z_COORD]),2.0)));
//
//				float rtmp = helper::computeCorr(synapsesPerPreTypeColumn1, synapsesPerPreTypeColumn2);
//				float dtmp = sqrt(pow((somaLocation2[X_COORD]-somaLocation1[X_COORD]),2.0)+pow((somaLocation2[Y_COORD]-somaLocation1[Y_COORD]),2.0)+pow((somaLocation2[Z_COORD]-somaLocation1[Z_COORD]),2.0));
//
//				int binpos = (int) dtmp/binsz;
//
//				//std::pair< float, int> binval (M,N);
//
//				// Write in csv file
//				CSVwriter << rtmp << "," << dtmp << "," << (*row2it)->cellType << "," << (*row2it)->column << std::endl;
//			}
//			//std::cout << i1/row1.size()*100 << "%" << std::endl;
//		}
//	}
//	else
//	{
//		std::cout << "ERROR! Writing " << Outputfname << " failed! " << std::endl;
//	}
//	CSVwriter.close();
}

/* Compute Correlation between Inputvectors of Combiniation of all cells of one Cell Type
 * Sort Correlation according to euclidean distance (x,y,z) between cells
 * inputfname = name of celltype (e.g. L34, SymLocal)
 * synapse = use Connectome (Number of Synapses; true) or Convergences (false) */
void outputCorrBinsItself(const char * inputfname, bool synapse)
{
	if (celltypeLabels2Int.count(std::string(inputfname)) > 0)
	{
		std::cout << "ERROR! Wrong Cell Type! " << inputfname << " does not exist!" << std::endl;
		return;
	}

	CellTable * connectome = new CellTable;

	if (synapse)
	{
		getConnectome(connectome);
	}
	else
	{
		getConvergence(connectome);
	}

	// Create Map with key (corresponding to bin position, with a pair: Mean,Number of samples
	int binsz = 50;

	std::list< unsigned int > columns;
	//columns.push_back(C2);
	std::list< unsigned int > cellTypes;
	cellTypes.push_back(celltypeLabels2Int[std::string(inputfname)]);

	std::vector< CellTableRow * > row1 = connectome->getPostColumnCelltypeRows(columns, cellTypes);
	// Initialize with 0s
	std::vector< bin * > valVector;
	for (int c=0; c!= 100; ++c)
	{
		bin * b = new bin;
		b->M = 0;
		b->N = 0;
		b->Q = 0;
		valVector.push_back(b);
	}

	// Get two rows
	std::list< unsigned int > cellTypes2;
	cellTypes2.push_back(celltypeLabels2Int[std::string(inputfname)]);
	std::vector< CellTableRow * > row2 = connectome->getPostColumnCelltypeRows(columns, cellTypes2);
	getCorrPerBin(valVector, row1, row2, binsz);

	// Output
	std::string s1 = int2CelltypeLabels[cellTypes.front()];
	std::string s2 = int2CelltypeLabels[cellTypes2.front()];

	std::string Outputfname = "/home/dudvary/Documents/Corr/" + s1 + "/MeanCorrValues_" + s1 + "-" + s2;
	if (synapse)
	{
		Outputfname = Outputfname + ".csv";
	}
	else
	{
		Outputfname = Outputfname + "_CP.csv";
	}

	std::ofstream CSVwriter;
	CSVwriter.open(Outputfname.c_str());

	if(!CSVwriter.fail())
	{
		CSVwriter << "bin,MeanR,N,SD" << std::endl;
		// Write in csv file
		int b = 0;
		for (std::vector< bin * >::iterator it = valVector.begin(); it!=valVector.end(); ++it, ++b)
		{
			float Q = (*it)->Q;
			int N = (*it)->N;
			CSVwriter << b*binsz+binsz/2 << "," << (*it)->M << "," << N << "," << sqrt(Q/(N-1)) << std::endl;
		}
	}
	else
	{
		std::cout << "ERROR! Writing " << Outputfname << " failed! " << std::endl;
	}
	CSVwriter.close();

	std::cout << s1 << "-" << s2 << std::endl;
}

/* read in Convergence Table csv and then (for verification) store it in new .csv file */
void writeConvergenceTableCSV()
{
	CellTable * convergence = new CellTable;
	getConvergence(convergence);

	const char * Outputfname = "/home/dudvary/Documents/test.csv";
	std::ofstream CSVwriter;
	CSVwriter.open(Outputfname);
	if(!CSVwriter.fail())
	{
		// Header
		CSVwriter <<  "CELLID,CELLTYPE,COLUMN,SOMA_X,SOMA_Y,SOMA_Z,INSIDE_COLUMN," ;
		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=convergence->header.begin(); headerit!=convergence->header.end(); ++headerit)
		{
			CSVwriter << int2ColumnLabels[headerit->first.first] << "-" << int2CelltypeLabels[headerit->first.second] << ",";
		}
		CSVwriter << std::endl;

		std::vector< CellTableRow * > rows = convergence->rows;

		for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
		{
			CSVwriter << (*rowit)->cellID << "," << int2CelltypeLabels[(*rowit)->cellType] << "," << int2ColumnLabels[(*rowit)->column] << ",";
			CSVwriter << (*rowit)->somaLocation[0] << "," << (*rowit)->somaLocation[1] << "," << (*rowit)->somaLocation[2] << "," << (*rowit)->insideColumn << ",";

			std::vector< float > synapsesPerPreTypeColumn = (*rowit)->synapsesPerPreTypeColumn;

			for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=convergence->header.begin(); headerit!=convergence->header.end(); ++headerit)
			{
				CSVwriter << synapsesPerPreTypeColumn.at(headerit->second) << ",";
			}

			CSVwriter << std::endl;
		}
	}
	else
	{
		std::cout << "ERROR! Writing " << Outputfname << " failed! " << std::endl;
	}
	CSVwriter.close();
}

/* Compute Mean and SD of Correlation for each bin of row1 and row2 combination */
void getCorrPerBin(std::vector< bin * > valVector, std::vector< CellTableRow * > row1, std::vector< CellTableRow * > row2, int binsz)
{
	std::set<int> cellIDlist;

	for (std::vector< CellTableRow * >::iterator row1it = row1.begin(); row1it != row1.end(); ++row1it)
	{
		std::vector< float > somaLocation1 = (*row1it)->somaLocation;
		std::list<double> synapsesPerPreTypeColumn1((*row1it)->synapsesPerPreTypeColumn.begin(),(*row1it)->synapsesPerPreTypeColumn.end());
		cellIDlist.insert((*row1it)->cellID);

		for (std::vector< CellTableRow * >::iterator row2it = row2.begin(); row2it != row2.end(); ++row2it)
		{
			if (cellIDlist.find((*row2it)->cellID) != cellIDlist.end())
			{
				continue;
			}

			std::vector< float > somaLocation2 = (*row2it)->somaLocation;
			std::list<double> synapsesPerPreTypeColumn2((*row2it)->synapsesPerPreTypeColumn.begin(),(*row2it)->synapsesPerPreTypeColumn.end());

			float rtmp = helper::computeCorr(synapsesPerPreTypeColumn1, synapsesPerPreTypeColumn2);
			float dtmp = sqrt(pow((somaLocation2[X_COORD]-somaLocation1[X_COORD]),2.0)+pow((somaLocation2[Y_COORD]-somaLocation1[Y_COORD]),2.0)+pow((somaLocation2[Z_COORD]-somaLocation1[Z_COORD]),2.0));
			int pos = (int) dtmp/binsz;

			// Store values in map
			if (valVector.size() < pos)
			{
				flush(std::cout << "Map size too small!" << std::endl);
				continue;
			}
			if (0 > pos)
			{
				flush(std::cout << "Something is wrong!" << std::endl);
				continue;
			}

			// Update Mean, SD, and Number of Samples,
			float tmp1 = rtmp - valVector.at(pos)->M;
			// Update Mean of correlation r: (Mean * N + t) / (N+1)
			valVector.at(pos)->M = (valVector.at(pos)->M * valVector.at(pos)->N + rtmp) / (valVector.at(pos)->N+1);
			// Update N (Number of Samples)
			valVector.at(pos)->N = valVector.at(pos)->N+1;
			float tmp2 = (rtmp - valVector.at(pos)->M);
			// Update SD (can be calculated from Q)
			valVector.at(pos)->Q = valVector.at(pos)->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
		}
	}
}

/* Function for computing Average Synapse Table per Postsynaptic Cell */
void computeSynapseMeanTable()
{
	CellTable * connectome = new CellTable;
	getConnectome(connectome);

	/* Compute Average Number of Synapses per Postsynaptic Cell Type */
	std::map< ColumnCellTypePair, unsigned int > header = connectome->header;
	std::list< unsigned int > columns;
	std::list< unsigned int > cellTypes;
	std::vector< float > synapsesMean;
	int sz = 0;

	// Go through each Postsynaptic Column
	for (int k = 0; k!=PostColumnList.size(); ++k)
	{
		columns.clear();
		columns.push_front(PostColumnList.at(k));

		// Go through each Postsynaptic Cell Type and compute its average number of synapses
		for (int i = 0; i!=PostCelltypeList.size(); ++i)
		{
			cellTypes.clear();
			cellTypes.push_front(PostCelltypeList.at(i));

			std::vector< CellTableRow * > selection = connectome->getPostColumnCelltypeRows(columns, cellTypes);

			if (selection.empty())
			{
				std::cout << "WARNING! " <<  int2CelltypeLabels[PostCelltypeList.at(i)] << " not in Table!" << std::endl;
				continue;
			}

			std::vector< float > synapsestmp = computeSynMean(selection);
			for (int j=0; j!=synapsestmp.size(); j++)
			{
				synapsesMean.push_back(synapsestmp.at(j));
			}

			sz = synapsestmp.size();
		}
	}

	/* Store in .csv Output */
	std::string Outputfname = "/home/dudvary/Documents/SynapsesPerCellTable_cpp_Mean.csv";
	writeSynapseMeanTable(Outputfname.c_str(), synapsesMean, header, sz);

	delete connectome;
}

/* Store in .csv Output
 * Outputfname: Output file name for .csv file
 * synapsesMean: vector array with mean synapses
 * header: for labels of header
 * sz: length of one row */
void writeSynapseMeanTable(const char * Outputfname, std::vector< float > synapsesMean, std::map< ColumnCellTypePair, unsigned int > header, int sz)
{
	std::ofstream CSVwriter;
	CSVwriter.open(Outputfname);
	if(!CSVwriter.fail())
	{
		// Header
		CSVwriter <<  "PreColumn,PreCellType" ;
		for (int k = 0; k!=PostColumnList.size(); ++k)
		{
			for (int i = 0; i!=PostCelltypeList.size(); ++i)
			{
				CSVwriter << "," << int2ColumnLabels[PostColumnList.at(k)] << "-" << int2CelltypeLabels[PostCelltypeList.at(i)];
			}
		}
		CSVwriter << std::endl;

		// Rows with Average Synapses
		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
		{
			CSVwriter << int2ColumnLabels[headerit->first.first] << "," << int2CelltypeLabels[headerit->first.second];

			for (int k = 0; k!=PostColumnList.size(); ++k)
			{
				for (int i = 0; i!=PostCelltypeList.size(); ++i)
				{
					int pos = (i+PostCelltypeList.size()*k)*sz;
					CSVwriter << "," << synapsesMean.at(headerit->second-8+pos);
				}
			}
			CSVwriter << std::endl;
		}
	}
	else
	{
		std::cout << "ERROR! Writing " << Outputfname << " failed! " << std::endl;
	}
	CSVwriter.close();
}

/* Compute Mean over Synapses for given selection
 * returns vector containing mean synapses per postsynaptic cell type*/
std::vector< float > computeSynMean(std::vector< CellTableRow * > selection)
{
	std::vector< CellTableRow * >::iterator row = selection.begin();
	std::vector< float > synapsesMean (((*row)->synapsesPerPreTypeColumn).size(),0.0);

	// Go through rows (individual cells)
	for (row=selection.begin(); row != selection.end(); ++row)
	{
		std::vector< float > synapsesPerPreTypeColumn = (*row)->synapsesPerPreTypeColumn;

		if (synapsesMean.size() != synapsesPerPreTypeColumn.size())
		{
			std::cout << "ERROR! Not same size!" << std::endl;
			continue;
		}

		// Go through synapses of each presynaptic column/celltype
		for (int i=0; i<synapsesMean.size(); i++)
		{
			synapsesMean.at(i) += synapsesPerPreTypeColumn.at(i);
		}
	}

	for (int i=0; i<synapsesMean.size(); i++)
	{
		synapsesMean.at(i) /= selection.size();
	}

	return synapsesMean;
}

/* Get Number of Cells in each Column
 * Returns Map containing ColumnCellTypePair and corresponding number of cells in nearest column
 * Not same function as for cell table (cell table only returns number of local types!) */
std::map< ColumnCellTypePair, unsigned int > getCellNumbers()
{
	const char * inputFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/data/nrCells.csv";
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers(inputFilename);
	return numCells;
}

/* Get Number of Cells in each Column
 * Returns Map containing ColumnCellTypePair and corresponding number of cells in nearest column
 * Not same function as for cell table (cell table only returns number of local types!) */
std::map< ColumnCellTypePair, unsigned int > getCellNumbers(const char* fname)
{
	std::ifstream inputStream(fname);
	std::map< ColumnCellTypePair, unsigned int > numCells;

	if(!inputStream.fail())
	{
		if(celltypeLabels2Int.size()==0 || ColumnLabels2Int.size()==0)
		{
			initializeConstants();
		}

		std::string currentLine;
		std::vector< int > colVector;

		bool headerFound = false;
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				if (currentLine.find("CELL TYPE/COLUMN TOTAL") != std::string::npos)
				{
					std::getline(inputStream, currentLine);

					std::size_t delim = currentLine.find("\t");

					while(delim != std::string::npos)
					{
						currentLine = currentLine.substr(delim+1);
						delim = currentLine.find("\t");
						colVector.push_back(ColumnLabels2Int[currentLine.substr(0,delim)]);
					}
					std::getline(inputStream, currentLine);
					headerFound = true;
				}

				if (currentLine.find("CELL TYPE/COLUMN INSIDE COLUMN") != std::string::npos)
				{
					break;
				}

				if (headerFound)
				{
					std::size_t delim = currentLine.find("\t");

					if (currentLine.find("TOTAL") != std::string::npos)
					{
						break;
					}

					if (celltypeLabels2Int.count(currentLine.substr(0,delim))==0)
					{
						std::cout << "Unknown Celltype: " << currentLine.substr(0,delim) << ";" << std::endl;
						break;
					}

					int celltypeID = celltypeLabels2Int[currentLine.substr(0,delim)];
					int i = 0;

					while(delim != std::string::npos)
					{
						if (i>=colVector.size())
							break;

						currentLine = currentLine.substr(delim+1);
						delim = currentLine.find("\t");
						std::string tmp = currentLine.substr(0,delim);
						unsigned int n = atoi(tmp.c_str());
						numCells.insert(std::pair< ColumnCellTypePair, unsigned int >(ColumnCellTypePair(colVector.at(i), celltypeID), n));
						i++;
					}

					if ((i)!=colVector.size())
					{
						flush(std::cout << "Header size (" << i << ") and Column Vector size (" << colVector.size() << ") do not match!" << std::endl);
						break;
					}
				}

			}
		}

		inputStream.close();

		/* Add SymLocalaxon and SymLocal (sum of SymLocals) */
		for (std::vector<int>::iterator it = colVector.begin() ; it != colVector.end(); ++it)
		{
			unsigned int numDend = numCells[(ColumnCellTypePair((*it), SymLocal1))] + numCells[(ColumnCellTypePair((*it), SymLocal2))] +
					numCells[(ColumnCellTypePair((*it), SymLocal3))] + numCells[(ColumnCellTypePair((*it), SymLocal4))] +
					numCells[(ColumnCellTypePair((*it), SymLocal5))] + numCells[(ColumnCellTypePair((*it), SymLocal6))];

			unsigned int numAxon = numCells[(ColumnCellTypePair((*it), SymLocal1axon))] + numCells[(ColumnCellTypePair((*it), SymLocal2axon))] +
					numCells[(ColumnCellTypePair((*it), SymLocal3axon))] + numCells[(ColumnCellTypePair((*it), SymLocal4axon))] +
					numCells[(ColumnCellTypePair((*it), SymLocal5axon))] + numCells[(ColumnCellTypePair((*it), SymLocal6axon))];

			numCells.insert(std::pair< ColumnCellTypePair, unsigned int >(ColumnCellTypePair((*it), SymLocal), numDend));
			numCells.insert(std::pair< ColumnCellTypePair, unsigned int >(ColumnCellTypePair((*it), SymLocalaxon), numAxon));
		}
	}
	else
	{
		std::cout << "ERROR! Reading nrCells table " << fname << " failed!" << std::endl;
	}

	return numCells;
}

int getTotalCellNumbers(std::map< ColumnCellTypePair, unsigned int > numCells, std::list< unsigned int > columnsPre, std::list< unsigned int > cellTypesPre)
{
	int totalCells = 0;
	for (std::list<unsigned int>::iterator ceit = cellTypesPre.begin(); ceit != cellTypesPre.end(); ++ceit)
	{
		for (std::list<unsigned int>::iterator coit = columnsPre.begin(); coit != columnsPre.end(); ++coit)
		{
			//std::cout << numCells[ColumnCellTypePair(*coit,*ceit)] << ";" << std::endl;
			totalCells = totalCells + numCells[ColumnCellTypePair(*coit,*ceit)];
		}
	}
	return totalCells;
}

/* ------------------------------------
 * Additional Initialization Functions
 * ------------------------------------
 */

/* Soma Positions of Presynaptic Cells located in C2 Column
 * to allow computation of Inter-Soma-Distance between Pre- and Postsynaptic cell
 */
void initalizeSomaPositions()
{
	if(celltypeLabels2SomaPos.size())
		celltypeLabels2SomaPos.clear();

	/* These are number of dendrites not axons! */
	/* C2 */
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L2axon),somaPosition(-72,442,411)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L34axon),somaPosition(-80,414,139)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L4pyaxon),somaPosition(-90,410,49)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L4spaxon),somaPosition(-91,408,24)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L4ssaxon),somaPosition(-89, 408 ,-38)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L5staxon),somaPosition(-98,375,-342)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L5ttaxon),somaPosition(-103,371,-435)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L6ccaxon),somaPosition(-109,350,-642)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L6ccinvaxon),somaPosition(-112,338,-779)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L6ctaxon),somaPosition(-120,329,-922)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,SymLocalaxon),somaPosition(-92,398,-124)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L1axon),somaPosition(-51,477,608)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L23Transaxon),somaPosition(-76,436,342)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L45Symaxon),somaPosition(-87,393,-196)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L45Peakaxon),somaPosition(-79,403,-83)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C2,L56Transaxon),somaPosition(-110,346,-620)));

	/* Remaining Columns in 3x3 Grid */
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L1axon),somaPosition(-474.67,832.124,598.593)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L2axon),somaPosition(-481.868,795.59,414.612)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L23Transaxon),somaPosition(-485.262,794.195,347.445)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L34axon),somaPosition(-492.115,765.562,167.822)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L45Peakaxon),somaPosition(-483.291,752.277,-54.9236)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L45Symaxon),somaPosition(-481.912,737.636,-153.89)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L4pyaxon),somaPosition(-481.693,761.247,69.1608)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L4spaxon),somaPosition(-502.015,746.315,55.7415)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L4ssaxon),somaPosition(-506.437,742.378,-13.505)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L56Transaxon),somaPosition(-509.78,688.374,-551.504)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L5staxon),somaPosition(-498.334,712.935,-290.351)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L5ttaxon),somaPosition(-499.997,710.386,-370.537)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L6ccaxon),somaPosition(-504.777,690.296,-552.633)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L6ccinvaxon),somaPosition(-508.315,676.074,-666.661)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,L6ctaxon),somaPosition(-513.216,665.847,-800.329)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B1,SymLocalaxon),somaPosition(-498.038,747.163,-64.7091)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L1axon),somaPosition(-90.7204,923.3,562.975)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L2axon),somaPosition(-113.169,867.188,392.291)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L23Transaxon),somaPosition(-123.96,849.92,297.57)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L34axon),somaPosition(-145.46,810.196,122.897)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L45Peakaxon),somaPosition(-158.694,775.253,-72.6116)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L45Symaxon),somaPosition(-166.898,748.801,-172.977)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L4pyaxon),somaPosition(-149.176,797.318,33.1878)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L4spaxon),somaPosition(-146.378,786.993,15.7714)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L4ssaxon),somaPosition(-152.059,775.094,-48.6073)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L56Transaxon),somaPosition(-207.072,676.432,-604.931)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L5staxon),somaPosition(-174.518,724.365,-341.662)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L5ttaxon),somaPosition(-187.596,709.643,-424.476)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L6ccaxon),somaPosition(-204.979,674.169,-621.368)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L6ccinvaxon),somaPosition(-215.187,653.349,-745.449)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,L6ctaxon),somaPosition(-223.542,628.505,-874.579)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B2,SymLocalaxon),somaPosition(-156.631,779.983,-55.833)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L1axon),somaPosition(233.521,896.117,524.326)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L2axon),somaPosition(233.835,863.023,340.995)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L23Transaxon),somaPosition(214.253,840.281,257.261)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L34axon),somaPosition(191.148,787.759,73.2667)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L45Peakaxon),somaPosition(177.068,756.191,-132.742)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L45Symaxon),somaPosition(168.636,733.377,-226.564)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L4pyaxon),somaPosition(178.723,767.754,-13.9058)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L4spaxon),somaPosition(193.469,777.172,-39.1041)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L4ssaxon),somaPosition(187.033,765.6,-93.7587)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L56Transaxon),somaPosition(95.5421,621.811,-689.83)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L5staxon),somaPosition(133.07,687.681,-396.992)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L5ttaxon),somaPosition(127.383,658.612,-485.81)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L6ccaxon),somaPosition(97.7243,620.074,-690.878)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L6ccinvaxon),somaPosition(78.4572,586.265,-829.777)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,L6ctaxon),somaPosition(61.8776,558.918,-958.072)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(B3,SymLocalaxon),somaPosition(166.004,752.589,-129.65)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L1axon),somaPosition(-497.975,391.558,635.023)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L2axon),somaPosition(-514.848,415.784,428.54)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L23Transaxon),somaPosition(-497.561,421.434,365.476)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L34axon),somaPosition(-498.758,413.717,173.791)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L45Peakaxon),somaPosition(-482.559,413.349,-64.0819)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L45Symaxon),somaPosition(-506.972,393.315,-166.148)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L4pyaxon),somaPosition(-497.719,409.418,78.1653)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L4spaxon),somaPosition(-512.758,406.548,56.4845)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L4ssaxon),somaPosition(-509.832,402.907,-5.75627)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L56Transaxon),somaPosition(-468.8,386.108,-577.688)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L5staxon),somaPosition(-485.141,394.836,-304.559)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L5ttaxon),somaPosition(-480.053,392.482,-387.328)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L6ccaxon),somaPosition(-470.169,382.171,-580.465)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L6ccinvaxon),somaPosition(-462.099,383.09,-703.325)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,L6ctaxon),somaPosition(-461.313,376.173,-842.341)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C1,SymLocalaxon),somaPosition(-497.488,406.064,-78.2332)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L1axon),somaPosition(392.63,522.242,557.183)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L2axon),somaPosition(389.859,468.021,357.04)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L23Transaxon),somaPosition(388.819,465.383,274.949)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L34axon),somaPosition(352.297,422.563,73.9589)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L45Peakaxon),somaPosition(340.415,385.803,-155.944)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L45Symaxon),somaPosition(316.51,374.995,-273.349)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L4pyaxon),somaPosition(348.333,407.677,-18.2404)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L4spaxon),somaPosition(349.537,406.872,-45.7528)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L4ssaxon),somaPosition(339.702,394.506,-101.236)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L56Transaxon),somaPosition(280.997,284.422,-718.165)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L5staxon),somaPosition(305.577,341.613,-409.25)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L5ttaxon),somaPosition(297.306,322.779,-505.623)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L6ccaxon),somaPosition(269.756,287.288,-725.12)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L6ccinvaxon),somaPosition(258.811,266.367,-865.561)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,L6ctaxon),somaPosition(242.791,241.482,-1012.03)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(C3,SymLocalaxon),somaPosition(325.804,383.398,-173.889)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L1axon),somaPosition(-455.613,12.6841,623.263)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L2axon),somaPosition(-422.337,26.5812,428.075)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L23Transaxon),somaPosition(-428.602,27.6627,358.004)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L34axon),somaPosition(-397.779,39.9939,156.644)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L45Peakaxon),somaPosition(-374.31,41.698,-115.293)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L45Symaxon),somaPosition(-400.394,38.6471,-184.837)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L4pyaxon),somaPosition(-391.282,43.4998,58.6822)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L4spaxon),somaPosition(-421.238,35.5006,41.2472)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L4ssaxon),somaPosition(-419.182,36.8277,-22.4172)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L56Transaxon),somaPosition(-369.53,50.3987,-621.067)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L5staxon),somaPosition(-383.619,43.9711,-338.044)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L5ttaxon),somaPosition(-384.013,54.4897,-420.588)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L6ccaxon),somaPosition(-371.357,51.7178,-627.432)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L6ccinvaxon),somaPosition(-371.319,52.2426,-750.961)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,L6ctaxon),somaPosition(-365.151,60.8253,-898.724)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D1,SymLocalaxon),somaPosition(-390.632,38.6198,-148.713)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L1axon),somaPosition(35.3985,-76.4903,619.542)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L2axon),somaPosition(7.63309,-21.279,416.579)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L23Transaxon),somaPosition(3.89785,-5.11769,337.288)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L34axon),somaPosition(7.24448,-7.56717,130.94)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L45Peakaxon),somaPosition(-3.87959,-17.1745,-114.551)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L45Symaxon),somaPosition(13.4235,-10.5902,-230.57)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L4pyaxon),somaPosition(12.6718,-12.9464,32.9645)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L4spaxon),somaPosition(3.83879,1.09866,16.1632)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L4ssaxon),somaPosition(2.15841,1.44619,-48.7626)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L56Transaxon),somaPosition(-4.33361,-9.76969,-676.358)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L5staxon),somaPosition(5.51313,-5.96642,-372.964)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L5ttaxon),somaPosition(8.39872,-8.29476,-471.132)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L6ccaxon),somaPosition(0.714012,-12.5874,-687.476)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L6ccinvaxon),somaPosition(-5.86441,-6.79541,-821.327)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,L6ctaxon),somaPosition(-3.37389,-10.409,-973.46)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D2,SymLocalaxon),somaPosition(3.48909,-5.9456,-187.761)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L1axon),somaPosition(452.237,11.3112,605.226)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L2axon),somaPosition(451.244,1.37625,396.179)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L23Transaxon),somaPosition(446.772,6.00185,295.102)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L34axon),somaPosition(437.294,-13.084,87.3785)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L45Peakaxon),somaPosition(433.306,-22.9179,-182.793)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L45Symaxon),somaPosition(426.122,-29.9964,-255.563)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L4pyaxon),somaPosition(435.575,-22.7169,-9.41535)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L4spaxon),somaPosition(416.374,-2.66745,-31.4017)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L4ssaxon),somaPosition(413.633,-4.72932,-92.0035)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L56Transaxon),somaPosition(387.667,-48.4738,-741.718)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L5staxon),somaPosition(409.591,-36.5382,-427.957)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L5ttaxon),somaPosition(400.297,-46.4944,-523.873)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L6ccaxon),somaPosition(393.329,-62.9091,-758.924)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L6ccinvaxon),somaPosition(375.678,-72.6315,-898.987)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,L6ctaxon),somaPosition(368.187,-70.8495,-1053.06)));
	celltypeLabels2SomaPos.insert(std::make_pair< ColumnCellTypePair, somaPosition>(ColumnCellTypePair(D3,SymLocalaxon),somaPosition(419.083,-22.2313,-191.71)));
}

/* Lists with all CellType and Column Indices */
void initializeLists()
{
	if(PostCelltypeList.size())
		PostCelltypeList.clear();
	PostCelltypeList.push_back(L2);
	PostCelltypeList.push_back(L34);
	PostCelltypeList.push_back(L4py);
	PostCelltypeList.push_back(L4sp);
	PostCelltypeList.push_back(L4ss);
	PostCelltypeList.push_back(L5st);
	PostCelltypeList.push_back(L5tt);
	PostCelltypeList.push_back(L6cc);
	PostCelltypeList.push_back(L6ccinv);
	PostCelltypeList.push_back(L6ct);
	PostCelltypeList.push_back(SymLocal);
	PostCelltypeList.push_back(L1);
	PostCelltypeList.push_back(L23Trans);
	PostCelltypeList.push_back(L45Sym);
	PostCelltypeList.push_back(L45Peak);
	PostCelltypeList.push_back(L56Trans);

	if(PreCelltypeList.size())
		PreCelltypeList.clear();
	PreCelltypeList.push_back(L2axon);
	PreCelltypeList.push_back(L34axon);
	PreCelltypeList.push_back(L4pyaxon);
	PreCelltypeList.push_back(L4spaxon);
	PreCelltypeList.push_back(L4ssaxon);
	PreCelltypeList.push_back(L5staxon);
	PreCelltypeList.push_back(L5ttaxon);
	PreCelltypeList.push_back(L6ccaxon);
	PreCelltypeList.push_back(L6ccinvaxon);
	PreCelltypeList.push_back(L6ctaxon);
	PreCelltypeList.push_back(SymLocalaxon);
	PreCelltypeList.push_back(L1axon);
	PreCelltypeList.push_back(L23Transaxon);
	PreCelltypeList.push_back(L45Symaxon);
	PreCelltypeList.push_back(L45Peakaxon);
	PreCelltypeList.push_back(L56Transaxon);

	if(PostColumnList.size())
		PostColumnList.clear();
//	PostColumnList.push_back(B1);
//	PostColumnList.push_back(B2);
//	PostColumnList.push_back(B3);
//	PostColumnList.push_back(C1);
//	PostColumnList.push_back(C2);
//	PostColumnList.push_back(C3);
//	PostColumnList.push_back(D1);
//	PostColumnList.push_back(D2);
//	PostColumnList.push_back(D3);
	PostColumnList.push_back(A1);
	PostColumnList.push_back(A2);
	PostColumnList.push_back(A3);
	PostColumnList.push_back(A4);
	PostColumnList.push_back(B1);
	PostColumnList.push_back(B2);
	PostColumnList.push_back(B3);
	PostColumnList.push_back(B4);
	PostColumnList.push_back(C1);
	PostColumnList.push_back(C2);
	PostColumnList.push_back(C3);
	PostColumnList.push_back(C4);
	PostColumnList.push_back(D1);
	PostColumnList.push_back(D2);
	PostColumnList.push_back(D3);
	PostColumnList.push_back(D4);
	PostColumnList.push_back(E1);
	PostColumnList.push_back(E2);
	PostColumnList.push_back(E3);
	PostColumnList.push_back(E4);
	PostColumnList.push_back(Alpha);
	PostColumnList.push_back(Beta);
	PostColumnList.push_back(Gamma);
	PostColumnList.push_back(Delta);

	if(PreColumnList.size())
		PreColumnList.clear();
	PreColumnList.push_back(A1);
	PreColumnList.push_back(A2);
	PreColumnList.push_back(A3);
	PreColumnList.push_back(A4);
	PreColumnList.push_back(B1);
	PreColumnList.push_back(B2);
	PreColumnList.push_back(B3);
	PreColumnList.push_back(B4);
	PreColumnList.push_back(C1);
	PreColumnList.push_back(C2);
	PreColumnList.push_back(C3);
	PreColumnList.push_back(C4);
	PreColumnList.push_back(D1);
	PreColumnList.push_back(D2);
	PreColumnList.push_back(D3);
	PreColumnList.push_back(D4);
	PreColumnList.push_back(E1);
	PreColumnList.push_back(E2);
	PreColumnList.push_back(E3);
	PreColumnList.push_back(E4);
	PreColumnList.push_back(Alpha);
	PreColumnList.push_back(Beta);
	PreColumnList.push_back(Gamma);
	PreColumnList.push_back(Delta);
}

/* Constant and correspnding name as char string */
void initializeConstants()
{
	if(int2CelltypeLabels.size())
		int2CelltypeLabels.clear();
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L2,"L2"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L34,"L34"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4py,"L4py"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4sp,"L4sp"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4ss,"L4ss"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5st,"L5st"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5tt,"L5tt"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6cc,"L6cc"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccinv,"L6ccinv"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ct,"L6ct"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal,"SymLocal"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal1,"SymLocal1"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal2,"SymLocal2"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal3,"SymLocal3"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal4,"SymLocal4"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal5,"SymLocal5"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal6,"SymLocal6"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L1,"L1"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L23Trans,"L23Trans"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Sym,"L45Sym"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Peak,"L45Peak"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L56Trans,"L56Trans"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L2axon,"L2axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L34axon,"L34axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4pyaxon,"L4pyaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4spaxon,"L4spaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4ssaxon,"L4ssaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5staxon,"L5staxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5ttaxon,"L5ttaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccaxon,"L6ccaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccinvaxon,"L6ccinvaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ctaxon,"L6ctaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocalaxon,"SymLocalaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal1axon,"SymLocal1axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal2axon,"SymLocal2axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal3axon,"SymLocal3axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal4axon,"SymLocal4axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal5axon,"SymLocal5axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal6axon,"SymLocal6axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L1axon,"L1axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L23Transaxon,"L23Transaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Symaxon,"L45Symaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Peakaxon,"L45Peakaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L56Transaxon,"L56Transaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(VPM,"VPM"));

	if(int2ColumnLabels.size())
		int2ColumnLabels.clear();
	int2ColumnLabels.insert(std::pair< int, const char * >(Alpha, "Alpha"));
	int2ColumnLabels.insert(std::pair< int, const char * >(A1, "A1"));
	int2ColumnLabels.insert(std::pair< int, const char * >(A2, "A2"));
	int2ColumnLabels.insert(std::pair< int, const char * >(A3, "A3"));
	int2ColumnLabels.insert(std::pair< int, const char * >(A4, "A4"));
	int2ColumnLabels.insert(std::pair< int, const char * >(Beta, "Beta"));
	int2ColumnLabels.insert(std::pair< int, const char * >(B1, "B1"));
	int2ColumnLabels.insert(std::pair< int, const char * >(B2, "B2"));
	int2ColumnLabels.insert(std::pair< int, const char * >(B3, "B3"));
	int2ColumnLabels.insert(std::pair< int, const char * >(B4, "B4"));
	int2ColumnLabels.insert(std::pair< int, const char * >(Gamma, "Gamma"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C1, "C1"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C2, "C2"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C3, "C3"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C4, "C4"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C5, "C5"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C6, "C6"));
	int2ColumnLabels.insert(std::pair< int, const char * >(Delta, "Delta"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D1, "D1"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D2, "D2"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D3, "D3"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D4, "D4"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D5, "D5"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D6, "D6"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E1, "E1"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E2, "E2"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E3, "E3"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E4, "E4"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E5, "E5"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E6, "E6"));

	if(celltypeLabels2Int.size())
		celltypeLabels2Int.clear();
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L2"),L2));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L34"),L34));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4py"),L4py));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4sp"),L4sp));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4ss"),L4ss));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5st"),L5st));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5tt"),L5tt));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6cc"),L6cc));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccinv"),L6ccinv));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ct"),L6ct));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal"),SymLocal));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal1"),SymLocal1));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal2"),SymLocal2));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal3"),SymLocal3));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal4"),SymLocal4));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal5"),SymLocal5));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal6"),SymLocal6));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L1"),L1));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L23Trans"),L23Trans));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Sym"),L45Sym));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Peak"),L45Peak));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L56Trans"),L56Trans));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L2axon"),L2axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L34axon"),L34axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4pyaxon"),L4pyaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4spaxon"),L4spaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4ssaxon"),L4ssaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5staxon"),L5staxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5ttaxon"),L5ttaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccaxon"),L6ccaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccinvaxon"),L6ccinvaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ctaxon"),L6ctaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocalaxon"),SymLocalaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal1axon"),SymLocal1axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal2axon"),SymLocal2axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal3axon"),SymLocal3axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal4axon"),SymLocal4axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal5axon"),SymLocal5axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal6axon"),SymLocal6axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L1axon"),L1axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L23Transaxon"),L23Transaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Symaxon"),L45Symaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Peakaxon"),L45Peakaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L56Transaxon"),L56Transaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("VPM"),VPM));

	if(ColumnLabels2Int.size())
		ColumnLabels2Int.clear();
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("Alpha"), Alpha));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("A1"), A1));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("A2"), A2));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("A3"), A3));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("A4"), A4));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("Beta"), Beta));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("B1"), B1));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("B2"), B2));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("B3"), B3));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("B4"), B4));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("Gamma"), Gamma));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C1"), C1));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C2"), C2));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C3"), C3));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C4"), C4));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C5"), C5));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C6"), C6));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("Delta"), Delta));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D1"), D1));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D2"), D2));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D3"), D3));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D4"), D4));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D5"), D5));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D6"), D6));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E1"), E1));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E2"), E2));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E3"), E3));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E4"), E4));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E5"), E5));
	ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E6"), E6));
}

///* Connection Probability vs. Distance of Postsynaptic Cell to vertical Column axis */
//int main()
//{
//	initializeConstants();
//	initializeLists();
//
//	std::string inputFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_3x3_final/convergence3x3.csv";
//	//std::string inputFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_3x3_final/divergence3x3.csv";
//
//	CellTable * connectome = new CellTable;
//	Reader * cellTableReader = new Reader(inputFilename.c_str(), inputFilename.c_str());
//	cellTableReader->readConvergenceTableCSV(connectome);
//
//	const char * inputFilename2 = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_3x3_final/nrCells.csv";
//  connectome->mergeLocal(inputFilename2);
//
//	// Extract CellType/Column Combination (postsynaptic)
//	std::list< unsigned int > columns;
//	std::list< unsigned int > cellTypes;
//	std::map< ColumnCellTypePair, unsigned int > header = connectome->header;
//
//	// Extract postsynaptic L4ss cells in C2
//	//columns.push_back(C2);
//	cellTypes.push_back(L4ss);
//
//	// Reference Column (C2)
//	BarrelField * BF = new BarrelField();
//	std::map< int, Column * > avgColumns = BF->avgColumns;
//	Column * C2col = avgColumns.find(C2)->second;
//
//	std::map< int, bin *> binMap;
//	double binsz = 50;
//
//	std::vector< CellTableRow * > rows = connectome->getPostColumnCelltypeRows(columns, cellTypes);
//	// Compute Connection Probability, go through all post-synaptic cells
//	for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
//	{
//		double t, projectedPt[3], somaLocation[3];
//		somaLocation[0] = (*rowit)->somaLocation[0]; // x
//		somaLocation[1] = (*rowit)->somaLocation[1]; // y
//		somaLocation[2] = (*rowit)->somaLocation[2]; // z
//
//		// Compute Distance to Vertical Column Axis
//		double dist = vtkLine::DistanceToLine(somaLocation, C2col->top, C2col->bottom, t, projectedPt);
//		dist = sqrt(dist);
//		int pos = (int) dist/binsz;
//
//		std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;
//
//		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
//		{
//			// Presynaptic Cells
//			if (headerit->first.first==C2 && headerit->first.second==L4ssaxon)
//			{
//				std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
//
//				// Add to histogram if not in histogram yet
//				if (binMap.count(pos) == 0)
//				{
//					bin * b = new bin;
//					b->M = (*Convit);
//					b->N = 1;
//					b->Q = 0;
//					binMap[pos] = b;
//				}
//				else // Update histogram
//				{
//					// Update Mean, SD, and Number of Samples,
//					float tmp1 = (*Convit) - binMap[pos]->M;
//					// Update Mean of correlation r: (Mean * N + t) / (N+1)
//					binMap[pos]->M = (binMap[pos]->M * binMap[pos]->N + (*Convit)) / (binMap[pos]->N+1);
//					// Update N (Number of Samples)
//					binMap[pos]->N = binMap[pos]->N+1;
//					float tmp2 = ((*Convit) - binMap[pos]->M);
//					// Update SD (can be calculated from Q)
//					binMap[pos]->Q = binMap[pos]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
//				}
//			}
//		}
//	}
//	for (std::map< int, bin *>::iterator it = binMap.begin(); it != binMap.end(); ++it)
//	{
//		bin * b = it->second;
//		int N = b->N;
//		float M = b->M;
//		float Q = b->Q;
//
//		std::cout << (it->first) * binsz+binsz/2 << "," << M << "," << N <<std::endl;
//	}
//}

///* Read in Convergences Table, output Average and SD table of Connection Probability for Celltype/Column Combinations (transpose!)
//int main()
//{
//	initializeConstants();
//	initializeLists();
//
//	/* L4ss vs L4sp onto L6cc/L6ccinv/L6ct */
//	std::string inputFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_3x3_final/convergenceApical3x3.csv";
//	CellTable * connectome = new CellTable;
//	Reader * cellTableReader = new Reader(inputFilename.c_str(), inputFilename.c_str());
//	cellTableReader->readConvergenceTableCSV(connectome);
//
//	const char * inputFilename2 = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_3x3_final/nrCells.csv";
//  connectome->mergeLocal(inputFilename2);
//
//	// Extract CellType/Column Combination (postsynaptic)
//	std::list< unsigned int > columns;
//	std::list< unsigned int > cellTypes;
//	std::map< ColumnCellTypePair, unsigned int > header = connectome->header;
//
//	// Extract postsynaptic L6 cells in C2 (in L6A)
//	columns.push_back(C2);
//	cellTypes.push_back(L6cc);
//	cellTypes.push_back(L6ccinv);
//	cellTypes.push_back(L6ct);
//
//	 std::vector< bin *> binList;
//	 // 0: C2-L4ss 1: C2-L4sp 2: C2-L4sp+C2-L4py
//	 for (int i = 0; i < 3; i++)
//	 {
//			bin * b = new bin;
//			b->M = 0;
//			b->N = 0;
//			b->Q = 0;
//			binList.push_back(b);
//	 }
//
//	BarrelField * BF = new BarrelField();
//	std::vector< CellTableRow * > rows = connectome->getPostColumnCelltypeRows(columns, cellTypes);
//	// Compute Connection Probability
//	for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
//	{
//		double somaLocation[3];
//		somaLocation[0] = (*rowit)->somaLocation[0];
//		somaLocation[1] = (*rowit)->somaLocation[1];
//		somaLocation[2] = (*rowit)->somaLocation[2];
//		double dist = BF->piaDistance(somaLocation);
//
//		// L6 1350-1850 -> L6A 1350-1600
////		if (dist>1600)
////		{
////			continue;
////		}
////		if (dist<1350)
////		{
////			continue;
////		}
//
//		std::vector< float > convergenceVal = (*rowit)->synapsesPerPreTypeColumn;
//
//		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
//		{
//			int pos = -1;
//
//			if (headerit->first.first==C2 && headerit->first.second==L4ssaxon)
//			{
//				pos = 0;
//			}
//
//			if (headerit->first.first==C2 && headerit->first.second==L4spaxon)
//			{
//				pos = 1;
//			}
//
//			if (headerit->first.first==C2 && headerit->first.second==L4pyaxon)
//			{
//				pos = 2;
//			}
//
//			if (pos != -1)
//			{
//				std::vector<float>::iterator Convit = convergenceVal.begin()+(headerit->second);
//				std::vector< bin *>::iterator Binit = binList.begin()+pos;
//
//				// Update Mean, SD, and Number of Samples,
//				float tmp1 = (*Convit) - (*Binit)->M;
//				// Update Mean of correlation r: (Mean * N + t) / (N+1)
//				(*Binit)->M = ((*Binit)->M * (*Binit)->N + (*Convit)) / ((*Binit)->N+1);
//				// Update N (Number of Samples)
//				(*Binit)->N = (*Binit)->N+1;
//				float tmp2 = ((*Convit) - (*Binit)->M);
//				// Update SD (can be calculated from Q)
//				(*Binit)->Q = (*Binit)->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
//
//				if(pos==1)
//				{
//					Binit = binList.begin()+2;
//					// Update Mean, SD, and Number of Samples,
//					float tmp1 = (*Convit) - (*Binit)->M;
//					// Update Mean of correlation r: (Mean * N + t) / (N+1)
//					(*Binit)->M = ((*Binit)->M * (*Binit)->N + (*Convit)) / ((*Binit)->N+1);
//					// Update N (Number of Samples)
//					(*Binit)->N = (*Binit)->N+1;
//					float tmp2 = ((*Convit) - (*Binit)->M);
//					// Update SD (can be calculated from Q)
//					(*Binit)->Q = (*Binit)->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
//				}
//			}
//		}
//	}
//
//	std::cout << "Postsynapse L6" << std::endl;
//	std::cout << "Presynapse,M,SD" << std::endl;
//	for (int i = 0; i < 3; i++)
//	{
//		int N = binList.at(i)->N;
//		float Q = binList.at(i)->Q;
//		float M = binList.at(i)->M;
//		std::cout << "L4," << M << "," << sqrt(Q/(N-1)) << std::endl;
//	}
//}

//int main()
//{
//	initializeConstants();
//	initializeLists();
//
//	/* Read in Convergence Table */
//	const char * inputFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_3x3_final/convergence3x3.csv";
//	CellTable * connectome = new CellTable;
//	Reader * cellTableReader = new Reader(inputFilename, inputFilename);
//	cellTableReader->readConvergenceTableCSV(connectome);
//	const char * inputFilename2 = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_3x3_final/nrCells.csv";
//  connectome->mergeLocal(inputFilename2);
//
//	// Extract CellType/Column Combination (postsynaptic)
//	std::list< unsigned int > columns;
//	std::list< unsigned int > cellTypes;
//	cellTypes.push_back(L4ss);
//	std::vector< CellTableRow * > rows = connectome->getPostColumnCelltypeRows(columns, cellTypes);
//	std::map< ColumnCellTypePair, unsigned int > header = connectome->header;
//
//	// Presynaptic CellType/Column Combination
//	ColumnCellTypePair pair = ColumnCellTypePair(C2, L4ssaxon);
//	//ColumnCellTypePair pair = ColumnCellTypePair(C3, L4ssaxon);
//	int headerpos = header.find(pair)->second;
//
//	// Reference Column (C2)
//	BarrelField * BF = new BarrelField();
//	std::map< int, Column * > avgColumns = BF->avgColumns;
//	Column * C2col = avgColumns.find(C2)->second;
//
//	std::map< int, bin *> binMap;
//	int binsz = 50;
//
//	// Compute Connection Probability
//	for (std::vector< CellTableRow * >::iterator rowit = rows.begin(); rowit != rows.end(); ++rowit)
//	{
//		double t, projectedPt[3], somaLocation[3];
//		somaLocation[0] = (*rowit)->somaLocation[0];
//		somaLocation[1] = (*rowit)->somaLocation[1];
//		somaLocation[2] = (*rowit)->somaLocation[2];
//
//		// Distance to C2 Vertical Column Axis
//		double dist = vtkLine::DistanceToLine(somaLocation, C2col->top, C2col->bottom, t, projectedPt);
//		dist = sqrt(dist);
//		int pos = (int) dist/binsz;
//
//		// Extract Convergence Value of Presynaptic Col/CellType Pair
//		float c = (*rowit)->synapsesPerPreTypeColumn.at(headerpos);
//
//		// Add to histogram if not in histogram yet
//		if (binMap.count(pos) == 0)
//		{
//			bin * b = new bin;
//			b->M = c;
//			b->N = 1;
//			b->Q = 0;
//			binMap[pos] = b;
//		}
//		else // Update histogram
//		{
//			// Update Mean, SD, and Number of Samples,
//			float tmp1 = c - binMap[pos]->M;
//			// Update Mean of correlation r: (Mean * N + t) / (N+1)
//			binMap[pos]->M = (binMap[pos]->M * binMap[pos]->N + c) / (binMap[pos]->N+1);
//			// Update N (Number of Samples)
//			binMap[pos]->N = binMap[pos]->N+1;
//			float tmp2 = (c - binMap[pos]->M);
//			// Update SD (can be calculated from Q)
//			binMap[pos]->Q = binMap[pos]->Q + tmp1*tmp2; // Q + (rtmp - M_prev)(rtmp - M)
//		}
//	}
//
//	for (std::map< int, bin *>::iterator it = binMap.begin(); it != binMap.end(); ++it)
//	{
//		bin * b = it->second;
//		int N = b->N;
//		float M = b->M;
//		float Q = b->Q;
//
//		//std::cout << (it->first) * binsz+binsz/2 << "," << M << "," << sqrt(Q/(N-1)) << "," << N << std::endl;
//		std::cout << (it->first) * binsz+binsz/2 << "," << M << "," << N <<std::endl;
//	}
//}
