/****************************************************************************/
/*                                                                          */
/* Program:   ConnectionMatrixVisualization                                 */
/*                                                                          */
/* File:      uniqueConnectionMatrixRows.cpp                                */
/*                                                                          */
/* Purpose:   load NeuroNet connection matrix and save only unique rows     */
/*            (i.e. individual presynaptic axons)                           */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max Planck Institute for Biological Cybernetics               */
/*            Spemannstr. 38-44                                             */
/*            72076 Tuebingen                                               */
/*            Germany                                                       */
/*                                                                          */
/* EMail:     robert.egger@tuebingen.mpg.de                                 */
/*                                                                          */
/* History:   16.06.2014                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "../../common/typedefs.h"
#include "../../common/amiraReader.h"
#include "matrixanalyzer.h"

#define DEBUG

SelectionType getPresynapticCellIDs(ConnectionMatrix * connectome);
SelectionType getPostsynapticCellIDs(ConnectionMatrix * connectome);
std::vector< float > getDenseConnectionMatrixRow(ConnectionMatrix * connectome, unsigned int presynapticID, SelectionType postsynapticIDs);
float absoluteDifference(std::vector< float > x, std::vector< float > y);
float innerProductNorm(std::vector< float > x, std::vector< float > y);
bool equalWithTolerance(std::vector< float > x, std::vector< float > y, float epsilon);
unsigned int findCelltype(unsigned int cellID, ConnectionMatrix * connectome);
std::map< unsigned int, unsigned int > readAxonDendriteIDFile(const char * filename);
void writeUniqueMatrixRows(ConnectionMatrix * connectome, std::map< unsigned int, SelectionType > uniqueIDsPerCelltype, const char * outputFilename);
void writeCompleteMatrixRows(ConnectionMatrix* connectome, std::map< unsigned int, SelectionType > presynapticCellIDsPerCellType, SelectionType postsynapticCellIDs,
							 std::map< unsigned int, unsigned int > axonDendriteCorrespondence, const char* outputFilename);
void writeConnectivitySummary(ConnectionMatrix* connectome, std::map< unsigned int, SelectionType > presynapticCellIDsPerCellType, SelectionType postsynapticCellIDs,
							 std::map< unsigned int, unsigned int > axonDendriteCorrespondence, const char* outputFilename);

void initializeConstants();
std::map< unsigned int, const char * > int2CelltypeLabels;
std::map< unsigned int, const char * > int2ColumnLabels;

int main( int argc , char * argv[])
{
	if(argc == 2)
	{
		initializeConstants();
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[1];
// 		float epsilon = atof(argv[2]);
		float epsilon = 1e-5;
		ConnectionMatrix * connectome = new ConnectionMatrix;
		Reader * matrixReader = new Reader(inputFilename);
		matrixReader->readConnectionMatrix(connectome);
		
		SelectionType presynapticIDs = getPresynapticCellIDs(connectome);
		SelectionType postsynapticIDs = getPostsynapticCellIDs(connectome);
// 		SelectionType presynapticIDs = connectome->getPreColumnCelltypeSelection(D2, VPM);
// 		SelectionType postsynapticIDs = connectome->getPostColumnCelltypeSelection(D2, L4ss);
		std::vector< std::vector< float > > matrixRows;
		std::vector< unsigned int > matrixRowIDs;
		std::cout << "Populating dense matrix rows:" << std::endl;
		std::cout << presynapticIDs.size() << " rows, " << postsynapticIDs.size() << " columns" << std::endl;
		for(int ii = 0; ii < presynapticIDs.size(); ++ii)
		{
			std::vector< float > row = getDenseConnectionMatrixRow(connectome, presynapticIDs[ii], postsynapticIDs);
			matrixRows.push_back(row);
			matrixRowIDs.push_back(presynapticIDs[ii]);
		}
		
		std::list< unsigned int > postCelltypeList;
		postCelltypeList.push_back(L34);
		postCelltypeList.push_back(L4py);
		postCelltypeList.push_back(L4sp);
		postCelltypeList.push_back(L4ss);
		postCelltypeList.push_back(L45Peak);
		postCelltypeList.push_back(L45Sym);
		SelectionType preIDsVPM = connectome->preTypeIDs[VPM];
		
		std::cout << "Cell type\tavg VPM innervation\tSTD VPM innervation\tnumber" << std::endl;
		std::list< unsigned int >::const_iterator postCelltypeListIt;
		for(postCelltypeListIt = postCelltypeList.begin(); postCelltypeListIt != postCelltypeList.end(); ++postCelltypeListIt)
		{
			SelectionType postIDsVPM = connectome->postTypeIDs[*postCelltypeListIt];
			std::vector< float > individualCellValues;
			float avgInnervation = 0, stdInnervation = 0;
			for(int ii = 0; ii < postIDsVPM.size(); ++ii)
			{
				float singleCellInnervation = 0;
				for(int jj = 0; jj < preIDsVPM.size(); ++jj)
				{
					MatrixIndexType innervationIndex(preIDsVPM[jj], postIDsVPM[ii]);
					if(connectome->matrix.find(innervationIndex) != connectome->matrix.end())
					{
						singleCellInnervation += connectome->matrix[innervationIndex];
					}
				}
				avgInnervation += singleCellInnervation;
				individualCellValues.push_back(singleCellInnervation);
			}
			avgInnervation = avgInnervation/postIDsVPM.size();
			for(int ii = 0; ii < individualCellValues.size(); ++ii)
			{
				stdInnervation += (avgInnervation - individualCellValues[ii])*(avgInnervation - individualCellValues[ii]);
			}
			stdInnervation = sqrt(stdInnervation/(postIDsVPM.size() - 1));
			std::cout << int2CelltypeLabels[*postCelltypeListIt] << "\t" << avgInnervation << "\t" << stdInnervation << "\t" << postIDsVPM.size() << std::endl;
			
// 			if(*postCelltypeListIt == L45Sym)
// 			{
// 				std::cout << "********************************" << std::endl;
// 				std::cout << "L45Sym individual values" << std::endl;
// 				std::cout << "post ID\tinnervation" << std::endl;
// 				for(int ii = 0; ii < postIDsVPM.size(); ++ii)
// 				{
// 					std::cout << postIDsVPM[ii] << "\t" << individualCellValues[ii] << std::endl;
// 				}
// 				std::cout << "********************************" << std::endl;
// 			}
		}
		
		std::cout << "Checking if " << matrixRows.size() << " rows are equal with tolerance " << epsilon << std::endl;
		std::list< unsigned int > duplicateRows;
		std::list< unsigned int > duplicateIDs;
		unsigned int equal = 0, notEqual = 0;
		for(int ii = 0; ii < matrixRows.size(); ++ii)
		{
			for(int jj = ii+1; jj < matrixRows.size(); ++jj)
			{
				if(equalWithTolerance(matrixRows[ii], matrixRows[jj], epsilon))
				{
					++equal;
					duplicateRows.push_back(jj);
					duplicateIDs.push_back(matrixRowIDs[jj]);
				}
				else
				{
					++notEqual;
				}
			}
		}
		duplicateIDs.sort();
		duplicateIDs.unique();
		duplicateRows.sort();
		duplicateRows.unique();
		
		std::cout << "Found " << duplicateIDs.size() << " duplicate IDs" << std::endl;
		std::cout << "Found " << duplicateRows.size() << " duplicate rows" << std::endl;
		std::cout << "Found " << equal << " equal row pairs; " << notEqual << " not equal row pairs" << std::endl;
		
		std::map< unsigned int, SelectionType > uniqueIDsPerCelltype;
		for(int ii = 0; ii < presynapticIDs.size(); ++ii)
		{
			unsigned int presynapticID = presynapticIDs[ii];
			if(std::find(duplicateIDs.begin(), duplicateIDs.end(), presynapticID) == duplicateIDs.end())
			{
				unsigned int cellType = findCelltype(presynapticID, connectome);
// 				unsigned int cellType = VPM;
				if(uniqueIDsPerCelltype.find(cellType) == uniqueIDsPerCelltype.end())
				{
					SelectionType typeUniqueIDs;
					typeUniqueIDs.push_back(presynapticID);
					uniqueIDsPerCelltype.insert(std::pair< unsigned int, SelectionType >(cellType, typeUniqueIDs));
				}
				else
				{
					uniqueIDsPerCelltype[cellType].push_back(presynapticID);
				}
			}
		}
#ifdef DEBUG
		std::map< unsigned int, SelectionType >::const_iterator uniqueIDsPerCelltypeIt;
		std::cout << "Cell type\tunique rows" << std::endl;
		for(uniqueIDsPerCelltypeIt = uniqueIDsPerCelltype.begin(); uniqueIDsPerCelltypeIt != uniqueIDsPerCelltype.end(); ++uniqueIDsPerCelltypeIt)
		{
			const char * cellTypeLabel = int2CelltypeLabels[uniqueIDsPerCelltypeIt->first];
			unsigned int nrOfUniqueRows = uniqueIDsPerCelltypeIt->second.size();
			std::cout << cellTypeLabel << "\t" << nrOfUniqueRows << std::endl;
			
			std::string outputFilename(inputFilename);
			outputFilename += "_";
			outputFilename += cellTypeLabel;
			outputFilename += "_all_unique.tif";
			connectome->writeConnectionMatrixAsImage(outputFilename.c_str(), uniqueIDsPerCelltypeIt->second, postsynapticIDs);
		}
// 		std::string outputFilename(inputFilename);
// 		outputFilename += "_all_all_unique.tif";
// 		SelectionType L4ssDendSelection = connectome->getPostColumnCelltypeSelection(D2, L4ss);
// 		connectome->writeConnectionMatrixAsImage(outputFilename.c_str(), uniqueIDsPerCelltype[VPM], L4ssDendSelection);
// 		connectome->writeConnectionMatrixAsImage(outputFilename.c_str(), uniqueIDsPerCelltype[VPM], postsynapticIDs);
#endif
		writeUniqueMatrixRows(connectome, uniqueIDsPerCelltype, outputFilename);
		
		delete matrixReader, delete connectome;
	}
	
	if(argc == 3)
	{
		initializeConstants();
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[1];
		const char * axonDendIDFilename = argv[2];
// 		float epsilon = atof(argv[2]);
		float epsilon = 1e-5;
		ConnectionMatrix * connectome = new ConnectionMatrix;
		Reader * matrixReader = new Reader(inputFilename);
		matrixReader->readConnectionMatrix(connectome);
		std::map< unsigned int, unsigned int > axonDendriteCorrespondence = readAxonDendriteIDFile(axonDendIDFilename);
#ifdef DEBUG
		std::cout << "axonDendriteCorrespondence.size(): " << axonDendriteCorrespondence.size() << std::endl;
		std::cout << "First element:" << std::endl;
		std::cout << axonDendriteCorrespondence.begin()->first << " --- " << axonDendriteCorrespondence.begin()->second << std::endl;
		std::cout << "Last element:" << std::endl;
		std::cout << axonDendriteCorrespondence.rbegin()->first << " --- " << axonDendriteCorrespondence.rbegin()->second << std::endl;
#endif
		
		std::map< unsigned int, unsigned int > cellTypeNumbers;
		cellTypeNumbers[L34] = 0;
		cellTypeNumbers[L4py] = 0;
		cellTypeNumbers[L4sp] = 0;
		cellTypeNumbers[L4ss] = 0;
		cellTypeNumbers[SymLocal] = 0;
		cellTypeNumbers[L23Trans] = 0;
		cellTypeNumbers[L45Peak] = 0;
		cellTypeNumbers[L45Sym] = 0;
		cellTypeNumbers[VPM] = connectome->preTypeNumbers[VPM];
		cellTypeNumbers[L34axon] = 0;
		cellTypeNumbers[L4pyaxon] = 0;
		cellTypeNumbers[L4spaxon] = 0;
		cellTypeNumbers[L4ssaxon] = 0;
		cellTypeNumbers[SymLocalaxon] = 0;
		cellTypeNumbers[L23Transaxon] = 0;
		cellTypeNumbers[L45Peakaxon] = 0;
		cellTypeNumbers[L45Symaxon] = 0;
		
		SelectionType presynapticIDs = getPresynapticCellIDs(connectome);
		SelectionType postsynapticIDs = getPostsynapticCellIDs(connectome);
		
		SelectionType isolatedPreIDs;
		SelectionType isolatedPostIDs;
		std::map< unsigned int, SelectionType > completePreIDs;
		SelectionType completePostIDs;
		for(int i = 0; i < presynapticIDs.size(); ++i)
		{
			unsigned int preID = presynapticIDs[i];
			if(axonDendriteCorrespondence.find(preID) != axonDendriteCorrespondence.end())
			{
				unsigned int correspondingID = axonDendriteCorrespondence[preID];
				if(std::find(postsynapticIDs.begin(), postsynapticIDs.end(), correspondingID) == postsynapticIDs.end())
				{
					isolatedPreIDs.push_back(preID);
				}
				else
				{
					unsigned int cellType = connectome->IDColumnCelltypeMap[preID].second;
					if(cellType >= SymLocal1 && cellType <= SymLocal6)
					{
						cellType = SymLocal;
					}
					if(cellType >= SymLocal1axon && cellType <= SymLocal6axon)
					{
						cellType = SymLocalaxon;
					}
					cellTypeNumbers[cellType]++;
					if(completePreIDs.find(cellType) != completePreIDs.end())
					{
						completePreIDs[cellType].push_back(preID);
					}
					else
					{
						SelectionType newIDVec;
						newIDVec.push_back(preID);
						completePreIDs[cellType] = newIDVec;
					}
				}
			}
			else
			{
				if(connectome->IDColumnCelltypeMap[preID].second != VPM)
				{
					isolatedPreIDs.push_back(preID);
				}
			}
		}
		for(int i = 0; i < postsynapticIDs.size(); ++i)
		{
			unsigned int postID = postsynapticIDs[i];
			std::map< unsigned int, unsigned int >::const_iterator axonDendriteCorrespondenceIt;
			for(axonDendriteCorrespondenceIt = axonDendriteCorrespondence.begin();
				axonDendriteCorrespondenceIt != axonDendriteCorrespondence.end(); ++axonDendriteCorrespondenceIt)
			{
				if(axonDendriteCorrespondenceIt->second == postID)
				{
					unsigned int correspondingID = axonDendriteCorrespondenceIt->first;
					if(std::find(presynapticIDs.begin(), presynapticIDs.end(), correspondingID) == presynapticIDs.end())
					{
						isolatedPostIDs.push_back(postID);
					}
					else
					{
						unsigned int cellType = connectome->IDColumnCelltypeMap[postID].second;
						if(cellType >= SymLocal1 && cellType <= SymLocal6)
						{
							cellType = SymLocal;
						}
						if(cellType >= SymLocal1axon && cellType <= SymLocal6axon)
						{
							cellType = SymLocalaxon;
						}
						cellTypeNumbers[cellType]++;
						completePostIDs.push_back(postID);
					}
				}
			}
		}
		completePreIDs[VPM] = connectome->preTypeIDs[VPM];
		
		std::cout << "Found " << isolatedPreIDs.size() << " isolated presynaptic cells" << std::endl;
#ifdef DEBUG
		std::cout << "ID\tColumn\tCell type" << std::endl;
		for(int i = 0; i < isolatedPreIDs.size(); ++i)
		{
			unsigned int preID = isolatedPreIDs[i];
			std::pair< unsigned int, unsigned int > columnCellType = connectome->IDColumnCelltypeMap[preID];
			unsigned int column = columnCellType.first;
			unsigned int cellType = columnCellType.second;
			std::cout << preID << "\t" << int2ColumnLabels[column] << "\t" << int2CelltypeLabels[cellType] << std::endl;
		}
#endif
		std::cout << "Found " << isolatedPostIDs.size() << " isolated postsynaptic cells" << std::endl;
#ifdef DEBUG
		std::cout << "ID\tColumn\tCell type" << std::endl;
		for(int i = 0; i < isolatedPostIDs.size(); ++i)
		{
			unsigned int postID = isolatedPostIDs[i];
			std::pair< unsigned int, unsigned int > columnCellType = connectome->IDColumnCelltypeMap[postID];
			unsigned int column = columnCellType.first;
			unsigned int cellType = columnCellType.second;
			std::cout << postID << "\t" << int2ColumnLabels[column] << "\t" << int2CelltypeLabels[cellType] << std::endl;
		}
#endif
		std::cout << "Complete matrix:" << std::endl;
		std::cout << "Cell type\tNr. of cells" << std::endl;
		std::map< unsigned int, unsigned int >::const_iterator cellTypeNumbersIt;
		for(cellTypeNumbersIt = cellTypeNumbers.begin(); cellTypeNumbersIt != cellTypeNumbers.end(); ++cellTypeNumbersIt)
		{
			std::cout << int2CelltypeLabels[cellTypeNumbersIt->first] << "\t" << cellTypeNumbersIt->second << std::endl;
		}
		
		std::string summaryOutStr(outputFilename);
		summaryOutStr += "_connectivity_summary.csv";
		writeConnectivitySummary(connectome, completePreIDs, completePostIDs, axonDendriteCorrespondence, summaryOutStr.c_str());
		
		std::string matrixOutStr(outputFilename);
		matrixOutStr += "_complete_innervation_matrix.csv";
		writeCompleteMatrixRows(connectome, completePreIDs, completePostIDs, axonDendriteCorrespondence, matrixOutStr.c_str());
// 		
// 		std::cout << "Checking if " << matrixRows.size() << " rows are equal with tolerance " << epsilon << std::endl;
// 		std::list< unsigned int > duplicateRows;
// 		std::list< unsigned int > duplicateIDs;
// 		unsigned int equal = 0, notEqual = 0;
// 		for(int ii = 0; ii < matrixRows.size(); ++ii)
// 		{
// 			for(int jj = ii+1; jj < matrixRows.size(); ++jj)
// 			{
// 				if(equalWithTolerance(matrixRows[ii], matrixRows[jj], epsilon))
// 				{
// 					++equal;
// 					duplicateRows.push_back(jj);
// 					duplicateIDs.push_back(matrixRowIDs[jj]);
// 				}
// 				else
// 				{
// 					++notEqual;
// 				}
// 			}
// 		}
// 		duplicateIDs.sort();
// 		duplicateIDs.unique();
// 		duplicateRows.sort();
// 		duplicateRows.unique();
// 		
// 		std::cout << "Found " << duplicateIDs.size() << " duplicate IDs" << std::endl;
// 		std::cout << "Found " << duplicateRows.size() << " duplicate rows" << std::endl;
// 		std::cout << "Found " << equal << " equal row pairs; " << notEqual << " not equal row pairs" << std::endl;
// 		
// 		std::map< unsigned int, SelectionType > uniqueIDsPerCelltype;
// 		for(int ii = 0; ii < presynapticIDs.size(); ++ii)
// 		{
// 			unsigned int presynapticID = presynapticIDs[ii];
// 			if(std::find(duplicateIDs.begin(), duplicateIDs.end(), presynapticID) == duplicateIDs.end())
// 			{
// 				unsigned int cellType = findCelltype(presynapticID, connectome);
// // 				unsigned int cellType = VPM;
// 				if(uniqueIDsPerCelltype.find(cellType) == uniqueIDsPerCelltype.end())
// 				{
// 					SelectionType typeUniqueIDs;
// 					typeUniqueIDs.push_back(presynapticID);
// 					uniqueIDsPerCelltype.insert(std::pair< unsigned int, SelectionType >(cellType, typeUniqueIDs));
// 				}
// 				else
// 				{
// 					uniqueIDsPerCelltype[cellType].push_back(presynapticID);
// 				}
// 			}
// 		}
// 		writeUniqueMatrixRows(connectome, uniqueIDsPerCelltype, outputFilename);
		
		delete matrixReader, delete connectome;
	}
	
	else if(argc != 2 && argc != 3)
	{
		std::cout << "Parameters: [inputFilename] [outputFilename]" << std::endl;
	}
	
	return 0;
}

SelectionType getPresynapticCellIDs(ConnectionMatrix* connectome)
{
	SelectionType presynapticCellIDs;
	std::map< unsigned int, SelectionType >::const_iterator preTypeIDsIt;
	for(preTypeIDsIt = connectome->preTypeIDs.begin(); preTypeIDsIt != connectome->preTypeIDs.end(); ++preTypeIDsIt)
	{
		SelectionType typeIDs = preTypeIDsIt->second;
		for(int ii = 0; ii < typeIDs.size(); ++ii)
		{
			presynapticCellIDs.push_back(typeIDs[ii]);
		}
	}
	return presynapticCellIDs;
}

SelectionType getPostsynapticCellIDs(ConnectionMatrix* connectome)
{
	SelectionType postsynapticCellIDs;
	std::map< unsigned int, SelectionType >::const_iterator postTypeIDsIt;
	for(postTypeIDsIt = connectome->postTypeIDs.begin(); postTypeIDsIt != connectome->postTypeIDs.end(); ++postTypeIDsIt)
	{
		SelectionType typeIDs = postTypeIDsIt->second;
		for(int ii = 0; ii < typeIDs.size(); ++ii)
		{
			postsynapticCellIDs.push_back(typeIDs[ii]);
		}
	}
	return postsynapticCellIDs;
}

std::vector< float > getDenseConnectionMatrixRow(ConnectionMatrix* connectome, unsigned int presynapticID, SelectionType postsynapticIDs)
{
	std::vector< float > innervationRow;
	for(int ii = 0; ii < postsynapticIDs.size(); ++ii)
	{
		float innervation = 0;
		unsigned int postsynapticID = postsynapticIDs[ii];
		std::map< MatrixIndexType, float >::const_iterator matrixIt = connectome->matrix.find(MatrixIndexType(presynapticID, postsynapticID));
		if(matrixIt != connectome->matrix.end())
		{
			innervation = matrixIt->second;
		}
		innervationRow.push_back(innervation);
	}
	return innervationRow;
}

float absoluteDifference(std::vector< float > x, std::vector< float > y)
{
	if(x.size() != y.size())
	{
		std::cout << "Error: Input vectors are of different dimension!" << std::endl;
		return -2.0;
	}
	
	float diff = 0;
	for(int ii = 0; ii < x.size(); ++ii)
	{
		diff += fabs(x[ii] - y[ii]);
	}
	if(x.size())
	{
		diff = diff/x.size();
	}
	return diff;
}

float innerProductNorm(std::vector< float > x, std::vector< float > y)
{
	if(x.size() != y.size())
	{
		std::cout << "Error: Input vectors are of different dimension!" << std::endl;
		return -2.0;
	}
	
	double inner = 0.0;
	double normx = 0.0, normy = 0.0;
	for(int ii = 0; ii < x.size(); ++ii)
	{
		double x_ = (double)(x[ii]);
		double y_ = (double)(y[ii]);
		inner += x_*y_;
		normx += x_*x_;
		normy += y_*y_;
	}
// 	normx = sqrt(normx);
// 	normy = sqrt(normy);
// 	if(normx && normy)
// 	{
// 		inner = inner/(normx*normy);
// 	}
// 	projection
	if(normx)
	{
		inner = inner/normx;
	}
	return (float)inner;
}

bool equalWithTolerance(std::vector< float > x, std::vector< float > y, float epsilon)
{
	if(x.size() != y.size())
	{
		std::cout << "Error: Input vectors are of different dimension!" << std::endl;
		return false;
	}
	for(int ii = 0; ii < x.size(); ++ii)
	{
// 		float diff = fabs(x[ii] - y[ii]);
// 		float norm = std::max(fabs(x[ii]), fabs(y[ii]));
// 		float rel = norm > 0 ? diff/norm : diff;
		if(fabs(x[ii] - y[ii]) > epsilon)
// 		if(rel > epsilon)
		{
			return false;
		}
	}
	return true;
}

unsigned int findCelltype(unsigned int cellID, ConnectionMatrix* connectome)
{
	unsigned int cellType = 0;
	std::map< unsigned int, SelectionType >::const_iterator preTypeIDsIt;
	for(preTypeIDsIt = connectome->preTypeIDs.begin(); preTypeIDsIt != connectome->preTypeIDs.end(); ++preTypeIDsIt)
	{
		SelectionType typeIDs = preTypeIDsIt->second;
		if(std::find(typeIDs.begin(), typeIDs.end(), cellID) != typeIDs.end())
		{
			cellType = preTypeIDsIt->first;
			break;
		}
	}
	return cellType;
}

std::map< unsigned int, unsigned int > readAxonDendriteIDFile(const char* filename)
{
	std::map< unsigned int, unsigned int > axonDendriteCorrespondence;
	std::ifstream inputStream(filename);
	if(!inputStream.fail())
	{
		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				if(currentLine.find("#") == 0)
				{
					continue;
				}
				
				size_t delim = currentLine.find("\t");
				std::string preStr = currentLine.substr(0, delim);
				std::string postStr = currentLine.substr(delim+1, currentLine.size()-delim-1);
				unsigned int preIndex = atoi(preStr.c_str());
				unsigned int postIndex = atoi(postStr.c_str());
				axonDendriteCorrespondence.insert(MatrixIndexType(preIndex, postIndex));
			}
		}
	}
	else
	{
		std::cout << "Error reading axon/dendrite correspondence file " << filename << std::endl;
	}
	
	return axonDendriteCorrespondence;
}

void writeUniqueMatrixRows( ConnectionMatrix* connectome, std::map< unsigned int, SelectionType > uniqueIDsPerCelltype, const char* outputFilename )
{
	SelectionType postsynapticCellIDs = getPostsynapticCellIDs(connectome);
	std::map< unsigned int, SelectionType >::const_iterator uniqueIDsPerCelltypeIt;
	
	std::string format = outputFilename;
	if(format.find(".csv") == std::string::npos)
	{
		format += ".csv";
	}
	std::ofstream OutputFile( format.c_str() );
	
	OutputFile << "CELLID\tCELLTYPE\tCOLUMN";
	for(uniqueIDsPerCelltypeIt = uniqueIDsPerCelltype.begin(); uniqueIDsPerCelltypeIt != uniqueIDsPerCelltype.end(); ++uniqueIDsPerCelltypeIt)
	{
		const char * cellType = int2CelltypeLabels[uniqueIDsPerCelltypeIt->first];
		SelectionType uniqueIDs = uniqueIDsPerCelltypeIt->second;
		for(int ii = 0; ii < uniqueIDs.size(); ++ii)
		{
			OutputFile << "\t" << cellType << "_" << uniqueIDs[ii];
		}
	}
	OutputFile << std::endl;
	
	for(int ii = 0; ii < postsynapticCellIDs.size(); ++ii)
	{
		unsigned int postCellID = postsynapticCellIDs[ii];
		const char * cellType = int2CelltypeLabels[connectome->IDColumnCelltypeMap[postCellID].second];
		const char * column = int2ColumnLabels[connectome->IDColumnCelltypeMap[postCellID].first];
		OutputFile << postCellID << "\t" << cellType << "\t" << column;
		for(uniqueIDsPerCelltypeIt = uniqueIDsPerCelltype.begin(); uniqueIDsPerCelltypeIt != uniqueIDsPerCelltype.end(); ++uniqueIDsPerCelltypeIt)
		{
			SelectionType uniqueIDs = uniqueIDsPerCelltypeIt->second;
			for(int jj = 0; jj < uniqueIDs.size(); ++jj)
			{
				unsigned int preCellID = uniqueIDs[jj];
				float innervation = 0;
				MatrixIndexType innervationIndex(preCellID, postCellID);
				if(connectome->matrix.find(innervationIndex) != connectome->matrix.end())
				{
					innervation = connectome->matrix[innervationIndex];
				}
				OutputFile << "\t" << innervation;
			}
		}
		OutputFile << std::endl;
	}
	OutputFile.close();
}

void writeCompleteMatrixRows( ConnectionMatrix* connectome, std::map< unsigned int, SelectionType > presynapticCellIDsPerCellType, SelectionType postsynapticCellIDs,
							  std::map< unsigned int, unsigned int > axonDendriteCorrespondence, const char* outputFilename )
{
	std::string format = outputFilename;
	if(format.find(".csv") == std::string::npos)
	{
		format += ".csv";
	}
	std::ofstream OutputFile( format.c_str() );
	
	OutputFile << "CELLID\tCELLTYPE\tCOLUMN";
	std::map< unsigned int, SelectionType >::const_iterator presynapticCellIDsPerCellTypeIt;
	for(presynapticCellIDsPerCellTypeIt = presynapticCellIDsPerCellType.begin();
		presynapticCellIDsPerCellTypeIt != presynapticCellIDsPerCellType.end(); ++presynapticCellIDsPerCellTypeIt)
	{
		const char * cellType = int2CelltypeLabels[presynapticCellIDsPerCellTypeIt->first];
		SelectionType presynapticCellIDs = presynapticCellIDsPerCellTypeIt->second;
		for(int ii = 0; ii < presynapticCellIDs.size(); ++ii)
		{
			OutputFile << "\t" << cellType << "_" << presynapticCellIDs[ii];
		}
	}
	OutputFile << std::endl;
	
	for(int ii = 0; ii < postsynapticCellIDs.size(); ++ii)
	{
		unsigned int postCellID = postsynapticCellIDs[ii];
		const char * cellType = int2CelltypeLabels[connectome->IDColumnCelltypeMap[postCellID].second];
		const char * column = int2ColumnLabels[connectome->IDColumnCelltypeMap[postCellID].first];
		OutputFile << postCellID << "\t" << cellType << "\t" << column;
		for(presynapticCellIDsPerCellTypeIt = presynapticCellIDsPerCellType.begin();
			presynapticCellIDsPerCellTypeIt != presynapticCellIDsPerCellType.end(); ++presynapticCellIDsPerCellTypeIt)
		{
			SelectionType presynapticCellIDs = presynapticCellIDsPerCellTypeIt->second;
			for(int jj = 0; jj < presynapticCellIDs.size(); ++jj)
			{
				unsigned int preCellID = presynapticCellIDs[jj];
				float innervation = 0;
				MatrixIndexType innervationIndex(preCellID, postCellID);
				if(connectome->matrix.find(innervationIndex) != connectome->matrix.end())
				{
					innervation = connectome->matrix[innervationIndex];
				}
				// set diagonals to 0 (no self-innervation)
				if(axonDendriteCorrespondence.find(preCellID) != axonDendriteCorrespondence.end())
					if(axonDendriteCorrespondence[preCellID] == postCellID)
					{
						innervation = 0;
					}
				OutputFile << "\t" << innervation;
			}
		}
		OutputFile << std::endl;
	}
	OutputFile.close();
}

void writeConnectivitySummary(ConnectionMatrix* connectome, std::map< unsigned int, SelectionType > presynapticCellIDsPerCellType, SelectionType postsynapticCellIDs,
							  std::map< unsigned int, unsigned int > axonDendriteCorrespondence, const char* outputFilename)
{
	std::list< unsigned int > excPreCellTypes;
	std::list< unsigned int > inhPreCellTypes;
	excPreCellTypes.push_back(VPM);
	excPreCellTypes.push_back(L2axon);
	excPreCellTypes.push_back(L34axon);
	excPreCellTypes.push_back(L4pyaxon);
	excPreCellTypes.push_back(L4spaxon);
	excPreCellTypes.push_back(L4ssaxon);
	excPreCellTypes.push_back(L5staxon);
	excPreCellTypes.push_back(L5ttaxon);
	excPreCellTypes.push_back(L6ccaxon);
	excPreCellTypes.push_back(L6ccinvaxon);
	excPreCellTypes.push_back(L6ctaxon);
	inhPreCellTypes.push_back(SymLocalaxon);
	inhPreCellTypes.push_back(L23Transaxon);
	inhPreCellTypes.push_back(L45Peakaxon);
	inhPreCellTypes.push_back(L45Symaxon);
	inhPreCellTypes.push_back(L56Transaxon);
	std::map< unsigned int, SelectionType > postsynapticCellIDsPerCellType;
	for(int i = 0; i < postsynapticCellIDs.size(); ++i)
	{
		unsigned int postID = postsynapticCellIDs[i];
		unsigned int cellType = connectome->IDColumnCelltypeMap[postID].second;
		if(cellType >= SymLocal1 && cellType <= SymLocal6)
		{
			cellType = SymLocal;
		}
		if(postsynapticCellIDsPerCellType.find(cellType) != postsynapticCellIDsPerCellType.end())
		{
			postsynapticCellIDsPerCellType[cellType].push_back(postID);
		}
		else
		{
			SelectionType newIDVec;
			newIDVec.push_back(postID);
			postsynapticCellIDsPerCellType[cellType] = newIDVec;
		}
	}
	
	std::string format = outputFilename;
	if(format.find(".csv") == std::string::npos)
	{
		format += ".csv";
	}
	std::ofstream OutputFile( format.c_str() );
	OutputFile << "Cell type\tquantity\tavg\tSTD\tmin\tmax" << std::endl;
	std::map< unsigned int, SelectionType >::const_iterator postsynapticCellIDsPerCellTypeIt;
	for(postsynapticCellIDsPerCellTypeIt = postsynapticCellIDsPerCellType.begin();
		postsynapticCellIDsPerCellTypeIt != postsynapticCellIDsPerCellType.end(); ++postsynapticCellIDsPerCellTypeIt)
	{
		unsigned int cellType = postsynapticCellIDsPerCellTypeIt->first;
		SelectionType postsynapticCellIDs = postsynapticCellIDsPerCellTypeIt->second;
		
		unsigned int nrOfCells = postsynapticCellIDs.size();
		std::vector< float > totalInnervationPerCell;
		std::vector< float > excInnervationPerCell;
		std::vector< float > inhInnervationPerCell;
		std::vector< float > EIRatioPerCell;
		
		float avgInnervation = 0, stdInnervation = 0, minInnervation = 0, maxInnervation = 0;
		float avgExcInnervation = 0, stdExcInnervation = 0, minExcInnervation = 0, maxExcInnervation = 0;
		float avgInhInnervation = 0, stdInhInnervation = 0, minInhInnervation = 0, maxInhInnervation = 0;
		float avgEIRatio = 0, stdEIRatio = 0, minEIRatio = 0, maxEIRatio = 0;
		for(int ii = 0; ii < nrOfCells; ++ii)
		{
			float singleCellInnervation = 0, singleCellExcInnervation = 0, singleCellInhInnervation = 0, singleCellEIRatio = -1;
			unsigned int postID = postsynapticCellIDs[ii];
			std::list< unsigned int >::const_iterator preCellTypesIt;
			for(preCellTypesIt = excPreCellTypes.begin(); preCellTypesIt != excPreCellTypes.end(); ++preCellTypesIt)
			{
				unsigned int preCellType = *preCellTypesIt;
				SelectionType presynapticCellIDs = presynapticCellIDsPerCellType[preCellType];
				for(int jj = 0; jj < presynapticCellIDs.size(); ++jj)
				{
					unsigned int preID = presynapticCellIDs[jj];
					MatrixIndexType innervationIndex(preID, postID);
					if(connectome->matrix.find(innervationIndex) != connectome->matrix.end())
					{
						float innervation = connectome->matrix[innervationIndex];
						// set diagonals to 0 (no self-innervation)
						if(axonDendriteCorrespondence.find(preID) != axonDendriteCorrespondence.end())
							if(axonDendriteCorrespondence[preID] == postID)
							{
								innervation = 0;
							}
						singleCellInnervation += innervation;
						singleCellExcInnervation += innervation;
					}
				}
			}
			for(preCellTypesIt = inhPreCellTypes.begin(); preCellTypesIt != inhPreCellTypes.end(); ++preCellTypesIt)
			{
				unsigned int preCellType = *preCellTypesIt;
				SelectionType presynapticCellIDs = presynapticCellIDsPerCellType[preCellType];
				for(int jj = 0; jj < presynapticCellIDs.size(); ++jj)
				{
					unsigned int preID = presynapticCellIDs[jj];
					MatrixIndexType innervationIndex(preID, postID);
					if(connectome->matrix.find(innervationIndex) != connectome->matrix.end())
					{
						float innervation = connectome->matrix[innervationIndex];
						// set diagonals to 0 (no self-innervation)
						if(axonDendriteCorrespondence.find(preID) != axonDendriteCorrespondence.end())
							if(axonDendriteCorrespondence[preID] == postID)
							{
								innervation = 0;
							}
						singleCellInnervation += innervation;
						singleCellInhInnervation += innervation;
					}
				}
			}
			singleCellEIRatio = singleCellExcInnervation/singleCellInhInnervation;
			avgInnervation += singleCellInnervation;
			avgExcInnervation += singleCellExcInnervation;
			avgInhInnervation += singleCellInhInnervation;
			avgEIRatio += singleCellEIRatio;
			totalInnervationPerCell.push_back(singleCellInnervation);
			excInnervationPerCell.push_back(singleCellExcInnervation);
			inhInnervationPerCell.push_back(singleCellInhInnervation);
			EIRatioPerCell.push_back(singleCellEIRatio);
		}
		
		avgInnervation = avgInnervation/nrOfCells;
		avgExcInnervation = avgExcInnervation/nrOfCells;
		avgInhInnervation = avgInhInnervation/nrOfCells;
		avgEIRatio = avgEIRatio/nrOfCells;
		
		for(int ii = 0; ii < totalInnervationPerCell.size(); ++ii)
		{
			stdInnervation += (avgInnervation - totalInnervationPerCell[ii])*(avgInnervation - totalInnervationPerCell[ii]);
			stdExcInnervation += (avgExcInnervation - excInnervationPerCell[ii])*(avgExcInnervation - excInnervationPerCell[ii]);
			stdInhInnervation += (avgInhInnervation - inhInnervationPerCell[ii])*(avgInhInnervation - inhInnervationPerCell[ii]);
			stdEIRatio += (avgEIRatio - EIRatioPerCell[ii])*(avgEIRatio - EIRatioPerCell[ii]);
		}
		stdInnervation = sqrt(stdInnervation/(nrOfCells - 1));
		stdExcInnervation = sqrt(stdExcInnervation/(nrOfCells - 1));
		stdInhInnervation = sqrt(stdInhInnervation/(nrOfCells - 1));
		stdEIRatio = sqrt(stdEIRatio/(nrOfCells - 1));
		
		maxInnervation = *std::max_element(totalInnervationPerCell.begin(), totalInnervationPerCell.end());
		maxExcInnervation = *std::max_element(excInnervationPerCell.begin(), excInnervationPerCell.end());
		maxInhInnervation = *std::max_element(inhInnervationPerCell.begin(), inhInnervationPerCell.end());
		maxEIRatio = *std::max_element(EIRatioPerCell.begin(), EIRatioPerCell.end());
		
		minInnervation = *std::min_element(totalInnervationPerCell.begin(), totalInnervationPerCell.end());
		minExcInnervation = *std::min_element(excInnervationPerCell.begin(), excInnervationPerCell.end());
		minInhInnervation = *std::min_element(inhInnervationPerCell.begin(), inhInnervationPerCell.end());
		minEIRatio = *std::min_element(EIRatioPerCell.begin(), EIRatioPerCell.end());
		
		OutputFile << int2CelltypeLabels[cellType] << "\tNumber\t" << nrOfCells << std::endl;
		OutputFile << int2CelltypeLabels[cellType] << "\tTotal input\t" << avgInnervation << "\t" << stdInnervation << "\t" << minInnervation << "\t" << maxInnervation << std::endl;
		OutputFile << int2CelltypeLabels[cellType] << "\tEx input\t" << avgExcInnervation << "\t" << stdExcInnervation << "\t" << minExcInnervation << "\t" << maxExcInnervation << std::endl;
		OutputFile << int2CelltypeLabels[cellType] << "\tInh input\t" << avgInhInnervation << "\t" << stdInhInnervation << "\t" << minInhInnervation << "\t" << maxInhInnervation << std::endl;
		OutputFile << int2CelltypeLabels[cellType] << "\tE-I ratio\t" << avgEIRatio << "\t" << stdEIRatio << "\t" << minEIRatio << "\t" << maxEIRatio << std::endl;
	}
	OutputFile.close();
}

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
}

