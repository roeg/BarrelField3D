/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      matrixanalyzer.cpp                                            */
/*                                                                          */
/* Purpose:   Implementation of matrixanalyzer.h                            */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max Planck Institute for Biological Cybernetics               */
/*            Spemannstr. 38-44                                             */
/*            72076 Tuebingen                                               */
/*            Germany                                                       */
/*                                                                          */
/* EMail:     robert.egger@tuebingen.mpg.de                                 */
/*                                                                          */
/* History:   17.06.2014                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#include "matrixanalyzer.h"

TripletMotif::TripletMotif(std::list< ConnectionType > connections)
{
	size = 3;
	for(unsigned int ii = 0; ii < size; ++ii)
	{
		for(unsigned int jj = 0; jj < size; ++jj)
		{
			this->connections[ii][jj] = false;
		}
	}
	std::list< ConnectionType >::const_iterator connectionIt;
	for(connectionIt = connections.begin(); connectionIt != connections.end(); ++connectionIt)
	{
		this->connections[connectionIt->first][connectionIt->second] = true;
	}
}

double TripletMotif::computeOccurrenceProbability(std::vector< std::vector< double > > innervation)
{
	if(!innervation.size())
	{
		throw std::out_of_range("Innervation matrix empty");
	}
	else if(innervation.size() != size || innervation[0].size() != size)
	{
		throw std::out_of_range("Motif and innervation matrix dimensions not matching");
	}
	
	double probability = 1;
	for(unsigned int ii = 0; ii < size; ++ii)
	{
		for(unsigned int jj = 0; jj < size; ++jj)
		{
			if(ii == jj)
			{
				continue;
			}
			double prob = 1 - exp(-1*innervation[ii][jj]);
			if(connections[ii][jj])
			{
				probability *= prob;
			}
			else
			{
				probability *= (1 - prob);
			}
		}
	}
	return probability;
}

unsigned int TripletMotif::computeOccurrencesWithStrength(std::vector< std::vector< unsigned int > > synapseRangeMatrix)
{
	if(!synapseRangeMatrix.size())
	{
		throw std::out_of_range("Synapse range matrix empty");
	}
	else if(synapseRangeMatrix.size() != size || synapseRangeMatrix[0].size() != size)
	{
		throw std::out_of_range("Motif and synapse range matrix dimensions not matching");
	}
	
	unsigned int totalNumber = 1;
	for(unsigned int ii = 0; ii < size; ++ii)
	{
		for(unsigned int jj = 0; jj < size; ++jj)
		{
			if(ii == jj)
			{
				continue;
			}
			if(connections[ii][jj])
			{
				totalNumber *= synapseRangeMatrix[ii][jj];
			}
		}
	}
	return totalNumber;
}


CellTriplet::CellTriplet(unsigned int preCell1, unsigned int preCell2, unsigned int preCell3, unsigned int postCell1, unsigned int postCell2, unsigned int postCell3)
{
	preCellIndex.push_back(preCell1);
	preCellIndex.push_back(preCell2);
	preCellIndex.push_back(preCell3);
	postCellIndex.push_back(postCell1);
	postCellIndex.push_back(postCell2);
	postCellIndex.push_back(postCell3);
}

void CellTriplet::setInnervationMatrix(ConnectionMatrix* connectome)
{
	for(int ii = 0; ii < 3; ++ii)
	{
		std::vector< double > emptyRow;
		this->innervation.push_back(emptyRow);
		for(int jj = 0; jj < 3; ++jj)
		{
			if(ii == jj)
			{
				this->innervation[ii].push_back(0);
				continue;
			}
			MatrixIndexType index(preCellIndex[ii], postCellIndex[jj]);
			double innervation = 0;
			if(connectome->matrix.find(index) != connectome->matrix.end())
			{
				innervation = connectome->matrix[index];
			}
			this->innervation[ii].push_back(innervation);
		}
	}
}


MatrixAnalyzer::MatrixAnalyzer()
{
	initializeConstants();
	this->connectome = NULL;
}

MatrixAnalyzer::MatrixAnalyzer(ConnectionMatrix* connectome)
{
	initializeConstants();
	this->connectome = connectome;
}

MatrixAnalyzer::~MatrixAnalyzer()
{

}

Profile* MatrixAnalyzer::computeInnervationHistogram(unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType, float binSize)
{
	Profile * hist = new Profile(binSize);
	SelectionType preSelection = connectome->getPreColumnCelltypeSelection(preColumn, preType);
	SelectionType postSelection = connectome->getPostColumnCelltypeSelection(postColumn, postType);
	
	std::cout << "Computing innervation histogram of " << preSelection.size()*postSelection.size() << " connections..." << std::endl;
	std::cout << "Pre: " << int2CelltypeLabels[preType] << " " << int2Labels[preColumn] << std::endl;
	std::cout << "Post: " << int2CelltypeLabels[postType] << " " << int2Labels[postColumn] << std::endl;
	
	for(int ii = 0; ii < preSelection.size(); ++ii)
	{
		for(int jj = 0; jj < postSelection.size(); ++jj)
		{
			MatrixIndexType index(preSelection[ii], postSelection[jj]);
			std::map< MatrixIndexType, float >::const_iterator matrixIt = connectome->matrix.find(index);
			if(matrixIt != connectome->matrix.end())
			{
				float innervation = connectome->matrix[index];
				unsigned int bin = (unsigned int)(innervation/binSize);
				hist->incrementBin(bin);
			}
			else
			{
				hist->incrementBin(0);
			}
		}
	}
	
	return hist;
}

Profile* MatrixAnalyzer::computeProbabilityHistogram(unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType, float binSize)
{
	Profile * hist = new Profile(binSize);
	SelectionType preSelection = connectome->getPreColumnCelltypeSelection(preColumn, preType);
	SelectionType postSelection = connectome->getPostColumnCelltypeSelection(postColumn, postType);
	
	std::cout << "Computing probability histogram of " << preSelection.size()*postSelection.size() << " connections..." << std::endl;
	std::cout << "Pre: " << int2CelltypeLabels[preType] << " " << int2Labels[preColumn] << std::endl;
	std::cout << "Post: " << int2CelltypeLabels[postType] << " " << int2Labels[postColumn] << std::endl;
	
	for(int ii = 0; ii < preSelection.size(); ++ii)
	{
		for(int jj = 0; jj < postSelection.size(); ++jj)
		{
			MatrixIndexType index(preSelection[ii], postSelection[jj]);
			std::map< MatrixIndexType, float >::const_iterator matrixIt = connectome->matrix.find(index);
			if(matrixIt != connectome->matrix.end())
			{
				float innervation = connectome->matrix[index];
				float prob = 1 - exp(-1*innervation);
				unsigned int bin = (unsigned int)(prob/binSize);
				hist->incrementBin(bin);
			}
			else
			{
				hist->incrementBin(0);
			}
		}
	}
	
	return hist;
}

Profile* MatrixAnalyzer::computeSynapseNumberHistogram(unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType)
{
	const double percentile = 0.99;
	Profile * hist = new Profile(1);
	SelectionType preSelection = connectome->getPreColumnCelltypeSelection(preColumn, preType);
	SelectionType postSelection = connectome->getPostColumnCelltypeSelection(postColumn, postType);
	
	std::cout << "Computing synapse number histogram of " << preSelection.size()*postSelection.size() << " connections..." << std::endl;
	std::cout << "Pre: " << int2CelltypeLabels[preType] << " " << int2Labels[preColumn] << std::endl;
	std::cout << "Post: " << int2CelltypeLabels[postType] << " " << int2Labels[postColumn] << std::endl;
	
	for(int ii = 0; ii < preSelection.size(); ++ii)
	{
		for(int jj = 0; jj < postSelection.size(); ++jj)
		{
			MatrixIndexType index(preSelection[ii], postSelection[jj]);
			std::map< MatrixIndexType, float >::const_iterator matrixIt = connectome->matrix.find(index);
			if(matrixIt != connectome->matrix.end())
			{
				float innervation = connectome->matrix[index];
				Profile * poissonHist = computePercentilePoissonHistogram(innervation, percentile);
				hist->addProfile(poissonHist);
				delete poissonHist;
			}
			else
			{
				hist->incrementBin(0);
			}
		}
	}
	
	return hist;
}

void MatrixAnalyzer::analyzeTripletMotifs(const char* outputFilename, unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType, unsigned int nrOfTriplets)
{
	std::vector< Profile* > tripletMotifVec = computeTripletMotifDistribution(preColumn, preType, postColumn, postType, nrOfTriplets);
	Profile * finalTriplet = tripletMotifVec[nrOfTriplets-1];
	Profile * tripletMotifConvergence = new Profile(1);
	for(int ii = 0; ii < nrOfTriplets-1; ++ii)
	{
		Profile * intermediateTriplet = tripletMotifVec[ii];
		double rmse = computeHistogramRMSE(intermediateTriplet, finalTriplet);
		tripletMotifConvergence->addSegment(rmse, ii);
	}
	
	std::string motifOutName(outputFilename);
	motifOutName += "_motifHist.csv";
	finalTriplet->writeProfile(motifOutName.c_str());
	std::string convergenceOutName(outputFilename);
	convergenceOutName += "_motifConvergence.csv";
	tripletMotifConvergence->writeProfile(convergenceOutName.c_str(), 1);
	
	delete tripletMotifConvergence;
	for(int ii = 0; ii < tripletMotifVec.size(); ++ii)
	{
		delete tripletMotifVec[ii];
	}
}

std::vector< Profile* > MatrixAnalyzer::computeTripletMotifDistribution(unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType, unsigned int nrOfTriplets)
{
	std::cout << "Analyzing triplet motif distribution for " << nrOfTriplets << " triplets..." << std::endl;
	std::map< unsigned int, std::list< TripletMotif * > > tripletMotifs = initializeNonRedundantTripletMotifs();
	SelectionType preSelection = connectome->getPreColumnCelltypeSelection(preColumn, preType);
	SelectionType postSelection = connectome->getPostColumnCelltypeSelection(postColumn, postType);
	std::list< CellTriplet * > triplets = initializeNonRedundantCellTriplets(preSelection, postSelection, nrOfTriplets);
	
// 	// simple test case:
// 	// hard-code innervation of manually selected cell triplet,
// 	// i.e. the one from the python code
// 	CellTriplet * sampleTriplet = new CellTriplet(0,0,0,0,0,0);
// 	std::vector< double > emptyRow1;
// 	sampleTriplet->innervation.push_back(emptyRow1);
// 	sampleTriplet->innervation[0].push_back(0.0);
// 	sampleTriplet->innervation[0].push_back(0.15);
// 	sampleTriplet->innervation[0].push_back(1.17);
// 	std::vector< double > emptyRow2;
// 	sampleTriplet->innervation.push_back(emptyRow2);
// 	sampleTriplet->innervation[1].push_back(1.14);
// 	sampleTriplet->innervation[1].push_back(0.0);
// 	sampleTriplet->innervation[1].push_back(1.78);
// 	std::vector< double > emptyRow3;
// 	sampleTriplet->innervation.push_back(emptyRow3);
// 	sampleTriplet->innervation[2].push_back(0.68);
// 	sampleTriplet->innervation[2].push_back(0.63);
// 	sampleTriplet->innervation[2].push_back(0.0);
// 	std::list< CellTriplet * > triplets;
// 	triplets.push_back(sampleTriplet);
	
	std::cout << "Computing probability of occurrence for " << triplets.size() << " triplets..." << std::endl;
	std::vector< Profile * > tripletMotifVec;
	Profile * motifHistogram = new Profile(1);
	std::list< CellTriplet * >::const_iterator tripletsIt;
	for(tripletsIt = triplets.begin(); tripletsIt != triplets.end(); ++tripletsIt)
	{
		CellTriplet * currentTriplet = *tripletsIt;
		for(int ii = 0; ii < tripletMotifs.size(); ++ii)
		{
			double motifProb = 0;
			std::list< TripletMotif * > motifList = tripletMotifs[ii];
			std::list< TripletMotif * >::const_iterator motifListIt;
			for(motifListIt = motifList.begin(); motifListIt != motifList.end(); ++motifListIt)
			{
				TripletMotif * currentMotif = *motifListIt;
				motifProb += currentMotif->computeOccurrenceProbability(currentTriplet->innervation);
			}
			motifHistogram->addSegment(motifProb, ii);
		}
		Profile * intermediateHist = new Profile(1);
		intermediateHist->addProfile(motifHistogram);
		tripletMotifVec.push_back(intermediateHist);
	}
	
	// normalize
	for(int ii = 0; ii < tripletMotifVec.size(); ++ii)
	{
		for(int jj = 0; jj < tripletMotifs.size(); ++jj)
		{
			Profile * tmpHist = tripletMotifVec[ii];
			double histValue = tmpHist->getProfile()->at(jj);
			histValue /= (double)(ii+1);
			tmpHist->getProfile()->at(jj) = histValue;
		}
	}
	
	for(int ii = 0; ii < tripletMotifs.size(); ++ii)
	{
		std::list< TripletMotif * > motifList = tripletMotifs[ii];
		std::list< TripletMotif * >::iterator motifListIt;
		for(motifListIt = motifList.begin(); motifListIt != motifList.end(); ++motifListIt)
		{
			delete *motifListIt;
		}
	}
	std::list< CellTriplet * >::iterator tripletsCleanIt;
	for(tripletsCleanIt = triplets.begin(); tripletsCleanIt != triplets.end(); ++tripletsCleanIt)
	{
		delete *tripletsCleanIt;
	}
	
	return tripletMotifVec;
}

void MatrixAnalyzer::writeConnectionMatrix(const char* outputFilename, unsigned int preColumn, unsigned int preType, unsigned int postColumn, unsigned int postType)
{
	SelectionType preSelection = connectome->getPreColumnCelltypeSelection(preColumn, preType);
	SelectionType postSelection = connectome->getPostColumnCelltypeSelection(postColumn, postType);
	
	std::cout << "Writing connection matrix with " << preSelection.size()*postSelection.size() << " entries as image..." << std::endl;
	std::cout << "Pre: " << int2CelltypeLabels[preType] << " " << int2Labels[preColumn] << std::endl;
	std::cout << "Post: " << int2CelltypeLabels[postType] << " " << int2Labels[postColumn] << std::endl;
	
	connectome->writeConnectionMatrixAsImage(outputFilename, preSelection, postSelection);
}

std::vector< unsigned int > MatrixAnalyzer::parseInputParameters(int argc, char* argv[], int argOffset)
{
	std::vector< unsigned int > parameters;
	for(int ii = argOffset; ii < argc; ++ii)
	{
		bool foundParameter = false;
		std::string paramString(argv[ii]);
		std::map< std::string, unsigned int >::const_iterator labelIt;
		for(labelIt = labels2Int.begin(); labelIt != labels2Int.end(); ++labelIt)
		{
			if(!foundParameter && labelIt->first.compare(paramString) == 0)
			{
				foundParameter = true;
				parameters.push_back(labelIt->second);
			}
		}
		if(foundParameter)
		{
			continue;
		}
		for(labelIt = celltypeLabels2Int.begin(); labelIt != celltypeLabels2Int.end(); ++labelIt)
		{
			if(!foundParameter && labelIt->first.compare(paramString) == 0)
			{
				foundParameter = true;
				parameters.push_back(labelIt->second);
			}
		}
		if(!foundParameter)
		{
			std::cout << "Error! Unknown parameter: " << argv[ii] << std::endl;
			break;
		}
	}
	
	return parameters;
}

void MatrixAnalyzer::writeSynapsesPerCellRows(const char* outputFilename, std::list< unsigned int > postColumns, std::list< unsigned int > postTypes, CellTable* table)
{
	std::vector< CellTableRow * > rows = table->getPostColumnCelltypeRows(postColumns, postTypes);
	std::cout << "Selected " << rows.size() << " rows..." << std::endl;
	std::string outputString(outputFilename);
	if(outputString.find(".csv") == std::string::npos)
	{
		outputString += ".csv";
	}
	std::ofstream OutputFile(outputString.c_str());
	OutputFile << "CELLID\tCELLTYPE\tCOLUMN\tSOMA_X\tSOMA_Y\tSOMA_Z\tINSIDE_COLUMN\tTOTAL";
	std::map< ColumnCellTypePair, unsigned int >::const_iterator headerIt;
	for(headerIt = table->header.begin(); headerIt != table->header.end(); ++headerIt)
	{
		OutputFile << "\t";
		OutputFile << int2Labels[headerIt->first.first];
		OutputFile << "_";
		OutputFile << int2CelltypeLabels[headerIt->first.second];
	}
	OutputFile << std::endl;
	
	for(int i = 0; i < rows.size(); ++i)
	{
		CellTableRow * row = rows[i];
		OutputFile << row->cellID << "\t";
		OutputFile << int2CelltypeLabels[row->cellType] << "\t";
		OutputFile << int2Labels[row->column] << "\t";
		OutputFile << row->somaLocation[0] << "\t";
		OutputFile << row->somaLocation[1] << "\t";
		OutputFile << row->somaLocation[2] << "\t";
		OutputFile << row->insideColumn << "\t";
		OutputFile << row->totalSynapses;
		for(headerIt = table->header.begin(); headerIt != table->header.end(); ++headerIt)
		{
			OutputFile << "\t";
			OutputFile << row->synapsesPerPreTypeColumn[headerIt->second];
		}
		OutputFile << std::endl;
	}
	OutputFile.close();
}

Profile* MatrixAnalyzer::computePercentilePoissonHistogram(float innervation, double percentile)
{
	if(percentile < 0 || percentile > 1)
	{
		throw std::out_of_range("Percentile out of range [0,1]");
	}
	Profile * hist = new Profile(1);
	double targetConnectionProb = 1 - exp(-1*innervation);
	double cumulativeConnectionProb = 0;
	unsigned int nSyn = 1;
	hist->addSegment(1-targetConnectionProb, 0);
	while(cumulativeConnectionProb < percentile*targetConnectionProb)
	{
		double binProb = pow(innervation, nSyn)/vtkMath::Factorial(nSyn)*exp(-1*innervation);
		hist->addSegment(binProb, nSyn);
		cumulativeConnectionProb += binProb;
		++nSyn;
	}
	
	return hist;
}

std::map< unsigned int, std::list< TripletMotif* > > MatrixAnalyzer::initializeNonRedundantTripletMotifs()
{
	std::list< std::list< ConnectionType > >::const_iterator connectionIt;
	
	// motifs for network analysis: three edges
	std::list< std::list< ConnectionType > > recurrentLoopConnections;
	std::list< std::list< ConnectionType > > directedLoopConnections;
	std::list< std::list< ConnectionType > > recurrentIncompleteLoopConnections;
	std::list< std::list< ConnectionType > > directedRecurrentLoopConnections;
	std::list< std::list< ConnectionType > > recurrentFeedForwardConvergentConnections;
	std::list< std::list< ConnectionType > > recurrentFeedForwardDivergentConnections;
	std::list< std::list< ConnectionType > > feedForwardConnections;
	std::list< TripletMotif * > recurrentLoopMotifs;
	std::list< TripletMotif * > directedLoopMotifs;
	std::list< TripletMotif * > recurrentIncompleteLoopMotifs;
	std::list< TripletMotif * > directedRecurrentLoopMotifs;
	std::list< TripletMotif * > recurrentFeedForwardConvergentMotifs;
	std::list< TripletMotif * > recurrentFeedForwardDivergentMotifs;
	std::list< TripletMotif * > feedForwardMotifs;
	
	std::list< ConnectionType > recurrentLoopConnections1;
	recurrentLoopConnections1.push_back(ConnectionType(0,1));
	recurrentLoopConnections1.push_back(ConnectionType(1,0));
	recurrentLoopConnections1.push_back(ConnectionType(1,2));
	recurrentLoopConnections1.push_back(ConnectionType(2,1));
	recurrentLoopConnections1.push_back(ConnectionType(0,2));
	recurrentLoopConnections1.push_back(ConnectionType(2,0));
	recurrentLoopConnections.push_back(recurrentLoopConnections1);
	std::list< ConnectionType > directedLoopConnections1;
	directedLoopConnections1.push_back(ConnectionType(0,1));
	directedLoopConnections1.push_back(ConnectionType(1,2));
	directedLoopConnections1.push_back(ConnectionType(2,0));
	std::list< ConnectionType > directedLoopConnections2;
	directedLoopConnections2.push_back(ConnectionType(1,0));
	directedLoopConnections2.push_back(ConnectionType(0,2));
	directedLoopConnections2.push_back(ConnectionType(2,1));
	directedLoopConnections.push_back(directedLoopConnections1);
	directedLoopConnections.push_back(directedLoopConnections2);
	std::list< ConnectionType > recurrentIncompleteLoopConnection1;
	recurrentIncompleteLoopConnection1.push_back(ConnectionType(0,1));
	recurrentIncompleteLoopConnection1.push_back(ConnectionType(1,2));
	recurrentIncompleteLoopConnection1.push_back(ConnectionType(2,1));
	recurrentIncompleteLoopConnection1.push_back(ConnectionType(0,2));
	recurrentIncompleteLoopConnection1.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentIncompleteLoopConnection2;
	recurrentIncompleteLoopConnection2.push_back(ConnectionType(1,0));
	recurrentIncompleteLoopConnection2.push_back(ConnectionType(1,2));
	recurrentIncompleteLoopConnection2.push_back(ConnectionType(2,1));
	recurrentIncompleteLoopConnection2.push_back(ConnectionType(0,2));
	recurrentIncompleteLoopConnection2.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentIncompleteLoopConnection3;
	recurrentIncompleteLoopConnection3.push_back(ConnectionType(0,1));
	recurrentIncompleteLoopConnection3.push_back(ConnectionType(1,0));
	recurrentIncompleteLoopConnection3.push_back(ConnectionType(1,2));
	recurrentIncompleteLoopConnection3.push_back(ConnectionType(0,2));
	recurrentIncompleteLoopConnection3.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentIncompleteLoopConnection4;
	recurrentIncompleteLoopConnection4.push_back(ConnectionType(0,1));
	recurrentIncompleteLoopConnection4.push_back(ConnectionType(1,0));
	recurrentIncompleteLoopConnection4.push_back(ConnectionType(2,1));
	recurrentIncompleteLoopConnection4.push_back(ConnectionType(0,2));
	recurrentIncompleteLoopConnection4.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentIncompleteLoopConnection5;
	recurrentIncompleteLoopConnection5.push_back(ConnectionType(0,1));
	recurrentIncompleteLoopConnection5.push_back(ConnectionType(1,0));
	recurrentIncompleteLoopConnection5.push_back(ConnectionType(1,2));
	recurrentIncompleteLoopConnection5.push_back(ConnectionType(2,1));
	recurrentIncompleteLoopConnection5.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentIncompleteLoopConnection6;
	recurrentIncompleteLoopConnection6.push_back(ConnectionType(0,1));
	recurrentIncompleteLoopConnection6.push_back(ConnectionType(1,0));
	recurrentIncompleteLoopConnection6.push_back(ConnectionType(1,2));
	recurrentIncompleteLoopConnection6.push_back(ConnectionType(2,1));
	recurrentIncompleteLoopConnection6.push_back(ConnectionType(0,2));
	recurrentIncompleteLoopConnections.push_back(recurrentIncompleteLoopConnection1);
	recurrentIncompleteLoopConnections.push_back(recurrentIncompleteLoopConnection2);
	recurrentIncompleteLoopConnections.push_back(recurrentIncompleteLoopConnection3);
	recurrentIncompleteLoopConnections.push_back(recurrentIncompleteLoopConnection4);
	recurrentIncompleteLoopConnections.push_back(recurrentIncompleteLoopConnection5);
	recurrentIncompleteLoopConnections.push_back(recurrentIncompleteLoopConnection6);
	std::list< ConnectionType > directedRecurrentLoopConnections1;
	directedRecurrentLoopConnections1.push_back(ConnectionType(0,1));
	directedRecurrentLoopConnections1.push_back(ConnectionType(1,2));
	directedRecurrentLoopConnections1.push_back(ConnectionType(0,2));
	directedRecurrentLoopConnections1.push_back(ConnectionType(2,0));
	std::list< ConnectionType > directedRecurrentLoopConnections2;
	directedRecurrentLoopConnections2.push_back(ConnectionType(1,0));
	directedRecurrentLoopConnections2.push_back(ConnectionType(2,1));
	directedRecurrentLoopConnections2.push_back(ConnectionType(0,2));
	directedRecurrentLoopConnections2.push_back(ConnectionType(2,0));
	std::list< ConnectionType > directedRecurrentLoopConnections3;
	directedRecurrentLoopConnections3.push_back(ConnectionType(0,1));
	directedRecurrentLoopConnections3.push_back(ConnectionType(1,0));
	directedRecurrentLoopConnections3.push_back(ConnectionType(1,2));
	directedRecurrentLoopConnections3.push_back(ConnectionType(2,0));
	std::list< ConnectionType > directedRecurrentLoopConnections4;
	directedRecurrentLoopConnections4.push_back(ConnectionType(0,1));
	directedRecurrentLoopConnections4.push_back(ConnectionType(1,0));
	directedRecurrentLoopConnections4.push_back(ConnectionType(2,1));
	directedRecurrentLoopConnections4.push_back(ConnectionType(0,2));
	std::list< ConnectionType > directedRecurrentLoopConnections5;
	directedRecurrentLoopConnections5.push_back(ConnectionType(0,1));
	directedRecurrentLoopConnections5.push_back(ConnectionType(1,2));
	directedRecurrentLoopConnections5.push_back(ConnectionType(2,1));
	directedRecurrentLoopConnections5.push_back(ConnectionType(2,0));
	std::list< ConnectionType > directedRecurrentLoopConnections6;
	directedRecurrentLoopConnections6.push_back(ConnectionType(1,0));
	directedRecurrentLoopConnections6.push_back(ConnectionType(1,2));
	directedRecurrentLoopConnections6.push_back(ConnectionType(2,1));
	directedRecurrentLoopConnections6.push_back(ConnectionType(0,2));
	directedRecurrentLoopConnections.push_back(directedRecurrentLoopConnections1);
	directedRecurrentLoopConnections.push_back(directedRecurrentLoopConnections2);
	directedRecurrentLoopConnections.push_back(directedRecurrentLoopConnections3);
	directedRecurrentLoopConnections.push_back(directedRecurrentLoopConnections4);
	directedRecurrentLoopConnections.push_back(directedRecurrentLoopConnections5);
	directedRecurrentLoopConnections.push_back(directedRecurrentLoopConnections6);
	std::list< ConnectionType > recurrentFeedForwardConvergentConnections1;
	recurrentFeedForwardConvergentConnections1.push_back(ConnectionType(1,0));
	recurrentFeedForwardConvergentConnections1.push_back(ConnectionType(1,2));
	recurrentFeedForwardConvergentConnections1.push_back(ConnectionType(2,1));
	recurrentFeedForwardConvergentConnections1.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentFeedForwardConvergentConnections2;
	recurrentFeedForwardConvergentConnections2.push_back(ConnectionType(0,1));
	recurrentFeedForwardConvergentConnections2.push_back(ConnectionType(2,1));
	recurrentFeedForwardConvergentConnections2.push_back(ConnectionType(0,2));
	recurrentFeedForwardConvergentConnections2.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentFeedForwardConvergentConnections3;
	recurrentFeedForwardConvergentConnections3.push_back(ConnectionType(0,2));
	recurrentFeedForwardConvergentConnections3.push_back(ConnectionType(1,2));
	recurrentFeedForwardConvergentConnections3.push_back(ConnectionType(0,1));
	recurrentFeedForwardConvergentConnections3.push_back(ConnectionType(1,0));
	recurrentFeedForwardConvergentConnections.push_back(recurrentFeedForwardConvergentConnections1);
	recurrentFeedForwardConvergentConnections.push_back(recurrentFeedForwardConvergentConnections2);
	recurrentFeedForwardConvergentConnections.push_back(recurrentFeedForwardConvergentConnections3);
	std::list< ConnectionType > recurrentFeedForwardDivergentConnections1;
	recurrentFeedForwardDivergentConnections1.push_back(ConnectionType(0,1));
	recurrentFeedForwardDivergentConnections1.push_back(ConnectionType(0,2));
	recurrentFeedForwardDivergentConnections1.push_back(ConnectionType(1,2));
	recurrentFeedForwardDivergentConnections1.push_back(ConnectionType(2,1));
	std::list< ConnectionType > recurrentFeedForwardDivergentConnections2;
	recurrentFeedForwardDivergentConnections2.push_back(ConnectionType(1,0));
	recurrentFeedForwardDivergentConnections2.push_back(ConnectionType(1,2));
	recurrentFeedForwardDivergentConnections2.push_back(ConnectionType(0,2));
	recurrentFeedForwardDivergentConnections2.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentFeedForwardDivergentConnections3;
	recurrentFeedForwardDivergentConnections3.push_back(ConnectionType(2,0));
	recurrentFeedForwardDivergentConnections3.push_back(ConnectionType(2,1));
	recurrentFeedForwardDivergentConnections3.push_back(ConnectionType(0,1));
	recurrentFeedForwardDivergentConnections3.push_back(ConnectionType(1,0));
	recurrentFeedForwardDivergentConnections.push_back(recurrentFeedForwardDivergentConnections1);
	recurrentFeedForwardDivergentConnections.push_back(recurrentFeedForwardDivergentConnections2);
	recurrentFeedForwardDivergentConnections.push_back(recurrentFeedForwardDivergentConnections3);
	std::list< ConnectionType > feedForwardConnections1;
	feedForwardConnections1.push_back(ConnectionType(0,1));
	feedForwardConnections1.push_back(ConnectionType(0,2));
	feedForwardConnections1.push_back(ConnectionType(1,2));
	std::list< ConnectionType > feedForwardConnections2;
	feedForwardConnections2.push_back(ConnectionType(0,1));
	feedForwardConnections2.push_back(ConnectionType(0,2));
	feedForwardConnections2.push_back(ConnectionType(2,1));
	std::list< ConnectionType > feedForwardConnections3;
	feedForwardConnections3.push_back(ConnectionType(1,0));
	feedForwardConnections3.push_back(ConnectionType(1,2));
	feedForwardConnections3.push_back(ConnectionType(0,2));
	std::list< ConnectionType > feedForwardConnections4;
	feedForwardConnections4.push_back(ConnectionType(1,0));
	feedForwardConnections4.push_back(ConnectionType(1,2));
	feedForwardConnections4.push_back(ConnectionType(2,0));
	std::list< ConnectionType > feedForwardConnections5;
	feedForwardConnections5.push_back(ConnectionType(2,0));
	feedForwardConnections5.push_back(ConnectionType(2,1));
	feedForwardConnections5.push_back(ConnectionType(0,1));
	std::list< ConnectionType > feedForwardConnections6;
	feedForwardConnections6.push_back(ConnectionType(2,0));
	feedForwardConnections6.push_back(ConnectionType(2,1));
	feedForwardConnections6.push_back(ConnectionType(1,0));
	feedForwardConnections.push_back(feedForwardConnections1);
	feedForwardConnections.push_back(feedForwardConnections2);
	feedForwardConnections.push_back(feedForwardConnections3);
	feedForwardConnections.push_back(feedForwardConnections4);
	feedForwardConnections.push_back(feedForwardConnections5);
	feedForwardConnections.push_back(feedForwardConnections6);
	
	for(connectionIt = recurrentLoopConnections.begin(); connectionIt != recurrentLoopConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		recurrentLoopMotifs.push_back(newMotif);
	}
	for(connectionIt = directedLoopConnections.begin(); connectionIt != directedLoopConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		directedLoopMotifs.push_back(newMotif);
	}
	for(connectionIt = recurrentIncompleteLoopConnections.begin(); connectionIt != recurrentIncompleteLoopConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		recurrentIncompleteLoopMotifs.push_back(newMotif);
	}
	for(connectionIt = directedRecurrentLoopConnections.begin(); connectionIt != directedRecurrentLoopConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		directedRecurrentLoopMotifs.push_back(newMotif);
	}
	for(connectionIt = recurrentFeedForwardConvergentConnections.begin(); connectionIt != recurrentFeedForwardConvergentConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		recurrentFeedForwardConvergentMotifs.push_back(newMotif);
	}
	for(connectionIt = recurrentFeedForwardDivergentConnections.begin(); connectionIt != recurrentFeedForwardDivergentConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		recurrentFeedForwardDivergentMotifs.push_back(newMotif);
	}
	for(connectionIt = feedForwardConnections.begin(); connectionIt != feedForwardConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		feedForwardMotifs.push_back(newMotif);
	}
	
	// motifs for network analysis: two edges
	std::list< std::list< ConnectionType > > recurrentIncompleteConnections;
	std::list< std::list< ConnectionType > > feedForwardIncompleteConnections;
	std::list< std::list< ConnectionType > > recurrentDivergentConnections;
	std::list< std::list< ConnectionType > > recurrentConvergentConnections;
	std::list< std::list< ConnectionType > > feedForwardConvergentConnections;
	std::list< std::list< ConnectionType > > feedForwardDivergentConnections;
	std::list< TripletMotif * > recurrentIncompleteMotifs;
	std::list< TripletMotif * > feedForwardIncompleteMotifs;
	std::list< TripletMotif * > recurrentDivergentMotifs;
	std::list< TripletMotif * > recurrentConvergentMotifs;
	std::list< TripletMotif * > feedForwardConvergentMotifs;
	std::list< TripletMotif * > feedForwardDivergentMotifs;
	
	std::list< ConnectionType > recurrentIncompleteConnections1;
	recurrentIncompleteConnections1.push_back(ConnectionType(0,1));
	recurrentIncompleteConnections1.push_back(ConnectionType(1,0));
	recurrentIncompleteConnections1.push_back(ConnectionType(0,2));
	recurrentIncompleteConnections1.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentIncompleteConnections2;
	recurrentIncompleteConnections2.push_back(ConnectionType(1,2));
	recurrentIncompleteConnections2.push_back(ConnectionType(2,1));
	recurrentIncompleteConnections2.push_back(ConnectionType(0,2));
	recurrentIncompleteConnections2.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentIncompleteConnections3;
	recurrentIncompleteConnections3.push_back(ConnectionType(0,1));
	recurrentIncompleteConnections3.push_back(ConnectionType(1,0));
	recurrentIncompleteConnections3.push_back(ConnectionType(1,2));
	recurrentIncompleteConnections3.push_back(ConnectionType(2,1));
	recurrentIncompleteConnections.push_back(recurrentIncompleteConnections1);
	recurrentIncompleteConnections.push_back(recurrentIncompleteConnections2);
	recurrentIncompleteConnections.push_back(recurrentIncompleteConnections3);
	std::list< ConnectionType > feedForwardIncompleteConnections1;
	feedForwardIncompleteConnections1.push_back(ConnectionType(0,1));
	feedForwardIncompleteConnections1.push_back(ConnectionType(1,2));
	std::list< ConnectionType > feedForwardIncompleteConnections2;
	feedForwardIncompleteConnections2.push_back(ConnectionType(2,1));
	feedForwardIncompleteConnections2.push_back(ConnectionType(1,0));
	std::list< ConnectionType > feedForwardIncompleteConnections3;
	feedForwardIncompleteConnections3.push_back(ConnectionType(1,2));
	feedForwardIncompleteConnections3.push_back(ConnectionType(2,0));
	std::list< ConnectionType > feedForwardIncompleteConnections4;
	feedForwardIncompleteConnections4.push_back(ConnectionType(0,2));
	feedForwardIncompleteConnections4.push_back(ConnectionType(2,1));
	std::list< ConnectionType > feedForwardIncompleteConnections5;
	feedForwardIncompleteConnections5.push_back(ConnectionType(2,0));
	feedForwardIncompleteConnections5.push_back(ConnectionType(0,1));
	std::list< ConnectionType > feedForwardIncompleteConnections6;
	feedForwardIncompleteConnections6.push_back(ConnectionType(1,0));
	feedForwardIncompleteConnections6.push_back(ConnectionType(0,2));
	feedForwardIncompleteConnections.push_back(feedForwardIncompleteConnections1);
	feedForwardIncompleteConnections.push_back(feedForwardIncompleteConnections2);
	feedForwardIncompleteConnections.push_back(feedForwardIncompleteConnections3);
	feedForwardIncompleteConnections.push_back(feedForwardIncompleteConnections4);
	feedForwardIncompleteConnections.push_back(feedForwardIncompleteConnections5);
	feedForwardIncompleteConnections.push_back(feedForwardIncompleteConnections6);
	std::list< ConnectionType > recurrentDivergentConnections1;
	recurrentDivergentConnections1.push_back(ConnectionType(0,1));
	recurrentDivergentConnections1.push_back(ConnectionType(1,0));
	recurrentDivergentConnections1.push_back(ConnectionType(0,2));
	std::list< ConnectionType > recurrentDivergentConnections2;
	recurrentDivergentConnections2.push_back(ConnectionType(0,1));
	recurrentDivergentConnections2.push_back(ConnectionType(1,0));
	recurrentDivergentConnections2.push_back(ConnectionType(1,2));
	std::list< ConnectionType > recurrentDivergentConnections3;
	recurrentDivergentConnections3.push_back(ConnectionType(1,2));
	recurrentDivergentConnections3.push_back(ConnectionType(2,1));
	recurrentDivergentConnections3.push_back(ConnectionType(1,0));
	std::list< ConnectionType > recurrentDivergentConnections4;
	recurrentDivergentConnections4.push_back(ConnectionType(1,2));
	recurrentDivergentConnections4.push_back(ConnectionType(2,1));
	recurrentDivergentConnections4.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentDivergentConnections5;
	recurrentDivergentConnections5.push_back(ConnectionType(0,2));
	recurrentDivergentConnections5.push_back(ConnectionType(2,0));
	recurrentDivergentConnections5.push_back(ConnectionType(0,1));
	std::list< ConnectionType > recurrentDivergentConnections6;
	recurrentDivergentConnections6.push_back(ConnectionType(0,2));
	recurrentDivergentConnections6.push_back(ConnectionType(2,0));
	recurrentDivergentConnections6.push_back(ConnectionType(2,1));
	recurrentDivergentConnections.push_back(recurrentDivergentConnections1);
	recurrentDivergentConnections.push_back(recurrentDivergentConnections2);
	recurrentDivergentConnections.push_back(recurrentDivergentConnections3);
	recurrentDivergentConnections.push_back(recurrentDivergentConnections4);
	recurrentDivergentConnections.push_back(recurrentDivergentConnections5);
	recurrentDivergentConnections.push_back(recurrentDivergentConnections6);
	std::list< ConnectionType > recurrentConvergentConnections1;
	recurrentConvergentConnections1.push_back(ConnectionType(0,1));
	recurrentConvergentConnections1.push_back(ConnectionType(1,0));
	recurrentConvergentConnections1.push_back(ConnectionType(2,0));
	std::list< ConnectionType > recurrentConvergentConnections2;
	recurrentConvergentConnections2.push_back(ConnectionType(0,1));
	recurrentConvergentConnections2.push_back(ConnectionType(1,0));
	recurrentConvergentConnections2.push_back(ConnectionType(2,1));
	std::list< ConnectionType > recurrentConvergentConnections3;
	recurrentConvergentConnections3.push_back(ConnectionType(1,2));
	recurrentConvergentConnections3.push_back(ConnectionType(2,1));
	recurrentConvergentConnections3.push_back(ConnectionType(0,1));
	std::list< ConnectionType > recurrentConvergentConnections4;
	recurrentConvergentConnections4.push_back(ConnectionType(1,2));
	recurrentConvergentConnections4.push_back(ConnectionType(2,1));
	recurrentConvergentConnections4.push_back(ConnectionType(0,2));
	std::list< ConnectionType > recurrentConvergentConnections5;
	recurrentConvergentConnections5.push_back(ConnectionType(0,2));
	recurrentConvergentConnections5.push_back(ConnectionType(2,0));
	recurrentConvergentConnections5.push_back(ConnectionType(1,0));
	std::list< ConnectionType > recurrentConvergentConnections6;
	recurrentConvergentConnections6.push_back(ConnectionType(0,2));
	recurrentConvergentConnections6.push_back(ConnectionType(2,0));
	recurrentConvergentConnections6.push_back(ConnectionType(1,2));
	recurrentConvergentConnections.push_back(recurrentConvergentConnections1);
	recurrentConvergentConnections.push_back(recurrentConvergentConnections2);
	recurrentConvergentConnections.push_back(recurrentConvergentConnections3);
	recurrentConvergentConnections.push_back(recurrentConvergentConnections4);
	recurrentConvergentConnections.push_back(recurrentConvergentConnections5);
	recurrentConvergentConnections.push_back(recurrentConvergentConnections6);
	std::list< ConnectionType > feedForwardConvergentConnections1;
	feedForwardConvergentConnections1.push_back(ConnectionType(1,0));
	feedForwardConvergentConnections1.push_back(ConnectionType(2,0));
	std::list< ConnectionType > feedForwardConvergentConnections2;
	feedForwardConvergentConnections2.push_back(ConnectionType(0,1));
	feedForwardConvergentConnections2.push_back(ConnectionType(2,1));
	std::list< ConnectionType > feedForwardConvergentConnections3;
	feedForwardConvergentConnections3.push_back(ConnectionType(0,2));
	feedForwardConvergentConnections3.push_back(ConnectionType(1,2));
	feedForwardConvergentConnections.push_back(feedForwardConvergentConnections1);
	feedForwardConvergentConnections.push_back(feedForwardConvergentConnections2);
	feedForwardConvergentConnections.push_back(feedForwardConvergentConnections3);
	std::list< ConnectionType > feedForwardDivergentConnections1;
	feedForwardDivergentConnections1.push_back(ConnectionType(0,1));
	feedForwardDivergentConnections1.push_back(ConnectionType(0,2));
	std::list< ConnectionType > feedForwardDivergentConnections2;
	feedForwardDivergentConnections2.push_back(ConnectionType(1,0));
	feedForwardDivergentConnections2.push_back(ConnectionType(1,2));
	std::list< ConnectionType > feedForwardDivergentConnections3;
	feedForwardDivergentConnections3.push_back(ConnectionType(2,0));
	feedForwardDivergentConnections3.push_back(ConnectionType(2,1));
	feedForwardDivergentConnections.push_back(feedForwardDivergentConnections1);
	feedForwardDivergentConnections.push_back(feedForwardDivergentConnections2);
	feedForwardDivergentConnections.push_back(feedForwardDivergentConnections3);
	
	for(connectionIt = recurrentIncompleteConnections.begin(); connectionIt != recurrentIncompleteConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		recurrentIncompleteMotifs.push_back(newMotif);
	}
	for(connectionIt = feedForwardIncompleteConnections.begin(); connectionIt != feedForwardIncompleteConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		feedForwardIncompleteMotifs.push_back(newMotif);
	}
	for(connectionIt = recurrentDivergentConnections.begin(); connectionIt != recurrentDivergentConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		recurrentDivergentMotifs.push_back(newMotif);
	}
	for(connectionIt = recurrentConvergentConnections.begin(); connectionIt != recurrentConvergentConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		recurrentConvergentMotifs.push_back(newMotif);
	}
	for(connectionIt = feedForwardConvergentConnections.begin(); connectionIt != feedForwardConvergentConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		feedForwardConvergentMotifs.push_back(newMotif);
	}
	for(connectionIt = feedForwardDivergentConnections.begin(); connectionIt != feedForwardDivergentConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		feedForwardDivergentMotifs.push_back(newMotif);
	}
	
	// motifs for network analysis: one edge
	std::list< std::list< ConnectionType > > recurrentSparseConnections;
	std::list< std::list< ConnectionType > > feedForwardSparseConnections;
	std::list< TripletMotif * > recurrentSparseMotifs;
	std::list< TripletMotif * > feedForwardSparseMotifs;
	
	std::list< ConnectionType > recurrentSparseConnections1;
	recurrentSparseConnections1.push_back(ConnectionType(0,1));
	recurrentSparseConnections1.push_back(ConnectionType(1,0));
	std::list< ConnectionType > recurrentSparseConnections2;
	recurrentSparseConnections2.push_back(ConnectionType(1,2));
	recurrentSparseConnections2.push_back(ConnectionType(2,1));
	std::list< ConnectionType > recurrentSparseConnections3;
	recurrentSparseConnections3.push_back(ConnectionType(0,2));
	recurrentSparseConnections3.push_back(ConnectionType(2,0));
	recurrentSparseConnections.push_back(recurrentSparseConnections1);
	recurrentSparseConnections.push_back(recurrentSparseConnections2);
	recurrentSparseConnections.push_back(recurrentSparseConnections3);
	std::list< ConnectionType > feedForwardSparseConnections1;
	feedForwardSparseConnections1.push_back(ConnectionType(0,1));
	std::list< ConnectionType > feedForwardSparseConnections2;
	feedForwardSparseConnections2.push_back(ConnectionType(1,0));
	std::list< ConnectionType > feedForwardSparseConnections3;
	feedForwardSparseConnections3.push_back(ConnectionType(1,2));
	std::list< ConnectionType > feedForwardSparseConnections4;
	feedForwardSparseConnections4.push_back(ConnectionType(2,1));
	std::list< ConnectionType > feedForwardSparseConnections5;
	feedForwardSparseConnections5.push_back(ConnectionType(0,2));
	std::list< ConnectionType > feedForwardSparseConnections6;
	feedForwardSparseConnections6.push_back(ConnectionType(2,0));
	feedForwardSparseConnections.push_back(feedForwardSparseConnections1);
	feedForwardSparseConnections.push_back(feedForwardSparseConnections2);
	feedForwardSparseConnections.push_back(feedForwardSparseConnections3);
	feedForwardSparseConnections.push_back(feedForwardSparseConnections4);
	feedForwardSparseConnections.push_back(feedForwardSparseConnections5);
	feedForwardSparseConnections.push_back(feedForwardSparseConnections6);
	
	for(connectionIt = recurrentSparseConnections.begin(); connectionIt != recurrentSparseConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		recurrentSparseMotifs.push_back(newMotif);
	}
	for(connectionIt = feedForwardSparseConnections.begin(); connectionIt != feedForwardSparseConnections.end(); ++connectionIt)
	{
		TripletMotif * newMotif = new TripletMotif(*connectionIt);
		feedForwardSparseMotifs.push_back(newMotif);
	}
	
	// motifs for network analysis: empty motif
	std::list< ConnectionType > emtpyConnections;
	std::list< TripletMotif * > emptyMotif;
	TripletMotif * emtpyMotifPtr = new TripletMotif(emtpyConnections);
	emptyMotif.push_back(emtpyMotifPtr);
	
	std::map< unsigned int, std::list< TripletMotif* > > nonRedundantTriplets;
	nonRedundantTriplets[0] = recurrentLoopMotifs;
	nonRedundantTriplets[1] = directedLoopMotifs;
	nonRedundantTriplets[2] = recurrentIncompleteLoopMotifs;
	nonRedundantTriplets[3] = directedRecurrentLoopMotifs;
	nonRedundantTriplets[4] = recurrentFeedForwardConvergentMotifs;
	nonRedundantTriplets[5] = recurrentFeedForwardDivergentMotifs;
	nonRedundantTriplets[6] = feedForwardMotifs;
	nonRedundantTriplets[7] = recurrentIncompleteMotifs;
	nonRedundantTriplets[8] = feedForwardIncompleteMotifs;
	nonRedundantTriplets[9] = recurrentDivergentMotifs;
	nonRedundantTriplets[10] = recurrentConvergentMotifs;
	nonRedundantTriplets[11] = feedForwardConvergentMotifs;
	nonRedundantTriplets[12] = feedForwardDivergentMotifs;
	nonRedundantTriplets[13] = recurrentSparseMotifs;
	nonRedundantTriplets[14] = feedForwardSparseMotifs;
	nonRedundantTriplets[15] = emptyMotif;
	
	return nonRedundantTriplets;
}

std::list< CellTriplet* > MatrixAnalyzer::initializeNonRedundantCellTriplets(SelectionType preSelection, SelectionType postSelection, unsigned int nrOfTriplets)
{
	std::cout << "Selecting " << nrOfTriplets << " cell triplets from connection matrix..." << std::endl;
// 	std::srand(1234567);
	std::srand(std::time(NULL));
	const unsigned int NMAX = std::min(preSelection.size(), postSelection.size());
	std::list< CellTriplet * > triplets;
	std::list< std::vector< unsigned int > > usedTriplets;
	while(triplets.size() < nrOfTriplets)
	{
		bool nonRedundantTriplet = true;
		unsigned int index1, index2, index3;
		unsigned int preCell1, preCell2, preCell3, postCell1, postCell2, postCell3;
		
		index1 = std::rand()%NMAX;
		index2 = std::rand()%NMAX;
		index3 = std::rand()%NMAX;
		
		if(index1 == index2 || index1 == index3 || index2 == index3)
		{
// 			std::cout << std::endl;
// 			std::cout << "index1 = " << index1 << std::endl;
// 			std::cout << "index2 = " << index2 << std::endl;
// 			std::cout << "index3 = " << index3 << std::endl;
// 			std::cout << "Duplicate index, repeating draw..." << std::endl;
			continue;
		}
		
		std::list< std::vector< unsigned int > >::const_iterator usedTripletsIt;
		for(usedTripletsIt = usedTriplets.begin(); usedTripletsIt != usedTriplets.end(); ++usedTripletsIt)
		{
			unsigned int matches = 0;
			if(index1 == usedTripletsIt->at(0) || index1 == usedTripletsIt->at(1) || index1 == usedTripletsIt->at(2))
			{
				++matches;
			}
			if(index2 == usedTripletsIt->at(0) || index2 == usedTripletsIt->at(1) || index2 == usedTripletsIt->at(2))
			{
				++matches;
			}
			if(index3 == usedTripletsIt->at(0) || index3 == usedTripletsIt->at(1) || index3 == usedTripletsIt->at(2))
			{
				++matches;
			}
			
			if(matches > 1)
			{
				nonRedundantTriplet = false;
				break;
			}
		}
		
		if(nonRedundantTriplet)
		{
			std::vector< unsigned int > usedTriplet;
			usedTriplet.push_back(index1);
			usedTriplet.push_back(index2);
			usedTriplet.push_back(index3);
			usedTriplets.push_back(usedTriplet);
			preCell1 = preSelection[index1];
			preCell2 = preSelection[index2];
			preCell3 = preSelection[index3];
			postCell1 = postSelection[index1];
			postCell2 = postSelection[index2];
			postCell3 = postSelection[index3];
			CellTriplet * newTriplet = new CellTriplet(preCell1, preCell2, preCell3, postCell1, postCell2, postCell3);
			newTriplet->setInnervationMatrix(connectome);
			triplets.push_back(newTriplet);
		}
		if(!(triplets.size()%((unsigned int)(0.1*nrOfTriplets))))
		{
			double percentage = (double)(triplets.size())/nrOfTriplets*100;
			std::flush(std::cout << percentage << "% done\r");
		}
	}
	std::cout << std::endl;
	return triplets;
}

double MatrixAnalyzer::computeHistogramRMSE(Profile* profile1, Profile* profile2)
{
	if(profile1->getProfile()->size() != profile2->getProfile()->size())
	{
		throw std::out_of_range("MatrixAnalyzer::computeHistogramRMSE: range of histograms profile1 and profile2 not matching");
	}
	double rmse = 0;
	for(int ii = 0; ii < profile1->getProfile()->size(); ++ii)
	{
		rmse += (profile1->getProfile()->at(ii) - profile2->getProfile()->at(ii))*(profile1->getProfile()->at(ii) - profile2->getProfile()->at(ii));
	}
	
	return sqrt(rmse/profile1->getProfile()->size());
}

void MatrixAnalyzer::initializeConstants()
{
	if(this->int2Labels.size())
		this->int2Labels.clear();
	this->int2Labels.insert(std::pair< unsigned int, const char * >(Alpha, "Alpha"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(A1, "A1"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(A2, "A2"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(A3, "A3"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(A4, "A4"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(Beta, "Beta"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(B1, "B1"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(B2, "B2"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(B3, "B3"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(B4, "B4"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(Gamma, "Gamma"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(C1, "C1"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(C2, "C2"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(C3, "C3"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(C4, "C4"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(C5, "C5"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(C6, "C6"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(Delta, "Delta"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(D1, "D1"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(D2, "D2"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(D3, "D3"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(D4, "D4"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(D5, "D5"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(D6, "D6"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(E1, "E1"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(E2, "E2"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(E3, "E3"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(E4, "E4"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(E5, "E5"));
	this->int2Labels.insert(std::pair< unsigned int, const char * >(E6, "E6"));
	if(this->labels2Int.size())
		this->labels2Int.clear();
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("Alpha"), Alpha));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("A1"), A1));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("A2"), A2));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("A3"), A3));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("A4"), A4));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("Beta"), Beta));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("B1"), B1));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("B2"), B2));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("B3"), B3));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("B4"), B4));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("Gamma"), Gamma));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("C1"), C1));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("C2"), C2));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("C3"), C3));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("C4"), C4));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("C5"), C5));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("C6"), C6));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("Delta"), Delta));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("D1"), D1));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("D2"), D2));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("D3"), D3));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("D4"), D4));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("D5"), D5));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("D6"), D6));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("E1"), E1));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("E2"), E2));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("E3"), E3));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("E4"), E4));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("E5"), E5));
	this->labels2Int.insert(std::pair< std::string, unsigned int >(std::string("E6"), E6));
	if(this->celltypeLabels2Int.size())
		this->celltypeLabels2Int.clear();
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L2"),L2));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L34"),L34));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4py"),L4py));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4sp"),L4sp));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4ss"),L4ss));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5st"),L5st));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5tt"),L5tt));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6cc"),L6cc));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccinv"),L6ccinv));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ct"),L6ct));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal"),SymLocal));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal1"),SymLocal1));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal2"),SymLocal2));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal3"),SymLocal3));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal4"),SymLocal4));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal5"),SymLocal5));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal6"),SymLocal6));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L1"),L1));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L23Trans"),L23Trans));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Sym"),L45Sym));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Peak"),L45Peak));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L56Trans"),L56Trans));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L2axon"),L2axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L34axon"),L34axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4pyaxon"),L4pyaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4spaxon"),L4spaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4ssaxon"),L4ssaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5staxon"),L5staxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5ttaxon"),L5ttaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccaxon"),L6ccaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccinvaxon"),L6ccinvaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ctaxon"),L6ctaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocalaxon"),SymLocalaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal1axon"),SymLocal1axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal2axon"),SymLocal2axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal3axon"),SymLocal3axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal4axon"),SymLocal4axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal5axon"),SymLocal5axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal6axon"),SymLocal6axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L1axon"),L1axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L23Transaxon"),L23Transaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Symaxon"),L45Symaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Peakaxon"),L45Peakaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L56Transaxon"),L56Transaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("VPM"),VPM));
	if(this->int2CelltypeLabels.size())
		this->int2CelltypeLabels.clear();
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L2,"L2"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L34,"L34"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4py,"L4py"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4sp,"L4sp"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4ss,"L4ss"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5st,"L5st"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5tt,"L5tt"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6cc,"L6cc"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccinv,"L6ccinv"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ct,"L6ct"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal,"SymLocal"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal1,"SymLocal1"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal2,"SymLocal2"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal3,"SymLocal3"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal4,"SymLocal4"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal5,"SymLocal5"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal6,"SymLocal6"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L1,"L1"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L23Trans,"L23Trans"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Sym,"L45Sym"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Peak,"L45Peak"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L56Trans,"L56Trans"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L2axon,"L2axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L34axon,"L34axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4pyaxon,"L4pyaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4spaxon,"L4spaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4ssaxon,"L4ssaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5staxon,"L5staxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5ttaxon,"L5ttaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccaxon,"L6ccaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccinvaxon,"L6ccinvaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ctaxon,"L6ctaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocalaxon,"SymLocalaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal1axon,"SymLocal1axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal2axon,"SymLocal2axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal3axon,"SymLocal3axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal4axon,"SymLocal4axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal5axon,"SymLocal5axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal6axon,"SymLocal6axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L1axon,"L1axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L23Transaxon,"L23Transaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Symaxon,"L45Symaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Peakaxon,"L45Peakaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L56Transaxon,"L56Transaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(VPM,"VPM"));
}

