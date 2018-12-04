#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "../../common/typedefs.h"
#include "../../common/basics.h"
#include "../../common/amiraReader.h"
#include "../../Interneuron/src/helper.h"
#include <set>
#include <utility>

int main(int argc, const char** argv)
{
	/* Truncate in vivo morphology
	 * Outputs truncated morphology
	 * Outputs
	 * filename,basallength_original,apicallength_original,axonlength_original,basallength_cut,apicallength_cut,axonlength_cut,SliceDistToSoma,SliceWidth,SlicingDimension
	 */
	if(argc == 6)
	{
		// arg[1] = Path/to/AmiraSpatialGraph
		// arg[2] = Distance of Soma to Slicing Edge (if negative, use maxSoma, otherwise use minSoma; usually between 10 and 50)
		// arg[3] = Slicing Width (usually 300 or 350)
		// arg[4] = Orientation (0: XCOORD; 1: YCOORD; 2: ZCOORD)
		// arg[5] = Outputpath/to/AmiraSpatialGraph
		// e.g. ./GenerateTruncatedCell "/nas1/Data_daniel/Network/L5/L5_InVitro/L5tt_ExampleCell.am" 10 300 0 "/nas1/Data_daniel/Network/L5/L5_InVitro/L5tt_ExampleCell_Cut.am"

		const char * Inputfilename = argv[1]; //"/nas1/Data_daniel/Network/L5/L5_InVitro/L5tt_ExampleCell.am";
		double distSomaSlice = atof(argv[2]);
		double sliceWdth = atof(argv[3]);
		int cut_COORD = atoi(argv[4]);
		std::string Outputfilename = argv[5];

		// Read in Spatial Graph
		AmiraSpatialGraph * sg = helper::getSpatialGraph(Inputfilename);

		double lenBasalDendrite = helper::lengthCell(sg, BasalDendrite);
		double lenApicalDendrite = helper::lengthCell(sg, ApicalDendrite);
		double lenAxon = helper::lengthCell(sg, Axon);

		// Get bounding box of soma
		double boundSoma[6];
		sg->getBoundingBox(Soma, boundSoma);

		// Compute slicing box given by distance to soma and slicing width
		// if distance to soma is negative, take maximum soma position [......maxSoma..]
		// if distance to soma is positive, take minimum soma position [..minSoma......]
		// enforces that the soma is always within the slicing box
		double s1;
		double s2;
		if (distSomaSlice<0)
		{	// [maxSoma-dist-sliceWdth maxSoma-dist]
			double maxSoma = boundSoma[cut_COORD*2+1];
			s1 = maxSoma-distSomaSlice-sliceWdth;
			s2 = maxSoma-distSomaSlice;
		}
		else
		{	// [minSoma-dist minSoma-dist+sliceWdth]
			double minSoma = boundSoma[cut_COORD*2];
			s1 = minSoma-distSomaSlice;
			s2 = minSoma-distSomaSlice+sliceWdth;
		}

		if ((boundSoma[cut_COORD*2]<s1) || (boundSoma[cut_COORD*2+1]>s2))
		{
			std::cout << "WARNING! Soma is truncated!" << std::endl;
			std::cout << "	Soma [min max]: [" << boundSoma[cut_COORD*2] << " " << boundSoma[cut_COORD*2+1] << "]" << std::endl;
			std::cout << "	Slice [min max]: [" << s1 << " " << s2 << "]" << std::endl;
		}

		// Copy Spatial Graph
		// sgOriginal 	: original in vivo morphology
		// sg 			: truncated morphology
		AmiraSpatialGraph * sgOriginal = new AmiraSpatialGraph;
		sgOriginal->mergeSpatialGraph(sg);

		// Truncate Spatial Graph
		double bounds[2] = {s1, s2};
		sg->truncateSpatialGraph(bounds,cut_COORD);

		// Shift Spatial Graphs in Slice Coordinate System. SliceDimension [0 300] LateralDimension [Soma at 0] Depth [Keep]
		double centerPt[3];
		helper::centerOfSpatialGraph(Soma, centerPt, sg);
		centerPt[cut_COORD] = -s1;
		if (cut_COORD==X_COORD)
		{
			centerPt[Y_COORD] = -centerPt[Y_COORD];
		}
		else
		{
			centerPt[X_COORD] = -centerPt[X_COORD];
		}
		centerPt[Z_COORD] = 0;
		TransformPointerType sub = TransformPointerType::New();
		sub->Translate(centerPt);

		sg->setTransformation(sub);
		sg->applyTransformation();

		sgOriginal->setTransformation(sub);
		sgOriginal->applyTransformation();

		// Save truncated Spatial Graph
		Reader * amWriterCut = new Reader(Outputfilename.c_str(), Outputfilename.c_str());
		amWriterCut->setSpatialGraph(sg);
		amWriterCut->writeSpatialGraphFile();

		// Save truncated Spatial Graph
		std::string OutputfilenameOriginal = Outputfilename + "_invivo";
		Reader * amWriter = new Reader(OutputfilenameOriginal.c_str(), OutputfilenameOriginal.c_str());
		amWriter->setSpatialGraph(sgOriginal);
		amWriter->writeSpatialGraphFile();

		// Compute Length of Truncated SpatialGraph and Display invivo and truncated length values, and parameters
		double lenBasalDendriteCut = helper::lengthCell(sg, BasalDendrite);
		double lenApicalDendriteCut = helper::lengthCell(sg, ApicalDendrite);
		double lenAxonCut = helper::lengthCell(sg, Axon);

		std::cout << Outputfilename << "," << lenBasalDendrite << "," << lenApicalDendrite << "," << lenAxon << ",";
		std::cout << lenBasalDendriteCut << "," << lenApicalDendriteCut << "," << lenAxonCut << ",";
		std::cout << distSomaSlice << "," << sliceWdth << "," << cut_COORD << std::endl;
	}
	else if (argc == 5)
	{
		/* Shifts L5tt Dendrites laterally
		 * Outputs shifted L5tt File */

		// arg[1] = Path/to/cell.txt
		// arg[2] = min lateral
		// arg[3] = max lateral
		// arg[4] = Outputpath
		// e.g. ./GenerateTruncatedCell "/nas1/Data_daniel/Network/L5/L5_InVitro/Cut_300/cell.txt" -175 175 "/nas1/Data_daniel/Network/L5/L5_InVitro/L5ttDendrites_am/"

		std::string pathToCellList = argv[1]; // "/nas1/Data_daniel/Network/L5/L5_InVitro/Cut_300/cell.txt";
		double minLateral = atof(argv[2]); // -175
		double maxLateral = atof(argv[3]); // 175
		std::string outputPath = argv[4]; // "/nas1/Data_daniel/Network/L5/L5_InVitro/L5ttDendrites/"
		double stepSize = 50;

		std::string pathToCell = helper::getRootFromPath(pathToCellList.c_str());
		std::ifstream inputStream(pathToCellList.c_str());
		if(!inputStream.fail())
		{
			std::string currentLine;

			while(!inputStream.eof())
			{
				getline(inputStream,currentLine);
				if (currentLine.empty())
					break;

				// check if filename was cut in x (_0.am) or y (_0_*.am)
				bool invivo = false;

				std::size_t found = currentLine.find(".am");
				std::string filename = currentLine.substr(0,found);

				std::size_t found2 = filename.find("_invivo");
				if (found2!=std::string::npos)
				{
					filename = currentLine.substr(0,found2);
					invivo = true;
				}

				unsigned foundY = filename.find_last_of("_");
				std::string strY = filename.substr(foundY+1);
				unsigned foundX = filename.find_last_of("_",foundY-1);
				std::string strX = filename.substr(foundX+1,foundY-foundX-1);

				//std::cout << filename << " " << atoi(strX.c_str()) << " " << atoi(strY.c_str()) << std::endl;

				int posX = atoi(strX.c_str());
				int posY = atoi(strY.c_str());

				int cut_COORD;
				if ((posY == 0) && (posX == 0))
				{
					std::cout << "Error! Two slicing dimensions (0) were found in filename!" << std::endl;
					std::cout << filename << " " << strX << " " << strY << std::endl;
					return 0;
				}
				if ((posY != 0) && (posX != 0))
				{
					std::cout << "Error! Slicing Dimension (0) was not found in filename!" << std::endl;
					std::cout << filename << " " << strX << " " << strY << std::endl;
					return 0;
				}

				std::string currentFilename = pathToCell + currentLine;

				// Read in Spatial Graph
				AmiraSpatialGraph * sg = helper::getSpatialGraph(currentFilename.c_str());

				// Only keep Basal and ApicalDendrites
				sg->removeLabel(Axon);

				for (double m = minLateral; m <= maxLateral; m+=stepSize)
				{
					AmiraSpatialGraph * sgShifted = new AmiraSpatialGraph;
					sgShifted->mergeSpatialGraph(sg);

					// shift Morphology in Slice Coordinate System
					double centerPt[3];
					helper::centerOfSpatialGraph(Soma, centerPt, sgShifted);
					std::string shift;

					if (posY == 0)
					{
						centerPt[Y_COORD] = m;
						centerPt[X_COORD] = 0;
						shift = "Y";
					}
					else
					{
						centerPt[X_COORD] = m;
						centerPt[Y_COORD] = 0;
						shift = "X";
					}
					centerPt[Z_COORD] = 0;
					TransformPointerType sub = TransformPointerType::New();
					sub->Translate(centerPt);
					sgShifted->setTransformation(sub);
					sgShifted->applyTransformation();

					// Store shifted Dendritic Morphology
					// path/to/currentName_Dist_Xpos
					std::ostringstream strs;
					strs << m;
					std::string Outputfilename = outputPath + filename + "_shift" + shift + "_" + strs.str();
					if (invivo)
					{
						Outputfilename += "_invivo";
					}

					Reader * amWriter = new Reader(Outputfilename.c_str(), Outputfilename.c_str());
					amWriter->setSpatialGraph(sgShifted);
					amWriter->writeSpatialGraphFile();

					delete amWriter;
					delete sgShifted;
				}
			}
		}
	}
	/* New Version for L4ss
	 * 1) Get Axon Morphology, construct slice around it, truncate Axon and save spatigal graph
	 * 2) Find all Dendrites that are within slice, truncate respective Dendrites and save spatial graphs
	 * 3) Output .csv file containing all parameters
	 * Input:
	 * 		-path/to/Axons -> cells.txt
	 * 		-path/to/Dendrites -> cells.txt
	 * 		-min Slicing box
	 * 		-slicing dimension
	 * 		-slicing width
	 * 		-max lateral
	 * 		-outputpath */
	else if (argc == 8)
	{
		// Assign Parameters
		double boundsSlice[2];
		std::string pathToAxonMorphologies = argv[1]; // "/nas1/Data_daniel/Network/L4ss/L4ss_InVitro/L4ssAxons/";
		std::string pathToDendriteMorphologies = argv[2]; // "/nas1/Data_daniel/Network/L4ss/L4ss_InVitro/L4ssDendrites/";
		boundsSlice[0] = atof(argv[3]); // min Slice
		int cut_COORD = atoi(argv[4]); // (0: XCOORD; 1: YCOORD; 2: ZCOORD)
		double sliceWdth = atof(argv[5]); // 300
		double maxLateral = atof(argv[6]); // 400
		std::string outputPath = argv[7]; // "/nas1/Data_daniel/Network/L4ss/L4ss_InVitro/SliceAxon1_10/"
		boundsSlice[1] = boundsSlice[0]+sliceWdth;

		// Get pathfilenames of Dendrites and Axons
		std::string pathtmp = pathToDendriteMorphologies + "cells.txt";
		std::vector<std::string> filenamesDendrite = helper::returnNames(pathtmp.c_str(), pathToDendriteMorphologies.c_str());
		if (filenamesDendrite.size()==0)
		{
			std::cout << "ERROR! No Dendrite Cell List found in " << pathtmp << "!" << std::endl;
			return 0;
		}

		pathtmp = pathToAxonMorphologies + "cells.txt";
		std::vector<std::string> filenamesAxon = helper::returnNames(pathtmp.c_str(), pathToAxonMorphologies.c_str());
		if (filenamesAxon.size()==0)
		{
			std::cout << "ERROR! No Axon Cell List found in " << pathtmp << "!" << std::endl;
			return 0;
		}

		if (cut_COORD==Z_COORD)
		{
			std::cout << "ERROR! Cutting Dimension along cortical axis is not defined yet!" << std::endl;
			return 0;
		}

		helper::mkDir(outputPath.c_str());
		std::string outputCSV = outputPath + "parameters.csv";
		std::ofstream outStream(outputCSV.c_str());
		outStream << "Name,BasalLength,ApicalLength,AxonLength,BasalLengthSlice,ApicalLengthSlice,AxonLengthSlice,TissueDepth,SliceWidth,CutDimension,SomaX,SomaY,SomaZ,SliceBound1,SliceBound2" << std::endl;

		if(!outStream.fail())
		{
			/*****************
			 * Truncate Axon *
			 *****************/
			double somaExtent[2] = {1E9, -1E9};
			int countAxons = 0;

			for (std::vector<std::string>::iterator nametmp = filenamesAxon.begin(); nametmp != filenamesAxon.end(); ++nametmp)
			{
				AmiraSpatialGraph * sg = helper::getSpatialGraph((*nametmp).c_str());
				// Soma Position and BB
				double posSoma[3];
				helper::centerOfSpatialGraph(Soma, posSoma, sg);
				double boundSoma[6];
				sg->getBoundingBox(Soma, boundSoma);

				// Check if Somata of Axon Morphology is within Slice
				if ((boundSoma[cut_COORD*2]<boundsSlice[0]) || (boundSoma[cut_COORD*2+1]>boundsSlice[1]))
				{
					delete sg;
					continue;
				}

				double lenAxon = helper::lengthCell(sg, Axon);

				if (lenAxon==0)
				{
					std::cout << "WARNING! SpatialGraph " << (*nametmp) << " does not have any axon labels!" << std::endl;
					delete sg;
					continue;
				}

				// Truncate Axon
				sg->truncateSpatialGraph(boundsSlice,cut_COORD);
				// Compute Length of Truncated Axon SpatialGraph
				double lenAxonCut = helper::lengthCell(sg, Axon);
				double distSomaSlice = posSoma[cut_COORD]-boundsSlice[0];

				std::string Outputfilename = outputPath + helper::getFilenameFromPath((*nametmp).c_str());
				std::size_t found = Outputfilename.find(".am");
				Outputfilename = Outputfilename.substr(0,found);

				Reader * amWriter = new Reader(Outputfilename.c_str(), Outputfilename.c_str());
				amWriter->setSpatialGraph(sg);
				amWriter->writeSpatialGraphFile();

				outStream << Outputfilename << ",";
				outStream << "0,0," << lenAxon << ",";
				outStream << "0,0," << lenAxonCut << ",";
				outStream << distSomaSlice << "," << sliceWdth << "," << cut_COORD << ",";
				outStream << posSoma[0] << "," << posSoma[1] << "," << posSoma[2] << ",";
				outStream << boundsSlice[0] << "," << boundsSlice[1] << std::endl;

				// Update somaExtent for maximal InterSomaDistance
				if (cut_COORD==X_COORD)
				{
					somaExtent[0] = std::min(somaExtent[0],posSoma[Y_COORD]);
					somaExtent[1] = std::max(somaExtent[1],posSoma[Y_COORD]);
				}
				else
				{
					somaExtent[0] = std::min(somaExtent[0],posSoma[X_COORD]);
					somaExtent[1] = std::max(somaExtent[1],posSoma[X_COORD]);
				}

				delete sg;
				delete amWriter;
				countAxons++;
			}

			if (countAxons==0)
			{
				outStream.close();
				std::cout << "WARNING! No Axon found in slice [" << boundsSlice[0] << " " << boundsSlice[1] << "]" << std::endl;
				//return 0;
			}

			/**********************
			 * Truncate Dendrites *
			 **********************/
			int countDendrites = 0;
			for (std::vector<std::string>::iterator nametmp = filenamesDendrite.begin(); nametmp != filenamesDendrite.end(); ++nametmp)
			{
				AmiraSpatialGraph * sg = helper::getSpatialGraph((*nametmp).c_str());
				// Soma Position and BB
				double posSoma[3];
				helper::centerOfSpatialGraph(Soma, posSoma, sg);
				double boundSoma[6];
				sg->getBoundingBox(Soma, boundSoma);

				// Check if Somata of Dendrite Morphology is within Slice
				if ((boundSoma[cut_COORD*2]<boundsSlice[0]) || (boundSoma[cut_COORD*2+1]>boundsSlice[1]))
				{
					delete sg;
					continue;
				}

				// Check if Somata of Dendrite is within range of Somata of Axon
				// If Somata are more than maxLateral apart, skip Dendrite
				double disttmp[2];
				if (cut_COORD==X_COORD)
				{
					disttmp[0] = fabs(somaExtent[0]-posSoma[Y_COORD]);
					disttmp[1] = fabs(somaExtent[1]-posSoma[Y_COORD]);
				}
				else
				{
					disttmp[0] = fabs(somaExtent[0]-posSoma[X_COORD]);
					disttmp[1] = fabs(somaExtent[1]-posSoma[X_COORD]);
				}

				if (disttmp[0]>maxLateral && disttmp[1]>maxLateral)
				{
					delete sg;
					continue;
				}

				double lenApical = helper::lengthCell(sg, ApicalDendrite);
				double lenBasal = helper::lengthCell(sg, BasalDendrite);

				if (lenBasal==0 && lenApical==0)
				{
					std::cout << "WARNING! SpatialGraph " << (*nametmp) << " does not have any dendrite labels!" << std::endl;
					delete sg;
					continue;
				}

				// Truncate Dendrite
				sg->truncateSpatialGraph(boundsSlice,cut_COORD);
				// Compute Length of Truncated Axon SpatialGraph
				double lenApicalCut = helper::lengthCell(sg, ApicalDendrite);
				double lenBasalCut = helper::lengthCell(sg, BasalDendrite);
				double distSomaSlice = posSoma[cut_COORD]-boundsSlice[0];

				std::string Outputfilename = outputPath + helper::getFilenameFromPath((*nametmp).c_str());
				std::size_t found = Outputfilename.find(".am");
				Outputfilename = Outputfilename.substr(0,found);

				Reader * amWriter = new Reader(Outputfilename.c_str(), Outputfilename.c_str());
				amWriter->setSpatialGraph(sg);
				amWriter->writeSpatialGraphFile();

				outStream << Outputfilename << ",";
				outStream << lenBasal << "," << lenApical << ",0,";
				outStream << lenBasalCut << "," << lenApicalCut << ",0,";
				outStream << distSomaSlice << "," << sliceWdth << "," << cut_COORD << ",";
				outStream << posSoma[0] << "," << posSoma[1] << "," << posSoma[2] << ",";
				outStream << boundsSlice[0] << "," << boundsSlice[1] << std::endl;

				delete sg;
				delete amWriter;
				countDendrites++;
			}

			if (countDendrites==0)
			{
				outStream.close();
				std::cout << "WARNING! No Dendrites found in slice [" << boundsSlice[0] << " " << boundsSlice[1] << "]" << std::endl;
				return 0;
			}

			outStream.close();
		}
		else
		{
			std::cout << "ERROR! Writing CSV file failed! Path: " << outputCSV << std::endl;
		}
	}
	/* Extracting Morphologies from SpatialGraphSetFile.bb
	 * Only implemented for L4ss, L5tt, and L5st so far */
	else if (argc == 2)
	{
		std::string CellTypeStr = argv[1];
		if (CellTypeStr.compare("L4ss") == 0)
		{
			const char * spatialGraphSetFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/Realization20151007_ascii/Morphologies.am";
			int CellType = L4ss;
			const char * outputpath = "/nas1/Data_daniel/Network/L4ss/L4ss_InVitro/L4ssDendrites/";
			double bounds[6];
			bounds[0] = -550;
			bounds[1] = 410;
			bounds[2] = -140;
			bounds[3] = 840;
			bounds[4] = -(std::numeric_limits<double>::infinity());
			bounds[5] = std::numeric_limits<double>::infinity();

			helper::extractMorphologies(spatialGraphSetFilename,CellType,outputpath,bounds);
		}
		else if (CellTypeStr.compare("L5tt") == 0)
		{
			const char * spatialGraphSetFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/Realization20151007_ascii/Morphologies.am";
			int CellType = L5tt;
			const char * outputpath = "/nas1/Data_daniel/Network/L5/L5_InVitroNetwork/L5ttDendrites/";
			double bounds[6];
			bounds[0] = -490;
			bounds[1] = 340;
			bounds[2] = -20;
			bounds[3] = 820;
			bounds[4] = -(std::numeric_limits<double>::infinity());
			bounds[5] = std::numeric_limits<double>::infinity();

			helper::extractMorphologies(spatialGraphSetFilename,CellType,outputpath,bounds);
		}
		else if (CellTypeStr.compare("L5st") == 0)
		{
			const char * spatialGraphSetFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/Realization20151007_ascii/Morphologies.am";
			int CellType = L5st;
			const char * outputpath = "/nas1/Data_daniel/Network/L5/L5_InVitroNetwork/L5stDendrites/";
			double bounds[6];
			bounds[0] = -490;
			bounds[1] = 340;
			bounds[2] = -20;
			bounds[3] = 820;
			bounds[4] = -(std::numeric_limits<double>::infinity());
			bounds[5] = std::numeric_limits<double>::infinity();

			helper::extractMorphologies(spatialGraphSetFilename,CellType,outputpath,bounds);
		}
		else if (CellTypeStr.compare("L2") == 0)
		{	// Extract only cells close to C2 barrel
			const char * spatialGraphSetFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/Realization20151007_ascii/Morphologies.am";
			int CellType = L2;
			const char * outputpath = "/nas1/Data_daniel/Network/L4-L23/L2Dendrites/";
			double bounds[6];
			bounds[0] = -87-200; // xCoord: C2: -87 C2rad: 181
			bounds[1] = -87+200;
			bounds[2] = 410-200; // yCoord: C2: 410 C2rad: 181
			bounds[3] = 410+200;
			bounds[4] = -(std::numeric_limits<double>::infinity());
			bounds[5] = std::numeric_limits<double>::infinity();

			helper::extractMorphologies(spatialGraphSetFilename,CellType,outputpath,bounds);
		}
		else if (CellTypeStr.compare("L34") == 0)
		{ 	// Extract only cells close to C2 barrel
			const char * spatialGraphSetFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/Realization20151007_ascii/Morphologies.am";
			int CellType = L34;
			const char * outputpath = "/nas1/Data_daniel/Network/L4-L23/L34Dendrites/";
			double bounds[6];
			bounds[0] = -87-200; // xCoord: C2: -87 C2rad: 181
			bounds[1] = -87+200;
			bounds[2] = 410-200; // yCoord: C2: 410 C2rad: 181
			bounds[3] = 410+200;
			bounds[4] = -(std::numeric_limits<double>::infinity());
			bounds[5] = std::numeric_limits<double>::infinity();

			helper::extractMorphologies(spatialGraphSetFilename,CellType,outputpath,bounds);
		}
		else if (CellTypeStr.compare("L4ssIN") == 0)
		{
			const char * spatialGraphSetFilename = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/Realization20151007_ascii/Morphologies.am";
			int CellType = L4ss;
			const char * outputpath = "/nas1/Data_daniel/Network/L4INSlice/L4ssDendrites/";
			double bounds[6];
			bounds[0] = -150;
			bounds[1] = 150;
			bounds[2] = -400;
			bounds[3] = 400;
			bounds[4] = -(std::numeric_limits<double>::infinity());
			bounds[5] = std::numeric_limits<double>::infinity();

			helper::extractMorphologies(spatialGraphSetFilename,CellType,outputpath,bounds);
		}
		else
		{
			std::cout << CellTypeStr << " not implemented yet!" << std::endl;
		}
	}
	else
	{
		std::cout << "Invalid number of inputs! " << std::endl;
	}
	return 0;
}
