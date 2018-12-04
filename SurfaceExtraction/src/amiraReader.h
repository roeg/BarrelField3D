/****************************************************************************/
/*                                                                          */
/* File:      amiraReader.h										                */
/*                                                                          */
/* Purpose:   class for basic functions and	common variables				*/
/*			  especially writing and reading functions						*/
/*                                                                          */
/* Author:    Marcel Oberlaender                                            */
/*            Max-Planck-Institute of Neurobiology                          */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/* Co-Author: Stefan Reissl                                                 */
/*            Max-Planck-Institute of Neurobiology                          */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/* Co-Author: Robert Egger                                                  */
/*            Max-Planck-Institute for Medical Research                     */
/*            Jahnstrasse 19                                                */
/*            D-69120 Heidelberg                                            */
/*                                                                          */
/* EMail:     oberlaender@neuro.mpg.de                                      */
/*                                                                          */
/* History:   03.01.2010                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#pragma once
#include "typedefs.h"
#include "basics.h"
#include "spatialGraph.h"

#ifndef AMIRAREADER
#define AMIRAREADER


class Reader : public Basics
{
	public:
		Reader(const char * filename);
		Reader(const char * filename, const char * outputFilename);
		~Reader();
		
		void readSpatialGraphFile();
		void writeSpatialGraphFile();
		
		void getTransformationsFromSpatialGraphFile();
		
		void writeAmiraSurfaceFile(PolyDataPointerType triangleData);
		PolyDataPointerType readAmiraSurfaceFile();
		
		AmiraSpatialGraph * getSpatialGraph() { return inputSpatialGraph; }
		void setSpatialGraph(AmiraSpatialGraph * outputSpatialGraph) { inputSpatialGraph = outputSpatialGraph; }
		
	private:
		std::string letters;
		std::string numbers;
		std::string signs;
		std::string otherChars;
		std::string whitespace;
		
		AmiraSpatialGraph * inputSpatialGraph;
// 		const char * inputFilename;						// filename of the input images
// 		const char * outputFilename;					// filename of CellClusterList.csv, image files (*.sr), landmark file and histo stuff 
		
		//******** basic input and output functions ********
		
};



#endif