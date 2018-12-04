/****************************************************************************/
/*                                                                          */
/* File:      cell_count.cpp                                                */
/*                                                                          */
/* Purpose:   this program seeks to extract the exact number ond position   */
/*            of neuronal cellbodies labelled with a fluorescent dye        */
/*            from confocal or 2-Photon stacks                              */
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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "typedefs.h"
#include "bloodvesselpattern.h"
/*old Values*/
float XYSAMPLING = 0;
float ZSAMPLING = 0;

/*float XYSAMPLING = 0.311;
float ZSAMPLING = 0.6;

float XYSAMPLING = 0.395879;
float ZSAMPLING = 1;*/

float averageSomaRadius = 6.0;
float zScale = 1.5;
unsigned long BINSIZE = 0;
unsigned long BOXSIZE = 150;

int WriteSpatialGraphFile(const char * output_file, AmiraContourGraphType * amira_bvp_graph);
AmiraContourGraphType * BloodVesselPatternToAmiraSpatialGraph(BloodVesselPattern * bvPattern);
std::list<std::vector<float> > SetContourVertices( std::list<std::list<std::vector<float> > > contour_edges);

int main( int argc , char * argv[])
{
	XYSAMPLING = 0.36075;
	ZSAMPLING = 1.0;
	if(argc == 3 || argc == 4 || argc == 5)
	{
		const char * inputFilename = argv[1];
		const char * outputFilename = argv[2];
		int threshold = -1;
		if(argc >= 4)
			XYSAMPLING = atof(argv[3]);
		if(argc == 5)
			threshold = atoi(argv[4]);
// 		
// 		// 		unsigned int start = atoi(argv[3]);
// 		// 		unsigned int stop  = atoi(argv[4]);
// 		unsigned int start = 1;
// 		unsigned int stop  = 1;
// 		int sliceNo = atoi(argv[3]);
// 		
// 		// 		int int1 = atoi(argv[5]);
// 		// 		int int2 = atoi(argv[6]);
// 		
// 		// 		float arg1 = atof(argv[7]);
// 		// 		float arg2 = atof(argv[8]);
// 		// 		float arg3 = atof(argv[9]);
// 		float arg1 = atof(argv[4]);
// 		float arg2 = atof(argv[5]);
// 		float arg3 = atof(argv[6]);
		
		BloodVesselPattern * vesselPattern = new BloodVesselPattern(inputFilename, outputFilename);
		vesselPattern->ExtractBloodVesselPattern(threshold);
		AmiraContourGraphType * amiraVesselPattern = BloodVesselPatternToAmiraSpatialGraph(vesselPattern);
		WriteSpatialGraphFile(outputFilename, amiraVesselPattern);
	}
	//getchar();
	return 0;
}

/************************************************************************************/
/*WriteSpatialGraphFile(spatial_graph, filename) writes an "am"-file                */
/************************************************************************************/

int WriteSpatialGraphFile(const char * output_file, AmiraContourGraphType * amira_bvp_graph)
{
// 	std::list<std::list<Compartment * > >::iterator edge_list_it;
// 	std::list<Compartment * >::iterator edge_it;
	std::list<std::vector<float> >::iterator contour_it;
	std::list<std::list<std::vector<float> > >::iterator edge_list_contour_it;
	
	int number_of_edge_points = 0;
	
// 	for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
// 	{
// 		number_of_edge_points += (*edge_list_it).size();
// 	}
// 	for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
// 	{
// 		number_of_edge_points += (*edge_list_contour_it).size();
// 	}
	for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
	{
		number_of_edge_points += (*edge_list_contour_it).size();
	}
	
	std::string format = output_file;
	format += ".am";
	
	#ifdef DEBUG
	std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
	//std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
	#endif
	
	std::ofstream NeuroMorphData( format.c_str() );
	
	NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	NeuroMorphData << "# This SpatialGraph file was created by the Neuron Registration Tool NeuroMap " << std::endl;
	NeuroMorphData << "# NeuroMap was programmed by Robert Egger," << std::endl;
	NeuroMorphData << "# Max-Planck-Florida Institute, Jupiter, Florida " << std::endl;
	
	NeuroMorphData << "define VERTEX " << amira_bvp_graph->vertice_list.size() << std::endl;
	NeuroMorphData << "define EDGE " << amira_bvp_graph->edge_list.size()  << std::endl;
	NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
	
	NeuroMorphData << "Parameters {GraphLabels {"                                	<<std::endl;
	NeuroMorphData << "        Neuron { "                                	<<std::endl;
	NeuroMorphData << "            Dendrite {"                           	<<std::endl;
	NeuroMorphData << "                ApicalDendrite {"                 	<<std::endl;
	NeuroMorphData << "                    Color 1 0.5 0.5,"          	<<std::endl;
	NeuroMorphData << "                    Id 4 }"                     	<<std::endl;
	NeuroMorphData << "                BasalDendrite {"         		<<std::endl;
	NeuroMorphData << "                    Color 0.8 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                    Id 5 }"				<<std::endl;
	NeuroMorphData << "                Color 1 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 3 }"				<<std::endl;
	NeuroMorphData << "            Axon {"					<<std::endl;
	NeuroMorphData << "                Color 0 0 1,"			<<std::endl;
	NeuroMorphData << "                Id 6 }"				<<std::endl;
	NeuroMorphData << "            Soma {"					<<std::endl;
	NeuroMorphData << "                Color 1 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 7 }"				<<std::endl;
	NeuroMorphData << "            Color 1 0 0,"				<<std::endl;
	NeuroMorphData << "            Id 2 }"					<<std::endl;
	NeuroMorphData << "        Landmark {"					<<std::endl;
	NeuroMorphData << "            Pia {"					<<std::endl;
	NeuroMorphData << "                Color 0 1 0.5,"			<<std::endl;
	NeuroMorphData << "                Id 9 }"				<<std::endl;
	NeuroMorphData << "            Vessel {"				<<std::endl;
	NeuroMorphData << "                Color 1 0.5 0,"			<<std::endl;
	NeuroMorphData << "                Id 10 }"				<<std::endl;
	NeuroMorphData << "            Barrel {"				<<std::endl;
	NeuroMorphData << "                aRow {"				<<std::endl;
	NeuroMorphData << "                    A1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 13 }"			<<std::endl;
	NeuroMorphData << "                    A2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 14 }"			<<std::endl;
	NeuroMorphData << "                    A3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 15 }"			<<std::endl;
	NeuroMorphData << "                    A4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 16 }"			<<std::endl;
	NeuroMorphData << "                Color 1 0.2 0.2,"			<<std::endl;
	NeuroMorphData << "                Id 12 }"				<<std::endl;
	NeuroMorphData << "                bRow {"				<<std::endl;
	NeuroMorphData << "                    B1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 18 }"			<<std::endl;
	NeuroMorphData << "                    B2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 19 }"			<<std::endl;
	NeuroMorphData << "                    B3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 20 }"			<<std::endl;
	NeuroMorphData << "                    B4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 21 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                    Id 17 }"				<<std::endl;
	NeuroMorphData << "                cRow {"				<<std::endl;
	NeuroMorphData << "                    C1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 23 }"			<<std::endl;
	NeuroMorphData << "                    C2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 24 }"			<<std::endl;
	NeuroMorphData << "                    C3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 25 }"			<<std::endl;
	NeuroMorphData << "                    C4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 26 }"			<<std::endl;
	NeuroMorphData << "                    C5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 27 }"			<<std::endl;
	NeuroMorphData << "                    C6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 28 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                    Id 22 }"				<<std::endl;
	NeuroMorphData << "                dRow {"				<<std::endl;
	NeuroMorphData << "                    D1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 30 }"			<<std::endl;
	NeuroMorphData << "                    D2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 31 }"			<<std::endl;
	NeuroMorphData << "                    D3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 32 }"			<<std::endl;
	NeuroMorphData << "                    D4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 33 }"			<<std::endl;
	NeuroMorphData << "                    D5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 34 }"			<<std::endl;
	NeuroMorphData << "                    D6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 35 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                    Id 29 }"				<<std::endl;
	NeuroMorphData << "                eRow {"				<<std::endl;
	NeuroMorphData << "                    E1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 37 }"			<<std::endl;
	NeuroMorphData << "                    E2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 38 }"			<<std::endl;
	NeuroMorphData << "                    E3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 39 }"			<<std::endl;
	NeuroMorphData << "                    E4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 40 }"			<<std::endl;
	NeuroMorphData << "                    E5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 41 }"			<<std::endl;
	NeuroMorphData << "                    E6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 42 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                    Id 36 }"				<<std::endl;
	NeuroMorphData << "                greekRow {"				<<std::endl;
	NeuroMorphData << "                    Alpha {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 44 }"			<<std::endl;
	NeuroMorphData << "                    Beta {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 45 }"			<<std::endl;
	NeuroMorphData << "                    Gamma {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 46 }"			<<std::endl;
	NeuroMorphData << "                    Delta {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 47 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                    Id 43 }"				<<std::endl;
	NeuroMorphData << "                Color 0 1 0,"			<<std::endl;
	NeuroMorphData << "                Id 11 }"				<<std::endl;
	NeuroMorphData << "            WhiteMatter {"				<<std::endl;
	NeuroMorphData << "                Color 0.5 1 0.75,"			<<std::endl;
	NeuroMorphData << "                Id 48 }"				<<std::endl;
	NeuroMorphData << "            OtherBarrels {"				<<std::endl;
	NeuroMorphData << "                Color 1 0 1,"			<<std::endl;
	NeuroMorphData << "                Id 49 }"				<<std::endl;
	NeuroMorphData << "            ZAxis {"					<<std::endl;
	NeuroMorphData << "                Color 0 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 50 }"				<<std::endl;
	NeuroMorphData << "            Color 0 1 1,"				<<std::endl;
	NeuroMorphData << "            Id 8 }"					<<std::endl;
	NeuroMorphData << "        Id 0,"					<<std::endl;
	NeuroMorphData << "        Color 0 0 0 }"				<<std::endl;
	NeuroMorphData << "ContentType \"HxSpatialGraph\" }"                 	<<std::endl;
	
	NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 " 		<< std::endl;
	NeuroMorphData << "VERTEX {int GraphLabels } @2 " 			<< std::endl;
	
	NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 " 		<< std::endl;
	NeuroMorphData << "EDGE { int NumEdgePoints } @4 " 			<< std::endl;
	NeuroMorphData << "EDGE { int GraphLabels } @5 " 			<< std::endl;
	
	NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 " 	<< std::endl;
	NeuroMorphData << "POINT { float Radius } @7 " 				<< std::endl;
	
	NeuroMorphData << "@1 # Vertices xyz coordinates" << std::endl;
	for(contour_it=amira_bvp_graph->vertice_list.begin(); contour_it!=amira_bvp_graph->vertice_list.end(); ++contour_it)
		NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD]  << " " << (*contour_it)[Z_COORD]  << std::endl;
	
	NeuroMorphData << "\n@2 # Vertex Graph Label" << std::endl;
	for(int i = 0; i < amira_bvp_graph->vertice_list.size(); i++)
	{
		NeuroMorphData << 10 << std::endl;
	}
	
	NeuroMorphData << "@3 # Edge Identifiers" << std::endl;
	int last_index = 0;
	for(int i=0; i < amira_bvp_graph->vertice_list.size(); i++)
	{
		NeuroMorphData << last_index+i << " " << last_index+i << std::endl;
	}
	
	NeuroMorphData << "@4 # Number of Points per Edge" << std::endl;
	for(edge_list_contour_it=amira_bvp_graph->edge_list.begin(); edge_list_contour_it!=amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
	{
		int number_of_contur_points = 0;
		
		number_of_contur_points = (*edge_list_contour_it).size();
		
		NeuroMorphData <<  number_of_contur_points <<std::endl;
	}
	
	NeuroMorphData << "\n@5 # Edge Graph Labels" << std::endl;
	for(edge_list_contour_it=amira_bvp_graph->edge_list.begin(); edge_list_contour_it!=amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
	{
		NeuroMorphData << 10 << std::endl;
	}
	
	NeuroMorphData << "@6 # Point xyz coordinates" << std::endl;
	for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
	{	  
		for(contour_it = (*edge_list_contour_it).begin(); contour_it != (*edge_list_contour_it).end(); ++contour_it)
		{
			NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
		}
		
	}
	
	NeuroMorphData << "@7 # Radius at Point" << std::endl;
	for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
	{	  
		for(contour_it = (*edge_list_contour_it).begin(); contour_it != (*edge_list_contour_it).end(); ++contour_it)
		{
			NeuroMorphData << 0 << std::endl;
		}
		
	}
	
	return 0;
};

/************************************************************************************/
/*BloodVesselPatternToAmiraSpatialGraph()                                           */
/************************************************************************************/
AmiraContourGraphType * BloodVesselPatternToAmiraSpatialGraph(BloodVesselPattern * bvPattern)
{
	std::list<std::list<std::vector<float> > > edges = bvPattern->GetAmiraEdges();
	std::list<std::vector<float> > vertice_list = SetContourVertices(edges);
	
	AmiraContourGraphType *tmp_contour = new AmiraContourGraphType(vertice_list, edges);
	return tmp_contour;
};

/************************************************************************************/
/*SetContourVertices(contour_edges)                                            */
/************************************************************************************/
std::list<std::vector<float> > SetContourVertices( std::list<std::list<std::vector<float> > > contour_edges)
{
	std::list<std::vector<float> > vertices;
	
	std::list<std::list<std::vector<float> > >::iterator edge_list_it;
	std::list<std::vector<float> >::iterator edge_it;
	
	for(edge_list_it = contour_edges.begin(); edge_list_it != contour_edges.end(); ++edge_list_it)
	{	  
		edge_it = (*edge_list_it).begin();
		vertices.push_back(*edge_it);
		#ifdef DEBUG
		for(edge_it = (*edge_list_it).begin(); edge_it != (*edge_list_it).end(); ++edge_it)
		{
			std::cout<< "Vessel: " << (*edge_it)[0] << " " << (*edge_it)[1] << " " << (*edge_it)[2] << std::endl;
		}
		#endif
	}
	
	return vertices;
};
