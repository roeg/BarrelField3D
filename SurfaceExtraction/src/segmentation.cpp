/****************************************************************************/
/*                                                                          */
/* File:      segmentation.cpp 						    */
/*                                                                          */
/* Purpose:   class for segmentation of either Pia/WM images (4x air)	    */
/*            or barrel field images(40x oil)                               */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*                                                                          */
/*                                                                          */
/* EMail: Robert.Egger@maxplanckflorida.org                                 */
/*                                                                          */
/* History:   22.12.2010                                                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/
#include "segmentation.h"

// #define DEBUG
// #define PIPELINE_DOC

AmiraContourGraphType::AmiraContourGraphStruct(std::list<std::vector<float> > _vertice_list, std::list<std::list<std::vector<float> > >  _edge_list):vertice_list(_vertice_list),edge_list(_edge_list)
{};
AmiraContourGraphType::AmiraContourGraphStruct()
{};

Segmentation::Segmentation(int start, int stop, const char * inputFilename, const char * outputFilename, const char * medProjFilename, int samplingRate, int slice, float arg1, float arg2, float arg3)
{
	originalImage = ImageType::New();
	inputImage = ImageType::New();
	medProjImage = ImageType::New();
	
	this->inputFilename = inputFilename;
	this->outputFilename = outputFilename;
	this->start = start;
	this->stop = stop;
	readImage(start, stop, inputFilename);
	inputImage->Update();
	maximum_region = inputImage->GetLargestPossibleRegion();
	
	originalImage->SetRegions(maximum_region);
	originalImage->Allocate();
	originalImage->FillBuffer(0);
	originalImage->Update();
	copy(originalImage, inputImage);
	
	readImage(1, 1, medProjFilename);
	medProjImage = inputImage;
	medProjImage->Update();
	
// 	read2DImage(inputFilename);
	this->XYDownSamplingRate = samplingRate;
	this->sliceNumber = slice;
	
	this->factor = arg1;
	this->factor2 = arg2;
	this->factor3 = arg3;
	#ifdef DEBUG
	std::string debugFilename(outputFilename);
	debugFilename += "_segmentation_debug.log";
	DebugLog.open(debugFilename.c_str());
	#endif
};

Segmentation::Segmentation(int start, int stop, const char * inputFilename, const char * outputFilename,float factor, float factor2, float factor3, int slice)
{
	originalImage = ImageType::New();
	inputImage = ImageType::New();
	
	this->inputFilename = inputFilename;
	this->outputFilename = outputFilename;
	this->start = start;
	this->stop = stop;
	readImage(start, stop, inputFilename);
	inputImage->Update();
	maximum_region = inputImage->GetLargestPossibleRegion();
	
	originalImage->SetRegions(maximum_region);
	originalImage->Allocate();
	originalImage->FillBuffer(0);
	originalImage->Update();
	copy(originalImage, inputImage);
	// 	read2DImage(inputFilename);
	this->factor = factor;
	this->factor2 = factor2;
	this->factor3 = factor3;
	this->sliceNumber = slice;
	this->XYDownSamplingRate = DOWNSAMPLING;
};

Segmentation::Segmentation( int start, int stop, int samplingRate, ImageType::Pointer segmentedStack, ImageType::Pointer preprocStack, BloodVesselPattern * vesselPattern, std::vector< std::vector< Contour * > > barrelContours, const char * outputFilename)
{
	this->start = start;
	this->stop = stop;
	this->XYDownSamplingRate = samplingRate;
	this->outputFilename = outputFilename;
	inputImage = ImageType::New();
	preprocImage = ImageType::New();
	segmentedImageStack = ImageType::New();
	preprocImageStack = ImageType::New();
	segmentedImageStack = segmentedStack;
	preprocImageStack = preprocStack;
	barrelZContours = barrelContours;
// 	bvPattern = new BloodVesselPattern;
	bvPattern = vesselPattern;
// 	std::flush(std::cout << "converting bvp to spatial graph..." << std::endl);
	BloodVesselPatternToAmiraSpatialGraph();
	computeBounds(segmentedImageStack);
	#ifdef DEBUG
	std::string debugFilename(outputFilename);
	debugFilename += "_segmentation_debug2.log";
	DebugLog.open(debugFilename.c_str());
	#endif
};

Segmentation::~Segmentation()
{//to be filled
	#ifdef DEBUG
	DebugLog.close();
	#endif
};

void Segmentation::piaWMExtraction(int WMbegin, bool SegmentPia)
{
	std:: cout << "Processing " << inputFilename << std::endl;
	if(SegmentPia)
	{
		std::flush(std::cout << "Calculating Pia contour...");
		GlobalThreshold(150);
		std::cout << "done!" << std::endl;
	}
	else if(sliceNumber < WMbegin)
	{
		amira_contour_graph = new AmiraContourGraphStruct();
	}

	bool WMSlice = 0;
	if(sliceNumber >= WMbegin)
	{
		std::flush(std::cout << "Calculating WM contour...");
		if(!SegmentPia)
			amira_contour_graph = new AmiraContourGraphStruct();
		getWMContour(1.0, 1.0);
		WMSlice = 1;
		std::cout << "done!" << std::endl;
	}
	
// // 	writeInputImage("_threshold.bmp");
	std::flush(std::cout << "Calculating BVP...");
	bvPattern = ExtractBloodVesselPattern(1, WMSlice);
	std::cout << "done!" << std::endl;
	BloodVesselPatternToAmiraSpatialGraph();
	WriteSpatialGraphFile(outputFilename);
};

void Segmentation::barrelFieldBvp()
{
	std:: cout << "Processing " << inputFilename << std::endl;
	std::flush(std::cout << "Calculating BVP...");
// 	inputImage = medProjImage;
// 	inputImage->Update();
	bvPattern = ExtractBloodVesselPattern(2, 0);
	std::cout << "done!" << std::endl;
	BloodVesselPatternToAmiraSpatialGraph();
	WriteBVPSpatialGraphFile(outputFilename);
};

void Segmentation::barrelExtraction(bool preProcessed)
{
	std:: cout << "Processing " << inputFilename << std::endl;
	std::flush(std::cout << "Calculating BVP...");
	inputImage = medProjImage;
	inputImage->Update();
	bvPattern = ExtractBloodVesselPattern(2, 0);
	std::cout << "done!" << std::endl;
	BloodVesselPatternToAmiraSpatialGraph();
	
	preprocImageStack = ImageType::New();
	preprocImageStack->SetRegions(originalImage->GetLargestPossibleRegion());
	preprocImageStack->Allocate();
	segmentedImageStack = ImageType::New();
	segmentedImageStack->SetRegions(originalImage->GetLargestPossibleRegion());
	segmentedImageStack->Allocate();
	
	ImageType::RegionType planeRegion;
	ImageType::RegionType planeImageRegion;
	ImageType::IndexType planeImageIndex;
	ImageType::IndexType planeIndex;
	ImageType::SizeType planeSize;
	planeImageIndex = originalImage->GetLargestPossibleRegion().GetIndex();
	planeIndex = originalImage->GetLargestPossibleRegion().GetIndex();
	planeSize = originalImage->GetLargestPossibleRegion().GetSize();
	planeSize[2] = 1;
	planeImageRegion.SetIndex(planeImageIndex);
	planeImageRegion.SetSize(planeSize);
	inputImage->SetRegions(planeImageRegion);
	inputImage->Allocate();
	preprocImage = ImageType::New();
	preprocImage->SetRegions(planeImageRegion);
	preprocImage->Allocate();
	fgRemovedImage = ImageType::New();
	fgRemovedImage->SetRegions(planeImageRegion);
	fgRemovedImage->Allocate();
	sigmoidImage = ImageType::New();
	sigmoidImage->SetRegions(planeImageRegion);
	sigmoidImage->Allocate();
	segmentedImage = ImageType::New();
	segmentedImage->SetRegions(planeImageRegion);
	segmentedImage->Allocate();
	std::flush(std::cout << "Segmenting barrels...");
	for(int z = 0; z < (stop - start + 1); ++z)
	{
		planeIndex[2] = z;
		planeRegion.SetIndex(planeIndex);
		planeRegion.SetSize(planeSize);
		copy(inputImage, originalImage, planeRegion);
		maximum_region = inputImage->GetLargestPossibleRegion();
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			barrelZContours[z][ii]->addAttribute(z);
		segmentBarrels(barrelZContours[z].size(), barrelZContours[z], z, preProcessed);
	}
	std::cout << "done!" << std::endl;
};

/************************************************************************************/
/*BloodVesselPatternToAmiraSpatialGraph()                                           */
/************************************************************************************/
void Segmentation::BloodVesselPatternToAmiraSpatialGraph()
{
	std::list<std::list<std::vector<float> > > edges = bvPattern->GetAmiraEdges();
	std::list<std::vector<float> > vertice_list = SetContourVertices(edges);
	
	AmiraContourGraphType *tmp_contour = new AmiraContourGraphType(vertice_list, edges);
	amira_bvp_graph = tmp_contour;
};

/************************************************************************************/
/*SetContourVertices(contour_edges)                                            */
/************************************************************************************/
std::list<std::vector<float> > Segmentation::SetContourVertices( std::list<std::list<std::vector<float> > > contour_edges)
{
	std::list<std::vector<float> > vertices;
	
	std::list<std::list<std::vector<float> > >::iterator edge_list_it;
	std::list<std::vector<float> >::iterator edge_it;
	
	for(edge_list_it = contour_edges.begin(); edge_list_it != contour_edges.end(); ++edge_list_it)
	{	  
		edge_it = (*edge_list_it).begin();
		if(edge_it != edge_list_it->end())
			vertices.push_back(*edge_it);
// 		#ifdef DEBUG
// 		for(edge_it = (*edge_list_it).begin(); edge_it != (*edge_list_it).end(); ++edge_it)
// 		{
// 			std::cout<< "Vessel: " << (*edge_it)[0] << " " << (*edge_it)[1] << " " << (*edge_it)[2] << std::endl;
// 		}
// 		#endif
	}
	
	return vertices;
};

/************************************************************************************/
/*WriteSpatialGraphFile(spatial_graph, filename) writes an "am"-file                */
/************************************************************************************/
int Segmentation::WriteSpatialGraphFile(const char * output_file)
{
	std::list<std::vector<float> >::iterator contour_it;
	std::list<std::list<std::vector<float> > >::iterator edge_list_contour_it;
	
	int number_of_edge_points = 0;
	int number_of_contours = 0;
	bool contours_empty = amira_contour_graph->vertice_list.empty();
	for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
	{
		number_of_edge_points += (*edge_list_contour_it).size();
		if((*edge_list_contour_it).size())
			++number_of_contours;
	}
	for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
	{
		number_of_edge_points += (*edge_list_contour_it).size();
	}
	
	std::string format = output_file;
	format += ".am";
	
// 	#ifdef DEBUG
// 	std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
// 	//std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
// 	#endif
	
	std::ofstream NeuroMorphData( format.c_str() );
	
	NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	NeuroMorphData << "# This SpatialGraph file was created by the Neuron Reconstruction Tool NeuroMorph " << std::endl;
	NeuroMorphData << "# NeuroMorph was programmed by Marcel Oberlaender and Philip J. Broser," << std::endl;
	NeuroMorphData << "# Max-Planck-Institute for Medical Research Heidelberg, Germany " << std::endl;
	
	NeuroMorphData << "define VERTEX " << amira_contour_graph->vertice_list.size() + amira_bvp_graph->vertice_list.size() << std::endl;
	NeuroMorphData << "define EDGE " << number_of_contours + amira_bvp_graph->edge_list.size() << std::endl;
	NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
// 	NeuroMorphData << "nVertices " << amira_contour_graph->vertice_list.size() + amira_bvp_graph->vertice_list.size() << std::endl;
// 	NeuroMorphData << "nEdges " <<  number_of_contours + amira_bvp_graph->edge_list.size()  << std::endl;
// 	NeuroMorphData << "define Graph " <<  amira_contour_graph->vertice_list.size() + amira_bvp_graph->vertice_list.size() + number_of_contours +  amira_bvp_graph->edge_list.size() << std::endl;
// 	NeuroMorphData << "define EdgePoints " << number_of_edge_points << std::endl;
	
	NeuroMorphData << "Parameters {GraphLabels {"                        	<<std::endl;
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
	NeuroMorphData << "            Color 0 1 1,"				<<std::endl;
	NeuroMorphData << "            Id 8 }"					<<std::endl;
	NeuroMorphData << "        Id 0,"					<<std::endl;
	NeuroMorphData << "        Color 0 0 0 }"				<<std::endl;
	NeuroMorphData << "ContentType \"HxSpatialGraph\" }"      	        <<std::endl;
	
	NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 " 		<< std::endl;
	NeuroMorphData << "VERTEX {int GraphLabels } @2 " 			<< std::endl;
	
	NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 " 		<< std::endl;
	NeuroMorphData << "EDGE { int NumEdgePoints } @4 " 			<< std::endl;
	NeuroMorphData << "EDGE { int GraphLabels } @5 " 			<< std::endl;
	
	NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 " 	<< std::endl;
	NeuroMorphData << "POINT { float Radius } @7 " 				<< std::endl;
	
	NeuroMorphData << "@1 # Vertices xyz coordinates" 			<< std::endl;
	for(contour_it=amira_contour_graph->vertice_list.begin(); contour_it!=amira_contour_graph->vertice_list.end(); ++contour_it)
		NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD]  << " " << (*contour_it)[Z_COORD]  << std::endl;
	for(contour_it=amira_bvp_graph->vertice_list.begin(); contour_it!=amira_bvp_graph->vertice_list.end(); ++contour_it)
		NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD]  << " " << (*contour_it)[Z_COORD]  << std::endl;
	
	NeuroMorphData << "@2 # Vertex Graph Labels" << std::endl;
	for(contour_it=amira_contour_graph->vertice_list.begin(); contour_it!=amira_contour_graph->vertice_list.end(); ++contour_it)
		NeuroMorphData << Pia << std::endl;
	for(contour_it=amira_bvp_graph->vertice_list.begin(); contour_it!=amira_bvp_graph->vertice_list.end(); ++contour_it)
		NeuroMorphData << Vessel << std::endl;
	
	NeuroMorphData << "@3 # Edge Identifiers" << std::endl;
	int last_index = amira_contour_graph->vertice_list.size();
	for(int i=0; i < amira_contour_graph->vertice_list.size(); i++)
	{
		NeuroMorphData << /*last_index+i+1*/ i << " " << /*last_index+i+1*/ i << std::endl;
	}
	for(int i=0; i < amira_bvp_graph->vertice_list.size(); i++)
	{
		NeuroMorphData << last_index+i << " " << last_index+i << std::endl;
	}
	
	NeuroMorphData << "@4 # Number of Points per Edge" << std::endl;
	for(edge_list_contour_it=amira_contour_graph->edge_list.begin(); edge_list_contour_it!=amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
	{
		int number_of_contur_points = 0;
		
		number_of_contur_points = (*edge_list_contour_it).size();
		
		NeuroMorphData <<  number_of_contur_points <<std::endl;
	}
	for(edge_list_contour_it=amira_bvp_graph->edge_list.begin(); edge_list_contour_it!=amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
	{
		int number_of_contur_points = 0;
		
		number_of_contur_points = (*edge_list_contour_it).size();
		
		NeuroMorphData <<  number_of_contur_points <<std::endl;
	}
	
	NeuroMorphData << "@5 # Edge Graph Labels" << std::endl;
	for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
		if(edge_list_contour_it->size())
			NeuroMorphData << Pia << std::endl;
	for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
		NeuroMorphData << Vessel << std::endl;
	
	NeuroMorphData << "@6 # Point xyz coordinates" << std::endl;
	for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
	{	  
		for(contour_it = (*edge_list_contour_it).begin(); contour_it != (*edge_list_contour_it).end(); ++contour_it)
		{
			NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
		}
		
	}
	for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
	{	  
		for(contour_it = (*edge_list_contour_it).begin(); contour_it != (*edge_list_contour_it).end(); ++contour_it)
		{
			NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
		}
		
	}
	
	NeuroMorphData << "@7 # Radius at Point" << std::endl;
	for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
	{	  
		for(contour_it = (*edge_list_contour_it).begin(); contour_it != (*edge_list_contour_it).end(); ++contour_it)
		{
			NeuroMorphData << 0 << std::endl;
		}
		
	}
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
/*WriteSpatialGraphFile(spatial_graph, filename) writes an "am"-file                */
/************************************************************************************/
int Segmentation::WriteBVPSpatialGraphFile(const char * output_file)
{
	std::list<std::vector<float> >::iterator contour_it;
	std::list<std::list<std::vector<float> > >::iterator edge_list_contour_it;
	
	int number_of_edge_points = 0;
	for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
	{
		number_of_edge_points += (*edge_list_contour_it).size();
	}
	
	std::string format = output_file;
	format += ".am";
	
	// 	#ifdef DEBUG
	// 	std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
	// 	//std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
	// 	#endif
	
	std::ofstream NeuroMorphData( format.c_str() );
	
	NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	NeuroMorphData << "# This SpatilaGraph file was created by the Neuron Reconstruction Tool NeuroMorph " << std::endl;
	NeuroMorphData << "# NeuroMorph was programmed by Marcel Oberlaender and Philip J. Broser," << std::endl;
	NeuroMorphData << "# Max-Planck-Institute for Medical Research Heidelberg, Germany " << std::endl;
	
	
	NeuroMorphData << "define VERTEX " << amira_bvp_graph->vertice_list.size() << std::endl;
	NeuroMorphData << "define EDGE " << amira_bvp_graph->edge_list.size() << std::endl;
	NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
	
	NeuroMorphData << "Parameters {GraphLabels {"                        	<<std::endl;
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
	NeuroMorphData << "            Color 0 1 1,"				<<std::endl;
	NeuroMorphData << "            Id 8 }"					<<std::endl;
	NeuroMorphData << "        Id 0,"					<<std::endl;
	NeuroMorphData << "        Color 0 0 0 }"				<<std::endl;
	NeuroMorphData << "ContentType \"HxSpatialGraph\" }"      	        <<std::endl;
	
	NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 " 		<< std::endl;
	NeuroMorphData << "VERTEX {int GraphLabels } @2 " 			<< std::endl;
	
	NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 " 		<< std::endl;
	NeuroMorphData << "EDGE { int NumEdgePoints } @4 " 			<< std::endl;
	NeuroMorphData << "EDGE { int GraphLabels } @5 " 			<< std::endl;
	
	NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 " 	<< std::endl;
	NeuroMorphData << "POINT { float Radius } @7 " 				<< std::endl;
	
	NeuroMorphData << "@1 # Vertices xyz coordinates" 			<< std::endl;
	for(contour_it=amira_bvp_graph->vertice_list.begin(); contour_it!=amira_bvp_graph->vertice_list.end(); ++contour_it)
		NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD]  << " " << (*contour_it)[Z_COORD]  << std::endl;
	
	NeuroMorphData << "@2 # Vertex Graph Labels" << std::endl;
	for(contour_it=amira_bvp_graph->vertice_list.begin(); contour_it!=amira_bvp_graph->vertice_list.end(); ++contour_it)
		NeuroMorphData << Vessel << std::endl;
	
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
	
	NeuroMorphData << "@5 # Edge Graph Labels" << std::endl;
	for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
		NeuroMorphData << Vessel << std::endl;
	
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
/*WriteSpatialGraphFile(spatial_graph, filename) writes an "am"-file                */
/*method for only writing barrels direct from contour vector                        */
/************************************************************************************/
int Segmentation::WriteBarrelSpatialGraphFile(const char * output_file)
{
	std::list<std::vector<float> >::iterator contour_it;
	std::list<std::vector<float> >::iterator edge_list_contour_it;
	
	int number_of_edge_points = 0;
	int number_of_contours = 0;
	
	for(int z = 0; z < barrelZContours.size(); ++z)
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			if(barrelZContours[z][ii]->getValid())
			{
				if(barrelZContours[z][ii]->edgeListPointer()->size())
				{
					++number_of_contours;
					number_of_edge_points += barrelZContours[z][ii]->edgeListPointer()->size() + 1;
				}
			}
	
	std::string format = output_file;
	format += ".am";
	
// 	#ifdef DEBUG
// 	std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
// 	//std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
// 	#endif
	
	std::ofstream NeuroMorphData( format.c_str() );
	
	NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	NeuroMorphData << "# This SpatilaGraph file was created by the Neuron Reconstruction Tool NeuroMorph " << std::endl;
	NeuroMorphData << "# NeuroMorph was programmed by Marcel Oberlaender and Philip J. Broser," << std::endl;
	NeuroMorphData << "# Max-Planck-Institute for Medical Research Heidelberg, Germany " << std::endl;
	
	NeuroMorphData << "nVertices " << number_of_contours << std::endl;
	NeuroMorphData << "nEdges " << number_of_contours << std::endl;
	NeuroMorphData << "define Graph " << 2*number_of_contours << std::endl;
	NeuroMorphData << "define EdgePoints " << number_of_edge_points << std::endl;
	
	NeuroMorphData << "Parameters {VertexLabels"<< std::endl;
	NeuroMorphData << "           {UndefinedVertexLabel{Color 1 1 1}"   	<<std::endl;
	NeuroMorphData << "            LowEndingVertexLabel{Color 0 0 1}"   	<<std::endl;
	NeuroMorphData << "            HighEndingVertexLabel{Color 0 1 0}"  	<<std::endl;
	NeuroMorphData << "            IntersecVertexLabel{Color 1 1 0}"     	<<std::endl; 
	NeuroMorphData << "            NormalEndingVertexLabel{Color 1 0 0}}"	<<std::endl; 
	NeuroMorphData << "            EdgeLabels"                           	<<std::endl;
	NeuroMorphData << "           {UndefinedEdgeLabel{Color 1 1 1}"     	<<std::endl;
	NeuroMorphData << "            BottomEdgeLabel{Color 0 0 1}"        	<<std::endl;
	NeuroMorphData << "            TopEdgeLabel{Color 0 1 0}"           	<<std::endl;
	NeuroMorphData << "            IntermediateEdgeLabel{Color 1 1 0}"  	<<std::endl;
	NeuroMorphData << "            TopToBottomEdgeLabel{Color 1 0 0}}"   	<<std::endl;
	NeuroMorphData << "    GraphLabels {"                                	<<std::endl;
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
	NeuroMorphData << "            Color 0 1 1,"				<<std::endl;
	NeuroMorphData << "            Id 8 }"					<<std::endl;
	NeuroMorphData << "        Id 0,"					<<std::endl;
	NeuroMorphData << "        Color 0 0 0 }"				<<std::endl;
	NeuroMorphData << "ContentType \"HxSpatialGraph\" }"      	        <<std::endl;
	
	NeuroMorphData << "Vertices { float[3] VertexCoordinates } = @1 " 	<< std::endl;
	NeuroMorphData << "Vertices {int VertexLabel } = @2 " 			<< std::endl;
	
	NeuroMorphData << "EdgeData { int[2] EdgeConnectivity} = @3 " 		<< std::endl;
	NeuroMorphData << "EdgeData { int NumEdgePoints} = @4 " 		<< std::endl;
	NeuroMorphData << "EdgeData { int EdgeLabel} = @5 " 			<< std::endl;
	
	NeuroMorphData << "EdgePoints { float[3] EdgePointCoordinates} = @6 " 	<< std::endl;
	NeuroMorphData << "EdgePoints { float Radius} = @7 " 			<< std::endl;
	
	NeuroMorphData << "Graph { int GraphLabel} = @8 " 			<< std::endl;
	
	NeuroMorphData << "@1 # Vertices xyz coordinates" 			<< std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			if(barrelZContours[z][ii]->getValid())
				if(barrelZContours[z][ii]->edgeListPointer()->size())
					NeuroMorphData << barrelZContours[z][ii]->edgeListPointer()->front()[X_COORD] << " " << barrelZContours[z][ii]->edgeListPointer()->front()[Y_COORD]  << " " << barrelZContours[z][ii]->edgeListPointer()->front()[Z_COORD]  << std::endl;
	
	NeuroMorphData << "@2 # Vertex Label" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			if(barrelZContours[z][ii]->getValid())
				if(barrelZContours[z][ii]->edgeListPointer()->size())
					NeuroMorphData << 0 << std::endl;
	
	NeuroMorphData << "@3 # Edge Identifiers" << std::endl;
	for(int i=0; i < number_of_contours; i++)
	{
		NeuroMorphData << i << " " << i << std::endl;
	}
	
	NeuroMorphData << "@4 # Number of Points per Edge" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			if(barrelZContours[z][ii]->getValid())
				if(barrelZContours[z][ii]->edgeListPointer()->size())
					NeuroMorphData << barrelZContours[z][ii]->edgeListPointer()->size() + 1 <<std::endl;
	
	NeuroMorphData << "@5 # EdgeLabels" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			if(barrelZContours[z][ii]->getValid())
				if(barrelZContours[z][ii]->edgeListPointer()->size())
					NeuroMorphData << 0 << std::endl;
	
	NeuroMorphData << "@6 # Point xyz coordinates" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			if(barrelZContours[z][ii]->getValid())
				if(barrelZContours[z][ii]->edgeListPointer()->size())
				{
					for(contour_it = barrelZContours[z][ii]->edgeListPointer()->begin(); contour_it != barrelZContours[z][ii]->edgeListPointer()->end(); ++contour_it)
						NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
					NeuroMorphData << barrelZContours[z][ii]->edgeListPointer()->front()[X_COORD] << " " << barrelZContours[z][ii]->edgeListPointer()->front()[Y_COORD]  << " " << barrelZContours[z][ii]->edgeListPointer()->front()[Z_COORD]  << std::endl;
				}
	
	NeuroMorphData << "@7 # Radius at Point" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			if(barrelZContours[z][ii]->getValid())
				if(barrelZContours[z][ii]->edgeListPointer()->size())
				{
					for(int jj = 0; jj < barrelZContours[z][ii]->edgeListPointer()->size() + 1; ++jj)
					{
						// 					NeuroMorphData << 1 << std::endl;
						if(barrelZContours[z][ii]->getOptimizeFlag())
							NeuroMorphData << (*(barrelZContours[z][ii]->attributesPointer()))[10] << std::endl;
						else
							NeuroMorphData << (*(barrelZContours[z][ii]->attributesPointer()))[6] << std::endl;
					}
				}
	
	//first, labels for all vertices
	//then, labels for all edges
	NeuroMorphData << "@8 # Graph Property" << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			if(barrelZContours[z][ii]->getValid())
				if(barrelZContours[z][ii]->edgeListPointer()->size())
				{
					if(barrelZContours[z][ii]->getInsideBarrel())
						NeuroMorphData << Barrel << std::endl;
					else
						NeuroMorphData << 0 << std::endl;
				}
	for(int z = 0; z < barrelZContours.size(); ++z)
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			if(barrelZContours[z][ii]->getValid())
				if(barrelZContours[z][ii]->edgeListPointer()->size())
				{
						if(barrelZContours[z][ii]->getInsideBarrel())
							NeuroMorphData << Barrel << std::endl;
						else
							NeuroMorphData << 0 << std::endl;
				}
	
	return 0;
};

BloodVesselPattern * Segmentation::ExtractBloodVesselPattern(int imageType, bool WM)
{
	
// 	#ifdef DEBUG
// 	std::cout << "Calc Max Image " << std::endl;
// 	#endif
	
// 	int maxZ=   input_image->GetLargestPossibleRegion().GetSize(2);
	BloodVesselPattern * vesselPattern= new BloodVesselPattern();
	ImageToCalcRescaleFilterType::Pointer image2calcFilter = ImageToCalcRescaleFilterType::New();
	
	//Calacualte Max Intensity Projection 
// 	int startZ= maxZ/4; 
// 	int stopZ= maxZ-startZ; //has to be dynamic!
	if(WM)
		copy(inputImage, wmBvpImage);
	else if(!WM && imageType == 1)
	{
		copy(inputImage, originalImage);
		bvpPreparation(1.0, 1.0);
	}
// 	writeInputImage("bvp_fg.tif");
	CalcImageType::Pointer tmpImage = CalcImageType::New();
// 	image2calcFilter->SetInput(originalImage);
	image2calcFilter->SetInput(inputImage);
	image2calcFilter->SetOutputMaximum(255.0);
	image2calcFilter->SetOutputMinimum(0.0);
	image2calcFilter->Update();
	tmpImage = image2calcFilter->GetOutput();
	tmpImage->Update();
	
	tmpImage = invertImage(tmpImage);
	
	CalcImageType::RegionType tmpRegion = tmpImage->GetLargestPossibleRegion();
	CalcImageType::SizeType tmpSize = tmpRegion.GetSize();
	CalcImageType::IndexType tmpIndex = tmpRegion.GetIndex();
	
	CalcImage2DType::Pointer projectionImage = CalcImage2DType::New();
	CalcImage2DType::RegionType calc2DRegion;
	CalcImage2DType::SizeType calc2DSize;
	CalcImage2DType::IndexType calc2DIndex;
	
	calc2DSize[0] = tmpSize[0];
	calc2DSize[1] = tmpSize[1];
	calc2DIndex[0] = tmpIndex[0];
	calc2DIndex[1] = tmpIndex[1];
	
	calc2DRegion.SetIndex(calc2DIndex);
	calc2DRegion.SetSize(calc2DSize);
	
	projectionImage->SetRegions(calc2DRegion);
	projectionImage->Allocate();
	projectionImage->Update();
	
	ConstIteratorType bgIter(inputImage, inputImage->GetLargestPossibleRegion());
	CalcIteratorType originalIter( tmpImage, tmpImage->GetLargestPossibleRegion() );
	CalcIterator2DType projectionIter( projectionImage, projectionImage->GetLargestPossibleRegion() );
	
	for( originalIter.GoToBegin(), projectionIter.GoToBegin()/*, bgIter.GoToBegin()*/ ; !originalIter.IsAtEnd() && !projectionIter.IsAtEnd()/* && !bgIter.IsAtEnd()*/ ; ++originalIter, ++projectionIter/*, ++bgIter*/ )
	{
		projectionIter.Set(originalIter.Get());
	}
	
	projectionImage->Update();
	
	//Blood Vessel Should be next to zero
	
// 	#ifdef DEBUG	
	//Write Temp Data
//  	std::cout << "Create DEBUG Image" << std::endl;
// 	CalcToImage2DRescaleFilterType::Pointer rescaler=CalcToImage2DRescaleFilterType::New();
// 	rescaler->SetOutputMinimum( 0 ); 
// 	rescaler->SetOutputMaximum( 255 ); 
// 	rescaler->SetInput(projectionImage);
// 	
// 	Single2DWriterType::Pointer calcWriter=Single2DWriterType::New();
// 	calcWriter->SetInput(rescaler->GetOutput());
// 	calcWriter->SetFileName(outputFilename);
// 	calcWriter->Update();
// 	#endif
	
	//std::cout << "Detect Vessel " << std::endl;
	vesselPattern->detectBloodVessel(projectionImage, imageType, WM, XYDownSamplingRate);
	
	
// 	#ifdef DEBUG
// 	vesselPattern-> plotBloodVesselPatternToImage(projectionImage);
// 	rescaler->SetInput(projectionImage);
// 	rescaler->Update();
// 	calcWriter->SetInput(rescaler->GetOutput());
// 	calcWriter->SetFileName("bloodvesselafterdetection.tif");
// 	calcWriter->Update();
// 	#endif
	
	
	
	return vesselPattern;
};

/****************************************************************************/
/*bvpPreparation() prepares the image similar to the pre-processing         */
/*for the WM detection algorithm in getWMContour()                          */
/****************************************************************************/
void Segmentation::bvpPreparation( float alpha_arg, float beta_arg)
{
	IteratorType2 inputIter(inputImage, maximum_region);
	
// 	std::cout << "Calculating ROI for BVP..." << std::endl;
	
	for(inputIter.GoToBegin(); !inputIter.IsAtEnd(); ++inputIter)
	{
		if(inputIter.Get() <= 150)
			inputIter.Set(255);
		else
			inputIter.Set(0);
	}
	
	ImageType::Pointer fgImage = MarkBackground(10000);
// 	writeBinaryImage(bgImage, "_simpleBG", 10);
	fgImage = invertImage(fgImage);
	regionGrowingObjectLabeling(fgImage, 1, 1);
	labeledToGreyscaleImage(fgImage, 1);
	binaryErosion(fgImage, 1, 100);
	regionGrowingObjectLabeling(fgImage, 1, 1);
	
// 	std::cout << "Adjusting intensity..." << std::endl;
	
	unsigned long fgpixel = 0, totalpixel = sizeOfObjects[0];
	float mean = 0, stddev = 0, fgfraction = 0;
	copy(inputImage, originalImage);
	
	unsigned int * fghist;
	fghist = foregroundHistogram(fgImage, inputImage);
	
	for(int ii = 0; ii < 256; ++ii)
	{
		fgpixel += fghist[ii];
		mean += (float)ii*(float)fghist[ii];
	}
	mean = mean/(float)fgpixel;
	
	for(int ii = 0; ii < 256; ++ii )
	{
		stddev += (float)fghist[ii]*((float)ii - mean)*((float)ii - mean);
	}
	stddev = stddev/(float)(fgpixel);
	stddev = std::sqrt(stddev);
	
	float alpha = alpha_arg*stddev;
	float beta = mean + beta_arg*stddev;
// 	std::cout << "alpha = " << alpha << "\t---\tbeta = " << beta << std::endl;
	
	IteratorType2 sigmoidIter(inputImage, maximum_region);
	ConstIteratorType fgIter2(fgImage, maximum_region);
	
	for(sigmoidIter.GoToBegin(), fgIter2.GoToBegin(); !sigmoidIter.IsAtEnd() && !fgIter2.IsAtEnd(); ++sigmoidIter, ++fgIter2)
	{
		if(fgIter2.Get())
		{
			float tmp = 255.0f/(1 + std::exp((beta - (float)sigmoidIter.Get())/alpha));
			sigmoidIter.Set((int)(tmp + 0.5));
		}
		else
			sigmoidIter.Set(0);
	}
	
	inputImage->Update();
};

/****************************************************************************/
/*PRIVATE METHODS                                                           */
/****************************************************************************/

/****************************************************************************/
/*Histogram() evokes computation of mean and returns standard-deviation     */
/****************************************************************************/

void Segmentation::Histogram(float * moments)
{
	HistogramSampleType::Pointer histo_sample = HistogramSampleType::New();
	histo_sample->SetMeasurementVectorSize( 1 );
	
	ConstIteratorType histo_it( inputImage, inputImage->GetRequestedRegion() );
	
	for( histo_it.GoToBegin(); !histo_it.IsAtEnd(); ++histo_it)
	{
		MeasurementHistogramType greyvalue;
		greyvalue = histo_it.Get();
		histo_sample->PushBack(greyvalue);
	}
	
	MeanAlgorithmType::Pointer mean_algorithm = MeanAlgorithmType::New();
	
	mean_algorithm->SetInputSample( histo_sample );
	mean_algorithm->Update();
	
	CovarianceAlgorithmType::Pointer covariance_algorithm = CovarianceAlgorithmType::New();
	
	covariance_algorithm->SetInputSample( histo_sample );
	covariance_algorithm->SetMean( mean_algorithm->GetOutput() );
	covariance_algorithm->Update();
	
	float mean = mean_algorithm->GetOutput()->GetElement(0);
	CovarianceAlgorithmType::OutputType covariance = *(covariance_algorithm->GetOutput());
	float variance = *(covariance.operator[](0));
	float standard_deviation = std::sqrt(variance);
	
	//return (mean); //1.5 default
// 	return (mean + (MAGNIFICATION * 1 * standard_deviation)); //1.5 default
	moments[0] = mean;
	moments[1] = standard_deviation;
};

/****************************************************************************/
/*Histogram() evokes computation of mean and returns standard-deviation     */
/****************************************************************************/
void Segmentation::Histogram(float * moments, ImageType::Pointer workingImage)
{
	HistogramSampleType::Pointer histo_sample = HistogramSampleType::New();
	histo_sample->SetMeasurementVectorSize( 1 );
	
	ConstIteratorType histo_it( workingImage, workingImage->GetRequestedRegion() );
	
	for( histo_it.GoToBegin(); !histo_it.IsAtEnd(); ++histo_it)
	{
		MeasurementHistogramType greyvalue;
		greyvalue = histo_it.Get();
		histo_sample->PushBack(greyvalue);
	}
	
	MeanAlgorithmType::Pointer mean_algorithm = MeanAlgorithmType::New();
	
	mean_algorithm->SetInputSample( histo_sample );
	mean_algorithm->Update();
	
	CovarianceAlgorithmType::Pointer covariance_algorithm = CovarianceAlgorithmType::New();
	
	covariance_algorithm->SetInputSample( histo_sample );
	covariance_algorithm->SetMean( mean_algorithm->GetOutput() );
	covariance_algorithm->Update();
	
	float mean = mean_algorithm->GetOutput()->GetElement(0);
	CovarianceAlgorithmType::OutputType covariance = *(covariance_algorithm->GetOutput());
	float variance = *(covariance.operator[](0));
	float standard_deviation = std::sqrt(variance);
	
	//return (mean); //1.5 default
	// 	return (mean + (MAGNIFICATION * 1 * standard_deviation)); //1.5 default
	moments[0] = mean;
	moments[1] = standard_deviation;
};

/****************************************************************************/
/*calculateRegionalHistograms() calculates the histograms                   */
/*of all regions passed in the values vector                                */
/****************************************************************************/
void Segmentation::calculateRegionalHistograms( int** histogram, std::list< std::list< unsigned char > > values, unsigned int divisions )
{
	std::list< std::list< unsigned char > >::iterator regionIter;
	unsigned int region = 0;
	for( regionIter = values.begin(); regionIter != values.end(); ++regionIter, ++region )
	{
		if(regionIter->size())
		{
			
			unsigned long pixelcount = regionIter->size();
			
			for(int i = 0; i < 256; ++i)
			{
				histogram[region][i] = 0;
			}
			
			std::list< unsigned char >::iterator pxValIter;
			
			for( pxValIter = regionIter->begin() ; pxValIter != regionIter->end() ; ++pxValIter )
			{
				++histogram[region][*pxValIter];
			}
// 			#ifdef DEBUG
// 			std::cout << "Region: " << region << " - no of pixels in histogram: " << pixelcount /*<< " - mean: " << mean << " +- " << standard_deviation*/ << std::endl;
// 			#endif
		}
	}
};

/****************************************************************************/
/*calculateRegionalHistograms() calculates the histograms                   */
/*of all regions passed in the values vector                                */
/****************************************************************************/
void Segmentation::calculateRegionalHistograms( int** histogram, std::list< unsigned char > * values, unsigned int divisions )
{
	for(unsigned int region = 0; region < divisions; ++region )
	{
		if(values[region].size())
		{
			
			unsigned long pixelcount = values[region].size();
			
			for(int i = 0; i < 256; ++i)
			{
				histogram[region][i] = 0;
			}
			
			std::list< unsigned char >::iterator pxValIter;
			
			for( pxValIter = values[region].begin() ; pxValIter != values[region].end() ; ++pxValIter )
			{
				++histogram[region][*pxValIter];
			}
// 			#ifdef DEBUG
// 			std::cout << "Region: " << region << " - no of pixels in histogram: " << pixelcount /*<< " - mean: " << mean << " +- " << standard_deviation*/ << std::endl;
// 			#endif
		}
	}
};

/****************************************************************************/
/*calculateRegionalHistograms() calculates the mean and standard deviation  */
/*of all regions passed in the values vector                                */
/****************************************************************************/
void Segmentation::calculateRegionalHistograms( float * means, float * sigmas, std::list< unsigned char > * values, unsigned int divisions )
{
	for( int region = 0 ; region < divisions ; ++region )
	{
		if(values[region].size())
		{
			unsigned long histogram[256];
			unsigned long pixelcount = values[region].size();
			unsigned long tmpmean = 0;
			float tmpsigma = 0;
			float mean = 0, standard_deviation = 0;
			
			for(int i = 0; i < 256; ++i)
			{
				histogram[i] = 0;
			}
			
			std::list< unsigned char >::iterator regionIter;
			
			for( regionIter = values[region].begin() ; regionIter != values[region].end() ; ++regionIter )
			{
				++histogram[*regionIter];
			}
			
			
			for(int i=0; i < 256; ++i)
			{
				tmpmean += i*histogram[i];
			}
			
			mean = (float)tmpmean/(float)pixelcount;
			
			for(int i=0; i < 256; ++i)
			{
				tmpsigma += histogram[i]*((float)i-mean)*((float)i-mean);
			}
			
			standard_deviation = std::sqrt(tmpsigma/(float)pixelcount);
// 			#ifdef DEBUG
// 			std::cout << "Region: " << region << " - no of pixels in histogram: " << pixelcount << " - mean: " << mean << " +- " << standard_deviation << std::endl;
// 			#endif
			means[region] = mean;
			sigmas[region] = standard_deviation;
		}
	}
};

/****************************************************************************/
/*calculateRegionalHistograms() calculates the mean and standard deviation  */
/*of all regions passed in the values vector                                */
/****************************************************************************/
void Segmentation::calculateRegionalHistograms( float * means, float * sigmas, std::vector< unsigned int > * values, unsigned int divisions )
{
	for( int region = 0 ; region < divisions ; ++region )
	{
		if(values[region].size())
		{
			
			unsigned long histogram[256];
			unsigned long pixelcount = values[region].size();
			unsigned long tmpmean = 0;
			float tmpsigma = 0;
			float mean = 0, standard_deviation = 0;
			
			for(int i = 0; i < 256; ++i)
			{
				histogram[i] = 0;
			}
			
			std::vector< unsigned int >::iterator regionIter;
			
			for( regionIter = values[region].begin() ; regionIter != values[region].end() ; ++regionIter )
			{
				++histogram[*regionIter];
			}
			
			
			for(int i=0; i < 256; ++i)
			{
				tmpmean += i*histogram[i];
			}
			
			mean = (float)tmpmean/(float)pixelcount;
			
			for(int i=0; i < 256; ++i)
			{
				tmpsigma += histogram[i]*((float)i-mean)*((float)i-mean);
			}
			
			standard_deviation = std::sqrt(tmpsigma/(float)pixelcount);
// #ifdef DEBUG
// 			std::cout << "Region: " << region << " - no of pixels in histogram: " << pixelcount << " - mean: " << mean << " +- " << standard_deviation << std::endl;
// #endif
			means[region] = mean;
			sigmas[region] = standard_deviation;
		}
	}
};

/****************************************************************************/
/*foregroundHistogram() calculates the histogram of the foreground region   */
/*background is assumed to be set to 0 in bgImage                           */
/****************************************************************************/
unsigned int * Segmentation::foregroundHistogram( ImageType::Pointer bgImage, ImageType::Pointer tmpImage )
{
	unsigned int * histogram = new unsigned int[256];
	
	for(int i = 0; i < 256; ++i)
	{
		histogram[i] = 0;
	}
	
	ConstIteratorType bgIter(bgImage, bgImage->GetLargestPossibleRegion());
	ConstIteratorType foregroundIter(tmpImage, tmpImage->GetLargestPossibleRegion());
	
	for( bgIter.GoToBegin(), foregroundIter.GoToBegin(); !bgIter.IsAtEnd() && !foregroundIter.IsAtEnd(); ++bgIter, ++foregroundIter )
	{
		if(bgIter.Get())
			++histogram[foregroundIter.Get()];
	}
	
	return histogram;
}

/****************************************************************************/
/*foregroundHistogram() calculates the histogram of the maximum region      */
/****************************************************************************/
int * Segmentation::completeHistogram( ImageType::Pointer tmpImage )
{
	int * histogram = new int[256];
	
	for(int i = 0; i < 256; ++i)
	{
		histogram[i] = 0;
	}
	
	ConstIteratorType foregroundIter(tmpImage, tmpImage->GetLargestPossibleRegion());
	for( foregroundIter.GoToBegin(); !foregroundIter.IsAtEnd(); ++foregroundIter )
	{
			++histogram[foregroundIter.Get()];
	}
// 	histogram[255] = 0;
	return histogram;
}

/****************************************************************************/
/*objectHistogram() calculates the histogram of a single object             */
/*background is assumed to be set to 0 and the object the only object       */
/*in the image                                                              */
/****************************************************************************/
unsigned int * Segmentation::objectHistogram( ImageType::Pointer tmpImage )
{
	unsigned int * histogram = new unsigned int[256];
	
	for(int i = 0; i < 256; ++i)
	{
		histogram[i] = 0;
	}
	
	ConstIteratorType foregroundIter(tmpImage, tmpImage->GetLargestPossibleRegion());
	
	for( foregroundIter.GoToBegin(); !foregroundIter.IsAtEnd(); ++foregroundIter )
	{
		if(foregroundIter.Get())
			++histogram[foregroundIter.Get()];
	}
	
	return histogram;
}

void Segmentation::sigmoidIntensityMapping(float a, float b, float p)
{
	unsigned int intensityMap[256];
	float beta = 0, alpha = 0, mean = 0, sigma = 0;
	float * moments =  new float[2];
	
	Histogram(moments, fgRemovedImage);
	mean = moments[0];
	sigma = moments[1];
	
	alpha = a*sigma;
	beta = mean + b*sigma;
	
	for(int ii = 0; ii < 256; ++ii)
	{
		double tmp = 255/pow( 1 + exp((beta - ii)/alpha) , p);
		intensityMap[ii] = (int)(tmp + 0.5);
		if(intensityMap[ii] == 255)
			intensityMap[ii] = 254;		//255 should still be fg marker!
	}
	
	sigmoidImage->FillBuffer(0);
	ConstIteratorType readIter(fgRemovedImage, maximum_region);
	IteratorType2 mapIter(sigmoidImage, maximum_region);
	for(readIter.GoToBegin(), mapIter.GoToBegin(); !readIter.IsAtEnd() && !mapIter.IsAtEnd(); ++readIter, ++mapIter)
		mapIter.Set(intensityMap[readIter.Get()]);
	sigmoidImage->Update();
};

/****************************************************************************/
/*voronoiRegionThreshold() decides foreground/background for bricks         */
/*depending on the local brick statistics.                                  */
/*this is the non-optimizier region growing version                         */
/****************************************************************************/
void Segmentation::voronoiRegionThreshold(unsigned int noOfRegions, std::vector< Contour * > barrelContours, int z)
{
	#ifdef DEBUG
	DebugLog << "Voronoi region growing: preparations" << std::endl;
	#endif
	
	std::list< unsigned char > * pixelValueLists = new std::list< unsigned char >[noOfRegions];
	std::list< unsigned char > * pxCircleList1 = new std::list< unsigned char >[noOfRegions];
	std::list< unsigned char > * pxCircleList2 = new std::list< unsigned char >[noOfRegions];
	std::list< unsigned char > * pxCircleList3 = new std::list< unsigned char >[noOfRegions];
	std::list< unsigned char > * pxOutCircleList1 = new std::list< unsigned char >[noOfRegions];
	std::list< unsigned char > * pxOutCircleList2 = new std::list< unsigned char >[noOfRegions];
	std::list< unsigned char > * pxOutCircleList3 = new std::list< unsigned char >[noOfRegions];
	std::list< ImageType::IndexType > * pixelIndexLists = new std::list< ImageType::IndexType >[noOfRegions];
	
	int ** histogram = new int * [noOfRegions];
	int ** circ1Histogram = new int * [noOfRegions];
	int ** circ2Histogram = new int * [noOfRegions];
	int ** circ3Histogram = new int * [noOfRegions];
	int ** outCirc1Histogram = new int * [noOfRegions];
	int ** outCirc2Histogram = new int * [noOfRegions];
	int ** outCirc3Histogram = new int * [noOfRegions];
	for(int ii = 0; ii < noOfRegions; ++ii)
	{
		histogram[ii] = new int[256];
		circ1Histogram[ii] = new int[256];
		circ2Histogram[ii] = new int[256];
		circ3Histogram[ii] = new int[256];
		outCirc1Histogram[ii] = new int[256];
		outCirc2Histogram[ii] = new int[256];
		outCirc3Histogram[ii] = new int[256];
	}
	unsigned int * threshold = new unsigned int[noOfRegions];
	float * regionMeans = new float[noOfRegions];
	float * circ1Means = new float[noOfRegions];
	float * circ2Means = new float[noOfRegions];
	float * circ3Means = new float[noOfRegions];
	float * outCirc1Means = new float[noOfRegions];
	float * outCirc2Means = new float[noOfRegions];
	float * outCirc3Means = new float[noOfRegions];
	float * regionSigmas = new float[noOfRegions];
	float * circ1Sigmas = new float[noOfRegions];
	float * circ2Sigmas = new float[noOfRegions];
	float * circ3Sigmas = new float[noOfRegions];
	float * outCirc1Sigmas = new float[noOfRegions];
	float * outCirc2Sigmas = new float[noOfRegions];
	float * outCirc3Sigmas = new float[noOfRegions];
	float * SNR = new float[noOfRegions];
	
	for(int i = 0 ; i < noOfRegions ; ++i)
	{
		regionMeans[i] = 0;
		circ1Means[i] = 0, circ2Means[i] = 0, circ3Means[i] = 0;
		outCirc1Means[i] = 0, outCirc2Means[i] = 0, outCirc3Means[i] = 0;
		regionSigmas[i] = 0;
		circ1Sigmas[i] = 0, circ2Sigmas[i] = 0, circ3Sigmas[i] = 0;
		outCirc1Sigmas[i] = 0, outCirc2Sigmas[i] = 0, outCirc3Sigmas[i] = 0;
	}
	
	ConstIteratorType pxValIter(inputImage, maximum_region);
	ConstCalcIteratorType voronoiIter(voronoiMap, maximum_region);
	ConstCalcIteratorType distanceIter(distanceMap, maximum_region);
	for(pxValIter.GoToBegin(), voronoiIter.GoToBegin(), distanceIter.GoToBegin();
	!pxValIter.IsAtEnd() && !voronoiIter.IsAtEnd() && !distanceIter.IsAtEnd();
	++pxValIter, ++voronoiIter, ++distanceIter)
	{
		if(voronoiIter.Get())
			if(pxValIter.Get() != 255)
			{
				int cell = (int)(voronoiIter.Get()+0.5)-1;
				pixelValueLists[cell].push_back(pxValIter.Get());
				pixelIndexLists[cell].push_back(pxValIter.GetIndex());
				float dist = distanceIter.Get();
				if(dist <= 45)
					pxCircleList1[cell].push_back(pxValIter.Get());
				if(dist > 45)
					pxOutCircleList1[cell].push_back(pxValIter.Get());
				if(dist <= 65)
					pxCircleList2[cell].push_back(pxValIter.Get());
				if(dist > 65)
					pxOutCircleList2[cell].push_back(pxValIter.Get());
				if(dist <= 85)
					pxCircleList3[cell].push_back(pxValIter.Get());
				if(dist > 85)
					pxOutCircleList3[cell].push_back(pxValIter.Get());
			}
	}
	
	calculateRegionalHistograms( histogram, pixelValueLists, noOfRegions );
	calculateRegionalHistograms( regionMeans, regionSigmas, pixelValueLists, noOfRegions );
	
	calculateRegionalHistograms(circ1Histogram, pxCircleList1, noOfRegions);
	calculateRegionalHistograms(circ1Means, circ1Sigmas, pxCircleList1, noOfRegions);
	
	calculateRegionalHistograms(circ2Histogram, pxCircleList2, noOfRegions);
	calculateRegionalHistograms(circ2Means, circ2Sigmas, pxCircleList2, noOfRegions);
	
	calculateRegionalHistograms(circ3Histogram, pxCircleList3, noOfRegions);
	calculateRegionalHistograms(circ3Means, circ3Sigmas, pxCircleList3, noOfRegions);
	
	calculateRegionalHistograms(outCirc1Histogram, pxOutCircleList1, noOfRegions);
	calculateRegionalHistograms(outCirc1Means, outCirc1Sigmas, pxOutCircleList1, noOfRegions);
	
	calculateRegionalHistograms(outCirc2Histogram, pxOutCircleList2, noOfRegions);
	calculateRegionalHistograms(outCirc2Means, outCirc2Sigmas, pxOutCircleList2, noOfRegions);
	
	calculateRegionalHistograms(outCirc3Histogram, pxOutCircleList3, noOfRegions);
	calculateRegionalHistograms(outCirc3Means, outCirc3Sigmas, pxOutCircleList3, noOfRegions);
	
	for(int regionID = 0; regionID < noOfRegions; ++regionID)
	{
		if(!barrelContours[regionID]->getValid())
			continue;
		
		threshold[regionID] = otsuThreshold(histogram[regionID]);
		
		barrelContours[regionID]->addAttribute(regionID+1);
		float circ1SNR = outCirc1Means[regionID]/circ1Means[regionID];
		float circ2SNR = outCirc2Means[regionID]/circ2Means[regionID];
		float circ3SNR = outCirc3Means[regionID]/circ3Means[regionID];
		barrelContours[regionID]->addAttribute(circ1SNR);
		barrelContours[regionID]->addAttribute(circ2SNR);
		barrelContours[regionID]->addAttribute(circ3SNR);
		
	}
	
	segmentedImage->FillBuffer(0);
	for(int ii = 0; ii < noOfRegions; ++ii)
	{
		if(!barrelContours[ii]->getValid())
			continue;
		
		std::list< unsigned char > pxValueList2;
		ConstIteratorType pxValIter2(sigmoidImage, maximum_region);
		for(pxValIter2.GoToBegin(), voronoiIter.GoToBegin(), distanceIter.GoToBegin();
		!pxValIter2.IsAtEnd() && !voronoiIter.IsAtEnd() && !distanceIter.IsAtEnd();
		++pxValIter2, ++voronoiIter, ++distanceIter)
		{
			int cell = (int)(voronoiIter.Get()+0.5)-1;
			if(cell == ii && pxValIter2.Get() != 255)
			{
				pxValueList2.push_back(pxValIter2.Get());
			}
		}
		int * sigmoidHisto = new int[256];
		float * sigmoidMean, * sigmoidSigma;
		sigmoidMean = new float, sigmoidSigma = new float;
		calculateRegionalHistograms(&sigmoidHisto, &pxValueList2, 1);
		calculateRegionalHistograms(sigmoidMean, sigmoidSigma, &pxValueList2, 1);
		threshold[ii] = otsuThreshold(sigmoidHisto);
		
		int maxLevel = threshold[ii];
		double minLevel = 0;
		for(int jj = 0; jj < 256; ++jj)
			if(sigmoidHisto[jj])
			{
				minLevel = jj;
				break;
			}
		barrelContours[ii]->addAttribute(maxLevel);
		
		int startLevel = std::max(0.5*maxLevel + 1.0*(*sigmoidMean - maxLevel) + minLevel, minLevel + 10);
		if(startLevel >= maxLevel)
			startLevel = maxLevel - 1;
		
		#ifdef DEBUG
		DebugLog << "Voronoi region growing: region growing" << std::endl;
		#endif
		
		maxLevelNeighborhoodRegionGrowing(segmentedImage, pixelIndexLists[ii], startLevel, maxLevel, 1);
		delete [] sigmoidHisto, delete sigmoidMean, delete sigmoidSigma;
		pxValueList2.clear();
	}
	
	for(int ii = 0; ii < noOfRegions; ++ii)
	{
		pixelValueLists[ii].clear();
		pxCircleList1[ii].clear(), pxCircleList2[ii].clear(), pxCircleList3[ii].clear();
		pxOutCircleList1[ii].clear(), pxOutCircleList2[ii].clear(), pxOutCircleList3[ii].clear();
		pixelIndexLists[ii].clear();
	}
	delete [] pixelValueLists, delete [] pixelIndexLists;
	delete [] pxCircleList1, delete [] pxCircleList2, delete [] pxCircleList3;
	delete [] pxOutCircleList1, delete [] pxOutCircleList2, delete [] pxOutCircleList3;
	for(int ii = 0; ii < noOfRegions; ++ii)
	{
		delete [] histogram[ii];
		delete [] circ1Histogram[ii], delete [] circ2Histogram[ii], delete [] circ3Histogram[ii];
		delete [] outCirc1Histogram[ii], delete [] outCirc2Histogram[ii], delete [] outCirc3Histogram[ii];
	}
	delete [] histogram, delete [] threshold, delete [] regionMeans, delete [] regionSigmas/*, delete [] regionAmplitudes*/, delete [] SNR;
	delete [] circ1Histogram, delete [] circ2Histogram, delete [] circ3Histogram;
	delete [] outCirc1Histogram, delete [] outCirc2Histogram, delete [] outCirc3Histogram;
	delete [] circ1Means, delete [] circ2Means, delete [] circ3Means;
	delete [] outCirc1Means, delete [] outCirc2Means, delete [] outCirc3Means;
	delete [] circ1Sigmas, delete [] circ2Sigmas, delete [] circ3Sigmas;
	delete [] outCirc1Sigmas, delete [] outCirc2Sigmas, delete [] outCirc3Sigmas;
};

/****************************************************************************/
/*voronoiRegionThreshold() decides foreground/background for bricks         */
/*depending on the local brick statistics.                                  */
/*this is the optimizier region growing version                         */
/****************************************************************************/
void Segmentation::optimizedVoronoiRegionThreshold(unsigned int noOfRegions, std::vector< Contour * > barrelContours, int z)
{
	#ifdef DEBUG
	DebugLog << "Optimized Voronoi region growing: preparations" << std::endl;
	#endif
	
	std::list< unsigned char > * pixelValueLists = new std::list< unsigned char >[noOfRegions];
	std::list< ImageType::IndexType > * pixelIndexLists = new std::list< ImageType::IndexType >[noOfRegions];
	
	int ** histogram = new int * [noOfRegions];
	for(int ii = 0; ii < noOfRegions; ++ii)
	{
		histogram[ii] = new int[256];
	}
	unsigned int * threshold = new unsigned int[noOfRegions];
	float * regionMeans = new float[noOfRegions];
	float * regionSigmas = new float[noOfRegions];
	float * SNR = new float[noOfRegions];
	
	for(int i = 0 ; i < noOfRegions ; ++i)
	{
		regionMeans[i] = 0;
		regionSigmas[i] = 0;
	}
	
	ConstIteratorType pxValIter(preprocImage, maximum_region);
	ConstCalcIteratorType voronoiIter(voronoiMap, maximum_region);
	ConstCalcIteratorType distanceIter(distanceMap, maximum_region);
	for(pxValIter.GoToBegin(), voronoiIter.GoToBegin(), distanceIter.GoToBegin();
	!pxValIter.IsAtEnd() && !voronoiIter.IsAtEnd() && !distanceIter.IsAtEnd();
	++pxValIter, ++voronoiIter, ++distanceIter)
	{
		if(voronoiIter.Get())
			if(pxValIter.Get() != 255)
			{
				int cell = (int)(voronoiIter.Get()+0.5)-1;
				pixelValueLists[cell].push_back(pxValIter.Get());
				pixelIndexLists[cell].push_back(pxValIter.GetIndex());
			}
	}
	
	calculateRegionalHistograms( histogram, pixelValueLists, noOfRegions );
	calculateRegionalHistograms( regionMeans, regionSigmas, pixelValueLists, noOfRegions );
	
	for(int regionID = 0; regionID < noOfRegions; ++regionID)
	{
		if(!barrelContours[regionID]->getValid())
			continue;
		threshold[regionID] = otsuThreshold(histogram[regionID]);
	}
	segmentedImage->FillBuffer(0);
	for(int ii = 0; ii < noOfRegions; ++ii)
	{
		if(!barrelContours[ii]->getValid())
			continue;
		
		if(barrelContours[ii]->getOptimizeFlag())
		{
			int maxLevel = threshold[ii];
			double minLevel = 0;
			for(int jj = 0; jj < 256; ++jj)
				if(histogram[ii][jj])
				{
					minLevel = jj;
					break;
				}
			int startLevel = lround(minLevel);
			if(startLevel >= maxLevel)
				startLevel = maxLevel - 1;
			
			#ifdef DEBUG
			DebugLog << "Optimized Voronoi region growing: region growing" << std::endl;
			#endif
			
			optimizedNeighborhoodRegionGrowing(segmentedImage, pixelIndexLists[ii], startLevel, maxLevel, barrelContours[ii]);
		}
	}
	
	for(int ii = 0; ii < noOfRegions; ++ii)
	{
		pixelValueLists[ii].clear();
		pixelIndexLists[ii].clear();
	}
	delete [] pixelValueLists, delete [] pixelIndexLists;
	for(int ii = 0; ii < noOfRegions; ++ii)
	{
		delete [] histogram[ii];
	}
	delete [] histogram, delete [] threshold, delete [] regionMeans, delete [] regionSigmas, delete [] SNR;
};

// starting from seeds @ startLevel, the intensity level is slowly raised
// and pixels neighboring already accepted pixels are also marked as foreground.
// algorithm stops at maxLevel.
void Segmentation::maxLevelNeighborhoodRegionGrowing(ImageType::Pointer segmentedImage, std::list< ImageType::IndexType > pxIndexList, int startLevel, int maxLevel, bool regular)
{
	ImageType::SizeType radius;
	radius[0] = 1, radius[1] = 1, radius[2] = 0;
	ConstNeighborhoodIteratorType pxValIter;
	
	if(regular)
		pxValIter = ConstNeighborhoodIteratorType(radius, sigmoidImage, maximum_region);
	else
		pxValIter = ConstNeighborhoodIteratorType(radius, preprocImage, maximum_region);
	
	SegNeighborhoodIteratorType segIter(radius, segmentedImage, maximum_region);
	std::list< ImageType::IndexType >::iterator pxIndexIter;
	for(pxIndexIter = pxIndexList.begin(); pxIndexIter != pxIndexList.end(); ++pxIndexIter)
	{
		pxValIter.SetLocation(*pxIndexIter);
		segIter.SetLocation(*pxIndexIter);
		if(pxValIter.GetCenterPixel() <= startLevel)
			segIter.SetCenterPixel(255);
	}
	for(int level = startLevel + 1; level <= maxLevel; ++level)
	{
		unsigned int scoreHistogram[9];
		for(int jj = 0; jj < 9; ++jj)
			scoreHistogram[jj] = 0;
		std::list< ImageType::IndexType > * scoreIndexList = new std::list< ImageType::IndexType >[9];
		for(pxIndexIter = pxIndexList.begin(); pxIndexIter != pxIndexList.end(); ++pxIndexIter)
		{
			pxValIter.SetLocation(*pxIndexIter);
			segIter.SetLocation(*pxIndexIter);
			if(pxValIter.GetCenterPixel() == level)
			{
				int score = 0;
				for(int jj = 0; jj < 8; ++jj)
				{
					if(segIter.GetPixel(neighborhood8[jj]))
						++score;
				}
				++scoreHistogram[score];
				scoreIndexList[score].push_back(pxValIter.GetIndex());
			}
		}
		unsigned int foreground = 0;
		for(int scoreLevel = 8; scoreLevel > 0; --scoreLevel)
		{
			if(scoreHistogram[scoreLevel])
			{
				std::list< ImageType::IndexType >::iterator levelSegIter;
				for(levelSegIter = scoreIndexList[scoreLevel].begin(); levelSegIter != scoreIndexList[scoreLevel].end(); ++levelSegIter)
				{
					segIter.SetLocation(*levelSegIter);
					segIter.SetCenterPixel(255);
				}
				foreground += scoreHistogram[scoreLevel];
				segmentedImage->Update();
			}
		}
		for(int jj = 0; jj < 9; ++jj)
			scoreIndexList[jj].clear();
		delete [] scoreIndexList;
	}
};

//
void Segmentation::optimizedNeighborhoodRegionGrowing(ImageType::Pointer segmentedImage, std::list< ImageType::IndexType > pxIndexList, int startLevel, int maxLevel, Contour * contour)
{
	#ifdef DEBUG
	DebugLog << "\toptimized neighborhood region growing for barrel ID = " << contour->getBarrelID() << std::endl;
	DebugLog << "\t\tbarrel valid = " << contour->getValid() << std::endl;
	DebugLog << "\t\tbarrel optimize flag = " << contour->getOptimizeFlag() << std::endl;
	DebugLog << std::endl;
	#endif
	ConstIndexIteratorType pxValIter(preprocImage, maximum_region);
	std::list< unsigned char > fgPxList;
	std::list< unsigned char > bgPxList;
	std::list< ImageType::IndexType > fgIndexList;
	std::list< ImageType::IndexType > bgIndexList;
	std::list< ImageType::IndexType >::iterator pxIndexIter;
	for(pxIndexIter = pxIndexList.begin(); pxIndexIter != pxIndexList.end(); ++pxIndexIter)
	{
		pxValIter.SetIndex(*pxIndexIter);
		if(pxValIter.Get() <= startLevel)
		{
			fgPxList.push_back(pxValIter.Get());
			fgIndexList.push_back(*pxIndexIter);
		}
		else
		{
			bgPxList.push_back(pxValIter.Get());
			bgIndexList.push_back(*pxIndexIter);
		}
	}
	
	float * levelCostFunction = new float[maxLevel-startLevel+1];
	float * levelSNR = new float[maxLevel-startLevel+1];
	float * levelRadius = new float[maxLevel-startLevel+1];
	float * levelDensity = new float[maxLevel-startLevel+1];
	float * fgMean, *fgSigma, *bgMean, *bgSigma;
	fgMean = new float, fgSigma = new float, bgMean = new float, bgSigma = new float;
	calculateRegionalHistograms(fgMean, fgSigma, &fgPxList, 1);
	calculateRegionalHistograms(bgMean, bgSigma, &bgPxList, 1);
	
	levelSNR[0] = (*bgMean)/(*fgMean);
	levelRadius[0] = averageSphereAroundPxIDs(fgIndexList);
	levelDensity[0] = (float)fgPxList.size()/(PI*levelRadius[0]*levelRadius[0]);	//Pi doesn't really matter...
	for(unsigned int level = startLevel + 1; level <= maxLevel; ++level)
	{
		std::list< unsigned char >::iterator bgLevelIt;
		for(pxIndexIter = pxIndexList.begin(); pxIndexIter != pxIndexList.end(); ++pxIndexIter)
		{
			pxValIter.SetIndex(*pxIndexIter);
			if(pxValIter.Get() == level)
			{
				fgPxList.push_back(pxValIter.Get());
				fgIndexList.push_back(*pxIndexIter);
			}
		}
		levelRadius[level-startLevel] = averageSphereAroundPxIDs(fgIndexList);
		levelDensity[level-startLevel] = (float)fgPxList.size()/(PI*levelRadius[level-startLevel]*levelRadius[level-startLevel]);
		levelSNR[level-startLevel] = avgSphereSNR(levelRadius[level-startLevel], &fgIndexList, &pxIndexList);
	}
	
	#ifdef PIPELINE_DOC
	std::string costFctFilename(outputFilename);
	char * barrelIDChar = new char[32];
	sprintf(barrelIDChar, "%02d", contour->getBarrelID());
	costFctFilename += "_opt_parameters_barrel_";
	costFctFilename += barrelIDChar;
	costFctFilename += ".csv";
	std::ofstream costFctFile(costFctFilename.c_str());
	costFctFile << "ID (level - start level)\tSNR\tRadius\tDensity\tcost fct" << std::endl;
	#endif
	
	//calculate DeltaSNR/DeltaR to determine best radius
	//best radius @ minimum (DeltaSNR/DeltaR negative)
	for(int ii = 0; ii < maxLevel-startLevel+1; ++ii)
	{
		if(levelDensity[ii] == 0)
			levelCostFunction[ii] = 0;
		else
			levelCostFunction[ii] = levelSNR[ii]/levelDensity[ii];
		
		#ifdef PIPELINE_DOC
		costFctFile << ii << "\t" << levelSNR[ii] << "\t" << levelRadius[ii] << "\t" << levelDensity[ii] << "\t" << levelCostFunction[ii] << std::endl;
		#endif
	}
	
	float max = 0, deltaThresh = 2;
	int optimalLevel = startLevel;
	int optimalID = 0;
	bool foundLocalMax = 0;
	while(!max && deltaThresh <= 3)
	{
		for(int ii = 1; ii < maxLevel-startLevel; ++ii)
			if(levelCostFunction[ii] > levelCostFunction[ii-1] && levelCostFunction[ii] > levelCostFunction[ii+1])
			{
				float deltaForward = levelCostFunction[ii] - levelCostFunction[ii+1];
				float deltaBack = levelCostFunction[ii] - levelCostFunction[ii-1];
				if(deltaForward < deltaThresh && deltaBack < deltaThresh)
				{
					if(levelCostFunction[ii] > max)
					{
						optimalID = ii;
						optimalLevel = startLevel + ii;
						max = levelCostFunction[ii];
						foundLocalMax = 1;
					}
				}
			}
		++deltaThresh;
	}
	if(!foundLocalMax)
	{
		std::cout << "Warning! Could not find local max for optimal radius/SNR value!" << std::endl;
		std::cout << "Setting radius to average radius..." << std::endl;
		#ifdef DEBUG
		DebugLog << "\t\tWarning! Could not find local max for optimal radius/SNR value!" << std::endl;
		DebugLog << "\t\tSetting radius to average radius..." << std::endl;
		DebugLog << std::endl;
		#endif
		float avgRadius = 0;
		for(int ii = 0; ii < maxLevel-startLevel+1; ++ii)
			avgRadius += levelRadius[ii];
		avgRadius = avgRadius/(float)(maxLevel-startLevel+1);
// 		std::cout << "avg radius = " << avgRadius << std::endl;
		float minDist = abs(levelRadius[0] - avgRadius);
		for(int ii = 0; ii < maxLevel-startLevel+1; ++ii)
			if(abs(levelRadius[ii] - avgRadius) < minDist)
			{
				optimalID = ii;
				optimalLevel = startLevel + ii;
			}
	}
	
	float optimalRadius = levelRadius[optimalID];
	if(optimalRadius < 20)
		optimalRadius = 20;
	
	#ifdef PIPELINE_DOC
	costFctFile << "found local max: " << foundLocalMax << std::endl;
	costFctFile << "optimal ID = " << optimalID << std::endl;
	costFctFile << "optimal level = " << optimalLevel << std::endl;
	costFctFile << "optimal radius = " << optimalRadius << std::endl;
	costFctFile.close();
	#endif
	
	ImageType::IndexType segmentationCenter;
	segmentationCenter.Fill(0);
	unsigned long segCircleCount = 0;
	std::list< ImageType::IndexType >::const_iterator fgPxIter;
	for(fgPxIter = fgIndexList.begin(); fgPxIter != fgIndexList.end(); ++fgPxIter)
	{
		ImageType::IndexType tmpIndex = *fgPxIter;
		pxValIter.SetIndex(tmpIndex);
		if(pxValIter.Get() <= optimalLevel)
		{
			segmentationCenter[0] += tmpIndex[0];
			segmentationCenter[1] += tmpIndex[1];
			++segCircleCount;
		}
	}
	segmentationCenter[0] = (float)segmentationCenter[0]/(float)segCircleCount + 0.5;
	segmentationCenter[1] = (float)segmentationCenter[1]/(float)segCircleCount + 0.5;
	
	fgPxList.clear();
	pxIndexIter = pxIndexList.begin();
	int optCircleMinLevel = optimalLevel;
	while(pxIndexIter != pxIndexList.end())
	{
		ImageType::IndexType tmpIndex = *pxIndexIter;
		float dist = sqrt((tmpIndex[0] - segmentationCenter[0])*(tmpIndex[0] - segmentationCenter[0]) + (tmpIndex[1] - segmentationCenter[1])*(tmpIndex[1] - segmentationCenter[1]));
		if(dist <= optimalRadius)
		{
			pxValIter.SetIndex(*pxIndexIter);
			PixelType tmpPxVal = pxValIter.Get();
			fgPxList.push_back(tmpPxVal);
			if(tmpPxVal < optCircleMinLevel)
				optCircleMinLevel = tmpPxVal;
		}
		++pxIndexIter;
	}
	
	int * segCircleHisto = new int[256];
	calculateRegionalHistograms(&segCircleHisto, &fgPxList, 1);
	unsigned int segCircleThreshold = otsuThreshold(segCircleHisto);
	contour->addAttribute(segCircleThreshold);
	
	startLevel = optCircleMinLevel;
	startLevel = (startLevel + segCircleThreshold)/2;
	if(startLevel >= segCircleThreshold)
		startLevel = segCircleThreshold - 1;
	bool regular = !(contour->getOptimizeFlag());
	#ifdef DEBUG
	DebugLog << "\t\tstarting actual region growing algorithm..." << std::endl;
	DebugLog << std::endl;
	#endif
	maxLevelNeighborhoodRegionGrowing(segmentedImage, pxIndexList, startLevel, segCircleThreshold, regular);
	
	delete [] levelCostFunction, delete [] levelSNR, delete [] levelRadius, delete [] levelDensity, delete [] segCircleHisto;
	delete fgMean, delete fgSigma, delete bgMean, delete bgSigma;
};

//calculate radius of sphere as average of bounding box dimensions
//around pixels supplied in pxIDs
float Segmentation::averageSphereAroundPxIDs(std::list< ImageType::IndexType > pxIDs)
{
	if(!pxIDs.size())
		return 0;
	
	ImageType::IndexType currIndex;
	std::list< ImageType::IndexType >::iterator pxIDIt;
	int minX, maxX, minY, maxY;
	pxIDIt = pxIDs.begin();
	currIndex = *pxIDIt;
	minX = maxX = currIndex[0];
	minY = maxY = currIndex[1];
	++pxIDIt;
	for( ; pxIDIt != pxIDs.end(); ++pxIDIt)
	{
		currIndex = *pxIDIt;
		if(currIndex[0] < minX)
			minX = currIndex[0];
		if(currIndex[0] > maxX)
			maxX = currIndex[0];
		if(currIndex[1] < minY)
			minY = currIndex[1];
		if(currIndex[1] > maxY)
			maxY = currIndex[1];
	}
	
	return (float)(abs(maxX-minX) + abs(maxY-minY))*0.25;	// 1/2 for average and 1/2 for radius
};

float Segmentation::avgSphereSNR(float radius, std::list< ImageType::IndexType > * pxIDs, std::list< ImageType::IndexType > * pxIndexList)
{
	if(!pxIDs->size())
		return 0;
	
	ImageType::IndexType segmentationCenter;
	segmentationCenter.Fill(0);
	unsigned long segCircleCount = 0;
	std::list< ImageType::IndexType >::const_iterator fgPxIter;
	for(fgPxIter = pxIDs->begin(); fgPxIter != pxIDs->end(); ++fgPxIter)
	{
		ImageType::IndexType tmpIndex = *fgPxIter;
		segmentationCenter[0] += tmpIndex[0];
		segmentationCenter[1] += tmpIndex[1];
		++segCircleCount;
	}
	segmentationCenter[0] = (float)segmentationCenter[0]/(float)segCircleCount + 0.5;
	segmentationCenter[1] = (float)segmentationCenter[1]/(float)segCircleCount + 0.5;
	
	std::list< unsigned char > * insideCircle = new std::list< unsigned char >;
	std::list< unsigned char > * outsideCircle = new std::list< unsigned char >;
	ConstIndexIteratorType pxValIter(preprocImage, maximum_region);
	std::list< ImageType::IndexType >::const_iterator fgIndexIter;
	for(fgIndexIter = pxIndexList->begin(); fgIndexIter != pxIndexList->end(); ++fgIndexIter)
	{
		ImageType::IndexType tmpIndex = *fgIndexIter;
		pxValIter.SetIndex(tmpIndex);
		float dist = sqrt((tmpIndex[0] - segmentationCenter[0])*(tmpIndex[0] - segmentationCenter[0]) + (tmpIndex[1] - segmentationCenter[1])*(tmpIndex[1] - segmentationCenter[1]));
		if(dist <= radius)
			insideCircle->push_back(pxValIter.Get());
		else
			outsideCircle->push_back(pxValIter.Get());
	}
	float * fgMean, * fgSigma, * bgMean, * bgSigma;
	fgMean = new float, fgSigma = new float, bgMean = new float, bgSigma = new float;
	if(insideCircle->size() && outsideCircle->size())
	{
		calculateRegionalHistograms(fgMean, fgSigma, insideCircle, 1);
		calculateRegionalHistograms(bgMean, bgSigma, outsideCircle, 1);
		if(*fgMean)
		{
			float in = *fgMean, out = *bgMean;
			insideCircle->clear();
			outsideCircle->clear();
			delete insideCircle, delete outsideCircle;
			delete fgMean, delete fgSigma, delete bgMean, delete bgSigma;
			return out/in;
		}
	}
	insideCircle->clear();
	outsideCircle->clear();
	delete insideCircle, delete outsideCircle;
	delete fgMean, delete fgSigma, delete bgMean, delete bgSigma;
	return 0;
};

/****************************************************************************/
/*Implementation of Otsu threshold for given histogram                      */
/****************************************************************************/
unsigned int Segmentation::otsuThreshold(int * histogram)
{
	unsigned int currentThreshold = 0;
	unsigned long totalNoOfPixels = 0;
	double weightk = 0, sigmab = 0, mean = 0, meank1 = 0, meank2 = 0;
	
	for(int ii = 0; ii < 256; ++ii)
		totalNoOfPixels += histogram[ii];
// 	std::cout << "totalNoOfPixels = " << totalNoOfPixels << std::endl;
	for(int ii = 0; ii < 256; ++ii)
	{
		mean += (double)ii*(double)histogram[ii];
	}
	mean /= (double)totalNoOfPixels;
	
	for(int k = 1; k < 255; ++k)	//sigmab(k==0) == sigmab(k==255) == 0
	{
		weightk = 0, meank1 = 0, meank2 = 0;	// sigmab = weightk*(1-weightk)*(mean1-mean2)^2
		for(int ii = 0; ii <= k; ++ii)
		{
			weightk += (double)histogram[ii];
			meank1 += (double)ii*(double)histogram[ii];
		}
		for(int ii = k + 1; ii < 255; ++ii)
		{
			meank2 += (double)ii*(double)histogram[ii];
		}
		weightk /= (double)totalNoOfPixels;
		meank1 /= (double)totalNoOfPixels;
		meank2 /= (double)totalNoOfPixels;
		if(weightk)
			meank1 /= weightk;
		else
			meank1 = 0;
		if((1 - weightk))
			meank2 /= (1 - weightk);
		else
			meank2 = 0;
		double tmp = weightk*(1 - weightk)*(meank1 - meank2)*(meank1 - meank2);
// 		std::cout << k << "\t" << tmp << std::endl;
		if(tmp > sigmab)
		{
			sigmab = tmp;
			currentThreshold = k;
		}
	}
	
	return currentThreshold;
};

/****************************************************************************/
/*read Amira landmark file with manual Barrel centroid markers              */
/****************************************************************************/
PointSetType::Pointer Segmentation::readBarrelMarkerFile(const char * markerFilename)
{
	std::ifstream inputStream(markerFilename);
	PointSetType::Pointer points = PointSetType::New();
	
	if(!inputStream.fail())
	{
		std::string currentLine;
		
		unsigned int currentIndex = 0;
		const unsigned int point = 1, cell = 2;
		unsigned int pointID = 0, cellID = 0;
		
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				std::string::size_type loc1, loc2, loc3;
				
				if(!currentIndex)
				{
					loc1 = currentLine.find("@1", 0);
					if(loc1 == 0)
					{
						
						currentIndex = point;
						continue;
// 						char * tmp = new char[currentLine.size() - 9];
// 						currentLine.copy(tmp, currentLine.size() - 9, 9);
// 						int noOfPoints = atoi(tmp);
// 						points->SetDataTypeToFloat();
// 						points->SetNumberOfPoints(noOfPoints);
// 						delete [] tmp;
					}
				}
				
				else if(currentIndex == point)
				{
					loc1 = currentLine.find_first_of("0123456789", 0);
					loc2 = currentLine.find_first_of("+-", 0);
					if(loc2 != std::string::npos)
						if(loc2 < loc1)
							loc1 = loc2;
					loc2 = currentLine.find_first_of(" \t", loc1 + 1);
					loc3 = currentLine.find_first_of(" \t", loc2 + 1);
					char * tmp1 = new char[loc2 - loc1];
					char * tmp2 = new char[loc3 - loc2 - 1];
// 					char * tmp3 = new char[currentLine.size() - loc3 - 1];
					currentLine.copy(tmp1, loc2 - loc1, loc1);
					currentLine.copy(tmp2, loc3 - loc2 - 1, loc2 + 1);
// 					currentLine.copy(tmp3, currentLine.size() - loc3 - 1, loc3 + 1);
					PointType tmpPoint;
					tmpPoint[0] = atof(tmp1), tmpPoint[1] = atof(tmp2);
					points->SetPoint(pointID, tmpPoint);
// 					float pointCoords[] = {atof(tmp1), atof(tmp2), atof(tmp3)};
// 					points->SetPoint(pointID, pointCoords);
					++pointID;
					delete [] tmp1;
					delete [] tmp2;
// 					delete [] tmp3;
				}
			}
		}
// 		surface->SetPoints(points);
// 		surface->Update();
	}
	inputStream.close();
	
// 	points->Print(std::cout);
// 	PointType printPoint;
// 	for(int ii = 0; ii < points->GetNumberOfPoints(); ++ii)
// 	{
// 		points->GetPoint(ii, &printPoint);
// 		std::cout << "Point " << ii << " = [" << printPoint[0] << "," << printPoint[1] << "]" << std::endl;
// 	}
	return points;
};

/****************************************************************************/
/*computes the 2D Voronoi diagram from the barrel markers                   */
/****************************************************************************/
void Segmentation::computeVoronoiDiagram(PointSetType::Pointer markers)
{
	ImageType::Pointer markerImage = ImageType::New();
	ImageType::RegionType markerRegion;
	ImageType::SizeType markerSize;
	ImageType::IndexType markerIndex;
	markerSize = maximum_region.GetSize();
	markerSize[2] = 1;
	markerIndex = maximum_region.GetIndex();
	markerRegion.SetIndex(markerIndex);
	markerRegion.SetSize(markerSize);
	markerImage->SetRegions(markerRegion);
	markerImage->Allocate();
	markerImage->FillBuffer(0);
	IndexIteratorType markerIter(markerImage, markerImage->GetLargestPossibleRegion());
	for(int ii = 0; ii < markers->GetNumberOfPoints(); ++ii)
	{
		PointType seedPoint;
		ImageType::IndexType seedIndex;
		markers->GetPoint(ii, &seedPoint);
		seedIndex[0] = (int)(seedPoint[0] + 0.5), seedIndex[1] = (int)(seedPoint[1] + 0.5), seedIndex[2] = maximum_region.GetIndex(2);
		markerIter.SetIndex(seedIndex);
		markerIter.Set(255);
	}
	markerImage->Update();
	
	DistanceMapImageFilterType::Pointer distanceFilter = DistanceMapImageFilterType::New();
	distanceFilter->InputIsBinaryOn();
	distanceFilter->SetInput(markerImage);
	distanceFilter->Update();
	
	voronoiMap = CalcImageType::New();
	distanceMap = CalcImageType::New();
	voronoiMap = distanceFilter->GetVoronoiMap();
	voronoiMap->Update();
	distanceMap = distanceFilter->GetDistanceMap();
	distanceMap->Update();
	
// 	CalcToImageRescaleFilterType::Pointer cToIFilter = CalcToImageRescaleFilterType::New();
// 	cToIFilter->SetOutputMinimum(1);
// 	cToIFilter->SetOutputMaximum(markers->GetNumberOfPoints());
// 	cToIFilter->SetInput(voronoiMap);
// 	cToIFilter->Update();
// 	writeBinaryImage(cToIFilter->GetOutput(), "_voronoi_map", 0);
};

//for each object, determine Voronoi cell it belongs to and whether it has to be split.
//mark (split) objects with label corresponding to the Voronoi cell they belong to
//Can potentially lead to problems in images where #markers>=256 !!!
void Segmentation::splitWrongBarrelObjects(CalcImageType::Pointer voronoiDiagram, unsigned int noOfCells, bool optimizedThreshold)
{
	ImageType::Pointer singleObjectImage = ImageType::New();
	singleObjectImage->SetRegions(maximum_region);
	singleObjectImage->Allocate();
	if(optimizedThreshold)
		regionGrowingObjectLabeling(inputImage, 1, 0);
	else
		regionGrowingObjectLabeling(inputImage, 5, 0);
	unsigned int * voronoiRegionSize = new unsigned int[noOfCells];
	for(int ii = 0; ii < noOfCells; ++ii)
		voronoiRegionSize[ii] = 0;
	ConstCalcIteratorType voronoiSizeIter(voronoiDiagram, voronoiDiagram->GetLargestPossibleRegion());
	for(voronoiSizeIter.GoToBegin(); !voronoiSizeIter.IsAtEnd(); ++voronoiSizeIter)
		if(voronoiSizeIter.Get())
			++voronoiRegionSize[(int)(voronoiSizeIter.Get() + 0.5) - 1];
	
// 	unsigned int * objectSplitFlag = new unsigned int[nrOfObjects];
// 	for(int ii = 0; ii < nrOfObjects; ++ii)
// 		objectSplitFlag[ii] = 0;
	
	#ifdef PIPELINE_DOC
	if(!optimizedThreshold)
	{
// 		ImageType::Pointer tmpWriteImage = ImageType::New();
// 		tmpWriteImage->SetRegions(inputImage->GetLargestPossibleRegion());
// 		tmpWriteImage->Allocate();
// 		copy(tmpWriteImage, inputImage);
// 		writeBinaryImage(tmpWriteImage, "_barrels_step11_a_split_segmented_objects", 0);
		writeInputImage("_barrels_step11_a_split_segmented_objects.tif");
	}
	else
		writeInputImage("_barrels_step11_a_opt_split_segmented_objects.tif");
	#endif
	
	for(int label = 1; label <=nrOfObjects; ++label)
	{
		singleObjectImage->FillBuffer(0);
		ConstCalcIteratorType voronoiIter(voronoiDiagram, voronoiDiagram->GetLargestPossibleRegion());
		IteratorType2 selectIter(singleObjectImage, singleObjectImage->GetLargestPossibleRegion());
		ConstObjectIteratorType objectIter(labeledImage, labeledImage->GetLargestPossibleRegion());
		for(voronoiIter.GoToBegin(), selectIter.GoToBegin(), objectIter.GoToBegin();
		!voronoiIter.IsAtEnd() && !selectIter.IsAtEnd() && !objectIter.IsAtEnd();
		++voronoiIter, ++selectIter, ++objectIter)
		{
			if(objectIter.Get() == label && voronoiIter.Get())
				selectIter.Set((int)(voronoiIter.Get() + 0.5));
		}
		singleObjectImage->Update();
		
		const unsigned int * objectHist = objectHistogram(singleObjectImage);
		unsigned int memberOfCells = 0;
		for(int ii = 0; ii < 256; ++ii)
			if(objectHist[ii])
				++memberOfCells;
// 		std::cout << "Object " << label << " is part of " << memberOfCells << " Voronoi cells" << std::endl;
		if(!memberOfCells)
			continue;	//Error! should be member of at least one cell
		if(memberOfCells == 1)
		{
			IteratorType2 splittingIter(inputImage, maximum_region);
			for(splittingIter.GoToBegin(), selectIter.GoToBegin(); !splittingIter.IsAtEnd() && !selectIter.IsAtEnd(); ++splittingIter, ++selectIter)
			{
				if(selectIter.Get())
					splittingIter.Set(selectIter.Get());
			}
			inputImage->Update();
			continue;
		}
		else
		{
			unsigned int * cellIDs = new unsigned int[memberOfCells];
			unsigned int currentID = 0;
			unsigned long totalObjectSize = 0;
			for(int ii = 0; ii < 256; ++ii)
				if(objectHist[ii])
				{
					cellIDs[currentID] = ii;
					++currentID;
					totalObjectSize += objectHist[ii];
				}
			for(int ii = 0; ii < memberOfCells; ++ii)
			{
// 				std::cout << "Checking cell " << cellIDs[ii] << ":" << std::endl;
				ImageType::Pointer tmpImage = ImageType::New();
				tmpImage->SetRegions(maximum_region);
				tmpImage->Allocate();
				tmpImage->FillBuffer(0);
				IteratorType2 tmpIter(tmpImage, maximum_region);
				for(objectIter.GoToBegin(), selectIter.GoToBegin(), tmpIter.GoToBegin();
				!objectIter.IsAtEnd() && !selectIter.IsAtEnd() && !tmpIter.IsAtEnd();
				++objectIter, ++selectIter, ++tmpIter)
				{
					if(objectIter.Get() == label && selectIter.Get() == cellIDs[ii])
							tmpIter.Set(255);
				}
				tmpImage->Update();
				
				ConnectedFilterType::Pointer labeller = ConnectedFilterType::New();
				RelabelType::Pointer relabeller = RelabelType::New();
				ObjectImageType::Pointer labeledImage2 = ObjectImageType::New();
				labeller->SetInput(tmpImage);
				labeller->Update();
				relabeller->SetInput(labeller->GetOutput());
				relabeller->SetMinimumObjectSize( 1 );
				labeledImage2 = relabeller->GetOutput();
				labeledImage2->Update();
				unsigned int splitNumber = relabeller->GetNumberOfObjects();
				const std::vector<unsigned long> tmpObjectSize = relabeller->GetSizeOfObjectsInPixels();
// 				std::cout << splitNumber << " parts in this cell with total size " << objectHist[cellIDs[ii]] << std::endl;
				bool * objectSplitFlag = new bool[splitNumber];
// 				std::cout << "splitting flags:" << std::endl;
				for(int jj = 0; jj < splitNumber; ++jj)
				{
					if((tmpObjectSize[jj] > 0.1*sizeOfObjects[label-1] || tmpObjectSize[jj] > 0.2*voronoiRegionSize[cellIDs[ii]-1]) && sizeOfObjects[label-1] > 1E04)
						objectSplitFlag[jj] = 1;
					else
						objectSplitFlag[jj] = 0;
				}
				
				ImageType::SizeType radius;
				radius[0] = 1, radius[1] = 1, radius[2] = 0;
				SegNeighborhoodIteratorType neighborIter(radius, singleObjectImage, singleObjectImage->GetLargestPossibleRegion());
				ConstObjectIteratorType tmpObjectIter(labeledImage2, labeledImage2->GetLargestPossibleRegion());
				ConstObjectIndexIteratorType tmpObjectIndexIter(labeledImage2, labeledImage2->GetLargestPossibleRegion());
				IteratorType2 splittingIter(inputImage, maximum_region);
				for(neighborIter.GoToBegin(), tmpObjectIter.GoToBegin(), splittingIter.GoToBegin();
				!neighborIter.IsAtEnd() && !tmpObjectIter.IsAtEnd() && !splittingIter.IsAtEnd();
				++neighborIter, ++tmpObjectIter, ++splittingIter)
				{
					if(tmpObjectIter.Get())
					{
						if(objectSplitFlag[tmpObjectIter.Get()-1])
						{
							bool border = 0;
							for(int jj = 0; jj < 4; ++jj)
							{
								unsigned int px = neighborIter.GetPixel(neighborhood4[jj]);
								if(px != cellIDs[ii] && px != 0)
								{
									border = 1;
									break;
								}
							}
							if(border)
								splittingIter.Set(cellIDs[ii]);
							else
								splittingIter.Set(255);
						}
						else
							splittingIter.Set(255);
					}
				}
				inputImage->Update();
				delete [] objectSplitFlag;
			}
			
			delete [] cellIDs;
		}
	}
	#ifdef PIPELINE_DOC
	if(!optimizedThreshold)
	{
		writeInputImage("_barrels_step11_b_split_segmented_objects.tif");
// 		copy(tmpWriteImage, inputImage);
// 		writeBinaryImage(tmpWriteImage, "_barrels_step11_b_split_segmented_objects", 0);
	}
	else
		writeInputImage("_barrels_step11_b_opt_split_segmented_objects.tif");
	#endif
// 	fillBinaryHoles();
// 	watershedObjectSplitting(objectSplitFlag);
// 	delete [] objectSplitFlag;
	
	ImageType::SizeType radius;
	radius[0] = 1, radius[1] = 1, radius[2] = 0;
	SegNeighborhoodIteratorType finalSplittingIter(radius, inputImage, maximum_region);
	for(finalSplittingIter.GoToBegin(); !finalSplittingIter.IsAtEnd(); ++finalSplittingIter)
	{
		unsigned int px = finalSplittingIter.GetCenterPixel();
		if(px)
			if(px != 255)
			{
				bool realBorder = 0;
				for(int ii = 0; ii < 4; ++ii)
				{
					unsigned int otherpx = finalSplittingIter.GetPixel(neighborhood4[ii]);
					if(otherpx && otherpx < 255 && otherpx != px)
						realBorder = 1;
				}
				if(!realBorder)
					finalSplittingIter.SetCenterPixel(255);
			}
	}
	for(finalSplittingIter.GoToBegin(); !finalSplittingIter.IsAtEnd(); ++finalSplittingIter)
		if(finalSplittingIter.GetCenterPixel() < 255)
			finalSplittingIter.SetCenterPixel(0);
	inputImage->Update();
	#ifdef PIPELINE_DOC
	if(!optimizedThreshold)
	{
		writeInputImage("_barrels_step11_c_split_segmented_objects.tif");
// 		copy(tmpWriteImage, inputImage);
// 		writeBinaryImage(tmpWriteImage, "_barrels_step11_c_split_segmented_objects", 0);
	}
	else
		writeInputImage("_barrels_step11_c_opt_split_segmented_objects.tif");
	#endif
	if(!optimizedThreshold)
	{
		binaryErosion(inputImage, 1, 1);	//erode so that holes can be closed w/o closing gaps between barrels
		binaryOpening(inputImage, 1, 3);
		fillBinaryHoles();
	}
	#ifdef PIPELINE_DOC
	if(!optimizedThreshold)
	{
// 		writeBinaryImage(inputImage, "_barrels_step11_b_fill_holes", 0);
		writeInputImage("_barrels_step11_d_fill_holes.tif");
// 		copy(tmpWriteImage, inputImage);
// 		writeBinaryImage(tmpWriteImage, "_barrels_step11_d_fill_holes", 0);
	}
	else
		writeInputImage("_barrels_step11_d_opt_fill_holes.tif");
	#endif
	
// // 	regionGrowingObjectLabeling(inputImage, 100, 0);
// // 	labeledToGreyscaleImage(inputImage, nrOfObjects);
// 	writeInputImage("_split_wrong_objects.tif");
};

// split flagged objects using the ITK watershed filter in (nr of objects) largest objects
// specified by objectSplitFlag
void Segmentation::watershedObjectSplitting(unsigned int * objectSplitFlag)
{
	WatershedFilterType::Pointer watershedFilter = WatershedFilterType::New();
	
	for(int label = 1; label <= nrOfObjects; ++label)
		if(objectSplitFlag[label-1])
		{
			long cellBounds[] = {boundingBox[1], boundingBox[0], boundingBox[3], boundingBox[2]};	//xmin, xmax, ymin, ymax
			long * offset = new long[6];
			ImageType::IndexType tmpIndex;
			ConstObjectIndexIteratorType labelIndexIter(labeledImage, labeledImage->GetLargestPossibleRegion());
			bool emptyCell = 1;
			for(labelIndexIter.GoToBegin(); !labelIndexIter.IsAtEnd(); ++labelIndexIter)
			{
				if(labelIndexIter.Get() == label)
				{
					emptyCell = 0;
					tmpIndex = labelIndexIter.GetIndex();
					if(tmpIndex[0] < cellBounds[0])
						cellBounds[0] = tmpIndex[0];
					if(tmpIndex[0] > cellBounds[1])
						cellBounds[1] = tmpIndex[0];
					if(tmpIndex[1] < cellBounds[2])
						cellBounds[2] = tmpIndex[1];
					if(tmpIndex[1] > cellBounds[3])
						cellBounds[3] = tmpIndex[1];
				}
			}
			if(emptyCell)
			{
				std::cout << "Error! Object with label " << label << " could not be split! Not found in labeledImage." << std::endl;
				continue;
			}
			offset[0] = std::max(cellBounds[0]-1, boundingBox[0]);
			offset[1] = std::min(cellBounds[1]+1, boundingBox[1]);
			offset[2] = std::max(cellBounds[2]-1, boundingBox[2]);
			offset[3] = std::min(cellBounds[3]+1, boundingBox[3]);
			offset[4] = maximum_region.GetIndex()[2];
			offset[5] = maximum_region.GetIndex()[2];
		}
};

//for all objects in one Voronoi cell, determine whether they are actually
//part of the barrel or not based on their distance to the marker
//and merge/delete them if they are/are not part of the same barrel
void Segmentation::mergeTrueBarrelObjects(CalcImageType::Pointer distanceMap, CalcImageType::Pointer voronoiDiagram, unsigned int noOfCells, std::vector< Contour * > barrelContours,int zCoord, bool optimizedThreshold)
{
	#ifdef DEBUG
	DebugLog << "in Segmentation::mergeTrueBarrelObjects()" << std::endl;
	DebugLog << "optimization mode = " << optimizedThreshold << std::endl;
	DebugLog << std::endl;
	#endif
	float threshold1;	//threshold for removal of distant objects	use pointers to adjust them individually to contour (opt/non-opt)
	float threshold2;	//threshold for removal of close objects
	if(optimizedThreshold)
	{
		threshold1 = 200;	//threshold for removal of distant objects
		threshold2 = 25;	//threshold for removal of close objects
		regionGrowingObjectLabeling(inputImage, 10, 0);
	}
	else
	{
		threshold1 = 200;	//threshold for removal of distant objects
		threshold2 = 30;	//threshold for removal of close objects
		regionGrowingObjectLabeling(inputImage, 100, 0);
	}
// 	std::list<std::list<std::vector<float> > > pxcoord_edges;
	
	float * minDistance1 = new float[nrOfObjects];
	int * partOfCell = new int[nrOfObjects];
	for(int ii = 0; ii < nrOfObjects; ++ii)
		minDistance1[ii] = 1E06;
	ConstCalcIteratorType distanceIter(distanceMap, distanceMap->GetLargestPossibleRegion());
	ConstCalcIteratorType voronoiIter(voronoiDiagram, voronoiDiagram->GetLargestPossibleRegion());
	ConstObjectIteratorType objectIter(labeledImage, labeledImage->GetLargestPossibleRegion());
	for(distanceIter.GoToBegin(), voronoiIter.GoToBegin(), objectIter.GoToBegin();
	!distanceIter.IsAtEnd() && !voronoiIter.IsAtEnd() && !objectIter.IsAtEnd();
	++distanceIter, ++voronoiIter, ++objectIter)
	{
		if(objectIter.Get() && voronoiIter.Get())
			if(distanceIter.Get() < minDistance1[objectIter.Get()-1])
			{
				minDistance1[objectIter.Get()-1] = distanceIter.Get();
				partOfCell[objectIter.Get()-1] = (int)(voronoiIter.Get() + 0.5);
			}
	}
	
// 	IteratorType2 labelIter(inputImage, maximum_region);
// 	for(objectIter.GoToBegin(), labelIter.GoToBegin(); !objectIter.IsAtEnd() && !labelIter.IsAtEnd(); ++objectIter, ++labelIter)
// 		if(objectIter.Get())
// 			labelIter.Set(partOfCell[objectIter.Get()-1]);
// 	
// 	inputImage->Update();
// 	writeInputImage("_cell_labeled_objects.tif");
	
	ObjectIteratorType deleteIter(labeledImage, labeledImage->GetLargestPossibleRegion());
	IteratorType2 labelIter(inputImage, maximum_region);
	for(deleteIter.GoToBegin(), labelIter.GoToBegin(); !deleteIter.IsAtEnd() && !labelIter.IsAtEnd(); ++deleteIter, ++labelIter)
		if(deleteIter.Get())
		{
			if(minDistance1[deleteIter.Get()-1] > threshold1)
			{
				deleteIter.Set(0);
				labelIter.Set(0);
			}
			else
				labelIter.Set(partOfCell[deleteIter.Get()-1]);
		}
// 	for(labelIter.GoToBegin(); !labelIter.IsAtEnd(); ++labelIter)
// 		if(labelIter.Get() == 255)
// 			labelIter.Set(0);
	labeledImage->Update();
	inputImage->Update();
	#ifdef PIPELINE_DOC
	if(!optimizedThreshold)
	{
		writeInputImage("_barrels_step12_a_delete_distant_objects.tif");
	// 	writeBinaryImage(inputImage, "_barrels_step12_a_delete_distant_objects", 0);
	}
	else
		writeInputImage("_barrels_step12_a_opt_delete_distant_objects.tif");
	#endif
// 	labeledToGreyscaleImage(inputImage, nrOfObjects);
// 	writeBinaryImage(inputImage, "_delete_distant_objects", zCoord);
// 	writeInputImage("_delete_distant_objects.tif");
	for(int voronoiID = 1; voronoiID <= noOfCells; ++voronoiID)
	{
		#ifdef DEBUG
		DebugLog << "\tvoronoi ID = " << voronoiID << std::endl;
		DebugLog << "\tBarrel ID = " << barrelContours[voronoiID-1]->getBarrelID() << std::endl;
		DebugLog << "\tBarrel valid = " << barrelContours[voronoiID-1]->getValid() << std::endl;
		DebugLog << std::endl;
		#endif
		if(barrelContours[voronoiID-1]->getValid())
		{
			long cellBounds[] = {boundingBox[1], boundingBox[0], boundingBox[3], boundingBox[2]};	//xmin, xmax, ymin, ymax
			long * offset = new long[6];
			ImageType::IndexType tmpIndex;
			ConstIndexIteratorType labelIndexIter(inputImage, maximum_region);
			bool emptyCell = 1;
			for(labelIndexIter.GoToBegin(); !labelIndexIter.IsAtEnd(); ++labelIndexIter)
			{
				if(labelIndexIter.Get() == voronoiID)
				{
					emptyCell = 0;
					tmpIndex = labelIndexIter.GetIndex();
					if(tmpIndex[0] < cellBounds[0])
						cellBounds[0] = tmpIndex[0];
					if(tmpIndex[0] > cellBounds[1])
						cellBounds[1] = tmpIndex[0];
					if(tmpIndex[1] < cellBounds[2])
						cellBounds[2] = tmpIndex[1];
					if(tmpIndex[1] > cellBounds[3])
						cellBounds[3] = tmpIndex[1];
				}
			}
			if(emptyCell)
				continue;
			offset[0] = std::max(cellBounds[0]-1, boundingBox[0]);
			offset[1] = std::min(cellBounds[1]+1, boundingBox[1]);
			offset[2] = std::max(cellBounds[2]-1, boundingBox[2]);
			offset[3] = std::min(cellBounds[3]+1, boundingBox[3]);
			offset[4] = maximum_region.GetIndex()[2];
			offset[5] = maximum_region.GetIndex()[2];
	// 		std::cout << "Bounds of region " << voronoiID << " are [" << cellBounds[0] << "," << cellBounds[1] << "] x [" << cellBounds[2] << "," << cellBounds[3] << "]" << std::endl;
	// 		std::flush(std::cout << "Offsets of region " << voronoiID << " are [" << offset[0] << "," << offset[1] << "] x [" << offset[2] << "," << offset[3] << "] x [" << offset[4] << "," << offset[5] << "]" << std::endl);
			ImageType::Pointer cellImage = ImageType::New();
			ImageType::RegionType regionOfInterest;
			ImageType::RegionType cropRegion;
			ImageType::IndexType roiIndex;
			ImageType::SizeType roiSize;
			roiIndex[0] = offset[0], roiIndex[1] = offset[2], roiIndex[2] = offset[4];
			roiSize[0] = offset[1] - offset[0] + 1, roiSize[1] = offset[3] - offset[2] + 1, roiSize[2] = offset[5] - offset[4] + 1;
			regionOfInterest.SetIndex(roiIndex), regionOfInterest.SetSize(roiSize);
	// 		regionOfInterest.Print(std::cout);
			cropRegion.SetIndex(maximum_region.GetIndex()), cropRegion.SetSize(roiSize);
			cellImage->SetRegions(cropRegion);
			cellImage->Allocate();
			cellImage->FillBuffer(0);
	// 		cellImage.Print(std::cout);
	// 		std::cout << "Writing all objects of region " << voronoiID << " to smaller region..." << std::endl;
			ConstIteratorType cellObjectIter(inputImage, regionOfInterest);
			IteratorType2 cellRegionIter(cellImage, cropRegion);
			for(cellObjectIter.GoToBegin(), cellRegionIter.GoToBegin(); !cellObjectIter.IsAtEnd() && !cellRegionIter.IsAtEnd(); ++cellObjectIter, ++cellRegionIter)
				if(cellObjectIter.Get() == voronoiID)
					cellRegionIter.Set(255);
			cellImage->Update();
			#ifdef PIPELINE_DOC
			if(!optimizedThreshold)
			{
				writeBinaryImage(cellImage, "_barrels_step12_object_image", voronoiID);
			}
			else
				writeBinaryImage(cellImage, "_barrels_step12_opt_object_image", voronoiID);
			#endif
			
			regionGrowingObjectLabeling(cellImage, 1, 0);
			
			if(nrOfObjects)
			{
				if(nrOfObjects == 1)
				{
					#ifdef DEBUG
					DebugLog << "\tstarting contour extraction for barrel ID = " << barrelContours[voronoiID-1]->getBarrelID() << std::endl;
					DebugLog << std::endl;
					#endif
					std::list< std::vector< float > > contour = barrelContourExtraction(cellImage, offset, 40, zCoord, barrelContours[voronoiID-1]);
					barrelContours[voronoiID-1]->setEdgeList(contour);
					#ifdef DEBUG
					DebugLog << "\tcompleted contour extraction for barrel ID = " << barrelContours[voronoiID-1]->getBarrelID() << std::endl;
					DebugLog << std::endl;
					#endif
					continue;
				}
				
	// 			std::flush(std::cout << nrOfObjects << " objects in region " << voronoiID << std::endl);
				float * minDistance2 = new float[nrOfObjects];
				for(int ii = 0; ii < nrOfObjects; ++ii)
					minDistance2[ii] = 1E06;
				
				ObjectIteratorType objectLabelIter(labeledImage, labeledImage->GetLargestPossibleRegion());
				ConstCalcIteratorType objectDistanceIter(distanceMap, regionOfInterest);
				for(objectLabelIter.GoToBegin(), objectDistanceIter.GoToBegin();
				!objectLabelIter.IsAtEnd() && !objectDistanceIter.IsAtEnd();
				++objectLabelIter, ++objectDistanceIter)
				{
					unsigned long tmpLabel = objectLabelIter.Get();
					float tmpDist = objectDistanceIter.Get();
					if(tmpDist < minDistance2[tmpLabel-1])
						minDistance2[tmpLabel-1] = tmpDist;
				}
				
				unsigned long closestObject = 0;
				for(unsigned long ii = 1; ii < nrOfObjects; ++ii)
					if(minDistance2[ii] < minDistance2[closestObject])
						closestObject = ii;
				++closestObject;	// closestObject == label of closest object
	// 			std::flush(std::cout << "closestObject = " << closestObject << std::endl);
	// 			for(int ii = 0; ii < nrOfObjects; ++ii)
	// 				std::flush(std::cout << "minDistance2[" << ii << "] = " << minDistance2[ii] << std::endl);
				
				float ** distances = new float *[nrOfObjects];
				for(int ii = 0; ii < nrOfObjects; ++ii)
				{
					distances[ii] = new float[nrOfObjects];
					for(int jj = 0; jj < nrOfObjects; ++jj)
						distances[ii][jj] = 1E06;
				}
				for(unsigned long thisLabel = 1; thisLabel <= nrOfObjects; ++thisLabel)
				{
	// 				std::flush(std::cout << "Calculating distances to object " << thisLabel << std::endl);
					
					ImageType::Pointer thisObjectImage = ImageType::New();
					CalcImageType::Pointer thisObjectDistImage = CalcImageType::New();
					DistanceMapImageFilterType::Pointer thisObjectDistTransform = DistanceMapImageFilterType::New();
					
					thisObjectImage->SetRegions(cellImage->GetLargestPossibleRegion());
					thisObjectImage->Allocate();
					thisObjectImage->FillBuffer(0);
					
					ConstObjectIteratorType thisLabelIter(labeledImage, labeledImage->GetLargestPossibleRegion());
					IteratorType2 pasteObjectIter(thisObjectImage, thisObjectImage->GetLargestPossibleRegion());
					for(thisLabelIter.GoToBegin(), pasteObjectIter.GoToBegin(); !thisLabelIter.IsAtEnd() && !pasteObjectIter.IsAtEnd(); ++thisLabelIter, ++pasteObjectIter)
						if(thisLabelIter.Get() == thisLabel)
							pasteObjectIter.Set(255);
					
					thisObjectDistTransform->InputIsBinaryOn();
					thisObjectDistTransform->SetInput(thisObjectImage);
					thisObjectDistTransform->Update();
					thisObjectDistImage = thisObjectDistTransform->GetDistanceMap();
					
					ConstCalcIteratorType thisObjectDistIter(thisObjectDistImage, thisObjectDistImage->GetLargestPossibleRegion());
	// 				std::flush(std::cout << "Calculating new minimum distances for all objects..." << std::endl);
					for(thisLabelIter.GoToBegin(), thisObjectDistIter.GoToBegin(); !thisLabelIter.IsAtEnd() && !thisObjectDistIter.IsAtEnd(); ++thisLabelIter, ++thisObjectDistIter)
					{
						unsigned long tmpLabel = thisLabelIter.Get();
						if(tmpLabel)
						{
							if(thisObjectDistIter.Get() < distances[thisLabel-1][tmpLabel-1])
								distances[thisLabel-1][tmpLabel-1] = thisObjectDistIter.Get();
						}
					}
				}
				
				int minSize = 40;
				for(int ii = 0; ii < nrOfObjects; ++ii)
				{
					if(minDistance2[ii] >  minSize)
						minSize = (int)(minDistance2[ii] + 0.5) + 1;
					
					for(int jj = ii; jj < nrOfObjects; ++jj)
						distances[jj][ii] = distances[ii][jj];	//enforce symmetry in case of numerical errors of dist transform
				}
				// decide do/do not delete object:
				// if the graph of all objects where the edges have weights d(i,j) with i,j e {label}
				// offers a connection from object i to closestObject without any weight on that path:
				// d(i,j) > threshold2, then keep object
				bool * objectDeleteFlag = new bool[nrOfObjects];
				GraphType distanceGraph;
				std::list< unsigned long > labelList;
				for(unsigned long thisLabel = 1; thisLabel <= nrOfObjects; ++thisLabel)
				{
					if(thisLabel != closestObject)
					{
						objectDeleteFlag[thisLabel-1] = true;
						labelList.push_back(thisLabel);
					}
					else
						objectDeleteFlag[thisLabel-1] = false;
					
					std::list< unsigned long > otherVertices;
					for(unsigned long nextLabel = 1; nextLabel <= nrOfObjects; ++nextLabel)
					{
						if(nextLabel == thisLabel)
							continue;
						if(distances[thisLabel-1][nextLabel-1] <= threshold2)
							otherVertices.push_back(nextLabel);
					}
					distanceGraph.insert(std::pair< unsigned long, std::list< unsigned long > >(thisLabel, otherVertices));
				}
				// use Dijkstra's algorithm to find path if it exists
				for(unsigned long thisLabel = 1; thisLabel <= nrOfObjects; ++thisLabel)
				{
					if(thisLabel == closestObject)
						continue;
					
					if(distanceGraph[thisLabel].size())
					{
						GraphType::iterator graphIter;
						std::list< unsigned long >::iterator contourIter;
						std::map< unsigned long, float > distance;
						std::list< unsigned long > tmpLabelList;
						tmpLabelList.push_back(thisLabel);
						for(contourIter = distanceGraph[thisLabel].begin(); contourIter != distanceGraph[thisLabel].end(); ++contourIter)
							tmpLabelList.push_back(*contourIter);
						
						for(graphIter = distanceGraph.begin(); graphIter != distanceGraph.end(); ++graphIter)
						{
							distance[graphIter->first] = 1E09;
						}
						distance[tmpLabelList.front()] = 0;
						
						std::list< unsigned long >::iterator vertexIter;
						std::list< unsigned long >::iterator neighborIter;
						
						while(!tmpLabelList.empty())
						{
							unsigned long tmpDist = 0;
							long nextVertex = 0;
							
							vertexIter = tmpLabelList.begin();
							nextVertex = *vertexIter;
							if(nextVertex == closestObject)
							{
								objectDeleteFlag[thisLabel-1] = false;
								break;
							}
							
							tmpLabelList.erase(vertexIter);
							
							for(neighborIter = distanceGraph[nextVertex].begin(); neighborIter != distanceGraph[nextVertex].end(); ++neighborIter)
							{
								float newDist = 0;
								newDist = distance[nextVertex] + 1;
								
								if(newDist < distance[*neighborIter])
								{
									distance[*neighborIter] = newDist;
									tmpLabelList.push_back(*neighborIter);
								}
							}
						}
					}
				}
				ObjectIteratorType deleteLabelIter(labeledImage, labeledImage->GetLargestPossibleRegion());
				for(deleteLabelIter.GoToBegin(); !deleteLabelIter.IsAtEnd(); ++deleteLabelIter)
				{
					unsigned long tmpLabel = deleteLabelIter.Get();
					if(tmpLabel)
						if(objectDeleteFlag[tmpLabel-1])
							deleteLabelIter.Set(0);
				}
				labeledImage->Update();
				labeledToGreyscaleImage(cellImage, nrOfObjects);
				cellImage->Update();
				#ifdef PIPELINE_DOC
				if(!optimizedThreshold)
				{
					writeBinaryImage(cellImage, "_barrels_step13_distant_objects_deleted", voronoiID);
				}
				else
					writeBinaryImage(cellImage, "_barrels_step13_opt_distant_objects_deleted", voronoiID);
				#endif
				#ifdef DEBUG
				DebugLog << "\tstarting contour extraction for barrel ID = " << barrelContours[voronoiID-1]->getBarrelID() << std::endl;
				DebugLog << std::endl;
				#endif
				std::list< std::vector< float > > contour = barrelContourExtraction(cellImage, offset, minSize, zCoord, barrelContours[voronoiID-1]);
				barrelContours[voronoiID-1]->setEdgeList(contour);
				#ifdef DEBUG
				DebugLog << "\tcompleted contour extraction for barrel ID = " << barrelContours[voronoiID-1]->getBarrelID() << std::endl;
				DebugLog << std::endl;
				#endif
			}
		}
	}
};

//calculate barrel contour from (possibly several) barrel objects
//use closing for computational efficiency instead of convex hull
std::list< std::vector<float> > Segmentation::barrelContourExtraction(ImageType::Pointer barrelImage, long * offset, unsigned int minSize, int zCoord, Contour * contour)
{
// 	std::flush(std::cout << "starting barrel contour extraction..." << std::endl);
// 	std::flush(std::cout << "setting up edgeImage and paddedImage..." << std::endl);
	ImageType::Pointer edgeImage = ImageType::New();
	ImageType::Pointer paddedImage = ImageType::New();
	ImageType::RegionType paddedRegion;
	ImageType::RegionType copyRegion;
	ImageType::IndexType paddingOffset;
	ImageType::SizeType paddedSize;
	ImageType::SizeType barrelSize = barrelImage->GetLargestPossibleRegion().GetSize();
	paddedSize[0] = barrelSize[0] + 2*(minSize + 1);
	paddedSize[1] = barrelSize[1] + 2*(minSize + 1);
	paddedSize[2] = barrelSize[2];
	paddingOffset[0] = minSize + 1;
	paddingOffset[1] = minSize + 1;
	paddingOffset[2] = barrelImage->GetLargestPossibleRegion().GetIndex()[2];
	offset[0] -= paddingOffset[0];
	offset[2] -= paddingOffset[1];
	int zOffset = zCoord - (stop - start)/2;
	
	copyRegion.SetIndex(paddingOffset);
	copyRegion.SetSize(barrelSize);
	paddedRegion.SetIndex(barrelImage->GetLargestPossibleRegion().GetIndex());
	paddedRegion.SetSize(paddedSize);
	
// 	copyRegion.Print(std::cout);
// 	paddedRegion.Print(std::cout);
	
	paddedImage->SetRegions(paddedRegion);
	paddedImage->Allocate();
	paddedImage->FillBuffer(0);
	edgeImage->SetRegions(paddedRegion);
	edgeImage->Allocate();
	edgeImage->FillBuffer(0);
// 	std::flush(std::cout << "Allimages/regions allocated!" << std::endl);
// 	barrelImage->Print(std::cout);
// 	paddedImage->Print(std::cout);
// 	edgeImage->Print(std::cout);
	ConstIteratorType copyIter(barrelImage, barrelImage->GetLargestPossibleRegion());
	IteratorType2 pasteIter(paddedImage, copyRegion);
	for(copyIter.GoToBegin(), pasteIter.GoToBegin(); !copyIter.IsAtEnd() && !pasteIter.IsAtEnd(); ++copyIter, ++pasteIter)
		pasteIter.Set(copyIter.Get());
	paddedImage->Update();
	// 	std::flush(std::cout << "Barrel image pasted to padded image, starting edge detection..." << std::endl);
	
	binaryClosing(paddedImage, 1, minSize);
	fillBinaryHoles();
// 	binaryConvexHull(paddedImage, minSize);
	
	#ifdef PIPELINE_DOC
	if(contour->getOptimizeFlag())
	{
		writeBinaryImage(paddedImage, "_final_barrel_object_opt_barrel", contour->getBarrelID());
	}
	else
	{
		writeBinaryImage(paddedImage, "_final_barrel_object_barrel", contour->getBarrelID());
	}
	#endif
	
	//calculate segmented object size
	ConnectedFilterType::Pointer labeller = ConnectedFilterType::New();
	RelabelType::Pointer relabeller = RelabelType::New();
	ObjectFilterType::Pointer objectFilter = ObjectFilterType::New();
	ObjectImageType::Pointer labeledImage2 = ObjectImageType::New();
	labeller->SetInput(paddedImage);
	labeller->SetFullyConnected(0);
	labeller->Update();
	relabeller->SetInput(labeller->GetOutput());
	relabeller->SetMinimumObjectSize( 1 );
	labeledImage2 = relabeller->GetOutput();
	labeledImage2->Update();
	unsigned int objectSize = 0;
	if(relabeller->GetSizeOfObjectsInPixels().size())
		objectSize = relabeller->GetSizeOfObjectsInPixels()[0];
	float radius = sqrt((float)objectSize/PI);
	contour->addAttribute(radius);
	
	//calculate SNR outside/inside segmented object
	std::list< unsigned char > * fgPxValList = new std::list< unsigned char >;
	std::list< unsigned char > * bgPxValList = new std::list< unsigned char >;
	ConstCalcIteratorType voronoiIter(voronoiMap, maximum_region);
	ConstIndexIteratorType inputIter(preprocImage, maximum_region);
	ConstIndexIteratorType barrelObjectIter(paddedImage, paddedRegion);
	unsigned long voronoiRegionSize = 0;
	for(voronoiIter.GoToBegin(), inputIter.GoToBegin(); !voronoiIter.IsAtEnd() && !inputIter.IsAtEnd(); ++voronoiIter, ++inputIter)
	{
		if(voronoiIter.Get() == (*(contour->attributesPointer()))[1] )
		{
			++voronoiRegionSize;
			ImageType::IndexType tmpIndex = voronoiIter.GetIndex();
			tmpIndex[0] -= offset[0], tmpIndex[1] -= offset[2];
// 			paddedRegion.Print(std::cout);
			if(isInBox(paddedRegion, tmpIndex))
			{
				barrelObjectIter.SetIndex(tmpIndex);
				if(barrelObjectIter.Get())
				{
					if(inputIter.Get() != 255)
						fgPxValList->push_back(inputIter.Get());
				}
				else
					if(inputIter.Get() != 255)
						bgPxValList->push_back(inputIter.Get());
			}
		}
	}
// 	std::flush(std::cout << "done with object!" << std::endl);
	float * fgMean = new float, * bgMean = new float;
	float * fgSigma = new float, * bgSigma = new float;
	calculateRegionalHistograms(fgMean, fgSigma, fgPxValList, 1);
	calculateRegionalHistograms(bgMean, bgSigma, bgPxValList, 1);
// 	std::cout << "bgMean = " << *bgMean << "\t---\tbgSize = " << bgPxValList->size() << std::endl;
// 	std::cout << "fgMean = " << *fgMean << "\t---\tfgSize = " << fgPxValList->size() << std::endl;
// 	std::cout << "voronoiRegionSize = " << voronoiRegionSize << std::endl;
// 	writeBinaryImage(paddedImage, "_paddedImage", 0);
	float SNR = 0;
	float sizeRatio = 0;
	if(bgPxValList->size() && fgPxValList->size())
	{
		if(*fgMean )
			SNR = (*bgMean)/(*fgMean);
		sizeRatio = (float)fgPxValList->size()/(float)(fgPxValList->size() + bgPxValList->size());
	}
	contour->addAttribute(SNR);
	contour->addAttribute(sizeRatio);
	delete fgMean, delete fgSigma, delete bgMean, delete bgSigma;
	fgPxValList->clear(), bgPxValList->clear();
	delete fgPxValList, delete bgPxValList;
	
	ShapedNeighborhoodIteratorType neighborIter(radius1, paddedImage, paddedImage->GetLargestPossibleRegion());
	ShapedNeighborhoodIteratorType::Iterator it;
	IteratorType2 edgeIter(edgeImage, edgeImage->GetLargestPossibleRegion());
// 	ImageType::IndexType tmpIndex;
	for(int ii = 0; ii < 4; ++ii)
		neighborIter.ActivateOffset(neighborhood4[ii]);
	
	for(edgeIter.GoToBegin(), neighborIter.GoToBegin(); !edgeIter.IsAtEnd() && !neighborIter.IsAtEnd(); ++edgeIter, ++neighborIter)
	{
		if(neighborIter.GetCenterPixel())
		{
// 			tmpIndex = neighborIter.GetIndex();
			bool neighbor = 0;
			for(it = neighborIter.Begin(); !it.IsAtEnd(); ++it)
			{
				if(!it.Get())
				{
					neighbor = 1;
					break;
				}
			}
			if(!neighbor)
				edgeIter.Set(0);
			else
				edgeIter.Set(255);
		}
		else
			edgeIter.Set(0);
	}
	edgeImage->Update();
	// in the extremely rare, but possible case that there are several contours after the closing,
	// keep only the largest contour object
// 	ConnectedFilterType::Pointer labeller = ConnectedFilterType::New();
// 	RelabelType::Pointer relabeller = RelabelType::New();
// 	ObjectFilterType::Pointer objectFilter = ObjectFilterType::New();
// 	ObjectImageType::Pointer labeledImage2 = ObjectImageType::New();
	labeller->SetInput(edgeImage);
	labeller->SetFullyConnected(1);
	labeller->Update();
	relabeller->SetInput(labeller->GetOutput());
	relabeller->SetMinimumObjectSize( 1 );
	labeledImage2 = relabeller->GetOutput();
	labeledImage2->Update();
// 	unsigned int splitNumber = relabeller->GetNumberOfObjects();
// 	const std::vector<unsigned long> tmpObjectSize = relabeller->GetSizeOfObjectsInPixels();
	objectFilter->SetInput( labeledImage2 );
	objectFilter->SetOutsideValue (0);
	objectFilter->SetInsideValue (255);
	objectFilter->SetLowerThreshold (1);
	objectFilter->SetUpperThreshold (1);
	edgeImage = objectFilter->GetOutput();
	edgeImage->Update();
// 	writeBinaryImage(edgeImage, "_barrel_contour", zCoord);
// 	std::flush(std::cout << "Edge detection finished!" << std::endl);
	std::list< std::vector< float > > barrelContour = contourExtraction(edgeImage);
	std::list< std::vector< float > >::iterator edgePointIt;
	for(edgePointIt = barrelContour.begin(); edgePointIt != barrelContour.end(); ++edgePointIt)
	{
		for(int jj = 0; jj < 3; ++jj)
			(*edgePointIt)[jj] += (float)offset[jj*2];
		
		(*edgePointIt)[2] += (float)zOffset;
	}
	
	return barrelContour;
};

/****************************************************************************/
/*marks background with -1 in calcImage for simple object recognition       */
/****************************************************************************/
ImageType::Pointer Segmentation::MarkBackground(unsigned int highThreshold)
{
	ImageType::Pointer tmpImage = invertImage(inputImage);
	
	regionGrowingObjectLabeling(tmpImage, highThreshold, 1);
	labeledToGreyscaleImage(tmpImage, nrOfObjects);
	
	return tmpImage;
};

ImageType::Pointer Segmentation::invertImage( ImageType::Pointer workingImage )
{
	ImageType::Pointer invertedImage = ImageType::New();
	
	ImageType::RegionType invertedRegion = workingImage->GetLargestPossibleRegion();
	
	invertedImage->SetRegions(invertedRegion);
	invertedImage->Allocate();
	invertedImage->Update();
	
	ConstIteratorType originalIter( workingImage, workingImage->GetLargestPossibleRegion() );
	IteratorType2 invertIter( invertedImage, invertedImage->GetLargestPossibleRegion() );
	
	for( originalIter.GoToBegin(), invertIter.GoToBegin() ; !originalIter.IsAtEnd() && !invertIter.IsAtEnd() ; ++originalIter, ++invertIter )
	{
		invertIter.Set(255 - originalIter.Get());
	}
	
	invertedImage->Update();
	return invertedImage;
};

CalcImageType::Pointer Segmentation::invertImage( CalcImageType::Pointer workingImage )
{
	CalcImageType::Pointer invertedImage = CalcImageType::New();
	
	CalcImageType::RegionType invertedRegion = workingImage->GetLargestPossibleRegion();
	
	invertedImage->SetRegions(invertedRegion);
	invertedImage->Allocate();
	invertedImage->Update();
	
	ConstCalcIteratorType originalIter( workingImage, workingImage->GetLargestPossibleRegion() );
	CalcIteratorType invertIter( invertedImage, invertedImage->GetLargestPossibleRegion() );
	
	for( originalIter.GoToBegin(), invertIter.GoToBegin() ; !originalIter.IsAtEnd() && !invertIter.IsAtEnd() ; ++originalIter, ++invertIter )
	{
		invertIter.Set(255 - originalIter.Get());
	}
	
	invertedImage->Update();
	return invertedImage;
};

/****************************************************************************/
/*medianFilter()                                                            */
/****************************************************************************/
void Segmentation::medianFilter( unsigned int radius )
{
	MedianFilterType::Pointer medianFilter = MedianFilterType::New();
	ImageType::SizeType mediansize;
	
	mediansize[0] = radius;
	mediansize[1] = radius;
	mediansize[2] = 0;
	
	medianFilter->SetRadius( mediansize );
	medianFilter->SetInput( inputImage );
	medianFilter->Update();
	
	inputImage = medianFilter->GetOutput();
	inputImage->Update();
};

/****************************************************************************/
/*copy( image1, image2 ) performs an iterator copy of image2 to image1      */
/****************************************************************************/
void Segmentation::copy( ImageType::Pointer workingImage1, ImageType::Pointer workingImage2 )
{
	ConstIteratorType reader( workingImage2, workingImage2->GetLargestPossibleRegion() );
	IteratorType2 writer( workingImage1, workingImage1->GetLargestPossibleRegion() );
	
	for( reader.GoToBegin(), writer.GoToBegin() ; !reader.IsAtEnd() && !writer.IsAtEnd() ; ++reader, ++writer )
	{
		writer.Set( reader.Get() );
	}
	
	workingImage1->Update();
};

/****************************************************************************/
/*copy( image1, image2 ) performs an iterator copy of image2 to image1      */
/****************************************************************************/
void Segmentation::copy( Image2DType::Pointer workingImage1, ImageType::Pointer workingImage2, ImageType::RegionType planeRegion )
{
	ConstIteratorType reader( workingImage2, planeRegion );
	Iterator2DType writer( workingImage1, workingImage1->GetLargestPossibleRegion() );
	
	for( reader.GoToBegin(), writer.GoToBegin() ; !reader.IsAtEnd() && !writer.IsAtEnd() ; ++reader, ++writer )
	{
		writer.Set( reader.Get() );
	}
	
	workingImage1->Update();
};

/****************************************************************************/
/*copy( image1, image2 ) performs an iterator copy of image2 to image1      */
/****************************************************************************/
void Segmentation::copy( ImageType::Pointer workingImage1, ImageType::Pointer workingImage2, ImageType::RegionType planeRegion )
{
	ConstIteratorType reader( workingImage2, planeRegion );
	IteratorType2 writer( workingImage1, workingImage1->GetLargestPossibleRegion() );
	
	for( reader.GoToBegin(), writer.GoToBegin() ; !reader.IsAtEnd() && !writer.IsAtEnd() ; ++reader, ++writer )
	{
		writer.Set( reader.Get() );
	}
	
	workingImage1->Update();
};

/****************************************************************************/
/*copy( image1, image2 ) performs an iterator copy of image2 to image1      */
/****************************************************************************/
void Segmentation::copy( ImageType::Pointer workingImage1, Image2DType::Pointer workingImage2, ImageType::RegionType planeRegion )
{
	Const2DIteratorType reader( workingImage2, workingImage2->GetLargestPossibleRegion() );
	IteratorType2 writer( workingImage1, planeRegion );
	
	for( reader.GoToBegin(), writer.GoToBegin() ; !reader.IsAtEnd() && !writer.IsAtEnd() ; ++reader, ++writer )
	{
		writer.Set( reader.Get() );
	}
	
	workingImage1->Update();
};

/****************************************************************************/
/*copy( image1, image2 ) performs an iterator copy of image2 to image1      */
/****************************************************************************/
void Segmentation::copy( ImageType::Pointer workingImage1, ImageType::Pointer workingImage2, ImageType::RegionType planeRegion, int second )
{
	ConstIteratorType reader( workingImage2, workingImage2->GetLargestPossibleRegion() );
	IteratorType2 writer( workingImage1, planeRegion );
	
	for( reader.GoToBegin(), writer.GoToBegin() ; !reader.IsAtEnd() && !writer.IsAtEnd() ; ++reader, ++writer )
	{
		writer.Set( reader.Get() );
	}
	
	workingImage1->Update();
};

/****************************************************************************/
/*GlobalThreshold() sets pixels below a grayvalue of 10 to 0                */
/****************************************************************************/

void Segmentation::GlobalThreshold(int thresh)
{
	int end_z = inputImage->GetLargestPossibleRegion().GetSize(2);
	
// 	int lower_threshold = (int)(mean + factor * std);
	int lower_threshold = thresh;
	
// 	std::cout<< "Lower Threshold is Mean (" << mean << ") + " << thresh << " times standard deviation (" << std << "): " << lower_threshold << std::endl;
	
	#pragma omp parallel for	
	for(int i = 0; i < end_z; i++)
	{
// 		ImageType::Pointer parallelImage = ImageType::New();
		ImageType::RegionType parallel_region;
		ImageType::IndexType parallel_index;
		ImageType::SizeType parallel_size;
		
		parallel_index[0] = 0;
		parallel_index[1] = 0;
		parallel_index[2] = i;
		
		parallel_size[0] = inputImage->GetLargestPossibleRegion().GetSize(0);
		parallel_size[1] = inputImage->GetLargestPossibleRegion().GetSize(1);
		parallel_size[2] = 1;
		
		parallel_region.SetIndex(parallel_index);
		parallel_region.SetSize(parallel_size);
// 		parallelImage->SetRegions(parallel_region);
// 		parallelImage->Allocate();
		
// 		ConstIteratorType inputIter(inputImage, parallel_region);
// 		IteratorType2 planeIter(parallelImage, parallel_region);
		IteratorType2 planeIter(inputImage, parallel_region);
		
		for(planeIter.GoToBegin(); !planeIter.IsAtEnd(); ++planeIter)
		{
			if(planeIter.Get() <= lower_threshold)
				planeIter.Set(255);
			
			else planeIter.Set(0);
		}
	}
	
	
	inputImage->Update();
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_pia_step1_threshold", 0);
	#endif
	
	binaryOpening(inputImage, 1, 40); // was: 10
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_pia_step2_opening_R40", 0);
	#endif
	
	regionGrowingObjectLabeling(inputImage, 1, 1);
	// 		relabelImage(labeledImage, 1);
	labeledToGreyscaleImage(inputImage, 1);	//use nrOfObjects for all objects instead
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_pia_step3_region_growing1", 0);
	#endif
	
	binaryClosing(inputImage, 1, 40);
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_pia_step4_closing_R40", 0);
	#endif
// 	writeBinaryImage(inputImage, "_threshold", lower_threshold);
	
// 	ImageType::Pointer bgImage = ImageType::New();
// 	bgImage->SetRegions(inputImage->GetLargestPossibleRegion());
// 	bgImage->Allocate();
// 	copy(bgImage, inputImage);	//for BVP, so ROI is already restricted to the slice!
	
// 	ImageType::Pointer tmpImage = ImageType::New();
	ImageType::Pointer bgImage = MarkBackground(10000);
// 	writeBinaryImage(bgImage, "_simpleBG", 10);
	bgImage = invertImage(bgImage);
	#ifdef PIPELINE_DOC
	writeBinaryImage(bgImage, "_pia_step5_bgImage", 0);
	#endif
	
	ShapedNeighborhoodIteratorType neighborIter(radius1, bgImage, bgImage->GetLargestPossibleRegion());
// 	ShapedNeighborhoodIteratorType neighborIter(radius1, inputImage, inputImage->GetLargestPossibleRegion());
	ShapedNeighborhoodIteratorType::Iterator it;
	IndexIteratorType edgeIter(inputImage, maximum_region);
	ImageType::IndexType tmpIndex;
	for(int ii = 0; ii < 4; ++ii)
		neighborIter.ActivateOffset(neighborhood4[ii]);
	
	for(edgeIter.GoToBegin(), neighborIter.GoToBegin(); !edgeIter.IsAtEnd() && !neighborIter.IsAtEnd(); ++edgeIter, ++neighborIter)
	{
		if(/*edgeIter.Get()*/neighborIter.GetCenterPixel())
		{
			tmpIndex = neighborIter.GetIndex();
			bool neighbor = 0;
			for(it = neighborIter.Begin(); !it.IsAtEnd(); ++it)
			{
				if(!it.Get())
				{
					neighbor = 1;
					break;
				}
			}
			
			if(!neighbor)
				edgeIter.Set(0);
			
			else
				edgeIter.Set(255);
		}
		
		else
			edgeIter.Set(0);
	}
	
	inputImage->Update();
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_pia_step6_edge1", 0);
	#endif
// 	writeInputImage("_edgetest.tif");
	//keep really only largest contour (in case of holes...)
	regionGrowingObjectLabeling(inputImage, 1, 1);
	labeledToGreyscaleImage(inputImage, 1);
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_pia_step7_edge2", 0);
	#endif
// 
	std::list<std::list<std::vector<float> > > pxcoord_edges;
	pxcoord_edges.push_back(contourExtraction(inputImage));
	contourSmoothing(pxcoord_edges.back(), 25);
	contourSampling(pxcoord_edges.back(), 100);
// 	std::cout << "writing amira mesh file..." << std::endl;
	std::list<std::list<std::vector<float> > > amira_edges = GetAmiraContours(pxcoord_edges);
	
	std::list<std::vector<float> > vertex_list = SetContourVertices(amira_edges);
	
	AmiraContourGraphType *tmp_contour = new AmiraContourGraphType(vertex_list, amira_edges);
	amira_contour_graph = tmp_contour;
	
// 	copy(inputImage, bgImage);	//for BVP, so ROI is already restricted to the slice!
};

/****************************************************************************/
/*extract contour of a single object binary image (assumes 2D image for     */
/* 4-connectivity                                                           */
/****************************************************************************/
std::list<std::vector<float > > Segmentation::contourExtraction( ImageType::Pointer workingImage )
{
	ShapedNeighborhoodIteratorType pruningIter(radius1, workingImage, workingImage->GetLargestPossibleRegion());
	ShapedNeighborhoodIteratorType::Iterator iter;
	ImageType::IndexType tmpIndex;
	for(int ii = 0; ii < 8; ++ii)
		pruningIter.ActivateOffset(neighborhood8[ii]);

	bool change = false;
	int count = 1;
	unsigned int neighbor[8];
	const unsigned int ordered_neighborhood[8] = {0,3,5,6,7,4,2,1};	//returns indices for clockwise iteration of px neighborhood
	do
	{
		change = false;
// 		std::cout << "in pruning iteration " << count << std::endl;
		++count;
		for(pruningIter.GoToBegin(); !pruningIter.IsAtEnd(); ++pruningIter)
		{
			if(pruningIter.GetCenterPixel())
			{
				int ii = 0;
				for(iter = pruningIter.Begin(); !iter.IsAtEnd(); ++iter, ++ii )
				{
					neighbor[ii] = iter.Get();
				}
				
				unsigned int lastpx = neighbor[ordered_neighborhood[7]], pxchanges = 0;
				
				for(int jj = 0; jj < 8; ++jj)
				{
					if(neighbor[ordered_neighborhood[jj]] && !lastpx) ++pxchanges;
					lastpx = neighbor[ordered_neighborhood[jj]];
				}
				
				if(pxchanges == 1 )
				{
					pruningIter.SetCenterPixel(0);
					change = true;
				}
			}
		}
	} while(change && count < 20);
	
	workingImage->Update();
// 	std::cout << "Writing pruned contour image..." <<std::endl;
// 	writeInputImage("_contour.tif");
	
	
	GraphType contourGraph = createGraphFromContourImage(workingImage);
	
	std::list< std::vector< int > > contourIndices;
	contourIndices = minimalDistanceGraph(contourGraph, workingImage, 10 );	//10 deprecated; downsampling seperately after smoothing through averaging
	
	std::list< std::vector< float > > contourIndicesFloat;
	std::list< std::vector< int > >::iterator it;
	for(it = contourIndices.begin(); it != contourIndices.end(); ++it)
	{
		std::vector< float > tmp;
		tmp.push_back((*it)[0]);
		tmp.push_back((*it)[1]);
		tmp.push_back((*it)[2]);
		contourIndicesFloat.push_back(tmp);
	}
	
	return contourIndicesFloat;
};

/****************************************************************************/
/*creates undirected Graph from binary image. neighbors are assigned with   */
/*respect to euclidean distance, i.e. first 4-, then 8 neighborhood is      */
/*checked for neighboring pixels                                            */
/****************************************************************************/
GraphType Segmentation::createGraphFromContourImage( ImageType::Pointer workingImage )
{
	GraphType newGraph;
	ShapedNeighborhoodIteratorType contourIter(radius1, workingImage, workingImage->GetLargestPossibleRegion());
	ShapedNeighborhoodIteratorType::Iterator iter;
	const ImageType::SizeType inputSize = workingImage->GetLargestPossibleRegion().GetSize();
	ImageType::IndexType tmpIndex;
	unsigned long tmpVertex = 0;
	unsigned int neighbors[8];
	const unsigned int neighbors4[4] = {1,3,4,6};
	const unsigned int neighbors8[4] = {0,2,5,7};
	for(int ii = 0; ii < 8; ++ii)
		contourIter.ActivateOffset(neighborhood8[ii]);
	
// 	std::cout << "creating graph form contour image..." << std::endl;
	for(contourIter.GoToBegin(); !contourIter.IsAtEnd(); ++contourIter)
	{
		if(contourIter.GetCenterPixel() == 255)
		{
			contourIter.SetCenterPixel(1);
			tmpIndex = contourIter.GetIndex();
			std::list< unsigned long > neighborVertices;
			tmpVertex = tmpIndex[1]*inputSize[0] + tmpIndex[0];
			
			int ii = 0;
			for(iter = contourIter.Begin(); !iter.IsAtEnd(); ++iter, ++ii)
			{
				neighbors[ii] = iter.Get();
			}
			
			if(neighbors[neighbors4[0]])
			{
				neighborVertices.push_back(tmpVertex-inputSize[0]);
				neighbors[neighbors8[0]] = 0;	//do not process 8-connected corners that also have a 4-connection with a 4-connected neighbor of tmpIndex
				neighbors[neighbors8[1]] = 0;
			}
			
			if(neighbors[neighbors4[1]])
			{
				neighborVertices.push_back(tmpVertex-1);
				neighbors[neighbors8[0]] = 0;
				neighbors[neighbors8[2]] = 0;
			}
			
			if(neighbors[neighbors4[2]])
			{
				neighborVertices.push_back(tmpVertex+1);
				neighbors[neighbors8[1]] = 0;
				neighbors[neighbors8[3]] = 0;
			}
			
			if(neighbors[neighbors4[3]])
			{
				neighborVertices.push_back(tmpVertex+inputSize[0]);
				neighbors[neighbors8[2]] = 0;
				neighbors[neighbors8[3]] = 0;
			}
			
			if(neighbors[neighbors8[0]])
			{
				neighborVertices.push_back(tmpVertex-inputSize[0]-1);
			}
			
			if(neighbors[neighbors8[1]])
			{
				neighborVertices.push_back(tmpVertex-inputSize[0]+1);
			}
			
			if(neighbors[neighbors8[2]])
			{
				neighborVertices.push_back(tmpVertex+inputSize[0]-1);
			}
			
			if(neighbors[neighbors8[3]])
			{
				neighborVertices.push_back(tmpVertex+inputSize[0]+1);
			}
			
			newGraph[tmpVertex] = neighborVertices;
		}
	}
	
	return newGraph;
};


/****************************************************************************/
/*Implementation of Dijkstra's algorithm for a contour graph                */
/*(i.e., a closed circle in the graph is assumed)                           */
/****************************************************************************/
std::list< std::vector< int > > Segmentation::minimalDistanceGraph( GraphType& graph, ImageType::Pointer workingImage, unsigned int samplingRate )
{	
	srand48(std::time(NULL));
	unsigned long firstVertex = 0, lastVertex = 0;
	unsigned long findRouteCount = 0;
	GraphType minGraph;
	do
	{
		GraphType mutableGraph(graph);
		GraphType::iterator graphIter;
		std::list< unsigned long > unvisitedVertices;
		std::vector< unsigned long > allVertices;
		std::map< unsigned long, unsigned long > shortestGraph;
		std::map< unsigned long, unsigned int > distance;
		
		for(graphIter = mutableGraph.begin(); graphIter != mutableGraph.end(); ++graphIter)
		{
			unvisitedVertices.push_back(graphIter->first);
			allVertices.push_back(graphIter->first);
			distance[graphIter->first] = 1E09;
		}
		//initialize moving direction along hte graph by arbitrarily selecting one of the first node's neighbors as end point
		//and disconnect it from the first node. can lead to irregular in case the contour is small and opposite edges 
		//8- or 4-connected
		unsigned long startvertex = 0;
		bool foundGoodStart = 0;
		unsigned long count = 0;
		while(!foundGoodStart && count < 1E06)
		{
			startvertex = (unsigned long)(drand48()*(double)allVertices.size());
			startvertex = allVertices[startvertex];
			
			firstVertex = startvertex;
			if(mutableGraph[firstVertex].size() == 2)
			{
				lastVertex = mutableGraph[firstVertex].back();
				mutableGraph[firstVertex].pop_back();
				distance[firstVertex] = 0;
				foundGoodStart = 1;
			}
			++count;
		}
		if(!foundGoodStart)
		{
			std::cout << "Warning! Could not find starting vertex with 2 neighbor vertices in graph!" << std::endl;
			#ifdef DEBUG
			DebugLog << "Warning! Could not find starting vertex with 2 neighbor vertices in graph!" << std::endl;
			DebugLog << "graph.size() = " << graph.size() << std::endl;
			DebugLog << std::endl;
			#endif
			std::list< std::vector< int > > emptyList;
			return emptyList;
		}
// 		if(foundGoodStart)
// 			std::cout << "found good starting vertex..." << std::endl;
		std::list< unsigned long >::iterator vertexIter;
		std::list< unsigned long >::iterator neighborIter;
		bool at_last_vertex = 0;
		
// 		std::cout << "starting minimal closed loop graph search..." << std::endl;
// 		std::cout << "# of unvisited nodes: " << unvisitedVertices.size() << std::endl;
		unsigned long mainLoopCount = 0;
		while(!unvisitedVertices.empty()/* && !at_last_vertex*/ && mainLoopCount < 1E08)
		{
			unsigned long tmpDist = 1E08;
			unsigned long nextVertex = 0;
			
			for(vertexIter = unvisitedVertices.begin(); vertexIter != unvisitedVertices.end(); ++vertexIter)
			{
				if(distance[*vertexIter] < tmpDist)
				{
					nextVertex = *vertexIter;
					tmpDist = distance[*vertexIter];
				}
			}
			
			unvisitedVertices.remove(nextVertex);
			
			for(neighborIter = mutableGraph[nextVertex].begin(); neighborIter != mutableGraph[nextVertex].end(); ++neighborIter)
			{
				if(std::find(unvisitedVertices.begin(), unvisitedVertices.end(), *neighborIter) != unvisitedVertices.end())
				{
					unsigned int newDist = distance[nextVertex] + 1;
					
					if(newDist < distance[*neighborIter])
					{
						distance[*neighborIter] = newDist;
						shortestGraph[*neighborIter] = nextVertex;
						
// 						if(*neighborIter == lastVertex) at_last_vertex = 1;
						if(*neighborIter == lastVertex)
						{
							at_last_vertex = 1;
							break;
						}
					}
				}
			}
			++mainLoopCount;
		}
		
// 		std::cout << "graph search finished..." << std::endl;
	// 	std::cout << "lastVertex: " << lastVertex << " - shortestGraph[lastVertex]: " << shortestGraph[lastVertex] << " - shortestGraph.size(): " << shortestGraph.size() << std::endl;
		
		unsigned long tmp2 = lastVertex;
		unsigned long shortestGraphCount = 0;
		while(!shortestGraph.empty() && shortestGraphCount < 1E06)
		{
			std::list< unsigned long > nextNeighbor;
			unsigned long tmp = shortestGraph[tmp2];
	// 		if(tmp == firstVertex)
	// 		{
	// 			std::cout << "lastVertex = " << lastVertex << " - firstVertex = " << firstVertex << std::endl;
	// 		}
			if(tmp2 == tmp)
			{
				std::cout << "Error! No closed loop in graph! Check edge image for discontinuities." << std::endl;
				break;
			}
			shortestGraph.erase(tmp2);
			nextNeighbor.push_back(tmp);
			minGraph[tmp2] = nextNeighbor;
			tmp2 = tmp;
			if(tmp2 == firstVertex) break;
			++shortestGraphCount;
		}
		if(shortestGraphCount == 1E06)
		{
			std::cout << "Warning! Could not find closed loop in graph!" << std::endl;
			#ifdef DEBUG
			DebugLog << "Warning! Could not find closed loop in shortestGraph!" << std::endl;
			DebugLog << "graph.size() = " << graph.size() << std::endl;
			DebugLog << "shortestGraph.size() = " << shortestGraph.size() << std::endl;
			DebugLog << std::endl;
			#endif
			std::list< std::vector< int > > emptyList;
			return emptyList;
		}
		
		unvisitedVertices.clear();
		shortestGraph.clear();
		distance.clear();
// 		std::cout << "size of minimal graph = " << minGraph.size() << std::endl;
// 		std::cout << "size of total graph = " << graph.size() << std::endl;
// 		std::cout << "minimal graph created..." << std::endl;
		++findRouteCount;
	} while(minGraph.size() <= 0.5*(float)graph.size() && findRouteCount < 1E04);

	std::list< std::vector< int > > indexList;
	unsigned long count = 0;
	unsigned int x_index = 0, y_index = 0, z_index = 0;
	const ImageType::SizeType inputSize = workingImage->GetLargestPossibleRegion().GetSize();
	const ImageType::IndexType inputIndex = workingImage->GetLargestPossibleRegion().GetIndex();
	unsigned long tmp = lastVertex, next = 0;

// 	std::cout << "graphIter->first = " << graphIter->first << " - tmp = " << tmp << " - graph[tmp].back() = " << graph[tmp].back() << std::endl;
// 	std::cout << "sampling graph with total size " << minGraph.size() << " with sampling rate " << samplingRate << "..." << std::endl;
	unsigned long reverseCount = 0;
	while(!minGraph.empty() && reverseCount < 1E06)
	{
		x_index = tmp%inputSize[0];
		y_index = tmp/inputSize[0];	//integer division!!!
		z_index = inputIndex[2];
		std::vector< int > tmpIndex;
		tmpIndex.push_back(x_index);
		tmpIndex.push_back(y_index);
		tmpIndex.push_back(z_index);
		indexList.push_back(tmpIndex);
		
		next = minGraph[tmp].back();
		// 		std::cout << "next: " << next << std::endl;
		if(tmp == next)
		{
			std::cout << "next: " << next << " - current: " << tmp << " - length of remaining graph: " << minGraph.size() << std::endl;
			break;
		}
		minGraph.erase(tmp);
		tmp = next;
		++reverseCount;
	}
	if(reverseCount == 1E06)
	{
		std::cout << "Warning! Could not find closed loop in minGraph!" << std::endl;
		#ifdef DEBUG
		DebugLog << "Warning! Could not find closed loop in minGraph!" << std::endl;
		DebugLog << "graph.size() = " << graph.size() << std::endl;
		DebugLog << "minGraph.size() = " << minGraph.size() << std::endl;
		DebugLog << std::endl;
		#endif
		std::list< std::vector< int > > emptyList;
		return emptyList;
	}

// 	std::cout << "index list complete!" << std::endl;
	return indexList;
};

/****************************************************************************/
/*Implementation of contour smoothing by averaging x- and y- values of each */
/*index with #(2*radius+1) neighbors going into the averaging               */
/*CAUTION: with large radius, natural curves will be systematically shrunk  */
/****************************************************************************/
void Segmentation::contourSmoothing( std::list< std::vector< float > >& contour, int radius )
{
	const unsigned int contourSize = contour.size();
	if(contourSize <= radius)
		radius = 1;
	std::vector< std::vector< float > > contourVec;
	std::list< std::vector< float > > neighborhood;
	std::list< std::vector< float > >::iterator contourIt;
	std::list< std::vector< float > >::iterator neighborIt;
	
	if(contourSize > 1)
	{
		for(contourIt = contour.begin(); contourIt != contour.end(); ++contourIt)
		{
			contourVec.push_back(*contourIt);
		}
		
		for(int ii = -1*radius; ii <= radius; ++ii)
		{
			if(ii >= 0)
				neighborhood.push_back(contourVec[ii]);
			else
				neighborhood.push_back(contourVec[contourSize - 1 + ii]);
		}
		
		int ii;
		for(ii = 0, contourIt = contour.begin(); ii < contourSize && contourIt != contour.end(); ++ii, ++contourIt)
		{
			float tmp_x = 0, tmp_y = 0;
			
			for(neighborIt = neighborhood.begin(); neighborIt != neighborhood.end(); ++neighborIt)
			{
				tmp_x += (*neighborIt)[0];
				tmp_y += (*neighborIt)[1];
			}
			
			tmp_x /= (float)(2*radius + 1);
			tmp_y /= (float)(2*radius + 1);
			
			(*contourIt)[0] = tmp_x;
			(*contourIt)[1] = tmp_y;
			
			neighborhood.pop_front();
			if(ii + radius < contourSize)
				neighborhood.push_back(contourVec[ii + radius]);
			else
				neighborhood.push_back(contourVec[ii + radius - contourSize]);
		}
	}
};

/****************************************************************************/
/*sampling of contour: only every 'samplingRate'th coordinate is kept       */
/****************************************************************************/
void Segmentation::contourSampling( std::list< std::vector< float > >& contour, int samplingRate )
{
	std::list< std::vector< float > > sampledContour;
	std::list< std::vector< float > >::iterator contourIt;
	int count = 1;
	if(contour.size() <= samplingRate)
		samplingRate = 1;
	
	for(contourIt = contour.begin(); contourIt != contour.end(); ++count)
	{
		if(count == samplingRate)
		{
			count = 0;
			++contourIt;
		}
		else
			contourIt = contour.erase(contourIt);
	}
};

/****************************************************************************/
/*in erronous case where the contour at plane z is empty, set it equal to   */
/*a neighboring (in z) contour                                              */
/****************************************************************************/
void Segmentation::fillEmptyBarrelContour(int barrelID, int z)
{
	int offset = 1;
	if(z >= 0 && z < barrelZContours.size())
	{
		while(!barrelZContours[z][barrelID]->edgeListPointer()->size())
		{
			if(z - offset >= 0 && z - offset < barrelZContours.size())
				if(barrelZContours[z-offset][barrelID]->getValid() && barrelZContours[z-offset][barrelID]->edgeListPointer()->size())
				{
					barrelZContours[z][barrelID]->replaceEdgeList(barrelZContours[z-offset][barrelID]->edgeListPointer(), offset);
					break;
				}
			if(z + offset >= 0 && z + offset < barrelZContours.size())
				if(barrelZContours[z+offset][barrelID]->getValid() && barrelZContours[z+offset][barrelID]->edgeListPointer()->size())
				{
					barrelZContours[z][barrelID]->replaceEdgeList(barrelZContours[z+offset][barrelID]->edgeListPointer(), -offset);
					break;
				}
			if((z - offset < 0) && (z + offset >= barrelZContours.size()))
				break;
			++offset;
		}
	}
};

/****************************************************************************/
/*getWMContour() extracts the WM contour from the foreground image          */
/****************************************************************************/
void Segmentation::getWMContour( float alpha_arg, float beta_arg )
{
	copy(inputImage, originalImage);
	IteratorType2 inputIter(inputImage, maximum_region);
	
// 	std::cout << "Calculating ROI for WM..." << std::endl;
	
	for(inputIter.GoToBegin(); !inputIter.IsAtEnd(); ++inputIter)
	{
		if(inputIter.Get() <= 150)
			inputIter.Set(255);
		else
			inputIter.Set(0);
	}
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_wm_step1_threshold", 0);
	#endif
	
	ImageType::Pointer fgImage = MarkBackground(10000);
// 	writeBinaryImage(bgImage, "_simpleBG", 10);
	fgImage = invertImage(fgImage);
	#ifdef PIPELINE_DOC
	writeBinaryImage(fgImage, "_wm_step2_fgImage", 0);
	#endif
	regionGrowingObjectLabeling(fgImage, 1, 1);
	labeledToGreyscaleImage(fgImage, 1);
	#ifdef PIPELINE_DOC
	writeBinaryImage(fgImage, "_wm_step3_fgImage_region_growing1", 0);
	#endif
	binaryErosion(fgImage, 1, 100);
	#ifdef PIPELINE_DOC
	writeBinaryImage(fgImage, "_wm_step4_fgImage_erosion_R100", 0);
	#endif
	regionGrowingObjectLabeling(fgImage, 1, 1);
	#ifdef PIPELINE_DOC
	writeBinaryImage(fgImage, "_wm_step5_fgImage_region_growing2", 0);
	#endif
	
// 	std::cout << "Adjusting WM intensity..." << std::endl;
	
	unsigned long fgpixel = 0, totalpixel = sizeOfObjects[0];
	float mean = 0, stddev = 0, fgfraction = 0;
	copy(inputImage, originalImage);
	
	unsigned int * fghist;
	fghist = foregroundHistogram(fgImage, inputImage);
	
	for(int ii = 0; ii < 256; ++ii)
	{
		fgpixel += fghist[ii];
		mean += (float)ii*(float)fghist[ii];
	}
	mean = mean/(float)fgpixel;
	
	for(int ii = 0; ii < 256; ++ii )
	{
		stddev += (float)fghist[ii]*((float)ii - mean)*((float)ii - mean);
	}
	stddev = stddev/(float)(fgpixel);
	stddev = std::sqrt(stddev);
	
	int roiThresh = (int)(mean + 0.15*stddev);
	
	ImageType::Pointer wmImage = ImageType::New();
	wmImage->SetRegions(maximum_region);
	wmImage->Allocate();
	ConstIteratorType histIter(inputImage, maximum_region);
	ConstIteratorType fgIter(fgImage, maximum_region);
	IteratorType2 wmIter(wmImage, maximum_region);
	for(histIter.GoToBegin(), fgIter.GoToBegin(), wmIter.GoToBegin(); !histIter.IsAtEnd() && !fgIter.IsAtEnd() && !wmIter.IsAtEnd(); ++histIter, ++fgIter, ++wmIter)
	{
		if(fgIter.Get() && histIter.Get() > roiThresh)
			wmIter.Set(255);
		else
			wmIter.Set(0);
	}
	wmImage->Update();
	#ifdef PIPELINE_DOC
	writeBinaryImage(wmImage, "_wm_step6_wmImage_threshold", 0);
	#endif
	regionGrowingObjectLabeling(wmImage, 1, 0);
	labeledToGreyscaleImage(wmImage, 1);
	#ifdef PIPELINE_DOC
	writeBinaryImage(wmImage, "_wm_step7_wmImage_region_growing1", 0);
	#endif
	
	unsigned int * hist;
	hist = foregroundHistogram(wmImage, inputImage);
	mean = 0, stddev = 0;
	fgpixel = 0;
	
	for(int ii = 0; ii < 256; ++ii)
	{
		fgpixel += hist[ii];
		mean += (float)ii*(float)hist[ii];
	}
	mean = mean/(float)fgpixel;
	
	for(int ii = 0; ii < 256; ++ii )
	{
		stddev += (float)hist[ii]*((float)ii - mean)*((float)ii - mean);
	}
	stddev = stddev/(float)(fgpixel);
	stddev = std::sqrt(stddev);
	
// 	std::cout << "S\thighpxno\ttotalpxno\tmean\tstddev\tmapped mean\tmapped stddev" << std::endl;
// 	std::cout << sliceNumber << "\t" << fgpixel << "\t" << totalpixel << "\t" << mean << "\t" << stddev;
// 	for(int ii = 0; ii < 256; ++ii)
// 	{
// 		totalpixel += hist[ii];
// // 		mean += (float)ii*(float)hist[ii];
// 	}
	fgfraction = (float)fgpixel/(float)totalpixel;
	delete [] hist;
	
	float alpha = alpha_arg*stddev;
	float beta = mean + beta_arg*fgfraction*stddev;
// 	std::cout << "alpha = " << alpha << "\t---\tbeta = " << beta << std::endl;
	
	IteratorType2 sigmoidIter(inputImage, maximum_region);
	ConstIteratorType fgIter2(fgImage, maximum_region);
	
	for(sigmoidIter.GoToBegin(), fgIter2.GoToBegin(); !sigmoidIter.IsAtEnd() && !fgIter2.IsAtEnd(); ++sigmoidIter, ++fgIter2)
	{
		if(fgIter2.Get())
		{
			float tmp = 255.0f/(1 + std::exp((beta - (float)sigmoidIter.Get())/alpha));
			sigmoidIter.Set((int)(tmp + 0.5));
		}
		else
			sigmoidIter.Set(0);
	}
	
	inputImage->Update();
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_wm_step8_sigmoid", 0);
	#endif
// 	writeInputImage("_fg_sigmoid_map.tif");
	unsigned int * hist2;
	hist2 = foregroundHistogram(wmImage, inputImage);
	mean = 0, stddev = 0;
	
	for(int ii = 0; ii < 256; ++ii)
	{
		mean += (float)ii*(float)hist[ii];
	}
	mean = mean/(float)fgpixel;
	
	for(int ii = 0; ii < 256; ++ii )
	{
		stddev += (float)hist[ii]*((float)ii - mean)*((float)ii - mean);
	}
	stddev = stddev/(float)(fgpixel);
	stddev = std::sqrt(stddev);
	
// 	std::cout << "highpxno\tmean\tstddev\ttotalpxno" << std::endl;
// 	std::cout << "\t" << mean << "\t" << stddev << std::endl;
	
	int seed =  (int)(mean + (1.5 + 1.0*fgfraction)*stddev + 0.5);
	int bg = (int)(mean - 0.5*stddev + 0.5);
	delete [] hist2;
// 	std::cout << "seed\tbg" << std::endl;
// 	std::cout << seed << "\t" << bg << std::endl;
	copy(fgImage, inputImage);
	wmImage = regionGrowingNeighborhoodThreshold(0.3, 20, 0.25, seed, bg );	//currently:sigma = 0.3, delta = 20, gamma = 0.25
	#ifdef PIPELINE_DOC
	writeBinaryImage(wmImage, "_wm_step11_region_growing_final", 0);
	#endif
	binaryClosing(wmImage, 1, 25);
	#ifdef PIPELINE_DOC
	writeBinaryImage(wmImage, "_wm_step12_closing_R25", 0);
	#endif
// 	writeBinaryImage(wmImage, "_closing", 0);
	createBvpImage(wmImage, fgImage);
	
	ShapedNeighborhoodIteratorType neighborIter(radius1, wmImage, wmImage->GetLargestPossibleRegion());
	IteratorType2 edgeIter(inputImage, inputImage->GetLargestPossibleRegion());
	ShapedNeighborhoodIteratorType::Iterator it;
	ImageType::IndexType tmpIndex;
	for(int ii = 0; ii < 4; ++ii)
		neighborIter.ActivateOffset(neighborhood4[ii]);
	
	for(neighborIter.GoToBegin(), edgeIter.GoToBegin(); !neighborIter.IsAtEnd() && !edgeIter.IsAtEnd(); ++neighborIter, ++edgeIter)
	{
		if(neighborIter.GetCenterPixel())
		{
			tmpIndex = neighborIter.GetIndex();
			bool neighbor = 0;
			for(it = neighborIter.Begin(); !it.IsAtEnd(); ++it)
			{
				if(!it.Get())
				{
					neighbor = 1;
					break;
				}
			}
			
			if(!neighbor)
				edgeIter.Set(0);
			
			else
				edgeIter.Set(255);
		}
		
		else
			edgeIter.Set(0);
	}
	
	inputImage->Update();
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_wm_step13_edge1", 0);
	#endif
	regionGrowingObjectLabeling(inputImage, 1, 1);
	labeledToGreyscaleImage(inputImage, 1);
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_wm_step14_edge2", 0);
	#endif
// 	writeInputImage("_edgetest.tif");
	
	std::list<std::list<std::vector<float> > > pxcoord_edges;
	pxcoord_edges.push_back(contourExtraction(inputImage));
	contourSmoothing(pxcoord_edges.back(), 25);
	contourSampling(pxcoord_edges.back(), 100);
	// 	std::cout << "writing amira mesh file..." << std::endl;
	std::list<std::list<std::vector<float> > > amira_edges = GetAmiraContours(pxcoord_edges);
	
	std::list<std::vector<float> > vertex_list = SetContourVertices(amira_edges);
	
	std::list<std::vector<float> >::iterator vertexIt;
	std::list<std::list<std::vector<float> > >::iterator edgesIt;
	for(vertexIt = vertex_list.begin(); vertexIt != vertex_list.end(); ++vertexIt)
	{
		amira_contour_graph->vertice_list.push_back(*vertexIt);
	}
	for(edgesIt = amira_edges.begin(); edgesIt != amira_edges.end(); ++edgesIt)
	{
		amira_contour_graph->edge_list.push_back(*edgesIt);
	}
};

/****************************************************************************/
/*region growing of WM from bright starting region; neighborhood of each    */
/*pixel determines whether it belongs to connected component or not         */
/****************************************************************************/
ImageType::Pointer Segmentation::regionGrowingNeighborhoodThreshold( float epsilon, float delta, float gamma, int highThresh, int lowThresh )
{
// 	std::cout << "Region growing preparation..." << std::endl;
	
	medianFilter(5);
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_wm_step9_median", 0);
	#endif
	ImageType::Pointer segmentedImage = ImageType::New();
	segmentedImage->SetRegions(maximum_region);
	segmentedImage->Allocate();
	segmentedImage->FillBuffer(0);
	ImageType::Pointer otsuImage = ImageType::New();
	otsuImage->SetRegions(maximum_region);
	otsuImage->Allocate();
	otsuImage->FillBuffer(0);
	std::list< ImageType::IndexType > queue;
	
	//ensure that seed is not blood vessel
	//but actually part of the largest bright object (i.e. WM)
	int * completeHist = completeHistogram(inputImage);
	unsigned int otsuThresh = otsuThreshold(completeHist);
	
	ConstIndexIteratorType threshIter(inputImage, maximum_region);
	IteratorType2 otsuIter(otsuImage, maximum_region);
	for(otsuIter.GoToBegin(), threshIter.GoToBegin(); !otsuIter.IsAtEnd() && !threshIter.IsAtEnd(); ++otsuIter, ++threshIter)
		if(threshIter.Get() > otsuThresh)
			otsuIter.Set(255);
	
	regionGrowingObjectLabeling(otsuImage, 1, 0);
	labeledToGreyscaleImage(otsuImage, 1);
	
	bool foundSeed = 0;
	int seedThresh = highThresh;
	while(!foundSeed)
	{
// 		std::flush(std::cout << "seedThresh = " << seedThresh << std::endl);
		IteratorType2 segIter(segmentedImage, maximum_region);
		ConstIteratorType seedThreshIter(inputImage, maximum_region);
		ConstIteratorType wmObjectIter(otsuImage, maximum_region);
		for(segIter.GoToBegin(), seedThreshIter.GoToBegin(), wmObjectIter.GoToBegin();
		!segIter.IsAtEnd() && !seedThreshIter.IsAtEnd() && !wmObjectIter.IsAtEnd();
		++segIter, ++seedThreshIter, ++wmObjectIter)
		{
			if(seedThreshIter.Get() > seedThresh && wmObjectIter.Get())
				segIter.Set(255);
		}
		segmentedImage->Update();
		
		regionGrowingObjectLabeling(segmentedImage, 1, 0);
		labeledToGreyscaleImage(segmentedImage, 1);
		if(nrOfObjects)
		{
// 			std::flush(std::cout << "nrOfObjects = " << nrOfObjects << std::endl);
// 			std::flush(std::cout << "sizeOfObjects[0] = " << sizeOfObjects[0] << std::endl);
			if(sizeOfObjects[0] > 1500)
				break;
		}
		--seedThresh;
		IteratorType2 segResetIter(segmentedImage, maximum_region);
		for(segResetIter.GoToBegin(); !segResetIter.IsAtEnd(); ++segResetIter)
			segResetIter.Set(0);
	}
	#ifdef PIPELINE_DOC
	writeBinaryImage(segmentedImage, "_wm_step10_seed", 0);
	#endif
	
// 	regionGrowingObjectLabeling(segmentedImage, 1, 0);
// 	labeledToGreyscaleImage(segmentedImage, 1);
// 	writeBinaryImage(segmentedImage, "_seed_region", 0);
	
	IndexIteratorType segIter2(segmentedImage, maximum_region);
	ImageType::SizeType planeRadius1;
	planeRadius1[0] = 1;
	planeRadius1[1] = 1;
	planeRadius1[2] = 0;
	SegNeighborhoodIteratorType neighborIter(planeRadius1, segmentedImage, maximum_region);
	int seedNeighbors[] = {1,3,5,7};
	
	for(segIter2.GoToBegin(), neighborIter.GoToBegin(); !segIter2.IsAtEnd() && !neighborIter.IsAtEnd(); ++segIter2, ++neighborIter)
	{
		if(segIter2.Get())
		{
// 			std::cout << "Checking index " << segIter2.GetIndex() << std::endl;
			for(int ii = 0; ii < 4; ++ii)
			{
// 				std::cout << "checking neighbor " << seedNeighbors[ii] << " @ index " << neighborIter.GetIndex(seedNeighbors[ii]) << std::endl;
				if(!neighborIter.GetPixel(seedNeighbors[ii]))
					queue.push_back(neighborIter.GetIndex(seedNeighbors[ii]));
			}
		}
	}
// 	std::cout << "queue.size() = " << queue.size() << std::endl;
	std::list< ImageType::IndexType >::iterator queueIt;
	for(queueIt = queue.begin(); queueIt != queue.end(); ++queueIt)
	{
		segIter2.SetIndex(queue.front());
		segIter2.Set(1);
	}
	
	ImageType::SizeType radius;
	radius[0] = 7;
	radius[1] = 7;
	radius[2] = 0;
	SegNeighborhoodIteratorType rgIter(radius, inputImage, maximum_region);
	int highThreshold = highThresh;
	int lowThreshold = lowThresh;
	int neighbors4[] = {-1*(int)(2*radius[0]+1), (int)(2*radius[1]+1), -1, 1};
// 	std::cout << "neighbors4: [" << neighbors4[0] << "," << neighbors4[1] << "," << neighbors4[2] << "," << neighbors4[3] << "]" << std::endl;
// 	segmentedImage->FillBuffer(0);

// 	std::flush(std::cout << "Starting region growing algorithm...");
	
	unsigned long count = 0;
	while(!queue.empty())
	{
		rgIter.SetLocation(queue.front());
		
		if(rgIter.GetCenterPixel() > highThreshold)
		{
			segIter2.SetIndex(queue.front());
			segIter2.Set(255);
			unsigned int centerOffset = rgIter.GetCenterNeighborhoodIndex();
// 			std::cout << "high threshold ok" << std::endl;
// 			std::cout << "centerOffset = " << centerOffset << std::endl;
// 			std::cout << "current index = " << queue.front() << std::endl;
			for(int ii = 0; ii < 4; ++ii)
			{
// 				std::cout << "neighbor index " << ii << " = " << rgIter.GetIndex(centerOffset + neighbors4[ii]) << std::endl;
				segIter2.SetIndex(rgIter.GetIndex(centerOffset + neighbors4[ii]));
				if(!segIter2.Get()/* && count < 1E06*/)
				{
					queue.push_back(segIter2.GetIndex());
					segIter2.Set(1);
				}
			}
		}
		else if(rgIter.GetCenterPixel() < lowThreshold)
		{
// 			std::cout << "low threshold not ok" << std::endl;
// 			std::cout << "current index = " << queue.front() << std::endl;
			segIter2.SetIndex(queue.front());
			segIter2.Set(128);
		}
		else if(rgIter.GetCenterPixel() >= lowThreshold && rgIter.GetCenterPixel() <= highThreshold)
		{
			float mean = 0, stddev = 0;
			float neighborhoodSize = (2*radius[0]+1)*(2*radius[1]+1)*(2*radius[2]+1);
			for(int ii = 0; ii < neighborhoodSize; ++ii)
			{
				mean += rgIter.GetPixel(ii);
			}
			
			mean /= neighborhoodSize;
			
			for(int ii = 0; ii < neighborhoodSize; ++ii)
			{
				stddev += (rgIter.GetPixel(ii) - mean)*(rgIter.GetPixel(ii) - mean);
			}
			stddev /= neighborhoodSize;
			stddev = std::sqrt(stddev);
			
			int neighborThresh = (int)(mean + epsilon*stddev/**stddev*/ + 0.5);
			unsigned int highPx = 0;
			
			for(int ii = 0; ii < neighborhoodSize; ++ii)
			{
				if(rgIter.GetPixel(ii) > neighborThresh)
					++highPx;
			}
			
			if(mean > lowThreshold + delta && highPx > gamma*neighborhoodSize)
			{
				segIter2.SetIndex(queue.front());
				segIter2.Set(255);
				unsigned int centerOffset = rgIter.GetCenterNeighborhoodIndex();
// 				std::cout << "medium threshold ok" << std::endl;
// 				std::cout << "centerOffset = " << centerOffset << std::endl;
// 				std::cout << "current index = " << queue.front() << std::endl;
				for(int ii = 0; ii < 4; ++ii)
				{
// 					std::cout << "neighbor index " << ii << " = " << rgIter.GetIndex(centerOffset + neighbors4[ii]) << std::endl;
					segIter2.SetIndex(rgIter.GetIndex(centerOffset + neighbors4[ii]));
					if(!segIter2.Get()/* && count < 1E06*/)
					{
						queue.push_back(segIter2.GetIndex());
						segIter2.Set(1);
					}
				}
			}
			else
			{
// 				std::cout << "medium threshold not ok" << std::endl;
// 				std::cout << "current index = " << queue.front() << std::endl;
				segIter2.SetIndex(queue.front());
				segIter2.Set(128);
			}
		}
		
		queue.pop_front();
		++count;
	}
	
// 	std::cout << " complete!" << std::endl;
// 	std::cout << "count = " << count << std::endl;
	
	segmentedImage->Update();
	
	for(segIter2.GoToBegin(); !segIter2.IsAtEnd(); ++segIter2)
	{
		if(segIter2.Get() != 255)
			segIter2.Set(0);
	}
	segmentedImage->Update();
// 	writeBinaryImage(segmentedImage, "_region_growing", 0);
	return segmentedImage;
};

void Segmentation::createBvpImage ( ImageType::Pointer wmImage, ImageType::Pointer fgImage )
{
	wmBvpImage = ImageType::New();
	wmBvpImage->SetRegions(wmImage->GetLargestPossibleRegion());
	wmBvpImage->Allocate();
	copy(wmBvpImage, wmImage);
	binaryDilation(wmBvpImage, 1, 100);
	
	ConstIteratorType readIter(fgImage, fgImage->GetLargestPossibleRegion());
	IteratorType2 setIter(wmBvpImage, wmBvpImage->GetLargestPossibleRegion());
	for(readIter.GoToBegin(), setIter.GoToBegin(); !readIter.IsAtEnd() && !setIter.IsAtEnd(); ++readIter, ++setIter)
	{
		if(!setIter.Get())
			setIter.Set(readIter.Get());
		else
			setIter.Set(0);
	}
	wmBvpImage->Update();
	#ifdef PIPELINE_DOC
	writeBinaryImage(wmBvpImage, "_bvp_image", 0);
	#endif
};

/*******************************************************/
/*binary3DOpening(image)                               */
/*******************************************************/
void Segmentation::binaryOpening( ImageType::Pointer workingImage, bool binaryBall, int size )
{
	ImageType::Pointer currImage = ImageType::New();
	currImage->SetRegions(workingImage->GetLargestPossibleRegion() );
	currImage->Allocate();
	currImage->FillBuffer(0);
	currImage->Update();
	
	BinaryErodeFilterBallType::Pointer erodeBallFilter = BinaryErodeFilterBallType::New();
	BinaryDilateFilterBallType::Pointer dilateBallFilter = BinaryDilateFilterBallType::New();
	StructuringElementBallType structuringBallElement;
	
	BinaryErodeFilterCrossType::Pointer erodeCrossFilter = BinaryErodeFilterCrossType::New();
	BinaryDilateFilterCrossType::Pointer dilateCrossFilter = BinaryDilateFilterCrossType::New();
	StructuringElementCrossType structuringCrossElement;
	
	if(binaryBall)
	{
// 		std::cout << "Opening with binary ball structuring element!" << std::endl;
		ImageType::SizeType radius;
		radius[0] = size;
		radius[1] = size;
		radius[2] = 0;
		structuringBallElement.SetRadius(radius);
		structuringBallElement.CreateStructuringElement();
		
		erodeBallFilter->SetKernel(structuringBallElement);
		erodeBallFilter->SetInput(workingImage);
		erodeBallFilter->SetErodeValue(255);
		erodeBallFilter->Update();
		
		dilateBallFilter->SetKernel(structuringBallElement);
		dilateBallFilter->SetInput(erodeBallFilter->GetOutput());
		dilateBallFilter->SetDilateValue(255);
		dilateBallFilter->Update();
		
		currImage = dilateBallFilter->GetOutput();
	}
	
	else // binary cross is DEFINED by ITK to have a radius of 1x1x1
	{
// 		std::cout << "Opening with binary cross structuring element!" << std::endl;
		ImageType::SizeType radius;
		radius[0] = size;
		radius[1] = size;
		radius[2] = 0;
		structuringCrossElement.SetRadius(radius);
		structuringCrossElement.CreateStructuringElement();
// 		std::cout << "structuring element radius: " << structuringCrossElement.GetRadius() << std::endl;
		
		erodeCrossFilter->SetKernel(structuringCrossElement);
		erodeCrossFilter->SetInput(workingImage);
		erodeCrossFilter->SetErodeValue(255);
		erodeCrossFilter->Update();
		
		dilateCrossFilter->SetKernel(structuringCrossElement);
		dilateCrossFilter->SetInput(erodeCrossFilter->GetOutput());
		dilateCrossFilter->SetDilateValue(255);
		dilateCrossFilter->Update();
		
		currImage = dilateCrossFilter->GetOutput();
	}
	
	currImage->Update();
	
	workingImage->FillBuffer(0);
	workingImage->Update();
	copy( workingImage, currImage );
	currImage=0;
};

/*******************************************************/
/*binary3DClosing(image)                               */
/*******************************************************/
void Segmentation::binaryClosing( ImageType::Pointer workingImage, bool binaryBall, int size )
{
	ImageType::Pointer currImage = ImageType::New();
	currImage->SetRegions(workingImage->GetLargestPossibleRegion() );
	currImage->Allocate();
	currImage->FillBuffer(0);
	currImage->Update();
	
	BinaryErodeFilterBallType::Pointer erodeBallFilter = BinaryErodeFilterBallType::New();
	BinaryDilateFilterBallType::Pointer dilateBallFilter = BinaryDilateFilterBallType::New();
	StructuringElementBallType structuringBallElement;
	
	BinaryErodeFilterCrossType::Pointer erodeCrossFilter = BinaryErodeFilterCrossType::New();
	BinaryDilateFilterCrossType::Pointer dilateCrossFilter = BinaryDilateFilterCrossType::New();
	StructuringElementCrossType structuringCrossElement;
	
	if(binaryBall)
	{
		// 		std::cout << "Opening with binary ball structuring element!" << std::endl;
		ImageType::SizeType radius;
		radius[0] = size;
		radius[1] = size;
		radius[2] = 0;
		structuringBallElement.SetRadius(radius);
		structuringBallElement.CreateStructuringElement();
		
		dilateBallFilter->SetKernel(structuringBallElement);
		dilateBallFilter->SetInput(workingImage);
		dilateBallFilter->SetDilateValue(255);
		dilateBallFilter->Update();
		
		erodeBallFilter->SetKernel(structuringBallElement);
		erodeBallFilter->SetInput(dilateBallFilter->GetOutput());
		erodeBallFilter->SetErodeValue(255);
		erodeBallFilter->Update();
		
		currImage = erodeBallFilter->GetOutput();
	}
	
	else // binary cross is DEFINED by ITK to have a radius of 1x1x1
	{
		// 		std::cout << "Opening with binary cross structuring element!" << std::endl;
		ImageType::SizeType radius;
		radius[0] = size;
		radius[1] = size;
		radius[2] = 0;
		structuringCrossElement.SetRadius(radius);
		structuringCrossElement.CreateStructuringElement();
		// 		std::cout << "structuring element radius: " << structuringCrossElement.GetRadius() << std::endl;
		
		dilateCrossFilter->SetKernel(structuringCrossElement);
		dilateCrossFilter->SetInput(workingImage);
		dilateCrossFilter->SetDilateValue(255);
		dilateCrossFilter->Update();
		
		erodeCrossFilter->SetKernel(structuringCrossElement);
		erodeCrossFilter->SetInput(dilateCrossFilter->GetOutput());
		erodeCrossFilter->SetErodeValue(255);
		erodeCrossFilter->Update();
		
		currImage = erodeCrossFilter->GetOutput();
	}
	
	currImage->Update();
	
	workingImage->FillBuffer(0);
	workingImage->Update();
	copy( workingImage, currImage );
	currImage=0;
};

/*******************************************************/
/*binaryErosion(image)                                 */
/*******************************************************/
void Segmentation::binaryErosion( ImageType::Pointer workingImage, bool binaryBall, int size )
{
	ImageType::Pointer currImage = ImageType::New();
	currImage->SetRegions(workingImage->GetLargestPossibleRegion() );
	currImage->Allocate();
	currImage->FillBuffer(0);
	currImage->Update();
	
	BinaryErodeFilterBallType::Pointer erodeBallFilter = BinaryErodeFilterBallType::New();
	StructuringElementBallType structuringBallElement;
	
	BinaryErodeFilterCrossType::Pointer erodeCrossFilter = BinaryErodeFilterCrossType::New();
	StructuringElementCrossType structuringCrossElement;
	
	if(binaryBall)
	{
		// 		std::cout << "Opening with binary ball structuring element!" << std::endl;
		ImageType::SizeType radius;
		radius[0] = size;
		radius[1] = size;
		radius[2] = 0;
		structuringBallElement.SetRadius(radius);
		structuringBallElement.CreateStructuringElement();
		
		erodeBallFilter->SetKernel(structuringBallElement);
		erodeBallFilter->SetInput(workingImage);
		erodeBallFilter->SetErodeValue(255);
		erodeBallFilter->Update();
		
		currImage = erodeBallFilter->GetOutput();
	}
	
	else // binary cross is DEFINED by ITK to have a radius of 1x1x1
	{
		// 		std::cout << "Opening with binary cross structuring element!" << std::endl;
		ImageType::SizeType radius;
		radius[0] = size;
		radius[1] = size;
		radius[2] = 0;
		structuringCrossElement.SetRadius(radius);
		structuringCrossElement.CreateStructuringElement();
		
		erodeCrossFilter->SetKernel(structuringCrossElement);
		erodeCrossFilter->SetInput(workingImage);
		erodeCrossFilter->SetErodeValue(255);
		erodeCrossFilter->Update();
		
		currImage = erodeCrossFilter->GetOutput();
	}
	
	currImage->Update();
	
	workingImage->FillBuffer(0);
	workingImage->Update();
	copy( workingImage, currImage );
	currImage=0;
};

/*******************************************************/
/*binaryDilation(image)                                */
/*******************************************************/
void Segmentation::binaryDilation( ImageType::Pointer workingImage, bool binaryBall, int size )
{
	ImageType::Pointer currImage = ImageType::New();
	currImage->SetRegions(workingImage->GetLargestPossibleRegion() );
	currImage->Allocate();
	currImage->FillBuffer(0);
	currImage->Update();
	
	BinaryDilateFilterBallType::Pointer dilateBallFilter = BinaryDilateFilterBallType::New();
	StructuringElementBallType structuringBallElement;
	
	BinaryDilateFilterCrossType::Pointer dilateCrossFilter = BinaryDilateFilterCrossType::New();
	StructuringElementCrossType structuringCrossElement;
	
	if(binaryBall)
	{
		// 		std::cout << "Opening with binary ball structuring element!" << std::endl;
		ImageType::SizeType radius;
		radius[0] = size;
		radius[1] = size;
		radius[2] = 0;
		structuringBallElement.SetRadius(radius);
		structuringBallElement.CreateStructuringElement();
		
		dilateBallFilter->SetKernel(structuringBallElement);
		dilateBallFilter->SetInput(workingImage);
		dilateBallFilter->SetDilateValue(255);
		dilateBallFilter->Update();
		
		currImage = dilateBallFilter->GetOutput();
	}
	
	else // binary cross is DEFINED by ITK to have a radius of 1x1x1
	{
		// 		std::cout << "Opening with binary cross structuring element!" << std::endl;
		ImageType::SizeType radius;
		radius[0] = size;
		radius[1] = size;
		radius[2] = 0;
		structuringCrossElement.SetRadius(radius);
		structuringCrossElement.CreateStructuringElement();
		
		dilateCrossFilter->SetKernel(structuringCrossElement);
		dilateCrossFilter->SetInput(workingImage);
		dilateCrossFilter->SetDilateValue(255);
		dilateCrossFilter->Update();
		
		currImage = dilateCrossFilter->GetOutput();
	}
	
	currImage->Update();
	
	workingImage->FillBuffer(0);
	workingImage->Update();
	copy( workingImage, currImage );
	currImage=0;
};

/****************************************************************************/
/*substractImages(a,b,c) gives image c with c=b-a                           */
/****************************************************************************/
void Segmentation::substractImages(ImageType::Pointer image1, ImageType::Pointer image2, ImageType::Pointer image3)
{
	ConstIteratorType it1( image1, image1->GetLargestPossibleRegion() );
	ConstIteratorType it2( image2, image2->GetLargestPossibleRegion() );
	IteratorType2 it3( image3, image3->GetLargestPossibleRegion() );
	
	for( it1.GoToBegin(), it2.GoToBegin(), it3.GoToBegin() ; !it1.IsAtEnd(); ++it1, ++it2, ++it3)
	{
		it3.Set( it2.Get()-it1.Get() );
	}
	
	image3->Update();
};

/****************************************************************************/
/*regionGrowingObjectLabeling()                                             */
/****************************************************************************/
void Segmentation::regionGrowingObjectLabeling( ImageType::Pointer workingImage, int minimum_size, bool neighborhood26 = 0 )
{
	ConnectedFilterType::Pointer labeller = ConnectedFilterType::New();
	labeller->SetInput(workingImage);
	if(neighborhood26)
	{
// 		std::cout << "Setting connectivity for objectImage to 26-connected!" << std::endl;
		labeller->FullyConnectedOn();
	}
	labeller->Update();
	
	RelabelType::Pointer relabeller = RelabelType::New();
	relabeller->SetInput(labeller->GetOutput());
	relabeller->SetMinimumObjectSize( minimum_size );
// 	relabeller->Update();
	
	labeledImage = ObjectImageType::New();
	labeledImage = relabeller->GetOutput();
// 	labeledImage = labeller->GetOutput();
	labeledImage->Update();
	
	nrOfObjects = relabeller->GetNumberOfObjects();
	//unsigned long nrOfOriginalObjects = relabeller->GetOriginalNumberOfObjects();
	
	
	const std::vector<unsigned long> tmpObjectSize = relabeller->GetSizeOfObjectsInPixels();
	
	sizeOfObjects = tmpObjectSize;
};

/****************************************************************************/
/*regionGrowingObjectLabeling()                                             */
/****************************************************************************/
void Segmentation::relabelImage( ObjectImageType::Pointer workingImage, int minimum_size )
{
	RelabelType::Pointer relabeller = RelabelType::New();
	relabeller->SetInput(workingImage);
	relabeller->SetMinimumObjectSize( minimum_size );
// 	relabeller->Update();
	
	labeledImage = relabeller->GetOutput();
	labeledImage->Update();
	
	nrOfObjects = relabeller->GetNumberOfObjects();
	//unsigned long nrOfOriginalObjects = relabeller->GetOriginalNumberOfObjects();
	
	
	const std::vector<unsigned long> tmpObjectSize = relabeller->GetSizeOfObjectsInPixels();
	
	sizeOfObjects = tmpObjectSize;
};

/****************************************************************************/
/*labeledToGreyscaleImage()                                                 */
/****************************************************************************/
void Segmentation::labeledToGreyscaleImage( ImageType::Pointer& workingImage, unsigned long upperThreshold )
{
	ObjectFilterType::Pointer objectFilter = ObjectFilterType::New();
	
	objectFilter->SetInput( labeledImage );
	objectFilter->SetOutsideValue ( 0 );
	objectFilter->SetInsideValue ( 255 );
	
	objectFilter->SetLowerThreshold ( 1 );
	objectFilter->SetUpperThreshold ( upperThreshold );
// 	objectFilter->Update();
	
	workingImage = objectFilter->GetOutput();
	workingImage->Update();
};

/****************************************************************************/
/*segmentation for extracting the barrel contours                           */
/****************************************************************************/
void Segmentation::segmentBarrels(unsigned int noOfPoints, std::vector< Contour * > barrelContours, int z, bool alreadyPreProcessed)
{
	#ifdef DEBUG
	DebugLog << "Plane " << z << " processing..." << std::endl;
	#endif
	if(!alreadyPreProcessed)
	{
		ImageType::Pointer holeImage = ImageType::New();
		holeImage->SetRegions(inputImage->GetLargestPossibleRegion());
		holeImage->Allocate();
		
		ConstIteratorType readIter(inputImage, inputImage->GetLargestPossibleRegion());
		IteratorType2 threshIter(holeImage, holeImage->GetLargestPossibleRegion());
		for(threshIter.GoToBegin(), readIter.GoToBegin(); !threshIter.IsAtEnd() && !readIter.IsAtEnd(); ++threshIter, ++readIter)
		{
			if(readIter.Get() > 250)
				threshIter.Set(255);
			else
				threshIter.Set(0);
		}
		holeImage->Update();
		#ifdef PIPELINE_DOC
		writeBinaryImage(holeImage, "_barrels_step1_holes", 0);
		#endif
		
		#ifdef DEBUG
		DebugLog << "Median filter running..." << std::endl;
		#endif
		medianFilter(5);
		#ifdef PIPELINE_DOC
		writeBinaryImage(inputImage, "_barrels_step2_median_R5", 0);
		#endif
		#ifdef DEBUG
		writeInputImage("_xy_median.tif");
		#endif
		
	// 	fillHoles(holeImage);
		
		morphologicalForegroundRemoval( 5, 0.5, 0 );
		
		#ifdef DEBUG
		DebugLog << "Running bilateral filter..." << std::endl;
		#endif
		
		input2DImage = Image2DType::New();
		Image2DType::RegionType plane2DRegion;
		Image2DType::SizeType plane2DSize;
		Image2DType::IndexType plane2DIndex;
		
		plane2DSize[0] = maximum_region.GetSize()[0];
		plane2DSize[1] = maximum_region.GetSize()[1];
		plane2DIndex[0] = maximum_region.GetIndex()[0];
		plane2DIndex[1] = maximum_region.GetIndex()[1];
		plane2DRegion.SetIndex(plane2DIndex);
		plane2DRegion.SetSize(plane2DSize);
		input2DImage->SetRegions(plane2DRegion);
		input2DImage->Allocate();
		copy(input2DImage, inputImage, maximum_region);
		
		ImageToCalc2DRescaleFilterType::Pointer iToCalcFilter = ImageToCalc2DRescaleFilterType::New();
		CalcToImage2DRescaleFilterType::Pointer calcToIFilter = CalcToImage2DRescaleFilterType::New();
		BilateralImage2DFilterType::Pointer bilateralFilter = BilateralImage2DFilterType::New();
		
		iToCalcFilter->SetOutputMinimum(0);
		iToCalcFilter->SetOutputMaximum(255);
		iToCalcFilter->SetInput(input2DImage);
		
		BilateralImage2DFilterType::ArrayType sigma2D;
		for(int ii = 0; ii < 2; ++ii)
			sigma2D[ii] = 5;
		bilateralFilter->SetDomainSigma(sigma2D);
		bilateralFilter->SetRangeSigma(20);
		bilateralFilter->SetInput(iToCalcFilter->GetOutput());
	// 	bilateralFilter->SetDebug(1);
		
		calcToIFilter->SetOutputMinimum(0);
		calcToIFilter->SetOutputMaximum(255);
		calcToIFilter->SetInput(bilateralFilter->GetOutput());
		
		input2DImage = calcToIFilter->GetOutput();
		input2DImage->Update();
		#ifdef DEBUG
		writeInput2DImage("_bilateral_filter.tif");
		#endif
	// 	std::flush(std::cout << "paste 2D to 3D image..." << std::endl);
		copy(inputImage, input2DImage, maximum_region);
		#ifdef PIPELINE_DOC
		writeBinaryImage(inputImage, "_barrels_step5_bilateral_filter", 0);
		#endif
		
		#ifdef DEBUG
		DebugLog << "filling holes in foreground image..." << std::endl;
		#endif
		
// 		fgRemovedImage = ImageType::New();
// 		fgRemovedImage->SetRegions(maximum_region);
// 		fgRemovedImage->Allocate();
		copy(fgRemovedImage, inputImage);
		ImageType::Pointer foregroundImage = ImageType::New();
		foregroundImage->SetRegions(maximum_region);
		foregroundImage->Allocate();
		IteratorType2 fgReadIter(inputImage, maximum_region);
		IteratorType2 fgWriteIter(foregroundImage, maximum_region);
		for(fgReadIter.GoToBegin(), fgWriteIter.GoToBegin(); !fgReadIter.IsAtEnd() && !fgWriteIter.IsAtEnd(); ++fgReadIter, ++fgWriteIter)
		{
			fgWriteIter.Set(fgReadIter.Get() >= 250 ? 255 : 0);
			if(fgReadIter.Get() >= 250)
				fgReadIter.Set(255);
		}
		
		foregroundImage->Update();
		#ifdef PIPELINE_DOC
		writeBinaryImage(foregroundImage, "_barrels_step6_foreground_image", 0);
		#endif
		fillHoles2(foregroundImage);
	// 	writeInputImage("_foreground_removed.tif");
	// 	ImageType::Pointer invertedImage = invertImage(foregroundImage);
	// 	foregroundHistogramEqualization(/*2, 2*/inputImage, inputImage);
		
	// 	sigmoidIntensityMapping(0.4, -0.85, 2);
	// 	writeInputImage("_fglinear.tif");
	// 	std::cout << "Calculating regional histograms..." << std::endl;
	}
// 	foregroundHistogramIntensityMapping(0.45, -0.3, 1);
	else
	{
// 		fgRemovedImage = ImageType::New();
// 		fgRemovedImage->SetRegions(maximum_region);
// 		fgRemovedImage->Allocate();
		#ifdef DEBUG
		DebugLog << "filling holes in foreground image..." << std::endl;
		#endif
		copy(fgRemovedImage, inputImage);
		ImageType::Pointer foregroundImage = ImageType::New();
		foregroundImage->SetRegions(maximum_region);
		foregroundImage->Allocate();
		IteratorType2 fgReadIter(inputImage, maximum_region);
		IteratorType2 fgWriteIter(foregroundImage, maximum_region);
		for(fgReadIter.GoToBegin(), fgWriteIter.GoToBegin(); !fgReadIter.IsAtEnd() && !fgWriteIter.IsAtEnd(); ++fgReadIter, ++fgWriteIter)
		{
			fgWriteIter.Set(fgReadIter.Get() >= 250 ? 255 : 0);
			if(fgReadIter.Get() >= 250)
				fgReadIter.Set(255);
		}
		
		foregroundImage->Update();
		fillHoles2(foregroundImage);
	}
	
	#ifdef DEBUG
	DebugLog << "Sigmoid intensity mapping..." << std::endl;
	#endif
	sigmoidIntensityMapping(0.4, -0.85, 2);
	#ifdef PIPELINE_DOC
	writeBinaryImage(sigmoidImage, "_barrels_step9_sigmoid", 0);
	#endif
// 	writeBinaryImage(sigmoidImage, "_sigmoid_image", z);
// 	preprocImage = ImageType::New();
// 	preprocImage->SetRegions(inputImage->GetLargestPossibleRegion());
// 	preprocImage->Allocate();
	copy(preprocImage, inputImage);
// 	std::cout << "nrOfPoints = " << noOfPoints << std::endl;
// 	std::cout << "barrelContours[z].size() = " << barrelContours.size() << std::endl;
	voronoiRegionThreshold(noOfPoints, barrelContours, z);
	#ifdef PIPELINE_DOC
	writeBinaryImage(segmentedImage, "_barrels_step10_voronoi_region_growing_final", 0);
	#endif
	
	ImageType::RegionType planeRegion;
// 	ImageType::RegionType planeImageRegion;
// 	ImageType::IndexType planeImageIndex;
	ImageType::IndexType planeIndex;
	ImageType::SizeType planeSize;
// 	planeImageIndex = originalImage->GetLargestPossibleRegion().GetIndex();
	planeIndex = originalImage->GetLargestPossibleRegion().GetIndex();
	planeIndex[2] = z;
	planeSize = originalImage->GetLargestPossibleRegion().GetSize();
	planeSize[2] = 1;
// 	planeImageRegion.SetIndex(planeImageIndex);
// 	planeImageRegion.SetSize(planeSize);
	planeRegion.SetIndex(planeIndex);
	planeRegion.SetSize(planeSize);
// 	copy(sigmoidImageStack, sigmoidImage, planeRegion, 1);
	copy(preprocImageStack, preprocImage, planeRegion, 1);
	copy(segmentedImageStack, segmentedImage, planeRegion, 1);
	
// 	voronoiRegionIdentification(noOfPoints);
// 	localRegionIdentification(3, 3);
// 	fillBinaryHoles();
};

/****************************************************************************/
/*remove cell bodies with morphological filtering                           */
/****************************************************************************/
void Segmentation::morphologicalForegroundRemoval( float rad, float arg2, float arg3 )
{
	ImageType::Pointer openImage = ImageType::New();
	openImage->SetRegions(inputImage->GetLargestPossibleRegion());
	openImage->Allocate();
	copy(openImage, inputImage);
	
	OpeningFilterType::Pointer grayscaleOpenFilter = OpeningFilterType::New();
	StructuringElementBallType binaryBall;
	ImageType::SizeType radius;
	radius[0] = (unsigned int)rad;
	radius[1] = (unsigned int)rad;
	radius[2] = 0;
	
	binaryBall.SetRadius(radius);
	binaryBall.CreateStructuringElement();
// 	std::flush(std::cout << "Structuring element size = " << binaryBall.GetSize() << std::endl);
// 	std::flush(std::cout << "Grayscale opening running...");
	grayscaleOpenFilter->SetKernel(binaryBall);
	grayscaleOpenFilter->SetInput(openImage);
	openImage = grayscaleOpenFilter->GetOutput();
	openImage->Update();
// 	std::flush(std::cout << "finished!" << std::endl);
	ImageType::Pointer tophatImage = ImageType::New();
	tophatImage->SetRegions(inputImage->GetLargestPossibleRegion());
	tophatImage->Allocate();
	tophatImage->FillBuffer(0);
// 	std::flush(std::cout << "Calculating top hat image..." << std::endl);
	substractImages(openImage, inputImage, tophatImage);
	#ifdef PIPELINE_DOC
	writeBinaryImage(tophatImage, "_barrels_step3_tophat_R5", 0);
	#endif
	
	float * tophatHist = new float[2];
	Histogram(tophatHist, tophatImage);
	
	int threshold = (int)(tophatHist[0] + arg2*tophatHist[1] + 0.5);
// 	std::cout << "Top Hat image threshold = " << threshold << std::endl;
	ConstIteratorType threshIter(tophatImage, tophatImage->GetLargestPossibleRegion());
	IteratorType2 markIter(inputImage, inputImage->GetLargestPossibleRegion());
	for(threshIter.GoToBegin(), markIter.GoToBegin(); !threshIter.IsAtEnd() && !markIter.IsAtEnd(); ++threshIter, ++markIter)
	{
		if(threshIter.Get() >= threshold)
			markIter.Set(255);
// 		else
// 			threshIter.Set(0);
	}
	inputImage->Update();
	#ifdef PIPELINE_DOC
	writeBinaryImage(inputImage, "_barrels_step4_fg_markers", 0);
	#endif
// 	tophatImage->Update();
	
// 	fillHoles2(tophatImage);
// 	writeInputImage("_tophat_filled2.tif");
	
// 	writeBinaryImage(openImage, "_open", 0);
};

/****************************************************************************/
/*fillHoles2() locates holes resulting from cropping other color channels   */
/*from the current image and fills the border pixels iteratively            */
/*with the mean of the surrounding pixels                                   */
/****************************************************************************/
void Segmentation::fillHoles2( ImageType::Pointer artifactImage )
{
// 	std::cout << "Filling holes!" << std::endl;
	
// 	setArtifactsToZero(inputImage);
	enlargeHoles( artifactImage );
	#ifdef PIPELINE_DOC
	writeBinaryImage(artifactImage, "_barrels_step7_enlarge_holes_foreground_image", 0);
	#endif
	
	regionGrowingObjectLabeling( artifactImage, 5, 0 );
// 	labeledToGreyscaleImage(artifactImage, nrOfObjects);
	#ifdef DEBUG
	writeBinaryImage(artifactImage, "_artifact_image_labelled", 0);
	#endif
	
	std::list< ImageType::IndexType > objectBorders;
	std::list< ImageType::IndexType > objects;
	
	// 		float * objectBorderMean = new float [nrOfObjects];
	std::list< unsigned int > borderPixelMean;
	
	ConstObjectIndexIteratorType objectIter(labeledImage, labeledImage->GetLargestPossibleRegion());
	IteratorType2 inputIter(fgRemovedImage, fgRemovedImage->GetLargestPossibleRegion());
	for(objectIter.GoToBegin(), inputIter.GoToBegin(); !objectIter.IsAtEnd() && !inputIter.IsAtEnd(); ++objectIter, ++inputIter)
	{
		if(objectIter.Get())
		{
			objects.push_back(objectIter.GetIndex());
			inputIter.Set(0);
		}
	}
	fgRemovedImage->Update();
	#ifdef DEBUG
	writeBinaryImage(fgRemovedImage, "_fg_removed_image", 0);
	#endif
	unsigned long remainingPx = objects.size();
	while(!objects.empty())
	{
		calculateBorderPixels( fgRemovedImage, objectBorders, objects, borderPixelMean );
		if(objects.size() == remainingPx)		// in case of 0 areas surrounding foreground objects -> leads to infinite loop otherwise
			break;
		remainingPx = objects.size();
	}
	#ifdef PIPELINE_DOC
	writeBinaryImage(fgRemovedImage, "_barrels_step8_fill_holes", 0);
	#endif
	
	objectBorders.clear();
	objects.clear();
	borderPixelMean.clear();
};

void Segmentation::enlargeHoles( ImageType::Pointer workingImage )
{
	SegNeighborhoodIteratorType::RadiusType borderradius;
	borderradius.Fill(1);
	SegNeighborhoodIteratorType borderIter( borderradius, workingImage, workingImage->GetLargestPossibleRegion() );
	
	IndexVectorType borderVector;
	
	//ObjectImageType::OffsetType bordercenter = {{0,0,0}};
	//ObjectImageType::OffsetType bordertop = {{0,0,1}};
	//ObjectImageType::OffsetType borderbottom = {{0,0,-1}};
	ObjectImageType::OffsetType borderright = {{0,1,0}};
	ObjectImageType::OffsetType borderleft = {{0,-1,0}};
	ObjectImageType::OffsetType borderfront = {{1,0,0}};
	ObjectImageType::OffsetType borderback = {{-1,0,0}};
	
	for( borderIter.GoToBegin() ; !borderIter.IsAtEnd() ; ++borderIter )
	{
		if(borderIter.GetCenterPixel())
		{
			bool noBorder = true;
			
			borderIter.GetPixel( borderright, noBorder );
			
			if(noBorder && !borderIter.GetPixel(borderright))
			{
				borderVector.push_back(borderIter.GetIndex(borderright));
			}
			
			
			borderIter.GetPixel( borderleft, noBorder );
			
			if(noBorder && !borderIter.GetPixel(borderleft))
			{
				borderVector.push_back(borderIter.GetIndex(borderleft));
			}
			
			
			borderIter.GetPixel( borderfront, noBorder );
			
			if(noBorder && !borderIter.GetPixel(borderfront))
			{
				borderVector.push_back(borderIter.GetIndex(borderfront));
			}
			
			borderIter.GetPixel( borderback, noBorder );
			
			if(noBorder && !borderIter.GetPixel(borderback))
			{
				borderVector.push_back(borderIter.GetIndex(borderback));
			}
		}
	}
	
	IndexIteratorType indexIter( workingImage, workingImage->GetLargestPossibleRegion() );
	
	for( int i = 0 ; i < borderVector.size(); ++i )
	{
		indexIter.SetIndex(borderVector[i]);
		indexIter.Set(255);
	}
	
	borderVector.clear();
};

void Segmentation::calculateBorderPixels( ImageType::Pointer workingImage, std::list< ImageType::IndexType >& objectBorders, std::list< ImageType::IndexType >& objects, std::list< unsigned int >& borderMean )
{
	#ifdef DEBUG
	std::flush(DebugLog << "in calculateBorderPixels():" << std::endl);
	std::flush(DebugLog << "objects.size() = " << objects.size() << std::endl);
	#endif
	
	ShapedNeighborhoodIteratorType::RadiusType radius;
	radius.Fill(1);
	ShapedNeighborhoodIteratorType neighborIter( radius, workingImage, workingImage->GetLargestPossibleRegion() );
	for(int ii = 0; ii < 8; ++ii)
		neighborIter.ActivateOffset(neighborhood8[ii]);
	ShapedNeighborhoodIteratorType::Iterator iter;
	
	std::list< ImageType::IndexType >::iterator objectIter;
	objectIter = objects.begin();
	while(objectIter != objects.end())
	{
		neighborIter.SetLocation(*objectIter);

		ImageType::IndexType debugIndex;
		debugIndex[0] = 12, debugIndex[1] = 0, debugIndex[2] = 0;
		
		float mean = 0;
		unsigned int borderpxlcount = 0;
		bool border = false;
		
		for(iter = neighborIter.Begin(); !iter.IsAtEnd(); ++iter)
		{
			if(iter.Get())
			{
				if(!border)
				{
					border = true;
				}
				mean += iter.Get();
				++borderpxlcount;
			}
		}
		
		if(border && borderpxlcount > 1)
		{
			objectBorders.push_back(*objectIter);
			mean /= (float)borderpxlcount;
			borderMean.push_back((unsigned int)(mean + 0.5));
// 			borderMean.push_back(borderpxlcount);
			objectIter = objects.erase(objectIter);
			
		}
		else
			++objectIter;
	}
	
	#ifdef DEBUG
	std::flush(DebugLog << "objectBorders.size() = " << objectBorders.size() << std::endl);
	#endif
	
	IndexIteratorType setIter(workingImage, workingImage->GetLargestPossibleRegion());
	std::list< ImageType::IndexType >::iterator indexIt;
	std::list< unsigned int >::iterator valueIt;
	for(indexIt = objectBorders.begin(), valueIt = borderMean.begin(); indexIt != objectBorders.end() && valueIt != borderMean.end(); ++indexIt, ++valueIt)\
	{
		setIter.SetIndex(*indexIt);
		setIter.Set(*valueIt);
// 		setIter.Set(255);
	}
	workingImage->Update();
	objectBorders.clear();
	borderMean.clear();
	
// 	#ifdef DEBUG
// 	std::cout << objectpxlcount << " object pixels, " << borderpxlcount << " border pixels." << std::endl;
// 	#endif
};

void Segmentation::fillBinaryHoles()
{
	ImageType::Pointer invertedImage = ImageType::New();
	invertedImage = invertImage(inputImage);
	regionGrowingObjectLabeling(invertedImage, 1, 0);
	ConstObjectIteratorType holeIter(labeledImage, labeledImage->GetLargestPossibleRegion());
	IteratorType2 fillIter(inputImage, inputImage->GetLargestPossibleRegion());
	for(holeIter.GoToBegin(), fillIter.GoToBegin(); !holeIter.IsAtEnd() && !fillIter.IsAtEnd(); ++holeIter, ++fillIter)
		if(holeIter.Get() > 1)
			fillIter.Set(255);
	
	inputImage->Update();
};

void Segmentation::fillBinaryHoles(ImageType::Pointer workingImage)
{
	ImageType::Pointer invertedImage = ImageType::New();
	invertedImage = invertImage(workingImage);
	
	ConnectedFilterType::Pointer labeller = ConnectedFilterType::New();
	labeller->SetInput(invertedImage);
	labeller->Update();
	
	RelabelType::Pointer relabeller = RelabelType::New();
	relabeller->SetInput(labeller->GetOutput());
	relabeller->SetMinimumObjectSize( 1 );
// 	relabeller->Update();
	
	ObjectImageType::Pointer thisLabeledImage = ObjectImageType::New();
	thisLabeledImage = relabeller->GetOutput();
	thisLabeledImage->Update();
	unsigned long labelThreshold = 0;
	unsigned long nrOfLabeledObjects = relabeller->GetNumberOfObjects();
	std::vector< unsigned long> labeledObjectSize = relabeller->GetSizeOfObjectsInPixels();
	for(int ii = 0; ii < nrOfLabeledObjects; ++ii)
		if(labeledObjectSize[ii] < 5E03)
		{
			labelThreshold = ii;
			break;
		}
	
	invertedImage = invertImage(invertedImage);
	ConstObjectIteratorType holeIter(thisLabeledImage, thisLabeledImage->GetLargestPossibleRegion());
	IteratorType2 fillIter(invertedImage, invertedImage->GetLargestPossibleRegion());
	for(holeIter.GoToBegin(), fillIter.GoToBegin(); !holeIter.IsAtEnd() && !fillIter.IsAtEnd(); ++holeIter, ++fillIter)
		if(holeIter.Get() > labelThreshold)
			fillIter.Set(255);
	
	copy(workingImage, invertedImage);
};

//interface for Section class
std::vector< std::vector< Contour * > > Segmentation::getRegularBarrelContours()
{
	ImageType::RegionType planeRegion;
	ImageType::RegionType planeImageRegion;
	ImageType::IndexType planeImageIndex;
	ImageType::IndexType planeIndex;
	ImageType::SizeType planeSize;
	planeImageIndex = originalImage->GetLargestPossibleRegion().GetIndex();
	planeIndex = originalImage->GetLargestPossibleRegion().GetIndex();
	planeSize = originalImage->GetLargestPossibleRegion().GetSize();
	planeSize[2] = 1;
	planeImageRegion.SetIndex(planeImageIndex);
	planeImageRegion.SetSize(planeSize);
	inputImage->SetRegions(planeImageRegion);
	inputImage->Allocate();
	preprocImage->SetRegions(planeImageRegion);
	preprocImage->Allocate();
	bool optimizedThreshold = 0;
	std::flush(std::cout << "Region growing from barrel markers...");
	for(int z = 0; z < (stop - start + 1); ++z)
	{
		planeIndex[2] = z;
		planeRegion.SetIndex(planeIndex);
		planeRegion.SetSize(planeSize);
		copy(inputImage, segmentedImageStack, planeRegion);
// 		#ifdef PIPELINE_DOC
// 		writeInputImage("_barrels_step11_0_input_image.tif");
// 		#endif
		copy(preprocImage, preprocImageStack, planeRegion);
		maximum_region = inputImage->GetLargestPossibleRegion();
		splitWrongBarrelObjects(voronoiMap,  barrelZContours[z].size(), optimizedThreshold);
// 		#ifdef PIPELINE_DOC
// 		writeInputImage("_barrels_step12_0_input_image.tif");
// 		#endif
		mergeTrueBarrelObjects(distanceMap, voronoiMap,  barrelZContours[z].size(), barrelZContours[z], z, optimizedThreshold);
	}
	std::cout << "done!" << std::endl;
	
	return barrelZContours;
};

//interface for Section class
std::vector< std::vector< Contour * > > Segmentation::getOptimizedBarrelContours()
{
	ImageType::RegionType planeRegion;
	ImageType::RegionType planeImageRegion;
	ImageType::IndexType planeImageIndex;
	ImageType::IndexType planeIndex;
	ImageType::SizeType planeSize;
	planeImageIndex = preprocImageStack->GetLargestPossibleRegion().GetIndex();
	planeIndex = preprocImageStack->GetLargestPossibleRegion().GetIndex();
	planeSize = preprocImageStack->GetLargestPossibleRegion().GetSize();
	planeSize[2] = 1;
	planeImageRegion.SetIndex(planeImageIndex);
	planeImageRegion.SetSize(planeSize);
	inputImage->SetRegions(planeImageRegion);
	inputImage->Allocate();
	preprocImage->SetRegions(planeImageRegion);
	preprocImage->Allocate();
	segmentedImage = ImageType::New();
	segmentedImage->SetRegions(planeImageRegion);
	segmentedImage->Allocate();
	std::flush(std::cout << "Optimizer region growing from barrel markers...");
	#ifdef DEBUG
	DebugLog << "Optimizer region growing from barrel markers" << std::endl;
	DebugLog << std::endl;
	#endif
	for(int z = 0; z < (stop - start + 1); ++z)
	{
		planeIndex[2] = z;
		planeRegion.SetIndex(planeIndex);
		planeRegion.SetSize(planeSize);
		copy(preprocImage, preprocImageStack, planeRegion);
		maximum_region = inputImage->GetLargestPossibleRegion();
		optimizedVoronoiRegionThreshold(barrelZContours[z].size(), barrelZContours[z], z);
		#ifdef PIPELINE_DOC
		writeBinaryImage(segmentedImage, "_barrels_step14_opt_voronoi_region_growing", 0);
		#endif
		copy(inputImage, segmentedImage);
		bool optimizedThreshold = 1;
		splitWrongBarrelObjects(voronoiMap, barrelZContours[z].size(), optimizedThreshold);
		mergeTrueBarrelObjects(distanceMap, voronoiMap, barrelZContours[z].size(), barrelZContours[z], z, optimizedThreshold);
	}
	std::cout << "done!" << std::endl;
	
	#ifdef DEBUG
	DebugLog << "Optimizer region growing from barrel markers completed!" << std::endl;
	DebugLog << std::endl;
	#endif
	
	return barrelZContours;
};

//interface for Section class
void Segmentation::writeBarrelContours(/*std::vector< Contour * > barrelContours*/)
{
// 	if(!(amira_bvp_graph->vertice_list.size()))
// 	{
// 		std::flush(std::cout << "converting bvp to spatial graph..." << std::endl);
// 		BloodVesselPatternToAmiraSpatialGraph();
// 	}
	#ifdef DEBUG
	DebugLog << "barrel vector passed into Segmentation object" << std::endl;
	DebugLog << "starting Segmentation::writeBarrelContours()" << std::endl;
	DebugLog << "barrelZContours.size() = " << barrelZContours.size() << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
	{
		DebugLog << "\tbarrelZContours[z].size() = " << barrelZContours[z].size() << std::endl;
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
		{
			DebugLog << "\t\tbarrelZContours[" << z << "][" << ii << "] edge list size = " << barrelZContours[z][ii]->edgeListPointer()->size() << std::endl;
			DebugLog << "\t\tbarrelZContours[" << z << "][" << ii << "] valid = " << barrelZContours[z][ii]->getValid() << std::endl;
			DebugLog << "\t\tbarrelZContours[" << z << "][" << ii << "] barrel ID = " << barrelZContours[z][ii]->getBarrelID() << std::endl;
		}
	}
	DebugLog << std::endl;
	#endif
	for(int z = 0; z < barrelZContours.size(); ++z)
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
			if(barrelZContours[z][ii]->getValid())
			{
	//			std::flush(std::cout << "Starting contour smoothing... contour size = " << barrelZContours[z]->getEdgeList()->size() << std::endl);
				if(!barrelZContours[z][ii]->edgeListPointer()->size())
				{
					fillEmptyBarrelContour(ii, z);
		// 			std::flush(std::cout << "replaced contour size = " << barrelZContours[z]->getEdgeList()->size() << std::endl);
				}
				if(barrelZContours[z][ii]->edgeListPointer()->size())
				{
					contourSmoothing(*(barrelZContours[z][ii]->edgeListPointer()), 5);
		// 			std::flush(std::cout << "Starting contour sampling... contour size = ");
					contourSampling(*(barrelZContours[z][ii]->edgeListPointer()), 10);
		// 			std::flush(std::cout << barrelZContours[z]->getEdgeList()->size() << std::endl);
				}
			}
	#ifdef DEBUG
	DebugLog << "barrel vector after contour filling, smoothing and sampling" << std::endl;
	DebugLog << "barrelZContours.size() = " << barrelZContours.size() << std::endl;
	for(int z = 0; z < barrelZContours.size(); ++z)
	{
		DebugLog << "\tbarrelZContours[" << z << "].size() = " << barrelZContours[z].size() << std::endl;
		for(int ii = 0; ii < barrelZContours[z].size(); ++ii)
		{
			DebugLog << "\t\tbarrelZContours[" << z << "][" << ii << "] edge list size = " << barrelZContours[z][ii]->edgeListPointer()->size() << std::endl;
			DebugLog << "\t\tbarrelZContours[" << z << "][" << ii << "] valid = " << barrelZContours[z][ii]->getValid() << std::endl;
			DebugLog << "\t\tbarrelZContours[" << z << "][" << ii << "] barrel ID = " << barrelZContours[z][ii]->getBarrelID() << std::endl;
		}
	}
	DebugLog << std::endl;
	#endif
// // 	std::flush(std::cout << "getting amira contours from barrel contours" << std::endl);
// 	// 	std::list<std::list<std::vector<float> > > amira_edges = GetAmiraContours(pxcoord_edges);
// 	std::list<std::list<std::vector<float> > > amira_edges = GetAmiraContours(barrelZContours, barrelZContours[0].size());
// 	
// 	std::list<std::vector<float> > vertex_list = SetContourVertices(amira_edges);
// 	
// 	AmiraContourGraphType *tmp_contour = new AmiraContourGraphType(vertex_list, amira_edges);
// 	amira_contour_graph = tmp_contour;
// 	
// 	#ifdef DEBUG
// 	DebugLog << "amira_edges" << std::endl;
// 	DebugLog << "amira_edges.size() = " << amira_edges.size() << std::endl;
// 	std::list<std::list<std::vector<float> > >::const_iterator amira_edges_it;
// 	int cnt = 0;
// 	for(amira_edges_it = amira_edges.begin(); amira_edges_it != amira_edges.end(); ++amira_edges_it, ++cnt)
// 		DebugLog << "\tamira_edges[" << cnt << "] edge list size = " << amira_edges_it->size() << std::endl;
// 	DebugLog << std::endl;
// 	#endif
	
	std::flush(std::cout << "writing attribute files" << std::endl);
	#ifdef DEBUG
	DebugLog << "writing attribute files" << std::endl;
	DebugLog << std::endl;
	#endif
	for(int ii = 0; ii < barrelZContours[0].size(); ++ii)
		if(barrelZContours[0][ii]->getValid())
		{
			std::string output(this->outputFilename);
			output += "_barrelID_%02d";
			output += "_attributes.csv";
			char * outChar = new char[256];
			sprintf(outChar, output.c_str(), barrelZContours[0][ii]->getBarrelID());
			// 	std::flush(std::cout << output.c_str() << std::endl);
			std::ofstream AttributeFile(outChar);
			AttributeFile << "\tPlane\tvoronoi ID\tcirc1 SNR\tcirc2 SNR\tcirc3 SNR\treg. threshold\tradius\treg. SNR\toptimize flag\tbarrelArea/VoronoiArea" << std::endl;
			// 	AttributeFile << "Plane\tregion mean" << std::endl;
			std::vector< float >::const_iterator attrIt;
			for(int z = 0; z < barrelZContours.size(); ++z)
			{
				// 		AttributeFile << z;
				int attrCount = 0;
				for(attrIt = barrelZContours[z][ii]->readAttributes(); attrIt != barrelZContours[z][ii]->attributesEnd(); ++attrIt, ++attrCount)
				{
					AttributeFile << "\t" << *attrIt;
					if(attrCount == 7)
						AttributeFile << "\t" << barrelZContours[z][ii]->getOptimizeFlag();
				}
				AttributeFile << std::endl;
			}
			AttributeFile.close();
		}
	
	std::flush(std::cout << "writing spatial graph file" << std::endl);
	#ifdef DEBUG
	DebugLog << "writing spatial graph file" << std::endl;
	DebugLog << std::endl;
	#endif
	WriteBarrelSpatialGraphFile(this->outputFilename);
// 	WriteSpatialGraphFile(this->outputFilename);
	
// 	CalcToImageRescaleFilterType::Pointer cToIFilter = CalcToImageRescaleFilterType::New();
// 	cToIFilter->SetOutputMinimum(1);
// 	cToIFilter->SetOutputMaximum(nrOfBarrelMarkerPoints);
// 	cToIFilter->SetInput(voronoiMap);
// 	cToIFilter->Update();
// 	writeBinaryImage(cToIFilter->GetOutput(), "_voronoi_map", 0);
// 	writeBinaryImage(preprocImage, "_preproc", 0);
};

//interface for Section class
void Segmentation::writeBarrelContours(std::vector< Contour * > barrelContours)
{
	#ifdef DEBUG
	DebugLog << "barrel vector passed into Segmentation object" << std::endl;
	DebugLog << "starting Segmentation::writeBarrelContours()" << std::endl;
	DebugLog << "barrelContours.size() = " << barrelContours.size() << std::endl;
	for(int z = 0; z < barrelContours.size(); ++z)
	{
		DebugLog << "\t\tbarrelContours[" << z << "] edge list size = " << barrelContours[z]->edgeListPointer()->size() << std::endl;
		DebugLog << "\t\tbarrelContours[" << z << "] valid = " << barrelContours[z]->getValid() << std::endl;
		DebugLog << "\t\tbarrelContours[" << z << "] barrel ID = " << barrelContours[z]->getBarrelID() << std::endl;
	}
	DebugLog << std::endl;
	#endif
	for(int z = 0; z < barrelContours.size(); ++z)
		if(barrelContours[z]->getValid())
		{
//			std::flush(std::cout << "Starting contour smoothing... contour size = " << barrelZContours[z]->getEdgeList()->size() << std::endl);
			if(!barrelContours[z]->edgeListPointer()->size())
			{
// 				fillEmptyBarrelContour(ii, z);
	// 			std::flush(std::cout << "replaced contour size = " << barrelZContours[z]->getEdgeList()->size() << std::endl);
			}
			if(barrelContours[z]->edgeListPointer()->size())
			{
				contourSmoothing(*(barrelContours[z]->edgeListPointer()), 5);
	// 			std::flush(std::cout << "Starting contour sampling... contour size = ");
				contourSampling(*(barrelContours[z]->edgeListPointer()), 10);
	// 			std::flush(std::cout << barrelZContours[z]->getEdgeList()->size() << std::endl);
			}
		}
	#ifdef DEBUG
	DebugLog << "barrel vector after contour filling, smoothing and sampling" << std::endl;
	DebugLog << "barrelContours.size() = " << barrelContours.size() << std::endl;
	for(int z = 0; z < barrelContours.size(); ++z)
	{
		DebugLog << "\t\tbarrelContours[" << z << "] edge list size = " << barrelContours[z]->edgeListPointer()->size() << std::endl;
		DebugLog << "\t\tbarrelContours[" << z << "] valid = " << barrelContours[z]->getValid() << std::endl;
		DebugLog << "\t\tbarrelContours[" << z << "] barrel ID = " << barrelContours[z]->getBarrelID() << std::endl;
	}
	DebugLog << std::endl;
	#endif
	
	std::flush(std::cout << "writing attribute files" << std::endl);
	#ifdef DEBUG
	DebugLog << "writing attribute files" << std::endl;
	DebugLog << std::endl;
	#endif
	for(int ii = 0; ii < barrelZContours[0].size(); ++ii)
		if(barrelZContours[0][ii]->getValid())
		{
			std::string output(this->outputFilename);
			output += "_barrelID_%02d";
			output += "_attributes.csv";
			char * outChar = new char[256];
			sprintf(outChar, output.c_str(), barrelZContours[0][ii]->getBarrelID());
			// 	std::flush(std::cout << output.c_str() << std::endl);
			std::ofstream AttributeFile(outChar);
			AttributeFile << "\tPlane\tvoronoi ID\tcirc1 SNR\tcirc2 SNR\tcirc3 SNR\treg. threshold\tradius\treg. SNR\toptimize flag\tbarrelArea/VoronoiArea" << std::endl;
			// 	AttributeFile << "Plane\tregion mean" << std::endl;
			std::vector< float >::const_iterator attrIt;
			for(int z = 0; z < barrelZContours.size(); ++z)
			{
				// 		AttributeFile << z;
				int attrCount = 0;
				for(attrIt = barrelZContours[z][ii]->readAttributes(); attrIt != barrelZContours[z][ii]->attributesEnd(); ++attrIt, ++attrCount)
				{
					AttributeFile << "\t" << *attrIt;
					if(attrCount == 7)
						AttributeFile << "\t" << barrelZContours[z][ii]->getOptimizeFlag();
				}
				AttributeFile << std::endl;
			}
			AttributeFile.close();
		}
	
	std::flush(std::cout << "writing spatial graph file" << std::endl);
	#ifdef DEBUG
	DebugLog << "writing spatial graph file" << std::endl;
	DebugLog << std::endl;
	#endif
	WriteBarrelSpatialGraphFile(this->outputFilename);
};

void Segmentation::writeVoronoiDiagram()
{
	CalcToImageRescaleFilterType::Pointer cToIFilter = CalcToImageRescaleFilterType::New();
	cToIFilter->SetOutputMinimum(0);
	cToIFilter->SetOutputMaximum(nrOfBarrelMarkerPoints);
	cToIFilter->SetInput(voronoiMap);
	cToIFilter->Update();
	writeBinaryImage(cToIFilter->GetOutput(), "_voronoi_map", 0);
};





















