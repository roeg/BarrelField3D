/****************************************************************************/
/*                                                                          */
/* File:      basics.cpp                                                    */
/*                                                                          */
/* Purpose:   class for basic functions and common variables                */
/*            especially writing and reading functions                      */
/*            extended functionality (classes Contour & BarrelMarker)       */
/*            added for barrel segmentation pipeline                        */
/*                                                                          */
/* Authors:   Marcel Oberlaender                                            */
/*            Max-Planck-Institute of Neurobiology                          */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/*            Stefan Reissl                                                 */
/*            Max-Planck-Institute of Neurobiology                          */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/*            Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   22.12.2010                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#include "basics.h"

Basics::Basics(void)
{
	inputFilename	= 0;
	outputFilename	= 0;
	inputImage=0;
	look_up_table = CreateLookUpTable();
	createNeighborhood();
	radius1[0] = 1;
	radius1[1] = 1;
	radius1[2] = 1;
}

Basics::~Basics(void)
{
	
}

void Basics::writeInputImage(const char * name)
{
	ImageType::Pointer image = ImageType::New();
	image = inputImage;
	
	Writer2DType::Pointer writer = Writer2DType::New();
	
	writer->SetInput( inputImage );
	
	NameGeneratorType::Pointer writerNameGen = NameGeneratorType::New();
	
	std::string format = outputFilename;
	format += name;
	
	writerNameGen->SetSeriesFormat( format.c_str() );
	
	ImageType::RegionType region = inputImage->GetLargestPossibleRegion();
	ImageType::IndexType start = region.GetIndex();
	ImageType::SizeType size = region.GetSize();
	
	const unsigned int firstSlice = start[2];
	const unsigned int lastSlice  = start[2] + size[2] -1;
	
	writerNameGen->SetStartIndex( firstSlice );
	writerNameGen->SetEndIndex( lastSlice );
	writerNameGen->SetIncrementIndex( 1 );
	
	writer->SetFileNames( writerNameGen->GetFileNames() );
	
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err )
	{
		std::cerr << "WriterExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}
};

void Basics::writeInput2DImage(const char * name)
{
// 	ImageType::Pointer image = ImageType::New();
// 	image = inputImage;
	
	Single2DWriterType::Pointer writer = Single2DWriterType::New();
	
	writer->SetInput( input2DImage );
	
	NameGeneratorType::Pointer writerNameGen = NameGeneratorType::New();
	
	std::string format = outputFilename;
	format += name;
	
// 	writerNameGen->SetSeriesFormat( format.c_str() );
	
// 	ImageType::RegionType region = inputImage->GetLargestPossibleRegion();
// 	ImageType::IndexType start = region.GetIndex();
// 	ImageType::SizeType size = region.GetSize();
// 	
// 	const unsigned int firstSlice = start[2];
// 	const unsigned int lastSlice  = start[2] + size[2] -1;
// 	
// 	writerNameGen->SetStartIndex( firstSlice );
// 	writerNameGen->SetEndIndex( lastSlice );
// 	writerNameGen->SetIncrementIndex( 1 );
/*	
	writer->SetFileName( writerNameGen->GetFileNames() );*/

	writer->SetFileName( format );
	
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err )
	{
		std::cerr << "WriterExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}
};

// ImageType::Pointer Basics::LoadCluster(CellCluster * cluster, const char * path)
// {
// 	ImageType::Pointer retImage= ImageType::New();
// 	std::stringstream iStream;
// 	iStream << cluster->GetID();
// 	unsigned char color;
// 	
// 	std::string filename = path;
// 	filename += "Cluster_NR_";
// 	filename += iStream.str();
// 	filename += ".sr";
// 	
// 	const char* in = filename.c_str();
// 	std::ifstream filereader(in);
// 	
// 	retImage->SetRegions(cluster->GetPureRegions());
// 	retImage->Allocate();
// 	retImage->FillBuffer(0);
// 	retImage->Update();
// 	
// 	IteratorType2 iter(retImage, retImage->GetLargestPossibleRegion());
// 	
// 	if(!filereader.fail())
// 	{
// 		for(iter.GoToBegin(); !iter.IsAtEnd() && !filereader.eof(); ++iter)
// 		{
// 			try
// 			{
// 				filereader >> color;
// 				
// 				if(color!=0)
// 					iter.Set(color);		
// 			}
// 			catch(...)
// 			{
// 				std::cout<<"something is wrong \n";
// 			}
// 			
// 		}
// 		
// 		filereader.close();
// 	}
// 	else
// 		std::cout <<"Filereader Error \n";
// 	
// 	return retImage;
// }

// ImageType::Pointer Basics::LoadClusterWithBox(CellCluster * cluster)
// {
// 	ImageType::Pointer retImage= ImageType::New();
// 	std::stringstream iStream;
// 	iStream << cluster->GetID();
// 	unsigned char color;
// 	ImageType::IndexType pos;
// 	
// 	std::string filename = outputFilename;
// 	filename += "Cluster_NR_";
// 	filename += iStream.str();
// 	filename += ".sr";
// 	
// 	const char* in = filename.c_str();
// 	std::ifstream filereader(in);
// 	
// 	retImage->SetRegions(cluster->GetRegionsWithBox());
// 	retImage->Allocate();
// 	retImage->FillBuffer(0);
// 	retImage->Update();
// 	
// 	IteratorType2 iter(retImage, cluster->GetPureRegions());
// 	
// 	//	int counter=0;
// 	
// 	if(!filereader.fail())
// 	{
// 		for(iter.GoToBegin(); !iter.IsAtEnd() && !filereader.eof(); ++iter)
// 		{
// 			try
// 			{
// 				pos=iter.GetIndex();
// 				
// 				pos[0]+=1;
// 				pos[1]+=1;
// 				pos[2]+=1;
// 				
// 				filereader >> color;
// 				
// 				if(color!=0)
// 					color=255;
// 				else
// 					color=0;
// 				
// 				retImage->SetPixel(pos,color);
// 			}
// 			catch(...)
// 			{
// 				std::cout<<"something is wrong \n";
// 			}
// 			
// 		}
// 		
// 		filereader.close();
// 	}
// 	
// 	return retImage;
// }


// void Basics::writeCellClusterList(CellClusterListType clusters, const char * file)
// {
// 	std::string filename = outputFilename;
// 	filename += file;//"CellClusterList.csv";
// 	const char* out = filename.c_str();
// 	std::ofstream filewriter( out, std::ios::in | std::ios::out | std::ios::trunc );
// 	
// 	CellClusterListType::iterator listIt;
// 	
// 	
// 	if( !filewriter.fail() )
// 	{	
// 		for(listIt = clusters.begin() ; listIt != clusters.end(); ++listIt)
// 		{
// 			filewriter << (float)listIt->GetID()  << " " << (float)listIt->GetSize() << " "
// 			<< (float)listIt->GetXMin() << " " << (float)listIt->GetXMax() << " "
// 			<< (float)listIt->GetYMin() << " " << (float)listIt->GetYMax() << " "
// 			<< (float)listIt->GetZMin() << " " << (float)listIt->GetZMax() << " "
// 			<< (float)listIt->getLandmark()[0] << " " << (float)listIt->getLandmark()[1] <<" " << (float)listIt->getLandmark()[2] << " "
// 			<< (float)listIt->GetStartPixel()[0] << " " << (float)listIt->GetStartPixel()[1] <<" " << (float)listIt->GetStartPixel()[2] << " \n";
// 		}
// 		filewriter.close();
// 	}
// }

bool Basics::IsNotAtImageBorder(float x, float y, float z)
{
	x=y=z=1;
	/*if(x<1 ||y<1 ||z<1)
	return false;
	
	if(x>(maximum_region.GetSize(0)-1) || y>(maximum_region.GetSize(1)-1) || z>(maximum_region.GetSize(2)-1))
	return false;*/
	
	return true;
}

void Basics::readImage( int start, int stop, const char * inputFilename )
{
	SeriesReaderType::Pointer seriesReader= SeriesReaderType::New();
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	
	nameGenerator->SetSeriesFormat(inputFilename);
	nameGenerator->SetStartIndex(start);
	nameGenerator->SetEndIndex(stop);
	nameGenerator->SetIncrementIndex( 1 );
	
	seriesReader->SetFileNames( nameGenerator->GetFileNames() );
	
	try
	{
		seriesReader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}
	
	inputImage=seriesReader->GetOutput();
	
	boundingBox = new long[6];
	boundingBox[0] = inputImage->GetLargestPossibleRegion().GetIndex()[0];
	boundingBox[1] = boundingBox[0] + inputImage->GetLargestPossibleRegion().GetSize()[0] - 1;
	boundingBox[2] = inputImage->GetLargestPossibleRegion().GetIndex()[1];
	boundingBox[3] = boundingBox[2] + inputImage->GetLargestPossibleRegion().GetSize()[1] - 1;
	boundingBox[4] = inputImage->GetLargestPossibleRegion().GetIndex()[2];
	boundingBox[5] = boundingBox[4] + inputImage->GetLargestPossibleRegion().GetSize()[2] - 1;
};

void Basics::read2DImage( const char * inputFilename )
{
	Reader2DType::Pointer planeReader= Reader2DType::New();
// 	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	
// 	nameGenerator->SetSeriesFormat(inputFilename);
// 	nameGenerator->SetStartIndex(start);
// 	nameGenerator->SetEndIndex(stop);
// 	nameGenerator->SetIncrementIndex( 1 );
	
// 	planeReader->SetFileName( nameGenerator->GetFileNames() );
	planeReader->SetFileName( inputFilename );
	
	try
	{
		planeReader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}
	
	input2DImage=planeReader->GetOutput();
	
	boundingBox = new long[4];
	boundingBox[0] = inputImage->GetLargestPossibleRegion().GetIndex()[0];
	boundingBox[1] = boundingBox[0] + input2DImage->GetLargestPossibleRegion().GetSize()[0] - 1;
	boundingBox[2] = inputImage->GetLargestPossibleRegion().GetIndex()[1];
	boundingBox[3] = boundingBox[2] + input2DImage->GetLargestPossibleRegion().GetSize()[1] - 1;
};

// bool Basics::isClusterValid(CellCluster *cluster, ImageType::RegionType region)
// {
// 	bool isValid = true;
// 	bool isProbValid = true;
// 	
// 	//cluster parameters
// 	//unsigned int id                       = cluster->GetID();
// 	unsigned long size                    = cluster->GetSize();
// 	ImageType::IndexType start            = cluster->GetStartPixel();
// 	unsigned int x_max                    = cluster->GetXMax();
// 	unsigned int x_min                    = cluster->GetXMin();
// 	unsigned int y_max                    = cluster->GetYMax();
// 	unsigned int y_min                    = cluster->GetYMin();
// 	unsigned int z_max                    = cluster->GetZMax();
// 	unsigned int z_min                    = cluster->GetZMin();
// 	std::vector<float>  probabilites      = cluster->getClusterProbabilities();
// 	std::vector<float> landmark           = cluster->getLandmark();
// 	//bool borderCell                       = cluster->isAtBorder();
// 	//bool isObsolete                       = cluster->IsObsolete();
// 	ImageType::RegionType boundingBox = cluster->GetPureRegions();
// 	
// 	//cluster constraints
// 	unsigned long imageSizeX  = region.GetSize(0);
// 	unsigned long imageSizeY  = region.GetSize(1);
// 	unsigned long imageSizeZ  = region.GetSize(2);
// 	unsigned long maxVolume   = imageSizeX * imageSizeY * imageSizeZ;
// 	
// 	if(size == 0 || size > maxVolume)
// 	{
// 		std::cout<< "Cluster not valid!!! False Size: " << size << std::endl;
// 		isValid = false;
// 	}
// 	if(start[0]>=imageSizeX || start[1] >= imageSizeY || start[2] >= imageSizeZ)
// 	{
// 		std::cout<< "Cluster not valid!!! False Start Pixel: " << start << std::endl;
// 		isValid = false;
// 	}
// 	if(x_min > x_max || y_min > y_max || z_min > z_max || x_min >= imageSizeX || x_max >= imageSizeX || y_min >= imageSizeY || y_max >= imageSizeY || z_min >= imageSizeZ || z_max >= imageSizeZ)
// 	{
// 		std::cout<< "Cluster not valid!!! False Max Values: " << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << std::endl;
// 		isValid = false;
// 	}
// 	if(probabilites.size()>0)
// 	{
// 		isProbValid = false;
// 		for(int i=0; i<probabilites.size();i++)
// 		{
// 			if(probabilites[i]!=0)
// 			{
// 				isProbValid = true;
// 				break;
// 			}
// 		}
// 		
// 		if(isProbValid==false)
// 		{
// 			std::cout<< "Cluster not valid in morphCellCount!!! Prob set to 1!" << std::endl;
// 			
// 			landmark[3]=1;
// 			probabilites[0]=1;
// 			
// 			cluster->setLandmark(landmark);
// 			cluster->setClusterProbabilities(probabilites);
// 			
// 			isProbValid=true;
// 		}
// 		
// 		if(isProbValid==false)
// 			isValid=false;
// 	}
// 	if(landmark.size()>0)
// 	{	
// 		if(landmark[0]<0 || landmark[0]>= imageSizeX || landmark[1]<0 || landmark[1]>= imageSizeY || landmark[2]<0 || landmark[2]>= imageSizeZ || landmark[3]<0 || landmark.size()<3)
// 		{
// 			std::cout<< "Cluster not valid!!! False Landmark Pos: " << landmark[0] << " " << landmark[1] << " " << landmark[2] << " " << landmark[3] << std::endl;
// 			isValid = false;
// 		}
// 	}
// 	if(boundingBox.GetSize(0)>imageSizeX || boundingBox.GetSize(1)>imageSizeY || boundingBox.GetSize(2) >imageSizeZ || boundingBox.GetSize(0)<1 || boundingBox.GetSize(1)<1 || boundingBox.GetSize(2)<1)
// 	{
// 		std::cout<< "Cluster not valid!!! False Bounding Box: " << boundingBox.GetSize() << std::endl;
// 		isValid = false;
// 	}
// 	return isValid;
// };

// CellClusterListType Basics::readCellClusterList( const char * _filename )
// {
// 	CellClusterListType retList;
// 	
// 	std::ifstream filereader(_filename);
// 	
// 	std::string line,data;
// 	std::string::size_type pos=0;
// 	float number[14];
// 	
// 	if( !filereader.fail())
// 	{
// 		while (std::getline(filereader, line))
// 		{
// 			unsigned int counter=0;
// 			CellCluster cluster;
// 			std::vector<float> landmark;
// 			ImageType::IndexType index;
// 			
// 			while(line.find(' ') != std::string::npos)
// 			{
// 				pos=line.find(' ');
// 				data=line.substr(0,pos);
// 				number[counter]=atof(line.c_str());
// 				line.erase(0,pos+1);
// 				counter++;
// 			}
// 			
// 			cluster.SetID(  (unsigned int) 	number[0]);
// 			cluster.SetSize((unsigned int) 	number[1]);
// 			
// 			cluster.SetXMin((unsigned int) 	number[2]);
// 			cluster.SetXMax((unsigned int) 	number[3]);
// 			
// 			cluster.SetYMin((unsigned int) 	number[4]);
// 			cluster.SetYMax((unsigned int) 	number[5]);
// 			
// 			cluster.SetZMin((unsigned int) 	number[6]);
// 			cluster.SetZMax((unsigned int) 	number[7]);
// 			
// 			landmark.push_back(number[ 8]);
// 			landmark.push_back(number[ 9]);
// 			landmark.push_back(number[10]);
// 			landmark.push_back(1);
// 			
// 			cluster.setLandmark(landmark);
// 			
// 			index[0]=(unsigned int)(number[ 8]+0.5);
// 			index[1]=(unsigned int)(number[ 9]+0.5);
// 			index[2]=(unsigned int)(number[10]+0.5);
// 			
// 			cluster.SetStartPixel(index);
// 			cluster.InitLocalHistogramm();
// 			
// 			
// 			if(isClusterValid(&cluster))// && cluster.GetID()==22)
// 			retList.push_back(cluster);	
// 		}
// 		
// 		filereader.close();
// 	}
// 	else
// 	{	
// 		std::cout<<"FilereaderError: "<< _filename << std::endl; 
// 	}
// 	
// 	retList.sort();
// 	
// 	return retList;
// }

// void Basics::DeleteCluster(CellCluster * cluster)
// {
// 	std::stringstream iStream;
// 	iStream << cluster->GetID();
// 	
// 	std::string filename = outputFilename; 
// 	filename += "Cluster_NR_";
// 	filename += iStream.str();
// 	filename += ".sr";
// 	
// 	const char* out = filename.c_str();
// 	
// 	std::remove(out);
// }

// void Basics::SaveCluster(CellCluster * cluster, ImageType::Pointer image)
// {
// 	std::stringstream iStream;
// 	iStream << cluster->GetID();
// 	unsigned char color;
// 	
// 	std::string filename = outputFilename; 
// 	filename += "Cluster_NR_";
// 	filename += iStream.str();
// 	filename += ".sr";
// 	
// 	const char* out = filename.c_str();
// 	
// 	std::ofstream filewriter( out, std::ios::in | std::ios::out | std::ios::trunc );
// 	ConstIteratorType iter(image, image->GetLargestPossibleRegion());
// 	
// 	if( !filewriter.fail() )
// 	{
// 		
// 		for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
// 		{
// 			color=iter.Get();
// 			filewriter << color;		
// 		}
// 		
// 		filewriter.close();
// 	}
// }

void Basics::writeBinaryImage(ImageType::Pointer binaryImage, unsigned int objectNumber)
{
	writeBinaryImage(binaryImage, "_Obj_Nr_", objectNumber);
}

void Basics::writeImagePlanes()
{
	Writer2DType::Pointer writer = Writer2DType::New();
	writer->SetInput( inputImage );
	
	NameGeneratorType::Pointer writerNameGen = NameGeneratorType::New();
	
	std::string format = outputFilename;
	format += "_%04d.";
	format += "tif";
	
	writerNameGen->SetSeriesFormat( format.c_str() );
	
	ImageType::RegionType region = inputImage->GetLargestPossibleRegion();
	ImageType::IndexType start = region.GetIndex();
	ImageType::SizeType size = region.GetSize();
	
	const unsigned int firstSlice = start[2];
	const unsigned int lastSlice  = start[2] + size[2] -1;
	
	writerNameGen->SetStartIndex( firstSlice );
	writerNameGen->SetEndIndex( lastSlice );
	writerNameGen->SetIncrementIndex( 1 );
	
	writer->SetFileNames( writerNameGen->GetFileNames() );
	
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err )
	{
		std::cerr << "WriterExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}
};

void Basics::writeBinaryImage(ImageType::Pointer binaryImage, const char * label, unsigned int objectNumber)
{
	Writer2DType::Pointer ImageWriter = Writer2DType::New();
	std::stringstream iStream;
	ImageWriter->SetInput( binaryImage );
	NameGeneratorType::Pointer ImageWriterNameGen = NameGeneratorType::New();
	
	iStream << objectNumber;
	std::string ImageFormat = this->outputFilename;
	ImageFormat += label;
	ImageFormat += iStream.str();
// 	ImageFormat += "_";
// 	ImageFormat += "%04d";
	ImageFormat += ".tif";	
	
	ImageWriterNameGen->SetSeriesFormat( ImageFormat.c_str() );
	
	ImageType::RegionType region = binaryImage->GetLargestPossibleRegion();
	
	ImageType::IndexType start = region.GetIndex();
	ImageType::SizeType size = region.GetSize();
	
	const unsigned int firstSlice = start[2];
	const unsigned int lastSlice  = start[2] + size[2] -1;
	ImageWriterNameGen->SetStartIndex( firstSlice );
	ImageWriterNameGen->SetEndIndex( lastSlice );
	ImageWriterNameGen->SetIncrementIndex( 1 );
	
	ImageWriter->SetFileNames(ImageWriterNameGen->GetFileNames() );
	
	try
	{
		ImageWriter->Update();
	}
	catch (itk::ExceptionObject & err )
	{
		std::cout << "binaryWriterExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
	}
}

void Basics::writeMeanProjectionImage(const char * name)
{
	char filename[1024];
	
	ImageType::RegionType 			imageRegion;
	ImageType::RegionType::IndexType 	imageIndex;
	ImageType::RegionType::SizeType 	imageSize;
	
	Image2DType::RegionType 		pRegion;
	Image2DType::RegionType::IndexType 	pIndex;
	Image2DType::RegionType::SizeType	pSize;
	
	unsigned int x_size = maximum_region.GetSize(0);
	unsigned int y_size = maximum_region.GetSize(1);
	unsigned int z_size = maximum_region.GetSize(2);
	
	//unsigned int length=x_size*y_size;
	unsigned int ** colorMap = new unsigned int * [x_size];
	
	for(int i=0;i<x_size;i++)
		colorMap[i]=new unsigned int[y_size];
	
	unsigned char color;
	
	for(int x=0; x< x_size; x++)
		for(int y=0; y< y_size; y++)
			colorMap[x][y]=0;
		
		pIndex[0]=0;
	pIndex[1]=0;
	
	pSize[0]=x_size;
	pSize[1]=y_size;
	
	pRegion.SetIndex(pIndex);
	pRegion.SetSize(pSize);
	
	Image2DType::Pointer projectionImage = Image2DType::New();
	projectionImage->SetRegions(pRegion);
	projectionImage->Allocate();
	projectionImage->FillBuffer(0);
	projectionImage->Update();
	
	for(unsigned int i=0; i< z_size; ++i)
	{
		std::cout << "Writing MeanProjectionImage plane " << i+1 << " of " << z_size << "               \r";
		imageIndex[0]=0;
		imageIndex[1]=0;
		imageIndex[2]=i;
		
		imageSize[0]=maximum_region.GetSize(0);
		imageSize[1]=maximum_region.GetSize(1);
		imageSize[2]=1;
		
		imageRegion.SetSize(imageSize);
		imageRegion.SetIndex(imageIndex);
		
		ConstIteratorType imageIt(inputImage, imageRegion);
		
		
		ImageType::IndexType pos;
		
		for(imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt)
		{
			color=imageIt.Get();
			pos=imageIt.GetIndex();
			colorMap[pos[0]][pos[1]]+=color;
		}
	}
	
	unsigned int maxColor=0;
	
	for(int x=0; x< x_size; x++)
	{
		for(int y=0; y< y_size; y++)
		{
			colorMap[x][y]/=z_size;
			
			if(colorMap[x][y]>maxColor)
				maxColor=colorMap[x][y];
		}
	}
	
	for(int x=0; x< x_size; x++)
	{
		for(int y=0; y< y_size; y++)
		{
			colorMap[x][y]*=255;
			colorMap[x][y]/=maxColor;
		}
	}
	
	Image2DType::RegionType::IndexType 	pixelIndex;
	for(int x=0; x< x_size; x++)
	{
		for(int y=0; y< y_size; y++)
		{
			pixelIndex[0]=x;
			pixelIndex[1]=y;
			
			color=colorMap[x][y];
			
			projectionImage->SetPixel(pixelIndex, color);
		}
	}
	
	Single2DWriterType::Pointer 	planeWriter 	= Single2DWriterType::New();
	
	strcpy(filename, outputFilename);
	strcat(filename, name);
	
	planeWriter->SetFileName(filename);
	planeWriter->SetInput(projectionImage);
	
	try
	{
		planeWriter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "WriterrExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}
	
	std::cout << "\tWriting MeanProjectionImage:\tdone                          \n";
}

void Basics::writeMaxProjectionImage(const char * name)
{
	char filename[1024];
	
	ImageType::RegionType 			imageRegion;
	ImageType::RegionType::IndexType 	imageIndex;
	ImageType::RegionType::SizeType 	imageSize;
	
	Image2DType::RegionType 		pRegion;
	Image2DType::RegionType::IndexType 	pIndex;
	Image2DType::RegionType::SizeType	pSize;
	
	unsigned int length=maximum_region.GetSize(2);
	
	pIndex[0]=0;
	pIndex[1]=0;
	
	pSize[0]=maximum_region.GetSize(0);
	pSize[1]=maximum_region.GetSize(1);
	
	pRegion.SetIndex(pIndex);
	pRegion.SetSize(pSize);
	
	Image2DType::Pointer projectionImage = Image2DType::New();
	projectionImage->SetRegions(pRegion);
	projectionImage->Allocate();
	projectionImage->FillBuffer(0);
	projectionImage->Update();
	
	for(unsigned int i=0; i< length; ++i)
	{
		std::cout << "Writing projectionImage plane " << i+1 << " of " << length << "                \r";
		imageIndex[0]=0;
		imageIndex[1]=0;
		imageIndex[2]=i;
		
		imageSize[0]=maximum_region.GetSize(0);
		imageSize[1]=maximum_region.GetSize(1);
		imageSize[2]=1;
		
		imageRegion.SetSize(imageSize);
		imageRegion.SetIndex(imageIndex);
		
		Iterator2DType	pIt(projectionImage, pRegion);
		ConstIteratorType imageIt(inputImage, imageRegion);
		
		unsigned char imageColor, pColor;
		
		for(pIt.GoToBegin(), imageIt.GoToBegin(); !pIt.IsAtEnd(); ++pIt,++imageIt)
		{
			imageColor=imageIt.Get();
			pColor=pIt.Get();
			
			if(imageColor>pColor)
				pIt.Set(imageColor);
		}
	}
	
	
	Single2DWriterType::Pointer 	planeWriter 	= Single2DWriterType::New();
	
	strcpy(filename, outputFilename);
	strcat(filename, name);
	
	planeWriter->SetFileName(filename);
	planeWriter->SetInput(projectionImage);
	
	try
	{
		planeWriter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "WriterExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}
	std::cout << "\tWriting projectionImage plane:\tdone                   \n";
}

NeighborhoodOffsetVectorType Basics::CreateLookUpTable()
{
	SegNeighborhoodIteratorType::OffsetType offset;
	NeighborhoodOffsetVectorType look_up_table(0);
	
	for(int z = -1; z <= 1 ; z++)
	{
		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				if( (x == 0) && (y == 0) && (z == 0) )
				{}
				else
				{
					offset[0]=x;
					offset[1]=y;
					offset[2]=z;
					look_up_table.push_back(offset);
				}
			}
		}
	}
	
	return look_up_table;
};

void Basics::createNeighborhood()
{
	ShapedNeighborhoodIteratorType::OffsetType offset4, offset8;
	
	for(int x=-1; x<=1; x++)
		for(int y=-1; y<=1; y++)
		{
			if( (x == 0) && (y == 0) )
				continue;
			else
			{
				offset8[0] = x;
				offset8[1] = y;
				offset8[2] = 0;
				neighborhood8.push_back(offset8);
			}
			if( (x && !y) || (!x && y) )
			{
				offset4[0] = x;
				offset4[1] = y;
				offset4[2] = 0;
				neighborhood4.push_back(offset4);
			}
		}
};

void Basics::computeBounds(ImageType::Pointer image)
{
	boundingBox = new long[6];
	boundingBox[0] = image->GetLargestPossibleRegion().GetIndex()[0];
	boundingBox[1] = boundingBox[0] + image->GetLargestPossibleRegion().GetSize()[0] - 1;
	boundingBox[2] = image->GetLargestPossibleRegion().GetIndex()[1];
	boundingBox[3] = boundingBox[2] + image->GetLargestPossibleRegion().GetSize()[1] - 1;
	boundingBox[4] = image->GetLargestPossibleRegion().GetIndex()[2];
	boundingBox[5] = boundingBox[4] + image->GetLargestPossibleRegion().GetSize()[2] - 1;
};

bool Basics::IsInBounds(ImageType::Pointer image, ImageType::IndexType _index)
{
	ImageType::RegionType region=image->GetLargestPossibleRegion();
	ImageType::IndexType index=region.GetIndex();
	ImageType::SizeType size=region.GetSize();
	#ifdef DEBUG
	std::cout << "Index:" << _index << "start:" << index << "size:" << size <<  std::endl;
	#endif
	if(_index[0]<index[0])		return false;
	if(_index[1]<index[1])		return false;
	if(_index[2]<index[2])		return false;
	
	if(_index[0]>size[0])		return false;
	if(_index[1]>size[1])		return false;
	if(_index[2]>size[2])		return false;
	
	return true;
};

std::list<std::list<std::vector<float> > > Basics::GetAmiraContours(std::list<std::list< std::vector< float > > > contours)
{
	std::list<std::list<std::vector<float> > > amira_edges;
	std::list<std::vector<float> > tmp_edge;
	
	std::list<std::list<std::vector<float> > >::iterator it;
	for(it=contours.begin();it!=contours.end();++it)
	{
		tmp_edge = GetAmiraContour(*it, 0);
		amira_edges.push_back(tmp_edge);
	}
	
	return amira_edges;
};

std::list<std::list<std::vector<float> > > Basics::GetAmiraContours(Contour * contours, int nrOfContours)
{
	std::list<std::list<std::vector<float> > > amira_edges;
	std::list<std::vector<float> > tmp_edge;
	
	for(int ii = 0; ii < nrOfContours; ++ii)
	{
		tmp_edge = GetAmiraContour(*(contours[ii].edgeListPointer()), 0);
		amira_edges.push_back(tmp_edge);
	}
	
	return amira_edges;
};

std::list<std::list<std::vector<float> > > Basics::GetAmiraContours(std::vector< Contour * > contours, int nrOfContours)
{
	std::list<std::list<std::vector<float> > > amira_edges;
	std::list<std::vector<float> > tmp_edge;
	
	for(int ii = 0; ii < nrOfContours; ++ii)
	{
		tmp_edge = GetAmiraContour(*(contours[ii]->edgeListPointer()), 0);
		amira_edges.push_back(tmp_edge);
	}
	
	return amira_edges;
};

std::list<std::list<std::vector<float> > > Basics::GetAmiraContours(std::vector< std::vector< Contour * > > contours, int nrOfContours)
{
	std::list<std::list<std::vector<float> > > amira_edges;
	
	for(int ii = 0; ii < nrOfContours; ++ii)
	{
		if(contours[0][ii]->getValid())
		{
			for(int z = 0; z < contours.size(); ++z)
			{
				std::list<std::vector<float> > tmp_edge;
				tmp_edge = GetAmiraContour(*(contours[z][ii]->edgeListPointer()), 1);
				amira_edges.push_back(tmp_edge);
			}
		}
	}
	
	return amira_edges;
};

//worldCoordinates: 1 if contour is already in world coordinates, 0 otherwise
std::list<std::vector<float> >  Basics::GetAmiraContour(std::list<std::vector<float> >& structure, bool worldCoordinates)
{
	std::list<std::vector<float> > tmp_edge;
	
	std::list< std::vector< float > >::const_iterator it = structure.begin();
	for (it; it != structure.end(); ++it)
	{
		std::vector<float> contour_point;
		if(!worldCoordinates)
		{
			contour_point.push_back((*it)[X_COORD]*XYSAMPLING*DOWNSAMPLING);
			contour_point.push_back((*it)[Y_COORD]*XYSAMPLING*DOWNSAMPLING);
			contour_point.push_back((*it)[Z_COORD]*ZSAMPLING);
		}
		else
		{
			contour_point.push_back((*it)[X_COORD]);
			contour_point.push_back((*it)[Y_COORD]);
			contour_point.push_back((*it)[Z_COORD]);
		}
		
		tmp_edge.push_back(contour_point);
	}
	
	std::vector<float> contour_point;
	it = structure.begin();
	if(it != structure.end())
	{
		if(!worldCoordinates)
		{
			contour_point.push_back((*it)[X_COORD]*XYSAMPLING*DOWNSAMPLING);
			contour_point.push_back((*it)[Y_COORD]*XYSAMPLING*DOWNSAMPLING);
			contour_point.push_back((*it)[Z_COORD]*ZSAMPLING);
		}
		else
		{
			contour_point.push_back((*it)[X_COORD]);
			contour_point.push_back((*it)[Y_COORD]);
			contour_point.push_back((*it)[Z_COORD]);
		}
		tmp_edge.push_back(contour_point);
	}
	
	return tmp_edge;
};



Contour::Contour()
{
	attributes.reserve(1);
};

Contour::Contour(std::list<std::vector<float> > _edge_list)
{
	this->edgeList = _edge_list;
};

Contour::Contour(Contour * otherContour)
{
	std::flush(std::cout << "assigning edge list..." << std::endl);
	if(otherContour->edgeList.size())
	{
		this->edgeList.assign(otherContour->edgeList.begin(), otherContour->edgeList.end());
// 		std::list< std::vector< float > >::const_iterator edgeIt;
// 		for(edgeIt = otherContour->edge_list.begin(); edgeIt != otherContour->edge_list.end(); ++edgeIt)
// 			this->edge_list.push_back(*edgeIt);
	}
	std::flush(std::cout << "assigning attribute list..." << std::endl);
	if(otherContour->attributes.size())
	{
		this->attributes.assign(otherContour->attributes.begin(), otherContour->attributes.end());
// 		std::vector< float >::const_iterator attrIt;
// 		for(attrIt = otherContour->attributes.begin(); attrIt != otherContour->attributes.end(); ++attrIt)
// 			this->attributes.push_back(*attrIt);
	}
	this->optimize = otherContour->getOptimizeFlag();
	this->validContour = otherContour->getValid();
	this->barrelID = otherContour->getBarrelID();
};

Contour::~Contour()
{
	this->edgeList.clear();
	this->attributes.clear();
}

void Contour::setEdgeList(std::list<std::vector<float> > * _edge_list)
{
	edgeList.clear();
	std::list<std::vector<float> >::const_iterator it;
	for(it = _edge_list->begin(); it != _edge_list->end(); ++it)
		edgeList.push_back(*it);
};


void Contour::replaceEdgeList(std::list<std::vector<float> > * _edge_list, float zOffset)
{
	if(edgeList.size())
		edgeList.clear();
	std::list<std::vector<float> >::const_iterator it;
	for(it = _edge_list->begin(); it != _edge_list->end(); ++it)
	{
		edgeList.push_back(*it);
		if(edgeList.rbegin()->size() > 2)
			(*(edgeList.rbegin()))[2] += zOffset;
	}
};

void Contour::prepareForOutput()
{
	contourSmoothing(5);
	contourSampling(10);
};

/****************************************************************************/
/*Implementation of contour smoothing by averaging x- and y- values of each */
/*index with #(2*radius+1) neighbors going into the averaging               */
/*CAUTION: with large radius, natural curves will be systematically shrunk  */
/****************************************************************************/
void Contour::contourSmoothing(int radius)
{
	const unsigned int contourSize = this->edgeList.size();
	if(contourSize <= radius)
		radius = 1;
	std::vector< std::vector< float > > contourVec;
	std::list< std::vector< float > > neighborhood;
	std::list< std::vector< float > >::iterator contourIt;
	std::list< std::vector< float > >::iterator neighborIt;
	
	if(contourSize > 1)
	{
		for(contourIt = this->edgeList.begin(); contourIt != this->edgeList.end(); ++contourIt)
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
		for(ii = 0, contourIt = this->edgeList.begin(); ii < contourSize && contourIt != this->edgeList.end(); ++ii, ++contourIt)
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
void Contour::contourSampling(int samplingRate)
{
	std::list< std::vector< float > >::iterator contourIt;
	int count = 1;
	if(this->edgeList.size() <= samplingRate)
		samplingRate = 1;
	
	for(contourIt = this->edgeList.begin(); contourIt != this->edgeList.end(); ++count)
	{
		if(count == samplingRate)
		{
			count = 0;
			++contourIt;
		}
		else
			contourIt = this->edgeList.erase(contourIt);
	}
};



BarrelMarker::BarrelMarker()
{
	//tbd
};
BarrelMarker::~BarrelMarker()
{
	//tbd
};

/****************************************************************************/
/*read Amira landmark file with manual Barrel centroid markers              */
/****************************************************************************/
void BarrelMarker::readBarrelMarkerFile(const char * markerFilename, ImageType::RegionType region)
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
	marker = PointSetType::New();
	marker = points;
	computeVoronoiDiagram(region);
};

/****************************************************************************/
/*computes the 2D Voronoi diagram from the barrel markers                   */
/****************************************************************************/
void BarrelMarker::computeVoronoiDiagram(ImageType::RegionType imageRegion)
{
	ImageType::Pointer markerImage = ImageType::New();
	ImageType::RegionType markerRegion;
	ImageType::SizeType markerSize;
	ImageType::IndexType markerIndex;
	markerSize = imageRegion.GetSize();
	markerSize[2] = 1;
	markerIndex = imageRegion.GetIndex();
	markerRegion.SetIndex(markerIndex);
	markerRegion.SetSize(markerSize);
	markerImage->SetRegions(markerRegion);
	markerImage->Allocate();
	markerImage->FillBuffer(0);
	IndexIteratorType markerIter(markerImage, markerImage->GetLargestPossibleRegion());
	for(int ii = 0; ii < marker->GetNumberOfPoints(); ++ii)
	{
		PointType seedPoint;
		ImageType::IndexType seedIndex;
		marker->GetPoint(ii, &seedPoint);
		seedIndex[0] = (int)(seedPoint[0] + 0.5), seedIndex[1] = (int)(seedPoint[1] + 0.5), seedIndex[2] = imageRegion.GetIndex(2);
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
	
// 	unsigned int noOfRegions = marker->GetNumberOfPoints();
// 	bool * borderCells = new bool[noOfRegions];
// 	for(int ii = 0; ii < noOfRegions; ++ii)
// 		borderCells[ii] = 0;
// 	ImageType::RegionType topRow, bottomRow, leftColumn, rightColumn;
// 	ImageType::SizeType rowSize, columnSize;
// 	ImageType::IndexType topLeftIndex, topRightIndex, bottomLeftIndex;
// 	rowSize[0] = voronoiMap->GetLargestPossibleRegion().GetSize()[0];
// 	rowSize[1] = 1;
// 	rowSize[2] = 1;
// 	columnSize[0] = 1;
// 	columnSize[1] = voronoiMap->GetLargestPossibleRegion().GetSize()[1];
// 	columnSize[2] = 1;
// 	topLeftIndex = voronoiMap->GetLargestPossibleRegion().GetIndex();
// 	topRightIndex = topLeftIndex;
// 	topRightIndex[0] += rowSize[0] - 1;
// 	bottomLeftIndex = topLeftIndex;
// 	bottomLeftIndex[1] += columnSize[1] - 1;
// 	topRow.SetIndex(topLeftIndex), topRow.SetSize(rowSize);
// 	bottomRow.SetIndex(bottomLeftIndex), bottomRow.SetSize(rowSize);
// 	leftColumn.SetIndex(topLeftIndex), leftColumn.SetSize(columnSize);
// 	rightColumn.SetIndex(topRightIndex), rightColumn.SetSize(columnSize);
// 	// 	leftColumn.Print(std::cout);
// 	// 	rightColumn.Print(std::cout);
// 	// 	topRow.Print(std::cout);
// 	// 	bottomRow.Print(std::cout);
// 	ConstCalcIteratorType borderIter1(voronoiMap, topRow);
// 	ConstCalcIteratorType borderIter2(voronoiMap, bottomRow);
// 	ConstCalcIteratorType borderIter3(voronoiMap, leftColumn);
// 	ConstCalcIteratorType borderIter4(voronoiMap, rightColumn);
// 	for(borderIter1.GoToBegin(); !borderIter1.IsAtEnd(); ++borderIter1)
// 		borderCells[(int)(borderIter1.Get()+0.5)-1] = 1;
// 	for(borderIter2.GoToBegin(); !borderIter2.IsAtEnd(); ++borderIter2)
// 		borderCells[(int)(borderIter2.Get()+0.5)-1] = 1;
// 	for(borderIter3.GoToBegin(); !borderIter3.IsAtEnd(); ++borderIter3)
// 		borderCells[(int)(borderIter3.Get()+0.5)-1] = 1;
// 	for(borderIter4.GoToBegin(); !borderIter4.IsAtEnd(); ++borderIter4)
// 		borderCells[(int)(borderIter4.Get()+0.5)-1] = 1;
// 	
// 	ConstCalcIteratorType distIter(distanceMap, distanceMap->GetLargestPossibleRegion());
// 	CalcIteratorType maxDistIter(voronoiMap, voronoiMap->GetLargestPossibleRegion());
// 	for(distIter.GoToBegin(), maxDistIter.GoToBegin(); !distIter.IsAtEnd() && !maxDistIter.IsAtEnd(); ++distIter, ++maxDistIter)
// 		if(borderCells[(int)(maxDistIter.Get()+0.5)-1])
// 			if(distIter.Get() > 200)
// 				maxDistIter.Set(0);
// 	voronoiMap->Update();
	
// 	CalcToImageRescaleFilterType::Pointer cToIFilter = CalcToImageRescaleFilterType::New();
// 	cToIFilter->SetOutputMinimum(1);
// 	cToIFilter->SetOutputMaximum(marker->GetNumberOfPoints());
// 	cToIFilter->SetInput(voronoiMap);
// 	cToIFilter->Update();
// 	writeBinaryImage(cToIFilter->GetOutput(), "_voronoi_map", 0);
};

