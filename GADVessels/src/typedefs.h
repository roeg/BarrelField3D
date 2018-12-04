/****************************************************************************/
/*                                                                          */
/* File:      typedefs.h                                                    */
/*                                                                          */
/* Purpose:   header file for all necessary inculdes and typedefs for files */
/*            associated with the NeuroMorph or CellCount projects          */
/*                                                                          */
/* Author:    Marcel Oberlaender                                            */
/*            Max-Planck-Institute for Neurobiologie                        */
/*            Am Kolpferspitz 18                                            */
/*            D-82152 Martinsried (Munich)                                  */
/*                                                                          */
/* Co-Author: Robert Egger                                                  */
/*            Max-Planck-Institute for Medical Research                     */
/*            Jahnstrasse 19                                                */
/*            D-69120 Heidelberg                                            */
/*                                                                          */
/* EMail:     regger@mpimf-heidelberg.mpg.de                                */
/*                                                                          */
/* History:   17.01.2008                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#ifndef TYPEDEF
#define TYPEDEF

//#define DEBUG


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.1415926

#define _CRT_SECURE_NO_DEPRECATE
#define _SCL_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include "string"

#include "itkImage.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <complex>
#include <utility>
#include <ctime>

// for Reading and Writing Images
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
// for simple Filters
#include "itkBinaryThresholdImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkImportImageFilter.h"
#include "itkMultiplyImageFilter.h"
//FFT Filters
#include "itkVnlFFTRealToComplexConjugateImageFilter.h"
#include "itkVnlFFTComplexConjugateToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"
// for Morphology Filters
#include "itkGrayscaleErodeImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkBinaryCrossStructuringElement.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkSubtractImageFilter.h"
#include "itkRGBPixel.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkNeighborhoodConnectedImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkReconstructionByDilationImageFilter.h"
#include "itkReconstructionByErosionImageFilter.h"
#include "itkWatershedImageFilter.h"
#include "itkWatershedSegmenter.h"
#include "itkWatershedMiniPipelineProgressCommand.h"
#include "itkWatershedSegmentTreeGenerator.h"
// for Reading and Writing of Stacks
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkTIFFImageIO.h"
#include "itkBMPImageIO.h"
#include "itkPNGImageIO.h"
#include "itkRGBPixel.h"
#include "itkLandmarkSpatialObject.h"
#include "itkSpatialObjectWriter.h"
#include "itkSpatialObjectReader.h"
#include "itkMetaLandmarkConverter.h"
// for Iterators
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
// for Histograms and Statistics...
#include "itkScalarImageToListAdaptor.h"
#include "itkListSampleToHistogramGenerator.h"
#include "itkCastImageFilter.h"
#include "itkJointDomainImageToListAdaptor.h"
#include "itkImageToListAdaptor.h"
#include "itkScalarToArrayCastImageFilter.h"
#include "itkSubsample.h"
#include "itkMeanCalculator.h"
#include "itkCovarianceCalculator.h"
#include "itkFixedArray.h"
// ... and Clusters
#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkMinimumDecisionRule.h"
#include "itkEuclideanDistance.h"
#include "itkSampleClassifier.h"


#include "itkPasteImageFilter.h"
#include "itkNoiseImageFilter.h"
//-------------------------------------------------------------------------------------------------------------------
//Type Definitions for 3D Grayvalue Images and Iterators
//-------------------------------------------------------------------------------------------------------------------
  
  typedef unsigned char						PixelType;
  typedef unsigned long						ObjectPixelType;
  typedef float 						CalcPixelType;

  typedef itk::Image< PixelType, 3 >				ImageType;
  typedef itk::Image< ObjectPixelType, 3 >			ObjectImageType;
  typedef itk::Image< PixelType, 2 >				Image2DType;  
  typedef itk::Image< ObjectPixelType, 2 >			ObjectImage2DType;
  typedef itk::Image< CalcPixelType, 3 >			CalcImageType;
  typedef itk::Image< CalcPixelType, 2 >			CalcImage2DType;
  
  typedef itk::ConstantBoundaryCondition< ImageType >		ConstantBoundaryConditionType;

  typedef itk::PasteImageFilter< ImageType >						PasteFilterType;
  typedef itk::PasteImageFilter< CalcImageType, CalcImageType, CalcImageType >		PasteCalcFilterType;
  typedef itk::PasteImageFilter< Image2DType, Image2DType, Image2DType >		Paste2DFilterType;
  typedef itk::PasteImageFilter< CalcImage2DType, CalcImage2DType, CalcImage2DType >	Paste2DCalcFilterType;
  typedef itk::NoiseImageFilter< ImageType, ImageType > NoiseFilterType;
  //typedef itk::PasteImageFilter< ImageType,Image2DType,ImageType >		Paste2Dto3DFilterType;

  typedef itk::ImageFileReader< ImageType >			ReaderType;
  typedef itk::ImageFileReader< Image2DType >			Reader2DType;
  typedef itk::ImageFileReader< ObjectImage2DType >		ObjectReader2DType;
  typedef itk::ImageFileReader< CalcImage2DType >		CalcReader2DType;
  typedef itk::ImageSeriesReader< ImageType >			SeriesReaderType;
  typedef itk::NumericSeriesFileNames				NameGeneratorType;
  typedef itk::ImageFileWriter< ImageType >			WriterType;
  typedef itk::ImageFileWriter< ObjectImageType >		ObjectWriterType;
  typedef itk::ImageSeriesReader< ImageType > 			SeriesReaderType;
  typedef itk::ImageSeriesWriter< ImageType, Image2DType >	Writer2DType;
  typedef itk::ImageFileWriter< Image2DType >			Single2DWriterType;
  typedef itk::SpatialObjectWriter< 3 >				LandmarkWriterType;
  typedef itk::SpatialObjectReader< 3 >				LandmarkReaderType;
  
  typedef std::complex< CalcPixelType >									ComplexType;
  typedef itk::VnlFFTRealToComplexConjugateImageFilter< CalcPixelType, 2 > 				ForwardFFT2DFilterType;
  typedef itk::VnlFFTComplexConjugateToRealImageFilter< CalcPixelType, 2 >				InverseFFT2DFilterType;
  typedef ForwardFFT2DFilterType::OutputImageType 							ComplexImage2DType;
  typedef itk::ComplexToImaginaryImageFilter< ComplexImage2DType, CalcImage2DType >			Imaginary2DFilterType;
  typedef itk::ComplexToRealImageFilter< ComplexImage2DType, CalcImage2DType >				Real2DFilterType;
  typedef itk::ComplexToModulusImageFilter< ComplexImage2DType, CalcImage2DType >			Modulus2DFilterType;
  typedef itk::MultiplyImageFilter< ComplexImage2DType, ComplexImage2DType, ComplexImage2DType >	MultiplyComplexImage2DFilterType;
  typedef itk::VnlFFTRealToComplexConjugateImageFilter< CalcPixelType, 3 > 				ForwardFFTFilterType;
  typedef itk::VnlFFTComplexConjugateToRealImageFilter< CalcPixelType, 3 >				InverseFFTFilterType;
  typedef ForwardFFTFilterType::OutputImageType 							ComplexImageType;
  typedef itk::ComplexToImaginaryImageFilter< ComplexImageType, CalcImageType >				ImaginaryFilterType;
  typedef itk::ComplexToRealImageFilter< ComplexImageType, CalcImageType >				RealFilterType;
  typedef itk::ComplexToModulusImageFilter< ComplexImageType, CalcImageType >				ModulusFilterType;
  typedef itk::MultiplyImageFilter< ComplexImageType, ComplexImageType, ComplexImageType >		MultiplyComplexImageFilterType;

  typedef itk::RelabelComponentImageFilter< ObjectImageType, ObjectImageType >		RelabelType;
  typedef itk::RelabelComponentImageFilter< ObjectImage2DType, ObjectImage2DType >	Relabel2DType;
  typedef itk::BinaryThresholdImageFilter< ObjectImageType, ImageType >			ObjectFilterType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType >			BinaryFilterType;
  typedef itk::ConnectedComponentImageFilter< ImageType, ObjectImageType >		ConnectedFilterType;
  typedef itk::ConnectedComponentImageFilter< ObjectImageType, ObjectImageType >	ObjectConnectedFilterType;
  typedef itk::ConnectedComponentImageFilter< Image2DType,ObjectImage2DType >		Connected2DFilterType;
  typedef itk::NeighborhoodConnectedImageFilter< ImageType, ImageType >			ConnectedNeighborhoodFilterType;
  
  typedef itk::ThresholdImageFilter< ImageType >		ThresholdFilter;

  typedef itk::ImageRegionConstIterator< ImageType >			ConstIteratorType;
  typedef itk::ImageRegionConstIteratorWithIndex< ImageType >		ConstIndexIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType >		IndexIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< Image2DType >		IndexIterator2DType;
  typedef itk::ImageRegionConstIterator< Image2DType >			Const2DIteratorType;
  typedef itk::ImageRegionIterator< ObjectImageType >			IteratorType;
  typedef itk::ImageRegionConstIterator< ObjectImageType >		ConstObjectIteratorType;
  typedef itk::ImageRegionConstIteratorWithIndex< ObjectImageType >	ConstObjectIndexIteratorType;
  typedef itk::ImageRegionConstIterator< CalcImageType >		ConstCalcIteratorType;
  typedef itk::ImageRegionIterator< CalcImageType >			CalcIteratorType;
  typedef itk::ImageRegionIterator< CalcImage2DType >			CalcIterator2DType;
  typedef itk::ImageRegionIteratorWithIndex< CalcImage2DType >		CalcIndexIterator2DType;
  typedef itk::ImageRegionIteratorWithIndex< CalcImageType >		CalcIndexIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< ComplexImage2DType >	ComplexIndexIterator2DType;
  typedef itk::ImageRegionConstIterator< ComplexImage2DType >		ConstComplexIterator2DType;
  typedef itk::ImageRegionIterator< ComplexImageType >			ComplexIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< ComplexImageType >		ComplexIndexIteratorType;
  typedef itk::ImageRegionConstIterator< ComplexImageType >		ConstComplexIteratorType;

  typedef itk::ConstNeighborhoodIterator< ImageType >					ConstNeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< ImageType>						IteratorType2;
  typedef itk::ImageRegionIterator< ObjectImageType>					ObjectIteratorType;
  typedef itk::NeighborhoodIterator< ImageType >					SegNeighborhoodIteratorType;
  typedef itk::NeighborhoodIterator< Image2DType >					SegNeighborhoodIterator2DType;
  typedef itk::NeighborhoodIterator< ObjectImageType >					NeighborhoodIteratorType;
  typedef itk::NeighborhoodIterator< CalcImageType >					CalcNeighborhoodIteratorType;
  typedef itk::ShapedNeighborhoodIterator< ImageType >					ShapedNeighborhoodIteratorType;
  typedef itk::ConstNeighborhoodIterator< ObjectImageType >				ConstObjectNeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< Image2DType >					Iterator2DType;
  typedef itk::ImageRegionConstIterator< Image2DType >					ConstIterator2DType;
  typedef itk::ImageRegionIterator< ObjectImage2DType >					ObjectIterator2DType;
  typedef itk::ImageRegionConstIteratorWithIndex<ObjectImage2DType>			ConstObjectIndexIterator2DType;

  typedef itk::RescaleIntensityImageFilter< ObjectImageType, ImageType > 	RescaleFilterType;
  typedef itk::RescaleIntensityImageFilter< CalcImageType,ImageType >   	RescalerCalcType; 
  typedef itk::RescaleIntensityImageFilter< ImageType, CalcImageType >		ImageToCalcRescaleFilterType;
  typedef itk::RescaleIntensityImageFilter< CalcImageType, ImageType >		CalcToImageRescaleFilterType;
  typedef itk::RescaleIntensityImageFilter< Image2DType, CalcImage2DType >	ImageToCalc2DRescaleFilterType;
  typedef itk::RescaleIntensityImageFilter< CalcImage2DType, Image2DType >	CalcToImage2DRescaleFilterType;
  typedef itk::SigmoidImageFilter< ImageType, ImageType >			SigmoidFilterType;
  typedef itk::RegionOfInterestImageFilter< ImageType, ImageType >		RegionFilterType;
  typedef itk::RegionOfInterestImageFilter< CalcImageType, CalcImageType >	CalcRegionFilterType;
  typedef itk::RegionOfInterestImageFilter< CalcImage2DType, CalcImage2DType >	CalcRegionFilter2DType;
  typedef itk::GradientMagnitudeImageFilter< CalcImageType, CalcImageType >	GradientFilterType;

  typedef itk::BinaryCrossStructuringElement< PixelType, 3 >					StructuringElementCrossType;
  typedef itk::BinaryErodeImageFilter< ImageType, ImageType, StructuringElementCrossType >	BinaryErodeFilterCrossType;
  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementCrossType >	BinaryDilateFilterCrossType;

  typedef itk::BinaryBallStructuringElement< PixelType, 3 >					StructuringElementBallType;
  typedef itk::BinaryErodeImageFilter< ImageType, ImageType, StructuringElementBallType >	BinaryErodeFilterBallType;
  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementBallType >	BinaryDilateFilterBallType;
//1 Ball, 0 Cross

  typedef itk::GrayscaleErodeImageFilter< ImageType, ImageType, StructuringElementBallType >  	ErodeFilterType;
  typedef itk::GrayscaleDilateImageFilter< ImageType, ImageType, StructuringElementBallType > 	DilateFilterType;
  typedef itk::GrayscaleMorphologicalOpeningImageFilter< ImageType, ImageType, StructuringElementBallType > OpeningFilterType;

  typedef itk::MedianImageFilter< ImageType, ImageType >					MedianFilterType;
  typedef itk::ReconstructionByDilationImageFilter< CalcImageType, CalcImageType >		ReconstructionByDilationFilterType;
  typedef itk::ReconstructionByErosionImageFilter< CalcImageType, CalcImageType >		ReconstructionByErosionFilterType;	

  typedef float						MeasurementVectorBase;
  typedef itk::Vector< MeasurementVectorBase, 4 >	MeasurementVectorType;
  typedef itk::Vector< unsigned char, 1 >		MeasurementHistogramType;
  typedef itk::Vector< unsigned long, 1 >		MeasurementObjectHistogramType;
  typedef itk::FixedArray< CalcPixelType, 3 >		CalcPixelArrayType;

  typedef itk::Vector< unsigned long, 3 >					MeasurementClusterVectorType;
  typedef itk::Statistics::ListSample< MeasurementClusterVectorType >		ClusterSampleType;
  typedef itk::Statistics::WeightedCentroidKdTreeGenerator< ClusterSampleType > TreeGeneratorType;
  typedef TreeGeneratorType::KdTreeType						TreeType;
  typedef itk::Statistics::KdTreeBasedKmeansEstimator< TreeType >		EstimatorType;

//KMEANS LABELING
  typedef itk::Statistics::ListSample< MeasurementClusterVectorType >		SampleType;
  typedef itk::Statistics::EuclideanDistance< MeasurementClusterVectorType >	MembershipFunctionType;
  typedef itk::MinimumDecisionRule DecisionRuleType;
  typedef itk::Statistics::SampleClassifier< ClusterSampleType > 		ClassifierType;

  typedef itk::Statistics::ListSample< MeasurementHistogramType >	HistogramSampleType;
  typedef itk::Statistics::MeanCalculator< HistogramSampleType >	MeanAlgorithmType;
  typedef itk::Statistics::CovarianceCalculator< HistogramSampleType >	CovarianceAlgorithmType;

  typedef itk::LandmarkSpatialObject< 3 >		LandmarkType;
  typedef itk::SpatialObjectPoint< 3 >			LandmarkPointType;
  typedef itk::MetaLandmarkConverter< 3 >		MetaLandmarkConverterType;

  typedef std::vector< ImageType::IndexType >				IndexVectorType;
  typedef std::vector< ObjectImageType::IndexType >			ObjectIndexVectorType;
  typedef std::vector< ImageType::SizeType >				SizeVectorType;
  typedef std::vector< unsigned long >					LongVectorType;
  typedef std::vector<SegNeighborhoodIteratorType::OffsetType>		NeighborhoodOffsetVectorType;
  typedef std::vector< ShapedNeighborhoodIteratorType::OffsetType >	ShapedNeighborhoodOffsetVectorType;
  typedef std::map< unsigned long, std::list< unsigned long > >		GraphType;
  typedef std::map< long, std::list< long > >				SignedGraphType;

  extern float XYSAMPLING;
  extern float ZSAMPLING;
  extern float averageSomaRadius;
  extern float zScale;
  extern unsigned long BINSIZE;
  extern unsigned long BOXSIZE;
  
  //VoxelFeatures of the Measurement Vector Type
  #define X_COORD		0
  #define Y_COORD		1
  #define Z_COORD		2
//168fines
#endif