#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 ) //warning C4996: 'strcpy': This function or variable may be unsafe.
#endif
/*  All header files needed */

// basic itk
#include "itkImage.h"
#include "itkVector.h"

// registration
#include "itkImageRegistrationMethod.h"
#include "itkSymmetricForcesDemonsRegistrationFilter.h"

// transform
#include "itkIdentityTransform.h"
#include "itkTranslationTransform.h"
#include "itkBSplineDeformableTransform.h"

// metric
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkSquaredDifferenceImageFilter.h" //??
// interpolator
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
// optimizer
#include "itkLBFGSOptimizer.h"
#include "itkGradientDescentOptimizer.h"
#include "itkConjugateGradientOptimizer.h"

// reader / writer
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"

// filter
#include "itkResampleImageFilter.h"
#include "itkBSplineResampleImageFunction.h"
#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkBSplineDecompositionImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkDivideByConstantImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkWarpVectorImageFilter.h"
#include "itkInverseDeformationFieldImageFilter.h"
#include "itkIterativeInverseDeformationFieldImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkLaplacianImageFilter.h"
#include "itkDeformationFieldJacobianDeterminantFilter.h"
#include "itkDisplacementFieldJacobianDeterminantFilter.h"


// for affine transformation  // 20100524
#include "itkTransform.h"
#include "itkAffineTransform.h"
#include "itkTransformToDeformationFieldSource.h"
#include "itkImageRegistrationMethod.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredTransformInitializer.h"


// for Diffeomorphic Demons
#include <itkCommand.h>
#include <itkDiffeomorphicDemonsRegistrationFilter.h>
#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include <itkFastSymmetricForcesDemonsRegistrationFilter.h>
#include <itkGridForwardWarpImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMultiResolutionPDEDeformableRegistration.h>
#include <itkTransformToDeformationFieldSource.h>
#include <itkVectorCentralDifferenceImageFunction.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include <itkWarpHarmonicEnergyCalculator.h>

// including itksys::SystemTools::MakeDirectory(char*)
#include <itksys/SystemTools.hxx>
#include <metaCommand.h>
// 
#include <string>
#include <vector>
#include <errno.h>
#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <sstream>  // add 20091216

#if defined(_WIN32) && !defined(__CYGWIN__)

#else

#include <sys/types.h>
#include <unistd.h>
#include <dlfcn.h>
#define DBL_MAX         1.7976931348623158e+308 /* max value */

#endif


// #include <windows.h>
//#include "apcluster.h"
//#include "mvcd.h"
//#include "cres.h"
//#include "matrixSHEN.h"

#define MAX_FILE_NAME_LENGTH 1024

// using namespace std;

// global type definitions
const    unsigned int   ImageDimension = 3;
typedef double			CoordinateRepType;
const   unsigned int	SpaceDimension = ImageDimension;
const	unsigned int	SplineOrder = 3;
const	unsigned int	lowRes = 16;
const	unsigned int	mediumRes = 32;
const	unsigned int	highRes = 64;
//const	unsigned int	lowRes = 8;
//const	unsigned int	mediumRes = 16;
//const	unsigned int	highRes = 32;
const	unsigned int	gridBorderLength = 3;

// basic data type
typedef unsigned char										CharPixelType;     // for image IO usage
typedef float												FloatPixelType;    // for 
//typedef INT32												Int32PixelType;
typedef int													IntPixelType;
typedef short												ShortPixelType;

typedef float												InternalPixelType; // for internal processing usage
typedef itk::Vector< InternalPixelType, ImageDimension >	VectorPixelType;
// basic image type
typedef itk::Image< CharPixelType, ImageDimension >			CharImageType;
//typedef itk::Image< Int32PixelType, ImageDimension >		Int32ImageType;
typedef itk::Image< IntPixelType, ImageDimension >			IntImageType;
typedef itk::Image< ShortPixelType, ImageDimension >		ShortImageType;
typedef itk::Image< FloatPixelType, ImageDimension >		FloatImageType;

typedef itk::Image< InternalPixelType, ImageDimension >		InternalImageType;
typedef itk::Image<  VectorPixelType, ImageDimension >		DeformationFieldType;
// basic iterator type
typedef itk::ImageRegionIterator< DeformationFieldType >	DeformationFieldIteratorType;
typedef itk::ImageRegionIterator< InternalImageType >		InternalImageIteratorType;
typedef itk::ImageRegionIterator< CharImageType >			CharImageIteratorType;
// basic image reader/writer related type
typedef itk::ImageFileReader< CharImageType >						CharImageReaderType;
typedef itk::ImageFileReader< InternalImageType >					InternalImageReaderType;
typedef itk::ImageFileWriter< InternalImageType >					InternalImageWriterType;
typedef itk::CastImageFilter< InternalImageType, CharImageType >	Internal2CharCastFilterType;
typedef itk::CastImageFilter< InternalImageType, FloatImageType>	Internal2FloatCastFilterType;
typedef itk::CastImageFilter< InternalImageType, IntImageType>		Internal2IntCastFilterType;
typedef itk::CastImageFilter< InternalImageType, ShortImageType>	Internal2ShortCastFilterType;

typedef itk::WarpImageFilter< InternalImageType, InternalImageType, DeformationFieldType  >     InternalWarpFilterType;
typedef itk::ImageFileWriter< CharImageType >						CharImageWriterType;
typedef itk::ImageFileWriter< IntImageType >						IntImageWriterType;
typedef itk::ImageFileWriter< FloatImageType >						FloatImageWriterType;
typedef itk::ImageFileWriter< ShortImageType >						ShortImageWriterType;

typedef itk::ImageFileReader< DeformationFieldType >				DeformationFieldReaderType;
typedef itk::ImageFileWriter< DeformationFieldType >				DeformationFieldWriterType;
// basic transform file reader/writer related type
typedef itk::TransformFileWriter						TransformFileWriterType;
typedef itk::TransformFileReader						TransformFileReaderType;
typedef itk::TransformFileReader::TransformListType		TransformListType;
/////////////////////////////////////////////////////////////////////////////
// metric type
typedef itk::MeanSquaresImageToImageMetric< InternalImageType, InternalImageType >  InternalMeanSquaresMetricType;
// interpolator type
typedef itk::LinearInterpolateImageFunction< InternalImageType, double >			InternalLinearInterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction< InternalImageType, double >	InternalNNInterpolatorType;
typedef itk::BSplineInterpolateImageFunction<InternalImageType, double, double>		InternalBSplineInterpolatorType;
// transform type
typedef itk::BSplineDeformableTransform <CoordinateRepType, SpaceDimension, SplineOrder> BSplineTransformType;
typedef itk::TranslationTransform< double, ImageDimension >		TranslationTransformType;
// optimizer type
typedef itk::LBFGSOptimizer										LBFGSOptimizerType;
typedef itk::GradientDescentOptimizer							GradientDescentOptimizerType;
typedef itk::ConjugateGradientOptimizer							ConjugateGradientOptimizerType;
// registration type
typedef itk::ImageRegistrationMethod <InternalImageType, InternalImageType>	ImageRegistrationType;

// affine transform type // 20100524
typedef itk::ImageRegistrationMethod <InternalImageType, InternalImageType>	InternalImageRegistrationType;
typedef itk::AffineTransform <CoordinateRepType, SpaceDimension> AffineTransformType;
typedef itk::TransformToDeformationFieldSource <DeformationFieldType, CoordinateRepType> Transform2DeformationFieldGeneratorType;
typedef itk::RegularStepGradientDescentOptimizer				RegularStepGradientDescentOptimizerType;
typedef itk::ImageRegistrationMethod<InternalImageType, InternalImageType >    ImageRegistrationType;
typedef itk::CenteredTransformInitializer<AffineTransformType, InternalImageType, InternalImageType >  AffineTransformInitializerType;

//////////////////////////////////////////////////////////////////////////////
// image filter type
typedef itk::ResampleImageFilter< InternalImageType, InternalImageType >    ResampleFilterType;
typedef itk::MeanImageFilter< InternalImageType, InternalImageType > InternalMeanFilterType;
typedef itk::HistogramMatchingImageFilter<InternalImageType,InternalImageType >   InternalHistMatchFilterType;

// operation on internal image
typedef itk::AddImageFilter< InternalImageType, InternalImageType, InternalImageType > AddInternalImageFilterType;
typedef itk::MultiplyByConstantImageFilter< InternalImageType, float, InternalImageType > MultiplyInternalImageFilterType;
typedef itk::DivideByConstantImageFilter< InternalImageType, float, InternalImageType >   DivideInternalImageFilterType;

////////////////////////////////////////////////////////////////////////////
// operation on deformation fields
typedef itk::WarpVectorImageFilter< DeformationFieldType, DeformationFieldType, DeformationFieldType>  WarpVectorFilterType;
typedef itk::IterativeInverseDeformationFieldImageFilter<DeformationFieldType, DeformationFieldType >  IterativeInverseDeformationFieldImageFilterType;
typedef itk::InverseDeformationFieldImageFilter<DeformationFieldType, DeformationFieldType >  InverseDeformationFieldImageFilterType;
typedef itk::AddImageFilter< DeformationFieldType, DeformationFieldType, DeformationFieldType > AddDeformationFieldFilterType;
typedef itk::MultiplyByConstantImageFilter< DeformationFieldType, float, DeformationFieldType > MultiplyDeformationFieldFilterType;
typedef itk::DivideByConstantImageFilter< DeformationFieldType, float, DeformationFieldType >   DivideDeformationFieldFilterType;
/////////////////////////////////////////////////////////////////////////////////
typedef itk::AddImageFilter< DeformationFieldType, DeformationFieldType, DeformationFieldType > AddImageFilterType;
typedef itk::MultiplyByConstantImageFilter< DeformationFieldType, float, DeformationFieldType > MultiplyImageFilterType;
typedef itk::DivideByConstantImageFilter< DeformationFieldType, float, DeformationFieldType >   DivideImageFilterType;

// sub-functions

// read and write deformation field
DeformationFieldType::Pointer ReadDeformationField(char* filename);
void ReadDeformationField(char* filename, DeformationFieldType::Pointer *deformationfield);
void ReadDeformationField(char* filename, DeformationFieldType::Pointer &deformationfield);
void WriteDeformationField(char* filename, DeformationFieldType::Pointer deformationfield);

// read and write image
InternalImageType::Pointer ReadImage(char *filename);
void ReadImage(char *filename, InternalImageType::Pointer *image);
void ReadImage(char *filename, InternalImageType::Pointer &image);
void ReadCharImage(char *filename, CharImageType::Pointer &image);
void WriteImage(char *filename, InternalImageType::Pointer image, char* outputType);
void WriteImage(char *filename, InternalImageType::Pointer image);
void WriteImageUCHAR(char *filename, InternalImageType::Pointer image);
void WriteImageINT(char *filename, InternalImageType::Pointer image);
void WriteImageSHORT(char *filename, InternalImageType::Pointer image);
void WriteImageFLOAT(char *filename, InternalImageType::Pointer image);

// read and write transform
void ReadTransform(char *filename, BSplineTransformType::Pointer &bsplineTransform);
void WriteTransform(char *filename, BSplineTransformType::Pointer bsplineTransform);

// compose deformation fields
void ComposeDeformationFields(DeformationFieldType::Pointer input, 
							  DeformationFieldType::Pointer deformationField, 
							  DeformationFieldType::Pointer &composedDeformationField);
void ComposeDeformationFieldsAndSave(char* inputDeformationFieldFileName,
									 char* deformationFieldFileName,
									 char* composedDeformationFieldFileName);
// compose bspline transforms
void ComposeBSplineTransforms(BSplineTransformType::Pointer input,
							  BSplineTransformType::Pointer tranform,
							  BSplineTransformType::Pointer &composedTransform);
void ComposeBSplineTransformsAndSave (char* inputTransformFileName,
									  char* transformFileName,
									  char* composedTransformFileName);
// inverse deformation field by DG's algorithm
void InverseDeformationFieldDG(DeformationFieldType::Pointer deformationField, 
							   DeformationFieldType::Pointer &deformationFieldInverse);
void InverseDeformationFieldDG3D(DeformationFieldType::Pointer deformationField, 
								 DeformationFieldType::Pointer &deformationFieldInverse);
void InverseDeformationFieldDG3D(float*** dfx, float*** dfy, float*** dfz, 
								 int x_size, int y_size, int z_size, 
								 float*** rdfx, float*** rdfy, float*** rdfz);
// apply deformation field on image and write deformed image
void ApplyDeformationFieldAndWriteWithFileNames(char* movingImageName, char* deformationFieldFileName, 
												char* deformedImageName, bool isLinearInterpolator);
void ApplyDeformationField(InternalImageType::Pointer movingImage, 
						   DeformationFieldType::Pointer deformationField,
						   InternalImageType::Pointer &deformedImage,
						   bool isLinearInterpolator);
// apply bspline transform on the moving image
void ApplyTransformAndWriteWithFileNames(char* movingImageName, char* transformFileName, 
										 char* deformedImageName, bool isLinearInterpolator);
void ApplyBSplineTransform(InternalImageType::Pointer movingImage,
						   BSplineTransformType::Pointer bsplineTransform,
						   InternalImageType::Pointer &deformedImage,
						   bool isLinearInterpolator);
// calculate weighted average deformation field
void calculateWeightedAverageDeformationField(DeformationFieldType::Pointer * allDeformationFields, double* weight, 
											  int totalNumber, DeformationFieldType::Pointer &weightedAveragedeformationField);
// Diffeomorphic Demons Registration
void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  char* deformedImageFileName, char* deformationFieldFileName);
void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  char* deformedImageFileName, char* deformationFieldFileName,
							  unsigned int * iter, int res);
void DiffeoDemonsRegistrationWithoutHistMatching(char* fixedImageFileName, char* movingImageFileName,
												 char* deformedImageFileName, char* deformationFieldFileName);
void DiffeoDemonsRegistrationWithoutHistMatching(char* fixedImageFileName, char* movingImageFileName,
												 char* deformedImageFileName, char* deformationFieldFileName,
												 unsigned int * iter, int res);
void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  char* deformedImageFileName, char* deformationFieldFileName, bool doHistMatch);
void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  DeformationFieldType::Pointer &deformationField, bool doHistMatch);
void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  DeformationFieldType::Pointer &deformationField, bool doHistMatch, float sigmaDef);

void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  char* deformedImageFileName, char* deformationFieldFileName, bool doHistMatch,
							  unsigned int * iter, int res);
void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFielFileName,
														 char* deformedImageFileName, char* deformationFieldFileName);
void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFielFileName,
														 char* deformedImageFileName, char* deformationFieldFileName,
														 unsigned int * iter, int res);
void DiffeoDemonsRegistrationWithInitialDeformationFieldWithoutHistMatching(char* fixedImageFileName, char* movingImageFileName, 
																			char* initialDeformationFielFileName,
																			char* deformedImageFileName, char* deformationFieldFileName);
void DiffeoDemonsRegistrationWithInitialDeformationFieldWithoutHistMatching(char* fixedImageFileName, char* movingImageFileName, 
																			char* initialDeformationFielFileName,
																			char* deformedImageFileName, char* deformationFieldFileName,
																			unsigned int * iter, int res);
void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFieldFileName,
														 char* deformedImageFileName, char* deformationFieldFileName,
														 bool doHistMatch);
void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFieldFileName,
														 char* deformedImageFileName, char* deformationFieldFileName,
														 bool doHistMatch, float sigmaDef);

void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFieldFileName,
														 DeformationFieldType::Pointer &deformationField, 
														 bool doHistMatch);
void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFieldFileName,
														 char* deformedImageFileName, char* deformationFieldFileName,
														 bool doHistMatch, unsigned int * iter, int res);
void DiffeoDemonsRegistration(InternalImageType::Pointer fixedImage, 
							  InternalImageType::Pointer movingImage,
							  InternalImageType::Pointer deformedImage,
							  DeformationFieldType::Pointer &deformationField);
void DiffeoDemonsRegistration(InternalImageType::Pointer fixedImage, 
							  InternalImageType::Pointer movingImage,
							  DeformationFieldType::Pointer &deformationField);
// BSpline registration
void BSplineRegistration(char* fixedImageFileName, char* movingImageFileName,
						 char* deformedImageFileName, char* deformationFieldFileName);
void BSplineRegistration(InternalImageType::Pointer fixedImage, InternalImageType::Pointer movingImage,
						 InternalImageType::Pointer &deformedImage, BSplineTransformType::Pointer &bsplineTransform);
void BSplineRegistrationWithInitialTransform(char* fixedImageFileName, char* movingImageFileName, 
											 char* initialDeformationFielFileName,
											 char* deformedImageFileName, char* deformationFieldFileName);
void BSplineRegistrationWithInitialTransform(InternalImageType::Pointer fixedImage, InternalImageType::Pointer movingImage,
											 BSplineTransformType::Pointer initialBSplineTransform, 
											 InternalImageType::Pointer &deformedImage, BSplineTransformType::Pointer &bsplineTransform);
// generate zero deformation field
void GenerateZeroDeformationField(DeformationFieldType::Pointer &deformationFieldZero, int imageSize);
void GenerateZeroDeformationField3D(DeformationFieldType::Pointer &deformationFieldZero, int imageSize_x, int imageSize_y, int imageSize_z)
;
// generate zero bspline transform
void GenerateZeroBSplineTransform(BSplineTransformType::Pointer &transformZero, InternalImageType::Pointer fixedImage, int gridSize);
// Other filter
// histogram matching
void HistogramMatching(InternalImageType::Pointer inputImage, 
					   InternalImageType::Pointer referenceImage, 
					   InternalImageType::Pointer &outputImage);
void HistogramMatchingWithNames(char* inputImageName, char* referenceImageName, char* outputImageName);

// calculate Average Images
void calculateMeanImage(char** imagefiles, int filenumber, InternalImageType::Pointer &meanImage);




// other auxiliary sub-functions
void myitoa(int num, char* str, int digit);
void bubbleSort(double* arr, int* index, int n);
void PairwiseDistanceAmongImages(char** imageFileNames, int totalNumber, double** distanceMatrix);
double calculateDistance(char* imageName1, char* imageName2);
double calculateDistanceMSD(char* imageName1, char* imageName2);
void SaveMatrix2File(double** matrix, int iSize, int jSize, char* martixFileName);
void ReadMatrixFile(double** matrix, int iSize, int jSize, char* martixFileName);
void CalculateIsomapKNN(double** distanceMatrix, int matrixSize, double** distanceIsomapMatrix, int K);
void FindMaxLocationInMatrix(double** distanceMatrix, int matrixSize, int &max, int &neighbor);
void FindMaxLocationInMatrix(double** distanceMatrix, int matrixSize, int* subGroupSize, int &max, int &neighbor);
void FindGraphMaxTree(double** distanceMatrixTotal, int matrixSize, int* tree);
void FindKNNGraphMaxTree(double** distanceMatrixTotal, int matrixSize, int kNN, int* tree);
void FindLeavesInTree(int* tree, int treeSize, bool* isLeave);
void FindRoot(int* tree, int treeSize, int &root);
void CalculateNodeSegmentsToRoot(int* tree, int treeSize, int* height);
void SaveArrayWithIndex(int* arr, int arrSize, char* filename);
void SaveArrayWithIndex(double* arr, int arrSize, char* filename);
void SaveArrayWithIndexSize(int* arr, int arrSize, char* filename);
void SaveArrayWithIndexSize(double* arr, int arrSize, char* filename);
void CalculateTotalPathLength(int* tree, int treeSize, double* edgeLength, double &totalLength);
void CalculateTotalPathLength(double* edgeLength, int* subTreeSize, int treeSize, double & totalLegnth);
void CalculateSubTreeSize(int* tree, int treeSize, int* subTreeSize);
void GetTreeHeight(int* tree, int treeSize, int &treeHeight);
void SaveTreeWithInfo(int* tree, int treeSize, char* filename);
void ReadTreeWithInfo(char* filename, int treeSize, int* tree);
 int callback(double *curr_a, double *curr_r, int N, int *curr_idx, int I, double curr_netsim, int iter);
/*
void APClustering(double** distanceMatrix, int groupSize, int &APclustersize, 
				  int* &APclustersubrep, int* &APclustersubsize, int** &APcluster);
void APClusteringAndSave(double** distanceMatrix, int groupSize, char* APFileName);
void APClusteringWithLibName(char* APLibName, double** distanceMatrix, int groupSize, int &APclustersize, 
				  int* &APclustersubrep, int* &APclustersubsize, int** &APcluster);
void APClusteringAndSaveWithLibName(char* APLibName, double** distanceMatrix, int groupSize, char* APFileName);
void SaveAPClustering(char* APClusteringFileName, int APclustersize, 
					  int* APclustersubrep, int* APclustersubsize, int** APcluster);
void ReadAPClustering(char* APClusteringFileName, int &APclustersize, 
					  int* &APclustersubrep, int* &APclustersubsize, int** &APcluster);
*/
void FindGraphMaxTreeWithInitialWeight(double** distanceMatrixTotal, int matrixSize, int* initialWeight, int* tree);
// float DistanceMetric

// void AverageDeformationFields

// void AverageDeformationFieldsWithWeights

// do Affinity Propagation Clustering
/*
void APClusteringUnix64(double** distanceMatrix, int groupSize, int &APclustersize, 
						int* &APclustersubrep, int* &APclustersubsize, int** &APcluster);
void APClusteringUnix64WithLibName(char* APLibName, double** distanceMatrix, int groupSize, int &APclustersize, 
						int* &APclustersubrep, int* &APclustersubsize, int** &APcluster);
*/
void MakeFileName(char* filename, char* cur_id, char* padding, int cur_index, char* filetype);
void MakeFileName(char* filename, char* cur_id, char* padding, char* index_string, char* filetype);
void MakeFileName(char* filename, char* id1, char* padding1, char* id2, char* padding2, int cur_index, char* filetype);
void MakeFileName(char* filename, char* id1, char* padding1, char* id2, char* padding2, char* index_string, char* filetype);

void MyPause();


// do affine registration 20100524
void AffineRegistration(char* fixedImageFileName, char* movingImageFileName, char* transformedFileName, char* affineTransformFileName);
void ComposeAffineAndDense(char* affineTransformationFileName, char* nonRigidDeformationFieldFileName, char* finalDeformationFieldFileName);
void PaddingImage(char* originalImageFileName, char* paddedImageFileName, int sizex, int sizey, int sizez);
void ResampleImageWithSizeSpacing(char* originalImageFileName, char* resampledImageFileName, int sizex, int sizey, int sizez, float spacingx, float spacingy, float spacingz);
