/*=========================================================================
  Copyright (c) Hongjun Jia

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

20101122	Create:  Multi-Atlas-Multi-Sample Label Fusion method 
=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
// #pragma warning ( disable : 4244 )
#endif
// for math
#include <vcl_iostream.h>
//#include <vul/vul_timer.h>
#include <vnl/vnl_random.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matlab_print.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_svd_economy.h>

// To include all related header files
#include "itkJia3D.h"
#define LABEL_MAX 255

// global bool variables to adjust the  procedure
bool isRunReg = true; //true;  // do registration to get pairwise deformation field
bool isConverged = false;  // the iterative update can stop
bool isEvaluate = false;  // if false, we do not know the ground-truth of labels
bool isDebug = false; //false;//true; // if true, print out more information
bool isCompressed = false;//true; //false; // if true, compress deformation field by *1000 and

//int numEigenVector = 2; // t
//int numSampleEachDirection = 4; // n  n^t total number of intermediate templates
int numEigenVector = 4; // t
int numSampleEachDirection = 4; // n  n^t total number of intermediate templates

// global variables
int atlas_size = 0;
unsigned int filenumber = 0;
int simulate_size = 0;
int sample_size = 0;
int allfilenumber = 0;
int allatlasnumber = 0;

int localPatchSize = 1; //(2r+1)*(2r+1)*(2r+1) is the volume of local patch
int iter = 1; // iter=1: out 001 using only atlases
              // iter=2: out 002 using all available images including atlases and other labeled images
int imx, imy, imz;// input image sizes

// demons registration parameters
int iterInResolutions[4][3]={{5,3,2},{10,5,5},{15,10,5},{20,15,10}};
int itereach = 2; // 
int itereach0 = 0;int itereach1 = 1;int itereach2 = 2;int itereach3 = 3;
double sigmaDef = 1.5;
double sigmaDef10 = 1.0;double sigmaDef15 = 1.5;double sigmaDef20 = 2.0;
double sigmaDef25 = 2.5;double sigmaDef30 = 3.0;double sigmaDef35 = 3.5;
bool doHistMatch = true;
bool doHistMatchTrue = true;bool doHistMatchFalse = false;
// 
DeformationFieldType::SpacingType df_spacing;
DeformationFieldType::DirectionType df_direction;
DeformationFieldType::PointType	df_origin;

typedef itk::Vector< ShortPixelType, ImageDimension >		ShortVectorPixelType;
typedef itk::Image<  ShortVectorPixelType, ImageDimension >	ShortDeformationFieldType;
typedef itk::ImageFileWriter< ShortDeformationFieldType >	ShortDeformationFieldWriterType;
typedef itk::CastImageFilter< DeformationFieldType, ShortDeformationFieldType >	Internal2ShortDeformationFieldCastFilterType;
typedef itk::CastImageFilter< ShortDeformationFieldType, DeformationFieldType >	Short2InternalDeformationFieldCastFilterType;

// return code: 
/* 101 -> no id files; 
   102 -> atlas_size <= 0, or atlas_size>=filenumber;
   103 -> some files do not exist; 
   104 -> not the same size or resolution;
   105 -> new built tree invalide
   */
using namespace std;

void DiffeoDemonsRegistrationWithParameters(char* fixedImageFileName, char* movingImageFileName, 
											char* deformedImageFileName, char* deformationFieldFileName, 
											double sigmaDef, bool doHistMatch, int* iterInResolutions);
void DiffeoDemonsRegistrationWithInitialWithParameters(char* fixedImageFileName, char* movingImageFileName, 
													   char* initDeformationFieldFileName, char* deformedImageFileName, char* deformationFieldFileName, 
													   double sigmaDef, bool doHistMatch, int* iterInResolutions);
void LabelFusion(char* curSampleImgName, char* outSampleSegName, char** allWarpedAtlasImgNames,
				 char** allDeformationFieldNames, char** allAtlasSegNames, int numOfAtlases);
void GaussianWeightedLabelFusion(InternalImageType::Pointer curSampleImgPtr, InternalImageType::Pointer outSampleSegPtr,
								 InternalImageType::Pointer* warpedImgPtrs, InternalImageType::Pointer* warpedSegPtrs,
								 int numOfAtlases);
void FindCenterBasedOnShortestTotalDistance(double** distanceMatrix, int matrixSize, int &center);
void FromEdgeSegment2Tree(int* tree, bool* treeEdgeUsed, int* v1, int* v2, int curRoot, int matrixSize);
void generateMSTFromMatrix( double** curDistance, int nnode, int* treeMST );
void generateMSTFromMatrixWithFixedRoot( double** curDistance, int nnode, int* treeMST, int root );
void TreeBasedRegistrationFast(int* tree, int tree_size, char** originalIntensityImageFileNames,
							   char** deformedImgFileNames, char** deformationFieldFileNames);
void TreeBasedRegistrationFastOniTree(int* itree_t, int itree_size_t, int atlas_size_t, int simulate_size_t, int sample_size_t,
									  char** originalIntensityImageFileNames,
									  char** deformedImgFileNames,
									  char** deformationFieldFileNames);
bool ValidateTree(int* tree, int tree_size);
void RemoveFile(char* filename);
void RemoveFile(char* filename, bool verbose);
void RemoveMoreFiles(char* filenames, bool verbose);
void DownResampleDeformationField(char* deformationFieldFileName, char* resampledDeformationFieldFileName, int sampleRate);
void DownResampleDeformationFieldNotValue(char* deformationFieldFileName, char* resampledDeformationFieldFileName, int sampleRate);
void UpResampleDeformationField(char* deformationFieldFileName, char* resampledDeformationFieldFileName, int sampleRate);
void UpResampleDeformationFieldNotValue(char* deformationFieldFileName, char* resampledDeformationFieldFileName, int sampleRate);
void LoadIntoArray(char* resampledDeformationFieldFileName, float* df_vector);
void SaveFromArray(char* deformationFieldFileName, float* df_vector, int sx, int sy, int sz);
void DoPCATraining(char** deformationFieldFileNames, int numFiles, char** allImgFileName, int root);
void DeCompressDeformationField(char* filename);
void CompressDeformationField(char* filename);
void DeCompressDeformationFieldFromShort(char* filename, DeformationFieldType::Pointer &deformationfield);
void CompressDeformationField2Short(char* filename);
void ComposeDeformationFieldsAndSaveCompressed(char* inputDeformationFieldFileName,
											   char* deformationFieldFileName,
											   char* composedDeformationFieldFileName);
void ApplyDeformationFieldAndWriteWithTypeWithFileNames(char* movingImageFileName, 
														char* deformationFieldFileName, char* deformedImageFileName, bool isLinear);

int main( int argc, char *argv[] )
{
	if (argc != 3)
	{
		std::cerr << "Missing Parameters of MultiAtlasBasedMultiImageSegmentation: " << std::endl;
		std::cerr << "Usage: MultiAtlasBasedMultiImageSegmentation all_ids.txt number_of_atlas" << std::endl;
		std::cerr << "All parameters are needed as defined below:" << std::endl;
		std::cerr << " all_ids.txt          : file name for all ids of atlases and new samples, e.g., " << std::endl;
		std::cerr << "                        id:AA, image file name:AA_cbq_000.hdr(img), segmentation file name (AA_seg_000.hdr(img)!" << std::endl;
		std::cerr << " number_of_atlas      : the number of atlases starting from the beginning of id list." << std::endl;
		//std::cerr << " fusion_type[0]		: fusion type, 0 for unweighted, 1 for weighted" << std::endl;
		//std::cerr << " neighbor_size[0]		: neighborhood region size, from -n to +n" << std::endl;

		std::cerr << "==========More descriptions about ...=======" << std::endl;
		std::cerr << "  " << std::endl;
		return EXIT_FAILURE;
	}
	atlas_size = atoi(argv[2]);
	///////////////////////////////
	// count all ids and assign names to all intensity and label images (atlases and samples)
	
	if (!itksys::SystemTools::FileExists(argv[1], true))
	{
		std::cerr << "No id list file: " << argv[1] << "!" << std::endl;
		return 101;//EXIT_FAILURE;
	}

	std::cerr << "Reading all ids ... ";
	std::ifstream idFile;
	idFile.open (argv[1]); //ios::binary for opening in binary mode
	idFile.seekg (0, ios::end);// get length of file:
	filenumber = (unsigned int)(idFile.tellg()/3 + 0.5);
	idFile.seekg (0, ios::beg);
	idFile.close();

	if ((atlas_size<=0) || (atlas_size>=filenumber))
	{
		std::cerr << "atlas_size should be within [1," << filenumber-1 << "]!" << std::endl;
		return 102; //EXIT_FAILURE;
	}

	sample_size = filenumber - atlas_size;
	simulate_size = 2*atlas_size; 
	//simulate_size = 2*(atlas_size-1);
	int numAllCombinations = (int)pow((float)numSampleEachDirection, numEigenVector);
	if (numAllCombinations < simulate_size)
	{		
		// std::cerr << "Please increase number of simulated images to be larger than " << numAllCombinations << "!" << std::endl;
		// return 109;//EXIT_FAILURE;
		simulate_size = numAllCombinations;
	}
	allatlasnumber = atlas_size + simulate_size;
	allfilenumber = atlas_size + sample_size + simulate_size;


	idFile.open (argv[1]);
	char** sub_ids = new char*[filenumber];
	char** imgfilenames = new char*[filenumber];
	char** segfilenames = new char*[filenumber];

	for (int i = 0; i< filenumber; i++)
	{
		//char temp_char;
		sub_ids[i] = new char[3];
		idFile >> sub_ids[i][0];
		idFile >> sub_ids[i][1];
		//idFile >> temp_char;
		sub_ids[i][2] = '\0';

		imgfilenames[i] = new char[MAX_FILE_NAME_LENGTH];
		segfilenames[i] = new char[MAX_FILE_NAME_LENGTH];
		MakeFileName(imgfilenames[i], sub_ids[i], "_cbq_", 0, "hdr");
		MakeFileName(segfilenames[i], sub_ids[i], "_seg_", 0, "hdr");
	}
	idFile.close();
	std::cerr << " Done! " << std::endl;

	// count all ids and assign names to all intensity and label images (atlases and samples)
	//////////////////////////////////////
	// do sanity check
	std::cerr << "Checking all files ... ";

	bool allFilesExist = true;
	for (int i = 0; i < filenumber; i++)
	{
		char* tempfileimg_img = new char[MAX_FILE_NAME_LENGTH];
		MakeFileName(tempfileimg_img, sub_ids[i], "_cbq_", 0, "img");
		if (!itksys::SystemTools::FileExists(imgfilenames[i], true))
			allFilesExist = false;
		if (!itksys::SystemTools::FileExists(tempfileimg_img, true))
			allFilesExist = false;
		delete[] tempfileimg_img;
	}
	if (isEvaluate)
		for (int i = 0; i < filenumber; i++) // for real applications: filenumber->atlas_size
		{
			char* tempfileseg_img = new char[MAX_FILE_NAME_LENGTH];
			MakeFileName(tempfileseg_img, sub_ids[i], "_seg_", 0, "img");
			if (!itksys::SystemTools::FileExists(segfilenames[i], true))
				allFilesExist = false;
			if (!itksys::SystemTools::FileExists(tempfileseg_img, true))
				allFilesExist = false;
			delete[] tempfileseg_img;
		}
	else
		for (int i = 0; i < atlas_size; i++) // for real applications: filenumber->atlas_size
		{
			char* tempfileseg_img = new char[MAX_FILE_NAME_LENGTH];
			MakeFileName(tempfileseg_img, sub_ids[i], "_seg_", 0, "img");
			if (!itksys::SystemTools::FileExists(segfilenames[i], true))
				allFilesExist = false;
			if (!itksys::SystemTools::FileExists(tempfileseg_img, true))
				allFilesExist = false;
			delete[] tempfileseg_img;
		}

	std::cerr << " Done! " << std::endl;

	if (!allFilesExist)
	{
		std::cerr << "Please check if all .hdr and .img file exist!" << std::endl;
		for (int i = 0; i < filenumber; i++) delete[] sub_ids[i];		delete[] sub_ids;
		for (int i = 0; i < filenumber; i++) delete[] imgfilenames[i];	delete[] imgfilenames;
		for (int i = 0; i < filenumber; i++) delete[] segfilenames[i];	delete[] segfilenames;
		return 103; //EXIT_FAILURE;
	}
	else
	{
		std::cerr << "All .hdr files and corresponding .img files exist" << std::endl;
	}
	// in the future, the size or resolution of each image should be checked as well.

	// do sanity check
	////////////////////////////////
	// read the size of images
	{
		InternalImageType::Pointer inputImage1 = 0;
		ReadImage(imgfilenames[0], inputImage1);
		InternalImageType::SizeType input_size = inputImage1->GetLargestPossibleRegion().GetSize();
		imx = input_size[0]; imy = input_size[1]; imz = input_size[2];
	}
	std::cerr << "Find root atlas ... " ;

	// read the size of images
	////////////////////////////////
	// build tree on atlases to find the root
	int root = -1;  // root node of tree
	{
		int tree_size = atlas_size;
		double** distanceMatrix = new double*[tree_size];
		for (int i = 0; i < tree_size; i++) distanceMatrix[i] = new double[tree_size];
		for (int i = 0; i < tree_size; i++)
			for (int j = 0; j < tree_size; j++)
				distanceMatrix[i][j] = 0.0;
		int* tree = new int[tree_size]; // each tree

		// fill distance matrix
		PairwiseDistanceAmongImages(imgfilenames, tree_size, distanceMatrix);
		if (isDebug)
		{
			char distFileName[MAX_FILE_NAME_LENGTH];
			MakeFileName(distFileName, "dist", "_", "all", "_", "atlases","txt");
			SaveMatrix2File(distanceMatrix, tree_size, tree_size, distFileName);
		}

		// build MST
		generateMSTFromMatrix(distanceMatrix, tree_size, tree);
		FindRoot(tree, tree_size, root);
		// save distance matrix and tree
		if (isDebug)
		{
			char treeFileName[MAX_FILE_NAME_LENGTH];
			MakeFileName(treeFileName, "tree", "_", "all", "_", "atlases","txt");
			SaveTreeWithInfo(tree, tree_size, treeFileName);
		}
		for (int i = 0; i< tree_size; i++) delete[] distanceMatrix[i]; delete[] distanceMatrix;
		delete[] tree;
	}
	std::cerr << "The root is " << root << ":" << sub_ids[root] << "! ";
	std::cerr << "Done!" << std::endl;
	// build tree on atlases to find the root
	//////////////////////////////
	// registration between root and other atlases
	std::cerr << "Run coarse registration ... ";
	for (int i = 0; i < atlas_size; i++)
	{
		std::cerr << i << ", ";
		if (i == root) continue;
		// if not the root image do registration
		{
			//prepare file names
			char fixedImageFileName[MAX_FILE_NAME_LENGTH];
			char movingImageFileName[MAX_FILE_NAME_LENGTH];
			char deformedImageFileName[MAX_FILE_NAME_LENGTH];
			char deformationFieldFileName[MAX_FILE_NAME_LENGTH];
			MakeFileName(fixedImageFileName, sub_ids[root], "_cbq_", 0, "hdr");
			MakeFileName(movingImageFileName, sub_ids[i], "_cbq_", 0, "hdr");
			MakeFileName(deformedImageFileName, sub_ids[root], "_", sub_ids[i], "_cbq_", "reg", "hdr");
			MakeFileName(deformationFieldFileName, sub_ids[i], "_deform_", 0, "mha");
			
			if (!itksys::SystemTools::FileExists(deformationFieldFileName, true))
			{
				// rough registration
				DiffeoDemonsRegistrationWithParameters(fixedImageFileName, movingImageFileName, 
					deformedImageFileName, deformationFieldFileName, 
					sigmaDef20, doHistMatchTrue, iterInResolutions[itereach0]);
				
				// compress it	
				if (isCompressed)
					CompressDeformationField2Short(deformationFieldFileName);
			}


			RemoveFile(deformedImageFileName);
			char deformedImageFileNameImg[MAX_FILE_NAME_LENGTH];
			MakeFileName(deformedImageFileNameImg, sub_ids[root], "_", sub_ids[i], "_cbq_", "reg", "img");
			RemoveFile(deformedImageFileNameImg);
		}
	}
	std::cerr << "Done!" << std::endl;
	// registration between root and other atlases
	///////////////////////////////
	// do PCA simulation
	std::cerr << "Build deformation field model ... ";
	{
		char** allDeformationFieldFileNames = new char*[atlas_size-1];
		for (int i = 0; i < atlas_size-1; i++)
			allDeformationFieldFileNames[i] = new char[MAX_FILE_NAME_LENGTH];
		int index = 0;
		for (int i = 0; i < atlas_size; i++)
		{
			if (i == root) continue;
			MakeFileName(allDeformationFieldFileNames[index], sub_ids[i], "_deform_", 0, "mha");
			index++;
		}
		DoPCATraining(allDeformationFieldFileNames, atlas_size-1, imgfilenames, root);
		for (int i = 0; i < atlas_size-1; i++) delete[] allDeformationFieldFileNames[i];
		delete[] allDeformationFieldFileNames;
	}
	// do PCA simulation
	/////////////////////////
	// generate simulated template with deformation field
	std::cerr << "Generate simulated templates ... ";
	char** simulateDeformationFieldFileNames = new char*[simulate_size];
	char** simulateTemplateFileNames = new char*[simulate_size];
	{
		// 
		for (int i = 0; i < simulate_size; i++)
		{
			std::cerr << i << ", ";
			simulateDeformationFieldFileNames[i] = new char[MAX_FILE_NAME_LENGTH];
			simulateTemplateFileNames[i] = new char[MAX_FILE_NAME_LENGTH];
			MakeFileName(simulateDeformationFieldFileNames[i], "simulated", "_deform_", i, "mha");
			MakeFileName(simulateTemplateFileNames[i], "simulated", "_cbq_", i, "hdr");
			//load simulated deformation field
			DeformationFieldType::Pointer df = 0;
			if (!isCompressed)
				ReadDeformationField(simulateDeformationFieldFileNames[i], df);
			else
				DeCompressDeformationFieldFromShort(simulateDeformationFieldFileNames[i], df);

			DeformationFieldType::Pointer invdf = DeformationFieldType::New();
			invdf -> SetRegions (df->GetLargestPossibleRegion());
			invdf -> SetSpacing (df->GetSpacing());
			invdf -> SetDirection (df->GetDirection());
			invdf -> SetOrigin (df->GetOrigin());
			invdf -> Allocate();

			InverseDeformationFieldDG3D(df, invdf);
			InternalImageType::Pointer rootImg = 0;
			ReadImage(imgfilenames[root], rootImg);
			InternalImageType::Pointer simImg = 0;
			ApplyDeformationField(rootImg, invdf, simImg, true);
			WriteImage(simulateTemplateFileNames[i], simImg);
		}

	}
	// delete model related deformation field
	if (!isDebug)
	{
		for (int i = 0; i < atlas_size; i++)
		{
			if (i == root) continue;
			{
				char deformationFieldFileName[MAX_FILE_NAME_LENGTH];
				MakeFileName(deformationFieldFileName, sub_ids[i], "_deform_", 0, "mha");
				RemoveFile(deformationFieldFileName);

				char subdeformationFieldFileName[MAX_FILE_NAME_LENGTH];
				strcpy(subdeformationFieldFileName, deformationFieldFileName);
				subdeformationFieldFileName[strlen(subdeformationFieldFileName)-4] = 0;
				strcat(subdeformationFieldFileName,"_sub.mha");
				RemoveFile(subdeformationFieldFileName);
			}
		}
	}
	std::cerr << "Done!" << std::endl;

	// generate simulated template with deformation field
	//////////////////////////////
	// build the combinative tree
	std::cerr << "Build combinative tree ... ";

	int tree_size = allatlasnumber;
	double** distanceMatrix = new double*[tree_size];
	for (int i = 0; i < tree_size; i++) distanceMatrix[i] = new double[tree_size];
	for (int i = 0; i < tree_size; i++)
		for (int j = 0; j < tree_size; j++)
			distanceMatrix[i][j] = 0.0;
	int* tree = new int[tree_size]; // each tree

	// fill distance matrix
	char** atlasfilenames = new char*[tree_size];
	for (int i = 0; i < tree_size; i++)
	{
		atlasfilenames[i] = new char[MAX_FILE_NAME_LENGTH];
		if (i < atlas_size)
			strcpy(atlasfilenames[i],imgfilenames[i]);
		else
			strcpy(atlasfilenames[i], simulateTemplateFileNames[i-atlas_size]);
	}
	PairwiseDistanceAmongImages(atlasfilenames, tree_size, distanceMatrix);

	// build MST
	generateMSTFromMatrixWithFixedRoot(distanceMatrix, tree_size, tree, root);

	// save distance matrix and tree
	if (isDebug)
	{
		char distFileName[MAX_FILE_NAME_LENGTH];
		MakeFileName(distFileName, "dist", "_", "all", "_", "atlases_and_simulated","txt");
		SaveMatrix2File(distanceMatrix, tree_size, tree_size, distFileName);
		char treeFileName[MAX_FILE_NAME_LENGTH];
		MakeFileName(treeFileName, "tree", "_", "all", "_", "atlases_and_simulated","txt");
		SaveTreeWithInfo(tree, tree_size, treeFileName);
	}

	// build the combinative tree
	/////////////////////////////////////
	// copy current tree into an incremental tree

    int itree_size = allfilenumber;
	int* itree = new int[itree_size];
	for (int i = 0; i < itree_size; i++) itree[i] = 0;
	for (int i = 0; i < tree_size; i++) itree[i] = tree[i];

	// copy current tree into an incremental tree
	//////////////////////////////
	// build the incremental tree
	char** allfilenames = new char*[allfilenumber];
	for (int i = 0; i < allfilenumber; i++)
	{
		allfilenames[i] = new char[MAX_FILE_NAME_LENGTH];
		if (i < atlas_size)
			strcpy(allfilenames[i],imgfilenames[i]);
		else if (i < allatlasnumber)
			strcpy(allfilenames[i], simulateTemplateFileNames[i-atlas_size]);
		else
			strcpy(allfilenames[i], imgfilenames[i-allatlasnumber+atlas_size]);
	}

	//std::cerr << "Calculating pairwise distance cross training and test images ...";
	double** cross_dist = new double* [sample_size];
	for (int i = 0; i < sample_size; i++)
	{
		cross_dist[i] = new double[allfilenumber];
		for (int j = 0; j < allfilenumber; j++)
			cross_dist[i][j] = calculateDistanceMSD(imgfilenames[i+atlas_size], allfilenames[j]);
	}
	if (isDebug)
		SaveMatrix2File(cross_dist, sample_size, allfilenumber, "sample2atlas_dist.txt");

	int itree_size_cur = tree_size;
	int sample_left_size = allfilenumber - itree_size_cur;

	int* atlas_cur = new int[allfilenumber];
	for (int i = 0; i < allfilenumber; i++) atlas_cur[i] = 0;
	for (int i = 0; i < allatlasnumber; i++) atlas_cur[i] = i;
	int* sample_cur = new int[sample_size];
	for (int i = 0; i < sample_size; i++) sample_cur[i] = i + allatlasnumber;

	for (int ii = 0; ii < sample_size; ii++)
	{
		//// create a (sample_left_size)by(itree_size_cur) matrix
		//int** dist_matrix_cur = new int*[sample_left_size];
		//for (int i = 0; i < sample_left_size; i++) dist_matrix_cur[i] = new int[itree_size_cur];
		//for (int i = 0; i < sample_left_size; i++)
		//{
		//	for (int j = 0; j < itree_size_cur; j++)
		//	{
		//		dist_matrix_cur[i][j] = cross_dist[sample_cur[i]-atlas_size][atlas_cur[j]];
		//	}
		//}
		// find the min for each row
		double* test2train_min_dist = new double[sample_left_size];
		int* test2train_min_index = new int[sample_left_size];
		{
			int* index_train = new int[itree_size_cur];
			double* dist_temp1 = new double[itree_size_cur];

			for (int j = 0; j < sample_left_size; j++)
			{
				for (int i = 0; i < itree_size_cur; i++)
				{
					index_train[i] = i;
					dist_temp1[i] = cross_dist[sample_cur[j]-allatlasnumber][atlas_cur[i]];
				}
				bubbleSort(dist_temp1, index_train, itree_size_cur);
				test2train_min_index[j] = atlas_cur[index_train[0]];
				test2train_min_dist[j] = dist_temp1[0];
			}
			delete[] dist_temp1;
			delete[] index_train;
		}
		// find the sample index to be attached
		int* index_test = new int[sample_left_size];
		double* dist_temp2 = new double[sample_left_size];
		for (int j = 0; j < sample_left_size; j++) 
		{
			index_test[j] = j;
			dist_temp2[j] = test2train_min_dist[j];
		}
		bubbleSort(dist_temp2, index_test, sample_left_size);

		// attach it to the tree 
		int sample_to_attach = sample_cur[index_test[0]];
		int best_match_atlas = test2train_min_index[index_test[0]];// parent node
		itree[sample_to_attach] = best_match_atlas;

		// include this sample to the list of current atlas
		atlas_cur[itree_size_cur] = sample_to_attach;
		// switch this sample to  the end of the list of sample_cur
		for (int i = index_test[0]; i < sample_size; i++) sample_cur[i] = sample_cur[i+1];
		sample_cur[sample_size-1] = sample_to_attach;

		delete[] dist_temp2;
		delete[] index_test;

		// the tree increases and sample left decreases
		itree_size_cur++;
		sample_left_size--;
		//delete
		//for (int i = 0; i < sample_left_size; i++) delete[] dist_matrix_cur[i];
		//delete[] dist_matrix_cur;
		delete[] test2train_min_dist;
		delete[] test2train_min_index;
	}
	// delete
	delete[] atlas_cur;
	delete[] sample_cur;
	for (int i = 0; i < sample_size; i++) delete[] cross_dist[i];	delete[] cross_dist;

	// check if new tree is actually a tree
	bool iTreeValid = false;
	iTreeValid = ValidateTree(itree, itree_size);
	if (!iTreeValid)
	{
		std::cerr << "The new built tree is NOT valid!" << std::endl;
		for (int i = 0; i < filenumber; i++) delete[] sub_ids[i];	delete[] sub_ids;
		for (int i = 0; i < filenumber; i++) delete[] imgfilenames[i];	delete[] imgfilenames;
		for (int i = 0; i < filenumber; i++) delete[] segfilenames[i];	delete[] segfilenames;
		for (int i = 0; i< tree_size; i++) delete[] distanceMatrix[i]; delete[] distanceMatrix;
		delete[] tree;
		delete[] itree;
		for (int i = 0; i < simulate_size; i++) delete[] simulateDeformationFieldFileNames[i];	delete[] simulateDeformationFieldFileNames;
		for (int i = 0; i < simulate_size; i++) delete[] simulateTemplateFileNames[i];	delete[] simulateTemplateFileNames;
		for (int i = 0; i < allatlasnumber; i++) delete[] atlasfilenames[i]; delete[] atlasfilenames;
		for (int i = 0; i < allfilenumber; i++) delete[] allfilenames[i]; delete[] allfilenames;

		return 105;// EXIT_FAILURE;
	}
	// save new tree
	if (isDebug)
	{
		SaveTreeWithInfo(itree, itree_size, "itree.txt");
		// std::cerr << "Tree is built up!" << std::endl;
	}
	std::cerr << "Done!" << std::endl;

	// build the incremental tree
	////////////////////////////////////
	// do registration on the tree
	// 20110505
	std::cerr << "Start registration ... ";
	{
		// registration to the root
		// prepare file names
		char** originalIntensityImageFileNames = new char*[filenumber];
		char** deformedImgFileNames = new char*[filenumber];
		char** deformedImgFileNamesImg = new char*[filenumber];
		char** deformedSegFileNames = new char*[filenumber];
		char** deformedSegFileNamesImg = new char*[filenumber];
		char** deformationFieldFileNames = new char*[filenumber];

		for (int i = 0; i < filenumber; i++)
		{
			originalIntensityImageFileNames[i] = new char[MAX_FILE_NAME_LENGTH];
			MakeFileName(originalIntensityImageFileNames[i], sub_ids[i],"_cbq_",0,"hdr");

			deformedImgFileNames[i] = new char[MAX_FILE_NAME_LENGTH];
			deformedImgFileNamesImg[i] = new char[MAX_FILE_NAME_LENGTH];
			deformedSegFileNames[i] = new char[MAX_FILE_NAME_LENGTH];
			deformedSegFileNamesImg[i] = new char[MAX_FILE_NAME_LENGTH];
			deformationFieldFileNames[i] = new char[MAX_FILE_NAME_LENGTH];

			MakeFileName(deformationFieldFileNames[i], sub_ids[root], "_", sub_ids[i], "_deform_", 0, "mha");
			MakeFileName(deformedImgFileNames[i], sub_ids[root], "_", sub_ids[i], "_cbq_", 0, "hdr");
			MakeFileName(deformedImgFileNamesImg[i], sub_ids[root], "_", sub_ids[i], "_cbq_", 0, "img");
			MakeFileName(deformedSegFileNames[i], sub_ids[root], "_", sub_ids[i], "_seg_", 0, "hdr");
			MakeFileName(deformedSegFileNamesImg[i], sub_ids[root], "_", sub_ids[i], "_seg_", 0, "img");

		}
		// do tree-based registration
		//TreeBasedRegistrationFast(itree, itree_size, originalIntensityImageFileNames,
		//	deformedImgFileNames, deformationFieldFileNames);
		TreeBasedRegistrationFastOniTree(itree, itree_size, atlas_size, simulate_size, sample_size,
			originalIntensityImageFileNames, deformedImgFileNames, deformationFieldFileNames);

		// warp the corresponding seg file
		if (isEvaluate)
		{
			for (int i = 0; i < filenumber; i++)
			{
				if (i ==  root)
					continue;
				//
				bool isWarpedSegExist = true;
				if (!itksys::SystemTools::FileExists(deformedSegFileNames[i], true)) 
					isWarpedSegExist = false;
				if (!itksys::SystemTools::FileExists(deformedSegFileNamesImg[i], true)) 
					isWarpedSegExist = false;

				if (!isWarpedSegExist)
				{
					//
					ApplyDeformationFieldAndWriteWithFileNames(segfilenames[i], 
						deformationFieldFileNames[i], deformedSegFileNames[i], false);
				}
				else
				{
					//
				}
			}
		}
		// delete file names
		for (int i = 0; i < filenumber; i++) delete[] originalIntensityImageFileNames[i]; 
		delete[] originalIntensityImageFileNames;
		for (int i = 0; i < filenumber; i++) delete[] deformedImgFileNames[i]; 
		delete[] deformedImgFileNames;
		for (int i = 0; i < filenumber; i++) delete[] deformedImgFileNamesImg[i]; 
		delete[] deformedImgFileNamesImg;
		for (int i = 0; i < filenumber; i++) delete[] deformedSegFileNames[i]; 
		delete[] deformedSegFileNames;
		for (int i = 0; i < filenumber; i++) delete[] deformedSegFileNamesImg[i]; 
		delete[] deformedSegFileNamesImg;
		for (int i = 0; i < filenumber; i++) delete[] deformationFieldFileNames[i]; 
		delete[] deformationFieldFileNames;
	}// end of do registration to the root

	//
	RemoveMoreFiles("simulated_*", false);
	
	{
		std::cerr << "Reverse deformation field ... ";
		// reverse each deformation field (from the root)
		for (int i = 0; i < filenumber; i++)
		{
			//if (isDebug)
			{
				std::cerr << i << ", ";
			}
			if (i == root) continue;
			// prepare file names
			char* originalImgImageFileName = new char[MAX_FILE_NAME_LENGTH];
			char* originalSegImageFileName = new char[MAX_FILE_NAME_LENGTH];
			MakeFileName(originalImgImageFileName, sub_ids[root],"_cbq_",0,"hdr");
			MakeFileName(originalSegImageFileName, sub_ids[root],"_seg_",0,"hdr");

			char* fixedImgImageFileName = new char[MAX_FILE_NAME_LENGTH];
			MakeFileName(fixedImgImageFileName, sub_ids[i],"_cbq_",0,"hdr");

			char* deformedImgFileName = new char[MAX_FILE_NAME_LENGTH];
			char* deformedImgFileNameImg = new char[MAX_FILE_NAME_LENGTH];
			char* deformedSegFileName = new char[MAX_FILE_NAME_LENGTH];
			char* deformedSegFileNameImg = new char[MAX_FILE_NAME_LENGTH];
			char* inversedDeformationFieldFileName = new char[MAX_FILE_NAME_LENGTH];

			MakeFileName(deformedImgFileName, sub_ids[i], "_", sub_ids[root], "_cbq_", 0, "hdr");
			MakeFileName(deformedImgFileNameImg, sub_ids[i], "_", sub_ids[root], "_cbq_", 0, "img");
			MakeFileName(deformedSegFileName, sub_ids[i], "_", sub_ids[root], "_seg_", 0, "hdr");
			MakeFileName(deformedSegFileNameImg, sub_ids[i], "_", sub_ids[root], "_seg_", 0, "img");
			MakeFileName(inversedDeformationFieldFileName, sub_ids[i], "_", sub_ids[root], "_deform_", 0, "mha");
			
			char* deformationFieldFileName = new char[MAX_FILE_NAME_LENGTH];
			MakeFileName(deformationFieldFileName, sub_ids[root], "_", sub_ids[i], "_deform_", 0, "mha");

			// if exist??
			bool isReversedExist = true;
			if (!itksys::SystemTools::FileExists(deformedImgFileName, true)) 	isReversedExist = false;
			if (!itksys::SystemTools::FileExists(deformedImgFileNameImg, true)) isReversedExist = false;
			//if (!itksys::SystemTools::FileExists(deformedSegFileName, true)) 	isReversedExist = false;
			//if (!itksys::SystemTools::FileExists(deformedSegFileNameImg, true)) isReversedExist = false;
			if (!itksys::SystemTools::FileExists(inversedDeformationFieldFileName, true)) 	isReversedExist = false;
			if (isReversedExist) 
			{
				//CompressDeformationField2Short(inversedDeformationFieldFileName);
				continue;
			}

			// reverse deformation field
			DeformationFieldType::Pointer deformationField = 0;
			if (!isCompressed)
				ReadDeformationField(deformationFieldFileName, deformationField);
			else
				DeCompressDeformationFieldFromShort(deformationFieldFileName, deformationField);
			DeformationFieldType::Pointer inversedDeformationField = DeformationFieldType::New();
			inversedDeformationField -> SetRegions (deformationField->GetLargestPossibleRegion());
			inversedDeformationField -> SetSpacing (deformationField->GetSpacing());
			inversedDeformationField -> SetDirection (deformationField->GetDirection());
			inversedDeformationField -> SetOrigin (deformationField->GetOrigin());
			inversedDeformationField -> Allocate();
			InverseDeformationFieldDG3D(deformationField, inversedDeformationField);
			// output reversed deformation field
			WriteDeformationField(inversedDeformationFieldFileName, inversedDeformationField);
			if (isCompressed)
				CompressDeformationField2Short(inversedDeformationFieldFileName);
			// //update
			DiffeoDemonsRegistrationWithInitialWithParameters(fixedImgImageFileName, originalImgImageFileName, 
				inversedDeformationFieldFileName, deformedImgFileName, inversedDeformationFieldFileName, 
				sigmaDef15, doHistMatch, iterInResolutions[itereach0]);
			// CompressDeformationField2Short(inversedDeformationFieldFileName);
			// apply deformation field on seg file
			if (isEvaluate)
			{
				ApplyDeformationFieldAndWriteWithFileNames(originalSegImageFileName, 
					inversedDeformationFieldFileName, deformedSegFileName, false);
			}

			// delete
			delete[] originalImgImageFileName;
			delete[] originalSegImageFileName;

			delete[] fixedImgImageFileName;

			delete[] deformedImgFileName;
			delete[] deformedImgFileNameImg;
			delete[] deformedSegFileName;
			delete[] deformedSegFileNameImg;

			delete[] deformationFieldFileName;
			delete[] inversedDeformationFieldFileName;
		}
		std::cout << "Done!" << std::endl;
	}// end of reverse each deformation field (from the root)
	// std::cerr << "done." << std::endl;
	
	// do registration on the tree
	/////////////////////////////////
	// with the help of root node, obtain the pairwise registration among all images
	std::cerr << "Generate all pairwise registration ..."; 
	for (int all_index = 0; all_index < filenumber; all_index++)
	{
		//if (isDebug)
			std::cerr << all_index << ": ";
		if (all_index == root) continue;

		for (int sample_index = atlas_size; sample_index < filenumber; sample_index++)
		{
			//if (isDebug)
				std::cerr << "[" << all_index << "," << sample_index << "], ";
			
			if (sample_index == all_index) continue;

			{
				// warp all_index(moving) to sample_index(fixed)
				char* fixedImgImageFileName = new char[MAX_FILE_NAME_LENGTH]; // index2
				char* movingImgImageFileName = new char[MAX_FILE_NAME_LENGTH];// index1
				char* movingSegImageFileName = new char[MAX_FILE_NAME_LENGTH];// index1
				MakeFileName(fixedImgImageFileName, sub_ids[sample_index],"_cbq_",0,"hdr");
				MakeFileName(movingImgImageFileName, sub_ids[all_index],"_cbq_",0,"hdr");
				MakeFileName(movingSegImageFileName, sub_ids[all_index],"_seg_",0,"hdr");

				char* deformedImgImageFileName = new char[MAX_FILE_NAME_LENGTH];
				char* deformedSegImageFileName = new char[MAX_FILE_NAME_LENGTH];
				char* deformedImgImageFileNameImg = new char[MAX_FILE_NAME_LENGTH];
				char* deformedSegImageFileNameImg = new char[MAX_FILE_NAME_LENGTH];
				char* deformationFieldFileName = new char[MAX_FILE_NAME_LENGTH];
				MakeFileName(deformedImgImageFileName, sub_ids[sample_index], "_", sub_ids[all_index], "_cbq_", 0, "hdr");
				MakeFileName(deformedSegImageFileName, sub_ids[sample_index], "_", sub_ids[all_index], "_seg_", 0, "hdr");
				MakeFileName(deformedImgImageFileNameImg, sub_ids[sample_index], "_", sub_ids[all_index], "_cbq_", 0, "img");
				MakeFileName(deformedSegImageFileNameImg, sub_ids[sample_index], "_", sub_ids[all_index], "_seg_", 0, "img");
				MakeFileName(deformationFieldFileName, sub_ids[sample_index], "_", sub_ids[all_index], "_deform_", 0, "mha");

				//if exist
				bool isComposeExist = true;
				if (!itksys::SystemTools::FileExists(deformedImgImageFileName, true)) 	isComposeExist = false;
				if (!itksys::SystemTools::FileExists(deformedImgImageFileNameImg, true)) isComposeExist = false;
				//if (!itksys::SystemTools::FileExists(deformedSegImageFileName, true)) 	isComposeExist = false;
				//if (!itksys::SystemTools::FileExists(deformedSegImageFileNameImg, true)) isComposeExist = false;
				if (!itksys::SystemTools::FileExists(deformationFieldFileName, true)) 	isComposeExist = false;
				if (isComposeExist)
				{
					//CompressDeformationField2Short(deformationFieldFileName);
					continue;
				}

				// compose
				char* inputDeformationFieldFileName = new char[MAX_FILE_NAME_LENGTH];
				char* secondDeformationFieldFileName = new char[MAX_FILE_NAME_LENGTH];
				MakeFileName(inputDeformationFieldFileName, sub_ids[root],"_",sub_ids[all_index],"_deform_",0,"mha");
				MakeFileName(secondDeformationFieldFileName, sub_ids[sample_index],"_",sub_ids[root],"_deform_",0,"mha");
				// output composed deformation field
				if (!isCompressed)
					ComposeDeformationFieldsAndSave(inputDeformationFieldFileName,
						secondDeformationFieldFileName, deformationFieldFileName);
				else
					ComposeDeformationFieldsAndSaveCompressed(inputDeformationFieldFileName,
					secondDeformationFieldFileName, deformationFieldFileName);

				//// apply composed deformation field
				//ApplyDeformationFieldAndWriteWithFileNames(movingImgImageFileName, 
				//	deformationFieldFileName, deformedImgImageFileName, true);
				// update
				DiffeoDemonsRegistrationWithInitialWithParameters(fixedImgImageFileName, movingImgImageFileName, 
					deformationFieldFileName, deformedImgImageFileName, deformationFieldFileName, 
					sigmaDef15, doHistMatch, iterInResolutions[itereach0]);
				//CompressDeformationField2Short(deformationFieldFileName);
				if (isEvaluate)
				{
					ApplyDeformationFieldAndWriteWithFileNames(movingSegImageFileName, 
						deformationFieldFileName, deformedSegImageFileName, false);
				}

				//delete
				delete[] fixedImgImageFileName;
				delete[] movingImgImageFileName;
				delete[] movingSegImageFileName;
				delete[] deformedImgImageFileName;
				delete[] deformedSegImageFileName;
				delete[] deformedImgImageFileNameImg;
				delete[] deformedSegImageFileNameImg;
				delete[] deformationFieldFileName;
				delete[] inputDeformationFieldFileName;
				delete[] secondDeformationFieldFileName;
			}
		}// end of sample_index
	}// end of all_index

	std::cout << "Done!" << std::endl;

	// with the help of root node, obtain the pairwise registration among all images
	////////////////////////////////
	// do all multi-atlas based segmentation with weighted label fusion

	std::cerr << "Start to process segmentation..." << std::endl;
	iter = 1;
	while (!isConverged)
	{
		//if (isDebug)
			std::cerr << "iteration "<< iter << ": ";
		// iter == 1, only using atlases; otherwise using all images
		// number of other image to be used
		int numOtherAtlases = 0;
		if (iter ==1 ) numOtherAtlases = atlas_size;
		else numOtherAtlases = filenumber - 1;

		// do label fusion for each new sample
		for (int sample_index = atlas_size;sample_index < filenumber; sample_index++)
		{
			//if (isDebug)
			//	std::cerr << sample_index << ":" << sub_ids[sample_index] << ", ";
			std::cerr << sample_index << ", ";
			// the sample to be processed
			char* curSampleImgName = new char[MAX_FILE_NAME_LENGTH];
			MakeFileName(curSampleImgName, sub_ids[sample_index], "_cbq_", 0, "hdr");
			// the output label image of current sample
			char* outSampleSegName = new char[MAX_FILE_NAME_LENGTH];
			MakeFileName(outSampleSegName, sub_ids[sample_index], "_seg_", iter, "hdr");
			

			// prepare all current atlases names
			char** allWarpedAtlasImgNames = new char*[numOtherAtlases];
			char** allDeformationFieldNames = new char*[numOtherAtlases];
			char** allAtlasSegNames = new char*[numOtherAtlases] ;
			for (int i = 0; i < numOtherAtlases; i++)
			{
				allWarpedAtlasImgNames[i] = new char[MAX_FILE_NAME_LENGTH];
				allDeformationFieldNames[i] = new char[MAX_FILE_NAME_LENGTH];
				allAtlasSegNames[i] = new char[MAX_FILE_NAME_LENGTH];
			}
			int atlas_index = 0;
			for (int i = 0; i < filenumber; i++)
			{
				if ((iter == 1)&&(i >= atlas_size))
					continue;
				if (i == sample_index)
					continue;
				// assign names
				MakeFileName(allWarpedAtlasImgNames[atlas_index], sub_ids[sample_index], "_", sub_ids[i], "_cbq_", 0, "hdr");
				MakeFileName(allDeformationFieldNames[atlas_index], sub_ids[sample_index], "_", sub_ids[i], "_deform_", 0, "mha");
				if (i < atlas_size)
					MakeFileName(allAtlasSegNames[atlas_index], sub_ids[i], "_seg_", 0, "hdr");
				else
					MakeFileName(allAtlasSegNames[atlas_index], sub_ids[i], "_seg_", iter-1, "hdr");
				// go to next atlas
				atlas_index++;
			}
			// do weighted label fusion
				
			LabelFusion(curSampleImgName, outSampleSegName, allWarpedAtlasImgNames,	allDeformationFieldNames, allAtlasSegNames, numOtherAtlases);

			// delete newed variables
			for (int i = 0; i < numOtherAtlases; i++) delete[] allWarpedAtlasImgNames[i];	delete[] allWarpedAtlasImgNames;
			for (int i = 0; i < numOtherAtlases; i++) delete[] allDeformationFieldNames[i];	delete[] allDeformationFieldNames;
			for (int i = 0; i < numOtherAtlases; i++) delete[] allAtlasSegNames[i];	delete[] allAtlasSegNames;
			delete[] curSampleImgName;
			delete[] outSampleSegName;
		}

		// to decide if converged
		if (iter >= 3)
			isConverged = true; // 

		// calculate overlap rate
		if (isEvaluate)
		{
			// pairwise segmentation accuracy

			// pairwise segmentation consistency
		}
		//if (isDebug)	
			std::cerr << std::endl;
		iter++;
	}

	std::cerr << "Done! "<< std::endl;
	// do all multi-atlas based segmentation with weighted label fusion
	/////////////////////////////////
	// 
	RemoveMoreFiles("*_*_cbq_000*", false);
	RemoveMoreFiles("*_*_deform_000*", false);
	// 
	/////////////////////////////////
	//

	for (int i = 0; i < filenumber; i++) delete[] sub_ids[i];	delete[] sub_ids;
	for (int i = 0; i < filenumber; i++) delete[] imgfilenames[i];	delete[] imgfilenames;
	for (int i = 0; i < filenumber; i++) delete[] segfilenames[i];	delete[] segfilenames;
	for (int i = 0; i< tree_size; i++) delete[] distanceMatrix[i]; delete[] distanceMatrix;
	delete[] tree;
	delete[] itree;
	for (int i = 0; i < simulate_size; i++) delete[] simulateDeformationFieldFileNames[i];	delete[] simulateDeformationFieldFileNames;
	for (int i = 0; i < simulate_size; i++) delete[] simulateTemplateFileNames[i];	delete[] simulateTemplateFileNames;
	for (int i = 0; i < tree_size; i++) delete[] atlasfilenames[i]; delete[] atlasfilenames;
	for (int i = 0; i < allfilenumber; i++) delete[] allfilenames[i]; delete[] allfilenames;

	return EXIT_SUCCESS;

}
void LabelFusion(char* curSampleImgName, char* outSampleSegName, char** allWarpedAtlasImgNames,
						 char** allDeformationFieldNames, char** allAtlasSegNames, int numOfAtlases)
{
	// create pointers for each images including: cur_image, our_seg, 
	// warpedImg(after histogram matching), warpedSeg (after warping)

	// cur sample image pointer 
	InternalImageType::Pointer curSampleImgPtr = 0;
	ReadImage(curSampleImgName, curSampleImgPtr);

	// output sample label image pointer
	InternalImageType::Pointer outSampleSegPtr = 0;
	ReadImage(curSampleImgName, outSampleSegPtr);

	// atlas related images pointers
	InternalImageType::Pointer* warpedImgPtrs = new InternalImageType::Pointer[numOfAtlases];
	InternalImageType::Pointer* warpedSegPtrs = new InternalImageType::Pointer[numOfAtlases];
	for (int i = 0; i < numOfAtlases; i++)
	{
		// histogram match each warped atlas to the current sample
		warpedImgPtrs[i] = 0;
		InternalImageType::Pointer curWarpedAtlasPtr = 0;
		ReadImage(allWarpedAtlasImgNames[i], curWarpedAtlasPtr);
		HistogramMatching(curWarpedAtlasPtr, curSampleImgPtr, warpedImgPtrs[i]);

		// load each deformation field
		DeformationFieldType::Pointer deformFieldPtr = 0;
		if (!isCompressed)
			ReadDeformationField(allDeformationFieldNames[i], deformFieldPtr);
		else
			DeCompressDeformationFieldFromShort(allDeformationFieldNames[i], deformFieldPtr);

		// warp each atlas label image to the current sample space
		InternalImageType::Pointer curAtlasSegPtr = 0;
		ReadImage(allAtlasSegNames[i], curAtlasSegPtr);
		ApplyDeformationField(curAtlasSegPtr, deformFieldPtr, warpedSegPtrs[i], false);
	}

	// weighted label fusion
	GaussianWeightedLabelFusion(curSampleImgPtr, outSampleSegPtr, warpedImgPtrs, warpedSegPtrs, numOfAtlases);

	// output segmentation image
	WriteImage(outSampleSegName, outSampleSegPtr);

	delete[] warpedImgPtrs;
	delete[] warpedSegPtrs;

	return;
}

void GaussianWeightedLabelFusion(InternalImageType::Pointer curSampleImgPtr, InternalImageType::Pointer outSampleSegPtr,
								 InternalImageType::Pointer* warpedImgPtrs, InternalImageType::Pointer* warpedSegPtrs,
								 int numOfAtlases)
{
	// introduce the iterators of each image
	InternalImageIteratorType sampleImgIt(curSampleImgPtr, curSampleImgPtr->GetLargestPossibleRegion());//GetRequestedRegion());
	InternalImageIteratorType sampleSegIt(outSampleSegPtr, outSampleSegPtr->GetLargestPossibleRegion());//GetRequestedRegion());
	
	InternalImageIteratorType* warpedImgIts = new InternalImageIteratorType[numOfAtlases];
	InternalImageIteratorType* warpedSegIts = new InternalImageIteratorType[numOfAtlases];
	for (int i = 0; i < numOfAtlases; i++)
	{
		InternalImageIteratorType it1( warpedImgPtrs[i], warpedImgPtrs[i]->GetLargestPossibleRegion());//GetRequestedRegion() );
		warpedImgIts[i] = it1;
		InternalImageIteratorType it2( warpedSegPtrs[i], warpedSegPtrs[i]->GetLargestPossibleRegion());//GetRequestedRegion() );
		warpedSegIts[i] = it2;
	}

	InternalImageIteratorType::IndexType idx;
	InternalImageIteratorType::IndexType idx1;

	int	temp1, temp2;
	double kernel_sigma = 1; // std of Gaussian approximation, set to the median of all distances 
	int* index = new int[numOfAtlases];
	double* sort1 = new double[numOfAtlases];
	int* label_index = new int[LABEL_MAX+1];

	// for each voxel do label fusion,(vx, vy, vz) is the current location
	for (int vz = 0; vz < imz; vz++)
	{
		//if (isDebug)
		//	std::cerr << vz << ", ";
		for (int vy = 0; vy < imy; vy++)
		{
			for (int vx = 0; vx < imx; vx++)
			{
				// process each voxel
				idx.SetElement(0,vx);
				idx.SetElement(1,vy);
				idx.SetElement(2,vz);

				// get labels from all atlases on this current location
				int* label_pool = new int[numOfAtlases];
				for (int i = 0; i < numOfAtlases; i++)
				{

					warpedSegIts[i].SetIndex(idx);
					label_pool[i] = warpedSegIts[i].Get();
				}

				// get the local patch differences
				double* mse = new double[numOfAtlases];
				for (int i = 0; i < numOfAtlases; i++) mse[i] = 0.0;
				for (int i = 0; i < numOfAtlases; i++)
				{
					double msec = 0.0;
					// compare with all other images at the local patch of current voxel
					// (nx,ny,nz) is the incremental position on (vx,vy,vz)
					for(int nz = 0-localPatchSize; nz <= localPatchSize; nz++)
					{
						for(int ny = 0-localPatchSize; ny <= localPatchSize; ny++)
						{
							for(int nx = 0-localPatchSize; nx <= localPatchSize; nx++)
							{
								if ((vx+nx>=0)&&(vx+nx<imx)&&(vy+ny>=0)&&(vy+ny<imy)&&(vz+nz>=0)&&(vz+nz<imz))
								{
									// within the valid range, and then get voxel intensity
									// 
									idx1.SetElement(0,vx+nx);
									idx1.SetElement(1,vy+ny);
									idx1.SetElement(2,vz+nz);

									sampleImgIt.SetIndex(idx1);
									warpedImgIts[i].SetIndex(idx1);

									temp1 = sampleImgIt.Get();
									temp2 = warpedImgIts[i].Get();
									
									// add together differences, squared sum
									msec += (temp1-temp2)*(temp1-temp2);
								}
							}
						}
					}
					mse[i] = sqrt(msec);
				}// end of for nz

				// sort all difference

				for (int i = 0; i < numOfAtlases; i++)
				{
					index[i] = i;
					sort1[i] = mse[i];
				}
				bubbleSort(sort1, index, numOfAtlases);
				kernel_sigma = sort1[(int)((numOfAtlases-1)/2)]+0.0001;
				
				// weight each label

				double* weight_sum =  new double[LABEL_MAX+1];
				for (int i = 0; i < LABEL_MAX+1; i++) weight_sum[i] = 0.0;
				for (int i = 0; i < numOfAtlases; i++)
				{
					// add-up all weights
					double weight = 0;
					weight = exp(0-(mse[i]*mse[i])/(2*kernel_sigma*kernel_sigma));
					weight_sum[label_pool[i]]+=weight;
				}
				
				// weighted label fusion
				for (int i = 0; i < LABEL_MAX+1; i++)	label_index[i] = i;
				
				bubbleSort(weight_sum, label_index, LABEL_MAX+1);
				int label_final;
				label_final = label_index[LABEL_MAX];

				// assign label to the segmentations
				sampleSegIt.SetIndex(idx);
				sampleSegIt.Set(label_final);
				

				delete[] mse;
				delete[] label_pool;
				delete[] weight_sum;
			}
		}
	}// end of for vz

	//
	delete[] warpedImgIts;
	delete[] warpedSegIts;
	delete[] sort1;
	delete[] index;
	delete[] label_index;

	return;
}


void DiffeoDemonsRegistrationWithParameters(char* fixedImageFileName, char* movingImageFileName, 
											char* deformedImageFileName, char* deformationFieldFileName, 
											double sigmaDef, bool doHistMatch, int* iterInResolutions)
{
	int res = 3;
	// read fixed and moving images
	InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
	InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

	fixedImageReader->SetFileName( fixedImageFileName );
	movingImageReader->SetFileName( movingImageFileName );

	// do histogram matching and some parameters need to set manually
	InternalHistMatchFilterType::Pointer matcher = InternalHistMatchFilterType::New();
	if (doHistMatch)
	{
		matcher->SetInput( movingImageReader->GetOutput() );
		matcher->SetReferenceImage( fixedImageReader->GetOutput() );
		matcher->SetNumberOfHistogramLevels( 1024 ); // 1024
		matcher->SetNumberOfMatchPoints( 7 );
		matcher->ThresholdAtMeanIntensityOn();
	}

	// Set up the demons filter
	typedef itk::PDEDeformableRegistrationFilter
		< InternalImageType, InternalImageType, DeformationFieldType>   BaseRegistrationFilterType;
	BaseRegistrationFilterType::Pointer filter;

	// s <- s o exp(u) (Diffeomorphic demons)
	typedef itk::DiffeomorphicDemonsRegistrationFilter
		< InternalImageType, InternalImageType, DeformationFieldType>
		ActualRegistrationFilterType;
	typedef ActualRegistrationFilterType::GradientType GradientType;

	ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();

	float maxStepLength = 2.0;
	actualfilter->SetMaximumUpdateStepLength( maxStepLength );
	unsigned int gradientType = 0;
	actualfilter->SetUseGradientType(static_cast<GradientType>(gradientType) );
	filter = actualfilter;

	// set up smoothing kernel for deformation field
	//float sigmaDef = 1.5;
	if (sigmaDef > 0.1)
	{
		filter->SmoothDeformationFieldOn();
		filter->SetStandardDeviations( sigmaDef );
	}
	else
	{
		filter->SmoothDeformationFieldOff();
	}

	// set up smoothing kernel for update field
	float sigmaUp = 0.0;
	if ( sigmaUp > 0.1 )
	{
		filter->SmoothUpdateFieldOn();
		filter->SetUpdateFieldStandardDeviations( sigmaUp );
	}
	else
	{
		filter->SmoothUpdateFieldOff();
	}

	//filter->SetIntensityDifferenceThreshold( 0.001 );

	// Set up the multi-resolution filter
	typedef itk::MultiResolutionPDEDeformableRegistration<
		InternalImageType, InternalImageType, DeformationFieldType, InternalPixelType >   MultiResRegistrationFilterType;
	MultiResRegistrationFilterType::Pointer multires = MultiResRegistrationFilterType::New();

	typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<DeformationFieldType,double> FieldInterpolatorType;

	FieldInterpolatorType::Pointer VectorInterpolator = FieldInterpolatorType::New();

#if ( ITK_VERSION_MAJOR > 3 ) || ( ITK_VERSION_MAJOR == 3 && ITK_VERSION_MINOR > 8 )
	multires->GetFieldExpander()->SetInterpolator(VectorInterpolator);
#endif
	std::vector<unsigned int> curNumIterations;
	// unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < res; i++)		curNumIterations.push_back(iterInResolutions[i]);

	multires->SetRegistrationFilter( filter );
	multires->SetNumberOfLevels( curNumIterations.size() );
	multires->SetNumberOfIterations( &curNumIterations[0] );
	multires->SetFixedImage( fixedImageReader->GetOutput() );
	if (doHistMatch)
	{
		multires->SetMovingImage( matcher->GetOutput() );
	}
	else
	{
		multires->SetMovingImage( movingImageReader->GetOutput() );
	}

	// Compute the deformation field
	try
	{
		multires->UpdateLargestPossibleRegion();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at demons registration between input images!" << std::endl; 
		std::cerr << excep << std::endl;
	}

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	ApplyDeformationFieldAndWriteWithTypeWithFileNames(movingImageFileName, 
		deformationFieldFileName, deformedImageFileName, true);

	//// write deformed image
	//InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	//InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	//warper->SetInput( movingImageReader->GetOutput() );
	//warper->SetInterpolator( interpolator );
	//warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	//warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	//warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	//warper->SetDeformationField( multires->GetOutput() );
	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	//if (isCompressed)
	//	CompressDeformationField2Short(deformationFieldFileName);
	return;
}

void DiffeoDemonsRegistrationWithInitialWithParameters(char* fixedImageFileName, char* movingImageFileName, 
													   char* initDeformationFieldFileName, char* deformedImageFileName, char* deformationFieldFileName, 
													   double sigmaDef, bool doHistMatch, int* iterInResolutions)
{
	int res = 3;
	// read initial deformation field file
	DeformationFieldType::Pointer initDeformationField = 0;
	if (!isCompressed)
		ReadDeformationField(initDeformationFieldFileName, initDeformationField);
	else
		DeCompressDeformationFieldFromShort(initDeformationFieldFileName, initDeformationField);
	
	// read fixed and moving images
	InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
	InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

	fixedImageReader->SetFileName( fixedImageFileName );
	movingImageReader->SetFileName( movingImageFileName );

	// do histogram matching and some parameters need to set manually
	InternalHistMatchFilterType::Pointer matcher = InternalHistMatchFilterType::New();
	if (doHistMatch)
	{
		matcher->SetInput( movingImageReader->GetOutput() );
		matcher->SetReferenceImage( fixedImageReader->GetOutput() );
		matcher->SetNumberOfHistogramLevels( 1024 ); // 1024
		matcher->SetNumberOfMatchPoints( 7 );
		matcher->ThresholdAtMeanIntensityOn();
	}

	// Set up the demons filter
	typedef itk::PDEDeformableRegistrationFilter
		< InternalImageType, InternalImageType, DeformationFieldType>   BaseRegistrationFilterType;
	BaseRegistrationFilterType::Pointer filter;

	// s <- s o exp(u) (Diffeomorphic demons)
	typedef itk::DiffeomorphicDemonsRegistrationFilter
		< InternalImageType, InternalImageType, DeformationFieldType>
		ActualRegistrationFilterType;
	typedef ActualRegistrationFilterType::GradientType GradientType;

	ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();

	float maxStepLength = 2.0;
	actualfilter->SetMaximumUpdateStepLength( maxStepLength );
	unsigned int gradientType = 0;
	actualfilter->SetUseGradientType(static_cast<GradientType>(gradientType) );
	filter = actualfilter;

	// set up smoothing kernel for deformation field
	//float sigmaDef = 1.5;
	if (sigmaDef > 0.1)
	{
		filter->SmoothDeformationFieldOn();
		filter->SetStandardDeviations( sigmaDef );
	}
	else
	{
		filter->SmoothDeformationFieldOff();
	}

	// set up smoothing kernel for update field
	float sigmaUp = 0.0;
	if ( sigmaUp > 0.1 )
	{
		filter->SmoothUpdateFieldOn();
		filter->SetUpdateFieldStandardDeviations( sigmaUp );
	}
	else
	{
		filter->SmoothUpdateFieldOff();
	}

	//filter->SetIntensityDifferenceThreshold( 0.001 );

	// Set up the multi-resolution filter
	typedef itk::MultiResolutionPDEDeformableRegistration<
		InternalImageType, InternalImageType, DeformationFieldType, InternalPixelType >   MultiResRegistrationFilterType;
	MultiResRegistrationFilterType::Pointer multires = MultiResRegistrationFilterType::New();

	typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<DeformationFieldType,double> FieldInterpolatorType;

	FieldInterpolatorType::Pointer VectorInterpolator = FieldInterpolatorType::New();

#if ( ITK_VERSION_MAJOR > 3 ) || ( ITK_VERSION_MAJOR == 3 && ITK_VERSION_MINOR > 8 )
	multires->GetFieldExpander()->SetInterpolator(VectorInterpolator);
#endif
	std::vector<unsigned int> curNumIterations;
	// unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < res; i++)		curNumIterations.push_back(iterInResolutions[i]);

	multires->SetRegistrationFilter( filter );
	multires->SetNumberOfLevels( curNumIterations.size() );
	multires->SetNumberOfIterations( &curNumIterations[0] );
	multires->SetFixedImage( fixedImageReader->GetOutput() );
	if (doHistMatch)
	{
		multires->SetMovingImage( matcher->GetOutput() );
	}
	else
	{
		multires->SetMovingImage( movingImageReader->GetOutput() );
	}

	// set initial
	multires->SetArbitraryInitialDeformationField( initDeformationField );

	// Compute the deformation field
	try
	{
		multires->UpdateLargestPossibleRegion();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at demons registration between input images!" << std::endl; 
		std::cerr << excep << std::endl;
	}

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	ApplyDeformationFieldAndWriteWithTypeWithFileNames(movingImageFileName, 
		deformationFieldFileName, deformedImageFileName, true);

	//// write deformed image
	//InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	//InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	//warper->SetInput( movingImageReader->GetOutput() );
	//warper->SetInterpolator( interpolator );
	//warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	//warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	//warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	//warper->SetDeformationField( multires->GetOutput() );
	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	//if (isCompressed)
	//	CompressDeformationField2Short(deformationFieldFileName);
	return;

}

void generateMSTFromMatrix( double** curDistance, int nnode, int* treeMST )
{

	// implement Prim algorithm
	int* tree_v1 = new int[nnode-1];
	int* tree_v2 = new int[nnode-1];
	double* tree_el = new double[nnode-1];

	double* close_lowcost = new double[nnode];
	int* close_vec = new int[nnode];

	for (int i = 1; i < nnode; i++)
	{
		close_vec[i] = 0;
		close_lowcost[i] = curDistance[i][0];
	}
	close_lowcost[0] = -1;

	for (int i =0; i < nnode-1; i++)
	{
		// find node in graph
		int k;
		for (k =1; k < nnode; k++)
			if (close_lowcost[k] != -1)
				break;

		for (int j = k+1; j < nnode; j++)
			if (close_lowcost[j] != -1)
				if (close_lowcost[j] < close_lowcost[k])
					k = j;

		// add to tree
		tree_el[i] = close_lowcost[k];
		tree_v1[i] = k;
		tree_v2[i] = close_vec[k];
		close_lowcost[k] = -1;

		// adjust the candidate node set
		for (int j = 1; j < nnode; j++)
			if (close_lowcost[j] > curDistance[j][k])
			{
				close_lowcost[j] = curDistance[j][k];
				close_vec[j] = k;
			}
	}

	// implement Prim algorithm
	/////////////////////////////////////
	// output tree with node and edge

	//for (int i = 0; i < nnode-1; i++)
	//	std::cout << "(" << tree_v1[i] << ", "	<< tree_v2[i] << ") = " << tree_el[i] << std::endl;

	// output tree with node and edge
	/////////////////////////////////////////////
	// rebuild the distance matrix based on connectivity on MST
	double** curDistanceTemp = new double*[nnode];
	for (int i = 0; i < nnode; i++)
	{
		curDistanceTemp[i] = new double[nnode];
		for (int j = 0; j < nnode; j++)
			if (i != j)
				curDistanceTemp[i][j] = DBL_MAX;
			else
				curDistanceTemp[i][j] = 0;
	}
	for (int i = 0; i < nnode-1; i++)
	{
		curDistanceTemp[tree_v1[i]][tree_v2[i]] = tree_el[i];
		curDistanceTemp[tree_v2[i]][tree_v1[i]] = tree_el[i];
	}

	// rebuild the distance matrix based on connectivity on MST
	/////////////////////////////////////////////
	// set root to the default on in this particular scenario
	// otherwise, find the root on MST based on the shortest total distance to all other nodes

	int root = nnode-1;

	// find root based on the tree path
	FindCenterBasedOnShortestTotalDistance(curDistanceTemp, nnode, root);
	// find root based on the original graph
	//FindCenterBasedOnShortestTotalDistance(curDistance, nnode, root);

	// set root to the default on in this particular scenario
	// otherwise, find the root on MST based on the shortest total distance to all other nodes
	/////////////////////////////////////
	// generate tree based on MST and root
	//int* treeMST = new int[nnode];
	bool* treeEdgeUsed = new bool[nnode-1];
	for (int i = 0; i < nnode-1; i++) treeEdgeUsed[i] = false;
	treeMST[root] = -1;

	FromEdgeSegment2Tree(treeMST, treeEdgeUsed, tree_v1, tree_v2, root, nnode);

	// check if the generated tree is valid
	bool isAllEdgeUsed = true;
	for (int i = 0; i < nnode-1; i++)
		if (!treeEdgeUsed[i])
			isAllEdgeUsed = false;
	if (isAllEdgeUsed)
	{
		// //SaveTreeWithInfo(treeMST, nnode, currentTreeFileName);
		//std::cout << "Tree is generated!@" << std::endl;
	}
	else
		std::cout << "Tree cannot be generated!@" << std::endl;
	// generate tree based on MST and root
	/////////////////////////////////////
	// clean up
	//delete[] treeMST;
	delete[] close_vec;
	delete[] close_lowcost;
	delete[] tree_el;
	delete[] tree_v1;
	delete[] tree_v2;
	for (int i = 0; i < nnode; i++) delete[] curDistanceTemp[i];
	delete[] curDistanceTemp;
	return;
}
void generateMSTFromMatrixWithFixedRoot( double** curDistance, int nnode, int* treeMST, int root )
{

	// implement Prim algorithm
	int* tree_v1 = new int[nnode-1];
	int* tree_v2 = new int[nnode-1];
	double* tree_el = new double[nnode-1];

	double* close_lowcost = new double[nnode];
	int* close_vec = new int[nnode];

	for (int i = 1; i < nnode; i++)
	{
		close_vec[i] = 0;
		close_lowcost[i] = curDistance[i][0];
	}
	close_lowcost[0] = -1;

	for (int i =0; i < nnode-1; i++)
	{
		// find node in graph
		int k;
		for (k =1; k < nnode; k++)
			if (close_lowcost[k] != -1)
				break;

		for (int j = k+1; j < nnode; j++)
			if (close_lowcost[j] != -1)
				if (close_lowcost[j] < close_lowcost[k])
					k = j;

		// add to tree
		tree_el[i] = close_lowcost[k];
		tree_v1[i] = k;
		tree_v2[i] = close_vec[k];
		close_lowcost[k] = -1;

		// adjust the candidate node set
		for (int j = 1; j < nnode; j++)
			if (close_lowcost[j] > curDistance[j][k])
			{
				close_lowcost[j] = curDistance[j][k];
				close_vec[j] = k;
			}
	}

	// implement Prim algorithm
	/////////////////////////////////////
	// output tree with node and edge

	//for (int i = 0; i < nnode-1; i++)
	//	std::cout << "(" << tree_v1[i] << ", "	<< tree_v2[i] << ") = " << tree_el[i] << std::endl;

	// output tree with node and edge
	/////////////////////////////////////////////
	// rebuild the distance matrix based on connectivity on MST

	double** curDistanceTemp = new double*[nnode];
	for (int i = 0; i < nnode; i++)
	{
		curDistanceTemp[i] = new double[nnode];
		for (int j = 0; j < nnode; j++)
			if (i != j)
				curDistanceTemp[i][j] = DBL_MAX;
			else
				curDistanceTemp[i][j] = 0;
	}
	for (int i = 0; i < nnode-1; i++)
	{
		curDistanceTemp[tree_v1[i]][tree_v2[i]] = tree_el[i];
		curDistanceTemp[tree_v2[i]][tree_v1[i]] = tree_el[i];
	}

	// rebuild the distance matrix based on connectivity on MST
	/////////////////////////////////////////////
	// set root to the default on in this particular scenario
	// otherwise, find the root on MST based on the shortest total distance to all other nodes

	//int root = nnode-1;

	//FindCenterBasedOnShortestTotalDistance(curDistance, nnode, root);

	// set root to the default on in this particular scenario
	// otherwise, find the root on MST based on the shortest total distance to all other nodes
	/////////////////////////////////////
	// generate tree based on MST and root
	//int* treeMST = new int[nnode];
	bool* treeEdgeUsed = new bool[nnode-1];
	for (int i = 0; i < nnode-1; i++) treeEdgeUsed[i] = false;
	treeMST[root] = -1;

	FromEdgeSegment2Tree(treeMST, treeEdgeUsed, tree_v1, tree_v2, root, nnode);

	// check if the generated tree is valid
	bool isAllEdgeUsed = true;
	for (int i = 0; i < nnode-1; i++)
		if (!treeEdgeUsed[i])
			isAllEdgeUsed = false;
	if (isAllEdgeUsed)
	{
		// //SaveTreeWithInfo(treeMST, nnode, currentTreeFileName);
		//std::cout << "Tree is generated!@" << std::endl;
	}
	else
		std::cout << "Tree cannot be generated!@" << std::endl;
	// generate tree based on MST and root
	/////////////////////////////////////
	// clean up
	//delete[] treeMST;
	delete[] close_vec;
	delete[] close_lowcost;
	delete[] tree_el;
	delete[] tree_v1;
	delete[] tree_v2;
	for (int i = 0; i < nnode; i++) delete[] curDistanceTemp[i];
	delete[] curDistanceTemp;
	return;
}

void FromEdgeSegment2Tree(int* tree, bool* treeEdgeUsed, int* v1, int* v2, int curRoot, int matrixSize)
{
	int child = -1;
	for (int i = 0; i < matrixSize-1; i++)
	{
		if ( (treeEdgeUsed[i] == false) && ((v1[i] == curRoot) || (v2[i] == curRoot)) )
		{
			child = v1[i]+v2[i]-curRoot;
			tree[child] = curRoot;
			treeEdgeUsed[i] = true;
			FromEdgeSegment2Tree(tree, treeEdgeUsed, v1, v2, child, matrixSize);
		}
	}
	return;
}
void FindCenterBasedOnShortestTotalDistance(double** distanceMatrix, int matrixSize, int &center)
{

	// input distance matrix should have very large value on the position defined by two un-connecting nodes

	////////////////////////
	double** distance = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)		distance[i] = new double[matrixSize];

	double** distance_temp = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)		distance_temp[i] = new double[matrixSize];

	// make a copy of distanceMatrix into distance and find maximum
	for (int i = 0; i < matrixSize; i++)
		for (int j = 0; j < matrixSize; j++)
			distance[i][j] = distanceMatrix[i][j];

	// make distance matrix be symmetric
	for (int i = 0; i < matrixSize; i++)
		for (int j = i+1; j < matrixSize; j++)
			if (distance[i][j] > distance[j][i])
				distance[i][j] = distance[j][i];
			else
				distance[j][i] = distance[i][j];
	// make a copy of the symmetric distance matrix
	for(int i = 0; i < matrixSize; i++)
		for(int j = 0; j < matrixSize; j++)
			distance_temp[i][j] = distance[i][j];

	// find the minimum distance between two points on kNN map
	int ind = 0;
	while(ind < matrixSize)
	{
		for (int i = 0; i < matrixSize; i++)
			for (int j = 0; j < matrixSize; j++)
				if (distance[i][ind]+distance[ind][j]<distance[i][j])
					distance_temp[i][j] = distance[i][ind]+distance[ind][j];
		for (int i = 0; i < matrixSize; i++)
			for (int j = 0; j < matrixSize; j++)
				distance[i][j] = distance_temp[i][j];
		ind = ind+1;
	}

	// calculate the sum of each row

	int* index = new int[matrixSize];

	double* distanceSum = new double[matrixSize];
	for (int i = 0; i < matrixSize; i++) distanceSum[i] = 0.0;
	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j< matrixSize; j++)
		{
			distanceSum[i] += distance[i][j];
		}
		index[i] = i;
	}
	bubbleSort(distanceSum, index, matrixSize);

	// center = index[matrixSize-1];
	center = index[0];
	// clean up
	delete[] index;
	delete[] distanceSum;

	for (int i = 0; i < matrixSize; i++)		delete[] distance[i];
	delete[] distance;
	for (int i = 0; i < matrixSize; i++)		delete[] distance_temp[i];
	delete[] distance_temp;

	return;
}

bool ValidateTree(int* tree, int tree_size)
{
	// to check if there is one and only one root
	int root;
	int root_num = 0;
	for (int i = 0; i < tree_size; i++)
	{
		if (tree[i] == -1)	
		{
			root = i;
			root_num++;
		}
		if ((tree[i] < -1) || (tree[i] >= tree_size))
			return false;
	}
	if (root_num != 1)
		return false;

	// to check if all node can trace back to the root
	for (int i = 0; i < tree_size; i++)
	{
		//
		bool isRooted = false;
		int max_step = tree_size;
		int j = 0;
		int cur_node = i;
		while (j < max_step)
		{
			if (cur_node == root) 
			{
				isRooted = true;
				break;
			}
			else
			{
				cur_node = tree[cur_node];
				j++;
			}
		}
		if (!isRooted)
			return false;
	}

	return true;
}
void TreeBasedRegistrationFast(int* tree, int tree_size, char** originalIntensityImageFileNames,
							   char** deformedImgFileNames, char** deformationFieldFileNames)
{
	int root = -1; // 
	FindRoot(tree, tree_size, root);

	// calculate the depth of each node
	int* node_depth = new int[tree_size];
	for (int i = 0; i < tree_size; i++)
	{
		int curnode = i;
		node_depth[i] = 0;
		bool isRoot = false;

		while (!isRoot)
		{
			if (tree[curnode]>=0)
			{
				node_depth[i]++;
				curnode = tree[curnode];
			}
			else
				isRoot = true;
		}
	}

	// sort all nodes by its depth
	int* index = new int[tree_size];
	{
		//
		double* arr = new double[tree_size];
		for (int i = 0; i < tree_size; i++)
		{
			arr[i] = node_depth[i];
			index[i] = i;
		}
		bubbleSort(arr, index, tree_size); // index[0] = root
		delete[] arr;
	}

	// start to register each image to the root node step by step
	for (int ii = 1; ii < tree_size; ii++) // starting from 1, since index[0] showing the root
	{
		if (isDebug)
			std::cerr << ii << ", ";
		int i = index[ii];
		int curnode = index[ii];
		int parentnode = tree[curnode];

		bool isRegistered = true;
		if (!itksys::SystemTools::FileExists(deformationFieldFileNames[curnode], true)) isRegistered = false;
		if (!itksys::SystemTools::FileExists(deformedImgFileNames[curnode], true)) isRegistered = false;
		char* deformedImgFileNameImg = new char[MAX_FILE_NAME_LENGTH];
		strcpy(deformedImgFileNameImg, deformedImgFileNames[curnode]);
		deformedImgFileNameImg[strlen(deformedImgFileNameImg)-3] = 'i';
		deformedImgFileNameImg[strlen(deformedImgFileNameImg)-2] = 'm';
		deformedImgFileNameImg[strlen(deformedImgFileNameImg)-1] = 'g';
		if (!itksys::SystemTools::FileExists(deformedImgFileNameImg, true)) isRegistered = false;

		// do registration
		if (isRegistered)
		{
			continue; // goto next one
		}

		if (parentnode == root) // direct registration without initial deformation
		{

			//DeformationFieldType::Pointer deformationField = 0;
			//DiffeoDemonsRegistration(originalIntensityImageFileNames[root], originalIntensityImageFileNames[i],
			//	deformationField, sigmaDef, doHistMatch, iterInResolutions);

			//WriteDeformationField(outputDeformationFieldFileName, deformationField);

			//ApplyDeformationFieldAndWriteWithFileNames(originalIntensityImageFileNames[i], 
			//	outputDeformationFieldFileName, outputImageName, true);
			DiffeoDemonsRegistrationWithParameters(originalIntensityImageFileNames[root], 
				originalIntensityImageFileNames[curnode], 
				deformedImgFileNames[curnode], deformationFieldFileNames[curnode], 
				sigmaDef, doHistMatch, iterInResolutions[itereach]);

		}
		else // registration with initial deformation
		{
			//
			//char initDeformationFieldFileName[MAX_FILE_NAME_LENGTH];
			//strcpy(initDeformationFieldFileName, originalIntensityImageFileNames[parentnode]);
			//initDeformationFieldFileName[strlen(initDeformationFieldFileName)-4] = 0;
			//strcat(initDeformationFieldFileName, "_deform.mha");

			//DiffeoDemonsRegistrationWithInitialDeformationField(originalIntensityImageFileNames[root], originalIntensityImageFileNames[i],
			//	initDeformationFieldFileName, outputImageName, outputDeformationFieldFileName, sigmaDef, doHistMatch, iterInResolutions);
			DiffeoDemonsRegistrationWithInitialWithParameters(originalIntensityImageFileNames[root], 
				originalIntensityImageFileNames[curnode], deformationFieldFileNames[parentnode],
				deformedImgFileNames[curnode], deformationFieldFileNames[curnode], 
				sigmaDef, doHistMatch, iterInResolutions[itereach]);

		}
	}
	std::cout << "Done!" << std::endl;

	delete[] node_depth;
	delete[] index;
}
void TreeBasedRegistrationFastOniTree(int* itree_t, int itree_size_t, int atlas_size_t, int simulate_size_t, int sample_size_t,
									  char** originalIntensityImageFileNames,
									  char** deformedImgFileNames,
									  char** deformationFieldFileNames)
{
	int root = -1; // 
	FindRoot(itree_t, itree_size_t, root);

	// calculate the depth of each node
	int* node_depth = new int[itree_size_t];
	for (int i = 0; i < itree_size_t; i++)
	{
		int curnode = i;
		node_depth[i] = 0;
		bool isRoot = false;

		while (!isRoot)
		{
			if (itree_t[curnode]>=0)
			{
				node_depth[i]++;
				curnode = itree_t[curnode];
			}
			else
				isRoot = true;
		}
	}

	// sort all nodes by its depth
	int* index = new int[itree_size_t];
	{
		//
		double* arr = new double[itree_size_t];
		for (int i = 0; i < itree_size_t; i++)
		{
			arr[i] = node_depth[i];
			index[i] = i;
		}
		bubbleSort(arr, index, itree_size_t); // index[0] = root
		delete[] arr;
	}

	// start to register each image to the root node step by step
	for (int ii = 1; ii < itree_size_t; ii++) // starting from 1, since index[0] showing the root
	{
		//if (isDebug)
			std::cerr << ii << ", ";
		int i = index[ii];
		int curnode = index[ii];
		int parentnode = itree_t[curnode];

		if ((curnode>=atlas_size_t) && (curnode < atlas_size_t+simulate_size_t))
		{
			// this is a simulated image, then pass
			continue;
		}

		int curnode_t = 0;

		if (curnode >= atlas_size_t+simulate_size_t)
			curnode_t = curnode - simulate_size_t;
		else
			curnode_t = curnode;
		int parentnode_t = 0;
		if (parentnode >= atlas_size_t+simulate_size_t)
			parentnode_t = parentnode - simulate_size_t;
		else
			parentnode_t = parentnode;

		bool isRegistered = true;
		{

			if (!itksys::SystemTools::FileExists(deformationFieldFileNames[curnode_t], true)) isRegistered = false;
			if (!itksys::SystemTools::FileExists(deformedImgFileNames[curnode_t], true)) isRegistered = false;
			
			char* deformedImgFileNameImg = new char[MAX_FILE_NAME_LENGTH];
			strcpy(deformedImgFileNameImg, deformedImgFileNames[curnode_t]);
			deformedImgFileNameImg[strlen(deformedImgFileNameImg)-3] = 'i';
			deformedImgFileNameImg[strlen(deformedImgFileNameImg)-2] = 'm';
			deformedImgFileNameImg[strlen(deformedImgFileNameImg)-1] = 'g';
			
			if (!itksys::SystemTools::FileExists(deformedImgFileNameImg, true)) isRegistered = false;

		}

		// do registration
		if (isRegistered)
		{
			//CompressDeformationField2Short(deformationFieldFileNames[curnode_t]);
			continue; // goto next one
		}

		if (parentnode == root) // direct registration without initial deformation
		{

			//DeformationFieldType::Pointer deformationField = 0;
			//DiffeoDemonsRegistration(originalIntensityImageFileNames[root], originalIntensityImageFileNames[i],
			//	deformationField, sigmaDef, doHistMatch, iterInResolutions);

			//WriteDeformationField(outputDeformationFieldFileName, deformationField);

			//ApplyDeformationFieldAndWriteWithFileNames(originalIntensityImageFileNames[i], 
			//	outputDeformationFieldFileName, outputImageName, true);
			DiffeoDemonsRegistrationWithParameters(originalIntensityImageFileNames[root], 
				originalIntensityImageFileNames[curnode_t], 
				deformedImgFileNames[curnode_t], deformationFieldFileNames[curnode_t], 
				sigmaDef, doHistMatch, iterInResolutions[itereach]);
			//if (isCompressed)
			//	CompressDeformationField2Short(deformationFieldFileNames[curnode_t]);

		}
		else // registration with initial deformation
		{
			//
			//char initDeformationFieldFileName[MAX_FILE_NAME_LENGTH];
			//strcpy(initDeformationFieldFileName, originalIntensityImageFileNames[parentnode]);
			//initDeformationFieldFileName[strlen(initDeformationFieldFileName)-4] = 0;
			//strcat(initDeformationFieldFileName, "_deform.mha");

			//DiffeoDemonsRegistrationWithInitialDeformationField(originalIntensityImageFileNames[root], originalIntensityImageFileNames[i],
			//	initDeformationFieldFileName, outputImageName, outputDeformationFieldFileName, sigmaDef, doHistMatch, iterInResolutions);
			
			if ((parentnode < atlas_size_t) || (parentnode >= atlas_size_t+simulate_size_t))
			{
				// parent is a real image
				DiffeoDemonsRegistrationWithInitialWithParameters(originalIntensityImageFileNames[root], 
					originalIntensityImageFileNames[curnode_t], deformationFieldFileNames[parentnode_t],
					deformedImgFileNames[curnode_t], deformationFieldFileNames[curnode_t], 
					sigmaDef, doHistMatch, iterInResolutions[itereach]);
				//if (isCompressed)
				//	CompressDeformationField2Short(deformationFieldFileNames[curnode_t]);
			}
			else
			{
				// parent is a simulated image
				char initDeformationFieldFileName[MAX_FILE_NAME_LENGTH];
				MakeFileName(initDeformationFieldFileName, "simulated", "_deform_", 
					parentnode-atlas_size_t, "mha");

				DiffeoDemonsRegistrationWithInitialWithParameters(originalIntensityImageFileNames[root], 
					originalIntensityImageFileNames[curnode_t], initDeformationFieldFileName,
					deformedImgFileNames[curnode_t], deformationFieldFileNames[curnode_t], 
					sigmaDef, doHistMatch, iterInResolutions[itereach]);
				//if (isCompressed)
				//	CompressDeformationField2Short(deformationFieldFileNames[curnode_t]);

			}

		}
	}
	std::cout << "Done!" << std::endl;

	delete[] node_depth;
	delete[] index;
}
void RemoveFile(char* filename)
{
	if (itksys::SystemTools::FileExists(filename, true))
	{
#if defined(_WIN32) && !defined(__CYGWIN__)
		std::ostringstream oss (std::ostringstream::out);
		oss << "del /F " << filename << std::endl;
		//std::cerr << oss.str().c_str();
		system(oss.str().c_str());

#else
		std::ostringstream oss (std::ostringstream::out);
		oss << "rm -f " << filename << std::endl;
		//std::cerr << oss.str().c_str();
		system(oss.str().c_str());
#endif
	}
	return;
}
void RemoveFile(char* filename, bool verbose)
{
	if (itksys::SystemTools::FileExists(filename, true))
	{
#if defined(_WIN32) && !defined(__CYGWIN__)
		std::ostringstream oss (std::ostringstream::out);
		oss << "del /F " << filename << std::endl;
		if (verbose)
			std::cerr << oss.str().c_str();
		system(oss.str().c_str());

#else
		std::ostringstream oss (std::ostringstream::out);
		oss << "rm -f " << filename << std::endl;
		if (verbose)
			std::cerr << oss.str().c_str();
		system(oss.str().c_str());
#endif
	}
	else
	{
		if (verbose)
			std::cerr << filename << " does not exist!" << std::endl;
	}
	return;
}
void RemoveMoreFiles(char* filenames, bool verbose)
{
#if defined(_WIN32) && !defined(__CYGWIN__)
		std::ostringstream oss (std::ostringstream::out);
		oss << "del /F " << filenames << std::endl;
		if (verbose)
			std::cerr << oss.str().c_str();
		system(oss.str().c_str());

#else
		std::ostringstream oss (std::ostringstream::out);
		oss << "rm -f " << filenames << std::endl;
		if (verbose)
			std::cerr << oss.str().c_str();
		system(oss.str().c_str());
#endif
	return;
}
void DoPCATraining(char** deformationFieldFileNames, int numFiles, char** allImgFileName, int root)
{
	bool doPCATraining = true;
	int sampleRate = 2;
	bool doResampleDeformationField = true;
	int size_x = 0;
	int size_y = 0;
	int size_z = 0;
	int size_xn = 0;
	int size_yn = 0;
	int size_zn = 0;
	int size_df = 0;  // size of deformation field
	int size_dfn = 0;  // size of sub-sampled deformation field
	//int numEigenVector = 4; // t
	//int numSampleEachDirection = 4; // n  n^t total number of intermediate templates
	
	//c = sqrt(2)*erfinv(h);
	//float h[] = {-0.6, -0.2, 0.2, 0.6};
	//float c[] = {-0.8416, -0.2533, 0.2533, 0.8416};

	float h3[] = {-0.5, 0.0, 0.5};
	float c3[] = {-0.6745, 0.0, 0.6745};

	float h4[] = {-0.6, -0.2, 0.2, 0.6};
	float c4[] = {-0.8416, -0.2533, 0.2533, 0.8416};

	float h5[] = {-0.667, -0.333, 0.0, 0.333, 0.667};
	float c5[] = {-0.9681, -0.4303, 0.0, 0.4303, 0.9681};

	float* c = NULL;
	float* h = NULL;
	if (numSampleEachDirection == 3)
	{
		c = new float[3];
		h = new float[3];
		for (int i = 0; i < numSampleEachDirection; i++)
		{
			c[i] = c3[i];
			h[i] = h3[i];
		}
	}
	else if (numSampleEachDirection == 4)
	{
		c = new float[4];
		h = new float[4];
		for (int i = 0; i < numSampleEachDirection; i++)
		{
			c[i] = c4[i];
			h[i] = h4[i];
		}
	}
	else if (numSampleEachDirection == 5)
	{
		c = new float[5];
		h = new float[5];
		for (int i = 0; i < numSampleEachDirection; i++)
		{
			c[i] = c5[i];
			h[i] = h5[i];
		}
	}
	else
	{
		std::cerr << "The number of samples on each eigen vector direction should be between 3 and 5!" << std::endl;
		exit(1);
	}
	// some global parameters
	// SVD   (Vt == V^T,transpose; Si = inv(S), inverse)
	// D = U S Vt   =>  DV = US => DVSi = U
	
	///////////////////////////////////////
	// initialize
	
	DeformationFieldType::Pointer curDeformationField = 0;
	if (!isCompressed)
		ReadDeformationField(deformationFieldFileNames[0], curDeformationField);
	else
		DeCompressDeformationFieldFromShort(deformationFieldFileNames[0], curDeformationField);

	DeformationFieldType::SizeType im_size = curDeformationField->GetLargestPossibleRegion().GetSize();
	size_x = im_size[0]; size_y = im_size[1]; size_z = im_size[2];
	size_xn = (size_x-1)/sampleRate+1;size_yn = (size_y-1)/sampleRate+1;size_zn = (size_z-1)/sampleRate+1;
	size_dfn = (size_xn)*(size_yn)*(size_zn)*3;
	size_df = size_x*size_y*size_z*3;

	df_spacing = curDeformationField->GetSpacing();
	df_direction = curDeformationField->GetDirection();
	df_origin = curDeformationField->GetOrigin();
	
	// initialize
	////////////////////////////////////////
	// do PCA training	
	////////////////////////////////////////
	// std::cerr << "Resample deformation fields ..." << std::endl;
	vnl_matrix<float> df_eigenvector(size_dfn, numEigenVector);
	vnl_vector<float> df_eigenvalues(numEigenVector);

	////////////////////////////////////////
	// down sample
	for (int i = 0; i < numFiles; i++)
		//for (int i = 1; i < 2; i++)
	{
		// std::cerr << i << ", ";
		// make file names
		char deformationFieldFileName[MAX_FILE_NAME_LENGTH];
		strcpy(deformationFieldFileName, deformationFieldFileNames[i]);

		char resampledDeformationFieldFileName[MAX_FILE_NAME_LENGTH];
		strcpy(resampledDeformationFieldFileName, deformationFieldFileNames[i]);
		resampledDeformationFieldFileName[strlen(resampledDeformationFieldFileName)-4] = 0;
		strcat(resampledDeformationFieldFileName,"_sub.mha");
		if (doResampleDeformationField)
		{
			DownResampleDeformationField(deformationFieldFileName, resampledDeformationFieldFileName, sampleRate);
		}
	}
	// std::cerr << "Done!" << std::endl;

	// down sample
	/////////////////////////////////
	
	vnl_matrix<float> df_mat(size_dfn, numFiles, 0.0);
	vnl_vector<float> df_mean(size_dfn, 0.0);

	if (doPCATraining)
	{
		// build all vectors from deformation fields

		for (int i = 0; i < numFiles; i++)
		{
			char resampledDeformationFieldFileName[MAX_FILE_NAME_LENGTH];
			strcpy(resampledDeformationFieldFileName, deformationFieldFileNames[i]);
			resampledDeformationFieldFileName[strlen(resampledDeformationFieldFileName)-4] = 0;
			strcat(resampledDeformationFieldFileName,"_sub.mha");
			float* df_arr = new float[size_dfn];
			LoadIntoArray(resampledDeformationFieldFileName, df_arr);
			vnl_vector<float> df_cur(size_dfn);
			df_cur.copy_in(df_arr);
			//for (int ii = 0; ii < size_dfn; ii++)
			//{
			//	df_cur[ii] = df_arr[ii];
			//	df_cur.data_block()[ii] = df_arr[ii];
			//}
			df_mean += df_cur;
			df_mat.set_column(i, df_cur);
			delete[] df_arr;
		}
		// calculate mean vector
		df_mean /= (float)(numFiles);

		// subtract mean vector from each vector
		for (int i = 0; i < numFiles; i++)
		{
			df_mat.set_column(i, df_mat.get_column(i)-df_mean);
		}

		// do SVD
		// df_mat.inplace_transpose();
		vnl_svd_economy<float> svd_e(df_mat);
		//vnl_matlab_print(vcl_cerr, svd_e.lambdas());
		//vcl_cerr << vcl_endl;
		//vnl_matlab_print(vcl_cerr, svd_e.V());
		//vcl_cerr << vcl_endl;
		//vcl_cerr << vcl_endl;
		int rows = svd_e.V().rows();
		int cols = svd_e.V().cols();

		// build eigen-vector matrix of original deformation field matrix

		for (int i = 0; i < numEigenVector; i++)
		{
			// U = D*V*Sinv
			df_eigenvalues.data_block()[i] = svd_e.lambdas().data_block()[i];
			df_eigenvector.set_column(i,df_mat*svd_e.V().get_column(i));
			//df_eigenvector.scale_column(i, 1/df_eigenvalues.data_block()[i]);
		}
		// calculate coefficients and generate all intermediate templates
		//char intermediateTemplateListName[MAX_FILE_NAME_LENGTH];
		//strcpy(intermediateTemplateListName, "simulated_image_list.txt");
		//std::ofstream intermediateFileNamesFile;
		//intermediateFileNamesFile.open (intermediateTemplateListName);
	}// end if doPCA
	if (1) // if (0) // if (1)
	{
		int* coeff = new int[numEigenVector];
		int numAllCombinations = (int)pow((float)numSampleEachDirection, numEigenVector);

		//for (int i = 0; i < 0; i++)
		for (int i = 0; i < numAllCombinations; i++)
		{
			std::cerr << i << ", ";
			//
			for (int j = 0; j < numEigenVector; j++) coeff[j] = 0;
			int num = i;
			int index = 0;
			while (num > 0)
			{
				coeff[index] = (num % numSampleEachDirection);
				num = (num - coeff[index])/numSampleEachDirection;
				index++;
			}
			// generate subsampled intermediate deformation field
			// create a new file list
			vnl_vector<float> df_intermediate_sub(size_dfn);
			df_intermediate_sub = df_mean;
			for (int j = 0; j < numEigenVector; j++)
			{
				//df_intermediate_sub += c[coeff[j]]*df_eigenvector.get_column(j)*df_eigenvalues.data_block()[j];
				df_intermediate_sub += c[coeff[j]]*df_eigenvector.get_column(j);
			}

			char index_string[5]; 		myitoa( i, index_string, 3 );
			char intermediateSubDeformationFieldFileName[] = "inter_deform_sub_000.mha";
			intermediateSubDeformationFieldFileName[strlen(intermediateSubDeformationFieldFileName)-7]=0;
			strcat(intermediateSubDeformationFieldFileName, index_string);
			strcat(intermediateSubDeformationFieldFileName, ".mha");
			char intermediateDeformationFieldFileName[] = "inter_deform_000.mha";
			intermediateDeformationFieldFileName[strlen(intermediateDeformationFieldFileName)-7]=0;
			strcat(intermediateDeformationFieldFileName, index_string);
			strcat(intermediateDeformationFieldFileName, ".mha");

			SaveFromArray(intermediateSubDeformationFieldFileName, df_intermediate_sub.data_block(), size_xn, size_yn,size_zn);
			// upSample the above deformation field and save
			if (doResampleDeformationField)
			{
				UpResampleDeformationField(intermediateSubDeformationFieldFileName,
					intermediateDeformationFieldFileName, sampleRate);
			}

			char intermediateTemplateFileName[] = "inter_template_000.hdr";
			intermediateTemplateFileName[strlen(intermediateTemplateFileName)-7]=0;
			strcat(intermediateTemplateFileName, index_string);
			strcat(intermediateTemplateFileName, ".hdr");

			// reverse deformation field
			char intermediateReversedDeformationFieldFileName[] = "inter_deform_reverse_000.mha";
			intermediateReversedDeformationFieldFileName[strlen(intermediateReversedDeformationFieldFileName)-7]=0;
			strcat(intermediateReversedDeformationFieldFileName, index_string);
			strcat(intermediateReversedDeformationFieldFileName, ".mha");
			DeformationFieldType::Pointer deformationField = 0;
			if (!isCompressed)
				ReadDeformationField(intermediateDeformationFieldFileName, deformationField);
			else
				DeCompressDeformationFieldFromShort(intermediateDeformationFieldFileName, deformationField);

			DeformationFieldType::Pointer reversedDeformationField = DeformationFieldType::New();
			reversedDeformationField -> SetRegions (deformationField->GetLargestPossibleRegion());
			reversedDeformationField -> SetSpacing (deformationField->GetSpacing());
			reversedDeformationField -> SetDirection (deformationField->GetDirection());
			reversedDeformationField -> SetOrigin (deformationField->GetOrigin());
			reversedDeformationField -> Allocate();
			InverseDeformationFieldDG3D(deformationField, reversedDeformationField);
			WriteDeformationField(intermediateReversedDeformationFieldFileName, reversedDeformationField);
			if (isCompressed)
				CompressDeformationField2Short(intermediateReversedDeformationFieldFileName);
			// apply intermediate deformation field to get intermediate template
			InternalImageType::Pointer templateImage = 0;
			ReadImage(allImgFileName[root], templateImage);
			InternalImageType::Pointer intermediateTemplateImage = 0;
			ApplyDeformationField(templateImage, reversedDeformationField, intermediateTemplateImage,true);
			// write image
			WriteImage(intermediateTemplateFileName, intermediateTemplateImage);
			
			//intermediateFileNamesFile << intermediateTemplateFileName << std::endl;
		}// end for numOfAllCombinations
		//intermediateFileNamesFile.close();
		delete[] coeff;
		
		// do selection
		{
			std::cerr << "Select the best ... ";
			// calculate pairwise distance between each atlas image and each simulated image
			double** dist = new double*[atlas_size];
			for (int i = 0; i < atlas_size; i++)
			{
				dist[i] = new double[numAllCombinations];
				for (int j = 0; j < numAllCombinations; j++)
				{
					dist[i][j] = 0.0;
				}
			}
			double distmax = 0.0;
			for (int i = 0; i < atlas_size; i++)
			{
				for (int j = 0; j < numAllCombinations; j++)
				{
					char index_string[5]; 		myitoa( j, index_string, 3 );
					char curInterTempFileName[] = "inter_template_000.hdr";
					curInterTempFileName[strlen(curInterTempFileName)-7]=0;
					strcat(curInterTempFileName, index_string);
					strcat(curInterTempFileName, ".hdr");

					dist[i][j] = calculateDistanceMSD(allImgFileName[i], curInterTempFileName);
					if (dist[i][j] > distmax)
					{
						distmax = dist[i][j];
					}
				}
			}
			
			// calculate the minimum distance between each simulated image and all atlas images
			double* dist_min_each_inter = new double[numAllCombinations];
			int* index = new int[numAllCombinations];
			
			for (int j = 0; j < numAllCombinations; j++)
			{
				index[j] = j;
				dist_min_each_inter[j] = distmax;
				for (int i = 0; i < atlas_size; i++)
				{
					if (dist[i][j] < dist_min_each_inter[j])
					{
						dist_min_each_inter[j] = dist[i][j];
					}
				}
			}
			// sort
			bubbleSort(dist_min_each_inter, index, numAllCombinations);

			// copy files to new names after selection
			for (int i = 0; i < simulate_size; i++)
			{
				int index_inter = index[i + (numAllCombinations-simulate_size)];
				
				char index_string[5]; 		myitoa( index_inter, index_string, 3 );
				char curInterTempFileName[] = "inter_deform_000.mha";
				curInterTempFileName[strlen(curInterTempFileName)-7]=0;
				strcat(curInterTempFileName, index_string);
				strcat(curInterTempFileName, ".mha");

				char outDFName[MAX_FILE_NAME_LENGTH];
				MakeFileName(outDFName, "simulated", "_deform_", i, "mha");

				itksys::SystemTools::CopyFileAlways(curInterTempFileName, outDFName);
			}
			// remove irrelevant files
			for (int i = 0; i < numAllCombinations; i++)
			{
				//
				char curInterTempFileName[MAX_FILE_NAME_LENGTH];
				char curInterTempImgFileName[MAX_FILE_NAME_LENGTH];
				char curInterTempDeformFileName[MAX_FILE_NAME_LENGTH];
				char curInterTempRevFileName[MAX_FILE_NAME_LENGTH];
				char curInterTempSubFileName[MAX_FILE_NAME_LENGTH];
				MakeFileName(curInterTempFileName, "inter", "_template_", i, "hdr");
				MakeFileName(curInterTempImgFileName, "inter", "_template_", i, "img");
				MakeFileName(curInterTempDeformFileName, "inter", "_deform_", i, "mha");
				MakeFileName(curInterTempRevFileName, "inter", "_deform_reverse_", i, "mha");
				MakeFileName(curInterTempSubFileName, "inter", "_deform_sub_", i, "mha");
				RemoveFile(curInterTempFileName);
				RemoveFile(curInterTempImgFileName);
				RemoveFile(curInterTempDeformFileName);
				RemoveFile(curInterTempRevFileName);
				RemoveFile(curInterTempSubFileName);
			}
			std::cerr << "Done!" << std::endl;
			// clear newed variables
			delete[] index;
			for (int i = 0; i < atlas_size; i++) delete[] dist[i];
			delete[] dist;
			delete[] dist_min_each_inter;
		}

	} 
	else
	{
		// simplify the process
		
		// first round
		{
			//
			for (int index = 0; index < numFiles; index++)
			{
				int index1 = index;
				int index2 = index+1;
				if (index2 == numFiles) index2 = 0;
				DeformationFieldType::Pointer df1 = 0;
				DeformationFieldType::Pointer df2 = 0;
				if (!isCompressed)
				{
					ReadDeformationField(deformationFieldFileNames[index1], df1);
					ReadDeformationField(deformationFieldFileNames[index2], df2);
				}
				else
				{
					DeCompressDeformationFieldFromShort(deformationFieldFileNames[index1], df1);
					DeCompressDeformationFieldFromShort(deformationFieldFileNames[index2], df2);
				}
				
				MultiplyImageFilterType::Pointer multiplyImageFilter1 = MultiplyImageFilterType::New();
				multiplyImageFilter1->SetInput(df1);
				multiplyImageFilter1->SetConstant(0.5);
				MultiplyImageFilterType::Pointer multiplyImageFilter2 = MultiplyImageFilterType::New();
				multiplyImageFilter2->SetInput(df2);
				multiplyImageFilter2->SetConstant(0.5);

				AddImageFilterType::Pointer addImageSum = AddImageFilterType::New();
				addImageSum->SetInput1(multiplyImageFilter1->GetOutput());
				addImageSum->SetInput2(multiplyImageFilter2->GetOutput());
				try
				{
					addImageSum->Update();
				}
				catch( itk::ExceptionObject & excep )
				{
					std::cerr << "Exception caught!" << std::endl;
					std::cerr << excep << std::endl;
				}
				
				DeformationFieldType::Pointer avgDF = 0;
				avgDF = addImageSum->GetOutput();
				avgDF -> DisconnectPipeline();
				
				char outDFName[MAX_FILE_NAME_LENGTH];
				MakeFileName(outDFName, "simulated", "_deform_", index,"mha");
				WriteDeformationField(outDFName, avgDF);
				if (isCompressed)
					CompressDeformationField2Short(outDFName);
			}
		}
		// second round
		{
			//
			for (int index = 0; index < numFiles; index++)
			{
				DeformationFieldType::Pointer df1 = 0;
				if (!isCompressed)
					ReadDeformationField(deformationFieldFileNames[index], df1);
				else
					DeCompressDeformationFieldFromShort(deformationFieldFileNames[index], df1);

				MultiplyImageFilterType::Pointer multiplyImageFilter1 = MultiplyImageFilterType::New();
				multiplyImageFilter1->SetInput(df1);
				multiplyImageFilter1->SetConstant(-1.0);
				
				try
				{
					multiplyImageFilter1->Update();
				}
				catch( itk::ExceptionObject & excep )
				{
					std::cerr << "Exception caught!" << std::endl;
					std::cerr << excep << std::endl;
				}

				DeformationFieldType::Pointer avgDF = 0;
				avgDF = multiplyImageFilter1->GetOutput();
				avgDF -> DisconnectPipeline();

				char outDFName[MAX_FILE_NAME_LENGTH];
				MakeFileName(outDFName, "simulated", "_deform_", index+numFiles,"mha");
				WriteDeformationField(outDFName, avgDF);
				if (isCompressed)
					CompressDeformationField2Short(outDFName);

			}
		}
	}
	// do PCA training
	////////////////////////////////////////
	delete[] c;
	delete[] h;
	return;
}

void DownResampleDeformationField(char* deformationFieldFileName, char* resampledDeformationFieldFileName, int sampleRate)
{
	// 
	int rx, ry, rz;
	rx = sampleRate;
	ry = sampleRate;
	rz = sampleRate;

	DeformationFieldType::Pointer originImage = 0;
	if (!isCompressed)
		ReadDeformationField(deformationFieldFileName, originImage);
	else
		DeCompressDeformationFieldFromShort(deformationFieldFileName, originImage);

	int im_x, im_y, im_z;
	int im_xn, im_yn, im_zn;
	DeformationFieldType::SizeType im_size = originImage->GetLargestPossibleRegion().GetSize();
	im_x = im_size[0]; im_xn = (im_x-1)/rx+1;
	im_y = im_size[1]; im_yn = (im_y-1)/ry+1;
	im_z = im_size[2]; im_zn = (im_z-1)/rz+1;

	// create an array to store original image
	float*** originImageRawX = new float**[im_z];
	float*** originImageRawY = new float**[im_z];
	float*** originImageRawZ = new float**[im_z];
	for (int k = 0; k < im_z; k++) 
	{
		originImageRawX[k] = new float*[im_y];
		originImageRawY[k] = new float*[im_y];
		originImageRawZ[k] = new float*[im_y];
		for (int j = 0; j < im_y; j++)
		{
			originImageRawX[k][j] = new float[im_x];
			originImageRawY[k][j] = new float[im_x];
			originImageRawZ[k][j] = new float[im_x];
			for (int i = 0; i < im_x; i++)
			{
				originImageRawX[k][j][i] = 0.0;
				originImageRawY[k][j][i] = 0.0;
				originImageRawZ[k][j][i] = 0.0;
			}
		}
	}

	// create an array to store sampled image
	float*** sampledImageRawX = new float**[im_zn];
	float*** sampledImageRawY = new float**[im_zn];
	float*** sampledImageRawZ = new float**[im_zn];
	for (int k = 0; k < im_zn; k++) 
	{
		sampledImageRawX[k] = new float*[im_yn];
		sampledImageRawY[k] = new float*[im_yn];
		sampledImageRawZ[k] = new float*[im_yn];
		for (int j = 0; j < im_yn; j++)
		{
			sampledImageRawX[k][j] = new float[im_xn];
			sampledImageRawY[k][j] = new float[im_xn];
			sampledImageRawZ[k][j] = new float[im_xn];
			for (int i = 0; i < im_xn; i++)
			{
				sampledImageRawX[k][j][i] = 0.0;
				sampledImageRawY[k][j][i] = 0.0;
				sampledImageRawZ[k][j][i] = 0.0;
			}
		}
	}

	// load original image
	DeformationFieldIteratorType itOrigin(originImage, originImage->GetLargestPossibleRegion() );
	VectorPixelType vectorPixel;
	DeformationFieldType::IndexType idx;
	int idx_x, idx_y, idx_z;

	for (itOrigin.GoToBegin(); !itOrigin.IsAtEnd(); ++itOrigin)
	{
		vectorPixel = itOrigin.Get();
		idx = itOrigin.GetIndex();
		idx_x = idx.GetElement(0);
		idx_y = idx.GetElement(1);
		idx_z = idx.GetElement(2);

		//pixel = itOrigin.Get();
		//idx = itOrigin.GetIndex();
		originImageRawX[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(0);
		originImageRawY[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(1);
		originImageRawZ[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(2);
	}


	// resample image
	for (int k = 0; k < im_zn; k++) 
	{
		for (int j = 0; j < im_yn; j++)
		{
			for (int i = 0; i < im_xn; i++)
			{
				sampledImageRawX[k][j][i] = originImageRawX[k*rz][j*ry][i*rx]/rx;
				sampledImageRawY[k][j][i] = originImageRawY[k*rz][j*ry][i*rx]/ry;
				sampledImageRawZ[k][j][i] = originImageRawZ[k*rz][j*ry][i*rx]/rz;
			}
		}
	}

	// create the resampled image
	DeformationFieldType::Pointer sampledImage = DeformationFieldType::New();
	DeformationFieldType::IndexType start;
	start[0] = 0;start[1] = 0;start[2] = 0;
	DeformationFieldType::SizeType size;
	size[0] = im_xn;size[1] = im_yn;size[2] = im_zn;
	DeformationFieldType::RegionType region;
	region.SetSize (size);
	region.SetIndex (start);

	sampledImage -> SetRegions (region);
	// sampledImage -> SetRegions (originImage->GetLargestPossibleRegion());
	sampledImage -> SetSpacing (originImage->GetSpacing());
	sampledImage -> SetDirection (originImage->GetDirection());
	sampledImage -> SetOrigin (originImage->GetOrigin());

	sampledImage -> Allocate();

	DeformationFieldIteratorType itSampled(sampledImage, sampledImage->GetLargestPossibleRegion() );
	for (itSampled.GoToBegin(); !itSampled.IsAtEnd(); ++itSampled)
	{
		idx = itSampled.GetIndex();
		vectorPixel.SetElement(0, sampledImageRawX[idx[2]][idx[1]][idx[0]]);
		vectorPixel.SetElement(1, sampledImageRawY[idx[2]][idx[1]][idx[0]]);
		vectorPixel.SetElement(2, sampledImageRawZ[idx[2]][idx[1]][idx[0]]);
		itSampled.Set(vectorPixel);
	}

	WriteDeformationField(resampledDeformationFieldFileName, sampledImage);
	if (isCompressed)
		CompressDeformationField2Short(resampledDeformationFieldFileName);
	// delete newed variables
	for (int k = 0; k < im_z; k++) 
	{
		for (int j = 0; j < im_y; j++)
		{
			delete[] originImageRawX[k][j];
			delete[] originImageRawY[k][j];
			delete[] originImageRawZ[k][j];
		}
		delete[] originImageRawX[k];
		delete[] originImageRawY[k];
		delete[] originImageRawZ[k];
	}
	delete[] originImageRawX;
	delete[] originImageRawY;
	delete[] originImageRawZ;
	for (int k = 0; k < im_zn; k++) 
	{
		for (int j = 0; j < im_yn; j++)
		{
			delete[] sampledImageRawX[k][j];
			delete[] sampledImageRawY[k][j];
			delete[] sampledImageRawZ[k][j];
		}
		delete[] sampledImageRawX[k];
		delete[] sampledImageRawY[k];
		delete[] sampledImageRawZ[k];
	}
	delete[] sampledImageRawX;
	delete[] sampledImageRawY;
	delete[] sampledImageRawZ;

	return;
}
void DownResampleDeformationFieldNotValue(char* deformationFieldFileName, char* resampledDeformationFieldFileName, int sampleRate)
{
	// 
	int rx, ry, rz;
	rx = sampleRate;
	ry = sampleRate;
	rz = sampleRate;

	DeformationFieldType::Pointer originImage = 0;
	ReadDeformationField(deformationFieldFileName, originImage);

	int im_x, im_y, im_z;
	int im_xn, im_yn, im_zn;
	DeformationFieldType::SizeType im_size = originImage->GetLargestPossibleRegion().GetSize();
	im_x = im_size[0]; im_xn = im_x/rx;
	im_y = im_size[1]; im_yn = im_y/ry;
	im_z = im_size[2]; im_zn = im_z/rz;

	// create an array to store original image
	float*** originImageRawX = new float**[im_z];
	float*** originImageRawY = new float**[im_z];
	float*** originImageRawZ = new float**[im_z];
	for (int k = 0; k < im_z; k++) 
	{
		originImageRawX[k] = new float*[im_y];
		originImageRawY[k] = new float*[im_y];
		originImageRawZ[k] = new float*[im_y];
		for (int j = 0; j < im_y; j++)
		{
			originImageRawX[k][j] = new float[im_x];
			originImageRawY[k][j] = new float[im_x];
			originImageRawZ[k][j] = new float[im_x];
			for (int i = 0; i < im_x; i++)
			{
				originImageRawX[k][j][i] = 0.0;
				originImageRawY[k][j][i] = 0.0;
				originImageRawZ[k][j][i] = 0.0;
			}
		}
	}

	// create an array to store sampled image
	float*** sampledImageRawX = new float**[im_zn];
	float*** sampledImageRawY = new float**[im_zn];
	float*** sampledImageRawZ = new float**[im_zn];
	for (int k = 0; k < im_zn; k++) 
	{
		sampledImageRawX[k] = new float*[im_yn];
		sampledImageRawY[k] = new float*[im_yn];
		sampledImageRawZ[k] = new float*[im_yn];
		for (int j = 0; j < im_yn; j++)
		{
			sampledImageRawX[k][j] = new float[im_xn];
			sampledImageRawY[k][j] = new float[im_xn];
			sampledImageRawZ[k][j] = new float[im_xn];
			for (int i = 0; i < im_xn; i++)
			{
				sampledImageRawX[k][j][i] = 0.0;
				sampledImageRawY[k][j][i] = 0.0;
				sampledImageRawZ[k][j][i] = 0.0;
			}
		}
	}

	// load original image
	DeformationFieldIteratorType itOrigin(originImage, originImage->GetLargestPossibleRegion() );
	VectorPixelType vectorPixel;
	DeformationFieldType::IndexType idx;
	int idx_x, idx_y, idx_z;

	for (itOrigin.GoToBegin(); !itOrigin.IsAtEnd(); ++itOrigin)
	{
		vectorPixel = itOrigin.Get();
		idx = itOrigin.GetIndex();
		idx_x = idx.GetElement(0);
		idx_y = idx.GetElement(1);
		idx_z = idx.GetElement(2);

		//pixel = itOrigin.Get();
		//idx = itOrigin.GetIndex();
		originImageRawX[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(0);
		originImageRawY[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(1);
		originImageRawZ[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(2);
	}


	// resample image
	for (int k = 0; k < im_zn; k++) 
	{
		for (int j = 0; j < im_yn; j++)
		{
			for (int i = 0; i < im_xn; i++)
			{
				//sampledImageRawX[k][j][i] = originImageRawX[k*rz][j*ry][i*rx]/rx;
				//sampledImageRawY[k][j][i] = originImageRawY[k*rz][j*ry][i*rx]/ry;
				//sampledImageRawZ[k][j][i] = originImageRawZ[k*rz][j*ry][i*rx]/rz;
				sampledImageRawX[k][j][i] = originImageRawX[k*rz][j*ry][i*rx];
				sampledImageRawY[k][j][i] = originImageRawY[k*rz][j*ry][i*rx];
				sampledImageRawZ[k][j][i] = originImageRawZ[k*rz][j*ry][i*rx];
			}
		}
	}

	// create the resampled image
	DeformationFieldType::Pointer sampledImage = DeformationFieldType::New();
	DeformationFieldType::IndexType start;
	start[0] = 0;start[1] = 0;start[2] = 0;
	DeformationFieldType::SizeType size;
	size[0] = im_xn;size[1] = im_yn;size[2] = im_zn;
	DeformationFieldType::RegionType region;
	region.SetSize (size);
	region.SetIndex (start);

	sampledImage -> SetRegions (region);
	// sampledImage -> SetRegions (originImage->GetLargestPossibleRegion());
	sampledImage -> SetSpacing (originImage->GetSpacing());
	sampledImage -> SetDirection (originImage->GetDirection());
	sampledImage -> SetOrigin (originImage->GetOrigin());

	sampledImage -> Allocate();

	DeformationFieldIteratorType itSampled(sampledImage, sampledImage->GetLargestPossibleRegion() );
	for (itSampled.GoToBegin(); !itSampled.IsAtEnd(); ++itSampled)
	{
		idx = itSampled.GetIndex();
		vectorPixel.SetElement(0, sampledImageRawX[idx[2]][idx[1]][idx[0]]);
		vectorPixel.SetElement(1, sampledImageRawY[idx[2]][idx[1]][idx[0]]);
		vectorPixel.SetElement(2, sampledImageRawZ[idx[2]][idx[1]][idx[0]]);
		itSampled.Set(vectorPixel);
	}

	WriteDeformationField(resampledDeformationFieldFileName, sampledImage);

	// delete newed variables
	for (int k = 0; k < im_z; k++) 
	{
		for (int j = 0; j < im_y; j++)
		{
			delete[] originImageRawX[k][j];
			delete[] originImageRawY[k][j];
			delete[] originImageRawZ[k][j];
		}
		delete[] originImageRawX[k];
		delete[] originImageRawY[k];
		delete[] originImageRawZ[k];
	}
	delete[] originImageRawX;
	delete[] originImageRawY;
	delete[] originImageRawZ;
	for (int k = 0; k < im_zn; k++) 
	{
		for (int j = 0; j < im_yn; j++)
		{
			delete[] sampledImageRawX[k][j];
			delete[] sampledImageRawY[k][j];
			delete[] sampledImageRawZ[k][j];
		}
		delete[] sampledImageRawX[k];
		delete[] sampledImageRawY[k];
		delete[] sampledImageRawZ[k];
	}
	delete[] sampledImageRawX;
	delete[] sampledImageRawY;
	delete[] sampledImageRawZ;

	return;
}
void UpResampleDeformationField(char* deformationFieldFileName, char* resampledDeformationFieldFileName, int sampleRate)
{
	// upsample
	int rx, ry, rz;
	rx = sampleRate;
	ry = sampleRate;
	rz = sampleRate;

	DeformationFieldType::Pointer originImage = 0;
	if (!isCompressed)
		ReadDeformationField(deformationFieldFileName, originImage);
	else
		DeCompressDeformationFieldFromShort(deformationFieldFileName, originImage);

	int im_x, im_y, im_z;
	int im_xn, im_yn, im_zn;
	DeformationFieldType::SizeType im_size = originImage->GetLargestPossibleRegion().GetSize();
	im_x = im_size[0];
	im_y = im_size[1];
	im_z = im_size[2];
	//im_xn = im_x*rx; im_yn = im_y*ry; im_zn = im_z*rz;
	im_xn = imx; im_yn = imy; im_zn = imz;

	// create an array to store original image
	float*** originImageRawX = new float**[im_z];
	float*** originImageRawY = new float**[im_z];
	float*** originImageRawZ = new float**[im_z];
	for (int k = 0; k < im_z; k++) 
	{
		originImageRawX[k] = new float*[im_y];
		originImageRawY[k] = new float*[im_y];
		originImageRawZ[k] = new float*[im_y];
		for (int j = 0; j < im_y; j++)
		{
			originImageRawX[k][j] = new float[im_x];
			originImageRawY[k][j] = new float[im_x];
			originImageRawZ[k][j] = new float[im_x];
			for (int i = 0; i < im_x; i++)
			{
				originImageRawX[k][j][i] = 0.0;
				originImageRawY[k][j][i] = 0.0;
				originImageRawZ[k][j][i] = 0.0;
			}
		}
	}

	// create an array to store sampled image
	float*** sampledImageRawX = new float**[im_zn];
	float*** sampledImageRawY = new float**[im_zn];
	float*** sampledImageRawZ = new float**[im_zn];
	for (int k = 0; k < im_zn; k++) 
	{
		sampledImageRawX[k] = new float*[im_yn];
		sampledImageRawY[k] = new float*[im_yn];
		sampledImageRawZ[k] = new float*[im_yn];
		for (int j = 0; j < im_yn; j++)
		{
			sampledImageRawX[k][j] = new float[im_xn];
			sampledImageRawY[k][j] = new float[im_xn];
			sampledImageRawZ[k][j] = new float[im_xn];
			for (int i = 0; i < im_xn; i++)
			{
				sampledImageRawX[k][j][i] = 0.0;
				sampledImageRawY[k][j][i] = 0.0;
				sampledImageRawZ[k][j][i] = 0.0;
			}
		}
	}

	// load original image
	DeformationFieldIteratorType itOrigin(originImage, originImage->GetLargestPossibleRegion() );
	VectorPixelType vectorPixel;
	DeformationFieldType::IndexType idx;
	int idx_x, idx_y, idx_z;

	for (itOrigin.GoToBegin(); !itOrigin.IsAtEnd(); ++itOrigin)
	{
		vectorPixel = itOrigin.Get();
		idx = itOrigin.GetIndex();
		idx_x = idx.GetElement(0);
		idx_y = idx.GetElement(1);
		idx_z = idx.GetElement(2);

		//pixel = itOrigin.Get();
		//idx = itOrigin.GetIndex();
		originImageRawX[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(0);
		originImageRawY[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(1);
		originImageRawZ[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(2);
	}


	// resample image
	for (int k = 0; k < im_zn; k++) 
	{
		for (int j = 0; j < im_yn; j++)
		{
			for (int i = 0; i < im_xn; i++)
			{
				sampledImageRawX[k][j][i] = originImageRawX[k/rz][j/ry][i/rx]*rx;
				sampledImageRawY[k][j][i] = originImageRawY[k/rz][j/ry][i/rx]*ry;
				sampledImageRawZ[k][j][i] = originImageRawZ[k/rz][j/ry][i/rx]*rz;
			}
		}
	}

	// create the resampled image
	DeformationFieldType::Pointer sampledImage = DeformationFieldType::New();
	DeformationFieldType::IndexType start;
	start[0] = 0;start[1] = 0;start[2] = 0;
	DeformationFieldType::SizeType size;
	size[0] = im_xn;size[1] = im_yn;size[2] = im_zn;
	DeformationFieldType::RegionType region;
	region.SetSize (size);
	region.SetIndex (start);

	sampledImage -> SetRegions (region);
	// sampledImage -> SetRegions (originImage->GetLargestPossibleRegion());
	sampledImage -> SetSpacing (originImage->GetSpacing());
	sampledImage -> SetDirection (originImage->GetDirection());
	sampledImage -> SetOrigin (originImage->GetOrigin());

	sampledImage -> Allocate();

	DeformationFieldIteratorType itSampled(sampledImage, sampledImage->GetLargestPossibleRegion() );
	for (itSampled.GoToBegin(); !itSampled.IsAtEnd(); ++itSampled)
	{
		idx = itSampled.GetIndex();
		vectorPixel.SetElement(0, sampledImageRawX[idx[2]][idx[1]][idx[0]]);
		vectorPixel.SetElement(1, sampledImageRawY[idx[2]][idx[1]][idx[0]]);
		vectorPixel.SetElement(2, sampledImageRawZ[idx[2]][idx[1]][idx[0]]);
		itSampled.Set(vectorPixel);
	}

	WriteDeformationField(resampledDeformationFieldFileName, sampledImage);

	if (isCompressed)
		CompressDeformationField2Short(resampledDeformationFieldFileName);
	// delete newed variables
	for (int k = 0; k < im_z; k++) 
	{
		for (int j = 0; j < im_y; j++)
		{
			delete[] originImageRawX[k][j];
			delete[] originImageRawY[k][j];
			delete[] originImageRawZ[k][j];
		}
		delete[] originImageRawX[k];
		delete[] originImageRawY[k];
		delete[] originImageRawZ[k];
	}
	delete[] originImageRawX;
	delete[] originImageRawY;
	delete[] originImageRawZ;
	for (int k = 0; k < im_zn; k++) 
	{
		for (int j = 0; j < im_yn; j++)
		{
			delete[] sampledImageRawX[k][j];
			delete[] sampledImageRawY[k][j];
			delete[] sampledImageRawZ[k][j];
		}
		delete[] sampledImageRawX[k];
		delete[] sampledImageRawY[k];
		delete[] sampledImageRawZ[k];
	}
	delete[] sampledImageRawX;
	delete[] sampledImageRawY;
	delete[] sampledImageRawZ;


	return;
}

void LoadIntoArray(char* resampledDeformationFieldFileName, float* df_vector)
{
	// read deformation field and load it into array
	DeformationFieldType::Pointer dfImage = 0;
	if (!isCompressed)
		ReadDeformationField(resampledDeformationFieldFileName, dfImage);
	else
		DeCompressDeformationFieldFromShort(resampledDeformationFieldFileName, dfImage);

	int im_x, im_y, im_z;
	DeformationFieldType::SizeType im_size = dfImage->GetLargestPossibleRegion().GetSize();
	im_x = im_size[0];
	im_y = im_size[1];
	im_z = im_size[2];

//	int dim = 3*im_x*im_y*im_z;

	// load original image
	DeformationFieldIteratorType itOrigin(dfImage, dfImage->GetLargestPossibleRegion() );
	VectorPixelType vectorPixel;
	//DeformationFieldType::IndexType idx;
	//int idx_x, idx_y, idx_z;
	int index = 0;
	for (itOrigin.GoToBegin(); !itOrigin.IsAtEnd(); ++itOrigin)
	{
		vectorPixel = itOrigin.Get();

		df_vector[index] = vectorPixel.GetElement(0);
		df_vector[index+1] = vectorPixel.GetElement(1);
		df_vector[index+2] = vectorPixel.GetElement(2);
		index += 3;
	}


	return;
}
void SaveFromArray(char* deformationFieldFileName, float* df_vector, int sx, int sy, int sz)
{
	// save deformation field array into file
	// create the resampled image
	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
	DeformationFieldType::IndexType start;
	start[0] = 0;start[1] = 0;start[2] = 0;
	DeformationFieldType::SizeType size;
	size[0] = sx;size[1] = sy;size[2] = sz;
	DeformationFieldType::RegionType region;
	region.SetSize (size);
	region.SetIndex (start);

	deformationField -> SetRegions (region);
	// sampledImage -> SetRegions (originImage->GetLargestPossibleRegion());
	//deformationField -> SetSpacing (originImage->GetSpacing());
	//deformationField -> SetDirection (originImage->GetDirection());
	//deformationField -> SetOrigin (originImage->GetOrigin());
	deformationField -> SetSpacing (df_spacing);
	deformationField -> SetDirection (df_direction);
	deformationField -> SetOrigin (df_origin);

	deformationField -> Allocate();

	DeformationFieldIteratorType itDF(deformationField, deformationField->GetLargestPossibleRegion() );
	VectorPixelType vectorPixel;
	int index = 0;
	for (itDF.GoToBegin(); !itDF.IsAtEnd(); ++itDF)
	{
		vectorPixel.SetElement(0, df_vector[index]);
		vectorPixel.SetElement(1, df_vector[index+1]);
		vectorPixel.SetElement(2, df_vector[index+2]);
		itDF.Set(vectorPixel);
		index += 3;
	}

	WriteDeformationField(deformationFieldFileName, deformationField);
	if (isCompressed)
		CompressDeformationField2Short(deformationFieldFileName);

	return;
}
void CompressDeformationField(char* filename)
{
	DeformationFieldType::Pointer deformationfield = 0;
	ReadDeformationField(filename, deformationfield);
	
	DeformationFieldWriterType::Pointer deformationFieldWriter = DeformationFieldWriterType::New();
	deformationFieldWriter->SetFileName( filename );
	deformationFieldWriter->SetInput (deformationfield);
	deformationFieldWriter->SetUseCompression( true );
	try
	{
		deformationFieldWriter->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << err << std::endl; 
		return;
	} 
	return;
}
void DeCompressDeformationField(char* filename)
{
	DeformationFieldType::Pointer deformationfield = 0;
	ReadDeformationField(filename, deformationfield);

	DeformationFieldWriterType::Pointer deformationFieldWriter = DeformationFieldWriterType::New();
	deformationFieldWriter->SetFileName( filename );
	deformationFieldWriter->SetInput (deformationfield);
	deformationFieldWriter->SetUseCompression( false );
	try
	{
		deformationFieldWriter->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << err << std::endl; 
		return;
	} 
	return;
}

void CompressDeformationField2Short(char* filename)
{
	DeformationFieldType::Pointer inputDeformationField = 0;
	ReadDeformationField(filename, inputDeformationField);

	// multiply it with 1000
	MultiplyDeformationFieldFilterType::Pointer multiplyfilter = MultiplyDeformationFieldFilterType::New();
	multiplyfilter->SetInput(inputDeformationField);
	multiplyfilter->SetConstant((float)1000);

	// save to short
	Internal2ShortDeformationFieldCastFilterType::Pointer caster = Internal2ShortDeformationFieldCastFilterType::New();
	caster->SetInput(multiplyfilter->GetOutput());

	ShortDeformationFieldWriterType::Pointer deformationFieldWriter = ShortDeformationFieldWriterType::New();
	deformationFieldWriter->SetFileName( filename );
	deformationFieldWriter->SetInput (caster->GetOutput());
	deformationFieldWriter->SetUseCompression( true );

	try
	{
		deformationFieldWriter->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << err << std::endl; 
		return;
	} 
	return;
}
void DeCompressDeformationFieldFromShort(char* filename, DeformationFieldType::Pointer &deformationfield)
{
	// load input deformation field
	DeformationFieldType::Pointer inputDeformationField = 0;
	ReadDeformationField(filename, inputDeformationField);

	// divide it with 1000
	DivideDeformationFieldFilterType::Pointer dividefilter = DivideDeformationFieldFilterType::New();
	dividefilter->SetInput(inputDeformationField);
	dividefilter->SetConstant((float)1000);

	try
	{
		dividefilter->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << err << std::endl; 
		return;
	} 
	deformationfield = dividefilter->GetOutput();
	deformationfield->DisconnectPipeline();
	return;
}

void ComposeDeformationFieldsAndSaveCompressed(char* inputDeformationFieldFileName,
									 char* deformationFieldFileName,
									 char* composedDeformationFieldFileName)
{
	DeformationFieldType::Pointer inputDeformationField = DeformationFieldType::New();
	if (!isCompressed)
		ReadDeformationField(inputDeformationFieldFileName, inputDeformationField);
	else
		DeCompressDeformationFieldFromShort(inputDeformationFieldFileName, inputDeformationField);

	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
	if (!isCompressed)
		ReadDeformationField(deformationFieldFileName, deformationField);
	else
		DeCompressDeformationFieldFromShort(deformationFieldFileName, deformationField);

	DeformationFieldType::Pointer composedDeformationField = DeformationFieldType::New();

	ComposeDeformationFields(inputDeformationField, deformationField, composedDeformationField);

	WriteDeformationField(composedDeformationFieldFileName, composedDeformationField);

	if (isCompressed)
		CompressDeformationField2Short(composedDeformationFieldFileName);
	return;
}

void ApplyDeformationFieldAndWriteWithTypeWithFileNames(char* movingImageFileName, 
														char* deformationFieldFileName, char* deformedImageFileName, bool isLinear)
{
	itk::ImageIOBase::Pointer imageIO;
	try
	{
		imageIO = itk::ImageIOFactory::CreateImageIO(movingImageFileName, itk::ImageIOFactory::ReadMode);
		if ( imageIO )
		{
			imageIO->SetFileName(movingImageFileName);
			imageIO->ReadImageInformation();
		}
		else
		{
			std::cout << "Could not read the image information of "<< movingImageFileName <<"." << std::endl;
			exit( EXIT_FAILURE );
		}
	}
	catch( itk::ExceptionObject& err )
	{
		std::cout << "Could not read the image information of "<< movingImageFileName <<"." << std::endl;
		std::cout << err << std::endl;
		exit( EXIT_FAILURE );
	}
	//{UNKNOWNCOMPONENTTYPE,UCHAR,CHAR,USHORT,SHORT,UINT,INT,ULONG,LONG,FLOAT,DOUBLE} IOComponentType;
	//                    0     1    2      3     4    5   6     7    8     9     10
	int input_type = imageIO->GetComponentType(); // 9:float, 10:double

	//
	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
	ReadDeformationField(deformationFieldFileName, deformationField);

	InternalImageType::Pointer movingImage = InternalImageType::New();
	ReadImage(movingImageFileName,movingImage);

	InternalImageType::Pointer deformedImage = InternalImageType::New();

	ApplyDeformationField(movingImage, deformationField, deformedImage, isLinear);

	if (input_type == 1) // UCHAR
	{		
		WriteImageUCHAR(deformedImageFileName, deformedImage);
	}
	else if (input_type == 4) // SHORT
	{
		WriteImageSHORT(deformedImageFileName, deformedImage);
	}
	else if (input_type == 6) // INT
	{
		WriteImageINT(deformedImageFileName, deformedImage);
	}
	else if (input_type == 9) // FLOAT
	{
		WriteImageFLOAT(deformedImageFileName, deformedImage);
	}
	else if (input_type == 10) // DOUBLE
	{
		WriteImageFLOAT(deformedImageFileName, deformedImage);
	}
	else
	{
		WriteImageFLOAT(deformedImageFileName, deformedImage);
	}

}
// end
// 
