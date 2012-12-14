

#include "itkJia3D.h"

// read and write deformation field
DeformationFieldType::Pointer ReadDeformationField(char* filename)
{
	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( filename );
	try
	{
		deformationFieldReader->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << err << std::endl; 
		return 0;
	} 
	return deformationFieldReader ->GetOutput();
}

void ReadDeformationField(char* filename, DeformationFieldType::Pointer *deformationfield)
{
	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( filename );
	try
	{
		deformationFieldReader->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << err << std::endl; 
		return;
	} 
	*deformationfield = deformationFieldReader->GetOutput();
	return;
}

void ReadDeformationField(char* filename, DeformationFieldType::Pointer &deformationfield)
{
	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( filename );
	try
	{
		deformationFieldReader->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << err << std::endl; 
		return;
	} 
	deformationfield = deformationFieldReader->GetOutput();
	return;
}

void WriteDeformationField(char* filename, DeformationFieldType::Pointer deformationfield)
{
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


// read and write image
InternalImageType::Pointer ReadImage(char *filename)
{
	InternalImageReaderType::Pointer internalImageReader = InternalImageReaderType::New();
	internalImageReader->SetFileName(filename);
	try
	{
		internalImageReader->Update();
	}
	catch(itk::ExceptionObject & err )
	{
		std::cerr << err << std::endl;
		return 0;
	}
	return internalImageReader->GetOutput();
}

void ReadImage(char *filename, InternalImageType::Pointer *image)
{
	InternalImageReaderType::Pointer internalImageReader = InternalImageReaderType::New();
	internalImageReader->SetFileName(filename);
	try
	{
		internalImageReader->Update();
	}
	catch(itk::ExceptionObject & err )
	{
		std::cerr << err << std::endl;
		return;
	}
	*image = internalImageReader->GetOutput();
	return;
}

void ReadImage(char *filename, InternalImageType::Pointer &image)
{
	InternalImageReaderType::Pointer internalImageReader = InternalImageReaderType::New();
	internalImageReader->SetFileName(filename);
	try
	{
		internalImageReader->Update();
	}
	catch(itk::ExceptionObject & err )
	{
		std::cerr << err << std::endl;
		return;
	}
	image = internalImageReader->GetOutput();
	return;
}

void ReadCharImage(char *filename, CharImageType::Pointer &image)
{
	CharImageReaderType::Pointer charImageReader = CharImageReaderType::New();
	charImageReader->SetFileName(filename);
	try
	{
		charImageReader->Update();
	}
	catch(itk::ExceptionObject & err )
	{
		std::cerr << err << std::endl;
		return;
	}
	image = charImageReader->GetOutput();
	return;
}
void WriteImage(char *filename, InternalImageType::Pointer image, char* outputType)
{

	//
	if (strcmp(outputType, "uchar") == 0)
	{
		WriteImageUCHAR(filename, image);
	}
	else if (strcmp(outputType, "short") == 0)
	{
		WriteImageSHORT(filename, image);
	}
	else if (strcmp(outputType, "int") == 0)
	{
		WriteImageINT(filename, image);
	}
	else if (strcmp(outputType, "float") == 0)
	{
		WriteImageFLOAT(filename, image);
	}
	else // default
	{
		WriteImage(filename, image);
	}
}
void WriteImage(char *filename, InternalImageType::Pointer image)
{
	Internal2CharCastFilterType::Pointer caster = Internal2CharCastFilterType::New();
	caster->SetInput(image);
	CharImageWriterType::Pointer writer = CharImageWriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(caster->GetOutput());
	writer->SetUseCompression( false );
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << err << std::endl;
		return;
	}
	return;
}

void WriteImageUCHAR(char *filename, InternalImageType::Pointer image)
{
	Internal2CharCastFilterType::Pointer caster = Internal2CharCastFilterType::New();
	caster->SetInput(image);
	CharImageWriterType::Pointer writer = CharImageWriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(caster->GetOutput());
	writer->SetUseCompression( false );
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << err << std::endl;
		return;
	}
	return;
}

void WriteImageINT(char *filename, InternalImageType::Pointer image)
{
	Internal2IntCastFilterType::Pointer caster = Internal2IntCastFilterType::New();
	caster->SetInput(image);
	IntImageWriterType::Pointer writer = IntImageWriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(caster->GetOutput());
	writer->SetUseCompression( false );
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << err << std::endl;
		return;
	}
	return;
}
void WriteImageSHORT(char *filename, InternalImageType::Pointer image)
{
	Internal2ShortCastFilterType::Pointer caster = Internal2ShortCastFilterType::New();
	caster->SetInput(image);
	ShortImageWriterType::Pointer writer = ShortImageWriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(caster->GetOutput());
	writer->SetUseCompression( false );
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << err << std::endl;
		return;
	}
	return;
}
void WriteImageFLOAT(char *filename, InternalImageType::Pointer image)
{
	Internal2FloatCastFilterType::Pointer caster = Internal2FloatCastFilterType::New();
	caster->SetInput(image);
	FloatImageWriterType::Pointer writer = FloatImageWriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(caster->GetOutput());
	writer->SetUseCompression( false );
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << err << std::endl;
		return;
	}
	return;
}
// read and write transform
void ReadTransform(char *filename, BSplineTransformType::Pointer &bsplineTransform)
{
	TransformFileReaderType::Pointer        transformFileReader = TransformFileReaderType::New();
	transformFileReader->SetFileName(filename);
	try
	{
		transformFileReader->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << err << std::endl;
		return;
	}

	TransformListType*   transformList = transformFileReader->GetTransformList();
	BSplineTransformType::Pointer bsplineTransform_loaded = BSplineTransformType::New();
	bsplineTransform_loaded->SetFixedParameters(transformList->front()->GetFixedParameters());

	BSplineTransformType::ParametersType bsplineParameters_loaded;
	bsplineParameters_loaded.set_size(bsplineTransform_loaded->GetNumberOfParameters());
	bsplineTransform_loaded->SetParameters(bsplineParameters_loaded);
	bsplineTransform_loaded->SetParametersByValue(transformList->front()->GetParameters());

	bsplineTransform = bsplineTransform_loaded;

	return;
}
void WriteTransform(char *filename, BSplineTransformType::Pointer bsplineTransform)
{
	TransformFileWriterType::Pointer  transformFileWriter = TransformFileWriterType::New();
	transformFileWriter->SetFileName(filename);
	transformFileWriter->SetPrecision(12);
	transformFileWriter->SetInput(bsplineTransform);

	try
	{
		transformFileWriter->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << err << std::endl;
		return;
	}
	return;
}
// compose deformation fields
void ComposeDeformationFields(DeformationFieldType::Pointer input, 
							  DeformationFieldType::Pointer deformationField, 
							  DeformationFieldType::Pointer &composedDeformationField)
{
	WarpVectorFilterType::Pointer  vectorWarper = WarpVectorFilterType::New();
	vectorWarper->SetInput( input );
	vectorWarper->SetDeformationField( deformationField );
	vectorWarper->SetOutputOrigin(deformationField->GetOrigin());
	vectorWarper->SetOutputSpacing(deformationField->GetSpacing());
	vectorWarper->SetOutputDirection(deformationField->GetDirection());

	AddImageFilterType::Pointer addImageSum = AddImageFilterType::New();
	addImageSum->SetInput1(vectorWarper->GetOutput());
	addImageSum->SetInput2(deformationField);
	try
	{
		addImageSum->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << err << std::endl;
		return;
	}
	composedDeformationField = addImageSum->GetOutput();
	composedDeformationField->DisconnectPipeline();
	return;

}

void ComposeDeformationFieldsAndSave(char* inputDeformationFieldFileName,
									 char* deformationFieldFileName,
									 char* composedDeformationFieldFileName)
{
	DeformationFieldType::Pointer inputDeformationField = DeformationFieldType::New();
	ReadDeformationField(inputDeformationFieldFileName, inputDeformationField);

	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
	ReadDeformationField(deformationFieldFileName, deformationField);

	DeformationFieldType::Pointer composedDeformationField = DeformationFieldType::New();

	ComposeDeformationFields(inputDeformationField, deformationField, composedDeformationField);

	WriteDeformationField(composedDeformationFieldFileName, composedDeformationField);

	return;
}
// compose bspline transforms
void ComposeBSplineTransforms(BSplineTransformType::Pointer input,
							  BSplineTransformType::Pointer tranform,
							  BSplineTransformType::Pointer &composedTransform)
{
	
	return;
}
void ComposeBSplineTransformsAndSave (char* inputTransformFileName,
									  char* transformFileName,
									  char* composedTransformFileName)
{
	return;
}
// Inverse deformation field by DG's algorithm
void InverseDeformationFieldDG(DeformationFieldType::Pointer deformationField, 
							   DeformationFieldType::Pointer &deformationFieldInverse)
{
	int SHIFT = 2;
	float OUTSIDE = 0;//100000.0;
	int samplenum = 5;
	float interval = 1.0/(2*samplenum+1);
	unsigned int i, j;
	int x, y;
	float ii, jj;
	float mdl_subvoxelx, mdl_subvoxely;
	float disp_subvoxelx, disp_subvoxely;
	mdl_subvoxelx = 0.0; mdl_subvoxely = 0.0; disp_subvoxelx = 0.0; disp_subvoxely = 0.0;
	DeformationFieldIteratorType dfNewIt ( deformationFieldInverse, deformationFieldInverse->GetRequestedRegion() );
	DeformationFieldIteratorType dfIt ( deformationField, deformationField->GetRequestedRegion() );

	// initialize two 2D float matrix to store deformation field
	unsigned int image_size = deformationField->GetRequestedRegion().GetSize()[0];
	float** dfx = new float*[image_size];
	for(i = 0; i<image_size; i++) dfx[i] = new float[image_size];
	float** dfy = new float*[image_size];
	for(i = 0; i<image_size; i++) dfy[i] = new float[image_size];

	// load deformationFieldBA into 2 2D matrix dfx and dfy
	VectorPixelType vectorPixel;
	for (i = 0, j = 0, dfIt.GoToBegin(); !dfIt.IsAtEnd(); ++dfIt)
	{
		vectorPixel = dfIt.Get();
		dfx[i][j] = vectorPixel.GetElement(1);
		dfy[i][j] = vectorPixel.GetElement(0);
		j++;
		if (j==image_size)
		{
			i++;
			j=0;
		}
		//if ((i==image_size-1)&&(j==image_size-1))
		//{
		//	std::cerr << "read end!" << std::endl;
		//}
	}

	// allocate some internal data matrices, weights matrix and	enlarged inverse df matrix with borders
	float** totalweights = new float*[image_size+2*SHIFT];
	for(i = 0; i<image_size+2*SHIFT; i++) totalweights[i] = new float[image_size+2*SHIFT];
	float** rdfbx = new float*[image_size+2*SHIFT];
	for(i = 0; i<image_size+2*SHIFT; i++) rdfbx[i] = new float[image_size+2*SHIFT];
	float** rdfby = new float*[image_size+2*SHIFT];
	for(i = 0; i<image_size+2*SHIFT; i++) rdfby[i] = new float[image_size+2*SHIFT];

	// initialize these matrices
	for (i = 0; i < image_size+2*SHIFT; i++)
	{	
		for (j = 0; j < image_size+2*SHIFT; j++)
		{	
			totalweights[i][j]=0.0;	
			rdfbx[i][j]=0.0;
			rdfby[i][j]=0.0;
		}
	}

	// estimating
	for (i = 0; i < image_size; i++)
	{
		// std::cerr << i << ", ";
		for (j = 0; j < image_size; j++)
		{
			// std::cerr << "(" << i << ", " << j << "), " ;
			for (x = -samplenum; x<=samplenum; x++)
			{
				for (y = -samplenum; y<=samplenum; y++)
				{
					mdl_subvoxelx = x*interval + i;
					mdl_subvoxely = y*interval + j;
					
					// call interpolateDisplacement
					{
						int ni,nj,nip1,njp1;
						float b,c,b1,c1;
						ni = (int)mdl_subvoxelx; nip1 = ni+1;
						nj = (int)mdl_subvoxely; njp1 = nj+1;
						if((ni>=0)&&(ni<(int)image_size-1)&&(nj>=0)&&(nj<(int)image_size-1))
						{
							b = mdl_subvoxelx-ni;b1 = 1.0-b;
							c = mdl_subvoxely-nj;c1 = 1.0-c;
							disp_subvoxelx = ( dfx[ni][nj]*(b1*c1) + 
											   dfx[nip1][nj]*(b*c1)+ 
											   dfx[ni][njp1]*(b1*c)+ 
											   dfx[nip1][njp1]*(b*c) )/( (b1*c1)+(b*c1)+(b1*c)+(b*c) ) ; 
							disp_subvoxely = ( dfy[ni][nj]*(b1*c1) + 
											   dfy[nip1][nj]*(b*c1)+ 
											   dfy[ni][njp1]*(b1*c)+ 
											   dfy[nip1][njp1]*(b*c) )/( (b1*c1)+(b*c1)+(b1*c)+(b*c) ) ; 

						}
						else if (((ni==(int)image_size-1)&&(nj>=0)&&(nj<(int)image_size-1))||((ni>=0)&&(ni<(int)image_size-1)&&(nj==(int)image_size-1)))
						{
							disp_subvoxelx = dfx[ni][nj];
							disp_subvoxely = dfy[ni][nj];
						}

					}

					ii = mdl_subvoxelx + disp_subvoxelx;
					jj = mdl_subvoxely + disp_subvoxely;
					
					// call iterativeEstimate
					{
						int ni,nj,nip1,njp1;
						float b,c,b1,c1,combined,weight;
						ii += SHIFT; jj += SHIFT;
						ni = (int)ii; nip1 = ni+1;
						nj = (int)jj; njp1 = nj+1;
						if(ni>=0 && ni<(int)image_size+2*SHIFT-1  &&  nj>=0 && nj<(int)image_size+2*SHIFT-1)
						{
							b = ii-ni ;        b1 = 1.0-b ;
							c = jj-nj ;        c1 = 1.0-c ;
							combined = (b1*c1)+(b*c1)+(b1*c)+(b*c) ;

							weight = (b1*c1)/combined ;
							totalweights[ni][nj]   += weight ;
							rdfbx[ni][nj] += disp_subvoxelx*weight ;
							rdfby[ni][nj] += disp_subvoxely*weight ;

							weight = (b*c1)/combined ;
							totalweights[nip1][nj]   += weight ;
							rdfbx[nip1][nj] += disp_subvoxelx*weight ;
							rdfby[nip1][nj] += disp_subvoxely*weight ;

							weight = (b1*c)/combined ;
							totalweights[ni][njp1]   += weight ;
							rdfbx[ni][njp1] += disp_subvoxelx*weight ;
							rdfby[ni][njp1] += disp_subvoxely*weight ;

							weight = (b*c)/combined ;
							totalweights[nip1][njp1]   += weight ;
							rdfbx[nip1][njp1] += disp_subvoxelx*weight ;
							rdfby[nip1][njp1] += disp_subvoxely*weight ;

						}

					}

				}// end for x
			}// end for y
		}// end for j
	}// end for i

	// allocate inverse deformation field
	float** rdfx = new float*[image_size];
	for(i = 0; i<image_size; i++) rdfx[i] = new float[image_size];
	float** rdfy = new float*[image_size];
	for(i = 0; i<image_size; i++) rdfy[i] = new float[image_size];

	// normalize the enlarged rdfb to rdf
	for (i = 0; i< image_size; i++)
	{
		for (j = 0; j < image_size; j++)
		{
			if (totalweights[i+SHIFT][j+SHIFT]>0)
			{
				rdfx[i][j] = rdfbx[i+SHIFT][j+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT]);
				rdfy[i][j] = rdfby[i+SHIFT][j+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT]);
			}
			else
			{
				rdfx[i][j] = OUTSIDE;
				rdfy[i][j] = OUTSIDE;
			}
		}
	}

	for (i = 0, j = 0, dfNewIt.GoToBegin(); !dfNewIt.IsAtEnd(); ++dfNewIt)
	{
		vectorPixel.SetElement(1,rdfx[i][j]);
		vectorPixel.SetElement(0,rdfy[i][j]);
		dfNewIt.Set(vectorPixel);
		j++;
		if (j==image_size)
		{
			i++;
			j=0;
		}
	}

	// free memory
	for(i = 0; i<image_size; i++) delete[] dfx[i]; delete[] dfx;
	for(i = 0; i<image_size; i++) delete[] dfy[i]; delete[] dfy;
	for(i = 0; i<image_size; i++) delete[] rdfx[i]; delete[] rdfx;
	for(i = 0; i<image_size; i++) delete[] rdfy[i]; delete[] rdfy;
	for(i = 0; i<image_size; i++) delete[] rdfbx[i]; delete[] rdfbx;
	for(i = 0; i<image_size; i++) delete[] rdfby[i]; delete[] rdfby;
	for(i = 0; i<image_size; i++) delete[] totalweights[i]; delete[] totalweights;

}

// apply deformation field on image and write deformed image
void ApplyDeformationFieldAndWriteWithFileNames(char* movingImageName, char* deformationFieldFileName, 
												char* deformedImageName, bool isLinearInterpolator)
{
	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
	ReadDeformationField(deformationFieldFileName, deformationField);

	InternalImageType::Pointer movingImage = InternalImageType::New();
	ReadImage(movingImageName,movingImage);

	InternalImageType::Pointer deformedImage = InternalImageType::New();

	ApplyDeformationField(movingImage, deformationField, deformedImage, isLinearInterpolator);

	WriteImage(deformedImageName, deformedImage);

}

void ApplyDeformationField(InternalImageType::Pointer movingImage, 
						   DeformationFieldType::Pointer deformationField,
						   InternalImageType::Pointer &deformedImage,
						   bool isLinearInterpolator)
{
	if (isLinearInterpolator)
	{
		InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();
		InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();

		warper->SetInput( movingImage );
		warper->SetInterpolator( interpolator );
		warper->SetOutputSpacing( movingImage->GetSpacing() );
		warper->SetOutputOrigin( movingImage->GetOrigin() );
		warper->SetOutputDirection( movingImage->GetDirection() );
		warper->SetDeformationField( deformationField );
		try
		{
			warper->Update();
		}
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << "ExceptionObject caught !" << std::endl; 
			std::cerr << err << std::endl; 
		} 
		deformedImage = warper->GetOutput();
		deformedImage->DisconnectPipeline();
	}
	else
	{
		InternalNNInterpolatorType::Pointer interpolator = InternalNNInterpolatorType::New();
		InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();

		warper->SetInput( movingImage );
		warper->SetInterpolator( interpolator );
		warper->SetOutputSpacing( movingImage->GetSpacing() );
		warper->SetOutputOrigin( movingImage->GetOrigin() );
		warper->SetOutputDirection( movingImage->GetDirection() );
		warper->SetDeformationField( deformationField );
		try
		{
			warper->Update();
		}
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << "ExceptionObject caught !" << std::endl; 
			std::cerr << err << std::endl; 
		} 
		deformedImage = warper->GetOutput();
		deformedImage->DisconnectPipeline();
	}

}

// apply bspline transform on the moving image
void ApplyTransformAndWriteWithFileNames(char* movingImageName, char* transformFileName, 
												char* deformedImageName, bool isLinearInterpolator)
{
	BSplineTransformType::Pointer bsplineTransform = BSplineTransformType::New();
	ReadTransform(transformFileName, bsplineTransform);

	InternalImageType::Pointer movingImage = InternalImageType::New();
	ReadImage(movingImageName,movingImage);

	InternalImageType::Pointer deformedImage = InternalImageType::New();

	ApplyBSplineTransform(movingImage, bsplineTransform, deformedImage, isLinearInterpolator);

	WriteImage(deformedImageName, deformedImage);

}

void ApplyBSplineTransform(InternalImageType::Pointer movingImage,
						   BSplineTransformType::Pointer bsplineTransform,
						   InternalImageType::Pointer &deformedImage,
						   bool isLinearInterpolator)
{
	InternalLinearInterpolatorType::Pointer LinearInterpolator = InternalLinearInterpolatorType::New();
	InternalNNInterpolatorType::Pointer nnInterpolator = InternalNNInterpolatorType::New();

	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetTransform( bsplineTransform );
	resampler->SetInput( movingImage );
	resampler->SetSize( movingImage->GetLargestPossibleRegion().GetSize() );
	resampler->SetOutputOrigin(  movingImage->GetOrigin() );
	resampler->SetOutputSpacing( movingImage->GetSpacing() );
	resampler->SetOutputDirection( movingImage->GetDirection() );
	if (isLinearInterpolator)
		resampler->SetInterpolator(LinearInterpolator);
	else
		resampler->SetInterpolator(nnInterpolator);
	resampler->SetDefaultPixelValue( 0 );

	try
	{    
		resampler->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << "ExceptionObject caught!" << std::endl; 
		std::cerr << err << std::endl; 
		return;
	} 
	deformedImage = resampler->GetOutput();
	return;
}
// calculate weighted average deformation field
void calculateWeightedAverageDeformationField(DeformationFieldType::Pointer * allDeformationFields, double* weight, 
											  int totalNumber, DeformationFieldType::Pointer &weightedAveragedeformationField)
{
	int j2 = 0;
	double weightSum = 0;
	int numberOfDeformationLoaded = 0;
	double adjust_ratio = 1.0;

	// 1st deformation field
	MultiplyImageFilterType::Pointer multiplyImageFilter1 = MultiplyImageFilterType::New();
	multiplyImageFilter1->SetInput(allDeformationFields[j2]);
	multiplyImageFilter1->SetConstant((float)weight[j2]);
	weightSum += weight[j2];
	numberOfDeformationLoaded += 1;

	DivideImageFilterType::Pointer divideImage = DivideImageFilterType::New();
	DeformationFieldType::Pointer sumImage = DeformationFieldType::New();

	j2+=1;
	if (j2 < totalNumber)
	{
		AddImageFilterType::Pointer addImageSum = AddImageFilterType::New();

		MultiplyImageFilterType::Pointer multiplyImageFilter2 = MultiplyImageFilterType::New();
		multiplyImageFilter2->SetInput(allDeformationFields[j2]);
		multiplyImageFilter2->SetConstant((float)weight[j2]);
		weightSum += weight[j2];
		numberOfDeformationLoaded += 1;

		addImageSum->SetInput1(multiplyImageFilter1->GetOutput());
		addImageSum->SetInput2(multiplyImageFilter2->GetOutput());
		try
		{
			addImageSum->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught at adding deformation field" << std::endl;
			std::cerr << excep << std::endl;
		}

		for (int j = j2+1; j < totalNumber; j++) 
		{
			MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
			multiplyImageFilter->SetInput(allDeformationFields[j]);
			multiplyImageFilter->SetConstant((float)weight[j]);
			weightSum += weight[j];
			numberOfDeformationLoaded += 1;

			addImageSum->SetInput1(addImageSum->GetOutput());
			addImageSum->SetInput2(multiplyImageFilter->GetOutput());
			try
			{
				addImageSum->Update();
			}
			catch( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught at adding deformation field" << std::endl;
				std::cerr << excep << std::endl;
			}
		}

		sumImage = addImageSum->GetOutput();
		sumImage -> DisconnectPipeline();
		divideImage->SetConstant((float)weightSum);
		divideImage->SetInput(sumImage);

	}
	else
	{
		divideImage->SetConstant((float)weightSum);
		divideImage->SetInput(multiplyImageFilter1->GetOutput());
	}

	MultiplyImageFilterType::Pointer multiplyImage = MultiplyImageFilterType::New();
	multiplyImage->SetConstant((float)adjust_ratio);
	multiplyImage->SetInput(divideImage->GetOutput());

	try
	{
		multiplyImage->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at dividing deformation field" << std::endl;
		std::cerr << excep << std::endl;
	}
	weightedAveragedeformationField = multiplyImage->GetOutput();
	weightedAveragedeformationField -> DisconnectPipeline();
	
	return;
}
// Diffeomorphic Demons Registration, some need to make changes to use sub-functions
void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  char* deformedImageFileName, char* deformationFieldFileName)
{
	// read fixed and moving images
	InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
	InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

	fixedImageReader->SetFileName( fixedImageFileName );
	movingImageReader->SetFileName( movingImageFileName );

	// do histogram matching and some parameters need to set manually
	InternalHistMatchFilterType::Pointer matcher = InternalHistMatchFilterType::New();
	matcher->SetInput( movingImageReader->GetOutput() );
	matcher->SetReferenceImage( fixedImageReader->GetOutput() );
	matcher->SetNumberOfHistogramLevels( 1024 ); // 1024
	matcher->SetNumberOfMatchPoints( 7 );
	matcher->ThresholdAtMeanIntensityOn();

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
	float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

	multires->SetRegistrationFilter( filter );
	multires->SetNumberOfLevels( curNumIterations.size() );
	multires->SetNumberOfIterations( &curNumIterations[0] );
	multires->SetFixedImage( fixedImageReader->GetOutput() );
	multires->SetMovingImage( matcher->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );
	
	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  char* deformedImageFileName, char* deformationFieldFileName,
							  unsigned int * iter, int res)
{
	// read fixed and moving images
	InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
	InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

	fixedImageReader->SetFileName( fixedImageFileName );
	movingImageReader->SetFileName( movingImageFileName );

	// do histogram matching and some parameters need to set manually
	InternalHistMatchFilterType::Pointer matcher = InternalHistMatchFilterType::New();
	matcher->SetInput( movingImageReader->GetOutput() );
	matcher->SetReferenceImage( fixedImageReader->GetOutput() );
	matcher->SetNumberOfHistogramLevels( 1024 ); // 1024
	matcher->SetNumberOfMatchPoints( 7 );
	matcher->ThresholdAtMeanIntensityOn();

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
	float sigmaDef = 1.5;
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
	for (int i = 0; i < res; i++)		curNumIterations.push_back(iter[i]);

	multires->SetRegistrationFilter( filter );
	multires->SetNumberOfLevels( curNumIterations.size() );
	multires->SetNumberOfIterations( &curNumIterations[0] );
	multires->SetFixedImage( fixedImageReader->GetOutput() );
	multires->SetMovingImage( matcher->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );

	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistrationWithoutHistMatching(char* fixedImageFileName, char* movingImageFileName,
							  char* deformedImageFileName, char* deformationFieldFileName)
{
	// read fixed and moving images
	InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
	InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

	fixedImageReader->SetFileName( fixedImageFileName );
	movingImageReader->SetFileName( movingImageFileName );

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
	float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

	multires->SetRegistrationFilter( filter );
	multires->SetNumberOfLevels( curNumIterations.size() );
	multires->SetNumberOfIterations( &curNumIterations[0] );
	multires->SetFixedImage( fixedImageReader->GetOutput() );
	multires->SetMovingImage( movingImageReader->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );
	
	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistrationWithoutHistMatching(char* fixedImageFileName, char* movingImageFileName,
												 char* deformedImageFileName, char* deformationFieldFileName,
												 unsigned int * iter, int res)
{
	// read fixed and moving images
	InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
	InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

	fixedImageReader->SetFileName( fixedImageFileName );
	movingImageReader->SetFileName( movingImageFileName );

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
	float sigmaDef = 1.5;
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
	for (int i = 0; i < res; i++)		curNumIterations.push_back(iter[i]);

	multires->SetRegistrationFilter( filter );
	multires->SetNumberOfLevels( curNumIterations.size() );
	multires->SetNumberOfIterations( &curNumIterations[0] );
	multires->SetFixedImage( fixedImageReader->GetOutput() );
	multires->SetMovingImage( movingImageReader->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );
	
	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  char* deformedImageFileName, char* deformationFieldFileName, bool doHistMatch)
{
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
	float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );
	
	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  DeformationFieldType::Pointer &deformationField, bool doHistMatch)
{
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
	float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

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

	deformationField = multires->GetOutput();

	//// write deformed image
	//InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	//InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	//warper->SetInput( movingImageReader->GetOutput() );
	//warper->SetInterpolator( interpolator );
	//warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	//warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	//warper->SetDeformationField( multires->GetOutput() );
	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	//// write deformation field into a file
	//DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	//fieldWriter->SetFileName( deformationFieldFileName );
	//fieldWriter->SetInput( multires->GetOutput() );
	//fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  DeformationFieldType::Pointer &deformationField, bool doHistMatch, float sigmaDef)
{
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
	// float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

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

	deformationField = multires->GetOutput();

	return;
}

void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName, char* preDeformationFieldFileName,
							  DeformationFieldType::Pointer &deformationField, bool doHistMatch, float sigmaDef)
{
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
	// float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

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

	// introduce previous deformation field file
	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( preDeformationFieldFileName );
	try 
	{
		deformationFieldReader->Update();  // do demons registration
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at reading deformation field!" << std::endl; 
		std::cerr << excep << std::endl;
	}
	multires->SetArbitraryInitialDeformationField( deformationFieldReader->GetOutput() );

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

	deformationField = multires->GetOutput();

	return;
}

void DiffeoDemonsRegistration(char* fixedImageFileName, char* movingImageFileName,
							  char* deformedImageFileName, char* deformationFieldFileName, bool doHistMatch,
							  unsigned int * iter, int res)
{
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
	float sigmaDef = 1.5;
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
	for (int i = 0; i < res; i++)		curNumIterations.push_back(iter[i]);

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );
	
	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFieldFileName,
														 char* deformedImageFileName, char* deformationFieldFileName)
{
	// read fixed and moving images
	InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
	InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

	fixedImageReader->SetFileName( fixedImageFileName );
	movingImageReader->SetFileName( movingImageFileName );

	// do histogram matching and some parameters need to set manually
	InternalHistMatchFilterType::Pointer matcher = InternalHistMatchFilterType::New();
	matcher->SetInput( movingImageReader->GetOutput() );
	matcher->SetReferenceImage( fixedImageReader->GetOutput() );
	matcher->SetNumberOfHistogramLevels( 1024 ); // 1024
	matcher->SetNumberOfMatchPoints( 7 );
	matcher->ThresholdAtMeanIntensityOn();

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
	float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

	multires->SetRegistrationFilter( filter );
	multires->SetNumberOfLevels( curNumIterations.size() );
	multires->SetNumberOfIterations( &curNumIterations[0] );
	multires->SetFixedImage( fixedImageReader->GetOutput() );
	multires->SetMovingImage( matcher->GetOutput() );

	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( initialDeformationFieldFileName );
	try 
	{
		deformationFieldReader->Update();  // do demons registration
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at reading deformation field!" << std::endl; 
		std::cerr << excep << std::endl;
	}
	multires->SetArbitraryInitialDeformationField( deformationFieldReader->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );

	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFieldFileName,
														 char* deformedImageFileName, char* deformationFieldFileName,
														 unsigned int * iter, int res)
{
	// read fixed and moving images
	InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
	InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

	fixedImageReader->SetFileName( fixedImageFileName );
	movingImageReader->SetFileName( movingImageFileName );

	// do histogram matching and some parameters need to set manually
	InternalHistMatchFilterType::Pointer matcher = InternalHistMatchFilterType::New();
	matcher->SetInput( movingImageReader->GetOutput() );
	matcher->SetReferenceImage( fixedImageReader->GetOutput() );
	matcher->SetNumberOfHistogramLevels( 1024 ); // 1024
	matcher->SetNumberOfMatchPoints( 7 );
	matcher->ThresholdAtMeanIntensityOn();

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
	float sigmaDef = 1.5;
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
	for (int i = 0; i < res; i++)		curNumIterations.push_back(iter[i]);

	multires->SetRegistrationFilter( filter );
	multires->SetNumberOfLevels( curNumIterations.size() );
	multires->SetNumberOfIterations( &curNumIterations[0] );
	multires->SetFixedImage( fixedImageReader->GetOutput() );
	multires->SetMovingImage( matcher->GetOutput() );

	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( initialDeformationFieldFileName );
	try 
	{
		deformationFieldReader->Update();  // do demons registration
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at reading deformation field!" << std::endl; 
		std::cerr << excep << std::endl;
	}
	multires->SetArbitraryInitialDeformationField( deformationFieldReader->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );

	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistrationWithInitialDeformationFieldWithoutHistMatching(char* fixedImageFileName, char* movingImageFileName, 
																			char* initialDeformationFieldFileName,
																			char* deformedImageFileName, char* deformationFieldFileName)
{
	// read fixed and moving images
	InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
	InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

	fixedImageReader->SetFileName( fixedImageFileName );
	movingImageReader->SetFileName( movingImageFileName );

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
	float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

	multires->SetRegistrationFilter( filter );
	multires->SetNumberOfLevels( curNumIterations.size() );
	multires->SetNumberOfIterations( &curNumIterations[0] );
	multires->SetFixedImage( fixedImageReader->GetOutput() );
	multires->SetMovingImage( movingImageReader->GetOutput() );

	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( initialDeformationFieldFileName );
	try 
	{
		deformationFieldReader->Update();  // do demons registration
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at reading deformation field!" << std::endl; 
		std::cerr << excep << std::endl;
	}
	multires->SetArbitraryInitialDeformationField( deformationFieldReader->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );

	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistrationWithInitialDeformationFieldWithoutHistMatching(char* fixedImageFileName, char* movingImageFileName, 
																			char* initialDeformationFieldFileName,
																			char* deformedImageFileName, char* deformationFieldFileName,
																			unsigned int * iter, int res)
{
	// read fixed and moving images
	InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
	InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

	fixedImageReader->SetFileName( fixedImageFileName );
	movingImageReader->SetFileName( movingImageFileName );

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
	float sigmaDef = 1.5;
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
	for (int i = 0; i < res; i++)		curNumIterations.push_back(iter[i]);

	multires->SetRegistrationFilter( filter );
	multires->SetNumberOfLevels( curNumIterations.size() );
	multires->SetNumberOfIterations( &curNumIterations[0] );
	multires->SetFixedImage( fixedImageReader->GetOutput() );
	multires->SetMovingImage( movingImageReader->GetOutput() );

	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( initialDeformationFieldFileName );
	try 
	{
		deformationFieldReader->Update();  // do demons registration
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at reading deformation field!" << std::endl; 
		std::cerr << excep << std::endl;
	}
	multires->SetArbitraryInitialDeformationField( deformationFieldReader->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );

	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFieldFileName,
														 char* deformedImageFileName, char* deformationFieldFileName,
														 bool doHistMatch)
{
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
	float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

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

	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( initialDeformationFieldFileName );
	try 
	{
		deformationFieldReader->Update();  // do demons registration
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at reading deformation field!" << std::endl; 
		std::cerr << excep << std::endl;
	}
	multires->SetArbitraryInitialDeformationField( deformationFieldReader->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );

	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFieldFileName,
														 char* deformedImageFileName, char* deformationFieldFileName,
														 bool doHistMatch, float sigmaDef)
{
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
	// float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

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

	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( initialDeformationFieldFileName );
	try 
	{
		deformationFieldReader->Update();  // do demons registration
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at reading deformation field!" << std::endl; 
		std::cerr << excep << std::endl;
	}
	multires->SetArbitraryInitialDeformationField( deformationFieldReader->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );

	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFieldFileName,
														 DeformationFieldType::Pointer &deformationField, 
														 bool doHistMatch)
{
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
	float sigmaDef = 1.5;
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
	unsigned int curNumOfIterations[] = {15,10,5};
	for (int i = 0; i < 3; i++)		curNumIterations.push_back(curNumOfIterations[i]);

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

	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( initialDeformationFieldFileName );
	try 
	{
		deformationFieldReader->Update();  // do demons registration
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at reading deformation field!" << std::endl; 
		std::cerr << excep << std::endl;
	}
	multires->SetArbitraryInitialDeformationField( deformationFieldReader->GetOutput() );

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

	deformationField = multires->GetOutput();

	//// write deformed image
	//InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	//InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	//warper->SetInput( movingImageReader->GetOutput() );
	//warper->SetInterpolator( interpolator );
	//warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	//warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	//warper->SetDeformationField( multires->GetOutput() );
	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	//// write deformation field into a file
	//DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	//fieldWriter->SetFileName( deformationFieldFileName );
	//fieldWriter->SetInput( multires->GetOutput() );
	//fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistrationWithInitialDeformationField(char* fixedImageFileName, char* movingImageFileName, 
														 char* initialDeformationFieldFileName,
														 char* deformedImageFileName, char* deformationFieldFileName,
														 bool doHistMatch, unsigned int * iter, int res)
{
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
	float sigmaDef = 1.5;
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
	for (int i = 0; i < res; i++)		curNumIterations.push_back(iter[i]);

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

	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( initialDeformationFieldFileName );
	try 
	{
		deformationFieldReader->Update();  // do demons registration
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught at reading deformation field!" << std::endl; 
		std::cerr << excep << std::endl;
	}
	multires->SetArbitraryInitialDeformationField( deformationFieldReader->GetOutput() );

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

	// write deformed image
	InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
	warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
	warper->SetDeformationField( multires->GetOutput() );

	warper->Update();
	InternalImageType::Pointer outputImage = warper->GetOutput();
	outputImage->DisconnectPipeline();
	WriteImage(deformedImageFileName, outputImage);

	//CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
	//Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
	//writer->SetFileName( deformedImageFileName );
	//caster->SetInput( warper->GetOutput() );
	//writer->SetInput( caster->GetOutput()   );
	//writer->Update();

	// write deformation field into a file
	DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
	fieldWriter->SetFileName( deformationFieldFileName );
	fieldWriter->SetInput( multires->GetOutput() );
	fieldWriter->Update();

	return;
}

void DiffeoDemonsRegistration(InternalImageType::Pointer fixedImage, 
							  InternalImageType::Pointer movingImage,
							  InternalImageType::Pointer deformedImage,
							  DeformationFieldType::Pointer &deformationField)
{
	return;
}
void DiffeoDemonsRegistration(InternalImageType::Pointer fixedImage, 
							  InternalImageType::Pointer movingImage,
							  DeformationFieldType::Pointer &deformationField)
{
	return;
}
// BSpline Registration
void BSplineRegistration(char* fixedImageFileName, char* movingImageFileName,
						 char* deformedImageFileName, char* transformFileName)
{
	InternalImageType::Pointer fixedImage = InternalImageType::New();
	ReadImage(fixedImageFileName, fixedImage);

	InternalImageType::Pointer movingImage = InternalImageType::New();
	ReadImage(movingImageFileName, movingImage);

	InternalImageType::Pointer matchedMovingImage = InternalImageType::New();
	HistogramMatching(movingImage, fixedImage, matchedMovingImage);

	InternalImageType::Pointer deformedImage = InternalImageType::New();

	BSplineTransformType::Pointer bsplineTransform = BSplineTransformType::New();

	BSplineRegistration(fixedImage, matchedMovingImage, deformedImage, bsplineTransform);

	WriteImage(deformedImageFileName, deformedImage);

	ReadTransform("transform.txt", bsplineTransform);

	WriteTransform(transformFileName, bsplineTransform);

	return;
}
void BSplineRegistration(InternalImageType::Pointer fixedImage, InternalImageType::Pointer movingImage,
						 InternalImageType::Pointer &deformedImage, BSplineTransformType::Pointer &bsplineTransform)
{
	InternalMeanSquaresMetricType::Pointer		    metric        = InternalMeanSquaresMetricType::New();
	LBFGSOptimizerType::Pointer						optimizer     = LBFGSOptimizerType::New();
	InternalLinearInterpolatorType::Pointer			interpolator  = InternalLinearInterpolatorType::New();
	ImageRegistrationType::Pointer					registration  = ImageRegistrationType::New();

	registration->SetMetric(        metric        );
	registration->SetOptimizer(     optimizer     );
	registration->SetInterpolator(  interpolator  );
	registration->SetFixedImage(  fixedImage   );
	registration->SetMovingImage(   movingImage  );

	InternalImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
	registration->SetFixedImageRegion( fixedRegion );


	BSplineTransformType::Pointer  transformLow = BSplineTransformType::New();
	registration->SetTransform( transformLow );

	BSplineTransformType::RegionType::SizeType   gridBorderSize;
	BSplineTransformType::RegionType::SizeType   gridLowSizeOnImage;
	BSplineTransformType::RegionType::SizeType   totalGridSize;

	gridBorderSize.Fill( gridBorderLength );    // Border for spline order gridBorderLength= 3 ( 1 lower, 2 upper )
	gridLowSizeOnImage.Fill( lowRes );
	totalGridSize = gridLowSizeOnImage + gridBorderSize;

	BSplineTransformType::RegionType bsplineRegion;
	bsplineRegion.SetSize( totalGridSize );

	BSplineTransformType::SpacingType spacingLow = fixedImage->GetSpacing();
	BSplineTransformType::OriginType originLow = fixedImage->GetOrigin();;

	InternalImageType::SizeType fixedImageSize = fixedRegion.GetSize();

	for(unsigned int r = 0; r < ImageDimension; r++)
	{ spacingLow[r] *= static_cast<double>(fixedImageSize[r] - 1) / static_cast<double>(gridLowSizeOnImage[r] - 1);}

	InternalImageType::DirectionType gridDirection = fixedImage->GetDirection();
	BSplineTransformType::SpacingType gridOriginOffsetLow = gridDirection * spacingLow;

	BSplineTransformType::OriginType gridOriginLow = originLow - gridOriginOffsetLow; 

	transformLow->SetGridSpacing( spacingLow );
	transformLow->SetGridOrigin( gridOriginLow );
	transformLow->SetGridRegion( bsplineRegion );
	transformLow->SetGridDirection( gridDirection );

	
	BSplineTransformType::ParametersType parametersLow( transformLow->GetNumberOfParameters() );
	parametersLow.Fill( 0.0 );
	transformLow->SetParameters( parametersLow );

	// start to register in low resolution
	registration->SetInitialTransformParameters( transformLow->GetParameters() );
	optimizer->SetGradientConvergenceTolerance( 0.05 );
	optimizer->SetLineSearchAccuracy( 0.9 );
	optimizer->SetDefaultStepLength( 1.5 );
	optimizer->TraceOff();
	optimizer->SetMaximumNumberOfFunctionEvaluations( 1000 );

	// std::cout << "Starting Registration with low resolution transform" << std::endl;

	try 
	{ 
		registration->StartRegistration(); 
	} 
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << "ExceptionObject caught at low resolution !" << std::endl; 
		std::cerr << err << std::endl; 
		return;
	} 

	// Finally we use the last transform parameters in order to resample the image.
	transformLow->SetParameters( registration->GetLastTransformParameters() );

	////////////////////////////////////////////////////////
	//  Once the registration has finished with the low resolution grid, we
	//  proceed to instantiate the medium resolution
	BSplineTransformType::Pointer  transformMedium = BSplineTransformType::New();
	BSplineTransformType::RegionType::SizeType   gridMediumSizeOnImage;
	gridMediumSizeOnImage.Fill( mediumRes );
	totalGridSize = gridMediumSizeOnImage + gridBorderSize;
	bsplineRegion.SetSize( totalGridSize );
	BSplineTransformType::SpacingType spacingMedium = fixedImage->GetSpacing();
	BSplineTransformType::OriginType  originMedium  = fixedImage->GetOrigin();;

	for(unsigned int rm = 0; rm < ImageDimension; rm++)
	{ spacingMedium[rm] *= static_cast<double>(fixedImageSize[rm] - 1) / static_cast<double>(gridMediumSizeOnImage[rm] - 1);}

	BSplineTransformType::SpacingType gridOriginOffsetMedium = gridDirection * spacingMedium;
	BSplineTransformType::OriginType gridOriginMedium = originMedium - gridOriginOffsetMedium; 

	transformMedium->SetGridSpacing( spacingMedium );
	transformMedium->SetGridOrigin( gridOriginMedium );
	transformMedium->SetGridRegion( bsplineRegion );
	transformMedium->SetGridDirection( gridDirection );
	BSplineTransformType::ParametersType parametersMedium( transformMedium->GetNumberOfParameters() );
	parametersMedium.Fill( 0.0 );

	int counter = 0;
	for (unsigned int k = 0; k < SpaceDimension; k++ )
	{
		typedef BSplineTransformType::ImageType ParametersImageType;
		typedef itk::ResampleImageFilter<ParametersImageType,ParametersImageType> ResamplerType;
		typedef itk::BSplineResampleImageFunction<ParametersImageType,double> FunctionType;
		typedef itk::IdentityTransform<double,SpaceDimension> IdentityTransformType;
		typedef itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType> DecompositionType;
		typedef itk::ImageRegionIterator<ParametersImageType> Iterator;

		ResamplerType::Pointer upsamplerMedium = ResamplerType::New();
		FunctionType::Pointer functionMedium = FunctionType::New();
		IdentityTransformType::Pointer identityMedium = IdentityTransformType::New();

		upsamplerMedium->SetInput( transformLow->GetCoefficientImage()[k] );
		upsamplerMedium->SetInterpolator( functionMedium );
		upsamplerMedium->SetTransform( identityMedium );
		upsamplerMedium->SetSize( transformMedium->GetGridRegion().GetSize() );
		upsamplerMedium->SetOutputSpacing( transformMedium->GetGridSpacing() );
		upsamplerMedium->SetOutputOrigin( transformMedium->GetGridOrigin() );
		upsamplerMedium->SetOutputDirection( fixedImage->GetDirection() );

		DecompositionType::Pointer decompositionMedium = DecompositionType::New();
		decompositionMedium->SetSplineOrder( SplineOrder );
		decompositionMedium->SetInput( upsamplerMedium->GetOutput() );
		try
		{
			decompositionMedium->Update();
		}
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << err << std::endl; 
			return;
		} 
		ParametersImageType::Pointer newCoefficientsMedium = decompositionMedium->GetOutput();

		// copy the coefficients into the parameter array
		Iterator it( newCoefficientsMedium, transformMedium->GetGridRegion() );
		double readdata = 0;
		while ( !it.IsAtEnd() )
		{ 
			readdata = it.Get();
			parametersMedium[ counter++ ] = it.Get();
			++it;
		}
	}

	// start registration on medium resolution  
	transformMedium->SetParameters( parametersMedium );
	// std::cout << "Starting Registration with medium resolution transform" << std::endl;
	registration->SetInitialTransformParameters( transformMedium->GetParameters() );
	registration->SetTransform( transformMedium );
	optimizer->SetGradientConvergenceTolerance( 0.01 );
	optimizer->SetDefaultStepLength( 1.0 );
	try 
	{ 
		registration->StartRegistration(); 
	} 
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << "ExceptionObject caught at medium resolution !" << std::endl; 
		std::cerr << err << std::endl; 
		return;
	} 

	// Finally we use the last transform parameters in order to resample the image.
	transformMedium->SetParameters( registration->GetLastTransformParameters() );

	////////////////////////////////////////////////
	//  Once the registration has finished with the medium resolution grid, we
	//  proceed to instantiate a higher resolution
	BSplineTransformType::Pointer  transformHigh = BSplineTransformType::New();
	BSplineTransformType::RegionType::SizeType   gridHighSizeOnImage;
	gridHighSizeOnImage.Fill( highRes );
	totalGridSize = gridHighSizeOnImage + gridBorderSize;
	bsplineRegion.SetSize( totalGridSize );
	BSplineTransformType::SpacingType spacingHigh = fixedImage->GetSpacing();
	BSplineTransformType::OriginType  originHigh  = fixedImage->GetOrigin();;

	for(unsigned int rh = 0; rh < ImageDimension; rh++)
	{ spacingHigh[rh] *= static_cast<double>(fixedImageSize[rh] - 1) / static_cast<double>(gridHighSizeOnImage[rh] - 1);}

	BSplineTransformType::SpacingType gridOriginOffsetHigh = gridDirection * spacingHigh;
	BSplineTransformType::OriginType gridOriginHigh = originHigh - gridOriginOffsetHigh; 

	transformHigh->SetGridSpacing( spacingHigh );
	transformHigh->SetGridOrigin( gridOriginHigh );
	transformHigh->SetGridRegion( bsplineRegion );
	transformHigh->SetGridDirection( gridDirection );
	BSplineTransformType::ParametersType parametersHigh( transformHigh->GetNumberOfParameters() );
	parametersHigh.Fill( 0.0 );

	//  Now we need to initialize the BSpline coefficients of the higher resolution
	//  transform. This is done by first computing the actual deformation field 
	//  at the higher resolution from the lower resolution BSpline coefficients. 
	//  Then a BSpline decomposition is done to obtain the BSpline coefficient of 
	//  the higher resolution transform.

	counter = 0;
	for (unsigned int k = 0; k < SpaceDimension; k++ )
	{
		typedef BSplineTransformType::ImageType ParametersImageType;
		typedef itk::ResampleImageFilter<ParametersImageType,ParametersImageType> ResamplerType;
		typedef itk::BSplineResampleImageFunction<ParametersImageType,double> FunctionType;
		typedef itk::IdentityTransform<double,SpaceDimension> IdentityTransformType;
		typedef itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType> DecompositionType;
		typedef itk::ImageRegionIterator<ParametersImageType> Iterator;

		ResamplerType::Pointer upsampler = ResamplerType::New();
		FunctionType::Pointer function = FunctionType::New();
		IdentityTransformType::Pointer identity = IdentityTransformType::New();

		upsampler->SetInput( transformMedium->GetCoefficientImage()[k] );
		upsampler->SetInterpolator( function );
		upsampler->SetTransform( identity );
		upsampler->SetSize( transformHigh->GetGridRegion().GetSize() );
		upsampler->SetOutputSpacing( transformHigh->GetGridSpacing() );
		upsampler->SetOutputOrigin( transformHigh->GetGridOrigin() );
		upsampler->SetOutputDirection( fixedImage->GetDirection() );

		DecompositionType::Pointer decompositionHigh = DecompositionType::New();
		decompositionHigh->SetSplineOrder( SplineOrder );
		decompositionHigh->SetInput( upsampler->GetOutput() );
		try
		{
			decompositionHigh->Update();
		}
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << err << std::endl; 
			return;
		} 

		ParametersImageType::Pointer newCoefficients = decompositionHigh->GetOutput();

		// copy the coefficients into the parameter array
		Iterator it( newCoefficients, transformHigh->GetGridRegion() );
		double readdata = 0;
		while ( !it.IsAtEnd() )
		{ 
			readdata = it.Get();
			parametersHigh[ counter++ ] = it.Get();
			++it;
		}
	}
	// start registration on high resolution  
	transformHigh->SetParameters( parametersHigh );
	// std::cout << "Starting Registration with high resolution transform" << std::endl;
	registration->SetInitialTransformParameters( transformHigh->GetParameters() );
	registration->SetTransform( transformHigh );
	optimizer->SetGradientConvergenceTolerance( 0.005 );
	optimizer->SetDefaultStepLength( 0.5 );
	try 
	{ 
		registration->StartRegistration(); 
	} 
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << "ExceptionObject caught at high resolution!" << std::endl; 
		std::cerr << err << std::endl; 
		return;
	} 
	// Finally we use the last transform parameters in order to resample the image.
	transformHigh->SetParameters( registration->GetLastTransformParameters() );

	bsplineTransform = transformHigh;

	ApplyBSplineTransform(movingImage, bsplineTransform, deformedImage, true);

	WriteTransform("transform.txt", bsplineTransform);

	return;
}
void BSplineRegistrationWithInitialTransform(char* fixedImageFileName, char* movingImageFileName, 
											 char* initialTransformFileName,
											 char* deformedImageFileName, char* transformFileName)
{
	InternalImageType::Pointer fixedImage = InternalImageType::New();
	ReadImage(fixedImageFileName, fixedImage);

	InternalImageType::Pointer movingImage = InternalImageType::New();
	ReadImage(movingImageFileName, movingImage);

	InternalImageType::Pointer matchedMovingImage = InternalImageType::New();
	HistogramMatching(movingImage, fixedImage, matchedMovingImage);

	BSplineTransformType::Pointer initialBSplineTransform = BSplineTransformType::New();
	ReadTransform(initialTransformFileName, initialBSplineTransform);

	InternalImageType::Pointer deformedImage = InternalImageType::New();

	BSplineTransformType::Pointer bsplineTransform = BSplineTransformType::New();

	BSplineRegistrationWithInitialTransform(fixedImage, matchedMovingImage, initialBSplineTransform, deformedImage, bsplineTransform);

	WriteImage(deformedImageFileName, deformedImage);

	ReadTransform("transform.txt", bsplineTransform);

	WriteTransform(transformFileName, bsplineTransform);

	return;
}
void BSplineRegistrationWithInitialTransform(InternalImageType::Pointer fixedImage, InternalImageType::Pointer movingImage,
											 BSplineTransformType::Pointer initialBSplineTransform, 
											 InternalImageType::Pointer &deformedImage, BSplineTransformType::Pointer &bsplineTransform)
{
	const unsigned int highRes = 64;
	const unsigned int gridBorderLength = 3;

	InternalMeanSquaresMetricType::Pointer		    metric        = InternalMeanSquaresMetricType::New();
	LBFGSOptimizerType::Pointer						optimizer     = LBFGSOptimizerType::New();
	InternalLinearInterpolatorType::Pointer			interpolator  = InternalLinearInterpolatorType::New();
	ImageRegistrationType::Pointer					registration  = ImageRegistrationType::New();

	registration->SetMetric(        metric        );
	registration->SetOptimizer(     optimizer     );
	registration->SetInterpolator(  interpolator  );
	registration->SetFixedImage(  fixedImage   );
	registration->SetMovingImage(   movingImage  );

	InternalImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
	registration->SetFixedImageRegion( fixedRegion );


	BSplineTransformType::Pointer  transformHigh = BSplineTransformType::New();
	registration->SetTransform( transformHigh );

	BSplineTransformType::RegionType::SizeType   gridBorderSize;
	BSplineTransformType::RegionType::SizeType   gridHighSizeOnImage;
	BSplineTransformType::RegionType::SizeType   totalGridSize;

	gridBorderSize.Fill( gridBorderLength );    // Border for spline order gridBorderLength= 3 ( 1 lower, 2 upper )
	gridHighSizeOnImage.Fill( highRes );
	totalGridSize = gridHighSizeOnImage + gridBorderSize;

	BSplineTransformType::RegionType bsplineRegion;
	bsplineRegion.SetSize( totalGridSize );

	BSplineTransformType::SpacingType spacingHigh = fixedImage->GetSpacing();
	BSplineTransformType::OriginType originHigh = fixedImage->GetOrigin();;

	InternalImageType::SizeType fixedImageSize = fixedRegion.GetSize();

	for(unsigned int r = 0; r < ImageDimension; r++)
	{ spacingHigh[r] *= static_cast<double>(fixedImageSize[r] - 1) / static_cast<double>(gridHighSizeOnImage[r] - 1);}

	InternalImageType::DirectionType gridDirection = fixedImage->GetDirection();
	BSplineTransformType::SpacingType gridOriginOffsetHigh = gridDirection * spacingHigh;

	BSplineTransformType::OriginType gridOriginHigh = originHigh - gridOriginOffsetHigh; 

	transformHigh->SetGridSpacing( spacingHigh );
	transformHigh->SetGridOrigin( gridOriginHigh );
	transformHigh->SetGridRegion( bsplineRegion );
	transformHigh->SetGridDirection( gridDirection );


	BSplineTransformType::ParametersType parametersHigh( transformHigh->GetNumberOfParameters() );
	parametersHigh.Fill( 0.0 );
	transformHigh->SetParameters( parametersHigh );

	// start to register in low resolution
	registration->SetInitialTransformParameters( initialBSplineTransform->GetParameters() );
	optimizer->SetGradientConvergenceTolerance( 0.005 );
	optimizer->SetLineSearchAccuracy( 0.9 );
	optimizer->SetDefaultStepLength( 0.5 );
	optimizer->TraceOff();
	optimizer->SetMaximumNumberOfFunctionEvaluations( 1000 );

	// std::cout << "Starting Registration with low resolution transform" << std::endl;

	try 
	{ 
		registration->StartRegistration(); 
	} 
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << "ExceptionObject caught at low resolution !" << std::endl; 
		std::cerr << err << std::endl; 
		return;
	} 

	// Finally we use the last transform parameters in order to resample the image.
	transformHigh->SetParameters( registration->GetLastTransformParameters() );

	bsplineTransform = transformHigh;

	ApplyBSplineTransform(movingImage, bsplineTransform, deformedImage, true);

	WriteTransform("transform.txt", bsplineTransform);

	return;
}
// generate zero deformation field
void GenerateZeroDeformationField(DeformationFieldType::Pointer &deformationFieldZero, int imageSize)
{
	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();

	DeformationFieldType::IndexType start;
	for (unsigned int i = 0; i < ImageDimension; i++) start[i] = 0;
	DeformationFieldType::SizeType size;
	for (unsigned int i = 0; i < ImageDimension; i++) size[i] = imageSize;
	DeformationFieldType::RegionType region;
	region.SetSize (size);
	region.SetIndex (start);

	deformationField -> SetRegions (region);
	deformationField -> Allocate();

	DeformationFieldIteratorType dfIt ( deformationField, deformationField->GetRequestedRegion() );
	VectorPixelType vectorPixel;
	for (dfIt.GoToBegin(); !dfIt.IsAtEnd(); ++dfIt)
	{
		vectorPixel.SetElement(0,0.0);
		vectorPixel.SetElement(1,0.0);
		dfIt.Set(vectorPixel);
	}

	deformationFieldZero = deformationField;
	return;
}
void GenerateZeroDeformationField3D(DeformationFieldType::Pointer &deformationFieldZero, int imageSize_x, int imageSize_y, int imageSize_z)
{
	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();

	DeformationFieldType::IndexType start;
	for (unsigned int i = 0; i < ImageDimension; i++) start[i] = 0;
	DeformationFieldType::SizeType size;
	// for (int i = 0; i < ImageDimension; i++) size[i] = imageSize;
	size[0] = imageSize_x;
	size[1] = imageSize_y;
	size[2] = imageSize_z;

	DeformationFieldType::RegionType region;
	region.SetSize (size);
	region.SetIndex (start);

	deformationField -> SetRegions (region);
	deformationField -> Allocate();

	DeformationFieldIteratorType dfIt ( deformationField, deformationField->GetRequestedRegion() );
	VectorPixelType vectorPixel;
	for (dfIt.GoToBegin(); !dfIt.IsAtEnd(); ++dfIt)
	{
		vectorPixel.SetElement(0,0.0);
		vectorPixel.SetElement(1,0.0);
		vectorPixel.SetElement(2,0.0);
		dfIt.Set(vectorPixel);
	}

	deformationFieldZero = deformationField;
	return;
}
// generate zero bspline transform
void GenerateZeroBSplineTransform(BSplineTransformType::Pointer &transformZero, InternalImageType::Pointer fixedImage, int gridSize)
{
	BSplineTransformType::Pointer transform = BSplineTransformType::New();

	BSplineTransformType::RegionType::SizeType   gridBorderSize;
	BSplineTransformType::RegionType::SizeType   gridHighSizeOnImage;
	BSplineTransformType::RegionType::SizeType   totalGridSize;

	gridBorderSize.Fill( gridBorderLength );    // Border for spline order gridBorderLength= 3 ( 1 lower, 2 upper )
	gridHighSizeOnImage.Fill( highRes );
	totalGridSize = gridHighSizeOnImage + gridBorderSize;

	BSplineTransformType::RegionType bsplineRegion;
	bsplineRegion.SetSize( totalGridSize );

	BSplineTransformType::SpacingType spacingHigh = fixedImage->GetSpacing();
	// for (int i = 0; i < ImageDimension; i++) spacingHigh[i] = 1.0;
	BSplineTransformType::OriginType originHigh = fixedImage->GetOrigin();;
	// for (int i = 0; i < ImageDimension; i++) originHigh[i] = 0.0;
	InternalImageType::SizeType fixedImageSize = fixedImage->GetBufferedRegion().GetSize();
	// for (int i = 0; i < ImageDimension; i++) fixedImageSize[i] = imageSize;

	for(unsigned int r = 0; r < ImageDimension; r++)
	{ spacingHigh[r] *= static_cast<double>(fixedImageSize[r] - 1) / static_cast<double>(gridHighSizeOnImage[r] - 1);}

	InternalImageType::DirectionType gridDirection = fixedImage->GetDirection();
	BSplineTransformType::SpacingType gridOriginOffsetHigh = gridDirection * spacingHigh;

	BSplineTransformType::OriginType gridOriginHigh = originHigh - gridOriginOffsetHigh; 

	transform->SetGridSpacing( spacingHigh );
	transform->SetGridOrigin( gridOriginHigh );
	transform->SetGridRegion( bsplineRegion );
	transform->SetGridDirection( gridDirection );


	BSplineTransformType::ParametersType parametersHigh( transform->GetNumberOfParameters() );
	parametersHigh.Fill( 0.0 );
	transform->SetParameters( parametersHigh );
	transform->SetIdentity();

	transformZero = transform;
	return;
}

// other filter
// histogram matching
void HistogramMatching(InternalImageType::Pointer inputImage, 
					   InternalImageType::Pointer referenceImage, 
					   InternalImageType::Pointer &outputImage)
{
	InternalHistMatchFilterType::Pointer matcher = InternalHistMatchFilterType::New();
	matcher->SetInput( inputImage );
	matcher->SetReferenceImage( referenceImage );
	matcher->SetNumberOfHistogramLevels( 1024 ); // 1024
	matcher->SetNumberOfMatchPoints( 7 );
	matcher->ThresholdAtMeanIntensityOn();
	try
	{
		matcher->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << err << std::endl;
		return;
	}
	outputImage = matcher->GetOutput();
	return;
}
void HistogramMatchingWithNames(char* inputImageName, char* referenceImageName, char* outputImageName)
{
	InternalImageType::Pointer inputImage = InternalImageType::New();
	ReadImage(inputImageName, inputImage);

	InternalImageType::Pointer referenceImage = InternalImageType::New();
	ReadImage(referenceImageName, referenceImage);

	InternalImageType::Pointer outputImage = InternalImageType::New();

	HistogramMatching(inputImage, referenceImage, outputImage);

	WriteImage(outputImageName, outputImage);
	return;
}
// calculate mean image
void calculateMeanImage(char** imagefiles, int filenumber, InternalImageType::Pointer &meanImage)
{
	if (filenumber == 1)
	{
		ReadImage(imagefiles[0], meanImage);
	}
	else
	{
		int numberOfDeformationLoaded = 0;

		InternalImageReaderType::Pointer imageReader1 = InternalImageReaderType::New();
		InternalImageReaderType::Pointer imageReader2 = InternalImageReaderType::New();
		imageReader1->SetFileName(imagefiles[0]);
		imageReader2->SetFileName(imagefiles[1]);

		AddInternalImageFilterType::Pointer addImageSum = AddInternalImageFilterType::New();
		addImageSum->SetInput1(imageReader1->GetOutput());
		addImageSum->SetInput2(imageReader2->GetOutput());
		numberOfDeformationLoaded += 2;
		try
		{
			addImageSum->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught at adding deformation field" << std::endl;
			std::cerr << excep << std::endl;
		}

		for (int j = 2; j < filenumber; j++) 
		{
			InternalImageReaderType::Pointer imageReader = InternalImageReaderType::New();
			imageReader->SetFileName(imagefiles[j]);

			addImageSum->SetInput1(addImageSum->GetOutput());
			addImageSum->SetInput2(imageReader->GetOutput());
			numberOfDeformationLoaded += 1;
			try
			{
				addImageSum->Update();
			}
			catch( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception caught at adding deformation field" << std::endl;
				std::cerr << excep << std::endl;
			}
		}
		//////////////////////////////
		InternalImageType::Pointer sumImage = InternalImageType::New();
		sumImage = addImageSum->GetOutput();
		sumImage -> DisconnectPipeline();
		//InternalImageIteratorType meanIter1 ( sumImage, sumImage->GetRequestedRegion() );
		//int index1 = 0;
		//float meanPixel1;
		//for (meanIter1.GoToBegin(); !meanIter1.IsAtEnd(); ++meanIter1)
		//{
		//	if (index1 == (78-1)*256+(133-1))
		//	{
		//		index1 = index1;
		//		meanPixel1 = meanIter1.Get();
		//	}
		//	// for (int i = 0; i < filenumber; i++) pixel[i] = Iter[i].Get();
		//	index1++;
		//}

		/////////////////////////////
		DivideInternalImageFilterType::Pointer divideImage = DivideInternalImageFilterType::New();
		const float factor = (float)filenumber;
		divideImage->SetConstant(factor);
		// divideImage->SetConstant(100.0);
		// float cs = divideImage->GetConstant();
		divideImage->SetInput(sumImage);
		try
		{
			divideImage->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught at dividing deformation field" << std::endl;
			std::cerr << excep << std::endl;
		}
		//////////////////////////////
		//InternalImageIteratorType meanIter2 ( divideImage->GetOutput(), divideImage->GetOutput()->GetRequestedRegion() );
		//int index2 = 0;
		//float meanPixel2;
		//for (meanIter2.GoToBegin(); !meanIter2.IsAtEnd(); ++meanIter2)
		//{
		//	if (index2 == (78-1)*256+(133-1))
		//	{
		//		index2 = index2;
		//		meanPixel2 = meanIter2.Get();
		//	}
		//	// for (int i = 0; i < filenumber; i++) pixel[i] = Iter[i].Get();
		//	index2++;
		//}

		/////////////////////////////

		meanImage = divideImage->GetOutput();
		meanImage -> DisconnectPipeline();

	}
	return;
}


// other auxiliary sub-functions
// convert an integer into string with specified digits
void myitoa(int num, char* str, int digit)
{
	// translate the last (1~3) digits of a number to string
	if (digit == 1)
	{
		str[0] = num%10 + 48;
		str[1] = 0;
	}
	else if (digit == 2)
	{
		str[0] = (num%100)/10 + 48;
		str[1] = num%10 + 48;
		str[2] = 0;
	}
	else if (digit == 3)
	{
		str[0] = (num%1000)/100 + 48;
		str[1] = (num%100)/10 + 48;
		str[2] = num%10 + 48;
		str[3] = 0;
	}
	else
	{
		str[0] = 0;
	}
	return;
}


// bubble sort on double array
void bubbleSort(double* arr, int* index, int n) 
{
	bool swapped = true;
	int j = 0;
	double tmp;
	int index_tmp;
	while (swapped) 
	{
		swapped = false;
		j++;
		for (int i = 0; i < n - j; i++) 
		{
			if (arr[i] > arr[i + 1]) 
			{
				tmp = arr[i];
				arr[i] = arr[i + 1];
				arr[i + 1] = tmp;

				index_tmp = index[i];
				index[i] = index[i+1];
				index[i+1] = index_tmp;

				swapped = true;
			}
		}
	}
}
// calculate pair-wise distance in an image population
void PairwiseDistanceAmongImages(char** imageFileNames, int totalNumber, double** distanceMatrix)
{
	for (int i = 0; i< totalNumber; i++)
	{
		//std::cerr << i << ", ";
		for (int j = i; j < totalNumber; j++)
		{
			if (j == i)
				distanceMatrix[i][j] = 0.0;
			else
			{
				double dist = 0.0;
				// dist = calculateDistance(imageFileNames[i], imageFileNames[j]);
				dist = calculateDistanceMSD(imageFileNames[i], imageFileNames[j]);
				distanceMatrix[i][j] = sqrt(dist);
				distanceMatrix[j][i] = sqrt(dist);
			}
		}
	}
	//std::cerr << "done!" << endl;
	return;
}
// calculate distance between two images
double calculateDistance(char* imageName1, char* imageName2)
{
	double dist;

	InternalImageReaderType::Pointer  imageReader1  = InternalImageReaderType::New();
	InternalImageReaderType::Pointer  imageReader2 = InternalImageReaderType::New();
	imageReader1->SetFileName( imageName1 );
	imageReader2->SetFileName( imageName2 );
	try 
	{
		imageReader1->Update();
		imageReader2->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	InternalImageType::ConstPointer image1 = imageReader1->GetOutput();
	InternalImageType::ConstPointer image2 = imageReader2->GetOutput();

	TranslationTransformType::Pointer transform = TranslationTransformType::New();
	InternalNNInterpolatorType::Pointer interpolator = InternalNNInterpolatorType::New();
	InternalMeanSquaresMetricType::Pointer metric = InternalMeanSquaresMetricType::New();
	transform->SetIdentity();
	metric->SetTransform( transform );
	metric->SetInterpolator( interpolator );
	metric->SetFixedImage(image1);
	metric->SetMovingImage(image2);
	metric->SetFixedImageRegion(  image1->GetBufferedRegion()  );

	try 
	{
		metric->Initialize();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	InternalMeanSquaresMetricType::TransformParametersType displacement( ImageDimension );
	for (unsigned int i = 0; i < ImageDimension; i++)
		displacement[i] = 0.0; 
	dist = metric->GetValue( displacement );

	return dist;
}

// calculate distance between two images by mean squared difference
double calculateDistanceMSD(char* imageName1, char* imageName2)
{
	double dist;

	InternalImageReaderType::Pointer  imageReader1  = InternalImageReaderType::New();
	InternalImageReaderType::Pointer  imageReader2 = InternalImageReaderType::New();
	imageReader1->SetFileName( imageName1 );
	imageReader2->SetFileName( imageName2 );
	try 
	{
		imageReader1->Update();
		imageReader2->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	InternalImageType::Pointer image1 = imageReader1->GetOutput();
	InternalImageType::Pointer image2 = imageReader2->GetOutput();
	
	InternalImageType::SizeType size;
	size = image1->GetLargestPossibleRegion().GetSize();
	InternalImageIteratorType it1(image1, image1->GetLargestPossibleRegion() );
	InternalImageIteratorType it2(image2, image2->GetLargestPossibleRegion() );
	dist = 0.0;
	float value1, value2;
	for (it1.GoToBegin(),it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2)
	{
		value1 = it1.Get();
		value2 = it2.Get();
		dist += (value1-value2)*(value1-value2);
	}
	dist = dist/(size[0]*size[1]*size[2]);
	return dist;
}


// write 2D double matrix file
void SaveMatrix2File(double** matrix, int iSize, int jSize, char* martixFileName)
{
	std::ofstream outfile;
	outfile.open (martixFileName);

	for(int i = 0; i < iSize; i++)
	{
		for(int j = 0; j < jSize; j++)
		{
				outfile <<  ' ' ;
				outfile << std::right << std::fixed << std::setw(8) << std::setprecision(2) << matrix[i][j];
		}
		outfile << std::endl;
	}
	outfile.close();
	return;
}
// read 2D double matrix file
void ReadMatrixFile(double** matrix, int iSize, int jSize, char* martixFileName)
{
	std::ifstream infile;
	double temp;
	infile.open (martixFileName);

	for(int i = 0; i < iSize; i++)
	{
		for(int j = 0; j < jSize; j++)
		{
			
			infile >> temp ;
			matrix[i][j] = temp;
		}
	}
	infile.close();
	return;
}
// calculate kNN Isomap pair-wise distance matrix
void CalculateIsomapKNN(double** distanceMatrix, int matrixSize, double** distanceIsomapMatrix, int K)
{
	double** distance = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)		distance[i] = new double[matrixSize];

	double** distance_temp = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)		distance_temp[i] = new double[matrixSize];

	int* index = new int[matrixSize];

	double* distance1 = new double[matrixSize];

	double INF = 0;
	double distMax = 0;

	// make a copy of distanceMatrix into distance and find maximum
	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			distance[i][j] = distanceMatrix[i][j];
			if (distance[i][j]>distMax)
				distMax = distance[i][j];
		}
	}
	// setup an infinite number
	INF = distMax*matrixSize*1000;
	// choose those pairs k-NN
	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			distance1[j] = distance[i][j];
			index[j] = j;
		}
		bubbleSort(distance1, index, matrixSize);
		for (int j = K+1; j < matrixSize; j++)
		{
			distance[i][index[j]] = INF;
		}
	}
	// make distance matrix be symmetric
	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = i+1; j < matrixSize; j++)
		{
			if (distance[i][j] > distance[j][i])
				distance[i][j] = distance[j][i];
			else
				distance[j][i] = distance[i][j];
		}
	}
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
	for (int i = 0; i < matrixSize; i++)
		for (int j = 0; j < matrixSize; j++)
			distanceIsomapMatrix[i][j] = distance[i][j];

	delete[] index;
	delete[] distance1;
	for (int i = 0; i < matrixSize; i++)		delete[] distance[i];
	delete[] distance;
	for (int i = 0; i < matrixSize; i++)		delete[] distance_temp[i];
	delete[] distance_temp;

	return;
}

// find two nodes in a graph whose connecting path has the largest distance
void FindMaxLocationInMatrix(double** distanceMatrix, int matrixSize, int &max, int &neighbor)
{
	int max1, max2;
	double dist1, dist2;
	int neighbor1, neighbor2;
	double** distance = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)		distance[i] = new double[matrixSize];

	double* distance_max = new double[matrixSize];

	int* index = new int[matrixSize];

	double* distance1 = new double[matrixSize];

	// double distMax = 0;

	// make a copy of distanceMatrix into distance and find maximum
	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			if (i > j)
				distance[i][j] = distanceMatrix[i][j];
			else
				distance[i][j] = 0.0;
		}
	}

	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			distance1[j] = distance[i][j];
			index[j] = j;
		}
		bubbleSort(distance1, index, matrixSize);
		distance_max[i] = distance1[matrixSize-1];
	}
	for (int i = 0; i < matrixSize; i++) index[i] = i;
	bubbleSort(distance_max, index, matrixSize);
	max1 = index[matrixSize-1];

	for (int j = 0; j < matrixSize; j++)
	{
		for (int i = 0; i < matrixSize; i++)
		{
			distance1[i] = distance[i][j];
			index[i] = i;
		}
		bubbleSort(distance1, index, matrixSize);
		distance_max[j] = distance1[matrixSize-1];
	}
	for (int i = 0; i < matrixSize; i++) index[i] = i;
	bubbleSort(distance_max, index, matrixSize);
	max2 = index[matrixSize-1];

	// choose between max1 and max2 based on whose nearest neighbor has larger distance
	for (int i = 0; i < matrixSize; i++)
	{
		distance1[i] = distanceMatrix[i][max1];
		index[i] = i;
	}
	bubbleSort(distance1, index, matrixSize);
	dist1 = distance1[1];
	neighbor1 = index[1];

	for (int i = 0; i < matrixSize; i++)
	{
		distance1[i] = distanceMatrix[i][max2];
		index[i] = i;
	}
	bubbleSort(distance1, index, matrixSize);
	dist2 = distance1[1];
	neighbor2 = index[1];

	max = (dist1>dist2)?max1:max2;
	neighbor = (dist1>dist2)?neighbor1:neighbor2;

	delete[] index;
	delete[] distance1;
	for (int i = 0; i < matrixSize; i++)		delete[] distance[i];
	delete[] distance;
	delete[] distance_max;
	
	return;
}

// find two nodes in a graph whose connecting path has the largest distance considering the subgroup size
void FindMaxLocationInMatrix(double** distanceMatrix, int matrixSize, int* subGroupSize, int &max, int &neighbor)
{
	int max1, max2;
	double dist1, dist2;
	int neighbor1, neighbor2;
	int subsize1, subsize2;

	double** distance = new double*[matrixSize];
	for (int i = 0; i < matrixSize; i++)		distance[i] = new double[matrixSize];

	double* distance_max = new double[matrixSize];

	int* index = new int[matrixSize];

	double* distance1 = new double[matrixSize];

	// double distMax = 0;

	// make a copy of distanceMatrix into distance and find maximum
	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			if (i > j)
				distance[i][j] = distanceMatrix[i][j];
			else
				distance[i][j] = 0.0;
		}
	}

	for (int i = 0; i < matrixSize; i++)
	{
		for (int j = 0; j < matrixSize; j++)
		{
			distance1[j] = distance[i][j];
			index[j] = j;
		}
		bubbleSort(distance1, index, matrixSize);
		distance_max[i] = distance1[matrixSize-1];
	}
	for (int i = 0; i < matrixSize; i++) index[i] = i;
	bubbleSort(distance_max, index, matrixSize);
	max1 = index[matrixSize-1];

	for (int j = 0; j < matrixSize; j++)
	{
		for (int i = 0; i < matrixSize; i++)
		{
			distance1[i] = distance[i][j];
			index[i] = i;
		}
		bubbleSort(distance1, index, matrixSize);
		distance_max[j] = distance1[matrixSize-1];
	}
	for (int i = 0; i < matrixSize; i++) index[i] = i;
	bubbleSort(distance_max, index, matrixSize);
	max2 = index[matrixSize-1];

	// choose between max1 and max2 based on whose nearest neighbor has larger distance
	for (int i = 0; i < matrixSize; i++)
	{
		distance1[i] = distanceMatrix[i][max1];
		index[i] = i;
	}
	bubbleSort(distance1, index, matrixSize);
	dist1 = distance1[1];
	neighbor1 = index[1];

	for (int i = 0; i < matrixSize; i++)
	{
		distance1[i] = distanceMatrix[i][max2];
		index[i] = i;
	}
	bubbleSort(distance1, index, matrixSize);
	dist2 = distance1[1];
	neighbor2 = index[1];
	subsize1 = subGroupSize[max1];
	subsize2 = subGroupSize[max2];
	if (subsize1 == subsize2)
	{
		max = (dist1>dist2)?max1:max2;
		neighbor = (dist1>dist2)?neighbor1:neighbor2;
	}
	else
	{
		max = (subsize1<subsize2)?max1:max2;
		neighbor = (subsize1<subsize2)?neighbor1:neighbor2;
	}

	delete[] index;
	delete[] distance1;
	for (int i = 0; i < matrixSize; i++)		delete[] distance[i];
	delete[] distance;
	delete[] distance_max;

	return;
}
// find a tree in a graph by merging one node on the longest path to its nearest neighbor
void FindGraphMaxTree(double** distanceMatrixTotal, int matrixSize, int* tree)
{
	/////////////////////////////////
	// initialize some variables

	int groupLevel = 1;
	bool isEnd = false;
	int idTemplate = -1;

	int groupSize = matrixSize;
	int** groupInfo = new int*[matrixSize];
	int*  subGroupSize = new int[matrixSize];
	int*  subGroupRep = new int[matrixSize];
	for (int i = 0; i < matrixSize; i++)
	{
		subGroupSize[i] = 1;
		groupInfo[i] = new int[matrixSize];
		for (int j = 0; j < matrixSize; j++) groupInfo[i][j] = i;
		subGroupRep[i] = groupInfo[i][0];
	}

	// int groupSizeNew = matrixSize;
	int** groupInfoNew = new int*[matrixSize];
	int*  subGroupSizeNew = new int[matrixSize];
	int*  subGroupRepNew = new int[matrixSize];
	for (int i = 0; i < matrixSize; i++)
	{
		subGroupSizeNew[i] = 1;
		groupInfoNew[i] = new int[matrixSize];
		for (int j = 0; j < matrixSize; j++) groupInfoNew[i][j] = i;
		subGroupRepNew[i] = groupInfoNew[i][0];
	}

	// start iteration
	while (!isEnd)
	{
		//////////////////////////////////////
		// extract current pair-wise distance matrix from distanceMatrixTotal
		double** distanceMatrix = new double*[groupSize];
		for (int i = 0; i < groupSize; i++) distanceMatrix[i] = new double[groupSize];
		for (int i = 0; i < groupSize; i++)
		{
			for (int j = 0; j < groupSize; j++)
			{
				int ii, jj;
				ii = subGroupRep[i];
				jj = subGroupRep[j];
				distanceMatrix[i][j] = distanceMatrixTotal[ii][jj];
				//distanceMatrix[i][j] = distanceMatrixTotal[subGroupRep[i]][subGroupRep[j]];
			}
		}

		// extract current pair-wise distance matrix from distanceMatrixTotal
		////////////////////////////////////
		// find the longest distance and locate the current moving subject and its neighbor
		// move idCur to idCurNeighbor

		int idCur = -1;
		int idCurNeighbor = -1; 

		if (groupSize == 1)  // add 20091221
		{
			tree[0] = -1;
		}
		else if (groupSize == 2)
		{
			// when only two subgroups, less move to more
			if (subGroupSize[0] > subGroupSize[1])
			{
				idCur = 1;
				idCurNeighbor = 0;
				idTemplate = subGroupRep[0];
			}
			else
			{
				idCur = 0;
				idCurNeighbor = 1;
				idTemplate = subGroupRep[1];
			}
			tree[subGroupRep[idCur]] = subGroupRep[idCurNeighbor];
			tree[subGroupRep[idCurNeighbor]] = -1;
		}
		else
		{
			FindMaxLocationInMatrix(distanceMatrix, groupSize, subGroupSize, idCur, idCurNeighbor);
			tree[subGroupRep[idCur]] = subGroupRep[idCurNeighbor];
		}
		// find the longest distance and locate the current moving subject and its neighbor
		//////////////////////////////
		// generate new group info
		int groupSizeNew = groupSize - 1;
		int index = 0;
		for (int i = 0; i < groupSize; i++)
		{
			if ((i != idCur) && (i != idCurNeighbor))
			{
				subGroupSizeNew[index] = subGroupSize[i];
				subGroupRepNew[index] = subGroupRep[i];
				groupInfoNew[index] = new int[subGroupSizeNew[index]];
				for (int j = 0; j < subGroupSizeNew[index]; j++)
				{
					groupInfoNew[index][j] = groupInfo[i][j];
				}
				index++;
			}
			else if (i == idCurNeighbor)
			{
				subGroupSizeNew[index] = subGroupSize[idCurNeighbor]+subGroupSize[idCur];
				subGroupRepNew[index] = subGroupRep[idCurNeighbor];
				groupInfoNew[index] = new int[subGroupSizeNew[index]];
				for (int j = 0; j < subGroupSize[idCurNeighbor]; j++)
				{
					groupInfoNew[index][j] = groupInfo[idCurNeighbor][j];
				}
				for (int j = 0; j < subGroupSize[idCur]; j++)
				{
					groupInfoNew[index][j+subGroupSize[idCurNeighbor]] = groupInfo[idCur][j];
				}
				index++;
			}
		}
		// generate new group info
		////////////////////////////////////
		// goto next level

		groupLevel ++;
		if (groupSize <= 2) // update on 20091221 from groupSize == 2
			isEnd = true;

		// goto next level
		//////////////////////////////////
		// delete newed variables for current level

		for (int i = 0; i < groupSize; i++) delete[] distanceMatrix[i];
		delete[] distanceMatrix;

		// generate group information for next level
		for (int i = 0; i < groupSizeNew; i++)
		{
			subGroupRep[i] = subGroupRepNew[i];
			subGroupSize[i] = subGroupSizeNew[i];
			groupInfo[i] = new int[subGroupSize[i]];
			for (int j = 0; j < subGroupSize[i]; j++)
				groupInfo[i][j] = groupInfoNew[i][j];
		}
		groupSize = groupSizeNew;


	} // end of while (!isEnd)


	// start iteration
	///////////////////////////////////
	//// output tree
	//std::ofstream outfile;
	//outfile.open ("tree_output.txt");
	//outfile << matrixSize << std::endl;
	//for(int i = 0; i < matrixSize; i++)
	//{
	//	outfile <<  ' ' ;
	//	outfile << right << setw(3) << i; 
	//	outfile << ' ';
	//	outfile << right << setw(3) << tree[i];
	//	outfile << std::endl;
	//}
	//outfile.close();
	// output tree
	////////////////////////////////////
	// delete newed variables
	for (int i = 0; i < matrixSize; i++) delete[] groupInfo[i];
	delete[] groupInfo;
	delete[] subGroupSize;
	delete[] subGroupRep;

	for (int i = 0; i < matrixSize; i++) delete[] groupInfoNew[i];
	delete[] groupInfoNew;
	delete[] subGroupSizeNew;
	delete[] subGroupRepNew;
	
	return;
}
// find a tree in a graph by iteratively updating the kNN-isomap
void FindKNNGraphMaxTree(double** distanceMatrixTotal, int matrixSize, int kNN, int* tree)
{
	/////////////////////////////////
	// initialize some variables

	int groupLevel = 1;
	bool isEnd = false;
	int idTemplate = -1;

	int groupSize = matrixSize;
	int** groupInfo = new int*[matrixSize];
	int*  subGroupSize = new int[matrixSize];
	int*  subGroupRep = new int[matrixSize];
	for (int i = 0; i < matrixSize; i++)
	{
		subGroupSize[i] = 1;
		groupInfo[i] = new int[matrixSize];
		for (int j = 0; j < matrixSize; j++) groupInfo[i][j] = i;
		subGroupRep[i] = groupInfo[i][0];
	}

	// int groupSizeNew = matrixSize;
	int** groupInfoNew = new int*[matrixSize];
	int*  subGroupSizeNew = new int[matrixSize];
	int*  subGroupRepNew = new int[matrixSize];
	for (int i = 0; i < matrixSize; i++)
	{
		subGroupSizeNew[i] = 1;
		groupInfoNew[i] = new int[matrixSize];
		for (int j = 0; j < matrixSize; j++) groupInfoNew[i][j] = i;
		subGroupRepNew[i] = groupInfoNew[i][0];
	}

	// start iteration
	while (!isEnd)
	{
		//////////////////////////////////////
		// extract current pair-wise distance matrix from distanceMatrixTotal
		double** distanceMatrix = new double*[groupSize];
		for (int i = 0; i < groupSize; i++) distanceMatrix[i] = new double[groupSize];
		for (int i = 0; i < groupSize; i++)
		{
			for (int j = 0; j < groupSize; j++)
			{
				int ii, jj;
				ii = subGroupRep[i];
				jj = subGroupRep[j];
				distanceMatrix[i][j] = distanceMatrixTotal[ii][jj];
				//distanceMatrix[i][j] = distanceMatrixTotal[subGroupRep[i]][subGroupRep[j]];
			}
		}

		if (kNN >= groupSize) kNN = groupSize-1;
		double ** distanceIsomapMatrix = new double*[groupSize];
		for (int i = 0;i < groupSize; i++) distanceIsomapMatrix[i] = new double[groupSize];
		CalculateIsomapKNN(distanceMatrix, groupSize, distanceIsomapMatrix, kNN);

		// extract current pair-wise distance matrix from distanceMatrixTotal
		////////////////////////////////////
		// find the longest distance and locate the current moving subject and its neighbor
		// move idCur to idCurNeighbor

		int idCur = -1;
		int idCurNeighbor = -1; 

		if (groupSize == 2)
		{
			// when only two subgroups, less move to more
			if (subGroupSize[0] > subGroupSize[1])
			{
				idCur = 1;
				idCurNeighbor = 0;
				idTemplate = subGroupRep[0];
			}
			else
			{
				idCur = 0;
				idCurNeighbor = 1;
				idTemplate = subGroupRep[1];
			}
			tree[subGroupRep[idCur]] = subGroupRep[idCurNeighbor];
			tree[subGroupRep[idCurNeighbor]] = -1;
		}
		else
		{
			FindMaxLocationInMatrix(distanceIsomapMatrix, groupSize, subGroupSize, idCur, idCurNeighbor);
			tree[subGroupRep[idCur]] = subGroupRep[idCurNeighbor];
		}
		// find the longest distance and locate the current moving subject and its neighbor
		//////////////////////////////
		// generate new group info
		int groupSizeNew = groupSize - 1;
		int index = 0;
		for (int i = 0; i < groupSize; i++)
		{
			if ((i != idCur) && (i != idCurNeighbor))
			{
				subGroupSizeNew[index] = subGroupSize[i];
				subGroupRepNew[index] = subGroupRep[i];
				groupInfoNew[index] = new int[subGroupSizeNew[index]];
				for (int j = 0; j < subGroupSizeNew[index]; j++)
				{
					groupInfoNew[index][j] = groupInfo[i][j];
				}
				index++;
			}
			else if (i == idCurNeighbor)
			{
				subGroupSizeNew[index] = subGroupSize[idCurNeighbor]+subGroupSize[idCur];
				subGroupRepNew[index] = subGroupRep[idCurNeighbor];
				groupInfoNew[index] = new int[subGroupSizeNew[index]];
				for (int j = 0; j < subGroupSize[idCurNeighbor]; j++)
				{
					groupInfoNew[index][j] = groupInfo[idCurNeighbor][j];
				}
				for (int j = 0; j < subGroupSize[idCur]; j++)
				{
					groupInfoNew[index][j+subGroupSize[idCurNeighbor]] = groupInfo[idCur][j];
				}
				index++;
			}
		}
		// generate new group info
		////////////////////////////////////
		// goto next level

		groupLevel ++;
		if (groupSize == 2)
			isEnd = true;

		// goto next level
		//////////////////////////////////
		// delete newed variables for current level

		for (int i = 0; i < groupSize; i++) delete[] distanceMatrix[i];
		delete[] distanceMatrix;
		for (int i = 0; i < groupSize; i++) delete[] distanceIsomapMatrix[i];
		delete[] distanceIsomapMatrix;

		// generate group information for next level
		for (int i = 0; i < groupSizeNew; i++)
		{
			subGroupRep[i] = subGroupRepNew[i];
			subGroupSize[i] = subGroupSizeNew[i];
			groupInfo[i] = new int[subGroupSize[i]];
			for (int j = 0; j < subGroupSize[i]; j++)
				groupInfo[i][j] = groupInfoNew[i][j];
		}
		groupSize = groupSizeNew;


	} // end of while (!isEnd)


	// start iteration
	///////////////////////////////////
	//// output tree
	//std::ofstream outfile;
	//outfile.open ("tree_output.txt");
	//outfile << matrixSize << std::endl;
	//for(int i = 0; i < matrixSize; i++)
	//{
	//	outfile <<  ' ' ;
	//	outfile << right << setw(3) << i; 
	//	outfile << ' ';
	//	outfile << right << setw(3) << tree[i];
	//	outfile << std::endl;
	//}
	//outfile.close();
	// output tree
	////////////////////////////////////
	// delete newed variables
	for (int i = 0; i < matrixSize; i++) delete[] groupInfo[i];
	delete[] groupInfo;
	delete[] subGroupSize;
	delete[] subGroupRep;

	for (int i = 0; i < matrixSize; i++) delete[] groupInfoNew[i];
	delete[] groupInfoNew;
	delete[] subGroupSizeNew;
	delete[] subGroupRepNew;

}
// find all leaves in a tree
void FindLeavesInTree(int* tree, int treeSize, bool* isLeave)
{
	for(int i = 0; i < treeSize; i++) isLeave[i] = false;
	for(int i = 0; i < treeSize; i++)
	{
		if (tree[i] >= 0)
			isLeave[tree[i]] = true;
	}
	return;
}
// calculate the height of each node in a tree
void CalculateNodeSegmentsToRoot(int* tree, int treeSize, int* height)
{
	for(int i = 0; i < treeSize; i++) height[i] = 0;
	for(int i = 0; i < treeSize; i++)
	{
		bool isAtRoot = false;
		int curNode = i;
		while (!isAtRoot)
		{
			if (tree[curNode] >= 0)
			{
				height[i]++;
				curNode = tree[curNode];
			}
			else
			{
				isAtRoot = true;
			}
		}
	}
	return;
}
// Find the root node in a tree
void FindRoot(int* tree, int treeSize, int &root)
{
	for(int i = 0; i < treeSize; i++)
	{
		if (tree[i] == -1)
		{
			root = i;
			return;
		}
	}
	root = -1;
	return;
}
// save 1D integer array with index
void SaveArrayWithIndex(int* arr, int arrSize, char* filename)
{
	std::ofstream outfile;
	outfile.open (filename);
	for(int i = 0; i < arrSize; i++)
	{
		outfile <<  ' ' ;
		outfile << std::right << std::setw(3) << i; 
		outfile << ' ';
		outfile << std::right << std::setw(3) << arr[i];
		outfile << std::endl;
	}
	outfile.close();

	return;
}
// save 1D integer array with index and size
void SaveArrayWithIndexSize(int* arr, int arrSize, char* filename)
{
	std::ofstream outfile;
	outfile.open (filename);
	outfile << arrSize << std::endl;
	for(int i = 0; i < arrSize; i++)
	{
		outfile <<  ' ' ;
		outfile << std::right << std::setw(3) << i; 
		outfile << ' ';
		outfile << std::right << std::setw(3) << arr[i];
		outfile << std::endl;
	}
	outfile.close();

	return;
}
// save 1D integer array with index
void SaveArrayWithIndex(double* arr, int arrSize, char* filename)
{
	std::ofstream outfile;
	outfile.open (filename);
	for(int i = 0; i < arrSize; i++)
	{
		outfile <<  ' ' ;
		outfile << std::right << std::setw(3) << i; 
		outfile << ' ';
		outfile << std::right << std::fixed << std::setw(8) << std::setprecision(2) << arr[i];
		outfile << std::endl;
	}
	outfile.close();

	return;
}
// save 1D integer array with index and size
void SaveArrayWithIndexSize(double* arr, int arrSize, char* filename)
{
	std::ofstream outfile;
	outfile.open (filename);
	outfile << arrSize << std::endl;
	for(int i = 0; i < arrSize; i++)
	{
		outfile <<  ' ' ;
		outfile << std::right << std::setw(3) << i; 
		outfile << ' ';
		outfile << std::right << std::fixed << std::setw(8) << std::setprecision(2) << arr[i];
		outfile << std::endl;
	}
	outfile.close();

	return;
}
// calculate the total path length from each leaf node to the root
void CalculateTotalPathLength(int* tree, int treeSize, double* edgeLength, double &totalLength)
{
	int* multiplierOnEdge = new int[treeSize];
	for (int i = 0; i < treeSize; i++)
		multiplierOnEdge[i] = 0;
	for(int i = 0; i < treeSize; i++)
	{
		bool isAtRoot = false;
		int curNode = i;
		multiplierOnEdge[i]++; // count for itself
		while (!isAtRoot)
		{
			if (tree[curNode] >= 0)
			{
				multiplierOnEdge[tree[curNode]]++; // count for each of its predecessors
				curNode = tree[curNode];
			}
			else
			{
				multiplierOnEdge[tree[curNode]]++;
				isAtRoot = true;
			}
		}
	}
	// SaveArrayWithIndexSize(multiplierOnEdge, treeSize, "multiplier.txt");
	totalLength = 0;
	for (int i = 0; i < treeSize; i++) totalLength += edgeLength[i]*multiplierOnEdge[i];
	return;
}
// calculate the total path length based on edge length and subtree size
void CalculateTotalPathLength(double* edgeLength, int* subTreeSize, int treeSize, double & totalLegnth)
{
	totalLegnth = 0.0;
	for (int i = 0; i < treeSize; i++)
		totalLegnth += edgeLength[i]*subTreeSize[i];
	return;
}
// calculate the sub-tree size for each node
void CalculateSubTreeSize(int* tree, int treeSize, int* subTreeSize)
{
	for (int i = 0; i < treeSize; i++)
		subTreeSize[i] = 0;
	for(int i = 0; i < treeSize; i++)
	{
		bool isAtRoot = false;
		int curNode = i;
		subTreeSize[curNode]++; // count for itself
		while (!isAtRoot)
		{
			if (tree[curNode] >= 0)
			{
				subTreeSize[tree[curNode]]++; // count for each of its predecessors
				curNode = tree[curNode];
			}
			else
			{
				// subTreeSize[tree[curNode]]++;
				isAtRoot = true;
			}
		}
	}
	return;
}
// get the height of the tree
void GetTreeHeight(int* tree, int treeSize, int &treeHeight)
{
	int* height = new int[treeSize];
	CalculateNodeSegmentsToRoot(tree, treeSize, height);
	treeHeight = 0;
	for (int i = 0; i < treeSize; i++)
		if (height[i] > treeHeight)
			treeHeight = height[i];
	return;
}
// save tree with root
void SaveTreeWithInfo(int* tree, int treeSize, char* filename)
{
	std::ofstream outfile;
	outfile.open (filename);
	outfile << "Tree Size : " << treeSize << std::endl;
	int root;
	FindRoot(tree, treeSize, root);
	outfile << "Tree Root : " << root << std::endl;
	int height;
	GetTreeHeight(tree, treeSize, height);
	outfile << "Tree Height : " << height << std::endl;

	for(int i = 0; i < treeSize; i++)
	{
		outfile <<  ' ' ;
		outfile << std::right << std::setw(3) << i; 
		outfile << ' ';
		outfile << std::right << std::fixed << std::setw(8) << std::setprecision(2) << tree[i];
		outfile << std::endl;
	}
	outfile.close();

	return;

}
// read tree with info
void ReadTreeWithInfo(char* filename, int treeSize, int* tree)
{
	int _treeSize;
	int _treeRoot;
	int _treeHeight;

	std::ifstream inTreeFile;
	inTreeFile.open(filename);
	//Tree Size :
	//Tree Root :
	//Tree Height :
	//char temp;
	inTreeFile.seekg (strlen("Tree Size : "), std::ios::beg);
	inTreeFile >> _treeSize;
	//inTreeFile >> temp; //T
	//inTreeFile >> temp; //r
	//inTreeFile >> temp; //e
	//inTreeFile >> temp; //e
	//inTreeFile >> temp; //R
	//inTreeFile >> temp; //o
	//inTreeFile >> temp; //o
	//inTreeFile >> temp; //t
	//inTreeFile >> temp; //:


	inTreeFile.seekg (strlen(" Tree Root : "), std::ios::cur);
	inTreeFile >> _treeRoot;
	//inTreeFile >> temp; //T
	//inTreeFile >> temp; //r
	//inTreeFile >> temp; //e
	//inTreeFile >> temp; //e
	//inTreeFile >> temp; //H
	//inTreeFile >> temp; //e
	//inTreeFile >> temp; //i
	//inTreeFile >> temp; //g
	//inTreeFile >> temp; //h
	//inTreeFile >> temp; //t
	//inTreeFile >> temp; //:


	inTreeFile.seekg (strlen(" Tree Height : "), std::ios::cur);
	inTreeFile >> _treeHeight;

	int temp_int;
	for (int i = 0; i < treeSize; i++)
	{
		inTreeFile >> temp_int;
		inTreeFile >> temp_int;
		tree[i] = temp_int;
	}
	inTreeFile.close();

	return;
}
// callback function for APClustering; returning non-zero will abort apcluster
int callback(double *curr_a, double *curr_r, int N, int *curr_idx, int I, double curr_netsim, int iter)
{
	if(curr_netsim==-HUGE_VAL) 
		printf("%d\t%s\n", iter, "NaN"); /* cannot calculate netsim */
	else 
		printf("%d\t%f\n", iter, curr_netsim);
	return(0); /* 0=continue apcluster */
}

// do Affinity Propagation Clustering

// find a tree in a graph by merging one node on the longest path to its nearest neighbor with initial weights
void FindGraphMaxTreeWithInitialWeight(double** distanceMatrixTotal, int matrixSize, int* initialWeight, int* tree)
{
	/////////////////////////////////
	// initialize some variables

	int groupLevel = 1;
	bool isEnd = false;
	int idTemplate = -1;

	int groupSize = matrixSize;
	int** groupInfo = new int*[matrixSize];
	int*  subGroupSize = new int[matrixSize];
	int*  subGroupRep = new int[matrixSize];
	for (int i = 0; i < matrixSize; i++)
	{
		subGroupSize[i] = initialWeight[i];  // original == 1
		groupInfo[i] = new int[matrixSize];
		for (int j = 0; j < matrixSize; j++) groupInfo[i][j] = i;
		subGroupRep[i] = groupInfo[i][0];
	}

	// int groupSizeNew = matrixSize;
	int** groupInfoNew = new int*[matrixSize];
	int*  subGroupSizeNew = new int[matrixSize];
	int*  subGroupRepNew = new int[matrixSize];
	for (int i = 0; i < matrixSize; i++)
	{
		subGroupSizeNew[i] = initialWeight[i];  // original == 1
		groupInfoNew[i] = new int[matrixSize];
		for (int j = 0; j < matrixSize; j++) groupInfoNew[i][j] = i;
		subGroupRepNew[i] = groupInfoNew[i][0];
	}

	// start iteration
	while (!isEnd)
	{
		//////////////////////////////////////
		// extract current pair-wise distance matrix from distanceMatrixTotal
		double** distanceMatrix = new double*[groupSize];
		for (int i = 0; i < groupSize; i++) distanceMatrix[i] = new double[groupSize];
		for (int i = 0; i < groupSize; i++)
		{
			for (int j = 0; j < groupSize; j++)
			{
				int ii, jj;
				ii = subGroupRep[i];
				jj = subGroupRep[j];
				distanceMatrix[i][j] = distanceMatrixTotal[ii][jj];
				//distanceMatrix[i][j] = distanceMatrixTotal[subGroupRep[i]][subGroupRep[j]];
			}
		}

		// extract current pair-wise distance matrix from distanceMatrixTotal
		////////////////////////////////////
		// find the longest distance and locate the current moving subject and its neighbor
		// move idCur to idCurNeighbor

		int idCur = -1;
		int idCurNeighbor = -1; 

		if (groupSize == 1)  // add 20091215
		{
			tree[0] = -1;
		}
		else if (groupSize == 2)
		{
			// when only two subgroups, less move to more
			if (subGroupSize[0] > subGroupSize[1])
			{
				idCur = 1;
				idCurNeighbor = 0;
				idTemplate = subGroupRep[0];
			}
			else
			{
				idCur = 0;
				idCurNeighbor = 1;
				idTemplate = subGroupRep[1];
			}
			tree[subGroupRep[idCur]] = subGroupRep[idCurNeighbor];
			tree[subGroupRep[idCurNeighbor]] = -1;
		}
		else
		{
			FindMaxLocationInMatrix(distanceMatrix, groupSize, subGroupSize, idCur, idCurNeighbor);
			tree[subGroupRep[idCur]] = subGroupRep[idCurNeighbor];
		}
		// find the longest distance and locate the current moving subject and its neighbor
		//////////////////////////////
		// generate new group info
		int groupSizeNew = groupSize - 1;
		int index = 0;
		for (int i = 0; i < groupSize; i++)
		{
			if ((i != idCur) && (i != idCurNeighbor))
			{
				subGroupSizeNew[index] = subGroupSize[i];
				subGroupRepNew[index] = subGroupRep[i];
				groupInfoNew[index] = new int[subGroupSizeNew[index]];
				for (int j = 0; j < subGroupSizeNew[index]; j++)
				{
					groupInfoNew[index][j] = groupInfo[i][j];
				}
				index++;
			}
			else if (i == idCurNeighbor)
			{
				subGroupSizeNew[index] = subGroupSize[idCurNeighbor]+subGroupSize[idCur];
				subGroupRepNew[index] = subGroupRep[idCurNeighbor];
				groupInfoNew[index] = new int[subGroupSizeNew[index]];
				for (int j = 0; j < subGroupSize[idCurNeighbor]; j++)
				{
					groupInfoNew[index][j] = groupInfo[idCurNeighbor][j];
				}
				for (int j = 0; j < subGroupSize[idCur]; j++)
				{
					groupInfoNew[index][j+subGroupSize[idCurNeighbor]] = groupInfo[idCur][j];
				}
				index++;
			}
		}
		// generate new group info
		////////////////////////////////////
		// goto next level

		groupLevel ++;
		if (groupSize <= 2)  // update on 20091215 from groupSize == 2
			isEnd = true;

		// goto next level
		//////////////////////////////////
		// delete newed variables for current level

		for (int i = 0; i < groupSize; i++) delete[] distanceMatrix[i];
		delete[] distanceMatrix;

		// generate group information for next level
		for (int i = 0; i < groupSizeNew; i++)
		{
			subGroupRep[i] = subGroupRepNew[i];
			subGroupSize[i] = subGroupSizeNew[i];
			groupInfo[i] = new int[subGroupSize[i]];
			for (int j = 0; j < subGroupSize[i]; j++)
				groupInfo[i][j] = groupInfoNew[i][j];
		}
		groupSize = groupSizeNew;


	} // end of while (!isEnd)


	// start iteration
	///////////////////////////////////
	//// output tree
	//std::ofstream outfile;
	//outfile.open ("tree_output.txt");
	//outfile << matrixSize << std::endl;
	//for(int i = 0; i < matrixSize; i++)
	//{
	//	outfile <<  ' ' ;
	//	outfile << right << setw(3) << i; 
	//	outfile << ' ';
	//	outfile << right << setw(3) << tree[i];
	//	outfile << std::endl;
	//}
	//outfile.close();
	// output tree
	////////////////////////////////////
	// delete newed variables
	for (int i = 0; i < matrixSize; i++) delete[] groupInfo[i];
	delete[] groupInfo;
	delete[] subGroupSize;
	delete[] subGroupRep;

	for (int i = 0; i < matrixSize; i++) delete[] groupInfoNew[i];
	delete[] groupInfoNew;
	delete[] subGroupSizeNew;
	delete[] subGroupRepNew;

	return;
}


// Inverse deformation field by DG's algorithm
void InverseDeformationFieldDG3D(DeformationFieldType::Pointer deformationField, 
							   DeformationFieldType::Pointer &deformationFieldInverse)
{
	int SHIFT = 2;
	float OUTSIDE = 0;//100000.0;//0;
	int samplenum = 1;
	float interval = 1.0/(2*samplenum+1);
	unsigned int i, j, k;
	int x, y, z;
	float ii, jj, kk;
	float mdl_subvoxelx, mdl_subvoxely, mdl_subvoxelz;
	float disp_subvoxelx, disp_subvoxely, disp_subvoxelz;
	mdl_subvoxelx = 0.0; mdl_subvoxely = 0.0; mdl_subvoxelz = 0.0; disp_subvoxelx = 0.0; disp_subvoxely = 0.0; disp_subvoxelz = 0.0;
	DeformationFieldIteratorType dfNewIt ( deformationFieldInverse, deformationFieldInverse->GetRequestedRegion() );
	DeformationFieldIteratorType dfIt ( deformationField, deformationField->GetRequestedRegion() );

	// initialize three 3D float matrix to store deformation field
	//unsigned int image_size = deformationField->GetRequestedRegion().GetSize()[0];
	unsigned int x_size = deformationField->GetRequestedRegion().GetSize()[0];
	unsigned int y_size = deformationField->GetRequestedRegion().GetSize()[1];
	unsigned int z_size = deformationField->GetRequestedRegion().GetSize()[2];
	
	float*** dfx = new float**[x_size];
	float*** dfy = new float**[x_size];
	float*** dfz = new float**[x_size];
	for(i = 0; i < x_size; i++) 
	{
		dfx[i] = new float*[y_size];
		dfy[i] = new float*[y_size];
		dfz[i] = new float*[y_size];
	}
	for(i = 0; i < x_size; i++) 
	{
		for (j = 0; j < y_size; j++)
		{
			dfx[i][j] = new float[z_size];
			dfy[i][j] = new float[z_size];
			dfz[i][j] = new float[z_size];
		}
	}

	// load deformationFieldBA into 3 3D matrix dfx and dfy and dfz
	VectorPixelType vectorPixel;
	DeformationFieldType::IndexType idx;
	//for (i = 0, j = 0, k = 0, dfIt.GoToBegin(); !dfIt.IsAtEnd(); ++dfIt)
	for (dfIt.GoToBegin(); !dfIt.IsAtEnd(); ++dfIt)
	{
		vectorPixel = dfIt.Get();
		idx = dfIt.GetIndex();
		dfx[idx[0]][idx[1]][idx[2]] = vectorPixel.GetElement(0);
		dfy[idx[0]][idx[1]][idx[2]] = vectorPixel.GetElement(1);
		dfz[idx[0]][idx[1]][idx[2]] = vectorPixel.GetElement(2);
		//dfx[i][j][k] = vectorPixel.GetElement(0);
		//dfy[i][j][k] = vectorPixel.GetElement(1);
		//dfz[i][j][k] = vectorPixel.GetElement(2);
		//k++;
		//if (k == z_size)
		//{
		//	j++;
		//	k = 0;
		//	if (j == image_size)
		//	{
		//		i++;
		//		j = 0;
		//	}
		//}
		//i++;
		//if (i == image_size)
		//{
		//	j++;
		//	i = 0;
		//	if (j == image_size)
		//	{
		//		k++;
		//		j = 0;
		//	}
		//}
		//if ((i == image_size-1) && (j == image_size-1) && (k == z_size-1))
		//{
		//	std::cerr << "read end!" << std::endl;
		//}
	}

	// allocate some internal data matrices, weights matrix and	enlarged inverse df matrix with borders
	float*** totalweights = new float**[x_size+2*SHIFT];
	float***	    rdfbx = new float**[x_size+2*SHIFT];
	float***	    rdfby = new float**[x_size+2*SHIFT];
	float***	    rdfbz = new float**[x_size+2*SHIFT];
	for(i = 0; i < x_size+2*SHIFT; i++)
	{
		totalweights[i] = new float*[y_size+2*SHIFT];
			   rdfbx[i] = new float*[y_size+2*SHIFT];
			   rdfby[i] = new float*[y_size+2*SHIFT];
			   rdfbz[i] = new float*[y_size+2*SHIFT];
	}
	for(i = 0; i < x_size+2*SHIFT; i++) 
	{
		for (j = 0; j < y_size+2*SHIFT; j++)
		{
			totalweights[i][j] = new float[z_size+2*SHIFT];
				   rdfbx[i][j] = new float[z_size+2*SHIFT];
				   rdfby[i][j] = new float[z_size+2*SHIFT];
				   rdfbz[i][j] = new float[z_size+2*SHIFT];
		}
	}

	// initialize these matrices
	for (i = 0; i < x_size+2*SHIFT; i++)
	{	
		for (j = 0; j < y_size+2*SHIFT; j++)
		{	
			for (k = 0; k < z_size+2*SHIFT; k++)
			{
				totalweights[i][j][k] = 0.0;	
					   rdfbx[i][j][k] = 0.0;
					   rdfby[i][j][k] = 0.0;
					   rdfbz[i][j][k] = 0.0;
			}
		}
	}

	// estimating
	for (i = 0; i < x_size; i++)
	{
		// std::cerr << i << ", ";
		for (j = 0; j < y_size; j++)
		{
			// std::cerr << "(" << i << ", " << j << "), " ;
			for (k = 0; k < z_size; k++)
			{
				for (x = -samplenum; x <= samplenum; x++)
				{
					for (y = -samplenum; y <= samplenum; y++)
					{
						for (z = -samplenum; z <= samplenum; z++)
						{
							mdl_subvoxelx = x*interval + i;
							mdl_subvoxely = y*interval + j;
							mdl_subvoxelz = z*interval + k;

							// call interpolateDisplacement
							{
								int ni, nj, nk, nip1, njp1, nkp1;
								float b, c, d, b1, c1, d1;
								ni = (int)mdl_subvoxelx; nip1 = ni + 1;
								nj = (int)mdl_subvoxely; njp1 = nj + 1;
								nk = (int)mdl_subvoxelz; nkp1 = nk + 1;

								if((ni>=0)&&(ni<(int)x_size-1)&&(nj>=0)&&(nj<(int)y_size-1)&&(nk>=0)&&(nk<(int)z_size-1))
								{
									b = mdl_subvoxelx - ni; b1 = 1.0-b;
									c = mdl_subvoxely - nj; c1 = 1.0-c;
									d = mdl_subvoxelz - nk; d1 = 1.0-d;

									disp_subvoxelx = (d1*( dfx[ni][nj][nk]*(b1*c1) + dfx[nip1][nj][nk]*(b*c1)+ 
										dfx[ni][njp1][nk]*(b1*c)+ dfx[nip1][njp1][nk]*(b*c) ) +
										d*(dfx[ni][nj][nkp1]*(b1*c1) + dfx[nip1][nj][nkp1]*(b*c1)+ 
										dfx[ni][njp1][nkp1]*(b1*c)+ dfx[nip1][njp1][nkp1]*(b*c)))
										/(d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c))) ; 
									disp_subvoxely = (d1*( dfy[ni][nj][nk]*(b1*c1) + dfy[nip1][nj][nk]*(b*c1)+ 
										dfy[ni][njp1][nk]*(b1*c)+ dfy[nip1][njp1][nk]*(b*c) ) +
										d*(dfy[ni][nj][nkp1]*(b1*c1) + dfy[nip1][nj][nkp1]*(b*c1)+ 
										dfy[ni][njp1][nkp1]*(b1*c)+ dfy[nip1][njp1][nkp1]*(b*c)))
										/(d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c))) ; 
									disp_subvoxelz = (d1*( dfz[ni][nj][nk]*(b1*c1) + dfz[nip1][nj][nk]*(b*c1)+ 
										dfz[ni][njp1][nk]*(b1*c)+ dfz[nip1][njp1][nk]*(b*c) ) +
										d*(dfz[ni][nj][nkp1]*(b1*c1) + dfz[nip1][nj][nkp1]*(b*c1)+ 
										dfz[ni][njp1][nkp1]*(b1*c)+ dfz[nip1][njp1][nkp1]*(b*c)))
										/(d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c))) ; 
									//disp_subvoxely = ( dfy[ni][nj]*(b1*c1) + 
									//	dfy[nip1][nj]*(b*c1)+ 
									//	dfy[ni][njp1]*(b1*c)+ 
									//	dfy[nip1][njp1]*(b*c) )/( (b1*c1)+(b*c1)+(b1*c)+(b*c) ) ; 

								}
								else if (((ni==(int)x_size-1)&&(nj>=0)&&(nj<(int)y_size-1)  &&(nk>=0)&&(nk<(int)z_size-1))
									   ||((ni>=0)&&(ni<(int)x_size-1) &&(nj==(int)y_size-1) &&(nk>=0)&&(nk<(int)z_size-1))
									   ||((ni>=0)&&(ni<(int)x_size-1) &&(nj>=0)&&(nj<(int)y_size-1)  &&(nk=(int)z_size-1)))
								{
									disp_subvoxelx = dfx[ni][nj][nk];
									disp_subvoxely = dfy[ni][nj][nk];
									disp_subvoxelz = dfz[ni][nj][nk];
								}

							}

							ii = mdl_subvoxelx + disp_subvoxelx;
							jj = mdl_subvoxely + disp_subvoxely;
							kk = mdl_subvoxelz + disp_subvoxelz;

							// call iterativeEstimate
							{
								int ni,nj,nk,nip1,njp1,nkp1;
								float b,c,d,b1,c1,d1,combined,weight;
								ii += SHIFT; jj += SHIFT; kk += SHIFT;
								ni = (int)ii; nip1 = ni+1;
								nj = (int)jj; njp1 = nj+1;
								nk = (int)kk; nkp1 = nk+1;
								if(ni>=0 && ni<(int)x_size+2*SHIFT-1  &&  nj>=0 && nj<(int)y_size+2*SHIFT-1 &&  nk>=0 && nk<(int)z_size+2*SHIFT-1)
								{
									b = ii-ni ;        b1 = 1.0-b ;
									c = jj-nj ;        c1 = 1.0-c ;
									d = kk-nk ;        d1 = 1.0-d ;
									combined = d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c));

									weight = d1*(b1*c1)/combined ;
									totalweights[ni][nj][nk] += weight ;
										   rdfbx[ni][nj][nk] += disp_subvoxelx*weight ;
										   rdfby[ni][nj][nk] += disp_subvoxely*weight ;
										   rdfbz[ni][nj][nk] += disp_subvoxelz*weight ;

									weight = d1*(b*c1)/combined ;
									totalweights[nip1][nj][nk] += weight ;
										   rdfbx[nip1][nj][nk] += disp_subvoxelx*weight ;
										   rdfby[nip1][nj][nk] += disp_subvoxely*weight ;
										   rdfbz[nip1][nj][nk] += disp_subvoxelz*weight ;

									weight = d1*(b1*c)/combined ;
									totalweights[ni][njp1][nk] += weight ;
										   rdfbx[ni][njp1][nk] += disp_subvoxelx*weight ;
										   rdfby[ni][njp1][nk] += disp_subvoxely*weight ;
										   rdfbz[ni][njp1][nk] += disp_subvoxelz*weight ;

									weight = d1*(b*c)/combined ;
									totalweights[nip1][njp1][nk]  += weight ;
										   rdfbx[nip1][njp1][nk] += disp_subvoxelx*weight ;
										   rdfby[nip1][njp1][nk] += disp_subvoxely*weight ;
										   rdfbz[nip1][njp1][nk] += disp_subvoxelz*weight ;

									weight = d*(b1*c1)/combined ;
									totalweights[ni][nj][nkp1] += weight ;
										   rdfbx[ni][nj][nkp1] += disp_subvoxelx*weight ;
										   rdfby[ni][nj][nkp1] += disp_subvoxely*weight ;
										   rdfbz[ni][nj][nkp1] += disp_subvoxelz*weight ;

									weight = d*(b*c1)/combined ;
									totalweights[nip1][nj][nkp1] += weight ;
										   rdfbx[nip1][nj][nkp1] += disp_subvoxelx*weight ;
										   rdfby[nip1][nj][nkp1] += disp_subvoxely*weight ;
										   rdfbz[nip1][nj][nkp1] += disp_subvoxelz*weight ;

									weight = d*(b1*c)/combined ;
									totalweights[ni][njp1][nkp1] += weight ;
										   rdfbx[ni][njp1][nkp1] += disp_subvoxelx*weight ;
										   rdfby[ni][njp1][nkp1] += disp_subvoxely*weight ;
										   rdfbz[ni][njp1][nkp1] += disp_subvoxelz*weight ;

									weight = d*(b*c)/combined ;
									totalweights[nip1][njp1][nkp1] += weight ;
										   rdfbx[nip1][njp1][nkp1] += disp_subvoxelx*weight ;
										   rdfby[nip1][njp1][nkp1] += disp_subvoxely*weight ;
										   rdfbz[nip1][njp1][nkp1] += disp_subvoxelz*weight ;
								}

							}
						}// end for z
					}// end for y
				}// end for z
			}// end for k
		}// end for j
	}// end for i

	// allocate inverse deformation field
	float*** rdfx = new float**[x_size];
	float*** rdfy = new float**[x_size];
	float*** rdfz = new float**[x_size];
	for(i = 0; i<x_size; i++)
	{
		rdfx[i] = new float*[y_size];
		rdfy[i] = new float*[y_size];
		rdfz[i] = new float*[y_size];
	}
	for(i = 0; i<x_size; i++) 
	{
		for (j = 0; j < y_size; j++)
		{
			rdfx[i][j] = new float[z_size];
			rdfy[i][j] = new float[z_size];
			rdfz[i][j] = new float[z_size];
		}
	}
	int count_outside = 0;
	// normalize the enlarged rdfb to rdf
	for (i = 0; i< x_size; i++)
	{
		for (j = 0; j < y_size; j++)
		{
			for (k = 0; k < z_size; k++)
			{
				if (totalweights[i+SHIFT][j+SHIFT][k+SHIFT]>0)
				{
					rdfx[i][j][k] = rdfbx[i+SHIFT][j+SHIFT][k+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT][k+SHIFT]);
					rdfy[i][j][k] = rdfby[i+SHIFT][j+SHIFT][k+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT][k+SHIFT]);
					rdfz[i][j][k] = rdfbz[i+SHIFT][j+SHIFT][k+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT][k+SHIFT]);
				}
				else
				{
					rdfx[i][j][k] = OUTSIDE;
					rdfy[i][j][k] = OUTSIDE;
					rdfz[i][j][k] = OUTSIDE;
					count_outside++;
				}
			}
		}
	}
	//for (i = 0, j = 0, k = 0, dfNewIt.GoToBegin(); !dfNewIt.IsAtEnd(); ++dfNewIt)
	for (dfNewIt.GoToBegin(); !dfNewIt.IsAtEnd(); ++dfNewIt)
	{
		idx = dfNewIt.GetIndex();
		vectorPixel.SetElement(0,rdfx[idx[0]][idx[1]][idx[2]]);
		vectorPixel.SetElement(1,rdfy[idx[0]][idx[1]][idx[2]]);
		vectorPixel.SetElement(2,rdfz[idx[0]][idx[1]][idx[2]]);
		dfNewIt.Set(vectorPixel);
		//k++;
		//if (k == z_size)
		//{
		//	j++;
		//	k = 0;
		//	if (j == image_size)
		//	{
		//		i++;
		//		j = 0;
		//	}
		//}
		//i++;
		//if (i == image_size)
		//{
		//	j++;
		//	i = 0;
		//	if (j == image_size)
		//	{
		//		k++;
		//		j = 0;
		//	}
		//}
		//j++;
		//if (j==image_size)
		//{
		//	i++;
		//	j=0;
		//}
	}

	// free memory
	for(i = 0; i<x_size; i++)
	{
		for(j = 0; j<y_size; j++)
		{
			delete[] dfx[i][j];			delete[] dfy[i][j];			delete[] dfz[i][j];
			delete[] rdfx[i][j];		delete[] rdfy[i][j];		delete[] rdfz[i][j];
		}
		delete[] dfx[i];		delete[] dfy[i];		delete[] dfz[i];
		delete[] rdfx[i];		delete[] rdfy[i];		delete[] rdfz[i];
	}
	for(i = 0; i<x_size+2*SHIFT; i++)
	{
		for(j = 0; j<y_size+2*SHIFT; j++)
		{
			delete[] rdfbx[i][j];			delete[] rdfby[i][j];			delete[] rdfbz[i][j];
			delete[] totalweights[i][j];
		}
		delete[] rdfbx[i];		delete[] rdfby[i];		delete[] rdfbz[i];
		delete[] totalweights[i]; 
	}
	delete[] dfx;	delete[] dfy;	delete[] dfz;
	delete[] rdfx;	delete[] rdfy;	delete[] rdfz;
	delete[] rdfbx;	delete[] rdfby;	delete[] rdfbz;
	delete[] totalweights;

}


void InverseDeformationFieldDG3D(float*** dfx, float*** dfy, float*** dfz, 
								 int x_size, int y_size, int z_size, 
								 float*** rdfx, float*** rdfy, float*** rdfz)
{
	int SHIFT = 2;
	float OUTSIDE = 0;//100000.0;//0;
	int samplenum = 1;
	float interval = 1.0/(2*samplenum+1);
	int i, j, k;
	int x, y, z;
	float ii, jj, kk;
	float mdl_subvoxelx, mdl_subvoxely, mdl_subvoxelz;
	float disp_subvoxelx, disp_subvoxely, disp_subvoxelz;
	mdl_subvoxelx = 0.0; mdl_subvoxely = 0.0; mdl_subvoxelz = 0.0; disp_subvoxelx = 0.0; disp_subvoxely = 0.0; disp_subvoxelz = 0.0;
/*
	DeformationFieldIteratorType dfNewIt ( deformationFieldInverse, deformationFieldInverse->GetRequestedRegion() );
	DeformationFieldIteratorType dfIt ( deformationField, deformationField->GetRequestedRegion() );

	// initialize three 3D float matrix to store deformation field
	unsigned int image_size = deformationField->GetRequestedRegion().GetSize()[0];
	unsigned int z_size =  deformationField->GetRequestedRegion().GetSize()[2];

	float*** dfx = new float**[image_size];
	float*** dfy = new float**[image_size];
	float*** dfz = new float**[image_size];
	for(i = 0; i < image_size; i++) 
	{
		dfx[i] = new float*[image_size];
		dfy[i] = new float*[image_size];
		dfz[i] = new float*[image_size];
	}
	for(i = 0; i < image_size; i++) 
	{
		for (j = 0; j < image_size; j++)
		{
			dfx[i][j] = new float[z_size];
			dfy[i][j] = new float[z_size];
			dfz[i][j] = new float[z_size];
		}
	}

	// load deformationFieldBA into 3 3D matrix dfx and dfy and dfz
	VectorPixelType vectorPixel;
	DeformationFieldType::IndexType idx;
	//for (i = 0, j = 0, k = 0, dfIt.GoToBegin(); !dfIt.IsAtEnd(); ++dfIt)
	for (dfIt.GoToBegin(); !dfIt.IsAtEnd(); ++dfIt)
	{
		vectorPixel = dfIt.Get();
		idx = dfIt.GetIndex();
		dfx[idx[0]][idx[1]][idx[2]] = vectorPixel.GetElement(0);
		dfy[idx[0]][idx[1]][idx[2]] = vectorPixel.GetElement(1);
		dfz[idx[0]][idx[1]][idx[2]] = vectorPixel.GetElement(2);
		//dfx[i][j][k] = vectorPixel.GetElement(0);
		//dfy[i][j][k] = vectorPixel.GetElement(1);
		//dfz[i][j][k] = vectorPixel.GetElement(2);
		//k++;
		//if (k == z_size)
		//{
		//	j++;
		//	k = 0;
		//	if (j == image_size)
		//	{
		//		i++;
		//		j = 0;
		//	}
		//}
		//i++;
		//if (i == image_size)
		//{
		//	j++;
		//	i = 0;
		//	if (j == image_size)
		//	{
		//		k++;
		//		j = 0;
		//	}
		//}
		//if ((i == image_size-1) && (j == image_size-1) && (k == z_size-1))
		//{
		//	std::cerr << "read end!" << std::endl;
		//}
	}
*/
	// allocate some internal data matrices, weights matrix and	enlarged inverse df matrix with borders
	float*** totalweights = new float**[x_size+2*SHIFT];
	float***	    rdfbx = new float**[x_size+2*SHIFT];
	float***	    rdfby = new float**[x_size+2*SHIFT];
	float***	    rdfbz = new float**[x_size+2*SHIFT];
	for(i = 0; i < x_size+2*SHIFT; i++)
	{
		totalweights[i] = new float*[y_size+2*SHIFT];
		rdfbx[i] = new float*[y_size+2*SHIFT];
		rdfby[i] = new float*[y_size+2*SHIFT];
		rdfbz[i] = new float*[y_size+2*SHIFT];
	}
	for(i = 0; i < x_size+2*SHIFT; i++) 
	{
		for (j = 0; j < y_size+2*SHIFT; j++)
		{
			totalweights[i][j] = new float[z_size+2*SHIFT];
			rdfbx[i][j] = new float[z_size+2*SHIFT];
			rdfby[i][j] = new float[z_size+2*SHIFT];
			rdfbz[i][j] = new float[z_size+2*SHIFT];
		}
	}

	// initialize these matrices
	for (i = 0; i < x_size+2*SHIFT; i++)
	{	
		for (j = 0; j < y_size+2*SHIFT; j++)
		{	
			for (k = 0; k < z_size+2*SHIFT; k++)
			{
				totalweights[i][j][k] = 0.0;	
				rdfbx[i][j][k] = 0.0;
				rdfby[i][j][k] = 0.0;
				rdfbz[i][j][k] = 0.0;
			}
		}
	}

	// estimating
	for (k = 0; k < z_size; k++)
	{
		std::cerr << k << ", ";
		for (j = 0; j < y_size; j++)
		{
			// std::cerr << "(" << i << ", " << j << "), " ;
			for (i = 0; i < x_size; i++)
			{
				for (z = -samplenum; z <= samplenum; z++)
				{
					for (y = -samplenum; y <= samplenum; y++)
					{
						for (x = -samplenum; x <= samplenum; x++)
						{
							mdl_subvoxelx = x*interval + i;
							mdl_subvoxely = y*interval + j;
							mdl_subvoxelz = z*interval + k;

							// call interpolateDisplacement
							{
								int ni, nj, nk, nip1, njp1, nkp1;
								float b, c, d, b1, c1, d1;
								ni = (int)mdl_subvoxelx; nip1 = ni + 1;
								nj = (int)mdl_subvoxely; njp1 = nj + 1;
								nk = (int)mdl_subvoxelz; nkp1 = nk + 1;

								if((ni>=0)&&(ni<x_size-1)&&(nj>=0)&&(nj<y_size-1)&&(nk>=0)&&(nk<z_size-1))
								{
									b = mdl_subvoxelx - ni; b1 = 1.0-b;
									c = mdl_subvoxely - nj; c1 = 1.0-c;
									d = mdl_subvoxelz - nk; d1 = 1.0-d;

									disp_subvoxelx = (d1*( dfx[ni][nj][nk]*(b1*c1) + dfx[nip1][nj][nk]*(b*c1)+ 
										dfx[ni][njp1][nk]*(b1*c)+ dfx[nip1][njp1][nk]*(b*c) ) +
										d*(dfx[ni][nj][nkp1]*(b1*c1) + dfx[nip1][nj][nkp1]*(b*c1)+ 
										dfx[ni][njp1][nkp1]*(b1*c)+ dfx[nip1][njp1][nkp1]*(b*c)))
										/(d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c))) ; 
									disp_subvoxely = (d1*( dfy[ni][nj][nk]*(b1*c1) + dfy[nip1][nj][nk]*(b*c1)+ 
										dfy[ni][njp1][nk]*(b1*c)+ dfy[nip1][njp1][nk]*(b*c) ) +
										d*(dfy[ni][nj][nkp1]*(b1*c1) + dfy[nip1][nj][nkp1]*(b*c1)+ 
										dfy[ni][njp1][nkp1]*(b1*c)+ dfy[nip1][njp1][nkp1]*(b*c)))
										/(d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c))) ; 
									disp_subvoxelz = (d1*( dfz[ni][nj][nk]*(b1*c1) + dfz[nip1][nj][nk]*(b*c1)+ 
										dfz[ni][njp1][nk]*(b1*c)+ dfz[nip1][njp1][nk]*(b*c) ) +
										d*(dfz[ni][nj][nkp1]*(b1*c1) + dfz[nip1][nj][nkp1]*(b*c1)+ 
										dfz[ni][njp1][nkp1]*(b1*c)+ dfz[nip1][njp1][nkp1]*(b*c)))
										/(d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c))) ; 
									//disp_subvoxely = ( dfy[ni][nj]*(b1*c1) + 
									//	dfy[nip1][nj]*(b*c1)+ 
									//	dfy[ni][njp1]*(b1*c)+ 
									//	dfy[nip1][njp1]*(b*c) )/( (b1*c1)+(b*c1)+(b1*c)+(b*c) ) ; 

								}
								else if (((ni==x_size-1)&&(nj>=0)&&(nj<y_size-1)  &&(nk>=0)&&(nk<z_size-1))
									||((ni>=0)&&(ni<x_size-1) &&(nj==y_size-1) &&(nk>=0)&&(nk<z_size-1))
									||((ni>=0)&&(ni<x_size-1) &&(nj>=0)&&(nj<y_size-1)  &&(nk=z_size-1)))
								{
									disp_subvoxelx = dfx[ni][nj][nk];
									disp_subvoxely = dfy[ni][nj][nk];
									disp_subvoxelz = dfz[ni][nj][nk];
								}

							}

							ii = mdl_subvoxelx + disp_subvoxelx;
							jj = mdl_subvoxely + disp_subvoxely;
							kk = mdl_subvoxelz + disp_subvoxelz;

							// call iterativeEstimate
							{
								int ni,nj,nk,nip1,njp1,nkp1;
								float b,c,d,b1,c1,d1,combined,weight;
								ii += SHIFT; jj += SHIFT; kk += SHIFT;
								ni = (int)ii; nip1 = ni+1;
								nj = (int)jj; njp1 = nj+1;
								nk = (int)kk; nkp1 = nk+1;
								if(ni>=0 && ni<x_size+2*SHIFT-1  &&  nj>=0 && nj<y_size+2*SHIFT-1 &&  nk>=0 && nk<z_size+2*SHIFT-1)
								{
									b = ii-ni ;        b1 = 1.0-b ;
									c = jj-nj ;        c1 = 1.0-c ;
									d = kk-nk ;        d1 = 1.0-d ;
									combined = d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c));

									weight = d1*(b1*c1)/combined ;
									totalweights[ni][nj][nk] += weight ;
									rdfbx[ni][nj][nk] += disp_subvoxelx*weight ;
									rdfby[ni][nj][nk] += disp_subvoxely*weight ;
									rdfbz[ni][nj][nk] += disp_subvoxelz*weight ;

									weight = d1*(b*c1)/combined ;
									totalweights[nip1][nj][nk] += weight ;
									rdfbx[nip1][nj][nk] += disp_subvoxelx*weight ;
									rdfby[nip1][nj][nk] += disp_subvoxely*weight ;
									rdfbz[nip1][nj][nk] += disp_subvoxelz*weight ;

									weight = d1*(b1*c)/combined ;
									totalweights[ni][njp1][nk] += weight ;
									rdfbx[ni][njp1][nk] += disp_subvoxelx*weight ;
									rdfby[ni][njp1][nk] += disp_subvoxely*weight ;
									rdfbz[ni][njp1][nk] += disp_subvoxelz*weight ;

									weight = d1*(b*c)/combined ;
									totalweights[nip1][njp1][nk]  += weight ;
									rdfbx[nip1][njp1][nk] += disp_subvoxelx*weight ;
									rdfby[nip1][njp1][nk] += disp_subvoxely*weight ;
									rdfbz[nip1][njp1][nk] += disp_subvoxelz*weight ;

									weight = d*(b1*c1)/combined ;
									totalweights[ni][nj][nkp1] += weight ;
									rdfbx[ni][nj][nkp1] += disp_subvoxelx*weight ;
									rdfby[ni][nj][nkp1] += disp_subvoxely*weight ;
									rdfbz[ni][nj][nkp1] += disp_subvoxelz*weight ;

									weight = d*(b*c1)/combined ;
									totalweights[nip1][nj][nkp1] += weight ;
									rdfbx[nip1][nj][nkp1] += disp_subvoxelx*weight ;
									rdfby[nip1][nj][nkp1] += disp_subvoxely*weight ;
									rdfbz[nip1][nj][nkp1] += disp_subvoxelz*weight ;

									weight = d*(b1*c)/combined ;
									totalweights[ni][njp1][nkp1] += weight ;
									rdfbx[ni][njp1][nkp1] += disp_subvoxelx*weight ;
									rdfby[ni][njp1][nkp1] += disp_subvoxely*weight ;
									rdfbz[ni][njp1][nkp1] += disp_subvoxelz*weight ;

									weight = d*(b*c)/combined ;
									totalweights[nip1][njp1][nkp1] += weight ;
									rdfbx[nip1][njp1][nkp1] += disp_subvoxelx*weight ;
									rdfby[nip1][njp1][nkp1] += disp_subvoxely*weight ;
									rdfbz[nip1][njp1][nkp1] += disp_subvoxelz*weight ;
								}

							}
						}// end for z
					}// end for y
				}// end for z
			}// end for k
		}// end for j
	}// end for i
/*
	// allocate inverse deformation field
	float*** rdfx = new float**[image_size];
	float*** rdfy = new float**[image_size];
	float*** rdfz = new float**[image_size];
	for(i = 0; i<image_size; i++)
	{
		rdfx[i] = new float*[image_size];
		rdfy[i] = new float*[image_size];
		rdfz[i] = new float*[image_size];
	}
	for(i = 0; i<image_size; i++) 
	{
		for (j = 0; j < image_size; j++)
		{
			rdfx[i][j] = new float[z_size];
			rdfy[i][j] = new float[z_size];
			rdfz[i][j] = new float[z_size];
		}
	}
*/
	int count_outside = 0;
	// normalize the enlarged rdfb to rdf
	for (i = 0; i< x_size; i++)
	{
		for (j = 0; j < y_size; j++)
		{
			for (k = 0; k < z_size; k++)
			{
				if (totalweights[i+SHIFT][j+SHIFT][k+SHIFT]>0)
				{
					rdfx[i][j][k] = rdfbx[i+SHIFT][j+SHIFT][k+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT][k+SHIFT]);
					rdfy[i][j][k] = rdfby[i+SHIFT][j+SHIFT][k+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT][k+SHIFT]);
					rdfz[i][j][k] = rdfbz[i+SHIFT][j+SHIFT][k+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT][k+SHIFT]);
				}
				else
				{
					rdfx[i][j][k] = OUTSIDE;
					rdfy[i][j][k] = OUTSIDE;
					rdfz[i][j][k] = OUTSIDE;
					count_outside++;
				}
			}
		}
	}
/*
	//for (i = 0, j = 0, k = 0, dfNewIt.GoToBegin(); !dfNewIt.IsAtEnd(); ++dfNewIt)
	for (dfNewIt.GoToBegin(); !dfNewIt.IsAtEnd(); ++dfNewIt)
	{
		idx = dfNewIt.GetIndex();
		vectorPixel.SetElement(0,rdfx[idx[0]][idx[1]][idx[2]]);
		vectorPixel.SetElement(1,rdfy[idx[0]][idx[1]][idx[2]]);
		vectorPixel.SetElement(2,rdfz[idx[0]][idx[1]][idx[2]]);
		dfNewIt.Set(vectorPixel);
		//k++;
		//if (k == z_size)
		//{
		//	j++;
		//	k = 0;
		//	if (j == image_size)
		//	{
		//		i++;
		//		j = 0;
		//	}
		//}
		//i++;
		//if (i == image_size)
		//{
		//	j++;
		//	i = 0;
		//	if (j == image_size)
		//	{
		//		k++;
		//		j = 0;
		//	}
		//}
		//j++;
		//if (j==image_size)
		//{
		//	i++;
		//	j=0;
		//}
	}
*/
	// free memory
	for(i = 0; i<x_size+2*SHIFT; i++)
	{
		for(j = 0; j<y_size+2*SHIFT; j++)
		{
			delete[] rdfbx[i][j];			delete[] rdfby[i][j];			delete[] rdfbz[i][j];
			delete[] totalweights[i][j];
		}
		delete[] rdfbx[i];		delete[] rdfby[i];		delete[] rdfbz[i];
		delete[] totalweights[i]; 
	}
	delete[] rdfbx;	delete[] rdfby;	delete[] rdfbz;
	delete[] totalweights;
/*
	for(i = 0; i<image_size; i++)
	{
		for(j = 0; j<image_size; j++)
		{
			delete[] dfx[i][j];			delete[] dfy[i][j];			delete[] dfz[i][j];
			delete[] rdfx[i][j];		delete[] rdfy[i][j];		delete[] rdfz[i][j];
		}
		delete[] dfx[i];		delete[] dfy[i];		delete[] dfz[i];
		delete[] rdfx[i];		delete[] rdfy[i];		delete[] rdfz[i];
	}
	for(i = 0; i<image_size+2*SHIFT; i++)
	{
		for(j = 0; j<image_size+2*SHIFT; j++)
		{
			delete[] rdfbx[i][j];			delete[] rdfby[i][j];			delete[] rdfbz[i][j];
			delete[] totalweights[i][j];
		}
		delete[] rdfbx[i];		delete[] rdfby[i];		delete[] rdfbz[i];
		delete[] totalweights[i]; 
	}
	delete[] dfx;	delete[] dfy;	delete[] dfz;
	delete[] rdfx;	delete[] rdfy;	delete[] rdfz;
	delete[] rdfbx;	delete[] rdfby;	delete[] rdfbz;
	delete[] totalweights;
*/
}


/*
// do Affinity Propagation Clustering on Unix64 system
void APClusteringUnix64(double** distanceMatrix, int groupSize, int &APclustersize, 
						int* &APclustersubrep, int* &APclustersubsize, int** &APcluster)
{
	//
	//std::cout << "groupsize: " << groupSize << std::endl;
	// for APcluster
	//double sij[]={-15.561,-1.861,-4.066,-9.293,-3.581,-3.093,-27.967,-32.147,-44.921,-47.967,-17.665,-15.628,-6.430,-36.111,-21.062,-32.094,-46.606,-65.053,-64.292,-42.983,-22.071,-11.619,-35.921,-51.256,-19.754,	-1.861,-15.561,-1.839,-8.387,-8.612,-4.909,-15.438,-19.376,-28.647,-31.421,-9.245,-6.721,-1.980,-22.467,-10.440,-22.388,-38.132,-52.074,-49.382,-31.872,-16.521,-6.795,-24.393,-36.228,-10.369,-4.066,-1.839,-15.561,-2.561,-7.097,-2.363,-17.266,-25.109,-31.390,-27.590,-5.415,-8.475,-6.370,-17.876,-10.448,-13.314,-23.992,-36.746,-36.043,-20.618,-7.652,-1.947,-15.996,-27.065,-15.957,-9.293,-8.387,-2.561,-15.561,-6.642,-2.241,-28.773,-41.141,-45.675,-34.572,-8.612,-18.211,-17.008,-22.646,-19.008,-11.184,-15.561,-29.072,-31.518,-16.679,-3.816,-2.318,-15.480,-27.314,-30.247,-3.581,-8.612,-7.097,-6.642,-15.561,-1.458,-44.952,-53.686,-66.515,-62.544,-24.441,-29.052,-18.491,-47.060,-34.200,-34.301,-42.021,-63.480,-66.777,-44.127,-20.513,-13.814,-40.627,-58.475,-37.880,-3.093,-4.909,-2.363,-2.241,-1.458,-15.561,-32.239,-41.536,-50.820,-45.330,-14.087,-19.384,-13.103,-32.147,-22.741,-21.832,-29.598,-47.020,-49.076,-30.043,-11.593,-6.297,-26.753,-41.495,-28.560,-27.967,-15.438,-17.266,-28.773,-44.952,-32.239,-15.561,-2.126,-2.110,-5.972,-7.313,-1.788,-9.199,-5.633,-1.042,-18.227,-38.981,-40.223,-31.372,-24.265,-24.331,-15.561,-14.764,-16.689,-2.962,-32.147,-19.376,-25.109,-41.141,-53.686,-41.536,-2.126,-15.561,-2.609,-13.515,-16.243,-4.613,-9.944,-14.491,-5.458,-32.592,-59.021,-60.772,-49.446,-40.740,-39.186,-26.203,-28.088,-30.241,-1.572,-44.921,-28.647,-31.390,-45.675,-66.515,-50.820,-2.110,-2.609,-15.561,-5.895,-15.692,-7.693,-18.816,-8.843,-5.778,-27.112,-50.933,-47.903,-36.342,-32.455,-37.470,-28.040,-21.515,-20.068,-7.021,-47.967,-31.421,-27.590,-34.572,-62.544,-45.330,-5.972,-13.515,-5.895,-15.561,-8.956,-11.047,-26.444,-1.281,-6.170,-11.085,-25.670,-20.974,-13.189,-12.655,-21.379,-19.102,-6.893,-4.432,-17.315,-17.665,-9.245,-5.415,-8.612,-24.441,-14.087,-7.313,-16.243,-15.692,-8.956,-15.561,-4.386,-10.533,-3.674,-2.907,-3.817,-14.577,-19.552,-16.461,-8.008,-4.993,-1.995,-3.791,-8.877,-12.546,-15.628,-6.721,-8.475,-18.211,-29.052,-19.384,-1.788,-4.613,-7.693,-11.047,-4.386,-15.561,-3.332,-7.760,-0.737,-15.948,-34.891,-40.428,-33.715,-23.112,-17.983,-8.972,-14.360,-19.640,-2.097,-6.430,-1.980,-6.370,-17.008,-18.491,-13.103,-9.199,-9.944,-18.816,-26.444,-10.533,-3.332,-15.561,-20.313,-7.075,-26.722,-47.318,-58.700,-52.987,-36.801,-23.903,-11.639,-26.800,-36.526,-3.661,-36.111,-22.467,-17.876,-22.646,-47.060,-32.147,-5.633,-14.491,-8.843,-1.281,-3.674,-7.760,-20.313,-15.561,-3.767,-5.337,-17.474,-16.075,-10.422,-7.417,-12.333,-10.623,-2.812,-2.958,-15.514,-21.062,-10.440,-10.448,-19.008,-34.200,-22.741,-1.042,-5.458,-5.778,-6.170,-2.907,-0.737,-7.075,-3.767,-15.561,-11.658,-28.921,-32.061,-25.402,-17.291,-15.521,-8.556,-9.602,-13.088,-4.295,-32.094,-22.388,-13.314,-11.184,-34.301,-21.832,-18.227,-32.592,-27.112,-11.085,-3.817,-15.948,-26.722,-5.337,-11.658,-15.561,-3.898,-6.270,-5.573,-0.837,-2.326,-5.129,-0.547,-3.830,-29.393,-46.606,-38.132,-23.992,-15.561,-42.021,-29.598,-38.981,-59.021,-50.933,-25.670,-14.577,-34.891,-47.318,-17.474,-28.921,-3.898,-15.561,-2.802,-5.848,-2.293,-4.558,-12.745,-6.267,-10.108,-54.020,-65.053,-52.074,-36.746,-29.072,-63.480,-47.020,-40.223,-60.772,-47.903,-20.974,-19.552,-40.428,-58.700,-16.075,-32.061,-6.270,-2.802,-15.561,-1.153,-2.547,-11.918,-21.792,-6.608,-6.127,-59.772,-64.292,-49.382,-36.043,-31.518,-66.777,-49.076,-31.372,-49.446,-36.342,-13.189,-16.461,-33.715,-52.987,-10.422,-25.402,-5.573,-5.848,-1.153,-15.561,-2.341,-13.500,-21.388,-4.453,-2.409,-50.490,-42.983,-31.872,-20.618,-16.679,-44.127,-30.043,-24.265,-40.740,-32.455,-12.655,-8.008,-23.112,-36.801,-7.417,-17.291,-0.837,-2.293,-2.547,-2.341,-15.561,-4.618,-9.909,-1.176,-2.902,-38.551,-22.071,-16.521,-7.652,-3.816,-20.513,-11.593,-24.331,-39.186,-37.470,-21.379,-4.993,-17.983,-23.903,-12.333,-15.521,-2.326,-4.558,-11.918,-13.500,-4.618,-15.561,-2.213,-4.956,-12.123,-32.179,-11.619,-6.795,-1.947,-2.318,-13.814,-6.297,-15.561,-26.203,-28.040,-19.102,-1.995,-8.972,-11.639,-10.623,-8.556,-5.129,-12.745,-21.792,-21.388,-9.909,-2.213,-15.561,-7.185,-15.493,-19.053,-35.921,-24.393,-15.996,-15.480,-40.627,-26.753,-14.764,-28.088,-21.515,-6.893,-3.791,-14.360,-26.800,-2.812,-9.602,-0.547,-6.267,-6.608,-4.453,-1.176,-4.956,-7.185,-15.561,-1.669,-26.659,-51.256,-36.228,-27.065,-27.314,-58.475,-41.495,-16.689,-30.241,-20.068,-4.432,-8.877,-19.640,-36.526,-2.958,-13.088,-3.830,-10.108,-6.127,-2.409,-2.902,-12.123,-15.493,-1.669,-15.561,-31.941,-19.754,-10.369,-15.957,-30.247,-37.880,-28.560,-2.962,-1.572,-7.021,-17.315,-12.546,-2.097,-3.661,-15.514,-4.295,-29.393,-54.020,-59.772,-50.490,-38.551,-32.179,-19.053,-26.659,-31.941,-15.561};
	//unsigned int i[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24};
	//unsigned int j[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
	//int m=25, N=625;	// 25 data points, 625 similarites & preferences 

	double* APsij = new double[groupSize*groupSize];
	double* APsijd = new double[groupSize*groupSize-groupSize];
	unsigned int* APi = new unsigned int[groupSize*groupSize];
	unsigned int* APj = new unsigned int[groupSize*groupSize];

	for (int i = 0; i < groupSize; i++)
	{
		for (int j = 0; j < groupSize; j++)
		{
			APi[i*groupSize+j] = i;
			APj[i*groupSize+j] = j;
		}
	}

	///////////////////////////////////////
	// run affinity propagation on the current pairwise distance matrix

	//  // use APCluster to do clustering and save results
	int APindexij = 0;
	for(int i = 0; i < groupSize; i++)
	{
		for (int j = 0;j < groupSize;j++) 
		{
			APsij[i*groupSize+j] = 0-distanceMatrix[i][j];
			if (i!=j)
			{
				APsijd[APindexij] = 0-distanceMatrix[i][j];
				APindexij++;
			}
		}
	}

	//   // find median of APsijd to set all APsij[i][i]
	int* APSortindexij = new int[groupSize*groupSize - groupSize];
	for (int i = 0; i < groupSize*groupSize - groupSize; i++) 
		APSortindexij[i]=i;
	bubbleSort(APsijd, APSortindexij, groupSize*groupSize - groupSize);
	double APmedian = (APsijd[(groupSize*groupSize - groupSize)/2]+APsijd[(groupSize*groupSize - groupSize)/2-1])/2;
	for (int i = 0; i < groupSize; i++) 
		APsij[i*groupSize+i]=APmedian;
	delete[] APSortindexij;

	//   // run AP
	int  (*apcluster32)(double*,unsigned int*, unsigned int*, unsigned int, int*, double*, APOPTIONS*); // function pointer 
	//typedef int (__stdcall *APcluster32)(double*, unsigned int*, unsigned int*, unsigned int, int*, double*, APOPTIONS*); 
	
	APOPTIONS AP_options={0};
	void *APdlh =  NULL;//HMODULE APdlh = NULL;
	char *APerror;
	int APiter, *APidx = 0; 
	double APnetsim = 0.0; 

	//
	//if (!(APdlh=LoadLibrary("apclusterwin.dll"))){printf("LoadLibrary() failed: %d\n", GetLastError());return;}
	//APcluster32 apcluster32_here = NULL;
	//apcluster32_here = (APcluster32)GetProcAddress(APdlh, "apcluster32");
	//if (!apcluster32_here) {printf("GetProcAddress() failed: %d\n", GetLastError()); return;}
	///

	if (!(APdlh=dlopen("./apclusterunix64.so", RTLD_LAZY))) {printf("dlh open error 64%s\n",dlerror()); return;}
	apcluster32 = (int (*)(double*,unsigned int*, unsigned int*, unsigned int, int*, double*, APOPTIONS*))dlsym(APdlh, "apcluster32");
	if((APerror = dlerror())!=NULL) { printf("%s\n",APerror); return;}

	AP_options.cbSize = sizeof(APOPTIONS);	AP_options.lambda = 0.9;	AP_options.minimum_iterations = 1;
	AP_options.converge_iterations = 200;	AP_options.nonoise = 0;		AP_options.maximum_iterations = 2000;
	AP_options.progressf=NULL;				AP_options.progress=NULL; //callback;	

	//APiter = (*apcluster32_here)(APsij, APi, APj, groupSize*groupSize, APidx=(int*)calloc(groupSize,sizeof(*APidx)), &APnetsim, &AP_options);
	APiter = (*apcluster32)(APsij, APi, APj, groupSize*groupSize, APidx=(int*)calloc(groupSize,sizeof(*APidx)), &APnetsim, &AP_options); // actual function call 
	//APiter = (*apcluster32)(sij, i, j, N, APidx=(int*)calloc(m,sizeof(*APidx)), &APnetsim, &AP_options); // actual function call 
	//std::cout << "err1";
	if(APiter>0){}	else {printf("Error code: %d\n", APiter);return;}
	//std::cout << "err2";
	//std::cout << APiter << std::endl; 
	//for(int APk = 0; APk < m; APk++) printf("%d\n", APidx[APk]+1);
	//for(int APk = 0; APk < groupSize; APk++) printf("%d\n", APidx[APk]+1);
	dlclose(APdlh);// FreeLibrary(APdlh); 

	//////
	//std::cout << "groupsize: " << groupSize << std::endl;
	//for (int i = 0; i < groupSize*groupSize; i++)
	//{
	//	std::cout << "(" << APi[i] << "," << APj[i] << ") " << APsij[i] << std::endl;
	//}

	//  // analyze AP output APidx to get the number of classes
	int* APindex = new int[groupSize];
	double* APoutput = new double[groupSize];

	for (int i = 0; i < groupSize; i++)
	{
		APindex[i]=i;
		APoutput[i]= (double)APidx[i];
	}
	bubbleSort(APoutput, APindex, groupSize);

	int APclassnum = 1;
	int APclassindexcur = APoutput[0];
	for (int i = 1; i<groupSize; i++) 
	{
		if (APclassindexcur != APoutput[i])
		{
			APclassnum++;
			APclassindexcur = APoutput[i];
		}
	}

	//  // get exemplars for each class, count class size and get the first (smallest) index in each class
	int* APclassrep = new int[APclassnum];
	for (int i = 0; i < APclassnum; i++) APclassrep[i] = 0;
	int* APclasssize = new int[APclassnum];
	for (int i = 0; i < APclassnum; i++) APclasssize[i] = 0;
	int* APclassfirst = new int[APclassnum];
	for (int i = 0; i < APclassnum; i++) APclassfirst[i] = 0;
	int* APclassstart = new int[APclassnum];
	for (int i = 0; i < APclassnum; i++) APclassstart[i] = 0;

	int APclassindex = 0;
	APclassstart[APclassindex] = 0;
	APclasssize[APclassindex] = 1;
	APclassrep[APclassindex] = (int)APoutput[0];
	APclassfirst[APclassindex] = APindex[0];
	for (int i = 1; i < groupSize; i++) 
	{
		if (APclassrep[APclassindex] != (int)APoutput[i])
		{
			APclassindex++;
			APclassstart[APclassindex] = i;
			APclasssize[APclassindex]++;
			APclassrep[APclassindex] = (int)APoutput[i];
			if (APclasssize[APclassindex] == 1)
				APclassfirst[APclassindex] = APindex[i];
			else if (APindex[i]<APclassfirst[APclassindex])
				APclassfirst[APclassindex] = APindex[i];
		}
		else
		{
			APclasssize[APclassindex]++;
			if (APindex[i]<APclassfirst[APclassindex])
				APclassfirst[APclassindex] = APindex[i];
		}
	}

	//  // test output
	//for (int i = 0; i < APclassnum; i++) printf("%d ", APclassrep[i]); printf("\n");
	//for (int i = 0; i < APclassnum; i++) printf("%d ", APclassfirst[i]); printf("\n");
	//for (int i = 0; i < APclassnum; i++) printf("%d ", APclasssize[i]); printf("\n");
	//for (int i = 0; i < APclassnum; i++) printf("%d ", APclassstart[i]); printf("\n");

	//  // sort all classes based on APclassfirst;
	int* APclassfirstindex = new int[APclassnum];
	for (int i = 0;i < APclassnum;i++) APclassfirstindex[i] = i;
	double* APclassfirstdouble = new double[APclassnum];
	for (int i = 0;i < APclassnum;i++) APclassfirstdouble[i] = (double)APclassfirst[i];
	bubbleSort(APclassfirstdouble, APclassfirstindex, APclassnum);

	//  // write into the data structure of final AP results
	APclustersize = APclassnum;
	APcluster = new int*[APclustersize];
	APclustersubsize = new int[APclustersize];
	APclustersubrep = new int[APclustersize];

	for (int i = 0; i < APclustersize; i++)
	{
		APcluster[i] = new int[APclasssize[APclassfirstindex[i]]];
		APclustersubsize[i] = APclasssize[APclassfirstindex[i]];
		APclustersubrep[i] = APclassrep[APclassfirstindex[i]];
	}
	for (int i = 0;i < APclustersize;i++)
		for (int j = 0; j < APclustersubsize[i];j++)
			APcluster[i][j] = (int)APindex[APclassstart[APclassfirstindex[i]]+j];

	delete[] APclassfirstindex;
	delete[] APclassfirstdouble;

	delete[] APclassstart;
	delete[] APclasssize;
	delete[] APclassfirst;
	delete[] APclassrep;

	delete[] APindex;
	delete[] APoutput;

	//  // delete APidx
	if(APidx) free(APidx);

	// run affinity propagation on the current pairwise distance matrix
	///////////////////////////////////////
	// delete newed variables
	delete[] APsij;
	delete[] APsijd;
	delete[] APi;
	delete[] APj;

	return;
}
// do Affinity Propagation Clustering on Unix64 system with specified lib/so file name

void APClusteringUnix64WithLibName(char* APLibName, double** distanceMatrix, int groupSize, int &APclustersize, 
						int* &APclustersubrep, int* &APclustersubsize, int** &APcluster)
{
	//
	//std::cout << "groupsize: " << groupSize << std::endl;
	// for APcluster
	//double sij[]={-15.561,-1.861,-4.066,-9.293,-3.581,-3.093,-27.967,-32.147,-44.921,-47.967,-17.665,-15.628,-6.430,-36.111,-21.062,-32.094,-46.606,-65.053,-64.292,-42.983,-22.071,-11.619,-35.921,-51.256,-19.754,	-1.861,-15.561,-1.839,-8.387,-8.612,-4.909,-15.438,-19.376,-28.647,-31.421,-9.245,-6.721,-1.980,-22.467,-10.440,-22.388,-38.132,-52.074,-49.382,-31.872,-16.521,-6.795,-24.393,-36.228,-10.369,-4.066,-1.839,-15.561,-2.561,-7.097,-2.363,-17.266,-25.109,-31.390,-27.590,-5.415,-8.475,-6.370,-17.876,-10.448,-13.314,-23.992,-36.746,-36.043,-20.618,-7.652,-1.947,-15.996,-27.065,-15.957,-9.293,-8.387,-2.561,-15.561,-6.642,-2.241,-28.773,-41.141,-45.675,-34.572,-8.612,-18.211,-17.008,-22.646,-19.008,-11.184,-15.561,-29.072,-31.518,-16.679,-3.816,-2.318,-15.480,-27.314,-30.247,-3.581,-8.612,-7.097,-6.642,-15.561,-1.458,-44.952,-53.686,-66.515,-62.544,-24.441,-29.052,-18.491,-47.060,-34.200,-34.301,-42.021,-63.480,-66.777,-44.127,-20.513,-13.814,-40.627,-58.475,-37.880,-3.093,-4.909,-2.363,-2.241,-1.458,-15.561,-32.239,-41.536,-50.820,-45.330,-14.087,-19.384,-13.103,-32.147,-22.741,-21.832,-29.598,-47.020,-49.076,-30.043,-11.593,-6.297,-26.753,-41.495,-28.560,-27.967,-15.438,-17.266,-28.773,-44.952,-32.239,-15.561,-2.126,-2.110,-5.972,-7.313,-1.788,-9.199,-5.633,-1.042,-18.227,-38.981,-40.223,-31.372,-24.265,-24.331,-15.561,-14.764,-16.689,-2.962,-32.147,-19.376,-25.109,-41.141,-53.686,-41.536,-2.126,-15.561,-2.609,-13.515,-16.243,-4.613,-9.944,-14.491,-5.458,-32.592,-59.021,-60.772,-49.446,-40.740,-39.186,-26.203,-28.088,-30.241,-1.572,-44.921,-28.647,-31.390,-45.675,-66.515,-50.820,-2.110,-2.609,-15.561,-5.895,-15.692,-7.693,-18.816,-8.843,-5.778,-27.112,-50.933,-47.903,-36.342,-32.455,-37.470,-28.040,-21.515,-20.068,-7.021,-47.967,-31.421,-27.590,-34.572,-62.544,-45.330,-5.972,-13.515,-5.895,-15.561,-8.956,-11.047,-26.444,-1.281,-6.170,-11.085,-25.670,-20.974,-13.189,-12.655,-21.379,-19.102,-6.893,-4.432,-17.315,-17.665,-9.245,-5.415,-8.612,-24.441,-14.087,-7.313,-16.243,-15.692,-8.956,-15.561,-4.386,-10.533,-3.674,-2.907,-3.817,-14.577,-19.552,-16.461,-8.008,-4.993,-1.995,-3.791,-8.877,-12.546,-15.628,-6.721,-8.475,-18.211,-29.052,-19.384,-1.788,-4.613,-7.693,-11.047,-4.386,-15.561,-3.332,-7.760,-0.737,-15.948,-34.891,-40.428,-33.715,-23.112,-17.983,-8.972,-14.360,-19.640,-2.097,-6.430,-1.980,-6.370,-17.008,-18.491,-13.103,-9.199,-9.944,-18.816,-26.444,-10.533,-3.332,-15.561,-20.313,-7.075,-26.722,-47.318,-58.700,-52.987,-36.801,-23.903,-11.639,-26.800,-36.526,-3.661,-36.111,-22.467,-17.876,-22.646,-47.060,-32.147,-5.633,-14.491,-8.843,-1.281,-3.674,-7.760,-20.313,-15.561,-3.767,-5.337,-17.474,-16.075,-10.422,-7.417,-12.333,-10.623,-2.812,-2.958,-15.514,-21.062,-10.440,-10.448,-19.008,-34.200,-22.741,-1.042,-5.458,-5.778,-6.170,-2.907,-0.737,-7.075,-3.767,-15.561,-11.658,-28.921,-32.061,-25.402,-17.291,-15.521,-8.556,-9.602,-13.088,-4.295,-32.094,-22.388,-13.314,-11.184,-34.301,-21.832,-18.227,-32.592,-27.112,-11.085,-3.817,-15.948,-26.722,-5.337,-11.658,-15.561,-3.898,-6.270,-5.573,-0.837,-2.326,-5.129,-0.547,-3.830,-29.393,-46.606,-38.132,-23.992,-15.561,-42.021,-29.598,-38.981,-59.021,-50.933,-25.670,-14.577,-34.891,-47.318,-17.474,-28.921,-3.898,-15.561,-2.802,-5.848,-2.293,-4.558,-12.745,-6.267,-10.108,-54.020,-65.053,-52.074,-36.746,-29.072,-63.480,-47.020,-40.223,-60.772,-47.903,-20.974,-19.552,-40.428,-58.700,-16.075,-32.061,-6.270,-2.802,-15.561,-1.153,-2.547,-11.918,-21.792,-6.608,-6.127,-59.772,-64.292,-49.382,-36.043,-31.518,-66.777,-49.076,-31.372,-49.446,-36.342,-13.189,-16.461,-33.715,-52.987,-10.422,-25.402,-5.573,-5.848,-1.153,-15.561,-2.341,-13.500,-21.388,-4.453,-2.409,-50.490,-42.983,-31.872,-20.618,-16.679,-44.127,-30.043,-24.265,-40.740,-32.455,-12.655,-8.008,-23.112,-36.801,-7.417,-17.291,-0.837,-2.293,-2.547,-2.341,-15.561,-4.618,-9.909,-1.176,-2.902,-38.551,-22.071,-16.521,-7.652,-3.816,-20.513,-11.593,-24.331,-39.186,-37.470,-21.379,-4.993,-17.983,-23.903,-12.333,-15.521,-2.326,-4.558,-11.918,-13.500,-4.618,-15.561,-2.213,-4.956,-12.123,-32.179,-11.619,-6.795,-1.947,-2.318,-13.814,-6.297,-15.561,-26.203,-28.040,-19.102,-1.995,-8.972,-11.639,-10.623,-8.556,-5.129,-12.745,-21.792,-21.388,-9.909,-2.213,-15.561,-7.185,-15.493,-19.053,-35.921,-24.393,-15.996,-15.480,-40.627,-26.753,-14.764,-28.088,-21.515,-6.893,-3.791,-14.360,-26.800,-2.812,-9.602,-0.547,-6.267,-6.608,-4.453,-1.176,-4.956,-7.185,-15.561,-1.669,-26.659,-51.256,-36.228,-27.065,-27.314,-58.475,-41.495,-16.689,-30.241,-20.068,-4.432,-8.877,-19.640,-36.526,-2.958,-13.088,-3.830,-10.108,-6.127,-2.409,-2.902,-12.123,-15.493,-1.669,-15.561,-31.941,-19.754,-10.369,-15.957,-30.247,-37.880,-28.560,-2.962,-1.572,-7.021,-17.315,-12.546,-2.097,-3.661,-15.514,-4.295,-29.393,-54.020,-59.772,-50.490,-38.551,-32.179,-19.053,-26.659,-31.941,-15.561};
	//unsigned int i[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24};
	//unsigned int j[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
	//int m=25, N=625;	// 25 data points, 625 similarites & preferences 

	double* APsij = new double[groupSize*groupSize];
	double* APsijd = new double[groupSize*groupSize-groupSize];
	unsigned int* APi = new unsigned int[groupSize*groupSize];
	unsigned int* APj = new unsigned int[groupSize*groupSize];

	for (int i = 0; i < groupSize; i++)
	{
		for (int j = 0; j < groupSize; j++)
		{
			APi[i*groupSize+j] = i;
			APj[i*groupSize+j] = j;
		}
	}

	///////////////////////////////////////
	// run affinity propagation on the current pairwise distance matrix

	//  // use APCluster to do clustering and save results
	int APindexij = 0;
	for(int i = 0; i < groupSize; i++)
	{
		for (int j = 0;j < groupSize;j++) 
		{
			APsij[i*groupSize+j] = 0-distanceMatrix[i][j];
			if (i!=j)
			{
				APsijd[APindexij] = 0-distanceMatrix[i][j];
				APindexij++;
			}
		}
	}

	//   // find median of APsijd to set all APsij[i][i]
	int* APSortindexij = new int[groupSize*groupSize - groupSize];
	for (int i = 0; i < groupSize*groupSize - groupSize; i++) 
		APSortindexij[i]=i;
	bubbleSort(APsijd, APSortindexij, groupSize*groupSize - groupSize);
	double APmedian = (APsijd[(groupSize*groupSize - groupSize)/2]+APsijd[(groupSize*groupSize - groupSize)/2-1])/2;
	for (int i = 0; i < groupSize; i++) 
		APsij[i*groupSize+i]=APmedian;
	delete[] APSortindexij;

	//   // run AP
	int  (*apcluster32)(double*,unsigned int*, unsigned int*, unsigned int, int*, double*, APOPTIONS*); // function pointer 
	//typedef int (__stdcall *APcluster32)(double*, unsigned int*, unsigned int*, unsigned int, int*, double*, APOPTIONS*); 

	APOPTIONS AP_options={0};
	void *APdlh =  NULL;//HMODULE APdlh = NULL;
	char *APerror;
	int APiter, *APidx = 0; 
	double APnetsim = 0.0; 

	//
	//if (!(APdlh=LoadLibrary("apclusterwin.dll"))){printf("LoadLibrary() failed: %d\n", GetLastError());return;}
	//APcluster32 apcluster32_here = NULL;
	//apcluster32_here = (APcluster32)GetProcAddress(APdlh, "apcluster32");
	//if (!apcluster32_here) {printf("GetProcAddress() failed: %d\n", GetLastError()); return;}
	///
	// "./apclusterunix64.so"
	if (!(APdlh=dlopen(APLibName, RTLD_LAZY))) {printf("dlh open error 64%s\n",dlerror()); return;}
	apcluster32 = (int (*)(double*,unsigned int*, unsigned int*, unsigned int, int*, double*, APOPTIONS*))dlsym(APdlh, "apcluster32");
	if((APerror = dlerror())!=NULL) { printf("%s\n",APerror); return;}

	AP_options.cbSize = sizeof(APOPTIONS);	AP_options.lambda = 0.9;	AP_options.minimum_iterations = 1;
	AP_options.converge_iterations = 200;	AP_options.nonoise = 0;		AP_options.maximum_iterations = 2000;
	AP_options.progressf=NULL;				AP_options.progress=NULL; //callback;	

	//APiter = (*apcluster32_here)(APsij, APi, APj, groupSize*groupSize, APidx=(int*)calloc(groupSize,sizeof(*APidx)), &APnetsim, &AP_options);
	APiter = (*apcluster32)(APsij, APi, APj, groupSize*groupSize, APidx=(int*)calloc(groupSize,sizeof(*APidx)), &APnetsim, &AP_options); // actual function call 
	//APiter = (*apcluster32)(sij, i, j, N, APidx=(int*)calloc(m,sizeof(*APidx)), &APnetsim, &AP_options); // actual function call 
	//std::cout << "err1";
	if(APiter>0){}	else {printf("Error code: %d\n", APiter);return;}
	//std::cout << "err2";
	//std::cout << APiter << std::endl; 
	//for(int APk = 0; APk < m; APk++) printf("%d\n", APidx[APk]+1);
	//for(int APk = 0; APk < groupSize; APk++) printf("%d\n", APidx[APk]+1);
	dlclose(APdlh);// FreeLibrary(APdlh); 

	//////
	//std::cout << "groupsize: " << groupSize << std::endl;
	//for (int i = 0; i < groupSize*groupSize; i++)
	//{
	//	std::cout << "(" << APi[i] << "," << APj[i] << ") " << APsij[i] << std::endl;
	//}

	//  // analyze AP output APidx to get the number of classes
	int* APindex = new int[groupSize];
	double* APoutput = new double[groupSize];

	for (int i = 0; i < groupSize; i++)
	{
		APindex[i]=i;
		APoutput[i]= (double)APidx[i];
	}
	bubbleSort(APoutput, APindex, groupSize);

	int APclassnum = 1;
	int APclassindexcur = APoutput[0];
	for (int i = 1; i<groupSize; i++) 
	{
		if (APclassindexcur != APoutput[i])
		{
			APclassnum++;
			APclassindexcur = APoutput[i];
		}
	}

	//  // get exemplars for each class, count class size and get the first (smallest) index in each class
	int* APclassrep = new int[APclassnum];
	for (int i = 0; i < APclassnum; i++) APclassrep[i] = 0;
	int* APclasssize = new int[APclassnum];
	for (int i = 0; i < APclassnum; i++) APclasssize[i] = 0;
	int* APclassfirst = new int[APclassnum];
	for (int i = 0; i < APclassnum; i++) APclassfirst[i] = 0;
	int* APclassstart = new int[APclassnum];
	for (int i = 0; i < APclassnum; i++) APclassstart[i] = 0;

	int APclassindex = 0;
	APclassstart[APclassindex] = 0;
	APclasssize[APclassindex] = 1;
	APclassrep[APclassindex] = (int)APoutput[0];
	APclassfirst[APclassindex] = APindex[0];
	for (int i = 1; i < groupSize; i++) 
	{
		if (APclassrep[APclassindex] != (int)APoutput[i])
		{
			APclassindex++;
			APclassstart[APclassindex] = i;
			APclasssize[APclassindex]++;
			APclassrep[APclassindex] = (int)APoutput[i];
			if (APclasssize[APclassindex] == 1)
				APclassfirst[APclassindex] = APindex[i];
			else if (APindex[i]<APclassfirst[APclassindex])
				APclassfirst[APclassindex] = APindex[i];
		}
		else
		{
			APclasssize[APclassindex]++;
			if (APindex[i]<APclassfirst[APclassindex])
				APclassfirst[APclassindex] = APindex[i];
		}
	}

	//  // test output
	//for (int i = 0; i < APclassnum; i++) printf("%d ", APclassrep[i]); printf("\n");
	//for (int i = 0; i < APclassnum; i++) printf("%d ", APclassfirst[i]); printf("\n");
	//for (int i = 0; i < APclassnum; i++) printf("%d ", APclasssize[i]); printf("\n");
	//for (int i = 0; i < APclassnum; i++) printf("%d ", APclassstart[i]); printf("\n");

	//  // sort all classes based on APclassfirst;
	int* APclassfirstindex = new int[APclassnum];
	for (int i = 0;i < APclassnum;i++) APclassfirstindex[i] = i;
	double* APclassfirstdouble = new double[APclassnum];
	for (int i = 0;i < APclassnum;i++) APclassfirstdouble[i] = (double)APclassfirst[i];
	bubbleSort(APclassfirstdouble, APclassfirstindex, APclassnum);

	//  // write into the data structure of final AP results
	APclustersize = APclassnum;
	APcluster = new int*[APclustersize];
	APclustersubsize = new int[APclustersize];
	APclustersubrep = new int[APclustersize];

	for (int i = 0; i < APclustersize; i++)
	{
		APcluster[i] = new int[APclasssize[APclassfirstindex[i]]];
		APclustersubsize[i] = APclasssize[APclassfirstindex[i]];
		APclustersubrep[i] = APclassrep[APclassfirstindex[i]];
	}
	for (int i = 0;i < APclustersize;i++)
		for (int j = 0; j < APclustersubsize[i];j++)
			APcluster[i][j] = (int)APindex[APclassstart[APclassfirstindex[i]]+j];

	delete[] APclassfirstindex;
	delete[] APclassfirstdouble;

	delete[] APclassstart;
	delete[] APclasssize;
	delete[] APclassfirst;
	delete[] APclassrep;

	delete[] APindex;
	delete[] APoutput;

	//  // delete APidx
	if(APidx) free(APidx);

	// run affinity propagation on the current pairwise distance matrix
	///////////////////////////////////////
	// delete newed variables
	delete[] APsij;
	delete[] APsijd;
	delete[] APi;
	delete[] APj;

	return;
}

*/
void MakeFileName(char* filename, char* cur_id, char* padding, int cur_index, char* filetype)
{
	// create on 20091228
	// to generate a file name like "{cur_id}{padding}{cur_index}.{type}"

	char index_string[5]; 

	if (cur_index >= 1000)
	{
		return;
	}
	else if (cur_index < 0)
	{
		index_string[0] = '\0';
	}
	else
		myitoa(cur_index, index_string, 3);

	int totallength = 0;
	totallength = strlen(cur_id) + strlen(padding) + strlen(index_string) 
		+ strlen(".") + strlen(filetype) + 1; // 1 for '\0'

	if (totallength > MAX_FILE_NAME_LENGTH)
	{
		return;
	}

	strcpy(filename, cur_id);
	strcat(filename, padding);
	strcat(filename, index_string);
	strcat(filename, ".");
	strcat(filename, filetype);
	return;
}

void MakeFileName(char* filename, char* cur_id, char* padding, char* index_string, char* filetype)
{
	// create on 20091228
	// to generate a file name like "{cur_id}{padding}{index_string}.{type}"
	int totallength = 0;
	totallength = strlen(cur_id) + strlen(padding) + strlen(index_string) 
		+ strlen(".") + strlen(filetype) + 1; // 1 for '\0'

	if (totallength > MAX_FILE_NAME_LENGTH)
	{
		return;
	}
	strcpy(filename, cur_id);
	strcat(filename, padding);
	strcat(filename, index_string);
	strcat(filename, ".");
	strcat(filename, filetype);
	return;
}
void MakeFileName(char* filename, char* id1, char* padding1, char* id2, char* padding2, int cur_index, char* filetype)
{
	// create on 20091228
	// to generate a file name like "{id1}{padding1}{id2}{padding2}{cur_index}.{type}"

	char index_string[5]; 

	if (cur_index >= 1000)
	{
		return;
	}
	else if (cur_index < 0)
	{
		index_string[0] = '\0';
	}
	else
		myitoa(cur_index, index_string, 3);

	int totallength = 0;
	totallength = strlen(id1) + strlen(padding1) + strlen(id2) + strlen(padding2) + strlen(index_string) 
		+ strlen(".") + strlen(filetype) + 1; // 1 for '\0'

	if (totallength > MAX_FILE_NAME_LENGTH)
	{
		return;
	}

	strcpy(filename, id1);
	strcat(filename, padding1);
	strcat(filename, id2);
	strcat(filename, padding2);
	strcat(filename, index_string);
	strcat(filename, ".");
	strcat(filename, filetype);
	return;
}

void MakeFileName(char* filename, char* id1, char* padding1, char* id2, char* padding2, char* index_string, char* filetype)
{
	// create on 20091228
	// to generate a file name like "{id1}{padding1}{id2}{padding2}{index_string}.{type}"

	int totallength = 0;
	totallength = strlen(id1) + strlen(padding1) + strlen(id2) + strlen(padding2) + strlen(index_string) 
		+ strlen(".") + strlen(filetype) + 1; // 1 for '\0'

	if (totallength > MAX_FILE_NAME_LENGTH)
	{
		return;
	}

	strcpy(filename, id1);
	strcat(filename, padding1);
	strcat(filename, id2);
	strcat(filename, padding2);
	strcat(filename, index_string);
	strcat(filename, ".");
	strcat(filename, filetype);
	return;
}

void MyPause()
{
	// create on 20091228
	// wait for any user input to continue

	//if (isDebugging)
	//{	
	std::cout << "Press any key to continue ...";
	getchar();
	std::cout << std::endl;
	//}

}


// do affine registration 20100524

void AffineRegistration(char* fixedImageFileName, char* movingImageFileName, char* transformedFileName, char* affineTransformFileName)
{
	// read images
	InternalImageType::Pointer fixedImage = 0;
	InternalImageType::Pointer movingImage = 0;

	ReadImage(fixedImageFileName, fixedImage);
	ReadImage(movingImageFileName, movingImage);

	// bool doAffine = false;
	AffineTransformType::Pointer finalTransform = AffineTransformType::New();
	if (!strcmp(fixedImageFileName, movingImageFileName))
	{
		// if fixed and moving are the same image
		finalTransform->SetIdentity();
		TransformFileWriterType::Pointer  transformFileWriter = TransformFileWriterType::New();
		transformFileWriter->SetFileName(affineTransformFileName);
		transformFileWriter->SetPrecision(12);
		transformFileWriter->SetInput(finalTransform);

		try
		{
			transformFileWriter->Update();
		}
		catch(itk::ExceptionObject & err)
		{
			std::cerr << err << std::endl;
			return;
		}

	}
	else
	{
		// affine
		InternalImageRegistrationType::Pointer				registration	= InternalImageRegistrationType::New();
		InternalMeanSquaresMetricType::Pointer				metric			= InternalMeanSquaresMetricType::New();

		RegularStepGradientDescentOptimizerType::Pointer	optimizer		= RegularStepGradientDescentOptimizerType::New();
		// ConjugateGradientOptimizerType::Pointer				optimizer		= ConjugateGradientOptimizerType::New();
		// GradientDescentOptimizerType::Pointer				optimizer		= GradientDescentOptimizerType::New();
		// LBFGSOptimizerType::Pointer								optimizer		= LBFGSOptimizerType::New();

		InternalLinearInterpolatorType::Pointer				interpolator	= InternalLinearInterpolatorType::New();
		AffineTransformType::Pointer						affineTransform = AffineTransformType::New();


		registration->SetMetric(        metric			);
		registration->SetOptimizer(     optimizer		);
		registration->SetInterpolator(  interpolator	);
		registration->SetTransform(		affineTransform	);

		registration->SetFixedImage(    fixedImage      );
		registration->SetMovingImage(   movingImage     );
		registration->SetFixedImageRegion(	fixedImage->GetBufferedRegion() );


		// initialize affine transformation
		AffineTransformInitializerType::Pointer initializer = AffineTransformInitializerType::New();
		initializer->SetTransform(   affineTransform );
		initializer->SetFixedImage(  fixedImage );
		initializer->SetMovingImage( movingImage );
		initializer->MomentsOn();
		initializer->InitializeTransform();

		registration->SetInitialTransformParameters( affineTransform->GetParameters() );

		double translationScale = 1.0 / 1000.0;
		RegularStepGradientDescentOptimizerType::ScalesType  optimizerScales( affineTransform->GetNumberOfParameters() );
		// ConjugateGradientOptimizerType::ScalesType  optimizerScales( affineTransform->GetNumberOfParameters() );
		// GradientDescentOptimizerType::ScalesType  optimizerScales( affineTransform->GetNumberOfParameters() );
		// LBFGSOptimizerType::ScalesType  optimizerScales( affineTransform->GetNumberOfParameters() );

		optimizerScales[0] =  1.0;	optimizerScales[1] =  1.0;	optimizerScales[2] =  1.0;
		optimizerScales[3] =  1.0;	optimizerScales[4] =  1.0;	optimizerScales[5] =  1.0;
		optimizerScales[6] =  1.0;	optimizerScales[7] =  1.0;	optimizerScales[8] =  1.0;
		optimizerScales[9]  =  translationScale;
		optimizerScales[10] =  translationScale;
		optimizerScales[11] =  translationScale;

		optimizer->SetScales( optimizerScales );

		double steplength = 4;

		unsigned int maxNumberOfIterations = 300;

		// for RegularStepGradientDescentOptimizerType
		optimizer->SetMaximumStepLength( steplength ); 
		optimizer->SetMinimumStepLength( 0.01 );
		optimizer->SetNumberOfIterations( maxNumberOfIterations );
		optimizer->MinimizeOn();

		// // for GradientDescentOptimizerType
		//optimizer->SetLearningRate(15.0);
		//optimizer->SetNumberOfIterations(200);
		//optimizer->MinimizeOn();

		//// for LBFGSOptimizerType
		//optimizer->SetGradientConvergenceTolerance( 0.05 );
		//optimizer->SetLineSearchAccuracy( 0.9 );
		//optimizer->SetDefaultStepLength( 1.5 );
		//optimizer->TraceOn();
		//optimizer->SetMaximumNumberOfFunctionEvaluations( 1000 );

		// // add observer

		//CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
		//optimizer->AddObserver( itk::IterationEvent(), observer );


		try 
		{ 
			registration->StartRegistration(); 
		} 
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << "ExceptionObject caught !" << std::endl; 
			std::cerr << err << std::endl; 
			return ;
		} 

		RegularStepGradientDescentOptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

		finalTransform->SetParameters( finalParameters );
		finalTransform->SetFixedParameters( affineTransform->GetFixedParameters() );

		// save tranform matrix
		// TransformFileReaderType::Pointer affineTransformWriter = TransformFileReaderType::New();
		TransformFileWriterType::Pointer  transformFileWriter = TransformFileWriterType::New();
		transformFileWriter->SetFileName(affineTransformFileName);
		transformFileWriter->SetPrecision(12);
		transformFileWriter->SetInput(finalTransform);

		try
		{
			transformFileWriter->Update();
		}
		catch(itk::ExceptionObject & err)
		{
			std::cerr << err << std::endl;
			return;
		}
		// std::cerr << "affine done!" << std::endl;
	}

	{
		// resample moving image to get warped image
		ResampleFilterType::Pointer resampler = ResampleFilterType::New();

		resampler->SetTransform( finalTransform );
		resampler->SetInput( movingImage );

		resampler->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
		resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
		resampler->SetOutputSpacing( fixedImage->GetSpacing() );
		resampler->SetOutputDirection( fixedImage->GetDirection() );
		resampler->SetDefaultPixelValue( 0 ); //100
		try
		{
			resampler->Update();
		}
		catch(itk::ExceptionObject & err)
		{
			std::cerr << err << std::endl;
			return;
		}

		InternalImageType::Pointer deformedImage = 0;
		deformedImage = resampler->GetOutput();
		deformedImage->DisconnectPipeline();

		WriteImage(transformedFileName, deformedImage);
		//WriteImage("AA_out.hdr", resampler->GetOutput());
	}
}

void ComposeAffineAndDense(char* affineTransformationFileName, char* nonRigidDeformationFieldFileName, char* finalDeformationFieldFileName)
{
	DeformationFieldType::Pointer nonRigidDeformationField = 0;
	ReadDeformationField(nonRigidDeformationFieldFileName, nonRigidDeformationField);

	AffineTransformType::Pointer finalTransform = AffineTransformType::New();
	DeformationFieldType::Pointer deformationField_Affine = 0;

	// read existing affine transform
	{
		TransformFileReaderType::Pointer        transformFileReader = TransformFileReaderType::New();
		transformFileReader->SetFileName(affineTransformationFileName);
		try
		{
			transformFileReader->Update();
		}
		catch(itk::ExceptionObject & err)
		{
			std::cerr << err << std::endl;
			return;
		}

		TransformListType*   transformList = transformFileReader->GetTransformList();
		// AffineTransformType::Pointer finalTransform = AffineTransformType::New();
		finalTransform->SetFixedParameters(transformList->front()->GetFixedParameters());

		AffineTransformType::ParametersType affineParameters;
		affineParameters.set_size(finalTransform->GetNumberOfParameters());
		finalTransform->SetParameters(affineParameters);
		finalTransform->SetParametersByValue(transformList->front()->GetParameters());
	}

	// affine transform to dense deformation field
	{
		// transform to deformation field
		Transform2DeformationFieldGeneratorType::IndexType indexT2D; indexT2D.Fill(0);
		Transform2DeformationFieldGeneratorType::Pointer defGenerator = Transform2DeformationFieldGeneratorType::New();
		defGenerator->SetOutputSize(nonRigidDeformationField->GetLargestPossibleRegion().GetSize());
		defGenerator->SetOutputSpacing (nonRigidDeformationField->GetSpacing());
		defGenerator->SetOutputOrigin (nonRigidDeformationField->GetOrigin());
		defGenerator->SetOutputDirection(nonRigidDeformationField->GetDirection());
		defGenerator->SetOutputIndex(indexT2D);

		defGenerator->SetTransform(finalTransform);

		try
		{
			defGenerator->Update();
		}
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << err << std::endl; 
			return;
		} 

		deformationField_Affine = defGenerator->GetOutput();
		deformationField_Affine->DisconnectPipeline();
	}
	// compose tow deformation field
	DeformationFieldType::Pointer composedDeformationField = 0;
	ComposeDeformationFields(deformationField_Affine, nonRigidDeformationField, composedDeformationField);
	WriteDeformationField(finalDeformationFieldFileName, composedDeformationField);

}
void PaddingImage(char* originalImageFileName, char* paddedImageFileName, int sizex, int sizey, int sizez)
{
	// padding the original image to a larger size
	InternalImageType::Pointer originImage = 0;
	ReadImage(originalImageFileName, originImage);
	InternalImageType::SizeType sizeOrigin = originImage->GetLargestPossibleRegion().GetSize();
	sizex = sizex>(int)sizeOrigin[0]?sizex:(int)sizeOrigin[0];
	sizey = sizey>(int)sizeOrigin[1]?sizey:(int)sizeOrigin[1];
	sizez = sizez>(int)sizeOrigin[2]?sizez:(int)sizeOrigin[2];

	// create the resampled image
	InternalImageType::Pointer newImage = InternalImageType::New();
	InternalImageType::IndexType start;
	start[0] = 0;start[1] = 0;start[2] = 0;
	InternalImageType::SizeType sizeNew;
	sizeNew[0] = sizex;sizeNew[1] = sizey;sizeNew[2] = sizez;
	InternalImageType::RegionType region;
	region.SetSize (sizeNew);
	region.SetIndex (start);

	newImage -> SetRegions (region);
	newImage -> SetSpacing (originImage->GetSpacing());
	newImage -> SetDirection (originImage->GetDirection());
	newImage -> SetOrigin (originImage->GetOrigin());
	newImage -> Allocate();
	
	InternalImageIteratorType itOrigin(originImage, originImage->GetLargestPossibleRegion() );
	InternalImageIteratorType itNew(newImage, newImage->GetLargestPossibleRegion() );
	InternalImageType::IndexType idxNew, idxOrigin;

	for (itNew.GoToBegin(), itOrigin.GoToBegin(); !itNew.IsAtEnd(); ++itNew)
	{
		idxNew = itNew.GetIndex();
		idxOrigin = itOrigin.GetIndex();
		if ( ( ((unsigned int)idxNew[0]>=(sizeNew[0]-sizeOrigin[0])/2) && ((unsigned int)idxNew[0]<(sizeNew[0]-sizeOrigin[0])/2+sizeOrigin[0]) )&& 
			 ( ((unsigned int)idxNew[1]>=(sizeNew[1]-sizeOrigin[1])/2) && ((unsigned int)idxNew[1]<(sizeNew[1]-sizeOrigin[1])/2+sizeOrigin[1]) ) && 
			 ( ((unsigned int)idxNew[2]>=(sizeNew[2]-sizeOrigin[2])/2) && ((unsigned int)idxNew[2]<(sizeNew[2]-sizeOrigin[2])/2+sizeOrigin[2]) ))
		{
			// within central region
			itNew.Set(itOrigin.Get());
			++itOrigin;
		}
		else
		{
			// out of central region
			itNew.Set(0.0);
		}
	}

	//
	WriteImage(paddedImageFileName, newImage);
}
void ResampleImageWithSizeSpacing(char* originalImageFileName, char* resampledImageFileName, int sizex, int sizey, int sizez, float spacingx, float spacingy, float spacingz)
{
	// read original image
	InternalImageType::Pointer originalImage = 0;
	ReadImage(originalImageFileName, originalImage);
	InternalImageType::SizeType sizeOriginal = originalImage->GetLargestPossibleRegion().GetSize();
	InternalImageType::SpacingType spacingOriginal = originalImage->GetSpacing();
//	InternalImageType::PointType originOriginal = originalImage->GetOrigin();
	InternalImageType::DirectionType directionOriginal = originalImage->GetDirection();
	//InternalImageType::IndexType index = originalImage->GetRequestedRegion().GetIndex();
	//InternalImageType::PointType out;
	//originalImage->TransformIndexToPhysicalPoint(index,out);

	// set parameters for new image
	InternalImageType::SpacingType spacingNew;
	spacingNew[0] = spacingx;spacingNew[1] = spacingy;spacingNew[2] = spacingz;
	InternalImageType::SizeType sizeNew;
	sizeNew[0] = sizex; sizeNew[1] = sizey; sizeNew[2] = sizez;
	InternalImageType::PointType originNew;
	originNew[0] = 0.0;originNew[1] = 0.0;originNew[2] = 0.0;
	//originNew[0] = 10.0;originNew[1] = 10.0;originNew[2] = 10.0;
	InternalImageType::DirectionType directionNew;
	directionNew = directionOriginal;

	// setup a filter
	ResampleFilterType::Pointer filter = ResampleFilterType::New();
	
	// set transform matrix, with only translation
	AffineTransformType::Pointer transform = AffineTransformType::New();
	transform->SetIdentity();
	AffineTransformType::OutputVectorType translation;
	translation[0] = 0 - (sizeNew[0]*spacingNew[0] - sizeOriginal[0]*spacingOriginal[0])/2;
	translation[1] = 0 - (sizeNew[1]*spacingNew[1] - sizeOriginal[1]*spacingOriginal[1])/2;
	translation[2] = 0 - (sizeNew[2]*spacingNew[2] - sizeOriginal[2]*spacingOriginal[2])/2;
	//	transform->Translate(translation);
	transform->Translate(directionOriginal*translation);
	filter->SetTransform(transform);

	// set interpolator
	InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();
	filter->SetInterpolator(interpolator);
	filter->SetDefaultPixelValue(0); // 100

	// set spacing using input
	filter->SetOutputSpacing(spacingNew);

	// set origin, always be 0 0 0 
	filter->SetOutputOrigin(originNew);

	// set size using input
	filter->SetSize(sizeNew);

	filter->SetOutputDirection(originalImage->GetDirection());

	// set input image
	filter->SetInput(originalImage);

	try
	{
		filter->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << err << std::endl; 
		exit( EXIT_FAILURE );
	} 

	InternalImageType::Pointer outputImage = 0;
	outputImage = filter->GetOutput();
	outputImage->DisconnectPipeline();

	WriteImage(resampledImageFileName, outputImage);

}



