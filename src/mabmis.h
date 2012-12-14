//
//  mabmis.h
//  
//
//  Created by Tanachat Nilanon on 12/14/12.
//
//

#ifndef ____mabmis__
#define ____mabmis__

#include <iostream>

#include "itkJia3D.h"

typedef itk::Vector< ShortPixelType, ImageDimension > ShortVectorPixelType;
typedef itk::Image< ShortVectorPixelType, ImageDimension > ShortDeformationFieldType;
typedef itk::ImageFileWriter< ShortDeformationFieldType > ShortDeformationFieldWriterType;
typedef itk::CastImageFilter< DeformationFieldType, ShortDeformationFieldType > InternalToShortDeformationFieldCastFilterType;
typedef itk::CastImageFilter< ShortDeformationFieldType, DeformationFieldType > ShortToInternalDeformationFieldCastFilterType;

#define EXIT_NOT_ENOUGH_PARAMETERS 101
#define EXIT_FILE_NOT_EXISTS 102

#endif /* defined(____mabmis__) */
