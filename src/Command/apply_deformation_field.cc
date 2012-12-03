#include <metacommand.h>

#include "itkJia3D.h"

struct Arguments {
    std::string input_image_file_name;
    std::string deformation_field_file_name;
    std::string output_image_file_name;
    unsigned int interpolation_method;
    
    friend std::ostream & operator << ( std::ostream &o, const Arguments &args ) {
        std::string interpolation_method_string;
        switch (args.interpolation_method) {
            case 0:
                interpolation_method_string = "Nearest neighbor";
                break;
                
            case 1:
                interpolation_method_string = "Linear";
                break;
                
            default:
                interpolation_method_string = "Unsupported";
                break;
        }
        
        return o
        << "Argument structure:" << std::endl
        << "  Input image file name: " << args.input_image_file_name << std::endl
        << "  Deformation field file name: " << args.deformation_field_file_name << std::endl
        << "  Output image file name: " << args.output_image_file_name << std::endl
        << "  Interpolation method: " << interpolation_method_string << std::endl;
    }
};

void HelpCallBack () {
    std::cout
    << std::endl
    << "gor@cs.unc.edu" << std::endl
    << "UNC BRIC IDEA Lab" << std::endl;
    exit(EXIT_FAILURE);
}

void ParseOptions ( int argc, char **argv, struct Arguments &args ) {
    // command line parser
    MetaCommand command;
    command.SetParseFailureOnUnrecognizedOption(true);
    command.SetHelpCallBack(HelpCallBack);
    
    command.SetAuthor("Tanachat Nilanon");
    
    command.SetOption("input_image_file_name", "i", true, "Name of input image");
    command.SetOptionLongTag("input_image_file_name", "input-image-file-name");
    command.AddOptionField("input_image_file_name", "filename", MetaCommand::STRING, true);
    
    command.SetOption("deformation_field_file_name", "f", true, "Name of deformation field");
    command.SetOptionLongTag("deformation_field_file_name", "deformation-field-file-name");
    command.AddOptionField("deformation_field_file_name", "filename", MetaCommand::STRING, true);
    
    command.SetOption("output_image_file_name", "o", true, "Name of output image");
    command.SetOptionLongTag("output_image_file_name", "output-image-file-name");
    command.AddOptionField("output_image_file_name", "filename", MetaCommand::STRING, true);
    
    command.SetOption("interpolation_method", "p", false, "Type of interpolation method; 0: Nearest neighbor, 1: Linear");
    command.SetOptionLongTag("interpolation_method", "interpolation-method");
    command.AddOptionField("interpolation_method", "type", MetaCommand::INT, true, "1");
    command.SetOptionRange("interpolation_method", "type", "0", "1");
    
    if (! command.Parse(argc, argv)) {
        exit(EXIT_FAILURE);
    }
    
    // store the parsed information
    args.input_image_file_name = command.GetValueAsString("input_image_file_name", "filename");
    args.deformation_field_file_name = command.GetValueAsString("deformation_field_file_name", "filename");
    args.output_image_file_name = command.GetValueAsString("output_image_file_name", "filename");
    args.interpolation_method = command.GetValueAsInt("interpolation_method", "type");
}

int main ( int argc, char *argv[] ) {
    struct Arguments args;
    ParseOptions(argc, argv, args);
    
    // read in image
    InternalImageType::Pointer input_image_pointer = 0;
    ReadImage(const_cast< char * >(args.input_image_file_name.c_str()), input_image_pointer);
    
    // read in deformation field
    DeformationFieldType::Pointer deformation_field_pointer = 0;
    ReadDeformationField(const_cast< char * >(args.deformation_field_file_name.c_str()), deformation_field_pointer);
    
    // apply deformation field
    InternalImageType::Pointer output_image_pointer = 0;
    ApplyDeformationField(input_image_pointer, deformation_field_pointer, output_image_pointer, args.interpolation_method == 1);
    
    // write out image
    WriteImage(const_cast< char * >(args.output_image_file_name.c_str()), output_image_pointer);
}

