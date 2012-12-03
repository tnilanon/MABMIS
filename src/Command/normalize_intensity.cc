#include <metacommand.h>

#include "itkJia3D.h"

namespace itk {
	namespace Functor {
		template< class TInput1, class TInput2 = TInput1, class TOutput = TInput1 >
		class DivideHack {
		public:
			DivideHack() {}
			~DivideHack() {}
            
			bool operator != ( const DivideHack & ) const
			{
				return false;
			}
            
			bool operator == ( const DivideHack & other ) const
			{
				return !( *this != other );
			}
            
			inline TOutput operator () ( const TInput1 & A, const TInput2 & B ) const
			{
				TOutput result;
				if(B < 1E-6)
					result = 0;
				else
					result = A / B;
				return static_cast< TOutput >(result);
			}
		};
	}
}

typedef itk::BinaryFunctorImageFilter< InternalImageType, InternalImageType, InternalImageType, itk::Functor::DivideHack< InternalImageType::PixelType > > DivideHackType;

struct Arguments {
    std::string input_fixed_image_file_name;
    std::string input_moving_image_file_name;
    std::string output_fixed_image_file_name;
    std::string output_moving_image_file_name;
    
    friend std::ostream & operator << ( std::ostream &o, const Arguments &args ) {
        return o
        << "Arguments structure:" << std::endl
        << "  Input fixed image file name: " << args.input_fixed_image_file_name << std::endl
        << "  Input moving image file name: " << args.input_moving_image_file_name << std::endl
        << "  Output fixed image file name: " << args.output_fixed_image_file_name << std::endl
        << "  Output moving image file name: " << args.output_moving_image_file_name << std::endl;
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
    
    command.SetOption("input_fixed_image_file_name", "i", true, "Name of input fixed image");
    command.SetOptionLongTag("input_fixed_image_file_name", "input-fixed-image-file-name");
    command.AddOptionField("input_fixed_image_file_name", "filename", MetaCommand::STRING, true);
    
    command.SetOption("input_moving_image_file_name", "I", true, "Name of input moving image");
    command.SetOptionLongTag("input_moving_image_file_name", "input-moving-image-file-name");
    command.AddOptionField("input_moving_image_file_name", "filename", MetaCommand::STRING, true);
    
    command.SetOption("output_fixed_image_file_name", "o", true, "Name of output fixed image");
    command.SetOptionLongTag("output_fixed_image_file_name", "output-fixed-image-file-name");
    command.AddOptionField("output_fixed_image_file_name", "filename", MetaCommand::STRING, true);
    
    command.SetOption("output_moving_image_file_name", "O", true, "Name of output moving image");
    command.SetOptionLongTag("output_moving_image_file_name", "output-moving-image-file-name");
    command.AddOptionField("output_moving_image_file_name", "filename", MetaCommand::STRING, true);
    
    if (! command.Parse(argc, argv)) {
        exit(EXIT_FAILURE);
    }
    
    // store the parsed information
    args.input_fixed_image_file_name = command.GetValueAsString("input_fixed_image_file_name", "filename");
    args.input_moving_image_file_name = command.GetValueAsString("input_moving_image_file_name", "filename");
    args.output_fixed_image_file_name = command.GetValueAsString("output_fixed_image_file_name", "filename");
    args.output_moving_image_file_name = command.GetValueAsString("output_moving_image_file_name", "filename");
}

int main( int argc, char *argv[] ) {
    const float kDenominator = 255;
    
    struct Arguments args;
    ParseOptions(argc, argv, args);
    
    // read in image
    InternalImageType::Pointer input_fixed_image_pointer = 0;
    InternalImageType::Pointer input_moving_image_pointer = 0;
    ReadImage(const_cast< char * >(args.input_fixed_image_file_name.c_str()), input_fixed_image_pointer);
    ReadImage(const_cast< char * >(args.input_moving_image_file_name.c_str()), input_moving_image_pointer);
    
    // histogram matching
    InternalImageType::Pointer histogram_matched_moving_image_pointer = 0;
    HistogramMatching(input_moving_image_pointer, input_fixed_image_pointer, histogram_matched_moving_image_pointer);
    
    // construct denominator field
    InternalImageType::Pointer denominator_image_pointer = InternalImageType::New();
    InternalImageType::RegionType image_region(input_fixed_image_pointer->GetLargestPossibleRegion().GetIndex(), input_fixed_image_pointer->GetLargestPossibleRegion().GetSize());
    denominator_image_pointer->SetRegions(image_region);
    denominator_image_pointer->Allocate();
    denominator_image_pointer->FillBuffer(kDenominator);
    
    // divide them and get output
    DivideHackType::Pointer divider = DivideHackType::New();
    divider->SetInput2(denominator_image_pointer);
    
    divider->SetInput1(input_fixed_image_pointer);
    divider->Update();
    InternalImageType::Pointer output_fixed_image_pointer = divider->GetOutput();
    output_fixed_image_pointer->DisconnectPipeline();
    
    divider->SetInput1(histogram_matched_moving_image_pointer);
    divider->Update();
    InternalImageType::Pointer output_moving_image_pointer = divider->GetOutput();
    output_moving_image_pointer->DisconnectPipeline();
    
    // write output
    WriteImage(const_cast< char * >(args.output_fixed_image_file_name.c_str()), output_fixed_image_pointer, const_cast< char * >("float"));
    WriteImage(const_cast< char * >(args.output_moving_image_file_name.c_str()), output_moving_image_pointer, const_cast< char * >("float"));
    
    return EXIT_SUCCESS;
}

