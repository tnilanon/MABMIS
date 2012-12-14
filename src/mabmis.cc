//
//  mabmis.cc
//  
//
//  Created by Tanachat Nilanon on 12/14/12.
//
//

#include <vcl_iostream.h>
#include <vnl/vnl_random.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matlab_print.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_svd_economy.h>

#include "mabmis.h"

// global boolean constants
const bool is_debug = true;

// global constants
//const int 

// global variables
int input_atlas_size;
int input_sample_size;
int simulate_atlas_size;

int main( int argc, char *argv[] ) {
    
    if (argc != 3) {
        std::cerr << "Missing parameters for Multi-Atlas-Based Multi-Image Segmentation: " << std::endl;
        std::cerr << "Usage: MABMIS atlas.txt sample.txt" << std::endl;
        std::cerr << std::endl;
        std::cerr << "  atlas.txt   : file name of the file containing atlases' prefixes" << std::endl;
        std::cerr << "  sample.txt  : file name of the file containing samples' prefixes" << std::endl;
        std::cerr << std::endl;
        std::cerr << "For example, if prefix is na01, there have to be na01_cbq_000.nii (intensity image). Furthermore, for atlases, there have to be na01_seg_000.nii (segmentation image)" << std::endl;
        return EXIT_NOT_ENOUGH_PARAMETERS;
    }
    
    ////////////////////////////////////////////////////////////////
    // read in all file prefixes
    
    if (! itksys::SystemTools::FileExists(argv[1], true)) {
        std::cerr << "File not exists: " << argv[1] << "!" << std::endl;
        return EXIT_FILE_NOT_EXISTS;
    }
    
    if (! itksys::SystemTools::FileExists(argv[2], true)) {
        std::cerr << "File not exists: " << argv[2] << "!" << std::endl;
        return EXIT_FILE_NOT_EXISTS;
    }
    
    std::ifstream file;
    
    std::cout << "Reading prefixes for atlases" << std::endl;
    file.open(argv[1]);
    file.seekg(0, std::ios::end);
    input_atlas_size = ceil(file.tellg() / 3);
    file.seekg(0, std::ios::beg);
    file.close();
    
    simulate_atlas_size = 2 * input_atlas_size;
    
    std::cout << "Reading prefixes for samples" << std::endl;
    file.open(argv[2]);
    file.seekg(0, std::ios::end);
    input_sample_size = ceil(file.tellg() / 3);
    file.seekg(0, std::ios::beg);
    file.close();
}

