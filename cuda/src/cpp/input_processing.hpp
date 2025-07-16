#ifndef INPUT_PROCESSING_H
#define INPUT_PROCESSING_H

#include <iostream>

// ---- processing the input parameters from file
void file_input_processing ( std::string file_name );

// ---- processing the input arguments from command line
void command_line_input_processing ( int argc, char* argv[], std::string argument_name );

void checking_and_printouts_input_parameters ();


#endif