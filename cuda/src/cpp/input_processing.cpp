#include "input_processing.hpp"
#include "namespace_input.hpp"

#include <getopt.h>

using namespace ns_input;

//-----------------------------------------------//
//------------- Function defintiion -------------//
//-----------------------------------------------//
void file_input_processing ( std::string file_name )
{
    std::string input_line;
    std::ifstream input_file(file_name);
    if ( input_file.is_open() )
    {   
        while ( getline ( input_file , input_line ) )
        {
            // std::cout << input_line << std::endl;                // print the unprocessed line to screen

            std::string string_delimiter = "_";                     // use "_" as delimiter in input file to avoid ambiguity
            size_t bgn_pos, end_pos;

            // first 'content' goes to input_name, which is used to decide how the following numbers should be processed
            bgn_pos = 0;
            end_pos = input_line.find( string_delimiter , bgn_pos );                
            std::string input_name = input_line.substr( bgn_pos , end_pos - bgn_pos );

            if ( strcmp( input_name.c_str(), "model" ) == 0 )
            {
                // for ( int i=0; i<ns_forward::N_dir; i++ )
                for ( const char& c_dir : {'X','Y'} )
                {
                    int i_dir = ns_forward::XY_to_01(c_dir);

                    bgn_pos = end_pos + string_delimiter.length();
                    end_pos = input_line.find( string_delimiter , bgn_pos );
                    prmt_M_sizes.at(i_dir) = strtol( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str(), nullptr, 10 );
                }
            }

            if ( strcmp( input_name.c_str(), "grids" ) == 0 )
            {
                // for ( int i=0; i<ns_forward::N_dir; i++ )
                for ( const char& c_dir : {'X','Y'} )
                {
                    int i_dir = ns_forward::XY_to_01(c_dir);

                    bgn_pos = end_pos + string_delimiter.length();
                    end_pos = input_line.find( string_delimiter , bgn_pos );
                    soln_M_sizes.at(i_dir) = strtol( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str(), nullptr, 10 );
                }
            }

            if ( strcmp( input_name.c_str(), "dxM" ) == 0 )
            {
                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                num_dx_prmt = strtol( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str(), nullptr, 10 );

                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                den_dx_prmt = strtol( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str(), nullptr, 10 );
            }

            if ( strcmp( input_name.c_str(), "dxS" ) == 0 )
            {
                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                num_dx_soln = strtol( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str(), nullptr, 10 );

                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                den_dx_soln = strtol( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str(), nullptr, 10 );
            }

            if ( strcmp( input_name.c_str(), "bdry" ) == 0 ) 
            {
                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                std::string B_type = input_line.substr( bgn_pos , end_pos - bgn_pos );

                { int i_dir = ns_forward::XY_to_01('X');  bdry_type_L.at(i_dir) = toupper( B_type[0] );  bdry_type_R.at(i_dir) = toupper( B_type[1] ); }
                { int i_dir = ns_forward::XY_to_01('Y');  bdry_type_L.at(i_dir) = toupper( B_type[2] );  bdry_type_R.at(i_dir) = toupper( B_type[3] ); }
            }

            if ( strcmp( input_name.c_str(), "Nt" ) == 0 ) 
            {
                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                Nt = strtol( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str(), nullptr, 10 );
            } 

            if ( strcmp( input_name.c_str(), "MAXdt" ) == 0 ) 
            {
                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                dt_max = std::stod( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str() ); 
            } 

            if ( strcmp( input_name.c_str(), "CFL" ) == 0 ) 
            {
                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                CFL_constant = std::stod( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str() ); 
            } 

            if ( strcmp( input_name.c_str(), "frequency" ) == 0 ) 
            {
                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                central_f = std::stod( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str() ); 
            } 

            if ( strcmp( input_name.c_str(), "delay" ) == 0 ) 
            {
                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                time_delay = std::stod( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str() ); 
            } 

            if ( strcmp( input_name.c_str(), "device" ) == 0 ) 
            {
                bgn_pos = end_pos + string_delimiter.length();
                end_pos = input_line.find( string_delimiter , bgn_pos );
                device_number = strtol( input_line.substr( bgn_pos , end_pos - bgn_pos ).c_str(), nullptr, 10 );
            } 

        }
        input_file.close();
    } // if ( input_file.is_open() )
    else 
    {
        printf( "Unable to open file: %s\n", file_name.c_str() ); fflush(stdout); exit(0); 
    }
    
} // file_input_processing()


void command_line_input_processing ( int argc, char* argv[], std::string argument_name )
{
    // The following struct option is defined in #include <getopt.h>
    static struct option long_options[] =
    {
        {"medium",   required_argument,  0,  0},
        {"B_type",   required_argument,  0,  0},
        {"Nt",       required_argument,  0,  0},
        {"dt_max",   required_argument,  0,  0},
        {"CFL",      required_argument,  0,  0},
        {"energy",   required_argument,  0,  0},
        {"f",        required_argument,  0,  0},
        {"t0",       required_argument,  0,  0},
        {"device",   required_argument,  0,  0},
        {0, 0, 0, 0}
    };

    optind = 1;              // 'default' (initial state) is 1; 
                             // 'reset' to 1 here in case this is a rescan
    while (1)
    {
        int option_index = -1;
        int opt = getopt_long_only (argc, argv, "", long_options, &option_index);        
        // optarg in the following is a 'system' defined (char *) pointer

        if (opt == 0)
        {
            // opt == 0: a valid command line option is captured
            // Next: if the argument_name passed in is "all" or the argument_name passed in
            //       matches the captured command line option name, we process the captured
            //       command line option (which involves comparing the captured option name 
            //       against the possible option names).
            if ( strcmp( argument_name.c_str(), "all"                           ) == 0 
              || strcmp( argument_name.c_str(), long_options[option_index].name ) == 0 )
            {
                if ( strcmp(long_options[option_index].name, "medium") == 0 ) { medium_name = optarg; }

                if ( strcmp(long_options[option_index].name, "B_type") == 0 ) 
                { 
                    { int i_dir = ns_forward::XY_to_01('X');  bdry_type_L.at(i_dir) = toupper( optarg[0] );  bdry_type_R.at(i_dir) = toupper( optarg[1] ); }
                    { int i_dir = ns_forward::XY_to_01('Y');  bdry_type_L.at(i_dir) = toupper( optarg[2] );  bdry_type_R.at(i_dir) = toupper( optarg[3] ); }
                }

                if ( strcmp(long_options[option_index].name, "Nt"    ) == 0 ) { Nt = int( strtol( optarg , nullptr, 10 ) ); } 

                if ( strcmp(long_options[option_index].name, "dt_max") == 0 ) { dt_max = std::stod( optarg ); } 

                if ( strcmp(long_options[option_index].name, "CFL")    == 0 ) { CFL_constant = std::stod( optarg ); } 

                if ( strcmp(long_options[option_index].name, "energy") == 0 ) { c_energy = toupper(optarg[0]); } 

                if ( strcmp(long_options[option_index].name, "f")      == 0 ) { central_f = std::stod( optarg ); } 

                if ( strcmp(long_options[option_index].name, "t0")     == 0 ) { time_delay = std::stod( optarg ); } 

                if ( strcmp(long_options[option_index].name, "device") == 0 ) { device_number = int( strtol( optarg , nullptr, 10 ) ); } 
            }
        }
        if ( opt == -1 ) break;
    }
}


void checking_and_printouts_input_parameters ()
{
    // check that the soln and prmt sizes match
    for ( int i_dir = 0; i_dir < ns_forward::N_dir; i_dir++ ) 
    {
        if ( prmt_M_sizes.at(i_dir) * num_dx_prmt * den_dx_soln 
          != soln_M_sizes.at(i_dir) * num_dx_soln * den_dx_prmt )
            { printf("Model size needs to match simulation size on direction %d.\n", i_dir); fflush(stdout); exit(0); }
    }
    
    // check the input character c_energy
    if ( c_energy != 'N' && c_energy != 'Y' )
        { printf("Input c_energy can only be N or Y.\n"); fflush(stdout); exit(0); }
    
    // print out some messages about the input parameters
    printf("prmt sizes: %d %d %f \n", prmt_M_sizes.at(0), prmt_M_sizes.at(1), static_cast<double>(num_dx_prmt) / den_dx_prmt );
    printf("soln sizes: %d %d %f \n", soln_M_sizes.at(0), soln_M_sizes.at(1), static_cast<double>(num_dx_soln) / den_dx_soln );

    printf( "Boundary types (X Y): %c%c %c%c\n", bdry_type_L.at( ns_forward::XY_to_01 ('X') ), bdry_type_R.at( ns_forward::XY_to_01 ('X') ) , 
                                                 bdry_type_L.at( ns_forward::XY_to_01 ('Y') ), bdry_type_R.at( ns_forward::XY_to_01 ('Y') ) );

    printf( "Boundary types (x y): %c%c %c%c\n", bdry_type_L.at(0), bdry_type_R.at(0), 
                                                 bdry_type_L.at(1), bdry_type_R.at(1) );    

    printf( "Nt (time steps): %d \n", Nt );

    printf( "CFL_constant: %f\n", CFL_constant );

    printf( "Calculate energy: %c\n", c_energy );

    printf( "\n" );
}