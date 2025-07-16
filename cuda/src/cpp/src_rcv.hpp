#ifndef SRC_RCV_H
#define SRC_RCV_H

#include "grid.hpp"

struct struct_src_input
{
	int                                       src_index;        // 'global' index, uniquely identifies this source
    std::array<char, ns_forward::N_dir     >  src_grid_type;    // grid type for this source
	std::array< int, ns_forward::N_dir * 2 >  src_location;     // 'global' location in 'physical' unit 
};                                                              // * 2 : 1 for num ; 1 for den


struct struct_rcv_input
{
    int                                       src_index;        // 'global' index, uniquely identifies the source this receiver belongs to

    int                                       rcv_index;        // index uniquely identifies this receiver within this source 
    std::array<char, ns_forward::N_dir     >  rcv_grid_type;    // grid type for this receiver
	std::array< int, ns_forward::N_dir * 2 >  rcv_location;     // 'global' location in 'physical' unit
};                                                              // * 2 : 1 for num ; 1 for den    
 

struct struct_src_forward
{
    bool is_src = false;

    int src_index = -1;
    std::array<char, ns_forward::N_dir> grid_type {'U','U'};
    
    int  I_SRC  = -1;    

    // the following three variables may not be needed
    int ix_src = -1;
    int iy_src = -1;
    
    // c stands for character ('S' and 'V' are valid options)
    char c_source_type = 'U'; 
    // Development note : c_source_type can be deduced from grid_type.
};


struct struct_rcv_forward
{
    int rcv_index = -1;

    int I_RCV = -1;
    
    int ix_rcv = -1;
    int iy_rcv = -1;
};


#endif