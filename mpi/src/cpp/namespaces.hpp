#ifndef NAMESPACES_H
#define NAMESPACES_H

#include <map>

#include "mpi.h"

#include "getopt.h"
#include "string.h"

#include "linear_memory.hpp"
#include "helper_mpi.hpp"

namespace namespace_mpi
{
    inline int rank = -(1<<30);
    inline int size = -(1<<30);
    
    inline MPI_Comm comm_cart = MPI_COMM_NULL;

    inline std::map<char , std::map<char , int> > map_rank
    { { 'x' , std::map<char , int> { {'L' , MPI_PROC_NULL } , {'R' , MPI_PROC_NULL } } } ,
      { 'y' , std::map<char , int> { {'L' , MPI_PROC_NULL } , {'R' , MPI_PROC_NULL } } } };

    inline int & rank_xL = map_rank.at('x').at('L');
    inline int & rank_xR = map_rank.at('x').at('R');

    inline int & rank_yL = map_rank.at('y').at('L');
    inline int & rank_yR = map_rank.at('y').at('R');

    inline int ndims = 2;
    inline int dims[2] = {0,0};  // need to init as {0,0} if want to call MPI_Dims_create
    inline std::map<char , int> map_dim { {'x',-1} , {'y',-1} };
// we can use reference as map value

    inline int periods[2] = {true,true};
    inline int reorder = true;

    inline int coords[2] = {-1,-1};
    inline std::map<char , int> map_coord { {'x',-1} , {'y',-1} };
// we can use reference as map value    
}
namespace ns_mpi = namespace_mpi;


namespace namespace_precision
{
#ifndef PRCSFLAG
    using precision = double;
#endif

#if PRCSFLAG == 64
    using precision = double;
#elif PRCSFLAG == 32    
    using precision = float;
#elif PRCSFLAG == 16        
    using precision = _Float16;
#endif

    MPI_Datatype precision_MPI_Datatype = -1;
}
namespace ns_prcs = namespace_precision;
// [2023/12/20]
// NOTE: If we use a macro, the type can be input at compile time by passing a
//       flag to make.


namespace namespace_config
{
    // [2023/11/16]
    // NOTE: inline or constexpr? Using "inline" we reserve the possibility of 
    //       changing them; using "constexpr" has the potential of faster code. 
    //       /* Leaning towards "inline" now. */

    inline std::string medium_name = "homogeneous";

    inline double min_velocity = 1.;
    // [2023/11/18]
    // NOTE: "min_velocity" should be determined before entering the simulation,
    //       which is possible. /* In other words, it should not be deduced from
    //       the input parameter fields - it may be too "late" in the process as
    //       we may need to use "min_velocity" to determine "dx" as well as "dt"
    //       early on in the simulation. */ /* Although, it might not have to be 
    //       specified as "early" as here. */

    inline double central_f  = 5.;
    inline double time_delay = 0.25;

    inline int base_ppw = 10;
    inline int ppw      = -1;

    inline long N_wavelength_x = 60;
    inline long N_wavelength_y = 60;

    inline long src_loc_i_base_ppw = 20 * base_ppw;
    inline long src_loc_j_base_ppw = 20 * base_ppw;

    inline long rcv_loc_i_base_ppw = 40 * base_ppw;
    inline long rcv_loc_j_base_ppw = 40 * base_ppw;

    // [2023/11/16] 
    // NOTE: Indices should be "long" rather than "int" becasue sooner or later int
    //       will run out of range. Better to have the consistency at the beginning.

    inline int multiple = 1;
    // [2023/12/21]
    // NOTE: If the "inline" keyword is forgotten. The code might still compile and 
    //       run, but misbehave - for example, passing in "-multiple 11" at command 
    //       line may not actually override the variable "multiple" above.

    inline double dx = 0./0.;
    inline double dt = 0./0.;
    // [2023/11/16]
    // NOTE: Both "dx" and "dt" should be double here. They should NOT be "downcast"
    //       until all operators (that involves them) have been prepared.

    inline int Nt = 60000;

    inline char char_cpst = '0';  // default to no compensation

    inline int it_debug = -1;


    // [2023/11/30]
    // NOTE: We should use the unaffixed version for "subdomain variable" names (and
    //       affix the "entire domain variable" names) since the subdomain ones will
    //       be used more often.

    inline std::map<char , int> map_Nel_entire { {'x', -1} , {'y', -1} };

    inline std::map<char , std::map<char , int> > map_bgn_entire
    { { 'x' , std::map<char , int> { {'N' , -1 } , {'M' , -1 } } } ,
      { 'y' , std::map<char , int> { {'N' , -1 } , {'M' , -1 } } } };

    inline std::map<char , std::map<char , int> > map_end_entire  // "end" is one past last
    { { 'x' , std::map<char , int> { {'N' , -1 } , {'M' , -1 } } } ,
      { 'y' , std::map<char , int> { {'N' , -1 } , {'M' , -1 } } } };

    // [2023/11/30]
    // NOTE: Padding is not involved in the "_entire" variables. However, duplication
    //       on the 'N' grid still needs to be accounted for.


    inline std::map<char , int> map_Nel { {'x', -1} , {'y', -1} };

    inline int & Nel_x = map_Nel.at('x');
    inline int & Nel_y = map_Nel.at('y');


    // {'x';'y'} -> {'M','N'} -> {'L','I','R'}
    inline std::map<char, std::map<char , std::map<char , int> > > map_len_LIR;


    inline std::map<char , std::map<char , int> > map_len
    { { 'x' , std::map<char , int> { {'N' , 0 } , {'M' , 0 } } } ,
      { 'y' , std::map<char , int> { {'N' , 0 } , {'M' , 0 } } } };

    inline int & Nx = map_len.at('x').at('N');
    inline int & Mx = map_len.at('x').at('M');

    inline int & Ny = map_len.at('y').at('N');
    inline int & My = map_len.at('y').at('M');


    inline std::map<char , std::map<char , int> > map_bgn
    { { 'x' , std::map<char , int> { {'N' , -1 } , {'M' , -1 } } } ,
      { 'y' , std::map<char , int> { {'N' , -1 } , {'M' , -1 } } } };

    inline std::map<char , std::map<char , int> > map_end
    { { 'x' , std::map<char , int> { {'N' , -1 } , {'M' , -1 } } } ,
      { 'y' , std::map<char , int> { {'N' , -1 } , {'M' , -1 } } } };


    inline double rho  = -1.;
    inline double vp   = -1.;
    inline double beta = -1.;


    std::map< std::array<char,2> , std::string > char2mesh
    {
        { {'N','M'} , "Vx" } ,
        { {'M','N'} , "Vy" } ,
        { {'M','M'} , "P"  }
    };

    std::map< std::string , std::array<char,2> > mesh2char
    {
        { "Vx" , {'N','M'} } ,
        { "Vy" , {'M','N'} } ,
        { "P"  , {'M','M'} }
    };
}
namespace ns_config = namespace_config;


namespace namespace_fields
{
    memory_resource<ns_prcs::precision,2> prmt_Vx;
    memory_resource<ns_prcs::precision,2> prmt_Vy;
    memory_resource<ns_prcs::precision,2> prmt_P ;

    memory_resource<ns_prcs::precision,2> soln_Vx;
    memory_resource<ns_prcs::precision,2> soln_Vy;
    memory_resource<ns_prcs::precision,2> soln_P ;

    memory_resource<ns_prcs::precision,2> rths_Vx;
    memory_resource<ns_prcs::precision,2> rths_Vy;
    memory_resource<ns_prcs::precision,2> rths_P ;


    std::map<std::string , memory_resource<ns_prcs::precision,2> *> mesh2prmt
    {
        { "Vx" , & prmt_Vx } ,
        { "Vy" , & prmt_Vy } ,
        { "P"  , & prmt_P  }
    }; 

    std::map<std::string , memory_resource<ns_prcs::precision,2> *> mesh2soln
    {
        { "Vx" , & soln_Vx } ,
        { "Vy" , & soln_Vy } ,
        { "P"  , & soln_P  }
    };

    std::map<std::string , memory_resource<ns_prcs::precision,2> *> mesh2rths
    {
        { "Vx" , & rths_Vx } ,
        { "Vy" , & rths_Vy } ,
        { "P"  , & rths_P  }
    };

    std::map<std::string , decltype(mesh2prmt) > field2mesh2resource
    { 
        {"prmt" , mesh2prmt} , 
        {"soln" , mesh2soln} , 
        {"rths" , mesh2rths} 
    };

    
    memory_resource<ns_prcs::precision> R_signal_Vx;
    memory_resource<ns_prcs::precision> R_signal_Vy;
    memory_resource<ns_prcs::precision> R_signal_P ;


    std::map<std::string , memory_resource<ns_prcs::precision> *> mesh2signal
    {
        { "Vx" , & R_signal_Vx } ,
        { "Vy" , & R_signal_Vy } ,
        { "P"  , & R_signal_P  }
    }; 


    memory_resource<double> R_energy_Vx;
    memory_resource<double> R_energy_Vy;
    memory_resource<double> R_energy_P ;

    std::map<std::string , memory_resource<double> *> mesh2energy
    {
        { "Vx" , & R_energy_Vx } ,
        { "Vy" , & R_energy_Vy } ,
        { "P"  , & R_energy_P  }
    }; 


    // [2024/02/01]
    // NOTE: The following fields are used to store the pre-calculated source value 
    //       at each time step. They are needed because we want to see if using the
    //       pre-stored version may make a difference for the observed variation in 
    //       elapsed time for every 1000 time steps. They aren't part of class_mesh 
    //       right now.
    memory_resource<double> R_source_Vx;
    memory_resource<double> R_source_Vy;
    memory_resource<double> R_source_P ;

    std::map<std::string , memory_resource<double> *> mesh2source
    {
        { "Vx" , & R_source_Vx } ,
        { "Vy" , & R_source_Vy } ,
        { "P"  , & R_source_P  }
    }; 


    memory_resource<double,2> weight_energy_Vx;
    memory_resource<double,2> weight_energy_Vy;
    memory_resource<double,2> weight_energy_P ;

    std::map<std::string , memory_resource<double,2> *> mesh2weight_energy
    {
        { "Vx" , & weight_energy_Vx } ,
        { "Vy" , & weight_energy_Vy } ,
        { "P"  , & weight_energy_P  }
    }; 


    memory_resource<double,2> weight_source_Vx;
    memory_resource<double,2> weight_source_Vy;
    memory_resource<double,2> weight_source_P ;


    // [2023/11/21]
    // NOTE: We don't need to declare separate buffers here because those extra 
    //       spaces are already built into the field vectors.
    //           If we do need to declare and organize the buffers, they should 
    //       adhere to fields, i.e., each field should have four (in 2D) buffers
    //       associated. (In other words, the buffers should not be associated 
    //       with meshes because multiple fields may be defined on the same mesh
    //       and they may not re-use the same buffers.)
    //           We could use an object
    //           std::map<"mesh_name", std::map<"field_name", run_time_vector>>
    //       to organize the fields. prmt_Vx, as an example, could be defined as 
    //       a reference to the entry.


    // [2023/11/17]
    // NOTE: We want these "resources" to be stand-alone, i.e., not as part of
    //       any classes. Classes that need convenient access to them can store
    //       pointers to them.

}
namespace ns_fields = namespace_fields;


class class_mesh 
{
public:
    std::string name = "0000";

    std::map<char,char> dir2type { {'x' , 'U'} , {'y' , 'U'} };

    char & T_x = dir2type.at('x');  // type of the grid 'N' or 'M'
    char & T_y = dir2type.at('y');  // type of the grid 'N' or 'M'

    bool bool_src = false;
    bool bool_rcv = false;

    long I_SRC = -1;
    long I_RCV = -1;

    // [2023/11/30]
    // NOTE: Currently, we assume that there is only 1 src and 1 rcv on each 
    //       (entire) grid. If src lands on the interface, it may need to be 
    //       duplicated (maybe as many as four times in 2D when it is on the 
    //       corner). If rcv lands on the interface, the duplicate should be 
    //       removed (no need to record and write twice).

    memory_accessor<ns_prcs::precision,2> hst_prmt;
    memory_accessor<ns_prcs::precision,2> hst_soln;
    memory_accessor<ns_prcs::precision,2> hst_rths;

    std::map< std::string , memory_accessor<ns_prcs::precision,2> & > hst_map_fields
    { { "prmt" , hst_prmt } , { "soln" , hst_soln } , { "rths" , hst_rths } };
    // rename hst_map_fields to hst_field2acc?

    memory_accessor<ns_prcs::precision> hst_R_signal;
    memory_accessor<double>             hst_R_energy;


    std::map<char , std::map<char , MPI_Datatype> > type_send
    { { 'x' , std::map<char , MPI_Datatype> { {'L' , -1 } , {'R' , -1 } } } ,
      { 'y' , std::map<char , MPI_Datatype> { {'L' , -1 } , {'R' , -1 } } } };

    std::map<char , std::map<char , MPI_Datatype> > type_recv
    { { 'x' , std::map<char , MPI_Datatype> { {'L' , -1 } , {'R' , -1 } } } ,
      { 'y' , std::map<char , MPI_Datatype> { {'L' , -1 } , {'R' , -1 } } } };
    // the two have relationships that needs to be asserted

    // [2023/11/28]
    // NOTE: We should call "MPI_Type_free ()" on these MPI_Datatype. They can't
    //       be called inside "dtor ()" because by the time "dtor ()" is called, 
    //       "MPI_Finalize ()" has already been called.
    //       
    //       "https://stackoverflow.com/a/57211527" seems to indicate that if we
    //       use assignment, e.g., 
    //              type_send.at('x').at('R') = type_recv.at('x').at('L'),
    //       then we only need to invoke "MPI_Type_free ()" on the original one, 
    //       i.e., type_recv.at('x').at('L'). If we use MPI_Type_dup () instead, 
    //       then we need to call MPI_Type_free () on both.
    

    std::map<char , std::map<char , int> > send_offset
    { { 'x' , std::map<char , int> { {'L' , -1 } , {'R' , -1 } } } ,
      { 'y' , std::map<char , int> { {'L' , -1 } , {'R' , -1 } } } };

    std::map<char , std::map<char , int> > recv_offset
    { { 'x' , std::map<char , int> { {'L' , -1 } , {'R' , -1 } } } ,
      { 'y' , std::map<char , int> { {'L' , -1 } , {'R' , -1 } } } };

    // [2023/11/27]
    // NOTE: *_offset are the indices added to the starting address of a memory. 
    //       For example, soln->ptr + send_offset.at('x').at('L') gives the bgn 
    //       of the send_xL buffer.


    // ---- send_bgn, send_len, recv_bgn, recv_len are for direct copying
    std::map< char , std::map< char , std::array<int,2> > > comm_len
    { { 'x' , std::map<char , std::array<int,2>> { {'L' , {-1,-1} } , {'R' , {-1,-1} } } } ,
      { 'y' , std::map<char , std::array<int,2>> { {'L' , {-1,-1} } , {'R' , {-1,-1} } } } };

    std::map< char , std::map< char , std::array<int,2> > > recv_bgn
    { { 'x' , std::map<char , std::array<int,2>> { {'L' , {-1,-1} } , {'R' , {-1,-1} } } } ,
      { 'y' , std::map<char , std::array<int,2>> { {'L' , {-1,-1} } , {'R' , {-1,-1} } } } };

    std::map< char , std::map< char , std::array<int,2> > > send_bgn
    { { 'x' , std::map<char , std::array<int,2>> { {'L' , {-1,-1} } , {'R' , {-1,-1} } } } ,
      { 'y' , std::map<char , std::array<int,2>> { {'L' , {-1,-1} } , {'R' , {-1,-1} } } } };


    class_mesh () {}

    void initialize ( char T_x, char T_y )
    {
        this->T_x = T_x;
        this->T_y = T_y;

        this->name = ns_config::char2mesh.at( {T_x,T_y} ); 

        this->hst_prmt.attach_memory( ns_fields::mesh2prmt.at(this->name)->acc );
        this->hst_soln.attach_memory( ns_fields::mesh2soln.at(this->name)->acc );
        this->hst_rths.attach_memory( ns_fields::mesh2rths.at(this->name)->acc );

        this->hst_R_signal.attach_memory( ns_fields::mesh2signal.at(this->name)->acc );
        this->hst_R_energy.attach_memory( ns_fields::mesh2energy.at(this->name)->acc );

        // ---- determine the (1D) indices for src and rcv 
        {
            using namespace ns_config;
            class_mpi_printf mpi_printf (MPI_COMM_WORLD);

            // ---- bgn and end global indices for current subdomain
            int x_bgn = map_bgn_entire.at('x').at(T_x);  int x_end = map_end_entire.at('x').at(T_x);
            int y_bgn = map_bgn_entire.at('y').at(T_y);  int y_end = map_end_entire.at('y').at(T_y);

            // [2024/07/02]
            // NOTE: If "map_bgn_entire" and "map_end_entire" are safe to use as above, "multiple" 
            //       should also be safe to use below.

            long src_loc_i = src_loc_i_base_ppw * multiple;  if ( T_x == 'M' ) { src_loc_i += (multiple / 2); }
            long src_loc_j = src_loc_j_base_ppw * multiple;  if ( T_y == 'M' ) { src_loc_j += (multiple / 2); }

            long rcv_loc_i = rcv_loc_i_base_ppw * multiple;  if ( T_x == 'M' ) { rcv_loc_i += (multiple / 2); }
            long rcv_loc_j = rcv_loc_j_base_ppw * multiple;  if ( T_y == 'M' ) { rcv_loc_j += (multiple / 2); }
            // [2024/07/04]
            // NOTE: The adjustment in indices for the half grid shift is still applied on the 'M' 
            //       grid, regardless of padding.

            int x_pad = map_bgn.at('x').at( this->T_x );
            int y_pad = map_bgn.at('y').at( this->T_y );
            // [2024/07/04]
            // NOTE: map_bgn.at(dir).at(T) is the same as map_len_LIR.at(dir).at(T).at('L'), i.e.,
            //       padding on this grid.


            if ( src_loc_i >= x_bgn && src_loc_i < x_end )
            if ( src_loc_j >= y_bgn && src_loc_j < y_end )
            {
                if ( T_x == 'M' && T_y == 'M' ) { bool_src = true; } 
                // ---- only set bool_src to true for P (we may want a better way to interface)
                
                if ( bool_src ) 
                {
                    this->I_SRC = this->hst_soln.convert_index ( x_pad + src_loc_i - x_bgn , 
                                                                 y_pad + src_loc_j - y_bgn );
                    printf("src %2s : %ld (%ld %ld) coords : (%d %d) \n", name.c_str() , I_SRC , 
                                                                          x_pad + src_loc_i - x_bgn , 
                                                                          y_pad + src_loc_j - y_bgn ,
                                                                          ns_mpi::coords[0] , ns_mpi::coords[1] );
                }
            }


            if ( rcv_loc_i >= x_bgn && rcv_loc_i < x_end )
            if ( rcv_loc_j >= y_bgn && rcv_loc_j < y_end )
            {
                bool_rcv = true;

                if ( ( T_x == 'N' && rcv_loc_i == x_bgn && ns_mpi::dims[0] != 1 ) 
                  || ( T_y == 'N' && rcv_loc_j == y_bgn && ns_mpi::dims[1] != 1 ) ) { bool_rcv = false; }
                // ---- to avoid duplicated recording

                if ( bool_rcv )
                { 
                    this->I_RCV = this->hst_soln.convert_index ( x_pad + rcv_loc_i - x_bgn , 
                                                                 y_pad + rcv_loc_j - y_bgn );
                    printf("rcv %2s : %ld (%ld %ld) coords : (%d %d) \n", name.c_str() , I_RCV , 
                                                                          x_pad + rcv_loc_i - x_bgn ,
                                                                          y_pad + rcv_loc_j - y_bgn ,
                                                                          ns_mpi::coords[0] , ns_mpi::coords[1] );
                }
            }
        }


        for ( char side : {'L','R'} )
        {
            using namespace ns_config;

            int block_count  = map_len_LIR.at('x').at(T_x).at(side);
            int block_length = map_len.at('y').at(T_y);
            int stride = ns_config::map_len.at('y').at(T_y);

            MPI_Type_vector(block_count, block_length, stride, ns_prcs::precision_MPI_Datatype, 
                            &type_recv.at('x').at(side));
            MPI_Type_commit(&type_recv.at('x').at(side));
        }

        for ( char side : {'L','R'} )
        {
            using namespace ns_config;

            int block_count  = map_len.at('x').at(T_x);
            int block_length = map_len_LIR.at('y').at(T_y).at(side);
            int stride = ns_config::map_len.at('y').at(T_y);

            MPI_Type_vector(block_count, block_length, stride, ns_prcs::precision_MPI_Datatype, 
                            &type_recv.at('y').at(side));
            MPI_Type_commit(&type_recv.at('y').at(side));
        }

        type_send.at('x').at('L') = type_recv.at('x').at('R');
        type_send.at('x').at('R') = type_recv.at('x').at('L');

        type_send.at('y').at('L') = type_recv.at('y').at('R');
        type_send.at('y').at('R') = type_recv.at('y').at('L');


        {
            using namespace ns_config;
            send_offset.at('x').at('L') = hst_soln.convert_index( map_end.at('x').at(T_x) - map_Nel.at('x') , 0 );
            send_offset.at('x').at('R') = hst_soln.convert_index(                           map_Nel.at('x') , 0 );

            send_offset.at('y').at('L') = hst_soln.convert_index( 0 , map_end.at('y').at(T_y) - map_Nel.at('y') );
            send_offset.at('y').at('R') = hst_soln.convert_index( 0 ,                           map_Nel.at('y') );
        }


        {
            using namespace ns_config;
            recv_offset.at('x').at('R') = hst_soln.convert_index( map_end.at('x').at(T_x) , 0 );
            recv_offset.at('x').at('L') = hst_soln.convert_index(                       0 , 0 );

            recv_offset.at('y').at('R') = hst_soln.convert_index( 0 , map_end.at('y').at(T_y) );
            recv_offset.at('y').at('L') = hst_soln.convert_index(                       0 , 0 );
        }

        {
            using namespace ns_config;

            recv_bgn.at('x').at('L') = { 0 , 0 };
            recv_bgn.at('x').at('R') = { map_end.at('x').at(T_x) , 0 };

            recv_bgn.at('y').at('L') = { 0 , 0 };
            recv_bgn.at('y').at('R') = { 0 , map_end.at('y').at(T_y) };


            send_bgn.at('x').at('L') = { map_end.at('x').at(T_x) - map_Nel.at('x') , 0 };
            send_bgn.at('x').at('R') = {                           map_Nel.at('x') , 0 };

            send_bgn.at('y').at('L') = { 0 , map_end.at('y').at(T_y) - map_Nel.at('y') };
            send_bgn.at('y').at('R') = { 0 ,                           map_Nel.at('y') };
            // 
            // [2023/12/19] NOTE: Be careful of the duplication on 'N' grid.


            comm_len.at('x').at('L') = { map_len_LIR.at('x').at(T_x).at('L') , map_len.at('y').at(T_y) };
            comm_len.at('x').at('R') = { map_len_LIR.at('x').at(T_x).at('R') , map_len.at('y').at(T_y) };

            comm_len.at('y').at('L') = { map_len.at('x').at(T_x) , map_len_LIR.at('y').at(T_y).at('L') };
            comm_len.at('y').at('R') = { map_len.at('x').at(T_x) , map_len_LIR.at('y').at(T_y).at('R') };
        }
    }
};


namespace namespace_meshes
{
    class_mesh mesh_Vx;
    class_mesh mesh_Vy;
    class_mesh mesh_P ;

    std::map<std::string , class_mesh *> name2mesh
    {
        {"Vx" , & mesh_Vx } ,
        {"Vy" , & mesh_Vy } ,
        {"P"  , & mesh_P  }
    };
}
namespace ns_meshes = namespace_meshes;


namespace namespace_helper_char
{
    constexpr char flip_side(char side)
    {
             if ( side == 'L' ) { return 'R'; }
        else if ( side == 'R' ) { return 'L'; }
        else throw std::logic_error("side can only be L or R");
    }

    constexpr char flip_type(char type)
    {
             if ( type == 'M' ) { return 'N'; }
        else if ( type == 'N' ) { return 'M'; }
        else throw std::logic_error("type can only be M or N");
    }
}
namespace ns_h_char = namespace_helper_char;


namespace namespace_action
{
    template<char dir>
    void comm_to_halo ( class_mesh & mesh )
    {
        memory_accessor<ns_prcs::precision,2> & soln = mesh.hst_soln;

        for ( char send_side : {'L','R'} )
        {
            char recv_side = ns_h_char::flip_side(send_side);

            void * sendbuf = soln.ptr + mesh.send_offset.at(dir).at(send_side);
            int sendcount = 1;
            MPI_Datatype sendtype = mesh.type_send.at(dir).at(send_side);
            int dest = ns_mpi::map_rank.at(dir).at(send_side);

            void * recvbuf = soln.ptr + mesh.recv_offset.at(dir).at(recv_side);
            int recvcount = 1;
            MPI_Datatype recvtype = mesh.type_recv.at(dir).at(recv_side);
            int orig = ns_mpi::map_rank.at(dir).at(recv_side);

            MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, 1, 
                         recvbuf, recvcount, recvtype, orig, 1, 
                         ns_mpi::comm_cart, MPI_STATUS_IGNORE);
        }
    }


    template<char dir, char type>
    void take_drvt ( int ix , int iy , 
                     class_mesh & drvt_mesh , class_mesh & soln_mesh )
    {
        static_assert(  dir == 'x' ||  dir == 'y' );
        static_assert( type == 'M' || type == 'N' );

        using namespace ns_meshes;

        constexpr int stencil_length = 4;
        constexpr std::array<double, 4> stencil = {1./24 , -9./8 , 9./8 , -1./24};
        constexpr std::array<int, 4> stencil_shift = type == 'N'  
                                                   ? std::array<int, 4> {-1, 0,1,2} 
                                                   : std::array<int, 4> {-2,-1,0,1};

        std::array<ns_prcs::precision, 4> stencil_dt_dx;
        for ( int iw = 0; iw<4; iw++ ) { stencil_dt_dx[iw] = stencil[iw] * ns_config::dt 
                                                                         / ns_config::dx; }

        static_assert(stencil_length == stencil.size());
        static_assert(stencil_length == stencil_shift.size());
        assert( type == drvt_mesh.dir2type.at(dir) );  // we should not assert here, too costly (unless static)

        ns_prcs::precision v_d = 0.;
        for ( long iw=0; iw<stencil_length; iw++ )
        {
            int ix_s = dir == 'x' ? ix + stencil_shift[iw] : ix;
            int iy_s = dir == 'y' ? iy + stencil_shift[iw] : iy;

            v_d += soln_mesh.hst_soln.at( ix_s , iy_s ) * stencil_dt_dx[iw];
        }
        long ind = drvt_mesh.hst_rths.convert_index(ix,iy);
        drvt_mesh.hst_rths.at(ind) += v_d * drvt_mesh.hst_prmt.at(ind);
        // [2023/11/28] NOTE: We pre-calculate ind so that it can be re-used by 
        //                    rths and prmt instead of calculating twice inside 
        //                    at ().
    }


    template<typename T>
    void sum_3op( T const a , T const b ,
                  T &     s , T &     t )
    {
          s = a + b;
        T z = s - a;
          t = b - z;
    }


    template<typename T>
    void sum_6op( T const a , T const b ,
                  T &     s , T &     t )
    {
            s =   a + b;
        T p_a =   s - b;
        T p_b =   s - p_a;

        T d_a =   a - p_a;
        T d_b =   b - p_b;

            t = d_a + d_b;
    } 


    inline void update_soln ( int ix , int iy , class_mesh & mesh )
    {
        char cpst = ns_config::char_cpst;

        memory_accessor<ns_prcs::precision,2> & soln = mesh.hst_soln;
        memory_accessor<ns_prcs::precision,2> & rths = mesh.hst_rths;        

        long ind = soln.convert_index(ix,iy);
        if ( cpst == '0' ) { soln.at(ind) += rths.at(ind); rths.at(ind) = 0.; }
        if ( cpst == '3' ) { sum_3op<ns_prcs::precision> ( soln.at(ind) , rths.at(ind) , soln.at(ind) , rths.at(ind) ); }
        if ( cpst == '6' ) { sum_6op<ns_prcs::precision> ( soln.at(ind) , rths.at(ind) , soln.at(ind) , rths.at(ind) ); }
    }
    // [2023/11/23]
    // NOTE: The three index calculations can be combined, we may even pass in
    //       the index directly as a single argument such that both take_drvt 
    //       and update_soln can reuse it.


    inline void cleanup_MPI ()
    {
        // ---- calll MPI_Type_free on send recv types
        for ( auto iter : ns_meshes::name2mesh )
        {
            class_mesh & mesh = * iter.second;

            for ( char  dir : {'x','y'} )
            for ( char side : {'L','R'} )
            {
                MPI_Type_free( &mesh.type_recv.at(dir).at(side) );
                MPI_Type_free( &mesh.type_send.at(dir).at(side) );
            }
        }
        // [2023/11/28]
        // NOTE: Cannot call these inside class dtor because dtor gets called 
        //       after MPI_Finalize ();
        //       
        //       I don't think we need to call MPI_Type_free () on type_send, 
        //       because it is assigned from type_recv, but calling it seems 
        //       to be ok, i.e., free "twice" seems to be alright.
    }
}
namespace ns_action = namespace_action;


namespace namespace_inputs
{
    inline void command_line_input_processing( int argc, char* argv[] )
    {
        using namespace ns_config;

        // The following struct option is defined in #include <getopt.h>
        static struct option long_options[] =
        {
            {"medium",    required_argument,  0,  0},
            {"Nt",        required_argument,  0,  0},
            {"multiple",  required_argument,  0,  0},
            {"procs",     required_argument,  0,  0},
            {"cpst",      required_argument,  0,  0},
            {0, 0, 0, 0}
        };

        optind = 1;              // 'default' (initial state) is 1; 
                                 // 'reset' to 1 here in case this is a rescan
        while (1)
        {
            int option_index = -1;
            int opt = getopt_long_only (argc, argv, "", long_options, &option_index);
            // optarg in the following is a "system" defined (char *) pointer

            if (opt == 0)  // opt == 0: a valid command line option is captured
            {
                if ( strcmp(long_options[option_index].name, "medium"  ) == 0 ) { medium_name = optarg; }

                if ( strcmp(long_options[option_index].name, "Nt"      ) == 0 ) { Nt = int( strtol( optarg , nullptr, 10 ) ); }

                if ( strcmp(long_options[option_index].name, "multiple") == 0 ) { multiple = int( strtol( optarg , nullptr, 10 ) ); }

                if ( strcmp(long_options[option_index].name, "procs"   ) == 0 ) 
                { 
                    std::stringstream ss (optarg);
                    std::string word;
                    int count = 0;
                    while ( ss >> word )
                    {
                        ns_mpi::dims[count] = int( strtol( word.c_str() , nullptr, 10 ) );
                        count++;
                    }
                    assert( count == ns_mpi::ndims && "input needs to match ndims (did you include numbers in quotes)" );
                }

                if ( strcmp(long_options[option_index].name, "cpst"    ) == 0 ) { char_cpst = optarg[0]; }
            }
            if ( opt == -1 ) break;
        }
    }
}
namespace ns_inputs = namespace_inputs;


namespace namespace_energy
{
    inline void prepare_energy_weight ()
    {
        for ( auto mesh : { "Vx" , "Vy" , "P" } )
        {
            char T_x = ns_config::mesh2char.at(mesh).at(0);
            char T_y = ns_config::mesh2char.at(mesh).at(1);

            long len_x = ns_config::map_len.at('x').at(T_x);
            long len_y = ns_config::map_len.at('y').at(T_y);

            auto & weight = ns_fields::mesh2weight_energy.at(mesh);
            weight->allocate_memory({len_x,len_y});

            const long bgn_x = ns_config::map_bgn.at('x').at(T_x);
            const long bgn_y = ns_config::map_bgn.at('y').at(T_y);

            const long end_x = ns_config::map_end.at('x').at(T_x);
            const long end_y = ns_config::map_end.at('y').at(T_y);

            using ns_config::dx;
            weight->acc.set_constant( dx * dx / 2. );

            for ( int ix = 0; ix < len_x; ix++ )
            for ( int iy = 0; iy < len_y; iy++ )
            {
                if ( ix < bgn_x || ix >= end_x || iy < bgn_y || iy >= end_y )
                    { ns_fields::mesh2weight_energy.at(mesh)->acc.at(ix,iy) = 0.; }
            }

            for ( int ix = bgn_x; ix < end_x; ix++ )
            if ( T_y == 'N' )
            {
                ns_fields::mesh2weight_energy.at(mesh)->acc.at(ix,bgn_y  ) /= 2.;
                ns_fields::mesh2weight_energy.at(mesh)->acc.at(ix,end_y-1) /= 2.;
            }

            for ( int iy = bgn_y; iy < end_y; iy++ )
            if ( T_x == 'N' )
            {
                ns_fields::mesh2weight_energy.at(mesh)->acc.at(bgn_x  ,iy) /= 2.;
                ns_fields::mesh2weight_energy.at(mesh)->acc.at(end_x-1,iy) /= 2.;
            }
        }
    }
}
namespace ns_energy = namespace_energy;
// [2024/01/28]
// NOTE: The above energy weights may work in the MPI environment without change? 


namespace namespace_experimental
{
    // place for experimental code
}
namespace ns_exp = namespace_experimental;


// [2023/11/28]
// NOTE: Keep the early namespaces "variable only" and place functions/actions 
//       in later namespaces, this would help with the inter-dependence issue. 


// [2025/07/15]
// NOTE: Prefer standalone (template) functions over class member functions.


#endif