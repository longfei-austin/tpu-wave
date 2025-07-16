#include <iostream>
#include <array>
#include <bitset>
#include <filesystem>
#include <algorithm>
#include <string>
#include <numeric>
#include <chrono>

#include "namespaces.hpp"

double tic,toc;
double Wtime ()
    { return std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::steady_clock::now().time_since_epoch() ).count() / 1000000.; }

int main(int argc, char* argv[]) 
{
#ifdef NDEBUG
    printf("NDEBUG is defined.\n");  // NDEBUG would turn off the "assert"s.
#endif

#if PRCSFLAG == 64
    printf("PRCSFLAG is 64 -> f64.\n");
#elif PRCSFLAG == 32    
    printf("PRCSFLAG is 32 -> f32.\n");
#elif PRCSFLAG == 16        
    printf("PRCSFLAG is 16 -> f16.\n");
#endif


tic = Wtime ();
    {
        using namespace ns_mpi;

        rank = 0;
        size = 1;

        if ( rank==0 ) { printf("multiple: %d Nt: %d\n", ns_config::multiple, ns_config::Nt); }


        map_dim.at('x') = dims[0];
        map_dim.at('y') = dims[1];

        assert( dims[0] * dims[1] == size );
        if ( rank == 0 ) { printf("dims: %d %d\n", dims[0], dims[1]); }


        map_coord.at('x') = coords[0];
        map_coord.at('y') = coords[1];

        printf("rank: %d ---- coords: (%d %d)\n", rank, coords[0], coords[1]);
        if (rank == 0) { printf("\n"); }


        // retrieve neighbors on x direction
        {
            ns_mpi::rank_xL = 0;
            ns_mpi::rank_xR = 0;
            printf("rank: %d ---- rank_xL: %d rank_xR: %d\n", rank, rank_xL, rank_xR);
        }

        // retrieve neighbors on y direction
        {
            ns_mpi::rank_yL = 0;
            ns_mpi::rank_yR = 0;
            printf("rank: %d ---- rank_xL: %d rank_xR: %d\n", rank, rank_yL, rank_yR);
        }
    }


    // ---- "prepare" the values for variables from ns_config
    {
        using namespace ns_config;
        ppw = base_ppw * multiple;


        dx = ( min_velocity / ( central_f * 2.5 ) ) / ppw;
        dt = 1e-4;
        // dt = 0.8 * (dx/1.) * (6./7.) / sqrt(2.);  // 1. is the max of vp below.
        // // [2024/06/20]
        // // NOTE: Relating dt to dx so that we can do large scale runs without 
        //          worrying about CFL. (6./7.) comes from the stencil; sqrt(2.) 
        //          comes from the dimension; 0.8 is the extra slack.

        rho  = 1.;
        vp   = 1.;
        beta = 1. / ( rho * vp * vp );


        map_Nel_entire.at('x') = N_wavelength_x * ppw;  // <- Nel_x
        map_Nel_entire.at('y') = N_wavelength_y * ppw;  // <- Nel_y


        for ( const char & dir : { 'x' , 'y' } )
        {
            const int dim_dir  = ns_mpi::map_dim.at(dir);
            const int quotient = map_Nel_entire.at(dir) / dim_dir;
            const int leftover = map_Nel_entire.at(dir) - dim_dir * quotient;

            std::vector<int> Nel_vec ( dim_dir , quotient );
            for ( int i = 0; i < leftover; i++ ) { Nel_vec.at(i) += 1; }

            map_Nel.at(dir) = Nel_vec.at( ns_mpi::map_coord.at(dir) );

            map_bgn_entire.at(dir).at('M') = std::accumulate( Nel_vec.begin(), 
                                                              Nel_vec.begin() + ns_mpi::map_coord.at(dir), 0 );
            map_end_entire.at(dir).at('M') = map_bgn_entire.at(dir).at('M') + map_Nel.at(dir);

            map_bgn_entire.at(dir).at('N') = map_bgn_entire.at(dir).at('M');

            map_end_entire.at(dir).at('N') = map_end_entire.at(dir).at('M') + 1;
        }


        for ( const char type : { 'M' , 'N' } )
        for ( const char dir  : { 'x' , 'y' } )
        {
            printf( "%c%c coords : (%d %d)  (bgn,end) : (% 4d % 4d) \n", dir, type, 
                    ns_mpi::coords[0], ns_mpi::coords[1], 
                    map_bgn_entire.at(dir).at(type), map_end_entire.at(dir).at(type) );
        }


        for ( char dir : {'x','y'} )
        {
            map_len_LIR[dir]['N']['I'] = map_Nel.at(dir) + 1;
            map_len_LIR[dir]['N']['L'] = 1;
            map_len_LIR[dir]['N']['R'] = 1;

            map_len_LIR[dir]['M']['I'] = map_Nel.at(dir);
            map_len_LIR[dir]['M']['L'] = 2;
            map_len_LIR[dir]['M']['R'] = 2;

            // [2025/07/15]
            // NOTE: 'L' and 'R' are the padded size (for communication).
        }


        for ( char dir : {'x','y'} )
        for ( char T : {'N' , 'M'} )
        {
            assert ( ( map_len_LIR.at(dir).at(T).at('I') >= map_len_LIR.at(dir).at(T).at('L') 
                    && map_len_LIR.at(dir).at(T).at('I') >= map_len_LIR.at(dir).at(T).at('R') )
                    && "There needs to be enough interior points to fill in neighbors halo.") ;
        }


        for ( char dir : {'x','y'} )
        for ( char T : {'N' , 'M'} )
        {
            map_len.at(dir).at(T) = map_len_LIR.at(dir).at(T).at('L')
                                  + map_len_LIR.at(dir).at(T).at('I')
                                  + map_len_LIR.at(dir).at(T).at('R');

            map_bgn.at(dir).at(T) = 0                     + map_len_LIR.at(dir).at(T).at('L');
            map_end.at(dir).at(T) = map_len.at(dir).at(T) - map_len_LIR.at(dir).at(T).at('R');
        }
    }
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n setup time : %5.4f\n", toc - tic); }


tic = Wtime ();
    {
        using namespace ns_config;

        for ( auto field : { "prmt" , "soln" , "rths" } )
        for ( auto mesh  : { "Vx" , "Vy" , "P" } )
        {
            std::array<char,2> mesh_type = ns_config::mesh2char.at(mesh);
            long len_x = map_len.at('x').at( mesh_type.at(0) );
            long len_y = map_len.at('y').at( mesh_type.at(1) );

            memory_resource<ns_prcs::precision,2> * res = ns_fields::field2mesh2resource.at(field).at(mesh);
            res->allocate_memory ( {len_x , len_y} );
        }

        ns_fields::prmt_Vx.acc.set_constant( 1. / rho  );
        ns_fields::prmt_Vy.acc.set_constant( 1. / rho  );
        ns_fields::prmt_P .acc.set_constant( 1. / beta );
    }
    // [2023/11/16]
    // NOTE: The above variables are declared in the namespace so that they can be 
    //       shared by different files. They are assigned values above because we 
    //       cannot do operations inside namespace. Their values should be assigned
    //       at the beginning of the main function (early in the execution), before 
    //       any other files may have needed them.
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n prepare field time: %5.4f\n", toc - tic); }


tic = Wtime ();
    {
        for ( auto mesh : { "Vx" , "Vy" , "P" } )
        { 
            ns_fields::mesh2signal.at(mesh)->allocate_memory (ns_config::Nt); 
            ns_fields::mesh2energy.at(mesh)->allocate_memory (ns_config::Nt); 
        }
    }
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n allocate signal time : %5.4f\n", toc - tic); }


tic = Wtime ();
    {
        ns_energy::prepare_energy_weight ();
    }
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n prepare energy time : %5.4f\n", toc - tic); }


tic = Wtime ();
    // ---- prepare the Grid class (they are convenient, not essential)
    {
        using namespace ns_meshes;
        using namespace ns_fields;

        mesh_Vx.initialize( 'N', 'M' );
        mesh_Vy.initialize( 'M', 'N' );
        mesh_P .initialize( 'M', 'M' );

        for ( const auto & iter : name2mesh )
            { assert( iter.first == iter.second->name ); }
    }
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n prepare mesh time : %5.4f\n", toc - tic); }


tic = Wtime (); toc = tic;

    for (int it=0; it<ns_config::Nt; it++)
    {
        if ( (it % 1000) == 0 ) 
        {   
            if ( ns_mpi::rank==0 ) { printf("\n %6d", it); fflush(stdout); }
            for ( auto name : { "Vx" , "Vy" , "P" } )
            { 
                class_mesh & mesh = * ns_meshes::name2mesh.at(name);
                if ( mesh.bool_rcv ) { printf(" % 16.15le   ", (double) mesh.hst_soln(mesh.I_RCV)); }
            }
            if ( ns_mpi::rank==0 ) { printf("\n"); fflush(stdout); }
            if ( ns_mpi::rank==0 ) 
                { double ttoc = Wtime (); printf("\n time for 1000 time steps : %5.4f\n", ttoc - toc); toc = ttoc; fflush(stdout); }
        }
        

        // update soln_Vx
        {
            ns_action::copy_to_halo<'x'> ( ns_meshes::mesh_P );

            using namespace ns_config;
            int bgn_x = map_bgn.at('x').at('N');    int bgn_y = map_bgn.at('y').at('M');
            int end_x = map_end.at('x').at('N');    int end_y = map_end.at('y').at('M');

            for ( long ix=bgn_x; ix<end_x; ix++ )
            for ( long iy=bgn_y; iy<end_y; iy++ )
            {
                ns_action::take_drvt<'x','N'> ( ix , iy , ns_meshes::mesh_Vx , ns_meshes::mesh_P );
                ns_action::update_soln ( ix , iy , ns_meshes::mesh_Vx );
            }
        }
        

        // update soln_Vy
        {
            ns_action::copy_to_halo<'y'> ( ns_meshes::mesh_P );

            using namespace ns_config;
            int bgn_x = map_bgn.at('x').at('M');    int bgn_y = map_bgn.at('y').at('N');
            int end_x = map_end.at('x').at('M');    int end_y = map_end.at('y').at('N');

            for ( long ix=bgn_x; ix<end_x; ix++ )
            for ( long iy=bgn_y; iy<end_y; iy++ )
            {
                ns_action::take_drvt<'y','N'> ( ix , iy , ns_meshes::mesh_Vy , ns_meshes::mesh_P );
                ns_action::update_soln ( ix , iy , ns_meshes::mesh_Vy );
            }
        }


        // update soln_P
        {
            ns_action::copy_to_halo<'x'> ( ns_meshes::mesh_Vx );
            ns_action::copy_to_halo<'y'> ( ns_meshes::mesh_Vy );

            using namespace ns_config;
            int bgn_x = map_bgn.at('x').at('M');    int bgn_y = map_bgn.at('y').at('M');
            int end_x = map_end.at('x').at('M');    int end_y = map_end.at('y').at('M');

            for ( long ix=bgn_x; ix<end_x; ix++ )
            for ( long iy=bgn_y; iy<end_y; iy++ )
            {
                ns_action::take_drvt <'x','M'> ( ix , iy , ns_meshes::mesh_P , ns_meshes::mesh_Vx );
                ns_action::take_drvt <'y','M'> ( ix , iy , ns_meshes::mesh_P , ns_meshes::mesh_Vy );
                ns_action::update_soln ( ix , iy , ns_meshes::mesh_P );
            }
            // [2023/11/18]
            // NOTE: For performance, it may actually be better if we separate the above loop 
            //       structure into two - one for 'x' dir drvt, one for 'y' dir drvt, so that 
            //       we may re-arrange the loop orders ("ix" and "iy") accordingly. 


            // apply source
            if ( ns_meshes::mesh_P.bool_src == true )
            {
                using namespace ns_meshes;

                double t, A;
                t = (it+1-1./2.)*dt - time_delay;          // No kidding, 1/2 would give zero
                A = (1-2*M_PI*M_PI*central_f*central_f*t*t)*exp(-M_PI*M_PI*central_f*central_f*t*t);

                mesh_P.hst_soln(mesh_P.I_SRC) += A * dt / (dx * dx);
            }
        }


        // record solution
        {
            for ( auto name : { "Vx" , "Vy" , "P" } )
            {
                class_mesh & mesh = * ns_meshes::name2mesh.at(name);
                if ( mesh.bool_rcv ) { mesh.hst_R_signal( it ) = mesh.hst_soln( mesh.I_RCV ); }
            }
        }


        // record energy
        {
            using namespace ns_config;
            
            for ( auto name : { "Vx" , "Vy" , "P" } )
            {
                class_mesh & mesh = * ns_meshes::name2mesh.at(name);
                memory_accessor<ns_prcs::precision,2> & soln = mesh.hst_soln;
                memory_accessor<ns_prcs::precision,2> & prmt = mesh.hst_prmt;
                memory_accessor<double,2> & weight = ns_fields::mesh2weight_energy.at(name)->acc;

                mesh.hst_R_energy(it) = 0.;
                for ( long i = 0; i < mesh.hst_soln.length; i++ )
                { 
                    mesh.hst_R_energy(it) += soln(i) * prmt(i) * soln(i) * weight(i);
                }
            }
        }
    } // for (int it=0; it<ns_config::Nt; it++)
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\nloop time : %5.4f\n", toc - tic); }


tic = Wtime ();
    // ---- output the solution and energy
    {
        using namespace ns_config;
        using namespace ns_fields;

        std::string dt_folder;
        {
            std::ostringstream ss; 
            ss << std::scientific << std::setprecision(6) << dt;
            dt_folder = "dt_" + ss.str();
            std::replace( dt_folder.begin(), dt_folder.end(), '.', 'p' );
        }

        std::string str_precision = "Unrecognized type";
        if ( std::is_same_v < ns_prcs::precision , double   > ) { str_precision = "f64"; }
        if ( std::is_same_v < ns_prcs::precision , float    > ) { str_precision = "f32"; }
        if ( std::is_same_v < ns_prcs::precision , _Float16 > ) { str_precision = "f16"; }
        if ( ns_mpi::rank==0 ) { printf( "\nPrecision type: %s .\n", str_precision.c_str() ); }

        std::string folder_name = "../result/acoustic/" + medium_name 
                                + "/src_" + std::to_string(src_loc_i_base_ppw) 
                                +     "_" + std::to_string(src_loc_j_base_ppw)
                                + "_rcv_" + std::to_string(rcv_loc_i_base_ppw) 
                                +     "_" + std::to_string(rcv_loc_j_base_ppw)
                                + "/base_size" + "_" + std::to_string(N_wavelength_x)
                                +              + "_" + std::to_string(N_wavelength_y)
                                + "/base_ppw_" + std::to_string(base_ppw)
                                + "/multiple_" + std::to_string(multiple)
                                + "/Nt_" + std::to_string(Nt)                                
                                + "/" + dt_folder
                                + "/" + str_precision
                                + "/" + "cpst_" + ns_config::char_cpst
                                + "/" + std::to_string(ns_mpi::size) 
                                      + "_" + std::to_string(ns_mpi::map_dim.at('x'))
                                      + "_" + std::to_string(ns_mpi::map_dim.at('y'));

        if ( ns_mpi::rank == 0 )
        {
            std::filesystem::create_directories( folder_name );
            std::cout << "Results will be stored in folder: \n" << folder_name << std::endl;
            // [2023/09/20] NOTE: create_directories () can take both char* and std::string.
        }
        

        for ( const std::string name : {"Vx", "Vy", "P"} )
        {
            if ( ns_meshes::name2mesh.at(name)->bool_rcv )
            {
                std::string file_name = folder_name + "/" + name; 
                FILE * fp = fopen( ( file_name + ".txt" ).c_str(), "w" );
                for ( int it = 0; it < Nt; it++ ) 
                    { fprintf( fp, "% 16.15e\n", (double) ns_meshes::name2mesh.at(name)->hst_R_signal.at(it) ); }
                fclose(fp);
            }
        }


        {
            using namespace ns_config;
            using namespace ns_meshes;
            memory_resource<double> res_E (Nt);
            memory_accessor<double> & E = res_E.acc;
            for ( int it=1; it<Nt; it++ )
            { 
                E(it) = mesh_Vx.hst_R_energy(it)      + mesh_Vy.hst_R_energy(it) 
                      + mesh_P. hst_R_energy(it) / 2. + mesh_P .hst_R_energy(it-1) / 2.; 
            }
            E(0) = mesh_Vx.hst_R_energy(0) + mesh_Vy.hst_R_energy(0) + mesh_P.hst_R_energy(0) / 2. + 0. / 2.;

            if ( ns_mpi::rank ==0 )
            {
                std::string file_name = folder_name + "/E";
                FILE * fp = fopen( ( file_name + ".txt" ).c_str(), "w" );
                for ( int it = 0; it < Nt; it++ )
                    { fprintf( fp, "% 16.15e\n", (double) E(it) ); }
                fclose(fp);
            }
        }

    }
toc = Wtime (); if ( ns_mpi::rank==0 ) 
{ printf("\n write output time : %5.4f\n", toc - tic); }

    return 0;
}