#include "forward.hpp"

using ns_input::Nt;

using ns_input::dx;
using ns_input::dt;

using ns_input::central_f;
using ns_input::time_delay;

using ns_input::bool_energy;


// [2023/06/19]
// NOTE: Didn't want to make this file tpp because of the using keywords above 
//       (e.g., using ns_input::Nt) - didn't want them to appear in the header. 
//       Instead, we'll use explicit template instantiation below.
// [2023/07/13]
// NOTE: It started to feel like there were lots of explicit instantiations.
template<int cpst_N , char cpst_S>
void Class_Forward_Specs::forward_simulation_periodic_y ()
{
    Class_Grid & Grid_SNN = * ( Map_Grid_pointers.at( {'N','N'} ) );
    Class_Grid & Grid_SMM = * ( Map_Grid_pointers.at( {'M','M'} ) );

    Class_Grid & Grid_SMN = * ( Map_Grid_pointers.at( {'M','N'} ) );
    Class_Grid & Grid_SNM = * ( Map_Grid_pointers.at( {'N','M'} ) );


    for (int it=0; it<Nt; it++) 
    {

        if ( bool_energy ) { ns_input::Record_E_p0.at(it) += Grid_SMM.energy_calculation(); }
        if ( bool_energy ) { ns_input::Record_E_p0.at(it) += Grid_SNN.energy_calculation(); }


        // ---- Vx
        Grid_SMN.calculate_derivative ( 'x' );
        Grid_SMN.secure_bdry ( 'x' );

        // Grid_SMN.calculate_derivative_periodic_modulo ( 'y' );
        Grid_SMN.calculate_derivative_periodic ( 'y' );

        Grid_SMN.update_cpst<cpst_N,cpst_S> ();


        // ---- Vy
        Grid_SNM.calculate_derivative ( 'x' );
        Grid_SNM.secure_bdry ( 'x' );

        // Grid_SNM.calculate_derivative_periodic_modulo ( 'y' );
        Grid_SNM.calculate_derivative_periodic ( 'y' );

        Grid_SNM.update_cpst<cpst_N,cpst_S> ();


        // Apply source
        if ( src_forward.is_src ) 
        {   
            if ( src_forward.c_source_type == 'V' ) // if source is on a velocity grid
            {
                double t, A;
                t = it*(double) dt - time_delay;
                A = (1-2*M_PI*M_PI*central_f*central_f*t*t)*exp(-M_PI*M_PI*central_f*central_f*t*t);

                Class_Grid & grid_src = * ( Map_Grid_pointers.at( src_forward.grid_type ) );

                for ( int i_field = 0; i_field < grid_src.N_soln; i_field++ )
                    { grid_src.Vec_soln.at(i_field).at(src_forward.I_SRC) 
                    = (double) grid_src.Vec_soln.at(i_field).at(src_forward.I_SRC) + A * (double) dt / ((double) dx * (double) dx); }
            }
        }


        if ( bool_energy ) { ns_input::Record_E_k.at(it) += Grid_SMN.energy_calculation (); }
        if ( bool_energy ) { ns_input::Record_E_k.at(it) += Grid_SNM.energy_calculation (); }


        // ---- Sxy
        Grid_SNN.calculate_derivative ( 'x' );
        Grid_SNN.secure_bdry ( 'x' );

        // Grid_SNN.calculate_derivative_periodic_modulo ( 'y' );
        Grid_SNN.calculate_derivative_periodic ( 'y' );

        Grid_SNN.update_cpst<cpst_N,cpst_S> ();


        // ---- {Sxx; Syy; Szz}
        Grid_SMM.calculate_derivative<true> ( 'x' );
        Grid_SMM.secure_bdry<true> ( 'x' );

        // Grid_SMM.calculate_derivative_periodic_modulo ( 'y' );
        Grid_SMM.calculate_derivative_periodic<true> ( 'y' );

        Grid_SMM.update_cpst<cpst_N,cpst_S> ();


        // Apply source
        if ( src_forward.is_src ) 
        {   
            if ( src_forward.c_source_type == 'S' ) // if source is on a stress grid
            {
                double t, A;
                t = (it+1./2.)*(double) dt - time_delay;       // No kidding, 1/2 must have given zero
                A = (1-2*M_PI*M_PI*central_f*central_f*t*t)*exp(-M_PI*M_PI*central_f*central_f*t*t);

                Class_Grid & grid_src = * ( Map_Grid_pointers.at( src_forward.grid_type ) );

                for ( int i_field = 0; i_field < grid_src.N_soln; i_field++ )
                    { grid_src.Vec_soln.at(i_field).at(src_forward.I_SRC) 
                    = (double) grid_src.Vec_soln.at(i_field).at(src_forward.I_SRC) + A * (double) dt / ((double) dx * (double) dx); }
            }
        }

        // Warning: Because of the strongly imposed periodic boundary condition, the first and last 
        //          grid points on y direction are supposed to be duplicates of each other. Therefore, 
        //          when applying the source, if the source happens to land on the periodic boundary, 
        //          it needs to be applied on both ends as below.
        // 
        //              grid_src.Vec_soln.at(i_field).at(src_forward.I_SRC + 80) += A * dt / (dx * dx);
        // 
        //              However, because of the special index mapping I am using, the last grid point
        //          was mapped to the first one when being used to calculate derivatives. Therefore, 
        //          although the results on the last grid point is incorrect, the simulation results
        //          elsewhere are not affected.
        //
        //              When calculating energy, if we adjust the energy parameters corresponding to 
        //          the first and last grid points using factors 1./2. and 1./2., the result would be 
        //          incorrect; same for factors 0. and 1.; but if using 1. and 0. instead, the result 
        //          would be correct.
        //              
        //              We should "forbid" the source to be placed on the last grid and use factors 1. 
        //          and 0. when calculating energy. This way, we don't need to test if the source is on 
        //          the strong boundary.


        if ( bool_energy ) { ns_input::Record_E_p1.at(it) += Grid_SMM.energy_calculation(); }
        if ( bool_energy ) { ns_input::Record_E_p1.at(it) += Grid_SNN.energy_calculation(); }


        // Store solution RESULT
        for ( const auto & iter_map : Map_grid_N_rcvs ) // go through the possible grid types
        {
            auto & grid_type = iter_map.first;

            if ( Map_grid_N_rcvs.at(grid_type) > 0 )     // if number of receiver on this grid is larger than 0
            {
                Class_Grid & grid_rcv = * ( Map_Grid_pointers.at(grid_type) );

                auto & Vec_struct_rcv = Map_grid_struct_rcv_forward.at(grid_type);
                auto & Vec_RESULT_rcv = Map_grid_RESULT_rcv.at(grid_type);

                for ( int i_rcv = 0; i_rcv < Map_grid_N_rcvs.at(grid_type); i_rcv++ )  // loop through the receivers
                {
                    int I_RCV              = Vec_struct_rcv.at(i_rcv).I_RCV;
                    auto & Solution_RESULT = Vec_RESULT_rcv.at(i_rcv);

                    for ( int i_field = 0; i_field < grid_rcv.N_soln; i_field++ )    // loop through the fields on this grid
                        { Solution_RESULT(i_field,it) = grid_rcv.Vec_soln.at(i_field).at(I_RCV); }
                }
            }
        }

        if ( (it % 1000) == 0 ) 
        { 
            printf( "\n" ); fflush(stdout); 
            Class_Grid & grid_src = * ( Map_Grid_pointers.at( src_forward.grid_type ) );
            for ( int i_field = 0; i_field < grid_src.N_soln; i_field++ )
            {
                printf( "Time step: %6d soln: % 12.11e \n", it, (double) grid_src.Vec_soln.at(i_field).at(src_forward.I_SRC) );
            }
        }
        // [2024/03/27]
        // NOTE: Made the above changes on the printing out the solution at the src location. Before, we always print out
        //       Sxx and Syy - this can lead to mismatch in nvcc_host_FDM_2D.exe and nvcc_device_FDM_2D.exe when the src
        //       is not applied on the {Sxx;Syy} grid. The reason is that the device code uses padding.

    }  // for (int it=0; it<Nt; it++) 
    printf( "\n" ); 

} // forward_simulation_periodic_y ()

template void Class_Forward_Specs::forward_simulation_periodic_y < 0,'0'> ();

template void Class_Forward_Specs::forward_simulation_periodic_y <-3,'C'> ();
template void Class_Forward_Specs::forward_simulation_periodic_y <-3,'R'> ();

template void Class_Forward_Specs::forward_simulation_periodic_y < 3,'C'> ();
template void Class_Forward_Specs::forward_simulation_periodic_y < 3,'R'> ();

template void Class_Forward_Specs::forward_simulation_periodic_y < 6,'C'> ();
template void Class_Forward_Specs::forward_simulation_periodic_y < 6,'R'> ();