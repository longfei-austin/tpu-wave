#include "grid.cuh"

__global__ void cuda_apply_source ( ns_type::cuda_precision * S , int ind , 
                                    ns_type::cuda_precision * R , int i_field ,
                                    double increment )
{
    R [ i_field ] += increment;
    device_slow2sum <ns_type::cuda_precision> ( S[ind] , R[i_field] , S[ind] , R[i_field] );
}

__global__ void cuda_apply_source ( ns_type::cuda_precision * S , int ind , double increment )
{
    S [ ind ] += increment;
}

__global__ void cuda_record_soln ( ns_type::cuda_precision * R , ns_type::cuda_precision * S , int I_RCV )
{
    * R = S [ I_RCV ];
}

__global__ void cuda_print_soln ( ns_type::cuda_precision * S , int ind , int it )
{
    printf( "Time step: %6d soln: % 12.11e \n", it, (double) S [ind] );
}


// [2023/07/13]
// NOTE: We can test the idea of using a global pointer to it using cuda_print_soln.
// [2023/07/18]
// NOTE: Actually, "it" is on host, we may still need to pass it in to the kernel
//       (__global__ function above). However, whichever host code (could be the 
//       main function) that launches the above kernel can have a (host) pointer 
//       pointing to "it". So we need at most pass "it" once.