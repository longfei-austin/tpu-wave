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

    ns_mpi::rank = 0;
    ns_mpi::size = 1;

tic = Wtime ();    
    long len_x = 600;
    long len_y = 600;
    long mem_len = len_x * len_y;
    
    constexpr long alignment = 32;    // 32 bytes; 
    long byte_size = ( ( mem_len*sizeof(ns_prcs::precision) + alignment - 1 ) / alignment ) * alignment;

    ns_prcs::precision * soln_ptr = nullptr;
    ns_prcs::precision * drvt_ptr = nullptr;

    {
        // Allocate memory
        using T = ns_prcs::precision;
        long L = mem_len;

        {
            T * ptr = static_cast<T *>( aligned_alloc(alignment, byte_size) );
            assert ( ptr != nullptr && "aligned_alloc returned nullptr" ); 
            std::iota(ptr, ptr+L, 0);
            soln_ptr = ptr;
        }

        {
            T * ptr = static_cast<T *>( aligned_alloc(alignment, byte_size) );
            assert ( ptr != nullptr && "aligned_alloc returned nullptr" ); 
            memset ( ptr, 0, byte_size );
            drvt_ptr = ptr;
        }
    }

    for (long ix = 0; ix < len_x; ix+=100)
    for (long iy = 0; iy < len_y; iy+=100)
    {
        long ind = ix * len_y + iy;
        printf( "% 4ld % 4ld % 8ld % e %e \n", ix, iy, ind, soln_ptr[ind], drvt_ptr[ind] );
    }
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n Allocate memory time : %5.4f\n", toc - tic); }


tic = Wtime ();
    {
        // Take derivatives (operating on strided data)
        constexpr int stencil_length = 4;
        constexpr std::array<double,stencil_length> stencil = {1./24 , -9./8 , 9./8 , -1./24};
        constexpr std::array<int,stencil_length> stencil_shift = std::array<int,4> {-2,-1,0,1};

        // x-derivative
        for ( long ix=2; ix<len_x-2; ix++ )
        for ( long iy=2; iy<len_y-2; iy++ )
        {
            long i_d = ix * len_y + iy;
            for ( long iw=0; iw<stencil_length; iw++ )
            {   
                long i_s = (ix + stencil_shift[iw]) * len_y + iy;
                drvt_ptr[i_d] += soln_ptr[i_s] * stencil[iw];
            }
        }
        for (long ix = 0; ix < len_x; ix+=100)
        for (long iy = 0; iy < len_y; iy+=100)
        {
            long ind = ix * len_y + iy;
            printf( "% 4ld % 4ld % 8ld % e %e \n", ix, iy, ind, soln_ptr[ind], drvt_ptr[ind] );
        }
        printf( "\n" );
        memset ( drvt_ptr, 0, byte_size );

        // y-derivative
        for ( long ix=2; ix<len_x-2; ix++ )
        for ( long iy=2; iy<len_y-2; iy++ )
        {
            long i_d = ix * len_y + iy;
            for ( long iw=0; iw<stencil_length; iw++ )
            {   
                long i_s = ix * len_y + (iy + stencil_shift[iw]);
                drvt_ptr[i_d] += soln_ptr[i_s] * stencil[iw];
            }
        }
        for (long ix = 0; ix < len_x; ix+=100)
        for (long iy = 0; iy < len_y; iy+=100)
        {
            long ind = ix * len_y + iy;
            printf( "% 4ld % 4ld % 8ld % e %e \n", ix, iy, ind, soln_ptr[ind], drvt_ptr[ind] );
        }
        printf( "\n" );
        memset ( drvt_ptr, 0, byte_size );
    }
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n Take derivatives time : %5.4f\n", toc - tic); }


tic = Wtime ();
    {
        // Apply source (change single entry)
        long ix_src = len_x / 2;
        long iy_src = len_y / 2;
        long I_SRC = ix_src * len_y + iy_src;

        printf("soln[I_SRC] before applying source: % e\n", soln_ptr[I_SRC]);
        for ( long it = 0; it < 1000; it++ ) {soln_ptr[I_SRC] += it; }
        printf("soln[I_SRC] after  applying source: % e\n", soln_ptr[I_SRC]);

        std::iota(soln_ptr, soln_ptr+mem_len, 0);
    }
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n Apply source time : %5.4f\n", toc - tic); }


tic = Wtime ();
    {
        // Record at receiver (copy single entry)
        long ix_src = len_x / 2;
        long iy_src = len_y / 2;
        long I_RCV = ix_src * len_y + iy_src;

        using T = ns_prcs::precision;
        constexpr long Nt = 1000;
        T * signal_ptr = static_cast<T *>( malloc( Nt*sizeof(T) ) );
        memset ( signal_ptr, 0, Nt*sizeof(T) );

        for ( long it = 0; it < 1000; it++ )
        {
            soln_ptr[I_RCV] += 1.;
            signal_ptr [it]  = soln_ptr[I_RCV];
        }

        for ( long it = 0; it < 1000; it+=10 )
            { printf("signal at I_RCV: % e\n", signal_ptr [it]); }

        std::iota(soln_ptr, soln_ptr+mem_len, 0);
    }
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n Record at receiver time : %5.4f\n", toc - tic); }


tic = Wtime ();
    {
        // Operate along the boundaries
        double sum = 0.;
        for (long ind = 0; ind < mem_len; ind++) { sum += soln_ptr[ind]; }

        printf("Sum before modifying boundary: % 16.15e % 16.15e\n", (mem_len-1)*mem_len/2.
                                                                   , sum );  // formula should qual accumulation 

        double c0 = soln_ptr[ 0          * len_y + 0];
        double c1 = soln_ptr[ 0          * len_y + (len_y - 1)];

        double c2 = soln_ptr[(len_x - 1) * len_y + 0];
        double c3 = soln_ptr[(len_x - 1) * len_y + (len_y - 1)];

        for (long ix = 0; ix < len_x; ix+=1)  // xL
            { long iy = 0;           long ind = ix * len_y + iy;  soln_ptr[ind] = 0; }

        for (long ix = 0; ix < len_x; ix+=1)  // xR
            { long iy = len_y -  1;  long ind = ix * len_y + iy;  soln_ptr[ind] = 0; }

        for (long iy = 0; iy < len_y; iy+=1)  // yL
            { long ix = 0;           long ind = ix * len_y + iy;  soln_ptr[ind] = 0; }

        for (long iy = 0; iy < len_y; iy+=1)  // yR
            { long ix = len_x - 1;   long ind = ix * len_y + iy;  soln_ptr[ind] = 0; }        

        sum = 0.;
        for (long ind = 0; ind < mem_len; ind++) { sum += soln_ptr[ind]; }

        printf("Sum after  modifying boundary: % 16.15e % 16.15e\n", (mem_len-1)*mem_len/2. 
                                                                   - (c0 + c1) * len_y / 2 
                                                                   - (c2 + c3) * len_y / 2 
                                                                   - (c0 + c2) * len_x / 2
                                                                   - (c1 + c3) * len_x / 2
                                                                   +  c0 + c1 + c2 + c3 
                                                                   , sum );  // formula should qual accumulation 
        std::iota(soln_ptr, soln_ptr+mem_len, 0);
    }
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n Boundary time : %5.4f\n", toc - tic); }


tic = Wtime ();
    {
        // calculate energy (aggregation)
        double E = 0.;
        for (long ix = 0; ix < len_x; ix+=1)
        for (long iy = 0; iy < len_y; iy+=1)
        {
            long ind = ix * len_y + iy;
            E += soln_ptr[ind] * soln_ptr[ind];
        }

        printf("Energy from accumulation and from formula:: % 16.15e % 16.15e\n", E, (mem_len-1)*mem_len*(2*mem_len-1)/6.);
        std::iota(soln_ptr, soln_ptr+mem_len, 0);
    }
toc = Wtime (); 
if ( ns_mpi::rank==0 ) { printf("\n Calculate energy time : %5.4f\n", toc - tic); }


    free(soln_ptr);
    free(drvt_ptr);

    return 0;
}