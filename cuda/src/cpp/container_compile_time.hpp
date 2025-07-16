#ifndef CONTAINER_COMPILE_TIME_H
#define CONTAINER_COMPILE_TIME_H

// #include <cstdlib>
#include <iostream>
#include <fstream>

// NOTE: PAY ATTENTION TO the empty initialization : array() in between 
//       the constructor declaration and the constructor body. 
//           Without this empty initialization, gcc compiles, clang and 
//       icc complain that the constructor does not produce constant 
//       expression. 
//           Not sure what is the precise language rule here. But this 
//       could make a good blog entry.

template<typename T, int L>
class compile_time_vector
{
    public:
        static constexpr int length = L;
        T array[L];
        // typedef T Array[L];    Array array;

        constexpr compile_time_vector( ) : array() 
            { for ( int i=0; i<L; i++ ) array[i] = static_cast<T>(0); }

        constexpr compile_time_vector( T const (&init_array)[L] ) : array()
            { for ( int i=0; i<L; i++ ) array[i] = init_array[i]; }

        // Question: should we make the following a constructor or a normal function?
        constexpr compile_time_vector( T const (&init_array)[L] , T const &a ) : array()
            { for ( int i=0; i<L; i++ ) array[i] = init_array[i] * a ; }

        constexpr compile_time_vector( compile_time_vector<T,L> const &init_compile_time_vector ) : array()
            { for ( int i=0; i<L; i++ ) array[i] = init_compile_time_vector.array[i]; }

        constexpr T operator[] (int i) const
            { return array[i]; }

        constexpr T operator() (int i) const
            { return array[i]; }

        constexpr T at (int const &i) const
        {
            if (i<0 || i>=L) printf("compile_time_vector.at(): entry needs to be within range." ); 
            return array[i]; 
        }

        constexpr void operator*= (T const &a)
            { for ( int i=0; i<L; i++ ) array[i] *= a; }

        constexpr void operator/= (T const &a)
            { for ( int i=0; i<L; i++ ) array[i] /= a; }

        constexpr void operator= ( compile_time_vector<T,L> const &init_compile_time_vector )
            { for ( int i=0; i<L; i++ ) array[i] = init_compile_time_vector.array[i]; }

        constexpr void operator= ( T const (&init_array)[L] )
            { for ( int i=0; i<L; i++ ) array[i] = init_array[i]; }
};


template<typename T, int M, int N>
class compile_time_matrix
{
    public:
        static constexpr int rows = M;
        static constexpr int cols = N;
        T array[M][N];

        constexpr compile_time_matrix( ) : array() 
        { 
            for ( int i=0; i<M; i++ ) 
                for ( int j=0; j<N; j++ ) { array[i][j] = static_cast<T>(0); }
        }

        constexpr compile_time_matrix( T const (&init_array)[M][N] ) : array()
        { 
            for ( int i=0; i<M; i++ ) 
                for ( int j=0; j<N; j++ ) { array[i][j] = init_array[i][j]; }
        }

        constexpr compile_time_matrix( T const (&init_array)[M][N] , T const &a ) : array()
        { 
            for ( int i=0; i<M; i++ ) 
                for ( int j=0; j<N; j++ ) { array[i][j] = init_array[i][j] * a; }
        }

        constexpr compile_time_matrix( compile_time_matrix<T,M,N> const &init_compile_time_matrix ) : array()
        { 
            for ( int i=0; i<M; i++ ) 
                for ( int j=0; j<N; j++ ) { array[i][j] = init_compile_time_matrix.array[i][j]; }
        }

        constexpr void operator*= (T const &a)
        { 
            for ( int i=0; i<M; i++ ) 
                for ( int j=0; j<N; j++ ) { array[i][j] *= a; }
        }

        constexpr void operator/= (T const &a)
        { 
            for ( int i=0; i<M; i++ ) 
                for ( int j=0; j<N; j++ ) { array[i][j] /= a; }
        }

        constexpr T operator() (int i , int j) const 
            { return array[i][j]; }

        constexpr T at (int &i , int &j) const
        { 
            if (i<0 || i>=M || j <0 || j>=N) printf("compile_time_matrix.at(): entry needs to be within range." ); 
            return array[i][j]; 
        }
};


#endif