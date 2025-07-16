#ifndef LINEAR_MEMORY_H
#define LINEAR_MEMORY_H

#include <iostream>
#include <vector>

#include <cassert>

#include <math.h>
#include <string.h>

template<typename T , int dim=1>
class memory_accessor
{
    public:
        long length = 0;
        T * ptr = nullptr;

        std::array<long,dim> dir_length {{-1}};
        std::array<long,dim> dir_stride {{-1}};        


        void attach_dimension ( std::array<long,dim> const & A_len ) 
        {
            this->dir_length = A_len;

            long deduced_length = 1;
            for ( int i = 0; i < dim; i++ ) { deduced_length *= dir_length.at(i); }
            assert ( length == deduced_length && "Error: length does not match" );

            for ( int i = 0; i < dim; i++ )
            {
                dir_stride.at(i) = 1;
                for ( int j = i + 1; j < dim; j++ ) 
                { 
                    dir_stride.at(i) *= dir_length.at(j); 
                }
            }
        }


    private:
        memory_accessor<T,dim> &
        private_attach_memory (T * p , long const L)
        {
            assert ( ptr == nullptr );  // NOTE: This check is not necessary, but may help eliminate some logic bugs.
            assert ( p   != nullptr );
            assert ( L > 0 );
            assert ( L * sizeof(T) < 1024.*1024.*1024. * 80. 
                     && "requested larger than 80GB; if intended, change the limit in this assert" );

            this->ptr = p;  this->length = L;

            return *this;
        }


    public:
        memory_accessor<T,dim> &
        attach_memory (T * p , std::array<long,dim> const & A_len)
        {
            long L = 1;  for ( unsigned int i = 0; i < A_len.size(); i++ ) { L *= A_len.at(i); }

            private_attach_memory (p , L);
            attach_dimension ( A_len );

            return *this; 
        }


        memory_accessor<T,dim> &
        attach_memory (T * p , long const L)
        { 
            private_attach_memory (p , L);
            attach_dimension ( {L} );

            return *this; 
        }


        memory_accessor<T,dim> &
        attach_memory (memory_accessor<T,dim> const & A)
            { attach_memory (A.ptr , A.dir_length);  return *this; }
        // [2023/12/13]
        // NOTE: If wanting attach_dimension to be called, need to pass A.dir_length
        //       instead of A.length.


        memory_accessor(){}  // Empty constructor in case we need to declare a memory_accessor before knowing its 
                             // size. Though with an empty body, it is still an user-defined default constructor. 
                             // Sometimes, its implicit generation may be inhibited for various reasons.
                             //     ( reference: https://en.cppreference.com/w/cpp/language/default_constructor )

        memory_accessor (T * p, std::array<long,dim> const & A_len)
            { attach_memory (p , A_len); }

        memory_accessor (T * p , long const L)
            { attach_memory (p , { L }); }


        // assignment operators
        void operator= ( std::initializer_list<T> il ) 
        { 
            if ( static_cast<long unsigned int>(length) != il.size() ) 
                { printf("Error: length - operator= - memory_accessor.\n"); fflush(stdout); exit(0); }

            long i = 0; for ( auto& e : il ) { this->at(i) = e; i++; }
        }

        void operator*= (T const &a) { for ( long i=0; i<length; i++ ) { ptr[i] *= a; } }
                                                       // { this->operator()(i) *= a; }

        void operator/= (T const &a) { for ( long i=0; i<length; i++ ) { ptr[i] /= a; } }
                                                       // { this->operator()(i) /= a; }

        T* data() { return ptr; }

        T& operator() (long const i) // read and write access
        { 
            assert ( ptr != nullptr && "Error: nullptr - operator() - memory_accessor." );
            assert ( i>=0 && i<length && "Out of bound access - memory_accessor." );

            return ptr[i]; 
        }


        T& at (long const i) // read and write access
        {
            assert ( ptr != nullptr && "Error: nullptr - operator() - memory_accessor." );
            assert ( i>=0 && i<length && "Out of bound access - memory_accessor." );

            return ptr[i];
        }


        T& at (long const i) const // read access
        {
            assert ( ptr != nullptr && "Error: nullptr - operator() - memory_accessor." );
            assert ( i>=0 && i<length && "Out of bound access - memory_accessor." );

            return ptr[i];
        }        


        long convert_index ( std::array<long,dim> const & A_ind ) 
        {
            long ind = 0;
            for ( int i = 0; i < dim; i++ ) { ind += A_ind.at(i) * dir_stride.at(i); }

            for ( int i = 0; i < dim; i++ ) 
                { assert( A_ind.at(i) >= 0 && A_ind.at(i) < dir_length.at(i) ); }
            assert ( ( ind >= 0 && ind < length ) && "Out of bound access - - memory_accessor" );

            return ind;
        }


        long convert_index ( long const & i0, long const & i1 ) const
        {
            static_assert ( dim == 2 , "This function is intended for a 2D array." );

            long ind = i0 * dir_stride[0]
                     + i1 * dir_stride[1];

            return ind;
        }
        // [2023/12/13]
        // NOTE: When debugging, we can wrap the two indices inside an {} so that the 
        //       convert_index ( A_ind ) version is called, which has bound checking.
        //       
        // [2023/12/14]
        // NOTE: Added const so that we may use it inside SYCL kernel. We may need a 
        //       non-const version somewhere down the road.


        long convert_index ( long const & i0, long const & i1, long const & i2 ) 
        {
            static_assert ( dim == 3, "This function is intended for a 3D array." );

            long ind = i0 * dir_stride[0]
                     + i1 * dir_stride[1]
                     + i2 * dir_stride[2];

            return ind;
        }


        T& at (std::array<long,dim> const & A_ind) // read and write access
        {
            assert ( ptr != nullptr && "Error: nullptr - at - memory_accessor" );

            long ind = convert_index (A_ind);            
            return ptr[ind];
        }


        // ---- an overloaded .at() that lets us retrieve the entry by two "long"s to 
        //      "mimic" the calling convention for matrices. /* dim need to be 2 */
        //
        T& at ( long const & i0, long const & i1 ) const
        {
            assert ( ptr != nullptr && "Error: nullptr - at - memory_accessor" );

            long ind = convert_index(i0,i1);
            return ptr[ind];
        }
        // [2023/12/14]
        // NOTE: Added const so that we may use it inside SYCL kernel. We may need a 
        //       non-const version somewhere down the road.


        // ---- an overloaded .at() that lets us retrieve the entry by three "long"s to 
        //      "mimic" the calling convention for a 3D tensor. /* dim need to be 3 */
        //
        T& at ( long const & i0, long const & i1, long const & i2 ) // read and write access
        {
            assert ( ptr != nullptr && "Error: nullptr - at - memory_accessor" );

            long ind = convert_index(i0,i1,i2);
            return ptr[ind];
        }


        void set_constant(T const & a) { for ( long i = 0; i < length; i++ ) { ptr[i] = a; } }


        T L2_norm () const
            { return sqrt( L2_norm_square () ); }  // Question: will sqrt from math.h work well with template?

        T L2_norm_square () const
        {
            T norm_square = static_cast<T>(0);
            for ( long i=0; i<length; i++ ) 
                { norm_square += this->ptr[i] * this->ptr[i]; }

            return norm_square;
        }

        void print ()
        {   
            for ( int ix = 0; ix < dir_length.at(0); ix++ )
            {
                for ( int iy = 0; iy < dir_length.at(1); iy++ )
                {
                    printf("%4.0f ", this->at(ix,iy));
                }
                printf("\n");
            }
            printf("\n");
        }

        void print_for_copy ( long const bgn , long const end ) const    // NOTE: "end" is NOT inclusive
        {   
            printf( "length : %ld \n" , end - bgn );
            printf( "{ { " );
            for ( long i=bgn; i<end; i++ ) 
            { 
                printf("% 16.15e" , static_cast< double > ( this->ptr[i] ) ); 
                if ( i < end-1 ) { printf(" , "); }
            }
            printf( " } };\n\n" );
        }

        T* bgn () { return ptr;          }
        T* end () { return ptr + length; } // one past the last
}; 
//class memory_accessor


template<typename T , int dim=1>
class memory_resource
{
    public:
        long length = 0;
        T * ptr = nullptr;

        std::string name = "undefined";
        // [2022/12/27]
        // NOTE: Give each memory_resource an std::string to store an optional name.
        //       This can be helpful for debugging or for ensuring that the intended 
        //       memory_resource is used.

        memory_accessor<T,dim> acc;


        memory_resource(){}  // Empty constructor in case we need to declare a memory_resource before knowing its 
                             // size.
                             //
                             // NOTE: Though with an empty body, it is still an user-defined default constructor. 
                             //       Sometimes, its implicit generation may be inhibited for various reasons.
                             //       ( reference: https://en.cppreference.com/w/cpp/language/default_constructor )

    private:
        void private_allocate (long const L)
        {
            assert ( ptr == nullptr && "ptr is already allocated" );
            assert ( L > 0 );
            assert ( (double) ( L * sizeof(T) ) < 1024.*1024.*1024. * 80. // use floating point to avoid integer overflow
                     && "requested size larger than 80GB; if intended, change the limit in the assert" );

            constexpr long alignment = 32;    // 32 bytes; 
            long byte_size = ( ( L*sizeof(T) + alignment - 1 ) / alignment ) * alignment;

            ptr = static_cast<T *>( aligned_alloc(alignment, byte_size) );
            assert ( ptr != nullptr && "aligned_alloc returned nullptr" ); 
            memset ( ptr, 0, byte_size );

            this->length = L;
            // NOTE: Assignment of L to this->length is here because sometimes we push_back an empty 
            //       memory_resource to std::vector and call allocate_memory explicitly afterwards.
            //           In this case, it would be "easy to forget" setting this->length accordingly
            //       if we don't include it here.
        }
        // NOTE: If empty constructor (currently the default constructor) is called, the above
        //       function needs to be called explicitly to allocate memory.


    public:
        void make_accessor (std::array<long,dim> const & A_len)
            { acc.attach_memory (this->ptr , A_len); }


        void allocate_memory (std::array<long,dim> const & A_len)
        {
            long L = 1;  for ( unsigned int i = 0; i < A_len.size(); i++ ) { L *= A_len.at(i); }
            private_allocate (L);

            make_accessor (A_len);
        }

        void allocate_memory (long const L)
            { private_allocate (L); make_accessor ( {L} ); }


        memory_resource (std::array<long,dim> const & A_len) // a user-supplied "non-speical" constructor
            { allocate_memory (A_len); }

        memory_resource (long const L)
            { allocate_memory ( {L} ); }


        void copy ( memory_resource<T> const & v )
        {
            assert( ( this != &v && this->ptr != v.ptr ) && "No self copy" );
            assert( this->ptr != nullptr && v.ptr != nullptr );
            assert( this->length > 0 && this->length == v.length );
            assert ( v.acc.ptr != nullptr && "input accessor not valid" );

            // NOTE: Syntactically, this copy member function should "work" for self assignment. Its
            //       behavior is (arguably) "logically" correct. However, such "self" assignment is 
            //       an indication that some error or misuse may have occurred from the "application"
            //       side. We add this (self copy) test to guard such misuse.

            for ( long i=0; i<length; i++ ) { ptr[i] = v.ptr[i]; } 
            this->acc = v.acc;
        }


        // copy ctor
        memory_resource( memory_resource<T> const & v ) : memory_resource ( v.acc.dir_length ) // : delegating to user-supplied ctor
            { copy (v); }

        // copy assignment operator
        memory_resource<T> & operator=( memory_resource<T> const & v )
            { copy (v); return *this; }
        // [2023/04/04] NOTE: Do we actually want to return void?


        // move ctor
        memory_resource( memory_resource<T> && v ) noexcept
        {
            this->length = v.length;  this->ptr = v.ptr;  this->acc = v.acc;

            v.length = 0;             v.ptr = nullptr;
        }


        // another move assignment operator
        void operator=( memory_resource<T> && v ) noexcept
        {
            // [2023/04/04]
            // NOTE: move assignment is only allowed for empty container; otherwise we 
            //       would need to destroy resource during the move assignment, which
            //       I feel should be discouraged.

            assert( ptr == nullptr && "move assignment is only allowed when moved-to vector is empty" );
            assert( this != &v && this->ptr != v.ptr );

            this->length = v.length;  this->ptr = v.ptr;  this->acc = v.acc;
            v.length = 0;             v.ptr = nullptr;  

            // [2023/04/04]
            // NOTE: Reiterate: the above implementation is only proper if "this" 
            //       memory_resource does not currently hold data.
        }


        // dtor
        ~memory_resource() { free(ptr); }
        // Questions: Is free(ptr) enough? What might be the corner cases that I haven't thought of?
}; 
// class memory_resource


#endif