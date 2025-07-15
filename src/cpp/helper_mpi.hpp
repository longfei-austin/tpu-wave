#ifndef HELPER_MPI_H
#define HELPER_MPI_H

#include <stdarg.h>  // for vprintf
#include <unistd.h>  // for usleep (unit: micro second)

#include "mpi.h"

// [2023/05/04]
// NOTE: The type MPI_Comm is actually int. Keep that in mind when writing
//       the interface. For example, it is alright to use = assignment and 
//       not have a dtor.

class class_mpi_printf
{
    public:
        int rank;
        int size;

        MPI_Comm comm;

        class_mpi_printf () { attach_comm ( MPI_COMM_WORLD ); }

        class_mpi_printf ( MPI_Comm input_comm )
            { attach_comm ( input_comm ); }
        
        void attach_comm ( MPI_Comm input_comm )
        {
            this->comm = input_comm;

            if ( this->comm != MPI_COMM_NULL )
            {
                MPI_Comm_rank( comm , & rank );
                MPI_Comm_size( comm , & size );
            }
        }

        long unsigned u_secs = 1<<5;

        void operator() (const char * format, ...)
        {
            if ( this->comm != MPI_COMM_NULL )
            {
                setbuf(stdout, NULL);  // supposed to disable buffering on stdout

                for ( int i_rank = 0; i_rank < size; i_rank++ )
                {
                    fflush(stdout); MPI_Barrier(comm);
                    if ( i_rank == rank ) 
                    {
                        va_list args;
                        va_start (args, format);
                        vprintf (format, args);
                        va_end (args);        
                    }
                    fflush(stdout); usleep(u_secs); MPI_Barrier(comm);
                }
                if ( rank == size - 1 ) { printf("\n"); }
                fflush(stdout); usleep(u_secs); MPI_Barrier(comm);
            }
        }
};

#endif