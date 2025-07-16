// #ifndef NAMESPACE_MPI_H
// #define NAMESPACE_MPI_H

// // [2022/12/25]
// // NOTE: Be AWARE that the file name the header guard are common enough that they 
// //       may be used by third party libraries.


// #include <vector>
// #include "src_rcv.hpp"
// #include "mpi.h"


// namespace namespace_mpi
// {
//     inline int mpi_rank  = -(1<<30);
//     // inline int num_procs = -(1<<30);

//     // inline int num_cmpt_procs = -(1<<30);

//     // [2023/03/07]
//     // NOTE: Should consider changing mpi_rank to mpi_rank_WORLD or mpi_rank_WL.

//     inline bool bool_proc_C = false;
//     inline bool bool_proc_M = false;


//     inline MPI_Comm comm_ND = MPI_COMM_NULL;

//     inline int  mpi_rank_ND  = -(1<<30);
//     inline int  num_procs_ND = -(1<<30);

//     inline int  num_cmpt_procs_ND = -(1<<30);


//     inline MPI_Comm comm_MG = MPI_COMM_NULL;

//     inline int  node_index = -(1<<30);
//     inline int  num_nodes  = -(1<<30);

//     // ---- simple functions for collecting elapsed times
//     inline double MPI_t_bgn;
//     inline double MPI_t_end;

//     inline void tic() { MPI_t_bgn = MPI_Wtime(); }

//     inline void toc(std::string message) 
//     {
//         MPI_t_end = MPI_Wtime();
//         printf("\n Elapsed time - %s : %16.15f. \n", message.c_str(), MPI_t_end - MPI_t_bgn );
//     }

//     inline void toc(double & accumulated_duration) 
//     {
//         MPI_t_end = MPI_Wtime();
//         accumulated_duration += MPI_t_end - MPI_t_bgn;
//     }
// }
// namespace ns_mpi = namespace_mpi;
// // 
// // [2022/12/25]
// // NOTE: inline varialbe is a c++-17 feature; it allows multiple definition 
// //       in multiple translation units; these variables have external linkage 
// //       by default; the linker will pick one definition and neglect the rest. 
// //           An alternative is to DECLARE extern variables in the header:
// //       
// //                    namespace namespace_mpi
// //                    {
// //                        extern int mpi_rank;
// //                    }
// //       and in exact one .cpp file, DEFINE it:
// //                    namespace namespace_mpi
// //                    {
// //                        int mpi_rank = -(1<<30);
// //                    }
// //       /* Initialize extern variable in header file will define it, will 
// //          still lead to violation of ODR rule. */
// // 
// //       In both the cases (inline and extern), we are using global variables.

// #endif