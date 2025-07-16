for fp in {64,32,16,}; do

	exe_name=cpu_acoustic
	prefix=cpu_acoustic
	postfix=

	make clean; 
    make ${exe_name} MACROFLAGS="-DPRCSFLAG=${fp}"; 			 
    mv ${exe_name}.exe ${exe_name}\_${fp}.exe;
	

	for c in {0,3,6,}; do
		for m in {$(seq 1 2 1),}; do
			for Nt in {60000,}; do
				for ppn in {12,}; do
					list=; 
					if [ $(($ppn%2)) -eq 1 ]
					then
						for l in {$(seq 1 1 ${ppn}),}; do list=${list}:1-51; done;
					else
						ppn_h=$(($ppn/2))
						for l in {$(seq 1 1 $ppn_h),}; do list=${list}:1-51; done;
						for l in {$(seq $(($ppn_h+1)) 1 $ppn),}; do list=${list}:53-103; done;
					fi

					job_file=job_fp_${fp}\_cpst_${c}\_multiple_${m}\_Nt_${Nt}_ppn_${ppn}\.pbs

					num_nodes=2
					np=$((${num_nodes} * ${ppn}))

					echo "
					#!/bin/bash -v
					#PBS -l select=${num_nodes}
					#PBS -A Performance
					#PBS -q prod
					#PBS -l walltime=2:00:00
                    #PBS -l filesystems=home_fs

					cd \${PBS_O_WORKDIR} 

					echo Jobid: \$PBS_JOBID
					echo Running on host \$HOSTNAME
					echo Running on nodes \$PBS_NODEFILE

					echo ${c} ${m} ${Nt} ${list}
					date; mpiexec -n ${np} -ppn ${ppn} --cpu-bind verbose,list${list} ./${exe_name}\_${fp}.exe -Nt ${Nt} -cpst ${c} -multiple ${m}; date

					cat \$PBS_JOBNAME.e\${PBS_JOBID%%.*} \$PBS_JOBNAME.o\${PBS_JOBID%%.*} > \
					\${PBS_JOBID%%.*}\_${prefix}\_fp_${fp}\_cpst_${c}\_multiple_${m}\_Nt_${Nt}\_np_${np}\_ppn_${ppn}\_${postfix}.out

					rm \$PBS_JOBNAME.e\${PBS_JOBID%%.*} \$PBS_JOBNAME.o\${PBS_JOBID%%.*}
					" > $job_file
					while [ `qstat -u longfei | wc -l` -eq 100 ]; do sleep 10; done;
					qsub $job_file
					# rm $job_file
				done; 
			done; 
		done; 
	done;
done;