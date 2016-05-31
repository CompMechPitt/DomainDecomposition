#!/bin/bash
DRIVER=$HOME/MNME_DDLSSC_FEM2D_CODE_SKELETON/Executables/INTEL/fem2d_poisson_driver.O

#Preconditioner to be tested
#num_precs=3
#lst_precs="identity#diagonal#neumann"
num_precs=2
lst_precs="diagonal#neumann"

#ny local per processor (i.e., global mesh has [n_x] * [P*n_y]  nodes)
num_local_ny=1
lst_local_ny="10"
nx=50

#Number of MPI tasks to be tested for the weak scaling test 
num_procs=6
lst_procs=1#2#4#8#16#32

id_prec=1
while [ $id_prec -le $num_precs ]
do
  prec=$(echo $lst_precs|cut -f$id_prec -d#)

  id_local_ny=1
  while [ $id_local_ny -le $num_local_ny ] #iterate over global problem sizes
  do
    k=$(echo $lst_local_ny|cut -f$id_local_ny -d#)
    id_proc=1
    while [ $id_proc -le $num_procs ]
    do
      proc=$(echo $lst_procs|cut -f$id_proc -d#) 
      let ny=$proc*$k
      let mult=$ny/$proc
  
      if [ $mult -gt 1 ] #execute if at least 2 rows of grid points per process 
      then
         echo $nx $ny $proc $prec
         dir="$nx"x"$ny"_"$proc"_"$prec"
         rm -Rf $dir
         mkdir $dir
         dir_abs=$(pwd)/$dir

         cat torque_template_weak_pure_mpi.sh | sed "s:LAUNCH_DIR:$dir_abs:g" |sed "s:PROGREXEC:$DRIVER:g"| \
                sed "s:PROCS:$proc:g"|sed "s:PREC:$prec:g"|sed "s:NX:$nx:g"|sed "s:NY:$ny:g" > $dir/"$dir"".sh"

         qsub -q iqDDM $dir/"$dir"".sh"
      fi

      let id_proc=id_proc+1
    done 
    let id_local_ny=id_local_ny+1
  done
  let id_prec=id_prec+1  
done

