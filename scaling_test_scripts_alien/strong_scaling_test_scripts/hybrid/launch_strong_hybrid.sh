#!/bin/bash
DRIVER=$HOME/MNME_DDLSSC_FEM2D_CODE_SKELETON/Executables/INTEL/fem2d_poisson_driver.Om

#Preconditioner to be tested
num_precs=1
lst_precs="identity#diagonal"

#Global problem sizes to be tested (n_x * n_y mesh nodes)
num_nx_ny=1
lst_nx="50"
lst_ny="50"

#Number of MPI tasks x OpenMP threads for each global problem size 
num_combs=5
lst_combs="1x1#1x2#1x4#1x8#1x16"

num_combs=5
lst_combs="2x1#2x2#2x4#2x8#2x16"

num_combs=4
lst_combs="4x1#4x2#4x4#4x8"

num_combs=3
lst_combs="8x1#8x2#8x4"

num_combs=2
lst_combs="16x1#16x2"

id_prec=1
while [ $id_prec -le $num_precs ]
do
  prec=$(echo $lst_precs|cut -f$id_prec -d#)

  id_nx_ny=1
  while [ $id_nx_ny -le $num_nx_ny ] #iterate over global problem sizes
  do
    nx=$(echo $lst_nx|cut -f$id_nx_ny -d#)
    ny=$(echo $lst_ny|cut -f$id_nx_ny -d#)

    id_comb=1
    while [ $id_comb -le $num_combs ]
    do
      comb=$(echo $lst_combs|cut -f$id_comb -d#) 
      procs=$(echo $comb | cut -f1 -d"x")
      let mult=$ny/$procs
  
      if [ $mult -gt 1 ] #execute if at least 2 rows of grid points per process 
      then
         echo $nx $ny $comb $prec
         dir="$nx"x"$ny"_"$comb"_"$prec"
         rm -Rf $dir
         mkdir $dir
         dir_abs=$(pwd)/$dir

         cat torque_template_strong_hybrid.sh | sed "s:LAUNCH_DIR:$dir_abs:g" |sed "s:PROGREXEC:$DRIVER:g"| \
                sed "s:COMB:$comb:g"|sed "s:PREC:$prec:g"|sed "s:NX:$nx:g"|sed "s:NY:$ny:g" > $dir/"$dir"".sh"

         qsub -q iqDDM $dir/"$dir"".sh"
      fi

      let id_comb=id_comb+1
    done 
    let id_nx_ny=id_nx_ny+1
  done
  let id_prec=id_prec+1  
done

