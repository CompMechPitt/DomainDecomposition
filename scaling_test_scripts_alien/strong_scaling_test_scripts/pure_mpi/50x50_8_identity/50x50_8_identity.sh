#!/bin/sh
#PBS -j oe
#PBS -l nodes=2:ppn=16
#PBS -N MNME-DDM 
#PBS -o 50x50_8_identity/50x50_8_identity.txt 

ulimit -m 7700000

cd /home/yuyang/DD_Final_Assignment/scaling_test_scripts_alien/strong_scaling_test_scripts/pure_mpi/50x50_8_identity

nodes=$(mpirun -np 32 hostname | sort | uniq)
if [ 8 -eq 1 ] 
then
 num_nodes=1
else
 num_nodes=2
fi
let PPN=8/$num_nodes

#####Generate hostfiles######
i=1
while [ $i -le $num_nodes ]
do
  j=1
  while [ $j -le $PPN ] 
  do
    host=`echo $nodes | sed s/" "" "*/#/g | cut -f$i -d#`
    echo $host >> hostfile
    let j=j+1
  done
  let i=i+1
done
############################

#####Generate rankfile######
core_pus="0#1#2#3#4#5#6#7#8#9#10#11#12#13#14#15"
core_pus="0#4#8#12#1#5#9#13#2#6#10#14#3#7#11#15"
i=1
task=0
while [ $i -le $num_nodes ]
do
  j=1
  while [ $j -le $PPN ]
  do
    host=`echo $nodes | sed s/" "" "*/#/g | cut -f$i -d#`
    echo "rank $task=$host slot=$(echo $core_pus|cut -f$j -d#)" >> rankfile
    let task=task+1
    let j=j+1
  done
  let i=i+1
done
############################

command="mpirun -np 8 --hostfile hostfile --rankfile rankfile /home/yuyang/MNME_DDLSSC_FEM2D_CODE_SKELETON/Executables/INTEL/fem2d_poisson_driver.O 50 50 0.0 1.0 0.0 1.0 identity"
#command="mpirun -np 8 --hostfile hostfile /home/yuyang/MNME_DDLSSC_FEM2D_CODE_SKELETON/Executables/INTEL/fem2d_poisson_driver.O 50 50 0.0 1.0 0.0 1.0 identity"
#command="mpirun -np 8 /home/yuyang/MNME_DDLSSC_FEM2D_CODE_SKELETON/Executables/INTEL/fem2d_poisson_driver.O 50 50 0.0 1.0 0.0 1.0 identity"
$command > stdout 2> stderr

# #PBS -N 50x50_8_identity
