#!/bin/sh
#PBS -j oe
#PBS -l nodes=2:ppn=16
#PBS -N MNME-DDM 
#PBS -o NXxNY_COMB_PREC/NXxNY_COMB_PREC.txt 

ulimit -m 7700000

cd LAUNCH_DIR

nodes=$(mpirun -np 32 hostname | sort | uniq)
PROCS=$(echo COMB | cut -f1 -d"x")

if [ $PROCS -eq 1 ] 
then
 num_nodes=1
 NODE1=`echo $nodes | sed s/" "" "*/#/g | cut -f1 -d#`
 cat ../template_COMB | sed s/NODE1/$NODE1/g > rankfile
else
 num_nodes=2
 NODE1=`echo $nodes | sed s/" "" "*/#/g | cut -f1 -d#`
 NODE2=`echo $nodes | sed s/" "" "*/#/g | cut -f2 -d#`
 cat ../template_COMB | sed s/NODE1/$NODE1/g | sed s/NODE2/$NODE2/g > rankfile 
fi

let PPN=$PROCS/$num_nodes
#####Generate hostfile######
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

THREADS=$(echo COMB | cut -f2 -d"x")
export OMP_NUM_THREADS=$THREADS
export MKL_NUM_THREADS=$THREADS
command="mpirun -np $PROCS -x OMP_NUM_THREADS -x MKL_NUM_THREADS --hostfile hostfile --rankfile rankfile PROGREXEC NX NY 0.0 1.0 0.0 1.0 PREC"
$command > stdout 2> stderr
# #PBS -N NXxNY_PROCS_PREC
