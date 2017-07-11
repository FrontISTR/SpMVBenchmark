#!/bin/bash


NCPU=`grep physical.id /proc/cpuinfo | sort -u | wc -l`
VCORE=`grep processor /proc/cpuinfo | wc -l`
RCORE=`cat /proc/cpuinfo | grep -P '^(physical|core) id\s*:' | paste - - | sort -u | wc -l`
CPUNAME=`cat /proc/cpuinfo |grep "model name"|head -n 1 | cut -d ':' -f 2`
CPUNAME=`cat /proc/cpuinfo |grep "model name"|head -n 1 | cut -d ':' -f 2`
CACHESIZE=`cat /proc/cpuinfo |grep "cache size"|head -n 1 | cut -d ':' -f 2`
NCOL=$1

echo $CPUNAME
echo $CACHESIZE
echo "NCOL: $NCOL"
echo "VectorLength: "$(( 9*$NCOL*$NCOL*$NCOL ))
echo "Vectorsize: "$(( 8*9*$NCOL*$NCOL*$NCOL/1024 ))"kB"
echo "Cores, GFLOPS"

function do_test () {
    
    export OMP_NUM_THREADS=$1
    echo -n "$1",
    time1=`./spmv33 $NCOL | grep "GFLOP"|awk '{print $2}'`
    echo $time1
}

for i in `seq 1 $RCORE`; do
    OMP=$i
    do_test $OMP
done

#for i in `seq 1 10`; do
#  OMP=$(( $i*$RCORE / 10 ))
#  if [ $OMP -le 4 ] ; then
#    continue
#  fi 
#  do_test $OMP
#done

