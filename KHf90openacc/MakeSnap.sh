#!/bin/bash

fstfile=`ls -1 ./meddata/Sc*.xss 2>/dev/null | head -1`
echo $fstfile
declare -i fstnum=`echo  ${fstfile##*/} | tr -cd '0123456789\n' |sed -e 's/^0\+\([0-9]\+\)$/\1/'`
echo $fstnum

lstfile=`ls -1 ./meddata/Sc*.xss 2>/dev/null | tail -1`
echo $lstfile
declare -i lstnum=`echo  ${lstfile##*/} | tr -cd '0123456789\n' |sed -e 's/^0\+\([0-9]\+\)$/\1/'`
echo $lstnum

for n in $(seq ${fstnum} ${lstnum}); do
    gnuplot -e ifnum=$n dn2dx.plt
done
